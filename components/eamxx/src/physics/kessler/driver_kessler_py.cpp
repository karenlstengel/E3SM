/*
 * driver_kessler_py.cpp
 *
 * Standalone driver for the JAX Kessler (kessler.py), called through
 * pybind11 the way KesslerMicrophysics::run_impl calls it, without the rest
 * of EAMxx. For finding where the time goes between EAMxx's fields and JAX,
 * and for trying device-resident alternatives.
 *
 * The "fields" are Kokkos views shaped like EAMxx's: Real (ncols, nlev_pad)
 * views, LayoutRight, with nlev padded to a multiple of SCREAM_PACK_SIZE, and
 * qv/qc/qr/qi as strided slices of one (ncols, 4, nlev_pad) tracers view.
 * They are handed to Python exactly as EAMxx does it:
 *   - PySession (share/core/eamxx_pysession.hpp, used as is) starts the
 *     interpreter and imports kessler.py with FPEs off;
 *   - make_py_array() is create_py_field() (share/field/field_pyutils.hpp)
 *     on a raw pointer: py::array(dtype, shape, byte strides, ptr,
 *     py::none()), i.e. a NumPy array viewing the memory, no copy;
 *   - the calls are kessler.init(latvap, P0, rhoqr, gravit, cpair, rair),
 *     kessler.run(...) and kessler.update(...) with run_impl's arguments.
 *
 * Backends (--backend, comma-separated):
 *   host    run_impl with py_backend=host: sync_to_host of the 14 fields
 *           (strided tracers through a contiguous helper, a fence per field,
 *           as Field::sync_to_host does), NumPy views of the host mirrors,
 *           run, update, sync_to_dev of the 12 outputs.
 *   device  run_impl with py_backend=device: kessler.set_backend("device"),
 *           NumPy views of the DEVICE pointers (create_py_field<Device>), no
 *           host syncs; qv/qc/qr go through contiguous device copies (run_impl's
 *           m_py_contig clones), packed before and unpacked after the calls.
 *           kessler.py imports the device pointers through
 *           __cuda_array_interface__ on a GPU and writes back with
 *           cudaMemcpy; on a CPU build they are host pointers and kessler.py
 *           uses its host path. (Before kessler.py had set_backend, this
 *           backend segfaulted on a GPU: kessler.py read the device pointers
 *           as host memory.)
 *   cai     GPU only. Each field goes to Python as an object with a
 *           __cuda_array_interface__ (strided tracers packed into contiguous
 *           device scratch first); driver_kessler_cai.py imports them with
 *           jnp.asarray (no host copy), runs the device-resident bridges and
 *           returns JAX arrays, whose device pointers the driver copies back
 *           into the fields on the device. No host transfers at all.
 *
 * Every backend starts from the same state (restored, untimed, before each
 * call), and every backend's outputs are compared bit for bit with the first
 * backend's. Times are per stage (seconds per call; the stages are the
 * GPTL timers run_impl has: sync_to_host, py_run, py_update, sync_to_dev).
 *
 * Environment (set before Python starts unless already set):
 *   JAX_PLATFORMS                  cuda on a GPU build, cpu otherwise
 *                                  (kessler.py itself defaults to cpu)
 *   XLA_PYTHON_CLIENT_PREALLOCATE  false, so JAX does not grab 75% of the
 *                                  GPU that Kokkos also allocates from
 *   KESSLER_PERF_LOG_PATH          ./driver_kessler_py_perf_log.csv
 *
 * Build and run on a compute node (see build_driver_kessler_py.sh):
 *   ./driver_kessler_py_gpu --grid ne16-gpu --backend host,cai --repeat 10
 *   ./driver_kessler_py_gpu --grid ne16-gpu --backend device
 *   ./driver_kessler_py_cpu --grid ne16-cpu --backend host,device
 *
 * --run-variant LIST runs each backend once per kessler::py_run_variant
 * value (kessler.set_run_variant, before warmup, as initialize_impl does).
 *
 * --passes reports the JAX kessler_run's per-column sub-cycle pass counts
 * and times instead (driver_kessler_passes.py); --variants times the
 * kessler_jax_experimental sub-cycle variants against it
 * (driver_kessler_variants.py).
 *
 * Checks: every run's outputs are compared bit for bit with the first
 * run's, at every step of a --nsteps chain; every later chain must repeat
 * the first chain's outputs bit for bit (determinism); and every output must
 * be finite, with qv/qc/qr/precl non-negative. check_kessler_py_bfb.sh runs
 * these over the grids, rain fractions, backends and variants.
 *
 * Exit code 0 only if every backend ran, matched the first one, repeated
 * itself, and passed the sanity checks.
 */

#include "share/core/eamxx_pysession.hpp"
#include "share/core/eamxx_types.hpp"
#include "share/physics/physics_constants.hpp"

#include <Kokkos_Core.hpp>
#include <pybind11/embed.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>
#include <vector>

#ifndef KESSLER_PY_DIR
#define KESSLER_PY_DIR "."
#endif

namespace py = pybind11;
using scream::Real;
using PC = scream::physics::Constants<Real>;

namespace {

using ExeSpace = Kokkos::DefaultExecutionSpace;
using DevMem   = ExeSpace::memory_space;
using View1d   = Kokkos::View<Real*,   Kokkos::LayoutRight, DevMem>;
using View2d   = Kokkos::View<Real**,  Kokkos::LayoutRight, DevMem>;
using View3d   = Kokkos::View<Real***, Kokkos::LayoutRight, DevMem>;
using Strided2 = Kokkos::View<Real**,  Kokkos::LayoutStride, DevMem>;

constexpr bool host_is_device =
  Kokkos::SpaceAccessibility<Kokkos::HostSpace, DevMem>::accessible;

// ---------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------
struct Options {
  int  ncols = 108;
  int  nz    = 58;
  Real dt    = 1800;
  Real rain_frac = 0.08;
  int  repeat = 5;
  int  warmup = 1;     // the first call compiles the JAX code
  int  nsteps = 1;     // calls chained from one state (each from the last one's output)
  bool per_call = false;
  bool py_warmup = false;
  bool passes = false;
  std::string variants;   // --variants: comma-separated, passed to Python
  std::vector<std::string> run_variants = {"current"};   // --run-variant
  std::vector<std::string> backends;
  std::string py_dir = KESSLER_PY_DIR;
};

void usage () {
  std::printf(
    "Usage: driver_kessler_py [options]\n"
    "  --grid G         ne16-cpu (108 cols), ne16-gpu (3457), ne60-cpu (1519),\n"
    "                   ne60-gpu (48601); 58 levels, dt 1800 s (one rank of the\n"
    "                   production layouts, as driver_kessler.cpp)\n"
    "  --ncols N, --nz N, --dt S\n"
    "  --rain-frac F    fraction of rainy columns (default 0.08)\n"
    "  --backend LIST   host, device, cai (default: host,device,cai on GPU,\n"
    "                   host,device on CPU). See the top of this file.\n"
    "  --repeat N       timed calls per backend (default 5)\n"
    "  --warmup N       untimed calls first (default 1; compiles JAX)\n"
    "  --nsteps N       chain N calls from one start state (each call starts from\n"
    "                   the previous one's output, as model steps do), then restart\n"
    "                   from the start state; every step is compared (default 1)\n"
    "  --per-call       print every call's stage times (warmup included)\n"
    "  --py-warmup      call kessler.warmup() after set_backend, as initialize_impl\n"
    "                   does with kessler::py_warmup=true; with --per-call, check\n"
    "                   for JIT cache misses after it\n"
    "  --passes         instead of the backends: per-column sub-cycle pass\n"
    "                   counts of the JAX kessler_run, and what restructuring its\n"
    "                   vmap(while_loop) would save (driver_kessler_passes.py)\n"
    "  --variants LIST  instead of the backends: time the kessler_jax_experimental\n"
    "                   sub-cycle variants (vmap,split,map,split_map,pallas,pallas_lanes or\n"
    "                   all) against the current kessler_run and compare outputs\n"
    "                   (driver_kessler_variants.py)\n"
    "  --run-variant LIST  kessler.set_run_variant for the host/device backends\n"
    "                   (kessler::py_run_variant): current, auto, pallas, split,\n"
    "                   map, split_map, vmap. Each backend runs once per variant;\n"
    "                   all are compared with the first (default: current)\n"
    "  --py-dir DIR     where kessler.py, bridge/, kessler_jax/ and\n"
    "                   driver_kessler_*.py are (default: the source dir)\n");
}

bool parse (int argc, char** argv, Options& o) {
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&]() -> std::string {
      if (i + 1 >= argc) { usage(); std::exit(2); }
      return argv[++i];
    };
    if      (a == "-h" || a == "--help") { usage(); std::exit(0); }
    else if (a == "--grid") {
      const std::string g = next();
      if      (g == "ne16-cpu") o.ncols = 108;
      else if (g == "ne16-gpu") o.ncols = 3457;
      else if (g == "ne60-cpu") o.ncols = 1519;
      else if (g == "ne60-gpu") o.ncols = 48601;
      else { std::printf("Unknown grid %s\n", g.c_str()); return false; }
      o.nz = 58; o.dt = 1800;
    }
    else if (a == "--ncols")     o.ncols = std::atoi(next().c_str());
    else if (a == "--nz")        o.nz = std::atoi(next().c_str());
    else if (a == "--dt")        o.dt = std::atof(next().c_str());
    else if (a == "--rain-frac") o.rain_frac = std::atof(next().c_str());
    else if (a == "--repeat")    o.repeat = std::max(1, std::atoi(next().c_str()));
    else if (a == "--warmup")    o.warmup = std::max(0, std::atoi(next().c_str()));
    else if (a == "--nsteps")    o.nsteps = std::max(1, std::atoi(next().c_str()));
    else if (a == "--py-dir")    o.py_dir = next();
    else if (a == "--per-call")  o.per_call = true;
    else if (a == "--py-warmup") o.py_warmup = true;
    else if (a == "--passes")    o.passes = true;
    else if (a == "--variants")  o.variants = next();
    else if (a == "--run-variant") {
      const std::string sv = next();
      o.run_variants.clear();
      for (size_t q = 0; q <= sv.size();) {
        const size_t e = std::min(sv.find(',', q), sv.size());
        o.run_variants.push_back(sv.substr(q, e - q));
        q = e + 1;
      }
    }
    else if (a == "--backend") {
      std::string s = next();
      o.backends.clear();
      for (size_t p = 0; p <= s.size();) {
        const size_t e = std::min(s.find(',', p), s.size());
        const std::string b = s.substr(p, e - p);
        if (b != "host" && b != "device" && b != "cai") {
          std::printf("Unknown backend '%s'\n", b.c_str()); return false;
        }
        if (b == "cai" && host_is_device) {
          std::printf("--backend cai needs a GPU build\n"); return false;
        }
        o.backends.push_back(b);
        p = e + 1;
      }
    }
    else if (a.rfind("--kokkos", 0) == 0) { /* handled by Kokkos */ }
    else { std::printf("Unknown option %s\n", a.c_str()); usage(); return false; }
  }
  if (o.backends.empty()) {
    o.backends = host_is_device ? std::vector<std::string>{"host", "device"}
                                : std::vector<std::string>{"host", "device", "cai"};
  }
  return true;
}

// ---------------------------------------------------------------------------
// Synthetic columns: the same profile as driver_kessler.cpp's build_profile
// (stretched 34 km grid, lapse-rate T, 85% RH below 10 km, cloud and rain in
// a fraction rain_frac of the columns with per-column amounts).
// ---------------------------------------------------------------------------
struct Profile { std::vector<Real> rho, z_mid, pk, theta, qv, qc, qr; };

Profile build_profile (const Options& o) {
  const int ncols = o.ncols, nz = o.nz;
  const Real ztop = 34000.0, p0 = PC::P0.value, rd = PC::Rair.value,
             cp = PC::Cpair.value, H = 7500.0;
  Profile p;
  for (auto* v : {&p.rho, &p.z_mid, &p.pk, &p.theta, &p.qv, &p.qc, &p.qr})
    v->assign(size_t(ncols) * nz, 0.0);
  std::vector<Real> z_int(nz + 1);
  for (int k = 0; k <= nz; ++k) z_int[k] = ztop * std::pow(Real(nz - k) / nz, 1.6);
  for (int c = 0; c < ncols; ++c) {
    const bool rainy = std::floor((c + 1) * o.rain_frac) > std::floor(c * o.rain_frac);
    const Real pert = std::sin(6.283185307179586 * c / std::max(ncols, 1));
    const Real amp  = 0.25 * std::pow(8.0, std::fmod(c * 0.6180339887498949, 1.0));
    for (int k = 0; k < nz; ++k) {
      const size_t i = size_t(c) * nz + k;
      const Real z  = 0.5 * (z_int[k] + z_int[k + 1]);
      const Real T  = std::max(300.0 + 2.0 * pert - 0.0065 * z, 210.0);
      const Real pr = p0 * std::exp(-z / H);
      const Real pk = std::pow(pr / p0, rd / cp);
      const Real qvs = 3.8 / (pr / 100.0) * std::exp(17.27 * (T - 273.0) / (T - 36.0));
      const Real rh  = (z < 10000.0) ? std::min(0.85 + 0.1 * pert, 0.99) : 1.0e-6;
      p.z_mid[i] = z;
      p.rho[i]   = pr / (rd * T);
      p.pk[i]    = pk;
      p.theta[i] = T / pk;
      p.qv[i]    = std::max(rh * qvs, 1.0e-12);
      if (rainy && z > 2000.0 && z < 6000.0) p.qc[i] = amp * 1.0e-3;
      if (rainy && z > 1000.0 && z < 4000.0) p.qr[i] = amp * 5.0e-4;
    }
  }
  return p;
}

// ---------------------------------------------------------------------------
// The fields, laid out like EAMxx's
// ---------------------------------------------------------------------------
// run_impl's fields, in its order. Tracers are components of `tracers`.
const char* const field_names[] = {
  "qv", "qc", "qr", "rho", "z_mid", "pk", "theta", "precl", "relhum",
  "phis", "T_mid", "T_mid_prev", "T_mid_tend", "st_energy"};
// run_impl's sync_to_dev list (everything but rho and z_mid)
const char* const output_names[] = {
  "qv", "qr", "qc", "pk", "theta", "precl", "relhum",
  "phis", "T_mid", "T_mid_prev", "T_mid_tend", "st_energy"};
const std::map<std::string, int> tracer_index = {{"qv", 0}, {"qc", 1}, {"qr", 2}};

struct Fields {
  int ncols, nz, npad;
  View3d tracers;                       // (ncols, 4, npad): qv, qc, qr, qi
  std::map<std::string, View2d> f2;     // (ncols, npad)
  std::map<std::string, View1d> f1;     // (ncols)
  View3d::HostMirror tracers_h;
  std::map<std::string, View2d::HostMirror> f2_h;
  std::map<std::string, View1d::HostMirror> f1_h;

  Fields (int nc, int nzz) : ncols(nc), nz(nzz),
    npad(((nzz + SCREAM_PACK_SIZE - 1) / SCREAM_PACK_SIZE) * SCREAM_PACK_SIZE)
  {
    tracers   = View3d("tracers", ncols, 4, npad);
    tracers_h = Kokkos::create_mirror_view(tracers);
    for (const char* n : field_names) {
      const std::string s = n;
      if (tracer_index.count(s)) continue;
      if (s == "precl" || s == "phis") {
        f1[s] = View1d(s, ncols);  f1_h[s] = Kokkos::create_mirror_view(f1[s]);
      } else {
        f2[s] = View2d(s, ncols, npad);  f2_h[s] = Kokkos::create_mirror_view(f2[s]);
      }
    }
  }
  bool is_1d (const std::string& n) const { return f1.count(n) > 0; }

  // The (ncols, nz) device/host view of a 2D field, as EAMxx's get_view
  // gives it: a strided slice of the padded allocation.
  Strided2 dev2 (const std::string& n) const {
    const auto lev = std::make_pair(0, nz);
    if (tracer_index.count(n)) return Kokkos::subview(tracers, Kokkos::ALL, tracer_index.at(n), lev);
    return Kokkos::subview(f2.at(n), Kokkos::ALL, lev);
  }
  Kokkos::View<Real**, Kokkos::LayoutStride, Kokkos::HostSpace> host2 (const std::string& n) const {
    const auto lev = std::make_pair(0, nz);
    if (tracer_index.count(n)) return Kokkos::subview(tracers_h, Kokkos::ALL, tracer_index.at(n), lev);
    return Kokkos::subview(f2_h.at(n), Kokkos::ALL, lev);
  }
};

void fill (Fields& F, const Profile& p) {
  Kokkos::deep_copy(F.tracers_h, 0);
  for (auto& kv : F.f2_h) Kokkos::deep_copy(kv.second, 0);
  for (auto& kv : F.f1_h) Kokkos::deep_copy(kv.second, 0);
  for (int c = 0; c < F.ncols; ++c) {
    for (int k = 0; k < F.npad; ++k) {
      const size_t i = size_t(c) * F.nz + std::min(k, F.nz - 1);  // padding repeats the last level
      F.tracers_h(c, 0, k) = p.qv[i];
      F.tracers_h(c, 1, k) = p.qc[i];
      F.tracers_h(c, 2, k) = p.qr[i];
      F.f2_h["rho"](c, k)   = p.rho[i];
      F.f2_h["z_mid"](c, k) = p.z_mid[i];
      F.f2_h["pk"](c, k)    = p.pk[i];
      F.f2_h["theta"](c, k) = p.theta[i];
      F.f2_h["T_mid"](c, k) = p.pk[i] * p.theta[i];
    }
  }
  Kokkos::deep_copy(F.tracers, F.tracers_h);
  for (auto& kv : F.f2) Kokkos::deep_copy(kv.second, F.f2_h[kv.first]);
  for (auto& kv : F.f1) Kokkos::deep_copy(kv.second, F.f1_h[kv.first]);
}

// All outputs on the host, (ncols * nz) per 2D field, for comparisons
std::map<std::string, std::vector<Real>> snapshot (const Fields& F) {
  std::map<std::string, std::vector<Real>> s;
  auto tr = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), F.tracers);
  for (const char* n : output_names) {
    const std::string name = n;
    auto& v = s[name];
    if (F.is_1d(name)) {
      auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), F.f1.at(name));
      v.assign(h.data(), h.data() + F.ncols);
    } else if (tracer_index.count(name)) {
      const int it = tracer_index.at(name);
      for (int c = 0; c < F.ncols; ++c)
        for (int k = 0; k < F.nz; ++k) v.push_back(tr(c, it, k));
    } else {
      auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), F.f2.at(name));
      for (int c = 0; c < F.ncols; ++c)
        for (int k = 0; k < F.nz; ++k) v.push_back(h(c, k));
    }
  }
  return s;
}

// ---------------------------------------------------------------------------
// Python objects for the fields
// ---------------------------------------------------------------------------
// create_py_field(): a NumPy array viewing (not copying) the memory
py::array make_py_array (const Real* ptr, std::vector<ssize_t> shape, std::vector<ssize_t> strides_elems) {
  std::vector<ssize_t> strides;
  for (auto s : strides_elems) strides.push_back(s * ssize_t(sizeof(Real)));
  return py::array(py::dtype::of<Real>(), shape, strides, ptr, py::none());
}

// The NumPy views of every field, over the host mirrors (host) or the device
// allocations (device), as m_py_fields_host / m_py_fields_dev hold them.
std::map<std::string, py::array> numpy_views (const Fields& F, const bool on_device) {
  std::map<std::string, py::array> m;
  for (const char* n : field_names) {
    const std::string s = n;
    if (F.is_1d(s)) {
      const Real* p = on_device ? F.f1.at(s).data() : F.f1_h.at(s).data();
      m[s] = make_py_array(p, {F.ncols}, {1});
    } else if (on_device) {
      const auto v = F.dev2(s);
      m[s] = make_py_array(v.data(), {F.ncols, F.nz}, {ssize_t(v.stride(0)), ssize_t(v.stride(1))});
    } else {
      const auto v = F.host2(s);
      m[s] = make_py_array(v.data(), {F.ncols, F.nz}, {ssize_t(v.stride(0)), ssize_t(v.stride(1))});
    }
  }
  return m;
}

// An object with a __cuda_array_interface__ (v3) for a C-contiguous device array
py::object cai_object (const Real* ptr, const std::vector<ssize_t>& shape) {
  py::dict cai;
  py::tuple shp(shape.size());
  for (size_t i = 0; i < shape.size(); ++i) shp[i] = py::int_(shape[i]);
  cai["shape"]   = shp;
  cai["typestr"] = std::string(sizeof(Real) == 8 ? "<f8" : "<f4");
  cai["data"]    = py::make_tuple(reinterpret_cast<std::uintptr_t>(ptr), false);
  cai["strides"] = py::none();   // C-contiguous
  cai["version"] = 3;
  // No "stream": the driver fences Kokkos before every call, so the data is ready.
  auto ns = py::module_::import("types").attr("SimpleNamespace")();
  py::setattr(ns, "__cuda_array_interface__", cai);
  return ns;
}

// Copy a JAX output (via its __cuda_array_interface__) into a device field
template <typename DstView>
void copy_from_cai (const py::handle& jax_array, const DstView& dst) {
  py::dict cai = jax_array.attr("__cuda_array_interface__");
  auto data = cai["data"].cast<py::tuple>();
  const Real* p = reinterpret_cast<const Real*>(data[0].cast<std::uintptr_t>());
  EKAT_REQUIRE_MSG(!cai.contains("strides") || cai["strides"].is_none(),
                   "Error! JAX output is not C-contiguous.\n");
  if constexpr (static_cast<int>(DstView::rank) == 1) {
    Kokkos::View<const Real*, DevMem, Kokkos::MemoryUnmanaged> src(p, dst.extent(0));
    Kokkos::deep_copy(dst, src);
  } else {
    Kokkos::View<const Real**, Kokkos::LayoutRight, DevMem, Kokkos::MemoryUnmanaged>
      src(p, dst.extent(0), dst.extent(1));
    Kokkos::deep_copy(dst, src);
  }
}

// ---------------------------------------------------------------------------
// One call of run_impl's Python section
// ---------------------------------------------------------------------------
using Times = std::map<std::string, double>;

struct Backend {
  std::string name;
  std::map<std::string, py::array> np;     // host/device: NumPy views (made once, like EAMxx)
  std::map<std::string, View2d> packed;    // cai: contiguous scratch for qv/qc/qr
  std::map<std::string, View2d> helper_d;  // host: Field's contiguous helpers for
  std::map<std::string, View2d::HostMirror> helper_h;  //   the strided tracers
};

// Field::sync_to_host / sync_to_dev for every listed field: contiguous fields
// copy their whole allocation; strided tracers go through a contiguous
// helper on each side. One fence per field (fence=true default).
void sync (const Fields& F, const std::vector<std::string>& names, const bool to_host,
           std::map<std::string, View2d>& helper_d, std::map<std::string, View2d::HostMirror>& helper_h) {
  if (host_is_device) return;   // Field::sync_to_host's early return
  const ExeSpace ex;
  for (const auto& n : names) {
    if (F.is_1d(n)) {
      if (to_host) Kokkos::deep_copy(ex, F.f1_h.at(n), F.f1.at(n));
      else         Kokkos::deep_copy(ex, F.f1.at(n), F.f1_h.at(n));
    } else if (tracer_index.count(n)) {
      if (!helper_d.count(n)) {
        helper_d[n] = View2d("helper_" + n, F.ncols, F.nz);
        helper_h[n] = Kokkos::create_mirror_view(helper_d[n]);
      }
      if (to_host) {
        Kokkos::deep_copy(ex, helper_d[n], F.dev2(n));
        Kokkos::deep_copy(ex, helper_h[n], helper_d[n]);
        ex.fence();
        Kokkos::deep_copy(F.host2(n), helper_h[n]);
      } else {
        Kokkos::deep_copy(helper_h[n], F.host2(n));
        Kokkos::deep_copy(ex, helper_d[n], helper_h[n]);
        Kokkos::deep_copy(ex, F.dev2(n), helper_d[n]);
      }
    } else {
      if (to_host) Kokkos::deep_copy(ex, F.f2_h.at(n), F.f2.at(n));
      else         Kokkos::deep_copy(ex, F.f2.at(n), F.f2_h.at(n));
    }
    Kokkos::fence();
  }
}

// driver_kessler_cai.last_times, as stage:part entries
void add_py_times (Times& t, const std::string& stage, const py::module_& m) {
  for (auto kv : m.attr("last_times").cast<py::dict>())
    t[stage + ":" + kv.first.cast<std::string>()] = kv.second.cast<double>();
}

Times call_once (Backend& b, Fields& F, const Options& o,
                 const py::module_& kessler, const py::module_& cai_mod) {
  Times t;
  const int ncols = F.ncols, nz = F.nz;
  const int lyr_surf = nz, lyr_toa = 1;   // as run_impl passes them (1-based)
  const double dt = o.dt;
  Kokkos::fence();

  if (b.name == "host" || b.name == "device") {
    auto& a = b.np;
    if (b.name == "device") {
      Kokkos::Timer tm;
      for (auto& kv : b.packed) Kokkos::deep_copy(kv.second, F.dev2(kv.first));
      Kokkos::fence();
      t["py_pack"] = tm.seconds();
    }
    if (b.name == "host") {
      Kokkos::Timer tm;
      std::vector<std::string> all(std::begin(field_names), std::end(field_names));
      sync(F, all, true, b.helper_d, b.helper_h);
      t["sync_to_host"] = tm.seconds();
    }
    Kokkos::Timer tm;
    kessler.attr("run")(ncols, nz, dt, lyr_surf, lyr_toa, a["rho"], a["z_mid"], a["pk"],
                        a["theta"], a["qv"], a["qc"], a["qr"], a["precl"], a["relhum"]);
    t["py_run"] = tm.seconds();
    tm.reset();
    kessler.attr("update")(ncols, nz, dt, a["z_mid"], a["pk"], a["theta"], a["phis"],
                           a["T_mid"], a["T_mid_prev"], a["T_mid_tend"], a["st_energy"]);
    t["py_update"] = tm.seconds();
    if (b.name == "host") {
      tm.reset();
      std::vector<std::string> outs(std::begin(output_names), std::end(output_names));
      sync(F, outs, false, b.helper_d, b.helper_h);
      t["sync_to_dev"] = tm.seconds();
    }
    if (b.name == "device") {
      tm.reset();
      for (auto& kv : b.packed) Kokkos::deep_copy(F.dev2(kv.first), kv.second);
      Kokkos::fence();
      t["py_unpack"] = tm.seconds();
    }
    return t;
  }

  // cai: pack the strided tracers into contiguous device scratch, hand every
  // field over as a __cuda_array_interface__ object, copy the JAX outputs
  // back into the fields on the device.
  Kokkos::Timer tm;
  for (const auto& kv : tracer_index) {
    if (!b.packed.count(kv.first)) b.packed[kv.first] = View2d("packed_" + kv.first, ncols, nz);
    Kokkos::deep_copy(b.packed[kv.first], F.dev2(kv.first));
  }
  Kokkos::fence();
  auto dev = [&](const std::string& n) -> py::object {
    if (F.is_1d(n)) return cai_object(F.f1.at(n).data(), {ncols});
    if (b.packed.count(n)) return cai_object(b.packed.at(n).data(), {ncols, nz});
    EKAT_REQUIRE_MSG(F.npad == nz, "Error! cai needs unpadded fields (Pack size 1).\n");
    return cai_object(F.f2.at(n).data(), {ncols, nz});
  };
  t["pack"] = tm.seconds();

  tm.reset();
  auto r = cai_mod.attr("run")(ncols, nz, dt, lyr_surf, lyr_toa, dev("rho"), dev("z_mid"),
                                    dev("pk"), dev("theta"), dev("qv"), dev("qc"), dev("qr"),
                                    dev("precl"), dev("relhum")).cast<py::tuple>();
  t["py_run"] = tm.seconds();
  add_py_times(t, "py_run", cai_mod);
  tm.reset();
  copy_from_cai(r[0], F.dev2("theta"));
  copy_from_cai(r[1], F.dev2("qv"));
  copy_from_cai(r[2], F.dev2("qc"));
  copy_from_cai(r[3], F.dev2("qr"));
  copy_from_cai(r[4], F.f1.at("precl"));
  copy_from_cai(r[5], F.dev2("relhum"));
  Kokkos::fence();
  t["copy_back_run"] = tm.seconds();

  tm.reset();
  auto u = cai_mod.attr("update")(ncols, nz, dt, dev("z_mid"), dev("pk"), dev("theta"),
                                       dev("phis"), dev("T_mid"), dev("T_mid_prev"),
                                       dev("T_mid_tend"), dev("st_energy")).cast<py::tuple>();
  t["py_update"] = tm.seconds();
  add_py_times(t, "py_update", cai_mod);
  tm.reset();
  copy_from_cai(u[0], F.dev2("T_mid"));
  copy_from_cai(u[1], F.dev2("T_mid_prev"));
  copy_from_cai(u[2], F.dev2("T_mid_tend"));
  copy_from_cai(u[3], F.dev2("st_energy"));
  Kokkos::fence();
  t["copy_back_update"] = tm.seconds();
  return t;
}

using Snapshot = std::map<std::string, std::vector<Real>>;

// Number of values that differ between two snapshots (bit for bit)
size_t ndiff (const Snapshot& a, const Snapshot& b) {
  size_t nd = 0;
  for (const auto& kv : a) {
    const auto& x = kv.second;
    const auto& y = b.at(kv.first);
    for (size_t i = 0; i < x.size(); ++i) nd += (x[i] != y[i]);
  }
  return nd;
}

// Physical sanity of one snapshot, as driver_kessler.cpp checks it: every
// output finite, and qv/qc/qr/precl non-negative. Returns the problems found.
std::vector<std::string> sanity (const Snapshot& s) {
  std::vector<std::string> bad;
  for (const auto& kv : s) {
    size_t nonfinite = 0, negative = 0;
    const bool nonneg = kv.first == "qv" || kv.first == "qc" || kv.first == "qr" || kv.first == "precl";
    for (const Real x : kv.second) {
      nonfinite += !std::isfinite(x);
      negative += nonneg && x < 0;
    }
    if (nonfinite) bad.push_back(kv.first + ": " + std::to_string(nonfinite) + " non-finite");
    if (negative)  bad.push_back(kv.first + ": " + std::to_string(negative) + " negative");
  }
  return bad;
}

bool compare (const std::string& a_name, const std::map<std::string, std::vector<Real>>& a,
              const std::string& b_name, const std::map<std::string, std::vector<Real>>& b,
              const std::string& what = "") {
  bool ok = true;
  std::printf("\n=== %s vs %s%s (bit for bit)\n", b_name.c_str(), a_name.c_str(), what.c_str());
  for (const auto& kv : a) {
    const auto& x = kv.second;
    const auto& y = b.at(kv.first);
    size_t nd = 0; Real mr = 0;
    for (size_t i = 0; i < x.size(); ++i) {
      if (x[i] == y[i]) continue;
      ++nd;
      const Real r = std::abs(x[i] - y[i]) / std::max(std::abs(x[i]), Real(1e-300));
      if (!(r <= mr)) mr = r;
    }
    if (nd) std::printf("  %-11s %zu of %zu differ (max rel %.3e)\n", kv.first.c_str(), nd, x.size(), double(mr));
    ok = ok && nd == 0;
  }
  std::printf("  %s\n", ok ? "MATCH" : "DIFFER");
  return ok;
}

} // anonymous namespace

int main (int argc, char** argv)
{
  Options o;
  if (!parse(argc, argv, o)) return 2;

  // Before Python (and so JAX) starts
  setenv("JAX_PLATFORMS", host_is_device ? "cpu" : "cuda", 0);
  setenv("XLA_PYTHON_CLIENT_PREALLOCATE", "false", 0);
  setenv("KESSLER_PERF_LOG_PATH", "./driver_kessler_py_perf_log.csv", 0);

  Kokkos::initialize(argc, argv);
  int rc = 0;
  {
    auto& session = scream::PySession::get();
    session.initialize();
    session.add_path(o.py_dir);
    {
      std::printf("driver_kessler_py: ncols=%d nz=%d (padded %d) dt=%g rain_frac=%g, Kokkos %s, "
                  "JAX_PLATFORMS=%s, py dir %s\n",
                  o.ncols, o.nz, ((o.nz + SCREAM_PACK_SIZE - 1) / SCREAM_PACK_SIZE) * SCREAM_PACK_SIZE, double(o.dt), double(o.rain_frac),
                  ExeSpace::name(), std::getenv("JAX_PLATFORMS"), o.py_dir.c_str());
      std::fflush(stdout);

      py::module_ kessler = session.safe_import("kessler");
      py::module_ cai_mod;
      if (std::find(o.backends.begin(), o.backends.end(), "cai") != o.backends.end())
        cai_mod = session.safe_import("driver_kessler_cai");
      // initialize_impl's call
      kessler.attr("init")(PC::LatVap.value, PC::P0.value, PC::RHOW.value,
                           PC::gravit.value, PC::Cpair.value, PC::Rair.value);
      std::printf("jax devices: %s\n",
                  py::str(py::module_::import("jax").attr("devices")()).cast<std::string>().c_str());

      // Count JAX traces and lowerings (each one is an in-memory jit cache
      // miss; an XLA compile or persistent-cache read follows it), by
      // function, so calls after a warmup can be checked for misses.
      py::dict jax_events;
      py::exec(R"(
import jax
events = []
def _listener(event, duration, **kw):
    if event.endswith(("jaxpr_trace_duration", "jaxpr_to_mlir_module_duration")):
        events.append(event.split("/")[-1].replace("_duration", "") + "[" + str(kw.get("fun_name", "?")) + "]")
jax.monitoring.register_event_duration_secs_listener(_listener)
)", jax_events);
      py::list jax_ev = jax_events["events"];
      auto take_events = [&]() {
        std::vector<std::string> v;
        for (auto e : jax_ev) v.push_back(e.cast<std::string>());
        jax_ev.attr("clear")();
        return v;
      };

      const Profile prof = build_profile(o);
      Fields F(o.ncols, o.nz);

      if (o.passes) {
        // The input state as kessler.run gets it with py_backend=host
        // (fill() writes the host mirrors first)
        fill(F, prof);
        auto a = numpy_views(F, false);
        py::module_ pm = session.safe_import("driver_kessler_passes");
        const bool same = pm.attr("report")(
          o.ncols, o.nz, double(o.dt), o.nz, 1, a["rho"], a["z_mid"], a["pk"], a["theta"],
          a["qv"], a["qc"], a["qr"], a["precl"], a["relhum"], o.repeat).cast<bool>();
        if (!same) rc = 1;
        o.backends.clear();
      }
      if (!o.variants.empty()) {
        fill(F, prof);
        auto a = numpy_views(F, false);
        py::list vl;
        if (o.variants == "all") {
          for (const char* v : {"vmap", "split", "map", "split_map", "pallas", "pallas_lanes"}) vl.append(v);
        } else {
          for (size_t q = 0; q <= o.variants.size();) {
            const size_t e = std::min(o.variants.find(',', q), o.variants.size());
            vl.append(o.variants.substr(q, e - q));
            q = e + 1;
          }
        }
        py::module_ vm = session.safe_import("driver_kessler_variants");
        const bool ran = vm.attr("compare")(
          o.ncols, o.nz, double(o.dt), o.nz, 1, a["rho"], a["z_mid"], a["pk"], a["theta"],
          a["qv"], a["qc"], a["qr"], a["precl"], a["relhum"], py::tuple(vl), o.repeat).cast<bool>();
        if (!ran) rc = 1;
        o.backends.clear();
      }

      std::vector<std::vector<Snapshot>> results;   // per run, per step of the first chain
      std::vector<std::map<std::string, std::vector<double>>> times;
      // (backend, run variant) pairs; cai bypasses kessler.run, so once only
      std::vector<std::pair<std::string, std::string>> runs;
      for (const auto& v : o.run_variants)
        for (const auto& name : o.backends)
          if (name != "cai" || v == o.run_variants[0]) runs.emplace_back(name, v);
      const bool tag = o.run_variants.size() > 1 || o.run_variants[0] != "current";
      std::vector<std::string> labels;
      for (const auto& nv : runs) {
        const std::string& name = nv.first;
        const std::string label = (tag && name != "cai") ? name + "/" + nv.second : name;
        labels.push_back(label);
        Backend b;
        b.name = name;
        if (name != "cai") b.np = numpy_views(F, name == "device");
        if (name == "device") {
          // run_impl's m_py_contig: contiguous copies of the strided tracers
          for (const auto& kv : tracer_index) {
            b.packed[kv.first] = View2d("packed_" + kv.first, F.ncols, F.nz);
            b.np[kv.first] = make_py_array(b.packed[kv.first].data(), {F.ncols, F.nz}, {F.nz, 1});
          }
        }
        // initialize_impl's set_backend call (cai bypasses kessler.run/update)
        kessler.attr("set_backend")(name == "device" ? "device" : "host");
        kessler.attr("set_run_variant")(nv.second);
        if (o.py_warmup && name != "cai") {
          // initialize_impl, after set_backend and set_run_variant:
          // warmup(m_num_cols, m_num_levs, lyr_surf, lyr_toa)
          take_events();
          Kokkos::Timer tw;
          kessler.attr("warmup")(o.ncols, o.nz, o.nz, 1);
          const auto ev = take_events();
          std::printf("  [%s] kessler.warmup: %.2f s, %zu traces/lowerings\n", label.c_str(), tw.seconds(), ev.size());
        }
        std::map<std::string, std::vector<double>> tb;
        std::vector<Snapshot> steps;   // outputs of each step of the first chain
        size_t timed_misses = 0, repeat_diffs = 0, repeat_checked = 0;
        take_events();
        for (int r = 0; r < o.warmup + o.repeat; ++r) {
          const int step = r % o.nsteps;
          if (step == 0) fill(F, prof);   // start state (untimed); later steps chain
          Kokkos::Timer total;
          Times t = call_once(b, F, o, kessler, cai_mod);
          const double tt = total.seconds();
          // Every call's outputs (untimed): the first chain's are kept, later
          // chains must repeat them bit for bit (catches races, e.g. a missing
          // barrier in a kernel, which would not differ on every call)
          if (r < o.nsteps) {
            steps.push_back(snapshot(F));
          } else {
            const size_t nd = ndiff(steps[step], snapshot(F));
            ++repeat_checked;
            if (nd) {
              ++repeat_diffs;
              std::printf("  [%s] call %d (step %d): %zu values differ from the first chain's step %d\n",
                          label.c_str(), r, step, nd, step);
            }
          }
          const auto ev = take_events();
          if (o.per_call) {
            std::printf("  [%s] call %d%s: total %.4e", label.c_str(), r, r < o.warmup ? " (warmup)" : "", tt);
            for (const auto& kv : t) std::printf(", %s %.4e", kv.first.c_str(), kv.second);
            std::printf(", jit misses %zu\n", ev.size());
            for (const auto& e : ev) std::printf("      %s\n", e.c_str());
          }
          if (r < o.warmup) continue;
          timed_misses += ev.size();
          for (const auto& kv : t) tb[kv.first].push_back(kv.second);
          tb["total"].push_back(tt);
        }
        std::printf("  [%s] JAX traces/lowerings during timed calls: %zu\n", label.c_str(), timed_misses);
        std::printf("  [%s] repeated calls identical to the first chain: %zu of %zu%s\n", label.c_str(),
                    repeat_checked - repeat_diffs, repeat_checked, repeat_diffs ? "  REPEAT DIFFERS" : "");
        if (repeat_diffs) rc = 1;
        if (results.empty())   // the first run: show that the chained steps evolve the state
          for (size_t k = 1; k < steps.size(); ++k)
            std::printf("  [%s] step %zu changed %zu values from step %zu\n", label.c_str(), k,
                        ndiff(steps[k - 1], steps[k]), k - 1);
        for (size_t k = 0; k < steps.size(); ++k) {
          const auto bad = sanity(steps[k]);
          for (const auto& m : bad) std::printf("  [%s] step %zu: SANITY %s\n", label.c_str(), k, m.c_str());
          if (!bad.empty()) rc = 1;
        }
        results.push_back(steps);
        times.push_back(tb);
      }

      // Every run vs the first, at every step of the chain
      size_t ncmp = 0, nbad = 0;
      for (size_t i = 1; i < results.size(); ++i)
        for (size_t k = 0; k < results[i].size(); ++k) {
          const std::string what = o.nsteps > 1 ? ", step " + std::to_string(k) : "";
          ++ncmp;
          if (!compare(labels[0], results[0][k], labels[i], results[i][k], what)) { rc = 1; ++nbad; }
        }
      if (ncmp)
        std::printf("\nbit-for-bit vs %s: %zu of %zu run/step comparisons match\n",
                    labels[0].c_str(), ncmp - nbad, ncmp);

      std::printf("\n--- seconds per call, %d timed call(s) after %d warmup (min / mean)\n",
                  o.repeat, o.warmup);
      for (size_t i = 0; i < times.size(); ++i) {
        std::printf("  %s\n", labels[i].c_str());
        for (const auto& kv : times[i]) {
          const auto& v = kv.second;
          double mean = 0; for (double x : v) mean += x; mean /= v.size();
          std::printf("    %-18s %.4e  %.4e\n", kv.first.c_str(),
                      *std::min_element(v.begin(), v.end()), mean);
        }
      }
      kessler.attr("finalize")();
    }
    session.finalize();
  }
  Kokkos::finalize();
  std::printf("\ndriver_kessler_py: %s\n", rc == 0 ? "ALL MATCH" : "DIFFERENCES (see above)");
  return rc;
}
