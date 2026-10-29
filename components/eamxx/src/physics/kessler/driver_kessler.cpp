/*
 * Standalone test driver for KesslerMicrophysicsFunctions::kessler_run.
 *
 * Built to chase the CPP_ne16np4_1day_gpu failure
 *   "Error! KESSLER: sub-cycling cannot converge (iter = 10, ...)"
 * where the per-column sub-cycle bookkeeping (dt0 / mask / time_counter)
 * held values kessler_run never writes (e.g. mask=6966.7), so no column
 * ever converged.
 *
 * The driver runs kessler_run outside of EAMxx on a synthetic set of
 * top-down columns shaped like the ne16np4 case (58 levels, ~3433
 * columns per rank, dt = 1800 s) and checks the result. It can build the
 * Workspace two ways:
 *
 *   managed : one Kokkos view per scratch array
 *   carved  : one contiguous allocation carved into Unmanaged views in
 *             exactly the order KesslerMicrophysics::init_buffers uses
 *
 * Before each call the Workspace can be "poisoned" with a garbage value,
 * so any read of scratch kessler_run did not initialise itself shows up
 * in the result. The device run is also compared against the same code
 * instantiated on the host.
 *
 * The device run is also repeated with qv/qc/qr passed the way EAMxx
 * passes them: as slices of one bundled (ncol, 4, npack) tracers array
 * ordered [qv, qr, qi, qc] (ekat::subview_1, what Field::get_view returns
 * for a tracer subfield), so each has column stride 4*npack instead of
 * npack. A kernel that ekat::scalarize()s one of these views -- which
 * drops the stride -- reads other tracers/columns. That was the
 * CPP_ne30np4_10day precl bug: the state was right, but precl came from
 * the wrong memory. Use --rain for this check to mean anything: with no
 * rain, surface qr is 0 everywhere and a misread precl is still 0.
 *
 * Checks applied to every run:
 *   - outputs finite; qv/qc/qr/precl non-negative
 *   - sub-cycling converged: mask == 0, time_counter == dt in every column
 *   - latent-heat balance: kessler_run only moves theta together with qv
 *     (theta += lv/(cp*pk)*x, qv -= x), so cp*pk*d(theta) + lv*d(qv) == 0
 *     at every point to round-off
 *   - relhum matches 100*qv/qvs recomputed on the host from the outputs
 *   - no rain (default, like the moist baroclinic wave IC at t=0): the
 *     subsaturated, cloud-free profile must come back unchanged (theta,
 *     qv bit-for-bit; qc = qr = precl = 0) -- exactly the case where the
 *     GPU build silently skipped the dt0/mask updates
 *   - --rain: precl > 0 in every rainy column and == 0 in every dry one
 * Comparisons:
 *   - device vs host reference (rtol 1e-8, atol 1e-14: GPU pow/exp and
 *     FMA contraction differ from the host in the last bits, and ~60 rain
 *     sub-cycles grow that to ~1e-9 relative on near-zero qc/qr)
 *   - managed vs carved workspace (bit-for-bit)
 *   - poisoned with --poison-val vs poisoned with NaN (bit-for-bit)
 *   - strided tracers vs contiguous qv/qc/qr (bit-for-bit)
 * Exit code 0 only if everything passes.
 *
 * Performance (see performance_suggestions.md, sections C and D):
 *   --impl takes a comma-separated list of implementations: prod (the
 *   production kessler_run, one team per column, item D1), hostloop (the
 *   previous host-driven sub-cycle loop, driver_kessler_hostloop.hpp), and
 *   the experimental variants in driver_kessler_variants.hpp: base (the
 *   variant code without I2 or I6, i.e. production before I2 was ported
 *   on 2026-10-28), i2 (skip pow on all-dry packs), i6 (per-level
 *   scratch in team scratch) and i2i6 (both); "both" is prod,hostloop
 *   and "all" is every one. Every run is done for each, and each
 *   non-prod implementation is compared to prod bit for bit.
 *   --repeat/--warmup time repeated calls on the same input
 *   (the state is restored, untimed, before each repeat), --kernel-timers
 *   adds a per-kernel breakdown (C3), and the hostloop runs print the
 *   per-column sub-cycle counts, i.e. how much work the host loop
 *   wasted (C1).
 *   --rain-frac sets the fraction of rainy columns (C2); --grid sets
 *   ncols/nz/dt to one rank of a production case. Typical use:
 *     OMP_NUM_THREADS=1 ./driver_kessler_cpu --perf --grid ne16-cpu --rain-frac 0.08
 *     ./driver_kessler_gpu --perf --grid ne16-gpu --rain-frac 0.08
 *
 * How to read the result:
 *   - device fails, host passes          -> bug in kessler_run's device path
 *   - carved fails, managed passes       -> Workspace carving/aliasing bug
 *   - poison comparison differs          -> kessler_run reads uninitialised scratch
 *   - strided comparison differs         -> kessler_run ignores a tracer view's
 *                                           stride (e.g. scalarizes qv/qc/qr)
 *   - everything passes here, but the
 *     model still fails                  -> the scratch is being clobbered or
 *                                           mis-set outside kessler_run
 *                                           (ATMBufferManager sharing, init order)
 *
 * Build: see build_driver_kessler.sh in this directory.
 * Usage: ./driver_kessler [--help] [options]
 */

#include "kessler_functions.hpp"
#include "kessler_functions_impl.hpp"  // needed for non-GPU / RDC builds
#include "driver_kessler_hostloop.hpp"
#include "driver_kessler_variants.hpp"
#include "driver_kessler_preprocess.hpp"

#include <Kokkos_Core.hpp>
#include <ekat_subview_utils.hpp>
#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace {

using scream::Real;
using PC = scream::physics::Constants<Real>;

using DevDevice  = scream::DefaultDevice;
using HostDevice = Kokkos::Device<Kokkos::DefaultHostExecutionSpace, Kokkos::HostSpace>;

using KMF_D = scream::kessler::KesslerMicrophysicsFunctions<Real, DevDevice>;
using KMF_H = scream::kessler::KesslerMicrophysicsFunctions<Real, HostDevice>;

// ---------------------------------------------------------------------------
// Options
// ---------------------------------------------------------------------------
// Implementations --impl can select (see the comment at the top).
enum class Impl { prod, hostloop, base, i2, i6, i2i6 };
const Impl all_impls[] = {Impl::prod, Impl::hostloop, Impl::base, Impl::i2, Impl::i6, Impl::i2i6};
const char* impl_name (const Impl i) {
  switch (i) {
    case Impl::prod:     return "prod";
    case Impl::hostloop: return "hostloop";
    case Impl::base:     return "base";
    case Impl::i2:       return "i2";
    case Impl::i6:       return "i6";
    case Impl::i2i6:     return "i2i6";
  }
  return "?";
}

struct Options {
  int  ncols      = 3433;   // ~ne16np4 columns per rank with 4 ranks
  int  nz         = 58;
  Real dt         = 1800.0; // ATM_NCPL=48
  int  nsteps     = 1;      // repeated kessler_run calls (like physics steps)
  bool managed    = true;
  bool carved     = true;
  bool host_ref   = true;
  bool poison     = true;
  Real poison_val = 6966.7; // value seen in mask in the failed GPU run
  Real rain_frac  = 0;      // fraction of columns with cloud/rain (sub-cycles)
  bool strided    = true;   // also run with qv/qc/qr as strided tracer slices
  std::vector<Impl> impls = {Impl::prod};  // --impl
  int  repeat     = 1;      // timed repeats of the nsteps calls
  int  warmup     = 0;      // untimed repeats before those
  bool kernel_timers = false;
  bool preprocess   = false;  // --preprocess: I4 test instead of kessler_run

  bool rain () const { return rain_frac > 0; }
};

void usage () {
  std::printf(
    "Usage: driver_kessler [options]\n"
    "  --ncols N         number of columns (default 3433)\n"
    "  --nz N            number of levels (default 58)\n"
    "  --dt S            physics time step in seconds (default 1800)\n"
    "  --grid G          ncols/nz/dt of one rank of a production case:\n"
    "                      ne16-cpu: 108 cols (13826 / 128 ranks), 58 levs, 1800 s\n"
    "                      ne16-gpu: 3457 cols (13826 / 4 ranks), 58 levs, 1800 s\n"
    "                      ne60-cpu: 1519 cols (194402 / 128 ranks), 58 levs, 1800 s\n"
    "                      ne60-gpu: 48601 cols (194402 / 4 ranks), 58 levs, 1800 s\n"
    "  --nsteps N        repeated kessler_run calls (default 1)\n"
    "  --ws MODE         workspace: managed | carved | both (default both)\n"
    "  --no-host-ref     skip the host reference run\n"
    "  --host-ref        do the host reference run (undoes --perf's --no-host-ref)\n"
    "  --no-poison       do not pre-fill the workspace with garbage (also skips\n"
    "                    the NaN-poison comparison run)\n"
    "  --poison-val X    garbage value for the workspace (default 6966.7)\n"
    "  --rain            same as --rain-frac 0.5\n"
    "  --rain-frac F     put cloud water and rain in a fraction F of the columns\n"
    "                    (spread evenly, each with its own rain amount) so\n"
    "                    sedimentation sub-cycles; 0.03-0.12 matches days 10-15\n"
    "                    of the dcmip2016 runs. kessler_run aborts after 100000\n"
    "                    sub-cycles.\n"
    "  --no-strided      skip the run with qv/qc/qr as strided slices of a\n"
    "                    tracers array (EAMxx's layout; combine with --rain,\n"
    "                    or a precl misread is invisible)\n"
    "  --impl LIST       comma-separated implementations (default prod):\n"
    "                    prod      production kessler_run\n"
    "                    hostloop  previous host-loop version\n"
    "                              (driver_kessler_hostloop.hpp)\n"
    "                    base      variant code without I2 or I6\n"
    "                              (= production before I2 was ported)\n"
    "                    i2        skip pow on all-dry packs\n"
    "                    i6        per-level scratch in team scratch\n"
    "                    i2i6      both (driver_kessler_variants.hpp)\n"
    "                    both = prod,hostloop; all = every one. Each non-prod\n"
    "                    implementation is compared to prod bit for bit.\n"
    "  --repeat N        time N repeats of the nsteps calls, restoring the\n"
    "                    input state before each (default 1)\n"
    "  --warmup N        untimed repeats before the timed ones (default 0)\n"
    "  --kernel-timers   per-kernel time breakdown of the timed repeats. Fences\n"
    "                    after every kernel, so on GPU take the total time from\n"
    "                    a run without it.\n"
    "  --preprocess      instead of kessler_run, compare and time run_impl's\n"
    "                    pre-processing, old vs new (I4, driver_kessler_preprocess.hpp);\n"
    "                    uses --grid/--ncols/--nz, --repeat, --warmup, --host-ref\n"
    "  --perf            timing setup: --ws managed --no-poison --no-strided\n"
    "                    --no-host-ref --impl both --repeat 5 --warmup 1\n"
    "                    (options after it override it)\n");
}

bool parse (int argc, char** argv, Options& opt) {
  for (int i = 1; i < argc; ++i) {
    const std::string a = argv[i];
    auto next = [&]() -> const char* {
      if (i + 1 >= argc) { std::printf("Missing value for %s\n", a.c_str()); std::exit(2); }
      return argv[++i];
    };
    if      (a == "--help" || a == "-h") { usage(); return false; }
    else if (a == "--ncols")       opt.ncols  = std::atoi(next());
    else if (a == "--nz")          opt.nz     = std::atoi(next());
    else if (a == "--dt")          opt.dt     = std::atof(next());
    else if (a == "--nsteps")      opt.nsteps = std::atoi(next());
    else if (a == "--ws") {
      const std::string m = next();
      opt.managed = (m == "managed" || m == "both");
      opt.carved  = (m == "carved"  || m == "both");
      if (!opt.managed && !opt.carved) { usage(); std::exit(2); }
    }
    else if (a == "--grid") {
      const std::string g = next();
      if      (g == "ne16-cpu") { opt.ncols = 108;  opt.nz = 58; opt.dt = 1800.0; }
      else if (g == "ne16-gpu") { opt.ncols = 3457; opt.nz = 58; opt.dt = 1800.0; }
      else if (g == "ne60-cpu") { opt.ncols = 1519; opt.nz = 58; opt.dt = 1800.0; }
      else if (g == "ne60-gpu") { opt.ncols = 48601; opt.nz = 58; opt.dt = 1800.0; }
      else { std::printf("Unknown grid: %s\n", g.c_str()); usage(); std::exit(2); }
    }
    else if (a == "--no-host-ref") opt.host_ref   = false;
    else if (a == "--host-ref")    opt.host_ref   = true;
    else if (a == "--no-poison")   opt.poison     = false;
    else if (a == "--poison-val")  opt.poison_val = std::atof(next());
    else if (a == "--rain")        opt.rain_frac  = 0.5;
    else if (a == "--rain-frac")   opt.rain_frac  = std::atof(next());
    else if (a == "--no-strided")  opt.strided    = false;
    else if (a == "--impl") {
      std::string m = next();
      if (m == "both") m = "prod,hostloop";
      if (m == "all")  m = "prod,hostloop,base,i2,i6,i2i6";
      opt.impls.clear();
      size_t pos = 0;
      while (pos <= m.size()) {
        const size_t end = std::min(m.find(',', pos), m.size());
        const std::string name = m.substr(pos, end - pos);
        bool found = false;
        for (const Impl im : all_impls) {
          if (name == impl_name(im)) {
            found = true;
            if (std::find(opt.impls.begin(), opt.impls.end(), im) == opt.impls.end())
              opt.impls.push_back(im);
          }
        }
        if (!found) { std::printf("Unknown --impl entry: '%s'\n", name.c_str()); usage(); std::exit(2); }
        pos = end + 1;
      }
    }
    else if (a == "--repeat")      opt.repeat     = std::max(1, std::atoi(next()));
    else if (a == "--warmup")      opt.warmup     = std::max(0, std::atoi(next()));
    else if (a == "--kernel-timers") opt.kernel_timers = true;
    else if (a == "--preprocess")  opt.preprocess = true;
    else if (a == "--perf") {
      opt.managed = true;  opt.carved = false;
      opt.poison  = false; opt.strided = false; opt.host_ref = false;
      opt.impls   = {Impl::prod, Impl::hostloop};
      opt.repeat  = 5;     opt.warmup = 1;
    }
    else if (a.rfind("--kokkos", 0) == 0) { /* handled by Kokkos */ }
    else { std::printf("Unknown option: %s\n", a.c_str()); usage(); std::exit(2); }
  }
  return true;
}

// ---------------------------------------------------------------------------
// Synthetic input profile (host, scalar, [col*nz + k], k=0 is TOA)
// ---------------------------------------------------------------------------
struct Profile {
  std::vector<Real> rho, z_mid, pk, theta, qv, qc, qr;
  std::vector<char> rainy;  // per column
};

Profile build_profile (const Options& opt) {
  const int ncols = opt.ncols, nz = opt.nz;
  const Real ztop = 34000.0;         // ~top of the 58-level grid in the failed run
  const Real p0   = PC::P0.value;    // exner reference pressure (Pa)
  const Real rd   = PC::Rair.value;
  const Real cp   = PC::Cpair.value;
  const Real H    = 7500.0;          // scale height for a simple p(z)

  Profile p;
  for (auto* v : {&p.rho, &p.z_mid, &p.pk, &p.theta, &p.qv, &p.qc, &p.qr}) {
    v->assign(size_t(ncols) * nz, 0.0);
  }

  std::vector<Real> z_int(nz + 1);
  for (int k = 0; k <= nz; ++k) {
    // Stretched grid: thin layers near the surface, thick aloft.
    const Real frac = Real(nz - k) / nz;
    z_int[k] = ztop * std::pow(frac, 1.6);
  }

  // A fraction rain_frac of the columns is rainy, spread evenly (column c
  // is rainy when floor((c+1)*F) > floor(c*F)). Each rainy column gets its
  // own rain/cloud amount, log-uniform in [0.25, 2] x the base amount, so
  // the per-column sub-cycle counts vary the way they do in the model.
  p.rainy.assign(ncols, 0);
  for (int c = 0; c < ncols; ++c) {
    p.rainy[c] = std::floor((c + 1) * opt.rain_frac) > std::floor(c * opt.rain_frac);
  }

  for (int c = 0; c < ncols; ++c) {
    const Real two_pi = 6.283185307179586;
    const Real pert = std::sin(two_pi * c / std::max(ncols, 1));
    const Real u    = std::fmod(c * 0.6180339887498949, 1.0);
    const Real amp  = 0.25 * std::pow(8.0, u);
    for (int k = 0; k < nz; ++k) {
      const size_t i = size_t(c) * nz + k;
      const Real z   = 0.5 * (z_int[k] + z_int[k + 1]);
      const Real T   = std::max(300.0 + 2.0 * pert - 0.0065 * z, 210.0);
      const Real pr  = p0 * std::exp(-z / H);
      const Real pk  = std::pow(pr / p0, rd / cp);

      // Kessler's own saturation formula: 3.8/p[hPa] * exp(...)
      const Real qvs = 3.8 / (pr / 100.0) * std::exp(17.27 * (T - 273.0) / (T - 36.0));
      const Real rh  = (z < 10000.0) ? std::min(0.85 + 0.1 * pert, 0.99) : 1.0e-6;

      p.z_mid[i] = z;
      p.rho[i]   = pr / (rd * T);        // kg/m^3, as the process interface computes
      p.pk[i]    = pk;
      p.theta[i] = T / pk;
      p.qv[i]    = std::max(rh * qvs, 1.0e-12);

      if (p.rainy[c]) {
        if (z > 2000.0 && z < 6000.0) p.qc[i] = amp * 1.0e-3;
        if (z > 1000.0 && z < 4000.0) p.qr[i] = amp * 5.0e-4;
      }
    }
  }
  return p;
}

// ---------------------------------------------------------------------------
// View helpers
// ---------------------------------------------------------------------------
template <typename KMF>
using PView2d = typename KMF::template view_2d<typename KMF::Pack>;
template <typename KMF>
using RView1d = typename KMF::template view_1d<Real>;

// Copy a [col*nz + k] host array into a packed (ncols, npack) view.
// Padding lanes (npack*Pack::n > nz) repeat the last valid level so the
// math on them stays benign.
template <typename KMF>
PView2d<KMF> make_2d (const std::string& name, const std::vector<Real>& src,
                      const int ncols, const int nz) {
  using Pack = typename KMF::Pack;
  const int npack = ekat::npack<Pack>(nz);
  PView2d<KMF> v(name, ncols, npack);
  auto h  = Kokkos::create_mirror_view(v);
  auto hs = ekat::scalarize(h);
  for (int c = 0; c < ncols; ++c) {
    for (int k = 0; k < npack * Pack::n; ++k) {
      hs(c, k) = src[size_t(c) * nz + std::min(k, nz - 1)];
    }
  }
  Kokkos::deep_copy(v, h);
  return v;
}

template <typename KMF>
std::vector<Real> to_vec (const PView2d<KMF>& v, const int nz) {
  auto h  = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
  auto hs = ekat::scalarize(h);
  const int ncols = h.extent(0);
  std::vector<Real> out(size_t(ncols) * nz);
  for (int c = 0; c < ncols; ++c)
    for (int k = 0; k < nz; ++k)
      out[size_t(c) * nz + k] = hs(c, k);
  return out;
}

template <typename KMF>
std::vector<Real> to_vec (const RView1d<KMF>& v) {
  auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
  return std::vector<Real>(h.data(), h.data() + h.extent(0));
}

// ---------------------------------------------------------------------------
// Workspace construction
// ---------------------------------------------------------------------------
template <typename KMF>
struct WorkspaceHolder {
  typename KMF::Workspace ws;
  RView1d<KMF> pool;  // owns the memory in carved mode
};

template <typename KMF>
WorkspaceHolder<KMF> make_workspace (const bool carved, const int ncols, const int nz) {
  using Pack = typename KMF::Pack;
  using Dev  = typename KMF::Device;
  const int npack = ekat::npack<Pack>(nz);

  WorkspaceHolder<KMF> h;
  auto& ws = h.ws;

  RView1d<KMF>* v1[] = {&ws.dt0, &ws.mask, &ws.time_counter, &ws.precl_acc};
  PView2d<KMF>* v2[] = {&ws.r, &ws.rhalf, &ws.velqr, &ws.sed, &ws.pc, &ws.f5};
  static_assert(KMF::Workspace::num_1d_scalar == 4, "update v1[] to match Workspace");
  static_assert(KMF::Workspace::num_2d_vector == 6, "update v2[] to match Workspace");

  if (!carved) {
    const char* n1[] = {"dt0", "mask", "time_counter", "precl_acc"};
    const char* n2[] = {"r", "rhalf", "velqr", "sed", "pc", "f5"};
    for (int i = 0; i < 4; ++i) *v1[i] = RView1d<KMF>(std::string(n1[i]), ncols);
    for (int i = 0; i < 6; ++i) *v2[i] = PView2d<KMF>(std::string(n2[i]), ncols, npack);
    return h;
  }

  // Same layout as KesslerMicrophysics::init_buffers: 4 x ncols Reals,
  // then 6 x (ncols, npack) Packs, back to back in one allocation.
  const size_t n_real = size_t(4) * ncols + size_t(6) * ncols * npack * Pack::n;
  h.pool = RView1d<KMF>("kessler_ws_pool", n_real);

  using u1d = Kokkos::View<Real*,  Kokkos::LayoutRight, Dev, Kokkos::MemoryUnmanaged>;
  using u2d = Kokkos::View<Pack**, Kokkos::LayoutRight, Dev, Kokkos::MemoryUnmanaged>;

  Real* mem = h.pool.data();
  for (int i = 0; i < 4; ++i) { *v1[i] = u1d(mem, ncols); mem += ncols; }
  Pack* pmem = reinterpret_cast<Pack*>(mem);
  for (int i = 0; i < 6; ++i) { *v2[i] = u2d(pmem, ncols, npack); pmem += size_t(ncols) * npack; }

  const size_t used = reinterpret_cast<Real*>(pmem) - h.pool.data();
  if (used != n_real) {
    std::printf("driver_kessler: carved workspace used %zu Reals, expected %zu\n", used, n_real);
    std::exit(3);
  }
  return h;
}

template <typename KMF>
void poison_workspace (const typename KMF::Workspace& ws, const Real val) {
  using Pack = typename KMF::Pack;
  // Explicit pointer arrays: nvcc's front end mis-parses a braced list in
  // a range-for inside a template.
  const RView1d<KMF>* v1[] = {&ws.dt0, &ws.mask, &ws.time_counter, &ws.precl_acc};
  const PView2d<KMF>* v2[] = {&ws.r, &ws.rhalf, &ws.velqr, &ws.sed, &ws.pc, &ws.f5};
  for (const auto* v : v1) Kokkos::deep_copy(*v, val);
  for (const auto* v : v2) Kokkos::deep_copy(*v, Pack(val));
  Kokkos::fence();
}

// ---------------------------------------------------------------------------
// Per-kernel timers, through Kokkos' tool callbacks (no Kokkos Tools library
// needed). Kernels are launched one at a time from the host, so one open
// slot is enough. The end callback fences, so asynchronous (GPU) kernels are
// timed to completion.
// ---------------------------------------------------------------------------
struct KernelTimes {
  bool on = false;
  std::string open_name;
  Kokkos::Timer open_timer;
  std::map<std::string, std::pair<double, long>> acc;  // name -> (seconds, count)
};
KernelTimes& ktimes () { static KernelTimes k; return k; }

void kt_begin (const char* name, const uint32_t, uint64_t* kid) {
  auto& k = ktimes();
  *kid = 0;
  k.open_name = name;
  k.open_timer.reset();
}
void kt_end (const uint64_t) {
  Kokkos::fence();
  auto& k = ktimes();
  if (!k.on) return;
  auto& e = k.acc[k.open_name];
  e.first += k.open_timer.seconds();
  ++e.second;
}
void install_kernel_timers () {
  Kokkos::Tools::Experimental::set_begin_parallel_for_callback(kt_begin);
  Kokkos::Tools::Experimental::set_end_parallel_for_callback(kt_end);
  Kokkos::Tools::Experimental::set_begin_parallel_reduce_callback(kt_begin);
  Kokkos::Tools::Experimental::set_end_parallel_reduce_callback(kt_end);
}

// ---------------------------------------------------------------------------
// One kessler_run experiment
// ---------------------------------------------------------------------------
struct RunResult {
  std::string label;
  bool threw = false;
  std::string err;
  std::vector<Real> theta, qv, qc, qr, relhum, precl, mask, time_counter;
  std::vector<double> times;  // seconds per kessler_run call, per timed repeat
  std::map<std::string, std::pair<double, long>> kernels;  // --kernel-timers
};

// Per-column sub-cycle counts (sorted) from the host loop's active-column
// count after each pass: n_active[i] columns needed more than i+1 passes.
std::vector<int> passes_from_active (const int ncols, const std::vector<int>& n_active) {
  std::vector<int> npass;
  int prev = ncols;
  for (size_t i = 0; i < n_active.size(); ++i) {
    npass.insert(npass.end(), prev - n_active[i], int(i) + 1);
    prev = n_active[i];
  }
  return npass;
}

// Sub-cycle passes per column for one call (C1). The host loop swept every
// column max(npass) times; prod does sum(npass) column-passes.
void report_passes (const int step, const std::vector<int>& npass, const bool histogram) {
  const int ncols = npass.size();
  long sum = 0; int mx = 0;
  for (int n : npass) { sum += n; mx = std::max(mx, n); }
  std::printf("  sub-cycles, step %d: max %d, mean %.2f; column-passes "
              "hostloop %ld, prod %ld, hostloop/prod %.2f\n",
              step, mx, double(sum) / std::max(ncols, 1), long(mx) * ncols, sum,
              double(mx) * ncols / std::max(sum, 1L));
  if (!histogram) return;
  // Buckets 1, 2-4, 5-16, 17-64, ... (powers of 4).
  std::vector<int> cnt;
  for (int n : npass) {
    int b = 0;
    for (long hi = 1; n > hi; hi *= 4) ++b;
    if (b >= int(cnt.size())) cnt.resize(b + 1, 0);
    ++cnt[b];
  }
  std::printf("  columns by sub-cycle count:");
  long lo = 1, hi = 1;
  for (size_t b = 0; b < cnt.size(); ++b) {
    if (cnt[b]) std::printf("  [%ld-%ld] %d", lo, hi, cnt[b]);
    lo = hi + 1; hi *= 4;
  }
  std::printf("\n");
}

template <typename KMF>
RunResult run_case (const Impl impl, const std::string& run_label, const Options& opt,
                    const Profile& prof, const bool carved,
                    const Real poison_val, const bool strided = false) {
  const int ncols = opt.ncols, nz = opt.nz;

  RunResult res;
  res.label = std::string(impl_name(impl)) + ":" + run_label;
  const std::string& label = res.label;

  auto rho    = make_2d<KMF>("rho",    prof.rho,   ncols, nz);
  auto z_mid  = make_2d<KMF>("z_mid",  prof.z_mid, ncols, nz);
  auto pk     = make_2d<KMF>("pk",     prof.pk,    ncols, nz);
  auto theta  = make_2d<KMF>("theta",  prof.theta, ncols, nz);
  auto qv     = make_2d<KMF>("qv",     prof.qv,    ncols, nz);
  auto qc     = make_2d<KMF>("qc",     prof.qc,    ncols, nz);
  auto qr     = make_2d<KMF>("qr",     prof.qr,    ncols, nz);
  auto relhum = make_2d<KMF>("relhum", std::vector<Real>(size_t(ncols) * nz, 0.0), ncols, nz);
  RView1d<KMF> precl("precl", ncols);

  // Strided mode: qv/qc/qr become slices of one (ncols, 4, npack) tracers
  // array ordered [qv, qr, qi, qc], like EAMxx's tracer bundle, so each has
  // column stride 4*npack. qi stays zero. Results are copied back into the
  // contiguous views afterwards (to_vec scalarizes, which needs contiguity).
  using Pack = typename KMF::Pack;
  using TracersView = Kokkos::View<Pack***, Kokkos::LayoutRight, typename KMF::Device>;
  TracersView tracers;
  PView2d<KMF> qv_in = qv, qc_in = qc, qr_in = qr;
  if (strided) {
    tracers = TracersView("tracers", ncols, 4, ekat::npack<Pack>(nz));
    Kokkos::deep_copy(tracers, Pack(0));
    Kokkos::deep_copy(Kokkos::subview(tracers, Kokkos::ALL, 0, Kokkos::ALL), qv);
    Kokkos::deep_copy(Kokkos::subview(tracers, Kokkos::ALL, 1, Kokkos::ALL), qr);
    Kokkos::deep_copy(Kokkos::subview(tracers, Kokkos::ALL, 3, Kokkos::ALL), qc);
    qv_in = ekat::subview_1(tracers, 0);
    qr_in = ekat::subview_1(tracers, 1);
    qc_in = ekat::subview_1(tracers, 3);
  }

  // Initial state, restored (untimed) before every repeat.
  auto theta0 = make_2d<KMF>("theta0", prof.theta, ncols, nz);
  auto qv0    = make_2d<KMF>("qv0",    prof.qv,    ncols, nz);
  auto qc0    = make_2d<KMF>("qc0",    prof.qc,    ncols, nz);
  auto qr0    = make_2d<KMF>("qr0",    prof.qr,    ncols, nz);
  auto restore = [&]() {
    Kokkos::deep_copy(theta, theta0);
    Kokkos::deep_copy(qv_in, qv0);
    Kokkos::deep_copy(qc_in, qc0);
    Kokkos::deep_copy(qr_in, qr0);
  };

  std::vector<int> n_active;  // hostloop: active columns after each pass

  auto wsh = make_workspace<KMF>(carved, ncols, nz);

  // Same arguments the process interface passes (see
  // eamxx_kessler_process_interface.cpp): top-down levels, pref in hPa.
  const int  lyr_surf = nz - 1;
  const int  lyr_toa  = 0;
  const Real rhoqr    = PC::RHOW.value;
  const Real latvap   = PC::LatVap.value;
  const Real pref     = PC::P0.value / 100.0;
  const Real cpair    = PC::Cpair.value;
  const Real rair     = PC::Rair.value;

  std::printf("\n=== %s: %d step(s) x %d repeat(s) (+%d warmup), ncols=%d nz=%d dt=%g, "
              "workspace=%s%s%s\n",
              label.c_str(), opt.nsteps, opt.repeat, opt.warmup, ncols, nz, double(opt.dt),
              carved ? "carved" : "managed",
              opt.poison ? (", poisoned with " + std::to_string(poison_val)).c_str() : "",
              strided ? (", q* = tracer slices, column stride " +
                         std::to_string(qr_in.stride(0)) + " (contiguous: " +
                         std::to_string(qr.stride(0)) + ")").c_str() : "");
  auto& kt = ktimes();
  kt.acc.clear();
  try {
    for (int rep = 0; rep < opt.warmup + opt.repeat; ++rep) {
      const bool timed = rep >= opt.warmup;
      restore();
      if (opt.poison) poison_workspace<KMF>(wsh.ws, poison_val);
      Kokkos::fence();
      double t = 0;
      for (int s = 0; s < opt.nsteps; ++s) {
        kt.on = timed && opt.kernel_timers;
        Kokkos::Timer timer;
        // Same arguments for every implementation; hostloop also reports
        // its per-pass active-column counts.
#define KESSLER_ARGS ncols, nz, opt.dt, lyr_surf, lyr_toa, rhoqr, latvap, pref, \
                     cpair, rair, rho, z_mid, pk, theta, qv_in, qc_in, qr_in, \
                     precl, relhum, wsh.ws
        switch (impl) {
          case Impl::prod:     KMF::kessler_run(KESSLER_ARGS); break;
          case Impl::hostloop: kessler_driver::kessler_run_hostloop<KMF>(KESSLER_ARGS, &n_active); break;
          case Impl::base:     kessler_driver::kessler_run_variant<KMF, false, false>(KESSLER_ARGS); break;
          case Impl::i2:       kessler_driver::kessler_run_variant<KMF, true,  false>(KESSLER_ARGS); break;
          case Impl::i6:       kessler_driver::kessler_run_variant<KMF, false, true >(KESSLER_ARGS); break;
          case Impl::i2i6:     kessler_driver::kessler_run_variant<KMF, true,  true >(KESSLER_ARGS); break;
        }
#undef KESSLER_ARGS
        Kokkos::fence();
        t += timer.seconds();
        kt.on = false;
        if (impl == Impl::hostloop && rep == 0 && opt.rain()) {
          report_passes(s, passes_from_active(ncols, n_active), s == 0);
        }
      }
      if (timed) res.times.push_back(t / opt.nsteps);
    }
  } catch (const std::exception& e) {
    kt.on = false;
    res.threw = true;
    res.err   = e.what();
  }
  res.kernels = kt.acc;

  if (strided) {
    Kokkos::deep_copy(qv, Kokkos::subview(tracers, Kokkos::ALL, 0, Kokkos::ALL));
    Kokkos::deep_copy(qr, Kokkos::subview(tracers, Kokkos::ALL, 1, Kokkos::ALL));
    Kokkos::deep_copy(qc, Kokkos::subview(tracers, Kokkos::ALL, 3, Kokkos::ALL));
  }

  res.theta        = to_vec<KMF>(theta, nz);
  res.qv           = to_vec<KMF>(qv, nz);
  res.qc           = to_vec<KMF>(qc, nz);
  res.qr           = to_vec<KMF>(qr, nz);
  res.relhum       = to_vec<KMF>(relhum, nz);
  res.precl        = to_vec<KMF>(precl);
  res.mask         = to_vec<KMF>(wsh.ws.mask);
  res.time_counter = to_vec<KMF>(wsh.ws.time_counter);
  return res;
}

// ---------------------------------------------------------------------------
// Checks and comparisons
// ---------------------------------------------------------------------------
// Saturation mixing ratio exactly as kessler_run computes it (Teten's
// formula, Klemp & Wilhelmson 1978 Eq 2.11), for the host-side relhum check.
Real kessler_qvs (const Real pk, const Real theta) {
  const Real xk   = PC::Cpair.value / PC::Rair.value;
  const Real pref = PC::P0.value / 100.0;
  const Real pc   = 3.8 / (std::pow(pk, xk) * pref);
  const Real T    = pk * theta;
  return pc * std::exp(17.27 * (T - 273.0) / (T - 36.0));
}

bool check_result (const RunResult& r, const Options& opt, const Profile& prof) {
  std::printf("\n--- checks: %s\n", r.label.c_str());
  if (r.threw) {
    std::printf("  FAIL: kessler_run threw:\n%s\n", r.err.c_str());
    return false;
  }

  bool ok = true;
  auto report = [&ok](const char* name, const bool pass, const std::string& detail) {
    std::printf("  [%s] %-18s %s\n", pass ? "ok  " : "FAIL", name, detail.c_str());
    ok = ok && pass;
  };
  auto fmt = [](const char* f, auto... args) {
    char buf[256]; std::snprintf(buf, sizeof(buf), f, args...); return std::string(buf);
  };
  auto count = [](const std::vector<Real>& v, auto pred) {
    int n = 0; for (Real x : v) if (pred(x)) ++n; return n;
  };
  auto nonfinite = [](Real x) { return !std::isfinite(x); };
  auto negative  = [](Real x) { return x < 0; };

  const int ncols = opt.ncols, nz = opt.nz;
  const Real cp = PC::Cpair.value, lv = PC::LatVap.value;

  // 1. Finite / non-negative outputs.
  struct F { const char* name; const std::vector<Real>* v; bool nonneg; };
  const F fields[] = {
    {"theta", &r.theta, false}, {"qv", &r.qv, true}, {"qc", &r.qc, true},
    {"qr", &r.qr, true}, {"relhum", &r.relhum, false}, {"precl", &r.precl, true},
  };
  for (const auto& f : fields) {
    const int nf = count(*f.v, nonfinite);
    const int ng = f.nonneg ? count(*f.v, negative) : 0;
    Real sum = 0; for (Real x : *f.v) sum += x;
    report(f.name, nf == 0 && ng == 0,
           fmt("sum=%-24.16e nonfinite=%d negative=%d", double(sum), nf, ng));
  }

  // 2. Sub-cycling converged in every column: mask == 0 and
  //    time_counter == dt (the bookkeeping that failed on GPU).
  {
    int bad_mask = 0, bad_time = 0;
    Real worst_mask = 0, worst_time = 0;
    for (size_t c = 0; c < r.mask.size(); ++c) {
      if (r.mask[c] != Real(0)) { ++bad_mask; worst_mask = r.mask[c]; }
      const Real terr = std::abs(r.time_counter[c] - opt.dt);
      if (!(terr <= Real(1e-5))) { ++bad_time; worst_time = std::max(worst_time, terr); }
    }
    report("converged", bad_mask == 0 && bad_time == 0,
           fmt("mask!=0 in %d cols (e.g. %g); |time_counter-dt|>1e-5 in %d cols (max %g)",
               bad_mask, double(worst_mask), bad_time, double(worst_time)));
  }

  // 3. Latent-heat balance: cp*pk*d(theta) + lv*d(qv) == 0 at every point.
  //    Normalised by cp*T (~3e5 J/kg); only the qv >= 0 clip can break it.
  {
    Real worst = 0; size_t iw = 0;
    for (size_t i = 0; i < r.theta.size(); ++i) {
      const Real res = cp * prof.pk[i] * (r.theta[i] - prof.theta[i])
                     + lv * (r.qv[i] - prof.qv[i]);
      const Real rel = std::abs(res) / (cp * prof.pk[i] * prof.theta[i]);
      if (!(rel <= worst)) { worst = rel; iw = i; }
    }
    report("latent heat", worst <= Real(1e-12),
           fmt("max |cp*pk*dtheta + lv*dqv|/(cp*T) = %.3e (col %zu, lev %zu)",
               double(worst), iw / nz, iw % nz));
  }

  // 4. relhum == 100*qv/qvs(pk, theta) from the returned state.
  {
    Real worst = 0;
    for (size_t i = 0; i < r.relhum.size(); ++i) {
      const Real expect = 100.0 * r.qv[i] / kessler_qvs(prof.pk[i], r.theta[i]);
      const Real rel = std::abs(r.relhum[i] - expect) / std::max(std::abs(expect), Real(1e-300));
      if (!(rel <= worst)) worst = rel;
    }
    report("relhum", worst <= Real(1e-10), fmt("max rel err vs host formula = %.3e", double(worst)));
  }

  if (!opt.rain()) {
    // 5a. Known answer: cloud-free, subsaturated columns are a no-op.
    int n_theta = 0, n_qv = 0, n_cond = 0, n_precl = 0;
    for (size_t i = 0; i < r.theta.size(); ++i) {
      if (r.theta[i] != prof.theta[i]) ++n_theta;
      if (r.qv[i]    != prof.qv[i])    ++n_qv;
      if (r.qc[i] != 0 || r.qr[i] != 0) ++n_cond;
    }
    for (Real p : r.precl) if (p != 0) ++n_precl;
    report("no-rain no-op", n_theta == 0 && n_qv == 0 && n_cond == 0 && n_precl == 0,
           fmt("changed: theta %d, qv %d pts; qc/qr != 0 at %d pts; precl != 0 in %d cols",
               n_theta, n_qv, n_cond, n_precl));
  } else {
    // 5b. Rain reaches the ground in rainy columns only. precl is the
    //     last step's rate, and later steps may have rained out, so the
    //     "every rainy column precipitates" half only holds for nsteps=1.
    int wet_dry = 0, dry_wet = 0;
    Real min_wet = std::numeric_limits<Real>::max();
    for (int c = 0; c < ncols; ++c) {
      const bool rainy = prof.rainy[c];
      if (rainy) {
        min_wet = std::min(min_wet, r.precl[c]);
        if (opt.nsteps == 1 && !(r.precl[c] > 0)) ++wet_dry;
      }
      else if (r.precl[c] != 0) ++dry_wet;
    }
    report("precl pattern", wet_dry == 0 && dry_wet == 0,
           fmt("rainy cols with precl<=0: %d (min precl %.3e); dry cols with precl!=0: %d",
               wet_dry, double(min_wet), dry_wet));
  }

  std::printf("  %s\n", ok ? "PASS" : "FAIL");
  return ok;
}

// Per-field comparison of two runs. A point matches when
//   |x - y| <= atol + rtol * max(|x|, |y|)
// (numpy.isclose-style), so near-zero mixing ratios are judged by atol
// rather than by a relative error that blows up as the values vanish.
// rtol = atol = 0 demands bit-for-bit equality. A NaN on either side never
// matches (every comparison with NaN is false).
bool compare (const RunResult& a, const RunResult& b, const Real rtol, const Real atol) {
  std::printf("\n--- compare: %s vs %s (rtol %g, atol %g)\n", a.label.c_str(), b.label.c_str(),
              double(rtol), double(atol));
  if (a.threw || b.threw) {
    std::printf("  skipped: a run threw\n");
    return false;
  }
  bool ok = true;
  struct P { const char* name; const std::vector<Real>* x; const std::vector<Real>* y; };
  const P pairs[] = {
    {"theta", &a.theta, &b.theta}, {"qv", &a.qv, &b.qv}, {"qc", &a.qc, &b.qc},
    {"qr", &a.qr, &b.qr}, {"relhum", &a.relhum, &b.relhum}, {"precl", &a.precl, &b.precl},
  };
  for (const auto& p : pairs) {
    Real max_abs = 0, max_rel = 0;
    size_t n_bad = 0;
    for (size_t i = 0; i < p.x->size(); ++i) {
      const Real x = (*p.x)[i], y = (*p.y)[i];
      if (x == y) continue;  // also covers +-0 and identical infinities
      // Written so a NaN on either side propagates into max_abs/max_rel
      // (std::max(a, NaN) would silently return a).
      const Real d = std::abs(x - y);
      const Real s = std::max(std::abs(x), std::abs(y));
      if (!(d <= max_abs)) max_abs = d;
      const Real rel = (s > 0) ? d / s : d;
      if (!(rel <= max_rel)) max_rel = rel;
      if (!(d <= atol + rtol * s)) ++n_bad;
    }
    const bool f_ok = (n_bad == 0);
    std::printf("  %-7s max_abs=%-12.4e max_rel=%-12.4e mismatched=%zu %s\n", p.name,
                double(max_abs), double(max_rel), n_bad, f_ok ? "" : "<-- differs");
    ok = ok && f_ok;
  }
  std::printf("  %s\n", ok ? "MATCH" : "DIFFER");
  return ok;
}

// ---------------------------------------------------------------------------
// --preprocess: run_impl's pre-processing, old vs new (I4). Inputs from the
// synthetic profile: p_mid = P0 * pk^(cp/R), T_mid = pk * theta,
// pseudo_density = rho * g * (layer thickness from z_mid), and qv as a strided
// slice of a tracers array, as EAMxx passes it. Compares pk, theta, dz, rho,
// z_mid and z_int on the real levels bit for bit, and times both versions.
// ---------------------------------------------------------------------------
template <typename KMF>
bool preprocess_test (const std::string& label, const Options& opt, const Profile& prof) {
  using Pack = typename KMF::Pack;
  using Dev  = typename KMF::Device;
  const int ncols = opt.ncols, nz = opt.nz;
  const Real cp = PC::Cpair.value, rd = PC::Rair.value, g = PC::gravit.value;

  std::vector<Real> p(size_t(ncols) * nz), T(p.size()), pdel(p.size());
  for (int c = 0; c < ncols; ++c) {
    for (int k = 0; k < nz; ++k) {
      const size_t i = size_t(c) * nz + k;
      p[i] = PC::P0.value * std::pow(prof.pk[i], cp / rd);
      T[i] = prof.pk[i] * prof.theta[i];
      const Real zu = (k > 0)      ? prof.z_mid[i - 1] : 2 * prof.z_mid[i] - prof.z_mid[i + 1];
      const Real zl = (k < nz - 1) ? prof.z_mid[i + 1] : 0;
      pdel[i] = prof.rho[i] * g * 0.5 * (zu - zl);
    }
  }
  auto p_mid = make_2d<KMF>("p_mid", p, ncols, nz);
  auto T_mid = make_2d<KMF>("T_mid", T, ncols, nz);
  auto pd    = make_2d<KMF>("pseudo_density", pdel, ncols, nz);
  auto qv    = make_2d<KMF>("qv", prof.qv, ncols, nz);
  using TracersView = Kokkos::View<Pack***, Kokkos::LayoutRight, Dev>;
  TracersView tracers("tracers", ncols, 4, ekat::npack<Pack>(nz));
  Kokkos::deep_copy(Kokkos::subview(tracers, Kokkos::ALL, 0, Kokkos::ALL), qv);
  PView2d<KMF> qv_s = ekat::subview_1(tracers, 0);

  const int np = ekat::npack<Pack>(nz), npi = ekat::npack<Pack>(nz + 1);
  auto outputs = [&](const std::string& tag) {
    kessler_driver::PreprocessViews<Dev> v;
    v.T_mid = T_mid; v.p_mid = p_mid; v.pseudo_density = pd; v.qv = qv_s;
    v.pk    = PView2d<KMF>("pk_"    + tag, ncols, np);
    v.theta = PView2d<KMF>("theta_" + tag, ncols, np);
    v.dz    = PView2d<KMF>("dz_"    + tag, ncols, np);
    v.rho   = PView2d<KMF>("rho_"   + tag, ncols, np);
    v.z_mid = PView2d<KMF>("z_mid_" + tag, ncols, np);
    v.z_int = PView2d<KMF>("z_int_" + tag, ncols, npi);
    return v;
  };
  const auto vo = outputs("old"), vn = outputs("new");

  auto time_it = [&](auto&& fn) {
    std::vector<double> t;
    for (int r = 0; r < opt.warmup + opt.repeat; ++r) {
      Kokkos::fence();
      Kokkos::Timer timer;
      fn();
      Kokkos::fence();
      if (r >= opt.warmup) t.push_back(timer.seconds());
    }
    return *std::min_element(t.begin(), t.end());
  };
  const double t_old = time_it([&]() { kessler_driver::preprocess_old<Dev>(ncols, nz, vo); });
  const double t_new = time_it([&]() { kessler_driver::preprocess_new<Dev>(ncols, nz, vn); });

  std::printf("\n=== preprocess (%s): ncols=%d nz=%d, qv column stride %d\n", label.c_str(), ncols, nz,
              int(qv_s.stride(0)));
  bool ok = true;
  struct F { const char* name; PView2d<KMF> a, b; int n; };
  const F fields[] = {{"pk", vo.pk, vn.pk, nz}, {"theta", vo.theta, vn.theta, nz},
                      {"dz", vo.dz, vn.dz, nz}, {"rho", vo.rho, vn.rho, nz},
                      {"z_mid", vo.z_mid, vn.z_mid, nz}, {"z_int", vo.z_int, vn.z_int, nz + 1}};
  for (const auto& f : fields) {
    const auto a = to_vec<KMF>(f.a, f.n), b = to_vec<KMF>(f.b, f.n);
    size_t n_diff = 0; Real max_rel = 0; int n_bad = 0;
    for (size_t i = 0; i < a.size(); ++i) {
      if (!std::isfinite(b[i])) ++n_bad;
      if (a[i] == b[i]) continue;
      ++n_diff;
      const Real rel = std::abs(a[i] - b[i]) / std::max(std::abs(a[i]), Real(1e-300));
      if (!(rel <= max_rel)) max_rel = rel;
    }
    std::printf("  %-6s differing values: %zu of %zu (max rel %.3e), non-finite: %d %s\n", f.name,
                n_diff, a.size(), double(max_rel), n_bad, (n_diff || n_bad) ? "<-- differs" : "");
    ok = ok && n_diff == 0 && n_bad == 0;
  }
  std::printf("  time per call (best of %d): old %.4e s, new %.4e s, old/new %.2fx\n  %s\n",
              opt.repeat, t_old, t_new, t_old / t_new, ok ? "MATCH" : "DIFFER");
  return ok;
}

} // anonymous namespace

// ---------------------------------------------------------------------------
int main (int argc, char** argv)
{
  MPI_Init(&argc, &argv);
  Kokkos::initialize(argc, argv);
  int rc = 0;
  {
    Options opt;
    if (!parse(argc, argv, opt)) {
      Kokkos::finalize();
      MPI_Finalize();
      return 0;
    }

    std::printf("driver_kessler: device exec space = %s (concurrency %d), host exec space = %s, "
                "Pack::n = %d\n",
                Kokkos::DefaultExecutionSpace::name(),
                int(Kokkos::DefaultExecutionSpace().concurrency()),
                Kokkos::DefaultHostExecutionSpace::name(), int(KMF_D::Pack::n));

    if (opt.kernel_timers) install_kernel_timers();

    const Profile prof = build_profile(opt);

    if (opt.preprocess) {
      bool ok = preprocess_test<KMF_D>("device", opt, prof);
      if (opt.host_ref) ok = preprocess_test<KMF_H>("host", opt, prof) && ok;
      std::printf("\ndriver_kessler --preprocess: %s\n", ok ? "ALL MATCH" : "DIFFERENCES (see above)");
      rc = ok ? 0 : 1;
      Kokkos::finalize();
      MPI_Finalize();
      return rc;
    }
    {
      int n_rainy = 0;
      for (char r : prof.rainy) n_rainy += r;
      std::printf("driver_kessler: %d of %d columns rainy (--rain-frac %g)\n",
                  n_rainy, opt.ncols, double(opt.rain_frac));
    }

    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    const bool ws_carved = !opt.managed;  // workspace mode for the NaN run

    const std::vector<Impl>& impls = opt.impls;

    std::vector<RunResult> runs;
    for (const Impl im : impls) {
      if (opt.managed)  runs.push_back(run_case<KMF_D>(im, "device/managed", opt, prof, false, opt.poison_val));
      if (opt.carved)   runs.push_back(run_case<KMF_D>(im, "device/carved",  opt, prof, true,  opt.poison_val));
      // Same device run with the scratch pre-filled with NaN instead: any
      // read of scratch kessler_run did not initialise changes the answer.
      if (opt.poison)   runs.push_back(run_case<KMF_D>(im, "device/nan-poison", opt, prof, ws_carved, nan));
      if (opt.host_ref) runs.push_back(run_case<KMF_H>(im, "host/managed",   opt, prof, false, opt.poison_val));
      // Same device run with qv/qc/qr as strided slices of a tracers array
      // (EAMxx's layout): any kernel that drops a tracer view's stride
      // changes the answer.
      if (opt.strided)  runs.push_back(run_case<KMF_D>(im, "device/strided-tracers", opt, prof,
                                                       ws_carved, opt.poison_val, true));
    }

    auto find = [&runs](const std::string& label) -> const RunResult* {
      for (const auto& r : runs) if (r.label == label) return &r;
      return nullptr;
    };

    bool all_ok = true;
    for (const auto& r : runs) all_ok = check_result(r, opt, prof) && all_ok;

    // Device vs host should agree to round-off; managed vs carved and
    // garbage- vs NaN-poisoned exactly.
    for (const Impl im : impls) {
      const std::string p = std::string(impl_name(im)) + ":";
      const RunResult* host    = find(p + "host/managed");
      const RunResult* managed = find(p + "device/managed");
      const RunResult* carved  = find(p + "device/carved");
      const RunResult* nanrun  = find(p + "device/nan-poison");
      const RunResult* strided = find(p + "device/strided-tracers");
      for (const auto* r : {managed, carved}) {
        if (host && r) all_ok = compare(*r, *host, Real(1e-8), Real(1e-14)) && all_ok;
      }
      if (managed && carved) all_ok = compare(*managed, *carved, Real(0), Real(0)) && all_ok;
      if (nanrun) {
        const RunResult* same_ws = ws_carved ? carved : managed;
        all_ok = compare(*same_ws, *nanrun, Real(0), Real(0)) && all_ok;
      }
      if (strided) {
        const RunResult* same_ws = ws_carved ? carved : managed;
        all_ok = compare(*same_ws, *strided, Real(0), Real(0)) && all_ok;
      }
    }
    // Every implementation does the same arithmetic per column as prod, so
    // each should match prod bit for bit on the same backend.
    auto run_part = [](const std::string& label) { return label.substr(label.find(':') + 1); };
    for (const auto& r : runs) {
      if (r.label.rfind("prod:", 0) == 0) continue;
      const RunResult* pr = find("prod:" + run_part(r.label));
      if (pr) all_ok = compare(r, *pr, Real(0), Real(0)) && all_ok;
    }

    // Timing summary: seconds per kessler_run call over the timed repeats.
    std::printf("\n--- timing (s per kessler_run call, %d timed repeat(s) x %d step(s), "
                "%d warmup)\n", opt.repeat, opt.nsteps, opt.warmup);
    std::printf("  %-32s %12s %12s %12s %14s\n", "run", "min", "mean", "max", "prod/this");
    for (const auto& r : runs) {
      if (r.threw || r.times.empty()) continue;
      const double mn = *std::min_element(r.times.begin(), r.times.end());
      const double mx = *std::max_element(r.times.begin(), r.times.end());
      double mean = 0; for (double t : r.times) mean += t; mean /= r.times.size();
      std::string speedup = "";
      if (r.label.rfind("prod:", 0) != 0) {
        const RunResult* pr = find("prod:" + run_part(r.label));
        if (pr && !pr->times.empty()) {
          char buf[32];
          std::snprintf(buf, sizeof(buf), "%.2fx",
                        *std::min_element(pr->times.begin(), pr->times.end()) / mn);
          speedup = buf;
        }
      }
      std::printf("  %-32s %12.4e %12.4e %12.4e %14s\n", r.label.c_str(), mn, mean, mx,
                  speedup.c_str());
    }
    if (opt.kernel_timers) {
      for (const auto& r : runs) {
        if (r.kernels.empty()) continue;
        double total = 0;
        for (const auto& k : r.kernels) total += k.second.first;
        std::printf("\n--- kernels: %s (summed over timed repeats; fenced)\n", r.label.c_str());
        std::vector<std::pair<std::string, std::pair<double, long>>> ks(r.kernels.begin(),
                                                                        r.kernels.end());
        std::sort(ks.begin(), ks.end(),
                  [](const auto& a, const auto& b) { return a.second.first > b.second.first; });
        for (const auto& k : ks) {
          std::printf("  %-28s %10.4e s %6.1f%%  %9ld calls  %10.3e s/call\n", k.first.c_str(),
                      k.second.first, 100 * k.second.first / std::max(total, 1e-300),
                      k.second.second, k.second.first / std::max(k.second.second, 1L));
        }
      }
    }

    std::printf("\ndriver_kessler: %s\n", all_ok ? "ALL PASS" : "FAILURES (see above)");
    rc = all_ok ? 0 : 1;
  }
  Kokkos::finalize();
  MPI_Finalize();
  return rc;
}
