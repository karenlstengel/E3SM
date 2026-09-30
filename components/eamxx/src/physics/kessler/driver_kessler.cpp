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
 * Exit code 0 only if everything passes.
 *
 * How to read the result:
 *   - device fails, host passes          -> bug in kessler_run's device path
 *   - carved fails, managed passes       -> Workspace carving/aliasing bug
 *   - poison comparison differs          -> kessler_run reads uninitialised scratch
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

#include <Kokkos_Core.hpp>
#include <mpi.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <exception>
#include <limits>
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
  bool rain       = false;  // add cloud/rain so sedimentation sub-cycles
};

void usage () {
  std::printf(
    "Usage: driver_kessler [options]\n"
    "  --ncols N         number of columns (default 3433)\n"
    "  --nz N            number of levels (default 58)\n"
    "  --dt S            physics time step in seconds (default 1800)\n"
    "  --nsteps N        repeated kessler_run calls (default 1)\n"
    "  --ws MODE         workspace: managed | carved | both (default both)\n"
    "  --no-host-ref     skip the host reference run\n"
    "  --no-poison       do not pre-fill the workspace with garbage (also skips\n"
    "                    the NaN-poison comparison run)\n"
    "  --poison-val X    garbage value for the workspace (default 6966.7)\n"
    "  --rain            add cloud water and rain so sedimentation sub-cycles\n"
    "                    (note: the TEMPORARY 'n_iter < 10' limit in\n"
    "                    kessler_run will trip if more than 10 sub-cycles\n"
    "                    are needed)\n");
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
    else if (a == "--no-host-ref") opt.host_ref   = false;
    else if (a == "--no-poison")   opt.poison     = false;
    else if (a == "--poison-val")  opt.poison_val = std::atof(next());
    else if (a == "--rain")        opt.rain       = true;
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

  for (int c = 0; c < ncols; ++c) {
    const Real two_pi = 6.283185307179586;
    const Real pert = std::sin(two_pi * c / std::max(ncols, 1));
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

      if (opt.rain && (c % 2 == 0)) {
        if (z > 2000.0 && z < 6000.0) p.qc[i] = 1.0e-3;
        if (z > 1000.0 && z < 4000.0) p.qr[i] = 5.0e-4;
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
// One kessler_run experiment
// ---------------------------------------------------------------------------
struct RunResult {
  std::string label;
  bool threw = false;
  std::string err;
  std::vector<Real> theta, qv, qc, qr, relhum, precl, mask, time_counter;
};

template <typename KMF>
RunResult run_case (const std::string& label, const Options& opt,
                    const Profile& prof, const bool carved,
                    const Real poison_val) {
  const int ncols = opt.ncols, nz = opt.nz;

  RunResult res;
  res.label = label;

  auto rho    = make_2d<KMF>("rho",    prof.rho,   ncols, nz);
  auto z_mid  = make_2d<KMF>("z_mid",  prof.z_mid, ncols, nz);
  auto pk     = make_2d<KMF>("pk",     prof.pk,    ncols, nz);
  auto theta  = make_2d<KMF>("theta",  prof.theta, ncols, nz);
  auto qv     = make_2d<KMF>("qv",     prof.qv,    ncols, nz);
  auto qc     = make_2d<KMF>("qc",     prof.qc,    ncols, nz);
  auto qr     = make_2d<KMF>("qr",     prof.qr,    ncols, nz);
  auto relhum = make_2d<KMF>("relhum", std::vector<Real>(size_t(ncols) * nz, 0.0), ncols, nz);
  RView1d<KMF> precl("precl", ncols);

  auto wsh = make_workspace<KMF>(carved, ncols, nz);
  if (opt.poison) poison_workspace<KMF>(wsh.ws, poison_val);

  // Same arguments the process interface passes (see
  // eamxx_kessler_process_interface.cpp): top-down levels, pref in hPa.
  const int  lyr_surf = nz - 1;
  const int  lyr_toa  = 0;
  const Real rhoqr    = PC::RHOW.value;
  const Real latvap   = PC::LatVap.value;
  const Real pref     = PC::P0.value / 100.0;
  const Real cpair    = PC::Cpair.value;
  const Real rair     = PC::Rair.value;

  std::printf("\n=== %s: %d step(s), ncols=%d nz=%d dt=%g, workspace=%s%s\n",
              label.c_str(), opt.nsteps, ncols, nz, double(opt.dt),
              carved ? "carved" : "managed",
              opt.poison ? (", poisoned with " + std::to_string(poison_val)).c_str() : "");
  try {
    for (int s = 0; s < opt.nsteps; ++s) {
      KMF::kessler_run(ncols, nz, opt.dt, lyr_surf, lyr_toa,
                       rhoqr, latvap, pref, cpair, rair,
                       rho, z_mid, pk, theta, qv, qc, qr,
                       precl, relhum, wsh.ws);
      Kokkos::fence();
    }
  } catch (const std::exception& e) {
    res.threw = true;
    res.err   = e.what();
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

  if (!opt.rain) {
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
      const bool rainy = (c % 2 == 0);
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

    std::printf("driver_kessler: device exec space = %s, host exec space = %s, Pack::n = %d\n",
                Kokkos::DefaultExecutionSpace::name(),
                Kokkos::DefaultHostExecutionSpace::name(), int(KMF_D::Pack::n));

    const Profile prof = build_profile(opt);

    const Real nan = std::numeric_limits<Real>::quiet_NaN();
    const bool ws_carved = !opt.managed;  // workspace mode for the NaN run

    std::vector<RunResult> runs;
    if (opt.managed)  runs.push_back(run_case<KMF_D>("device/managed", opt, prof, false, opt.poison_val));
    if (opt.carved)   runs.push_back(run_case<KMF_D>("device/carved",  opt, prof, true,  opt.poison_val));
    // Same device run with the scratch pre-filled with NaN instead: any
    // read of scratch kessler_run did not initialise changes the answer.
    if (opt.poison)   runs.push_back(run_case<KMF_D>("device/nan-poison", opt, prof, ws_carved, nan));
    if (opt.host_ref) runs.push_back(run_case<KMF_H>("host/managed",   opt, prof, false, opt.poison_val));

    auto find = [&runs](const std::string& label) -> const RunResult* {
      for (const auto& r : runs) if (r.label == label) return &r;
      return nullptr;
    };

    bool all_ok = true;
    for (const auto& r : runs) all_ok = check_result(r, opt, prof) && all_ok;

    // Device vs host should agree to round-off; managed vs carved and
    // garbage- vs NaN-poisoned exactly.
    const RunResult* host    = find("host/managed");
    const RunResult* managed = find("device/managed");
    const RunResult* carved  = find("device/carved");
    const RunResult* nanrun  = find("device/nan-poison");
    for (const auto* r : {managed, carved}) {
      if (host && r) all_ok = compare(*r, *host, Real(1e-8), Real(1e-14)) && all_ok;
    }
    if (managed && carved) all_ok = compare(*managed, *carved, Real(0), Real(0)) && all_ok;
    if (nanrun) {
      const RunResult* same_ws = ws_carved ? carved : managed;
      all_ok = compare(*same_ws, *nanrun, Real(0), Real(0)) && all_ok;
    }

    std::printf("\ndriver_kessler: %s\n", all_ok ? "ALL PASS" : "FAILURES (see above)");
    rc = all_ok ? 0 : 1;
  }
  Kokkos::finalize();
  MPI_Finalize();
  return rc;
}
