/*
 * driver_bridge_transpose.cpp
 *
 * Standalone check and timing of the Fortran bridge's C++ <-> Fortran layout
 * copies (params_helpers / params_computed ::transpose in
 * kessler_functions.hpp) against their previous version
 * (driver_bridge_transpose_ref.hpp, namespace kessler_ref), without the rest
 * of EAMxx. Not part of the CMake build; build_driver_bridge_transpose.sh
 * compiles it with the flags and Kokkos of an f90 case build.
 *
 * Per call, as kessler_eamxx_bridge_run does it:
 *   old: init (zero-fill all Fortran arrays), c2f of everything (plus host
 *        copies of everything on the GPU -> CPU build), Fortran, f2c of
 *        everything
 *   new: c2f of what the Fortran reads, Fortran, f2c of what it writes (on
 *        the GPU -> CPU build, one host copy per direction plus one 1-D array)
 * "Fortran" is simulated: every array the Fortran writes gets a value derived
 * from its inputs, the same for both versions. Checks, every call: the
 * Fortran-layout arrays the Fortran reads are equal after c2f, and every
 * EAMxx-layout array is equal after f2c, bit for bit.
 *
 * Usage: driver_bridge_transpose [--ncols N] [--nlev N] [--repeat N]
 */
#include "kessler_functions.hpp"
#include "driver_bridge_transpose_ref.hpp"

#include <Kokkos_Core.hpp>

#include <cstdio>
#include <cstdlib>
#include <string>
#include <vector>

// The old transposes time their host copies with EAMxx's GPTL timers, which
// this standalone driver does not link
namespace scream {
void start_timer (const std::string&) {}
void stop_timer (const std::string&) {}
}

using scream::Real;
using NewKMF = scream::kessler::KesslerMicrophysicsFunctions<Real, scream::DefaultDevice>;
using OldKMF = scream::kessler_ref::KesslerMicrophysicsFunctions<Real, scream::DefaultDevice>;
using Pack   = NewKMF::Pack;
using ExeSpace = NewKMF::KT::ExeSpace;
using DevMem   = ExeSpace::memory_space;
constexpr bool host_fortran =
#if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
  true;
#else
  false;
#endif

// The EAMxx-layout arrays, shared by both versions (the fields)
struct Fields {
  NewKMF::view_2d<Pack> rho, pk, z_mid, theta, qv, qc, qr, relhum, temp_prev, temp, temp_tend, st_energy;
  NewKMF::view_1d<Real> phis, precl;
  Fields (int nc, int np) {
    for (auto* v : {&rho, &pk, &z_mid, &theta, &qv, &qc, &qr, &relhum, &temp_prev, &temp, &temp_tend, &st_energy})
      *v = NewKMF::view_2d<Pack>("f", nc, np);
    phis = NewKMF::view_1d<Real>("phis", nc);
    precl = NewKMF::view_1d<Real>("precl", nc);
  }
  std::vector<NewKMF::view_2d<Pack>> all2 () const {
    return {rho, pk, z_mid, theta, qv, qc, qr, relhum, temp_prev, temp, temp_tend, st_energy};
  }
};

// Point a version's structs at the fields and carve its Fortran-layout arrays
// from one device block, in the order init_buffers uses (new order for both)
template <typename PH, typename PC>
void setup (PH& ph, PC& pc, Fields& F, Real* fb, int nc, int nl) {
  ph.rho = F.rho; ph.pk = F.pk; ph.z_mid = F.z_mid; ph.phis = F.phis;
  pc.theta = F.theta; pc.qv = F.qv; pc.qc = F.qc; pc.qr = F.qr; pc.relhum = F.relhum;
  pc.temp_prev = F.temp_prev; pc.temp = F.temp; pc.temp_tend = F.temp_tend; pc.st_energy = F.st_energy;
  pc.precl = F.precl;
  const size_t n2 = size_t(nc) * nl;
  using V2 = decltype(ph.f_rho);
  using V1 = decltype(ph.f_phis);
  V2* f2[12] = {&ph.f_rho, &ph.f_pk, &ph.f_z_mid, &pc.f_theta, &pc.f_qv, &pc.f_qc, &pc.f_qr, &pc.f_temp,
                &pc.f_relhum, &pc.f_temp_prev, &pc.f_temp_tend, &pc.f_st_energy};
  for (auto* v : f2) { *v = V2(fb, nc, nl); fb += n2; }
  pc.f_precl = V1(fb, nc); fb += nc;
  ph.f_phis  = V1(fb, nc);
}

// What the Fortran does to the Fortran-layout arrays, simulated (on the
// device, or on the host copies on the GPU -> CPU build)
template <typename PH, typename PC>
void fake_fortran (PH& ph, PC& pc, int nc, int nl) {
  if constexpr (host_fortran) {
    for (int k = 0; k < nl; ++k)
      for (int i = 0; i < nc; ++i) {
        pc.h_theta(i,k)     = pc.h_theta(i,k) * 1.001 + ph.h_pk(i,k);
        pc.h_qv(i,k)        = pc.h_qv(i,k) * 0.999;
        pc.h_qc(i,k)        = pc.h_qc(i,k) + 1.0e-6;
        pc.h_qr(i,k)        = pc.h_qr(i,k) * 0.5 + ph.h_rho(i,k);
        pc.h_relhum(i,k)    = pc.h_qv(i,k) * 100.0;
        pc.h_temp_prev(i,k) = pc.h_temp(i,k);
        pc.h_temp_tend(i,k) = pc.h_theta(i,k) - pc.h_temp(i,k);
        pc.h_temp(i,k)      = pc.h_temp(i,k) + pc.h_temp_tend(i,k);
        pc.h_st_energy(i,k) = pc.h_temp(i,k) * 1004.0 + ph.h_z_mid(i,k) + ph.h_phis(i);
      }
    for (int i = 0; i < nc; ++i) pc.h_precl(i) = pc.h_qr(i, nl - 1);
  } else {
    auto t = pc.f_theta; auto qv = pc.f_qv; auto qc = pc.f_qc; auto qr = pc.f_qr; auto rh = pc.f_relhum;
    auto tp = pc.f_temp_prev; auto tt = pc.f_temp_tend; auto T = pc.f_temp; auto se = pc.f_st_energy;
    auto pr = pc.f_precl; auto rho = ph.f_rho; auto pk = ph.f_pk; auto z = ph.f_z_mid; auto ps = ph.f_phis;
    Kokkos::parallel_for("fake fortran", Kokkos::RangePolicy<ExeSpace>(0, nc), KOKKOS_LAMBDA (const int i) {
      for (int k = 0; k < nl; ++k) {
        t(i,k)  = t(i,k) * 1.001 + pk(i,k);
        qv(i,k) = qv(i,k) * 0.999;
        qc(i,k) = qc(i,k) + 1.0e-6;
        qr(i,k) = qr(i,k) * 0.5 + rho(i,k);
        rh(i,k) = qv(i,k) * 100.0;
        tp(i,k) = T(i,k);
        tt(i,k) = t(i,k) - T(i,k);
        T(i,k)  = T(i,k) + tt(i,k);
        se(i,k) = T(i,k) * 1004.0 + z(i,k) + ps(i);
      }
      pr(i) = qr(i, nl - 1);
    });
  }
}

// All values of a view (of Reals or Packs) on the host
template <typename V>
std::vector<Real> to_vec (const V& v) {
  auto h = Kokkos::create_mirror_view_and_copy(Kokkos::HostSpace(), v);
  using T = typename decltype(h)::value_type;
  const Real* p = reinterpret_cast<const Real*>(h.data());
  return std::vector<Real>(p, p + h.size() * (sizeof(T) / sizeof(Real)));
}

void fill (Fields& F, int nc, int np, int seed) {
  int f = 0;
  for (auto v : F.all2()) {
    const Real base = 1.0 + f++ + 0.01 * seed;
    Kokkos::parallel_for("fill", Kokkos::RangePolicy<ExeSpace>(0, nc), KOKKOS_LAMBDA (const int i) {
      for (int p = 0; p < np; ++p)
        for (int l = 0; l < Pack::n; ++l) v(i, p)[l] = base + 0.001 * i + 1.0e-6 * (p * Pack::n + l);
    });
  }
  auto phis = F.phis; auto precl = F.precl;
  Kokkos::parallel_for("fill1", Kokkos::RangePolicy<ExeSpace>(0, nc), KOKKOS_LAMBDA (const int i) {
    phis(i) = 0.5 * i; precl(i) = -1.0;
  });
  Kokkos::fence();
}

#if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
// kessler_eamxx_bridge.cpp's copy_range
void copy_range (Real* dst, const Real* src, size_t n, bool to_host) {
  using DevRange  = Kokkos::View<Real*, DevMem, Kokkos::MemoryUnmanaged>;
  using HostRange = Kokkos::View<Real*, Kokkos::HostSpace, Kokkos::MemoryUnmanaged>;
  if (to_host) Kokkos::deep_copy(HostRange(dst, n), DevRange(const_cast<Real*>(src), n));
  else         Kokkos::deep_copy(DevRange(dst, n), HostRange(const_cast<Real*>(src), n));
}
#endif

int main (int argc, char** argv) {
  int nc = 48601, nl = 58, repeat = 20;
  for (int a = 1; a + 1 < argc; a += 2) {
    const std::string k = argv[a];
    if (k == "--ncols") nc = std::atoi(argv[a+1]);
    else if (k == "--nlev") nl = std::atoi(argv[a+1]);
    else if (k == "--repeat") repeat = std::atoi(argv[a+1]);
  }
  Kokkos::initialize(argc, argv);
  int bad = 0;
  {
    const int np = (nl + Pack::n - 1) / Pack::n;
    const size_t n2 = size_t(nc) * nl;
    std::printf("driver_bridge_transpose: ncols=%d nlev=%d, Kokkos %s, Pack::n=%d, Fortran on the %s\n",
                nc, nl, ExeSpace::name(), int(Pack::n), host_fortran ? "host (GPU -> CPU bridge)" : "same space");
    Fields Fo(nc, np), Fn(nc, np);
    Kokkos::View<Real*, DevMem> fb_old("fb_old", 12 * n2 + 2 * nc), fb_new("fb_new", 12 * n2 + 2 * nc);
    OldKMF::params_helpers oph; OldKMF::params_computed opc;
    NewKMF::params_helpers nph; NewKMF::params_computed npc;
    setup(oph, opc, Fo, fb_old.data(), nc, nl);
    setup(nph, npc, Fn, fb_new.data(), nc, nl);
#if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
    // old: one managed host mirror per array; new: one host block, same order as the device
    for (auto* v : {&oph.h_rho, &oph.h_pk, &oph.h_z_mid}) *v = decltype(oph.h_rho)("h", nc, nl);
    for (auto* v : {&opc.h_theta, &opc.h_qv, &opc.h_qc, &opc.h_qr, &opc.h_relhum, &opc.h_temp_prev,
                    &opc.h_temp, &opc.h_temp_tend, &opc.h_st_energy}) *v = decltype(opc.h_theta)("h", nc, nl);
    oph.h_phis = decltype(oph.h_phis)("h", nc); opc.h_precl = decltype(opc.h_precl)("h", nc);
    Kokkos::View<Real*, Kokkos::HostSpace> hb("hb", 12 * n2 + 2 * nc);
    {
      Real* h = hb.data();
      NewKMF::view_2dh<Real>* h2[12] = {&nph.h_rho, &nph.h_pk, &nph.h_z_mid, &npc.h_theta, &npc.h_qv, &npc.h_qc,
                                        &npc.h_qr, &npc.h_temp, &npc.h_relhum, &npc.h_temp_prev, &npc.h_temp_tend,
                                        &npc.h_st_energy};
      for (auto* v : h2) { *v = NewKMF::view_2dh<Real>(h, nc, nl); h += n2; }
      npc.h_precl = NewKMF::view_1dh<Real>(h, nc); h += nc;
      nph.h_phis  = NewKMF::view_1dh<Real>(h, nc);
    }
#endif
    using TD = ekat::TransposeDirection;
    // The 9 arrays the Fortran reads, as it sees them (host copies on the GPU -> CPU build)
    auto fortran_inputs = [&] (auto& ph, auto& pc) {
      std::vector<std::vector<Real>> v;
#if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
      auto hv = [] (const auto& a) { return std::vector<Real>(a.data(), a.data() + a.size()); };
      for (auto x : {hv(ph.h_rho), hv(ph.h_pk), hv(ph.h_z_mid), hv(ph.h_phis), hv(pc.h_theta), hv(pc.h_qv),
                     hv(pc.h_qc), hv(pc.h_qr), hv(pc.h_temp)}) v.push_back(x);
#else
      for (auto x : {to_vec(ph.f_rho), to_vec(ph.f_pk), to_vec(ph.f_z_mid), to_vec(ph.f_phis), to_vec(pc.f_theta),
                     to_vec(pc.f_qv), to_vec(pc.f_qc), to_vec(pc.f_qr), to_vec(pc.f_temp)}) v.push_back(x);
#endif
      return v;
    };
    std::vector<std::vector<Real>> old_in;
    double t_old[3] = {1e30, 1e30, 1e30}, t_new[2] = {1e30, 1e30};
    for (int r = 0; r < repeat; ++r) {
      fill(Fo, nc, np, r); fill(Fn, nc, np, r);
      Kokkos::Timer tm;
      // old
      oph.init(nc, nl); opc.init(nc, nl); Kokkos::fence();
      t_old[0] = std::min(t_old[0], tm.seconds()); tm.reset();
      oph.transpose<TD::c2f>(nc, nl); opc.transpose<TD::c2f>(nc, nl); Kokkos::fence();
      t_old[1] = std::min(t_old[1], tm.seconds());
      if (r == 0) old_in = fortran_inputs(oph, opc);
      fake_fortran(oph, opc, nc, nl); Kokkos::fence(); tm.reset();
      oph.transpose<TD::f2c>(nc, nl); opc.transpose<TD::f2c>(nc, nl); Kokkos::fence();
      t_old[2] = std::min(t_old[2], tm.seconds());
      // new
      tm.reset();
      nph.transpose<TD::c2f>(nc, nl); npc.transpose<TD::c2f>(nc, nl);
#if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
      copy_range(nph.h_rho.data(), nph.f_rho.data(), 8 * n2, true);
      copy_range(nph.h_phis.data(), nph.f_phis.data(), nc, true);
#endif
      Kokkos::fence();
      t_new[0] = std::min(t_new[0], tm.seconds());
      // c2f check: what the Fortran reads must match
      if (r == 0) {
        const char* in_names[9] = {"rho", "pk", "z_mid", "phis", "theta", "qv", "qc", "qr", "temp"};
        auto new_in = fortran_inputs(nph, npc);
        for (int f = 0; f < 9; ++f)
          if (new_in[f] != old_in[f]) { ++bad; std::printf("  c2f %s differs\n", in_names[f]); }
      }
      fake_fortran(nph, npc, nc, nl); Kokkos::fence(); tm.reset();
#if defined(EAMXX_ENABLE_GPU) && !defined(EAMXX_ENABLE_OPENACC)
      copy_range(npc.f_theta.data(), npc.h_theta.data(), 9 * n2, false);
      copy_range(npc.f_precl.data(), npc.h_precl.data(), nc, false);
#endif
      npc.transpose<TD::f2c>(nc, nl); Kokkos::fence();
      t_new[1] = std::min(t_new[1], tm.seconds());
      // f2c check: every EAMxx-layout array must match
      auto a = Fo.all2(), b = Fn.all2();
      const char* names[12] = {"rho", "pk", "z_mid", "theta", "qv", "qc", "qr", "relhum", "temp_prev", "temp",
                               "temp_tend", "st_energy"};
      for (int f = 0; f < 12; ++f)
        if (to_vec(a[f]) != to_vec(b[f])) { ++bad; std::printf("  call %d: f2c %s differs\n", r, names[f]); }
      if (to_vec(Fo.precl) != to_vec(Fn.precl)) { ++bad; std::printf("  call %d: f2c precl differs\n", r); }
      if (to_vec(Fo.phis) != to_vec(Fn.phis)) { ++bad; std::printf("  call %d: phis differs\n", r); }
    }
    std::printf("  s/call (min of %d)   old: init %.3e  c2f %.3e  f2c %.3e  total %.3e\n", repeat,
                t_old[0], t_old[1], t_old[2], t_old[0] + t_old[1] + t_old[2]);
    std::printf("                       new:              c2f %.3e  f2c %.3e  total %.3e   (%.2fx)\n",
                t_new[0], t_new[1], t_new[0] + t_new[1],
                (t_old[0] + t_old[1] + t_old[2]) / (t_new[0] + t_new[1]));
    std::printf("  %s\n", bad ? "DIFFER" : "MATCH: every array bit for bit, every call");
  }
  Kokkos::finalize();
  return bad ? 1 : 0;
}
