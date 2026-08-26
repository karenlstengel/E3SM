/*
 * Standalone driver for the Kessler Kokkos microphysics.
 *
 * Mirrors the logic in driver_kessler.F90:
 *  - Allocates host arrays, initialises them with the same Box-Muller
 *    scaled random values used in the Fortran driver.
 *  - Deep-copies data to device, runs kessler_run and the three
 *    kessler_update routines, then deep-copies results back to host.
 *  - Prints the sum of each output field so results can be compared
 *    with the Fortran reference.
 *
 * Build via the CMakeLists.txt in this directory (driver_kessler target).
 * Usage: ./driver_kessler
 */

#include "kessler_functions.hpp"
#include "kessler_functions_impl.hpp"  // needed for non-GPU ETI builds

#include <Kokkos_Core.hpp>

#include <cmath>
#include <cstdio>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

using Real = double;
using Device = Kokkos::DefaultExecutionSpace;
using KF = scream::kessler::KesslerFunctions<Real, Device>;

// Host view aliases
using HView1d = Kokkos::View<Real*,   Kokkos::HostSpace>;
using HView2d = Kokkos::View<Real**, Kokkos::LayoutRight, Kokkos::HostSpace>;

// Device view aliases (managed)
using DView1d = KF::view_1d<Real>;
using DView2d = KF::view_2d<Real>;

// ---------------------------------------------------------------------------
// Helper: sum all elements of a device 1D view (for printing)
// ---------------------------------------------------------------------------
static Real sum1d (const DView1d& v) {
  Real s = 0;
  Kokkos::parallel_reduce(v.extent(0),
    KOKKOS_LAMBDA(int i, Real& acc) { acc += v(i); }, s);
  return s;
}

// Helper: sum all elements of a device 2D view
static Real sum2d (const DView2d& v) {
  const int n0 = v.extent(0), n1 = v.extent(1);
  Real s = 0;
  Kokkos::parallel_reduce(n0 * n1,
    KOKKOS_LAMBDA(int idx, Real& acc) {
      acc += v(idx / n1, idx % n1);
    }, s);
  return s;
}

// ---------------------------------------------------------------------------
int main (int argc, char** argv)
{
  Kokkos::initialize(argc, argv);
  {
    // ------------------------------------------------------------------
    // Grid dimensions (match the Fortran driver)
    // ------------------------------------------------------------------
    const int ncols  = 128;
    const int nz     = 56;
    const Real dt    = 60.0;     // seconds
    const int lyr_surf = 0;      // 0-based C++ index: surface at k=0
    const int lyr_toa  = nz - 1; // top of atmosphere at k=nz-1

    // ------------------------------------------------------------------
    // Kessler physical constants (match the Fortran driver)
    // ------------------------------------------------------------------
    KF::KesslerData kd;
    kd.lv    = 2.5e6;    // J kg-1
    kd.pref  = 1000.0;   // hPa  (Fortran passes 100000 Pa, /100 in init)
    kd.rhoqr = 1000.0;   // kg m-3
    kd.cpair = 1004.0;   // J kg-1 K-1
    kd.rair  = 287.0;    // J kg-1 K-1

    // Persistent scratch workspace, allocated once and reused across calls
    KF::KesslerWorkspace workspace;
    workspace.init(ncols, nz);

    // ------------------------------------------------------------------
    // Build per-column scaling array: Normal(1, 0.1) via Box-Muller
    // (uses a fixed seed for reproducibility, matching the Fortran driver)
    // ------------------------------------------------------------------
    std::mt19937 rng(42);
    std::uniform_real_distribution<Real> udist(0.0, 1.0);
    std::vector<Real> arr(ncols);
    for (int i = 0; i < ncols; ++i) {
      const Real u1 = udist(rng);
      const Real u2 = udist(rng);
      const Real z  = std::sqrt(-2.0 * std::log(u1)) *
                       std::cos(2.0 * M_PI * u2);
      arr[i] = 1.0 + 0.1 * z;
    }

    // ------------------------------------------------------------------
    // Allocate and initialise host arrays
    // ------------------------------------------------------------------
    HView2d h_rho   ("h_rho",    ncols, nz);
    HView2d h_z     ("h_z",      ncols, nz);
    HView2d h_pk    ("h_pk",     ncols, nz);
    HView2d h_theta ("h_theta",  ncols, nz);
    HView2d h_qv    ("h_qv",     ncols, nz);
    HView2d h_qc    ("h_qc",     ncols, nz);
    HView2d h_qr    ("h_qr",     ncols, nz);
    HView2d h_relhum("h_relhum", ncols, nz);
    HView1d h_precl ("h_precl",  ncols);

    // kessler_update arrays
    HView2d h_temp     ("h_temp",     ncols, nz);
    HView2d h_exner    ("h_exner",    ncols, nz);
    HView2d h_zm       ("h_zm",       ncols, nz);
    HView2d h_temp_prev("h_temp_prev",ncols, nz);
    HView2d h_ttend_t  ("h_ttend_t",  ncols, nz);
    HView2d h_st_energy("h_st_energy",ncols, nz);
    HView1d h_phis     ("h_phis",     ncols);

    for (int i = 0; i < ncols; ++i) {
      const Real a = arr[i];
      for (int k = 0; k < nz; ++k) {
        // Fortran: z(i,k) = arr(i) * 100*(k-1)  →  C++ k=0..nz-1
        const Real z_ik = a * (100.0 * k);

        h_z     (i, k) = z_ik;
        h_rho   (i, k) = a * 1.2 * std::exp(-z_ik / 8000.0);
        h_pk    (i, k) = a * 1.0;
        h_theta (i, k) = a * (300.0 - 0.006 * z_ik);
        h_qv    (i, k) = a * 0.010;
        h_qc    (i, k) = a * 0.01;
        h_qr    (i, k) = a * 0.01;
        h_temp  (i, k) = a * 287.4;
        h_exner (i, k) = a * 0.97;
        h_zm    (i, k) = z_ik;  // same as z for this driver
      }
      h_phis(i) = a * 0.1;
    }

    // ------------------------------------------------------------------
    // Deep-copy host arrays to device
    // ------------------------------------------------------------------
    DView2d rho    = Kokkos::create_mirror_view_and_copy(Device(), h_rho);
    DView2d z      = Kokkos::create_mirror_view_and_copy(Device(), h_z);
    DView2d pk     = Kokkos::create_mirror_view_and_copy(Device(), h_pk);
    DView2d theta  = Kokkos::create_mirror_view_and_copy(Device(), h_theta);
    DView2d qv     = Kokkos::create_mirror_view_and_copy(Device(), h_qv);
    DView2d qc     = Kokkos::create_mirror_view_and_copy(Device(), h_qc);
    DView2d qr     = Kokkos::create_mirror_view_and_copy(Device(), h_qr);
    DView2d relhum ("relhum", ncols, nz);
    DView1d precl  ("precl",  ncols);

    DView2d temp      = Kokkos::create_mirror_view_and_copy(Device(), h_temp);
    DView2d exner     = Kokkos::create_mirror_view_and_copy(Device(), h_exner);
    DView2d zm        = Kokkos::create_mirror_view_and_copy(Device(), h_zm);
    DView2d temp_prev ("temp_prev", ncols, nz);
    DView2d ttend_t   ("ttend_t",   ncols, nz);
    DView2d st_energy ("st_energy", ncols, nz);
    DView1d phis      = Kokkos::create_mirror_view_and_copy(Device(), h_phis);

    // ------------------------------------------------------------------
    // Run kessler microphysics
    // ------------------------------------------------------------------
    KF::kessler_run(ncols, nz, dt, lyr_surf, lyr_toa, kd,
                    workspace, rho, z, pk,
                    theta, qv, qc, qr,
                    precl, relhum);

    // ------------------------------------------------------------------
    // kessler_update_timestep_init
    // (Fortran driver had a bug and omitted ncols/nz; fixed here)
    // ------------------------------------------------------------------
    KF::kessler_update_timestep_init(ncols, nz, temp, temp_prev, ttend_t);

    // ------------------------------------------------------------------
    // kessler_update_run
    // ------------------------------------------------------------------
    KF::kessler_update_run(ncols, nz, dt, theta, exner, temp_prev, ttend_t);

    // ------------------------------------------------------------------
    // kessler_update_timestep_final
    // (Fortran driver omitted ncols after nz; fixed here)
    // ------------------------------------------------------------------
    const Real gravit = 9.80616;  // m s-2 (standard EAMxx value)
    KF::kessler_update_timestep_final(ncols, nz, gravit,
                                      kd.cpair, temp, zm, phis, st_energy);

    // ------------------------------------------------------------------
    // Print field sums (matches the Fortran driver output format)
    // ------------------------------------------------------------------
    std::printf("theta:      %20.10e\n", sum2d(theta));
    std::printf("qv:         %20.10e\n", sum2d(qv));
    std::printf("qc:         %20.10e\n", sum2d(qc));
    std::printf("qr:         %20.10e\n", sum2d(qr));
    std::printf("precl:      %20.10e\n", sum1d(precl));
    std::printf("relhum:     %20.10e\n", sum2d(relhum));
    std::printf("temp_prev:  %20.10e\n", sum2d(temp_prev));
    std::printf("ttend_t:    %20.10e\n", sum2d(ttend_t));
    std::printf("st_energy:  %20.10e\n", sum2d(st_energy));

  } // Kokkos scope
  Kokkos::finalize();
  return 0;
}
