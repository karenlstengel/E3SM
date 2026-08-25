#ifndef DCMIP2016_FUNCTIONS_IMPL_HPP
#define DCMIP2016_FUNCTIONS_IMPL_HPP

// Template function definitions for BaroclinicWaveFunctions<ScalarT, DeviceT>.
//
// This header is included by dcmip2016.cpp (for explicit template instantiation)
// and should not be included elsewhere — use dcmip2016_functions.hpp instead.
//
// All intermediate arithmetic in the analytic formulas is done in double
// precision (via the Params constants and local double variables) so that the
// secant-method iteration converges regardless of whether ScalarT is float or
// double.  Results are cast to Scalar only at the final store.

#include "analytic_conditions/dcmip2016/dcmip2016_functions.hpp"

#include <ekat_pack_kokkos.hpp>
#include <ekat_subview_utils.hpp>

namespace scream {
namespace dcmip2016 {

// ==========================================================================
// eval_pressure_temperature
//
// Port of Fortran subroutine evaluate_pressure_temperature() from
// dcmip2016-baroclinic.F90.
//
// Given geometric altitude z (m) and (lon, lat) in radians, returns the
// hydrostatic pressure p (Pa) and temperature T (K) of the background state.
// Purely algebraic — GPU-safe (KOKKOS_FUNCTION).
// ==========================================================================
template <typename S, typename D>
KOKKOS_FUNCTION
void BaroclinicWaveFunctions<S,D>::eval_pressure_temperature(
    int deep, double X, double a, double Rd, double lon, double lat, double z,
    double& p, double& T)
{
  using P = Params;

  const double aref    = a / X;
  const double T0      = 0.5 * (P::T0E + P::T0P);
  const double constA  = 1.0 / P::lapse;
  const double constB  = (T0 - P::T0P) / (T0 * P::T0P);
  const double constC  = 0.5 * (P::K + 2.0) * (P::T0E - P::T0P) / (P::T0E * P::T0P);
  const double constH  = Rd * T0 / P::g;
  const double scaledZ = z / (P::B * constH);

  // Tau functions: weighting factors for the hydrostatic integration
  const double tau1    = constA * P::lapse / T0 * Kokkos::exp(P::lapse * z / T0)
                       + constB * (1.0 - 2.0*scaledZ*scaledZ) * Kokkos::exp(-scaledZ*scaledZ);
  const double tau2    = constC * (1.0 - 2.0*scaledZ*scaledZ) * Kokkos::exp(-scaledZ*scaledZ);

  // Vertically integrated tau functions
  const double inttau1 = constA * (Kokkos::exp(P::lapse * z / T0) - 1.0)
                       + constB * z * Kokkos::exp(-scaledZ*scaledZ);
  const double inttau2 = constC * z * Kokkos::exp(-scaledZ*scaledZ);

  // Radius ratio: 1 for shallow atmosphere, (z+a)/a for deep
  const double rratio  = (deep == 0) ? 1.0 : (z + aref) / aref;

  // Latitudinal structure factor (geopotential integral)
  const double inttermT = Kokkos::pow(rratio * Kokkos::cos(lat), P::K)
                        - (P::K / (P::K + 2.0))
                          * Kokkos::pow(rratio * Kokkos::cos(lat), P::K + 2.0);

  // Temperature from the inverse of the tau product
  T = 1.0 / (rratio * rratio * (tau1 - tau2 * inttermT));

  // Hydrostatic pressure from vertical integration
  p = P::p0 * Kokkos::exp(-P::g / Rd * (inttau1 - inttau2 * inttermT));
}

// ==========================================================================
// eval_z_from_p
//
// Port of Fortran subroutine evaluate_z_temperature() from
// dcmip2016-baroclinic.F90.
//
// Given pressure p_tgt (Pa), finds altitude z (m) and temperature T (K)
// using the secant method on eval_pressure_temperature.
//
// The sequential for-loop is valid inside a KOKKOS_FUNCTION: each Kokkos
// thread (one per column) runs the loop independently.  There are no
// synchronization primitives inside the loop, so GPU threads do not
// interfere with each other.
//
// All arithmetic is in double precision so the iteration converges even when
// ScalarT = float.
// ==========================================================================
template <typename S, typename D>
KOKKOS_FUNCTION
void BaroclinicWaveFunctions<S,D>::eval_z_from_p(
    int deep, double X, double a, double Rd, double lon, double lat, double p_tgt,
    double& z, double& T)
{
  // Initial bracket: sea level and 10 km
  double z0 = 0.0,     p0_tmp, T_tmp;
  double z1 = 10000.0, p1_tmp;
  eval_pressure_temperature(deep, X, a, Rd, lon, lat, z0, p0_tmp, T_tmp);
  eval_pressure_temperature(deep, X, a, Rd, lon, lat, z1, p1_tmp, T_tmp);

  // Secant iteration (mirrors the Fortran DO loop up to 1000 steps)
  double z2 = z1, p2_tmp = p1_tmp;
  for (int ix = 0; ix < 1000; ++ix) {
    z2 = z1 - (p1_tmp - p_tgt) * (z1 - z0) / (p1_tmp - p0_tmp);
    eval_pressure_temperature(deep, X, a, Rd, lon, lat, z2, p2_tmp, T_tmp);

    // Converge when relative pressure error is below 1e-13
    if (Kokkos::abs((p2_tmp - p_tgt) / p_tgt) < 1.e-13) break;

    z0 = z1; p0_tmp = p1_tmp;
    z1 = z2; p1_tmp = p2_tmp;
  }

  z = z2;
  // Final T evaluation at the converged altitude
  eval_pressure_temperature(deep, X, a, Rd, lon, lat, z, p2_tmp, T);
}

// ==========================================================================
// eval_exponential
//
// Port of Fortran function evaluate_exponential() from dcmip2016-baroclinic.F90.
//
// Returns the zonal-wind perturbation u' for perturbation type pertt=0.
// Algebraic — GPU-safe.
// ==========================================================================
template <typename S, typename D>
KOKKOS_FUNCTION
double BaroclinicWaveFunctions<S,D>::eval_exponential(double lon, double lat, double z)
{
  using P = Params;

  // Non-dimensionalized great-circle distance from perturbation center
  const double gc = (1.0 / P::pertexpr)
    * Kokkos::acos(Kokkos::sin(P::pertlat) * Kokkos::sin(lat)
                 + Kokkos::cos(P::pertlat) * Kokkos::cos(lat)
                   * Kokkos::cos(lon - P::pertlon));

  // Smooth polynomial vertical taper: 1 at z=0, 0 at z≥pertz
  const double taper = (z < P::pertz)
    ? 1.0 - 3.0*(z*z)/(P::pertz*P::pertz)
          + 2.0*(z*z*z)/(P::pertz*P::pertz*P::pertz)
    : 0.0;

  return (gc < 1.0) ? P::pertup * taper * Kokkos::exp(-gc * gc) : 0.0;
}

// ==========================================================================
// eval_streamfunction
//
// Port of Fortran function evaluate_streamfunction() from
// dcmip2016-baroclinic.F90.
//
// Returns stream-function value for perturbation type pertt=1.
// u' and v' are derived from this by numerical derivatives in wave_at_point.
// Algebraic — GPU-safe.
// ==========================================================================
template <typename S, typename D>
KOKKOS_FUNCTION
double BaroclinicWaveFunctions<S,D>::eval_streamfunction(double lon, double lat, double z)
{
  using P = Params;

  const double gc = (1.0 / P::pertr)
    * Kokkos::acos(Kokkos::sin(P::pertlat) * Kokkos::sin(lat)
                 + Kokkos::cos(P::pertlat) * Kokkos::cos(lat)
                   * Kokkos::cos(lon - P::pertlon));

  const double taper = (z < P::pertz)
    ? 1.0 - 3.0*(z*z)/(P::pertz*P::pertz)
          + 2.0*(z*z*z)/(P::pertz*P::pertz*P::pertz)
    : 0.0;

  const double cospert = (gc < 1.0) ? Kokkos::cos(0.5 * P::pi * gc) : 0.0;

  return -P::pertu0 * P::pertr * taper * cospert * cospert * cospert * cospert;
}

// ==========================================================================
// wave_at_point
//
// Port of Fortran subroutine baroclinic_wave_test() with zcoords=0.
// Computes the complete DCMIP2016 Test 1 state at a single grid point.
//
// Calls eval_z_from_p (secant iteration) — valid on GPU because each
// thread handles one column and runs this function for each level
// sequentially (see main() below).
//
// Parameters:
//   deep  : 0 = shallow atmosphere (standard for test 1)
//   moist : 1 = include moisture, 0 = dry
//   pertt : 0 = exponential perturbation, 1 = stream-function perturbation
//   X     : Earth scaling factor (1.0 = full Earth)
//   lon   : longitude in radians
//   lat   : latitude in radians
//   p     : mid-level pressure in Pa
// ==========================================================================
template <typename S, typename D>
KOKKOS_FUNCTION
typename BaroclinicWaveFunctions<S,D>::State
BaroclinicWaveFunctions<S,D>::wave_at_point(
    int deep, int moist, int pertt, double X, double a, double Rd, double Rvap,
    double lon, double lat, double p)
{
  using P = Params;
  State s;

  // ---- 1. Invert p -> z, get dry temperature ----
  double z, T;
  eval_z_from_p(deep, X, a, Rd, lon, lat, p, z, T);

  // ---- 2. Background zonal wind from thermal-wind balance ----
  const double aref       = a / X;
  const double omegaref   = P::omega * X;
  const double T0         = 0.5 * (P::T0E + P::T0P);
  const double constH     = Rd * T0 / P::g;
  const double constC     = 0.5 * (P::K + 2.0) * (P::T0E - P::T0P) / (P::T0E * P::T0P);
  const double scaledZ    = z / (P::B * constH);
  const double inttau2    = constC * z * Kokkos::exp(-scaledZ * scaledZ);
  const double rratio     = (deep == 0) ? 1.0 : (z + aref) / aref;

  // Thermal-wind integrand (Ullrich et al. 2015, eq. 7)
  const double inttermU   = Kokkos::pow(rratio * Kokkos::cos(lat), P::K - 1.0)
                          - Kokkos::pow(rratio * Kokkos::cos(lat), P::K + 1.0);
  const double bigU       = P::g / aref * P::K * inttau2 * inttermU * T;
  const double rcoslat    = (deep == 0) ? aref * Kokkos::cos(lat)
                                        : (z + aref) * Kokkos::cos(lat);
  const double omrcl      = omegaref * rcoslat;

  double u = -omrcl + Kokkos::sqrt(omrcl * omrcl + rcoslat * bigU);
  double v = 0.0;

  // ---- 3. Add wind perturbation to seed baroclinic instability ----
  if (pertt == 0) {
    u += eval_exponential(lon, lat, z);
  } else if (pertt == 1) {
    const double dx = P::dxepsilon;
    u -= 1.0 / (2.0 * dx)
       * (eval_streamfunction(lon, lat + dx, z)
        - eval_streamfunction(lon, lat - dx, z));
    v += 1.0 / (2.0 * dx * Kokkos::cos(lat))
       * (eval_streamfunction(lon + dx, lat, z)
        - eval_streamfunction(lon - dx, lat, z));
  }

  s.u    = Scalar(u);
  s.v    = Scalar(v);
  s.ps   = Scalar(P::p0);
  s.phis = Scalar(0.0);

  // ---- 4. Water-vapor specific humidity ----
  double qv = 0.0;
  if (moist == 1) {
    const double eta = p / P::p0;
    if (eta > P::moisttr) {
      qv = P::moistq0
         * Kokkos::exp(-Kokkos::pow(lat / P::moistqlat, 4))
         * Kokkos::exp(-Kokkos::pow((eta - 1.0) * P::p0 / P::moistqp, 2));
    } else {
      qv = P::moistqs;
    }
    // Recover actual T from virtual T: T_v = T*(1 + Mvap*qv)  =>  T = T_v/(1+Mvap*qv)
    const double Mvap = Rvap / Rd - 1.0;   // ~0.608
    T = T / (1.0 + Mvap * qv);
  }
  s.T  = Scalar(T);
  s.qv = Scalar(qv);

  return s;
}

// ==========================================================================
// main
//
// Kernel launcher: fills T_mid, horiz_winds, ps, phis, qv, qc, qr for all
// columns and levels using the DCMIP2016 Test 1 analytic formulas.
//
// Kernel structure:
//   - Outer: RangePolicy over ncols — one Kokkos thread per column.
//   - Inner: sequential for-loop over nlevs inside each thread.
//
// The level loop is sequential (not TeamVectorRange) because wave_at_point
// contains a secant iteration that must complete before the next level.
// Column independence ensures full GPU parallelism between columns.
//
// Views:
//   lat_deg / lon_deg : geometry data, degrees, scalar views (Real*)
//   hyam / hybm       : hybrid-pressure coefficients, scalar views
//   T_mid / horiz_winds / qv / qc / qr : Pack views from the field manager
//   ps / phis         : scalar 2D (ncols) views
// ==========================================================================
template <typename S, typename D>
void BaroclinicWaveFunctions<S,D>::main(
    int ncols, int nlevs,
    int deep, int moist, int pertt, Scalar X, Scalar a, Scalar Rd, Scalar Rvap,
    const view_1d<const Scalar>& lat_deg,
    const view_1d<const Scalar>& lon_deg,
    const view_1d<const Scalar>& hyam,
    const view_1d<const Scalar>& hybm,
    const view_2d<Pack>& T_mid,
    const view_3d<Pack>& horiz_winds,
    const view_1d<Scalar>& ps,
    const view_1d<Scalar>& phis,
    const view_2d<Pack>& qv,
    const view_2d<Pack>& qc,
    const view_2d<Pack>& qr)
{
  using ExeSpace = typename KT::ExeSpace;

  // Capture as plain values so the KOKKOS_LAMBDA can copy them
  const double X_d    = double(X);
  const double a_d    = double(a);
  const double Rd_d   = double(Rd);
  const double Rvap_d = double(Rvap);
  const int    deep_  = deep;
  const int    moist_ = moist;
  const int    pertt_ = pertt;
  const double p0     = Params::p0;
  const double pi     = Params::pi;
  const int    nlevs_ = nlevs;

  Kokkos::parallel_for(
    "dcmip2016_baroclinic_wave_ic",
    Kokkos::RangePolicy<ExeSpace>(0, ncols),
    KOKKOS_LAMBDA(const int icol) {

      // Convert geometry from degrees to radians (Params::pi / 180.0)
      const double lat = double(lat_deg(icol)) * pi / 180.0;
      const double lon = double(lon_deg(icol)) * pi / 180.0;

      // Surface fields: uniform for the flat-surface baroclinic wave
      ps(icol)   = Scalar(p0);
      phis(icol) = Scalar(0.0);

      // Scalarize column slices of the Pack views so we can index by level.
      // ekat::scalarize converts uview_1d<Pack> (nlev_packs) to
      // uview_1d<Scalar> (nlev_packs * Pack::n), and we only access k < nlevs.
      auto T_col  = ekat::scalarize(ekat::subview(T_mid, icol));
      // horiz_winds is (ncols, 2, nlev_packs): subview on icol and component
      auto u_col  = ekat::scalarize(Kokkos::subview(horiz_winds, icol, 0, Kokkos::ALL()));
      auto v_col  = ekat::scalarize(Kokkos::subview(horiz_winds, icol, 1, Kokkos::ALL()));
      auto qv_col = ekat::scalarize(ekat::subview(qv, icol));
      auto qc_col = ekat::scalarize(ekat::subview(qc, icol));
      auto qr_col = ekat::scalarize(ekat::subview(qr, icol));

      // Sequential level loop.  Each call to wave_at_point runs the secant
      // iteration to convergence before moving to the next level.
      for (int k = 0; k < nlevs_; ++k) {
        // Mid-level pressure from hybrid coordinate formula.
        // ps = p0 everywhere (no orography), so p_mid = (hyam + hybm)*p0.
        const double p_mid = (double(hyam(k)) + double(hybm(k))) * p0;

        const auto s = wave_at_point(deep_, moist_, pertt_, X_d, a_d, Rd_d, Rvap_d,
                                      lon, lat, p_mid);

        T_col(k)  = s.T;
        u_col(k)  = s.u;
        v_col(k)  = s.v;
        qv_col(k) = s.qv;
        qc_col(k) = Scalar(0.0);  // No initial cloud liquid
        qr_col(k) = Scalar(0.0);  // No initial rain
      }
    });

  Kokkos::fence();
}

} // namespace dcmip2016
} // namespace scream

#endif // DCMIP2016_FUNCTIONS_IMPL_HPP
