#ifndef DCMIP2016_FUNCTIONS_HPP
#define DCMIP2016_FUNCTIONS_HPP

// Declarations for the DCMIP2016 Test 1 analytic initial-condition functions.
//
// Original Fortran: components/homme/src/test_src/dcmip2016-baroclinic.F90
// Reference: Ullrich, Melvin, Staniforth & Jablonowski (2015), QJRMS 141(686):
//   "A proposed baroclinic wave test case for deep and shallow-atmosphere
//    dynamical cores", doi:10.1002/qj.2544
//
// Design follows the EAMxx physics convention:
//   - This header contains only the struct layout and function declarations.
//   - Template function definitions live in dcmip2016_functions_impl.hpp.
//   - Explicit instantiation for <Real, DefaultDevice> is in dcmip2016.cpp.
//   - The process interface calls BaroclinicWaveFunctions::main(), which
//     launches a Kokkos kernel on whatever device EAMxx is configured for.

#include "share/core/eamxx_types.hpp"

#include <ekat_pack.hpp>
#include <ekat_pack_kokkos.hpp>

namespace scream {
namespace dcmip2016 {

// ==========================================================================
// BaroclinicWaveFunctions<ScalarT, DeviceT>
//
// ScalarT : the model's floating-point scalar type (Real = double or float)
// DeviceT : the Kokkos device (DefaultDevice on GPU or CPU builds)
//
// All computation-heavy functions are KOKKOS_FUNCTION so they run on the
// same device as the rest of EAMxx.  The sequential secant iteration inside
// eval_z_from_p is valid inside a KOKKOS_FUNCTION because each Kokkos
// thread runs it independently (one thread per column).
// ==========================================================================
template <typename ScalarT, typename DeviceT>
struct BaroclinicWaveFunctions
{
  // ---- Type aliases (follows the EAMxx physics convention) ----
  using Scalar = ScalarT;
  using Device = DeviceT;

  // Pack: SIMD type used for level-vectorized operations in EAMxx.
  // Note: the analytic IC formulas are scalar per level (no level
  // vectorization is possible because the secant iteration is sequential),
  // but we accept Pack views in main() so the interface matches EAMxx
  // field-manager conventions.  We scalarize inside the Kokkos kernel.
  using Pack = ekat::Pack<Scalar, SCREAM_PACK_SIZE>;

  using KT          = KokkosTypes<Device>;
  using MemberType  = typename KT::MemberType;

  template <typename S>
  using view_1d = typename KT::template view_1d<S>;
  template <typename S>
  using view_2d = typename KT::template view_2d<S>;
  template <typename S>
  using view_3d = typename KT::template view_3d<S>;
  template <typename S>
  using uview_1d = typename ekat::template Unmanaged<view_1d<S>>;

  // ==========================================================================
  // Params: physical and test-case constants.
  //
  // Using double (not ScalarT) for the constants so that the analytic formula
  // intermediate computations are always performed in double precision.
  // This ensures the secant iteration converges regardless of whether
  // ScalarT is float or double.  Results are cast back to Scalar on output.
  //
  // Values match dcmip2016-baroclinic.F90 and EAMxx's physics_constants.hpp.
  // ==========================================================================
  struct Params {
    // Mathematical
    static constexpr double pi      = 3.14159265358979323846;

    // Earth reference values
    static constexpr double a       = 6.376e6;     // Radius (m)
    static constexpr double g       = 9.80616;     // Gravity (m/s^2)
    static constexpr double Rd      = 287.042;     // Dry-air gas constant (J/kg/K)
    static constexpr double cp      = 1004.64;     // Specific heat Cp (J/kg/K)
    static constexpr double Rvap    = 461.505;     // Water-vapor gas constant (J/kg/K)
    static constexpr double omega   = 7.292e-5;    // Earth's rotation rate (rad/s)
    static constexpr double p0      = 100000.0;    // Reference surface pressure (Pa)
    static constexpr double kappa   = Rd / cp;     // Poisson exponent

    // Mvap used in virtual-temperature correction: T_v = T*(1 + Mvap*qv)
    static constexpr double Mvap    = Rvap / Rd - 1.0;   // ≈ 0.608

    // Background-state parameters (Ullrich et al. 2015, Table 1)
    static constexpr double T0E     = 310.0;   // Equatorial surface temperature (K)
    static constexpr double T0P     = 240.0;   // Polar surface temperature (K)
    static constexpr double B       = 2.0;     // Jet half-width parameter
    static constexpr double K       = 3.0;     // Jet power-law exponent
    static constexpr double lapse   = 0.005;   // Lapse-rate parameter (K/m)

    // Exponential perturbation (pertt = 0): compact Gaussian in zonal wind
    static constexpr double pertup    = 1.0;              // Max speed (m/s)
    static constexpr double pertexpr  = 0.1;              // Radius (non-dim)
    static constexpr double pertlon   = pi / 9.0;         // Center longitude (rad) ~20°
    static constexpr double pertlat   = 2.0 * pi / 9.0;  // Center latitude  (rad) ~40°
    static constexpr double pertz     = 15000.0;          // Vertical taper height (m)

    // Stream-function perturbation (pertt = 1): cosine-bell
    static constexpr double pertu0    = 0.5;       // Max speed (m/s)
    static constexpr double pertr     = 1.0 / 6.0; // Radius (non-dim)

    // Finite-difference step for numerical gradient of stream function
    static constexpr double dxepsilon = 1.e-5;

    // Moisture parameters
    static constexpr double moistqlat  = 2.0 * pi / 9.0; // Latitudinal width (rad)
    static constexpr double moistqp    = 34000.0;         // Vertical pressure width (Pa)
    static constexpr double moisttr    = 0.1;             // Cutoff pressure ratio eta
    static constexpr double moistqs    = 1.e-12;          // Specific humidity above cutoff
    static constexpr double moistq0    = 0.018;           // Max surface specific humidity
  };

  // ==========================================================================
  // State: atmospheric state at a single grid point.
  // Returned by wave_at_point() and stored into field views by main().
  // ==========================================================================
  struct State {
    Scalar u;     // Zonal wind (m/s)
    Scalar v;     // Meridional wind (m/s)
    Scalar T;     // Temperature (K), virtual-T corrected when moist
    Scalar ps;    // Surface pressure (Pa) — uniform p0 (no orography)
    Scalar phis;  // Surface geopotential (m^2/s^2) — always 0
    Scalar qv;    // Water-vapor specific humidity (kg/kg); 0 if dry
  };

  // ==========================================================================
  // Function declarations
  // All KOKKOS_FUNCTION methods are implemented in dcmip2016_functions_impl.hpp
  // and explicitly instantiated in dcmip2016.cpp.
  // ==========================================================================

  // Given altitude z (m), return hydrostatic pressure p (Pa) and temperature
  // T (K) at position (lon, lat) in radians.  Purely algebraic — GPU-safe.
  KOKKOS_FUNCTION
  static void eval_pressure_temperature(
      int deep, double X, double lon, double lat, double z,
      double& p, double& T);

  // Given pressure p_tgt (Pa), find altitude z (m) and temperature T (K) via
  // a secant-method iteration on eval_pressure_temperature.
  // The sequential loop is valid inside a KOKKOS_FUNCTION because each Kokkos
  // thread runs it independently (one thread per column in main()).
  KOKKOS_FUNCTION
  static void eval_z_from_p(
      int deep, double X, double lon, double lat, double p_tgt,
      double& z, double& T);

  // Zonal-wind perturbation for the exponential type (pertt = 0).  GPU-safe.
  KOKKOS_FUNCTION
  static double eval_exponential(double lon, double lat, double z);

  // Stream-function value for the cosine-bell perturbation (pertt = 1).  GPU-safe.
  KOKKOS_FUNCTION
  static double eval_streamfunction(double lon, double lat, double z);

  // Compute the complete DCMIP2016 Test 1 state at one point.
  // Calls eval_z_from_p — valid on GPU (sequential loop per thread).
  KOKKOS_FUNCTION
  static State wave_at_point(
      int deep, int moist, int pertt, double X,
      double lon, double lat, double p);

  // Kernel launcher: fills all output fields for every column and level.
  // Launched on DeviceT via a Kokkos RangePolicy (one thread per column).
  static void main(
      int ncols, int nlevs,
      int deep, int moist, int pertt, Scalar X,
      const view_1d<const Scalar>& lat_deg,   // geometry, degrees
      const view_1d<const Scalar>& lon_deg,   // geometry, degrees
      const view_1d<const Scalar>& hyam,      // hybrid-pressure A coefficients
      const view_1d<const Scalar>& hybm,      // hybrid-pressure B coefficients
      const view_2d<Pack>& T_mid,             // (ncols, nlev_packs) temperature
      const view_3d<Pack>& horiz_winds,       // (ncols, 2, nlev_packs) u=0, v=1
      const view_1d<Scalar>& ps,              // (ncols) surface pressure
      const view_1d<Scalar>& phis,            // (ncols) surface geopotential
      const view_2d<Pack>& qv,               // (ncols, nlev_packs) water vapor
      const view_2d<Pack>& qc,               // (ncols, nlev_packs) cloud liquid
      const view_2d<Pack>& qr);              // (ncols, nlev_packs) rain water

}; // struct BaroclinicWaveFunctions

} // namespace dcmip2016
} // namespace scream

#endif // DCMIP2016_FUNCTIONS_HPP
