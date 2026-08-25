#ifndef SCREAM_REGISTER_ANALYTIC_CONDITIONS_PROCESSES_HPP
#define SCREAM_REGISTER_ANALYTIC_CONDITIONS_PROCESSES_HPP

#include "share/atm_process/atmosphere_process.hpp"

// Only include headers and register processes for libs that have been linked in

#ifdef EAMXX_HAS_DCMIP2016_IC
#include "analytic_conditions/dcmip2016/eamxx_dcmip2016_ic_process_interface.hpp"
#endif

namespace scream {

inline void register_analytic_conditions () {
  auto& proc_factory = AtmosphereProcessFactory::instance();
#ifdef EAMXX_HAS_DCMIP2016_IC
  // Register the DCMIP2016 Test 1 analytic IC process.
  // Activated when the dcmip2016_ic CMake target is linked (sets EAMXX_HAS_DCMIP2016_IC).
  proc_factory.register_product("dcmip2016_baroclinic_wave_ic",
    &create_atmosphere_process<DCMIP2016BaroclinicIC>);
#endif

  // If no physics was enabled, silence compile warning about unused var
  (void) proc_factory;
}

} // namespace scream

#endif // SCREAM_REGISTER_ANALYTIC_CONDITIONS_PROCESSES_HPP
