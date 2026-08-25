# analytic_conditions: idealized test-case infrastructure for EAMxx.
#
# Unlike physics parameterizations, processes here exist solely to set up
# idealized initial or boundary conditions for testing dynamical cores and
# physics schemes.  They are not intended for production climate simulations.
#
# Each subdirectory builds a library and sets a corresponding
# EAMXX_HAS_<NAME> compile definition.  The process is then registered with
# the AtmosphereProcessFactory in src/physics/register_physics.hpp so that
# it can be named in atm_procs_list in a YAML configuration file.

add_subdirectory(dcmip2016)
