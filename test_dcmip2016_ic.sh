#!/bin/bash
set -e

date

user=kstengel
scratch=/glade/derecho/scratch/$user/E3SM

####################################################################
# Machine, compset, etc.
####################################################################
CCSMROOT=$scratch/E3SM
# Standalone SCREAM compset (all-stub lnd/ice/ocn/rof/glc/wav) so we don't
# need land/ocean input data and don't have to fight the coupler -- we only
# care about the analytic IC process's output at t=0.
COMPSET=F2000-SCREAM-SA
RESOLUTION=ne4_ne4          # pure-GLL physics grid: matches the IC process's
                             # "no remap needed" grid-selection path
DYCORE=theta-l_kokkos
MACH=derecho
MYCOMPILER=nvidia
QUEUE_NAME=main

CASE_NAME="dcmip2016_ic_test"
CASE_ROOT="$scratch/e3sm_test/${CASE_NAME}"
CASE_SCRIPTS_DIR=${CASE_ROOT}/case
CASE_BUILD_DIR=${CASE_ROOT}/build
CASE_RUN_DIR=${CASE_ROOT}/run
export NETCDF_PATH=$NETCDF

####################################################################
# Create a new case
####################################################################
rm -rf $CASE_ROOT

cd $CCSMROOT/cime/scripts

./create_newcase --case ${CASE_NAME} --output-root ${CASE_ROOT} --script-root ${CASE_SCRIPTS_DIR} \
               --handle-preexisting-dirs u --compset ${COMPSET} --res ${RESOLUTION} --machine ${MACH} \
               --compiler ${MYCOMPILER} --project NTDD0004 --walltime "00:20:00" --verbose -q ${QUEUE_NAME}

####################################################################
# Configure & Compile
####################################################################
cd $CASE_SCRIPTS_DIR

./xmlchange EXEROOT=${CASE_BUILD_DIR}
./xmlchange RUNDIR=${CASE_RUN_DIR}

./xmlchange DEBUG=TRUE
./xmlchange NTASKS=16
./xmlchange NTHRDS=1
./xmlchange ROOTPE=0

./case.setup

./xmlchange CAM_TARGET=$DYCORE

# --- Wire in the analytic IC process ---
# "kessler" from the README's example atm_procs_list isn't actually a
# registered EAMxx process yet, so it's dropped.
#
# homme is dropped too, deliberately: EAMxx's output manager only ever
# writes a stream after atm_process_group->run(dt) has already executed and
# m_current_ts has been advanced (AtmosphereDriver::run(dt) in
# atmosphere_driver.cpp calls the process group's run() before out_mgr.run());
# there is no reachable path that writes the output stream at the bare,
# unmodified t=0 state. So if homme were in the list, the "instant, 1 nstep"
# output below would actually capture the state AFTER homme's first
# dynamics step, not the pure analytic IC.
#
# With only dcmip2016_baroclinic_wave_ic in the list, its run_impl is a
# documented no-op (all work happens once in initialize_impl), so nothing
# touches the fields between initialization and the output write -- the
# snapshot written after "step 1" is bit-for-bit the analytic IC.
./atmchange eamxx::atm_procs_list=dcmip2016_baroclinic_wave_ic

# phis=0 is analytically correct here, but must still be listed so any grid
# instance the IC process doesn't cover also gets initialized (see
# analytic_conditions/README.md).
./atmchange initial_conditions::phis=0.0

./atmchange atm_log_level=debug

# --- Output: dump every field the IC process fills, at t=0 (nstep 1) ---
cat > ${CASE_SCRIPTS_DIR}/dcmip2016_ic_output.yaml << 'EOF'
%YAML 1.1
---
filename_prefix: dcmip2016_ic
averaging_type: instant
fields:
  physics_gll:
    field_names:
      - T_mid
      - horiz_winds
      - ps
      - phis
      - qv
      - qc
      - qr
output_control:
  frequency: 1
  frequency_units: nsteps
...
EOF
./atmchange output_yaml_files+=${CASE_SCRIPTS_DIR}/dcmip2016_ic_output.yaml

./case.build

####################################################################
# Run E3SM
####################################################################
cd $CASE_SCRIPTS_DIR

./xmlchange RUN_TYPE="startup"
./xmlchange RUN_STARTDATE='0001-01-01'
./xmlchange RESUBMIT='0'
./xmlchange CONTINUE_RUN='FALSE'
./xmlchange STOP_N='1',STOP_OPTION='nsteps'
./xmlchange JOB_WALLCLOCK_TIME='00:10:00'
./xmlchange JOB_QUEUE=$QUEUE_NAME
./xmlchange BUDGETS=TRUE

# case.submit left commented out on purpose -- run this build/verify pass
# first, then submit separately once confirmed.
# ./case.submit
echo "Build and case config complete. Case dir: ${CASE_SCRIPTS_DIR}"
