#!/bin/bash

date

user=kstengel
scratch=/glade/derecho/scratch/$user/E3SM
# scratch=/glade/campaign/cisl/asap/$user/

####################################################################
# Machine, compset, etc.
####################################################################
CCSMROOT=${scratch}/E3SM
# CCSMROOT=/glade/derecho/scratch/$user/E3SM/E3SM
COMPSET=F2000-SCREAMv1-KESSLER
RESOLUTION=ne30_ne30 #ne30pg2_ne30pg2,ne4pg2_ne4pg2
DYCORE=theta-l_kokkos
MACH=derecho
MYCOMPILER=nvidiagpu
QUEUE_NAME=main

# CASE_NAME="${COMPSET}.${RESOLUTION}.${MACH}.${MYCOMPILER}.${DYCORE}"
CASE_NAME="ne30np4_5day_gpuOA"
CASE_ROOT="$scratch/e3sm_test/JAX_v_Fortran_perf/${CASE_NAME}"
CASE_SCRIPTS_DIR=${CASE_ROOT}/case
CASE_BUILD_DIR=${CASE_ROOT}/build
CASE_RUN_DIR=${CASE_ROOT}/run
CASE_ARCHIVE_DIR=${CASE_ROOT}/archive
export NETCDF_PATH=$NETCDF
export KESSLER_PERF_LOG_PATH=${CASE_SCRIPTS_DIR}/kessler_perf_log.csv

####################################################################
# Create a new case
####################################################################
rm -rf $CASE_ROOT

cd $CCSMROOT/cime/scripts

./create_newcase --case ${CASE_NAME} --output-root ${CASE_ROOT} --script-root ${CASE_SCRIPTS_DIR} \
               --handle-preexisting-dirs u --compset ${COMPSET} --res ${RESOLUTION} --machine ${MACH} \
               --compiler ${MYCOMPILER} --project NTDD0004 --walltime "00:59:00" --verbose -q ${QUEUE_NAME} \
               --user-mods-dir ${CCSMROOT}/components/eamxx//cime_config/testdefs/testmods_dirs/eamxx/L58-kessler

# ${CCSMROOT}/components/eamxx//cime_config/testdefs/testmods_dirs/eamxx/output/preset/2
####################################################################
# Configure & Compile
####################################################################
cd $CASE_SCRIPTS_DIR

./xmlchange EXEROOT=${CASE_BUILD_DIR}
./xmlchange RUNDIR=${CASE_RUN_DIR}

./xmlchange DEBUG=FALSE

./xmlchange NTASKS=4
./xmlchange NTHRDS=1
./xmlchange NGPUS_PER_NODE=4
./xmlchange GPU_TYPE=a100 # NVIDIA A100 GPUs in Derecho
./xmlchange OPENACC_GPU_OFFLOAD=TRUE # TRUE for with OpenACC
./xmlchange USE_OPENACC_BACKEND=TRUE
./xmlchange OPENMP_GPU_OFFLOAD=FALSE
./xmlchange KOKKOS_GPU_OFFLOAD=TRUE
./xmlchange OVERSUBSCRIBE_GPU=FALSE
./xmlchange ROOTPE='0'
./xmlchange DOUT_S=false

./case.setup

./xmlchange CAM_TARGET=$DYCORE
./xmlchange GMAKE_J='32'

./xmlchange ATM_NCPL=48 # 30 min time step, 48 time steps per day, daily output

./atmchange atm_log_level=info #debug
# ./atmchange physics::atm_procs_list=mac_aero_mic # this removes the rrtmgp physics
# ./atmchange mac_aero_mic::atm_procs_list=kessler #kessler
./atmchange save_field_manager_content=true
./atmchange output_yaml_files+=/glade/derecho/scratch/kstengel/E3SM/E3SM/output_control_JAX.yml
./atmchange initial_conditions::filename=/glade/derecho/scratch/kstengel/inputdata/atm/scream/init/FKESSLER_NE30NP4.cam.i.moist_baroclinic_wave_dcmip2016.nc
./atmchange enable_fine_grain_timers=false

# use below to match to stormspeed
./atmchange ctl_nl::dt_tracer_factor=6
./atmchange ctl_nl::hypervis_subcycle_q=6
./atmchange ctl_nl::se_ftype=2
./atmchange ctl_nl::se_nsplit=2
./atmchange ctl_nl::statefreq=488
./atmchange ctl_nl::transport_alg=12


# -------------------------------------------
./atmquery --listall
./case.build

# ####################################################################
# Run E3SM
# ####################################################################
cd $CASE_SCRIPTS_DIR

./xmlchange RUN_TYPE="startup"
if [[ $COMPSET == *"F20TR"* ]]; then
   ./xmlchange RUN_STARTDATE='1850-01-01'
elif [[ $COMPSET == "FMTHIST" || $COMPSET == "FLTHIST" ]]; then
   ./xmlchange RUN_STARTDATE='2001-01-01'
else
   ./xmlchange RUN_STARTDATE='0001-01-01'
fi
./xmlchange RESUBMIT='0'
./xmlchange CONTINUE_RUN='FALSE'
./xmlchange STOP_N='5',STOP_OPTION='ndays' # note that we need to run this for 10 days to see anything interesting
./xmlchange JOB_WALLCLOCK_TIME='07:00:00'
./xmlchange JOB_QUEUE=$QUEUE_NAME
./xmlchange BUDGETS=TRUE

if [[ $DYCORE == "theta-l_kokkos" ]]; then
cat << EOF >> user_nl_elm
   check_finidat_year_consistency = .false.
   check_dynpft_consistency = .false.
   create_crop_landunit = .false.
EOF
fi

./case.submit
