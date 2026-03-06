#!/bin/bash

date

user=kstengel
scratch=/glade/derecho/scratch/$user/E3SM

####################################################################
# Machine, compset, etc.
####################################################################
CCSMROOT=$scratch/E3SM
COMPSET=F2000-SCREAMv1-AQP1 #F20TR-SCREAMv1
RESOLUTION=ne4_ne4
DYCORE=theta-l_kokkos
MACH=derecho
MYCOMPILER=nvidiagpu
QUEUE_NAME=main

# CASE_NAME="${COMPSET}.${RESOLUTION}.${MACH}.${MYCOMPILER}.${DYCORE}"
CASE_NAME="kessler_ne4_short_eamxx_gpu_test"
CASE_ROOT="$scratch/e3sm_test/${CASE_NAME}"
CASE_SCRIPTS_DIR=${CASE_ROOT}/case
CASE_BUILD_DIR=${CASE_ROOT}/build
CASE_RUN_DIR=${CASE_ROOT}/run
CASE_ARCHIVE_DIR=${CASE_ROOT}/archive
export NETCDF_PATH=$NETCDF

####################################################################
# Create a new case 
####################################################################
rm -rf $CASE_ROOT

cd $CCSMROOT/cime/scripts

./create_newcase --case ${CASE_NAME} --output-root ${CASE_ROOT} --script-root ${CASE_SCRIPTS_DIR} \
               --handle-preexisting-dirs u --compset ${COMPSET} --res ${RESOLUTION} --machine ${MACH} \
               --compiler ${MYCOMPILER} --project NTDD0004 --walltime "00:59:00" --verbose -q ${QUEUE_NAME} \
               --user-mods-dir ${CCSMROOT}/components/eamxx//cime_config/testdefs/testmods_dirs/eamxx/output/preset/2 ${CCSMROOT}/components/eamxx//cime_config/testdefs/testmods_dirs/eamxx/L72 

####################################################################
# Configure & Compile
####################################################################
cd $CASE_SCRIPTS_DIR

./xmlchange EXEROOT=${CASE_BUILD_DIR}
./xmlchange RUNDIR=${CASE_RUN_DIR}

./xmlchange DEBUG=TRUE

./xmlchange NTASKS=4
./xmlchange NTHRDS=1
./xmlchange NGPUS_PER_NODE=4
./xmlchange GPU_TYPE=a100 # NVIDIA A100 GPUs in Derecho
./xmlchange OPENACC_GPU_OFFLOAD=FALSE
./xmlchange OPENMP_GPU_OFFLOAD=FALSE
./xmlchange KOKKOS_GPU_OFFLOAD=TRUE
./xmlchange OVERSUBSCRIBE_GPU=FALSE
./xmlchange ROOTPE='0'
./xmlchange DOUT_S=false
# ./xmlchange DOUT_S_ROOT=${CASE_ARCHIVE_DIR}

./case.setup

./xmlchange CAM_TARGET=$DYCORE
./xmlchange GMAKE_J='32'

# Turn off all other pyhsics except Kessler
./atmchange atm_log_level=debug
./atmchange physics::atm_procs_list=mac_aero_mic # this removes the rrtmgp physics
./atmchange mac_aero_mic::atm_procs_list=kessler #kessler
./atmchange save_field_manager_content=true
./atmchange output_yaml_files+=/glade/derecho/scratch/kstengel/E3SM/E3SM/output_control.yml
   
./case.build 
#####################################################################
# Run E3SM
#####################################################################
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
./xmlchange STOP_N='15',STOP_OPTION='ndays'
./xmlchange JOB_WALLCLOCK_TIME='00:05:00'
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
