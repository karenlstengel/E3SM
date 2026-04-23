#!/bin/bash
#PBS -N hommeOnly
#PBS -A NTDD0004
#PBS -j oe
#PBS -k eod
#PBS -q develop
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100
#PBS -M kstengel@ucar.edu
#PBS -m e

# -l select=1:ncpus=128:mem=200GB
# ./e3sm_short_stengel_fortran_nvidiacpu.sh


# -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100
./e3sm_short_stengel_fortran_nvidiagpu.sh