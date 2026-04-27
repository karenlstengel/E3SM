#!/bin/bash
#PBS -N IC_eamxx
#PBS -A NTDD0004
#PBS -j oe
#PBS -k eod
#PBS -q develop
#PBS -l walltime=00:30:00
#PBS -l select=1:ncpus=128:mem=200GB
#PBS -M kstengel@ucar.edu
#PBS -m e

# -l select=1:ncpus=128:mem=200GB
# ./e3sm_kessler_short_cpu.sh
./eamxx_kessler_analyticIC_cpu.sh


# -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100
# ./e3sm_kessler_short_gpu.sh
