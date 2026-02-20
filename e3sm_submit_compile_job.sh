#!/bin/bash
#PBS -N eamxx_compile
#PBS -A NTDD0004
#PBS -j oe
#PBS -k eod
#PBS -q main
#PBS -l walltime=00:20:00
#PBS -l select=1:ncpus=128:mem=200GB
#PBS -M kstengel@ucar.edu
#PBS -m e

./e3sm_kessler_short_cpu.sh