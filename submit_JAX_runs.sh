#!/bin/bash

# 1 and 5 day runs to get performance data
qsub -N JAX_1_cpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=128:mem=200GB  -M kstengel@ucar.edu -m e JAX_kessler_perf_1_cpu.sh
qsub -N JAX_1_cpu_c  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=128:mem=200GB  -M kstengel@ucar.edu -m e JAX_kessler_perf_1_cpu_comp.sh

qsub -N JAX_1_gpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e JAX_kessler_perf_1_gpu.sh
qsub -N JAX_1_gpu_c  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e JAX_kessler_perf_1_gpu_comp.sh

qsub -N JAX_1_gpucpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e JAX_kessler_perf_1_gpucpu.sh
qsub -N JAX_1_gpucpu_c  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e JAX_kessler_perf_1_gpucpu_comp.sh


qsub -N JAX_5_cpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=128:mem=200GB  -M kstengel@ucar.edu -m e JAX_kessler_perf_5_cpu.sh
qsub -N JAX_5_cpu_c  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=128:mem=200GB  -M kstengel@ucar.edu -m e JAX_kessler_perf_5_cpu_comp.sh

qsub -N JAX_5_gpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e JAX_kessler_perf_5_gpu.sh
qsub -N JAX_5_gpu_c  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e JAX_kessler_perf_5_gpu_comp.sh

qsub -N JAX_5_gpucpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e JAX_kessler_perf_5_gpucpu.sh
qsub -N JAX_5_gpucpu_c  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e JAX_kessler_perf_5_gpucpu_comp.sh

# 10 day runs to error data 

qsub -N JAX_10_cpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=128:mem=200GB  -M kstengel@ucar.edu -m e JAX_kessler_perf_10_cpu.sh

qsub -N JAX_10_gpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e JAX_kessler_perf_10_gpu.sh
