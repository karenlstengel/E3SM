#!/bin/bash

# -l select=1:ncpus=128:mem=200GB
qsub -N perf_1_cpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=128:mem=200GB  -M kstengel@ucar.edu -m e kessler_perf_1_cpu.sh

qsub -N perf_1_gpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e kessler_perf_1_gpu.sh

qsub -N perf_1_gpucpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e kessler_perf_1_gpucpu.sh


qsub -N perf_5_cpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=128:mem=200GB  -M kstengel@ucar.edu -m e kessler_perf_5_cpu.sh

qsub -N perf_5_gpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e kessler_perf_5_gpu.sh

qsub -N perf_5_gpucpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e kessler_perf_5_gpucpu.sh


# -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100
# qsub -N perf_10_cpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=128:mem=200GB  -M kstengel@ucar.edu -m e kessler_perf_10_cpu.sh

# qsub -N perf_10_gpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e kessler_perf_10_gpu.sh

# qsub -N perf_10_gpucpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e kessler_perf_10_gpucpu.sh

# OpenACC backend runs
# qsub -N perf_1_gpuOA  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e kessler_perf_1_gpuOA.sh

# qsub -N perf_10_gpuOA  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e kessler_perf_10_gpuOA.sh
