#!/bin/bash

# -l select=1:ncpus=128:mem=200GB
qsub -N claude_1_cpu  -A NTDD0004 -j oe -k eod -q main  -l walltime=00:30:00 -l select=1:ncpus=128:mem=200GB  -M kstengel@ucar.edu -m e kessler_claude_perf_1_cpu.sh

qsub -N claude_1_gpu  -A NTDD0004 -j oe -k eod -q main  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e kessler_claude_perf_1_gpu.sh

# -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100
qsub -N claude_10_cpu  -A NTDD0004 -j oe -k eod -q main  -l walltime=00:30:00 -l select=1:ncpus=128:mem=200GB  -M kstengel@ucar.edu -m e kessler_claude_perf_10_cpu.sh

qsub -N claude_10_gpu  -A NTDD0004 -j oe -k eod -q main  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e kessler_claude_perf_10_gpu.sh
