#!/bin/bash

# -l select=1:ncpus=128:mem=200GB

# qsub -N CPP_1_ne30_cpu  -A NTDD0004 -j oe -k eod -q main  -l walltime=00:30:00 -l select=1:ncpus=128:mem=200GB  -M kstengel@ucar.edu -m e kessler_claude_perf_1day_ne30_cpu.sh

# qsub -N CPP_1_ne30_gpu  -A NTDD0004 -j oe -k eod -q main  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e kessler_claude_perf_1day_ne30_gpu.sh

# qsub -N CPP_1_ne16_cpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=128:mem=200GB  -M kstengel@ucar.edu -m e kessler_claude_perf_1day_ne16_cpu.sh

qsub -N CPP_1_ne16_gpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e kessler_claude_perf_1day_ne16_gpu.sh

# qsub -N CPP_1_ne60_cpu  -A NTDD0004 -j oe -k eod -q main  -l walltime=00:30:00 -l select=1:ncpus=128:mem=200GB  -M kstengel@ucar.edu -m e kessler_claude_perf_1day_ne60_cpu.sh

# qsub -N CPP_1_ne60_gpu  -A NTDD0004 -j oe -k eod -q main  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e kessler_claude_perf_1day_ne60_gpu.sh

# # 5 day runs 
# qsub -N CPP_5_ne30_cpu  -A NTDD0004 -j oe -k eod -q main  -l walltime=00:30:00 -l select=1:ncpus=128:mem=200GB  -M kstengel@ucar.edu -m e kessler_claude_perf_5day_ne30_cpu.sh

# qsub -N CPP_5_ne30_gpu  -A NTDD0004 -j oe -k eod -q main  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e kessler_claude_perf_5day_ne30_gpu.sh

# qsub -N CPP_5_ne16_cpu  -A NTDD0004 -j oe -k eod -q main  -l walltime=00:30:00 -l select=1:ncpus=128:mem=200GB  -M kstengel@ucar.edu -m e kessler_claude_perf_5day_ne16_cpu.sh

# qsub -N CPP_5_ne16_gpu  -A NTDD0004 -j oe -k eod -q main  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e kessler_claude_perf_5day_ne16_gpu.sh

# qsub -N CPP_5_ne60_cpu  -A NTDD0004 -j oe -k eod -q main  -l walltime=00:30:00 -l select=1:ncpus=128:mem=200GB  -M kstengel@ucar.edu -m e kessler_claude_perf_5day_ne60_cpu.sh

# qsub -N CPP_5_ne60_gpu  -A NTDD0004 -j oe -k eod -q main  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e kessler_claude_perf_5day_ne60_gpu.sh

# 10 day runs for error checking
# -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100
# qsub -N CPP_10_cpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=128:mem=200GB  -M kstengel@ucar.edu -m e kessler_claude_perf_10_cpu.sh

# qsub -N CPP_10_gpu  -A NTDD0004 -j oe -k eod -q develop  -l walltime=00:30:00 -l select=1:ncpus=64:mem=480GB:ngpus=4:gpu_type=a100  -M kstengel@ucar.edu -m e kessler_claude_perf_10_gpu.sh
