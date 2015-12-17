#!/bin/sh

#qsub -q qexp -l select=1:ncpus=16:mpiprocs=1:cpu_freq=24,walltime=01:00:00 \
#     -l cpu_turbo_boost=0 -v NB_NODES=1,NB_PROCESS_PER_NODE=1,NB_CORES_PER_NODE=16 \
#     ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=1:ncpus=16:mpiprocs=1:cpu_freq=24,walltime=02:00:00\
#     -l cpu_turbo_boost=0 -v NB_NODES=1,NB_PROCESS_PER_NODE=1,NB_CORES_PER_NODE=16 \
#     ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=1:ncpus=16:mpiprocs=4:cpu_freq=24,walltime=01:00:00\
#     -l cpu_turbo_boost=0 -v NB_NODES=1,NB_PROCESS_PER_NODE=4,NB_CORES_PER_NODE=16 \
#     ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=1:ncpus=16:mpiprocs=8:cpu_freq=24,walltime=00:45:00\
#     -l cpu_turbo_boost=0 -v NB_NODES=1,NB_PROCESS_PER_NODE=8,NB_CORES_PER_NODE=16 \
#     ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=1:ncpus=16:mpiprocs=16:cpu_freq=24,walltime=00:30:00\
#     -l cpu_turbo_boost=0 -v NB_NODES=1,NB_PROCESS_PER_NODE=16,NB_CORES_PER_NODE=16 \
#     ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=4:ncpus=16:mpiprocs=1:cpu_freq=24,walltime=01:00:00\
#     -l cpu_turbo_boost=0 -v NB_NODES=4,NB_PROCESS_PER_NODE=1,NB_CORES_PER_NODE=16 \
#     ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=4:ncpus=16:mpiprocs=2:cpu_freq=24,walltime=00:45:00\
#     -l cpu_turbo_boost=0 -v NB_NODES=4,NB_PROCESS_PER_NODE=2,NB_CORES_PER_NODE=16 \
#     ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=4:ncpus=16:mpiprocs=16:cpu_freq=24,walltime=00:15:00\
#     -l cpu_turbo_boost=0 -v NB_NODES=4,NB_PROCESS_PER_NODE=16,NB_CORES_PER_NODE=16 \
#     ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=8:ncpus=16:mpiprocs=1:cpu_freq=24,walltime=00:45:00\
#     -l cpu_turbo_boost=0 -v NB_NODES=8,NB_PROCESS_PER_NODE=1,NB_CORES_PER_NODE=16 \
#     ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=8:ncpus=16:mpiprocs=2:cpu_freq=24,walltime=00:30:00\
#     -l cpu_turbo_boost=0 -v NB_NODES=8,NB_PROCESS_PER_NODE=2,NB_CORES_PER_NODE=16 \
#     ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=8:ncpus=16:mpiprocs=16:cpu_freq=24,walltime=00:15:00\
#     -l cpu_turbo_boost=0 -v NB_NODES=8,NB_PROCESS_PER_NODE=16,NB_CORES_PER_NODE=16 \
#     ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=16:ncpus=16:mpiprocs=1:cpu_freq=24,walltime=00:30:00\
#     -l cpu_turbo_boost=0 -v NB_NODES=16,NB_PROCESS_PER_NODE=16,NB_CORES_PER_NODE=16 \
#     ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=32:ncpus=16:mpiprocs=1:cpu_freq=24,walltime=00:30:00\
#     -l cpu_turbo_boost=0 -v NB_NODES=32,NB_PROCESS_PER_NODE=16,NB_CORES_PER_NODE=16 \
#     ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=64:ncpus=16:mpiprocs=1:cpu_freq=24,walltime=00:30:00\
#     -l cpu_turbo_boost=0 -v NB_NODES=64,NB_PROCESS_PER_NODE=16,NB_CORES_PER_NODE=16 \
#     ./measures_miniFEM.sh
