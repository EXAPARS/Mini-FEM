#!/bin/sh

#qsub -q qexp -l select=1:ncpus=24:mpiprocs=1,walltime=01:00:00 \
#     -l cpu_turbo_boost=0 -v NB_NODES=1,NB_PROCESS_PER_NODE=1,NB_CORES_PER_NODE=24 \
#     ./measures_miniFEM.sh

qsub -A DD-15-6 -q qprod -l select=4:ncpus=24:mpiprocs=2,walltime=01:00:00 \
     -l cpu_turbo_boost=0 -v NB_NODES=4,NB_PROCESS_PER_NODE=2,NB_CORES_PER_NODE=24 \
     ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=8:ncpus=24:mpiprocs=1,walltime=00:30:00 \
#     -l cpu_turbo_boost=0 -v NB_NODES=8,NB_PROCESS_PER_NODE=1,NB_CORES_PER_NODE=24 \
#     ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=16:ncpus=24:mpiprocs=16,walltime=00:30:00 \
#     -l cpu_turbo_boost=0 -v NB_NODES=16,NB_PROCESS_PER_NODE=16,NB_CORES_PER_NODE=24 \
#     ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=32:ncpus=24:mpiprocs=16,walltime=00:30:00 \
#     -l cpu_turbo_boost=0 -v NB_NODES=32,NB_PROCESS_PER_NODE=16,NB_CORES_PER_NODE=24 \
#     ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=64:ncpus=24:mpiprocs=16,walltime=00:30:00 \
#     -l cpu_turbo_boost=0 -v NB_NODES=64,NB_PROCESS_PER_NODE=16,NB_CORES_PER_NODE=24 \
#     ./measures_miniFEM.sh
