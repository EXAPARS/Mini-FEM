#!/bin/sh

#qsub -q qexp -l select=1:ncpus=24:mpiprocs=24,walltime=01:00:00 \
#     -l cpu_turbo_boost=0 -v NB_NODES=1,MAX_CORES=24 ./measures_miniFEM.sh

qsub -A DD-15-6 -q qprod -l select=1:ncpus=24:mpiprocs=24,walltime=01:00:00 \
     -l cpu_turbo_boost=0 -v NB_NODES=1,MAX_CORES=24 ./measures_miniFEM.sh

qsub -A DD-15-6 -q qprod -l select=4:ncpus=24:mpiprocs=96,walltime=00:30:00 \
     -l cpu_turbo_boost=0 -v NB_NODES=4,MAX_CORES=96 ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=8:ncpus=24:mpiprocs=1,walltime=00:30:00 \
#     -l cpu_turbo_boost=0 -v NB_NODES=8,MAX_CORES=192 ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=12:ncpus=24:mpiprocs=1,walltime=00:30:00 \
#     -l cpu_turbo_boost=0 -v NB_NODES=12,MAX_CORES=288 ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=16:ncpus=24:mpiprocs=1,walltime=00:30:00 \
#     -l cpu_turbo_boost=0 -v NB_NODES=16,MAX_CORES=384 ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=32:ncpus=24:mpiprocs=1,walltime=00:30:00 \
#     -l cpu_turbo_boost=0 -v NB_NODES=32,MAX_CORES=768 ./measures_miniFEM.sh

#qsub -A DD-15-6 -q qprod -l select=64:ncpus=24:mpiprocs=1,walltime=01:00:00 \
#     -l cpu_turbo_boost=0 -v NB_NODES=64,MAX_CORES=1536 ./measures_miniFEM.sh
