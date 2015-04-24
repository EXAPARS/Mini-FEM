#!/bin/sh

#qsub -q qexp -l select=1:ncpus=4:mpiprocs=4:cpu_freq=24,walltime=01:00:00 \
#     -v NB_NODES=1,MAX_CORES=4 ./mesures_miniFEM_anselm.sh

#qsub -q qexp -l select=1:ncpus=16:mpiprocs=16:cpu_freq=24,walltime=01:00:00 \
#     -v NB_NODES=1,MAX_CORES=16 ./mesures_miniFEM_anselm.sh

#qsub -A DD-15-6 -q qprod -l select=1:ncpus=16:cpu_freq=24,walltime=01:00:00 \
#     -v NB_NODES=1,MAX_CORES=16 ./mesures_miniFEM_anselm.sh

qsub -A DD-15-6 -q qprod -l select=4:ncpus=16:cpu_freq=24,walltime=00:30:00 \
     -v NB_NODES=4,MAX_CORES=64 ./mesures_miniFEM_anselm.sh

#qsub -A DD-15-6 -q qprod -l select=8:ncpus=16:cpu_freq=24,walltime=00:30:00 \
#     -v NB_NODES=8,MAX_CORES=128 ./mesures_miniFEM_anselm.sh
#
#qsub -A DD-15-6 -q qprod -l select=12:ncpus=16:cpu_freq=24,walltime=00:30:00 \
#     -v NB_NODES=12,MAX_CORES=192 ./mesures_miniFEM_anselm.sh
#
#qsub -A DD-15-6 -q qprod -l select=16:ncpus=16:cpu_freq=24,walltime=00:30:00 \
#     -v NB_NODES=16,MAX_CORES=256 ./mesures_miniFEM_anselm.sh
#
#qsub -A DD-15-6 -q qprod -l select=32:ncpus=16:cpu_freq=24,walltime=00:30:00 \
#     -v NB_NODES=32,MAX_CORES=512 ./mesures_miniFEM_anselm.sh
