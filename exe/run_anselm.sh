#!/bin/sh

qsub -A IT4I-3-4 -q qprod -l select=1:ncpus=16:cpu_freq=24,walltime=01:00:00 \
     -v NB_NODES=1,MAX_CORES=16 ./mesures_miniApp_anselm.sh

#qsub -A IT4I-3-4 -q qprod -l select=4:ncpus=16:cpu_freq=24,walltime=00:30:00 \
#     -v NB_NODES=4,MAX_CORES=64 ./mesures_miniApp_anselm.sh
#
#qsub -A IT4I-3-4 -q qprod -l select=8:ncpus=16:cpu_freq=24,walltime=00:30:00 \
#     -v NB_NODES=8,MAX_CORES=128 ./mesures_miniApp_anselm.sh
#
#qsub -A IT4I-3-4 -q qprod -l select=12:ncpus=16:cpu_freq=24,walltime=00:30:00 \
#     -v NB_NODES=12,MAX_CORES=192 ./mesures_miniApp_anselm.sh
#
#qsub -A IT4I-3-4 -q qprod -l select=16:ncpus=16:cpu_freq=24,walltime=00:30:00 \
#     -v NB_NODES=16,MAX_CORES=256 ./mesures_miniApp_anselm.sh
#
#qsub -A IT4I-3-4 -q qprod -l select=32:ncpus=16:cpu_freq=24,walltime=00:30:00 \
#     -v NB_NODES=32,MAX_CORES=512 ./mesures_miniApp_anselm.sh
