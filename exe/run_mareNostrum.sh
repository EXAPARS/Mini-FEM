#!/bin/sh

#BSUB -W 01:30
#BSUB -n 8
#BSUB -R "span[ptile=1]"
#BSUB -x

#BSUB -oo output_%J.out
#BSUB -eo output_%J.err
#BSUB -J miniFEM

export NB_NODES=8
export NB_PROCESS_PER_NODE=1
export NB_CORES_PER_NODE=16

./measures_miniFEM.sh
