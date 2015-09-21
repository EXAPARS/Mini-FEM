#!/bin/sh
cd /home/lthebaul/lthebaul/Mini-FEM/exe
export LD_LIBRARY_PATH=/apps/all/ncurses/5.9-intel-2015b/lib:/apps/all/imkl/11.2.3.187-iimpi-7.3.5-GNU-5.1.0-2.25/mkl/lib/intel64:/apps/all/imkl/11.2.3.187-iimpi-7.3.5-GNU-5.1.0-2.25/lib/intel64:/apps/all/impi/5.0.3.048-iccifort-2015.3.187-GNU-5.1.0-2.25/lib64:/apps/all/ifort/2015.3.187-GNU-5.1.0-2.25/lib/intel64:/apps/all/ifort/2015.3.187-GNU-5.1.0-2.25/lib:/apps/all/icc/2015.3.187-GNU-5.1.0-2.25/lib/intel64:/apps/all/icc/2015.3.187-GNU-5.1.0-2.25/lib:/apps/all/binutils/2.25-GCC-5.1.0-binutils-2.25/lib:/apps/all/GCC/5.1.0-binutils-2.25/lib/gcc/x86_64-unknown-linux-gnu/5.1.0:/apps/all/GCC/5.1.0-binutils-2.25/lib64:/apps/all/GCC/5.1.0-binutils-2.25/lib
export CILK_NWORKERS=1
./bin/miniFEM_DC_GASPI_CILK_MultithreadedComm_TreeCreation LM6 lap 2
