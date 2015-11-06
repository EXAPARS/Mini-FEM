#!/bin/sh
cd /home/lthebault/Dassault/Mini-FEM/exe
export PATH=$PATH:/home/lthebault/Programs/GPI-2/bin:/home/lthebault/Programs/CilkTools/bin
export LD_LIBRARY_PATH=/opt/intel/impi/4.1.2.040/intel64/lib:/opt/intel/composer_xe_2013_sp1.3.174/compiler/lib/intel64:/opt/intel/composer_xe_2013_sp1.3.174/mpirt/lib/intel64:/opt/intel/composer_xe_2013_sp1.3.174/compiler/lib/intel64:/opt/intel/composer_xe_2013_sp1.3.174/mkl/lib/intel64:/usr/local/lib
export CILK_NWORKERS=1
export DISPLAY=localhost:10.0
#./bin/miniFEM_DC_GASPI_CILK_BulkSynchronous_TreeCreation LM6 lap 1
./bin/miniFEM_DC_GASPI_CILK_MultithreadedComm_TreeCreation_Debug LM6 lap 1
#xterm -e gdb --ex run --args ./bin/miniFEM_DC_GASPI_CILK_MultithreadedComm_TreeCreation_Debug LM6 lap 1
