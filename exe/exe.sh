#!/bin/sh
cd /home/lthebault/Dassault/Mini-FEM/exe
export PATH=$PATH:/home/lthebault/Programs/GPI-2/bin:/home/lthebault/Programs/CilkTools/bin
export LD_LIBRARY_PATH=/opt/intel/impi/4.1.2.040/intel64/lib:/opt/intel/composer_xe_2013_sp1.3.174/compiler/lib/intel64:/opt/intel/composer_xe_2013_sp1.3.174/mpirt/lib/intel64:/opt/intel/composer_xe_2013_sp1.3.174/compiler/lib/intel64:/opt/intel/composer_xe_2013_sp1.3.174/mkl/lib/intel64:/usr/local/lib
export CILK_NWORKERS=4
#export DISPLAY=localhost:10.0
#./bin/miniFEM_REF_BulkSynchronous_GASPI_CILK LM6 ela 15
#./bin/miniFEM_DC_BulkSynchronous_TreeCreation_GASPI_CILK LM6 ela 15
#./bin/miniFEM_DC_BulkSynchronous_GASPI_CILK LM6 ela 15
#./bin/miniFEM_DC_MultithreadedComm_TreeCreation_GASPI_CILK LM6 ela 5
./bin/miniFEM_DC_MultithreadedComm_GASPI_CILK FGN1 ela 5
#xterm -e gdb --ex run --args ./bin/miniFEM_DC_MultithreadedComm_TreeCreation_Debug_GASPI_CILK EIB ela 1
#xterm -e gdb --ex run --args ./bin/miniFEM_DC_MultithreadedComm_Debug_GASPI_CILK EIB ela 1
