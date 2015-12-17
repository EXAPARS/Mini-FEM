Mini-FEM Proto-Application
==========================

Overview
--------

Mini-FEM is a proto-application reproducing the assembly step of FEM applications
working on 3D unstructured meshes.
Mini-FEM is parallelized at distributed memory level with MPI domain decomposition.
It can be in addition parallelized at shared memory level with either:
  - A mesh coloring approach based on the OpenMP runtime,
  - A Divide & Conquer (D&C) approach based on the Cilk Plus runtime,
  - A hybrid D&C + coloring approach also based on the Cilk Plus runtime.

How to compile
--------------

Mini-FEM requires CMake 2.8.10 or newer.
The build directory is "Mini-FEM/build".

A single command is required to compile a new binary:
    ./iMake $VERSION [$VECTOR_LENGTH] [tree]

- The $VERSION variable can be either:
    - "ref" to build a pure MPI version,
    - "coloring" to build a hybrid MPI+coloring version,
    - "DC" to build a hybrid MPI+D&C version,
    - or "dc-vec" to build a hybrid version using MPI, D&C, and coloring.

The 2 last options requires to have the DC-lib.
The path to the DC-lib can be set at the beginning of the iMake file.

- The $VECTOR_LENGTH variable must be specified when using the D&C Vec version.
  It can be either SSE, AVX or MIC depending on the target architecture.

- The tree option is used to create a new D&C tree and new permutation functions.
  If not specified, the application will try to read the existing tree and permutations.
  The created tree and permutations are stored in "Mini-FEM/data/$USE_CASE/DC_tree".
  They are associated to the code version (D&C or D&C Vec), to the partition size,
  and to the number of MPI processes. The files name can be read this way:
  $VERSION_$PARTITION_SIZE_$NB_PROCESS_$PROCESS_RANK

If you want to build binaries for Xeon Phi, you need to use the MIC option whatever
the code version.

How to execute
--------------

The execution directory is "Mini-FEM/exe".

Here is the command to execute a binary:
    mpirun -np $NB_PROCESS ./bin/$BINARY $USE_CASE $OPERATOR $NB_ITERATIONS

- The $BINARY variable refers to previously compiled binaries located in "Mini-FEM/bin".
- The $USE_CASE variable can be either "LM6" or "EIB".
    - LM6 is composed of around 27,500 nodes and 150,000 elements.
    - EIB is composed of around 1,000,000 nodes and 6,000,000 elements.
- The $OPERATOR variable can be either "ela" or "lap" corresponding resp. to elasticity
  and laplacian operators.
- The $NB_ITERATIONS variable corresponds to the number of iterations that you want to
  be executed.

To vary the number of threads when using the mesh coloring version, you need to use the
"OMP_NUM_THREADS" environment variable. When using one of the D&C version, you need to
use the "CILK_NWORKERS" environment variable.

For instance, in order to execute 10 iterations of the D&C version with 1 MPI process
and 4 Cilk threads using the EIB use case and the elasticity operator, the command is:
    export CILK_NWORKERS=4
    mpirun -np 1 ./bin/miniFEM_DC EIB ela 10

How to read the results
-----------------------

The proto-application output is composed of:
 - A summary of the parameters of execution.
 - The execution time of each prerequisite steps.
 - The execution time, in RDTSC cycles, of the matrix assembly and the
   preconditioner creation for each iteration.
 - A numerical checking.
