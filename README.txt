PPoPP'15 Artifact Evaluation
============================

Archive contents
================
 - The current README.txt
 - The accepted version of the paper: "accepted_paper.pdf".
 - The last version of the paper: "camera_ready_draft.pdf" which includes additional
   experiments on 4 Xeon Phi (KNC).
 - An artifact directory containing our proto-application, a CDE package, and an access
   to a test machine.


Confidentiality caution
=======================
The proto-application shared in the artifact directory is representative of an
industrial application coming from Dassault Aviation. This proto-application is under
licensing and it still contains Dassault Aviation property code.
It must not be shared outside the committee.

Among the five versions presented in the publication, four versions (Ref Improved,
Coloring, D&C, and D&C Hybrid) come from our proto-application and can be freely
experimented by the committee members. However, the Ref version coming from Dassault
Aviation is not in the present artifact package.

Moreover, the FGN mesh used for the 4 KNC experiment is a Dassault Aviation proprietary
mesh and it is also not in this package. Only the small LM6 mesh used for the locality
experiments and the EIB mesh used for the SandyBridge and the single KNC experiments.


Artifact experiments
====================
We provide 3 different ways to evaluate our artifact:
 (1) - An access to a 4 cores E31245 Xeon machine containing the proto-application with
       all its dependancies preconfigured and the EIB and LM6 test cases.
 (2) - The proto-application itself containing the sources and the build process.
 (3) - A CDE package containing a binary of the four versions: Ref Improved, Coloring,
       D&C, and D&C Hybrid. 

We do not provide a VM since it would not be representative of the performance of the
different versions.

How to use the option (1)
-------------------------
Please connect via SSH to "majorque.prism.uvsq.fr" using the following login and password:
Login: ppopp_aec
PWD: 2015_ppopp_aec_pwd

A binary of each of the four available code versions (Ref Improved, Coloring, D&C, and
D&C Hybrid) has already been built.
You can directly execute them in the "MiniFEM/exe" directory.

Here is the command to execute a binary:
    mpirun -np $NB_PROCESS ./bin/$BINARY $TEST_CASE $OPERATOR $NB_ITERATIONS

- The $BINARY variable can be one of the four binaries in the ./bin directory.
- The $TEST_CASE variable can be either "LM6" or "EIB".
    - LM6 is composed of around 27,500 nodes and 150,000 elements.
    - EIB is composed of around 1,000,000 nodes and 6,000,000 elements.
- The $OPERATOR variable can be either "ela" or "lap" corresponding resp. to elasticity
  and laplacian operators. In the publication we only used the elasticity operator.
- The $NB_ITERATIONS variable corresponds to the number of iterations that you want to
  be executed.

To vary the number of Cilk threads, use the "CILK_NWORKERS" environment variable.

For instance, to execute 10 iterations of the D&C version with 1 MPI process and 4 Cilk
threads using the EIB test case and the elasticity operator, the commands would be:
    export CILK_NWORKERS=4
    mpirun -np 1 ./bin/miniFEM_DC EIB ela 10

How to use the option (2)
-------------------------
If you want to compile yourself new binaries of the proto-application, you can access
the "majorque.prism.uvsq.fr" machine (see previous section for further details).

The build directory is "MiniFEM/build".
Here, you need to use the following command:
    ./iMake $VERSION [$VECTOR_LENGTH] [tree]

- The $VERSION variable can be either:
    - "ref" to build the Ref Improved version,
    - "coloring" to build the coloring version,
    - "DC" to build the D&C version,
    - or "hybrid" to build the D&C Hybrid version.
- The $VECTOR_LENGTH variable must be specified with the D&C Hybrid version.
  It can be either SSE, AVX or MIC depending on the target architecture.
- The tree option is used to create a new D&C tree and new permutation functions.
  If not specified, the application will try to read the existing tree and permutations.
  The created tree and permutations are stored in the "MiniFEM/data/$TEST_CASE/DC_tree"
  directory. They are associated to the code version (D&C or D&C Hybrid), to the
  partition size, and to the number of MPI processes.
  The files name can be read this way: $VERSION_$PARTITION_SIZE_$NB_PROCESS_$PROCESS_RANK.

If you plan to build binaries for Xeon Phi (whatever the code version), you need to use
the MIC option too.

You can also vary the D&C partitions size by modifying the "DC.h" header file. It is
located in "DC/include/". The value corresponds to the "MAX_ELEM_PER_PART" define
located at line 8.

Once you have build your new binary, you can execute it in the "MiniFEM/exe" directory
as explained in the previous section.

How to use the option (3)
-------------------------
To use the CDE packages, once your are in the given artifact directory, please go to:
"MiniFEM/exe/cde-package/cde-root/home/lthebault/Dassault/MiniFEM/exe"

Here you can find the CDE scripts to execute each of the four available code versions.
You can execute the scripts with the following command:
    mpirun -np $NB_PROCESS ./$SCRIPT.cde $TEST_CASE $OPERATOR $NB_ITERATIONS 

Where $SCRIPT.cde is one of the four available scripts.
The other variables are explained in the option (1) section.

Please note that we successfully run these scripts on different machines BUT the
performance using the CDE package were really unstable and poor.
Moreover the CILK_NWORKERS environment variable used to set to use the number of Cilk
threads has no effect.

How to read the results
-----------------------
The proto-application output is in 4 parts:
 - The first part summarizes the chosen parameter of execution.
 - The second part gives the execution time of all the preliminary steps.
 - The third part is the one which interest us. It gives the execution time (in RDTSC
   cycles) of the of the matrix assembly and the preconditioner creation for each
   iteration. The preconditioner creation step  also includes the MPI halos
   communication. Thus, the assembly step is achieved once both of these two steps are
   done. The measures made in the publication correspond to the sum of these two steps.
 - The last part is the numerical checking.
