#!/bin/bash

# Set the parameters
EXE_DIR=$HOME/lthebaul/Mini-FEM/exe
MACHINE_FILE=$EXE_DIR/GASPI_machine_file
EXE_FILE=$EXE_DIR/GASPI_exe.sh
TEST_CASE=EIB
VECTOR_LENGTH=AVX
NB_ITERATIONS=1

# Go to the appropriate directory, exit on failure
cd $EXE_DIR || exit

for DISTRI in 'XMPI' #'GASPI' 'XMPI'
do
    # Set the environment
    if [ $DISTRI == "GASPI" ]; then
        module load PrgEnv-intel/14.0.1 gpi2/1.1.1
        cat $PBS_NODEFILE | cut -d'.' -f1 > $MACHINE_FILE
    else
        module load PrgEnv-intel/14.0.1 impi/4.1.1.036
    fi

    for SHARED in 'CILK' #'OMP'
    do
        for VERSION in 'DC_TreeCreation' #'REF' 'DC' #'DC_HYBRID' 'COLORING_OMP'
        do
            echo -e "\n$VERSION version using $DISTRI and $SHARED"
            BINARY=$EXE_DIR/bin/miniFEM_$DISTRI\_$SHARED\_$VERSION

            for OPERATOR in 'ela' 'lap'
            do
                # Create the GASPI execution script
                if [ $DISTRI == "GASPI" ]; then
                    echo "#!/bin/sh" > $EXE_FILE
                    echo "cd $EXE_DIR" >> $EXE_FILE
                    echo "export LD_LIBRARY_PATH=$LD_LIBRARY_PATH" >> $EXE_FILE
                    echo "$BINARY $TEST_CASE $OPERATOR $NB_ITERATIONS" >> $EXE_FILE
                    chmod +x $EXE_FILE
                fi

                for PART_SIZE in 200 #50 200 500
                do
                    export elemPerPart=$PART_SIZE
               	    echo -e "\n$OPERATOR operator, $PART_SIZE elements max per partition"

                    for NB_PROCESS in 32 64 #1 4 8 12 16
                    do
                        for NB_THREADS in 1 #2 4 8 12 16
                        do
                            let "nbCores=$NB_PROCESS*$NB_THREADS"
                            if [ "$nbCores" -gt "$MAX_CORES" ]; then break; fi
                            if [ $VERSION == 'REF' ] && [ $NB_THREADS -gt 1 ]; then break; fi

                            echo "$NB_PROCESS processus, $NB_THREADS threads"
                            if [ $SHARED == "CILK" ]; then
                                export CILK_NWORKERS=$NB_THREADS
                            else
                                export OMP_NUM_THREADS=$NB_THREADS
                            fi

                            OUTPUT_FILE=$EXE_DIR/stdout_$TEST_CASE\_$NB_NODES\_$DISTRI\_$SHARED\_$VERSION\_$OPERATOR\_$PART_SIZE\_$NB_PROCESS\_$NB_THREADS

                            if [ $DISTRI == "GASPI" ]; then
                                gaspi_run -m $MACHINE_FILE $EXE_FILE > $OUTPUT_FILE
                            else
       	                        mpirun -np $NB_PROCESS $BINARY $TEST_CASE $OPERATOR $NB_ITERATIONS > $OUTPUT_FILE
                            fi

       	                    matrixASM=0
       	                    precond=0
                            firstASM=1
                            firstPrecond=1
       	                    while read line
       	                    do
       	                    	echo $line | grep -q "Matrix assembly"
       	                    	if [ $? == 0 ]; then
       	                    		tmpVar=`echo $line | cut -d " " -f5`
                                    echo $tmpVar >> ASM_$TEST_CASE\_$NB_NODES\_$DISTRI\_$SHARED\_$VERSION\_$OPERATOR\_$PART_SIZE\_$NB_PROCESS\_$NB_THREADS
                                    if [ $firstASM == 0 ]; then
       	                                let "matrixASM=$matrixASM+$tmpVar"
                                    fi
                                    firstASM=0
       	                    	fi
       	                    	echo $line | grep -q "Preconditioner"
       	                    	if [ $? == 0 ]; then
       	                    		tmpVar=`echo $line | cut -d " " -f4`
                                    echo $tmpVar >> Prec_$TEST_CASE\_$NB_NODES\_$DISTRI\_$SHARED\_$VERSION\_$OPERATOR\_$PART_SIZE\_$NB_PROCESS\_$NB_THREADS
                                    if [ $firstPrecond == 0 ]; then
       	                    		    let "precond=$precond+$tmpVar"
                                    fi
                                    firstPrecond=0
       	                    	fi
       	                    done < $OUTPUT_FILE
       	                    let "matrixASM=$matrixASM/($NB_ITERATIONS-1)"
       	                    let "precond=$precond/($NB_ITERATIONS-1)"
                            let "total=$matrixASM+$precond"

       	                    echo "Matrix assembly         : $matrixASM cycles"
       	                    echo "Preconditioner creation : $precond cycles"
       	                    echo "Total                   : $total cycles"
                        done
                    done
                done
            done
        done
    done
    if [ $DISTRI == "GASPI" ]; then
        rm $MACHINE_FILE $EXE_FILE
    fi
done > $EXE_DIR/miniFEM_$TEST_CASE\_$NB_NODES.log

exit
