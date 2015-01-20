#!/bin/sh

# Set the parameters
MAX_CORES=240
EXE_DIR=$HOME/lthebaul/MiniFEM/exe
TEST_CASE=EIB
VECTOR_LENGTH=AVX
NB_ITERATIONS=50

# Go to the appropriate directory, exit on failure
cd $EXE_DIR || exit

# Load the required modules
module load PrgEnv-intel/13.5.192

for OPERATOR in 'ela' #'lap'
do
    for PART_SIZE in 50 200 500
    do
        export elemPerPart=$PART_SIZE
   	    echo "$TEST_CASE $OPERATOR, $PART_SIZE elements max per partition"

        for VERSION in 'REF' 'COLORING_OMP' 'DC' 'DC_HYBRID'
        do
            BINARY=./bin/miniFEM_$VERSION
            OUTPUT_FILE=./stdout_$VERSION\_$NB_NODES
            echo -e "\n$VERSION"

            for NB_PROCESS in 1 4 8 12 16 32
            do
                for NB_THREADS in 1 2 4 8 12 16
                do
                    let "nbCores=$NB_PROCESS*$NB_THREADS"
                    if [ "$nbCores" -gt "$MAX_CORES" ]; then break; fi
                    if [ $VERSION == 'REF' ] && [ $NB_THREADS -gt 1 ]; then break; fi

                    export CILK_NWORKERS=$NB_THREADS
                    export OMP_NUM_THREADS=$NB_THREADS
                    echo "$NB_PROCESS processus, $NB_THREADS threads"

       	            mpirun -np $NB_PROCESS $BINARY $TEST_CASE $OPERATOR \
                               $NB_ITERATIONS > $OUTPUT_FILE

       	            matrixASM=0
       	            precond=0
                    firstASM=1
                    firstPrecond=1
       	            while read line
       	            do
       	            	echo $line | grep -q "Matrix assembly"
       	            	if [ $? == 0 ]; then
       	            		tmpVar=`echo $line | cut -d " " -f5`
                            echo $tmpVar >> miniApp_$TEST_CASE\_$NB_NODES\_$OPERATOR\_$PART_SIZE\_$VERSION\_$NB_PROCESS\_$NB_THREADS.ASM
                            if [ $firstASM == 0 ]; then
       	            		    let "matrixASM=$matrixASM+$tmpVar"
                            fi
                            firstASM=0
       	            	fi
       	            	echo $line | grep -q "Preconditioner"
       	            	if [ $? == 0 ]; then
       	            		tmpVar=`echo $line | cut -d " " -f4`
                            echo $tmpVar >> miniApp_$TEST_CASE\_$NB_NODES\_$OPERATOR\_$PART_SIZE\_$VERSION\_$NB_PROCESS\_$NB_THREADS.prec
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
        echo
    done
done > miniFEM_$TEST_CASE\_$NB_NODES.log

exit
