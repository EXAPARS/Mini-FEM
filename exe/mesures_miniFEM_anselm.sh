#!/bin/sh

# Set the parameters
MAX_CORES=240
EXE_DIR=$HOME/lthebaul/Mini-FEM/exe
MACHINE_FILE=$EXE_DIR/machine_file
EXE_FILE=$EXE_DIR/env_file.sh
TEST_CASE=LM6
VECTOR_LENGTH=AVX
NB_ITERATIONS=5

# Go to the appropriate directory, exit on failure
cd $EXE_DIR || exit

# Load the required modules
module load gpi2/1.1.1 PrgEnv-intel/14.0.1

# Create the machine and the environment file
cat $PBS_NODEFILE | cut -d'.' -f1 > $MACHINE_FILE
echo "#!/bin/sh" > $EXE_FILE
echo "export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:$LD_LIBRARY_PATH" >> $EXE_FILE
chmod +x $EXE_FILE

for OPERATOR in 'ela' #'lap'
do
    for PART_SIZE in 200 #50 200 500
    do
        export elemPerPart=$PART_SIZE
   	    echo "$TEST_CASE $OPERATOR, $PART_SIZE elements max per partition"

        for VERSION in 'REF' #'COLORING_OMP' 'DC' 'DC_HYBRID'
        do
            BINARY=$EXE_DIR/bin/miniFEM_GASPI_CILK_$VERSION
            OUTPUT_FILE=$EXE_DIR/stdout_$VERSION\_$NB_NODES
            echo "$BINARY $TEST_CASE $OPERATOR $NB_ITERATIONS" >> $EXE_FILE
            echo -e "\n$VERSION"

            for NB_PROCESS in 1 #4 8 12 16 32
            do
                for NB_THREADS in 1 #2 4 8 12 16
                do
                    let "nbCores=$NB_PROCESS*$NB_THREADS"
                    if [ "$nbCores" -gt "$MAX_CORES" ]; then break; fi
                    if [ $VERSION == 'REF' ] && [ $NB_THREADS -gt 1 ]; then break; fi

                    export CILK_NWORKERS=$NB_THREADS
                    export OMP_NUM_THREADS=$NB_THREADS
                    echo "$NB_PROCESS processus, $NB_THREADS threads"

       	            #mpirun -np $NB_PROCESS $BINARY $TEST_CASE $OPERATOR $NB_ITERATIONS > $OUTPUT_FILE
                    gaspi_run -m $MACHINE_FILE $EXE_FILE  > $OUTPUT_FILE

       	            matrixASM=0
       	            precond=0
                    firstASM=1
                    firstPrecond=1
       	            while read line
       	            do
       	            	echo $line | grep -q "Matrix assembly"
       	            	if [ $? == 0 ]; then
       	            		tmpVar=`echo $line | cut -d " " -f5`
                            echo $tmpVar >> miniFEM_$TEST_CASE\_$NB_NODES\_$OPERATOR\_$PART_SIZE\_$VERSION\_$NB_PROCESS\_$NB_THREADS.ASM
                            if [ $firstASM == 0 ]; then
       	                        let "matrixASM=$matrixASM+$tmpVar"
                            fi
                            firstASM=0
       	            	fi
       	            	echo $line | grep -q "Preconditioner"
       	            	if [ $? == 0 ]; then
       	            		tmpVar=`echo $line | cut -d " " -f4`
                            echo $tmpVar >> miniFEM_$TEST_CASE\_$NB_NODES\_$OPERATOR\_$PART_SIZE\_$VERSION\_$NB_PROCESS\_$NB_THREADS.prec
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
done > $EXE_DIR/miniFEM_$TEST_CASE\_$NB_NODES.log

rm $MACHINE_FILE $EXE_FILE

exit
