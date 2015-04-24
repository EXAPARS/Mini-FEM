/*  Copyright 2014 - UVSQ, Dassault Aviation
    Authors list: Loïc Thébault, Eric Petit

    This file is part of Mini-FEM.

    Mini-FEM is free software: you can redistribute it and/or modify it under the terms
    of the GNU Lesser General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later version.

    Mini-FEM is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License along with
    Mini-FEM. If not, see <http://www.gnu.org/licenses/>. */

#ifdef XMPI
    #include <mpi.h>
#elif GASPI
    #include <GASPI.h>
    #include "GASPI_handler.h"
#endif
#include <iostream>
#include <iomanip>
#include <DC.h>

#include "globals.h"
#include "FEM.h"
#include "matrix.h"
#include "coloring.h"
#include "IO.h"

// External Fortran functions
extern "C" {
	void dqmrd4_ (int *nbNodes, int *boundNodesCode, int *nbBoundNodes,
				  int *boundNodesList, int *error);
	void e_essbcm_(int *dimNode, int *nbNodes, int *nbBoundNodes,
				   int *boundNodesList, int *boundNodesCode, int *checkBounds);
}

// Global variables
string meshName, operatorName;
int *colorToElem = nullptr;
int nbTotalColors;
//int MAX_ELEM_PER_PART = strtol (getenv ("elemPerPart"), nullptr, 0);

// Help message
void help () {
	cerr << "Please specify:\n"
		 << " 1. The test case: LM6, EIB or FGN.\n"
		 << " 2. The operator: lap or ela.\n"
		 << " 3. The number of iterations.\n";
}

// Check arguments (test case, operator & number of iterations)
void check_args (int argCount, char **argValue, int *nbIter, int rank)
{
    if (argCount < 4) {
        if (rank == 0) help ();
        exit (EXIT_FAILURE);
    }
    meshName = argValue[1];
    if (meshName.compare ("LM6") && meshName.compare ("EIB") &&
        meshName.compare ("FGN")) {
        if (rank == 0) {
		    cerr << "Incorrect argument \"" << meshName << "\".\n";
		    help ();
        }
		exit (EXIT_FAILURE);
	}
	operatorName = argValue[2];
	if (operatorName.compare ("lap") && operatorName.compare ("ela")) {
        if (rank == 0) {
		    cerr << "Incorrect argument \"" << operatorName << "\".\n";
		    help ();
        }
		exit (EXIT_FAILURE);
	}
	*nbIter = strtol (argValue[3], nullptr, 0);
	if (*nbIter < 1) {
        if (rank == 0) cerr << "Number of iterations must be at least 1.\n";
		exit (EXIT_FAILURE);
	}
}

int main (int argCount, char **argValue)
{
	// Process initialization
	int nbBlocks = 0, rank = 0;
    #ifdef XMPI
        MPI_Init (&argCount, &argValue);
        MPI_Comm_size (MPI_COMM_WORLD, &nbBlocks);
        MPI_Comm_rank (MPI_COMM_WORLD, &rank);
    #elif GASPI
        gaspi_proc_init (GASPI_BLOCK);
        gaspi_proc_num ((gaspi_rank_t*)&nbBlocks);
        gaspi_proc_rank ((gaspi_rank_t*)&rank);
    #endif

    // Arguments initialization
    int nbIter;
    check_args (argCount, argValue, &nbIter, rank);
    if (rank == 0) {
        cout << "\t\t* Mini-App *\n\n"
        	 << "Test case              : \"" << meshName << "\"\n"
        	 << "Operator               : \"" << operatorName << "\"\n"
             << "Elements per partition :  "  << MAX_ELEM_PER_PART << "\n"
        	 << "Iterations             :  "  << nbIter << "\n\n"
             << scientific << setprecision (1);
    }

	// Declarations
	DC_timer timer;
    index_t nodeToElem;
	double *coord, *nodeToNodeValue, *prec;
	double t1, t2;
	int *nodeToNodeRow, *nodeToNodeColumn, *elemToNode, *intfIndex, *intfNodes,
        *dispList, *neighborList, *boundNodesCode, *boundNodesList,
        *checkBounds, *elemToEdge = nullptr;
	int nbElem, nbNodes, nbEdges, nbIntf, nbIntfNodes, nbDispNodes,
        nbBoundNodes, operatorDim, operatorID, error;

    // Set the operator dimension & ID
    if (!operatorName.compare ("lap")) {
        operatorDim = 1;
        operatorID  = 0;
    }
    else {
        operatorDim = DIM_NODE * DIM_NODE;
        operatorID  = 1;
    }

    // Get the input data from DefMesh
	if (rank == 0) {
        cout << "Reading input data...                ";
        timer.time_start ();
        t1 = DC_get_time ();
    }
	read_input_data (&coord, &elemToNode, &neighborList, &intfIndex, &intfNodes,
                     &dispList, &boundNodesCode, &nbElem, &nbNodes, &nbEdges, &nbIntf,
                     &nbIntfNodes, &nbDispNodes, &nbBoundNodes, nbBlocks, rank);
    if (rank == 0) {
        t2 = DC_get_time ();
        timer.time_stop ();
        cout << "done  (" << t2 - t1 << " seconds)\n";
        cout << "done  (" << timer.get_avg_time () << " seconds)\n";
        cout << "done  (" << timer.get_avg_cycles () << " cycles)\n";
    }

    // D&C version
    #if defined (DC) || defined (DC_HYBRID)
        // Set the path to the D&C tree and permutations
        #ifdef DC_HYBRID
            string treePath = (string)DATA_PATH + "/" + meshName + "/DC_tree/Hybrid_"
        #else
            string treePath = (string)DATA_PATH + "/" + meshName + "/DC_tree/DC_"
        #endif
                              + to_string ((long long)MAX_ELEM_PER_PART) + "_"
                              + to_string ((long long)nbBlocks) + "_"
                              + to_string ((long long)rank);

        // Creation of the D&C tree and permutations
        #ifdef TREE_CREATION
            if (rank == 0) {
                cout << "Creation of the D&C tree...          ";
                t1 = DC_get_time ();
            }
            DC_create_tree (elemToNode, nbElem, DIM_ELEM, nbNodes);
            if (rank == 0) {
                t2 = DC_get_time ();
            	cout << "done  (" << t2 - t1 << " seconds)\n";
            }

        // Reading of the D&C tree and permutations
        #else
            if (rank == 0) {
                cout << "Reading the D&C tree...              ";
                t1 = DC_get_time ();
            }
            DC_read_tree (treePath, nbElem, nbNodes);
            if (rank == 0) {
                t2 = DC_get_time ();
            	cout << "done  (" << t2 - t1 << " seconds)\n";
            }
        #endif

        // Apply permutations
        if (rank == 0) {
            cout << "Applying permutation...              ";
            t1 = DC_get_time ();
        }
        DC_permute_double_2d_array (coord, nbNodes, DIM_NODE);
        #ifndef TREE_CREATION
            DC_permute_int_2d_array (elemToNode, nullptr, nbElem, DIM_ELEM, 0);
        #endif
        DC_renumber_int_array (elemToNode, nbElem * DIM_ELEM, true);
        DC_renumber_int_array (intfNodes, nbIntfNodes, true);
        DC_renumber_int_array (dispList, nbDispNodes, true);
        DC_permute_int_1d_array (boundNodesCode, nbNodes);
        if (rank == 0) {
            t2 = DC_get_time ();
        	cout << "done  (" << t2 - t1 << " seconds)\n";
        }
    // Mesh coloring version
    #elif COLORING
        // Create the coloring
        if (rank == 0) {
            cout << "Coloring of the mesh...              ";
            t1 = DC_get_time ();
        }
        int *colorPerm = new int [nbElem];
        coloring_creation (elemToNode, colorPerm, nbElem, nbNodes);
        if (rank == 0) {
            t2 = DC_get_time ();
        	cout << "done  (" << t2 - t1 << " seconds)\n";
        }

        // Apply the element permutation
        if (rank == 0) {
            cout << "Applying permutation...              ";
            t1 = DC_get_time ();
        }
        DC_permute_int_2d_array (elemToNode, colorPerm, nbElem, DIM_ELEM, 0);
        delete[] colorPerm;
        if (rank == 0) {
            t2 = DC_get_time ();
        	cout << "done  (" << t2 - t1 << " seconds)\n";
        }
    #endif
	delete[] dispList;

    // Create the CSR matrix
	if (rank == 0) {
	    cout << "Creating CSR matrix...               ";
	    t1 = DC_get_time ();
    }
    nodeToElem.index = new int [nbNodes + 1];
    nodeToElem.value = new int [nbElem * DIM_ELEM];
    nodeToNodeRow    = new int [nbNodes + 1];
    nodeToNodeColumn = new int [nbEdges];
    DC_create_nodeToElem (nodeToElem, elemToNode, nbElem, DIM_ELEM, nbNodes);
    create_nodeToNode (nodeToNodeRow, nodeToNodeColumn, nodeToElem, elemToNode,
                       nbNodes);
    delete[] nodeToElem.value, delete[] nodeToElem.index;
	if (rank == 0) {
        t2 = DC_get_time ();
    	cout << "done  (" << t2 - t1 << " seconds)\n";
    }

    // Finalize and store the D&C tree
    #if (defined (DC) || defined (DC_HYBRID)) && defined (TREE_CREATION)
        if (rank == 0) {
            cout << "Finalizing the D&C tree...           ";
            t1 = DC_get_time ();
        }
        DC_finalize_tree (nodeToNodeRow, elemToNode, nbNodes);
        if (rank == 0) {
            t2 = DC_get_time ();
            cout << "done  (" << t2 - t1 << " seconds)\n"
                 << "Storing the D&C tree...              ";
            t1 = DC_get_time ();
        }
        DC_store_tree (treePath, nbElem, nbNodes);
        if (rank == 0) {
            t2 = DC_get_time ();
            cout << "done  (" << t2 - t1 << " seconds)\n";
        }
    #endif

    // Compute the index of each edge of each element
    #ifdef OPTIMIZED
        if (rank == 0) {
            cout << "Computing edges index...             ";
            t1 = DC_get_time ();
        }
        elemToEdge = new int [nbElem * VALUES_PER_ELEM];
        create_elemToEdge (nodeToNodeRow, nodeToNodeColumn, elemToNode, elemToEdge,
                          nbElem);
        if (rank == 0) {
            t2 = DC_get_time ();
            cout << "done  (" << t2 - t1 << " seconds)\n";
        }
    #endif

    // Compute the boundary conditions
    if (rank == 0) {
        cout << "Computing boundary conditions...     ";
        t1 = DC_get_time ();
    }
    int dimNode = DIM_NODE;
    boundNodesList = new int [nbBoundNodes];
    checkBounds    = new int [nbNodes * DIM_NODE];
    dqmrd4_ (&nbNodes, boundNodesCode, &nbBoundNodes, boundNodesList, &error);
    e_essbcm_ (&dimNode, &nbNodes, &nbBoundNodes, boundNodesList,
    		   boundNodesCode, checkBounds);
    delete[] boundNodesList, delete[] boundNodesCode;
    if (rank == 0) {
        t2 = DC_get_time ();
        cout << "done  (" << t2 - t1 << " seconds)\n";
    }

    // Initialization of the GASPI library
    #ifdef GASPI
        if (rank == 0) {
            cout << "Initializing GASPI lib...            ";
            t1 = DC_get_time ();
        }
        double *srcSegment = NULL, *destSegment = NULL;
	    int *destOffset = new int [nbIntf];
        gaspi_size_t segmentSize = nbIntfNodes * operatorDim * sizeof (double);
        gaspi_segment_id_t srcSegmentID, destSegmentID;
        gaspi_queue_id_t queueID;
        GASPI_init (&srcSegment, &destSegment, segmentSize, &srcSegmentID,
                    &destSegmentID, &queueID, rank);
        GASPI_offset_exchange (destOffset, intfIndex, neighborList, nbIntf, nbBlocks,
                               rank, operatorDim, destSegmentID, queueID);
        if (rank == 0) {
            t2 = DC_get_time ();
            cout << "done  (" << t2 - t1 << " seconds)\n";
        }
    #endif

    // Main loop with assembly, solver & update
    if (rank == 0) cout << "\nMain FEM loop...\n";
    nodeToNodeValue = new double [nbEdges * operatorDim];
    prec            = new double [nbNodes * operatorDim];
    FEM_loop (prec, coord, nodeToNodeValue, nodeToNodeRow, nodeToNodeColumn,
              elemToNode, elemToEdge, intfIndex, intfNodes, neighborList, checkBounds,
              nbElem, nbNodes, nbEdges, nbIntf, nbIntfNodes, nbIter, nbBlocks, rank,
    #ifdef XMPI
              operatorDim, operatorID);
    #elif GASPI
              operatorDim, operatorID, srcSegment, destSegment, destOffset,
              srcSegmentID, destSegmentID, queueID);
    #endif
    delete[] checkBounds, delete[] nodeToNodeColumn, delete[] nodeToNodeRow;
    delete[] neighborList, delete[] intfNodes, delete[] intfIndex;
    delete[] coord, delete[] elemToNode;
    #ifdef OPTIMIZED
        delete[] elemToEdge; 
    #endif

    // Check results on nodeToNodeValue & prec arrays
    if (rank == 0) cout << "done\n\nChecking results...\n" << setprecision (3);
    check_assembly (prec, nodeToNodeValue, nbEdges, nbNodes, operatorDim, nbBlocks,
                    rank);
    delete[] prec, delete[] nodeToNodeValue;
    if (rank == 0) cout << "done\n";

    #ifdef XMPI
	    MPI_Finalize ();
    #elif GASPI
        GASPI_finalize (destOffset, rank, srcSegmentID, destSegmentID, queueID);
        gaspi_proc_term (GASPI_BLOCK);
    #endif

	return EXIT_SUCCESS;
}
