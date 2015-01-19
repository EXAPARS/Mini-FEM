/*  Copyright 2014 - UVSQ, Dassault Aviation
    Authors list: Loïc Thébault, Eric Petit

    This file is part of MiniFEM.

    MiniFEM is free software: you can redistribute it and/or modify it under the terms
    of the GNU Lesser General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later version.

    MiniFEM is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License along with
    MiniFEM. If not, see <http://www.gnu.org/licenses/>. */

#include <cilk/cilk.h>

#include "globals.h"
#include "preconditioner.h"

// Create preconditioner for elasticity operator
void preconditioner_ela (double *prec, double *buffer, double *nodeToNodeValue,
                         int *nodeToNodeRow, int *nodeToNodeColumn,
                         int *intfIndex, int *intfNodes, int *neighborList,
                         int *checkBounds, int nbNodes, int nbBlocks,
                         int nbIntf, int nbIntfNodes, int mpiRank)
{
	// Copy matrix diagonal into preconditioner
    int operatorDim = DIM_NODE * DIM_NODE;
    #ifdef REF
        for (int i = 0; i < nbNodes; i++) {
    #elif COLORING
        //#pragma omp parallel for
        cilk_for (int i = 0; i < nbNodes; i++) {
    #elif DC
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbNodes; i++) {
        #else
            cilk_for (int i = 0; i < nbNodes; i++) {
        #endif
    #endif
       	    for (int j = nodeToNodeRow[i]; j < nodeToNodeRow[i+1]; j++) {
	    	    if (nodeToNodeColumn[j]-1 == i) {
		    	    for (int k = 0; k < operatorDim; k++) {
			    	    prec[i*operatorDim+k] = nodeToNodeValue[j*operatorDim+k];
    			    }   
       				break;
	       		}
	        }
        }

	// MPI communications
    int dimNode = DIM_NODE;
	ela_comm_mpi_ (&dimNode, &nbNodes, prec, &nbBlocks, buffer, &nbIntf,
                   &nbIntfNodes, neighborList, intfIndex, intfNodes);

	// Inversion of preconditioner
    int error;

    #ifdef REF
    	for (int i = 1; i <= nbNodes; i++) {
    #elif COLORING
        //#pragma omp parallel for
        cilk_for (int i = 1; i <= nbNodes; i++) {
    #elif DC
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 1; i <= nbNodes; i++) {
        #else
            cilk_for (int i = 1; i <= nbNodes; i++) {
        #endif
    #endif
	       	int curNode = i;
            ela_invert_prec_ (&dimNode, &nbNodes, nodeToNodeRow, nodeToNodeColumn,
                              prec, &error, checkBounds, &curNode);
	    }
}

// Create preconditioner for laplacian operator
void preconditioner_lap (double *prec, double *buffer, double *nodeToNodeValue,
                         int *nodeToNodeRow, int *nodeToNodeColumn, int *intfIndex,
                         int *intfNodes, int *neighborList, int nbNodes, int nbBlocks,
                         int nbIntf, int nbIntfNodes, int mpiRank)
{
	// Copy matrix diagonal into preconditioner
    #ifdef REF
    	for (int i = 0; i < nbNodes; i++) {
    #elif COLORING
        //#pragma omp parallel for
    	cilk_for (int i = 0; i < nbNodes; i++) {
    #elif DC
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbNodes; i++) {
        #else
    	    cilk_for (int i = 0; i < nbNodes; i++) {
        #endif
    #endif
		    for (int j = nodeToNodeRow[i]; j < nodeToNodeRow[i+1]; j++) {
			    if (nodeToNodeColumn[j]-1 == i) {
				    prec[i] = nodeToNodeValue[j];
				    break;
			    }
		    }
	    }

	// MPI communications
	lap_comm_mpi_ (&nbNodes, prec, &nbBlocks, buffer, &nbIntf, &nbIntfNodes,
                   neighborList, intfIndex, intfNodes);

	// Inversion of preconditioner
    #ifdef REF
    	for (int i = 0; i < nbNodes; i++) {
    #elif COLORING
        //#pragma omp parallel for
    	cilk_for (int i = 0; i < nbNodes; i++) {
    #elif DC
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbNodes; i++) {
        #else
    	    cilk_for (int i = 0; i < nbNodes; i++) {
        #endif
    #endif
		    prec[i] = 1.0 / prec[i];
	    }
}

// Call the appropriate function to create the preconditioner
void preconditioner (double *prec, double *buffer, double *nodeToNodeValue,
                     int *nodeToNodeRow, int *nodeToNodeColumn, int *intfIndex,
                     int *intfNodes, int *neighborList, int *checkBounds,
                     int nbNodes, int nbBlocks, int nbIntf, int nbIntfNodes,
                     int operatorDim, int operatorID, int mpiRank)
{
    // Preconditioner reset
    #ifdef REF
        for (int i = 0; i < nbNodes * operatorDim; i++) prec[i] = 0;
    #elif COLORING
        //#pragma omp parallel for
        cilk_for (int i = 0; i < nbNodes * operatorDim; i++) prec[i] = 0;
    #elif DC
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbNodes * operatorDim; i++) prec[i] = 0;
        #else
            cilk_for (int i = 0; i < nbNodes * operatorDim; i++) prec[i] = 0;
        #endif
    #endif

    if (operatorID == 0) {
        preconditioner_lap (prec, buffer, nodeToNodeValue, nodeToNodeRow,
                            nodeToNodeColumn, intfIndex, intfNodes,
                            neighborList, nbNodes, nbBlocks, nbIntf,
                            nbIntfNodes, mpiRank);
    }
    else {
        preconditioner_ela (prec, buffer, nodeToNodeValue, nodeToNodeRow,
                            nodeToNodeColumn, intfIndex, intfNodes,
                            neighborList, checkBounds, nbNodes, nbBlocks,
                            nbIntf, nbIntfNodes, mpiRank);
    }
}
