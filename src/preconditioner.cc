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

#ifdef CILK
    #include <cilk/cilk.h>
#endif

#include "globals.h"
#include "preconditioner.h"
#include "halo.h"

// Create preconditioner for elasticity operator
void preconditioner_ela (double *prec, double *buffer, double *nodeToNodeValue,
                         int *nodeToNodeRow, int *nodeToNodeColumn, int *intfIndex,
                         int *intfNodes, int *neighborList, int *checkBounds,
                         int nbNodes, int nbBlocks, int nbIntf, int nbIntfNodes,
                         int operatorDim, int operatorID, int rank)
{
	// Copy matrix diagonal into preconditioner
    #if defined (REF)
        for (int i = 0; i < nbNodes; i++) {
    #elif defined (DC) || defined (COLORING)
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbNodes; i++) {
        #elif CILK
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

    // Distributed communications
    if (nbBlocks > 1) {
        #ifdef XMPI
        #elif GASPI
        #endif
    }

    // Inversion of preconditioner
    int dimNode = DIM_NODE, error;
    #if defined (REF)
    	for (int i = 1; i <= nbNodes; i++) {
    #elif defined (DC) || defined (COLORING)
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 1; i <= nbNodes; i++) {
        #elif CILK
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
                         int nbIntf, int nbIntfNodes, int operatorDim, int operatorID,
                         int rank)
{
	// Copy matrix diagonal into preconditioner
    #if defined (REF)
    	for (int i = 0; i < nbNodes; i++) {
    #elif defined (DC) || defined (COLORING)
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbNodes; i++) {
        #elif CILK
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

	// Distributed communications
    if (nbBlocks > 1) {
        #ifdef XMPI
        #elif GASPI
        #endif
    }

	// Inversion of preconditioner
    #if defined (REF)
    	for (int i = 0; i < nbNodes; i++) {
    #elif defined (DC) || defined (COLORING)
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbNodes; i++) {
        #elif CILK
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
                     int operatorDim, int operatorID, int rank)
{
    // Preconditioner reset
    #if defined (REF)
        for (int i = 0; i < nbNodes * operatorDim; i++) prec[i] = 0;
    #elif defined (DC) || defined (COLORING)
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbNodes * operatorDim; i++) prec[i] = 0;
        #elif CILK
            cilk_for (int i = 0; i < nbNodes * operatorDim; i++) prec[i] = 0;
        #endif
    #endif

    if (operatorID == 0) {
        preconditioner_lap (prec, buffer, nodeToNodeValue, nodeToNodeRow,
                            nodeToNodeColumn, intfIndex, intfNodes, neighborList,
                            nbNodes, nbBlocks, nbIntf, nbIntfNodes, operatorDim,
                            operatorID, rank);
    }
    else {
        preconditioner_ela (prec, buffer, nodeToNodeValue, nodeToNodeRow,
                            nodeToNodeColumn, intfIndex, intfNodes, neighborList,
                            checkBounds, nbNodes, nbBlocks, nbIntf, nbIntfNodes,
                            operatorDim, operatorID, rank);
    }
}
