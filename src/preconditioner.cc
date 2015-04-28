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

// Inversion of the preconditioner
void prec_inversion (double *prec, int *nodeToNodeRow, int *nodeToNodeColumn,
                     int *checkBounds, int nbNodes, int operatorID)
{
    #ifdef REF
        for (int i = 0; i < nbNodes; i++) {
    #else
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbNodes; i++) {
        #elif CILK
            cilk_for (int i = 0; i < nbNodes; i++) {
        #endif
    #endif
        // Laplacian operator
        if (operatorID == 0) {
            prec[i] = 1.0 / prec[i];
        }
        // Elasticity operator
        else {
            int dimNode = DIM_NODE, curNode = i + 1, error;
            ela_invert_prec_ (&dimNode, &nbNodes, nodeToNodeRow, nodeToNodeColumn,
                              prec, &error, checkBounds, &curNode);
        }
    }
}

// Reset & initialization of the preconditioner
void prec_init (double *prec, double *nodeToNodeValue, int *nodeToNodeRow,
                int *nodeToNodeColumn, int nbNodes, int operatorDim, int operatorID)
{
    // Preconditioner reset
    #ifdef REF
        for (int i = 0; i < nbNodes * operatorDim; i++) prec[i] = 0;
    #else
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbNodes * operatorDim; i++) prec[i] = 0;
        #elif CILK
            cilk_for (int i = 0; i < nbNodes * operatorDim; i++) prec[i] = 0;
        #endif
    #endif

    // Copy matrix diagonal into preconditioner
    #ifdef REF
        for (int i = 0; i < nbNodes; i++) {
    #else
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbNodes; i++) {
        #elif CILK
            cilk_for (int i = 0; i < nbNodes; i++) {
        #endif
    #endif
        for (int j = nodeToNodeRow[i]; j < nodeToNodeRow[i+1]; j++) {
            if (nodeToNodeColumn[j]-1 == i) {
                if (operatorID == 0) {
                    prec[i] = nodeToNodeValue[j];
                }
                else {
                    for (int k = 0; k < operatorDim; k++) {
                        prec[i*operatorDim+k] = nodeToNodeValue[j*operatorDim+k];
                    }   
                }
                break;
            }
        }
    }
}
