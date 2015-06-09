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
#include "assembly.h"

#ifdef DC_VEC
// Vectorially compute the elements coefficient
inline void elem_coef_vec (double elemCoef[DIM_ELEM][DIM_NODE][VEC_SIZE],
                           double *coord, int *elemToNode, int elem)
{
    double xa[VEC_SIZE], xb[VEC_SIZE], xc[VEC_SIZE];
    double ya[VEC_SIZE], yb[VEC_SIZE], yc[VEC_SIZE];
    double za[VEC_SIZE], zb[VEC_SIZE], zc[VEC_SIZE];
    double vol[VEC_SIZE];

    // Compute elements coefficient in vectorial
    for (int i = 0; i < DIM_ELEM; i++) {
        int node[VEC_SIZE];
        for (int j = 0; j < VEC_SIZE; j++) {
    	    node[j] = elemToNode[(elem+j)*DIM_ELEM+i] - 1;
        }
        for (int j = 0; j < DIM_NODE; j++) {
            elemCoef[i][j][:] = coord[node[:]*DIM_NODE+j];
        }
    }

    xa[:] = elemCoef[0][0][:] - elemCoef[3][0][:];
    xb[:] = elemCoef[2][0][:] - elemCoef[3][0][:];
    xc[:] = elemCoef[1][0][:] - elemCoef[3][0][:];
    ya[:] = elemCoef[0][1][:] - elemCoef[3][1][:];
    yb[:] = elemCoef[2][1][:] - elemCoef[3][1][:];
    yc[:] = elemCoef[1][1][:] - elemCoef[3][1][:];
    za[:] = elemCoef[0][2][:] - elemCoef[3][2][:];
    zb[:] = elemCoef[2][2][:] - elemCoef[3][2][:];
    zc[:] = elemCoef[1][2][:] - elemCoef[3][2][:];
    elemCoef[0][0][:] = yb[:] * zc[:] - yc[:] * zb[:];
    elemCoef[0][1][:] = zb[:] * xc[:] - zc[:] * xb[:];
    elemCoef[0][2][:] = xb[:] * yc[:] - xc[:] * yb[:];
    elemCoef[1][0][:] = ya[:] * zb[:] - yb[:] * za[:];
    elemCoef[1][1][:] = za[:] * xb[:] - zb[:] * xa[:];
    elemCoef[1][2][:] = xa[:] * yb[:] - xb[:] * ya[:];
    elemCoef[2][0][:] = yc[:] * za[:] - ya[:] * zc[:];
    elemCoef[2][1][:] = zc[:] * xa[:] - za[:] * xc[:];
    elemCoef[2][2][:] = xc[:] * ya[:] - xa[:] * yc[:];
    elemCoef[3][0][:] = - (elemCoef[0][0][:] + elemCoef[1][0][:] +
                           elemCoef[2][0][:]);
    elemCoef[3][1][:] = - (elemCoef[0][1][:] + elemCoef[1][1][:] +
                           elemCoef[2][1][:]);
    elemCoef[3][2][:] = - (elemCoef[0][2][:] + elemCoef[1][2][:] +
                           elemCoef[2][2][:]);
    vol[:] = xa[:] * elemCoef[0][0][:] + ya[:] * elemCoef[0][1][:] +
             za[:] * elemCoef[0][2][:];

    for (int i = 0; i < DIM_ELEM; i++) {
        for (int j = 0; j < DIM_NODE; j++) {
            elemCoef[i][j][:] *= (1. / vol[:]);
        }
    }
}
#endif

// Sequentially compute the elements coefficient
inline void elem_coef_seq (double elemCoef[DIM_ELEM][DIM_NODE], double *coord,
                           int *elemToNode, int elem)
{
    double xa, xb, xc, ya, yb, yc, za, zb, zc, vol;

    for (int i = 0; i < DIM_ELEM; i++) {
        int node = elemToNode[elem*DIM_ELEM+i] - 1;
        for (int j = 0; j < DIM_NODE; j++) {
            elemCoef[i][j] = coord[node*DIM_NODE+j];
        }
    }

    xa = elemCoef[0][0] - elemCoef[3][0];
    xb = elemCoef[2][0] - elemCoef[3][0];
    xc = elemCoef[1][0] - elemCoef[3][0];
    ya = elemCoef[0][1] - elemCoef[3][1];
    yb = elemCoef[2][1] - elemCoef[3][1];
    yc = elemCoef[1][1] - elemCoef[3][1];
    za = elemCoef[0][2] - elemCoef[3][2];
    zb = elemCoef[2][2] - elemCoef[3][2];
    zc = elemCoef[1][2] - elemCoef[3][2];
    elemCoef[0][0] = yb * zc - yc * zb;
    elemCoef[0][1] = zb * xc - zc * xb;
    elemCoef[0][2] = xb * yc - xc * yb;
    elemCoef[1][0] = ya * zb - yb * za;
    elemCoef[1][1] = za * xb - zb * xa;
    elemCoef[1][2] = xa * yb - xb * ya;
    elemCoef[2][0] = yc * za - ya * zc;
    elemCoef[2][1] = zc * xa - za * xc;
    elemCoef[2][2] = xc * ya - xa * yc;
    elemCoef[3][0] = - (elemCoef[0][0] + elemCoef[1][0] + elemCoef[2][0]);
    elemCoef[3][1] = - (elemCoef[0][1] + elemCoef[1][1] + elemCoef[2][1]);
    elemCoef[3][2] = - (elemCoef[0][2] + elemCoef[1][2] + elemCoef[2][2]);
    vol = xa * elemCoef[0][0] + ya * elemCoef[0][1] + za * elemCoef[0][2];

    elemCoef[:][:] *= (1. / vol);
}

#ifdef DC_VEC
// Vectorial version of elasticity assembly on a given element interval
void assembly_ela_vec (void *userArgs, DCargs_t *DCargs)
{
    // Get user arguments
    userArgs_t *tmpArgs = (userArgs_t*)userArgs;
    double *coord           = tmpArgs->coord,
           *nodeToNodeValue = tmpArgs->nodeToNodeValue;
    int *nodeToNodeRow      = tmpArgs->nodeToNodeRow,
        *nodeToNodeColumn   = tmpArgs->nodeToNodeColumn,
        *elemToNode         = tmpArgs->elemToNode,
        *elemToEdge         = tmpArgs->elemToEdge;
    int operatorDim         = tmpArgs->operatorDim;

    // Get D&C arguments
    int firstElem = DCargs->firstElem,
        lastElem  = DCargs->lastElem;

    // If leaf is not a separator, reset locally the CSR matrix
    if (DCargs->isSep == 0) {
        int firstEdge = DCargs->firstEdge * operatorDim;
        int nbEdges   = (DCargs->lastEdge + 1) * operatorDim - firstEdge;
        nodeToNodeValue[firstEdge:nbEdges] = 0;
    }

    // For each block of VEC_SIZE elements in the interval
    for (int elem = firstElem; elem <= lastElem; elem += VEC_SIZE) {

        // Compute the element coefficient
        double elemCoef[DIM_ELEM][DIM_NODE][VEC_SIZE];
        elem_coef_vec (elemCoef, coord, elemToNode, elem);

        // Optimized version with precomputed edge index
        #ifdef OPTIMIZED
            int ctr = 0;
            // For each edge of current elements
            for (int i = 0; i < DIM_ELEM; i++) {
                for (int j = 0; j < DIM_ELEM; j++) {
                    double edgeCoef[9][VEC_SIZE];
                    int index[VEC_SIZE];
                    // Get edges index
                    for (int k = 0; k < VEC_SIZE; k++) {
                        index[k] = elemToEdge[(elem+k)*VALUES_PER_ELEM+ctr];
                    }
                    // Compute edges coefficients
                    edgeCoef[0][:] = elemCoef[i][0][:] * elemCoef[j][0][:] * 2.25
                                   + elemCoef[i][1][:] * elemCoef[j][1][:]
                                   + elemCoef[i][2][:] * elemCoef[j][2][:];
                    edgeCoef[1][:] = elemCoef[i][0][:] * elemCoef[j][1][:] * 1.25;
                    edgeCoef[2][:] = elemCoef[i][0][:] * elemCoef[j][2][:] * 1.25;
                    edgeCoef[3][:] = elemCoef[i][1][:] * elemCoef[j][0][:] * 1.25;
                    edgeCoef[4][:] = elemCoef[i][0][:] * elemCoef[j][0][:]
                                   + elemCoef[i][1][:] * elemCoef[j][1][:] * 2.25
                                   + elemCoef[i][2][:] * elemCoef[j][2][:];
                    edgeCoef[5][:] = elemCoef[i][1][:] * elemCoef[j][2][:] * 1.25;
                    edgeCoef[6][:] = elemCoef[i][2][:] * elemCoef[j][0][:] * 1.25;
                    edgeCoef[7][:] = elemCoef[i][2][:] * elemCoef[j][1][:] * 1.25;
                    edgeCoef[8][:] = elemCoef[i][0][:] * elemCoef[j][0][:]
                                   + elemCoef[i][1][:] * elemCoef[j][1][:]
                                   + elemCoef[i][2][:] * elemCoef[j][2][:] * 2.25;

                    // Update edges values
                    for (int k = 0; k < VEC_SIZE; k++) {
                        nodeToNodeValue[index[k]*operatorDim+0] += edgeCoef[0][k];
                        nodeToNodeValue[index[k]*operatorDim+1] += edgeCoef[1][k];
                        nodeToNodeValue[index[k]*operatorDim+2] += edgeCoef[2][k];
                        nodeToNodeValue[index[k]*operatorDim+3] += edgeCoef[3][k];
                        nodeToNodeValue[index[k]*operatorDim+4] += edgeCoef[4][k];
                        nodeToNodeValue[index[k]*operatorDim+5] += edgeCoef[5][k];
                        nodeToNodeValue[index[k]*operatorDim+6] += edgeCoef[6][k];
                        nodeToNodeValue[index[k]*operatorDim+7] += edgeCoef[7][k];
                        nodeToNodeValue[index[k]*operatorDim+8] += edgeCoef[8][k];
                    }
                    ctr++;
                }
            }
        #else
            // For each edge of current element
            for (int j = 0; j < VEC_SIZE; j++) {
                for (int k = 0; k < DIM_ELEM; k++) {
                  	int node1 = elemToNode[(elem+j)*DIM_ELEM+k] - 1;
                	for (int l = 0; l < DIM_ELEM; l++) {
                  		int node2 = elemToNode[(elem+j)*DIM_ELEM+l];
                		for (int m = nodeToNodeRow[node1]; m < nodeToNodeRow[node1+1];
                             m++) {
                			if (nodeToNodeColumn[m] == node2) {
                                // Add element contribution
                                nodeToNodeValue[m*operatorDim+0]
                                    += elemCoef[k][0][j] * elemCoef[l][0][j] * 2.25
                                     + elemCoef[k][1][j] * elemCoef[l][1][j]
                                     + elemCoef[k][2][j] * elemCoef[l][2][j];
                                nodeToNodeValue[m*operatorDim+1]
                                    += elemCoef[k][0][j] * elemCoef[l][1][j] * 1.25;
                                nodeToNodeValue[m*operatorDim+2]
                                    += elemCoef[k][0][j] * elemCoef[l][2][j] * 1.25;
                                nodeToNodeValue[m*operatorDim+3]
                                    += elemCoef[k][1][j] * elemCoef[l][0][j] * 1.25;
                                nodeToNodeValue[m*operatorDim+4]
                                    += elemCoef[k][0][j] * elemCoef[l][0][j]
                                     + elemCoef[k][1][j] * elemCoef[l][1][j] * 2.25
                                     + elemCoef[k][2][j] * elemCoef[l][2][j];
                                nodeToNodeValue[m*operatorDim+5]
                                    += elemCoef[k][1][j] * elemCoef[l][2][j] * 1.25;
                                nodeToNodeValue[m*operatorDim+6]
                                    += elemCoef[k][2][j] * elemCoef[l][0][j] * 1.25;
                                nodeToNodeValue[m*operatorDim+7]
                                    += elemCoef[k][2][j] * elemCoef[l][1][j] * 1.25;
                                nodeToNodeValue[m*operatorDim+8]
                                    += elemCoef[k][0][j] * elemCoef[l][0][j]
                                     + elemCoef[k][1][j] * elemCoef[l][1][j]
                                     + elemCoef[k][2][j] * elemCoef[l][2][j] * 2.25;
                                break;
                			}
                		}
                	}
                }
            }
        #endif
    }
}

// Vectorial version of laplacian assembly on a given element interval
void assembly_lap_vec (void *userArgs, DCargs_t *DCargs)
{
    // Get user arguments
    userArgs_t *tmpArgs = (userArgs_t*)userArgs;
    double *coord           = tmpArgs->coord,
           *nodeToNodeValue = tmpArgs->nodeToNodeValue;
    int *nodeToNodeRow      = tmpArgs->nodeToNodeRow,
        *nodeToNodeColumn   = tmpArgs->nodeToNodeColumn,
        *elemToNode         = tmpArgs->elemToNode,
        *elemToEdge         = tmpArgs->elemToEdge;
    int operatorDim         = tmpArgs->operatorDim;

    // Get D&C arguments
    int firstElem = DCargs->firstElem,
        lastElem  = DCargs->lastElem;

    // If leaf is not a separator, reset locally the CSR matrix
    if (DCargs->isSep == 0) {
        int firstEdge = DCargs->firstEdge * operatorDim;
        int nbEdges   = (DCargs->lastEdge + 1) * operatorDim - firstEdge;
        nodeToNodeValue[firstEdge:nbEdges] = 0;
    }

    // For each block of VEC_SIZE elements in the interval
    for (int elem = firstElem; elem <= lastElem; elem += VEC_SIZE) {

        // Compute the element coefficient
        double elemCoef[DIM_ELEM][DIM_NODE][VEC_SIZE];
        elem_coef_vec (elemCoef, coord, elemToNode, elem);

        // Optimized version with precomputed edge index
        #ifdef OPTIMIZED
            int ctr = 0;
            // For each edge of current elements
            for (int i = 0; i < DIM_ELEM; i++) {
                for (int j = 0; j < DIM_ELEM; j++) {
                    double edgeCoef[VEC_SIZE];
                    int index[VEC_SIZE];
                    // Get edges index
                    for (int k = 0; k < VEC_SIZE; k++) {
                        index[k] = elemToEdge[(elem+k)*VALUES_PER_ELEM+ctr];
                    }
                    // Compute edges coefficient
                    edgeCoef[:] = (elemCoef[i][0][:] * elemCoef[j][0][:] +
                                   elemCoef[i][1][:] * elemCoef[j][1][:] +
                                   elemCoef[i][2][:] * elemCoef[j][2][:]);
                    // Update edges value
                    for (int k = 0; k < VEC_SIZE; k++) {
                        nodeToNodeValue[index[k]] += edgeCoef[k];
                    }
                    ctr++;
                }
            }
        #else
            // For each edge of current element
            for (int j = 0; j < VEC_SIZE; j++) {
                for (int k = 0; k < DIM_ELEM; k++) {
                    int node1 = elemToNode[(elem+j)*DIM_ELEM+k] - 1;
                	for (int l = 0; l < DIM_ELEM; l++) {
                  		int node2 = elemToNode[(elem+j)*DIM_ELEM+l];
                		for (int m = nodeToNodeRow[node1]; m < nodeToNodeRow[node1+1];
                             m++) {
                			if (nodeToNodeColumn[m] == node2) {
                                // Add element contribution
                		    	nodeToNodeValue[m]
                                    += (elemCoef[k][0][j] * elemCoef[l][0][j]
                                     +  elemCoef[k][1][j] * elemCoef[l][1][j]
                                     +  elemCoef[k][2][j] * elemCoef[l][2][j]);
                                break;
                			}
                		}
                	}
                }
            }
        #endif
    }
}
#endif

// Sequential version of elasticity assembly on a given element interval
#if defined (DC) || defined (DC_VEC)
void assembly_ela_seq (void *userArgs, DCargs_t *DCargs)
#else
void assembly_ela_seq (void *userArgs, int firstElem, int lastElem)
#endif
{
    // Get user arguments
    userArgs_t *tmpArgs = (userArgs_t*)userArgs;
    double *coord           = tmpArgs->coord,
           *nodeToNodeValue = tmpArgs->nodeToNodeValue,
           *prec            = tmpArgs->prec;
    int *nodeToNodeRow      = tmpArgs->nodeToNodeRow,
        *nodeToNodeColumn   = tmpArgs->nodeToNodeColumn,
        *elemToNode         = tmpArgs->elemToNode,
        *elemToEdge         = tmpArgs->elemToEdge;
    int operatorDim         = tmpArgs->operatorDim;

    #if defined (DC) || defined (DC_VEC)
        // Get D&C arguments
        int firstElem = DCargs->firstElem,
            lastElem  = DCargs->lastElem;
    #endif
    #ifdef DC
        // If leaf is not a separator, reset locally the CSR matrix
        if (DCargs->isSep == 0) {
            int firstEdge = DCargs->firstEdge * operatorDim;
            int nbEdges   = (DCargs->lastEdge + 1) * operatorDim - firstEdge;
            nodeToNodeValue[firstEdge:nbEdges] = 0;
        }
    #endif

    #ifdef COLORING
        // For each element of the interval in parallel
        #ifdef OMP
            //#pragma omp parallel for  // Disable because of a bug
            for (int elem = firstElem; elem <= lastElem; elem++) {
        #elif CILK
            cilk_for (int elem = firstElem; elem <= lastElem; elem++) {
        #endif
    #else
        // For each element of the interval in sequential
        for (int elem = firstElem; elem <= lastElem; elem++) {
    #endif

        // Compute the element coefficient
        double elemCoef[DIM_ELEM][DIM_NODE];
        elem_coef_seq (elemCoef, coord, elemToNode, elem);

        // Optimized version with precomputed edge index
        #ifdef OPTIMIZED
            int ctr = 0;
            // For each edge of current element
            for (int i = 0; i < DIM_ELEM; i++) {
                for (int j = 0; j < DIM_ELEM; j++) {
                    // Get edge index & update its values
                    int index = elemToEdge[elem*VALUES_PER_ELEM+ctr];
                    nodeToNodeValue[index*operatorDim+0]
                        += elemCoef[i][0] * elemCoef[j][0] * 2.25
                         + elemCoef[i][1] * elemCoef[j][1]
                         + elemCoef[i][2] * elemCoef[j][2];
                    nodeToNodeValue[index*operatorDim+1]
                        += elemCoef[i][0] * elemCoef[j][1] * 1.25;
                    nodeToNodeValue[index*operatorDim+2]
                        += elemCoef[i][0] * elemCoef[j][2] * 1.25;
                    nodeToNodeValue[index*operatorDim+3]
                        += elemCoef[i][1] * elemCoef[j][0] * 1.25;
                    nodeToNodeValue[index*operatorDim+4]
                        += elemCoef[i][0] * elemCoef[j][0]
                         + elemCoef[i][1] * elemCoef[j][1] * 2.25
                         + elemCoef[i][2] * elemCoef[j][2];
                    nodeToNodeValue[index*operatorDim+5]
                        += elemCoef[i][1] * elemCoef[j][2] * 1.25;
                    nodeToNodeValue[index*operatorDim+6]
                        += elemCoef[i][2] * elemCoef[j][0] * 1.25;
                    nodeToNodeValue[index*operatorDim+7]
                        += elemCoef[i][2] * elemCoef[j][1] * 1.25;
                    nodeToNodeValue[index*operatorDim+8]
                        += elemCoef[i][0] * elemCoef[j][0]
                         + elemCoef[i][1] * elemCoef[j][1]
                         + elemCoef[i][2] * elemCoef[j][2] * 2.25;
                    ctr++;
                }
            }
        #else
            // For each edge of current element
            for (int j = 0; j < DIM_ELEM; j++) {
              	int node1 = elemToNode[elem*DIM_ELEM+j] - 1;
            	for (int k = 0; k < DIM_ELEM; k++) {
              		int node2 = elemToNode[elem*DIM_ELEM+k];
            		for (int l = nodeToNodeRow[node1];
                             l < nodeToNodeRow[node1+1]; l++) {
            			if (nodeToNodeColumn[l] == node2) {
                            // Add element contribution
                            nodeToNodeValue[l*operatorDim+0]
                                += elemCoef[j][0] * elemCoef[k][0] * 2.25
                                 + elemCoef[j][1] * elemCoef[k][1]
                                 + elemCoef[j][2] * elemCoef[k][2];
                            nodeToNodeValue[l*operatorDim+1]
                                += elemCoef[j][0] * elemCoef[k][1] * 1.25;
                            nodeToNodeValue[l*operatorDim+2]
                                += elemCoef[j][0] * elemCoef[k][2] * 1.25;
                            nodeToNodeValue[l*operatorDim+3]
                                += elemCoef[j][1] * elemCoef[k][0] * 1.25;
                            nodeToNodeValue[l*operatorDim+4]
                                += elemCoef[j][0] * elemCoef[k][0]
                                 + elemCoef[j][1] * elemCoef[k][1] * 2.25
                                 + elemCoef[j][2] * elemCoef[k][2];
                            nodeToNodeValue[l*operatorDim+5]
                                += elemCoef[j][1] * elemCoef[k][2] * 1.25;
                            nodeToNodeValue[l*operatorDim+6]
                                += elemCoef[j][2] * elemCoef[k][0] * 1.25;
                            nodeToNodeValue[l*operatorDim+7]
                                += elemCoef[j][2] * elemCoef[k][1] * 1.25;
                            nodeToNodeValue[l*operatorDim+8]
                                += elemCoef[j][0] * elemCoef[k][0]
                                 + elemCoef[j][1] * elemCoef[k][1]
                                 + elemCoef[j][2] * elemCoef[k][2] * 2.25;
            				break;
            			}
            		}
            	}
            }
        #endif
    }

    #ifndef BULK_SYNCHRONOUS
        // If leaf is not a separator, reset locally the preconditioner
        #ifdef DC
            if (DCargs->isSep == 0) {
                int firstNode = DCargs->firstNode * operatorDim;
                int nbNodes   = (DCargs->lastNode + 1) * operatorDim - firstNode;
                prec[firstNode:nbNodes] = 0;
            }
        #endif

        // Compute partially the preconditioner on each node owned by current leaf
        #if defined (DC) || defined (DC_VEC)
            for (int i = 0; i < DCargs->nbOwnedNodes; i++) {
                int node = DCargs->ownedNodes[i];
                for (int j = nodeToNodeRow[node]; j < nodeToNodeRow[node+1]; j++) {
                    if (nodeToNodeColumn[j]-1 == node) {
                        for (int k = 0; k < operatorDim; k++) {
                            prec[node*operatorDim+k] = nodeToNodeValue[j*operatorDim+k];
                        }
                        break;
                    }
                }
            }
        #endif
    #endif
}

// Sequential version of laplacian assembly on a given element interval
#if defined (DC) || defined (DC_VEC)
void assembly_lap_seq (void *userArgs, DCargs_t *DCargs)
#else
void assembly_lap_seq (void *userArgs, int firstElem, int lastElem)
#endif
{
    // Get user arguments
    userArgs_t *tmpArgs = (userArgs_t*)userArgs;
    double *coord           = tmpArgs->coord,
           *nodeToNodeValue = tmpArgs->nodeToNodeValue,
           *prec            = tmpArgs->prec;
    int *nodeToNodeRow      = tmpArgs->nodeToNodeRow,
        *nodeToNodeColumn   = tmpArgs->nodeToNodeColumn,
        *elemToNode         = tmpArgs->elemToNode,
        *elemToEdge         = tmpArgs->elemToEdge;
    int operatorDim         = tmpArgs->operatorDim;

    // Get D&C arguments
    #if defined (DC) || defined (DC_VEC)
        int firstElem = DCargs->firstElem,
            lastElem  = DCargs->lastElem;
    #endif

    // If leaf is not a separator, reset locally the CSR matrix
    #ifdef DC
        if (DCargs->isSep == 0) {
            int firstEdge = DCargs->firstEdge * operatorDim;
            int nbEdges = (DCargs->lastEdge + 1) * operatorDim - firstEdge;
            nodeToNodeValue[firstEdge:nbEdges] = 0;
        }
    #endif

    // For each element of the interval
    #ifdef COLORING
        #ifdef OMP
            //#pragma omp parallel for  // Disabled because of a bug
            for (int elem = firstElem; elem <= lastElem; elem++) {
        #elif CILK
            cilk_for (int elem = firstElem; elem <= lastElem; elem++) {
        #endif
    #else
        for (int elem = firstElem; elem <= lastElem; elem++) {
    #endif

        // Compute the element coefficient
        double elemCoef[DIM_ELEM][DIM_NODE];
        elem_coef_seq (elemCoef, coord, elemToNode, elem);

        // Optimized version with precomputed edge index
        #ifdef OPTIMIZED
            int ctr = 0;
            // For each edge of current element
            for (int i = 0; i < DIM_ELEM; i++) {
                for (int j = 0; j < DIM_ELEM; j++) {
                    // Get edge index & update its value
                    int index = elemToEdge[elem*VALUES_PER_ELEM+ctr];
                    nodeToNodeValue[index] += (elemCoef[i][0] * elemCoef[j][0] +
                                               elemCoef[i][1] * elemCoef[j][1] +
                                               elemCoef[i][2] * elemCoef[j][2]);
                    ctr++;
                }
            }
        #else
            // For each edge of current element
            for (int j = 0; j < DIM_ELEM; j++) {
              	int node1 = elemToNode[elem*DIM_ELEM+j] - 1;
            	for (int k = 0; k < DIM_ELEM; k++) {
              		int node2 = elemToNode[elem*DIM_ELEM+k] - 1;
            		for (int l = nodeToNodeRow[node1];
                             l < nodeToNodeRow[node1+1]; l++) {
            			if (nodeToNodeColumn[l] == (node2 + 1)) {
                            // Add element contribution
            				nodeToNodeValue[l] += (elemCoef[j][0] * elemCoef[k][0] +
                                                   elemCoef[j][1] * elemCoef[k][1] +
                                                   elemCoef[j][2] * elemCoef[k][2]);
            				break;
            			}
            		}
            	}
            }
        #endif
    }

    #ifndef BULK_SYNCHRONOUS
        // If leaf is not a separator, reset locally the preconditioner
        #ifdef DC
            if (DCargs->isSep == 0) {
                int firstNode = DCargs->firstNode * operatorDim;
                int nbNodes   = (DCargs->lastNode + 1) * operatorDim - firstNode;
                prec[firstNode:nbNodes] = 0;
            }
        #endif

        // Compute partially the preconditioner on each node owned by current leaf
        #if defined (DC) || defined (DC_VEC)
            for (int i = 0; i < DCargs->nbOwnedNodes; i++) {
                int node = DCargs->ownedNodes[i];
                for (int j = nodeToNodeRow[node]; j < nodeToNodeRow[node+1]; j++) {
                    if (nodeToNodeColumn[j]-1 == node) {
                        prec[node] = nodeToNodeValue[j];
                        break;
                    }
                }
            }
        #endif
    #endif
}

#ifdef COLORING
// Iterate over the colors & execute the assembly step on the elements of a same color
// in parallel
void coloring_assembly (userArgs_t *userArgs, int operatorID)
{
    // For each color
    for (int color = 0; color < nbTotalColors; color++) {

        // Get the interval of elements of the current color
        int firstElem = colorToElem[color];
        int lastElem  = colorToElem[color+1] - 1;

        // Call assembly function using laplacian operator
        if (operatorID == 0) {
            assembly_lap_seq (userArgs, firstElem, lastElem);
        }
        // Using elasticity operator
        else {
            assembly_ela_seq (userArgs, firstElem, lastElem);
        }
    }
}
#endif

// Call the appropriate function to perform the assembly step
void assembly (double *coord, double *nodeToNodeValue, double *prec, int *nodeToNodeRow,
               int *nodeToNodeColumn, int *elemToNode, int *elemToEdge, int nbElem,
               int nbEdges, int operatorDim, int operatorID)
{
    // Create the structure containing all the arguments needed for ASM
    userArgs_t userArgs = {
        coord, nodeToNodeValue, prec,
        nodeToNodeRow, nodeToNodeColumn, elemToNode, elemToEdge,
        operatorDim
    };

    #ifdef REF
        // Sequential reset of CSR matrix
        for (int i = 0; i < nbEdges * operatorDim; i++) {
            nodeToNodeValue[i] = 0;
        }
        // Sequential assembly using laplacian operator
        if (operatorID == 0) {
            assembly_lap_seq (&userArgs, 0, nbElem-1);
        }
        // Using elasticity operator
        else {
            assembly_ela_seq (&userArgs, 0, nbElem-1);
        }
    #elif COLORING
        // Parallel reset of CSR matrix
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbEdges * operatorDim; i++) {
                nodeToNodeValue[i] = 0;
            }
        #elif CILK
            cilk_for (int i = 0; i < nbEdges * operatorDim; i++) {
                nodeToNodeValue[i] = 0;
            }
        #endif
        // Coloring parallel assembly
        coloring_assembly (&userArgs, operatorID);
    #else
        // D&C parallel assembly using laplacian operator
        if (operatorID == 0) {
            #ifdef DC_VEC
                DC_tree_traversal (assembly_lap_seq, assembly_lap_vec, &userArgs);
            #else
                DC_tree_traversal (assembly_lap_seq, nullptr, &userArgs);
            #endif
        }
        // Using elasticity operator
        else {
            #ifdef DC_VEC
                DC_tree_traversal (assembly_ela_seq, assembly_ela_vec, &userArgs);
            #else
                DC_tree_traversal (assembly_ela_seq, nullptr, &userArgs);
            #endif
        }
    #endif
}
