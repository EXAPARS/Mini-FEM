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

#ifndef IO_H
#define IO_H

// Read reference norm of matrix & preconditioner arrays
void read_ref_assembly (double *refMatrixNorm, double *refPrecNorm, int nbBlocks,
                        int rank);

// Store reference norm of matrix & preconditioner arrays
extern "C"
void store_ref_assembly_ (double *refMatrix, double *refPrec, int *nbEdges,
                          int *nbNodes, int *operatorDim, int *nbBlocks,
                          int *rank);

// Read input data from DefMesh
void read_input_data (double **coord, int **elemToNode, int **neighborsList,
                      int **intfIndex, int **intfNodes, int **boundNodesCode,
                      int *nbElem, int *nbNodes, int *nbEdges, int *nbIntf,
                      int *nbIntfNodes, int *nbBoundNodes, int nbBlocks, int rank);

// Store necessary data from DefMesh
extern "C"
void store_input_data_ (double *coord, int *elemToNode, int *neighborsList,
                        int *intfIndex, int *intfNodes, int *boundNodesCode,
                        int *nbElem, int *dimElem, int *nbNodes, int *dimNode,
                        int *nbEdges, int *nbIntf, int *nbIntfNodes, int *nbBoundNodes,
                        int *nbBlocks, int *rank);

#endif
