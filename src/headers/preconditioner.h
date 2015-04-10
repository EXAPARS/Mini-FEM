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

#ifndef PRECOND_H
#define PRECOND_H

// External preconditioner Fortran functions
extern "C"
void ela_invert_prec_ (int *dimNode, int *nbNodes, int *nodeToNodeIndex,
                       int *nodeToNodeValue, double *prec, int *error,
                       int *checkBounds, int *curNode);

// Create preconditioner for elasticity operator
void preconditioner_ela (double *prec, double *nodeToNodeValue, int *nodeToNodeRow,
                         int *nodeToNodeColumn, int *intfIndex, int *intfNodes,
                         int *neighborList, int *checkBounds, int nbNodes,
                         int nbBlocks, int nbIntf, int nbIntfNodes, int operatorDim,
                         int operatorID, int rank);

// Create preconditioner for laplacian operator
void preconditioner_lap (double *prec, double *nodeToNodeValue, int *nodeToNodeRow,
                         int *nodeToNodeColumn, int *intfIndex, int *intfNodes,
                         int *neighborList, int nbNodes, int nbBlocks, int nbIntf,
                         int nbIntfNodes, int operatorDim, int operatorID, int rank);

// Call the appropriate function to create the preconditioner
void preconditioner (double *prec, double *nodeToNodeValue, int *nodeToNodeRow,
                     int *nodeToNodeColumn, int *intfIndex, int *intfNodes,
                     int *neighborList, int *checkBounds, int nbNodes, int nbBlocks,
                     int nbIntf, int nbIntfNodes, int operatorDim, int operatorID,
#ifdef GASPI
                     int rank, gaspi_pointer_t srcSegmentPtr,
                     gaspi_pointer_t destSegmentPtr,
                     const gaspi_segment_id_t srcSegmentID,
                     const gaspi_segment_id_t destSegmentID,
                     const gaspi_queue_id_t queueID);
#else
                     int rank);
#endif

#endif
