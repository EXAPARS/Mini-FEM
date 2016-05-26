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

#ifndef FEM_H
#define FEM_H

#ifdef GASPI
    #include <GASPI.h>
#endif
#include <DC.h>

// Return the euclidean norm of given array
double compute_double_norm (double *tab, int size);

// Check if current results match to the reference version
void check_results (double *prec, double *nodeToNodeValue, int nbEdges, int nbNodes,
                    int operatorDim, int nbBlocks, int rank);

// Get the average measures from all ranks and keep the max
void get_average_cycles (DC_timer &ASMtimer, DC_timer &precInitTimer,
                         DC_timer &haloTimer, DC_timer &precInverTimer,
                         int nbBlocks, int rank);

// Main loop iterating over the 3 main steps of FEM applications
void FEM_loop (double *prec, double *coord, double *nodeToNodeValue,
               int *nodeToNodeRow, int *nodeToNodeColumn, int *elemToNode,
               int *elemToEdge, int *intfIndex, int *intfNodes, int *neighborsList,
               int *checkBounds, int nbElem, int nbNodes, int nbEdges, int nbIntf,
               int nbIntfNodes, int nbIter, int nbBlocks, int rank, int operatorDim,
#ifdef XMPI
               int operatorID);
#elif GASPI
               int operatorID, int nbMaxComm, int nbNotifications,
               double *srcDataSegment, double *dstDataSegment, int *srcOffsetSegment,
               int *dstOffsetSegment, int *intfDstIndex,
               gaspi_segment_id_t srcDataSegmentID,
               gaspi_segment_id_t dstDataSegmentID,
               gaspi_segment_id_t srcOffsetSegmentID,
               gaspi_segment_id_t dstOffsetSegmentID, gaspi_queue_id_t queueID);
#endif
#endif
