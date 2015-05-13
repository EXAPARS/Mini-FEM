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

#ifndef GASPI_HANDLER_H
#define GASPI_HANDLER_H

#ifdef GASPI

#include <GASPI.h>

// Free the destination offset array, flush the GASPI queue & free the segments
void GASPI_finalize (int *destOffset, int nbBlocks, int rank,
                     gaspi_segment_id_t srcSegmentID, gaspi_segment_id_t destSegmentID,
                     gaspi_queue_id_t queueID);

// Get the adjacent domains destination offset
void GASPI_offset_exchange (int *destOffset, int *intfIndex, int *neighborList,
                            int nbIntf, int nbBlocks, int rank, int operatorDim,
                            gaspi_segment_id_t destSegmentID, gaspi_queue_id_t queueID);

// Initialization of the GASPI segments & creation of the segment pointers
void GASPI_init (double **srcSegment, double **destSegment, int **destOffset,
                 int nbIntf, int nbBlocks, int rank, gaspi_size_t segmentSize,
                 gaspi_segment_id_t *srcSegmentID, gaspi_segment_id_t *destSegmentID,
                 gaspi_queue_id_t *queueID);

#endif
#endif
