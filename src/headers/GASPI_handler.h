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

#define COMM_SIZE 100

// Free the destination offset array, flush the GASPI queue & free the segments
void GASPI_finalize (int *intfDstIndex, int nbBlocks,
                     gaspi_segment_id_t srcDataSegmentID,
                     gaspi_segment_id_t dstDataSegmentID,
                     gaspi_segment_id_t srcOffsetSegmentID,
                     gaspi_segment_id_t dstOffsetSegmentID,
                     gaspi_queue_id_t queueID);

// Waits until given queue is empty if it's at least half full
void GASPI_wait_for_queue_half_full (gaspi_queue_id_t queueID, int rank);

// Get the max number of communications
void GASPI_max_nb_communications (int *nbDCcomm, int *globalMax, int nbIntf,
                                  int nbBlocks);

// Get the number of notifications coming from adjacent domains
void GASPI_nb_notifications_exchange (int *intfIndex, int *neighborsList,
                                      int *nbNotifications, int nbIntf, int nbBlocks,
                                      int rank, gaspi_segment_id_t destOffsetSegmentID,
                                      gaspi_queue_id_t queueID);

// Get the adjacent domains destination offset
void GASPI_offset_exchange (int *intfDstIndex, int *intfIndex, int *neighborsList,
                            int nbIntf, int nbBlocks, int rank,
                            gaspi_segment_id_t dstOffsetSegmentID,
                            gaspi_queue_id_t queueID);

// Initialization of the GASPI segments & creation of the segment pointers
void GASPI_init (double **srcDataSegment, double **dstDataSegment,
                 int **srcOffsetSegment, int **dstOffsetSegment,
                 int **intfDstIndex, int nbIntf, int nbIntfNodes, int nbBlocks,
                 int operatorDim, gaspi_segment_id_t *srcDataSegmentID,
                 gaspi_segment_id_t *dstDataSegmentID,
                 gaspi_segment_id_t *srcOffsetSegmentID,
                 gaspi_segment_id_t *dstOffsetSegmentID, gaspi_queue_id_t *queueID);

#endif
#endif
