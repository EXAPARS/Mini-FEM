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

#ifdef GASPI

#include <iostream>
#include "globals.h"
#include "GASPI_handler.h"

// Free the destination offset array, flush the GASPI queue & free the segments
void GASPI_finalize (int *destOffset, int nbBlocks, int rank,
                     gaspi_segment_id_t srcSegmentID, gaspi_segment_id_t destSegmentID,
                     gaspi_queue_id_t queueID)
{
    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;

    SUCCESS_OR_DIE (gaspi_wait (queueID, GASPI_BLOCK));
    SUCCESS_OR_DIE (gaspi_barrier (GASPI_GROUP_ALL, GASPI_BLOCK));
    SUCCESS_OR_DIE (gaspi_segment_delete (srcSegmentID));
    SUCCESS_OR_DIE (gaspi_segment_delete (destSegmentID));
    SUCCESS_OR_DIE (gaspi_proc_term (GASPI_BLOCK));
    delete[] destOffset;
}

// Get the adjacent domains destination offset
void GASPI_offset_exchange (int *destOffset, int *intfIndex, int *neighborList,
                            int nbIntf, int nbBlocks, int rank, int operatorDim,
                            gaspi_segment_id_t destSegmentID, gaspi_queue_id_t queueID)
{
    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;

    // For each interface, send local offset to adjacent domain
    for (int i = 0; i < nbIntf; i++) {
        // The +1 is required since a notification value cannot be equal to 0...
        int localOffset = intfIndex[i] * operatorDim * sizeof (double) + 1;
        int dest = neighborList[i] - 1;
        SUCCESS_OR_DIE (gaspi_notify (destSegmentID, dest, rank, localOffset, queueID,
                                      GASPI_BLOCK));
    }

    // For each interface, receive the destination offset from adjacent domain
    for (int i = 0; i < nbIntf; i++) {
        int source = neighborList[i] - 1;
        gaspi_notification_t recvNotifyValue;
        gaspi_notification_id_t recvNotifyID;

        SUCCESS_OR_DIE (gaspi_notify_waitsome (destSegmentID, source, 1, &recvNotifyID,
                                               GASPI_BLOCK));
        SUCCESS_OR_DIE (gaspi_notify_reset (destSegmentID, recvNotifyID,
                                            &recvNotifyValue));

        // Remove the +1 of the local offset
        destOffset[i] = recvNotifyValue - 1;
    }
}

// Initialization of the GASPI segments & creation of the segment pointers
void GASPI_init (double **srcSegment, double **destSegment, int **destOffset,
                 int nbIntf, int nbBlocks, int rank, gaspi_size_t segmentSize,
                 gaspi_segment_id_t *srcSegmentID, gaspi_segment_id_t *destSegmentID,
                 gaspi_queue_id_t *queueID)
{
    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;

    gaspi_pointer_t srcSegmentPtr = NULL, destSegmentPtr = NULL;
    *destOffset = new int [nbIntf];
    *srcSegmentID  = 0;
    *destSegmentID = 1;
    *queueID = 0;

    SUCCESS_OR_DIE (gaspi_segment_create (*srcSegmentID, segmentSize, GASPI_GROUP_ALL,
                                          GASPI_BLOCK, GASPI_ALLOC_DEFAULT));
    SUCCESS_OR_DIE (gaspi_segment_create (*destSegmentID, segmentSize, GASPI_GROUP_ALL,
                                          GASPI_BLOCK, GASPI_ALLOC_DEFAULT));
    SUCCESS_OR_DIE (gaspi_segment_ptr (*srcSegmentID, &srcSegmentPtr));
    SUCCESS_OR_DIE (gaspi_segment_ptr (*destSegmentID, &destSegmentPtr));

    *srcSegment  = (double*)srcSegmentPtr;
    *destSegment = (double*)destSegmentPtr;
}

#endif
