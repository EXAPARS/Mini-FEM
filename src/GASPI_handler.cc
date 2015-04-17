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
#include "GASPI_handler.h"

using namespace std;

// Free the destination offset array, flush the GASPI queue & free the segments
void GASPI_finalize (int *destOffset, int rank, gaspi_segment_id_t srcSegmentID,
                     gaspi_segment_id_t destSegmentID, gaspi_queue_id_t queueID)
{
    gaspi_return_t check;
    delete[] destOffset;

    gaspi_wait (queueID, GASPI_BLOCK);
    gaspi_barrier (GASPI_GROUP_ALL, GASPI_BLOCK);
    gaspi_segment_delete (srcSegmentID);
    if (check != GASPI_SUCCESS) {
        cerr << "Error at source segment free from rank " << rank << endl;
        exit (EXIT_FAILURE);
    }
    gaspi_segment_delete (destSegmentID);
    if (check != GASPI_SUCCESS) {
        cerr << "Error at destination segment free from rank " << rank << endl;
        exit (EXIT_FAILURE);
    }
}

// Get the adjacent domains destination offset
void GASPI_offset_exchange (int *destOffset, int *intfIndex, int *neighborList,
                            int nbIntf, int nbBlocks, int rank, int operatorDim,
                            gaspi_segment_id_t destSegmentID, gaspi_queue_id_t queueID)
{
    // For each interface, send local offset to adjacent domain
    for (int i = 0; i < nbIntf; i++) {
        // The +1 is required since a notification value cannot be equal to 0...
        int localOffset = intfIndex[i] * operatorDim * sizeof (double) + 1;
        int dest = neighborList[i] - 1;
        gaspi_return_t check = gaspi_notify (destSegmentID, dest, rank, localOffset,
                                             queueID, GASPI_BLOCK);
        if (check != GASPI_SUCCESS) {
            cerr << "Notify error during offset exchange from rank " << rank << endl;
            exit (EXIT_FAILURE);
        }
    }

    // For each interface, receive the destination offset from adjacent domain
    for (int i = 0; i < nbIntf; i++) {
        int source = neighborList[i] - 1;
        gaspi_notification_t recvNotifyValue;
        gaspi_notification_id_t recvNotifyID;
        gaspi_return_t check;

        check = gaspi_notify_waitsome (destSegmentID, source, 1, &recvNotifyID,
                                       GASPI_BLOCK);
        if (check != GASPI_SUCCESS) {
            cerr << "Wait error during offset exchange from rank " << rank << endl;
            exit (EXIT_FAILURE);
        }
        check = gaspi_notify_reset (destSegmentID, recvNotifyID, &recvNotifyValue);
        if (check != GASPI_SUCCESS) {
            cerr << "Reset error during offset exchange from rank " << rank << endl;
            exit (EXIT_FAILURE);
        }

        // Remove the +1 of the local offset
        destOffset[i] = recvNotifyValue - 1;
    }
}

// Initialization of the GASPI segments & creation of the segment pointers
void GASPI_init (double **srcSegment, double **destSegment, gaspi_size_t segmentSize,
                 gaspi_segment_id_t *srcSegmentID, gaspi_segment_id_t *destSegmentID,
                 gaspi_queue_id_t *queueID, int rank)
{
    gaspi_pointer_t srcSegmentPtr = NULL, destSegmentPtr = NULL;
    gaspi_return_t check;

    *srcSegmentID  = 0;
    *destSegmentID = 1;
    *queueID = 0;

    check = gaspi_segment_create (*srcSegmentID, segmentSize, GASPI_GROUP_ALL,
                                  GASPI_BLOCK, GASPI_ALLOC_DEFAULT);
    if (check != GASPI_SUCCESS) {
        cerr << "Error at source segment creation from rank " << rank << endl;
        exit (EXIT_FAILURE);
    }
    check = gaspi_segment_create (*destSegmentID, segmentSize, GASPI_GROUP_ALL,
                                  GASPI_BLOCK, GASPI_ALLOC_DEFAULT);
    if (check != GASPI_SUCCESS) {
        cerr << "Error at destination segment creation from rank " << rank << endl;
        exit (EXIT_FAILURE);
    }
    check = gaspi_segment_ptr (*srcSegmentID, &srcSegmentPtr);
    if (check != GASPI_SUCCESS) {
        cerr << "Error at source segment pointer from rank " << rank << endl;
        exit (EXIT_FAILURE);
    }
    check = gaspi_segment_ptr (*destSegmentID, &destSegmentPtr);
    if (check != GASPI_SUCCESS) {
        cerr << "Error at destination segment pointer from rank " << rank << endl;
        exit (EXIT_FAILURE);
    }

    *srcSegment  = (double*)srcSegmentPtr;
    *destSegment = (double*)destSegmentPtr;
}

#endif
