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

#include <cstdlib>
#include <pthread.h>
#include <cmath>

#include "globals.h"
#include "GASPI_handler.h"

int *segmentPtr = nullptr, *commPtr = nullptr;
pthread_mutex_t *segmentMutex = nullptr;
//int COMM_SIZE = strtol (getenv ("commSize"), nullptr, 0);

// Free the destination offset array, flush the GASPI queue & free the segments
void GASPI_finalize (int *intfDstIndex, int nbBlocks,
                     gaspi_segment_id_t srcDataSegmentID,
                     gaspi_segment_id_t dstDataSegmentID,
                     gaspi_segment_id_t srcOffsetSegmentID,
                     gaspi_segment_id_t dstOffsetSegmentID,
                     gaspi_queue_id_t queueID)
{
    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;
    delete[] segmentMutex, delete[] commPtr, delete[] segmentPtr;
    delete[] intfDstIndex;

    SUCCESS_OR_DIE (gaspi_wait (queueID, GASPI_BLOCK));
    SUCCESS_OR_DIE (gaspi_barrier (GASPI_GROUP_ALL, GASPI_BLOCK));
    SUCCESS_OR_DIE (gaspi_segment_delete (srcDataSegmentID));
    SUCCESS_OR_DIE (gaspi_segment_delete (dstDataSegmentID));
    SUCCESS_OR_DIE (gaspi_segment_delete (srcOffsetSegmentID));
    SUCCESS_OR_DIE (gaspi_segment_delete (dstOffsetSegmentID));
    SUCCESS_OR_DIE (gaspi_proc_term (GASPI_BLOCK));
}

// Waits until given queue is empty if it's at least half full
void GASPI_wait_for_queue_half_full (gaspi_queue_id_t queueID, int rank)
{
    gaspi_number_t queueSizeMax;
    gaspi_number_t queueSize;

    SUCCESS_OR_DIE (gaspi_queue_size_max (&queueSizeMax));
    SUCCESS_OR_DIE (gaspi_queue_size (queueID, &queueSize));

    if (queueSize > queueSizeMax) {
        fprintf (stderr, "Rank %d has exceeded its queue capacity: %d / %d\n",
                 rank, queueSize, queueSizeMax);
        exit (EXIT_FAILURE);
    }
    else if (queueSize >= queueSizeMax/2) {
        SUCCESS_OR_DIE (gaspi_wait (queueID, GASPI_BLOCK));
    }
}

// Get the max number of communications
void GASPI_max_nb_communications (int *nbDCcomm, int *globalMax, int nbIntf,
                                  int nbBlocks)
{
    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;

    // Get the max number of communications of local process
    int localMax = 0;
    for (int i = 0; i < nbIntf; i++) {
        if (nbDCcomm[i] > localMax) localMax = nbDCcomm[i];
    }

    // Get the max number of communications for all process
    SUCCESS_OR_DIE (gaspi_allreduce (&localMax, globalMax, 1, GASPI_OP_MAX,
                                     GASPI_TYPE_INT, GASPI_GROUP_ALL, GASPI_BLOCK));
}

// Get the number of notifications coming from adjacent domains
void GASPI_nb_notifications_exchange (int *intfIndex, int *neighborsList,
                                      int *nbNotifications, int nbIntf, int nbBlocks,
                                      int rank, gaspi_segment_id_t dstOffsetSegmentID,
                                      gaspi_queue_id_t queueID)
{
    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;

    // For each interface, send the number of notifications
    for (int i = 0; i < nbIntf; i++) {
        int begin  = intfIndex[i],
            end    = intfIndex[i+1],
            size   = end - begin,
            nbComm = ceil ((float)size / COMM_SIZE);

        // The +1 is required since a notification value cannot be equal to 0...
        SUCCESS_OR_DIE (gaspi_notify (dstOffsetSegmentID, neighborsList[i]-1, rank,
                                      nbComm+1, queueID, GASPI_BLOCK));
    }

    // For each interface, receive the number of notifications coming from adjacent
    // domain
    for (int i = 0; i < nbIntf; i++) {
        gaspi_notification_t notifyValue;
        gaspi_notification_id_t notifyID;
        SUCCESS_OR_DIE (gaspi_notify_waitsome (dstOffsetSegmentID, neighborsList[i]-1,
                                               1, &notifyID, GASPI_BLOCK));
        SUCCESS_OR_DIE (gaspi_notify_reset (dstOffsetSegmentID, notifyID,
                                            &notifyValue));

        // Remove the +1 of the local offset
        (*nbNotifications) += (notifyValue - 1);
    }

    // Ensure that number of notifications are received by all
    SUCCESS_OR_DIE (gaspi_barrier (GASPI_GROUP_ALL, GASPI_BLOCK));
}

// Get the adjacent domains destination offset
void GASPI_offset_exchange (int *intfDstIndex, int *intfIndex, int *neighborsList,
                            int nbIntf, int nbBlocks, int rank,
                            gaspi_segment_id_t dstOffsetSegmentID,
                            gaspi_queue_id_t queueID)
{
    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;

    // For each interface, send local offset to adjacent domain
    for (int i = 0; i < nbIntf; i++) {
        // The +1 is required since a notification value cannot be equal to 0...
        SUCCESS_OR_DIE (gaspi_notify (dstOffsetSegmentID, neighborsList[i]-1, rank,
                                      intfIndex[i]+1, queueID, GASPI_BLOCK));
    }

    // For each interface, receive the destination offset from adjacent domain
    for (int i = 0; i < nbIntf; i++) {
        gaspi_notification_t notifyValue;
        gaspi_notification_id_t notifyID;
        SUCCESS_OR_DIE (gaspi_notify_waitsome (dstOffsetSegmentID, neighborsList[i]-1,
                                               1, &notifyID, GASPI_BLOCK));
        SUCCESS_OR_DIE (gaspi_notify_reset (dstOffsetSegmentID, notifyID,
                                            &notifyValue));

        // Remove the +1 of the local offset
        intfDstIndex[i] = notifyValue - 1;
    }

    // Ensure that offsets are received by all
    SUCCESS_OR_DIE (gaspi_barrier (GASPI_GROUP_ALL, GASPI_BLOCK));
}

// Initialization of the GASPI segments & creation of the segment pointers
void GASPI_init (double **srcDataSegment, double **dstDataSegment,
                 int **srcOffsetSegment, int **dstOffsetSegment,
                 int **intfDstIndex, int nbIntf, int nbIntfNodes, int nbBlocks,
                 int operatorDim, gaspi_segment_id_t *srcDataSegmentID,
                 gaspi_segment_id_t *dstDataSegmentID,
                 gaspi_segment_id_t *srcOffsetSegmentID,
                 gaspi_segment_id_t *dstOffsetSegmentID, gaspi_queue_id_t *queueID)
{
    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;

    *intfDstIndex = new int [nbIntf];
    segmentPtr    = new int [nbIntf] ();
    commPtr       = new int [nbIntf] ();
    segmentMutex  = new pthread_mutex_t [nbIntf];
    for (int i = 0; i < nbIntf; i++) {
        segmentMutex[i] = PTHREAD_MUTEX_INITIALIZER;
    }

    gaspi_pointer_t srcDataSegmentPtr = NULL, dstDataSegmentPtr = NULL,
                  srcOffsetSegmentPtr = NULL, dstOffsetSegmentPtr = NULL;
    gaspi_size_t dataSegmentSize = nbIntfNodes * operatorDim * sizeof (double),
               offsetSegmentSize = nbIntfNodes * sizeof (int);

    *srcDataSegmentID   = 0;
    *dstDataSegmentID   = 1;
    *srcOffsetSegmentID = 3;
    *dstOffsetSegmentID = 4;
    *queueID            = 0;

    SUCCESS_OR_DIE (gaspi_segment_create (*srcDataSegmentID, dataSegmentSize,
                                          GASPI_GROUP_ALL, GASPI_BLOCK,
                                          GASPI_ALLOC_DEFAULT));
    SUCCESS_OR_DIE (gaspi_segment_create (*dstDataSegmentID, dataSegmentSize,
                                          GASPI_GROUP_ALL, GASPI_BLOCK,
                                          GASPI_ALLOC_DEFAULT));
    SUCCESS_OR_DIE (gaspi_segment_create (*srcOffsetSegmentID, offsetSegmentSize,
                                          GASPI_GROUP_ALL, GASPI_BLOCK,
                                          GASPI_ALLOC_DEFAULT));
    SUCCESS_OR_DIE (gaspi_segment_create (*dstOffsetSegmentID, offsetSegmentSize,
                                          GASPI_GROUP_ALL, GASPI_BLOCK,
                                          GASPI_ALLOC_DEFAULT));
    SUCCESS_OR_DIE (gaspi_segment_ptr (*srcDataSegmentID, &srcDataSegmentPtr));
    SUCCESS_OR_DIE (gaspi_segment_ptr (*dstDataSegmentID, &dstDataSegmentPtr));
    SUCCESS_OR_DIE (gaspi_segment_ptr (*srcOffsetSegmentID, &srcOffsetSegmentPtr));
    SUCCESS_OR_DIE (gaspi_segment_ptr (*dstOffsetSegmentID, &dstOffsetSegmentPtr));

    *srcDataSegment   = (double*)srcDataSegmentPtr;
    *dstDataSegment   = (double*)dstDataSegmentPtr;
    *srcOffsetSegment = (int*)srcOffsetSegmentPtr;
    *dstOffsetSegment = (int*)dstOffsetSegmentPtr;
}

#endif
