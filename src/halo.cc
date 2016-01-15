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

#ifdef XMPI
    #include <mpi.h>
#endif
#ifdef CILK
    #include <cilk/cilk.h>
#endif
#include <pthread.h>
#include <iostream>

#include "globals.h"
#include "GASPI_handler.h"
#include "halo.h"

#ifdef MULTITHREADED_COMM
    extern int *segmentPtr, *commPtr;
    extern pthread_mutex_t *segmentMutex;
#endif

#ifdef XMPI

// Halo exchange between MPI ranks
void MPI_halo_exchange (double *prec, int *intfIndex, int *intfNodes,
                        int *neighborsList, int nbBlocks, int nbIntf, int nbIntfNodes,
                        int operatorDim, int rank)
{
    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;

    // Initialize communication buffers
    MPI_Request Req[nbIntf];
    double *bufferSend = new double [nbIntfNodes*operatorDim],
           *bufferRecv = new double [nbIntfNodes*operatorDim];

    // Initialize reception from adjacent domains
    for (int i = 0; i < nbIntf; i++) {
        int begin  = intfIndex[i],
            end    = intfIndex[i+1],
            size   = (end - begin) * operatorDim,
            source = neighborsList[i] - 1,
            tag    = neighborsList[i] + 100;
        MPI_Irecv (&(bufferRecv[begin*operatorDim]), size, MPI_DOUBLE, source, tag,
                   MPI_COMM_WORLD, &(Req[i]));
    }

    // Initialize send buffer
    #ifdef REF
        for (int i = 0; i < nbIntf; i++) {
            for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
    #else
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbIntf; i++) {
                #pragma omp parallel for
                for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
        #elif CILK
            cilk_for (int i = 0; i < nbIntf; i++) {
                cilk_for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
        #endif
    #endif
            int tmpNode = intfNodes[j] - 1;
            for (int k = 0; k < operatorDim; k++) {
                bufferSend[j*operatorDim+k] = prec[tmpNode*operatorDim+k];
            }
        }
    }

    // Sending local data to adjacent domains
    for (int i = 0; i < nbIntf; i++) {
        int begin = intfIndex[i],
            end   = intfIndex[i+1],
            size  = (end - begin) * operatorDim,
            dst  = neighborsList[i] - 1,
            tag   = rank + 101;
        MPI_Send (&(bufferSend[begin*operatorDim]), size, MPI_DOUBLE, dst, tag,
                  MPI_COMM_WORLD);
    }

    // Wait for incoming data
    MPI_Waitall (nbIntf, Req, MPI_STATUSES_IGNORE);

    // Assembling local and incoming data
    #ifdef REF
        for (int i = 0; i < nbIntf; i++) {
            for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
    #else
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbIntf; i++) {
                #pragma omp parallel for
                for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
        #elif CILK
            cilk_for (int i = 0; i < nbIntf; i++) {
                cilk_for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
        #endif
    #endif
            int tmpNode = intfNodes[j] - 1;
            for (int k = 0; k < operatorDim; k++) {
                prec[tmpNode*operatorDim+k] += bufferRecv[j*operatorDim+k];
            }
        }
    }

    // Free communication buffers
    delete[] bufferRecv, delete[] bufferSend;
}

#elif GASPI

// Halo exchange between GASPI ranks
void GASPI_halo_exchange (double *prec, double *srcDataSegment,
                          double *dstDataSegment, int *intfIndex, int *intfNodes,
                          int *neighborsList, int *intfDstIndex, int nbBlocks,
                          int nbIntf, int operatorDim, int rank, int iter,
                          const gaspi_segment_id_t srcDataSegmentID,
                          const gaspi_segment_id_t dstDataSegmentID,
                          const gaspi_queue_id_t queueID)
{
    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;

    // For each interface
    #ifdef REF
        for (int i = 0; i < nbIntf; i++) {
    #else
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbIntf; i++) {
        #elif CILK
            cilk_for (int i = 0; i < nbIntf; i++) {
        #endif
    #endif
        int begin       = intfIndex[i],
            end         = intfIndex[i+1],
            size        = (end - begin)    * operatorDim * sizeof (double),
            localOffset = begin            * operatorDim * sizeof (double),
            dstOffset   = intfDstIndex[i] * operatorDim * sizeof (double),
            neighbor    = neighborsList[i] - 1;
        gaspi_notification_id_t notifyID = rank;

        // Initialize source segment
        #ifdef REF
            for (int j = begin; j < end; j++) {
        #else
            #ifdef OMP
                #pragma omp parallel for
                for (int j = begin; j < end; j++) {
            #elif CILK
                cilk_for (int j = begin; j < end; j++) {
            #endif
        #endif
            int tmpNode = intfNodes[j] - 1;
            for (int k = 0; k < operatorDim; k++) {
                srcDataSegment[j*operatorDim+k] = prec[tmpNode*operatorDim+k];
            }
        }

        // Send local data to adjacent domain
        SUCCESS_OR_DIE (gaspi_write_notify (srcDataSegmentID, localOffset, neighbor,
                                            dstDataSegmentID, dstOffset, size,
                                            notifyID, rank+1, queueID, GASPI_BLOCK));
    }

    // For each interface
    #ifdef REF
        for (int i = 0; i < nbIntf; i++) {
    #else
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbIntf; i++) {
        #elif CILK
            cilk_for (int i = 0; i < nbIntf; i++) {
        #endif
    #endif
        gaspi_notification_t notifyValue;
        int recvIntf;

        // Wait & reset the first incoming notification
        while (1) {
            gaspi_notification_id_t notifyID;
            SUCCESS_OR_DIE (gaspi_notify_waitsome (dstDataSegmentID, 0, nbBlocks,
                                                   &notifyID, GASPI_BLOCK));
            SUCCESS_OR_DIE (gaspi_notify_reset (dstDataSegmentID, notifyID,
                                                &notifyValue));
            if (notifyValue) break;
        }

        // Look for the interface associated with the received notification
        for (recvIntf = 0; recvIntf < nbIntf; recvIntf++) {
            if ((neighborsList[recvIntf]-1) == (notifyValue-1)) break;
        }

        // Assemble local and incoming data
        #ifdef REF
            for (int j = intfIndex[recvIntf]; j < intfIndex[recvIntf+1]; j++) {
        #else
            #ifdef OMP
                #pragma omp parallel for
                for (int j = intfIndex[recvIntf]; j < intfIndex[recvIntf+1]; j++) {
            #elif CILK
                cilk_for (int j = intfIndex[recvIntf]; j < intfIndex[recvIntf+1]; j++){
            #endif
        #endif
            int tmpNode = intfNodes[j] - 1;
            for (int k = 0; k < operatorDim; k++) {
                prec[tmpNode*operatorDim+k] += dstDataSegment[j*operatorDim+k];
            }
        }
    }
}

#endif
#ifdef MULTITHREADED_COMM

// Wait for multithreaded GASPI notifications
void GASPI_multithreaded_wait (double *prec, double *dstDataSegment, int *intfNodes,
                               int *dstOffsetSegment, int nbNotifications,
                               int nbBlocks, int operatorDim, int rank,
                               gaspi_segment_id_t dstOffsetSegmentID)
{
    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;

    // For each incoming notification
    #ifdef REF
        for (int i = 0; i < nbNotifications; i++) {
    #else
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbNotifications; i++) {
        #elif CILK
            cilk_for (int i = 0; i < nbNotifications; i++) {
        #endif
    #endif
        gaspi_notification_id_t notifyID;
        gaspi_notification_t notifyValue;

        // Wait & reset
        while (1) {
            SUCCESS_OR_DIE (gaspi_notify_waitsome (dstOffsetSegmentID, 0, 65536,
                                                   &notifyID, GASPI_BLOCK));
            SUCCESS_OR_DIE (gaspi_notify_reset (dstOffsetSegmentID, notifyID,
                                                &notifyValue));
            if (notifyValue) break;
        }

        // Assemble local and incoming data
        int begin = (uint16_t)notifyValue,
            end   = begin + (notifyValue >> 16);
        for (int j = begin; j < end; j++) {
            int dst = intfNodes[dstOffsetSegment[j]] - 1;
            for (int k = 0; k < operatorDim; k++) {
                prec[dst*operatorDim+k] += dstDataSegment[j*operatorDim+k];
            }
        }
    }
}

// Send initialized parts of the preconditioner
void GASPI_multithreaded_send (void *userCommArgs, DCcommArgs_t *DCcommArgs)
{
    // Get user arguments
    userCommArgs_t *commArgs = (userCommArgs_t*)userCommArgs;
    double *prec           = commArgs->prec,
           *srcDataSegment = commArgs->srcDataSegment;
    int  *srcOffsetSegment = commArgs->srcOffsetSegment,
         *neighborsList    = commArgs->neighborsList,
         *intfIndex        = commArgs->intfIndex,
         *intfDstIndex     = commArgs->intfDstIndex;
    int nbBlocks           = commArgs->nbBlocks,
        nbIntf             = commArgs->nbIntf,
        nbMaxComm          = commArgs->nbMaxComm,
        operatorDim        = commArgs->operatorDim,
        rank               = commArgs->rank;
    const gaspi_segment_id_t srcDataSegmentID   = commArgs->srcDataSegmentID,
                             dstDataSegmentID   = commArgs->dstDataSegmentID,
                             srcOffsetSegmentID = commArgs->srcOffsetSegmentID,
                             dstOffsetSegmentID = commArgs->dstOffsetSegmentID;
    const gaspi_queue_id_t queueID = commArgs->queueID;

    // Get D&C arguments
    int *intfDCindex  = DCcommArgs->intfIndex,
        *intfDCnodes  = DCcommArgs->intfNodes,
        *intfDCdst    = DCcommArgs->intfDst,
        *commID       = DCcommArgs->commID;

    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;

    // For each interface
    for (int i = 0; i < nbIntf; i++) {

        // Go to the next one if current one is empty
        int begin    = intfDCindex[i],
            end      = intfDCindex[i+1],
            dataSize = end - begin;
        if (dataSize <= 0) continue;

        // Init variables
        int mySegmentPtr, myCommPtr, srcSegmentPtr, dataCommSize, offsetCommSize,
            neighbor = neighborsList[i] - 1;
        bool handleComm = false, lastComm = false;

        // Get the current pointers position and update them
        pthread_mutex_lock (&(segmentMutex[i]));
        mySegmentPtr   = segmentPtr[i],
        myCommPtr      = commPtr[i];
        segmentPtr[i] += dataSize;
        if (segmentPtr[i] >= (myCommPtr + COMM_SIZE)) {
            commPtr[i] += COMM_SIZE;
            handleComm = true;
        }
        pthread_mutex_unlock (&(segmentMutex[i]));
        srcSegmentPtr = intfIndex[i] + mySegmentPtr;

        // Check if it's the last comm
        if ((mySegmentPtr + dataSize) == (intfIndex[i+1] - intfIndex[i])) {
            handleComm = true;
            lastComm   = true;
        }

        // Initialize source segment
        for (int j = 0; j < dataSize; j++) {
            int tmpNode = intfDCnodes[begin+j] - 1;
            for (int k = 0; k < operatorDim; k++) {
                srcDataSegment[(srcSegmentPtr+j)*operatorDim+k] =
                          prec[tmpNode*operatorDim+k];
            }
            srcOffsetSegment[srcSegmentPtr+j] = intfDCdst[begin+j];
        }

        // If current D&C node handles a comm
        if (handleComm) {

            // If there is more than 1 comm per DC node, queue capacity may be exceeded
            int commSize = (mySegmentPtr + dataSize - myCommPtr);
            if ((commSize / COMM_SIZE) > 1) {
                cerr << "Error: communication max size is too small.\n";
                exit (EXIT_FAILURE);
            }

            // Init comm variables
            int srcCommPtr = intfIndex[i]    + myCommPtr,
                dstCommPtr = intfDstIndex[i] + myCommPtr,
                srcDataSegmentPtr   = srcCommPtr * operatorDim * sizeof (double),
                dstDataSegmentPtr   = dstCommPtr * operatorDim * sizeof (double),
                srcOffsetSegmentPtr = srcCommPtr               * sizeof (int),
                dstOffsetSegmentPtr = dstCommPtr               * sizeof (int);
            gaspi_notification_t notifyValue;
            gaspi_notification_id_t notifyID = rank * nbMaxComm + commID[i];

            // If it's a last comm
            if (lastComm) {
                dataCommSize   = commSize * operatorDim * sizeof (double);
                offsetCommSize = commSize               * sizeof (int);
                notifyValue    = commSize << 16 | dstCommPtr;
            }
            else {
                dataCommSize   = COMM_SIZE * operatorDim * sizeof (double);
                offsetCommSize = COMM_SIZE               * sizeof (int);
                notifyValue    = COMM_SIZE << 16 | dstCommPtr;
            }

            // Send local data to adjacent domain
            SUCCESS_OR_DIE (gaspi_write (srcDataSegmentID, srcDataSegmentPtr, neighbor,
                                         dstDataSegmentID, dstDataSegmentPtr,
                                         dataCommSize, queueID, GASPI_BLOCK));

            // Send data destination and notification to adjacent domain
            SUCCESS_OR_DIE (gaspi_write_notify (srcOffsetSegmentID,
                                                srcOffsetSegmentPtr, neighbor,
                                                dstOffsetSegmentID,
                                                dstOffsetSegmentPtr, offsetCommSize,
                                                notifyID, notifyValue, queueID,
                                                GASPI_BLOCK));
        }
    }
}

#endif
