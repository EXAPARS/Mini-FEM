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
#include <iostream>

#include "globals.h"
#include "halo.h"

#ifdef XMPI

// Halo exchange between MPI ranks
void MPI_halo_exchange (double *prec, int *intfIndex, int *intfNodes,
                        int *neighborsList, int nbBlocks, int nbIntf, int nbIntfNodes,
                        int operatorDim, int rank)
{
    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;

    // Initialize communication buffers
    double *bufferSend = new double [nbIntfNodes*operatorDim],
           *bufferRecv = new double [nbIntfNodes*operatorDim];

    // Initialize reception from adjacent domains
    for (int i = 0; i < nbIntf; i++) {
        int node1  = intfIndex[i],
            node2  = intfIndex[i+1],
            size   = (node2 - node1) * operatorDim,
            source = neighborsList[i] - 1,
            tag    = neighborsList[i] + 100;
        MPI_Irecv (&(bufferRecv[node1*operatorDim]), size, MPI_DOUBLE, source, tag,
                   MPI_COMM_WORLD, &(neighborsList[2*nbIntf+i]));
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
        int node1 = intfIndex[i],
            node2 = intfIndex[i+1],
            size  = (node2 - node1) * operatorDim,
            dest  = neighborsList[i] - 1,
            tag   = rank + 101;
        MPI_Send (&(bufferSend[node1*operatorDim]), size, MPI_DOUBLE, dest, tag,
                  MPI_COMM_WORLD);
    }

    // Waiting incoming data
    for (int i = 0; i < nbIntf; i++) {
        MPI_Wait (&(neighborsList[2*nbIntf+i]), MPI_STATUS_IGNORE);
    }

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
                          double *destDataSegment, int *intfIndex, int *intfNodes,
                          int *neighborsList, int *intfDestOffsets, int nbBlocks,
                          int nbIntf, int operatorDim, int rank, int iter,
                          const gaspi_segment_id_t segment1,
                          const gaspi_segment_id_t segment2,
                          const gaspi_queue_id_t queueID)
{
    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;

    // Double buffering flip/flop
    gaspi_segment_id_t srcDataSegmentID, destDataSegmentID;
    if ((segment1 % 2) == 0) {
        srcDataSegmentID  = segment1;
        destDataSegmentID = segment2;
    }
    else {
        srcDataSegmentID  = segment2;
        destDataSegmentID = segment1;
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
        int node1       = intfIndex[i],
            node2       = intfIndex[i+1],
            size        = (node2 - node1)    * operatorDim * sizeof (double),
            localOffset = node1              * operatorDim * sizeof (double),
            destOffset  = intfDestOffsets[i] * operatorDim * sizeof (double),
            neighbor    = neighborsList[i] - 1;
        gaspi_notification_id_t sendNotifyID = iter * nbBlocks + rank;

        // Initialize source segment
        #ifdef REF
            for (int j = node1; j < node2; j++) {
        #else
            #ifdef OMP
                #pragma omp parallel for
                for (int j = node1; j < node2; j++) {
            #elif CILK
                cilk_for (int j = node1; j < node2; j++) {
            #endif
        #endif
            int tmpNode = intfNodes[j] - 1;
            for (int k = 0; k < operatorDim; k++) {
                srcDataSegment[j*operatorDim+k] = prec[tmpNode*operatorDim+k];
            }
        }

        // Send local data to adjacent domain
        SUCCESS_OR_DIE (gaspi_write_notify (srcDataSegmentID, localOffset, neighbor,
                                            destDataSegmentID, destOffset, size,
                                            sendNotifyID, rank+1, queueID,
                                            GASPI_BLOCK));
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
        gaspi_notification_t recvNotifyValue;
        int recvIntf;

        // Wait & reset the first incoming notification
        while (1) {
            gaspi_notification_id_t recvNotifyID;
            SUCCESS_OR_DIE (gaspi_notify_waitsome (destDataSegmentID, iter*nbBlocks,
                                                   nbBlocks, &recvNotifyID,
                                                   GASPI_BLOCK));
            SUCCESS_OR_DIE (gaspi_notify_reset (destDataSegmentID, recvNotifyID,
                                                &recvNotifyValue));
            if (recvNotifyValue) break;
        }

        // Look for the interface associated with the received notification
        for (recvIntf = 0; recvIntf < nbIntf; recvIntf++) {
            if ((neighborsList[recvIntf]-1) == (recvNotifyValue-1)) break;
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
                prec[tmpNode*operatorDim+k] += destDataSegment[j*operatorDim+k];
            }
        }
    }
}

#endif
#ifdef MULTITHREADED_COMM

// Wait for multithreaded GASPI notifications
void GASPI_multithreaded_wait (int nbBlocks)
{
    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;
}

// Send initialized parts of the preconditioner
void GASPI_multithreaded_send (void *userCommArgs, DCcommArgs_t *DCcommArgs)
{
    // Get user arguments
    userCommArgs_t *tmpCommArgs = (userCommArgs_t*)userCommArgs;
    double *prec           = tmpCommArgs->prec,
           *srcDataSegment = tmpCommArgs->srcDataSegment;
    int  *srcOffsetSegment = tmpCommArgs->srcOffsetSegment,
         *neighborsList    = tmpCommArgs->neighborsList;
    int nbBlocks           = tmpCommArgs->nbBlocks,
        nbIntf             = tmpCommArgs->nbIntf,
        operatorDim        = tmpCommArgs->operatorDim,
        rank               = tmpCommArgs->rank,
        iter               = tmpCommArgs->iter;
    const gaspi_segment_id_t srcDataSegmentID    = tmpCommArgs->srcDataSegmentID,
                             destDataSegmentID   = tmpCommArgs->destDataSegmentID,
                             srcOffsetSegmentID  = tmpCommArgs->srcOffsetSegmentID,
                             destOffsetSegmentID = tmpCommArgs->destOffsetSegmentID;
    const gaspi_queue_id_t queueID = tmpCommArgs->queueID;

    // Get D&C arguments
    int *intfIndex = DCcommArgs->intfIndex,
        *intfNodes = DCcommArgs->intfNodes,
        *intfDest  = DCcommArgs->intfDest;
    int intfOffset = DCcommArgs->intfOffset;

    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;

    // For each interface
    for (int i = 0; i < nbIntf; i++) {
        int node1 = intfIndex[i],
            node2 = intfIndex[i+1],
            dataSegmentSize     = (node2 - node1)      * operatorDim * sizeof (double),
            dataSegmentOffset   = (intfOffset + node1) * operatorDim * sizeof (double),
            offsetSegmentSize   = (node2 - node1)      * sizeof (int),
            offsetSegmentOffset = (intfOffset + node1) * sizeof (int),
            neighbor = neighborsList[i] - 1;
        gaspi_notification_id_t sendNotifyID = iter * nbBlocks + rank;

        // Initialize source segment
        for (int j = node1; j < node2; j++) {
            int tmpNode = intfNodes[j] - 1;
            for (int k = 0; k < operatorDim; k++) {
                srcDataSegment[(intfOffset+j)*operatorDim+k] =
                          prec[tmpNode*operatorDim+k];
            }
            srcOffsetSegment[intfOffset+j] = intfDest[j];
        }
/*
        // Send local data to adjacent domain
        SUCCESS_OR_DIE (gaspi_write_notify (srcDataSegmentID, dataSegmentOffset,
                                            neighbor, destDataSegmentID,
                                            dataSegmentOffset, dataSegmentSize,
                                            sendNotifyID, rank+1, queueID,
                                            GASPI_BLOCK));

        // Send data destination to adjacent domain
        SUCCESS_OR_DIE (gaspi_write_notify (srcOffsetSegmentID, offsetSegmentOffset,
                                            neighbor, destOffsetSegmentID,
                                            offsetSegmentOffset, offsetSegmentSize,
                                            sendNotifyID, rank+1, queueID,
                                            GASPI_BLOCK));
*/
    }
}

#endif
