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

#include "globals.h"
#include "halo.h"

#ifdef XMPI

// Halo exchange between MPI ranks
void MPI_halo_exchange (double *prec, int *intfIndex, int *intfNodes,
                        int *neighborsList, int nbBlocks, int nbIntf, int nbIntfNodes,
                        int operatorDim, int rank)
{
    MPI_Request Req[nbIntf];

    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;

    // Initialize communication buffers
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
            dest  = neighborsList[i] - 1,
            tag   = rank + 101;
        MPI_Send (&(bufferSend[begin*operatorDim]), size, MPI_DOUBLE, dest, tag,
                  MPI_COMM_WORLD);
    }

    // Wait for incoming data
    MPI_Waitall(nbIntf, Req, MPI_STATUSES_IGNORE);
    //    for (int i = 0; i < nbIntf; i++) {
    //  MPI_Wait (&(Req[i]), MPI_STATUS_IGNORE);
    //}

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
                          int *neighborsList, int *intfDestIndex, int nbBlocks,
                          int nbIntf, int operatorDim, int rank, int iter,
                          const gaspi_segment_id_t srcDataSegmentID,
                          const gaspi_segment_id_t destDataSegmentID,
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
            destOffset  = intfDestIndex[i] * operatorDim * sizeof (double),
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
                                            destDataSegmentID, destOffset, size,
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
            SUCCESS_OR_DIE (gaspi_notify_waitsome (destDataSegmentID, 0, nbBlocks,
                                                   &notifyID, GASPI_BLOCK));
            SUCCESS_OR_DIE (gaspi_notify_reset (destDataSegmentID, notifyID,
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
                prec[tmpNode*operatorDim+k] += destDataSegment[j*operatorDim+k];
            }
        }
    }
}

#endif
#ifdef MULTITHREADED_COMM

// Wait for multithreaded GASPI notifications
void GASPI_multithreaded_wait (double *prec, double *destDataSegment, int *intfNodes,
                               int *destOffsetSegment, int nbNotifications,
                               int nbBlocks, int operatorDim, int iter,
                               gaspi_segment_id_t destOffsetSegmentID, int rank)
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
            SUCCESS_OR_DIE (gaspi_notify_waitsome (destOffsetSegmentID, 0, 65536,
                                                   &notifyID, GASPI_BLOCK));
            SUCCESS_OR_DIE (gaspi_notify_reset (destOffsetSegmentID, notifyID,
                                                &notifyValue));
            if (notifyValue) break;
        }

        // Assemble local and incoming data
        int begin = (uint16_t)notifyValue,
            end   = begin + (notifyValue >> 16);
        for (int j = begin; j < end; j++) {
            int dest = intfNodes[destOffsetSegment[j]] - 1;
            for (int k = 0; k < operatorDim; k++) {
                prec[dest*operatorDim+k] += destDataSegment[j*operatorDim+k];
            }
        }
    }
}

// Send initialized parts of the preconditioner
void GASPI_multithreaded_send (void *userCommArgs, DCcommArgs_t *DCcommArgs)
{
    // Get user arguments
    userCommArgs_t *tmpCommArgs = (userCommArgs_t*)userCommArgs;
    double *prec           = tmpCommArgs->prec,
           *srcDataSegment = tmpCommArgs->srcDataSegment;
    int  *srcOffsetSegment = tmpCommArgs->srcOffsetSegment,
         *neighborsList    = tmpCommArgs->neighborsList,
         *intfIndex        = tmpCommArgs->intfIndex,
         *intfDestIndex    = tmpCommArgs->intfDestIndex;
    int nbBlocks           = tmpCommArgs->nbBlocks,
        nbIntf             = tmpCommArgs->nbIntf,
        nbMaxComm          = tmpCommArgs->nbMaxComm,
        operatorDim        = tmpCommArgs->operatorDim,
        rank               = tmpCommArgs->rank;
//        iter               = tmpCommArgs->iter;
    const gaspi_segment_id_t srcDataSegmentID    = tmpCommArgs->srcDataSegmentID,
                             destDataSegmentID   = tmpCommArgs->destDataSegmentID,
                             srcOffsetSegmentID  = tmpCommArgs->srcOffsetSegmentID,
                             destOffsetSegmentID = tmpCommArgs->destOffsetSegmentID;
    const gaspi_queue_id_t queueID = tmpCommArgs->queueID;

    // Get D&C arguments
    int *intfDCindex  = DCcommArgs->intfIndex,
        *intfDCnodes  = DCcommArgs->intfNodes,
        *intfDCdest   = DCcommArgs->intfDest,
        *intfDCoffset = DCcommArgs->intfOffset,
        *commID       = DCcommArgs->commID;

    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;

    // For each interface
    for (int i = 0; i < nbIntf; i++) {

        // Go to the next one if current one is empty
        int begin = intfDCindex[i],
            end   = intfDCindex[i+1],
            size  = end - begin;
        if (size <= 0) continue;

        int srcOffset               = intfIndex[i]     + intfDCoffset[i],
            destOffset              = intfDestIndex[i] + intfDCoffset[i],
            srcDataSegmentOffset    = srcOffset  * operatorDim * sizeof (double),
            destDataSegmentOffset   = destOffset * operatorDim * sizeof (double),
            dataSegmentSize         = size       * operatorDim * sizeof (double),
            srcOffsetSegmentOffset  = srcOffset                * sizeof (int),
            destOffsetSegmentOffset = destOffset               * sizeof (int),
            offsetSegmentSize       = size                     * sizeof (int),
            neighbor                = neighborsList[i] - 1;
        //gaspi_notification_id_t notifyID = iter * nbBlocks * nbMaxComm +
        //                                   rank * nbMaxComm + commID[i];
        gaspi_notification_id_t notifyID = rank * nbMaxComm + commID[i];
        gaspi_notification_t notifyValue = size << 16 | destOffset;

        // Initialize source segment
        for (int j = 0; j < size; j++) {
            int tmpNode = intfDCnodes[begin+j] - 1;
            for (int k = 0; k < operatorDim; k++) {
                srcDataSegment[(srcOffset+j)*operatorDim+k] =
                          prec[tmpNode*operatorDim+k];
            }
            srcOffsetSegment[srcOffset+j] = intfDCdest[begin+j];
        }

        // Send local data to adjacent domain
        SUCCESS_OR_DIE (gaspi_write (srcDataSegmentID, srcDataSegmentOffset, neighbor,
                                     destDataSegmentID, destDataSegmentOffset,
                                     dataSegmentSize, queueID, GASPI_BLOCK));

        // Send data destination and notification to adjacent domain
        SUCCESS_OR_DIE (gaspi_write_notify (srcOffsetSegmentID, srcOffsetSegmentOffset,
                                            neighbor, destOffsetSegmentID,
                                            destOffsetSegmentOffset, offsetSegmentSize,
                                            notifyID, notifyValue, queueID,
                                            GASPI_BLOCK));
    }
}

#endif
