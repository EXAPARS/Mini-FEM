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
#include <iostream>

#include "globals.h"
#include "halo.h"

#ifdef XMPI

// Halo exchange between MPI ranks
void MPI_halo_exchange (double *prec, int *intfIndex, int *intfNodes,
                        int *neighborList, int nbNodes, int nbIntf, int nbIntfNodes,
                        int operatorDim, int operatorID, int rank)
{
    // Initialize communication buffers
    double *bufferSend = new double [nbIntfNodes*operatorDim];
    double *bufferRecv = new double [nbIntfNodes*operatorDim];

    // Initializing reception from adjacent domains
    for (int i = 0; i < nbIntf; i++) {
        int node1  = intfIndex[i];
        int node2  = intfIndex[i+1];
        int size   = (node2 - node1) * operatorDim;
        int source = neighborList[i] - 1;
        int tag    = neighborList[i] + 100;
        MPI_Irecv (&(bufferRecv[node1*operatorDim]), size, MPI_DOUBLE, source, tag,
                   MPI_COMM_WORLD, &(neighborList[2*nbIntf+i]));
    }

    // Buffering local data with laplacian operator
    if (operatorID == 0) {
        for (int i = 0; i < nbIntf; i++) {
            for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
                int tmpNode = intfNodes[j] - 1;
                for (int k = 0; k < operatorDim; k++) {
                    bufferSend[j*operatorDim+k] = prec[k*operatorDim+tmpNode];
                }
            }
        }
    }
    // Elasticity operator
    else {
        for (int i = 0; i < nbIntf; i++) {
            for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
                int tmpNode = intfNodes[j] - 1;
                for (int k = 0; k < operatorDim; k++) {
                    bufferSend[j*operatorDim+k] = prec[tmpNode*operatorDim+k];
                }
            }
        }
    }

    // Sending local data to adjacent domains
    for (int i = 0; i < nbIntf; i++) {
        int node1 = intfIndex[i];
        int node2 = intfIndex[i+1];
        int size  = (node2 - node1) * operatorDim;
        int dest  = neighborList[i] - 1;
        int tag   = rank + 101;
        MPI_Send (&(bufferSend[node1*operatorDim]), size, MPI_DOUBLE, dest, tag,
                  MPI_COMM_WORLD);
    }

    // Waiting incoming data
    for (int i = 0; i < nbIntf; i++) {
        MPI_Wait (&(neighborList[2*nbIntf+i]), MPI_STATUS_IGNORE);
    }

    // Assembling local and incoming data with laplacian operator
    if (operatorID == 0) {
        for (int i = 0; i < nbIntf; i++) {
            for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
                int tmpNode = intfNodes[j] - 1;
                for (int k = 0; k < operatorDim; k++) {
                    prec[k*operatorDim+tmpNode] += bufferRecv[j*operatorDim+k];
                }
            }
        }
    }
    // Elasticity operator
    else {
        for (int i = 0; i < nbIntf; i++) {
            for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
                int tmpNode = intfNodes[j] - 1;
                for (int k = 0; k < operatorDim; k++) {
                    prec[tmpNode*operatorDim+k] += bufferRecv[j*operatorDim+k];
                }
            }
        }
    }

    // Free communication buffers
    delete[] bufferRecv, delete[] bufferSend;
}

#elif GASPI

// Halo exchange between GASPI ranks
void GASPI_halo_exchange (double *prec, int *intfIndex, int *intfNodes,
                          int *neighborList, int nbNodes, int nbBlocks, int nbIntf,
                          int nbIntfNodes, int operatorDim, int operatorID, int rank,
                          gaspi_pointer_t srcSegmentPtr, gaspi_pointer_t destSegmentPtr,
                          const gaspi_segment_id_t srcSegmentID,
                          const gaspi_segment_id_t destSegmentID,
                          const gaspi_queue_id_t queueID)
{
    double *srcSegment  = (double*)srcSegmentPtr;
    double *destSegment = (double*)destSegmentPtr;
    int nbNotifiesLeft = nbIntf;

    // For each interface, sends local offset to adjacent domain
    for (int i = 0; i < nbIntf; i++) {
        // The +1 is required since a notification value cannot be equal to 0...
        int localOffset = intfIndex[i] * operatorDim * sizeof (double) + 1;
        int dest = neighborList[i] - 1;
        gaspi_return_t check = gaspi_notify (destSegmentID, dest, nbBlocks+rank,
                                             localOffset, queueID, GASPI_BLOCK);
//        if (check != GASPI_SUCCESS) {
//            cerr << "Notify error during offset exchange from rank " << rank << endl;
//            exit (EXIT_FAILURE);
//        }
    }

    // For each interface
    for (int i = 0; i < nbIntf; i++) {
        int node1       = intfIndex[i];
        int node2       = intfIndex[i+1];
        int localOffset = node1 * operatorDim * sizeof (double);
        int size        = (node2 - node1) * operatorDim * sizeof (double);
        int neighbor    = neighborList[i] - 1;
        gaspi_notification_t destOffset;
        gaspi_notification_id_t recvNotifyID;
        gaspi_return_t check;

        // Initializes source segment with laplacian operator
        if (operatorID == 0) {
            for (int j = node1; j < node2; j++) {
                int tmpNode = intfNodes[j] - 1;
                for (int k = 0; k < operatorDim; k++) {
                    srcSegment[j*operatorDim+k] = prec[k*operatorDim+tmpNode];
                }
            }
        }
        // Elasticity operator
        else {
            for (int j = node1; j < node2; j++) {
                int tmpNode = intfNodes[j] - 1;
                for (int k = 0; k < operatorDim; k++) {
                    srcSegment[j*operatorDim+k] = prec[tmpNode*operatorDim+k];
                }
            }
        }

        // Receives the destination offset from adjacent domain
        check = gaspi_notify_waitsome (destSegmentID, nbBlocks+neighbor, 1,
                                       &recvNotifyID, GASPI_BLOCK);
//        if (check != GASPI_SUCCESS) {
//            cerr << "Wait error during offset exchange from rank " << rank << endl;
//            exit (EXIT_FAILURE);
//        }
        check = gaspi_notify_reset (destSegmentID, recvNotifyID, &destOffset);
//        if (check != GASPI_SUCCESS) {
//            cerr << "Reset error during offset exchange from rank " << rank << endl;
//            exit (EXIT_FAILURE);
//        }
        destOffset--; // Remove the +1 of the local offset

        // Sends local data to adjacent domain
        check = gaspi_write_notify (srcSegmentID, localOffset, neighbor, destSegmentID,
                                    destOffset, size, rank, 1, queueID, GASPI_BLOCK);
//        if (check != GASPI_SUCCESS) {
//            cerr << "Write error during data exchange from rank " << rank << endl;
//            exit (EXIT_FAILURE);
//        }
    }

    // For each notification from adjacent domains
    while (nbNotifiesLeft > 0) {
        int recvIntf;

//    for (int i = 0; i < nbIntf; i++) {
//        int source = neighborList[i] - 1;

        gaspi_notification_t recvNotifyValue;
        gaspi_notification_id_t recvNotifyID;
        gaspi_return_t check;

        // Wait & reset the first incoming notification
        check = gaspi_notify_waitsome (destSegmentID, 0, nbBlocks, &recvNotifyID,
                                       GASPI_BLOCK);
//        check = gaspi_notify_waitsome (destSegmentID, source, 1, &recvNotifyID,
//                                       GASPI_BLOCK);
//        if (check != GASPI_SUCCESS) {
//            cerr << "Wait error during data exchange from rank " << rank << endl;
//            exit (EXIT_FAILURE);
//        }
        check = gaspi_notify_reset (destSegmentID, recvNotifyID, &recvNotifyValue);
//        if (check != GASPI_SUCCESS) {
//            cerr << "Reset error during data exchange from rank " << rank << endl;
//            exit (EXIT_FAILURE);
//        }
        nbNotifiesLeft--;

        // Look for the interface corresponding to the source rank
        for (recvIntf = 0; recvIntf < nbIntf; recvIntf++) {
            if ((neighborList[recvIntf]-1) == recvNotifyID) break;
        }

        // Assembling local and incoming data with laplacian operator
        if (operatorID == 0) {
            for (int j = intfIndex[recvIntf]; j < intfIndex[recvIntf+1]; j++) {
            //for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
                int tmpNode = intfNodes[j] - 1;
                for (int k = 0; k < operatorDim; k++) {
                    prec[k*operatorDim+tmpNode] += destSegment[j*operatorDim+k];
                }
            }
        }
        // Elasticity operator
        else {
            for (int j = intfIndex[recvIntf]; j < intfIndex[recvIntf+1]; j++) {
            //for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
                int tmpNode = intfNodes[j] - 1;
                for (int k = 0; k < operatorDim; k++) {
                    prec[tmpNode*operatorDim+k] += destSegment[j*operatorDim+k];
                }
            }
        }
    }
}

#endif
