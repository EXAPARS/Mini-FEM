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
        MPI_Irecv (&(bufferRecv[node1*operatorDim]), size, MPI_DOUBLE_PRECISION,
                   source, tag, MPI_COMM_WORLD, &(neighborList[2*nbIntf+i]));
    }

    // Buffering local data with laplacian operator
    if (operatorID == 0) {
        for (int i = 0; i < nbIntf; i++) {
            for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
                for (int k = 0; k < operatorDim; k++) {
                    int tmpNode = intfNodes[j] - 1;
                    bufferSend[j*operatorDim+k] = prec[k*operatorDim+tmpNode];
                }
            }
        }
    }
    // Elasticity operator
    else {
        for (int i = 0; i < nbIntf; i++) {
            for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
                for (int k = 0; k < operatorDim; k++) {
                    int tmpNode = intfNodes[j] - 1;
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
        MPI_Send (&(bufferSend[node1*operatorDim]), size, MPI_DOUBLE_PRECISION,
                  dest, tag, MPI_COMM_WORLD);
    }

    // Waiting incoming data
    for (int i = 0; i < nbIntf; i++) {
        MPI_Wait (&(neighborList[2*nbIntf+i]), MPI_STATUS_IGNORE);
    }

    // Assembling local and incoming data with laplacian operator
    if (operatorID == 0) {
        for (int i = 0; i < nbIntf; i++) {
            for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
                for (int k = 0; k < operatorDim; k++) {
                    int tmpNode = intfNodes[j] - 1;
                    prec[k*operatorDim+tmpNode] += bufferRecv[j*operatorDim+k];
                }
            }
        }
    }
    // Elasticity operator
    else {
        for (int i = 0; i < nbIntf; i++) {
            for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
                for (int k = 0; k < operatorDim; k++) {
                    int tmpNode = intfNodes[j] - 1;
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

    // For each interface
    for (int i = 0; i < nbIntf; i++) {
        int node1  = intfIndex[i];
        int node2  = intfIndex[i+1];
        int size   = (node2 - node1) * operatorDim;
        int offset = node1 * operatorDim * sizeof (double);
        int dest   = neighborList[i] - 1;
        gaspi_notification_id_t sendNotifyID = rank;
        gaspi_return_t check;

        // Buffering local data with laplacian operator
        if (operatorID == 0) {
            for (int j = node1; j < node2; j++) {
                for (int k = 0; k < operatorDim; k++) {
                    int tmpNode = intfNodes[j] - 1;
                    srcSegment[j*operatorDim+k] = prec[k*operatorDim+tmpNode];
                }
            }
        }
        // Elasticity operator
        else {
            for (int j = node1; j < node2; j++) {
                for (int k = 0; k < operatorDim; k++) {
                    int tmpNode = intfNodes[j] - 1;
                    srcSegment[j*operatorDim+k] = prec[tmpNode*operatorDim+k];
                }
            }
        }

        // Sending local data to adjacent domains
        check = gaspi_write_notify (srcSegmentID, offset, dest, destSegmentID, offset,
                                    size, sendNotifyID, i+1, queueID, GASPI_BLOCK);
        if (check != GASPI_SUCCESS) {
            cerr << "Error at gaspi_write_notify from rank " << rank << endl;
            exit (EXIT_FAILURE);
        }

        gaspi_printf ("%d : %d sends notif %d to %d containing nodes %d to %d\n", i,
                      rank, sendNotifyID, dest, node1, node2);
    }

    gaspi_barrier (GASPI_GROUP_ALL, GASPI_BLOCK);
    if (rank == 0) {
        gaspi_printf ("Sends done\n");
    }

    // For each notification from adjacent domains
    while (nbNotifiesLeft > 0) {
        gaspi_notification_t recvNotifyValue;
        gaspi_notification_id_t recvNotifyID;
        gaspi_return_t check;

        // Wait & reset the first incoming notification
        check = gaspi_notify_waitsome (destSegmentID, 0, nbBlocks, &recvNotifyID,
                                       GASPI_BLOCK);
        if (check != GASPI_SUCCESS) {
            cerr << "Error at gaspi_notify_waitsome from rank " << rank << endl;
            exit (EXIT_FAILURE);
        }
        check = gaspi_notify_reset (destSegmentID, recvNotifyID, &recvNotifyValue);
        if (check != GASPI_SUCCESS) {
            cerr << "Error at gaspi_notify_reset from rank " << rank << endl;
            exit (EXIT_FAILURE);
        }
        nbNotifiesLeft--;

        gaspi_printf ("%d : %d receives & resets notif %d containing nodes %d to %d\n",
                      nbNotifiesLeft, rank, recvNotifyID, intfIndex[recvNotifyValue-1],
                      intfIndex[recvNotifyValue]);

        // Assembling local and incoming data with laplacian operator
        if (operatorID == 0) {
            for (int j = intfIndex[recvNotifyValue-1]; j < intfIndex[recvNotifyValue];
                 j++) {
                for (int k = 0; k < operatorDim; k++) {
                    int tmpNode = intfNodes[j] - 1;
                    prec[k*operatorDim+tmpNode] += destSegment[j*operatorDim+k];
                }
            }
        }
        // Elasticity operator
        else {
            for (int j = intfIndex[recvNotifyValue-1]; j < intfIndex[recvNotifyValue];
                 j++) {
                for (int k = 0; k < operatorDim; k++) {
                    int tmpNode = intfNodes[j] - 1;
                    prec[tmpNode*operatorDim+k] += destSegment[j*operatorDim+k];
                }
            }
        }
    }
}

#endif
