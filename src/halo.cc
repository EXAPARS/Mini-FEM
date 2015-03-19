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
#elif GASPI
    #include <GASPI.h>
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
    int source, dest, size, tag, node1, node2 = 0;

    // Initializing reception from adjacent domains
    for (int i = 0; i < nbIntf; i++) {
        node1  = node2;
        node2  = intfIndex[i+1];
        size   = (node2 - node1) * operatorDim;
        source = neighborList[i] - 1;
        tag    = neighborList[i] + 100;
        MPI_Irecv (&(bufferRecv[node1*operatorDim]), size, MPI_DOUBLE_PRECISION,
                   source, tag, MPI_COMM_WORLD, &(neighborList[2*nbIntf+i]));
    }

    // Buffering local data with laplacian operator
    node2 = 0;
    if (operatorID == 0) {
        for (int i = 0; i < nbIntf; i++) {
            node1 = node2;
            node2 = intfIndex[i+1];
            for (int j = node1; j < node2; j++) {
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
            node1 = node2;
            node2 = intfIndex[i+1];
            for (int j = node1; j < node2; j++) {
                for (int k = 0; k < operatorDim; k++) {
                    int tmpNode = intfNodes[j] - 1;
                    bufferSend[j*operatorDim+k] = prec[tmpNode*operatorDim+k];
                }
            }
        }
    }

    // Sending local data to adjacent domains
    node2 = 0;
    for (int i = 0; i < nbIntf; i++) {
        node1 = node2;
        node2 = intfIndex[i+1];
        size  = (node2 - node1) * operatorDim;
        dest  = neighborList[i] - 1;
        tag   = rank + 101;
        MPI_Send (&(bufferSend[node1*operatorDim]), size, MPI_DOUBLE_PRECISION,
                  dest, tag, MPI_COMM_WORLD);
    }

    // Waiting incoming data
    for (int i = 0; i < nbIntf; i++) {
        MPI_Wait (&(neighborList[2*nbIntf+i]), MPI_STATUS_IGNORE);
    }

    // Assembling local and incoming data with laplacian operator
    node2 = 0;
    if (operatorID == 0) {
        for (int i = 0; i < nbIntf; i++) {
            node1 = node2;
            node2 = intfIndex[i+1];
            for (int j = node1; j < node2; j++) {
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
            node1 = node2;
            node2 = intfIndex[i+1];
            for (int j = node1; j < node2; j++) {
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
                          int *neighborList, int nbNodes, int nbIntf, int nbIntfNodes,
                          int operatorDim, int operatorID, int rank)
{
    gaspi_pointer_t *srcSegment, *destSegment = NULL;
    int dest, size, offset, node1, node2 = 0;
    const gaspi_segment_id_t srcID = 0, destID = 1;
    gaspi_size_t segmentSize = nbIntfNodes * operatorDim;
    gaspi_notification_id_t nbNotifiesLeft = nbIntf;

    // Creation of the GASPI segments
    gaspi_segment_create (srcID, segmentSize, GASPI_GROUP_ALL, GASPI_BLOCK,
                          GASPI_ALLOC_DEFAULT);
    gaspi_segment_create (destID, segmentSize, GASPI_GROUP_ALL, GASPI_BLOCK,
                          GASPI_ALLOC_DEFAULT);
    gaspi_segment_ptr (srcID, srcSegment);
    gaspi_segment_ptr (destID, destSegment);

    // For each interface
    for (int i = 0; i < nbIntf; i++) {
        node1  = node2;
        node2  = intfIndex[i+1];
        size   = (node2 - node1) * operatorDim;
        offset = node1 * operatorDim * sizeof (double);
        dest   = neighborList[i] - 1;

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
        gaspi_write_notify (srcID, offset, dest, destID, offset, size, rank, 42, 0,
                            GASPI_BLOCK);
    }

    // For each notification from adjacent domains
    while (nbNotifiesLeft > 0) {
        gaspi_notification_id_t recvID, recvValue;
        gaspi_return_t check;

        // Wait & reset the first incoming notification
        check = gaspi_notify_waitsome (destID, 0, nbIntf, &recvID, GASPI_BLOCK);
        if (check == GASPI_ERROR || check == GASPI_TIMEOUT) {
            cerr << "Error at gaspi_notify_waitsome\n";
            exit (EXIT_FAILURE);
        }
        check = gaspi_notify_reset (destID, recvID, &recvValue);
        if (check == GASPI_ERROR) {
            cerr << "Error at gaspi_notify_reset\n";
            exit (EXIT_FAILURE);
        }
        nbNotifiesLeft--;

        // Assembling local and incoming data with laplacian operator
        node1 = intfIndex[recvID];
        node2 = intfIndex[recvID+1];
        if (operatorID == 0) {
            for (int j = node1; j < node2; j++) {
                for (int k = 0; k < operatorDim; k++) {
                    int tmpNode = intfNodes[j] - 1;
                    prec[k*operatorDim+tmpNode] += destSegment[j*operatorDim+k];
                }
            }
        }
        // Elasticity operator
        else {
            for (int j = node1; j < node2; j++) {
                for (int k = 0; k < operatorDim; k++) {
                    int tmpNode = intfNodes[j] - 1;
                    prec[tmpNode*operatorDim+k] += destSegment[j*operatorDim+k];
                }
            }
        }
    }

    // Free GASPI segments
    gaspi_segment_delete (srcID);
    gaspi_segment_delete (destID);
}

#endif
