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
                        int *neighborList, int nbNodes, int nbBlocks, int nbIntf,
                        int nbIntfNodes, int operatorDim, int operatorID, int rank)
{
    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;

    // Initialize communication buffers
    double *bufferSend = new double [nbIntfNodes*operatorDim];
    double *bufferRecv = new double [nbIntfNodes*operatorDim];

    // Initialize reception from adjacent domains
    for (int i = 0; i < nbIntf; i++) {
        int node1  = intfIndex[i];
        int node2  = intfIndex[i+1];
        int size   = (node2 - node1) * operatorDim;
        int source = neighborList[i] - 1;
        int tag    = neighborList[i] + 100;
        MPI_Irecv (&(bufferRecv[node1*operatorDim]), size, MPI_DOUBLE, source, tag,
                   MPI_COMM_WORLD, &(neighborList[2*nbIntf+i]));
    }

    // Initialize send buffer with laplacian operator
    if (operatorID == 0) {
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
                    bufferSend[j*operatorDim+k] = prec[k*operatorDim+tmpNode];
                }
            }
        }
    }
    // Elasticity operator
    else {
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
                    prec[k*operatorDim+tmpNode] += bufferRecv[j*operatorDim+k];
                }
            }
        }
    }
    // Elasticity operator
    else {
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
    }

    // Free communication buffers
    delete[] bufferRecv, delete[] bufferSend;
}

#elif GASPI

// Halo exchange between GASPI ranks
void GASPI_halo_exchange (double *prec, double *srcSegment, double *destSegment,
                          int *intfIndex, int *intfNodes, int *neighborList,
                          int *destOffset, int nbNodes, int nbBlocks, int nbIntf,
                          int nbIntfNodes, int operatorDim, int operatorID, int rank,
                          int iter, const gaspi_segment_id_t segment1,
                          const gaspi_segment_id_t segment2,
                          const gaspi_queue_id_t queueID)
{
    // If there is only one domain, do nothing
    if (nbBlocks < 2) return;

    // Double buffering flip/flop
    gaspi_segment_id_t srcSegmentID, destSegmentID;
    if ((segment1 % 2) == 0) {
        srcSegmentID  = segment1;
        destSegmentID = segment2;
    }
    else {
        srcSegmentID  = segment2;
        destSegmentID = segment1;
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
        int node1       = intfIndex[i];
        int node2       = intfIndex[i+1];
        int localOffset = node1 * operatorDim * sizeof (double);
        int size        = (node2 - node1) * operatorDim * sizeof (double);
        int neighbor    = neighborList[i] - 1;
        gaspi_notification_id_t sendNotifyID = iter * nbBlocks + rank;

        // Initialize source segment with laplacian operator
        if (operatorID == 0) {
            #ifdef REF
                for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
            #else
                #ifdef OMP
                    #pragma omp parallel for
                    for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
                #elif CILK
                    cilk_for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
                #endif
            #endif
                int tmpNode = intfNodes[j] - 1;
                for (int k = 0; k < operatorDim; k++) {
                    srcSegment[j*operatorDim+k] = prec[k*operatorDim+tmpNode];
                }
            }
        }
        // Elasticity operator
        else {
            #ifdef REF
                for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
            #else
                #ifdef OMP
                    #pragma omp parallel for
                    for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
                #elif CILK
                    cilk_for (int j = intfIndex[i]; j < intfIndex[i+1]; j++) {
                #endif
            #endif
                int tmpNode = intfNodes[j] - 1;
                for (int k = 0; k < operatorDim; k++) {
                    srcSegment[j*operatorDim+k] = prec[tmpNode*operatorDim+k];
                }
            }
        }

        // Send local data to adjacent domain
        gaspi_write_notify (srcSegmentID, localOffset, neighbor, destSegmentID,
                            destOffset[i], size, sendNotifyID, rank+1, queueID,
                            GASPI_BLOCK);
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
            gaspi_notify_waitsome (destSegmentID, iter*nbBlocks, nbBlocks,
                                   &recvNotifyID, GASPI_BLOCK);
            gaspi_notify_reset (destSegmentID, recvNotifyID, &recvNotifyValue);
            if (recvNotifyValue) break;
        }

        // Look for the interface associated with the received notification
        for (recvIntf = 0; recvIntf < nbIntf; recvIntf++) {
            if ((neighborList[recvIntf]-1) == (recvNotifyValue-1)) break;
        }

        // Assemble local and incoming data with laplacian operator
        if (operatorID == 0) {
            #ifdef REF
                for (int j = intfIndex[recvIntf]; j < intfIndex[recvIntf+1]; j++) {
            #else
                #ifdef OMP
                    #pragma omp parallel for
                    for (int j = intfIndex[recvIntf]; j < intfIndex[recvIntf+1]; j++) {
                #elif CILK
                    cilk_for (int j = intfIndex[recvIntf]; j < intfIndex[recvIntf+1];
                              j++) {
                #endif
            #endif
                int tmpNode = intfNodes[j] - 1;
                for (int k = 0; k < operatorDim; k++) {
                    prec[k*operatorDim+tmpNode] += destSegment[j*operatorDim+k];
                }
            }
        }
        // Elasticity operator
        else {
            #ifdef REF
                for (int j = intfIndex[recvIntf]; j < intfIndex[recvIntf+1]; j++) {
            #else
                #ifdef OMP
                    #pragma omp parallel for
                    for (int j = intfIndex[recvIntf]; j < intfIndex[recvIntf+1]; j++) {
                #elif CILK
                    cilk_for (int j = intfIndex[recvIntf]; j < intfIndex[recvIntf+1];
                              j++) {
                #endif
            #endif
                int tmpNode = intfNodes[j] - 1;
                for (int k = 0; k < operatorDim; k++) {
                    prec[tmpNode*operatorDim+k] += destSegment[j*operatorDim+k];
                }
            }
        }
    }
}

#endif
