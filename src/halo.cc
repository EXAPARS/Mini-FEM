/*  Copyright 2014 - UVSQ, Dassault Aviation
    Authors list: Loïc Thébault, Eric Petit

    This file is part of MiniFEM.

    MiniFEM is free software: you can redistribute it and/or modify it under the terms
    of the GNU Lesser General Public License as published by the Free Software
    Foundation, either version 3 of the License, or (at your option) any later version.

    MiniFEM is distributed in the hope that it will be useful, but WITHOUT ANY
    WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
    PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License along with
    MiniFEM. If not, see <http://www.gnu.org/licenses/>. */

#include <mpi.h>
#include <stdio.h>

#include "globals.h"
#include "halo.h"

// Halo exchange between MPI domains (only for elasticity operator)
void halo_exchange (double *prec, int *intfIndex, int *intfNodes, int *neighborList,
                    int nbNodes, int nbIntf, int nbIntfNodes, int operatorDim,
                    int mpiRank)
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

    // Buffering local data
    node2 = 0;
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

    // Sending local data to adjacent domains
    node2 = 0;
    for (int i = 0; i < nbIntf; i++) {
        node1 = node2;
        node2 = intfIndex[i+1];
        size  = (node2 - node1) * operatorDim;
        dest  = neighborList[i] - 1;
        tag   = mpiRank + 101;
        MPI_Send (&(bufferSend[node1*operatorDim]), size, MPI_DOUBLE_PRECISION, dest,
                  tag, MPI_COMM_WORLD);
    }

    // Waiting incoming data
    for (int i = 0; i < nbIntf; i++) {
        MPI_Wait (&(neighborList[2*nbIntf+i]), MPI_STATUS_IGNORE);
    }

    // Assembling local and incoming data
    node2 = 0;
    for (int i = 0; i < nbIntf; i++) {
        node1  = node2;
        node2  = intfIndex[i+1];
        for (int j = node1; j < node2; j++) {
            for (int k = 0; k < operatorDim; k++) {
                int tmpNode = intfNodes[j] - 1;
                prec[tmpNode*operatorDim+k] += bufferRecv[j*operatorDim+k];
            }
        }
    }

    // Free communication buffers
    delete[] bufferRecv, delete[] bufferSend;
}
