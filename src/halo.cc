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

// Halo exchange between MPI domains
void halo_exchange (double *prec, double *buffer, int *intfIndex, int *intfNodes,
                    int *neighborList, int nbNodes, int nbIntf, int nbIntfNodes,
                    int operatorDim, int mpiRank)
{
    double *bufferSend = new double [nbIntfNodes*DIM_NODE*operatorDim];
    double *bufferRecv = new double [nbIntfNodes*DIM_NODE*operatorDim];
    int source, dest, size, tag, node1, node2 = 0;

    // Initializing reception from adjacent domains
    printf ("C++\n");
    for (int i = 0; i < nbIntf; i++) {
        node1  = node2;
        node2  = intfIndex[i+1];
        size   = (node2 - node1) * operatorDim;
        source = neighborList[i] - 1;
        tag    = neighborList[i] + 99;
        MPI_Irecv (&(bufferRecv[node1]), size,
                   MPI_DOUBLE_PRECISION, source, tag, MPI_COMM_WORLD,
                   &(neighborList[2*nbIntf+i]));
    }

    // Buffering local data
    node2 = 0;
    for (int i = 0; i < nbIntf; i++) {
        node1  = node2;
        node2  = intfIndex[i+1];
        for (int j = 0; j < operatorDim; j++) {
            for (int k = node1; k < node2; k++) {
                int tmpNode = intfNodes[k] - 1;
                bufferSend[k*operatorDim+j] = prec[tmpNode*operatorDim+j];
                if (i == 0 && j == 0 && k == node1) {
                    printf ("%d : prec[%d*%d+%d] = %f\n", mpiRank, tmpNode,
                            operatorDim, j, prec[tmpNode*operatorDim+j]);
                    for (int z = 0; z < operatorDim*DIM_NODE*nbNodes; z++)
                        if (prec[z] >= 2.1389718111366E-003 &&
                            prec[z] <= 2.1389718111367E-003)
                            printf ("%d : valeur %f atteinte pour z = %d\n", mpiRank,
                                    prec[z], z); //68086
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
        tag   = mpiRank + 100;
//        printf ("%d :\t%d\t%d\t%d\t%d\t%d\n", mpiRank, node1, node2, size, dest, tag);
        MPI_Send (&(bufferSend[node1]), size, MPI_DOUBLE_PRECISION, dest, tag,
                  MPI_COMM_WORLD);
    }

    // Waiting incoming data
    node2 = 0;
    for (int i = 0; i < nbIntf; i++) {
        node1 = node2;
        node2 = intfIndex[i+1];
        MPI_Wait (&(neighborList[2*nbIntf+i]), MPI_STATUS_IGNORE);
    }

    // Assembling local and incoming data
    node2 = 0;
    for (int i = 0; i < nbIntf; i++) {
        node1  = node2;
        node2  = intfIndex[i+1];
        for (int j = 0; j < operatorDim; j++) {
            for (int k = node1; k < node2; k++) {
                prec[intfNodes[k]*operatorDim+j] +=
                    bufferRecv[node1*operatorDim+j];
            }
        }
    }
}
