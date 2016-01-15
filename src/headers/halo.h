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

#ifndef HALO_H
#define HALO_H

#ifdef GASPI
    #include <GASPI.h>
#endif
#include <DC.h>

#ifdef MULTITHREADED_COMM
    // Structure containing the user arguments passed to multithreaded communication
    // function
    typedef struct userCommArgs_s {
        double *prec, *srcDataSegment;
        int *srcOffsetSegment, *neighborsList, *intfIndex, *intfDstIndex;
        int nbBlocks, nbIntf, nbMaxComm, operatorDim, rank, iter;
        const gaspi_segment_id_t srcDataSegmentID, dstDataSegmentID,
                                 srcOffsetSegmentID, dstOffsetSegmentID;
        const gaspi_queue_id_t queueID;
    } userCommArgs_t;
#endif

#ifdef XMPI

// Halo exchange between MPI ranks
void MPI_halo_exchange (double *prec, int *intfIndex, int *intfNodes,
                        int *neighborsList, int nbBlocks, int nbIntf, int nbIntfNodes,
                        int operatorDim, int rank);

#elif GASPI

// Halo exchange between GASPI ranks
void GASPI_halo_exchange (double *prec, double *srcDataSegment,
                          double *dstDataSegment, int *intfIndex, int *intfNodes,
                          int *neighborsList, int *intfDstIndex, int nbBlocks,
                          int nbIntf, int operatorDim, int rank, int iter,
                          const gaspi_segment_id_t srcDataSegmentID,
                          const gaspi_segment_id_t dstDataSegmentID,
                          const gaspi_queue_id_t queueID);

#endif
#ifdef MULTITHREADED_COMM

// Wait for multithreaded GASPI notifications
void GASPI_multithreaded_wait (double *prec, double *dstDataSegment, int *intfNodes,
                               int *dstOffsetSegment, int nbNotifications,
                               int nbBlocks, int operatorDim, int rank,
                               gaspi_segment_id_t dstOffsetSegmentID);

// Send initialized parts of the preconditioner
void GASPI_multithreaded_send (void *userCommArgs, DCcommArgs_t *DCcommArgs);

#endif
#endif
