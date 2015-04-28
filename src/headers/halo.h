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

#ifdef XMPI
// Halo exchange between MPI ranks
void MPI_halo_exchange (double *prec, int *intfIndex, int *intfNodes,
                        int *neighborList, int nbNodes, int nbBlocks, int nbIntf,
                        int nbIntfNodes, int operatorDim, int operatorID, int rank);

#elif GASPI

// Halo exchange between GASPI ranks
void GASPI_halo_exchange (double *prec, double *srcSegment, double *destSegment,
                          int *intfIndex, int *intfNodes, int *neighborList,
                          int *destOffset, int nbNodes, int nbBlocks, int nbIntf,
                          int nbIntfNodes, int operatorDim, int operatorID, int rank,
                          const gaspi_segment_id_t srcSegmentID,
                          const gaspi_segment_id_t destSegmentID,
                          const gaspi_queue_id_t queueID);
#endif

#endif
