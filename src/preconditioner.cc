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
#include <DC.h>
#include <iostream>

#include "globals.h"
#include "halo.h"
#include "preconditioner.h"

// Create preconditioner for elasticity operator
void preconditioner_ela (double *prec, double *nodeToNodeValue, int *nodeToNodeRow,
                         int *nodeToNodeColumn, int *intfIndex, int *intfNodes,
                         int *neighborList, int *checkBounds, int nbNodes,
                         int nbBlocks, int nbIntf, int nbIntfNodes, int operatorDim,
                         int operatorID, int rank
#ifdef XMPI
                         )
#elif GASPI
                         , double *srcSegment, double *destSegment, int *destOffset,
                         gaspi_segment_id_t srcSegmentID,
                         gaspi_segment_id_t destSegmentID,
                         gaspi_queue_id_t queueID)
#endif
{
    // Copy matrix diagonal into preconditioner
    #ifdef REF
        for (int i = 0; i < nbNodes; i++) {
    #else
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbNodes; i++) {
        #elif CILK
            cilk_for (int i = 0; i < nbNodes; i++) {
        #endif
    #endif
        for (int j = nodeToNodeRow[i]; j < nodeToNodeRow[i+1]; j++) {
            if (nodeToNodeColumn[j]-1 == i) {
        	    for (int k = 0; k < operatorDim; k++) {
            	    prec[i*operatorDim+k] = nodeToNodeValue[j*operatorDim+k];
        	    }   
        		break;
        	}
        }
    }

//    if (rank == 0) {
//        cout << "   Prec initialization     : " << globalElapsed << " cycles\n";
//    }

    // Distributed communications
    if (nbBlocks > 1) {
        #ifdef XMPI
            MPI_halo_exchange (prec, intfIndex, intfNodes, neighborList, nbNodes,
                               nbIntf, nbIntfNodes, operatorDim, operatorID, rank);
        #elif GASPI
            GASPI_halo_exchange (prec, intfIndex, intfNodes, neighborList, nbNodes,
                                 nbBlocks, nbIntf, nbIntfNodes, operatorDim,
                                 operatorID, rank, srcSegment, destSegment, destOffset,
                                 srcSegmentID, destSegmentID, queueID);
        #endif
    }

//    if (rank == 0) {
//        cout << "   Halo exchange           : " << globalElapsed << " cycles\n";
//    }

    // Inversion of preconditioner
    int dimNode = DIM_NODE, error;
    #ifdef REF
    	for (int i = 1; i <= nbNodes; i++) {
    #else
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 1; i <= nbNodes; i++) {
        #elif CILK
            cilk_for (int i = 1; i <= nbNodes; i++) {
        #endif
    #endif
        int curNode = i;
        ela_invert_prec_ (&dimNode, &nbNodes, nodeToNodeRow, nodeToNodeColumn,
                          prec, &error, checkBounds, &curNode);
    }

//    if (rank == 0) {
//        cout << "   Prec inversion          : " << globalElapsed << " cycles\n";
//    }
}

// Create preconditioner for laplacian operator
void preconditioner_lap (double *prec, double *nodeToNodeValue, int *nodeToNodeRow,
                         int *nodeToNodeColumn, int *intfIndex, int *intfNodes,
                         int *neighborList, int nbNodes, int nbBlocks, int nbIntf,
                         int nbIntfNodes, int operatorDim, int operatorID, int rank
#ifdef XMPI
                         )
#elif GASPI
                         , double *srcSegment, double *destSegment, int *destOffset,
                         gaspi_segment_id_t srcSegmentID,
                         gaspi_segment_id_t destSegmentID,
                         gaspi_queue_id_t queueID)
#endif
{
	// Copy matrix diagonal into preconditioner
    #ifdef REF
    	for (int i = 0; i < nbNodes; i++) {
    #else
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbNodes; i++) {
        #elif CILK
    	    cilk_for (int i = 0; i < nbNodes; i++) {
        #endif
    #endif
        for (int j = nodeToNodeRow[i]; j < nodeToNodeRow[i+1]; j++) {
    	    if (nodeToNodeColumn[j]-1 == i) {
    		    prec[i] = nodeToNodeValue[j];
    		    break;
    	    }
        }
    }

//    if (rank == 0) {
//        cout << "   Prec initialization     : " << globalElapsed << " cycles\n";
//    }

	// Distributed communications
    if (nbBlocks > 1) {
        #ifdef XMPI
            MPI_halo_exchange (prec, intfIndex, intfNodes, neighborList, nbNodes,
                               nbIntf, nbIntfNodes, operatorDim, operatorID, rank);
        #elif GASPI
            GASPI_halo_exchange (prec, intfIndex, intfNodes, neighborList, nbNodes,
                                 nbBlocks, nbIntf, nbIntfNodes, operatorDim,
                                 operatorID, rank, srcSegment, destSegment, destOffset,
                                 srcSegmentID, destSegmentID, queueID);
        #endif
    }

//    if (rank == 0) {
//        cout << "   Halo exchange           : " << globalElapsed << " cycles\n";
//    }

	// Inversion of preconditioner
    #ifdef REF
    	for (int i = 0; i < nbNodes; i++) {
    #else
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbNodes; i++) {
        #elif CILK
    	    cilk_for (int i = 0; i < nbNodes; i++) {
        #endif
    #endif
	    prec[i] = 1.0 / prec[i];
	}

//    if (rank == 0) {
//        cout << "   Prec inversion          : " << globalElapsed << " cycles\n";
//    }
}

// Call the appropriate function to create the preconditioner
void preconditioner (double *prec, double *nodeToNodeValue, int *nodeToNodeRow,
                     int *nodeToNodeColumn, int *intfIndex, int *intfNodes,
                     int *neighborList, int *checkBounds, int nbNodes, int nbBlocks,
                     int nbIntf, int nbIntfNodes, int operatorDim, int operatorID,
#ifdef XMPI
                     int rank)
#elif GASPI
                     int rank, double *srcSegment, double *destSegment,
                     int *destOffset, gaspi_segment_id_t srcSegmentID,
                     gaspi_segment_id_t destSegmentID, gaspi_queue_id_t queueID)
#endif
{
    // Preconditioner reset
    #ifdef REF
        for (int i = 0; i < nbNodes * operatorDim; i++) prec[i] = 0;
    #else
        #ifdef OMP
            #pragma omp parallel for
            for (int i = 0; i < nbNodes * operatorDim; i++) prec[i] = 0;
        #elif CILK
            cilk_for (int i = 0; i < nbNodes * operatorDim; i++) prec[i] = 0;
        #endif
    #endif

    if (operatorID == 0) {
        preconditioner_lap (prec, nodeToNodeValue, nodeToNodeRow, nodeToNodeColumn,
                            intfIndex, intfNodes, neighborList, nbNodes, nbBlocks,
                            nbIntf, nbIntfNodes, operatorDim, operatorID, rank
        #ifdef XMPI
                            );
        #elif GASPI
                            , srcSegment, destSegment, destOffset, srcSegmentID,
                            destSegmentID, queueID);
        #endif
    }
    else {
        preconditioner_ela (prec, nodeToNodeValue, nodeToNodeRow, nodeToNodeColumn,
                            intfIndex, intfNodes, neighborList, checkBounds, nbNodes,
                            nbBlocks, nbIntf, nbIntfNodes, operatorDim, operatorID,
        #ifdef XMPI
                            rank);
        #elif GASPI
                            rank, srcSegment, destSegment, destOffset, srcSegmentID,
                            destSegmentID, queueID);
        #endif
    }
}
