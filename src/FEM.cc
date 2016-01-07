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
#include <fstream>
#include <iomanip>
#include <string.h>
#include <cmath>
#ifdef CILKVIEW
    #include <cilkview.h>
#elif VTUNE
    #include <ittnotify.h>
#elif PAPI
    #include <papi.h>
#elif INTEL_PCM
    #include <cpucounters.h>
#endif

#include "globals.h"
#include "IO.h"
#include "GASPI_handler.h"
#include "halo.h"
#include "preconditioner.h"
#include "assembly.h"
#include "FEM.h"

// Return the euclidean norm of given array
double compute_double_norm (double *tab, int size)
{
    double norm = 0;
    for (int i = 0; i < size; i++) {
        norm += pow (tab[i], 2);
    }
    norm = sqrt (norm);
    return norm;
}

// Check if current results match to the reference version
void check_results (double *prec, double *nodeToNodeValue, int nbEdges, int nbNodes,
                    int operatorDim, int nbBlocks, int rank)
{
    double refMatrixNorm, refPrecNorm, MatrixNorm, precNorm;
    read_ref_assembly (&refMatrixNorm, &refPrecNorm, nbBlocks, rank);
    MatrixNorm = compute_double_norm (nodeToNodeValue, nbEdges*operatorDim);
    precNorm   = compute_double_norm (prec, nbNodes*operatorDim);

    // Store results in a file
    string fileName = "numerical_results_" + to_string ((long long)rank);
    ofstream resultFile (fileName, ios::out | ios::trunc);
    resultFile << "Numerical stability of rank " << rank << endl
               << "----------------------------------------------" << endl
               << "  Matrix -> reference norm : " << refMatrixNorm << endl
               << "              current norm : " << MatrixNorm << endl
               << "                difference : " << abs (refMatrixNorm - MatrixNorm) /
                                                          refMatrixNorm << endl << endl
               << "    Prec -> reference norm : " << refPrecNorm << endl
               << "              current norm : " << precNorm << endl
               << "                difference : " << abs (refPrecNorm - precNorm) /
                                                          refPrecNorm << endl
               << "----------------------------------------------" << endl;
    resultFile.close ();

    // Display results of rank 0
    if (rank == 0) {
        cout << "Numerical stability of rank 0" << endl
             << "----------------------------------------------" << endl
             << "  Matrix -> reference norm : " << refMatrixNorm << endl
             << "              current norm : " << MatrixNorm << endl
             << "                difference : " << abs (refMatrixNorm - MatrixNorm) /
                                                        refMatrixNorm << endl << endl
             << "    Prec -> reference norm : " << refPrecNorm << endl
             << "              current norm : " << precNorm << endl
             << "                difference : " << abs (refPrecNorm - precNorm) /
                                                        refPrecNorm << endl
             << "----------------------------------------------" << endl
             << "(see numerical_results files for all ranks)" << endl << endl;
    }
}

// Get the average measures from all ranks and keep the max
void get_average_cycles (DC_timer &ASMtimer, DC_timer &precInitTimer,
                         DC_timer &haloTimer, DC_timer &precInverTimer,
                         int nbBlocks, int rank)
{
    uint64_t localCycles[4], globalCycles[4];
    localCycles[0] = ASMtimer.get_avg_cycles ();
    localCycles[1] = precInitTimer.get_avg_cycles ();
    localCycles[2] = haloTimer.get_avg_cycles ();
    localCycles[3] = precInverTimer.get_avg_cycles ();

    if (nbBlocks > 1) {
        #ifdef XMPI
            MPI_Reduce (localCycles, globalCycles, 4, MPI_UINT64_T, MPI_MAX, 0,
                        MPI_COMM_WORLD);
        #elif GASPI
            SUCCESS_OR_DIE (gaspi_allreduce (localCycles, globalCycles, 4,
                                             GASPI_OP_MAX, GASPI_TYPE_ULONG,
                                             GASPI_GROUP_ALL, GASPI_BLOCK));
        #endif
    }
    else {
        memcpy (globalCycles, localCycles, 4 * sizeof (uint64_t));
    }

    if (rank == 0) {
        cout << "Average cycles\n";
        cout << "----------------------------------------------\n";
        cout << "  Matrix assembly               : " << globalCycles[0] << endl;
        cout << "  Preconditioner initialization : " << globalCycles[1] << endl;
        cout << "  Halo exchange                 : " << globalCycles[2] << endl;
        cout << "  Preconditioner inversion      : " << globalCycles[3] << endl;
        cout << "  Total                         : " << globalCycles[0]
                  + globalCycles[1] + globalCycles[2] + globalCycles[3] << endl;
        cout << "----------------------------------------------\n\n";
    }
}

// Main loop iterating over the 3 main steps of FEM applications
void FEM_loop (double *prec, double *coord, double *nodeToNodeValue,
               int *nodeToNodeRow, int *nodeToNodeColumn, int *elemToNode,
               int *elemToEdge, int *intfIndex, int *intfNodes, int *neighborsList,
               int *checkBounds, int nbElem, int nbNodes, int nbEdges, int nbIntf,
               int nbIntfNodes, int nbIter, int nbBlocks, int rank, int operatorDim,
#ifdef XMPI
               int operatorID)
#elif GASPI
               int operatorID, int nbMaxComm, int nbNotifications,
               double *srcDataSegment, double *destDataSegment, int *srcOffsetSegment,
               int *destOffsetSegment, int *intfDestIndex,
               gaspi_segment_id_t srcDataSegmentID,
               gaspi_segment_id_t destDataSegmentID,
               gaspi_segment_id_t srcOffsetSegmentID,
               gaspi_segment_id_t destOffsetSegmentID, gaspi_queue_id_t queueID)
#endif
{
    DC_timer ASMtimer, precInitTimer, haloTimer, precInverTimer;

    #ifdef VTUNE
    	__itt_pause ();
    #elif CILKVIEW
    	cilkview_data_t cilkviewData;
    	__cilkview_query (cilkviewData);
    #elif PAPI
        long long papiL3TCM;
        int event_set = PAPI_NULL, event_code = PAPI_L3_TCM;
        PAPI_library_init (PAPI_VER_CURRENT);
        PAPI_create_eventset (&event_set);
        PAPI_add_event (event_set, event_code);
        PAPI_start (event_set);
    #elif INTEL_PCM
        PCM *m = PCM::getInstance ();
        m->program (PCM::DEFAULT_EVENTS, nullptr);
        SystemCounterState pcmStart = getSystemCounterState ();
    #endif

    // Main FEM loop
    for (int iter = 0; iter < nbIter; iter++) {

        // Matrix assembly for bulk synchronous version + preconditioner initialization
        // and halo sending for multithreaded version
        if (rank == 0) cout << iter << ". Matrix assembly...                ";
        if (nbIter == 1 || iter > 0) ASMtimer.start_cycles ();
        assembly (coord, nodeToNodeValue, nodeToNodeRow, nodeToNodeColumn, elemToNode,
                  elemToEdge, nbElem, nbEdges, operatorDim, operatorID
        #ifdef MULTITHREADED_COMM
                  , prec, srcDataSegment, srcOffsetSegment, neighborsList, intfIndex,
                  intfDestIndex, nbBlocks, nbIntf, nbMaxComm, rank, iter,
                  srcDataSegmentID, destDataSegmentID, srcOffsetSegmentID,
                  destOffsetSegmentID, queueID
        #endif
                  );
        if (nbIter == 1 || iter > 0) ASMtimer.stop_cycles ();
        if (rank == 0) cout << "done\n";

        // Preconditioner initialization
        #ifdef BULK_SYNCHRONOUS
            if (rank == 0) cout << "   Preconditioner initialization...  ";
            if (nbIter == 1 || iter > 0) precInitTimer.start_cycles ();
            prec_init (prec, nodeToNodeValue, nodeToNodeRow, nodeToNodeColumn,
                       nbNodes, operatorDim);
            if (nbIter == 1 || iter > 0) precInitTimer.stop_cycles ();
            if (rank == 0) cout << "done\n";
        #endif

        // Distributed communications
        if (rank == 0) cout << "   Halo exchange...                  ";
        if (nbIter == 1 || iter > 0) haloTimer.start_cycles ();
        #ifdef MULTITHREADED_COMM
            // Wait for multithreaded GASPI notifications sent during assembly step
            GASPI_multithreaded_wait (prec, destDataSegment, intfNodes,
                                      destOffsetSegment, nbNotifications, nbBlocks,
                                      operatorDim, iter, destOffsetSegmentID, rank);
        #else
            // Halo exchange
            #ifdef XMPI
                MPI_halo_exchange (prec, intfIndex, intfNodes, neighborsList, nbBlocks,
                                   nbIntf, nbIntfNodes, operatorDim, rank);
            #elif GASPI
                GASPI_halo_exchange (prec, srcDataSegment, destDataSegment, intfIndex,
                                     intfNodes, neighborsList, intfDestIndex,
                                     nbBlocks, nbIntf, operatorDim, rank, iter,
                                     srcDataSegmentID, destDataSegmentID, queueID);
            #endif
        #endif
        if (nbIter == 1 || iter > 0) haloTimer.stop_cycles ();
        if (rank == 0) cout << "done\n";

        // Preconditioner inversion
        if (rank == 0) cout << "   Preconditioner inversion...       ";
        if (nbIter == 1 || iter > 0) precInverTimer.start_cycles ();
        prec_inversion (prec, nodeToNodeRow, nodeToNodeColumn, checkBounds, nbNodes,
                        operatorID);
        if (nbIter == 1 || iter > 0) precInverTimer.stop_cycles ();
        if (rank == 0) cout << "done\n\n";

        #ifdef GASPI
            // If queue is half full, wait for previous comms
            GASPI_wait_for_queue_half_full (queueID, rank);

            // Double buffering flip/flop
            gaspi_segment_id_t tmpSegmentID;
            tmpSegmentID        = srcDataSegmentID;
            srcDataSegmentID    = destDataSegmentID;
            destDataSegmentID   = tmpSegmentID;
            tmpSegmentID        = srcOffsetSegmentID;
            srcOffsetSegmentID  = destOffsetSegmentID;
            destOffsetSegmentID = tmpSegmentID;
        #endif
    }

    // Print the average measures
    get_average_cycles (ASMtimer, precInitTimer, haloTimer, precInverTimer, nbBlocks,
                        rank);

    #ifdef VTUNE
    	__itt_resume ();
    #elif CILKVIEW
        __cilkview_do_report (&cilkviewData, nullptr, "MiniFEM: ASM + Prec",
                              CV_REPORT_WRITE_TO_LOG);
    #elif PAPI
        PAPI_stop (event_set, &papiL3TCM);
        if (rank == 0) {
            cout << "\nPAPI counters\n  L3 cache misses : " << papiL3TCM << "\n\n";
        }
    #elif INTEL_PCM
        SystemCounterState pcmStop = getSystemCounterState ();
        cout << "\nPCM counters"
             << "\n  Cycles              : " << getCycles (pcmStart, pcmStop)
             << "\n  IPC                 : " << getIPC (pcmStart, pcmStop)
             << "\n  L2 cache misses     : " << getL2CacheMisses (pcmStart, pcmStop)
             << "\n  L3 cache misses     : " << getL3CacheMisses (pcmStart, pcmStop)
             << "\n  Bytes read from MC  : " << getBytesReadFromMC (pcmStart, pcmStop)
             << "\n  Bytes written to MC : " << getBytesWrittenToMC (pcmStart, pcmStop)
             << endl;
        m->cleanup ();
    #endif
}
