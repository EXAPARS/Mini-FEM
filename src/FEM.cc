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
#include <iomanip>
#include <cmath>
#include <DC.h>
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
#include "preconditioner.h"
#include "assembly.h"
#include "FEM.h"

using namespace std;

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

// Check if current assembly results match to the reference version
void check_assembly (double *prec, double *nodeToNodeValue, int nbEdges, int nbNodes,
                     int operatorDim, int nbBlocks, int rank)
{
    double refMatrixNorm, refPrecNorm, MatrixNorm, precNorm;
    read_ref_assembly (&refMatrixNorm, &refPrecNorm, nbBlocks, rank);
    MatrixNorm = compute_double_norm (nodeToNodeValue, nbEdges*operatorDim);
    precNorm   = compute_double_norm (prec, nbNodes*operatorDim);

    if (rank == 0) {
        cout << "Numerical stability\n" << setprecision (3);
        cout << "-----------------------------------------\n";
        cout << "  Matrix -> reference norm : " << refMatrixNorm << endl
             << "              current norm : " << MatrixNorm << endl
             << "                difference : " << abs (refMatrixNorm-MatrixNorm)
             << endl << endl;
        cout << "    Prec -> reference norm : " << refPrecNorm << endl
             << "              current norm : " << precNorm << endl
             << "                difference : " << abs (refPrecNorm - precNorm)
             << endl;
        cout << "-----------------------------------------\n";
    }
}

// Main loop iterating over the 3 main steps of FEM applications
void FEM_loop (double *prec, double *coord, double *nodeToNodeValue,
               int *nodeToNodeRow, int *nodeToNodeColumn, int *elemToNode,
               int *elemToEdge, int *intfIndex, int *intfNodes, int *neighborList,
               int *checkBounds, int nbElem, int nbNodes, int nbEdges, int nbIntf,
               int nbIntfNodes, int nbIter, int nbBlocks, int rank, int operatorDim,
#ifdef XMPI
               int operatorID)
#elif GASPI
               int operatorID, double *srcSegment, double *destSegment,
               int *destOffset, gaspi_segment_id_t srcSegmentID,
               gaspi_segment_id_t destSegmentID, gaspi_queue_id_t queueID)
#endif
{
    DC_timer ASMtimer, precTimer;
    uint64_t globalASMcycles, globalPrecCycles, localASMcycles, localPrecCycles;

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
    for (int iter = 1; iter <= nbIter; iter++) {
        
        // Matrix assembly
        if (nbIter == 1 || iter > 1) ASMtimer.start_cycles ();
        assembly (coord, nodeToNodeValue, nodeToNodeRow, nodeToNodeColumn, elemToNode,
                  elemToEdge, nbElem, nbEdges, operatorDim, operatorID);
        if (nbIter == 1 || iter > 1) ASMtimer.stop_cycles ();
        if (rank == 0) cout << iter << ". Matrix assembly                   done\n";

        // Preconditioner creation
        if (nbIter == 1 || iter > 1) precTimer.start_cycles ();
        preconditioner (prec, nodeToNodeValue, nodeToNodeRow, nodeToNodeColumn,
                        intfIndex, intfNodes, neighborList, checkBounds, nbNodes,
                        nbBlocks, nbIntf, nbIntfNodes, operatorDim, operatorID, rank
        #ifdef XMPI
                        );
        #elif GASPI
                        , srcSegment, destSegment, destOffset, srcSegmentID,
                        destSegmentID, queueID);
        #endif
        if (nbIter == 1 || iter > 1) precTimer.stop_cycles ();
        if (rank == 0) {
            cout << "   Preconditioner creation           done\n";
            if (iter != nbIter) cout << endl;
            else                cout << "done\n\n";
        }
    }

    // Get the max ASM & prec measures from all ranks
    localASMcycles  = ASMtimer.get_avg_cycles ();
    localPrecCycles = precTimer.get_avg_cycles ();
    #ifdef XMPI
        MPI_Reduce (&localASMcycles, &globalASMcycles, 1, MPI_UINT64_T, MPI_MAX, 0,
                    MPI_COMM_WORLD);
        MPI_Reduce (&localPrecCycles, &globalPrecCycles, 1, MPI_UINT64_T, MPI_MAX, 0,
                    MPI_COMM_WORLD);
    #elif GASPI
        gaspi_allreduce (&localASMcycles, &globalASMcycles, 1, GASPI_OP_MAX,
                         GASPI_TYPE_ULONG, GASPI_GROUP_ALL, GASPI_BLOCK);
        gaspi_allreduce (&localPrecCycles, &globalPrecCycles, 1, GASPI_OP_MAX,
                         GASPI_TYPE_ULONG, GASPI_GROUP_ALL, GASPI_BLOCK);
    #endif
    if (rank == 0) {
        cout << "Average cycles\n";
        cout << "-----------------------------------------\n";
        cout << "  Matrix assembly         : " << globalASMcycles << endl;
        cout << "  Preconditioner creation : " << globalPrecCycles << endl;
        cout << "-----------------------------------------\n\n";
    }

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
