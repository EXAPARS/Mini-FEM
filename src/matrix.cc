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

#include <cilk/cilk.h>

#include "globals.h"
#include "matrix.h"

// Create elem to edge array giving the index of each edge of each element
void create_elemToEdge (int *nodeToNodeRow, int *nodeToNodeColumn, int *elemToNode,
                        int *elemToEdge, int nbElem)
{
    // For each element
    #ifdef OMP
        #pragma omp parallel for
        for (int i = 0; i < nbElem; i++) {
    #else
        cilk_for (int i = 0; i < nbElem; i++) {
    #endif
        int ctr = 0;
        // For each edge of current element
        for (int j = 0; j < DIM_ELEM; j++) {
            int node1 = elemToNode[i*DIM_ELEM+j] - 1;
            for (int k = 0; k < DIM_ELEM; k++) {
                int node2 = elemToNode[i*DIM_ELEM+k] - 1;
                // Get the index of current edge from nodeToNode
                for (int l = nodeToNodeRow[node1];
                         l < nodeToNodeRow[node1+1]; l++) {
                    if (nodeToNodeColumn[l] == (node2 + 1)) {
                        elemToEdge[i*VALUES_PER_ELEM+ctr] = l;
                        ctr++;
                        break;
                    }
                }
            }
        }
    }
}

// Create node to node arrays from node to element and element to node
void create_nodeToNode (int *nodeToNodeRow, int *nodeToNodeColumn,
                        index_t &nodeToElem, int *elemToNode, int nbNodes)
{
    int nodeToNodeCtr = 0;

	// For each node
    for (int i = 0; i < nbNodes; i++) {
        int nbNeighbors = 0,
            nbMaxNeighbors = (nodeToElem.index[i+1] - nodeToElem.index[i]) * DIM_ELEM;
        int *neighbors = new int [nbMaxNeighbors];
        nodeToNodeRow[i] = nodeToNodeCtr;
		// For each neighbor element of current node
        for (int j = nodeToElem.index[i]; j < nodeToElem.index[i+1]; j++) {
            int elemNeighbor = nodeToElem.value[j];
		    // For each node of neighbor element
            for (int k = 0; k < DIM_ELEM; k++) {
                int nodeNeighbor = elemToNode[elemNeighbor*DIM_ELEM+k];
                bool isNew = true;
                // Check if current node neighbor is already stored
                for (int l = 0; l < nbNeighbors; l++) {
                    if (nodeNeighbor == neighbors[l]) {
                        isNew = false;
                    }
                }
				// If current neighbor is not in the list, add it
                if (isNew) {
                    nodeToNodeColumn[nodeToNodeCtr] = nodeNeighbor;
                    neighbors[nbNeighbors] = nodeNeighbor;
                    nodeToNodeCtr++;
                    nbNeighbors++;
                }
            }
        }
        delete[] neighbors;
    }
    nodeToNodeRow[nbNodes] = nodeToNodeCtr;
}
