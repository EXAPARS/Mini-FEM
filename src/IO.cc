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

#include <iostream>
#include <iomanip>
#include <fstream>

#include "globals.h"
#include "FEM.h"
#include "IO.h"

// Read reference norm of matrix & preconditioner arrays
void read_ref_assembly (double *refMatrixNorm, double *refPrecNorm, int nbBlocks,
                        int rank)
{
	string fileName = (string)DATA_PATH + "/" + meshName + "/checkings/" + operatorName
                      + "_" + to_string ((long long)nbBlocks) + "_"
                      + to_string ((long long)rank);
    ifstream refASM (fileName, ios::in);
    if (!refASM.is_open ()) {
        cerr << "Error: cannot read reference checking: " << fileName << ".\n";
        exit (EXIT_FAILURE);
    }
    refASM >> *refMatrixNorm >> *refPrecNorm;
    refASM.close ();
}

// Store reference norm of matrix & preconditioner arrays 
void store_ref_assembly_ (double *refMatrix, double *refPrec, int *nbEdges,
                          int *nbNodes, int *operatorDim, int *nbBlocks, int *rank)
{
    double refMatrixNorm, refPrecNorm;
    refMatrixNorm = compute_double_norm (refMatrix, (*nbEdges)*(*operatorDim));
    refPrecNorm   = compute_double_norm (refPrec, (*nbNodes)*(*operatorDim));

    string fileName = "checkings_" + to_string ((long long)*nbBlocks) + "_" +
                                     to_string ((long long)*rank);
    ofstream refASM (fileName, ios::out | ios::trunc);
    if (!refASM.is_open ()) {
        cerr << "Error: cannot store reference checking.\n";
        exit (EXIT_FAILURE);
    }
    refASM << setprecision(17) << refMatrixNorm << endl << refPrecNorm << endl;
    refASM.close ();
}

// Read input data from DefMesh
void read_input_data (double **coord, int **elemToNode, int **neighborsList,
                      int **intfIndex, int **intfNodes, int **boundNodesCode,
                      int *nbElem, int *nbNodes, int *nbEdges, int *nbIntf,
                      int *nbIntfNodes, int *nbBoundNodes, int nbBlocks, int rank)
{
    string fileName = (string)DATA_PATH + "/" + meshName + "/inputs/" + operatorName
                      + "_" + to_string ((long long)nbBlocks) + "_"
                      + to_string ((long long)rank);
    ifstream inputFile (fileName, ios::in | ios::binary);
    if (!inputFile.is_open ()) {
        cerr << "Error: cannot read input data: " << fileName << "\n";
        exit (EXIT_FAILURE);
    }

    inputFile.read ((char*)nbElem,       sizeof (int));
    inputFile.read ((char*)nbNodes,      sizeof (int));
    inputFile.read ((char*)nbEdges,      sizeof (int));
    inputFile.read ((char*)nbIntf,       sizeof (int));
    inputFile.read ((char*)nbIntfNodes,  sizeof (int));
    inputFile.read ((char*)nbBoundNodes, sizeof (int));

    *coord          = new double [(*nbNodes) * DIM_NODE];
    *elemToNode     = new int    [(*nbElem)  * DIM_ELEM];
    *neighborsList  = new int    [max (*nbIntf,1) * 3];
    *intfIndex      = new int    [(*nbIntf) + 1];
    *intfNodes      = new int    [*nbIntfNodes];
    *boundNodesCode = new int    [*nbNodes];

    inputFile.read ((char*)*coord,        (*nbNodes) * DIM_NODE * sizeof (double));
    inputFile.read ((char*)*elemToNode,    (*nbElem) * DIM_ELEM * sizeof (int));
    inputFile.read ((char*)*neighborsList, (max(*nbIntf,1) * 3) * sizeof (int));
    inputFile.read ((char*)*intfIndex,        ((*nbIntf) + 1)   * sizeof (int));
    inputFile.read ((char*)*intfNodes,         (*nbIntfNodes)   * sizeof (int));
    inputFile.read ((char*)*boundNodesCode,    (*nbNodes)       * sizeof (int));
    inputFile.close ();
}

// Store necessary data from DefMesh
void store_input_data_ (double *coord, int *elemToNode, int *neighborsList,
                        int *intfIndex, int *intfNodes, int *boundNodesCode,
                        int *nbElem, int *dimElem, int *nbNodes, int *dimNode,
                        int *nbEdges, int *nbIntf, int *nbIntfNodes, int *nbBoundNodes,
                        int *nbBlocks, int *rank)
{
    string fileName = "inputs_" + to_string ((long long)*nbBlocks) + "_" +
                                  to_string ((long long)*rank);
    ofstream inputFile (fileName, ios::out | ios::trunc | ios::binary);
    if (!inputFile.is_open ()) {
        cerr << "Error: cannot store input data.\n";
        exit (EXIT_FAILURE);
    }

    inputFile.write ((char*)nbElem,       sizeof (int));
    inputFile.write ((char*)nbNodes,      sizeof (int));
    inputFile.write ((char*)nbEdges,      sizeof (int));
    inputFile.write ((char*)nbIntf,	      sizeof (int));
    inputFile.write ((char*)nbIntfNodes,  sizeof (int));
    inputFile.write ((char*)nbBoundNodes, sizeof (int));

    inputFile.write ((char*)coord,      (*nbNodes) * (*dimNode) * sizeof (double));
    inputFile.write ((char*)elemToNode,  (*nbElem) * (*dimElem) * sizeof (int));
    inputFile.write ((char*)neighborsList, (max(*nbIntf,1) * 3) * sizeof (int));
    inputFile.write ((char*)intfIndex,        ((*nbIntf) + 1)   * sizeof (int));
    inputFile.write ((char*)intfNodes,         (*nbIntfNodes)   * sizeof (int));
    inputFile.write ((char*)boundNodesCode,    (*nbNodes)       * sizeof (int));
    inputFile.close ();
}
