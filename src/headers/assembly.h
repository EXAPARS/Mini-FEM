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

#ifndef ASSEMBLY_H
#define ASSEMBLY_H

#include <DC.h>

// Structure containing the user arguments passed to ASM function
typedef struct {
    double *coord, *nodeToNodeValue;
    int *nodeToNodeRow, *nodeToNodeColumn, *elemToNode, *elemToEdge;
    int operatorDim;
} userArgs_t;

#ifdef DC_VEC
// Vectorial version of elasticity assembly on a given element interval
void assembly_ela_vec (void *userArgs, DCargs_t *DCargs);

// Vectorial version of laplacian assembly on a given element interval
void assembly_lap_vec (void *userArgs, DCargs_t *DCargs);
#endif

// Sequential version of elasticity assembly on a given element interval
#if defined (DC) || defined (DC_VEC)
void assembly_ela_seq (void *userArgs, DCargs_t *DCargs);
#else
void assembly_ela_seq (void *userArgs, int firstElem, int lastElem);
#endif

// Sequential version of laplacian assembly on a given element interval
#if defined (DC) || defined (DC_VEC)
void assembly_lap_seq (void *userArgs, DCargs_t *DCargs);
#else
void assembly_lap_seq (void *userArgs, int firstElem, int lastElem);
#endif

#ifdef COLORING
// Iterate over the colors & execute the assembly step on the elements of a same color
// in parallel
void coloring_assembly (userArgs_t *userArgs, int operatorID);
#endif

// Call the appropriate function to perform the assembly step
void assembly (double *coord, double *nodeToNodeValue, int *nodeToNodeRow,
               int *nodeToNodeColumn, int *elemToNode, int *elemToEdge, int nbElem,
               int nbEdges, int operatorDim, int operatorID);

#endif
