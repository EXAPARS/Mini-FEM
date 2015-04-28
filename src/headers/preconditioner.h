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

#ifndef PRECONDITIONER_H
#define PRECONDITIONER_H

// External elasticity preconditioner Fortran functions
extern "C"
void ela_invert_prec_ (int *dimNode, int *nbNodes, int *nodeToNodeIndex,
                       int *nodeToNodeValue, double *prec, int *error,
                       int *checkBounds, int *curNode);

// Inversion of the preconditioner
void prec_inversion (double *prec, int *nodeToNodeRow, int *nodeToNodeColumn,
                     int *checkBounds, int nbNodes, int operatorID);

// Reset & initialization of the preconditioner
void prec_init (double *prec, double *nodeToNodeValue, int *nodeToNodeRow,
                int *nodeToNodeColumn, int nbNodes, int operatorDim, int operatorID);

#endif
