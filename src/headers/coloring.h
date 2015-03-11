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

#ifndef COLORING_H
#define COLORING_H

#include <DC.h>

// Fill the index of elements per color
void fill_color_index (int *colors, int *colorPart, int nbElem, int nbColors,
                       int offset);

// Assign a color to the elements of a given leaf following the longest colors strategy
// & return the number of colors
int create_longest_color_part (int *colorPart, list_t *elemToElem, int nbElem);

// Mesh coloring of the whole mesh (elemToNode)
void coloring_creation (int *elemToNode, int *colorPerm, int nbElem, int nbNodes);

#endif
