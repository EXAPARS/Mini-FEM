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

#include "globals.h"
#include "coloring.h"

#define MAX_COLOR 128

using namespace std;

// Fill the index of elements per color
void fill_color_index (int *colors, int *colorPart, int nbElem, int nbColors,
                       int offset)
{
    int *colorCount = new int [nbColors] ();

    // Count the number of elements of each color
    for (int i = 0; i < nbElem; i++) {
        colorCount[colorPart[i]]++;
    }
    // Create the coloring index
    colors[0] = offset;
    for (int i = 1; i <= nbColors; i++) {
        colors[i] = colors[i-1] + colorCount[i-1];
    }
    delete[] colorCount;
}

// Assign a color to the elements of a given leaf following the longest colors strategy
// & return the number of colors
int create_longest_color_part (int *colorPart, list_t *elemToElem, int nbElem)
{
	__uint128_t *elemToColor = new __uint128_t [nbElem] ();
	int nbColors = 0;

    // For each element of local interval
	for (int i = 0; i < nbElem; i++) {
		__uint128_t mask = 1, neighborColor = 0;
		int color = 0;

        // Get the color of all neigbor elements
        for (int j = 0; j < elemToElem[i].size; j++) {
            int neighbor = elemToElem[i].list[j];
            neighborColor |= elemToColor[neighbor];
		}
        // Get the first free color (position of the first 0 bit)
		while (neighborColor & mask) {
		    neighborColor = neighborColor >> 1;
		    color++;
		}
		if (color >= MAX_COLOR) {
		    cerr << "Error: Not enough colors !\n";
		    exit (EXIT_FAILURE);
        }
        // Assign the first free color to current element
		elemToColor[i] = (mask << color);
        colorPart[i] = color;

        // Compute the total number of colors
        if (color > nbColors) nbColors = color;
	}
    nbColors++;

	delete[] elemToColor;
    return nbColors;
}

// Mesh coloring of the whole mesh (elemToNode)
void coloring_creation (int *elemToNode, int *colorPerm, int nbElem, int nbNodes)
{
    // List the neighbor elements of each node
    index_t nodeToElem;
    nodeToElem.index = new int [nbNodes + 1];
    nodeToElem.value = new int [nbElem * DIM_ELEM];
    DC_create_nodeToElem (nodeToElem, elemToNode, nbElem, DIM_ELEM, nbNodes);

    // List the neighbor elements of each element
    list_t *elemToElem = new list_t [nbElem];
    DC_create_elemToElem (elemToElem, nodeToElem, elemToNode, 0, nbElem-1, DIM_ELEM);
    delete[] nodeToElem.value, delete[] nodeToElem.index;

    // Assign a color to each element
    int *colorPart = new int [nbElem];
    nbTotalColors = create_longest_color_part (colorPart, elemToElem, nbElem);
    delete[] elemToElem;

    // Fill the index of elements per color
    colorToElem = new int [nbTotalColors+1];
    fill_color_index (colorToElem, colorPart, nbElem, nbTotalColors, 0);

    // Create a permutation array to sort the elements per color
    DC_create_permutation (colorPerm, colorPart, nbElem, MAX_COLOR);
    delete[] colorPart;
}
