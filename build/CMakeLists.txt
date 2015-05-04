#  Copyright 2014 - UVSQ, Dassault Aviation
#  Authors list: Loïc Thébault, Eric Petit
#
#  This file is part of Mini-FEM.
#
#  Mini-FEM is free software: you can redistribute it and/or modify it under the terms
#  of the GNU Lesser General Public License as published by the Free Software
#  Foundation, either version 3 of the License, or (at your option) any later version.
#
#  Mini-FEM is distributed in the hope that it will be useful, but WITHOUT ANY
#  WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A
#  PARTICULAR PURPOSE. See the GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License along with
#  Mini-FEM. If not, see <http://www.gnu.org/licenses/>.

cmake_minimum_required (VERSION 2.6)

# Compilers
set (EXECUTABLE_OUTPUT_PATH ../bin)
set (CMAKE_Fortran_COMPILER ifort)
set (CMAKE_CXX_COMPILER icpc)

# Project
project (Mini-FEM)
enable_language (Fortran)

# Flags
set (flags "-O2 -mkl:sequential -DDATA_PATH=\\\"${DATA_PATH}\\\"")
add_definitions (-D${DISTRI} -D${SHARED} -D${VERSION})
if (${DISTRI} STREQUAL "XMPI")
    find_package (MPI)
elseif (${DISTRI} STREQUAL "GASPI")
    set (flags "${flags} -libverbs -lGPI2")
endif (${DISTRI} STREQUAL "XMPI")
if (${SHARED} STREQUAL "CILK")
	set (flags "${flags} -lcilkrts")
elseif (${SHARED} STREQUAL "OMP")
    set (flags "${flags} -openmp")
endif (${SHARED} STREQUAL "CILK")
if (${ARCHI} STREQUAL "SSE")
    set (flags "${flags} -msse")
    add_definitions (-DVEC_SIZE=2)
elseif (${ARCHI} STREQUAL "AVX")
    set (flags "${flags} -mavx")
    add_definitions (-DVEC_SIZE=4)
elseif (${ARCHI} STREQUAL "MIC")
    set (flags "${flags} -mmic")
    add_definitions (-DVEC_SIZE=8)
endif (${ARCHI} STREQUAL "SSE")
if (${TREE})
    add_definitions (-DTREE_CREATION)
endif (${TREE})
if (${OPTIMIZED})
	add_definitions (-DOPTIMIZED)
endif (${OPTIMIZED})
if (${DEBUG})
    set (flags "${flags} -g -O0")
endif (${DEBUG})
if (${VERBOSE})
    set (CMAKE_VERBOSE_MAKEFILE 1)
    set (flags "${flags} -vec-report6")
endif (${VERBOSE})
if (${MATT})
	set (flags "${flags} -finstrument-functions")
endif (${MATT})
if (${PAPI})
	set (flags "${flags} -lpapi -DPAPI")
endif (${PAPI})
if (${VTUNE})
	add_definitions (-DVTUNE)
endif (${VTUNE})
if (${CILKVIEW})
	add_definitions (-DCILKVIEW)
endif (${CILKVIEW})
set (CMAKE_CXX_FLAGS "${flags} -std=c++11")
set (CMAKE_Fortran_FLAGS "${flags} -fpp -convert big_endian")
set (CMAKE_EXE_LINKER_FLAGS "${flags} -nofor-main -lstdc++")

# Headers
include_directories (../src/headers)
include_directories (${DC_INCLUDE})
if (${DISTRI} STREQUAL "XMPI")
    include_directories (${MPI_INCLUDE_PATH})
endif (${DISTRI} STREQUAL "XMPI")
if (${VTUNE})
    include_directories (${VTUNE_INCLUDE})
endif (${VTUNE})
if (${CILKVIEW})
    include_directories (${CILKVIEW_INCLUDE})
endif (${CILKVIEW})

# Sources
file (GLOB_RECURSE source_list ../src/*)

# Executable
set (exec miniFEM_${VERSION}_${DISTRI}_${SHARED})
if (${VERSION} STREQUAL "DC_HYBRID" OR ${ARCHI} STREQUAL "MIC")
    set (exec ${exec}_${ARCHI})
endif (${VERSION} STREQUAL "DC_HYBRID" OR ${ARCHI} STREQUAL "MIC")
if (${TREE})
    set (exec ${exec}_TreeCreation)
endif (${TREE})
if (${OPTIMIZED})
    set (exec ${exec}_Opti)
endif (${OPTIMIZED})
if (${DEBUG})
    set (exec ${exec}_Debug)
endif (${DEBUG})
if (${MATT})
    set (exec ${exec}_MATT)
endif (${MATT})
if (${PAPI})
    set (exec ${exec}_PAPI)
endif (${PAPI})
if (${VTUNE})
    set (exec ${exec}_VTune)
endif (${VTUNE})
if (${CILKVIEW})
    set (exec ${exec}_CilkView)
endif (${CILKVIEW})
add_executable (
    ${exec}
    ${source_list}
)

# Link
set_property (TARGET ${exec} PROPERTY LINKER_LANGUAGE Fortran)
target_link_libraries (
    ${exec}
    ${DC_LIBRARIES}
)
if (${DISTRI} STREQUAL "XMPI")
    target_link_libraries (
        ${exec}
        ${MPI_LIBRARIES}
    )
elseif (${DISTRI} STREQUAL "GASPI")
    target_link_libraries (
        ${exec}
        ${GASPI_LIBRARIES}
    )
endif (${DISTRI} STREQUAL "XMPI")
if (${TREE})
    target_link_libraries (
        ${exec}
        ${METIS_LIBRARIES}
    )
endif (${TREE})
if (${VTUNE})
    target_link_libraries (
        ${exec}
        ${VTUNE_LIBRARIES}
    )
endif (${VTUNE})
