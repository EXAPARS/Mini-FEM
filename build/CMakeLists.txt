cmake_minimum_required (VERSION 2.6)

# Compilers
set (EXECUTABLE_OUTPUT_PATH ../bin)
set (CMAKE_Fortran_COMPILER ifort)
set (CMAKE_CXX_COMPILER icpc)
set (CMAKE_C_COMPILER icc)

# Project
project (MiniFEM)
enable_language (Fortran)
find_package (MPI)

# Flags
set (flags "-O2 -mkl:sequential -D${VERSION} -DDATA_PATH=\\\"${DATA_PATH}\\\"")
if (${HYBRID})
    add_definitions (-DHYBRID)
endif (${HYBRID})
if (${ARCHI} STREQUAL "SSE")
    set (flags "${flags} -msse -DVEC_SIZE=2")
elseif (${ARCHI} STREQUAL "AVX")
    set (flags "${flags} -mavx -DVEC_SIZE=4")
elseif (${ARCHI} STREQUAL "MIC")
    set (flags "${flags} -mmic -DVEC_SIZE=8")
endif (${ARCHI} STREQUAL "SSE")
if (${TREE})
    add_definitions (-DTREE_CREATION)
endif (${TREE})
if (${OPTIMIZED})
	add_definitions (-DOPTIMIZED)
endif (${OPTIMIZED})
if (${OMP})
	set (flags "${flags} -openmp")
    add_definitions (-DOMP)
endif (${OMP})
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
set (CMAKE_C_FLAGS ${flags})
set (CMAKE_CXX_FLAGS "${flags} -std=c++11")
set (CMAKE_Fortran_FLAGS "${flags} -fpp -convert big_endian")
set (CMAKE_EXE_LINKER_FLAGS "${flags} -nofor-main -lstdc++ -lcilkrts")

# Headers
include_directories (../src/headers)
include_directories (${MPI_INCLUDE_PATH})
include_directories (${DC_INCLUDE})
if (${VTUNE})
	include_directories (${VTUNE_INCLUDE})
endif (${VTUNE})
if (${CILKVIEW})
	include_directories (${CILKVIEW_INCLUDE})
endif (${CILKVIEW})

# Sources
file (GLOB_RECURSE source_list ../src/*)

# Executable
set (exec miniFEM_${VERSION})
if (${HYBRID})
        set (exec ${exec}_Hybrid_${ARCHI})
endif (${HYBRID})
if (${TREE})
	set (exec ${exec}_TreeCreation)
endif (${TREE})
if (${OPTIMIZED})
	set (exec ${exec}_Opti)
endif (${OPTIMIZED})
if (${OMP})
	set (exec ${exec}_OpenMP)
endif (${OMP})
if (${CILK} AND NOT ${REF})
    set (exec ${exec}_Cilk)
endif (${CILK} AND NOT ${REF})
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
if (${TREE})
    target_link_libraries (
        ${exec}
        ${MPI_LIBRARIES}
        ${DC_LIBRARIES}
        ${METIS_LIBRARIES}
        ${VTUNE_LIBRARIES}
    )
else (${TREE})
    target_link_libraries (
        ${exec}
        ${MPI_LIBRARIES}
        ${DC_LIBRARIES}
        ${VTUNE_LIBRARIES}
    )
endif (${TREE})
