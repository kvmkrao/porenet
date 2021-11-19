# =====================================================================
# \file CMakeLists.txt
# \author V Kotteda 
# \date 2021-02-04
# =====================================================================

# my cmake file for teton-project
cmake_minimum_required(VERSION 3.0) 

# set paths
SET(PETSC_DIR /home/vkotteda/Software/libraries/petsc-gcc7)
SET(PETSC_ARCH linux-gnu-release-cuda)
#SET(MPI_DIR /home/vkotteda/Software/libraries/mpi/gcc7/aware-c/) 


# set default c++ compilers
#SET(CMAKE_CXX_COMPILER ${MPI_DIR}/bin/mpicxx "-g -Wall -fno-omit-frame-pointer ")
#SET(CMAKE_C_COMPILER ${MPI_DIR}/bin/mpicc)

# set default flags for release and debug mode
# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# =====================================================================
# Project settings
# =====================================================================

MESSAGE("")
MESSAGE("====================================")
MESSAGE("Begin configuration")
MESSAGE("====================================")

# start generator
project(sinphase)

# other cmake setting
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_COLOR_MAKEFILE ON)

#find_package(MPI REQUIRED)

# =====================================================================
#SET src directory path 
SET(SRC ${CMAKE_SOURCE_DIR})

#set list of files to compile 
set(SOURCE_EXE
	sphase_unsteady_incomp.cpp
        write_vtk.cpp	
	linear_solver_petsc.cpp)

#include header files 
include_directories(".")
include_directories(${PETSC_DIR}/include
	${PETSC_DIR}/${PETSC_ARCH}/include
        /usr/local/include
        )

#For the static library:
set(PROJECT_LINK_LIBS libpetsc.so) 
LINK_DIRECTORIES(${PETSC_DIR}/${PETSC_ARCH}/lib/)

#set(PROJECT_LINK_LIBS libmpi.so)
#LINK_DIRECTORIES(${MPI_DIR}/lib/) 

ADD_EXECUTABLE(sinphase ${SOURCE_EXE})

TARGET_LINK_LIBRARIES(sinphase
	-lpetsc )

# =====================================================================
# Print information
# =====================================================================

MESSAGE("")
MESSAGE("====================================")
MESSAGE("Config Information:")
MESSAGE("====================================")
MESSAGE("")

MESSAGE("Build type: " ${CMAKE_BUILD_TYPE})
MESSAGE("Build shared library: " ${BUILD_SHARED_LIBS})
MESSAGE("Doxygen available: " ${DOXYGEN_FOUND})
MESSAGE("Installation path: " ${CMAKE_INSTALL_PREFIX})
MESSAGE("")

MESSAGE("Compilers:")
MESSAGE("\tCXX compiler: " ${CMAKE_CXX_COMPILER})
MESSAGE("\tC compiler: " ${CMAKE_C_COMPILER})
MESSAGE("")

MESSAGE("PETSc:")
MESSAGE("\tPETSc dir: " ${PETSC_DIR})
MESSAGE("")

