# =====================================================================
# \file CMakeLists.txt
# \author V Kotteda 
# \date 2021-02-04
# =====================================================================

# my cmake file for teton-project
cmake_minimum_required(VERSION 3.0) 
#SET(CMAKE_SYSTEM_NAME Generic)
# =====================================================================
# Cache variables
# =====================================================================

# set paths
SET(PETSC_DIR /home/vkotteda/Software/libraries/petsc-gcc7)
SET(PETSC_ARCH linux-gnu-release-cuda)
SET(BOOST_DIR /home/vkotteda/Software/libraries/boost/gcc7/install/)
SET(METIS_DIR /home/vkotteda/Software/libraries/parmetis/install/) 
SET(MPI_DIR /home/vkotteda/Software/libraries/mpi/gcc7/aware-c/) 
#SET(AMGX_DIR /home/exx/Downloads/AMGX-main)


# set default c++ compilers
SET(CMAKE_CXX_COMPILER ${MPI_DIR}/bin/mpicxx "-Wall -fno-omit-frame-pointer")
SET(CMAKE_C_COMPILER ${MPI_DIR}/bin/mpicc)

# set default flags for release and debug mode
# specify the C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)
SET(CMAKE_CXX_FLAGS -DUSE_DOUBLES=1)

# =====================================================================
# Project settings
# =====================================================================

MESSAGE("")
MESSAGE("====================================")
MESSAGE("Begin configuration")
MESSAGE("====================================")

# start generator
project(teton)

# other cmake setting
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_COLOR_MAKEFILE ON)

#find_package(MPI REQUIRED)

# =====================================================================
#SET src directory path 
SET(SRC ${CMAKE_SOURCE_DIR})

#set list of files to compile 
set(SOURCE_EXE
	single_phase3d_openpnm_network.cpp
	linear_solver_petsc.cpp)

#include header files 
include_directories(".")
#INCLUDE(GNUInstallDirs)
include_directories(${PETSC_DIR}/include
	${PETSC_DIR}/${PETSC_ARCH}/include
	${BOOST_DIR}/include
        /usr/local/include
        src/
        src/operations/
        src/pore_network/
        src/misc/
        src/simulation/)

#For the static library:
set(PROJECT_LINK_LIBS libmetis.a)
link_directories(${METIS_DIR}/lib/)

set(PROJECT_LINK_LIBS libboost_system.a  libboost_filesystem.a)
link_directories(${BOOST_DIR}/lib/)

set(PROJECT_LINK_LIBS libpetsc.so) 
LINK_DIRECTORIES(${PETSC_DIR}/${PETSC_ARCH}/lib/)

set(PROJECT_LINK_LIBS libmpi.so)
LINK_DIRECTORIES(${MPI_DIR}/lib/) 

#set(PROJECT_LINK_LIBS libamgx.a     libamgxsh.so)
#LINK_DIRECTORIES(/home/vkotteda/Software/libraries/Nvidia-AmgX/lib/) 

#added for Amgx 
#LINK_DIRECTORIES(
#    ${CUDA_LIBRARY_DIRS}
#    ${AMGX_LIBRARY_DIRS})

ADD_EXECUTABLE(teton ${SOURCE_EXE})

TARGET_LINK_LIBRARIES(teton
	-lpetsc 
	-lmetis 
	-lmpi
	-lboost_filesystem -lboost_system)

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

