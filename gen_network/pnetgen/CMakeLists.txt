# =====================================================================
# \file CMakeLists.txt
# \author V Kotteda 
# \date 2021-02-04
# =====================================================================

cmake_minimum_required(VERSION 3.0) 
# set paths
#SET(MPI_DIR /home/vkotteda/Software/libraries/mpi/gcc7/aware-c/) 

# set default c++ compilers
#SET(CMAKE_CXX_COMPILER ${MPI_DIR}/bin/mpicxx)
#SET(CMAKE_C_COMPILER ${MPI_DIR}/bin/mpicc)

# set default flags for release and debug mode
# specify the C++ standard
#set(CMAKE_CXX_STANDARD 17)
#set(CMAKE_CXX_STANDARD_REQUIRED True)
enable_language (Fortran)
# =====================================================================
set(CMAKE_Fortran_FLAGS "-g -O3") 
# Project settings
# =====================================================================

MESSAGE("")
MESSAGE("====================================")
MESSAGE("Begin configuration")
MESSAGE("====================================")

# start generator
project(pnetgen)

# other cmake setting
SET(CMAKE_VERBOSE_MAKEFILE ON)
SET(CMAKE_COLOR_MAKEFILE ON)

#find_package(MPI REQUIRED)
# =====================================================================
#SET src directory path 
SET(SRC ${CMAKE_SOURCE_DIR}/src)

#set list of files to compile 
set(SOURCE_EXE
    gen_network_pore.f90
    ${SRC}/sdot.f   
    ${SRC}/spofa.f 
    ${SRC}/advnst.f
    ${SRC}/genbet.f
    ${SRC}/genchi.f
    ${SRC}/genexp.f
    ${SRC}/genf.f
    ${SRC}/gengam.f
    ${SRC}/genmn.f
    ${SRC}/genmul.f
    ${SRC}/gennch.f
    ${SRC}/gennf.f
    ${SRC}/gennor.f
    ${SRC}/genprm.f
    ${SRC}/genunf.f
    ${SRC}/getcgn.f
    ${SRC}/getsd.f
    ${SRC}/ignbin.f
    ${SRC}/ignlgi.f
    ${SRC}/ignnbn.f
    ${SRC}/ignpoi.f
    ${SRC}/ignuin.f
    ${SRC}/initgn.f
    ${SRC}/inrgcm.f
    ${SRC}/lennob.f
    ${SRC}/mltmod.f
    ${SRC}/phrtsd.f
    ${SRC}/qrgnin.f
    ${SRC}/ranf.f
    ${SRC}/setall.f
    ${SRC}/setant.f
    ${SRC}/setgmn.f
    ${SRC}/setsd.f
    ${SRC}/sexpo.f
    ${SRC}/sgamma.f
    ${SRC}/snorm.f ) 

include_directories(".")
#INCLUDE(GNUInstallDirs)
include_directories(/usr/local/include
        src) 

ADD_EXECUTABLE(pnetgen ${SOURCE_EXE})

TARGET_LINK_LIBRARIES(pnetgen)
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
MESSAGE("")

MESSAGE("Compilers:")
MESSAGE("\tFortran compiler: " ${CMAKE_Fortran_COMPILER})
MESSAGE("")


