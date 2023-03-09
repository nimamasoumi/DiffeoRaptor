#-----------------------------------------------------------------------------
# Top level CMakeLists.txt file for PyCA
CMAKE_MINIMUM_REQUIRED(VERSION 2.6)
PROJECT(DIFF)

SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

if(COMMAND cmake_policy)
  cmake_policy(SET CMP0003 NEW)
endif(COMMAND cmake_policy)

INCLUDE_DIRECTORIES(
  ${CMAKE_CURRENT_SOURCE_DIR}/
  )

#-----------------------------------------------------------------------------
# set cmake module path
SET(DIFF_CMAKE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/CMake")
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${DIFF_CMAKE_DIR}")

# looks for PyCA bin directory
FIND_PATH(PyCA_DIR PyCAConfig.cmake)

IF(PyCA_DIR)
  INCLUDE(${PyCA_DIR}/PyCAConfig.cmake)
  INCLUDE_DIRECTORIES(${PyCA_INCLUDE_DIRECTORIES})
ELSE(PyCA_DIR)
  MESSAGE("PyCA_DIR (binary dir for compiled PyCA library) required but missing")
ENDIF(PyCA_DIR)

# == Setup PyCA Include Directories (files missing in PyCA_DIR): ==
set(PyCA_SOURCE "" CACHE FILEPATH "PyCA library path")
LINK_DIRECTORIES( ${PyCA_SOURCE} )
set(PyCA_SOURCE "" CACHE FILEPATH "PyCA header file path")
LIST(APPEND PyCA_INCLUDE_DIRECTORIES
     ${PyCA_SOURCE}/Code/Cxx/inc/types
     ${PyCA_SOURCE}/Code/Cxx/inc/math
     ${PyCA_SOURCE}/Code/Cxx/inc/alg
     ${PyCA_SOURCE}/Code/Cxx/inc/io
     ${PyCA_SOURCE}/Code/Cxx/src/math
)
INCLUDE_DIRECTORIES(${PyCA_INCLUDE_DIRECTORIES})

#-----------------------------------------------------------------------------
# Setup FFTW
# * include directories are added automatically
# * applications must add link libraries ${FFTWF_LIB} or ${FFTWF_THREADS_LIB}
MARK_AS_ADVANCED(USE_FFTWF)
OPTION(USE_FFTWF "Use single precision FFTW if found" ON)
MARK_AS_ADVANCED(USE_FFTWD)
OPTION(USE_FFTWD "Use single precision FFTW if found" ON)
FIND_PACKAGE(FFTW REQUIRED)
INCLUDE_DIRECTORIES(${FFTW_INCLUDE})

#-----------------------------------------------------------------------------
# Setup ITK
# * include directories are added automatically
# * applications must add appropriate libraries to link libraries
#MARK_AS_ADVANCED(USE_ITK)
#OPTION(USE_ITK "Use ITK if found" ON)
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  INCLUDE (${ITK_USE_FILE})
ELSE(ITK_FOUND)
  MESSAGE("Could not find ITK directory.")
ENDIF(ITK_FOUND)

#-----------------------------------------------------------------------------
# == Using Armadillo: ==
FIND_PACKAGE(Armadillo REQUIRED)
IF (Armadillo_FOUND)
    INCLUDE_DIRECTORIES(${ARMADILLO_INCLUDE_DIRS})
    MESSAGE(${ARMADILLO_INCLUDE_DIRS})
ENDIF()

#-----------------------------------------------------------------------------
# == Optional subdirectories to build ==
OPTION(BUILD_TESTING "Build tests?" ON)

# == Setup Include Directories ==
LIST(APPEND Flash_INCLUDE_DIRECTORIES
  ${Flash_SOURCE_DIR}/Cxx/inc/Base
  ${Flash_SOURCE_DIR}/Cxx/inc/Alg
  ${Flash_SOURCE_DIR}/Cxx/inc/
)

INCLUDE_DIRECTORIES(${Flash_INCLUDE_DIRECTORIES})

# == What to Build ==
SUBDIRS(${DIFF_SOURCE_DIR}/Cxx/src/)

IF(BUILD_TESTING)
  SUBDIRS(${DIFF_SOURCE_DIR}/Testing)
ENDIF(BUILD_TESTING)
