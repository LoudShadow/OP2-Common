cmake_minimum_required(VERSION 3.0)

project(LULESH CXX)

option(WITH_MPI    "Build LULESH with MPI"          FALSE)
option(WITH_OPENMP "Build LULESH with OpenMP"       FALSE)
option(WITH_SILO   "Build LULESH with silo support" FALSE)
option(WITH_OP2 "Build LULESH with OP2" TRUE)
option(OP2_SEQ "Build with OP2 seq" FALSE)
option(WITH_HDF5 "Build LULESH with HDF5" FALSE)

if(WITH_HDF5)
  find_path(
    HDF5_INCLUDE_DIR
    NAMES hdf5.h
    HINTS /dcs/pg21/u2148990/Documents/diss/op2_deps/hdf5/include
    )
  find_library(
    HDF5_LIBRARY
    NAMES hdf5
    HINTS /dcs/pg21/u2148990/Documents/diss/op2_deps/hdf5/lib)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(HDF5 DEFAULT_MSG
    HDF5_LIBRARY
    HDF5_INCLUDE_DIR)
  if(HDF5_FOUND)
    include_directories(${HDF5_INCLUDE_DIR})
    list(APPEND LULESH_EXTERNAL_LIBS ${HDF5_LIBRARY})
  endif()
endif()
if (WITH_MPI)
  find_package(MPI REQUIRED)
  include_directories(${MPI_C_INCLUDE_PATH} ${MPI_CXX_INCLUDE_PATH})
  add_definitions("-DUSE_MPI=1")
  list(APPEND LULESH_EXTERNAL_LIBS ${MPI_C_LIBRARIES} ${MPI_CXX_LIBRARIES})
else()
  add_definitions("-DUSE_MPI=0")
endif()

if (WITH_OPENMP)
  find_package(OpenMP REQUIRED)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

if (WITH_SILO)
  find_path(SILO_INCLUDE_DIR silo.h
    HINTS ${SILO_DIR}/include)
  find_library(SILO_LIBRARY
    NAMES siloh5
    HINTS ${SILO_DIR}/lib)

  include(FindPackageHandleStandardArgs)
  find_package_handle_standard_args(SILO DEFAULT_MSG
    SILO_LIBRARY
    SILO_INCLUDE_DIR)

  if (SILO_FOUND)
    add_definitions("-DVIZ_MESH")
    include_directories(${SILO_INCLUDE_DIR})
    #   Note: silo needs to be built as a dynamic lib, otherwise
    # there are additional dependencies (hdf5) which we don't know.
    # This would be fixed by silo providing a CMake package.
    list(APPEND LULESH_EXTERNAL_LIBS ${SILO_LIBRARY})
  endif()
endif()
if(WITH_OP2)
  if(OP2_SEQ)
    find_path(
      OP2_INCLUDE_DIR
      NAMES op_seq.h
      HINTS $ENV{OP2_INSTALL_PATH}/op2/include
      )
    find_library(
      OP2_LIBRARY
      NAMES op2_seq
      HINTS $ENV{OP2_INSTALL_PATH}/op2/lib/)
    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(OP2 DEFAULT_MSG
      OP2_LIBRARY
      OP2_INCLUDE_DIR)
  else()
    find_path(
      OP2_INCLUDE_DIR
      NAMES op_lib_mpi.h
      HINTS $ENV{OP2_INSTALL_PATH}/op2/include
      )
    find_library(
      OP2_LIBRARY
      NAMES op2_mpi
      HINTS $ENV{OP2_INSTALL_PATH}/op2/lib/)

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(OP2 DEFAULT_MSG
      OP2_LIBRARY
      OP2_INCLUDE_DIR)
  endif()
  if(OP2_FOUND)
    include_directories(${OP2_INCLUDE_DIR})
    list(APPEND LULESH_EXTERNAL_LIBS ${OP2_LIBRARY})
  endif()
endif()

set(LULESH_SOURCES
  # lulesh-comm.cc
  # lulesh-init.cc
  # lulesh-util.cc
  # lulesh-viz.cc
  lulesh.cpp)

set(LULESH_EXEC lulesh2.0)

add_executable(${LULESH_EXEC} ${LULESH_SOURCES})
target_link_libraries(${LULESH_EXEC} ${LULESH_EXTERNAL_LIBS})
