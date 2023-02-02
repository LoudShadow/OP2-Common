# OP2-LULESH

This is the README for the OP2 version of LULESH 2.0

The original source code can be found under https://github.com/LLNL/LULESH and more information on https://asc.llnl.gov/codes/proxy-apps/lulesh

### Before compiling

This repository requires a working version of the OP2 repository found under https://github.com/OP-DSL/OP2-Common
Both the makefile and cmake files require the environment variable OP2_INSTALL_PATH to be set to the root of the OP2 directory.

### Structure

At the moment the structure of this repository is still a work in progress and will be changed in future versions.
The main part of the application can be found in the lulesh.cpp file.
In the original application the following files still served a purpose but have mostly been consolidated in lulesh.cpp.

- lulesh-init.cc
- lulesh-util.cc
- lulesh-comm.cc
- lulesh_tuple.h
- lulesh-viz.cc
- lulesh.h

### Running LULESH

To compile LULESH with OP2 the location for the OP2 installation has to be changed in the Makefile.
Runtime arguments can be found in lulesh.cpp, with balance and regions having no impact at the moment.

### Further Development

LULESH is originally built using Cmake and has not been changed when developing the OP2 developer serial version. When following the "Developing an OP2 application" guide from teh OP2 DSL page, the section on how to include the different headers does not exaclty apply as cmake is quite different from a normal make system. However, the CMakeLists.txt has been adjusted to allow for compiling of the sequential and distributed header files. To do so the following can be done:

- Create a build directory (mkdir build) and cd into it (cd build)
- use "cmake -DCMAKE_BUILD_TYPE=Release .." to build the project
- Additional CMake variables can be supplied

  CMAKE_CXX_COMPILER Path to the C++ compiler (usually found automaticlly)
  MPI_CXX_COMPILER Path to the MPI C++ compiler (usually found automaticlly)

  WITH_MPI=On|Off Build with MPI (Default: On)
  WITH_OPENMP=On|Off Build with OpenMP support (Default: Off)
  WITH_SILO=On|Off Build with support for SILO. (Default: Off).
  WITH_OP2=On|Off Build with OP2 (Default: On and should obviously stay on)
  OP2_SEQ=On|Off Build with op2_seq header for developer serial. If off it will build with op_lib_mpi.h and allows developer distrubted version to be built

  SILO_DIR Path to SILO library (only needed when WITH_SILO is "On")

When developing the distributed version the MPI flag should stay on as some code for setting the variable which holds the number of ranks is used for initialisation and only set when USE_MPI flag is true. This might not be required at a later point.

Beyond that the options for silo doesnt do much as it has been uncommented in the application code.

When compiling the application the CUDA version seems to create some problems where it does not recognize any variable that is defined with #define and should be replaced by the compiler. This should be looked into in future versions.
