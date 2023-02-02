/*
  This is a Version 2.0 MPI + OpenMP implementation of LULESH

                 Copyright (c) 2010-2013.
      Lawrence Livermore National Security, LLC.
Produced at the Lawrence Livermore National Laboratory.
                  LLNL-CODE-461231
                All rights reserved.

This file is part of LULESH, Version 2.0.
Please also read this link -- http://www.opensource.org/licenses/index.php

//////////////
DIFFERENCES BETWEEN THIS VERSION (2.x) AND EARLIER VERSIONS:
* Addition of regions to make work more representative of multi-material codes
* Default size of each domain is 30^3 (27000 elem) instead of 45^3. This is
  more representative of our actual working set sizes
* Single source distribution supports pure serial, pure OpenMP, MPI-only,
  and MPI+OpenMP
* Addition of ability to visualize the mesh using VisIt
  https://wci.llnl.gov/codes/visit/download.html
* Various command line options (see ./lulesh2.0 -h)
 -q              : quiet mode - suppress stdout
 -i <iterations> : number of cycles to run
 -s <size>       : length of cube mesh along side
 -r <numregions> : Number of distinct regions (def: 11)
 -b <balance>    : Load balance between regions of a domain (def: 1)
 -c <cost>       : Extra cost of more expensive regions (def: 1)
 -f <filepieces> : Number of file parts for viz output (def: np/9)
 -p              : Print out progress
 -v              : Output viz file (requires compiling with -DVIZ_MESH
 -h              : This message

 printf("Usage: %s [opts]\n", execname);
      printf(" where [opts] is one or more of:\n");
      printf(" -q              : quiet mode - suppress all stdout\n");
      printf(" -i <iterations> : number of cycles to run\n");
      printf(" -s <size>       : length of cube mesh along side\n");
      printf(" -r <numregions> : Number of distinct regions (def: 11)\n");
      printf(" -b <balance>    : Load balance between regions of a domain (def: 1)\n");
      printf(" -c <cost>       : Extra cost of more expensive regions (def: 1)\n");
      printf(" -f <numfiles>   : Number of files to split viz dump into (def: (np+10)/9)\n");
      printf(" -p              : Print out progress\n");
      printf(" -v              : Output viz file (requires compiling with -DVIZ_MESH\n");
      printf(" -h              : This message\n");
      printf("\n\n");

*Notable changes in LULESH 2.0

* Split functionality into different files
lulesh.cc - where most (all?) of the timed functionality lies
lulesh-comm.cc - MPI functionality
lulesh-init.cc - Setup code
lulesh-viz.cc  - Support for visualization option
lulesh-util.cc - Non-timed functions
*
* The concept of "regions" was added, although every region is the same ideal
*    gas material, and the same sedov blast wave problem is still the only
*    problem its hardcoded to solve.
* Regions allow two things important to making this proxy app more representative:
*   Four of the LULESH routines are now performed on a region-by-region basis,
*     making the memory access patterns non-unit stride
*   Artificial load imbalances can be easily introduced that could impact
*     parallelization strategies.
* The load balance flag changes region assignment.  Region number is raised to
*   the power entered for assignment probability.  Most likely regions changes
*   with MPI process id.
* The cost flag raises the cost of ~45% of the regions to evaluate EOS by the
*   entered multiple. The cost of 5% is 10x the entered multiple.
* MPI and OpenMP were added, and coalesced into a single version of the source
*   that can support serial builds, MPI-only, OpenMP-only, and MPI+OpenMP
* Added support to write plot files using "poor mans parallel I/O" when linked
*   with the silo library, which in turn can be read by VisIt.
* Enabled variable timestep calculation by default (courant condition), which
*   results in an additional reduction.
* Default domain (mesh) size reduced from 45^3 to 30^3
* Command line options to allow numerous test cases without needing to recompile
* Performance optimizations and code cleanup beyond LULESH 1.0
* Added a "Figure of Merit" calculation (elements solved per microsecond) and
*   output in support of using LULESH 2.0 for the 2017 CORAL procurement
*
* Possible Differences in Final Release (other changes possible)
*
* High Level mesh structure to allow data structure transformations
* Different default parameters
* Minor code performance changes and cleanup

TODO in future versions
* Add reader for (truly) unstructured meshes, probably serial only
* CMake based build system

//////////////

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions
are met:

   * Redistributions of source code must retain the above copyright
     notice, this list of conditions and the disclaimer below.

   * Redistributions in binary form must reproduce the above copyright
     notice, this list of conditions and the disclaimer (as noted below)
     in the documentation and/or other materials provided with the
     distribution.

   * Neither the name of the LLNS/LLNL nor the names of its contributors
     may be used to endorse or promote products derived from this software
     without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL LAWRENCE LIVERMORE NATIONAL SECURITY, LLC,
THE U.S. DEPARTMENT OF ENERGY OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT,
INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY
OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE,
EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


Additional BSD Notice

1. This notice is required to be provided under our contract with the U.S.
   Department of Energy (DOE). This work was produced at Lawrence Livermore
   National Laboratory under Contract No. DE-AC52-07NA27344 with the DOE.

2. Neither the United States Government nor Lawrence Livermore National
   Security, LLC nor any of their employees, makes any warranty, express
   or implied, or assumes any liability or responsibility for the accuracy,
   completeness, or usefulness of any information, apparatus, product, or
   process disclosed, or represents that its use would not infringe
   privately-owned rights.

3. Also, reference herein to any specific commercial products, process, or
   services by trade name, trademark, manufacturer or otherwise does not
   necessarily constitute or imply its endorsement, recommendation, or
   favoring by the United States Government or Lawrence Livermore National
   Security, LLC. The views and opinions of authors expressed herein do not
   necessarily state or reflect those of the United States Government or
   Lawrence Livermore National Security, LLC, and shall not be used for
   advertising or product endorsement purposes.

*/

#include <climits>
#include <vector>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <sys/time.h>
#include <iostream>
#include <unistd.h>
#include <iomanip>

#include <op_seq.h>
// #include <hdf5.h>
#include <mpi.h>

#include "lulesh.h"
// #include "const.h"
int myRank;
struct cmdLineOpts opts;
double  m_hgcoef = double(3.0);            // hourglass control
// Cutoffs (treat as constants)
double  m_e_cut = double(1.0e-7);             // energy tolerance
double  m_p_cut = double(1.0e-7) ;             // pressure tolerance
double  m_q_cut = double(1.0e-7) ;             // q tolerance
double  m_v_cut = double(1.0e-10) ;             // relative volume tolerance
double  m_u_cut = double(1.0e-7) ;             // velocity tolerance

// Other constants (usually setable, but hardcoded in this proxy app)

double  m_ss4o3 = (double(4.0)/double(3.0));
double  m_qstop = double(1.0e+12);             // excessive q indicator
double  m_monoq_max_slope = double(1.0);
double  m_monoq_limiter_mult = double(2.0);
double  m_qlc_monoq = double(0.5);         // linear term coef for q
double  m_qqc_monoq = (double(2.0)/double(3.0));         // quadratic term coef for q
double  m_qqc = double(2.0);
double  m_eosvmax = double(1.0e+9);
double  m_eosvmin = double(1.0e-9);
double  m_pmin = double(0.);              // pressure floor
double  m_emin = double(-1.0e+15);              // energy floor
double  m_dvovmax = double(0.1);           // maximum allowable volume change
double  m_refdens = double(1.0);           // reference density

//Consts not defined originally, but for CUDA version
double  m_qqc2 = double(64.0) * m_qqc * m_qqc ;
double  m_ptiny = double(1.e-36) ;
double  m_gamma_t[4*8];
double  m_twelfth = double(1.0)/double(12.0);
double  m_sixth = double(1.0) / double(6.0);
double  m_c1s = double(2.0)/double(3.0);
double  m_ssc_thresh = double(.1111111e-36);
double  m_ssc_low = double(.3333333e-18);

// Identifiers for free ndoes or boundary nodes
#define FREE_NODE 0x00
#define BOUNDARY_NODE 0x01
#define FILE_NAME_PATH "file_out.h5"

// Stuff needed for boundary conditions
// 2 BCs on each of 6 hexahedral faces (12 bits)
#define XI_M        0x00007
#define XI_M_SYMM   0x00001
#define XI_M_FREE   0x00002
#define XI_M_COMM   0x00004

#define XI_P        0x00038
#define XI_P_SYMM   0x00008
#define XI_P_FREE   0x00010
#define XI_P_COMM   0x00020

#define ETA_M       0x001c0
#define ETA_M_SYMM  0x00040
#define ETA_M_FREE  0x00080
#define ETA_M_COMM  0x00100

#define ETA_P       0x00e00
#define ETA_P_SYMM  0x00200
#define ETA_P_FREE  0x00400
#define ETA_P_COMM  0x00800

#define ZETA_M      0x07000
#define ZETA_M_SYMM 0x01000
#define ZETA_M_FREE 0x02000
#define ZETA_M_COMM 0x04000

#define ZETA_P      0x38000
#define ZETA_P_SYMM 0x08000
#define ZETA_P_FREE 0x10000
#define ZETA_P_COMM 0x20000

         /* -t partiotioning library <S,PK,PG,PKG,K> */
#define Partition_PK  4
#define Partition_S   1
#define Partition_PG  2
#define Partition_PGK 3
#define Partition_K   0

#define Creation_Parallel 0
#define Creation_Root 1

#define Hide_Time 0
#define Show_Time 1

// inline double SQRT(double  arg) { return sqrt(arg) ; }
// inline double CBRT(double  arg) { return cbrt(arg) ; }
// inline double FABS(double  arg) { return fabs(arg) ; }
// enum { VolumeError = -1, QStopError = -2 } ;
// #define MAX(a, b) ( ((a) > (b)) ? (a) : (b))

#include "FBHourglassForceForElems.h"
#include "CalcVolumeDerivatives.h"
#include "CheckForNegativeElementVolume.h"
#include "IntegrateStressForElemsLoop.h"
#include "setForceToZero.h"
#include "CalcAccelForNodes.h"
#include "CalcVeloForNodes.h"
#include "CalcPosForNodes.h"
#include "CalcKinematicsForElem.h"
#include "CalcLagrangeElemRemaining.h"
#include "CalcMonotonicQGradientsForElem.h"
#include "CalcHalfStepBVC.h"
#include "CalcPHalfstep.h"
#include "CalcBVC.h"
#include "updateVolumesForElem.h"
#include "NoExcessiveArtificialViscosity.h"
#include "CalcMonotonicQRegionForElem.h"
#include "ALE3DRelevantCheck.h"
#include "CopyEOSValsIntoArray.h"
#include "CalcHalfSteps.h"
#include "CheckEOSLowerBound.h"
#include "CheckEOSUpperBound.h"
#include "CalcEOSWork.h"
#include "CopyTempEOSVarsBack.h"
#include "CopyVelocityToTempArray.h"
#include "ApplyLowerBoundToVelocity.h"
#include "ApplyUpperBoundToVelocity.h"
#include "CalcNewE.h"
#include "CalcNewEStep2.h"
#include "CalcNewEStep3.h"
#include "CalcNewEStep4.h"
#include "CalcQNew.h"
#include "CalcPNew.h"
#include "CalcSoundSpeedForElem.h"
#include "BoundaryX.h"
#include "BoundaryY.h"
#include "BoundaryZ.h"
#include "CalcHydroConstraint.h"
#include "initStressTerms.h"
#include "CalcCourantConstraint.h"
// Define arrays and constants
// Should maybe be moved to a header file to avoid clutter or to the main function to not be global
// Element-centered


static inline
void logVec(char *funcName, op_dat value){
   std::cout << std::scientific << std::setprecision(6);

   printf("%s\n",funcName);
   printf("\n");
   for (int i =0 ; i<64;i++){
      std::cout << i << "     " << std::right << (double)((double *)value->data)[i] << "\n";
   }
   printf("\n");
}


// Region information
int    m_numReg ;
int    m_cost; //imbalance cost
int *m_regElemSize ;   // Size of region sets
int *m_regNumList ;    // Region number per domain element
int **m_regElemlist ;  // region indexset

int* nodelist;

int* lxim; /* element connectivity across each face */
int* lxip;
int* letam;
int* letap;
int* lzetam;
int* lzetap;

int* elemBC;            /* symmetry/free-surface flags for each elem face */

double* dxx,*dyy,*dzz ;  /* principal strains -- temporary */

double* delv_xi, *delv_eta, *delv_zeta ;    /* velocity gradient -- temporary */

double* delx_xi, *delx_eta, *delx_zeta ;    /* coordinate gradient -- temporary */

double* e;
double* p;

double* q;
double* ql;
double* qq;

double* v;     /* relative volume */
double* vnew;  /* reference volume */
double* volo;  /* new relative volume -- temporary */
double* delv;  /* m_vnew - m_v */
double* m_vdov;  /* volume derivative over volume */

double* arealg; /* characteristic length of an element */

double* ss;

double* elemMass;

/* Node-centered */
double* x,*y,*z,*loc;                  /* coordinates */

double* xd,*yd,*zd;                 /* velocities */

double* xdd,*ydd,*zdd;              /* accelerations */

double* m_fx,*m_fy,*m_fz;                 /* forces */

double* nodalMass;          /* mass */

int* symmX,*symmY,*symmZ;              /* symmetry plane nodesets */

int* t_symmX,*t_symmY,*t_symmZ;
op_dat p_t_symmX,p_t_symmY,p_t_symmZ;

/* Other Constants */

// Variables to keep track of timestep, simulation time, and cycle
double  m_dtcourant ;         // courant constraint
double  m_dthydro ;           // volume change constraint
int   m_cycle ;             // iteration count for simulation
double  m_dtfixed ;           // fixed time increment
double  m_time ;              // current time
double  m_deltatime ;         // variable time increment
double  m_deltatimemultlb ;
double  m_deltatimemultub ;
double  m_dtmax ;             // maximum allowable time increment
double  m_stoptime ;          // end time for simulation

int sizeX ;
int sizeY ;
int sizeZ ;
int m_numElem ;
int m_numNode ;


//! global Data //
int g_sizeX ;
int g_sizeY ;
int g_sizeZ ;
int g_numElem ;
int g_numNode ;

// Region information
int *g_m_regElemSize ;   // Size of region sets
int *g_m_regNumList ;    // Region number per domain element
int **g_m_regElemlist ;  // region indexset

int* g_nodelist;

int* g_lxim; /* element connectivity across each face */
int* g_lxip;
int* g_letam;
int* g_letap;
int* g_lzetam;
int* g_lzetap;

int* g_elemBC;            /* symmetry/free-surface flags for each elem face */

double* g_e;
double* g_p;

double* g_q;

double* g_v;     /* relative volume */
double* g_volo;  /* new relative volume -- temporary */

double* g_ss;

double* g_elemMass;

/* Node-centered */
double* g_x,*g_y,*g_z,*g_loc;                  /* coordinates */

double* g_xd,*g_yd,*g_zd;                 /* velocities */

double* g_xdd,*g_ydd,*g_zdd;              /* accelerations */

double* g_nodalMass;          /* mass */

int* g_t_symmX,*g_t_symmY,*g_t_symmZ;

//! global Data end //






// Used in setup
int rowMin, rowMax;
int colMin, colMax;
int planeMin, planeMax ;

// OP2 Variables
op_set nodes;
op_set elems;
op_set symmetry;
op_set temp_vols;

op_map p_nodelist;
op_map p_symmX, p_symmY, p_symmZ;
op_map p_lxim;
op_map p_lxip;
op_map p_letam;
op_map p_letap;
op_map p_lzetam;
op_map p_lzetap;

op_dat p_fx, p_fy, p_fz;

op_dat p_xd, p_yd, p_zd;

op_dat p_xdd, p_ydd, p_zdd;

op_dat p_x, p_y, p_z;
op_dat p_loc;

op_dat p_qq;
op_dat p_ql;

op_dat p_e;

op_dat p_p;
op_dat p_q;
op_dat p_v;
op_dat p_volo;
op_dat p_arealg;
op_dat p_delv;
op_dat p_vnew;
op_dat p_vdov;

op_dat p_ss;
op_dat p_elemMass;
op_dat p_nodalMass;
op_dat p_elemBC;

op_dat p_delv_xi ;    /* velocity gradient -- temporary */
op_dat p_delv_eta ;
op_dat p_delv_zeta ;

op_dat p_delx_xi ;    /* coordinate gradient -- temporary */
op_dat p_delx_eta ;
op_dat p_delx_zeta ;

//Temporary OP2 Vars
double *sigxx, *sigyy, *sigzz;
op_dat p_sigxx, p_sigyy, p_sigzz;
double *determ;
op_dat p_determ;

double* dvdx,*dvdy, *dvdz;
double* x8n, *y8n, *z8n;
op_dat p_dvdx, p_dvdy, p_dvdz;
op_dat p_x8n, p_y8n, p_z8n;
op_dat p_dxx, p_dyy, p_dzz;

double* vnewc;
op_dat p_vnewc;


//EOS Temp variables
double *e_old;
double *delvc;
double *p_old;
double *q_old;
double *compression;
double *compHalfStep;
double *qq_old;
double *ql_old;
double *work;
double *p_new;
double *e_new;
double *q_new;
double *bvc;
double *pbvc;

double *pHalfStep;
op_dat p_pHalfStep;

op_dat p_e_old;
op_dat p_delvc;
op_dat p_p_old;
op_dat p_q_old;
op_dat p_compression;
op_dat p_compHalfStep;
op_dat p_qq_old;
op_dat p_ql_old;
op_dat p_work;
op_dat p_p_new;
op_dat p_e_new;
op_dat p_q_new;
op_dat p_bvc;
op_dat p_pbvc;

static inline
void allocateElems(){

   nodelist = (int*) malloc(m_numElem*8 * sizeof(int));

   lxim = (int*) malloc(m_numElem * sizeof(int));
   lxip = (int*) malloc(m_numElem * sizeof(int));
   letam = (int*) malloc(m_numElem * sizeof(int));
   letap = (int*) malloc(m_numElem * sizeof(int));
   lzetam = (int*) malloc(m_numElem * sizeof(int));
   lzetap = (int*) malloc(m_numElem * sizeof(int));

   elemBC = (int*) malloc(m_numElem * sizeof(int));

   e = (double*) malloc(m_numElem * sizeof(double));
   p = (double*) malloc(m_numElem * sizeof(double));

   q = (double*) malloc(m_numElem * sizeof(double));
   ql = (double*) malloc(m_numElem * sizeof(double));
   qq = (double*) malloc(m_numElem * sizeof(double));

   v = (double*) malloc(m_numElem * sizeof(double));

   volo = (double*) malloc(m_numElem * sizeof(double));
   delv = (double*) malloc(m_numElem * sizeof(double));
   m_vdov = (double*) malloc(m_numElem * sizeof(double));

   arealg = (double*) malloc(m_numElem * sizeof(double));

   ss = (double*) malloc(m_numElem * sizeof(double));

   elemMass = (double*) malloc(m_numElem * sizeof(double));

   vnew = (double*) malloc(m_numElem * sizeof(double));
   vnewc = (double*) malloc(m_numElem * sizeof(double));

   // sigxx = (double*) malloc(m_numElem* 3 * sizeof(double));
   sigxx = (double*) malloc(m_numElem* sizeof(double));
   sigyy = (double*) malloc(m_numElem* sizeof(double));
   sigzz = (double*) malloc(m_numElem* sizeof(double));
   determ = (double*) malloc(m_numElem * sizeof(double));

   dvdx = (double*) malloc(m_numElem * 8 * sizeof(double));
   dvdy = (double*) malloc(m_numElem * 8 * sizeof(double));
   dvdz = (double*) malloc(m_numElem * 8 * sizeof(double));
   x8n = (double*) malloc(m_numElem * 8 * sizeof(double));
   y8n = (double*) malloc(m_numElem * 8 * sizeof(double));
   z8n = (double*) malloc(m_numElem * 8 * sizeof(double));

   dxx = (double*) malloc(m_numElem * sizeof(double));
   dyy = (double*) malloc(m_numElem * sizeof(double));
   dzz = (double*) malloc(m_numElem * sizeof(double));

   delx_xi = (double*) malloc(m_numElem * sizeof(double));
   delx_eta = (double*) malloc(m_numElem * sizeof(double));
   delx_zeta = (double*) malloc(m_numElem * sizeof(double));

   delv_xi = (double*) malloc(m_numElem * sizeof(double));
   delv_eta = (double*) malloc(m_numElem * sizeof(double));
   delv_zeta = (double*) malloc(m_numElem * sizeof(double));

   //EOS Temp Vars
   e_old = (double*) malloc(m_numElem *sizeof(double));
   delvc = (double*) malloc(m_numElem *sizeof(double));
   p_old = (double*) malloc(m_numElem *sizeof(double));
   q_old = (double*) malloc(m_numElem *sizeof(double));
   compression = (double*) malloc(m_numElem *sizeof(double));
   compHalfStep = (double*) malloc(m_numElem *sizeof(double));
   qq_old = (double*) malloc(m_numElem *sizeof(double));
   ql_old = (double*) malloc(m_numElem *sizeof(double));
   work = (double*) malloc(m_numElem *sizeof(double));
   p_new = (double*) malloc(m_numElem *sizeof(double));
   e_new = (double*) malloc(m_numElem *sizeof(double));
   q_new = (double*) malloc(m_numElem *sizeof(double));
   bvc = (double*) malloc(m_numElem *sizeof(double));
   pbvc = (double*) malloc(m_numElem *sizeof(double));
   pHalfStep = (double*) malloc(m_numElem *sizeof(double));
};

static inline
void allocateGlobalElems(){
   g_nodelist = (int*) malloc(g_numElem*8 * sizeof(int));

   g_lxim = (int*) malloc(g_numElem * sizeof(int));
   g_lxip = (int*) malloc(g_numElem * sizeof(int));
   g_letam = (int*) malloc(g_numElem * sizeof(int));
   g_letap = (int*) malloc(g_numElem * sizeof(int));
   g_lzetam = (int*) malloc(g_numElem * sizeof(int));
   g_lzetap = (int*) malloc(g_numElem * sizeof(int));

   g_elemBC = (int*) malloc(g_numElem * sizeof(int));

   g_e = (double*) malloc(g_numElem * sizeof(double));
   g_p = (double*) malloc(g_numElem * sizeof(double));

   g_q = (double*) malloc(g_numElem * sizeof(double));

   g_v = (double*) malloc(g_numElem * sizeof(double));

   g_volo = (double*) malloc(g_numElem * sizeof(double));

   g_ss = (double*) malloc(g_numElem * sizeof(double));

   g_elemMass = (double*) malloc(g_numElem * sizeof(double));
};

static inline
void allocateNodes(){

   x = (double*) malloc(m_numNode * sizeof(double)); // Coordinates
   y = (double*) malloc(m_numNode * sizeof(double));
   z = (double*) malloc(m_numNode * sizeof(double));
   loc = (double*) malloc(3 * m_numNode * sizeof(double)); // Coordinates


   xd = (double*) malloc(m_numNode * sizeof(double)); // Velocities
   yd = (double*) malloc(m_numNode * sizeof(double));
   zd = (double*) malloc(m_numNode * sizeof(double));

   xdd = (double*) malloc(m_numNode * sizeof(double)); //Accelerations
   // xdd = (double*) malloc(m_numNode * 3 * sizeof(double)); //Accelerations
   ydd = (double*) malloc(m_numNode * sizeof(double));
   zdd = (double*) malloc(m_numNode * sizeof(double));

   m_fx = (double*) malloc(m_numNode * sizeof(double)); // Forces
   // m_fx = (double*) malloc(m_numNode * 3 * sizeof(double)); // Forces
   m_fy = (double*) malloc(m_numNode * sizeof(double));
   m_fz = (double*) malloc(m_numNode * sizeof(double));

   t_symmX = (int*) malloc(m_numNode * sizeof(int));
   t_symmY = (int*) malloc(m_numNode * sizeof(int));
   t_symmZ = (int*) malloc(m_numNode * sizeof(int));

   nodalMass = (double*) malloc(m_numNode * sizeof(double));  // mass
}

static inline
void allocateGlobalNodes(){

   g_x = (double*) malloc(g_numNode * sizeof(double)); // Coordinates
   g_y = (double*) malloc(g_numNode * sizeof(double));
   g_z = (double*) malloc(g_numNode * sizeof(double));
   g_loc = (double*) malloc(3 * g_numNode * sizeof(double));

   g_xd = (double*) malloc(g_numNode * sizeof(double)); // Velocities
   g_yd = (double*) malloc(g_numNode * sizeof(double));
   g_zd = (double*) malloc(g_numNode * sizeof(double));

   g_xdd = (double*) malloc(g_numNode * sizeof(double)); //Accelerations
   g_ydd = (double*) malloc(g_numNode * sizeof(double));
   g_zdd = (double*) malloc(g_numNode * sizeof(double));

   g_t_symmX = (int*) malloc(g_numNode * sizeof(int));
   g_t_symmY = (int*) malloc(g_numNode * sizeof(int));
   g_t_symmZ = (int*) malloc(g_numNode * sizeof(int));

   g_nodalMass = (double*) malloc(g_numNode * sizeof(double));  // mass
}


static inline
void initialiseGlobal(int nr){
   m_numReg = nr;

   m_dtfixed = double(-1.0e-6) ; // Negative means use courant condition
   m_stoptime  = double(1.0e-2); // *double(edgeElems*tp/45.0) ;

   // Initial conditions
   m_deltatimemultlb = double(1.1) ;
   m_deltatimemultub = double(1.2) ;
   m_dtcourant = double(1.0e+20) ;
   m_dthydro   = double(1.0e+20) ;
   m_dtmax     = double(1.0e-2) ;
   m_time    = double(0.) ;
   m_cycle   = int(0) ;
}


static int compute_local_size(int global_size, int mpi_comm_size,
                              int mpi_rank) {
  int local_size = global_size / mpi_comm_size;
  int remainder = (int)fmod(global_size, mpi_comm_size);

  if (mpi_rank < remainder) {
    local_size = local_size + 1;
  }
  return local_size;
}

//Given a cube size and an aray index caculate its x,y,z
static inline
void get_loc_from_index(int dx,int dy,int dz,int location,
      int *x, int *y, int *z){

   *x = location % dx ;
   *y = (location / dx) % dy ;
   *z = location / (dx*dy) ;
}

static inline
int get_rank_from_index(int index, int global_size, int number_of_ranks){
   int currentAccumulated =0;
   int local_size = global_size / number_of_ranks;
   int remainder = global_size % number_of_ranks;

   int rank=-1;
   while (currentAccumulated<=index){
      rank+=1;
      currentAccumulated+=local_size;
      if (rank < remainder) {
         currentAccumulated = currentAccumulated + 1;
      }
   }
   return rank;
}

static inline
Real_t volumeAt(int elemNumber,int meshEdgeElems,int edgeNodes,int edgeElems){
   int curr_z,curr_y,curr_x;

   get_loc_from_index(edgeElems,edgeElems,edgeElems,elemNumber, &curr_x, &curr_y, &curr_z);
   int nidx = (curr_z*edgeElems*edgeElems)
         + (curr_y*edgeElems)
         + (curr_x)

         + (curr_y)
         + (curr_z*edgeElems)

         + (curr_z*edgeNodes);
   int nodeNumbers[8];
   nodeNumbers[0] = nidx                                       ;
   nodeNumbers[1] = nidx                                   + 1 ;
   nodeNumbers[2] = nidx                       + edgeNodes + 1 ;
   nodeNumbers[3] = nidx                       + edgeNodes     ;
   nodeNumbers[4] = nidx + edgeNodes*edgeNodes                 ;
   nodeNumbers[5] = nidx + edgeNodes*edgeNodes             + 1 ;
   nodeNumbers[6] = nidx + edgeNodes*edgeNodes + edgeNodes + 1 ;
   nodeNumbers[7] = nidx + edgeNodes*edgeNodes + edgeNodes     ;

   double x_local[8], y_local[8], z_local[8] ;

   for( int lnode=0 ; lnode<8 ; ++lnode )
   {
      int gnode = nodeNumbers[lnode];
      get_loc_from_index(edgeNodes,edgeNodes,edgeNodes,gnode, &curr_x, &curr_y, &curr_z);

      x_local[lnode] = double(1.125)*double(curr_x)/double(meshEdgeElems);
      y_local[lnode] = double(1.125)*double(curr_y)/double(meshEdgeElems);
      z_local[lnode] = double(1.125)*double(curr_z)/double(meshEdgeElems);
   }
   // volume calculations
   return  CalcElemVolume(x_local, y_local, z_local );
}

static inline
int calc_Single_Elem_Offset_to_node(int loc_x,int loc_y,int loc_z,int edgeNodes){
   return (loc_x  ) + (loc_y  )*edgeNodes + (loc_z  )*edgeNodes*edgeNodes -loc_y     -(loc_z  ) -2*(loc_z  * edgeNodes);
}


//initialise using the static data
static inline
void initialise(int myRank,
               int nx, int tp, int nr, int balance, int cost, Int8_t numRanks){

   g_numElem = nx*nx*nx;
   g_numNode = (nx+1)*(nx+1)*(nx+1);
   g_sizeX = nx ;
   g_sizeY = nx ;
   g_sizeZ = nx ;

   m_numElem =  compute_local_size(g_numElem,(Int8_t)numRanks,myRank);
   m_numNode =  compute_local_size(g_numNode,(Int8_t)numRanks,myRank);

   int starting_m_numElem =0;
   int starting_m_numNode =0;
   for (int i =0; i<myRank; i++){
      starting_m_numElem +=  compute_local_size(g_numElem,(Int8_t)numRanks,i);
      starting_m_numNode +=  compute_local_size(g_numNode,(Int8_t)numRanks,i);
   }

   int colLoc =0;
   int rowLoc =0;
   int planeLoc=0;


   int edgeElems = nx ;
   int edgeNodes = edgeElems+1 ;

   //! Setup Comm Buffer, Should not be necessary in final app
   rowMin = (rowLoc == 0)        ? 0 : 1;
   rowMax = (rowLoc == tp-1)     ? 0 : 1;
   colMin = (colLoc == 0)        ? 0 : 1;
   colMax = (colLoc == tp-1)     ? 0 : 1;
   planeMin = (planeLoc == 0)    ? 0 : 1;
   planeMax = (planeLoc == tp-1) ? 0 : 1;

   //All orgional
   for(int i=0; i<m_numElem;++i){
      p[i] = double(0.0);
      e[i] = double(0.0);
      q[i] = double(0.0);
      ss[i] = double(0.0);
   }

   for(int i=0; i<m_numElem;++i){
      v[i] = double(1.0);
   }

   for (int i = 0; i<m_numNode;++i){
      xd[i] = double(0.0);
      yd[i] = double(0.0);
      zd[i] = double(0.0);
   }

   for (int i = 0; i<m_numNode;++i){
   // for (int i = 0; i<m_numNode*3;++i){
      xdd[i] = double(0.0);
      ydd[i] = double(0.0);
      zdd[i] = double(0.0);
   }

   for (int i=0; i<m_numNode; ++i) {
      nodalMass[i] = double(0.0) ;
   }

   //!Build Mesh function !HERE

   // Happy with this
   int meshEdgeElems = tp * nx;

   int curr_z,curr_y,curr_x;
   get_loc_from_index(edgeNodes,edgeNodes,edgeNodes,starting_m_numNode, &curr_x, &curr_y, &curr_z);


   int nidx = 0;

   while (nidx < m_numNode && curr_z<edgeNodes){
      double tz = double(1.125)*double(curr_z)/double(meshEdgeElems) ;
      while (nidx < m_numNode && curr_y<edgeNodes){
         double ty = double(1.125)*double(curr_y)/double(meshEdgeElems);
         while (nidx < m_numNode && curr_x<edgeNodes){
            double tx = double(1.125)*double(curr_x)/double(meshEdgeElems);
            x[nidx] = tx ;
            y[nidx] = ty ;
            z[nidx] = tz ;
            ++nidx;
            ++curr_x;
         }
         curr_x=0;
         ++curr_y;
      }
      curr_y=0;
      ++curr_z;
   }

   //Happy with this for now
   int zidx = 0 ;
   get_loc_from_index(edgeElems,edgeElems,edgeElems,starting_m_numElem, &curr_x, &curr_y, &curr_z);
   nidx = (curr_z*edgeElems*edgeElems)
         + (curr_y*edgeElems)
         + (curr_x)

         + (curr_y)
         + (curr_z*edgeElems)

         + (curr_z*edgeNodes);

   while (zidx < m_numElem && curr_z<edgeElems){
      while (zidx < m_numElem && curr_y<edgeElems){
         while (zidx < m_numElem && curr_x<edgeElems){
            int *localNode = &nodelist[zidx*int(8)] ;

            localNode[0] = nidx                                       ;
            localNode[1] = nidx                                   + 1 ;
            localNode[2] = nidx                       + edgeNodes + 1 ;
            localNode[3] = nidx                       + edgeNodes     ;
            localNode[4] = nidx + edgeNodes*edgeNodes                 ;
            localNode[5] = nidx + edgeNodes*edgeNodes             + 1 ;
            localNode[6] = nidx + edgeNodes*edgeNodes + edgeNodes + 1 ;
            localNode[7] = nidx + edgeNodes*edgeNodes + edgeNodes     ;

            ++zidx;
            ++nidx ;
            ++curr_x;
         }
         curr_x=0;
         ++curr_y;
         ++nidx ;
      }
      curr_y=0;
      ++curr_z;
      nidx += edgeNodes ;
   }

   //! End Build Mesh Function

   //! Start Create Region Sets
   //regions are just not used
   srand(0);
   // int myRank = 0;

   m_numReg = nr;
   m_regElemSize = (int*) malloc(m_numReg * sizeof(int));
   m_regElemlist = (int**) malloc(m_numReg * sizeof(int));

   //! End Create Region Sets

   //! Setup Symmetry Planes Function !HERE

   for (int i = 0; i<m_numNode;++i){
      t_symmX[i] |= FREE_NODE;
      t_symmY[i] |= FREE_NODE;
      t_symmZ[i] |= FREE_NODE;
   }

   //new Version
   nidx = 0 ;
   for (int i=0; i<edgeNodes; ++i) {
      int planeInc = i*edgeNodes*edgeNodes ;
      int rowInc   = i*edgeNodes ;
      for (int j=0; j<edgeNodes; ++j) {
         if (planeLoc == 0) {
            if ( rowInc + j>=starting_m_numNode && rowInc + j<starting_m_numNode+m_numNode ){
               t_symmZ[rowInc + j-starting_m_numNode] |= BOUNDARY_NODE;
            }
         }
         if (rowLoc == 0) {
            if ( planeInc + j>=starting_m_numNode && planeInc + j<starting_m_numNode+m_numNode ){
               t_symmY[planeInc+j-starting_m_numNode] |= BOUNDARY_NODE;
            }
         }
         if (colLoc == 0) {
            if ( planeInc + j*edgeNodes>=starting_m_numNode && planeInc + j*edgeNodes<starting_m_numNode+m_numNode ){
               t_symmX[planeInc + j*edgeNodes-starting_m_numNode] |= BOUNDARY_NODE;
            }
         }
         ++nidx ;
      }
   }
   //! End Setup Symmetry Plnes Function

   //! Setup Elem Connectivity Function
   for (int i=0; i<m_numElem; ++i) {
      lxim[i] = std::max(0,i+starting_m_numElem-1);
      lxip[i] = std::min(g_numElem-1  ,i+starting_m_numElem+1) ;
   }

   for (int i=0; i<m_numElem; ++i ){
      if (i+starting_m_numElem < edgeElems){
         letam[i] = i+starting_m_numElem;
      }else{
         letam[i] = i+starting_m_numElem-edgeElems ;
      }

      if (i+starting_m_numElem >= g_numElem-edgeElems){
         letap[i] = i+starting_m_numElem;
      }else{
         letap[i] = i+starting_m_numElem+edgeElems ;
      }
   }

   for (int i=0; i<m_numElem; ++i ){
      if (i+starting_m_numElem < edgeElems*edgeElems){
         lzetam[i] = i+starting_m_numElem;
      }else{
         lzetam[i] = i+starting_m_numElem-edgeElems*edgeElems ;
      }

      if (i+starting_m_numElem >= g_numElem-edgeElems*edgeElems){
         lzetap[i] = i+starting_m_numElem;
      }else{
         lzetap[i] = i+starting_m_numElem+edgeElems*edgeElems ;
      }
   }

   //! End Selem Connectivity Function


   //! Setup Boundary Conditions Function !HERE
   // set up boundary condition information
   int ghostIdx[6] ;  // offsets to ghost locations

   //New version
   for (int i=0; i<m_numElem; ++i) {
      elemBC[i] = int(0) ;
   }
   for (int i=0; i<6; ++i) {
      ghostIdx[i] = INT_MIN ;
   }

   int pidx = g_numElem ;
   if (planeMin != 0) {
      ghostIdx[0] = pidx ;
      pidx += g_sizeX*g_sizeY ;
   }
   if (planeMax != 0) {
      ghostIdx[1] = pidx ;
      pidx += g_sizeX*g_sizeY ;
   }
   if (rowMin != 0) {
      ghostIdx[2] = pidx ;
      pidx += g_sizeX*g_sizeZ ;
   }
   if (rowMax != 0) {
      ghostIdx[3] = pidx ;
      pidx += g_sizeX*g_sizeZ ;
   }
   if (colMin != 0) {
      ghostIdx[4] = pidx ;
      pidx += g_sizeY*g_sizeZ ;
   }
   if (colMax != 0) {
      ghostIdx[5] = pidx ;
   }

   for (int i=0; i<edgeElems; ++i) {

      int planeInc = i*edgeElems*edgeElems ;
      int rowInc   = i*edgeElems ;
      for (int j=0; j<edgeElems; ++j) {
         if (planeLoc == 0) {
            if ( (rowInc+j) >=starting_m_numElem && (rowInc+j)<starting_m_numElem+m_numElem ){
               elemBC[rowInc+j-starting_m_numElem] |= ZETA_M_SYMM ;
            }
         }
         else {
            if ( (rowInc+j) >=starting_m_numElem && (rowInc+j)<starting_m_numElem+m_numElem ){
               elemBC[rowInc+j-starting_m_numElem] |= ZETA_M_COMM ;
               lzetam[rowInc+j-starting_m_numElem] = ghostIdx[0] + rowInc + j ;
            }
         }
         if (planeLoc == tp-1) {
            if ( (rowInc+j+g_numElem-(edgeElems*edgeElems)) >=starting_m_numElem && (rowInc+j+g_numElem-(edgeElems*edgeElems))<starting_m_numElem+m_numElem ){
               elemBC[rowInc+j+g_numElem-(edgeElems*edgeElems)-starting_m_numElem] |= ZETA_P_FREE;
            }
         }
         else {
            if ( (rowInc+j+g_numElem-edgeElems*edgeElems) >=starting_m_numElem && (rowInc+j+g_numElem-edgeElems*edgeElems)<starting_m_numElem+m_numElem ){
               elemBC[rowInc+j+g_numElem-edgeElems*edgeElems-starting_m_numElem] |= ZETA_P_COMM ;
               lzetap[rowInc+j+g_numElem-edgeElems*edgeElems-starting_m_numElem] = ghostIdx[1] + rowInc + j ;
            }
         }
         if (rowLoc == 0) {
            if ( (planeInc+j) >=starting_m_numElem && (planeInc+j)<starting_m_numElem+m_numElem ){
               elemBC[planeInc+j-starting_m_numElem] |= ETA_M_SYMM ;
            }
         }
         else {
            if ( (planeInc+j) >=starting_m_numElem && (planeInc+j)<starting_m_numElem+m_numElem ){
               elemBC[planeInc+j-starting_m_numElem] |= ETA_M_COMM ;
               letam[planeInc+j-starting_m_numElem] = ghostIdx[2] + rowInc + j ;
            }
         }
         if (rowLoc == tp-1) {
            if ( (planeInc+j+edgeElems*edgeElems-edgeElems) >=starting_m_numElem && (planeInc+j+edgeElems*edgeElems-edgeElems)<starting_m_numElem+m_numElem ){
               elemBC[planeInc+j+edgeElems*edgeElems-edgeElems-starting_m_numElem] |= ETA_P_FREE ;
            }
         }
         else {
            if ( (planeInc+j+edgeElems*edgeElems-edgeElems) >=starting_m_numElem && (planeInc+j+edgeElems*edgeElems-edgeElems)<starting_m_numElem+m_numElem ){
               elemBC[planeInc+j+edgeElems*edgeElems-edgeElems-starting_m_numElem] |=  ETA_P_COMM ;
               letap[planeInc+j+edgeElems*edgeElems-edgeElems-starting_m_numElem] = ghostIdx[3] +  rowInc + j ;
            }
         }
         if (colLoc == 0) {
            if ( (planeInc+j*edgeElems) >=starting_m_numElem && (planeInc+j*edgeElems)<starting_m_numElem+m_numElem ){
               elemBC[planeInc+j*edgeElems-starting_m_numElem] |= XI_M_SYMM ;
            }
         }
         else {
            if ( (planeInc+j*edgeElems) >=starting_m_numElem && (planeInc+j*edgeElems)<starting_m_numElem+m_numElem ){
               elemBC[planeInc+j*edgeElems-starting_m_numElem] |= XI_M_COMM ;
               lxim[planeInc+j*edgeElems-starting_m_numElem] = ghostIdx[4] + rowInc + j ;
            }
         }
         if (colLoc == tp-1) {
            if ( (planeInc+j*edgeElems+edgeElems-1) >=starting_m_numElem && (planeInc+j*edgeElems+edgeElems-1)<starting_m_numElem+m_numElem ){
               elemBC[planeInc+j*edgeElems+edgeElems-1-starting_m_numElem] |= XI_P_FREE ;
            }
         }
         else {
            if ( (planeInc+j*edgeElems+edgeElems-1) >=starting_m_numElem && (planeInc+j*edgeElems+edgeElems-1)<starting_m_numElem+m_numElem ){
               elemBC[planeInc+j*edgeElems+edgeElems-1-starting_m_numElem] |= XI_P_COMM ;
               lxip[planeInc+j*edgeElems+edgeElems-1-starting_m_numElem] = ghostIdx[5] + rowInc + j ;
            }
         }
      }
   }


   //! End SetupBC Function

   // Setup defaults

   // These can be changed (requires recompile) if you want to run
   // with a fixed timestep, or to a different end time, but it's
   // probably easier/better to just run a fixed number of timesteps
   // using the -i flag in 2.x

   m_dtfixed = double(-1.0e-6) ; // Negative means use courant condition
   m_stoptime  = double(1.0e-2); // *double(edgeElems*tp/45.0) ;

   // Initial conditions
   m_deltatimemultlb = double(1.1) ;
   m_deltatimemultub = double(1.2) ;
   m_dtcourant = double(1.0e+20) ;
   m_dthydro   = double(1.0e+20) ;
   m_dtmax     = double(1.0e-2) ;
   m_time    = double(0.) ;
   m_cycle   = int(0) ;

   //! Difficult
   // initialize field data
   for (int i=0; i<m_numElem; ++i) {

      double x_local[8], y_local[8], z_local[8] ;
      int *m_elemToNode = &nodelist[i*int(8)] ;

      for( int lnode=0 ; lnode<8 ; ++lnode )
      {
        int gnode = m_elemToNode[lnode];
        get_loc_from_index(edgeNodes,edgeNodes,edgeNodes,gnode, &curr_x, &curr_y, &curr_z);

        x_local[lnode] = double(1.125)*double(curr_x)/double(meshEdgeElems);
        y_local[lnode] = double(1.125)*double(curr_y)/double(meshEdgeElems);
        z_local[lnode] = double(1.125)*double(curr_z)/double(meshEdgeElems);
      }

      // volume calculations
      double volume = CalcElemVolume(x_local, y_local, z_local);

      volo[i] = volume;
      elemMass[i] = volume ;
   }
   for (int i=0; i<m_numNode; ++i) {
      int *m_elemToNode = &nodelist[i*int(8)];
      get_loc_from_index(edgeNodes,edgeNodes,edgeNodes,i+starting_m_numNode, &curr_x, &curr_y, &curr_z);

      nodalMass[i] =0;
      if (curr_x<edgeElems && curr_y<edgeElems && curr_z<edgeElems)
      nodalMass[i] += volumeAt((curr_x  ) + (curr_y  )*edgeElems + (curr_z  )*edgeElems*edgeElems   ,meshEdgeElems,edgeNodes,edgeElems)/double(8.0) ;
      if (curr_x<edgeElems && curr_y<edgeElems && curr_z>0)
      nodalMass[i] += volumeAt((curr_x  ) + (curr_y  )*edgeElems + (curr_z-1)*edgeElems*edgeElems   ,meshEdgeElems,edgeNodes,edgeElems)/double(8.0) ;
      if (curr_x<edgeElems && curr_y>0         && curr_z<edgeElems)
      nodalMass[i] += volumeAt((curr_x  ) + (curr_y-1)*edgeElems + (curr_z  )*edgeElems*edgeElems   ,meshEdgeElems,edgeNodes,edgeElems)/double(8.0) ;
      if (curr_x<edgeElems && curr_y>0         && curr_z>0)
      nodalMass[i] += volumeAt((curr_x  ) + (curr_y-1)*edgeElems + (curr_z-1)*edgeElems*edgeElems   ,meshEdgeElems,edgeNodes,edgeElems)/double(8.0) ;
      if (curr_x>0         && curr_y<edgeElems && curr_z<edgeElems)
      nodalMass[i] += volumeAt((curr_x-1) + (curr_y  )*edgeElems + (curr_z  )*edgeElems*edgeElems   ,meshEdgeElems,edgeNodes,edgeElems)/double(8.0) ;
      if (curr_x>0         && curr_y<edgeElems && curr_z>0)
      nodalMass[i] += volumeAt((curr_x-1) + (curr_y  )*edgeElems + (curr_z-1)*edgeElems*edgeElems   ,meshEdgeElems,edgeNodes,edgeElems)/double(8.0) ;
      if (curr_x>0         && curr_y>0         && curr_z<edgeElems)
      nodalMass[i] += volumeAt((curr_x-1) + (curr_y-1)*edgeElems + (curr_z  )*edgeElems*edgeElems   ,meshEdgeElems,edgeNodes,edgeElems)/double(8.0) ;
      if (curr_x>0         && curr_y>0         && curr_z>0)
      nodalMass[i] += volumeAt((curr_x-1) + (curr_y-1)*edgeElems + (curr_z-1)*edgeElems*edgeElems   ,meshEdgeElems,edgeNodes,edgeElems)/double(8.0) ;

      // nodalMass[i] =0;
      // x   y    z -1
      // x   y    z -1
      // x   y -1 z -1
      // x   y -1 z -1
      // x-1 y    z -1
      // x-1 y    z -1
      // x-1 y -1 z -1
      // x-1 y -1 z -1
   }

   op_printf("Voume at 0: %e, Nodal Mass: %e\n", volo[0], nodalMass[0]);

   // deposit initial energy
   // An energy of 3.948746e+7 is correct for a problem with
   // 45 zones along a side - we need to scale it
   if (myRank==0){
      const double ebase = double(3.948746e+7);
      double scale = (nx * tp)/double(45.0);
      double einit = ebase*scale*scale*scale;
      if (rowLoc + colLoc + planeLoc == 0) {
         // Dump into the first zone (which we know is in the corner)
         // of the domain that sits at the origin
         e[0] = einit;
      }
      m_deltatime = (double(.5)*cbrt(volo[0]))/sqrt(double(2.0)*einit);
   }

   for (int i = 0; i < m_numNode; i++)
   {
      loc[3*i]=x[i];
      loc[3*i +1]=y[i];
      loc[3*i +2]=z[i];
   }
}

//initialise using the static data
static inline
void initialiseSingular(int colLoc,
               int rowLoc, int planeLoc,
               int nx, int tp, int nr, int balance, int cost, Int8_t numRanks){
   printf("%d %d %d %d %d %d\n",colLoc,rowLoc,planeLoc,nx,tp,nr);

   int edgeElems = nx ;
   int edgeNodes = edgeElems+1 ;

   // sizeX = edgeElems ;
   // sizeY = edgeElems ;
   // sizeZ = edgeElems ;
   // m_numElem = edgeElems * edgeElems * edgeElems;
   // m_numNode = edgeNodes * edgeNodes * edgeNodes;

   g_numElem = nx*nx*nx;
   g_numNode = (nx+1)*(nx+1)*(nx+1);
   g_sizeX = nx ;
   g_sizeY = nx ;
   g_sizeZ = nx ;

   printf("%d %d %d %d %d\n",g_sizeX,g_sizeY,g_sizeZ,g_numElem,g_numNode);


   g_m_regNumList = (int*) malloc(g_numElem * sizeof(int));

   printf("HERE NOW 3.1\n");
   //! Setup Comm Buffer, Should not be necessary in final app
   rowMin = (rowLoc == 0)        ? 0 : 1;
   rowMax = (rowLoc == tp-1)     ? 0 : 1;
   colMin = (colLoc == 0)        ? 0 : 1;
   colMax = (colLoc == tp-1)     ? 0 : 1;
   planeMin = (planeLoc == 0)    ? 0 : 1;
   planeMax = (planeLoc == tp-1) ? 0 : 1;

   for(int i=0; i<g_numElem;++i){
      g_p[i] = double(0.0);
      g_e[i] = double(0.0);
      g_q[i] = double(0.0);
      g_ss[i] = double(0.0);
   }

   for(int i=0; i<g_numElem;++i){
      g_v[i] = double(1.0);
   }

   for (int i = 0; i<g_numNode;++i){
      g_xd[i] = double(0.0);
      g_yd[i] = double(0.0);
      g_zd[i] = double(0.0);
   }

   for (int i = 0; i<g_numNode;++i){
   // for (int i = 0; i<m_numNode*3;++i){
      g_xdd[i] = double(0.0);
      g_ydd[i] = double(0.0);
      g_zdd[i] = double(0.0);
   }

   for (int i=0; i<g_numNode; ++i) {
      g_nodalMass[i] = double(0.0) ;
   }

   printf("HERE NOW 3.2\n");
   //!Build Mesh function !HERE
   int meshEdgeElems = tp * nx;
   // op_printf("Building Mesh with planeloc: %d, col: %d, rowLoc: %d, meshEdge: %d\n", planeLoc, colLoc, rowLoc, meshEdgeElems);
   int nidx = 0;
   double tz = double(1.125)*double(planeLoc*nx)/double(meshEdgeElems) ;
   for (int plane=0; plane<edgeNodes; ++plane) {
      double ty = double(1.125)*double(rowLoc*nx)/double(meshEdgeElems) ;
      for (int row=0; row<edgeNodes; ++row) {
         double tx = double(1.125)*double(colLoc*nx)/double(meshEdgeElems) ;
         for (int col=0; col<edgeNodes; ++col) {
         g_x[nidx] = tx ;
         g_y[nidx] = ty ;
         g_z[nidx] = tz ;
         ++nidx ;
         // tx += ds ; // may accumulate roundoff...
         tx = double(1.125)*double(colLoc*nx+col+1)/double(meshEdgeElems) ;
         }
         // ty += ds ;  // may accumulate roundoff...
         ty = double(1.125)*double(rowLoc*nx+row+1)/double(meshEdgeElems) ;
      }
      // tz += ds ;  // may accumulate roundoff...
      tz = double(1.125)*double(planeLoc*nx+plane+1)/double(meshEdgeElems) ;
   }




   // embed hexehedral elements in nodal point lattice
   int zidx = 0 ;
   nidx = 0 ;
   printf("HERE NOW 3.3\n");

   for (int plane=0; plane<edgeElems; ++plane) {
      for (int row=0; row<edgeElems; ++row) {
         for (int col=0; col<edgeElems; ++col) {
         int *localNode = &g_nodelist[zidx*int(8)] ;

         localNode[0] = nidx                                       ;
         localNode[1] = nidx                                   + 1 ;
         localNode[2] = nidx                       + edgeNodes + 1 ;
         localNode[3] = nidx                       + edgeNodes     ;
         localNode[4] = nidx + edgeNodes*edgeNodes                 ;
         localNode[5] = nidx + edgeNodes*edgeNodes             + 1 ;
         localNode[6] = nidx + edgeNodes*edgeNodes + edgeNodes + 1 ;
         localNode[7] = nidx + edgeNodes*edgeNodes + edgeNodes     ;
         ++zidx ;
         ++nidx ;
         }
         ++nidx ;
      }
      nidx += edgeNodes ;
   }

   //! End Build Mesh Function
   printf("HERE NOW 3.4\n");

   //! Start Create Region Sets
   srand(0);
   int myRank = 0;

   m_numReg = nr;
   g_m_regElemSize = (int*) malloc(m_numReg * sizeof(int));
   g_m_regElemlist = (int**) malloc(m_numReg * sizeof(int));
   int nextIndex = 0;
   //if we only have one region just fill it
   // Fill out the regNumList with material numbers, which are always
   // the region index plus one
   if(m_numReg == 1){
      printf("HERE NOW 3.4.1\n");

      while(nextIndex < g_numElem){
         g_m_regNumList[nextIndex] = 1;
         nextIndex++;
      }
      g_m_regElemSize[0] = 0;
   } else {//If we have more than one region distribute the elements.
      printf("HERE NOW 3.4.2\n");

      int regionNum;
      int regionVar;
      int lastReg = -1;
      int binSize;
      int elements;
      int runto = 0;
      int costDenominator = 0;
      int* g_regBinEnd = (int*) malloc(m_numReg * sizeof(int));
      //Determine the relative weights of all the regions.  This is based off the -b flag.  Balance is the value passed into b.
      for(int i=0; i<m_numReg;++i){
         g_m_regElemSize[i] = 0;
         costDenominator += pow((i+1), balance);//Total sum of all regions weights
         g_regBinEnd[i] = costDenominator;  //Chance of hitting a given region is (regBinEnd[i] - regBinEdn[i-1])/costDenominator
      }
      //Until all elements are assigned
      while (nextIndex < m_numElem) {
         //pick the region
         regionVar = rand() % costDenominator;
         int i = 0;
         while(regionVar >= g_regBinEnd[i]) i++;
         //rotate the regions based on MPI rank.  Rotation is Rank % m_numRegions this makes each domain have a different region with
         //the highest representation
         regionNum = ((i+myRank)% m_numReg) +1;
         while(regionNum == lastReg){
            regionVar = rand() % costDenominator;
            i = 0;
            while(regionVar >= g_regBinEnd[i]) i++;
            regionNum = ((i + myRank) % m_numReg) + 1;
         }
         //Pick the bin size of the region and determine the number of elements.
         binSize = rand() % 1000;
         if(binSize < 773) {
            elements = rand() % 15 + 1;
         }
         else if(binSize < 937) {
            elements = rand() % 16 + 16;
         }
         else if(binSize < 970) {
            elements = rand() % 32 + 32;
         }
         else if(binSize < 974) {
            elements = rand() % 64 + 64;
         }
         else if(binSize < 978) {
            elements = rand() % 128 + 128;
         }
         else if(binSize < 981) {
            elements = rand() % 256 + 256;
         }
         else
            elements = rand() % 1537 + 512;
         runto = elements + nextIndex;
         //Store the elements.  If we hit the end before we run out of elements then just stop.
         while (nextIndex < runto && nextIndex < m_numElem) {
            g_m_regNumList[nextIndex] = regionNum;
            nextIndex++;
         }
         lastReg = regionNum;
      }
      delete [] g_regBinEnd;
   }
   // Convert m_regNumList to region index sets
   // First, count size of each region
      printf("HERE NOW 3.4.3\n");

   for (int i=0 ; i<g_numElem ; ++i) {
      int r = g_m_regNumList[i]-1; // region index == regnum-1
      g_m_regElemSize[r]++;
   }
   // Second, allocate each region index set
   for (int i=0 ; i<m_numReg ; ++i) {
      g_m_regElemlist[i] = (int*) malloc(g_m_regElemSize[i]*sizeof(int));
      g_m_regElemSize[i] = 0;
   }
         printf("HERE NOW 3.4.4\n");

   // Third, fill index sets
   for (int i=0 ; i<g_numElem ; ++i) {

      int r = g_m_regNumList[i]-1;       // region index == regnum-1
      int regndx = g_m_regElemSize[r]++; // Note increment

      g_m_regElemlist[r][regndx] = i;
   }



   //! End Create Region Sets
   printf("HERE NOW 3.5\n");

   //! Setup Symmetry Planes Function !HERE
   for (int i = 0; i<g_numNode;++i){
   // for (int i = 0; i<m_numNode*3;++i){
      g_t_symmX[i] |= FREE_NODE;
      g_t_symmY[i] |= FREE_NODE;
      g_t_symmZ[i] |= FREE_NODE;
   }
   nidx = 0 ;
   for (int i=0; i<edgeNodes; ++i) {
      int planeInc = i*edgeNodes*edgeNodes ;
      int rowInc   = i*edgeNodes ;
      for (int j=0; j<edgeNodes; ++j) {
         if (planeLoc == 0) {
            // symmZ[nidx] = rowInc   + j ;
            g_t_symmZ[rowInc + j] |= BOUNDARY_NODE;
         }
         if (rowLoc == 0) {
            // symmY[nidx] = planeInc + j ;
            g_t_symmY[planeInc+j] |= BOUNDARY_NODE;
         }
         if (colLoc == 0) {
            // symmX[nidx] = planeInc + j*edgeNodes ;
            g_t_symmX[planeInc + j*edgeNodes] |= BOUNDARY_NODE;
         }
         ++nidx ;
      }
   }

   //! End Setup Symmetry Plnes Function
   printf("HERE NOW 3.6\n");

   //! Setup Elem Connectivity Function
   g_lxim[0] = 0 ;
   for (int i=1; i<g_numElem; ++i) {
      g_lxim[i]   = i-1 ;
      g_lxip[i-1] = i ;
   }
   g_lxip[g_numElem-1] = g_numElem-1 ;

   for (int i=0; i<edgeElems; ++i) {
      g_letam[i] = i ;
      g_letap[g_numElem-edgeElems+i] = g_numElem-edgeElems+i ;
   }
   for (int i=edgeElems; i<g_numElem; ++i) {
      g_letam[i] = i-edgeElems ;
      g_letap[i-edgeElems] = i ;
   }

   for (int i=0; i<edgeElems*edgeElems; ++i) {
      g_lzetam[i] = i ;
      g_lzetap[g_numElem-edgeElems*edgeElems+i] = g_numElem-edgeElems*edgeElems+i ;
   }
   for (int i=edgeElems*edgeElems; i<g_numElem; ++i) {
      g_lzetam[i] = i - edgeElems*edgeElems ;
      g_lzetap[i-edgeElems*edgeElems] = i ;
   }
   //! End Selem Connectivity Function
   printf("HERE NOW 3.7\n");

   //! Setup Boundary Conditions Function !HERE
   // set up boundary condition information
   int ghostIdx[6] ;  // offsets to ghost locations
   for (int i=0; i<g_numElem; ++i) {
      g_elemBC[i] = int(0) ;
   }
   for (int i=0; i<6; ++i) {
      ghostIdx[i] = INT_MIN ;
   }

  int pidx = g_numElem ;
  if (planeMin != 0) {
    ghostIdx[0] = pidx ;
    pidx += g_sizeX*g_sizeY ;
  }

  if (planeMax != 0) {
    ghostIdx[1] = pidx ;
    pidx += g_sizeX*g_sizeY ;
  }

  if (rowMin != 0) {
    ghostIdx[2] = pidx ;
    pidx += g_sizeX*g_sizeZ ;
  }

  if (rowMax != 0) {
    ghostIdx[3] = pidx ;
    pidx += g_sizeX*g_sizeZ ;
  }

  if (colMin != 0) {
    ghostIdx[4] = pidx ;
    pidx += g_sizeY*g_sizeZ ;
  }

  if (colMax != 0) {
    ghostIdx[5] = pidx ;
  }

  // symmetry plane or free surface BCs
  for (int i=0; i<edgeElems; ++i) {

   int planeInc = i*edgeElems*edgeElems ;
   int rowInc   = i*edgeElems ;
   for (int j=0; j<edgeElems; ++j) {
      if (planeLoc == 0) {
         g_elemBC[rowInc+j] |= ZETA_M_SYMM ;
      }
      else {
         g_elemBC[rowInc+j] |= ZETA_M_COMM ;
         g_lzetam[rowInc+j] = ghostIdx[0] + rowInc + j ;
      }
      if (planeLoc == tp-1) {
         g_elemBC[rowInc+j+m_numElem-edgeElems*edgeElems] |= ZETA_P_FREE;
      }
      else {
         g_elemBC[rowInc+j+m_numElem-edgeElems*edgeElems] |= ZETA_P_COMM ;
         g_lzetap[rowInc+j+m_numElem-edgeElems*edgeElems] = ghostIdx[1] + rowInc + j ;
      }
      if (rowLoc == 0) {
         g_elemBC[planeInc+j] |= ETA_M_SYMM ;
      }
      else {
         g_elemBC[planeInc+j] |= ETA_M_COMM ;
         g_letam[planeInc+j] = ghostIdx[2] + rowInc + j ;
      }
      if (rowLoc == tp-1) {
         g_elemBC[planeInc+j+edgeElems*edgeElems-edgeElems] |= ETA_P_FREE ;
      }
      else {
         g_elemBC[planeInc+j+edgeElems*edgeElems-edgeElems] |=  ETA_P_COMM ;
         g_letap[planeInc+j+edgeElems*edgeElems-edgeElems] = ghostIdx[3] +  rowInc + j ;
      }
      if (colLoc == 0) {
         g_elemBC[planeInc+j*edgeElems] |= XI_M_SYMM ;
      }
      else {
         g_elemBC[planeInc+j*edgeElems] |= XI_M_COMM ;
         g_lxim[planeInc+j*edgeElems] = ghostIdx[4] + rowInc + j ;
      }

      if (colLoc == tp-1) {
         g_elemBC[planeInc+j*edgeElems+edgeElems-1] |= XI_P_FREE ;
      }
      else {
         g_elemBC[planeInc+j*edgeElems+edgeElems-1] |= XI_P_COMM ;
         g_lxip[planeInc+j*edgeElems+edgeElems-1] = ghostIdx[5] + rowInc + j ;
      }
    }
   }
   //! End SetupBC Function
   printf("HERE NOW 3.8\n");

   // Setup defaults

   // These can be changed (requires recompile) if you want to run
   // with a fixed timestep, or to a different end time, but it's
   // probably easier/better to just run a fixed number of timesteps
   // using the -i flag in 2.x

   m_dtfixed = double(-1.0e-6) ; // Negative means use courant condition
   m_stoptime  = double(1.0e-2); // *double(edgeElems*tp/45.0) ;

   // Initial conditions
   m_deltatimemultlb = double(1.1) ;
   m_deltatimemultub = double(1.2) ;
   m_dtcourant = double(1.0e+20) ;
   m_dthydro   = double(1.0e+20) ;
   m_dtmax     = double(1.0e-2) ;
   m_time    = double(0.) ;
   m_cycle   = int(0) ;

   printf("HERE NOW 3.9\n");


   double mytotal =0;
   for (int i = 0; i < edgeNodes*edgeNodes*edgeNodes; i++)
   {
      mytotal+=g_x[i];
   }
   printf("curstep %e %d\n",mytotal,edgeNodes);



   // initialize field data
   for (int i=0; i<g_numElem; ++i) {

      double x_local[8], y_local[8], z_local[8] ;
      int *g_elemToNode = &g_nodelist[i*int(8)] ;
      for( int lnode=0 ; lnode<8 ; ++lnode )
      {
        int gnode = g_elemToNode[lnode];
        x_local[lnode] = g_x[gnode];
        y_local[lnode] = g_y[gnode];
        z_local[lnode] = g_z[gnode];
      }

      // volume calculations
      double volume = CalcElemVolume(x_local, y_local, z_local );

      g_volo[i] = volume ;
      g_elemMass[i] = volume ;
      for (int j=0; j<8; ++j) {
         int idx = g_elemToNode[j] ;
         g_nodalMass[idx] += volume / double(8.0) ;
      }
   }
   printf("HERE NOW 3.10\n");

   op_printf("Voume at 0: %e, Nodal Mass: %e\n", g_volo[0], g_nodalMass[0]);
   printf("HERE NOW 3.11\n");

   // deposit initial energy
   // An energy of 3.948746e+7 is correct for a problem with
   // 45 zones along a side - we need to scale it
   const double ebase = double(3.948746e+7);
   double scale = (nx * tp)/double(45.0);
   double einit = ebase*scale*scale*scale;
   if (rowLoc + colLoc + planeLoc == 0) {
      // Dump into the first zone (which we know is in the corner)
      // of the domain that sits at the origin
      g_e[0] = einit;
   }
   printf("HERE NOW 3.12\n");

   m_deltatime = (double(.5)*cbrt(g_volo[0]))/sqrt(double(2.0)*einit);

   for (int i = 0; i < edgeNodes*edgeNodes*edgeNodes; i++)
   {
      g_loc[3*i]=g_x[i];
      g_loc[3*i +1]=g_y[i];
      g_loc[3*i +2]=g_z[i];
   }

}

static void scatter_double_array(double *g_array, double *l_array,
                                 int comm_size, int g_size, int l_size,
                                 int elem_size) {
  int *sendcnts = (int *)malloc(comm_size * sizeof(int));
  int *displs = (int *)malloc(comm_size * sizeof(int));
  int disp = 0;

  for (int i = 0; i < comm_size; i++) {
    sendcnts[i] = elem_size * compute_local_size(g_size, comm_size, i);
  }
  for (int i = 0; i < comm_size; i++) {
    displs[i] = disp;
    disp = disp + sendcnts[i];
  }

  MPI_Scatterv(g_array, sendcnts, displs, MPI_DOUBLE, l_array,
               l_size * elem_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  free(sendcnts);
  free(displs);
}

static void scatter_int_array(int *g_array, int *l_array, int comm_size,
                              int g_size, int l_size, int elem_size) {
  int *sendcnts = (int *)malloc(comm_size * sizeof(int));
  int *displs = (int *)malloc(comm_size * sizeof(int));
  int disp = 0;

  for (int i = 0; i < comm_size; i++) {
    sendcnts[i] = elem_size * compute_local_size(g_size, comm_size, i);
  }
  for (int i = 0; i < comm_size; i++) {
    displs[i] = disp;
    disp = disp + sendcnts[i];
  }

  MPI_Scatterv(g_array, sendcnts, displs, MPI_INT, l_array, l_size * elem_size,
               MPI_INT, 0, MPI_COMM_WORLD);

  free(sendcnts);
  free(displs);
}

static inline
void distributeGlobalElems(int comm_size, int g_numElem , int numElem, int g_numNode, int numNode){
   //Element Information

   // g_nodelist = (int*) malloc(g_numElem*8 * sizeof(int));
   // g_lxim = (int*) malloc(g_numElem * sizeof(int));
   // g_lxip = (int*) malloc(g_numElem * sizeof(int));
   // g_letam = (int*) malloc(g_numElem * sizeof(int));
   // g_letap = (int*) malloc(g_numElem * sizeof(int));
   // g_lzetam = (int*) malloc(g_numElem * sizeof(int));
   // g_lzetap = (int*) malloc(g_numElem * sizeof(int));
   // g_elemBC = (int*) malloc(g_numElem * sizeof(int));
   scatter_int_array(g_nodelist,nodelist,comm_size,g_numElem,numElem,8);
   scatter_int_array(g_lxim,lxim,comm_size,g_numElem,numElem,1);
   scatter_int_array(g_lxip,lxip,comm_size,g_numElem,numElem,1);
   scatter_int_array(g_letam,letam,comm_size,g_numElem,numElem,1);
   scatter_int_array(g_letap,letap,comm_size,g_numElem,numElem,1);
   scatter_int_array(g_lzetam,lzetam,comm_size,g_numElem,numElem,1);
   scatter_int_array(g_lzetap,lzetap,comm_size,g_numElem,numElem,1);
   scatter_int_array(g_elemBC,elemBC,comm_size,g_numElem,numElem,1);
   free(g_nodelist);
   free(g_lxim);
   free(g_lxip);
   free(g_letam);
   free(g_letap);
   free(g_lzetam);
   free(g_lzetap);
   free(g_elemBC);

   // g_e = (double*) malloc(g_numElem * sizeof(double));
   // g_p = (double*) malloc(g_numElem * sizeof(double));
   // g_q = (double*) malloc(g_numElem * sizeof(double));
   // g_v = (double*) malloc(g_numElem * sizeof(double));
   // g_volo = (double*) malloc(g_numElem * sizeof(double));
   // g_ss = (double*) malloc(g_numElem * sizeof(double));
   // g_elemMass = (double*) malloc(g_numElem * sizeof(double));
   scatter_double_array(g_e,e,comm_size,g_numElem,numElem,1);
   scatter_double_array(g_p,p,comm_size,g_numElem,numElem,1);
   scatter_double_array(g_q,q,comm_size,g_numElem,numElem,1);
   scatter_double_array(g_v,v,comm_size,g_numElem,numElem,1);
   scatter_double_array(g_volo,volo,comm_size,g_numElem,numElem,1);
   scatter_double_array(g_ss,ss,comm_size,g_numElem,numElem,1);
   scatter_double_array(g_elemMass,elemMass,comm_size,g_numElem,numElem,1);
   free(g_e);
   free(g_p);
   free(g_q);
   free(g_v);
   free(g_volo);
   free(g_ss);
   free(g_elemMass);

   //Node Information
   // g_t_symmX = (int*) malloc(g_numNode * sizeof(int));
   // g_t_symmY = (int*) malloc(g_numNode * sizeof(int));
   // g_t_symmZ = (int*) malloc(g_numNode * sizeof(int));

   scatter_int_array(g_t_symmX,t_symmX,comm_size,g_numNode,numNode,1);
   scatter_int_array(g_t_symmY,t_symmY,comm_size,g_numNode,numNode,1);
   scatter_int_array(g_t_symmZ,t_symmZ,comm_size,g_numNode,numNode,1);
   free(g_t_symmX);
   free(g_t_symmY);
   free(g_t_symmZ);

   // g_x = (double*) malloc(g_numNode * sizeof(double)); // Coordinates
   // g_y = (double*) malloc(g_numNode * sizeof(double));
   // g_z = (double*) malloc(g_numNode * sizeof(double));

   // g_xd = (double*) malloc(g_numNode * sizeof(double)); // Velocities
   // g_yd = (double*) malloc(g_numNode * sizeof(double));
   // g_zd = (double*) malloc(g_numNode * sizeof(double));

   // g_xdd = (double*) malloc(g_numNode * sizeof(double)); //Accelerations
   // g_ydd = (double*) malloc(g_numNode * sizeof(double));
   // g_zdd = (double*) malloc(g_numNode * sizeof(double));

   // g_nodalMass = (double*) malloc(g_numNode * sizeof(double));  // mass
   scatter_double_array(g_x,x,comm_size,g_numNode,numNode,1);
   scatter_double_array(g_y,y,comm_size,g_numNode,numNode,1);
   scatter_double_array(g_z,z,comm_size,g_numNode,numNode,1);
   scatter_double_array(g_loc,loc,comm_size,g_numNode,numNode,3);
   free(g_x);
   free(g_y);
   free(g_z);
   free(g_loc);

   scatter_double_array(g_xd,xd,comm_size,g_numNode,numNode,1);
   scatter_double_array(g_yd,yd,comm_size,g_numNode,numNode,1);
   scatter_double_array(g_zd,zd,comm_size,g_numNode,numNode,1);
   scatter_double_array(g_xdd,xdd,comm_size,g_numNode,numNode,1);
   scatter_double_array(g_ydd,ydd,comm_size,g_numNode,numNode,1);
   scatter_double_array(g_zdd,zdd,comm_size,g_numNode,numNode,1);
   scatter_double_array(g_nodalMass,nodalMass,comm_size,g_numNode,numNode,1);
   free(g_xd);
   free(g_yd);
   free(g_zd);
   free(g_xdd);
   free(g_ydd);
   free(g_zdd);
   free(g_nodalMass);
}

static inline
void initOp2Vars(){
   nodes = op_decl_set(m_numNode, "nodes");
   elems = op_decl_set(m_numElem, "elems");
   // temp_vols = op_decl_set(m_numElem*8, "tempVols");

   p_nodelist = op_decl_map(elems, nodes, 8, nodelist, "nodelist");

   p_lxim = op_decl_map(elems, elems, 1, lxim, "lxim");
   p_lxip = op_decl_map(elems, elems, 1, lxip, "lxip");
   p_letam = op_decl_map(elems, elems, 1, letam, "letam");
   p_letap = op_decl_map(elems, elems, 1, letap, "letap");
   p_lzetam = op_decl_map(elems, elems, 1, lzetam, "lzetam");
   p_lzetap = op_decl_map(elems, elems, 1, lzetap, "lzetap");

   //Node Centred
   p_t_symmX = op_decl_dat(nodes, 1, "int", t_symmX, "t_symmX");
   p_t_symmY = op_decl_dat(nodes, 1, "int", t_symmY, "t_symmY");
   p_t_symmZ = op_decl_dat(nodes, 1, "int", t_symmZ, "t_symmZ");

   p_x = op_decl_dat(nodes, 1, "double", x, "p_x");
   p_y = op_decl_dat(nodes, 1, "double", y, "p_y");
   p_z = op_decl_dat(nodes, 1, "double", z, "p_z");
   p_loc = op_decl_dat(nodes, 3, "double", loc, "p_loc");
   p_xd = op_decl_dat(nodes, 1, "double", xd, "p_xd");
   p_yd = op_decl_dat(nodes, 1, "double", yd, "p_yd");
   p_zd = op_decl_dat(nodes, 1, "double", zd, "p_zd");
   p_xdd = op_decl_dat(nodes, 1, "double", xdd, "p_xdd");
   p_ydd = op_decl_dat(nodes, 1, "double", ydd, "p_ydd");
   p_zdd = op_decl_dat(nodes, 1, "double", zdd, "p_zdd");
   p_fx = op_decl_dat(nodes, 1, "double", m_fx, "p_fx");
   p_fy = op_decl_dat(nodes, 1, "double", m_fy, "p_fy");
   p_fz = op_decl_dat(nodes, 1, "double", m_fz, "p_fz");

   p_nodalMass = op_decl_dat(nodes, 1, "double", nodalMass, "p_nodalMass");
   //Elem Centred
   p_e = op_decl_dat(elems, 1, "double", e, "p_e");
   p_p = op_decl_dat(elems, 1, "double", p, "p_p");
   p_q = op_decl_dat(elems, 1, "double", q, "p_q");
   p_ql = op_decl_dat(elems, 1, "double", ql, "p_ql");
   p_qq = op_decl_dat(elems, 1, "double", qq, "p_qq");
   p_v = op_decl_dat(elems, 1, "double", v, "p_v");
   p_volo = op_decl_dat(elems, 1, "double", volo, "p_volo");
   p_delv = op_decl_dat(elems, 1, "double", delv, "p_delv");
   p_vdov = op_decl_dat(elems, 1, "double", m_vdov, "p_vdov");
   p_arealg = op_decl_dat(elems, 1, "double", arealg, "p_arealg");

   p_dxx = op_decl_dat(elems, 1, "double", dxx, "p_dxx");
   p_dyy = op_decl_dat(elems, 1, "double", dyy, "p_dyy");
   p_dzz = op_decl_dat(elems, 1, "double", dzz, "p_dzz");

   p_ss = op_decl_dat(elems, 1, "double", ss, "p_ss");
   p_elemMass = op_decl_dat(elems, 1, "double", elemMass, "p_elemMass");
   p_vnew = op_decl_dat(elems, 1, "double", vnew, "p_vnew");
   p_vnewc = op_decl_dat(elems, 1, "double", vnewc, "p_vnewc");

   //Temporary
   p_sigxx = op_decl_dat(elems, 1, "double", sigxx, "p_sigxx");
   p_sigyy = op_decl_dat(elems, 1, "double", sigyy, "p_sigyy");
   p_sigzz = op_decl_dat(elems, 1, "double", sigzz, "p_sigzz");
   p_determ = op_decl_dat(elems, 1, "double", determ, "p_determ");

   p_dvdx = op_decl_dat(elems, 8, "double",dvdx, "dvdx");
   p_dvdy = op_decl_dat(elems, 8, "double",dvdy, "dvdy");
   p_dvdz = op_decl_dat(elems, 8, "double",dvdz, "dvdz");
   p_x8n = op_decl_dat(elems, 8, "double",x8n, "x8n");
   p_y8n = op_decl_dat(elems, 8, "double",y8n, "y8n");
   p_z8n = op_decl_dat(elems, 8, "double",z8n, "z8n");

   p_delv_xi = op_decl_dat(elems, 1, "double", delv_xi, "p_delv_xi");
   p_delv_eta = op_decl_dat(elems, 1, "double", delv_eta, "p_delv_eta");
   p_delv_zeta = op_decl_dat(elems, 1, "double", delv_zeta, "p_delv_zeta");

   p_delx_xi = op_decl_dat(elems, 1, "double", delx_xi, "p_delx_xi");
   p_delx_eta = op_decl_dat(elems, 1, "double", delx_eta, "p_delx_eta");
   p_delx_zeta = op_decl_dat(elems, 1, "double", delx_zeta, "p_delx_zeta");

   p_elemBC = op_decl_dat(elems, 1, "int", elemBC, "p_elemBC");

   //EOS temp variables
   p_e_old = op_decl_dat(elems, 1, "double", e_old, "e_old");
   p_delvc = op_decl_dat(elems, 1, "double", delvc, "delvc");
   p_p_old = op_decl_dat(elems, 1, "double", p_old, "p_old");
   p_q_old = op_decl_dat(elems, 1, "double", q_old, "q_old");
   p_compression = op_decl_dat(elems, 1, "double", compression, "compression");
   p_compHalfStep = op_decl_dat(elems, 1, "double", compHalfStep, "compHalfStep");
   p_qq_old = op_decl_dat(elems, 1, "double", qq_old, "qq_old");
   p_ql_old = op_decl_dat(elems, 1, "double", ql_old, "ql_old");
   p_work = op_decl_dat(elems, 1, "double", work, "work");
   p_p_new = op_decl_dat(elems, 1, "double", p_new, "p_new");
   p_e_new = op_decl_dat(elems, 1, "double", e_new, "e_new");
   p_q_new = op_decl_dat(elems, 1, "double", q_new, "q_new");
   p_bvc = op_decl_dat(elems, 1, "double", bvc, "bvc"); ;
   p_pbvc = op_decl_dat(elems, 1, "double", pbvc, "pbvc");
   p_pHalfStep = op_decl_dat(elems, 1, "double", pHalfStep, "pHalfStep");
}


// This Function can be used when only reading from an HDF5 file
// It replaces the initialise function
static inline
void readOp2VarsFromHDF5(char* file){
   nodes = op_decl_set_hdf5(file, "nodes");
   elems = op_decl_set_hdf5(file, "elems");

   m_numElem = elems->size;
   m_numNode = nodes->size;

   p_nodelist = op_decl_map_hdf5(elems, nodes, 8, file, "nodelist");

   p_lxim = op_decl_map_hdf5(elems, elems, 1, file, "lxim");
   p_lxip = op_decl_map_hdf5(elems, elems, 1, file, "lxip");
   p_letam = op_decl_map_hdf5(elems, elems, 1, file, "letam");
   p_letap = op_decl_map_hdf5(elems, elems, 1, file, "letap");
   p_lzetam = op_decl_map_hdf5(elems, elems, 1, file, "lzetam");
   p_lzetap = op_decl_map_hdf5(elems, elems, 1, file, "lzetap");

   //Node Centred
   p_t_symmX = op_decl_dat_hdf5(nodes, 1, "int", file, "t_symmX");
   p_t_symmY = op_decl_dat_hdf5(nodes, 1, "int", file, "t_symmY");
   p_t_symmZ = op_decl_dat_hdf5(nodes, 1, "int", file, "t_symmZ");

   p_x = op_decl_dat_hdf5(nodes, 1, "double", file, "p_x");
   p_y = op_decl_dat_hdf5(nodes, 1, "double", file, "p_y");
   p_z = op_decl_dat_hdf5(nodes, 1, "double", file, "p_z");
   p_xd = op_decl_dat_hdf5(nodes, 1, "double", file, "p_xd");
   p_yd = op_decl_dat_hdf5(nodes, 1, "double", file, "p_yd");
   p_zd = op_decl_dat_hdf5(nodes, 1, "double", file, "p_zd");
   p_xdd = op_decl_dat_hdf5(nodes, 1, "double", file, "p_xdd");
   p_ydd = op_decl_dat_hdf5(nodes, 1, "double", file, "p_ydd");
   p_zdd = op_decl_dat_hdf5(nodes, 1, "double", file, "p_zdd");
   p_fx = op_decl_dat_hdf5(nodes, 1, "double", file, "p_fx");
   p_fy = op_decl_dat_hdf5(nodes, 1, "double", file, "p_fy");
   p_fz = op_decl_dat_hdf5(nodes, 1, "double", file, "p_fz");

   p_nodalMass = op_decl_dat_hdf5(nodes, 1, "double", file, "p_nodalMass");
   //Elem Centred
   p_e = op_decl_dat_hdf5(elems, 1, "double", file, "p_e");
   p_p = op_decl_dat_hdf5(elems, 1, "double", file, "p_p");
   p_q = op_decl_dat_hdf5(elems, 1, "double", file, "p_q");
   p_ql = op_decl_dat_hdf5(elems, 1, "double", file, "p_ql");
   p_qq = op_decl_dat_hdf5(elems, 1, "double", file, "p_qq");
   p_v = op_decl_dat_hdf5(elems, 1, "double", file, "p_v");
   p_volo = op_decl_dat_hdf5(elems, 1, "double", file, "p_volo");
   p_delv = op_decl_dat_hdf5(elems, 1, "double", file, "p_delv");
   p_vdov = op_decl_dat_hdf5(elems, 1, "double", file, "p_vdov");
   p_arealg = op_decl_dat_hdf5(elems, 1, "double", file, "p_arealg");

   p_dxx = op_decl_dat_hdf5(elems, 1, "double", file, "p_dxx");
   p_dyy = op_decl_dat_hdf5(elems, 1, "double", file, "p_dyy");
   p_dzz = op_decl_dat_hdf5(elems, 1, "double", file, "p_dzz");

   p_ss = op_decl_dat_hdf5(elems, 1, "double", file, "p_ss");
   p_elemMass = op_decl_dat_hdf5(elems, 1, "double", file, "p_elemMass");
   p_vnew = op_decl_dat_hdf5(elems, 1, "double", file, "p_vnew");
   p_vnewc = op_decl_dat_hdf5(elems, 1, "double", file, "p_vnewc");

   //Temporary
   // p_sigxx = op_decl_dat_hdf5(elems, 3, "double", file, "p_sigxx");
   p_sigxx = op_decl_dat_hdf5(elems, 1, "double", file, "p_sigxx");
   p_sigyy = op_decl_dat_hdf5(elems, 1, "double", file, "p_sigyy");
   p_sigzz = op_decl_dat_hdf5(elems, 1, "double", file, "p_sigzz");
   p_determ = op_decl_dat_hdf5(elems, 1, "double", file, "p_determ");

   p_dvdx = op_decl_dat_hdf5(elems, 8, "double", file, "dvdx");
   p_dvdy = op_decl_dat_hdf5(elems, 8, "double", file, "dvdy");
   p_dvdz = op_decl_dat_hdf5(elems, 8, "double", file, "dvdz");
   p_x8n = op_decl_dat_hdf5(elems, 8, "double", file, "x8n");
   p_y8n = op_decl_dat_hdf5(elems, 8, "double", file, "y8n");
   p_z8n = op_decl_dat_hdf5(elems, 8, "double", file, "z8n");

   p_delv_xi = op_decl_dat_hdf5(elems, 1, "double", file, "p_delv_xi");
   p_delv_eta = op_decl_dat_hdf5(elems, 1, "double", file, "p_delv_eta");
   p_delv_zeta = op_decl_dat_hdf5(elems, 1, "double", file, "p_delv_zeta");

   p_delx_xi = op_decl_dat_hdf5(elems, 1, "double", file, "p_delx_xi");
   p_delx_eta = op_decl_dat_hdf5(elems, 1, "double", file, "p_delx_eta");
   p_delx_zeta = op_decl_dat_hdf5(elems, 1, "double", file, "p_delx_zeta");

   p_elemBC = op_decl_dat_hdf5(elems, 1, "int", file, "p_elemBC");

   //EOS temp variables
   p_e_old = op_decl_dat_hdf5(elems, 1, "double", file, "e_old");
   p_delvc = op_decl_dat_hdf5(elems, 1, "double", file, "delvc");
   p_p_old = op_decl_dat_hdf5(elems, 1, "double", file, "p_old");
   p_q_old = op_decl_dat_hdf5(elems, 1, "double", file, "q_old");
   p_compression = op_decl_dat_hdf5(elems, 1, "double", file, "compression");
   p_compHalfStep = op_decl_dat_hdf5(elems, 1, "double", file, "compHalfStep");
   p_qq_old = op_decl_dat_hdf5(elems, 1, "double", file, "qq_old");
   p_ql_old = op_decl_dat_hdf5(elems, 1, "double", file, "ql_old");
   p_work = op_decl_dat_hdf5(elems, 1, "double", file, "work");
   p_p_new = op_decl_dat_hdf5(elems, 1, "double", file, "p_new");
   p_e_new = op_decl_dat_hdf5(elems, 1, "double", file, "e_new");
   p_q_new = op_decl_dat_hdf5(elems, 1, "double", file, "q_new");
   p_bvc = op_decl_dat_hdf5(elems, 1, "double", file, "bvc"); ;
   p_pbvc = op_decl_dat_hdf5(elems, 1, "double", file, "pbvc");
   p_pHalfStep = op_decl_dat_hdf5(elems, 1, "double", file, "pHalfStep");

   op_get_const_hdf5("deltatime", 1, "double", (char *)&m_deltatime,  file);
}

// This function was used as a hybrid for reading initialised variables and initalising the ones that are set to 0 itself
static inline
void readAndInitVars(char* file){
   nodes = op_decl_set_hdf5(file, "nodes");
   elems = op_decl_set_hdf5(file, "elems");

   m_numElem = elems->size;
   m_numNode = nodes->size;

   // Allocate Elems
   // e = (double*) malloc(m_numElem * sizeof(double));
   p = (double*) malloc(m_numElem * sizeof(double));
   q = (double*) malloc(m_numElem * sizeof(double));
   ss = (double*) malloc(m_numElem * sizeof(double));
   v = (double*) malloc(m_numElem * sizeof(double));
   delv = (double*) malloc(m_numElem * sizeof(double));
   m_vdov = (double*) malloc(m_numElem * sizeof(double));
   arealg = (double*) malloc(m_numElem * sizeof(double));
   vnew = (double*) malloc(m_numElem * sizeof(double));
   vnewc = (double*) malloc(m_numElem * sizeof(double));
   // Allocate Nodes
   xd = (double*) malloc(m_numNode * sizeof(double)); // Velocities
   yd = (double*) malloc(m_numNode * sizeof(double));
   zd = (double*) malloc(m_numNode * sizeof(double));
   xdd = (double*) malloc(m_numNode * sizeof(double)); //Accelerations
   ydd = (double*) malloc(m_numNode * sizeof(double));
   zdd = (double*) malloc(m_numNode * sizeof(double));
   m_fx = (double*) malloc(m_numNode * sizeof(double)); // Forces
   m_fy = (double*) malloc(m_numNode * sizeof(double));
   m_fz = (double*) malloc(m_numNode * sizeof(double));
   ql = (double*) malloc(m_numElem * sizeof(double));
   qq = (double*) malloc(m_numElem * sizeof(double));
   //Temporary Vars
   sigxx = (double*) malloc(m_numElem*  sizeof(double));
   sigyy = (double*) malloc(m_numElem*  sizeof(double));
   sigzz = (double*) malloc(m_numElem*  sizeof(double));
   determ = (double*) malloc(m_numElem * sizeof(double));

   dvdx = (double*) malloc(m_numElem * 8 * sizeof(double));
   dvdy = (double*) malloc(m_numElem * 8 * sizeof(double));
   dvdz = (double*) malloc(m_numElem * 8 * sizeof(double));
   x8n = (double*) malloc(m_numElem * 8 * sizeof(double));
   y8n = (double*) malloc(m_numElem * 8 * sizeof(double));
   z8n = (double*) malloc(m_numElem * 8 * sizeof(double));

   dxx = (double*) malloc(m_numElem * sizeof(double));
   dyy = (double*) malloc(m_numElem * sizeof(double));
   dzz = (double*) malloc(m_numElem * sizeof(double));

   delx_xi = (double*) malloc(m_numElem * sizeof(double));
   delx_eta = (double*) malloc(m_numElem * sizeof(double));
   delx_zeta = (double*) malloc(m_numElem * sizeof(double));

   delv_xi = (double*) malloc(m_numElem * sizeof(double));
   delv_eta = (double*) malloc(m_numElem * sizeof(double));
   delv_zeta = (double*) malloc(m_numElem * sizeof(double));

   //EOS Temp Vars
   e_old = (double*) malloc(m_numElem *sizeof(double));
   delvc = (double*) malloc(m_numElem *sizeof(double));
   p_old = (double*) malloc(m_numElem *sizeof(double));
   q_old = (double*) malloc(m_numElem *sizeof(double));
   compression = (double*) malloc(m_numElem *sizeof(double));
   compHalfStep = (double*) malloc(m_numElem *sizeof(double));
   qq_old = (double*) malloc(m_numElem *sizeof(double));
   ql_old = (double*) malloc(m_numElem *sizeof(double));
   work = (double*) malloc(m_numElem *sizeof(double));
   p_new = (double*) malloc(m_numElem *sizeof(double));
   e_new = (double*) malloc(m_numElem *sizeof(double));
   q_new = (double*) malloc(m_numElem *sizeof(double));
   bvc = (double*) malloc(m_numElem *sizeof(double));
   pbvc = (double*) malloc(m_numElem *sizeof(double));
   pHalfStep = (double*) malloc(m_numElem *sizeof(double));


   for(int i=0; i<m_numElem;++i){
      // p[i] = double(0.0);
      e[i] = double(0.0);
      q[i] = double(0.0);
      ss[i] = double(0.0);
   }
   // p_e = op_decl_dat(elems, 1, "double", e, "p_e");
   p_p = op_decl_dat(elems, 1, "double", p, "p_p");
   p_q = op_decl_dat(elems, 1, "double", q, "p_q");
   p_ss = op_decl_dat(elems, 1, "double", ss, "p_ss");

   for(int i=0; i<m_numElem;++i){
      v[i] = double(1.0);
   }
   p_v = op_decl_dat(elems, 1, "double", v, "p_v");

   for (int i = 0; i<m_numNode;++i){
      xd[i] = double(0.0);
      yd[i] = double(0.0);
      zd[i] = double(0.0);
   }
   p_xd = op_decl_dat(nodes, 1, "double", xd, "p_xd");
   p_yd = op_decl_dat(nodes, 1, "double", yd, "p_yd");
   p_zd = op_decl_dat(nodes, 1, "double", zd, "p_zd");

   for (int i = 0; i<m_numNode;++i){
   // for (int i = 0; i<m_numNode*3;++i){
      xdd[i] = double(0.0);
      ydd[i] = double(0.0);
      zdd[i] = double(0.0);
   }
   p_xdd = op_decl_dat(nodes, 1, "double", xdd, "p_xdd");
   p_ydd = op_decl_dat(nodes, 1, "double", ydd, "p_ydd");
   p_zdd = op_decl_dat(nodes, 1, "double", zdd, "p_zdd");

   p_ql = op_decl_dat(elems, 1, "double", ql, "p_ql");
   p_qq = op_decl_dat(elems, 1, "double", qq, "p_qq");
   p_delv = op_decl_dat(elems, 1, "double", delv, "p_delv");
   p_vdov = op_decl_dat(elems, 1, "double", m_vdov, "p_vdov");
   p_arealg = op_decl_dat(elems, 1, "double", arealg, "p_arealg");
   p_vnew = op_decl_dat(elems, 1, "double", vnew, "p_vnew");
   p_vnewc = op_decl_dat(elems, 1, "double", vnewc, "p_vnewc");
   p_fx = op_decl_dat(nodes, 1, "double", m_fx, "p_fx");
   p_fy = op_decl_dat(nodes, 1, "double", m_fy, "p_fy");
   p_fz = op_decl_dat(nodes, 1, "double", m_fz, "p_fz");

   //Temporary Vars
   p_sigxx = op_decl_dat(elems, 1, "double", sigxx, "p_sigxx");
   p_sigyy = op_decl_dat(elems, 1, "double", sigyy, "p_sigyy");
   p_sigzz = op_decl_dat(elems, 1, "double", sigzz, "p_sigzz");
   p_determ = op_decl_dat(elems, 1, "double", determ, "p_determ");

   p_dxx = op_decl_dat(elems, 1, "double", dxx, "p_dxx");
   p_dyy = op_decl_dat(elems, 1, "double", dyy, "p_dyy");
   p_dzz = op_decl_dat(elems, 1, "double", dzz, "p_dzz");

   p_dvdx = op_decl_dat(elems, 8, "double",dvdx, "dvdx");
   p_dvdy = op_decl_dat(elems, 8, "double",dvdy, "dvdy");
   p_dvdz = op_decl_dat(elems, 8, "double",dvdz, "dvdz");
   p_x8n = op_decl_dat(elems, 8, "double",x8n, "x8n");
   p_y8n = op_decl_dat(elems, 8, "double",y8n, "y8n");
   p_z8n = op_decl_dat(elems, 8, "double",z8n, "z8n");

   p_delv_xi = op_decl_dat(elems, 1, "double", delv_xi, "p_delv_xi");
   p_delv_eta = op_decl_dat(elems, 1, "double", delv_eta, "p_delv_eta");
   p_delv_zeta = op_decl_dat(elems, 1, "double", delv_zeta, "p_delv_zeta");

   p_delx_xi = op_decl_dat(elems, 1, "double", delx_xi, "p_delx_xi");
   p_delx_eta = op_decl_dat(elems, 1, "double", delx_eta, "p_delx_eta");
   p_delx_zeta = op_decl_dat(elems, 1, "double", delx_zeta, "p_delx_zeta");
   //EOS temp variables
   p_e_old = op_decl_dat(elems, 1, "double", e_old, "e_old");
   p_delvc = op_decl_dat(elems, 1, "double", delvc, "delvc");
   p_p_old = op_decl_dat(elems, 1, "double", p_old, "p_old");
   p_q_old = op_decl_dat(elems, 1, "double", q_old, "q_old");
   p_compression = op_decl_dat(elems, 1, "double", compression, "compression");
   p_compHalfStep = op_decl_dat(elems, 1, "double", compHalfStep, "compHalfStep");
   p_qq_old = op_decl_dat(elems, 1, "double", qq_old, "qq_old");
   p_ql_old = op_decl_dat(elems, 1, "double", ql_old, "ql_old");
   p_work = op_decl_dat(elems, 1, "double", work, "work");
   p_p_new = op_decl_dat(elems, 1, "double", p_new, "p_new");
   p_e_new = op_decl_dat(elems, 1, "double", e_new, "e_new");
   p_q_new = op_decl_dat(elems, 1, "double", q_new, "q_new");
   p_bvc = op_decl_dat(elems, 1, "double", bvc, "bvc"); ;
   p_pbvc = op_decl_dat(elems, 1, "double", pbvc, "pbvc");
   p_pHalfStep = op_decl_dat(elems, 1, "double", pHalfStep, "pHalfStep");

   //HDF5 Read
   p_nodelist = op_decl_map_hdf5(elems, nodes, 8, file, "nodelist");
   p_lxim = op_decl_map_hdf5(elems, elems, 1, file, "lxim");
   p_lxip = op_decl_map_hdf5(elems, elems, 1, file, "lxip");
   p_letam = op_decl_map_hdf5(elems, elems, 1, file, "letam");
   p_letap = op_decl_map_hdf5(elems, elems, 1, file, "letap");
   p_lzetam = op_decl_map_hdf5(elems, elems, 1, file, "lzetam");
   p_lzetap = op_decl_map_hdf5(elems, elems, 1, file, "lzetap");

   p_x = op_decl_dat_hdf5(nodes, 1, "double", file, "p_x");
   p_y = op_decl_dat_hdf5(nodes, 1, "double", file, "p_y");
   p_z = op_decl_dat_hdf5(nodes, 1, "double", file, "p_z");
   p_nodalMass = op_decl_dat_hdf5(nodes, 1, "double", file, "p_nodalMass");
   p_volo = op_decl_dat_hdf5(elems, 1, "double", file, "p_volo");
   p_e = op_decl_dat_hdf5(elems, 1, "double", file, "p_e");
   p_elemMass = op_decl_dat_hdf5(elems, 1, "double", file, "p_elemMass");
   p_elemBC = op_decl_dat_hdf5(elems, 1, "int", file, "p_elemBC");
   p_t_symmX = op_decl_dat_hdf5(nodes, 1, "int", file, "t_symmX");
   p_t_symmY = op_decl_dat_hdf5(nodes, 1, "int", file, "t_symmY");
   p_t_symmZ = op_decl_dat_hdf5(nodes, 1, "int", file, "t_symmZ");
   op_get_const_hdf5("deltatime", 1, "double", (char *)&m_deltatime,  file);

      //! Start Create Region Sets
   srand(0);
   int myRank = 0;

   m_numReg = 1;
   m_regElemSize = (int*) malloc(m_numReg * sizeof(int));
   m_regElemlist = (int**) malloc(m_numReg * sizeof(int));
   int nextIndex = 0;
   //if we only have one region just fill it
   // Fill out the regNumList with material numbers, which are always
   // the region index plus one
   if(m_numReg == 1){
      m_regNumList = (int*) malloc(m_numElem * sizeof(int));
      while(nextIndex < m_numElem){
         m_regNumList[nextIndex] = 1;
         nextIndex++;
      }
      m_regElemSize[0] = 0;
   } else {//If we have more than one region distribute the elements.
      int regionNum;
      int regionVar;
      int lastReg = -1;
      int binSize;
      int elements;
      int runto = 0;
      int costDenominator = 0;
      int* regBinEnd = (int*) malloc(m_numReg * sizeof(int));
      //Determine the relative weights of all the regions.  This is based off the -b flag.  Balance is the value passed into b.
      for(int i=0; i<m_numReg;++i){
         m_regElemSize[i] = 0;
         costDenominator += pow((i+1), 1);//Total sum of all regions weights
         regBinEnd[i] = costDenominator;  //Chance of hitting a given region is (regBinEnd[i] - regBinEdn[i-1])/costDenominator
      }
      //Until all elements are assigned
      while (nextIndex < m_numElem) {
         //pick the region
         regionVar = rand() % costDenominator;
         int i = 0;
         while(regionVar >= regBinEnd[i]) i++;
         //rotate the regions based on MPI rank.  Rotation is Rank % m_numRegions this makes each domain have a different region with
         //the highest representation
         regionNum = ((i+myRank)% m_numReg) +1;
         while(regionNum == lastReg){
            regionVar = rand() % costDenominator;
            i = 0;
            while(regionVar >= regBinEnd[i]) i++;
            regionNum = ((i + myRank) % m_numReg) + 1;
         }
         //Pick the bin size of the region and determine the number of elements.
         binSize = rand() % 1000;
         if(binSize < 773) {
            elements = rand() % 15 + 1;
         }
         else if(binSize < 937) {
            elements = rand() % 16 + 16;
         }
         else if(binSize < 970) {
            elements = rand() % 32 + 32;
         }
         else if(binSize < 974) {
            elements = rand() % 64 + 64;
         }
         else if(binSize < 978) {
            elements = rand() % 128 + 128;
         }
         else if(binSize < 981) {
            elements = rand() % 256 + 256;
         }
         else
            elements = rand() % 1537 + 512;
         runto = elements + nextIndex;
         //Store the elements.  If we hit the end before we run out of elements then just stop.
         while (nextIndex < runto && nextIndex < m_numElem) {
            m_regNumList[nextIndex] = regionNum;
            nextIndex++;
         }
         lastReg = regionNum;
      }
      delete [] regBinEnd;
   }
      // Convert m_regNumList to region index sets
   // First, count size of each region
   for (int i=0 ; i<m_numElem ; ++i) {
      int r = m_regNumList[i]-1; // region index == regnum-1
      m_regElemSize[r]++;
   }
   // Second, allocate each region index set
   for (int i=0 ; i<m_numReg ; ++i) {
      m_regElemlist[i] = (int*) malloc(m_regElemSize[i]*sizeof(int));
      m_regElemSize[i] = 0;
   }
   // Third, fill index sets
   for (int i=0 ; i<m_numElem ; ++i) {
      int r = m_regNumList[i]-1;       // region index == regnum-1
      int regndx = m_regElemSize[r]++; // Note increment
      m_regElemlist[r][regndx] = i;
   }
   //! End Create Region Sets
   // free(e);
   free(p);
   free(q);
   free(ss);
   free(v);
   free(delv);
   free(m_vdov);
   free(arealg);
   free(vnew);
   free(vnewc);
   free(xd);
   free(yd);
   free(zd);
   free(xdd);
   free(ydd);
   free(zdd);
   free(m_fx);
   free(m_fy);
   free(m_fz);
   free(ql);
   free(qq);
   free(sigxx);
   free(sigyy);
   free(sigzz);
   free(determ);
   free(dvdx);
   free(dvdy);
   free(dvdz);
   free(x8n);
   free(y8n);
   free(z8n);
   free(dxx);
   free(dyy);
   free(dzz);
   free(delx_xi);
   free(delx_eta);
   free(delx_zeta);
   free(delv_xi);
   free(delv_eta);
   free(delv_zeta);
   free(e_old);
   free(delvc);
   free(p_old);
   free(q_old);
   free(compression);
   free(compHalfStep);
   free(qq_old);
   free(ql_old);
   free(work);
   free(p_new);
   free(e_new);
   free(q_new);
   free(bvc);
   free(pbvc);
   free(pHalfStep);


}

/* Work Routines */

static inline
void TimeIncrement(int myrank)
{
   double targetdt = m_stoptime - m_time ;
   if ((m_dtfixed <= double(0.0)) && (m_cycle != int(0))) {
      double ratio ;
      double olddt = m_deltatime ;

      /* This will require a reduction in parallel */
      /* Reduction is done in CalcHydroConstraint and CalcCourantConstraint for OP2*/
      double gnewdt = double(1.0e+20) ;
      double newdt ;
      if (m_dtcourant < gnewdt) {
         gnewdt = m_dtcourant / double(2.0) ;
      }
      if (m_dthydro < gnewdt) {
         gnewdt = m_dthydro * double(2.0) / double(3.0) ;
      }

      newdt = gnewdt;

      ratio = newdt / olddt ;
      if (ratio >= double(1.0)) {
         if (ratio < m_deltatimemultlb) {
            newdt = olddt ;
         }
         else if (ratio > m_deltatimemultub) {
            newdt = olddt*m_deltatimemultub ;
         }
      }
      if (newdt > m_dtmax) {
         newdt = m_dtmax ;
      }
      m_deltatime = newdt ;

   }

   /* TRY TO PREVENT VERY SMALL SCALING ON THE NEXT CYCLE */
   if ((targetdt > m_deltatime) &&
       (targetdt < (double(4.0) * m_deltatime / double(3.0))) ) {
      targetdt = double(2.0) * m_deltatime / double(3.0) ;
   }

   if (targetdt < m_deltatime) {
      m_deltatime = targetdt ;
   }

   m_time += m_deltatime ;

   ++m_cycle ;
}



static inline
void InitStressTermsForElems()
{
   // pull in the stresses appropriate to the hydro integration
   // std::cout<<"init stress at " << myRank << "\n";
   op_par_loop(initStressTerms, "initStressTerms", elems,
               op_arg_dat(p_sigxx, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_sigyy, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_sigzz, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_p, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_q, -1, OP_ID, 1, "double", OP_READ));

   // MPI_Barrier(MPI_COMM_WORLD);
}

/******************************************/

static inline
void IntegrateStressForElems()
{
//   if (numthreads > 1) {
//      fx_elem = Allocate<double>(numElem8) ;
//      fy_elem = Allocate<double>(numElem8) ;
//      fz_elem = Allocate<double>(numElem8) ;
//   }
  // loop over all elements
   op_par_loop(IntegrateStressForElemsLoop, "IntegrateStressForElemsLoop", elems,
               op_arg_dat(p_x, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_x, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_x, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_y, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_y, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_y, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_z, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_z, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_z, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_fx, 0, p_nodelist, 1, "double", OP_INC), op_arg_dat(p_fx, 1, p_nodelist, 1, "double", OP_INC), op_arg_dat(p_fx, 2, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fx, 3, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fx, 4, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fx, 5, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fx, 6, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fx, 7, p_nodelist, 1, "double", OP_INC),
               op_arg_dat(p_fy, 0, p_nodelist, 1, "double", OP_INC), op_arg_dat(p_fy, 1, p_nodelist, 1, "double", OP_INC), op_arg_dat(p_fy, 2, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fy, 3, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fy, 4, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fy, 5, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fy, 6, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fy, 7, p_nodelist, 1, "double", OP_INC),
               op_arg_dat(p_fz, 0, p_nodelist, 1, "double", OP_INC), op_arg_dat(p_fz, 1, p_nodelist, 1, "double", OP_INC), op_arg_dat(p_fz, 2, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fz, 3, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fz, 4, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fz, 5, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fz, 6, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fz, 7, p_nodelist, 1, "double", OP_INC),
               op_arg_dat(p_determ, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_sigxx, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_sigyy, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_sigzz, -1, OP_ID, 1, "double", OP_READ)
               );


}

/******************************************/

static inline
void CalcFBHourglassForceForElems()
{
   /*************************************************
    *
    *     FUNCTION: Calculates the Flanagan-Belytschko anti-hourglass
    *               force.
    *
    *************************************************/

/*************************************************/
/*    compute the hourglass modes */
   op_par_loop(FBHourglassForceForElems, "CalcFBHourglassForceForElems", elems,
               op_arg_dat(p_xd, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_xd, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_xd, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_yd, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_yd, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_yd, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_zd, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_zd, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_zd, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_fx, 0, p_nodelist, 1, "double", OP_INC), op_arg_dat(p_fx, 1, p_nodelist, 1, "double", OP_INC), op_arg_dat(p_fx, 2, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fx, 3, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fx, 4, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fx, 5, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fx, 6, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fx, 7, p_nodelist, 1, "double", OP_INC),
               op_arg_dat(p_fy, 0, p_nodelist, 1, "double", OP_INC), op_arg_dat(p_fy, 1, p_nodelist, 1, "double", OP_INC), op_arg_dat(p_fy, 2, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fy, 3, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fy, 4, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fy, 5, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fy, 6, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fy, 7, p_nodelist, 1, "double", OP_INC),
               op_arg_dat(p_fz, 0, p_nodelist, 1, "double", OP_INC), op_arg_dat(p_fz, 1, p_nodelist, 1, "double", OP_INC), op_arg_dat(p_fz, 2, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fz, 3, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fz, 4, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fz, 5, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fz, 6, p_nodelist, 1, "double", OP_INC),op_arg_dat(p_fz, 7, p_nodelist, 1, "double", OP_INC),
               op_arg_dat(p_dvdx, -1, OP_ID, 8, "double", OP_READ),
               op_arg_dat(p_dvdy, -1, OP_ID, 8, "double", OP_READ),
               op_arg_dat(p_dvdz, -1, OP_ID, 8, "double", OP_READ),
               op_arg_dat(p_x8n, -1, OP_ID, 8, "double", OP_READ),
               op_arg_dat(p_y8n, -1, OP_ID, 8, "double", OP_READ),
               op_arg_dat(p_z8n, -1, OP_ID, 8, "double", OP_READ),
               op_arg_dat(p_determ, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_ss, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_elemMass, -1, OP_ID, 1, "double", OP_READ)
   );


   // MPI_Barrier(MPI_COMM_WORLD);
}

/******************************************/

static inline
void CalcHourglassControlForElems()
{
   op_par_loop(CalcVolumeDerivatives, "CalcVolumeDerivatives", elems,
               op_arg_dat(p_x, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_x, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_x, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_y, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_y, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_y, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_z, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_z, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_z, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_dvdx, -1, OP_ID, 8, "double", OP_WRITE),
               op_arg_dat(p_dvdy, -1, OP_ID, 8, "double", OP_WRITE),
               op_arg_dat(p_dvdz, -1, OP_ID, 8, "double", OP_WRITE),
               op_arg_dat(p_x8n, -1, OP_ID, 8, "double", OP_WRITE),
               op_arg_dat(p_y8n, -1, OP_ID, 8, "double", OP_WRITE),
               op_arg_dat(p_z8n, -1, OP_ID, 8, "double", OP_WRITE),
               op_arg_dat(p_v, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_determ, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_volo, -1, OP_ID, 1, "double", OP_READ)
               );
   if ( m_hgcoef > double(0.) ) {
      CalcFBHourglassForceForElems() ;
   }

   return ;
}

/******************************************/

static inline
void CalcVolumeForceForElems()
{
   if (m_numElem != 0) {
      /* Sum contributions to total stress tensor */
      // op_print("Init stress");
      InitStressTermsForElems();

      // call elemlib stress integration loop to produce nodal forces from
      // material stresses.
      // std::cout<<"Integrate stress at " << myRank << "\n";
      IntegrateStressForElems() ;

      // check for negative element volume
      op_par_loop(CheckForNegativeElementVolume, "CheckForNegativeElementVolume", elems,
                  op_arg_dat(p_determ, -1, OP_ID, 1, "double", OP_READ));

      CalcHourglassControlForElems() ;
   }
}

/******************************************/

static inline void CalcForceForNodes()
{
   // op_print("Force to zero");
   op_par_loop(setForceToZero, "setForceToZero", nodes,
               op_arg_dat(p_fx, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_fy, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_fz, -1, OP_ID, 1, "double", OP_WRITE));


  /* Calcforce calls partial, force, hourq */

  CalcVolumeForceForElems() ;

}

/******************************************/

static inline
void CalcAccelerationForNodes()
{
   op_par_loop(CalcAccelForNodes, "CalcAccelForNodes", nodes,
               op_arg_dat(p_xdd, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_ydd, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_zdd, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_fx, -1, OP_ID, 1, "double", OP_READ), op_arg_dat(p_fy, -1, OP_ID, 1, "double", OP_READ), op_arg_dat(p_fz, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_nodalMass, -1, OP_ID, 1, "double", OP_READ)
               );

}

/******************************************/

static inline
void ApplyAccelerationBoundaryConditionsForNodes()
{
   //Possible improvement would be to check if the current rank has a boundary node at pos 0
   op_par_loop(BoundaryX, "BoundaryX", nodes,
               op_arg_dat(p_xdd, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_t_symmX, -1, OP_ID, 1, "int", OP_READ)
               );


   op_par_loop(BoundaryY, "BoundaryY", nodes,
               op_arg_dat(p_ydd, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_t_symmY, -1, OP_ID, 1, "int", OP_READ)
               );


   op_par_loop(BoundaryZ, "BoundaryZ", nodes,
               op_arg_dat(p_zdd, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_t_symmZ, -1, OP_ID, 1, "int", OP_READ)
               );
}

/******************************************/

static inline
void CalcVelocityForNodes()
{
   op_par_loop(CalcVeloForNodes, "CalcVeloForNodes", nodes,
               op_arg_dat(p_xd, -1, OP_ID, 1, "double", OP_RW), op_arg_dat(p_yd, -1, OP_ID, 1, "double", OP_RW), op_arg_dat(p_zd, -1, OP_ID, 1, "double", OP_RW),
               op_arg_dat(p_xdd, -1, OP_ID, 1, "double", OP_READ), op_arg_dat(p_ydd, -1, OP_ID, 1, "double", OP_READ), op_arg_dat(p_zdd, -1, OP_ID, 1, "double", OP_READ),
               op_arg_gbl(&m_deltatime, 1, "double", OP_READ)
               );
}

/******************************************/

static inline
void CalcPositionForNodes()
{
   op_par_loop(CalcPosForNodes, "CalcPosForNodes", nodes,
               op_arg_dat(p_x, -1, OP_ID, 1, "double", OP_INC), op_arg_dat(p_y, -1, OP_ID, 1, "double", OP_INC), op_arg_dat(p_z, -1, OP_ID, 1, "double", OP_INC),
               op_arg_dat(p_xd, -1, OP_ID, 1, "double", OP_READ), op_arg_dat(p_yd, -1, OP_ID, 1, "double", OP_READ), op_arg_dat(p_zd, -1, OP_ID, 1, "double", OP_READ),
               op_arg_gbl(&m_deltatime, 1, "double", OP_READ)
               );
}

/******************************************/

static inline
void LagrangeNodal()
{
  /* time of boundary condition evaluation is beginning of step for force and
   * acceleration boundary conditions. */
//   op_print("Force For Nodes");
   CalcForceForNodes();

   // std::cout << "Acceleration at " << myRank << "\n";
   CalcAccelerationForNodes();
   // std::cout << "BCs at " << myRank << "\n";
   ApplyAccelerationBoundaryConditionsForNodes();

   // std::cout << "Velocity at " << myRank << "\n";
   CalcVelocityForNodes() ;
   // std::cout << "Calc Pos at " << myRank << "\n";
   CalcPositionForNodes();
   return;
}

/******************************************/

static inline
Real_t CalcElemVolume( const Real_t x0, const Real_t x1,
               const Real_t x2, const Real_t x3,
               const Real_t x4, const Real_t x5,
               const Real_t x6, const Real_t x7,
               const Real_t y0, const Real_t y1,
               const Real_t y2, const Real_t y3,
               const Real_t y4, const Real_t y5,
               const Real_t y6, const Real_t y7,
               const Real_t z0, const Real_t z1,
               const Real_t z2, const Real_t z3,
               const Real_t z4, const Real_t z5,
               const Real_t z6, const Real_t z7 )
{
  Real_t twelveth = Real_t(1.0)/Real_t(12.0);

  Real_t dx61 = x6 - x1;
  Real_t dy61 = y6 - y1;
  Real_t dz61 = z6 - z1;

  Real_t dx70 = x7 - x0;
  Real_t dy70 = y7 - y0;
  Real_t dz70 = z7 - z0;

  Real_t dx63 = x6 - x3;
  Real_t dy63 = y6 - y3;
  Real_t dz63 = z6 - z3;

  Real_t dx20 = x2 - x0;
  Real_t dy20 = y2 - y0;
  Real_t dz20 = z2 - z0;

  Real_t dx50 = x5 - x0;
  Real_t dy50 = y5 - y0;
  Real_t dz50 = z5 - z0;

  Real_t dx64 = x6 - x4;
  Real_t dy64 = y6 - y4;
  Real_t dz64 = z6 - z4;

  Real_t dx31 = x3 - x1;
  Real_t dy31 = y3 - y1;
  Real_t dz31 = z3 - z1;

  Real_t dx72 = x7 - x2;
  Real_t dy72 = y7 - y2;
  Real_t dz72 = z7 - z2;

  Real_t dx43 = x4 - x3;
  Real_t dy43 = y4 - y3;
  Real_t dz43 = z4 - z3;

  Real_t dx57 = x5 - x7;
  Real_t dy57 = y5 - y7;
  Real_t dz57 = z5 - z7;

  Real_t dx14 = x1 - x4;
  Real_t dy14 = y1 - y4;
  Real_t dz14 = z1 - z4;

  Real_t dx25 = x2 - x5;
  Real_t dy25 = y2 - y5;
  Real_t dz25 = z2 - z5;

#define TRIPLE_PRODUCT(x1, y1, z1, x2, y2, z2, x3, y3, z3) \
   ((x1)*((y2)*(z3) - (z2)*(y3)) + (x2)*((z1)*(y3) - (y1)*(z3)) + (x3)*((y1)*(z2) - (z1)*(y2)))

  Real_t volume =
    TRIPLE_PRODUCT(dx31 + dx72, dx63, dx20,
       dy31 + dy72, dy63, dy20,
       dz31 + dz72, dz63, dz20) +
    TRIPLE_PRODUCT(dx43 + dx57, dx64, dx70,
       dy43 + dy57, dy64, dy70,
       dz43 + dz57, dz64, dz70) +
    TRIPLE_PRODUCT(dx14 + dx25, dx61, dx50,
       dy14 + dy25, dy61, dy50,
       dz14 + dz25, dz61, dz50);

#undef TRIPLE_PRODUCT

  volume *= twelveth;

  return volume ;
}

/******************************************/

//inline
Real_t CalcElemVolume( const Real_t x[8], const Real_t y[8], const Real_t z[8] )
{
return CalcElemVolume( x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7],
                       y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7],
                       z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7]);
}

//static inline
void CalcKinematicsForElems()
{
   op_par_loop(CalcKinematicsForElem, "CalcKinematicsForElem", elems,
               op_arg_dat(p_x, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_x, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_x, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_y, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_y, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_y, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_z, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_z, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_z, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_xd, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_xd, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_xd, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_yd, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_yd, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_yd, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_zd, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_zd, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_zd, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_dxx, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_dyy, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_dzz, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_vnew, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_volo, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_delv, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_v, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_arealg, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_gbl(&m_deltatime, 1, "double", OP_READ)
   );
}

/******************************************/

static inline
void CalcLagrangeElements()
{
   if (m_numElem > 0) {
      CalcKinematicsForElems() ;
      //TODO review should be ok
      op_par_loop(CalcLagrangeElemRemaining, "CalcLagrangeElemRemaining", elems,
                  op_arg_dat(p_dxx, -1, OP_ID, 1, "double", OP_RW), op_arg_dat(p_dyy, -1, OP_ID, 1, "double", OP_RW), op_arg_dat(p_dzz, -1, OP_ID, 1, "double", OP_RW),
                  op_arg_dat(p_vdov, -1, OP_ID, 1, "double", OP_WRITE),
                  op_arg_dat(p_vnew, -1, OP_ID, 1, "double", OP_READ)
                  );
      // element loop to do some stuff not included in the elemlib function.
   }
}

/******************************************/

static inline
void CalcMonotonicQGradientsForElems()
{
   op_par_loop(CalcMonotonicQGradientsForElem, "CalcMonotonicQGradientsForElem", elems,
               op_arg_dat(p_x, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_x, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_x, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_x, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_y, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_y, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_y, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_y, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_z, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_z, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_z, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_z, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_xd, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_xd, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_xd, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_xd, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_yd, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_yd, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_yd, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_yd, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_zd, 0, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_zd, 1, p_nodelist, 1, "double", OP_READ), op_arg_dat(p_zd, 2, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 3, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 4, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 5, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 6, p_nodelist, 1, "double", OP_READ),op_arg_dat(p_zd, 7, p_nodelist, 1, "double", OP_READ),
               op_arg_dat(p_volo, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_vnew, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_delx_zeta, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_delv_zeta, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_delv_xi, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_delx_xi, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_delx_eta, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_delv_eta, -1, OP_ID, 1, "double", OP_WRITE)
               );
}

/******************************************/

static inline
void CalcMonotonicQRegionForElems()
{
   op_par_loop(CalcMonotonicQRegionForElem, "CalcMonotonicQRegionForElem",elems,
               op_arg_dat(p_delv_xi, -1, OP_ID, 1, "double", OP_READ), op_arg_dat(p_delv_xi, 0, p_lxim, 1, "double", OP_READ),op_arg_dat(p_delv_xi, 0, p_lxip, 1, "double", OP_READ),
               op_arg_dat(p_delv_eta, -1, OP_ID, 1, "double", OP_READ), op_arg_dat(p_delv_eta, 0, p_letam, 1, "double", OP_READ),op_arg_dat(p_delv_eta, 0, p_letap, 1, "double", OP_READ),
               op_arg_dat(p_delv_zeta, -1, OP_ID, 1, "double", OP_READ), op_arg_dat(p_delv_zeta, 0, p_lzetam, 1, "double", OP_READ),op_arg_dat(p_delv_zeta, 0, p_lzetap, 1, "double", OP_READ),
               op_arg_dat(p_delx_xi, -1, OP_ID, 1, "double", OP_READ), op_arg_dat(p_delx_eta, -1, OP_ID, 1, "double", OP_READ), op_arg_dat(p_delx_zeta, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_elemBC, -1, OP_ID, 1, "int", OP_READ),
               op_arg_dat(p_vdov, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_qq, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_ql, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_elemMass, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_volo, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_vnew, -1, OP_ID, 1, "double", OP_READ)
               );

}

/******************************************/

static inline
void CalcMonotonicQForElems()
{
   //
   // calculate the monotonic q for all regions
   //
   // The OP2 version does not support multiple regions yet
   // Code is left here for future reference
   // for (int r=0 ; r<domain.numReg() ; ++r) {
   // for (int r=0 ; r<m_numReg ; ++r) {
      // if (domain.regElemSize(r) > 0) {
      // if (m_regElemSize[r] > 0) {
      if(m_numElem > 0){
         // CalcMonotonicQRegionForElems(domain, r, ptiny) ;
         CalcMonotonicQRegionForElems() ;
      }
      // }
}


/******************************************/

static inline
void CalcQForElems()
{
   //
   // MONOTONIC Q option
   //
   if (m_numElem != 0) {

      /* Calculate velocity gradients */
      CalcMonotonicQGradientsForElems();

      CalcMonotonicQForElems();
      // Free up memory
      // domain.DeallocateGradients();
      // DeallocGradients();

      /* Don't allow excessive artificial viscosity */
      op_par_loop(NoExcessiveArtificialViscosity, "NoExcessiveArtificialViscosity", elems,
                  op_arg_dat(p_q, -1, OP_ID, 1, "double", OP_READ));
   }
}

/******************************************/

static inline
void CalcPressureForElemsHalfstep()
{
   op_par_loop(CalcHalfStepBVC, "CalcHalfStepBVC", elems,
               op_arg_dat(p_bvc, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_compHalfStep, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_pbvc, -1, OP_ID, 1, "double", OP_WRITE));

   //NOTE changed p_pHalfStep to write review
   op_par_loop(CalcPHalfstep, "CalcPHalfstep", elems,
               op_arg_dat(p_pHalfStep, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_bvc, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_e_new, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_vnewc, -1, OP_ID, 1, "double", OP_READ)
               );
}

static inline
void CalcPressureForElems()
{
   op_par_loop(CalcBVC, "CalcBVC", elems,
               op_arg_dat(p_bvc, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_compression, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_pbvc, -1, OP_ID, 1, "double", OP_WRITE));


   op_par_loop(CalcPNew, "CalcPNew", elems,
               op_arg_dat(p_p_new, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_bvc, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_e_new, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_vnewc, -1, OP_ID, 1, "double", OP_READ)
               );
}

/******************************************/

static inline
void CalcEnergyForElems(double* p_new, double* e_new, double* q_new,
                        double* bvc, double* pbvc,
                        double* p_old, double* e_old, double* q_old,
                        double* compression, double* compHalfStep,
                        double* vnewc, double* work, double* delvc, double pmin,
                        double p_cut, double  e_cut, double q_cut, double emin,
                        double* qq_old, double* ql_old,
                        double rho0,
                        double eosvmax,
                        int length, int *regElemList)
{
   // double *pHalfStep = Allocate<double>(length) ;


   op_par_loop(CalcNewE, "CalcNewE", elems,
               op_arg_dat(p_e_new, -1, OP_ID, 1, "double", OP_WRITE),
               op_arg_dat(p_e_old, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_delvc, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_p_old, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_q_old, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_work, -1, OP_ID, 1, "double", OP_READ)
   );

   CalcPressureForElemsHalfstep();

   //NOTE p_e_new may be an INC
   op_par_loop(CalcNewEStep2, "CalcNewEStep2", elems,
               op_arg_dat(p_compHalfStep, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_delvc, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_q_new, -1, OP_ID, 1, "double", OP_RW),
               op_arg_dat(p_pbvc, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_e_new, -1, OP_ID, 1, "double", OP_RW),
               op_arg_dat(p_bvc, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_pHalfStep, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_ql_old, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_qq_old, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_p_old, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_q_old, -1, OP_ID, 1, "double", OP_READ)
   );

   //NOTE p_e_new may be an INC
   op_par_loop(CalcNewEStep3, "CalcNewEStep3", elems,
               op_arg_dat(p_e_new, -1, OP_ID, 1, "double", OP_RW),
               op_arg_dat(p_work, -1, OP_ID, 1, "double", OP_READ)
   );

   CalcPressureForElems();

   //NOTE p_e_new may be an INC
   op_par_loop(CalcNewEStep4, "CalcNewEStep4", elems,
               op_arg_dat(p_delvc, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_pbvc, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_e_new, -1, OP_ID, 1, "double", OP_RW),
               op_arg_dat(p_vnewc, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_bvc, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_p_new, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_ql_old, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_qq_old, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_p_old, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_q_old, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_q_new, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_pHalfStep, -1, OP_ID, 1, "double", OP_READ)
   );

   CalcPressureForElems();
   //NOTE p_q_new could probably be a write
   op_par_loop(CalcQNew, "CalcQNew", elems,
               op_arg_dat(p_delvc, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_pbvc, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_e_new, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_vnewc, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_bvc, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_p_new, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_q_new, -1, OP_ID, 1, "double", OP_RW),
               op_arg_dat(p_ql_old, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_qq_old, -1, OP_ID, 1, "double", OP_READ)
               );


   return ;
}

/******************************************/

static inline
void CalcSoundSpeedForElems()
{
   op_par_loop(CalcSoundSpeedForElem, "CalcSoundSpeedForElem", elems,
               op_arg_dat(p_pbvc, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_e_new, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_vnewc, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_bvc, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_p_new, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_ss, -1, OP_ID, 1, "double", OP_WRITE)
               );
}

/******************************************/

static inline
void EvalEOSForElems(double *vnewc,
                     int numElemReg, int *regElemList, int rep)
{

   double  e_cut = m_e_cut ;
   double  p_cut = m_p_cut ;
   double  ss4o3 = m_ss4o3 ;
   double  q_cut = m_q_cut ;

   double eosvmax = m_eosvmax ;
   double eosvmin = m_eosvmin ;
   double pmin    = m_pmin ;
   double emin    = m_emin ;
   double rho0    = m_refdens ;

   // These temporaries will be of different size for
   // each call (due to different sized region element
   // lists)

   //loop to add load imbalance based on region number
   for(int j = 0; j < rep; j++) {
      /* compress data, minimal set */

         op_par_loop(CopyEOSValsIntoArray, "CopyEOSValsIntoArray", elems,
                     op_arg_dat(p_e_old, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_e, -1, OP_ID, 1, "double", OP_READ),
                     op_arg_dat(p_delvc, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_delv, -1, OP_ID, 1, "double", OP_READ),
                     op_arg_dat(p_p_old, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_p, -1, OP_ID, 1, "double", OP_READ),
                     op_arg_dat(p_q_old, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_q, -1, OP_ID, 1, "double", OP_READ),
                     op_arg_dat(p_qq_old, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_qq, -1, OP_ID, 1, "double", OP_READ),
                     op_arg_dat(p_ql_old, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_ql, -1, OP_ID, 1, "double", OP_READ));

         op_par_loop(CalcHalfSteps, "CalcHalfSteps", elems,
                     op_arg_dat(p_compression, -1, OP_ID, 1, "double", OP_WRITE),
                     op_arg_dat(p_vnewc, -1, OP_ID, 1, "double", OP_READ),
                     op_arg_dat(p_delvc, -1, OP_ID, 1, "double", OP_READ),
                     op_arg_dat(p_compHalfStep, -1, OP_ID, 1, "double", OP_WRITE)
         );


      /* Check for v > eosvmax or v < eosvmin */
         if ( eosvmin != double(0.) ) {
            op_par_loop(CheckEOSLowerBound, "CheckEOSLowerBound", elems,
                        op_arg_dat(p_vnewc, -1, OP_ID, 1, "double", OP_READ),
                        op_arg_dat(p_compHalfStep, -1, OP_ID, 1, "double", OP_WRITE),
                        op_arg_dat(p_compression, -1, OP_ID, 1, "double", OP_READ)
            );

         }
         if ( eosvmax != double(0.) ) {
            op_par_loop(CheckEOSUpperBound, "CheckEOSUpperBound", elems,
                        op_arg_dat(p_vnewc, -1, OP_ID, 1, "double", OP_READ),
                        op_arg_dat(p_compHalfStep, -1, OP_ID, 1, "double", OP_WRITE),
                        op_arg_dat(p_compression, -1, OP_ID, 1, "double", OP_WRITE),
                        op_arg_dat(p_p_old, -1, OP_ID, 1, "double", OP_WRITE)
            );

         }
         op_par_loop(CalcEOSWork, "CalcEOSWork", elems,
                     op_arg_dat(p_work, -1, OP_ID, 1 , "double", OP_WRITE));

      CalcEnergyForElems(p_new, e_new, q_new, bvc, pbvc,
                         p_old, e_old,  q_old, compression, compHalfStep,
                         vnewc, work,  delvc, pmin,
                         p_cut, e_cut, q_cut, emin,
                         qq_old, ql_old, rho0, eosvmax,
                         numElemReg, regElemList);
   }

   op_par_loop(CopyTempEOSVarsBack, "CopyTempEOSVarsBack", elems,
               op_arg_dat(p_p, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_p_new, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_e, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_e_new, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_q, -1, OP_ID, 1, "double", OP_WRITE), op_arg_dat(p_q_new, -1, OP_ID, 1, "double", OP_READ)
   );

   CalcSoundSpeedForElems() ;


}

/******************************************/

static inline
void ApplyMaterialPropertiesForElems()
{
   if (m_numElem != 0) {

      // MPI_Barrier(MPI_COMM_WORLD);
      // op_print("CopyVelo");
      op_par_loop(CopyVelocityToTempArray, "CopyVelocityToTempArray", elems,
                  op_arg_dat(p_vnewc, -1, OP_ID, 1, "double", OP_WRITE),
                  op_arg_dat(p_vnew, -1, OP_ID, 1, "double", OP_READ));

       // Bound the updated relative volumes with eosvmin/max
      // MPI_Barrier(MPI_COMM_WORLD);
      // op_print("Lower Bound");
       if (m_eosvmin != double(0.)) {
         op_par_loop(ApplyLowerBoundToVelocity, "ApplyLowerBoundToVelocity", elems,
                  op_arg_dat(p_vnewc, -1, OP_ID, 1, "double", OP_WRITE)
                  );

       }
      // MPI_Barrier(MPI_COMM_WORLD);
      // op_print("Upper bound");
       if (m_eosvmax != double(0.)) {
         op_par_loop(ApplyUpperBoundToVelocity, "ApplyUpperBoundToVelocity", elems,
                  op_arg_dat(p_vnewc, -1, OP_ID, 1, "double", OP_WRITE)
                  );
       }

       // This check may not make perfect sense in LULESH, but
       // it's representative of something in the full code -
       // just leave it in, please
      // MPI_Barrier(MPI_COMM_WORLD);
      // op_print("Ale3d");
      op_par_loop(ALE3DRelevantCheck, "ALE3DRelevantCheck", elems,
                  op_arg_dat(p_v, -1, OP_ID, 1, "double", OP_READ)
                  );



   //  for (int r=0 ; r<m_numReg ; r++) {
   //     int numElemReg = m_regElemSize[r];
   //     int *regElemList = m_regElemlist[r];
   //     int rep;
   //     //Determine load imbalance for this region
   //     //round down the number with lowest cost
   //    //  if(r < domain.numReg()/2)
   //     if(r < m_numReg/2)
   //  rep = 1;
   //     //you don't get an expensive region unless you at least have 5 regions
   //     else if(r < (m_numReg - (m_numReg+15)/20))
   //       rep = 1 + m_cost;
   //     //very expensive regions
   //     else
   //  rep = 10 * (1+ domain.cost());
   //  rep = 10 * (1+ m_cost);
   //     EvalEOSForElems(domain, vnewc, numElemReg, regElemList, rep);
   //    MPI_Barrier(MPI_COMM_WORLD);
   //    op_print("Eval EOS");
       EvalEOSForElems(vnewc, m_regElemSize[0], m_regElemlist[0], 1);
   //  }

   //  Release(&vnewc) ;
  }
}

/******************************************/

static inline
void UpdateVolumesForElems()
{
   if (m_numElem != 0) {
      op_par_loop(updateVolumesForElem, "UpdateVolumesForElems", elems,
               op_arg_dat(p_vnew, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_v, -1, OP_ID, 1, "double", OP_WRITE)
               );

   }
   return ;
}

/******************************************/

static inline
void LagrangeElements()
{
   // op_print("Calc Lagrange Elements");
  CalcLagrangeElements() ;

  /* Calculate Q.  (Monotonic q option requires communication) */
   //   MPI_Barrier(MPI_COMM_WORLD);
   //   op_print("Q For Elems");
  CalcQForElems() ;

   // MPI_Barrier(MPI_COMM_WORLD);
   // op_print("Material props");
  ApplyMaterialPropertiesForElems() ;

// MPI_Barrier(MPI_COMM_WORLD);
//    op_print("Update Vols for Elems");
  UpdateVolumesForElems() ;
}

/******************************************/

static inline
void CalcCourantConstraintForElems()
{

   op_par_loop(CalcCourantConstraint, "CalcCourantConstraint", elems,
               op_arg_dat(p_ss, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_vdov, -1, OP_ID, 1, "double", OP_READ),
               op_arg_dat(p_arealg, -1, OP_ID, 1, "double", OP_READ),
               op_arg_gbl(&m_dtcourant, 1, "double", OP_MIN)
   );

   return ;

}

/******************************************/

static inline
void CalcHydroConstraintForElems()
{

   op_par_loop(CalcHydroConstraint, "CalcHydroConstraint", elems,
               op_arg_dat(p_vdov, -1, OP_ID, 1, "double", OP_READ),
               op_arg_gbl(&m_dthydro, 1, "double", OP_MIN)
   );
   return ;
}

/******************************************/

static inline
void CalcTimeConstraintsForElems() {

   // Initialize conditions to a very large value
   m_dtcourant = 1.0e+20;
   m_dthydro = 1.0e+20;

   // for (int r=0 ; r < domain.numReg() ; ++r) {
   for (int r=0 ; r < m_numReg ; ++r) {
      /* evaluate time constraint */
      CalcCourantConstraintForElems() ;

      /* check hydro constraint */
      CalcHydroConstraintForElems() ;


   }
}

/******************************************/

static inline
void LagrangeLeapFrog()
{
   /* calculate nodal forces, accelerations, velocities, positions, with
    * applied boundary conditions and slide surface considerations */
   // op_print("Nodal");
   LagrangeNodal();

   // op_dump_dat()

   /* calculate element quantities (i.e. velocity gradient & q), and update
    * material states */
   // MPI_Barrier(MPI_COMM_WORLD);
   // op_print("Elements");
   LagrangeElements();

   CalcTimeConstraintsForElems();
}


/******************************************/

/* Helper function for converting strings to ints, with error checking */
template<typename IntT>
int StrToInt(const char *token, IntT *retVal)
{
   const char *c ;
   char *endptr ;
   const int decimal_base = 10 ;

   if (token == NULL)
      return 0 ;

   c = token ;
   *retVal = strtol(c, &endptr, decimal_base) ;
   if((endptr != c) && ((*endptr == ' ') || (*endptr == '\0')))
      return 1 ;
   else
      return 0 ;
}

static void PrintCommandLineOptions(char *execname, int myRank)
{
   if (myRank == 0) {

      printf("Usage: %s [opts]\n", execname);
      printf(" where [opts] is one or more of:\n");
      printf(" -q              : quiet mode - suppress all stdout\n");
      printf(" -i <iterations> : number of cycles to run\n");
      printf(" -s <size>       : length of cube mesh along side\n");
      printf(" -r <numregions> : Number of distinct regions (def: 11)\n");
      printf(" -b <balance>    : Load balance between regions of a domain (def: 1)\n");
      printf(" -c <cost>       : Extra cost of more expensive regions (def: 1)\n");
      printf(" -f <numfiles>   : Number of files to split viz dump into (def: (np+10)/9)\n");
      printf(" -p              : Print out progress\n");
      printf(" -v              : Output viz file (requires compiling with -DVIZ_MESH\n");
      printf(" -h              : This message\n");
      printf("\n\n");
   }
}

static void ParseError(const char *message, int myRank)
{
   if (myRank == 0) {
      printf("%s\n", message);
      MPI_Abort(MPI_COMM_WORLD, -1);

   }
}

void ParseCommandLineOptions(int argc, char *argv[],
                             int myRank, struct cmdLineOpts *opts)
{
   if(argc > 1) {
      int i = 1;

      while(i < argc) {
         int ok;
         /* -i <iterations> */
         if(strcmp(argv[i], "-i") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -i", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->its));
            if(!ok) {
               ParseError("Parse Error on option -i integer value required after argument\n", myRank);
            }
            i+=2;
         }
         /* -s <size, sidelength> */
         else if(strcmp(argv[i], "-s") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -s\n", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->nx));
            if(!ok) {
               ParseError("Parse Error on option -s integer value required after argument\n", myRank);
            }
            i+=2;
         }
    /* -r <numregions> */
         else if (strcmp(argv[i], "-r") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -r\n", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->numReg));
            if (!ok) {
               ParseError("Parse Error on option -r integer value required after argument\n", myRank);
            }
            i+=2;
         }
    /* -f <numfilepieces> */
         else if (strcmp(argv[i], "-f") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -f\n", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->numFiles));
            if (!ok) {
               ParseError("Parse Error on option -f integer value required after argument\n", myRank);
            }
            i+=2;
         }
         /* -p */
         else if (strcmp(argv[i], "-p") == 0) {
            opts->showProg = 1;
            i++;
         }
         /* -q */
         else if (strcmp(argv[i], "-q") == 0) {
            opts->quiet = 1;
            i++;
         }
         else if (strcmp(argv[i], "-b") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -b\n", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->balance));
            if (!ok) {
               ParseError("Parse Error on option -b integer value required after argument\n", myRank);
            }
            i+=2;
         }
         else if (strcmp(argv[i], "-c") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -c\n", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->cost));
            if (!ok) {
               ParseError("Parse Error on option -c integer value required after argument\n", myRank);
            }
            i+=2;
         }
         /* -v */
         else if (strcmp(argv[i], "-v") == 0) {
            #if VIZ_MESH
            opts->viz = 1;
            #else
            ParseError("Use of -v requires compiling with -DVIZ_MESH\n", myRank);
            #endif
            i++;
         }
         /* -h */
         else if (strcmp(argv[i], "-h") == 0) {
            PrintCommandLineOptions(argv[0], myRank);
            MPI_Abort(MPI_COMM_WORLD, 0);

         }
         else if(strcmp(argv[i], "OP_NO_REALLOC" ) == 0){
            i++;
         }
         else if(strcmp(argv[i], "-l") == 0){
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -l\n", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->itert));
            if(!ok) {
               ParseError("Parse Error on option -l integer value required after argument\n", myRank);
            }
            i+=2;
         }
            else if(strcmp(argv[i], "-z") == 0){
            if (i+1 >= argc) {
               ParseError("Missing integer argument to -z\n", myRank);
            }
            ok = StrToInt(argv[i+1], &(opts->fnm));
            if(!ok) {
               ParseError("Parse Error on option -z integer value required after argument\n", myRank);
            }
            i+=2;
         }

         /* -t partiotioning library <S,PK,PG,PKG,K> */
         else if(strcmp(argv[i], "-t") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing letter argumnet to -t\n", myRank);
            }
            if ( strcmp(argv[i+1], "S") == 0 ){
               opts->partition = Partition_S;
            }
            if ( strcmp(argv[i+1], "PK") == 0 ){
               opts->partition = Partition_PK;
            }
            if ( strcmp(argv[i+1], "PG") == 0 ){
               opts->partition = Partition_PG;
            }
            if ( strcmp(argv[i+1], "PGK") == 0 ){
               opts->partition = Partition_PGK;
            }
            if ( strcmp(argv[i+1], "K") == 0 ){
               opts->partition = Partition_K;
            }
            i+=2;
         }

         /* -t initialiseation method <S,PK,PG,PKG,K> */
         else if(strcmp(argv[i], "-d") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing letter argumnet to -i\n", myRank);
            }
            if ( strcmp(argv[i+1], "S") == 0 ){
               opts->creation = Creation_Root;
            }
            if ( strcmp(argv[i+1], "P") == 0 ){
               opts->creation = Creation_Parallel;
            }
            i+=2;
         }

         /* -t initialiseation method <S,PK,PG,PKG,K> */
         else if(strcmp(argv[i], "-tm") == 0) {
            if (i+1 >= argc) {
               ParseError("Missing letter argumnet to -i\n", myRank);
            }
            if ( strcmp(argv[i+1], "S") == 0 ){
               opts->time = Show_Time;
            }
            if ( strcmp(argv[i+1], "H") == 0 ){
               opts->time = Hide_Time;
            }
            i+=2;
         }


         /* -f <numfilepieces> */
         else {
            char msg[80];
            i++;
            // PrintCommandLineOptions(argv[0], myRank);
            // sprintf(msg, "ERROR: Unknown command line argument: %s\n", argv[i]);
            // ParseError(msg, myRank);
         }
      }
   }
}

void VerifyAndWriteFinalOutput(double elapsed_time,
                               int cycle,
                               int nx,
                               int numRanks, double* e)
{
   // GrindTime1 only takes a single domain into account, and is thus a good way to measure
   // processor speed indepdendent of MPI parallelism.
   // GrindTime2 takes into account speedups from MPI parallelism.
   // Cast to 64-bit integer to avoid overflows.
   Int8_t nx8 = nx;
   double grindTime1 = ((elapsed_time*1e6)/cycle)/(nx8*nx8*nx8);
   double grindTime2 = ((elapsed_time*1e6)/cycle)/(nx8*nx8*nx8*numRanks);

   int ElemId = 0;
   std::cout << "Run completed:\n";
   std::cout << "   Problem size        =  " << nx       << "\n";
   std::cout << "   MPI tasks           =  " << numRanks << "\n";
   std::cout << "   Iteration count     =  " << cycle << "\n";
   std::cout << "   Final Origin Energy =  ";
   std::cout << std::scientific << std::setprecision(6);
   // std::cout << std::setw(12) << locDom.e(ElemId) << "\n";
   std::cout << std::setw(12) << e[ElemId] << "\n";

   double   MaxAbsDiff = double(0.0);
   double TotalAbsDiff = double(0.0);
   double   MaxRelDiff = double(0.0);

   for (int j=0; j<nx; ++j) {
      for (int k=j+1; k<nx; ++k) {
         double AbsDiff = FABS(e[j*nx+k]-e[k*nx+j]);
         TotalAbsDiff  += AbsDiff;

         if (MaxAbsDiff <AbsDiff) MaxAbsDiff = AbsDiff;

         double RelDiff = AbsDiff / e[k*nx+j];

         if (MaxRelDiff <RelDiff)  MaxRelDiff = RelDiff;
      }
   }

   // Quick symmetry check
   std::cout << "   Testing Plane 0 of Energy Array on rank 0:\n";
   std::cout << "        MaxAbsDiff   = " << std::setw(12) << MaxAbsDiff   << "\n";
   std::cout << "        TotalAbsDiff = " << std::setw(12) << TotalAbsDiff << "\n";
   std::cout << "        MaxRelDiff   = " << std::setw(12) << MaxRelDiff   << "\n";

   // Timing information
   std::cout.unsetf(std::ios_base::floatfield);
   std::cout << std::setprecision(4);
   std::cout << "\nElapsed time         = " << std::setw(10) << elapsed_time << " (s)\n";
   std::cout << std::setprecision(8);
   std::cout << "Grind time (us/z/c)  = "  << std::setw(10) << grindTime1 << " (per dom)  ("
             << std::setw(10) << elapsed_time << " overall)\n";
   std::cout << "FOM                  = " << std::setw(10) << 1000.0/grindTime2 << " (z/s)\n\n";

   return ;
}

void InitMeshDecomp(int numRanks, int myRank,
                    int *col, int *row, int *plane, int *side)
{
   int testProcs;
   int dx, dy, dz;
   int myDom;

   // Assume cube processor layout for now
   testProcs = int(cbrt(double(numRanks))+0.5) ;
   if (testProcs*testProcs*testProcs != numRanks) {
      printf("Num processors must be a cube of an integer (1, 8, 27, ...)\n") ;
      MPI_Abort(MPI_COMM_WORLD, -1) ;
   }
   if (sizeof(double) != 4 && sizeof(double) != 8) {
      printf("MPI operations only support float and double right now...\n");
      MPI_Abort(MPI_COMM_WORLD, -1) ;
   }
   if (MAX_FIELDS_PER_MPI_COMM > CACHE_COHERENCE_PAD_REAL) {
      printf("corner element comm buffers too small.  Fix code.\n") ;
      MPI_Abort(MPI_COMM_WORLD, -1) ;
   }

   dx = testProcs ;
   dy = testProcs ;
   dz = testProcs ;

   // temporary test
   if (dx*dy*dz != numRanks) {
      printf("error -- must have as many domains as procs\n") ;
      MPI_Abort(MPI_COMM_WORLD, -1) ;
   }
   int remainder = dx*dy*dz % numRanks ;
   if (myRank < remainder) {
      myDom = myRank*( 1+ (dx*dy*dz / numRanks)) ;
   }
   else {
      myDom = remainder*( 1+ (dx*dy*dz / numRanks)) +
         (myRank - remainder)*(dx*dy*dz/numRanks) ;
   }

   *col = myDom % dx ;
   *row = (myDom / dx) % dy ;
   *plane = myDom / (dx*dy) ;
   *side = testProcs;

   return;
}


void writeFileADHFJ(int myRank,int m_numNodes,double *x,double *y, double *z, char *file ){
   int fileLength = strlen(file);
   int numLength=1;

   if (myRank<=0){
      numLength=1;
   }else{
      numLength=(int)(ceil(  log10(myRank+1)  ));
   }

   char *newName = (char*)malloc( (numLength+fileLength+2)*sizeof(char) );
   strcpy(newName, file);
   sprintf(&newName[fileLength],"%d",myRank);
   newName[fileLength+numLength] ='\0';

   printf(newName);
   FILE* ptr = fopen(newName,"w");
   for (int i = 0; i < m_numNodes; i++)
   {
      /* code */
      fprintf(ptr,"%d, %.6f, %.6f, %.6f\n",myRank,x[i],y[i],z[i]);
   }
   fclose(ptr);
   free(newName);
}


/******************************************/

int main(int argc, char *argv[])
{
   op_init(argc, argv, 1);

   int numRanks ;

// std::cout << "CHENCK THAT STUFF IS PRINTING\n";
//    std::cout << "USING MPI HERE JUST TO BE SURE\n";

   MPI_Comm_size(MPI_COMM_WORLD, &numRanks) ;
   MPI_Comm_rank(MPI_COMM_WORLD, &myRank) ;

   /* Set defaults that can be overridden by command line opts */
   opts.its = 9999999; // Iterations
   opts.nx  = 30; //Size
   opts.numReg = 1;
   opts.numFiles = (int)(numRanks+10)/9;
   opts.showProg = 0;
   opts.quiet = 0;
   opts.viz = 0;
   opts.balance = 1;
   opts.cost = 1;
   opts.itert = 1;
   m_numReg = 1;


   ParseCommandLineOptions(argc, argv, myRank, &opts);

   if ((myRank == 0) && (opts.quiet == 0)) {
      std::cout << "Running problem size " << opts.nx << "^3 per domain until completion\n";
      std::cout << "Num processors: "      << numRanks << "\n";
      #if _OPENMP
      std::cout << "Num threads: " << omp_get_max_threads() << "\n";
      #endif
      std::cout << "Total number of elements: " << (opts.nx*opts.nx*opts.nx) << " \n\n";
      std::cout << "To run other sizes, use -s <integer>.\n";
      std::cout << "To run a fixed number of iterations, use -i <integer>.\n";
      std::cout << "To run a more or less balanced region set, use -b <integer>.\n";
      std::cout << "To change the relative costs of regions, use -c <integer>.\n";
      std::cout << "To print out progress, use -p\n";
      std::cout << "To write an output file for VisIt, use -v\n";
      std::cout << "See help (-h) for more options\n\n";
   }

   // Set up the mesh and decompose. Assumes regular cubes for now
   int col, row, plane, side;
   // InitMeshDecomp(numRanks, myRank, &col, &row, &plane, &side);


   //locDom
   int edgeElems = opts.nx;
   int edgeNodes = edgeElems+1;

   sizeX = edgeElems ;
   sizeY = edgeElems ;
   sizeZ = edgeElems ;
   m_numElem = edgeElems * edgeElems * edgeElems;
   m_numNode = edgeNodes * edgeNodes * edgeNodes;


   //generate the global size
   g_numElem = opts.nx*opts.nx*opts.nx;
   g_numNode = (opts.nx+1)*(opts.nx+1)*(opts.nx+1);

   //Compute the local size
   m_numElem =  compute_local_size(g_numElem,(Int8_t)numRanks,myRank);
   m_numNode =  compute_local_size(g_numNode,(Int8_t)numRanks,myRank);

   allocateElems();
   allocateNodes();
   if (opts.creation == Creation_Parallel){
      initialise(myRank,opts.nx,1,opts.numReg,opts.balance, opts.cost,(Int8_t)numRanks);
   }else if (opts.creation == Creation_Root){
      if (myRank == 0){
         allocateGlobalElems();
         allocateGlobalNodes();

         initialiseSingular(0,0,0,opts.nx,1,opts.numReg,opts.balance, opts.cost,(Int8_t)numRanks);
      }
      distributeGlobalElems(numRanks,g_numElem,m_numElem,g_numNode,m_numNode);
   }
   MPI_Bcast(&m_deltatime,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

   // writeFileADHFJ(myRank,m_numNode,x,y,z,"/home/joseph/3rdYear/testingFolder/partitioning/Before/Node");



   MPI_Barrier(MPI_COMM_WORLD);

   // (int comm_size, int g_numElem , int numElem, int g_numNode, int numNode)
   // distributeGlobalElems(numRanks,g_numElem,m_numElem,g_numNode,m_numNode);

   MPI_Barrier(MPI_COMM_WORLD);


   m_dtfixed = double(-1.0e-6) ; // Negative means use courant condition
   m_stoptime  = double(1.0e-2); // *double(edgeElems*tp/45.0) ;

   // Initial conditions
   m_deltatimemultlb = double(1.1) ;
   m_deltatimemultub = double(1.2) ;
   m_dtcourant = double(1.0e+20) ;
   m_dthydro   = double(1.0e+20) ;
   m_dtmax     = double(1.0e-2) ;
   m_time    = double(0.) ;
   m_cycle   = int(0) ;

   m_gamma_t[0] = double( 1.);
   m_gamma_t[1] = double( 1.);
   m_gamma_t[2] = double(-1.);
   m_gamma_t[3] = double(-1.);
   m_gamma_t[4] = double(-1.);
   m_gamma_t[5] = double(-1.);
   m_gamma_t[6] = double( 1.);
   m_gamma_t[7] = double( 1.);
   m_gamma_t[8] = double( 1.);
   m_gamma_t[9] = double(-1.);
   m_gamma_t[10] = double(-1.);
   m_gamma_t[11] = double( 1.);
   m_gamma_t[12] = double(-1.);
   m_gamma_t[13] = double( 1.);
   m_gamma_t[14] = double( 1.);
   m_gamma_t[15] = double(-1.);
   m_gamma_t[16] = double( 1.);
   m_gamma_t[17] = double(-1.);
   m_gamma_t[18] = double( 1.);
   m_gamma_t[19] = double(-1.);
   m_gamma_t[20] = double( 1.);
   m_gamma_t[21] = double(-1.);
   m_gamma_t[22] = double( 1.);
   m_gamma_t[23] = double(-1.);
   m_gamma_t[24] = double(-1.);
   m_gamma_t[25] = double( 1.);
   m_gamma_t[26] = double(-1.);
   m_gamma_t[27] = double( 1.);
   m_gamma_t[28] = double( 1.);
   m_gamma_t[29] = double(-1.);
   m_gamma_t[30] = double( 1.);
   m_gamma_t[31] = double(-1.);

   //Declare Constants
   op_decl_const(1, "double", &m_e_cut);
   op_decl_const(1, "double", &m_p_cut);
   op_decl_const(1, "double", &m_q_cut);
   op_decl_const(1, "double", &m_v_cut);
   op_decl_const(1, "double", &m_u_cut);

   op_decl_const(1, "double", &m_hgcoef);
   op_decl_const(1, "double", &m_ss4o3);
   op_decl_const(1, "double", &m_qstop);
   op_decl_const(1, "double", &m_monoq_max_slope);
   op_decl_const(1, "double", &m_monoq_limiter_mult);
   op_decl_const(1, "double", &m_qlc_monoq);
   op_decl_const(1, "double", &m_qqc_monoq);
   op_decl_const(1, "double", &m_qqc);
   op_decl_const(1, "double", &m_eosvmax);
   op_decl_const(1, "double", &m_eosvmin);
   op_decl_const(1, "double", &m_pmin);
   op_decl_const(1, "double", &m_emin);
   op_decl_const(1, "double", &m_dvovmax);
   op_decl_const(1, "double", &m_refdens);
   op_decl_const(1, "double", &m_qqc2);

   op_decl_const(1, "double", &m_ptiny);
   op_decl_const(32, "double", m_gamma_t);
   op_decl_const(1, "double", &m_twelfth);
   op_decl_const(1, "double", &m_sixth);
   op_decl_const(1, "double", &m_c1s);
   op_decl_const(1, "double", &m_ssc_thresh);
   op_decl_const(1, "double", &m_ssc_low);

   char file[] = FILE_NAME_PATH;
   // readAndInitVars(file);
   // readOp2VarsFromHDF5(file);

   initOp2Vars();
   int siz = op_get_size(elems);
   std::cout <<  "Global Size: " << siz << ", Local Elems: " << elems->size << ", Local Nodes: "<< nodes->size << "\n";
   MPI_Barrier(MPI_COMM_WORLD);
   if (myRank==0){printf("\n\n");}
   op_diagnostic_output();

   switch (opts.partition)
   {
      case Partition_S:
         op_partition("PTSCOTCH", "KWAY", nodes, p_nodelist, p_loc);
         break;
      case Partition_PK:
         op_partition("PARMETIS", "KWAY", nodes, p_nodelist, p_loc);
         break;
      case Partition_PG:
         op_partition("PARMETIS", "GEOM", nodes, p_nodelist, p_loc);
         break;
      case Partition_PGK:
         op_partition("PARMETIS", "GEOMKWAY", nodes, p_nodelist, p_loc);
         break;
      case Partition_K:
         op_partition("KAHIP", "KWAY", nodes, p_nodelist, p_loc);
         break;
   }

   // op_dump_to_hdf5("file_out.h5");
   // op_write_const_hdf5("deltatime", 1, "double", (char*)&m_deltatime, FILE_NAME_PATH);

   MPI_Barrier(MPI_COMM_WORLD);
   if (myRank==0){printf("\n\n");}

   if (myRank==0){
      printf("size %d\n",p_x->set->size);
      printf("exec_size x%d\n",p_x->set->exec_size);
      printf("nonexec_size x%d\n",p_x->set->nonexec_size);

      printf("exec_size y%d\n",p_y->set->exec_size);
      printf("nonexec_size y %d\n",p_y->set->nonexec_size);

      printf("exec_size z %d\n",p_z->set->exec_size);
      printf("nonexec_size z%d\n",p_z->set->nonexec_size);
      printf("m_nodes %d \n",m_numNode);
   }



   // SHOW THE CUBE AS DIVIDED
   // writeFileADHFJ(myRank,
   //    p_x->set->size+p_x->set->exec_size+p_x->set->nonexec_size,
   //    (double *)p_x->data,
   //    (double *)p_y->data,
   //    (double *)p_z->data,"/home/joseph/3rdYear/testingFolder/partitioning/After/Node");

   // writeFileADHFJ(myRank,p_x->set->nonexec_size+p_x->set->exec_size,(double *)p_x->data,(double *)p_y->data,(double *)p_z->data,"/home/joseph/3rdYear/testingFolder/partitioning/All/Node");
   // double *outputX = (double*) malloc(m_numNode*sizeof(double));
   // double *outputY = (double*) malloc(m_numNode*sizeof(double));
   // double *outputZ = (double*) malloc(m_numNode*sizeof(double));
   // op_fetch_data(p_x, outputX);
   // op_fetch_data(p_y, outputY);
   // op_fetch_data(p_z, outputZ);
   // writeFileADHFJ(myRank,m_numNode,outputX,outputY,outputZ,"/home/joseph/3rdYear/testingFolder/partitioning/After/Node");
   // free(outputX);
   // free(outputY);
   // free(outputZ);

   MPI_Barrier(MPI_COMM_WORLD);


   // BEGIN timestep to solution */
   double cpu_t1, cpu_t2, wall_t1, wall_t2;
   op_timers(&cpu_t1, &wall_t1);

   // !Main Loop
   while((m_time < m_stoptime) && (m_cycle < opts.its)) {
      TimeIncrement(myRank) ;
      LagrangeLeapFrog() ;

      if ((opts.showProg != 0) && (opts.quiet == 0) && (myRank==0)) {
         std::cout << "rank = " << myRank       << ", "
                   << "cycle = " << m_cycle       << ", "
                   << std::scientific
                   << "time = " << double(m_time) << ", "
                   << "dt="     << double(m_deltatime) << "\n";
         std::cout.unsetf(std::ios_base::floatfield);
      }
   }

   op_timers(&cpu_t2, &wall_t2);
   double walltime = wall_t2 - wall_t1;

   // Write out final viz file */
   // Currently no support for visualisation
   // if (opts.viz) {
      // DumpToVisit(*locDom, opts.numFiles, myRank, numRanks) ;
   // }

   // writeFileADHFJ(myRank,
   //    p_x->set->size+p_x->set->exec_size+p_x->set->nonexec_size,
   //    (double *)p_x->data,
   //    (double *)p_y->data,
   //    (double *)p_z->data,"/home/joseph/3rdYear/testingFolder/partitioning/All/Node");

   double *verify_e = (double*) malloc(m_numElem*sizeof(double));
   op_fetch_data(p_e, verify_e);

   if ((myRank == 0) && (opts.quiet == 0)) {
      VerifyAndWriteFinalOutput(walltime, m_cycle, opts.nx, numRanks, verify_e);
   }

   if (opts.time){op_timing_output();}



   op_exit();

   return 0 ;
}
