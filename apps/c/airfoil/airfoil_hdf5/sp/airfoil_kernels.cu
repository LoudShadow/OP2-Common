//
// auto-generated by op2.py
//

//header
#include "op_lib_cpp.h"
#include "op_cuda_rt_support.h"
#include "op_cuda_reduction.h"

//global constants
#ifndef MAX_CONST_SIZE
#define MAX_CONST_SIZE 128
#endif

__constant__ float gam;
__constant__ float gm1;
__constant__ float cfl;
__constant__ float eps;
__constant__ float mach;
__constant__ float alpha;
__constant__ float qinf[4];

void op_decl_const_char(int dim, char const *type,
int size, char *dat, char const *name){
  if (!strcmp(name,"gam")) {
    cutilSafeCall(cudaMemcpyToSymbol(gam, dat, dim*size));
  }
  else
  if (!strcmp(name,"gm1")) {
    cutilSafeCall(cudaMemcpyToSymbol(gm1, dat, dim*size));
  }
  else
  if (!strcmp(name,"cfl")) {
    cutilSafeCall(cudaMemcpyToSymbol(cfl, dat, dim*size));
  }
  else
  if (!strcmp(name,"eps")) {
    cutilSafeCall(cudaMemcpyToSymbol(eps, dat, dim*size));
  }
  else
  if (!strcmp(name,"mach")) {
    cutilSafeCall(cudaMemcpyToSymbol(mach, dat, dim*size));
  }
  else
  if (!strcmp(name,"alpha")) {
    cutilSafeCall(cudaMemcpyToSymbol(alpha, dat, dim*size));
  }
  else
  if (!strcmp(name,"qinf")) {
    cutilSafeCall(cudaMemcpyToSymbol(qinf, dat, dim*size));
  }
  else
  {
    printf("error: unknown const name\n"); exit(1);
  }
}

//user kernel files
#include "save_soln_kernel.cu"
#include "adt_calc_kernel.cu"
#include "res_calc_kernel.cu"
#include "bres_calc_kernel.cu"
#include "update_kernel.cu"
