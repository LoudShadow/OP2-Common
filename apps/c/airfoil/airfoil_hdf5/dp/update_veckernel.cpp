//
// auto-generated by op2.py
//

//user function
inline void update(const double *qold, double *q, double *res, const double *adt, double *rms){
  double del, adti,rmsl;

  rmsl = 0.0f;
  adti = 1.0f/(*adt);

  for (int n=0; n<4; n++) {
    del    = adti*res[n];
    q[n]   = qold[n] - del;
    res[n] = 0.0f;
    rmsl  += del*del;
  }
  *rms += rmsl;
}


// host stub function
void op_par_loop_update(char const *name, op_set set,
  op_arg arg0,
  op_arg arg1,
  op_arg arg2,
  op_arg arg3,
  op_arg arg4){

  int nargs = 5;
  op_arg args[5];

  args[0] = arg0;
  args[1] = arg1;
  args[2] = arg2;
  args[3] = arg3;
  args[4] = arg4;
  //create aligned pointers for dats
  ALIGNED_double const double * __restrict__ ptr0 = (double *) arg0.data;
  __assume_aligned(ptr0,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr1 = (double *) arg1.data;
  __assume_aligned(ptr1,double_ALIGN);
  ALIGNED_double       double * __restrict__ ptr2 = (double *) arg2.data;
  __assume_aligned(ptr2,double_ALIGN);
  ALIGNED_double const double * __restrict__ ptr3 = (double *) arg3.data;
  __assume_aligned(ptr3,double_ALIGN);

  // initialise timers
  double cpu_t1, cpu_t2, wall_t1, wall_t2;
  op_timing_realloc(4);
  op_timers_core(&cpu_t1, &wall_t1);


  if (OP_diags>2) {
    printf(" kernel routine w/o indirection:  update");
  }

  int exec_size = op_mpi_halo_exchanges(set, nargs, args);

  if (exec_size >0) {

    #ifdef VECTORIZE
    #pragma novector
    for ( int n=0; n<(exec_size/SIMD_VEC)*SIMD_VEC; n+=SIMD_VEC ){
      double dat4[SIMD_VEC] = {0.0};
      #pragma simd
      for ( int i=0; i<SIMD_VEC; i++ ){
        update(
          &(ptr0)[4 * (n+i)],
          &(ptr1)[4 * (n+i)],
          &(ptr2)[4 * (n+i)],
          &(ptr3)[1 * (n+i)],
          &dat4[i]);
      }
      for ( int i=0; i<SIMD_VEC; i++ ){
        *(double*)arg4.data += dat4[i];
      }
    }
    //remainder
    for ( int n=(exec_size/SIMD_VEC)*SIMD_VEC; n<exec_size; n++ ){
    #else
    for ( int n=0; n<exec_size; n++ ){
    #endif
      update(
        &(ptr0)[4*n],
        &(ptr1)[4*n],
        &(ptr2)[4*n],
        &(ptr3)[1*n],
        (double*)arg4.data);
    }
  }

  // combine reduction data
  op_mpi_reduce(&arg4,(double*)arg4.data);
  op_mpi_set_dirtybit(nargs, args);

  // update kernel record
  op_timers_core(&cpu_t2, &wall_t2);
  OP_kernels[4].name      = name;
  OP_kernels[4].count    += 1;
  OP_kernels[4].time     += wall_t2 - wall_t1;
  OP_kernels[4].transfer += (float)set->size * arg0.size;
  OP_kernels[4].transfer += (float)set->size * arg1.size * 2.0f;
  OP_kernels[4].transfer += (float)set->size * arg2.size * 2.0f;
  OP_kernels[4].transfer += (float)set->size * arg3.size;
}
