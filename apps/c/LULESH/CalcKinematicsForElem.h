// static inline
// double CalcElemVolumeSingleArg( const double x0, const double x1,
//                const double x2, const double x3,
//                const double x4, const double x5,
//                const double x6, const double x7,
//                const double y0, const double y1,
//                const double y2, const double y3,
//                const double y4, const double y5,
//                const double y6, const double y7,
//                const double z0, const double z1,
//                const double z2, const double z3,
//                const double z4, const double z5,
//                const double z6, const double z7 )
// {
//   double twelveth = double(1.0)/double(12.0);

//   double dx61 = x6 - x1;
//   double dy61 = y6 - y1;
//   double dz61 = z6 - z1;

//   double dx70 = x7 - x0;
//   double dy70 = y7 - y0;
//   double dz70 = z7 - z0;

//   double dx63 = x6 - x3;
//   double dy63 = y6 - y3;
//   double dz63 = z6 - z3;

//   double dx20 = x2 - x0;
//   double dy20 = y2 - y0;
//   double dz20 = z2 - z0;

//   double dx50 = x5 - x0;
//   double dy50 = y5 - y0;
//   double dz50 = z5 - z0;

//   double dx64 = x6 - x4;
//   double dy64 = y6 - y4;
//   double dz64 = z6 - z4;

//   double dx31 = x3 - x1;
//   double dy31 = y3 - y1;
//   double dz31 = z3 - z1;

//   double dx72 = x7 - x2;
//   double dy72 = y7 - y2;
//   double dz72 = z7 - z2;

//   double dx43 = x4 - x3;
//   double dy43 = y4 - y3;
//   double dz43 = z4 - z3;

//   double dx57 = x5 - x7;
//   double dy57 = y5 - y7;
//   double dz57 = z5 - z7;

//   double dx14 = x1 - x4;
//   double dy14 = y1 - y4;
//   double dz14 = z1 - z4;

//   double dx25 = x2 - x5;
//   double dy25 = y2 - y5;
//   double dz25 = z2 - z5;

// #define TRIPLE_PRODUCT(x1, y1, z1, x2, y2, z2, x3, y3, z3) \
//    ((x1)*((y2)*(z3) - (z2)*(y3)) + (x2)*((z1)*(y3) - (y1)*(z3)) + (x3)*((y1)*(z2) - (z1)*(y2)))

//   double volume =
//     TRIPLE_PRODUCT(dx31 + dx72, dx63, dx20,
//        dy31 + dy72, dy63, dy20,
//        dz31 + dz72, dz63, dz20) +
//     TRIPLE_PRODUCT(dx43 + dx57, dx64, dx70,
//        dy43 + dy57, dy64, dy70,
//        dz43 + dz57, dz64, dz70) +
//     TRIPLE_PRODUCT(dx14 + dx25, dx61, dx50,
//        dy14 + dy25, dy61, dy50,
//        dz14 + dz25, dz61, dz50);

// #undef TRIPLE_PRODUCT

//   volume *= twelveth;

//   return volume ;
// }

// double CalcElemVolumes( const double x[8], const double y[8], const double z[8] ){
// return CalcElemVolumeSingleArg( x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7],
//                        y[0], y[1], y[2], y[3], y[4], y[5], y[6], y[7],
//                        z[0], z[1], z[2], z[3], z[4], z[5], z[6], z[7]);
// }
// #include "const.h"

// static inline //OP_FUN_PREFIX
// double AreaFacet( const double x0, const double x1,
//                  const double x2, const double x3,
//                  const double y0, const double y1,
//                  const double y2, const double y3,
//                  const double z0, const double z1,
//                  const double z2, const double z3)
// {
//    double fx = (x2 - x0) - (x3 - x1);
//    double fy = (y2 - y0) - (y3 - y1);
//    double fz = (z2 - z0) - (z3 - z1);
//    double gx = (x2 - x0) + (x3 - x1);
//    double gy = (y2 - y0) + (y3 - y1);
//    double gz = (z2 - z0) + (z3 - z1);
//    double area =
//       (fx * fx + fy * fy + fz * fz) *
//       (gx * gx + gy * gy + gz * gz) -
//       (fx * gx + fy * gy + fz * gz) *
//       (fx * gx + fy * gy + fz * gz);
//    return area ;
// }

// static inline
// double CalcElemCharacteristicLengtht( const double x[8],
//                                      const double y[8],
//                                      const double z[8],
//                                      const double volume)
// {
//    double a, charLength = double(0.0);

//    a = AreaFacet(x[0],x[1],x[2],x[3],
//                 y[0],y[1],y[2],y[3],
//                 z[0],z[1],z[2],z[3]) ;
//    charLength = std::fmax(a,charLength) ;

//    a = AreaFacet(x[4],x[5],x[6],x[7],
//                 y[4],y[5],y[6],y[7],
//                 z[4],z[5],z[6],z[7]) ;
//    charLength = std::fmax(a,charLength) ;

//    a = AreaFacet(x[0],x[1],x[5],x[4],
//                 y[0],y[1],y[5],y[4],
//                 z[0],z[1],z[5],z[4]) ;
//    charLength = std::fmax(a,charLength) ;

//    a = AreaFacet(x[1],x[2],x[6],x[5],
//                 y[1],y[2],y[6],y[5],
//                 z[1],z[2],z[6],z[5]) ;
//    charLength = std::fmax(a,charLength) ;

//    a = AreaFacet(x[2],x[3],x[7],x[6],
//                 y[2],y[3],y[7],y[6],
//                 z[2],z[3],z[7],z[6]) ;
//    charLength = std::fmax(a,charLength) ;

//    a = AreaFacet(x[3],x[0],x[4],x[7],
//                 y[3],y[0],y[4],y[7],
//                 z[3],z[0],z[4],z[7]) ;
//    charLength = std::fmax(a,charLength) ;

//    charLength = double(4.0) * volume / sqrt(charLength);

//    return charLength;
// }

// static inline
// void CalcElemShapeFunctionDerivativest( double const x[],
//                                        double const y[],
//                                        double const z[],
//                                        double b[][8],
//                                        double* const volume )
// {
//   const double x0 = x[0] ;   const double x1 = x[1] ;
//   const double x2 = x[2] ;   const double x3 = x[3] ;
//   const double x4 = x[4] ;   const double x5 = x[5] ;
//   const double x6 = x[6] ;   const double x7 = x[7] ;

//   const double y0 = y[0] ;   const double y1 = y[1] ;
//   const double y2 = y[2] ;   const double y3 = y[3] ;
//   const double y4 = y[4] ;   const double y5 = y[5] ;
//   const double y6 = y[6] ;   const double y7 = y[7] ;

//   const double z0 = z[0] ;   const double z1 = z[1] ;
//   const double z2 = z[2] ;   const double z3 = z[3] ;
//   const double z4 = z[4] ;   const double z5 = z[5] ;
//   const double z6 = z[6] ;   const double z7 = z[7] ;

//   double fjxxi, fjxet, fjxze;
//   double fjyxi, fjyet, fjyze;
//   double fjzxi, fjzet, fjzze;
//   double cjxxi, cjxet, cjxze;
//   double cjyxi, cjyet, cjyze;
//   double cjzxi, cjzet, cjzze;

//   fjxxi = double(.125) * ( (x6-x0) + (x5-x3) - (x7-x1) - (x4-x2) );
//   fjxet = double(.125) * ( (x6-x0) - (x5-x3) + (x7-x1) - (x4-x2) );
//   fjxze = double(.125) * ( (x6-x0) + (x5-x3) + (x7-x1) + (x4-x2) );

//   fjyxi = double(.125) * ( (y6-y0) + (y5-y3) - (y7-y1) - (y4-y2) );
//   fjyet = double(.125) * ( (y6-y0) - (y5-y3) + (y7-y1) - (y4-y2) );
//   fjyze = double(.125) * ( (y6-y0) + (y5-y3) + (y7-y1) + (y4-y2) );

//   fjzxi = double(.125) * ( (z6-z0) + (z5-z3) - (z7-z1) - (z4-z2) );
//   fjzet = double(.125) * ( (z6-z0) - (z5-z3) + (z7-z1) - (z4-z2) );
//   fjzze = double(.125) * ( (z6-z0) + (z5-z3) + (z7-z1) + (z4-z2) );

//   /* compute cofactors */
//   cjxxi =    (fjyet * fjzze) - (fjzet * fjyze);
//   cjxet =  - (fjyxi * fjzze) + (fjzxi * fjyze);
//   cjxze =    (fjyxi * fjzet) - (fjzxi * fjyet);

//   cjyxi =  - (fjxet * fjzze) + (fjzet * fjxze);
//   cjyet =    (fjxxi * fjzze) - (fjzxi * fjxze);
//   cjyze =  - (fjxxi * fjzet) + (fjzxi * fjxet);

//   cjzxi =    (fjxet * fjyze) - (fjyet * fjxze);
//   cjzet =  - (fjxxi * fjyze) + (fjyxi * fjxze);
//   cjzze =    (fjxxi * fjyet) - (fjyxi * fjxet);

//   /* calculate partials :
//      this need only be done for l = 0,1,2,3   since , by symmetry ,
//      (6,7,4,5) = - (0,1,2,3) .
//   */
//   b[0][0] =   -  cjxxi  -  cjxet  -  cjxze;
//   b[0][1] =      cjxxi  -  cjxet  -  cjxze;
//   b[0][2] =      cjxxi  +  cjxet  -  cjxze;
//   b[0][3] =   -  cjxxi  +  cjxet  -  cjxze;
//   b[0][4] = -b[0][2];
//   b[0][5] = -b[0][3];
//   b[0][6] = -b[0][0];
//   b[0][7] = -b[0][1];

//   b[1][0] =   -  cjyxi  -  cjyet  -  cjyze;
//   b[1][1] =      cjyxi  -  cjyet  -  cjyze;
//   b[1][2] =      cjyxi  +  cjyet  -  cjyze;
//   b[1][3] =   -  cjyxi  +  cjyet  -  cjyze;
//   b[1][4] = -b[1][2];
//   b[1][5] = -b[1][3];
//   b[1][6] = -b[1][0];
//   b[1][7] = -b[1][1];

//   b[2][0] =   -  cjzxi  -  cjzet  -  cjzze;
//   b[2][1] =      cjzxi  -  cjzet  -  cjzze;
//   b[2][2] =      cjzxi  +  cjzet  -  cjzze;
//   b[2][3] =   -  cjzxi  +  cjzet  -  cjzze;
//   b[2][4] = -b[2][2];
//   b[2][5] = -b[2][3];
//   b[2][6] = -b[2][0];
//   b[2][7] = -b[2][1];

//   /* calculate jacobian determinant (volume) */
//   *volume = double(8.) * ( fjxet * cjxet + fjyet * cjyet + fjzet * cjzet);
// }

// static inline
// void CalcElemVelocityGradientt( const double* const xvel,
//                                 const double* const yvel,
//                                 const double* const zvel,
//                                 const double b[][8],
//                                 const double detJ,
//                                 double* const d )
// {
//   const double inv_detJ = double(1.0) / detJ ;
//   double dyddx, dxddy, dzddx, dxddz, dzddy, dyddz;
//   const double* const pfx = b[0];
//   const double* const pfy = b[1];
//   const double* const pfz = b[2];

//   d[0] = inv_detJ * ( pfx[0] * (xvel[0]-xvel[6])
//                      + pfx[1] * (xvel[1]-xvel[7])
//                      + pfx[2] * (xvel[2]-xvel[4])
//                      + pfx[3] * (xvel[3]-xvel[5]) );

//   d[1] = inv_detJ * ( pfy[0] * (yvel[0]-yvel[6])
//                      + pfy[1] * (yvel[1]-yvel[7])
//                      + pfy[2] * (yvel[2]-yvel[4])
//                      + pfy[3] * (yvel[3]-yvel[5]) );

//   d[2] = inv_detJ * ( pfz[0] * (zvel[0]-zvel[6])
//                      + pfz[1] * (zvel[1]-zvel[7])
//                      + pfz[2] * (zvel[2]-zvel[4])
//                      + pfz[3] * (zvel[3]-zvel[5]) );

//   dyddx  = inv_detJ * ( pfx[0] * (yvel[0]-yvel[6])
//                       + pfx[1] * (yvel[1]-yvel[7])
//                       + pfx[2] * (yvel[2]-yvel[4])
//                       + pfx[3] * (yvel[3]-yvel[5]) );

//   dxddy  = inv_detJ * ( pfy[0] * (xvel[0]-xvel[6])
//                       + pfy[1] * (xvel[1]-xvel[7])
//                       + pfy[2] * (xvel[2]-xvel[4])
//                       + pfy[3] * (xvel[3]-xvel[5]) );

//   dzddx  = inv_detJ * ( pfx[0] * (zvel[0]-zvel[6])
//                       + pfx[1] * (zvel[1]-zvel[7])
//                       + pfx[2] * (zvel[2]-zvel[4])
//                       + pfx[3] * (zvel[3]-zvel[5]) );

//   dxddz  = inv_detJ * ( pfz[0] * (xvel[0]-xvel[6])
//                       + pfz[1] * (xvel[1]-xvel[7])
//                       + pfz[2] * (xvel[2]-xvel[4])
//                       + pfz[3] * (xvel[3]-xvel[5]) );

//   dzddy  = inv_detJ * ( pfy[0] * (zvel[0]-zvel[6])
//                       + pfy[1] * (zvel[1]-zvel[7])
//                       + pfy[2] * (zvel[2]-zvel[4])
//                       + pfy[3] * (zvel[3]-zvel[5]) );

//   dyddz  = inv_detJ * ( pfz[0] * (yvel[0]-yvel[6])
//                       + pfz[1] * (yvel[1]-yvel[7])
//                       + pfz[2] * (yvel[2]-yvel[4])
//                       + pfz[3] * (yvel[3]-yvel[5]) );
//   d[5]  = double( .5) * ( dxddy + dyddx );
//   d[4]  = double( .5) * ( dxddz + dzddx );
//   d[3]  = double( .5) * ( dzddy + dyddz );
// }

inline void CalcKinematicsForElem(
                                const double *p_x0, const double *p_x1, const double *p_x2, const double *p_x3, const double *p_x4, const double *p_x5, const double *p_x6, const double *p_x7,
                                const double *p_y0, const double *p_y1, const double *p_y2, const double *p_y3, const double *p_y4, const double *p_y5, const double *p_y6, const double *p_y7,
                                const double *p_z0, const double *p_z1, const double *p_z2, const double *p_z3, const double *p_z4, const double *p_z5, const double *p_z6, const double *p_z7,
                                const double *p_xd0, const double *p_xd1, const double *p_xd2, const double *p_xd3, const double *p_xd4, const double *p_xd5, const double *p_xd6, const double *p_xd7,
                                const double *p_yd0, const double *p_yd1, const double *p_yd2, const double *p_yd3, const double *p_yd4, const double *p_yd5, const double *p_yd6, const double *p_yd7,
                                const double *p_zd0, const double *p_zd1, const double *p_zd2, const double *p_zd3, const double *p_zd4, const double *p_zd5, const double *p_zd6, const double *p_zd7,
                                double *dxx, double *dyy, double *dzz,
                                double *vnew,
                                const double *volo,
                                double *delv,
                                const double *v,
                                double *arealg,
                                const double *deltaTime

){

   double B[3][8] ; /** shape function derivatives */
   double D[6] ;
   double x_local[8] ;
   double y_local[8] ;
   double z_local[8] ;
   double detJ = double(0.0) ;

   double volume ;
   double relativeVolume ;

   x_local[0] = p_x0[0];
   x_local[1] = p_x1[0];
   x_local[2] = p_x2[0];
   x_local[3] = p_x3[0];
   x_local[4] = p_x4[0];
   x_local[5] = p_x5[0];
   x_local[6] = p_x6[0];
   x_local[7] = p_x7[0];

   y_local[0] = p_y0[0];
   y_local[1] = p_y1[0];
   y_local[2] = p_y2[0];
   y_local[3] = p_y3[0];
   y_local[4] = p_y4[0];
   y_local[5] = p_y5[0];
   y_local[6] = p_y6[0];
   y_local[7] = p_y7[0];

   z_local[0] = p_z0[0];
   z_local[1] = p_z1[0];
   z_local[2] = p_z2[0];
   z_local[3] = p_z3[0];
   z_local[4] = p_z4[0];
   z_local[5] = p_z5[0];
   z_local[6] = p_z6[0];
   z_local[7] = p_z7[0];

   // Calc Elem Volume

   double dx61 = p_x6[0] - p_x1[0];
   double dy61 = p_y6[0] - p_y1[0];
   double dz61 = p_z6[0] - p_z1[0];

   double dx70 = p_x7[0] - p_x0[0];
   double dy70 = p_y7[0] - p_y0[0];
   double dz70 = p_z7[0] - p_z0[0];

   double dx63 = p_x6[0] - p_x3[0];
   double dy63 = p_y6[0] - p_y3[0];
   double dz63 = p_z6[0] - p_z3[0];

   double dx20 = p_x2[0] - p_x0[0];
   double dy20 = p_y2[0] - p_y0[0];
   double dz20 = p_z2[0] - p_z0[0];

   double dx50 = p_x5[0] - p_x0[0];
   double dy50 = p_y5[0] - p_y0[0];
   double dz50 = p_z5[0] - p_z0[0];

   double dx64 = p_x6[0] - p_x4[0];
   double dy64 = p_y6[0] - p_y4[0];
   double dz64 = p_z6[0] - p_z4[0];

   double dx31 = p_x3[0] - p_x1[0];
   double dy31 = p_y3[0] - p_y1[0];
   double dz31 = p_z3[0] - p_z1[0];

   double dx72 = p_x7[0] - p_x2[0];
   double dy72 = p_y7[0] - p_y2[0];
   double dz72 = p_z7[0] - p_z2[0];

   double dx43 = p_x4[0] - p_x3[0];
   double dy43 = p_y4[0] - p_y3[0];
   double dz43 = p_z4[0] - p_z3[0];

   double dx57 = p_x5[0] - p_x7[0];
   double dy57 = p_y5[0] - p_y7[0];
   double dz57 = p_z5[0] - p_z7[0];

   double dx14 = p_x1[0] - p_x4[0];
   double dy14 = p_y1[0] - p_y4[0];
   double dz14 = p_z1[0] - p_z4[0];

   double dx25 = p_x2[0] - p_x5[0];
   double dy25 = p_y2[0] - p_y5[0];
   double dz25 = p_z2[0] - p_z5[0];

   #define TRIPLE_PRODUCT(x1, y1, z1, x2, y2, z2, x3, y3, z3) \
      ((x1)*((y2)*(z3) - (z2)*(y3)) + (x2)*((z1)*(y3) - (y1)*(z3)) + (x3)*((y1)*(z2) - (z1)*(y2)))

   double temp_volume =
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
   temp_volume *= m_twelfth;

   volume = temp_volume;
   // End Calc Elem Volume


   // volume = CalcElemVolumes(x_local, y_local, z_local );
   relativeVolume = volume / volo[0] ;
   vnew[0] = relativeVolume ;
   delv[0] = relativeVolume - v[0] ;


   // arealg[0] = CalcElemCharacteristicLengtht(x_local, y_local, z_local,
   //                                        volume);
   // Start CalcElemCharacteristicLength function
   double a, charLength = double(0.0);
   double fx,fy,fz,gx,gy,gz,area;
   // a = AreaFacet(x_local[0],x_local[1],x_local[2],x_local[3],
   //             y_local[0],y_local[1],y_local[2],y_local[3],
   //             z_local[0],z_local[1],z_local[2],z_local[3]) ;
   fx = (x_local[2] - x_local[0]) - (x_local[3] - x_local[1]);
   fy = (y_local[2] - y_local[0]) - (y_local[3] - y_local[1]);
   fz = (z_local[2] - z_local[0]) - (z_local[3] - z_local[1]);
   gx = (x_local[2] - x_local[0]) + (x_local[3] - x_local[1]);
   gy = (y_local[2] - y_local[0]) + (y_local[3] - y_local[1]);
   gz = (z_local[2] - z_local[0]) + (z_local[3] - z_local[1]);
   area =
      (fx * fx + fy * fy + fz * fz) *
      (gx * gx + gy * gy + gz * gz) -
      (fx * gx + fy * gy + fz * gz) *
      (fx * gx + fy * gy + fz * gz);
   a = area;
   charLength = std::fmax(a,charLength) ;

   // a = AreaFacet(x_local[4],x_local[5],x_local[6],x_local[7],
   //              y_local[4],y_local[5],y_local[6],y_local[7],
   //              z_local[4],z_local[5],z_local[6],z_local[7]) ;
   fx = (x_local[6] - x_local[4]) - (x_local[7] - x_local[5]);
   fy = (y_local[6] - y_local[4]) - (y_local[7] - y_local[5]);
   fz = (z_local[6] - z_local[4]) - (z_local[7] - z_local[5]);
   gx = (x_local[6] - x_local[4]) + (x_local[7] - x_local[5]);
   gy = (y_local[6] - y_local[4]) + (y_local[7] - y_local[5]);
   gz = (z_local[6] - z_local[4]) + (z_local[7] - z_local[5]);
   area =
      (fx * fx + fy * fy + fz * fz) *
      (gx * gx + gy * gy + gz * gz) -
      (fx * gx + fy * gy + fz * gz) *
      (fx * gx + fy * gy + fz * gz);
   a = area;
   charLength = std::fmax(a,charLength) ;

   // a = AreaFacet(x_local[0],x_local[1],x_local[5],x_local[4],
   //              y_local[0],y_local[1],y_local[5],y_local[4],
   //              z_local[0],z_local[1],z_local[5],z_local[4]) ;
   fx = (x_local[5] - x_local[0]) - (x_local[4] - x_local[1]);
   fy = (y_local[5] - y_local[0]) - (y_local[4] - y_local[1]);
   fz = (z_local[5] - z_local[0]) - (z_local[4] - z_local[1]);
   gx = (x_local[5] - x_local[0]) + (x_local[4] - x_local[1]);
   gy = (y_local[5] - y_local[0]) + (y_local[4] - y_local[1]);
   gz = (z_local[5] - z_local[0]) + (z_local[4] - z_local[1]);
   area =
      (fx * fx + fy * fy + fz * fz) *
      (gx * gx + gy * gy + gz * gz) -
      (fx * gx + fy * gy + fz * gz) *
      (fx * gx + fy * gy + fz * gz);
   a = area;
   charLength = std::fmax(a,charLength) ;

   // a = AreaFacet(x_local[1],x_local[2],x_local[6],x_local[5],
   //              y_local[1],y_local[2],y_local[6],y_local[5],
   //              z_local[1],z_local[2],z_local[6],z_local[5]) ;
   fx = (x_local[6] - x_local[1]) - (x_local[5] - x_local[2]);
   fy = (y_local[6] - y_local[1]) - (y_local[5] - y_local[2]);
   fz = (z_local[6] - z_local[1]) - (z_local[5] - z_local[2]);
   gx = (x_local[6] - x_local[1]) + (x_local[5] - x_local[2]);
   gy = (y_local[6] - y_local[1]) + (y_local[5] - y_local[2]);
   gz = (z_local[6] - z_local[1]) + (z_local[5] - z_local[2]);
   area =
      (fx * fx + fy * fy + fz * fz) *
      (gx * gx + gy * gy + gz * gz) -
      (fx * gx + fy * gy + fz * gz) *
      (fx * gx + fy * gy + fz * gz);
   a = area;
   charLength = std::fmax(a,charLength) ;

   // a = AreaFacet(x_local[2],x_local[3],x_local[7],x_local[6],
   //              y_local[2],y_local[3],y_local[7],y_local[6],
   //              z_local[2],z_local[3],z_local[7],z_local[6]) ;
   fx = (x_local[7] - x_local[2]) - (x_local[6] - x_local[3]);
   fy = (y_local[7] - y_local[2]) - (y_local[6] - y_local[3]);
   fz = (z_local[7] - z_local[2]) - (z_local[6] - z_local[3]);
   gx = (x_local[7] - x_local[2]) + (x_local[6] - x_local[3]);
   gy = (y_local[7] - y_local[2]) + (y_local[6] - y_local[3]);
   gz = (z_local[7] - z_local[2]) + (z_local[6] - z_local[3]);
   area =
      (fx * fx + fy * fy + fz * fz) *
      (gx * gx + gy * gy + gz * gz) -
      (fx * gx + fy * gy + fz * gz) *
      (fx * gx + fy * gy + fz * gz);
   a = area;
   charLength = std::fmax(a,charLength) ;

   // a = AreaFacet(x_local[3],x_local[0],x_local[4],x_local[7],
   //              y_local[3],y_local[0],y_local[4],y_local[7],
   //              z_local[3],z_local[0],z_local[4],z_local[7]) ;
   fx = (x_local[4] - x_local[3]) - (x_local[7] - x_local[0]);
   fy = (y_local[4] - y_local[3]) - (y_local[7] - y_local[0]);
   fz = (z_local[4] - z_local[3]) - (z_local[7] - z_local[0]);
   gx = (x_local[4] - x_local[3]) + (x_local[7] - x_local[0]);
   gy = (y_local[4] - y_local[3]) + (y_local[7] - y_local[0]);
   gz = (z_local[4] - z_local[3]) + (z_local[7] - z_local[0]);
   area =
      (fx * fx + fy * fy + fz * fz) *
      (gx * gx + gy * gy + gz * gz) -
      (fx * gx + fy * gy + fz * gz) *
      (fx * gx + fy * gy + fz * gz);
   a = area;
   charLength = std::fmax(a,charLength) ;

   charLength = double(4.0) * volume / sqrt(charLength);

   //End Function
   arealg[0] = charLength;

   double dt2 = double(0.5) * (*deltaTime);

   x_local[0] -= dt2 * p_xd0[0];
   y_local[0] -= dt2 * p_yd0[0];
   z_local[0] -= dt2 * p_zd0[0];

   x_local[1] -= dt2 * p_xd1[0];
   y_local[1] -= dt2 * p_yd1[0];
   z_local[1] -= dt2 * p_zd1[0];

   x_local[2] -= dt2 * p_xd2[0];
   y_local[2] -= dt2 * p_yd2[0];
   z_local[2] -= dt2 * p_zd2[0];

   x_local[3] -= dt2 * p_xd3[0];
   y_local[3] -= dt2 * p_yd3[0];
   z_local[3] -= dt2 * p_zd3[0];

   x_local[4] -= dt2 * p_xd4[0];
   y_local[4] -= dt2 * p_yd4[0];
   z_local[4] -= dt2 * p_zd4[0];

   x_local[5] -= dt2 * p_xd5[0];
   y_local[5] -= dt2 * p_yd5[0];
   z_local[5] -= dt2 * p_zd5[0];

   x_local[6] -= dt2 * p_xd6[0];
   y_local[6] -= dt2 * p_yd6[0];
   z_local[6] -= dt2 * p_zd6[0];

   x_local[7] -= dt2 * p_xd7[0];
   y_local[7] -= dt2 * p_yd7[0];
   z_local[7] -= dt2 * p_zd7[0];

   // CalcElemShapeFunctionDerivativest( x_local, y_local, z_local,
   //                                  B, &detJ );
   // Start CalcElemShapeFunctionDerivativest function
//      const double x0 = x_local[0] ;   const double x1 = x_local[1] ;
//   const double x2 = x_local[2] ;   const double x3 = x_local[3] ;
//   const double x4 = x_local[4] ;   const double x5 = x_local[5] ;
//   const double x6 = x_local[6] ;   const double x7 = x_local[7] ;

//   const double y0 = y_local[0] ;   const double y1 = y_local[1] ;
//   const double y2 = y_local[2] ;   const double y3 = y_local[3] ;
//   const double y4 = y_local[4] ;   const double y5 = y_local[5] ;
//   const double y6 = y_local[6] ;   const double y7 = y_local[7] ;

//   const double z0 = z_local[0] ;   const double z1 = z_local[1] ;
//   const double z2 = z_local[2] ;   const double z3 = z_local[3] ;
//   const double z4 = z_local[4] ;   const double z5 = z_local[5] ;
//   const double z6 = z_local[6] ;   const double z7 = z_local[7] ;

  double fjxxi, fjxet, fjxze;
  double fjyxi, fjyet, fjyze;
  double fjzxi, fjzet, fjzze;
  double cjxxi, cjxet, cjxze;
  double cjyxi, cjyet, cjyze;
  double cjzxi, cjzet, cjzze;

   fjxxi = double(.125) * ( (x_local[6]-x_local[0]) + (x_local[5]-x_local[3]) - (x_local[7]-x_local[1]) - (x_local[4]-x_local[2]) );
   fjxet = double(.125) * ( (x_local[6]-x_local[0]) - (x_local[5]-x_local[3]) + (x_local[7]-x_local[1]) - (x_local[4]-x_local[2]) );
   fjxze = double(.125) * ( (x_local[6]-x_local[0]) + (x_local[5]-x_local[3]) + (x_local[7]-x_local[1]) + (x_local[4]-x_local[2]) );

   fjyxi = double(.125) * ( (y_local[6]-y_local[0]) + (y_local[5]-y_local[3]) - (y_local[7]-y_local[1]) - (y_local[4]-y_local[2]) );
   fjyet = double(.125) * ( (y_local[6]-y_local[0]) - (y_local[5]-y_local[3]) + (y_local[7]-y_local[1]) - (y_local[4]-y_local[2]) );
   fjyze = double(.125) * ( (y_local[6]-y_local[0]) + (y_local[5]-y_local[3]) + (y_local[7]-y_local[1]) + (y_local[4]-y_local[2]) );

   fjzxi = double(.125) * ( (z_local[6]-z_local[0]) + (z_local[5]-z_local[3]) - (z_local[7]-z_local[1]) - (z_local[4]-z_local[2]) );
   fjzet = double(.125) * ( (z_local[6]-z_local[0]) - (z_local[5]-z_local[3]) + (z_local[7]-z_local[1]) - (z_local[4]-z_local[2]) );
   fjzze = double(.125) * ( (z_local[6]-z_local[0]) + (z_local[5]-z_local[3]) + (z_local[7]-z_local[1]) + (z_local[4]-z_local[2]) );


  /* compute cofactors */
  cjxxi =    (fjyet * fjzze) - (fjzet * fjyze);
  cjxet =  - (fjyxi * fjzze) + (fjzxi * fjyze);
  cjxze =    (fjyxi * fjzet) - (fjzxi * fjyet);

  cjyxi =  - (fjxet * fjzze) + (fjzet * fjxze);
  cjyet =    (fjxxi * fjzze) - (fjzxi * fjxze);
  cjyze =  - (fjxxi * fjzet) + (fjzxi * fjxet);

  cjzxi =    (fjxet * fjyze) - (fjyet * fjxze);
  cjzet =  - (fjxxi * fjyze) + (fjyxi * fjxze);
  cjzze =    (fjxxi * fjyet) - (fjyxi * fjxet);

  /* calculate partials :
     this need only be done for l = 0,1,2,3   since , by symmetry ,
     (6,7,4,5) = - (0,1,2,3) .
  */
  B[0][0] =   -  cjxxi  -  cjxet  -  cjxze;
  B[0][1] =      cjxxi  -  cjxet  -  cjxze;
  B[0][2] =      cjxxi  +  cjxet  -  cjxze;
  B[0][3] =   -  cjxxi  +  cjxet  -  cjxze;
  B[0][4] = -B[0][2];
  B[0][5] = -B[0][3];
  B[0][6] = -B[0][0];
  B[0][7] = -B[0][1];

  B[1][0] =   -  cjyxi  -  cjyet  -  cjyze;
  B[1][1] =      cjyxi  -  cjyet  -  cjyze;
  B[1][2] =      cjyxi  +  cjyet  -  cjyze;
  B[1][3] =   -  cjyxi  +  cjyet  -  cjyze;
  B[1][4] = -B[1][2];
  B[1][5] = -B[1][3];
  B[1][6] = -B[1][0];
  B[1][7] = -B[1][1];

  B[2][0] =   -  cjzxi  -  cjzet  -  cjzze;
  B[2][1] =      cjzxi  -  cjzet  -  cjzze;
  B[2][2] =      cjzxi  +  cjzet  -  cjzze;
  B[2][3] =   -  cjzxi  +  cjzet  -  cjzze;
  B[2][4] = -B[2][2];
  B[2][5] = -B[2][3];
  B[2][6] = -B[2][0];
  B[2][7] = -B[2][1];

  /* calculate jacobian determinant (volume) */
  detJ = double(8.) * ( fjxet * cjxet + fjyet * cjyet + fjzet * cjzet);


   // CalcElemVelocityGradientt( xd_local, yd_local, zd_local,
   //                            B, detJ, D );
   const double inv_detJ = double(1.0) / detJ ;
   double dyddx, dxddy, dzddx, dxddz, dzddy, dyddz;
   const double* const pfx = B[0];
   const double* const pfy = B[1];
   const double* const pfz = B[2];

   D[0] = inv_detJ * ( pfx[0] * (p_xd0[0]-p_xd6[0])
                        + pfx[1] * (p_xd1[0]-p_xd7[0])
                        + pfx[2] * (p_xd2[0]-p_xd4[0])
                        + pfx[3] * (p_xd3[0]-p_xd5[0]) );

   D[1] = inv_detJ * ( pfy[0] * (p_yd0[0]-p_yd6[0])
                        + pfy[1] * (p_yd1[0]-p_yd7[0])
                        + pfy[2] * (p_yd2[0]-p_yd4[0])
                        + pfy[3] * (p_yd3[0]-p_yd5[0]) );

   D[2] = inv_detJ * ( pfz[0] * (p_zd0[0]-p_zd6[0])
                        + pfz[1] * (p_zd1[0]-p_zd7[0])
                        + pfz[2] * (p_zd2[0]-p_zd4[0])
                        + pfz[3] * (p_zd3[0]-p_zd5[0]) );

   dyddx  = inv_detJ * ( pfx[0] * (p_yd0[0]-p_yd6[0])
                        + pfx[1] * (p_yd1[0]-p_yd7[0])
                        + pfx[2] * (p_yd2[0]-p_yd4[0])
                        + pfx[3] * (p_yd3[0]-p_yd5[0]) );

   dxddy  = inv_detJ * ( pfy[0] * (p_xd0[0]-p_xd6[0])
                        + pfy[1] * (p_xd1[0]-p_xd7[0])
                        + pfy[2] * (p_xd2[0]-p_xd4[0])
                        + pfy[3] * (p_xd3[0]-p_xd5[0]) );

   dzddx  = inv_detJ * ( pfx[0] * (p_zd0[0]-p_zd6[0])
                        + pfx[1] * (p_zd1[0]-p_zd7[0])
                        + pfx[2] * (p_zd2[0]-p_zd4[0])
                        + pfx[3] * (p_zd3[0]-p_zd5[0]) );

   dxddz  = inv_detJ * ( pfz[0] * (p_xd0[0]-p_xd6[0])
                        + pfz[1] * (p_xd1[0]-p_xd7[0])
                        + pfz[2] * (p_xd2[0]-p_xd4[0])
                        + pfz[3] * (p_xd3[0]-p_xd5[0]) );

   dzddy  = inv_detJ * ( pfy[0] * (p_zd0[0]-p_zd6[0])
                        + pfy[1] * (p_zd1[0]-p_zd7[0])
                        + pfy[2] * (p_zd2[0]-p_zd4[0])
                        + pfy[3] * (p_zd3[0]-p_zd5[0]) );

   dyddz  = inv_detJ * ( pfz[0] * (p_yd0[0]-p_yd6[0])
                        + pfz[1] * (p_yd1[0]-p_yd7[0])
                        + pfz[2] * (p_yd2[0]-p_yd4[0])
                        + pfz[3] * (p_yd3[0]-p_yd5[0]) );
   D[5]  = double( .5) * ( dxddy + dyddx );
   D[4]  = double( .5) * ( dxddz + dzddx );
   D[3]  = double( .5) * ( dzddy + dyddz );
   //ENd Function
   dxx[0] = D[0];
   dyy[0] = D[1];
   dzz[0] = D[2];
}