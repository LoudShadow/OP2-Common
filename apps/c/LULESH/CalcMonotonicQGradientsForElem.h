inline void CalcMonotonicQGradientsForElem(
    const double *p_x0, const double *p_x1, const double *p_x2, const double *p_x3, const double *p_x4, const double *p_x5, const double *p_x6, const double *p_x7,
    const double *p_y0, const double *p_y1, const double *p_y2, const double *p_y3, const double *p_y4, const double *p_y5, const double *p_y6, const double *p_y7,
    const double *p_z0, const double *p_z1, const double *p_z2, const double *p_z3, const double *p_z4, const double *p_z5, const double *p_z6, const double *p_z7,
    const double *p_xd0, const double *p_xd1, const double *p_xd2, const double *p_xd3, const double *p_xd4, const double *p_xd5, const double *p_xd6, const double *p_xd7,
    const double *p_yd0, const double *p_yd1, const double *p_yd2, const double *p_yd3, const double *p_yd4, const double *p_yd5, const double *p_yd6, const double *p_yd7,
    const double *p_zd0, const double *p_zd1, const double *p_zd2, const double *p_zd3, const double *p_zd4, const double *p_zd5, const double *p_zd6, const double *p_zd7,
    const double *volo,
    const double *vnew,
    double *delx_zeta,
    double *delv_zeta,
    double *delv_xi,
    double *delx_xi,
    double *delx_eta,
    double *delv_eta
){
    double ax,ay,az;
    double dxv,dyv,dzv;

    double vol = volo[0]*vnew[0];
    double norm = double(1.0) / ( vol + m_ptiny );

    double dxj = double(-0.25)*((p_x0[0]+p_x1[0]+p_x5[0]+p_x4[0]) - (p_x3[0]+p_x2[0]+p_x6[0]+p_x7[0])) ;
    double dyj = double(-0.25)*((p_y0[0]+p_y1[0]+p_y5[0]+p_y4[0]) - (p_y3[0]+p_y2[0]+p_y6[0]+p_y7[0])) ;
    double dzj = double(-0.25)*((p_z0[0]+p_z1[0]+p_z5[0]+p_z4[0]) - (p_z3[0]+p_z2[0]+p_z6[0]+p_z7[0])) ;

    double dxi = double( 0.25)*((p_x1[0]+p_x2[0]+p_x6[0]+p_x5[0]) - (p_x0[0]+p_x3[0]+p_x7[0]+p_x4[0])) ;
    double dyi = double( 0.25)*((p_y1[0]+p_y2[0]+p_y6[0]+p_y5[0]) - (p_y0[0]+p_y3[0]+p_y7[0]+p_y4[0])) ;
    double dzi = double( 0.25)*((p_z1[0]+p_z2[0]+p_z6[0]+p_z5[0]) - (p_z0[0]+p_z3[0]+p_z7[0]+p_z4[0])) ;

    double dxk = double( 0.25)*((p_x4[0]+p_x5[0]+p_x6[0]+p_x7[0]) - (p_x0[0]+p_x1[0]+p_x2[0]+p_x3[0])) ;
    double dyk = double( 0.25)*((p_y4[0]+p_y5[0]+p_y6[0]+p_y7[0]) - (p_y0[0]+p_y1[0]+p_y2[0]+p_y3[0])) ;
    double dzk = double( 0.25)*((p_z4[0]+p_z5[0]+p_z6[0]+p_z7[0]) - (p_z0[0]+p_z1[0]+p_z2[0]+p_z3[0])) ;

    /* find delvk and delxk ( i cross j ) */
    ax = dyi*dzj - dzi*dyj ;
    ay = dzi*dxj - dxi*dzj ;
    az = dxi*dyj - dyi*dxj ;

    delx_zeta[0] = vol / sqrt(ax*ax + ay*ay + az*az + m_ptiny) ;

    ax *= norm ;
    ay *= norm ;
    az *= norm ;

    dxv = double(0.25)*((p_xd4[0]+p_xd5[0]+p_xd6[0]+p_xd7[0]) - (p_xd0[0]+p_xd1[0]+p_xd2[0]+p_xd3[0])) ;
    dyv = double(0.25)*((p_yd4[0]+p_yd5[0]+p_yd6[0]+p_yd7[0]) - (p_yd0[0]+p_yd1[0]+p_yd2[0]+p_yd3[0])) ;
    dzv = double(0.25)*((p_zd4[0]+p_zd5[0]+p_zd6[0]+p_zd7[0]) - (p_zd0[0]+p_zd1[0]+p_zd2[0]+p_zd3[0])) ;

    delv_zeta[0] = ax*dxv + ay*dyv + az*dzv ;

    /* find delxi and delvi ( j cross k ) */
    ax = dyj*dzk - dzj*dyk ;
    ay = dzj*dxk - dxj*dzk ;
    az = dxj*dyk - dyj*dxk ;

    delx_xi[0] = vol / sqrt(ax*ax + ay*ay + az*az + m_ptiny) ;

    ax *= norm ;
    ay *= norm ;
    az *= norm ;

    dxv = double(0.25)*((p_xd1[0]+p_xd2[0]+p_xd6[0]+p_xd5[0]) - (p_xd0[0]+p_xd3[0]+p_xd7[0]+p_xd4[0])) ;
    dyv = double(0.25)*((p_yd1[0]+p_yd2[0]+p_yd6[0]+p_yd5[0]) - (p_yd0[0]+p_yd3[0]+p_yd7[0]+p_yd4[0])) ;
    dzv = double(0.25)*((p_zd1[0]+p_zd2[0]+p_zd6[0]+p_zd5[0]) - (p_zd0[0]+p_zd3[0]+p_zd7[0]+p_zd4[0])) ;

    // domain.delv_xi(i) = ax*dxv + ay*dyv + az*dzv ;
    delv_xi[0] = ax*dxv + ay*dyv + az*dzv ;

    /* find delxj and delvj ( k cross i ) */

    ax = dyk*dzi - dzk*dyi ;
    ay = dzk*dxi - dxk*dzi ;
    az = dxk*dyi - dyk*dxi ;

    // domain.delx_eta(i) = vol / sqrt(ax*ax + ay*ay + az*az + m_ptiny) ;
    delx_eta[0] = vol / sqrt(ax*ax + ay*ay + az*az + m_ptiny) ;

    ax *= norm ;
    ay *= norm ;
    az *= norm ;

    dxv = double(-0.25)*((p_xd0[0]+p_xd1[0]+p_xd5[0]+p_xd4[0]) - (p_xd3[0]+p_xd2[0]+p_xd6[0]+p_xd7[0])) ;
    dyv = double(-0.25)*((p_yd0[0]+p_yd1[0]+p_yd5[0]+p_yd4[0]) - (p_yd3[0]+p_yd2[0]+p_yd6[0]+p_yd7[0])) ;
    dzv = double(-0.25)*((p_zd0[0]+p_zd1[0]+p_zd5[0]+p_zd4[0]) - (p_zd3[0]+p_zd2[0]+p_zd6[0]+p_zd7[0])) ;

    // domain.delv_eta(i) = ax*dxv + ay*dyv + az*dzv ;
    delv_eta[0] = ax*dxv + ay*dyv + az*dzv ;

}