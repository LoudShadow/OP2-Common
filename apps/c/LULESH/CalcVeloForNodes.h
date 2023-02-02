inline void CalcVeloForNodes(
                            double *xd, double *yd, double *zd,
                            const double *xdd, const double *ydd, const double *zdd,
                            const double *dt
){
    //Review need for OP_inc
    double xdtmp, ydtmp, zdtmp ;

    xdtmp = xd[0] + xdd[0] * (*dt) ;
    if( fabs(xdtmp) < m_u_cut ) xdtmp = double(0.0);
    xd[0] = xdtmp ;

    ydtmp = yd[0] + ydd[0] * (*dt) ;
    if( fabs(ydtmp) < m_u_cut ) ydtmp = double(0.0);
    yd[0] = ydtmp ;

    zdtmp = zd[0] + zdd[0] * (*dt) ;
    if( fabs(zdtmp) < m_u_cut ) zdtmp = double(0.0);

    zd[0] = zdtmp ;
}