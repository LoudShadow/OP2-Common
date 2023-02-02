inline void CalcPosForNodes(
                            double *x, double *y, double *z,
                            const double *xd, const double *yd, const double *zd,
                            const double *dt
){
    x[0] += xd[0] * (*dt) ;
    y[0] += yd[0] * (*dt) ;
    z[0] += zd[0] * (*dt) ;
}