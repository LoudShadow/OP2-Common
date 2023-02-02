inline void initStressTerms(double *sigxx,double *sigyy,double *sigzz, const double *p, const double *q)
{
    sigxx[0] = -p[0] - q[0];
    sigyy[0] = -p[0] - q[0];
    sigzz[0] = -p[0] - q[0];
}