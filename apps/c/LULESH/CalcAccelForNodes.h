inline void CalcAccelForNodes(
                            double *xdd, double *ydd, double *zdd,
                            const double *fx, const double *fy, const double *fz,
                            const double *nodalMass
                            ){
    xdd[0] = fx[0] / nodalMass[0];
    ydd[0] = fy[0] / nodalMass[0];
    zdd[0] = fz[0] / nodalMass[0];
}