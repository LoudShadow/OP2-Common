inline void CalcLagrangeElemRemaining(
                                    double *dxx, double *dyy, double *dzz,
                                    double *m_vdov,
                                    const double *vnew

){
    // calc strain rate and apply as constraint (only done in FB element)
    double vdov = dxx[0] + dyy[0] + dzz[0] ;
    double vdovthird = vdov/double(3.0) ;

    // make the rate of deformation tensor deviatoric
    m_vdov[0] = vdov ;

    dxx[0] -= vdovthird ;
    dyy[0] -= vdovthird ;
    dzz[0] -= vdovthird ;

// See if any volumes are negative, and take appropriate action.
    if (vnew[0] <= double(0.0))
    {
        // exit(-1); // Volume Error
    }
}