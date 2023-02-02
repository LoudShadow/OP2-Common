inline void CalcCourantConstraint(const double *ss, const double *vdov, const double *arealg, double *dtcourant){

    double dtf = ss[0] * ss[0];

    if ( vdov[0] < double(0.) ) {

        dtf = dtf
            + m_qqc2 * arealg[0] * arealg[0]
            * vdov[0] * vdov[0] ;
    }

    dtf = sqrt(dtf) ;

    dtf = arealg[0] / dtf ;

/* determine minimum timestep with its corresponding elem */
    if (vdov[0] != double(0.)) {
        if ( dtf < (*dtcourant) ) {
        (*dtcourant) = dtf ;
        // courant_elem = indx ;
        }
    }
}