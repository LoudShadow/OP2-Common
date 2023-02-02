inline void CalcHydroConstraint(const double *vdov, double *dthydro){
    double dtdvov = m_dvovmax / (fabs(vdov[0])+double(1.e-20)) ;

    if ( (*dthydro) > dtdvov ) {
            (*dthydro) = dtdvov ;
            // hydro_elem = indx ;
    }
}