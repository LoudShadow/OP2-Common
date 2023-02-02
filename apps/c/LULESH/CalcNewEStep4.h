inline void CalcNewEStep4(const double *delvc, const double *pbvc, double *e_new, const double *vnewc, const double *bvc, const double *p_new,
                            const double *ql_old, const double *qq_old, const double *p_old, const double *q_old, const double *q_new,
                            const double *pHalfStep){
    // int ielem = regElemList[0];
    double q_tilde ;

    if (delvc[0] > double(0.)) {
        q_tilde = double(0.) ;
    }
    else {
        double ssc = ( pbvc[0] * e_new[0]
                + vnewc[0] * vnewc[0] * bvc[0] * p_new[0] ) / m_refdens ;

        if ( ssc <= m_ssc_thresh ) {
        ssc = m_ssc_low ;
        } else {
        ssc = sqrt(ssc) ;
        }

        q_tilde = (ssc*ql_old[0] + qq_old[0]) ;
    }

    e_new[0] = e_new[0] - (  double(7.0)*(p_old[0]     + q_old[0])
                            - double(8.0)*(pHalfStep[0] + q_new[0])
                            + (p_new[0] + q_tilde)) * delvc[0]*m_sixth ;

    if (fabs(e_new[0]) < m_e_cut) {
        e_new[0] = double(0.)  ;
    }
    if (     e_new[0]  < m_emin ) {
        e_new[0] = m_emin ;
    }
}