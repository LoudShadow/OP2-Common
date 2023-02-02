inline void CalcNewEStep2(const double *compHalfStep, const double *delvc, double *q_new, const double *pbvc, double *e_new, const double *bvc,
                            const double *pHalfStep, const double *ql_old, const double *qq_old, const double *p_old, const double *q_old
){
    double vhalf = double(1.) / (double(1.) + compHalfStep[0]) ;

    if ( delvc[0] > double(0.) ) {
        q_new[0] /* = qq_old[0] = ql_old[0] */ = double(0.) ;
    }
    else {
        double ssc = ( pbvc[0] * e_new[0]
                + vhalf * vhalf * bvc[0] * pHalfStep[0] ) / m_refdens ;

        if ( ssc <= m_ssc_thresh ) {
        ssc = m_ssc_low ;
        } else {
        ssc = sqrt(ssc) ;
        }

        q_new[0] = (ssc*ql_old[0] + qq_old[0]) ;
    }

    e_new[0] = e_new[0] + double(0.5) * delvc[0]
        * (  double(3.0)*(p_old[0]     + q_old[0])
            - double(4.0)*(pHalfStep[0] + q_new[0])) ;
}