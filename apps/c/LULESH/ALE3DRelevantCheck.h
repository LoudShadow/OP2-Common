inline void ALE3DRelevantCheck(const double *v){
    double vc = v[0] ;

    if (m_eosvmin != double(0.)) {
        if (vc < m_eosvmin)
        vc = m_eosvmin ;
    }
    if (m_eosvmax != double(0.)) {
        if (vc > m_eosvmax)
        vc = m_eosvmax ;
    }
    if (vc <= 0.) {
        // exit(-1); // Volume Error
    }
}
