inline void CalcPNew(double *p_new, const double *bvc, const double *e_old, const double *vnewc){
    p_new[0] = bvc[0] * e_old[0] ;

    if    (fabs(p_new[0]) <  m_p_cut   )
        p_new[0] = double(0.0) ;

    if    ( vnewc[0] >= m_eosvmax ) /* impossible condition here? */
        p_new[0] = double(0.0) ;

    if    (p_new[0]       <  m_pmin)
        p_new[0]   = m_pmin ;
}