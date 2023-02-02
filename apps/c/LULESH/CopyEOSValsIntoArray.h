inline void CopyEOSValsIntoArray(
                                double *e_old, const double *e,
                                double *delvc, const double *delv,
                                double *p_old, const double *p,
                                double *q_old, const double *q,
                                double *qq_old, const double *qq,
                                double *ql_old, const double *ql){
    e_old[0] = e[0] ;
    delvc[0] = delv[0] ;

    p_old[0] = p[0] ;
    q_old[0] = q[0] ;

    qq_old[0] = qq[0] ;
    ql_old[0] = ql[0] ;
}