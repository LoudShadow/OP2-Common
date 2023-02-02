inline void CopyTempEOSVarsBack(
                                double *p, const double *p_new,
                                double *e, const double *e_new,
                                double *q, const double *q_new){
    p[0] = p_new[0];
    e[0] = e_new[0];
    q[0] = q_new[0];
}