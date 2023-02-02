inline void CalcNewE(double *e_new, const double *e_old, const double *delvc, const double *p_old, const double *q_old, const double *work){
    e_new[0] = e_old[0] - double(0.5) * delvc[0] * (p_old[0] + q_old[0])
        + double(0.5) * work[0];

    if (e_new[0]  < m_emin ) {
        e_new[0] = m_emin ;
    }
}