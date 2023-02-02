inline void CalcNewEStep3(double *e_new, const double *work){
    e_new[0] += double(0.5) * work[0];

    if (fabs(e_new[0]) < m_e_cut) {
        e_new[0] = double(0.)  ;
    }
    if (     e_new[0]  < m_emin ) {
        e_new[0] = m_emin ;
    }
}
