inline void CalcHalfStepBVC(double *bvc, const double *compression, double *pbvc){
    bvc[0] = m_c1s * (compression[0] + double(1.));
    pbvc[0] = m_c1s;
}
