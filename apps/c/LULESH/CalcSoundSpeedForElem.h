inline void CalcSoundSpeedForElem(const double *pbvc, const double *enewc, const double *vnewc, const double *bvc, const double *pnewc, double *ss){
    double ssTmp = (pbvc[0] * enewc[0] + vnewc[0] * vnewc[0] *
                bvc[0] * pnewc[0]) / m_refdens;
    if (ssTmp <= m_ssc_thresh) {
        ssTmp = m_ssc_low;
    }
    else {
        ssTmp = sqrt(ssTmp);
    }
    ss[0] = ssTmp ;
}