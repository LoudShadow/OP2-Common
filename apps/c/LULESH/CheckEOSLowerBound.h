inline void CheckEOSLowerBound(const double *vnewc, double *compHalfStep, const double *compression){
    if (vnewc[0] <= m_eosvmin) { /* impossible due to calling func? */
        compHalfStep[0] = compression[0] ;
    }
}