inline void ApplyLowerBoundToVelocity(double *vnewc){
    if (vnewc[0] < (m_eosvmin))
        vnewc[0] = (m_eosvmin) ;
}