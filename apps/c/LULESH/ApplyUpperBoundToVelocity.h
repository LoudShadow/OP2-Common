inline void ApplyUpperBoundToVelocity(double *vnewc){
    if (vnewc[0] > (m_eosvmax))
        vnewc[0] = (m_eosvmax) ;
}