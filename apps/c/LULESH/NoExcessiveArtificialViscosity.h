inline void NoExcessiveArtificialViscosity(const double *q){
    if ( q[0] > m_qstop ) {
        // exit(-2); // Q Stop Error
    }
}