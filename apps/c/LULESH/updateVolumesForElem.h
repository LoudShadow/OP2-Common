inline void updateVolumesForElem(const double *vnew, double *v){
    double tmpV = vnew[0] ;

    if ( fabs(tmpV - double(1.0)) < m_v_cut )
    tmpV = double(1.0) ;

    // domain.v(i) = tmpV ;
    v[0] = tmpV ;
}