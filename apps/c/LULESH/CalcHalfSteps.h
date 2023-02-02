inline void CalcHalfSteps(double *compression, const double *vnewc, const double* delvc, double *compHalfStep){
    double vchalf ;
    compression[0] = double(1.) / vnewc[0] - double(1.);
    vchalf = vnewc[0] - delvc[0] * double(.5);
    compHalfStep[0] = double(1.) / vchalf - double(1.);
}