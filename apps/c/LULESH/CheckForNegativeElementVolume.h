inline void CheckForNegativeElementVolume(const double *determ){
    if(determ[0] <= 0.0){
        // exit(-1); // Volume Error
    }
}