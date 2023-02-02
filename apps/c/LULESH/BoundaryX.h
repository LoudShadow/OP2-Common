#define FREE_NODE 0x00
#define BOUNDARY_NODE 0x01
inline void BoundaryX(double *xdd, const int *symmX){
    if(symmX[0] & 0x01){ //BOUNDARY NODE
        xdd[0] = double(0.0);
    }
}