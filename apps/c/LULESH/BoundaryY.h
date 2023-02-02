#define FREE_NODE 0x00
#define BOUNDARY_NODE 0x01
inline void BoundaryY(double *ydd, const int *symmY){
    if(symmY[0] & 0x01){ //BOUNDARY NODE
        ydd[0] = double(0.0);
    }
}
