#define FREE_NODE 0x00
#define BOUNDARY_NODE 0x01
inline void BoundaryZ(double *zdd, const int *symmZ){
    if(symmZ[0] & 0x01){ //BOUNDARY NODE
        zdd[0] = double(0.0);
    }
}