#include "point.h"

typedef std::vector<std::vector<point>> matrixVec;
typedef std::vector<std::vector<double>> matrixScalar;

struct Mesh{
    size_t NI, NJ;
    matrixVec coords, cellCentre, cellXFaceCentre, cellXFaceVector,
                                    cellYFaceCentre, cellYFaceVector;
    matrixScalar cellVolume;
    Mesh(char *fileName);
    Mesh();
    ~Mesh();
};
