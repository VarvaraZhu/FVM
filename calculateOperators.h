#include "coord.h"
#include <vector>

void calculateGradient(const matrixScalar& scalarField, matrixVec& gradient,\
    const matrixScalar& cellVolumes, const matrixVec& cellCenter,\
        const matrixVec& iFaceCenter, const matrixVec& iFaceVector,\
            const matrixVec& jFaceCenter, const matrixVec& jFaceVector);
