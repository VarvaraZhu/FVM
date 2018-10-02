#pragma once

#include "coord.h"
#include <string>

void calculateGeom(const matrixVec& mesh,\
   matrixScalar& cellVolumes, matrixVec& cellCenter, matrixVec& iFaceCenter, \
   matrixVec& iFaceVector, matrixVec& jFaceCenter, matrixVec& jFaceVector);
