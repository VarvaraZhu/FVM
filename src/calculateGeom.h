#pragma once

#include "coord.h"
#include <string>


void calculateMetric(const matrixVec& mesh,\
                      matrixScalar& cellVolumes, matrixVec& cellCenter, matrixVec& iFaceCenter, \
                        matrixVec& iFaceVector, matrixVec& jFaceCenter, matrixVec& jFaceVector);
