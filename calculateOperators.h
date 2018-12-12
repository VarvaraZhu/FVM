#pragma once

#include "coord.h"

#include <vector>
#include <string>

void calculateGradientGauss(const matrixScalar& Pressure, matrixVec& gradient, \
                              const matrixScalar& cellVolumes, const matrixVec& cellCenter, \
                                const matrixVec& iFaceCenter, const matrixVec& iFaceVector, \
                                  const matrixVec& jFaceCenter, const matrixVec& jFaceVector);

void calculateGradientIterGauss(const matrixScalar& Pressure, matrixVec& gradient, \
                                  const matrixScalar& cellVolumes, const matrixVec& cellCenter, \
                                    const matrixVec& iFaceCenter, const matrixVec& iFaceVector, \
                                      const matrixVec& jFaceCenter, const matrixVec& jFaceVector);

void calculateGradientOLS(const matrixScalar& Pressure, matrixVec& gradient, \
                            const matrixScalar& cellVolumes, const matrixVec& cellCenter, \
                              const matrixVec& iFaceCenter, const matrixVec& iFaceVector, \
                                const matrixVec& jFaceCenter, const matrixVec& jFaceVector);

void calculateDivergence(const matrixVec& Velocity, matrixScalar& divergence, \
                          const matrixScalar& cellVolumes, const matrixVec& cellCenter, \
                            const matrixVec& iFaceCenter, const matrixVec& iFaceVector, \
                              const matrixVec& jFaceCenter, const matrixVec& jFaceVector);

void calculatePVDivergence(const std::string& scheme, const matrixScalar& Pressure, const matrixVec& Velocity, matrixScalar& divergence, \
                            const matrixScalar& cellVolumes, const matrixVec& cellCenter, \
                              const matrixVec& iFaceCenter, const matrixVec& iFaceVector, \
                                const matrixVec& jFaceCenter, const matrixVec& jFaceVector, \
                                  const matrixVec& gradP = matrixVec());

void calculateLaplacian(const matrixScalar& Pressure, matrixScalar& laplacian,
 const matrixVec& gradP, \
                          const matrixScalar& cellVolumes,
 const matrixVec& cellCenter, \
                            const matrixVec& iFaceCenter,
 const matrixVec& iFaceVector, \
                              const matrixVec& jFaceCenter, const matrixVec& jFaceVector, bool &correction, size_t order);