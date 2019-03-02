#pragma once
#include "coord.h"

void ouputToTecplot(std::string outputFileName, const matrixVec& mesh, matrixScalar& Pressure, const matrixVec& V, const matrixVec& pV, \
  const matrixVec& gradP, const matrixScalar& divV, const matrixScalar& divPV, const matrixScalar& lapP, \
  const matrixVec& gradPExact, const matrixScalar& divVExact, const matrixScalar& divPVExact, const matrixScalar& lapPExact, \
  const matrixScalar& gradPError, const matrixScalar& divVError, const matrixScalar& divPVError, const matrixScalar& lapPError);