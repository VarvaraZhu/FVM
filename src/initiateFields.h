#pragma once

#include "coord.h"

void InitiatePressureField(matrixScalar &P, const matrixVec& cellCenter, const bool& linear);
void InitiateVelocityField(matrixVec &V, const matrixVec& cellCenter, const bool& linear);
void InitiatePVField(matrixVec &PV, const matrixScalar &P, const matrixVec& V);

void CalculateExactGradientField(matrixVec &gradPExact, const matrixVec& cellCenter, const bool& linear);
void CalculateExactLaplacianField(matrixScalar &lapPExact, const matrixVec& cellCenter, const bool& linear);
void CalculateExactDivergenceField(matrixScalar &divVExact, const matrixVec& cellCenter, const bool& linear);
void CalculateExactPVDivergenceField(matrixScalar &divPVExact, const matrixVec& cellCenter, const bool& linear);

void CalculateErrorField(const matrixVec &gradP, const matrixVec &gradPExact, matrixScalar &Error);
void CalculateErrorField(const matrixScalar &divV, const matrixScalar &divVExact, matrixScalar &Error);