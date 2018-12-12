#pragma once

#include "coord.h"

double Pressure(const coord &a, const bool & linear);
coord gradientExact(const coord &a, const bool& linear);
double laplacianExact(const coord &a, const bool& linear);
double error(const coord &a, const coord &b);

coord Velocity(const coord &a, const bool& linear);
double divergenceExact(const coord &a, const bool& linear);
double error(const double &a, const double &b);

double divergencePVExact(const coord &a, const bool& linear);

double linearInterpolation(const double& x1, const double& x2, const double& d1, const double& d2);
coord linearInterpolation(const double& x1, const double& x2, const coord& d1, const coord& d2);

double matrixMax(const matrixScalar& A);