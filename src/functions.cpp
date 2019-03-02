#include "functions.h"

#include <cmath>
#include <cstddef>

//Fuction defining Pressure field
double Pressure(const coord &a, const bool & linear){
  //return (a.x + a.y + 2 * a.x * a.x + a.y * a.y);
  if (linear == 1)
    return (a.x + a.y);
  else
    return sin(a.x + a.y);
}

//Fuction defining Velocity field
coord Velocity(const coord &a, const bool& linear) {
  if (linear == 1)
    return coord(a.x + 1, a.y + 1);
  else
    return coord(1 + cos(a.x), 1 + sin(a.y));
}

//Fuction defining exact Pressure gradient field
coord gradientExact(const coord &a, const bool& linear) {
  if (linear == 1)
    return coord(1, 1);
  else
    return coord(cos(a.x + a.y), cos(a.x + a.y));
}

//Fuction defining exact Pressure laplacian field
double laplacianExact(const coord &a, const bool& linear) {
  //return 6.0;
  if (linear == 1)
    return 0.;
  else
    return (-2 *sin(a.x + a.y));
}

//Fuction defining exact Velocity divergence field
double divergenceExact(const coord &a, const bool& linear) {
  if (linear == 1)
    return 2;
  else
    return (cos(a.y) - sin(a.x));
}

//Fuction defining exact Velocity x Pressure divergence field
double divergencePVExact(const coord &a, const bool& linear) {
  if (linear == 1)
    return (3 * a.x + 3 * a.y + 2.0);
  else
    return (sin(a.x + 2 * a.y) + 2 * cos(a.x + a.y) + cos(2 * a.x + a.y));
}

//Error function for scalar variable
double error(const coord &a, const coord &b) {
  return norm2(a - b) / norm2(b);
}

//Error function for vector variable
double error(const double &a, const double &b) {
  return std::abs(a - b) / std::abs(b);
}

//linear interpolation of scalar variable
double linearInterpolation(const double& x1, const double& x2, const double& d1, const double& d2){
  return(d1 * (x2 / (x1 + x2)) + d2 * (x1 / (x1 + x2)));
}

//linear interpolation of vector variable
coord linearInterpolation(const double& x1, const double& x2, const coord& d1, const coord& d2) {
  return(d1 * (x2 / (x1 + x2)) + d2 * (x1 / (x1 + x2)));
}

double max(const double& a, const double& b) {
  if (a > b) return a;
  else return b;
}

double matrixMax(const matrixScalar& A) {
  double res = A[1][1];

  for (size_t i = 1; i != A.size() - 1; i++)
    for (size_t j = 1; j != A[0].size() - 1; j++)
      res = max(res, A[i][j]);

  return res;
}
