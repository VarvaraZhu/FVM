#include "functions.h"

double Pressure(const coord &a){
  return (a.x + a.y);
}

double linearInterpolation(const double& x1, const double& x2, const double& d1, const double& d2){
  return((d1 * x1 + d2 * x2) / (d1 + d2));
}
