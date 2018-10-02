#pragma once

#include <vector>
#include <cmath>

struct coord{
  double x, y;

  coord(): x(0), y(0){}
  coord(const double x, const double y): x(x), y(y){}
  coord(const coord &a): x(a.x), y(a.y){}

  coord operator=(const coord &a){
    if(&a == this)
      return *this;
    this->x = a.x;
    this->y = a.y;
    return *this;
  }

  double operator*(const coord &a) const{
    return (this->x * a.x + this->y * a.y);

  }
  coord operator-(const coord &a) const{
    return coord(this->x - a.x, this->y - a.y);
  }
  coord operator+(const coord &a) const{
    return coord(this->x + a.x, this->y + a.y);
  }
  coord operator*(const double &d) const{
    return coord(d * this->x, d * this->y);
  }

};

inline double norm2(const coord& a){
  return pow((a * a), 0.5);
}

using matrixScalar = std::vector<std::vector<double>>;
using matrixVec = std::vector<std::vector<coord>>;
