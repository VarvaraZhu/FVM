struct point{
  double x, y;
  point(double x, double y): x(x), y(y) {};
  point(point &other): x(other.x), y(other.y){};
  ~point(){};
};
