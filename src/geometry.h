#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

#include "typedefs.h"
#include <cmath>

class Point{
public: 
  inline Point(double _x, double _y, double _z)
    : x(_x), y(_y), z(_z)
  {}

  inline Point(const Point& p)
    : x(p.x), y(p.y), z(p.z)
  {}

  double x, y, z;
};

class Ray {
public:
  inline Ray(const Point& start, double _vx, double _vy, double _vz)
    : x(start.x), y(start.y), z(start.z), vx(_vx), vy(_vy), vz(_vz)
  {}

  cdouble x, y, z, vx, vy, vz;
};
	
class Plane {
public:
 Plane(const Ray& normal_vect)
   : normal_x(normal_vect.x), normal_y(normal_vect.y), normal_z(normal_vect.z),
    normal_vx(normal_vect.vx), normal_vy(normal_vect.vy), normal_vz(normal_vect.vz)
    {}
  
  inline bool onFrontSide(const Point& p) const {
    return  (p.x - normal_x)*normal_vx + (p.y - normal_y)*normal_vy + (p.z - normal_z)*normal_vz >= 0;
  }

  cdouble normal_x, normal_y, normal_z;
  cdouble normal_vx, normal_vy, normal_vz;
};

class Box {
public:
  inline Box(double _xl, double _xu, double _yl, double _yu, double _zl, double _zu)
    : xl(_xl),xu(_xu),yl(_yl),yu(_yu),zl(_zl),zu(_zu)
  {}

  inline bool inBox(double x, double y, double z) const {
    return (x <= xu && x >= xl && y >= yl && y <= yu && z >= zl && z <= zu);
  }
    
  inline bool intersects(const Box& b) const {
    return ! ( xl >= b.xu || xu <= b.xl
	       || yl >= b.yu || yu <= b.yl 
	       || zl >= b.zu || zu <= b.zl);
  }
  
  inline Point center() const {
    return Point(0.5 * (xl + xu), 0.5 * (yl + yu), 0.5 * (zl + zu));
  }

  inline Point centeredCorner() const {
    return Point(0.5 * (xu - xl), 0.5 * (yu - yl), 0.5 * (zu - zl));
  }

  cdouble xl, xu, yl, yu, zl, zu;
};

inline double dist2(const Point& p1, const Point& p2) {
  
  double xd = p1.x - p2.x;
  double yd = p1.y - p2.y;
  double zd = p1.z - p2.z;

  return xd*xd + yd*yd + zd*zd;
}

inline double dist(const Point& p1, const Point& p2) {
  return sqrt(dist2(p1,p2));
}

#endif /* _GEOMETRY_H_ */
