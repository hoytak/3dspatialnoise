#ifndef _GEOMETRY_H_
#define _GEOMETRY_H_

class Point{
public: 
    inline Point(double x, double y, double z);
    double x, y, z;
};

class Ray {
public:
    inline Ray(const Point& start, double vx, double vy, double vz);
    
    double x, y, z, vx, vy, vz;
};
	
class Plane {
public:
    inline Plane(const Ray& normal_vect);
    
    inline bool onFrontSide(const Point& p) const;

    double normal_x, normal_y, normal_z;
    double normal_vx, normal_vy, normal_vz;
};

class Box {
public:
    inline Box(double xl, double xu, double yl, double yu, double zl, double zu);

    inline bool inBox(double x, double y, double z) const;
    
    inline bool intersects(const Box& b) const;

    double xl, xu, yl, yu, zl, zu;
};

/**************************************************/
/*  Inline functions of the above.                */

inline bool Plane::onFrontSide(const Point& p) const
{
    return  (p.x - normal_x)*normal_vx + (p.y - normal_y)*normal_vy + (p.z - normal_z)*normal_vz >= 0;
}

inline Point::Point(double _x, double _y, double _z)
    : x(_x), y(_y), z(_z)
{
}

inline Ray::Ray(const Point& start, double _vx, double _vy, double _vz)
    : x(start.x), y(start.y), z(start.z), vx(_vx), vy(_vy), vz(_vz)
{
}

inline Plane::Plane(const Ray& normal_vect)
    : normal_x(normal_vect.x), normal_y(normal_vect.y), normal_z(normal_vect.z),
      normal_vx(normal_vect.vx), normal_vy(normal_vect.vy), normal_vz(normal_vect.vz)
{
}

inline Box::Box(double _xl, double _xu, double _yl, double _yu, double _zl, double _zu)
	   : xl(_xl),xu(_xu),yl(_yl),yu(_yu),zl(_zl),zu(_zu)
{
}

inline bool Box::inBox(double x, double y, double z) const
{
    return (x <= xu && x >= xl && y >= yl && y <= yu && z >= zl && z <= zu);
}

inline bool Box::intersects(const Box& b) const
{
    return ! ( xl >= b.xu || xu <= b.xl
	       || yl >= b.yu || yu <= b.yl 
	       || zl >= b.zu || zu <= b.zl);
}

#endif /* _GEOMETRY_H_ */
