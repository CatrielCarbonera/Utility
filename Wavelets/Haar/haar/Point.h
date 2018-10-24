#pragma once
#include <cmath>

struct point_type
{
    point_type() {}
    point_type(double x, double y, double z) : x(x), y(y), z(z) {}
    point_type(double x, double y) : x(x), y(y), z(0) {}
    double x, y, z;
};


inline point_type operator*(double b, point_type const & a) { return point_type(a.x * b, a.y * b, a.z * b); }

inline point_type & operator+=(point_type & a, point_type const & b)
{
    a.x += b.x;
    a.y += b.y;
    a.z += b.z;
    return a;
}

inline point_type & operator*=(point_type & a, double b)
{
    a.x *= b;
    a.y *= b;
    a.z *= b;
    return a;
}

inline point_type & operator-=(point_type & a, point_type const &b)
{
    a.x -= b.x;
    a.y -= b.y;
    a.z -= b.z;
    return a;
}

inline point_type operator+(point_type const & a, point_type const & b) { return point_type(a.x + b.x, a.y + b.y, a.z + b.z); }

inline point_type operator-(point_type const & a, point_type const & b) { return point_type(a.x - b.x, a.y - b.y, a.z - b.z); }

inline double dot(point_type const & a, point_type const & b)
{
    return a.x * b.x + a.y * b.y + a.z * b.z;
}

inline double length_sq(point_type const & pt) 
{
    return dot(pt, pt);
}

inline double length(point_type const & a) { return std::sqrt(length_sq(a)); }
