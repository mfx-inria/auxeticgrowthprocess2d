#ifndef VEC2_H
#define VEC2_H

template<typename T>
class Vec2 {
 public:
    Vec2(T xx, T yy) : x(xx), y(yy) {}
    Vec2 operator + (const Vec2 &v) const { return Vec2(x + v.x, y + v.y); }
    Vec2 operator - (const Vec2 &v) const { return Vec2(x - v.x, y - v.y); }
    Vec2 operator * (const T &r) const { return Vec2(x * r, y * r); }
    Vec2 operator / (const T &r) const { return Vec2(x / r, y / r); }
    T dotProduct(const Vec2<T> &v) const { return x * v.x + y * v.y; }
    Vec2& operator /= (const T &r) { x /= r, y /= r; return *this; }
    T length() const { return sqrt(x * x + y * y); }
    T length_squared() const { return x * x + y * y; }
    T& operator[] (int i) { return (&x)[i]; }
    friend std::ostream& operator << (std::ostream &s, const Vec2<T> &v) { return s << '[' << v.x << ' ' << v.y << ']';}
    T x, y;
};
typedef Vec2<double> Vec2d;

Vec2d get_coordinate(unsigned int x, unsigned int y, unsigned int image_size) {
    Vec2d coordinate(static_cast<double>(x) + 0.5, static_cast<double>(y) + 0.5);
    coordinate /= static_cast<double>(image_size);
    return coordinate;
}

#endif
