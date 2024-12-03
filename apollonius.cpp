/**
 * Copyright (c) 2023-present Merlot.Rain
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#include "apollonius.h"
#include <assert.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <set>
#include <array>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define APO_TOLERANCE (1.0e-9)

#define APO_UNKOWN_TYPE (-1)
#define APO_POINT_TYPE (0)
#define APO_LINE_TYPE (1)
#define APO_CIRCLE_TYPE (2)

#define APO_IS_POINT(shp) ((shp).type == APO_POINT_TYPE)
#define APO_IS_LINE(shp) ((shp).type == APO_LINE_TYPE)
#define APO_IS_CIRCLE(shp) ((shp).type == APO_CIRCLE_TYPE)
#define APO_IS_NULL(shp) ((shp).type == APO_UNKOWN_TYPE)

#define APO_IS_NORMALE_NUMBER(x) (!std::isnan(x) && !std::isinf(x))

enum Side
{
  Left,
  Right,
  BothSides,
};

struct APO_Point
{
  double x;
  double y;
  bool valid;

  APO_Point();
  APO_Point(double x, double y, bool valid_in = true);
  void set(double vx, double vy);
  void setPolar(double radius, double angle);
  bool isValid() const;

  APO_Point operator+(const APO_Point& rhs) const;
  APO_Point operator-(const APO_Point& rhs) const;
  APO_Point operator*(double scale) const;
  APO_Point operator/(double scale) const;
  APO_Point& operator+=(const APO_Point& rhs);
  APO_Point& operator-=(const APO_Point& rhs);
  APO_Point& operator*=(double scale);
  APO_Point& operator/=(double scale);

  double getDistanceTo(const APO_Point& point, bool limited = true) const;
  bool equalsFuzzy(const APO_Point& point) const;
  void setAngle(double a);
  double getAngle() const;
  double getAngleTo(const APO_Point& point) const;
  double getMagnitude() const;
  void setMagnitude(double m);

  static APO_Point getAverage(const APO_Point& v1, const APO_Point& v2);
  static double getCrossProduct(const APO_Point& v1, const APO_Point& v2);
  static double getDotProduct(const APO_Point& v1, const APO_Point& v2);
  static APO_Point createPolar(double radius, double angle);
  static std::vector<APO_Point>
  getUnique(const std::vector<APO_Point>& vectors);
  APO_Point& rotate(double rotation);
  APO_Point& rotate(double rotation, const APO_Point& center);

  static APO_Point invalid;
  static APO_Point undefined;

  friend bool
  operator==(const APO_Point& lhs, const APO_Point& rhs)
  {
    if (lhs.valid == true && rhs.valid == true)
    {
      return lhs.x == rhs.x && lhs.y == rhs.y;
    }
    else if (lhs.valid == false && rhs.valid == false)
    {
      return true;
    }
    return false;
 }

  friend bool operator!=(const APO_Point& lhs, const APO_Point& rhs)
  {
    return !(lhs == rhs);
  }
};

struct APO_Line
{
  APO_Point begin_point;
  APO_Point end_point;

  APO_Line();
  APO_Line(const APO_Point& begin, const APO_Point& end);
  APO_Line(double x1, double y1, double x2, double y2);
  APO_Line(const APO_Point& begin, double angle, double distance);

  double getAngle() const;
  void setLength(double l, bool fromStart = true);
  double getLength() const;
  APO_Point getClosestPoint(const APO_Point& p, bool limited) const;
  APO_Point getMiddlePoint() const;
  bool rotate(double rotation, const APO_Point& center);
  bool move(const APO_Point& offset);
  double getDistanceTo(const APO_Point& point, bool limited = true) const;
  void reverse();
  bool moveTo(const APO_Point& dest);
  std::vector<APO_Line> getOffset(double distance, double num,
                                  Side side) const;
  APO_Point getVectorTo(const APO_Point& point, bool limited) const;
  Side getSideOfPoint(const APO_Point& point) const;

  static APO_Line undefined;
};

struct APO_Circle
{
  APO_Point center;
  double radius;

  APO_Circle();
  APO_Circle(double center_x, double center_y, double radius);
  APO_Circle(const APO_Point& center, double radius);

  static APO_Circle createFrom3Points(const APO_Point& p1, const APO_Point& p2,
                                      const APO_Point& p3);
  static APO_Circle createFrom2Points(const APO_Point& p1,
                                      const APO_Point& p2);

  bool containsPoint(const APO_Point& p) const;
};

struct APO_Shape
{
  int type;
  union
  {
    APO_Point point;
    APO_Line line;
    APO_Circle circle;
  };

  APO_Shape();
  APO_Shape(const APO_Point& p);
  APO_Shape(const APO_Line& l);
  APO_Shape(const APO_Circle& c);
};

template<class T>
bool isUndefined(const T& t)
{
  return false;
}

template<>
bool
isUndefined<APO_Point>(const APO_Point& t)
{
  return std::isnan(t.x) && std::isnan(t.y) && !t.valid;
}

template<>
bool 
isUndefined<APO_Line>(const APO_Line& t)
{
  return isUndefined<APO_Point>(t.begin_point) && isUndefined<APO_Point>(t.end_point);
}

template<>
bool isUndefined<APO_Circle>(const APO_Circle& t)
{
  return isUndefined<APO_Point>(t.center) && std::isnan(t.radius);
}

template<class T>
T undefined()
{
  return T::undefined;
}

/* ------------------------- static solve functions ------------------------- */

static bool getInverseShape(const APO_Shape& shp,
                            const APO_Shape& inversionCircle,
                            APO_Shape& inversed);

template <class T>
static std::vector<APO_Shape>
getInverseShapes(const std::vector<T>& shapes,
                 const APO_Shape& inversionCircle);

static std::vector<APO_Line>
getCircleTangentsThroughPoint(const APO_Circle& circle, const APO_Point& p);

static std::vector<APO_Line> getAllTangents(const APO_Shape& shp1,
                                            const APO_Shape& shp2);

static std::vector<APO_Circle>
getSolutions(const std::vector<APO_Shape>& shps);

static std::vector<APO_Circle> getSolution(const APO_Shape& shp1,
                                           const APO_Shape& shp2,
                                           const APO_Shape& shp3);

static std::vector<APO_Circle> getSolutionFromPPP(const APO_Point& point1,
                                                  const APO_Point& point2,
                                                  const APO_Point& point3);

static std::vector<APO_Circle> getSolutionFromPPC(const APO_Point& point1,
                                                  const APO_Point& point2,
                                                  const APO_Circle& circle);

static std::vector<APO_Circle> getSolutionFromPPL(const APO_Point& point1,
                                                  const APO_Point& point2,
                                                  const APO_Line& line);

static std::vector<APO_Circle> getSolutionFromPCC(const APO_Point& point,
                                                  const APO_Circle& circle1,
                                                  const APO_Circle& circle2);

static std::vector<APO_Circle> getSolutionFromPLL(const APO_Point& point,
                                                  const APO_Line& line1,
                                                  const APO_Line& line2);

static std::vector<APO_Circle> getSolutionFromPLC(const APO_Point& point,
                                                  const APO_Line& line,
                                                  const APO_Circle& circle);

static std::vector<APO_Circle> getSolutionFromLLL(const APO_Line& line1,
                                                  const APO_Line& line2,
                                                  const APO_Line& line3);

static std::vector<APO_Circle> getSolutionFromLLC(const APO_Line& line1,
                                                  const APO_Line& line2,
                                                  const APO_Circle& circle);

static std::vector<APO_Circle> getSolutionFromLCC(const APO_Line& line,
                                                  const APO_Circle& circle1,
                                                  const APO_Circle& circle2);

static std::vector<APO_Circle> getSolutionFromCCC(const APO_Circle& circle1,
                                                  const APO_Circle& circle2,
                                                  const APO_Circle& circle3);

template<class T, class E>
static std::vector<APO_Point> getIntersectionPoints(const T& t, const E& e, bool limited = true);
template <>
static std::vector<APO_Point> getIntersectionPoints<APO_Circle, APO_Circle>(
    const APO_Circle& t, const APO_Circle& e, bool limited);
template <>
static std::vector<APO_Point>
getIntersectionPoints<APO_Line, APO_Line>(const APO_Line& t, const APO_Line& e,
                                          bool limited);
template <>
static std::vector<APO_Point>
getIntersectionPoints<APO_Line, APO_Circle>(const APO_Line& t,
                                            const APO_Circle& e, bool limited);                   


template<class T, class E>
static bool isIntersectWith(const T& t, const E& e, bool limited);

static bool pointIsOnLine(const APO_Point& point, const APO_Line& line,
                          bool limited);
static bool pointIsOnCircle(const APO_Point& point, const APO_Circle& circle);

static bool fuzzyCompare(double a, double b);
static double getAngleDifference(double a1, double a2);
static std::vector<APO_Line> getAngleBisectors(const APO_Line& line1,
                                               const APO_Line& line2);
static APO_Point getCommonIntersectionPoint(const APO_Circle& c1,
                                            const APO_Circle& c2,
                                            const APO_Circle& c3);
static APO_Point getPowerCenter(const APO_Circle& c1,
                                            const APO_Circle& c2,
                                            const APO_Circle& c3);
static std::vector<APO_Line> getSimilarityAxes(const APO_Circle& c1,
                                               const APO_Circle& c2,
                                               const APO_Circle& c3);
static APO_Point getPole(const APO_Circle& circle, const APO_Line& polarLine);                                             
static std::vector<APO_Circle> getSolutionsCCCAlt(const APO_Circle& circle1,
                                                  const APO_Circle& circle2,
                                                  const APO_Circle& circle3);
template<class T>
static std::vector<T>
verify(const std::vector<T>& candidates, const APO_Shape& shape1,
       const APO_Shape& shape2, const APO_Shape& shape3);

static std::vector<APO_Circle>
getCircles2TR(const APO_Line shp1, const APO_Line& shp2, double radius,
              const APO_Point& pos = APO_Point());

template<class T>
static std::vector<T>
removeDuplicates(const std::vector<T>& shaps);


/* ------------------------- template function impls ------------------------ */

template <class T>
std::vector<APO_Shape>
getInverseShapes(const std::vector<T>& shapes,
                 const APO_Shape& inversionCircle)
{
  std::vector<APO_Shape> res;
  for (auto&& shp : shapes)
  {
    APO_Shape inversed;
    if (getInverseShape(shp, inversionCircle, inversed))
    {
      res.push_back(inversed);
    }
  }
  return res;
}

template <class T, class E>
std::vector<APO_Point>
getIntersectionPoints(const T& t, const E& e, bool limited)
{
  return std::vector<APO_Point>();
}

template <>
std::vector<APO_Point> getIntersectionPoints<APO_Circle, APO_Circle>(
    const APO_Circle& circle1, const APO_Circle& circle2, bool limited)
{

    std::vector<APO_Point> res;

    double r1 = circle1.radius;
    double r2 = circle2.radius;
    if (r1 < r2)
    {
      return getIntersectionPoints<APO_Circle, APO_Circle>(circle2, circle1);
    }

    APO_Point c1 = circle1.center;
    APO_Point c2 = circle2.center;

    APO_Point u  = c2 - c1;
    double u_mag = u.getMagnitude();

    // concentric
    if (u_mag < APO_TOLERANCE)
    {
      return res;
    }

    double tol = (r1 + r2) / 200000;

    // the two circles (almost) touch externally / internally in one point
    // (tangent):
    if (fabs(u_mag - (r1 + r2)) < tol || fabs(u_mag - fabs(r1 - r2)) < tol)
    {
      u.setMagnitude(r1);
      res.push_back(c1 + u);
      return res;
    }

    APO_Point v(u.y, -u.x);

    double s, t1, t2, term;
    s = 1.0 / 2.0 * ((r1 * r1 - r2 * r2) / (pow(u_mag, 2.0)) + 1.0);

    term = (r1 * r1) / (pow(u_mag, 2.0)) - s * s;

    // no intersection:
    if (term < 0.0)
    {
      return res;
    }

    // one or two intersections:
    t1 = sqrt(term);
    t2 = -sqrt(term);

    APO_Point sol1 = c1 + u * s + v * t1;
    APO_Point sol2 = c1 + u * s + v * t2;

    if (fabs(sol1.x - sol2.x) < tol && fabs(sol1.y - sol2.y) < tol)
    {
      res.push_back(sol1);
    }
    else
    {
      res.push_back(sol1);
      res.push_back(sol2);
    }
    return res;
}

template <>
std::vector<APO_Point>
getIntersectionPoints<APO_Line, APO_Line>(const APO_Line& line1, const APO_Line& line2,
                                          bool limited)
{
    std::vector<APO_Point> res;
    double a1 = line1.end_point.y - line1.begin_point.y;
    double b1 = line1.begin_point.x - line1.end_point.x;
    double c1 = a1 * line1.begin_point.x + b1 * line1.begin_point.y;

    double a2 = line2.end_point.y - line2.begin_point.y;
    double b2 = line2.begin_point.x - line2.end_point.x;
    double c2 = a2 * line2.begin_point.x + b2 * line2.begin_point.y;

    double det = a1 * b2 - a2 * b1;
    if (fabs(det) < 1.0e-6)
    {
      return res;
    }
    else
    {
      APO_Point v((b2 * c1 - b1 * c2) / det, (a1 * c2 - a2 * c1) / det);

      if ((!limited || 0 == pointIsOnLine(v, line1, limited))
          && (!limited || 0 == pointIsOnLine(v, line2, limited)))
      {
        res.push_back(v);
      }
    }
    return res;
}

template <>
std::vector<APO_Point>
getIntersectionPoints<APO_Line, APO_Circle>(const APO_Line& line,
                                            const APO_Circle& circle, bool limited)
{

    std::vector<APO_Point> res;

    APO_Point vLineCenter = line.getVectorTo(circle.center, false);
    double dist           = vLineCenter.getMagnitude();

    // special case: arc almost touches line (tangent with tiny gap or tiny
    // overlap):
    if (fuzzyCompare(dist, circle.radius))
    {
      APO_Point sol = circle.center - vLineCenter;
      if (!limited || pointIsOnLine(sol, line, true))
      {
        res.push_back(sol);
      }
      return res;
    }

    APO_Point p = line.begin_point;
    APO_Point d = line.end_point - line.begin_point;
    if (d.getMagnitude() < 1.0e-6)
    {
      return res;
    }

    APO_Point delta = p - circle.center;

    // root term:
    double term = std::pow(APO_Point::getDotProduct(d, delta), 2.0)
                  - std::pow(d.getMagnitude(), 2.0)
                        * (std::pow(delta.getMagnitude(), 2.0)
                           - std::pow(circle.radius, 2.0));

    // no intersection:
    if (term < 0.0)
    {
      return res;
    }

    // one or two intersections:
    double t1 = (-APO_Point::getDotProduct(d, delta) + sqrt(term))
                / std::pow(d.getMagnitude(), 2.0);
    double t2;
    bool tangent = false;

    // only one intersection:
    if (fabs(term) < APO_TOLERANCE)
    {
      t2      = t1;
      tangent = true;
    }

    // two intersections
    else
    {
      t2 = (-APO_Point::getDotProduct(d, delta) - sqrt(term))
           / std::pow(d.getMagnitude(), 2.0);
    }

    APO_Point sol1;
    APO_Point sol2(0, 0, false);

    sol1 = p + d * t1;

    if (!tangent)
    {
      sol2 = p + d * t2;
    }

    if (!limited || pointIsOnLine(sol1, line, true))
    {
      res.push_back(sol1);
    }
    if (sol2.isValid())
    {
      if (!limited || pointIsOnLine(sol2, line, true))
      {
        res.push_back(sol2);
      }
    }
    // ret.setTangent(tangent);

    // tangent with two intersections very close to each other:
    if (res.size() == 2 && res[0].equalsFuzzy(res[1]))
    {
      res.pop_back();
    }

    return res;

}

template <class T, class E>
bool
isIntersectWith(const T& t, const E& e, bool limited)
{
  auto&& ips = getIntersectionPoints<T, E>(t, e, limited);
  return !ips.empty();
}

template <class T>
inline std::vector<T>
verify(const std::vector<T>& candidates, const APO_Shape& shape1,
       const APO_Shape& shape2, const APO_Shape& shape3)
{
  return std::vector<T>();
}

template <class T>
std::vector<T>
removeDuplicates(const std::vector<T>& shaps)
{
  auto cmp = [&](const APO_Shape& a, const APO_Shape& b)
  {
    if (APO_IS_LINE(a) && APO_IS_LINE(b))
    {
      return a.line.begin_point.equalsFuzzy(b.line.begin_point)
             && b.line.end_point.equalsFuzzy(b.line.end_point);
    }
    if (APO_IS_CIRCLE(a) && APO_IS_CIRCLE(b))
    {
      return a.circle.center.equalsFuzzy(b.circle.center)
             && fuzzyCompare(a.circle.radius, b.circle.radius);
    }
    return false;
  };
  

  std::set<T, decltype(cmp)> sets(cmp);
  for (auto&& s : shaps)
  {
    sets.insert(s);
  }

  std::vector<T> res;
  for (auto it = sets.begin(); it != sets.end(); ++it)
  {
    res.push_back(*it);
  }

  return res;
}

/* --------------------------------- solves --------------------------------- */

/// http://www.geometer.org/mathcircles/inversion.pdf
bool
getInverseShape(const APO_Shape& shp, const APO_Shape& inversionCircle,
                APO_Shape& inversed)
{
  // inverse point
  // https://mathworld.wolfram.com/InversePoints.html
  //
  // Points, also called polar reciprocals, which are transformed into each
  // other through inversion about a given inversion circle C (or inversion
  // sphere). The points P and P^' are inverse points with respect to the
  // inversion circle if
  // OP * OP' = r^2
  // In this case, P^' is called the inversion pole and the line L through P
  // and perpendicular to OP is called the polar. In the above figure, the
  // quantity r^2 is called the circle power of the point P relative to the
  // circle C.
  if (APO_IS_POINT(shp))
  {
    double r         = inversionCircle.circle.radius;
    APO_Point center = inversionCircle.circle.center;
    double d         = shp.point.getDistanceTo(center);
    if (fabs(d) < APO_TOLERANCE)
    {
      inversed = APO_Shape(shp.point);
      return true;
    }

    double d_inverse = pow(r, 2) / d;
    inversed         = APO_Shape(
        APO_Point(center.x + (shp.point.x - center.x) * d_inverse / d,
                          center.y + (shp.point.y - center.y) * d_inverse / d));
    return true;
  }
  // line inverse
  if (APO_IS_LINE(shp))
  {
    APO_Point center = inversionCircle.circle.center;
    // A "line" that passes through O is inverted to itself. Note, of course
    // that the individual points of the "line" are inverted to other points
    // on the "line" except for the two points where it passes through k.
    if (pointIsOnLine(center, shp.line, false))
    {
      inversed = APO_Shape(shp.line);
      return true;
    }
    else
    {
      // Every "line" that does not pass through O is inverted to a circle
      // (no quotes: a real circle) that passes through O.
      APO_Line s;
      s.begin_point = center;
      s.end_point   = shp.line.getClosestPoint(center, false);
      // intersection_points from shp.line and s
      std::vector<APO_Point> ips = getIntersectionPoints(shp.line, s, false);
      if (ips.size() == 1)
      {
        APO_Point p = ips[0];
        APO_Shape pinverse;
        if (getInverseShape(p, inversionCircle, pinverse))
        {
          inversed = APO_Shape(APO_Circle::createFrom2Points(center, p));
          return true;
        }
      }
    }
  }

  if (APO_IS_CIRCLE(shp))
  {
    /// concentric circles
    APO_Circle circle = shp.circle;
    if (circle.center.equalsFuzzy(inversionCircle.circle.center))
    {
      APO_Shape inversed_point_shp;
      getInverseShape(
          APO_Point(circle.center.x + circle.radius, circle.center.y),
          inversionCircle, inversed_point_shp);
      double radius = circle.center.x - inversed_point_shp.point.x;
      if (radius < 0)
      {
        radius = fabs(radius);
      }
      inversed = APO_Shape(APO_Circle(circle.center, radius));
      return true;
    }
    else if (pointIsOnCircle(inversionCircle.circle.center, circle))
    {
      APO_Line s(inversionCircle.circle.center, circle.center);
      std::vector<APO_Point> ips = getIntersectionPoints(s, circle, false);
      if (ips.size() > 0)
      {
        return false;
      }

      APO_Point p = ips[0];
      if (p.equalsFuzzy(inversionCircle.circle.center))
      {
        if (ips.size() < 2)
        {
          return false;
        }
        p = ips[1];
      }

      APO_Shape pinverse;
      if (!getInverseShape(p, inversionCircle, pinverse))
      {
        return false;
      }
      inversed = APO_Shape(
          APO_Line(pinverse.point, s.getAngle() + M_PI / 2.0, 1.0));
      return true;
    }
    else
    {
      APO_Line l(inversionCircle.circle.center, circle.center);
      std::vector<APO_Point> ips = getIntersectionPoints(l, circle, false);
      if (ips.empty())
      {
        return false;
      }

      APO_Point p1 = ips[0];
      APO_Point p2 = ips[1];

      APO_Shape p1inverse, p2inverse;
      if (!getInverseShape(p1, inversionCircle, p1inverse))
      {
        return false;
      }
      if (!getInverseShape(p2, inversionCircle, p2inverse))
      {
        return false;
      }
      inversed = APO_Shape(
          APO_Circle::createFrom2Points(p1inverse.point, p2inverse.point));
      return true;
    }
  }

  return false;
}

std::vector<APO_Line>
getCircleTangentsThroughPoint(const APO_Circle& circle, const APO_Point& p)
{
  std::vector<APO_Line> res;
  // used when creating tangential circles to two parallel lines and point:
  if (fabs(circle.radius) < APO_TOLERANCE)
  {
    res.push_back(APO_Line(p, circle.center));
    res.push_back(APO_Line(p, circle.center));
    return res;
  }
  // point on the circle line (produces error):
  else if (pointIsOnCircle(p, circle))
  {
    APO_Line s(p, circle.center);
    res.push_back(APO_Line(p, s.getAngle() + M_PI / 2.0, 1.0));
    return res;
  }
  // pointis inside the circle:
  else if (circle.containsPoint(p))
  {
    return res;
  }
  // point outside circle:
  else
  {
    APO_Circle circle2 = APO_Circle::createFrom2Points(p, circle.center);
    std::vector<APO_Point> touching_points = getIntersectionPoints(circle2, circle);
    if (touching_points.size() == 1)
    {
      res.push_back(APO_Line(p, touching_points[0]));
    }
    else if (touching_points.size() == 2)
    {
      res.push_back(APO_Line(p, touching_points[0]));
      res.push_back(APO_Line(p, touching_points[1]));
    }
  }
  return res;
}

std::vector<APO_Line>
getAllTangents(const APO_Shape& shp1, const APO_Shape& shp2)
{
  // Collection of tangent lines to return
  std::vector<APO_Line> tangents;
  // Distance between center points, difference of radii, sum of radii
  double dC1C2, rDiff, rSum;
  // Tangent being constructed
  APO_Line tangent;
  // Angles
  double a1, a2, at;
  // Offset points
  APO_Point offset1, offset2;
  // Line shape
  APO_Line line;
  // Intersection points
  std::vector<APO_Point> ips;

  // Most common usage: Handle 2 circle shapes first:
  if (APO_IS_CIRCLE(shp1) && APO_IS_CIRCLE(shp2))
  {
    // Normalized circle shape clones, dummy to order c1-c2 on size
    APO_Circle c1, c2;
    // Radii of circle shapes
    double c1Radius, c2Radius;
    // Centers of circle shapes
    APO_Point c1Center, c2Center;

    // Validate first circle shape:
    if (std::isnan(shp1.circle.radius))
    {
      return tangents; // Empty, invalid radius
    }
    else if (fabs(shp1.circle.radius) < APO_TOLERANCE)
    {
      // Handle as zero sized 2D circle:
      c1        = shp1.circle;
      c1.radius = 0.0;
    }
    else
    {
      // Handle as normalized 2D circle:
      c1        = shp1.circle;
      c1.radius = fabs(shp1.circle.radius);
    }

    // Validate second circle shape:
    if (std::isnan(shp2.circle.radius))
    {
      return tangents; // Empty, invalid radius
    }
    else if (fabs(shp2.circle.radius) < APO_TOLERANCE)
    {
      // Handle as zero sized 2D circle:
      c2        = shp2.circle;
      c2.radius = 0.0;
    }
    else
    {
      // Handle as normalized 2D circle:
      c2        = shp2.circle;
      c2.radius = fabs(shp2.circle.radius);
    }

    // Ensure that c1 is the smaller circle:
    // Does not swap equal sized circles
    if (c1.radius > c2.radius)
    {
      std::swap(c1, c2);
    }

    // With 2 valid circle shapes;
    c1Radius = c1.radius;
    c1Center = c1.center;
    c2Radius = c2.radius;
    c2Center = c2.center;
    // Not expecting NaN with 2 valid circle shapes:
    dC1C2 = c1Center.getDistanceTo(c2Center);

    // Reject (almost) concentric circles:
    if (dC1C2 < 1e-6)
    {
      return tangents; // Empty, no solutions
    }

    // (Almost) internally touching circles:
    if (fuzzyCompare(dC1C2 + c1Radius, c2Radius))
    { // APO_TOLERANCE
      tangent = APO_Line(c2Center, c1Center);
      // With 2 radii larger than zero:
      if (c1Radius > 0.0)
      {
        tangent.setLength((dC1C2 + c1Radius) * 2, true); // fromStart
      }
      // Handle point on circle here instead of externally touching:
      // Ensuring that the single valid tangent is the first
      else
      {
        tangent.setLength(c2Radius * 2, true); // fromStart
      }
      tangent.rotate(M_PI / 2, tangent.getMiddlePoint());

      // First and final solution:
      tangents.push_back(tangent);
      return tangents; // One single solution
    }

    // Exclude other nested circles:
    if (dC1C2 + c1Radius < c2Radius)
    {
      return tangents; // Empty, no solutions
    }

    // Include external tangents:
    rDiff = c2Radius - c1Radius;
    if (dC1C2 > rDiff)
    {
      a1      = c1Center.getAngleTo(c2Center);
      a2      = std::asin(rDiff / dC1C2);
      offset1 = APO_Point();
      offset2 = APO_Point();

      // First solution:
      at = a1 + a2 + M_PI / 2.0;
      offset1.setPolar(c1Radius, at);
      offset2.setPolar(c2Radius, at);
      tangents.push_back(APO_Line(c1Center + offset1, c2Center + offset2));

      // Second solution, exclude for R1=R2=zero:
      if (c2Radius < APO_TOLERANCE)
      {
        // No second solution
      }
      else
      {
        at = a1 - a2 - M_PI / 2.0;
        offset1.setPolar(c1Radius, at);
        offset2.setPolar(c2Radius, at);
        tangents.push_back(APO_Line(c1Center + offset1, c2Center + offset2));
      }
    }
    // No external tangents:
    else
    {
      // No tangents
    }

    // (Almost) externally touching circles:
    rSum = c2Radius + c1Radius;
    if (fuzzyCompare(dC1C2, rSum))
    { // APO_TOLERANCE
      tangent = APO_Line(c2Center, c1Center);
      tangent.setLength(c2Radius * 2, true); // fromStart
      tangent.rotate(M_PI / 2, tangent.getMiddlePoint());

      // Third and final solution:
      tangents.push_back(tangent);
      // only one solution

      return tangents;
    }

    // Include internal tangents but only for radii larger than zero:
    if (dC1C2 > rSum && c1Radius > 0.0)
    {
      a1      = c1Center.getAngleTo(c2Center);
      a2      = std::asin(rSum / dC1C2);
      offset1 = APO_Point();
      offset2 = APO_Point();

      // Third solution:
      at = a1 + a2 + M_PI / 2.0;
      offset1.setPolar(c1Radius, at);
      offset2.setPolar(c2Radius, at);
      tangents.push_back(APO_Line(c1Center - offset1, c2Center + offset2));

      // Fourth solution:
      at = a1 - a2 - M_PI / 2.0;
      offset1.setPolar(c1Radius, at);
      offset2.setPolar(c2Radius, at);
      tangents.push_back(APO_Line(c1Center - offset1, c2Center + offset2));
    }
    // No internal tangents:
    else
    {
      // No tangents
    }

    return tangents;
  } // End 2 circles

  // With 2 line shapes (Circles with infinite radii):
  else if (APO_IS_LINE(shp1) && APO_IS_LINE(shp2))
  {
    APO_Line l1, l2;
    l1 = APO_Line(shp1.line.begin_point, shp1.line.end_point);
    l2 = APO_Line(shp2.line.begin_point, shp2.line.end_point);

    // The angle of a near zero-length line is zero by default (APO_Point::getAngle())
    // Exclude solutions for a line with almost no length:
    if (l1.getLength() <= 1.0e-6 || l2.getLength() <= 1.0e-6)
    {
      return tangents; // Empty, not processable line(s)
    }

    // Diversify on crossing or not:
    if (!isIntersectWith(shp1.line, shp2.line, false))
    {
      dC1C2 = l1.getDistanceTo(l2.begin_point, false); // unlimited
      // May return NaN, comparing with NaN is always false
      if (std::isnan(dC1C2) || dC1C2 > APO_TOLERANCE)
      {
        return tangents; // Empty, incorrect data or parallel at a distance
      }
      else
      {
        // Include one representation of itself when collinear:
        ips.push_back(l1.begin_point);
        ips.push_back(l1.end_point);
        ips.push_back(l2.begin_point);
        ips.push_back(l2.end_point);
        std::sort(ips.begin(), ips.end(),
                  [](const APO_Point& v1, const APO_Point& v2)
                  { return v1.y > v2.y || (v1.y == v2.y && v1.x < v2.x); });
        tangents.push_back(APO_Line(ips[0], ips[3]));
      }
    }
    // With an intersection point:
    else
    {
      a1            = (l1.getAngle() + l2.getAngle()) / 2;
      a2            = a1 + M_PI / 2;
      double length = (l1.getLength() + l2.getLength()) / 2;
      tangents.push_back(APO_Line(ips[0], a1, length));
      tangents.push_back(APO_Line(ips[0], a2, length));
    }
    return tangents; // No or 1-2 solution(s)
  } // End line-line

  APO_Shape c1, c2;
  APO_Shape item1, item2;

  // With a circle and a line shape:
  if (APO_IS_CIRCLE(shp1) && APO_IS_LINE(shp2))
  {

    line = APO_Line(shp2.line.begin_point, shp2.line.end_point);

    // The angle of a near zero-length line is zero by default (APO_Point::getAngle())
    // Exclude solutions for a line with almost no length:
    if (line.getLength() <= 1.0e-6)
    {
      return tangents; // Empty, not processable line
    }

    // Handle first shape as normalized 2D circle:
    c2               = shp1;
    c2.circle.radius = fabs(shp1.circle.radius);

  } // End circle-line

  // With a line and a circle shape:
  else if (APO_IS_LINE(shp1) && APO_IS_CIRCLE(shp2))
  {

    line = shp1.line;

    // The angle of a near zero-length line is zero by default (APO_Point::getAngle())
    // Exclude solutions for a line with almost no length:
    if (line.getLength() <= 1.0e-6)
    {
      return tangents; // Empty, not processable line
    }

    // Handle second shape as normalized 2D circle:
    c2               = shp2;
    c2.circle.radius = fabs(shp2.circle.radius);

  } // End line-circle

  // Least common usage: Support RPoint shapes as zero sized APO_Circle shapes:
  // Circle shapes are not guaranteed to be valid circles
  // Shapes are further validated and handled by a recursive call
  else if (APO_IS_POINT(shp1) || APO_IS_POINT(shp2))
  {
    // Ifso, convert first point into an APO_Circle:
    if (APO_IS_POINT(shp1))
    {
      item1 = APO_Shape(APO_Circle(shp1.point, 0.0));
    }
    else
    {
      item1 = shp1;
    }

    // Ifso, convert second point into an APO_Circle:
    if (APO_IS_POINT(shp2))
    {
      item2 = APO_Shape(APO_Circle(shp2.point, 0.0));
    }
    else
    {
      item2 = shp2;
    }
    // Handle the occurrence of points as circles:
    return getAllTangents(item1, item2);
  }
  // With any unsupported shape:
  else
  {
    return tangents; // Empty, incorrect data
  }

  // With validated line and circle shapes:
  // Solutions are tangent to the circle and in special tangent to the line at infinity
  double c2Radius    = c2.circle.radius;
  APO_Point c2Center = c2.circle.center;
  a1                 = line.getAngle();
  // Arbitrary tangent length
  double c1Radius = std::max(c2Radius, line.getLength() / 2);

  // Define a perpendicular diameter:
  APO_Point ip = line.getClosestPoint(c2Center, false); // unlimited
  if (ip.equalsFuzzy(c2Center))
  {                                                 // APO_TOLERANCE
    ips = getIntersectionPoints(line, c2.circle, false);// unlimited
    line = APO_Line(c2Center, ips[0]);
    line.rotate(M_PI / 2, c2Center);
  }
  else
  {
    line = APO_Line(c2Center, ip);
  }
  line.setLength(c2Radius, true); // fromStart
  line.reverse();
  line.setLength(c2Radius * 2, true); // fromStart

  // Include first parallels tangent at circle:
  tangents.push_back(APO_Line(line.begin_point, a1, c1Radius));
  // Include a second solution for R2>zero:
  if (c2Radius > 0.0)
  {
    a2 = a1 + M_PI;
    tangents.push_back(APO_Line(line.end_point, a2, c1Radius));
  }
  return tangents; // 1-2 solutions
}

std::vector<APO_Circle>
getSolutions(const std::vector<APO_Shape>& shps)
{
  assert(shps.size() == 3);

  return getSolution(shps[0], shps[1], shps[2]);
}

std::vector<APO_Circle>
getSolution(const APO_Shape& shp1, const APO_Shape& shp2,
            const APO_Shape& shp3)
{

  std::vector<APO_Point> pointObjs;
  std::vector<APO_Line> lineObjs;
  std::vector<APO_Circle> circleObjs;

#define OBJECT_CLASSIFICATION(shp)                                            \
  switch (shp.type)                                                           \
  {                                                                           \
  case APO_POINT_TYPE:                                                        \
    pointObjs.push_back(shp.point);                                           \
    break;                                                                    \
  case APO_LINE_TYPE:                                                         \
    lineObjs.push_back(shp.line);                                             \
    break;                                                                    \
  case APO_CIRCLE_TYPE:                                                       \
    circleObjs.push_back(shp.circle);                                         \
    break;                                                                    \
  };

  OBJECT_CLASSIFICATION(shp1);
  OBJECT_CLASSIFICATION(shp2);
  OBJECT_CLASSIFICATION(shp3);

  if (pointObjs.size() == 3)
  {

    return getSolutionFromPPP(pointObjs[0], pointObjs[1], pointObjs[2]);
  }
  else if (pointObjs.size() == 2)
  {
    if (circleObjs.size() == 1)
    {
      return getSolutionFromPPC(pointObjs[0], pointObjs[1], circleObjs[0]);
    }
    else if (lineObjs.size() == 1)
    {
      return getSolutionFromPPL(pointObjs[0], pointObjs[1], lineObjs[0]);
    }
  }
  else if (pointObjs.size() == 1)
  {
    if (circleObjs.size() == 2)
    {
      return getSolutionFromPCC(pointObjs[0], circleObjs[0], circleObjs[1]);
    }
    else if (lineObjs.size() == 2)
    {
      return getSolutionFromPLL(pointObjs[0], lineObjs[0], lineObjs[1]);
    }
    else if (circleObjs.size() == 1 && lineObjs.size() == 1)
    {
      return getSolutionFromPLC(pointObjs[0], lineObjs[0], circleObjs[0]);
    }
  }
  else if (pointObjs.size() == 0)
  {
    if (lineObjs.size() == 3)
    {
      return getSolutionFromLLL(lineObjs[0], lineObjs[1], lineObjs[2]);
    }
    else if (lineObjs.size() == 2 && circleObjs.size() == 1)
    {
      return getSolutionFromLLC(lineObjs[0], lineObjs[1], circleObjs[0]);
    }
    else if (lineObjs.size() == 1 && circleObjs.size() == 2)
    {
      return getSolutionFromLCC(lineObjs[0], circleObjs[0], circleObjs[1]);
    }
    else if (circleObjs.size() == 3)
    {
      return getSolutionFromCCC(circleObjs[0], circleObjs[1], circleObjs[2]);
    }
  }

#undef OBJECT_CLASSIFICATION
  return std::vector<APO_Circle>();
}

std::vector<APO_Circle>
getSolutionFromPPP(const APO_Point& point1, const APO_Point& point2,
                   const APO_Point& point3)
{
  return { APO_Circle::createFrom3Points(point1, point2, point3) };
}

std::vector<APO_Circle>
getSolutionFromPPC(const APO_Point& point1, const APO_Point& point2,
                   const APO_Circle& circle)
{
  auto ppc = [](const APO_Point& p1, const APO_Point& p2,
                const APO_Circle& circle) -> std::vector<APO_Circle>
  {
    APO_Line l1(circle.center, p1);
    APO_Point m = APO_Point::getAverage(p1, p2);
    APO_Line l2(m, p2.getAngleTo(p1) + M_PI / 2.0, 1.0);

    auto&& ips = getIntersectionPoints(l1, l2, false);
    if (ips.size() != 1)
    {
      return std::vector<APO_Circle>();
    }
    return { APO_Circle(ips[0], ips[0].getDistanceTo(p1)) };
  };

  std::vector<APO_Circle> res;

  if (0 == pointIsOnCircle(point1, circle)
      && 0 == pointIsOnCircle(point2, circle))
  {
    return { APO_Circle(circle.center, circle.radius) };
  }

  APO_Point pOnCircle, pOther;
  if (0 == pointIsOnCircle(point1, circle))
  {
    pOnCircle = point1;
    pOther    = point2;
    return ppc(pOnCircle, pOther, circle);
  }
  if (0 == pointIsOnCircle(point2, circle))
  {
    pOnCircle = point2;
    pOther    = point1;
    return ppc(pOnCircle, pOther, circle);
  }

  APO_Shape inversion_cirlce;
  inversion_cirlce.type          = APO_CIRCLE_TYPE;
  inversion_cirlce.circle.center = pOnCircle;
  inversion_cirlce.circle.radius = 10.0;
  APO_Shape circle_inverse;
  APO_Shape point2_inverse;
  if (!getInverseShape(circle, inversion_cirlce, circle_inverse)
      || !getInverseShape(point2, inversion_cirlce, point2_inverse))
  {
    return res;
  }

  std::vector<APO_Line> lines = getCircleTangentsThroughPoint(
      circle_inverse.circle, point2_inverse.point);
  for (int i = 0; i < lines.size(); ++i)
  {
    APO_Shape res_circle;
    getInverseShape(lines[i], inversion_cirlce, res_circle);
    res.push_back(
        APO_Circle(res_circle.circle.center, res_circle.circle.radius));
  }
  return res;
}

std::vector<APO_Circle>
getSolutionFromPPL(const APO_Point& point1, const APO_Point& point2,
                   const APO_Line& line)
{
  auto ppl = [](const APO_Point& point1, const APO_Point& point2,
                const APO_Line& line)
  {
    std::vector<APO_Circle> res;
    APO_Shape inversion_cirlce(APO_Circle(point1, 10));
    APO_Shape line_inverse;
    APO_Shape point2_inverse;
    if (!getInverseShape(line, inversion_cirlce, line_inverse)
        || !getInverseShape(point2, inversion_cirlce, point2_inverse))
    {
      return res;
    }

    std::vector<APO_Line> lines = getCircleTangentsThroughPoint(
        line_inverse.circle, point2_inverse.point);
    for (int i = 0; i < lines.size(); ++i)
    {
      APO_Shape res_circle;
      getInverseShape(lines[i], inversion_cirlce, res_circle);
      res.push_back(
          APO_Circle(res_circle.circle.center, res_circle.circle.radius));
    }
    return res;
  };

  bool swapPoint = false;
  if (pointIsOnLine(point1, line, true))
  {
    swapPoint = true;
    if (pointIsOnLine(point2, line, true))
    {
      return std::vector<APO_Circle>();
    }
  }

  if (swapPoint)
  {
    return ppl(point2, point1, line);
  }
  else
  {
    return ppl(point1, point2, line);
  }
}

std::vector<APO_Circle>
getSolutionFromPCC(const APO_Point& point, const APO_Circle& circle1,
                   const APO_Circle& circle2)
{
  std::vector<APO_Circle> res;

  // relative sized inversion circle:
  double r_inv = circle1.radius;
  if (circle2.radius > circle1.radius)
  {
    r_inv = circle2.radius;
  }
  APO_Shape inversionCircle(APO_Circle(point, r_inv));

  // construct inversion shape:
  APO_Shape c1Inverse;
  APO_Shape c2Inverse;
  if (!getInverseShape(circle1, inversionCircle, c1Inverse)
      || !getInverseShape(circle2, inversionCircle, c2Inverse))
  {
    return res;
  }

  // Get all tangent shapes for given inversion shapes:
  // Exploits an enhanced algorithm by CVH
  auto&& tangents = getAllTangents(c1Inverse, c2Inverse);

  // Return the re-inversion of all tangents (0-4 solutions):
  auto&& inversed = getInverseShapes<APO_Line>(tangents, inversionCircle);
  for (auto&& ii : inversed)
  {
    if (APO_IS_CIRCLE(ii))
      res.push_back(ii.circle);
  }
  return res;
}

std::vector<APO_Circle>
getSolutionFromPLL(const APO_Point& point, const APO_Line& line1,
                   const APO_Line& line2)
{
  std::vector<APO_Circle> res;

  // intersection between two lines line1, line2:
  std::vector<APO_Point> ipsLL = getIntersectionPoints(line1, line2, false);
  std::vector<APO_Circle> circles;
  std::vector<APO_Point> centers;

  std::vector<APO_Line> bisectorLines = getAngleBisectors(line1, line2);

  bool onLine1 = pointIsOnLine(point, line1, false);
  bool onLine2 = pointIsOnLine(point, line2, false);

  bool onBisector = false;
  for (size_t k = 0; k < bisectorLines.size(); k++)
  {
    if (pointIsOnLine(point, bisectorLines[k], false))
    {
      onBisector = true;
      break;
    }
  }

  // lines are parallel:
  if (ipsLL.size() == 0)
  {
    // middle line:
    APO_Point s         = line1.begin_point;
    APO_Point p         = line2.getClosestPoint(s, false);
    APO_Point center    = APO_Point::getAverage(s, p);
    APO_Line middleLine = line1;
    middleLine.move(center - middleLine.begin_point);
    // circle with radius c-s around point:
    APO_Circle circle = APO_Circle(point, center.getDistanceTo(s));
    // intersections between circle and middle line are candidates:
    centers = getIntersectionPoints(middleLine, circle, false);
  }

  // point is on line1 or line2:
  else if (onLine1 || onLine2)
  {
    APO_Line line = onLine1 ? line1 : line2;
    APO_Line orthoLine(point, line.getAngle() + M_PI / 2, 1.0);
    for (size_t k = 0; k < bisectorLines.size(); k++)
    {
      auto&& bisectorLine = bisectorLines[k];
      auto&& ips = getIntersectionPoints(bisectorLine, orthoLine, false);
      if (ips.size() != 1)
        continue;
      centers.push_back(ips[0]);
    }
  }

  else
  {
    std::vector<APO_Point> centerCandidates;

    // point on bisector:
    if (onBisector)
    {
      if (ipsLL.size() != 1)
      {
        return res;
      }
      // distance from point to line1 (radius of circle around point, tangential to lines):
      double rp = line1.getDistanceTo(point, false);
      // distance from intersection line1/line2 to point:
      double dp = ipsLL[0].getDistanceTo(point);
      // distances from intersection line1/line2 to intersection of bisector line with circle around point, touching line1, line2:
      double dc1 = dp + rp;
      double dc2 = dp - rp;
      // factors to scale circle to reach results:
      double f1 = dp / dc1;
      double f2 = dp / dc2;
      // radius of solution
      double r1 = rp * f1;
      double r2 = rp * f2;

      for (size_t k = 0; k < bisectorLines.size(); ++k)
      {
        auto&& bisectorLine = bisectorLines[k];
        double a            = bisectorLine.getAngle();
        centerCandidates.push_back(ipsLL[0]
                                   + (APO_Point::createPolar(dp + r2, a)));
        centerCandidates.push_back(ipsLL[0]
                                   + (APO_Point::createPolar(dp - r1, a)));
        centerCandidates.push_back(
            ipsLL[0] + (APO_Point::createPolar(dp + r2, a + M_PI)));
        centerCandidates.push_back(
            ipsLL[0] + (APO_Point::createPolar(dp - r1, a + M_PI)));
      }
    }

    // circle C tangential to two lines with center E, radius 10:
    else
    {
      circles = getCircles2TR(line1, line2, 10.0);
      if (circles.size() == 0)
      {
        return res;
      }

      for (size_t i = 0; i < circles.size(); ++i)
      {
        APO_Circle circle = circles[i];
        APO_Point e       = circle.center;

        // line L from intersection between the two lines to the point:
        // center of solution is on this line
        APO_Line line;
        if (ipsLL.size() == 0)
        {
          // lines parallel:
          line = line1;
          line.move(point - line1.begin_point);
        }
        else
        {
          line = APO_Line(ipsLL[0], point);
        }

        // intersections between line L and circle C -> G, H:
        std::vector<APO_Point> ipsLC
            = getIntersectionPoints(line, circle, false);
        if (ipsLC.size() != 2)
        {
          continue;
        }

        APO_Point g = ipsLC[0];
        APO_Point h = ipsLC[1];

        // two lines L1, L2 with same angle as EG, EH through point:
        APO_Line l1(e, g);
        APO_Line l2(e, h);
        l1.move(point - l1.begin_point);
        l2.move(point - l2.begin_point);

        // intersection of angle bisector and lines L1, L2 are centers of candidates:
        for (size_t k = 0; k < bisectorLines.size(); ++k)
        {
          auto&& bisectorLine = bisectorLines[k];
          auto&& t1           = getIntersectionPoints(l1, bisectorLine, false);
          centerCandidates.insert(centerCandidates.end(), t1.begin(),
                                  t1.end());
          auto&& t2 = getIntersectionPoints(l2, bisectorLine, false);
          centerCandidates.insert(centerCandidates.end(), t2.begin(),
                                  t2.end());
        }
      }
    }

    for (size_t c = 0; c < centerCandidates.size(); ++c)
    {
      APO_Point centerCandidate = centerCandidates[c];
      double dLine1             = line1.getDistanceTo(centerCandidate, false);
      double dLine2             = line1.getDistanceTo(centerCandidate, false);
      double dPoint             = point.getDistanceTo(centerCandidate, false);
      if (fuzzyCompare(dLine1, dPoint) && fuzzyCompare(dLine2, dPoint))
      {
        centers.push_back(centerCandidate);
      }
    }
    centers = APO_Point::getUnique(centers);
  }

  for (size_t c = 0; c < centers.size(); ++c)
  {
    double r = centers[c].getDistanceTo(point);
    if (fuzzyCompare(r, 0.0))
    {
      continue;
    }
    res.push_back(APO_Circle(centers[c], r));
  }
  return res;
}

std::vector<APO_Circle>
getSolutionFromPLC(const APO_Point& point, const APO_Line& line,
                   const APO_Circle& circle)
{
  APO_Point a = point;
  APO_Point c = circle.center;
  APO_Point f = line.getClosestPoint(c, false);
  APO_Line ortho(c, f);

  APO_Point cf = f - c;
  APO_Point ce = cf;
  ce.setMagnitude(circle.radius);
  // intersections of orthogonal through c to line with circle:
  APO_Point d1              = c + ce;
  APO_Point e1              = c - ce;
  std::vector<APO_Point> ds = { d1, e1 };
  std::vector<APO_Point> es = { e1, d1 };

  std::vector<APO_Point> centerCandidates;
  std::vector<APO_Point> ips;

  for (size_t i = 0; i < 2; i++)
  {
    APO_Point d = ds[i];
    APO_Point e = es[i];

    // special case:
    // a is on orthogonal cf:
    if (ortho.getDistanceTo(a, false) < APO_TOLERANCE)
    {
      APO_Line da_(d, f);
      da_.rotate(M_PI / 4, d);
      APO_Line par = line;
      par.moveTo(a);

      ips = getIntersectionPoints(da_, par, false);
      if (ips.size() != 1)
        continue;

      APO_Point a_ = ips[0];
      
      ips = getIntersectionPoints(da_, line, false);
      if (ips.size() != 1)
        continue;

      APO_Point f_ = ips[0];

      APO_Line a_e(a_, e);
      APO_Line f_p = a_e;
      f_p.moveTo(f_);

      ips = getIntersectionPoints(f_p, ortho, false);
      if (ips.size() != 1)
        continue;

      APO_Point p = ips[0];

      APO_Point apm = APO_Point::getAverage(a, p);
      double m      = line.getDistanceTo(apm, false);
      APO_Circle circ(a, m);
      APO_Line par2 = line;
      par2.moveTo(apm);

      ips = getIntersectionPoints(
          par2, circ, false);
      centerCandidates.insert(centerCandidates.end(), ips.begin(), ips.end());
    }
    // special case:
    // point is on line:
    else if (pointIsOnLine(a, line, false))
    {
      // similarity axis:
      APO_Line ea(e, a);
      ips = getIntersectionPoints(ea, circle, false);
      if (ips.size() != 2)
      {
        continue;
      }
      APO_Point i1 = ips[0];
      if (i1.equalsFuzzy(e))
      {
        i1 = ips[1];
      }
      APO_Line ci1(c, i1);

      // ortho through a:
      APO_Line orthoA(a, line.getAngle() + M_PI / 2, 1.0);

      ips = getIntersectionPoints(orthoA, ci1, false);
      if (ips.size() != 1)
        continue;
      centerCandidates.push_back(ips[0]);
    }
    else
    {
      APO_Line da(d, a);
      ips = getIntersectionPoints(da, line, false);
      if (ips.size() != 1)
        continue;
      APO_Point m = ips[0];

      APO_Circle efa = APO_Circle::createFrom3Points(e, f, a);
      auto&& ps = getIntersectionPoints(da, efa, false);
      if (ps.size() < 1)
          continue;

      APO_Point p = ps[0];
      if (p.equalsFuzzy(a) && ps.size() > 1)
      {
        p = ps[1];
      }

      APO_Line ap(a, p);
      APO_Line apOrtho = ap;
      apOrtho.rotate(M_PI / 2, ap.getMiddlePoint());

      auto&& tangents = getCircleTangentsThroughPoint(efa, m);
      for (size_t k = 0; k < tangents.size(); ++k)
      {
        auto&& tangent = tangents[k];

        double mt    = tangent.getLength();
        APO_Point mu = line.end_point - line.begin_point;
        mu.setMagnitude(mt);
        APO_Point u    = m + mu;
        APO_Point v    = m - mu;
        APO_Line orthU = line;
        orthU.rotate(M_PI / 2, u);
        APO_Line orthV = line;
        orthV.rotate(M_PI / 2, v);

        auto&& pps = getIntersectionPoints(orthU, apOrtho, false);
        centerCandidates.insert(centerCandidates.end(), pps.begin(),
                                pps.end());
        pps = getIntersectionPoints(orthV, apOrtho, false);
        centerCandidates.insert(centerCandidates.end(), pps.begin(),
                                pps.end());
      }
    }
  }

  std::vector<APO_Circle> res;
  for (auto&& centerCandidate : centerCandidates)
  {
    res.push_back(
        APO_Circle(centerCandidate, centerCandidate.getDistanceTo(a)));
  }
  return res;
}

std::vector<APO_Circle>
getSolutionFromLLL(const APO_Line& line1, const APO_Line& line2,
                   const APO_Line& line3)
{
  /*
  situations:
  0: all lines are parallel (no solutions)
  3: lines 2 and 3 are parallel
  5: lines 1 and 3 are parallel
  6: lines 1 and 2 are parallel
  7: none of the lines are parallel (4 solutions)
  */

  int situation = 0;
  std::vector<APO_Line> angleBisectors1, angleBisectors2, angleBisectors3;
  if (isIntersectWith(line1, line2, false))
  {
    situation += 1;
  }
  if (isIntersectWith(line1, line3, false))
  {
    situation += 2;
  }
  if (isIntersectWith(line2, line3, false))
  {
    situation += 4;
  }

  if (situation == 3 || situation == 5 || situation == 7)
  {
    angleBisectors1 = getAngleBisectors(line1, line2);
  }
  if (situation == 3 || situation == 6 || situation == 7)
  {
    angleBisectors2 = getAngleBisectors(line1, line3);
  }
  if (situation == 5 || situation == 6)
  {
    angleBisectors3 = getAngleBisectors(line2, line3);
  }

  std::vector<APO_Point> centerPoints;

#define LLL_CONCAT(l1, l2)                                                    \
  ips = getIntersectionPoints(l1, l2, false);                                 \
  centerPoints.insert(centerPoints.end(), ips.begin(), ips.end());

  std::vector<APO_Point> ips;
  if (situation == 3 || situation == 7)
  {
    LLL_CONCAT(angleBisectors1[0], angleBisectors2[0])
    LLL_CONCAT(angleBisectors1[0], angleBisectors2[1])
    LLL_CONCAT(angleBisectors1[1], angleBisectors2[0])
    LLL_CONCAT(angleBisectors1[1], angleBisectors2[1])
  }
  else if (situation == 5)
  {
    LLL_CONCAT(angleBisectors1[0], angleBisectors3[0])
    LLL_CONCAT(angleBisectors1[0], angleBisectors3[1])
    LLL_CONCAT(angleBisectors1[1], angleBisectors3[0])
    LLL_CONCAT(angleBisectors1[1], angleBisectors3[1])
  }
  else if (situation == 6)
  {
    LLL_CONCAT(angleBisectors2[0], angleBisectors3[0])
    LLL_CONCAT(angleBisectors2[0], angleBisectors3[1])
    LLL_CONCAT(angleBisectors2[1], angleBisectors3[0])
    LLL_CONCAT(angleBisectors2[1], angleBisectors3[1])
  }

#undef LLL_CONCAT

  std::vector<APO_Circle> ret;
  for (size_t i = 0; i < centerPoints.size(); ++i)
  {
    APO_Point c   = centerPoints[i];
    double radius = line1.getDistanceTo(c, false);
    ret.push_back(APO_Circle(c, radius));
  }
  return ret;
}

std::vector<APO_Circle>
getSolutionFromLLC(const APO_Line& line1, const APO_Line& line2,
                   const APO_Circle& circle)
{
  std::vector<APO_Line> parallels1
      = line1.getOffset(circle.radius, 1, BothSides);
  std::vector<APO_Line> parallels2
      = line2.getOffset(circle.radius, 1, BothSides);

  std::vector<APO_Shape> arr1;
  std::vector<APO_Shape> arr2;
  std::vector<APO_Shape> arr3;
  std::vector<APO_Shape> arr4;

  arr1.push_back(APO_Shape(circle.center));
  arr2.push_back(APO_Shape(circle.center));
  arr3.push_back(APO_Shape(circle.center));
  arr4.push_back(APO_Shape(circle.center));

  arr1.push_back(APO_Shape(parallels1[0]));
  arr2.push_back(APO_Shape(parallels1[0]));
  arr3.push_back(APO_Shape(parallels1[1]));
  arr4.push_back(APO_Shape(parallels1[1]));

  arr1.push_back(APO_Shape(parallels2[0]));
  arr2.push_back(APO_Shape(parallels2[1]));
  arr3.push_back(APO_Shape(parallels2[0]));
  arr4.push_back(APO_Shape(parallels2[1]));

  std::vector<APO_Circle> cArr1 = getSolution(arr1[0], arr1[1], arr1[2]);
  std::vector<APO_Circle> cArr2 = getSolution(arr2[0], arr2[1], arr2[2]);
  std::vector<APO_Circle> cArr3 = getSolution(arr3[0], arr3[1], arr3[2]);
  std::vector<APO_Circle> cArr4 = getSolution(arr4[0], arr4[1], arr4[2]);

  std::vector<APO_Circle> ret;

  auto llc_inner = [&](std::vector<APO_Circle>& cArr)
  {
    for (size_t k = 0; k < cArr.size(); ++k)
    {
      auto&& obj           = cArr[k];
      APO_Circle tmpCircle = obj;
      APO_Circle obj1      = obj;
      APO_Circle obj2      = obj;
      obj1.radius += circle.radius;
      obj2.radius = fabs(obj2.radius - circle.radius);

      double d1 = obj1.center.getDistanceTo(circle.center);
      if ((fuzzyCompare(d1, obj1.radius + circle.radius)
           || fuzzyCompare(d1, fabs(obj1.radius - circle.radius)))
          && fuzzyCompare(line1.getDistanceTo(obj1.center, false), obj1.radius)
          && fuzzyCompare(line2.getDistanceTo(obj1.center, false),
                          obj1.radius))
      {

        ret.push_back(obj1);
      }

      double d2 = obj2.center.getDistanceTo(circle.center);
      if ((fuzzyCompare(d2, obj2.radius + circle.radius)
           || fuzzyCompare(d2, fabs(obj2.radius - circle.radius)))
          && fuzzyCompare(line1.getDistanceTo(obj2.center, false), obj2.radius)
          && fuzzyCompare(line2.getDistanceTo(obj2.center, false),
                          obj2.radius))
      {

        ret.push_back(obj2);
      }
    }
  };

  llc_inner(cArr1);
  llc_inner(cArr2);
  llc_inner(cArr3);
  llc_inner(cArr4);
  return ret;
}

std::vector<APO_Circle>
getSolutionFromLCC(const APO_Line& line, const APO_Circle& circle1,
                   const APO_Circle& circle2)
{
  auto lCC_Circle_Create
      = [&](std::vector<APO_Circle>& ret, const std::vector<APO_Circle>& cArr)
  {
    for (size_t k = 0; k < cArr.size(); ++k)
    {
      auto&& obj = cArr[k];

      APO_Circle tmpCircle = obj;
      APO_Circle obj1      = obj;
      APO_Circle obj2      = obj;
      obj1.radius += circle1.radius;
      obj2.radius = fabs(obj2.radius - circle1.radius);

      double obj1c1 = obj1.center.getDistanceTo(circle1.center);
      double obj1c2 = obj1.center.getDistanceTo(circle2.center);
      if ((fuzzyCompare(obj1c1, circle1.radius + obj1.radius)
           || fuzzyCompare(obj1c1, fabs(circle1.radius - obj1.radius)))
          && (fuzzyCompare(obj1c2, circle2.radius + obj1.radius)
              || fuzzyCompare(obj1c2, fabs(circle2.radius - obj1.radius)))
          && fuzzyCompare(line.getDistanceTo(obj1.center, false), obj1.radius))
      {

        ret.push_back(obj1);
      }

      double obj2c1 = obj2.center.getDistanceTo(circle1.center);
      double obj2c2 = obj2.center.getDistanceTo(circle2.center);
      if ((fuzzyCompare(obj2c1, circle1.radius + obj2.radius)
           || fuzzyCompare(obj2c1, fabs(circle1.radius - obj2.radius)))
          && (fuzzyCompare(obj2c2, circle2.radius + obj2.radius)
              || fuzzyCompare(obj2c2, fabs(circle2.radius - obj2.radius)))
          && fuzzyCompare(line.getDistanceTo(obj2.center, false), obj2.radius))
      {

        ret.push_back(obj2);
      }
    }
  };

  auto LCC =
      [lCC_Circle_Create](const APO_Line& line, const APO_Circle& circle1,
                          const APO_Circle& circle2) -> std::vector<APO_Circle>
  {
    // find solutions for tangent to:
    // center point of smaller circle,
    // concentric circles for larger circle with distance = radius of smaller circle
    // parallels to line with distance = radius of smaller circle
    std::vector<APO_Shape> arr1;
    std::vector<APO_Shape> arr2;
    std::vector<APO_Shape> arr3;
    std::vector<APO_Shape> arr4;

    arr1.push_back(APO_Shape(circle1.center));
    arr2.push_back(APO_Shape(circle1.center));
    arr3.push_back(APO_Shape(circle1.center));
    arr4.push_back(APO_Shape(circle1.center));

    APO_Circle circle21 = circle2;
    circle21.radius += circle1.radius;
    arr1.push_back(circle21);
    arr3.push_back(circle21);

    if (fuzzyCompare(circle1.radius, circle2.radius))
    {
      arr2.push_back(APO_Shape(circle2.center));
      arr3.push_back(APO_Shape(circle2.center));
    }
    else
    {
      APO_Shape circle22    = circle2;
      circle22.circle.radius = fabs(circle22.circle.radius - circle1.radius);
      arr2.push_back(circle22);
      arr4.push_back(circle22);
    }

    std::vector<APO_Line> parallels
        = line.getOffset(circle1.radius, 1, BothSides);
    arr1.push_back(APO_Shape(parallels[0]));
    arr2.push_back(APO_Shape(parallels[0]));
    arr3.push_back(APO_Shape(parallels[1]));
    arr4.push_back(APO_Shape(parallels[1]));

    auto&& cArr1 = getSolution(arr1[0], arr1[1], arr1[2]);
    auto&& cArr2 = getSolution(arr2[0], arr2[1], arr2[2]);
    auto&& cArr3 = getSolution(arr3[0], arr3[1], arr3[2]);
    auto&& cArr4 = getSolution(arr4[0], arr4[1], arr4[2]);

    std::vector<APO_Circle> ret;
    lCC_Circle_Create(ret, cArr1);
    lCC_Circle_Create(ret, cArr2);
    lCC_Circle_Create(ret, cArr3);
    lCC_Circle_Create(ret, cArr4);

    return ret;
  };

  if (circle1.radius > circle2.radius)
  {
    return LCC(line, circle2, circle1);
  }
  else
  {
    return LCC(line, circle1, circle2);
  }
}

std::vector<APO_Circle>
getSolutionFromCCC(const APO_Circle& c1, const APO_Circle& c2,
                   const APO_Circle& c3)
{
  std::vector<APO_Circle> ret;

  APO_Circle circle1 = c1;
  APO_Circle circle2 = c2;
  APO_Circle circle3 = c3;
  APO_Circle locus;
  bool realLocus = false;
  double rDiff = std::numeric_limits<double>::quiet_NaN();

  // special case: at least two circles are concentric: no solution:
  // # Enhanced by CVH #
  // Verify if circle 1 is concentric with circle 2:
  if (c1.center.equalsFuzzy(c2.center))
  {
    if (c1.center.equalsFuzzy(c3.center))
    {
      // Failed: 3 concentric circles
      return std::vector<APO_Circle>(); // No or infinite solutions
    }
    rDiff     = c1.radius - c2.radius;
    locus     = c1;
    circle1   = c3;
    realLocus = true;
  }
  // Verify if circle 1 is concentric with circle 3:
  else if (c1.center.equalsFuzzy(c3.center))
  {
    if (c1.center.equalsFuzzy(c2.center))
    {
      // Failed: 3 concentric circles
      return std::vector<APO_Circle>(); // No or infinite solutions
    }
    rDiff     = c1.radius - c3.radius;
    locus     = c1;
    circle1   = c2;
    realLocus = true;
  }
  // Verify if circle 2 is concentric with circle 3:
  else if (c2.center.equalsFuzzy(c3.center))
  {
    if (c2.center.equalsFuzzy(c1.center))
    {
      // Failed: 3 concentric circles
      return std::vector<APO_Circle>(); // No or infinite solutions
    }
    rDiff   = c2.radius - c3.radius;
    locus   = c2;
    circle1   = c1;
    realLocus = true;
  }

  // Handle two concentric circles:
  if (realLocus)
  {
    if (std::isnan(rDiff))
    {
      return std::vector<APO_Circle>(); // No solutions
    }
    if (fuzzyCompare(rDiff, 0.0))
    {
      return std::vector<APO_Circle>(); // Infinite solutions
    }

    // Get intersections of two concentric with the other circle and the locus:
    locus.radius = locus.radius - rDiff / 2;
    circle1.radius = circle1.radius + fabs(rDiff / 2);
    std::vector<APO_Point> ips = getIntersectionPoints(locus, circle1);
    circle1.radius = circle1.radius - fabs(rDiff);
    // Avoid inversion through the center with negative radii:
    if (circle1.radius > 0)
    {
      auto&& lips = getIntersectionPoints(locus, circle1);
      ips.insert(ips.end(), lips.begin(), lips.end());
    }

    // If any, create tangent circles at intersections;
    if (ips.size() == 0)
    {
      return std::vector<APO_Circle>(); // No solutions
    }
    for (size_t i = 0; i < ips.size(); ++i)
    {
      ret.push_back(APO_Circle(ips[i], fabs(rDiff / 2)));
    }

    ret = removeDuplicates(ret); // Most likely none
    return ret;
  }

  // special case: three circles of equal size:
  if (fuzzyCompare(c1.radius, c2.radius)
      && fuzzyCompare(c1.radius, c3.radius))
  {
    // add outer and inner circles to result:
    auto&& sol = APO_Circle::createFrom3Points(c1.center, c2.center, c3.center);
    if (sol.center.isValid())
    {
      APO_Circle sol1 = sol;
      APO_Circle sol2 = sol;
      sol1.radius = sol1.radius + c1.radius;
      sol2.radius = fabs(sol2.radius - c1.radius);
      ret.push_back(sol1);
      ret.push_back(sol2);
    }
  }
  // circle1 is always the smallest:
  else
  {
    if (c2.radius <= c1.radius && c2.radius <= c3.radius)
    {
      circle1 = c2;
      circle2 = c1;
      circle3 = c3;
    }

    if (c3.radius <= c1.radius && c3.radius <= c2.radius)
    {
      circle1 = c3;
      circle2 = c1;
      circle3 = c2;
    }
  }

  // special case: three circles intersect in one point:
  APO_Point commonIP
      = getCommonIntersectionPoint(circle1, circle2, circle3);
  if (commonIP.isValid())
  {
    APO_Circle inversionCircle(commonIP, 10);
    auto&& shapesInverse =  getInverseShapes<APO_Circle>(
        { circle1, circle2, circle3 },
                                 inversionCircle);

    if (APO_IS_LINE(shapesInverse[0]) && APO_IS_LINE(shapesInverse[1])
        && APO_IS_LINE(shapesInverse[2]))
    {

      auto&& circlesTouching = getSolutions(shapesInverse);
      auto&& inversed = getInverseShapes(circlesTouching, inversionCircle);
      for (auto&& ii : inversed)
      {
        if (APO_IS_CIRCLE(ii))
          ret.push_back(ii.circle);
      }
    }

    return ret;
  }


  // special case: each circle intersects the other two,
  // at least one intersects through two points:
  std::vector<APO_Point> ips12= getIntersectionPoints(circle1, circle2);
  std::vector<APO_Point> ips13= getIntersectionPoints(circle1, circle3);
  std::vector<APO_Point> ips23= getIntersectionPoints(circle2, circle3);

  size_t nIps12 = ips12.size();
  size_t nIps13 = ips13.size();
  size_t nIps23 = ips23.size();

  if (nIps12 > 0 && nIps13 > 0 && nIps23 > 0
      && (nIps12 == 2 || nIps13 == 2 || nIps23 == 2))
  {
    std::vector<APO_Circle> inversionCircles;

    if (ips12.size() == 2)
    {
      double r = ips12[0].getDistanceTo(ips12[1]);
      inversionCircles.push_back(APO_Circle(ips12[0], r));
      inversionCircles.push_back(APO_Circle(ips12[1], r));
    }
    if (ips13.size() == 2)
    {
      double r = ips13[0].getDistanceTo(ips13[1]);
      inversionCircles.push_back(APO_Circle(ips13[0], r));
      inversionCircles.push_back(APO_Circle(ips13[1], r));
    }
    if (ips23.size() == 2)
    {
      double r = ips23[0].getDistanceTo(ips23[1]);
      inversionCircles.push_back(APO_Circle(ips23[0], r));
      inversionCircles.push_back(APO_Circle(ips23[1], r));
    }

    for (size_t i = 0; i < inversionCircles.size(); ++i)
    {
      APO_Shape circle1Inverse;
      getInverseShape(circle1, inversionCircles[i], circle1Inverse);
      APO_Shape circle2Inverse;
      getInverseShape(circle2, inversionCircles[i], circle2Inverse);
      APO_Shape circle3Inverse;
      getInverseShape(circle3, inversionCircles[i], circle3Inverse);

      auto&& iSol = getSolution(circle1Inverse, circle2Inverse, circle3Inverse);
      auto&& sol  = getInverseShapes(iSol, inversionCircles[i]);
      for (auto&& _soli : sol)
      {
        if (APO_IS_CIRCLE(_soli))
          ret.push_back(_soli.circle);
      }
    }

    ret = removeDuplicates(ret);

    return ret;
  }

  APO_Point powerCenter = getPowerCenter(circle1, circle2, circle3);
  //constructionShapes.push(new RPoint(powerCenter));

  if (!powerCenter.isValid())
  {
    return ret;
  }

  auto&& similarityAxes = getSimilarityAxes(circle1, circle2, circle3);
  for (size_t i = 0; i < similarityAxes.size(); i++)
  {
    //constructionShapes.push(similarityAxes[i]);
    APO_Point p, pp, q, qq, r, rr;

    APO_Point pole1 = getPole(circle1, similarityAxes[i]);
    APO_Point pole2 = getPole(circle2, similarityAxes[i]);
    APO_Point pole3 = getPole(circle3, similarityAxes[i]);
    if (!pole1.isValid() || !pole2.isValid() || !pole3.isValid())
    {
      continue;
    }

    APO_Line ray1(powerCenter, pole1);
    APO_Line ray2(powerCenter, pole2);
    APO_Line ray3(powerCenter, pole3);

    auto&& ips1 = getIntersectionPoints(ray1, circle1, false);
    auto&& ips2 = getIntersectionPoints(ray2, circle2, false);
    auto&& ips3 = getIntersectionPoints(ray3, circle3, false);

    double gotPoints = false;
    if (circle1.containsPoint(powerCenter)
        || circle2.containsPoint(powerCenter)
        || circle3.containsPoint(powerCenter))
    {
      std::vector<APO_Point> ipsRight, ipsLeft;

      auto _ipssProc =
          [&](const std::vector<APO_Point> & ips)
      {
        for (size_t n = 0; n < ips.size(); ++n)
        {
          auto&& ip = ips[n];
          if (similarityAxes[i].getSideOfPoint(ip) == Side::Right)
          {
            ipsRight.push_back(ip);
          }
          else
          {
            ipsLeft.push_back(ip);
          }
        }
      };

      _ipssProc(ips1);
      _ipssProc(ips2);
      _ipssProc(ips3);

      if (ipsRight.size() == 3 && ipsLeft.size() == 3)
      {
        p         = ipsRight[0];
        q         = ipsRight[1];
        r         = ipsRight[2];
        pp        = ipsLeft[0];
        qq        = ipsLeft[1];
        rr        = ipsLeft[2];
        gotPoints = true;
      }
    }

    if (!gotPoints)
    {
      std::sort(ips1.begin(), ips1.end(),
                [&](const APO_Point& v1, const APO_Point& v2) {
                  return powerCenter.getDistanceTo(v1)
                         < powerCenter.getDistanceTo(v2);
                });
      std::sort(ips2.begin(), ips2.end(),
                [&](const APO_Point& v1, const APO_Point& v2) {
                  return powerCenter.getDistanceTo(v1)
                         < powerCenter.getDistanceTo(v2);
                });
      std::sort(ips3.begin(), ips3.end(),
                [&](const APO_Point& v1, const APO_Point& v2) {
                  return powerCenter.getDistanceTo(v1)
                         < powerCenter.getDistanceTo(v2);
                });              

      if (ips1.size() != 2 || ips2.size() != 2 || ips3.size() != 2)
      {
        continue;
      }

      // alpha: +
      if (i == 0 || i == 3)
      {
        p  = ips1[0];
        pp = ips1[1];
      }
      // alpha: -
      else
      {
        p  = ips1[1];
        pp = ips1[0];
      }   

      // beta: +
      if (i == 0 || i == 2)
      {
        q  = ips2[0];
        qq = ips2[1];
      }
      // beta: -
      else
      {
        q  = ips2[1];
        qq = ips2[0];
      } 

      // gamma: +
      if (i == 0 || i == 1)
      {
        r  = ips3[0];
        rr = ips3[1];
      }
      // gamma: -
      else
      {
        r  = ips3[1];
        rr = ips3[0];
      }
    }

    ret.push_back(APO_Circle::createFrom3Points(p, q, r));
    ret.push_back(APO_Circle::createFrom3Points(pp, qq, rr));
  }

  auto&& cccAltRet = getSolutionsCCCAlt(c1, c2, c3);
  ret.insert(ret.end(), cccAltRet.begin(), cccAltRet.end());
  ret = removeDuplicates(ret);
  ret = verify(ret, c1, c2, c3);
  return ret;
}

bool
pointIsOnLine(const APO_Point& point, const APO_Line& line, bool limited)
{
  APO_Point vt = line.getClosestPoint(point, limited);
  if (!vt.isValid())
      return false;

  double vtm = vt.getMagnitude();
  return vtm < 1.0e-4 ? true : false;
}
bool
pointIsOnCircle(const APO_Point& point, const APO_Circle& circle)
{
  double d    = std::numeric_limits<double>::quiet_NaN();
  APO_Point v = point - circle.center;

  // point is at the center of the circle, infinite solutions:
  if (v.getMagnitude() < APO_TOLERANCE)
  {
    return false;
  }

  APO_Point vt = APO_Point::createPolar(v.getMagnitude() - circle.radius,
                                        v.getAngle());

  if (vt.isValid())
  {
    d = vt.getMagnitude();
  }

  if (std::isnan(d))
  {
    return false;
  }
  return d < APO_TOLERANCE;
}

bool
fuzzyCompare(double a, double b)
{
  return fabs(a-b) < APO_TOLERANCE;
}

double
getAngleDifference(double a1, double a2)
{
    double ret;

    if (a1 >= a2) {
        a2 += 2 * M_PI;
    }
    ret = a2 - a1;

    if (ret >= 2 * M_PI) {
        ret = 0.0;
    }

    return ret;
}

std::vector<APO_Line>
getAngleBisectors(const APO_Line& line1, const APO_Line& line2)
{
  double angle1 = (line1.getAngle() + line2.getAngle()) / 2;
  double angle2 = angle1 + M_PI / 2;

  auto&& ips = getIntersectionPoints(line1, line2, false);
  if (ips.empty())
  {
    return std::vector<APO_Line>();
  }
  auto&& point = ips[0];
  return { APO_Line(point, angle1, 1.0), APO_Line(point, angle2, 1.0) };
}

APO_Point
getCommonIntersectionPoint(const APO_Circle& c1, const APO_Circle& c2,
                           const APO_Circle& c3)
{
  auto&& ips1 = getIntersectionPoints(c1, c2, false);
  auto&& ips2 = getIntersectionPoints(c1, c3, false);

  if (ips1.size() != 2 || ips2.size() != 2)
    return APO_Point(0, 0, false);

  auto p1_1 = ips1[0];
  auto p1_2 = ips1[1];
  auto p2_1 = ips2[0];
  auto p2_2 = ips2[1];

  if (p1_1.equalsFuzzy(p2_1) || p1_1.equalsFuzzy(p2_2))
  {
    return p1_1;
  }
  else if (p1_2.equalsFuzzy(p2_1) || p1_2.equalsFuzzy(p2_2))
  {
    return p1_2;
  }
  else
  {
    return APO_Point(0, 0, false);
  }
}

APO_Point
getPowerCenter(const APO_Circle& c1, const APO_Circle& c2,
               const APO_Circle& c3)
{
  auto getRadicalAxis
      = [](const APO_Circle& c1, const APO_Circle& c2) -> APO_Line {
    if (!c1.center.isValid() || !c2.center.isValid() || std::isnan(c1.radius)
        || !std::isnan(c2.radius))
    {
      return APO_Line(); // Failed, invalid data
    }

    // Define a vector from center to center in 2D, its orientation and length, reject concentric:
    APO_Point v = c2.center- c1.center;
    double d = v.getMagnitude();
    if (d < APO_TOLERANCE)
    {
      return APO_Line(); // Failed, concentric circles
    }
    double dir = v.getAngle();

    // Define distance to center of first circle:
    double d1 = (d + c1.radius / d * c1.radius - c2.radius / d * c2.radius) / 2;

    // Define position on central axis:
    v.setMagnitude(d1);
    v += c1.center;

    // Define an offset vector, construct and return the radical axis:
    APO_Point offset = APO_Point::createPolar(50, dir + M_PI / 2);
    return APO_Line(v - offset, v+offset);
  };

  APO_Line radicalAxis1 = getRadicalAxis(c1, c2);
  APO_Line radicalAxis2 = getRadicalAxis(c1, c3);
  auto&& ips = getIntersectionPoints(radicalAxis1, radicalAxis2, false);
  if (ips.size() == 0)
  {
    return APO_Point(0,0,false);
  }
  return ips[0];
}

std::vector<APO_Line>
getSimilarityAxes(const APO_Circle& c1, const APO_Circle& c2,
                  const APO_Circle& c3)
{
  auto getHomotheticCenters
      = [](const APO_Circle& c1,
           const APO_Circle& c2) -> std::array<APO_Point, 2>
  {
    std::array<APO_Point, 2> ret
        = { APO_Point::undefined, APO_Point::undefined };
    bool ret_ist = false;
    APO_Point pos;
    APO_Point v1, v2;
    APO_Line axis;
    double a = 0.0;
    std::vector<APO_Point> ips;
    // Line connecting the offsets of centers
    APO_Line line;    

    /* # Remark by CVH #
    With some exceptions any pair of circles has 2 well defined homothetic centers

    Exceptions based on: https://en.wikipedia.org/wiki/Homothetic_center#Special_cases
    r1 == 0 && r2 != 0 -> ext = cc1 = int
    r1 != 0 && r2 == 0 -> ext = cc2 = int
    r1 == 0 && r2 == 0 -> ext = Null = int
    r1 == r2 && cc1 != cc2 -> ext = @Inf; int = halfway or Null when r1 == r2 == 0
    r1 == r2 && cc1 == cc2 -> ext = Null; int = cc1=cc2 or Null when r1 == r2 == 0
    r1 != r2 && cc1 != cc2 -> ext = Math; int = Math or both cc1 when r1 == 0 or both cc2 when r2 == 0
    r1 != r2 && cc1 == cc2 -> ext = cc1=cc2 = int or both cc1 when r1 == 0 or both cc2 when r2 == 0

    v1 & v2 are best seen as offsets upwards and/or downwards (Wiki figure 9)
    Both upwards for the external and only one downwards for the internal homothetic center
    Of interest can be that v1 is always considered to be upwards regarding the center of c1
    Rotating the vertical plane around the central axis we can do the math in XY using standard QCAD resources
    */

    // All in 2D, normalized, reject when invalid:
    APO_Point center1 = c1.center;
    APO_Point center2 = c2.center;
    double radius1    = fabs(c1.radius);
    double radius2    = fabs(c2.radius);
    if (!center1.isValid() || !center2.isValid() || std::isnan(radius1)
        || std::isnan(radius1))
    {
      return { APO_Point::invalid,
               APO_Point::invalid }; // Failed, invalid data
    }

    // Define flags:
    bool sameC  = center1.equalsFuzzy(center2);   // APO_TOLERANCE
    bool sameR  = fuzzyCompare(radius1, radius2); // APO_TOLERANCE
    bool zeroR1 = radius1 < APO_TOLERANCE;
    bool zeroR2 = radius2 < APO_TOLERANCE;

    // Force a zero radius if nearly zero:
    if (zeroR1)
    {
      radius1 = 0.0;
    }
    if (zeroR2)
    {
      radius2 = 0.0;
    }

    // Handle limit situations (Math won't solve most of these):
    if (sameR)
    {
      ret[0] = APO_Point::invalid; // Infinite or None when cc1 == cc2
      if (zeroR1 && zeroR2)
      {
        ret[1] = APO_Point::invalid; // None
      }
      else
      {
        ret[1]
            = APO_Point::getAverage(center1, center2); // cc1 or cc2 or halfway
      }
    }
    else
    {
      if (zeroR1)
      {
        ret[0] = center1; // Both cc1
        ret[1] = center1;
      }
      else if (zeroR2)
      {
        ret[0] = center2; // Both cc2
        ret[1] = center2;
      }
      else if (sameC)
      {
        pos    = APO_Point::getAverage(center1, center2); // cc1 or cc2
        ret[0] = pos;
        ret[1] = pos;
      }
    }

    // Shared Math for when one is still undefined at this point:
    if (isUndefined(ret[0]) || isUndefined(ret[1]))
    {
      axis = APO_Line(center1, center2);
      a    = axis.getAngle() + M_PI / 2;
      v1   = APO_Point();
      v1.setPolar(radius1, a);
      v2 = APO_Point();
    }

    // When still undefined get the external homothetic center by an intersection:
    if (isUndefined(ret[0]))
    {
      v2.setPolar(radius2, a);
      line = APO_Line(v1 + center1, v2 + center2);
      ips  = getIntersectionPoints(axis, line, false);
      if (ips.size() == 1 && ips[0].isValid())
      {
        ret[0] = ips[0];
      }
      else
      {
        ret[0] = APO_Point::invalid;
      }
    }

    // When still undefined get the internal homothetic center by an intersection:
    if (isUndefined(ret[1]))
    {
      v2.setPolar(radius2, a + M_PI);
      line = APO_Line(v1 + center1, v2 + center2);
      ips  = getIntersectionPoints(axis, line, false);
      if (ips.size() == 1 && ips[0].isValid())
      {
        ret[1] = ips[0];
      }
      else
      {
        ret[1] = APO_Point::invalid;
      }
    }

    // Return both homothetic centers:
    return ret;
  };

  std::vector<APO_Line> ret;
  // External and internal homothetic centers
  std::vector<APO_Point> hExtC, hIntC;
  double d12, d13, d23; // Intermediate distances

  // Get the homothetic centers for each pair of circles:
  // Both invalid when a shape is not a circle
  auto&& hCenters12 = getHomotheticCenters(c1, c2);
  auto&& hCenters13 = getHomotheticCenters(c1, c3);
  auto&& hCenters23 = getHomotheticCenters(c2, c3);

  // External and internal homothetic centers, in 2D by default:
  hExtC = { hCenters12[0], hCenters13[0], hCenters23[0] };
  hIntC = { hCenters12[1], hCenters13[1], hCenters23[1] };

  // Remove invalid external homothetic centers:
  hExtC.erase(std::remove(hExtC.begin(), hExtC.end(),
                          APO_Point::invalid),
              hExtC.end());

  // One solution is a line connecting all valid external homothetic centers
  // With at least 2 valid external homothetic centers:
  if (hExtC.size() == 2)
  {
    ret.push_back(APO_Line(hExtC[0], hExtC[1]));
  }
  // With all 3 valid external homothetic centers:
  else if (hExtC.size() > 2)
  {
    // Get all 3 intermediate distances:
    d12 = hExtC[0].getDistanceTo(hExtC[1]);
    d13 = hExtC[0].getDistanceTo(hExtC[2]);
    d23 = hExtC[1].getDistanceTo(hExtC[2]);

    // Comparing with NaN is always false
    // Avoid false assumptions:
    if (std::isnan(d12) || std::isnan(d13) || std::isnan(d23))
    {
      return ret; // Failed, at least one distance is invalid
    }

    // Diversify on the longest segment, 1 out 3:
    if (d13 > d12)
    {
      if (fuzzyCompare(d23, d13))
      {
        ret.push_back(
            APO_Line(APO_Point::getAverage(hExtC[0], hExtC[1]), hExtC[2]));
      }
      else if (d23 > d13)
      {
        ret.push_back(APO_Line(hExtC[1], hExtC[2]));
      }
      else
      {
        ret.push_back(APO_Line(hExtC[0], hExtC[2]));
      }
    }
    else
    {
      if (fuzzyCompare(d23, d12))
      {
        ret.push_back(
            APO_Line(APO_Point::getAverage(hExtC[2], hExtC[0]), hExtC[1]));
      }
      else if (d23 > d12)
      {
        ret.push_back(APO_Line(hExtC[1], hExtC[2]));
      }
      else
      {
        ret.push_back(APO_Line(hExtC[0], hExtC[1]));
      }
    }
  }

  // Each valid external homothetic centers is crossed a second time
  // This line emerges from one of the other internal homothetic centers
  // Process each valid external homothetic center:
  for (size_t i = 0; i < hExtC.size(); i++)
  {
    // Get intermediate distances with internal homothetic centers, excluding it's own :
    APO_Point p1 = hIntC[(i + 1) % 3];
    APO_Point p2 = hIntC[(i + 2) % 3];
    double d1    = hExtC[i].getDistanceTo(p1);
    double d2    = hExtC[i].getDistanceTo(p2);

    // Comparing with NaN is always false
    // Avoid false assumptions:
    if (std::isnan(d1) || std::isnan(d2))
    {
      continue; // Failed, at least one distance is invalid
    }

    // Diversify on the longest segment:
    if (fuzzyCompare(d1, d2))
    {
      ret.push_back(APO_Line(hExtC[i], APO_Point::getAverage(p1, p2)));
    }
    else if (d1 > d2)
    {
      ret.push_back(APO_Line(hExtC[i], p1));
    }
    else
    {
      ret.push_back(APO_Line(hExtC[i], p2));
    }
  }

  // Return the collected results (0-4):
  return ret;
}

APO_Point
getPole(const APO_Circle& circle, const APO_Line& polarLine)
{
  auto r      = circle.radius;
  auto center = circle.center;

  APO_Point p = polarLine.getClosestPoint(center, false);
  double op   = center.getDistanceTo(p);
  if (fabs(op) < APO_TOLERANCE)
  {
    return undefined<APO_Point>();
  }

  double opInverse = (r * r) / op;

  APO_Point v = p - center;
  v.setMagnitude(opInverse);
  return center + v;
}

std::vector<APO_Circle>
getSolutionsCCCAlt(const APO_Circle& circle1, const APO_Circle& circle2,
                   const APO_Circle& circle3)
{
  return std::vector<APO_Circle>();
}

std::vector<APO_Circle>
getCircles2TR(const APO_Line shp1, const APO_Line& shp2, double radius,
              const APO_Point& pos)
{
  return std::vector<APO_Circle>();
}

/* ------------------------ APO_Point function impls ------------------------ */

APO_Point APO_Point::invalid = APO_Point(0, 0, false);
APO_Point APO_Point::undefined
    = APO_Point(std::numeric_limits<double>::quiet_NaN(),
                std::numeric_limits<double>::quiet_NaN(), false);

APO_Point::APO_Point() : x(0.0), y(0.0), valid(true) {}

APO_Point::APO_Point(double vx, double vy, bool valid_in)
{
  x     = fabs(vx) < APO_TOLERANCE ? 0.0 : vx;
  y     = fabs(vy) < APO_TOLERANCE ? 0.0 : vy;
  valid = valid_in && APO_IS_NORMALE_NUMBER(x) && APO_IS_NORMALE_NUMBER(y);
}

void
APO_Point::set(double vx, double vy)
{
  x     = fabs(vx) < APO_TOLERANCE ? 0.0 : vx;
  y     = fabs(vy) < APO_TOLERANCE ? 0.0 : vy;
  valid = true;
}

bool
APO_Point::isValid() const
{
  return valid;
}

void
APO_Point::setPolar(double radius, double angle)
{
  x     = radius * cos(angle);
  y     = radius * sin(angle);
  x     = fabs(x) < APO_TOLERANCE ? 0.0 : x;
  y     = fabs(y) < APO_TOLERANCE ? 0.0 : y;
  valid = APO_IS_NORMALE_NUMBER(radius) && APO_IS_NORMALE_NUMBER(angle);
}

void
APO_Point::setAngle(double a)
{
  double m = getMagnitude();
  setPolar(m, a);
}

double
APO_Point::getAngle() const
{
  double ret = 0.0;
  double m   = getMagnitude();

  if (m > 1.0e-6)
  {
    double dp = getDotProduct(*this, APO_Point(1.0, 0.0));
    if (dp / m >= 1.0)
    {
      ret = 0.0;
    }
    else if (dp / m < -1.0)
    {
      ret = M_PI;
    }
    else
    {
      ret = acos(dp / m);
    }
    if (y < 0.0)
    {
      ret = 2 * M_PI - ret;
    }
  }
  return ret;
}

double
APO_Point::getAngleTo(const APO_Point& v) const
{
  if (!valid || !v.valid)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  else
  {
    return (v - *this).getAngle();
  }
}

void
APO_Point::setMagnitude(double m)
{
  double a = getAngle();
  setPolar(m, a);
}

double
APO_Point::getMagnitude() const
{
  if (!valid)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  return sqrt(x * x + y * y);
}

double
APO_Point::getDistanceTo(const APO_Point& v, bool limited) const
{
  if (!valid || !v.valid)
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
  else
  {
    return (*this - v).getMagnitude();
  }
}

bool
APO_Point::equalsFuzzy(const APO_Point& point) const
{
  return fuzzyCompare(x, point.x) && fuzzyCompare(y, point.y);
}

APO_Point
APO_Point::operator+(const APO_Point& v) const
{
  return APO_Point(x + v.x, y + v.y, valid && v.valid);
}

APO_Point
APO_Point::operator-(const APO_Point& v) const
{
  return APO_Point(x - v.x, y - v.y, valid && v.valid);
}

APO_Point
APO_Point::operator*(double s) const
{
  return APO_Point(x * s, y * s, valid);
}

APO_Point
APO_Point::operator/(double s) const
{
  return APO_Point(x / s, y / s, valid);
}

double
APO_Point::getDotProduct(const APO_Point& v1, const APO_Point& v2)
{
  return v1.x * v2.x + v1.y * v2.y;
}

APO_Point&
APO_Point::operator+=(const APO_Point& v)
{
  x += v.x;
  y += v.y;
  valid = valid && v.valid;
  return *this;
}

APO_Point&
APO_Point::operator-=(const APO_Point& v)
{
  x -= v.x;
  y -= v.y;
  valid = valid && v.valid;
  return *this;
}

APO_Point&
APO_Point::operator*=(double s)
{
  x *= s;
  y *= s;
  return *this;
}

APO_Point&
APO_Point::operator/=(double s)
{
  x /= s;
  y /= s;
  return *this;
}

APO_Point
APO_Point::getAverage(const APO_Point& v1, const APO_Point& v2)
{
  return (v1 + v2) / 2.0;
}

double
APO_Point::getCrossProduct(const APO_Point& v1, const APO_Point& v2)
{
  return v1.x * v2.y - v1.y * v2.x;
}

APO_Point
APO_Point::createPolar(double radius, double angle)
{
  APO_Point ret;
  ret.setPolar(radius, angle);
  return ret;
}

std::vector<APO_Point>
APO_Point::getUnique(const std::vector<APO_Point>& vectors)
{

  auto findFirstFuzzy
      = [](const std::vector<APO_Point>& vectors, const APO_Point& v)
  {
    for (int i = 0; i < vectors.size(); i++)
    {
      if (v.equalsFuzzy(vectors[i]))
      {
        return i;
      }
    }

    return -1;
  };

  auto containsFuzzy
      = [&](const std::vector<APO_Point>& vectors, const APO_Point& v)
  { return findFirstFuzzy(vectors, v) != -1; };

  std::vector<APO_Point> ret;
  for (int i = 0; i < vectors.size(); i++)
  {
    if (!containsFuzzy(ret, vectors[i]))
    {
      ret.push_back(vectors[i]);
    }
  }
  return ret;
}

APO_Point&
APO_Point::rotate(double rotation)
{
  if (!valid)
  {
    return *this;
  }

  double r = getMagnitude();
  double a = getAngle() + rotation;

  x = cos(a) * r;
  y = sin(a) * r;

  return *this;
}

APO_Point&
APO_Point::rotate(double rotation, const APO_Point& center)
{
  *this = center + (*this - center).rotate(rotation);
  return *this;
}

/* ------------------------- APO_Line function impls ------------------------ */

APO_Line APO_Line::undefined = APO_Line(APO_Point::undefined, APO_Point::undefined);

APO_Line::APO_Line() : begin_point(0, 0), end_point(0, 0) {}

APO_Line::APO_Line(const APO_Point& begin, const APO_Point& end)
    : begin_point(begin), end_point(end)
{
}

APO_Line::APO_Line(double x1, double y1, double x2, double y2) {}

APO_Line::APO_Line(const APO_Point& begin, double angle, double distance)
{
  end_point = begin_point + APO_Point::createPolar(distance, angle);
}

double
APO_Line::getAngle() const
{
  return begin_point.getAngleTo(end_point);
}

void
APO_Line::setLength(double l, bool fromStart)
{
  if (fromStart)
  {
    end_point = begin_point + APO_Point::createPolar(l, getAngle());
  }
  else
  {
    begin_point = end_point - APO_Point::createPolar(l, getAngle());
  }
}

double
APO_Line::getLength() const
{
  return begin_point.getDistanceTo(end_point);
}

APO_Point
APO_Line::getClosestPoint(const APO_Point& p, bool limited) const
{
  return APO_Point();
}

APO_Point
APO_Line::getMiddlePoint() const
{
  return APO_Point::getAverage(begin_point, end_point);
}

bool
APO_Line::rotate(double rotation, const APO_Point& center)
{
  begin_point.rotate(rotation, center);
  end_point.rotate(rotation, center);
  return true;
}

bool
APO_Line::move(const APO_Point& offset)
{
  if (!offset.isValid() || offset.getMagnitude() < APO_TOLERANCE)
    return false;

  begin_point += offset;
  end_point += offset;
  return true;
}

double
APO_Line::getDistanceTo(const APO_Point& point, bool limited) const
{
  APO_Point v = getVectorTo(point, limited);
  if (v.isValid())
    return v.getMagnitude();
  
  return std::numeric_limits<double>::quiet_NaN();
}

void
APO_Line::reverse()
{
  std::swap(begin_point, end_point);
}

bool
APO_Line::moveTo(const APO_Point& dest)
{
  APO_Point offset = dest - begin_point;
  return move(offset);
}

std::vector<APO_Line>
APO_Line::getOffset(double distance, double num, Side side) const
{
  std::vector<APO_Line> ret;

  std::vector<Side> sides;
  if (side == BothSides)
  {
    sides.push_back(Left);
    sides.push_back(Right);
  }
  else
  {
    sides.push_back(side);
  }

  for (size_t i = 0; i < sides.size(); ++i)
  {
    Side side_ = sides[i];
    double a       = 0.0;
    if (side_ == Left)
    {
      a = begin_point.getAngleTo(end_point) + M_PI / 2.0;
    }
    else
    {
      a = end_point.getAngleTo(begin_point) + M_PI / 2.0;
    }

    APO_Point distanceV;
    for (size_t n = 1; n <= num; ++n)
    {
      distanceV.setPolar(distance * n, a);
      APO_Line parallel = *this;
      parallel.move(distanceV);
      ret.push_back(parallel);
    }
  }
  return ret;
}

APO_Point
APO_Line::getVectorTo(const APO_Point& point, bool limited) const
{
  APO_Point ae = end_point - begin_point;
  APO_Point ap = point - begin_point;

  if (ae.getMagnitude() < 1.0e-6)
  {
    return APO_Point(0, 0, false);
  }

  if (ap.getMagnitude() < 1.0e-6)
  {
    // distance to start point is very small:
    return APO_Point(0, 0, false);
  }

  double b = APO_Point::getDotProduct(ap, ae) / APO_Point::getDotProduct(ae, ae);

  if (limited && (b < 0 || b > 1.0))
  {
    // orthogonal to line does not cross line, use distance to end point:
    APO_Point ret = begin_point.getDistanceTo(point)
                    > end_point.getDistanceTo(point) ? end_point : begin_point;
    return ret;
  }

  APO_Point closestPoint = begin_point + ae * b;

  return point - closestPoint;
}

Side
APO_Line::getSideOfPoint(const APO_Point& point) const
{
  double entityAngle  = getAngle();
  double angleToCoord = begin_point.getAngleTo(point);
  double angleDiff    = getAngleDifference(entityAngle, angleToCoord);

  if (angleDiff < M_PI)
  {
    return Side::Left;
  }
  else
  {
    return Side::Right;
  }
}

/* ------------------------ APO_Circle function impls ----------------------- */

APO_Circle::APO_Circle() : center(APO_Point::invalid), radius(0.0) {}

APO_Circle::APO_Circle(double center_x, double center_y, double radius)
    : center(center_x, center_y), radius(radius)
{
}

APO_Circle::APO_Circle(const APO_Point& center, double radius)
    : center(center), radius(radius)
{
}

APO_Circle
APO_Circle::createFrom3Points(const APO_Point& p1, const APO_Point& p2,
                              const APO_Point& p3)
{
  // intersection of two middle lines

  // middle points between first two points:
  APO_Point mp1 = APO_Point::getAverage(p1, p2);
  double a1 = p1.getAngleTo(p2) + M_PI / 2.0;
  // direction from middle point to center:
  APO_Point dir1 = APO_Point::createPolar(1.0, a1);

  // middle points between last two points:
  APO_Point mp2 = APO_Point::getAverage(p2, p3);
  double a2 = p2.getAngleTo(p3) + M_PI / 2.0;
  // direction from middle point to center:
  APO_Point dir2 = APO_Point::createPolar(1.0, a2);

  APO_Line midLine1(mp1, mp1 + dir1);
  APO_Line midLine2(mp2, mp2 + dir2);

  auto&& ips = getIntersectionPoints(midLine1, midLine2, false);
  if (ips.size() != 1)
    return APO_Circle();

  APO_Point center = ips[0];
  double radius = center.getDistanceTo(p3);

  return APO_Circle(center, radius);
}

APO_Circle
APO_Circle::createFrom2Points(const APO_Point& p1, const APO_Point& p2)
{
  APO_Point center  = (p1 + p2) / 2.0;
  double radius = p1.getDistanceTo(p2) / 2.0;
  return APO_Circle(center, radius);
}

bool
APO_Circle::containsPoint(const APO_Point& p) const
{
  return p.getDistanceTo(center) < radius;
}

/* ------------------------ APO_Shape function impls ----------------------- */

APO_Shape::APO_Shape() { type = APO_UNKOWN_TYPE; }

APO_Shape::APO_Shape(const APO_Point& p)
{
  point = p;
  type  = APO_POINT_TYPE;
}

APO_Shape::APO_Shape(const APO_Line& l)
{
  line = l;
  type = APO_LINE_TYPE;
}

APO_Shape::APO_Shape(const APO_Circle& c)
{
  circle = c;
  type   = APO_CIRCLE_TYPE;
}

/* ------------------------------- capi impls ------------------------------- */

struct apo_object
{
  std::vector<APO_Shape> objects;
};

struct apo_solution
{
  std::vector<APO_Circle> circles;
};

apo_t*
apo_init()
{
  apo_t* apo = new apo_object;
  return apo;
}

void
apo_destroy(apo_t* apo)
{
  delete apo;
}

int
apo_add_point(apo_t* apo, double x, double y)
{
  assert(apo);
  if (apo->objects.size() > 3)
  {
    return false;
  }
  else
  {
    APO_Shape obj(APO_Point(x, y));
    apo->objects.push_back(obj);
    return true;
  }
}

int
apo_add_line(apo_t* apo, double x1, double y1, double x2, double y2)
{
  assert(apo);
  if (apo->objects.size() > 3)
  {
    return false;
  }
  else
  {
    APO_Shape obj(APO_Line(x1, y1, x2, y2));
    apo->objects.push_back(obj);
    return true;
  }
}

int
apo_add_circle(apo_t* apo, double cx, double cy, double r)
{
  assert(apo);
  if (apo->objects.size() > 3)
  {
    return false;
  }
  else
  {
    APO_Shape obj(APO_Circle(cx, cy, r));
    apo->objects.push_back(obj);
    return true;
  }
}

int
apo_solve(apo_t* apo, apo_solution_t** solution)
{
  assert(apo);
  if (apo->objects.size() != 3)
    return false;

  *solution = new apo_solution();

  std::vector<APO_Circle> res
      = getSolution(apo->objects[0], apo->objects[1], apo->objects[2]);
  if (res.empty())
    return false;

  for (auto&& c : res)
  {
    (*solution)->circles.emplace_back(c);
  }

  return true;
}

void
apo_solution_destroy(apo_solution_t* solution)
{
  assert(solution);
  delete solution;
}

unsigned int
apo_solution_get_count(const apo_solution_t* solution)
{
  assert(solution);
  return (unsigned int) solution->circles.size();
}

void
apo_solution_get_circle(const apo_solution_t* solution, unsigned int idx,
                        double* cx, double* cy, double* radius)
{
  assert(solution);
  if (idx < solution->circles.size())
  {
    auto&& c = solution->circles.at(idx);
    *cx      = c.center.x;
    *cy      = c.center.y;
    *radius  = c.radius;
  }
}
