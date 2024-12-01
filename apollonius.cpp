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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define APO_TOLERANCE (1.0e-9)

#define APO_UNKOWN_TYPE (-1)
#define APO_POINT_TYPE (0)
#define APO_LINE_TYPE (1)
#define APO_CIRCLE_TYPE (2)

#define APO_IS_POINT(obj) ((obj).type == APO_POINT_TYPE)
#define APO_IS_LINE(obj) ((obj).type == APO_LINE_TYPE)
#define APO_IS_CIRCLE(obj) ((obj).type == APO_CIRCLE_TYPE)
#define APO_IS_NULL(obj) ((obj).type == APO_UNKOWN_TYPE)

#define NumberIsNormal(x) (!std::isnan(x) && !std::isinf(x))

#define APO_EMPTY_RES std::vector<APO_Circle>()

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

struct APO_Object
{
  int type;
  union
  {
    APO_Point point;
    APO_Line line;
    APO_Circle circle;
  };

  APO_Object();
  APO_Object(const APO_Point& p);
  APO_Object(const APO_Line& l);
  APO_Object(const APO_Circle& c);
};

/* ------------------------- static solve functions ------------------------- */

static bool getInverseShape(const APO_Object& shp,
                            const APO_Object& inversionCircle,
                            APO_Object& inversed);

template <class T>
static std::vector<APO_Circle>
getInverseShapes(const std::vector<T>& shapes,
                 const APO_Object& inversionCircle);

static std::vector<APO_Line>
getCircleTangentsThroughPoint(const APO_Circle& circle, const APO_Point& p);

static std::vector<APO_Line> getAllTangents(const APO_Object& obj1,
                                            const APO_Object& obj2);

static std::vector<APO_Circle>
getSolutions(const std::vector<APO_Object>& objs);

static std::vector<APO_Circle> getSolution(const APO_Object& obj1,
                                           const APO_Object& obj2,
                                           const APO_Object& obj3);

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

static std::vector<APO_Circle> getSolutionFromCCC(const APO_Object& circle1,
                                                  const APO_Object& circle2,
                                                  const APO_Object& circle3);

static bool getIntersectionLC(const APO_Line& line, const APO_Circle& circle,
                              bool limited, std::vector<APO_Point>& points);
static bool getIntersectionCC(const APO_Circle& circle1,
                              const APO_Circle& circle2,
                              std::vector<APO_Point>& points);
static bool getIntersectionLL(const APO_Line& line1, const APO_Line& line2,
                              bool limited, APO_Point& points);
static bool getIntersectsWithLL(const APO_Line& line1, const APO_Line& line2,
                                bool limited);

static bool pointIsOnLine(const APO_Point& point, const APO_Line& line,
                          bool limited);
static bool pointIsOnCircle(const APO_Point& point, const APO_Circle& circle);

static bool fuzzyCompare(double a, double b);
static void sortPointsLRTB(std::vector<APO_Point>& points);
static std::vector<APO_Line> getAngleBisectors(const APO_Line& line1,
                                               const APO_Line& line2);

static std::vector<APO_Circle>
getCircles2TR(const APO_Line shp1, const APO_Line& shp2, double radius,
              const APO_Point& pos = APO_Point());

/* --------------------------------- solves --------------------------------- */

/// http://www.geometer.org/mathcircles/inversion.pdf
bool
getInverseShape(const APO_Object& shp, const APO_Object& inversionCircle,
                APO_Object& inversed)
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
      inversed = APO_Object(shp.point);
      return true;
    }

    double d_inverse = pow(r, 2) / d;
    inversed         = APO_Object(
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
      inversed = APO_Object(shp.line);
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
      APO_Point p;
      if (getIntersectionLL(shp.line, s, false, p))
      {
        APO_Object pinverse;
        if (getInverseShape(p, inversionCircle, pinverse))
        {
          inversed = APO_Object(APO_Circle::createFrom2Points(center, p));
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
      APO_Object inversed_point_obj;
      getInverseShape(
          APO_Point(circle.center.x + circle.radius, circle.center.y),
          inversionCircle, inversed_point_obj);
      double radius = circle.center.x - inversed_point_obj.point.x;
      if (radius < 0)
      {
        radius = fabs(radius);
      }
      inversed = APO_Object(APO_Circle(circle.center, radius));
      return true;
    }
    else if (pointIsOnCircle(inversionCircle.circle.center, circle))
    {
      APO_Line s(inversionCircle.circle.center, circle.center);
      std::vector<APO_Point> ips;
      if (getIntersectionLC(s, circle, false, ips))
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

      APO_Object pinverse;
      if (!getInverseShape(p, inversionCircle, pinverse))
      {
        return false;
      }
      inversed = APO_Object(
          APO_Line(pinverse.point, s.getAngle() + M_PI / 2.0, 1.0));
      return true;
    }
    else
    {
      APO_Line l(inversionCircle.circle.center, circle.center);
      std::vector<APO_Point> ips;
      if (!getIntersectionLC(l, circle, false, ips))
      {
        return false;
      }

      APO_Point p1 = ips[0];
      APO_Point p2 = ips[1];

      APO_Object p1inverse, p2inverse;
      if (!getInverseShape(p1, inversionCircle, p1inverse))
      {
        return false;
      }
      if (!getInverseShape(p2, inversionCircle, p2inverse))
      {
        return false;
      }
      inversed = APO_Object(
          APO_Circle::createFrom2Points(p1inverse.point, p2inverse.point));
      return true;
    }
  }

  return false;
}

template <class T>
std::vector<APO_Circle>
getInverseShapes(const std::vector<T>& shapes,
                 const APO_Object& inversionCircle)
{
  std::vector<APO_Circle> res;
  for (auto&& shp : shapes)
  {
    APO_Object inversed;
    if (getInverseShape(shp, inversionCircle, inversed))
    {
      if (APO_IS_CIRCLE(inversed))
        res.push_back(inversed.circle);
    }
  }
  return res;
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
    std::vector<APO_Point> touching_points;
    getIntersectionCC(circle2, circle, touching_points);
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
getAllTangents(const APO_Object& obj1, const APO_Object& obj2)
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
  if (APO_IS_CIRCLE(obj1) && APO_IS_CIRCLE(obj2))
  {
    // Normalized circle shape clones, dummy to order c1-c2 on size
    APO_Circle c1, c2;
    // Radii of circle shapes
    double c1Radius, c2Radius;
    // Centers of circle shapes
    APO_Point c1Center, c2Center;

    // Validate first circle shape:
    if (std::isnan(obj1.circle.radius))
    {
      return tangents; // Empty, invalid radius
    }
    else if (fabs(obj1.circle.radius) < APO_TOLERANCE)
    {
      // Handle as zero sized 2D circle:
      c1        = obj1.circle;
      c1.radius = 0.0;
    }
    else
    {
      // Handle as normalized 2D circle:
      c1        = obj1.circle;
      c1.radius = fabs(obj1.circle.radius);
    }

    // Validate second circle shape:
    if (std::isnan(obj2.circle.radius))
    {
      return tangents; // Empty, invalid radius
    }
    else if (fabs(obj2.circle.radius) < APO_TOLERANCE)
    {
      // Handle as zero sized 2D circle:
      c2        = obj2.circle;
      c2.radius = 0.0;
    }
    else
    {
      // Handle as normalized 2D circle:
      c2        = obj2.circle;
      c2.radius = fabs(obj2.circle.radius);
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
  else if (APO_IS_LINE(obj1) && APO_IS_LINE(obj2))
  {
    APO_Line l1, l2;
    // Handle as 2D RLine shapes:
    // This would convert RXLine and RRay shapes
    l1 = APO_Line(obj1.line.begin_point, obj1.line.end_point);
    l2 = APO_Line(obj2.line.begin_point, obj2.line.end_point);

    // The angle of a near zero-length line is zero by default (RVector::getAngle())
    // Exclude solutions for a line with almost no length:
    if (l1.getLength() <= 1.0e-6 || l2.getLength() <= 1.0e-6)
    {
      return tangents; // Empty, not processable line(s)
    }

    // Diversify on crossing or not:
    // RLine.isParallel(...) may fail (FS#2495)
    APO_Point ip;
    if (!getIntersectionLL(obj1.line, obj2.line, false, ip))
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
        sortPointsLRTB(ips);
        tangents.push_back(APO_Line(ips[0], ips[3]));
      }
    }
    // With an intersection point:
    else
    {
      a1            = (l1.getAngle() + l2.getAngle()) / 2;
      a2            = a1 + M_PI / 2;
      double length = (l1.getLength() + l2.getLength()) / 2;
      tangents.push_back(APO_Line(ip, a1, length));
      tangents.push_back(APO_Line(ip, a2, length));
    }
    return tangents; // No or 1-2 solution(s)
  } // End line-line

  APO_Object c1, c2;
  APO_Object item1, item2;

  // With a circle and a line shape:
  if (APO_IS_CIRCLE(obj1) && APO_IS_LINE(obj2))
  {

    // Handle second shape as 2D RLine:
    // This would convert RXLine and RRay shapes
    line = APO_Line(obj2.line.begin_point, obj2.line.end_point);

    // The angle of a near zero-length line is zero by default (RVector::getAngle())
    // Exclude solutions for a line with almost no length:
    if (line.getLength() <= 1.0e-6)
    {
      return tangents; // Empty, not processable line
    }

    // Handle first shape as normalized 2D circle:
    c2               = obj1;
    c2.circle.radius = fabs(obj1.circle.radius);

  } // End circle-line

  // With a line and a circle shape:
  else if (APO_IS_LINE(obj1) && APO_IS_CIRCLE(obj2))
  {

    // Handle first shape as 2D RLine:
    // This would convert RXLine and RRay shapes
    line = obj1.line;

    // The angle of a near zero-length line is zero by default (RVector::getAngle())
    // Exclude solutions for a line with almost no length:
    if (line.getLength() <= 1.0e-6)
    {
      return tangents; // Empty, not processable line
    }

    // Handle second shape as normalized 2D circle:
    c2               = obj2;
    c2.circle.radius = fabs(obj2.circle.radius);

  } // End line-circle

  // Least common usage: Support RPoint shapes as zero sized RCircle shapes:
  // Circle shapes are not guaranteed to be valid circles
  // Shapes are further validated and handled by a recursive call
  else if (APO_IS_POINT(obj1) || APO_IS_POINT(obj2))
  {
    // Ifso, convert first point into an RCircle:
    if (APO_IS_POINT(obj1))
    {
      item1 = APO_Object(APO_Circle(obj1.point, 0.0));
    }
    else
    {
      item1 = obj1;
    }

    // Ifso, convert second point into an RCircle:
    if (APO_IS_POINT(obj2))
    {
      item2 = APO_Object(APO_Circle(obj2.point, 0.0));
    }
    else
    {
      item2 = obj2;
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
    getIntersectionLC(line, c2.circle, false, ips); // unlimited
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
getSolutions(const std::vector<APO_Object>& objs)
{
  assert(objs.size() == 3);

  return getSolution(objs[0], objs[1], objs[2]);
}

std::vector<APO_Circle>
getSolution(const APO_Object& obj1, const APO_Object& obj2,
            const APO_Object& obj3)
{

  std::vector<APO_Point> pointObjs;
  std::vector<APO_Line> lineObjs;
  std::vector<APO_Circle> circleObjs;

#define OBJECT_CLASSIFICATION(obj)                                            \
  switch (obj.type)                                                           \
  {                                                                           \
  case APO_POINT_TYPE:                                                        \
    pointObjs.push_back(obj.point);                                           \
    break;                                                                    \
  case APO_LINE_TYPE:                                                         \
    lineObjs.push_back(obj.line);                                             \
    break;                                                                    \
  case APO_CIRCLE_TYPE:                                                       \
    circleObjs.push_back(obj.circle);                                         \
    break;                                                                    \
  };

  OBJECT_CLASSIFICATION(obj1);
  OBJECT_CLASSIFICATION(obj2);
  OBJECT_CLASSIFICATION(obj3);

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
  return APO_EMPTY_RES;
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

    APO_Point ip;
    if (getIntersectionLL(l1, l2, 0, ip))
    {
      return APO_EMPTY_RES;
    }
    return { APO_Circle(ip, ip.getDistanceTo(p1)) };
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

  APO_Object inversion_cirlce;
  inversion_cirlce.type          = APO_CIRCLE_TYPE;
  inversion_cirlce.circle.center = pOnCircle;
  inversion_cirlce.circle.radius = 10.0;
  APO_Object circle_inverse;
  APO_Object point2_inverse;
  if (!getInverseShape(circle, inversion_cirlce, circle_inverse)
      || !getInverseShape(point2, inversion_cirlce, point2_inverse))
  {
    return res;
  }

  std::vector<APO_Line> lines = getCircleTangentsThroughPoint(
      circle_inverse.circle, point2_inverse.point);
  for (int i = 0; i < lines.size(); ++i)
  {
    APO_Object res_circle;
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
    APO_Object inversion_cirlce(APO_Circle(point1, 10));
    APO_Object line_inverse;
    APO_Object point2_inverse;
    if (!getInverseShape(line, inversion_cirlce, line_inverse)
        || !getInverseShape(point2, inversion_cirlce, point2_inverse))
    {
      return res;
    }

    std::vector<APO_Line> lines = getCircleTangentsThroughPoint(
        line_inverse.circle, point2_inverse.point);
    for (int i = 0; i < lines.size(); ++i)
    {
      APO_Object res_circle;
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
      return APO_EMPTY_RES;
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
  APO_Object inversionCircle(APO_Circle(point, r_inv));

  // construct inversion shape:
  APO_Object c1Inverse;
  APO_Object c2Inverse;
  if (!getInverseShape(circle1, inversionCircle, c1Inverse)
      || !getInverseShape(circle2, inversionCircle, c2Inverse))
  {
    return res;
  }

  // Get all tangent shapes for given inversion shapes:
  // Exploits an enhanced algorithm by CVH
  auto&& tangents = getAllTangents(c1Inverse, c2Inverse);

  // Return the re-inversion of all tangents (0-4 solutions):
  return getInverseShapes<APO_Line>(tangents, inversionCircle);
}

std::vector<APO_Circle>
getSolutionFromPLL(const APO_Point& point, const APO_Line& line1,
                   const APO_Line& line2)
{
  std::vector<APO_Circle> res;

  // intersection between two lines line1, line2:
  APO_Point ipLL;
  bool resLL = getIntersectionLL(line1, line2, false, ipLL);
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
  if (!resLL)
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
    getIntersectionLC(middleLine, circle, false, centers);
  }

  // point is on line1 or line2:
  else if (onLine1 || onLine2)
  {
    APO_Line line = onLine1 ? line1 : line2;
    APO_Line orthoLine(point, line.getAngle() + M_PI / 2, 1.0);
    for (size_t k = 0; k < bisectorLines.size(); k++)
    {
      auto&& bisectorLine = bisectorLines[k];
      if (getIntersectionLL(bisectorLine, orthoLine, false, ipLL))
        centers.push_back(ipLL);
    }
  }

  else
  {
    std::vector<APO_Point> centerCandidates;

    // point on bisector:
    if (onBisector)
    {
      if (!resLL)
      {
        return res;
      }
      // distance from point to line1 (radius of circle around point, tangential to lines):
      double rp = line1.getDistanceTo(point, false);
      // distance from intersection line1/line2 to point:
      double dp = ipLL.getDistanceTo(point);
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
        centerCandidates.push_back(ipLL
                                   + (APO_Point::createPolar(dp + r2, a)));
        centerCandidates.push_back(ipLL
                                   + (APO_Point::createPolar(dp - r1, a)));
        centerCandidates.push_back(
            ipLL + (APO_Point::createPolar(dp + r2, a + M_PI)));
        centerCandidates.push_back(
            ipLL + (APO_Point::createPolar(dp - r1, a + M_PI)));
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
        if (!resLL)
        {
          // lines parallel:
          line = line1;
          line.move(point - line1.begin_point);
        }
        else
        {
          line = APO_Line(ipLL, point);
        }

        // intersections between line L and circle C -> G, H:
        std::vector<APO_Point> ipsLC;
        getIntersectionLC(line, circle, false, ipsLC);
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
          APO_Point tmp;
          if (getIntersectionLL(l1, bisectorLine, false, tmp))
            centerCandidates.push_back(tmp);
          if (getIntersectionLL(l2, bisectorLine, false, tmp))
            centerCandidates.push_back(tmp);
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

      APO_Point a_;
      if (!getIntersectionLL(da_, par, false, a_))
      {
        continue;
      }

      APO_Point f_;
      if (!getIntersectionLL(da_, line, false, f_))
      {
        continue;
      }

      APO_Line a_e(a_, e);
      APO_Line f_p = a_e;
      f_p.moveTo(f_);

      APO_Point p;
      if (!getIntersectionLL(f_p, ortho, false, p))
      {
        continue;
      }

      p = ips[0];

      APO_Point apm = APO_Point::getAverage(a, p);
      double m      = line.getDistanceTo(apm, false);
      APO_Circle circ(a, m);
      APO_Line par2 = line;
      par2.moveTo(apm);

      if (getIntersectionLC(par2, circ, false, ips))
      {
        centerCandidates.insert(centerCandidates.end(), ips.begin(),
                                ips.end());
      }
    }
    // special case:
    // point is on line:
    else if (pointIsOnLine(a, line, false))
    {
      // similarity axis:
      APO_Line ea(e, a);
      getIntersectionLC(ea, circle, false, ips);
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

      APO_Point p;
      if (!getIntersectionLL(orthoA, ci1, false, p))
        continue;
      centerCandidates.push_back(p);
    }
    else
    {
      APO_Line da(d, a);
      APO_Point m;
      if (!getIntersectionLL(da, line, false, m))
        continue;

      APO_Circle efa = APO_Circle::createFrom3Points(e, f, a);
      std::vector<APO_Point> ps;
      if (!getIntersectionLC(da, efa, false, ps))
      {
        continue;
      }

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

        APO_Point pps;
        if (getIntersectionLL(orthU, tangent, false, pps))
        {
          centerCandidates.push_back(pps);
        }
        if (getIntersectionLL(orthV, tangent, false, pps))
        {
          centerCandidates.push_back(pps);
        }
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
  if (getIntersectsWithLL(line1, line2, false))
  {
    situation += 1;
  }
  if (getIntersectsWithLL(line1, line3, false))
  {
    situation += 2;
  }
  if (getIntersectsWithLL(line2, line3, false))
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
  if (getIntersectionLL(l1, l2, false, p))                                    \
    centerPoints.push_back(p);

  APO_Point p;
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
  double radius;
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

  std::vector<APO_Object> arr1;
  std::vector<APO_Object> arr2;
  std::vector<APO_Object> arr3;
  std::vector<APO_Object> arr4;

  arr1.push_back(APO_Object(circle.center));
  arr2.push_back(APO_Object(circle.center));
  arr3.push_back(APO_Object(circle.center));
  arr4.push_back(APO_Object(circle.center));

  arr1.push_back(APO_Object(parallels1[0]));
  arr2.push_back(APO_Object(parallels1[0]));
  arr3.push_back(APO_Object(parallels1[1]));
  arr4.push_back(APO_Object(parallels1[1]));

  arr1.push_back(APO_Object(parallels2[0]));
  arr2.push_back(APO_Object(parallels2[1]));
  arr3.push_back(APO_Object(parallels2[0]));
  arr4.push_back(APO_Object(parallels2[1]));

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
    std::vector<APO_Object> arr1;
    std::vector<APO_Object> arr2;
    std::vector<APO_Object> arr3;
    std::vector<APO_Object> arr4;

    arr1.push_back(APO_Object(circle1.center));
    arr2.push_back(APO_Object(circle1.center));
    arr3.push_back(APO_Object(circle1.center));
    arr4.push_back(APO_Object(circle1.center));

    APO_Circle circle21 = circle2;
    circle21.radius += circle1.radius;
    arr1.push_back(circle21);
    arr3.push_back(circle21);

    if (fuzzyCompare(circle1.radius, circle2.radius))
    {
      arr2.push_back(APO_Object(circle2.center));
      arr3.push_back(APO_Object(circle2.center));
    }
    else
    {
      APO_Object circle22    = circle2;
      circle22.circle.radius = fabs(circle22.circle.radius - circle1.radius);
      arr2.push_back(circle22);
      arr4.push_back(circle22);
    }

    std::vector<APO_Line> parallels
        = line.getOffset(circle1.radius, 1, BothSides);
    arr1.push_back(APO_Object(parallels[0]));
    arr2.push_back(APO_Object(parallels[0]));
    arr3.push_back(APO_Object(parallels[1]));
    arr4.push_back(APO_Object(parallels[1]));

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
getSolutionFromCCC(const APO_Object& c1, const APO_Object& c2,
                   const APO_Object& c3)
{
  return APO_EMPTY_RES;
}

bool
getIntersectionLC(const APO_Line& line, const APO_Circle& circle, bool limited,
                  std::vector<APO_Point>& points)
{
  return false;
}

bool
getIntersectionCC(const APO_Circle& circle1, const APO_Circle& circle2,
                  std::vector<APO_Point>& points)
{
  return false;
}

bool
getIntersectionLL(const APO_Line& line1, const APO_Line& line2, bool limited,
                  APO_Point& points)
{
  return false;
}

bool
pointIsOnLine(const APO_Point& point, const APO_Line& line, bool limited)
{
  return false;
}
bool
pointIsOnCircle(const APO_Point& point, const APO_Circle& circle)
{
  return false;
}

bool
fuzzyCompare(double a, double b)
{
  return false;
}

void
sortPointsLRTB(std::vector<APO_Point>& points)
{
}

/* ------------------------ APO_Point function impls ------------------------ */

APO_Point::APO_Point() : x(0.0), y(0.0), valid(true) {}

APO_Point::APO_Point(double vx, double vy, bool valid_in)
{
  x     = fabs(vx) < APO_TOLERANCE ? 0.0 : vx;
  y     = fabs(vy) < APO_TOLERANCE ? 0.0 : vy;
  valid = valid_in && NumberIsNormal(x) && NumberIsNormal(y);
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
  valid = NumberIsNormal(radius) && NumberIsNormal(angle);
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
APO_Point::getDistanceTo(const APO_Point& v, bool limited = true) const
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
  return false;
}

double
APO_Line::getDistanceTo(const APO_Point& point, bool limited) const
{
  return 0.0;
}

void
APO_Line::reverse()
{
  std::swap(begin_point, end_point);
}

bool
APO_Line::moveTo(const APO_Point& dest)
{
  return false;
}

std::vector<APO_Line>
APO_Line::getOffset(double distance, double num, Side side) const
{
  return std::vector<APO_Line>();
}

/* ------------------------ APO_Circle function impls ----------------------- */

APO_Circle::APO_Circle() : center(0.0, 0.0), radius(0.0) {}

APO_Circle::APO_Circle(double center_x, double center_y, double radius)
{
  center = APO_Point(center_x, center_y);
  radius = radius;
}

APO_Circle::APO_Circle(const APO_Point& center, double radius)
    : center(center), radius(radius)
{
}

APO_Circle
APO_Circle::createFrom3Points(const APO_Point& p1, const APO_Point& p2,
                              const APO_Point& p3)
{
  return APO_Circle();
}

APO_Circle
APO_Circle::createFrom2Points(const APO_Point& p1, const APO_Point& p2)
{
  return APO_Circle();
}

bool
APO_Circle::containsPoint(const APO_Point& p) const
{
  return false;
}

/* ------------------------ APO_Object function impls ----------------------- */

APO_Object::APO_Object() { type = APO_UNKOWN_TYPE; }

APO_Object::APO_Object(const APO_Point& p)
{
  point = p;
  type  = APO_POINT_TYPE;
}

APO_Object::APO_Object(const APO_Line& l)
{
  line = l;
  type = APO_LINE_TYPE;
}

APO_Object::APO_Object(const APO_Circle& c)
{
  circle = c;
  type   = APO_CIRCLE_TYPE;
}

/* ------------------------------- capi impls ------------------------------- */

struct apo_object
{
  std::vector<APO_Object> objects;
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
    APO_Object obj(APO_Point(x, y));
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
    APO_Object obj(APO_Line(x1, y1, x2, y2));
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
    APO_Object obj(APO_Circle(cx, cy, r));
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
  return solution->circles.size();
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