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

#define APO_POINT_TYPE (0)
#define APO_LINE_TYPE (1)
#define APO_CIRCLE_TYPE (2)

#define APO_IS_POINT(obj) ((obj).type == APO_POINT_TYPE)
#define APO_IS_LINE(obj) ((obj).type == APO_LINE_TYPE)
#define APO_IS_CIRCLE(obj) ((obj).type == APO_CIRCLE_TYPE)

struct APO_Point
{
  double x;
  double y;

  APO_Point();
  APO_Point(double x, double y);
  void setPolar(double radius, double angle);

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
  static APO_Point getAverage(const std::vector<APO_Point>& vectors);

  static double getCrossProduct(const APO_Point& v1, const APO_Point& v2);
  static double getDotProduct(const APO_Point& v1, const APO_Point& v2);
  static APO_Point createPolar(double radius, double angle);
  static std::vector<APO_Point>
  getUnique(const std::vector<APO_Point>& vectors, double tol = APO_TOLERANCE);
};

struct APO_Line
{
  APO_Point begin_point;
  APO_Point end_point;

  APO_Line();
  APO_Line(const APO_Point& begin, const APO_Point& end);
  APO_Line(double x1, double y1, double x2, double y2);
  APO_Line(const APO_Point& begin, double angle, double length);

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

static std::vector<APO_Circle> getSolutionFromPPP(const APO_Object& point1,
                                                  const APO_Object& point2,
                                                  const APO_Object& point3);

static std::vector<APO_Circle> getSolutionFromPPC(const APO_Object& point1,
                                                  const APO_Object& point2,
                                                  const APO_Object& circle);

static std::vector<APO_Circle> getSolutionFromPPL(const APO_Object& point1,
                                                  const APO_Object& point2,
                                                  const APO_Object& line);

static std::vector<APO_Circle> getSolutionFromPCC(const APO_Object& point,
                                                  const APO_Object& circle1,
                                                  const APO_Object& circle2);

static std::vector<APO_Circle> getSolutionFromPLL(const APO_Object& point,
                                                  const APO_Object& line1,
                                                  const APO_Object& line2);

static std::vector<APO_Circle> getSolutionFromPLC(const APO_Object& point,
                                                  const APO_Object& line,
                                                  const APO_Object& circle);

static std::vector<APO_Circle> getSolutionFromLLL(const APO_Object& line1,
                                                  const APO_Object& line2,
                                                  const APO_Object& line3);

static std::vector<APO_Circle> getSolutionFromLLC(const APO_Object& line1,
                                                  const APO_Object& line2,
                                                  const APO_Object& circle);

static std::vector<APO_Circle> getSolutionFromLCC(const APO_Object& line,
                                                  const APO_Object& circle2,
                                                  const APO_Object& circle3);

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
        if (getInverseShape(APO_Object(p), inversionCircle, pinverse))
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
      getInverseShape(APO_Object(APO_Point(circle.center.x + circle.radius,
                                           circle.center.y)),
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
      if (!getInverseShape(APO_Object(p), inversionCircle, pinverse))
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
      if (!getInverseShape(APO_Object(p1), inversionCircle, p1inverse))
      {
        return false;
      }
      if (!getInverseShape(APO_Object(p2), inversionCircle, p2inverse))
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
    if (getInverseShape(APO_Object(shp), inversionCircle, inversed))
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
getSolutionFromPPP(const APO_Object& point1, const APO_Object& point2,
                   const APO_Object& point3)
{
  if (!APO_IS_POINT(point1) || !APO_IS_POINT(point2) || !APO_IS_POINT(point3))
  {
    return std::vector<APO_Circle>();
  }

  return { APO_Circle::createFrom3Points(point1.point, point2.point,
                                         point3.point) };
}

std::vector<APO_Circle>
getSolutionFromPPC(const APO_Object& point1, const APO_Object& point2,
                   const APO_Object& circle)
{
  auto ppc = [](const APO_Point& p1, const APO_Point& p2,
                const APO_Object& circle) -> std::vector<APO_Circle>
  {
    APO_Line l1(circle.circle.center, p1);
    APO_Point m = APO_Point::getAverage(p1, p2);
    APO_Line l2(m, p2.getAngleTo(p1) + M_PI / 2.0, 1.0);

    APO_Point ip;
    if (getIntersectionLL(l1, l2, 0, ip))
    {
      return std::vector<APO_Circle>();
    }
    return { APO_Circle(ip, ip.getDistanceTo(p1)) };
  };

  std::vector<APO_Circle> res;

  if (!APO_IS_POINT(point1) || !APO_IS_POINT(point2) || !APO_IS_CIRCLE(circle))
  {
    return res;
  }

  if (0 == pointIsOnCircle(point1.point, circle.circle)
      && 0 == pointIsOnCircle(point2.point, circle.circle))
  {
    return { APO_Circle(circle.circle.center, circle.circle.radius) };
  }

  APO_Point pOnCircle, pOther;
  if (0 == pointIsOnCircle(point1.point, circle.circle))
  {
    pOnCircle = point1.point;
    pOther    = point2.point;
    return ppc(pOnCircle, pOther, circle);
  }
  if (0 == pointIsOnCircle(point2.point, circle.circle))
  {
    pOnCircle = point2.point;
    pOther    = point1.point;
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
    getInverseShape(APO_Object(lines[i]), inversion_cirlce, res_circle);
    res.push_back(
        APO_Circle(res_circle.circle.center, res_circle.circle.radius));
  }
  return res;
}

std::vector<APO_Circle>
getSolutionFromPPL(const APO_Object& point1, const APO_Object& point2,
                   const APO_Object& line)
{
  auto ppl = [](const APO_Object& point1, const APO_Object& point2,
                const APO_Object& line)
  {
    std::vector<APO_Circle> res;
    APO_Object inversion_cirlce(APO_Circle(point1.point, 10));
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
      getInverseShape(APO_Object(lines[i]), inversion_cirlce, res_circle);
      res.push_back(
          APO_Circle(res_circle.circle.center, res_circle.circle.radius));
    }
    return res;
  };

  if (!APO_IS_POINT(point1) || !APO_IS_POINT(point2) || !APO_IS_LINE(line))
    return std::vector<APO_Circle>();

  bool swapPoint = false;
  if (pointIsOnLine(point1.point, line.line, true))
  {
    swapPoint = true;
    if (pointIsOnLine(point2.point, line.line, true))
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
getSolutionFromPCC(const APO_Object& point, const APO_Object& circle1,
                   const APO_Object& circle2)
{
  std::vector<APO_Circle> res;
  if (!APO_IS_POINT(point) || !APO_IS_CIRCLE(circle1)
      || !APO_IS_CIRCLE(circle2))
    return res;

  // relative sized inversion circle:
  double r_inv = circle1.circle.radius;
  if (circle2.circle.radius > circle1.circle.radius)
  {
    r_inv = circle2.circle.radius;
  }
  APO_Object inversionCircle(APO_Circle(point.point, r_inv));

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
getSolutionFromPLL(const APO_Object& point, const APO_Object& line1,
                   const APO_Object& line2)
{
  std::vector<APO_Circle> res;
  if (!APO_IS_POINT(point) || !APO_IS_LINE(line1) || !APO_IS_LINE(line2))
  {
    return res;
  }

  // intersection between two lines line1, line2:
  APO_Point ipLL;
  bool resLL = getIntersectionLL(line1.line, line2.line, false, ipLL);
  std::vector<APO_Circle> circles;
  std::vector<APO_Point> centers;

  std::vector<APO_Line> bisectorLines
      = getAngleBisectors(line1.line, line2.line);

  bool onLine1 = pointIsOnLine(point.point, line1.line, false);
  bool onLine2 = pointIsOnLine(point.point, line2.line, false);

  bool onBisector = false;
  for (size_t k = 0; k < bisectorLines.size(); k++)
  {
    if (pointIsOnLine(point.point, bisectorLines[k], false))
    {
      onBisector = true;
      break;
    }
  }

  // lines are parallel:
  if (!resLL)
  {
    // middle line:
    APO_Point s         = line1.line.begin_point;
    APO_Point p         = line2.line.getClosestPoint(s, false);
    APO_Point center    = APO_Point::getAverage(s, p);
    APO_Line middleLine = line1.line;
    middleLine.move(center - middleLine.begin_point);
    // circle with radius c-s around point:
    APO_Circle circle = APO_Circle(point.point, center.getDistanceTo(s));
    // intersections between circle and middle line are candidates:
    getIntersectionLC(middleLine, circle, false, centers);
  }

  // point is on line1 or line2:
  else if (onLine1 || onLine2)
  {
    APO_Line line = onLine1 ? line1.line : line2.line;
    APO_Line orthoLine(point.point, line.getAngle() + M_PI / 2, 1.0);
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
      double rp = line1.line.getDistanceTo(point.point, false);
      // distance from intersection line1/line2 to point:
      double dp = ipLL.getDistanceTo(point.point);
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
      circles = getCircles2TR(line1.line, line2.line, 10.0);
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
          line = line1.line;
          line.move(point.point - line1.line.begin_point);
        }
        else
        {
          line = APO_Line(ipLL, point.point);
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
        l1.move(point.point - l1.begin_point);
        l2.move(point.point - l2.begin_point);

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
      double dLine1 = line1.line.getDistanceTo(centerCandidate, false);
      double dLine2 = line1.line.getDistanceTo(centerCandidate, false);
      double dPoint = point.point.getDistanceTo(centerCandidate, false);
      if (fuzzyCompare(dLine1, dPoint) && fuzzyCompare(dLine2, dPoint))
      {
        centers.push_back(centerCandidate);
      }
    }
    centers = APO_Point::getUnique(centers);
  }

  for (size_t c = 0; c < centers.size(); ++c)
  {
    double r = centers[c].getDistanceTo(point.point);
    if (fuzzyCompare(r, 0.0))
    {
      continue;
    }
    res.push_back(APO_Circle(centers[c], r));
  }
  return res;
}

std::vector<APO_Circle>
getSolutionFromPLC(const APO_Object& point, const APO_Object& line,
                   const APO_Object& circle)
{
  if (!APO_IS_POINT(point) || !APO_IS_LINE(line) || !APO_IS_CIRCLE(circle))
    return std::vector<APO_Circle>();

  APO_Point a = point.point;
  APO_Point c = circle.circle.center;
  APO_Point f = line.line.getClosestPoint(c, false);
  APO_Line ortho(c, f);

  APO_Point cf = f - c;
  APO_Point ce = cf;
  ce.setMagnitude(circle.circle.radius);
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
      APO_Line par = line.line;
      par.moveTo(a);

      APO_Point a_;
      if (!getIntersectionLL(da_, par, false, a_))
      {
        continue;
      }

      APO_Point f_;
      if (!getIntersectionLL(da_, line.line, false, f_))
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
      double m      = line.line.getDistanceTo(apm, false);
      APO_Circle circ(a, m);
      APO_Line par2 = line.line;
      par2.moveTo(apm);

      if (getIntersectionLC(par2, circ, false, ips))
      {
        centerCandidates.insert(centerCandidates.end(), ips.begin(),
                                ips.end());
      }
    }
    // special case:
    // point is on line:
    else if (pointIsOnLine(a, line.line, false))
    {
      // similarity axis:
      APO_Line ea(e, a);
      getIntersectionLC(ea, circle.circle, false, ips);
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
      APO_Line orthoA(a, line.line.getAngle() + M_PI / 2, 1.0);

      APO_Point p;
      if (!getIntersectionLL(orthoA, ci1, false, p))
        continue;
      centerCandidates.push_back(p);
    }
    else
    {
      APO_Line da(d, a);
      APO_Point m;
      if (!getIntersectionLL(da, line.line, false, m))
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
        APO_Point mu = line.line.end_point - line.line.begin_point;
        mu.setMagnitude(mt);
        APO_Point u    = m + mu;
        APO_Point v    = m - mu;
        APO_Line orthU = line.line;
        orthU.rotate(M_PI / 2, u);
        APO_Line orthV = line.line;
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
getSolutionFromLLL(const APO_Object& line1, const APO_Object& line2,
                   const APO_Object& line3)
{
  return std::vector<APO_Circle>();
}

std::vector<APO_Circle>
getSolutionFromLLC(const APO_Object& line1, const APO_Object& line2,
                   const APO_Object& circle)
{
  return std::vector<APO_Circle>();
}

std::vector<APO_Circle>
getSolutionFromLCC(const APO_Object& line, const APO_Object& circle2,
                   const APO_Object& circle3)
{
  return std::vector<APO_Circle>();
}

std::vector<APO_Circle>
getSolutionFromCCC(const APO_Object& circle1, const APO_Object& circle2,
                   const APO_Object& circle3)
{
  return std::vector<APO_Circle>();
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

  std::vector<APO_Object> pointObjs;
  std::vector<APO_Object> lineObjs;
  std::vector<APO_Object> circleObjs;

#define OBJECT_CLASSIFICATION(obj)                                            \
  switch (obj.type)                                                           \
  {                                                                           \
  case APO_POINT_TYPE:                                                        \
    pointObjs.push_back(obj);                                                 \
    break;                                                                    \
  case APO_LINE_TYPE:                                                         \
    lineObjs.push_back(obj);                                                  \
    break;                                                                    \
  case APO_CIRCLE_TYPE:                                                       \
    circleObjs.push_back(obj);                                                \
    break;                                                                    \
  };

  OBJECT_CLASSIFICATION(apo->objects[0]);
  OBJECT_CLASSIFICATION(apo->objects[1]);
  OBJECT_CLASSIFICATION(apo->objects[2]);

  std::vector<APO_Circle> res;
  if (pointObjs.size() == 3)
  {

    res = getSolutionFromPPP(pointObjs[0], pointObjs[1], pointObjs[2]);
  }
  else if (pointObjs.size() == 2)
  {
    if (circleObjs.size() == 1)
    {
      res = getSolutionFromPPC(pointObjs[0], pointObjs[1], circleObjs[0]);
    }
    else if (lineObjs.size() == 1)
    {
      res = getSolutionFromPPL(pointObjs[0], pointObjs[1], lineObjs[0]);
    }
  }
  else if (pointObjs.size() == 1)
  {
    if (circleObjs.size() == 2)
    {
      res = getSolutionFromPCC(pointObjs[0], circleObjs[0], circleObjs[1]);
    }
    else if (lineObjs.size() == 2)
    {
      res = getSolutionFromPLL(pointObjs[0], lineObjs[0], lineObjs[1]);
    }
    else if (circleObjs.size() == 1 && lineObjs.size() == 1)
    {
      res = getSolutionFromPLC(pointObjs[0], lineObjs[0], circleObjs[0]);
    }
  }
  else if (pointObjs.size() == 0)
  {
    if (lineObjs.size() == 3)
    {
      res = getSolutionFromLLL(lineObjs[0], lineObjs[1], lineObjs[2]);
    }
    else if (lineObjs.size() == 2 && circleObjs.size() == 1)
    {
      res = getSolutionFromLLC(lineObjs[0], lineObjs[1], circleObjs[0]);
    }
    else if (lineObjs.size() == 1 && circleObjs.size() == 2)
    {
      res = getSolutionFromLCC(lineObjs[0], circleObjs[0], circleObjs[1]);
    }
    else if (circleObjs.size() == 3)
    {
      res = getSolutionFromCCC(circleObjs[0], circleObjs[1], circleObjs[2]);
    }
  }

#undef OBJECT_CLASSIFICATION

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