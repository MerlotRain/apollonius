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

#ifndef APOLLONIUS_H
#define APOLLONIUS_H

#ifdef __cplusplus
extern "C"
{
#endif

#ifdef APO_BUILD_SHARED
#ifdef _WIN32
#ifdef APO_LIBRARY
#define APO_API __declspec(dllexport)
#else
#define APO_API __declspec(dllimport)
#endif
#else /* Unix */
#define APO_API __attribute__((visibility("default")))
#endif
#else
#define APO_API
#endif

  /**
   * @typedef apo_t
   * @brief Represents the Apollonius problem object.
   * This structure contains the state and data necessary for solving the problem.
   */
  typedef struct apo_object apo_t;

  /**
   * @typedef apo_solution_t
   * @brief Represents the solution to the Apollonius problem, containing the resulting circles.
   */
  typedef struct apo_solution apo_solution_t;

  /**
   * @brief Initializes a new Apollonius problem instance.
   * @return Pointer to a newly allocated Apollonius problem object, or NULL on failure.
   */
  APO_API apo_t* apo_init();

  /**
   * @brief Destroys an Apollonius problem object and frees associated resources.
   * @param apo Pointer to the Apollonius problem object to destroy.
   */
  APO_API void apo_destroy(apo_t* apo);

  /**
   * @brief Adds a point to the Apollonius problem.
   * @param apo Pointer to the Apollonius problem object.
   * @param x X-coordinate of the point.
   * @param y Y-coordinate of the point.
   * @return 1 if the point was added successfully, 0 otherwise.
   */
  APO_API int apo_add_point(apo_t* apo, double x, double y);

  /**
   * @brief Adds a line to the Apollonius problem.
   * @param apo Pointer to the Apollonius problem object.
   * @param x1 X-coordinate of the first point on the line.
   * @param y1 Y-coordinate of the first point on the line.
   * @param x2 X-coordinate of the second point on the line.
   * @param y2 Y-coordinate of the second point on the line.
   * @return 1 if the line was added successfully, 0 otherwise.
   */
  APO_API int apo_add_line(apo_t* apo, double x1, double y1, double x2,
                           double y2);

  /**
   * @brief Adds a circle to the Apollonius problem.
   * @param apo Pointer to the Apollonius problem object.
   * @param cx X-coordinate of the circle center.
   * @param cy Y-coordinate of the circle center.
   * @param r Radius of the circle.
   * @return 1 if the circle was added successfully, 0 otherwise.
   */
  APO_API int apo_add_circle(apo_t* apo, double cx, double cy, double r);

  /**
   * @brief Solves the Apollonius problem and computes the resulting circles.
   * @param apo Pointer to the Apollonius problem object.
   * @param solution Pointer to a variable that will receive the solution object.
   * @return 1 if the solution was computed successfully, 0 otherwise.
   */
  APO_API int apo_solve(apo_t* apo, apo_solution_t** solution);

  /**
   * @brief Destroys an Apollonius solution object and frees associated resources.
   * @param solution Pointer to the Apollonius solution object to destroy.
   */
  APO_API void apo_solution_destroy(apo_solution_t* solution);

  /**
   * @brief Gets the number of circles in the solution.
   * @param solution Pointer to the Apollonius solution object.
   * @return Number of circles in the solution.
   */
  APO_API unsigned int apo_solution_get_count(const apo_solution_t* solution);

  /**
   * @brief Retrieves the details of a circle in the solution.
   * @param solution Pointer to the Apollonius solution object.
   * @param idx Index of the circle (0-based).
   * @param cx Pointer to a variable to receive the X-coordinate of the circle center.
   * @param cy Pointer to a variable to receive the Y-coordinate of the circle center.
   * @param radius Pointer to a variable to receive the radius of the circle.
   */
  APO_API void apo_solution_get_circle(const apo_solution_t* solution,
                                       unsigned int idx, double* cx,
                                       double* cy, double* radius);

#ifdef __cplusplus
}
#endif

#endif