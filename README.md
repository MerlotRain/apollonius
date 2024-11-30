# Apollonius Circle Problem Library

This library provides a C/C++ API to solve the Apollonius Circle Problem, which involves finding the circle tangent to three given geometric objects (points, lines, or circles). It includes functions for adding geometric objects to the problem, solving for the resulting circles, and retrieving the details of the solution.


## Table of Contents

- [Apollonius Circle Problem Library](#apollonius-circle-problem-library)
  - [Table of Contents](#table-of-contents)
  - [Installation](#installation)
      - [Windows Specific Installation](#windows-specific-installation)
  - [Usage](#usage)
      - [Initialization](#initialization)
      - [Adding Points, Lines, and Circles](#adding-points-lines-and-circles)
      - [Solving the Problem](#solving-the-problem)
      - [Retrieving the Solution](#retrieving-the-solution)
      - [Cleaning Up](#cleaning-up)
  - [API Reference](#api-reference)
        - [apo_t* apo_init()](#apo_t-apo_init)
        - [void apo_destroy(apo_t* apo)](#void-apo_destroyapo_t-apo)
        - [int apo_add_point(apo_t* apo, double x, double y)](#int-apo_add_pointapo_t-apo-double-x-double-y)
        - [int apo_add_line(apo_t* apo, double x1, double y1, double x2, double y2)](#int-apo_add_lineapo_t-apo-double-x1-double-y1-double-x2-double-y2)
        - [int apo_add_circle(apo_t* apo, double cx, double cy, double r)](#int-apo_add_circleapo_t-apo-double-cx-double-cy-double-r)
        - [int apo_solve(apo_t* apo, apo_solution_t** solution)](#int-apo_solveapo_t-apo-apo_solution_t-solution)
        - [void apo_solution_destroy(apo_solution_t* solution)](#void-apo_solution_destroyapo_solution_t-solution)
        - [unsigned int apo_solution_get_count(const apo_solution_t* solution)](#unsigned-int-apo_solution_get_countconst-apo_solution_t-solution)
        - [void apo_solution_get_circle(const apo_solution_t* solution, unsigned int idx, double* cx, double* cy, double* radius)](#void-apo_solution_get_circleconst-apo_solution_t-solution-unsigned-int-idx-double-cx-double-cy-double-radius)
  - [License](#license)

## Installation

1. Clone the repository:

``` bash
git clone https://github.com/yourusername/apollonius.git
cd apollonius
```
2. Build the library using a CMake build system (optional for shared/static library):

``` bash
mkdir build
cd build
cmake ..
make
```

3. (Optional) Install the library:

``` bash
sudo make install
```

#### Windows Specific Installation

If you are on Windows and building as a shared library, ensure that you have Visual Studio installed and use the appropriate CMake generator.

## Usage

#### Initialization

To begin using the Apollonius library, initialize a new Apollonius problem instance using apo_init():

``` c
apo_t* apo = apo_init();
if (!apo) {
    // Handle error
}
``` 
#### Adding Points, Lines, and Circles

You can add points, lines, and circles to the problem. The library allows adding:

   * Points: Represented by (x, y) coordinates.

   * Lines: Defined by two points, (x1, y1) and (x2, y2).

   * Circles: Defined by a center (cx, cy) and a radius r.

Example for adding a point, a line, and a circle:

``` c
apo_add_point(apo, 0.0, 0.0);  // Add a point at (0, 0)
apo_add_line(apo, 0.0, 0.0, 1.0, 1.0);  // Add a line from (0, 0) to (1, 1)
apo_add_circle(apo, 2.0, 2.0, 1.0);  // Add a circle with center (2, 2) and radius 1
```
#### Solving the Problem

Once you have added the necessary objects (points, lines, and circles), you can solve for the resulting Apollonius circles using apo_solve():

```c
apo_solution_t* solution = NULL;
int result = apo_solve(apo, &solution);
if (result == 0 || solution == NULL) {
    // Handle error
}
```

#### Retrieving the Solution

After solving the problem, you can retrieve the resulting circles. Use the following functions to get the count of circles and details of each circle (center coordinates and radius):

```c
unsigned int count = apo_solution_get_count(solution);
for (unsigned int i = 0; i < count; i++) {
    double cx, cy, radius;
    apo_solution_get_circle(solution, i, &cx, &cy, &radius);
    // Do something with the circle data
}
```
#### Cleaning Up

Don't forget to free the memory for the Apollonius problem and solution objects when done:

```c
apo_destroy(apo);  // Destroy the Apollonius problem object
apo_solution_destroy(solution);  // Destroy the solution object
```

## API Reference

##### apo_t* apo_init()

Initializes a new Apollonius problem instance.
   
   * Returns: Pointer to a new apo_t object or NULL if allocation fails.

##### void apo_destroy(apo_t* apo)

Destroys an Apollonius problem object and frees associated resources.

   * Parameters: apo — Pointer to the Apollonius problem object to destroy.

##### int apo_add_point(apo_t* apo, double x, double y)

Adds a point to the Apollonius problem.

   * Parameters:

     * apo: Pointer to the Apollonius problem object.

     * x: X-coordinate of the point.

     * y: Y-coordinate of the point.

     * Returns: 1 if the point was added successfully, 0 otherwise.

##### int apo_add_line(apo_t* apo, double x1, double y1, double x2, double y2)

Adds a line to the Apollonius problem.

   * Parameters:

     * apo: Pointer to the Apollonius problem object.

     * x1, y1: Coordinates of the first point of the line.

     * x2, y2: Coordinates of the second point of the line.

     * Returns: 1 if the line was added successfully, 0 otherwise.

##### int apo_add_circle(apo_t* apo, double cx, double cy, double r)

Adds a circle to the Apollonius problem.

   * Parameters:

     * apo: Pointer to the Apollonius problem object.

     * cx, cy: Coordinates of the circle's center.

     * r: Radius of the circle.

     * Returns: 1 if the circle was added successfully, 0 otherwise.

##### int apo_solve(apo_t* apo, apo_solution_t** solution)

Solves the Apollonius problem and computes the resulting circles.

   * Parameters:

     * apo: Pointer to the Apollonius problem object.

     * solution: Pointer to a variable that will receive the solution object.

     * Returns: 1 if the solution was computed successfully, 0 otherwise.

##### void apo_solution_destroy(apo_solution_t* solution)

Destroys an Apollonius solution object and frees associated resources.

   * Parameters: solution — Pointer to the solution object to destroy.

##### unsigned int apo_solution_get_count(const apo_solution_t* solution)

Gets the number of circles in the solution.

   * Parameters: solution — Pointer to the Apollonius solution object.

     * Returns: Number of circles in the solution.

##### void apo_solution_get_circle(const apo_solution_t* solution, unsigned int idx, double* cx, double* cy, double* radius)

Retrieves the details of a circle in the solution.

   * Parameters:

     * solution: Pointer to the Apollonius solution object.

     * idx: Index of the circle (0-based).

     * cx, cy: Pointers to variables that will receive the circle's center coordinates.

     * radius: Pointer to a variable that will receive the circle's radius.

## License

This project is licensed under the MIT License - see the LICENSE file for details.
