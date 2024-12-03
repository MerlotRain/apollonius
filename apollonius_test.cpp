#include "apollonius.h"
#include <iostream>

#define APO_INIT                                                              \
  apo_t* apo = apo_init();                                                    \
  if (!apo)                                                                   \
    return;

#define APO_SOLVE                                                             \
  apo_solution* solution = nullptr;                                           \
  apo_solve(apo, &solution);                                                  \
  if (!solution)                                                              \
  {                                                                           \
    apo_destroy(apo);                                                         \
    return;                                                                   \
  }                                                                           \
  unsigned int size = apo_solution_get_count(solution);                       \
  std::cout << "solution size is : " << size << std::endl;                    \
  for (unsigned int i = 0; i < size; ++i)                                     \
  {                                                                           \
    double cx, cy, radius;                                                    \
    apo_solution_get_circle(solution, i, &cx, &cy, &radius);                  \
    std::cout << "circle " << i << " is (" << cx << ", " << cy                \
              << "), radius is " << radius << std::endl;                      \
  }                                                                           \
  apo_solution_destroy(solution);                                             \
  apo_destroy(apo);


int
test_ppp()
{
  {
    APO_INIT
    apo_add_point(apo, 0, 0);
    apo_add_point(apo, 5, 5);
    apo_add_point(apo, 10, 0);
    APO_SOLVE
  }

  {
    APO_INIT
    apo_add_point(apo, 0, 0);
    apo_add_point(apo, 0, 0);
    apo_add_point(apo, 0, 0);
    APO_SOLVE
  }
}

void
test_ppl()
{
  APO_INIT
  apo_add_point(apo, 20, 50);
  apo_add_point(apo, 40, 50);
  apo_add_line(apo, 10, 60, 50, 60);
  APO_SOLVE
}

void
test_ppc()
{
  {
    APO_INIT
    apo_add_point(apo, 30, 60);
    apo_add_point(apo, 50, 50);
    apo_add_circle(apo, 30, 40, 10);
    APO_SOLVE
  }
}

void
test_lll()
{
  APO_INIT
  apo_add_line(apo, 0, 10, 10, 20);
  apo_add_line(apo, 20, 20, 30, 10);
  apo_add_line(apo, 10, 0, 20, 0);
  APO_SOLVE
  //-20.35533906,25,25
  //50.35533906,25,25
  //15,-60.35533906,60.35533906
  //15,10.35533906,10.35533906
}

void
test_llc()
{
  APO_INIT
  apo_add_circle(apo, 30, 10, 10);
  apo_add_line(apo, 10, 20, 20, 30);
  apo_add_line(apo, 40, 30, 50, 20);
  APO_SOLVE
  // 30,28.28427125,8.28427125
  // 30,16.56854249,16.56854249
  // 30,-28.28427125,48.28427125
  // 30,-96.56854249,96.56854249
}

int
main()
{
  test_ppp();
  return 0;
}
  