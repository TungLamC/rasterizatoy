#include <variant>
#include <iostream>

#include "rasterizatoy.hpp"

using namespace std;
using namespace rasterizatoy;

typedef uintptr_t gsize;

Matrix4D fuck()
{
  return {
    1, 2, 3, 4,
    5, 6, 7, 8,
    5, 6, 7, 8,
    5, 6, 7, 8,
  };
}

int main()
{
  cout << fuck();
  return 0;
}

