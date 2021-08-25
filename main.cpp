#include <variant>
#include <iostream>

#include "rasterizatoy.hpp"

using namespace std;
using namespace rasterizatoy;

void fuck(std::function<void(int)> f)
{
  int position = 0;
  f(position);
}

int main()
{
  using U = BufferLayout<integer, decimal>;
  using A = BufferLayout<integer, decimal>;
  using V = BufferLayout<integer, decimal>;
  return 0;
}

