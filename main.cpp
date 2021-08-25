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
//  int a = 2;
//  std::function<void(int)> f = [&](auto b) { std::cout << a << ", " << b << std::endl; };
//  using U = BufferLayout<integer, decimal>;
//  using A = BufferLayout<integer, decimal>;
//  using V = BufferLayout<integer, decimal>;
  DEFINE_SHADER_LAYOUT(U, A, V);
  return 0;
}

