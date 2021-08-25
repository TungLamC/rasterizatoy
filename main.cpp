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
  using V = BufferLayout<integer, decimal>;
  using A = BufferLayout<integer, decimal>;
  DEFINE_BUFFER_LAYOUT(V, A);
  rasterizater::set_vertex_shader([&](std::tuple<integer, decimal>& varyings, const std::tuple<integer, decimal>& attributes) {
    return Vector4D{1.0, 2.0, 3.0, 4.0};
  });
  return 0;
}

