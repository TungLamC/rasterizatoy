#include "../rasterizatoy.hpp"

using namespace rasterizatoy;

int main()
{
  struct { Vector4D pos; Vector4D color; } vs_input[3] = {
    { { -0.5, -0.5, 0.00, 1}, {255, 0, 0, 1} },
    { { +0.5, -0.5, 0.00, 1}, {0, 255, 0, 1} },
    { { +0.0, +0.5, 0.00, 1}, {0, 0, 255, 1} },
  };
  Window window(800, 600);
  using Varying = VARYING_LAYOUT(Vector4D);
  auto raster = Rasterizater<Varying>(&window);
  raster.set_vertex_shader([&](uint32_t index, Varying& varying) -> Vector4D {
    Vector4D& color = std::get<0>(varying);
    color = vs_input[index].color;
    return vs_input[index].pos;
  });
  raster.set_fragment_shader([&](const Varying& varying) -> Vector4D {
    return std::get<0>(varying);
  });
  while (!window.should_close())
  {
    raster.clear(200, 150, 0);
    raster.draw_call(3);
    raster.swap_buffer();
  }
  return 0;
}