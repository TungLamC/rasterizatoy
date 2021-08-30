#include "../rasterizatoy.hpp"

using namespace rasterizatoy;

struct Vertex { Vector4D position; Vector4D color; };

inline static float tick = 0;

int main()
{
  RectArray<decimal> rect(3, 2);
  for (auto y = 0; y < 2; y++)
  {
    for (auto x = 0; x < 3; x++)
      std::cout << rect[x][y] << " ";
    std::cout << std::endl;
  }
  return 0;
//  std::vector<Vertex> vertices {
//    { { -0.5, -0.5, 0.00, 1}, {255, 0, 0, 1} },
//    { { +0.5, -0.5, 0.00, 1}, {0, 255, 0, 1} },
//    { { +0.0, +0.5, 0.00, 1}, {0, 0, 255, 1} },
//  };
//  Window window(800, 600);
//  using Varying = VARYING_LAYOUT(Vector4D);
//  auto raster = Rasterizater<Varying>(&window);
//  raster.set_vertex_shader([&](uint32_t index, Varying& varying) -> Vector4D {
//    decimal radius = 5;
//    decimal x = std::sin(radians(tick)) * radius;
//    decimal z = std::cos(radians(tick)) * radius;
//    Matrix4D view = look_at({x, 2, z}, {0, 0, 0}, {0, 1, 0});
//    Matrix4D projection = perspective(radians(45.f), (float)800 / (float)600, 10, 100);
//
//    Vector4D& color = std::get<0>(varying);
//    color = vertices[index].color;
//    return projection * view * vertices[index].position;
//  });
//  raster.set_fragment_shader([&](const Varying& varying) -> Vector4D {
//    return std::get<0>(varying);
//  });
//  while (!window.should_close())
//  {
//    ++tick;
//    raster.clear(0, 255, 0);
//    raster.draw_call(3);
//    raster.swap_buffer();
//  }
}