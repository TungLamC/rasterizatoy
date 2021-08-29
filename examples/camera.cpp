#ifndef RASTERIZATOY_CAMERA_CPP_H
#define RASTERIZATOY_CAMERA_CPP_H

#include "../rasterizatoy.hpp"

using namespace rasterizatoy;

inline static float tick = 0;

int main()
{
  struct { Vector4D position; Vector4D color; } vertices[36] = {
    {{-0.5f, -0.5f, -0.5f, +1.0f}, {255, 0, 0, 1}},
    {{+0.5f, +0.5f, -0.5f, +1.0f}, {0, 255, 0, 1}},
    {{+0.5f, -0.5f, -0.5f, +1.0f}, {0, 0, 255, 1}},
    {{+0.5f, +0.5f, -0.5f, +1.0f}, {255, 0, 0, 1}},
    {{-0.5f, -0.5f, -0.5f, +1.0f}, {0, 255, 0, 1}},
    {{-0.5f, +0.5f, -0.5f, +1.0f}, {0, 0, 255, 1}},

    {{-0.5f, -0.5f, +0.5f, +1.0f}, {255, 0, 0, 1}},
    {{+0.5f, -0.5f, +0.5f, +1.0f}, {0, 255, 0, 1}},
    {{+0.5f, +0.5f, +0.5f, +1.0f}, {0, 0, 255, 1}},
    {{+0.5f, +0.5f, +0.5f, +1.0f}, {255, 0, 0, 1}},
    {{-0.5f, +0.5f, +0.5f, +1.0f}, {0, 255, 0, 1}},
    {{-0.5f, -0.5f, +0.5f, +1.0f}, {0, 0, 255, 1}},

    {{-0.5f, +0.5f, +0.5f, +1.0f}, {255, 0, 0, 1}},
    {{-0.5f, +0.5f, -0.5f, +1.0f}, {0, 255, 0, 1}},
    {{-0.5f, -0.5f, -0.5f, +1.0f}, {0, 0, 255, 1}},
    {{-0.5f, -0.5f, -0.5f, +1.0f}, {255, 0, 0, 1}},
    {{-0.5f, -0.5f, +0.5f, +1.0f}, {0, 255, 0, 1}},
    {{-0.5f, +0.5f, +0.5f, +1.0f}, {0, 0, 255, 1}},

    {{+0.5f, +0.5f, +0.5f, +1.0f}, {255, 0, 0, 1}},
    {{+0.5f, -0.5f, -0.5f, +1.0f}, {0, 255, 0, 1}},
    {{+0.5f, +0.5f, -0.5f, +1.0f}, {0, 0, 255, 1}},
    {{+0.5f, -0.5f, -0.5f, +1.0f}, {255, 0, 0, 1}},
    {{+0.5f, +0.5f, +0.5f, +1.0f}, {0, 255, 0, 1}},
    {{+0.5f, -0.5f, +0.5f, +1.0f}, {0, 0, 255, 1}},

    {{-0.5f, -0.5f, -0.5f, +1.0f}, {255, 0, 0, 1}},
    {{+0.5f, -0.5f, -0.5f, +1.0f}, {0, 255, 0, 1}},
    {{+0.5f, -0.5f, +0.5f, +1.0f}, {0, 0, 255, 1}},
    {{+0.5f, -0.5f, +0.5f, +1.0f}, {255, 0, 0, 1}},
    {{-0.5f, -0.5f, +0.5f, +1.0f}, {0, 255, 0, 1}},
    {{-0.5f, -0.5f, -0.5f, +1.0f}, {0, 0, 255, 1}},

    {{-0.5f, +0.5f, -0.5f, +1.0f}, {255, 0, 0, 1}},
    {{+0.5f, +0.5f, +0.5f, +1.0f}, {0, 255, 0, 1}},
    {{+0.5f, +0.5f, -0.5f, +1.0f}, {0, 0, 255, 1}},
    {{+0.5f, +0.5f, +0.5f, +1.0f}, {255, 0, 0, 1}},
    {{-0.5f, +0.5f, -0.5f, +1.0f}, {0, 255, 0, 1}},
    {{-0.5f, +0.5f, +0.5f, +1.0f}, {0, 0, 255, 1}},
  };
  Window window(800, 600);
  using Varying = VARYING_LAYOUT(Vector4D);
  auto rasterizater = Rasterizater<Varying>(&window);
  rasterizater.set_vertex_shader([&](uint32_t index, Varying& varying) -> Vector4D {
    decimal radius = 5;
    decimal x = std::sin(radians(tick)) * radius;
    decimal z = std::cos(radians(tick)) * radius;
    Matrix4D view = look_at({x, 2, z}, {0, 0, 0}, {0, 1, 0});
    Matrix4D projection = perspective(radians(45.f), (float)800 / (float)600, 10, 100);
    std::get<0>(varying) = vertices[index].color;
    return projection * view * vertices[index].position;
  });
  rasterizater.set_fragment_shader([&](const Varying& varying) {
    return std::get<0>(varying);
  });
  while (!window.should_close())
  {
    ++tick;
    rasterizater.clear(0, 0, 0);
    rasterizater.draw_call(36);
    rasterizater.swap_buffer();
  }
  return 0;
}

#endif //RASTERIZATOY_CAMERA_CPP_H