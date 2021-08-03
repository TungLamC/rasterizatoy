
#ifndef RASTERIZATOY_CAMERA_CPP_H
#define RASTERIZATOY_CAMERA_CPP_H

#include "../rasterizatoy.hpp"

using namespace rasterizatoy;

struct CameraShader: Shader
{
public:
  virtual void vertex_shader(Vertex& vertex) override
  {
    static float tick = 0;
    ++tick;
    real_t radius = 5;
    real_t x = std::sin(radians(tick)) * radius;
    real_t z = std::cos(radians(tick)) * radius;
    auto view = look_at<real_t>({x, 0, z}, {0, 0, 0}, {0, 1, 0});
    auto projection = perspective<real_t>(45, (float)1280 / (float)720, 0.1, 100);
    vertex.position = view * vertex.position;
  }
};

int main()
{
  Primitive primitive{{{-0.5, -0.5, +0.0}}, {{+0.5, -0.5, +0.0}}, {{+0.0, +0.5, +0.0}}};

  Window window{1280, 720};
  rasterizater::set_current_context(&window);
  rasterizater::set_shader(new CameraShader());
  rasterizater::input_primitives(std::vector<Primitive>{primitive});

  while (!window.should_close())
  {
    rasterizater::clear({0, 0, 0});
    rasterizater::draw_call();
    rasterizater::swap_buffer();
  }
}

#endif //RASTERIZATOY_CAMERA_CPP_H