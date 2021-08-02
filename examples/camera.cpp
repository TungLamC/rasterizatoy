#ifndef RASTERIZATOY_CAMERA_CPP_H
#define RASTERIZATOY_CAMERA_CPP_H

#include <chrono>
#include "../rasterizatoy.hpp"

using namespace std;
using namespace std::chrono;
using namespace rasterizatoy;

struct CameraShader: Shader
{
public:
  virtual void vertex_shader(Vertex& vertex) override
  {
    auto now = duration_cast<seconds>(system_clock::now().time_since_epoch());
    real_t radius = 10;
    real_t x = std::sin(now.count()) * radius;
    real_t z = std::cos(now.count()) * radius;
    auto view = look_at<real_t>({x, 0, z}, {0, 0, 0}, {0, 1, 0});
    auto projection = perspective<real_t>(45, 1280 / 720, 0.1, 100);
    vertex.position = projection * view * vertex.position;
  }
};

int main()
{
  Primitive primitive
  {
    {{-0.5, -0.5, +0.0}/*, {255, 255, 255}*/},
    {{+0.5, -0.5, +1.0}/*, {255, 255, 255}*/}, 
    {{+0.0, +0.5, +0.0}/*, {255, 255, 255}*/}
  };

  Window window{1280, 720};
  rasterizater::set_current_context(&window);
  rasterizater::input_primitives(std::vector<Primitive>{primitive});
  rasterizater::set_shader(new CameraShader());

  while (!window.should_close())
  {
    rasterizater::clear({0, 0, 0});
    rasterizater::draw_call();
    rasterizater::swap_buffer();
  }
}

#endif //RASTERIZATOY_CAMERA_CPP_H
