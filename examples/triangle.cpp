#include "../rasterizatoy.hpp"

using namespace rasterizatoy;

int main()
{
  Vertex vertex1({-0.5, -0.5, +0.0, +1.0});
  Vertex vertex2({+0.5, -0.5, +0.0, +1.0});
  Vertex vertex3({+0.0, +0.5, +0.0, +1.0});
  Primitive primitive;
  primitive.vertices[0] = vertex1;
  primitive.vertices[1] = vertex2;
  primitive.vertices[2] = vertex3;

  Window window(1280, 720);  
  rasterizater::set_current_context(&window);
  rasterizater::input_primitives(std::vector<Primitive>{primitive});
  rasterizater::set_shader(new Shader());
  while (!window.should_close())
  {
    //    rasterizater::window_->draw_line(640, 540, 960, 180, {255, 0, 0});
    rasterizater::draw_call();
    rasterizater::swap_buffer();
  }
  return 0;
}