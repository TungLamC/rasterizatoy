#include "../rasterizatoy.hpp"

using namespace rasterizatoy;

int main()
{
  Primitive primitive{
    {{-0.5, -0.5, +0.0}, {255, 255, 255}}, {{+0.5, -0.5, +1.0}, {255, 255, 255}}, {{+0.0, +0.5, +0.0}, {255, 255, 255}}
  };

  Window window(1280, 720);
  rasterizater::set_current_context(&window);
  rasterizater::input_primitives(std::vector<Primitive>{primitive});
  rasterizater::set_shader(new Shader());
  while (!window.should_close())
  {
    rasterizater::clear({0, 0, 0});
    std::cout << "------------------------" << std::endl;
    rasterizater::draw_call();
    std::cout << "++++++++++++++++++++++++" << std::endl;
    rasterizater::swap_buffer();
  }
  return 0;
}