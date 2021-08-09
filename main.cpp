#include "rasterizatoy.hpp"

using namespace rasterizatoy;

int main()
{
  Window window(1280, 720);
  while (window.should_close())
    window.swap_buffer();
  return 0;
}