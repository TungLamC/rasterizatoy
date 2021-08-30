#include "../rasterizatoy.hpp"

using namespace rasterizatoy;

static inline decimal tick = 0;

class Texture: public RectArray<Vector4D>
{
public:
  using RectArray<Vector4D>::RectArray;

  inline Vector4D sample2d(decimal u, decimal v)
  {
    u = std::clamp(u, 0.0, 1.0); v = std::clamp(v, 0.0, 1.0);
    decimal x = (width() - 1) * u;
    decimal y = (height() - 1) * v;
    return this->value_[x][y];
  }
};

int main()
{
  Texture texutre(64, 64);
  for (auto y = 0; y < 64; y++)
  {
    for (auto x = 0; x < 64; x++)
    {
      auto i = x / 8;
      auto j = y / 8;
      if ((i + j) % 2 == 0)
        texutre[x][y] = Vector4D{0, 0, 0, 1};
      else
        texutre[x][y] = Vector4D{255, 255, 255, 1};
    }
  }

  struct { Vector4D position; Vector2D texcoord; } vertices[6] = {
    { { -0.5, -0.5, 0.0, 1.0 }, { 0, 0 } },
    { { +0.5, -0.5, 0.0, 1.0 }, { 1, 0 } },
    { { +0.5, +0.5, 0.0, 1.0 }, { 1, 1 } },
//    { { -0.5, -0.5, 0.0, 1.0 }, { 0, 0 } },
//    { { +0.5, +0.5, 0.0, 1.0 }, { 1, 1 } },
//    { { -0.5, +0.5, 0.0, 1.0 }, { 0, 1 } },
  };
  Window window(800, 600);
  using Varying = VARYING_LAYOUT(Vector2D);
  auto rasterizater = Rasterizater<Varying>(&window);
  rasterizater.set_vertex_shader([&](uint32_t index, Varying& varying) {
    decimal radius = 5;
    decimal x = std::sin(radians(tick)) * radius;
    decimal z = std::cos(radians(tick)) * radius;
    Matrix4D view = look_at({x, 2, z}, {0, 0, 0}, {0, 1, 0});
    Matrix4D projection = perspective(radians(45.f), (float)800 / (float)600, 10, 100);
    Vector2D& textcoords = LOCATION(0, varying);
    textcoords = vertices[index].texcoord;
    return projection * view * vertices[index].position;
  });
  rasterizater.set_fragment_shader([&](const Varying& varying) {
    const Vector2D& textcoords = LOCATION(0, varying);
    Vector4D color = texutre.sample2d(textcoords.x, textcoords.y);
    return color;
  });
  while (!window.should_close())
  {
    ++tick;
    rasterizater.clear(0, 0, 0);
    rasterizater.draw_call(6);
    rasterizater.swap_buffer();
  }
  return 0;
}