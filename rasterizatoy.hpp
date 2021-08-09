#ifndef RASTERIZATOY_HPP
#define RASTERIZATOY_HPP

#include "cimg.h"

#include <map>
#include <cmath>
#include <vector>
#include <string>
#include <variant>
#include <cassert>
#include <type_traits>

namespace rasterizatoy
{
using namespace cimg_library;

using integer = int;
using decimal = float;

template<size_t N, typename T> struct Vector;
template<size_t ROW, size_t COLUMN, typename T> class Matrix;
template<typename T> using Vector2 = Vector<2, T>;
template<typename T> using Vector3 = Vector<3, T>;
template<typename T> using Vector4 = Vector<4, T>;
template<typename T> using Matrix4 = Matrix<4, 4, T>;

using Vector2I = Vector<2, integer>;
using Vector3I = Vector<3, integer>;
using Vector4I = Vector<4, integer>;
using Vector2D = Vector<2, decimal>;
using Vector3D = Vector<3, decimal>;
using Vector4D = Vector<4, decimal>;
using Matrix4D = Matrix<4, 4, decimal>;

using ShaderPass = std::variant<Vector2D, Vector3D, Vector4D, Matrix4D>;

using RGB  = Vector3<uint8_t>;
using RGBA = Vector4<uint8_t>;
}

namespace rasterizatoy
{
namespace derived
{
template<size_t DIMENSION, typename T> struct VectorN { T components[DIMENSION]; };
template<typename T> struct Vector2 { union { struct { T x, y; }; struct { T u, v; }; T components[2]; }; };
template<typename T> struct Vector3 { union { struct { T x, y, z; }; struct { T r, g, b; }; T components[3]; }; };
template<typename T> struct Vector4 { union { struct { T x, y, z, w; }; struct { T r, g, b, a; }; T components[4]; }; };
}

template<size_t DIMENSION, typename T>
struct Vector:
  std::conditional_t<DIMENSION == 2, derived::Vector2<T>,
  std::conditional_t<DIMENSION == 3, derived::Vector3<T>,
  std::conditional_t<DIMENSION == 4, derived::Vector4<T>,
  derived::VectorN<DIMENSION, T>>>>
{
  inline Vector() = default;
  inline Vector(const T* pointer) { for (size_t i = 0; i < DIMENSION; i++) this->components[i] = pointer[i]; }
  inline Vector(const Vector<DIMENSION, T>& other) { for (size_t i = 0; i < DIMENSION; i++) this->components[i] = other.components[i]; }
  inline Vector(const std::initializer_list<T>& initializer) { for (size_t i = 0; i < DIMENSION; i++) this->components[i] = *(initializer.begin() + i); }

  inline T& operator[](size_t index) { assert(index < DIMENSION); return this->components[index]; }
  inline const T& operator[](size_t index) const { assert(index < DIMENSION); return this->components[index]; }

  T components[DIMENSION];
};

template<size_t ROW, size_t COLUMN, typename T>
struct Matrix
{
};
}

namespace rasterizatoy
{
struct Vertex
{
  inline Vertex(const Vector4D& position = {0, 0, 0, 1})
    : rhw{position.w}, position{position}, viewport{0, 0} { }

  decimal rhw;
  Vector4D position;
  Vector2D viewport;
};

struct Primitive
{
  Vertex vertices[3];
};

struct Shader
{
  virtual ShaderPass vertex_shader(const Vertex& input) = 0;
  virtual RGBA fragment_shader(const ShaderPass& input) = 0;
};

class Window
{
public:
  inline Window(uint32_t width, uint32_t height, const char* title = "rasterizatoy")
    : bitmap_(width, height, 1, 3, 0), display_(bitmap_, title) { }

  inline void swap_buffer() { bitmap_.display(display_); }

  inline uint32_t width() const { return display_.width(); }
  inline uint32_t height() const { return display_.height(); }
  inline bool should_close() const { return display_.is_closed(); }

  inline void set_pixel(uint32_t x, uint32_t y, const Vector3D& color)
  {
    y = height() - y - 1;
    uint8_t r = color.r < 0 ? 0 : color.r > 255 ? 255 : static_cast<uint8_t>(std::lround(color.r));
    uint8_t g = color.g < 0 ? 0 : color.g > 255 ? 255 : static_cast<uint8_t>(std::lround(color.g));
    uint8_t b = color.b < 0 ? 0 : color.b > 255 ? 255 : static_cast<uint8_t>(std::lround(color.b));
    bitmap_(x, y, 0) = r;
    bitmap_(x, y, 1) = g;
    bitmap_(x, y, 2) = b;
  }

private:
  CImg<uint8_t> bitmap_;
  CImgDisplay   display_;
};

class rasterizater
{
public:
  inline static void draw_call()
  {
    if (window_ == nullptr || shader_ == nullptr) return;

    for (Primitive primitive : primitives_)
    {
      for (Vertex& vertex : primitive.vertices)
      {
        shader_->vertex_shader(vertex);
      }
    }
  }

private:
  inline static Window* window_ = nullptr;
  inline static Shader* shader_ = nullptr;
  inline static std::vector<Primitive> primitives_ = {};
};
}

#endif //RASTERIZATOY_HPP
