#ifndef RASTERIZATOY_HPP
#define RASTERIZATOY_HPP

#include "cimg.h"
#include <vector>
#include <cassert>
#include <utility>
#include <initializer_list>

namespace rasterizatoy
{
using namespace cimg_library;

using real_t = double;

template<size_t N, typename T> class Vector;
template<size_t ROW, size_t COLUMN, typename T> class Matrix;

template<typename T> using Vector2 = Vector<2, T>;
template<typename T> using Vector3 = Vector<3, T>;
template<typename T> using Vector4 = Vector<4, T>;
template<typename T> using Matrix4 = Matrix<4, 4, T>;

using ColorRGB    = Vector3<uint8_t>;
using ColorBuffer = CImg<uint8_t>*;
}

//-------------------------------------------------------- math --------------------------------------------------------
namespace rasterizatoy
{
namespace derived
{
template<size_t N, typename T> struct VectorN { T components[N]; };
template<typename T> struct Vector2 { union { struct { T x, y; }; struct { T u, v; }; T components[2]; }; };
template<typename T> struct Vector3 { union { struct { T x, y, z; }; struct { T r, g, b; }; T components[3]; }; };
template<typename T> struct Vector4 { union { struct { T x, y, z, w; }; struct { T r, g, b, a; }; T components[4]; }; };
}

template<size_t N, typename T>
class Vector: public
  std::conditional_t<N == 2, derived::Vector2<T>,
  std::conditional_t<N == 3, derived::Vector3<T>,
  std::conditional_t<N == 4, derived::Vector4<T>,
  derived::VectorN<N, T>>>>
{
public:
  inline Vector() = default;
  inline Vector(const T* pointer) { for (size_t i = 0; i < N; i++) this->components[i] = pointer[i]; }
  inline Vector(const Vector<N, T>& other) { for (size_t i = 0; i < N; i++) this->components[i] = other.components[i]; }
  inline Vector(const std::initializer_list<T>& initializer) { for (size_t i = 0; i < N; i++) this->components[i] = *(initializer.begin() + i); }

  inline T& operator[](size_t index) { assert(index < N); return this->components[index]; }
  inline const T& operator[](size_t index) const { assert(index < N); return this->components[index]; }

  inline Vector<N, T> normalize() { return *this / length(); }
  inline Vector<N, T> normalize(std::in_place_t) { return *this /= length(); }
  inline T length(bool square = false) const { T result = 0; for (size_t i = 0; i < N; i++) result += this->components[i] * this->components[i]; return result; }
};

template<size_t N, typename T>
inline Vector<N, T> operator+(const Vector<N, T>& vector)
{
  return vector;
}

template<size_t N, typename T>
inline Vector<N, T> operator-(const Vector<N, T>& vector)
{
  Vector<N, T> result;
  for (size_t i = 0; i < N; i++) result[i] = -vector[i];
  return result;
}

template<size_t N, typename T>
inline Vector<N, T> operator==(const Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  for (size_t i = 0; i < N; i++) if (lhs[i] != rhs[i]) return false;
  return true;
}

template<size_t N, typename T>
inline Vector<N, T> operator!=(const Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  return !(lhs == rhs);
}

template<size_t N, typename T>
inline Vector<N, T> operator+(const Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  Vector<N, T> result;
  for (size_t i = 0; i < N; i++) result[i] = lhs[i] + rhs[i];
  return result;
}

template<size_t N, typename T>
inline Vector<N, T>& operator+=(Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  for (size_t i = 0; i < N; i++) lhs[i] = lhs[i] + rhs[i];
  return lhs;
}

template<size_t N, typename T>
inline Vector<N, T> operator-(const Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  return lhs + (-rhs);
}

template<size_t N, typename T>
inline Vector<N, T>& operator-=(Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  return lhs += (-rhs);
}

template<size_t N, typename T>
inline Vector<N, T> operator*(const Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  Vector<N, T> result;
  for (size_t i = 0; i < N; i++) result[i] = lhs[i] * rhs[i];
  return result;
}

template<size_t N, typename T>
inline Vector<N, T> operator*(const Vector<N, T>& vector, T scalar)
{
  Vector<N, T> result;
  for (size_t i = 0; i < N; i++) result[i] = vector[i] * scalar;
  return result;
}

template<size_t N, typename T>
inline Vector<N, T> operator*(T scalar, const Vector<N, T>& vector)
{
  return vector * scalar;
}

template<size_t N, typename T>
inline Vector<N, T>& operator*=(Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  for (size_t i = 0; i < N; i++) lhs[i] = lhs[i] * rhs[i];
  return lhs;
}

template<size_t N, typename T>
inline Vector<N, T>& operator*=(Vector<N, T>& vector, T scalar)
{
  for (size_t i = 0; i < N; i++) vector[i] = vector[i] * scalar;
  return vector;
}

template<size_t N, typename T>
inline Vector<N, T> operator/(const Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  Vector<N, T> result;
  for (size_t i = 0; i < N; i++) result[i] = lhs[i] / rhs[i];
  return result;
}

template<size_t N, typename T>
inline Vector<N, T> operator/(const Vector<N, T>& vector, T scalar)
{
  Vector<N, T> result;
  for (size_t i = 0; i < N; i++) result[i] = vector[i] / scalar;
  return result;
}

template<size_t N, typename T>
inline Vector<N, T> operator/( T scalar, const Vector<N, T>& vector)
{
  Vector<N, T> result;
  for (size_t i = 0; i < N; i++) result[i] = scalar / vector;
  return result;
}

template<size_t N, typename T>
inline Vector<N, T>& operator/=(Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  for (size_t i = 0; i < N; i++) lhs[i] = lhs[i] / rhs[i];
  return lhs;
}

template<size_t N, typename T>
inline Vector<N, T>& operator/=(const Vector<N, T>& vector, T scalar)
{
  for (size_t i = 0; i < N; i++) vector[i] = vector[i] / scalar;
  return vector;
}

template<size_t N, typename T>
inline T dot(const Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  T result = 0;
  for (size_t i = 0; i < N; i++) result += lhs[i] * rhs[i];
  return result;
}

template<size_t N, typename T>
inline std::conditional_t<N == 2, T, Vector<N, T>> cross(const Vector<N, T>& lhs, const Vector<N, T>& rhs) 
{
  static_assert(N == 2 | N == 3 || N == 4);
  if constexpr (N == 2)
    return lhs.x * rhs.y - lhs.y * rhs.x;
  if constexpr (N == 3)
    return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x};
  if constexpr (N == 4)
    return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x, lhs.w};
}

template<size_t N, typename T>
inline std::ostream& operator<<(std::ostream& stream, const Vector<N, T>& vector)
{
  stream << "[";
  for (size_t i = 0; i < N; i++) stream << vector[i] << (i < N - 1 ? ", " : "");
  stream << "]";
  return stream;
}

template<size_t ROW, size_t COLUMN, typename T>
class Matrix
{
public:
  inline Matrix() = default;

  inline Matrix(const std::initializer_list<Vector<COLUMN, T>>& initializer)
  {
    for (size_t row = 0; row < ROW; row++)
      this->operator[](row) = *(initializer.begin() + row);
  }

public:
  inline Vector<COLUMN, T>& operator[](size_t row) { assert(row < ROW); return value[row]; }
  inline const Vector<COLUMN, T>& operator[](size_t row) const { assert(row < ROW); return value[row]; }

  inline Vector<ROW, T> column_at(size_t column) const
  {
    Vector<ROW, T> result;
    for (size_t row = 0; row < ROW; row++)
      result[row] = value[row][column];
    return result;
  }

private:
  Vector<ROW, Vector<COLUMN, T>> value;
};

template<size_t ROW, size_t COLUMN, typename T>
inline bool operator==(const Matrix<ROW, COLUMN, T>& lhs, const Matrix<ROW, COLUMN, T>& rhs)
{
  for (size_t row = 0; row < ROW; row++)
    for (size_t column = 0; column < COLUMN; column++)
      if (lhs.matrix_[row][column] != rhs.matrix_[row][column]) return false;
      return true;
}

template<size_t ROW, size_t COLUMN, typename T>
inline bool operator!=(const Matrix<ROW, COLUMN, T>& lhs, const Matrix<ROW, COLUMN, T>& rhs)
{
  return !(lhs == rhs);
}

template<size_t ROW, size_t COLUMN, typename T>
inline Matrix<ROW, COLUMN, T> operator+(const Matrix<ROW, COLUMN, T>& matrix)
{
  return matrix;
}

template<size_t ROW, size_t COLUMN, typename T>
inline Matrix<ROW, COLUMN, T> operator-(const Matrix<ROW, COLUMN, T>& matrix)
{
  Matrix<ROW, COLUMN, T> result;
  for (size_t row = 0; row < ROW; row++)
    for (size_t column = 0; column < COLUMN; column++)
      result[row][column] = -matrix[row][column];
  return result;
}

template<size_t ROW, size_t COLUMN, typename T>
inline Matrix<ROW, COLUMN, T> operator+(const Matrix<ROW, COLUMN, T>& lhs, const Matrix<ROW, COLUMN, T>& rhs)
{
  Matrix<ROW, COLUMN, T> result;
  for (size_t row = 0; row < ROW; row++)
    for (size_t column = 0; column < COLUMN; column++)
      result[row][column] = lhs[row][column] + rhs[row][column];
  return result;
}

template<size_t ROW, size_t COLUMN, typename T>
inline Matrix<ROW, COLUMN, T> operator-(const Matrix<ROW, COLUMN, T>& lhs, const Matrix<ROW, COLUMN, T>& rhs)
{
  return lhs + (-rhs);
}

template<size_t ROW, size_t COLUMN, typename T>
inline Vector<ROW, T> operator*(const Matrix<ROW, COLUMN, T>& lhs, const Vector<COLUMN, T>& rhs)
{
  Vector<ROW, T> result;
  for (size_t i = 0; i < ROW; i++)
    result[i] = dot(lhs[i], rhs);
  return result;
}

template<size_t ROW, size_t LHS_COLUMN, size_t RHS_COLUMN, typename T>
inline Matrix<ROW, RHS_COLUMN, T> operator*(const Matrix<ROW, LHS_COLUMN, T>& lhs, const Matrix<LHS_COLUMN, RHS_COLUMN, T>& rhs)
{
  Matrix<ROW, RHS_COLUMN, T> result;
  for (size_t row = 0; row < ROW; row++)
    for (size_t column = 0; column < RHS_COLUMN; column++)
      result[row][column] = dot(lhs[row], rhs.column_at(column));
    return result;
}

template<typename T>
inline Matrix4<T> look_at(const Vector3<T>& eye, const Vector3<T>& center, const Vector3<T>& world_up)
{
  Vector3<T> z_axis = (eye - center).normalize(std::in_place);
  Vector3<T> x_axis = cross(world_up, z_axis).normalize(std::in_place);
  Vector3<T> y_axis = cross(z_axis, x_axis).normalize(std::in_place);
  return {
    {x_axis.x, x_axis.y, x_axis.z, -dot(x_axis, eye)},
    {y_axis.x, y_axis.y, y_axis.z, -dot(y_axis, eye)},
    {z_axis.x, z_axis.y, z_axis.z, -dot(z_axis, eye)},
    {       0,        0,        0,                 1}
  };
}

template<typename T>
inline Matrix4<T> perspective(T fovy, T aspect, T z_near, T z_far)
{
  fovy = fovy * static_cast<T>(0.01745329251994329576923690768489);
  return {
    {1 / (aspect * std::tan(fovy / 2)) , 0, 0, 0},
    {0, 1 / (std::tan(fovy / 2)), 0, 0},
    {0, 0, -(z_far + z_near) / (z_far - z_near), -(2 * z_far * z_near) / (z_far - z_near)},
    {0, 0, -1, 0}
  };
}

template<size_t ROW, size_t COLUMN, typename T>
inline std::ostream& operator<<(std::ostream& stream, const Matrix<ROW, COLUMN, T>& matrix)
{
  for (size_t row = 0; row < ROW; row++) stream << matrix[row] << std::endl;
  return stream;
}
}

//------------------------------------------------------- render -------------------------------------------------------
namespace rasterizatoy
{
struct Vertex
{
  Vertex(Vector4<real_t> position = {0, 0, 0, 0}): position(position) { }
  real_t rhw;
  Vector4<real_t> position;
  Vector2<size_t> viewport;
};

struct Primitive { Vertex vertices[3]; };

struct Shader
{
  virtual void vertex_shader(Vertex& vertex) { }
};

class Window
{
public:
  friend class rasterizater;

  inline Window(size_t width, size_t height, const char* title = "rasterizatoy")
    : bitmap_(width, height, 1, 3, 0), display_(bitmap_, title) { }

  void swap_buffer() { bitmap_.display(display_); }

  inline void close() { display_.close(); }
  inline size_t width() const { return display_.width(); }
  inline size_t height() const { return display_.height(); }
  inline bool should_close() const { return display_.is_closed(); }

  inline void clear(const ColorRGB& color) { cimg_forXY(bitmap_, x, y) { put_pixel(x, y, color); } }

  inline void put_pixel(size_t x, size_t y, const ColorRGB& color)
  {
    y = height() - y - 1;
    bitmap_(x, y, 0) = color.r;
    bitmap_(x, y, 1) = color.g;
    bitmap_(x, y, 2) = color.b;
  }

  void draw_line(size_t x1, size_t y1, size_t x2, size_t y2, const ColorRGB& color)
  {
    if (x1 == x2 && y1 == y2) { put_pixel(x1, y1, color); return; }
    if (x1 == x2) { for (size_t y = y1; y != y2; y += (y1 <= y2 ? 1 : -1)) put_pixel(x1, y, color); return; }
    if (y1 == y2) { for (size_t x = x1; x != x2; x += (x1 <= x2 ? 1 : -1)) put_pixel(x, y1, color); return; }

    size_t x, y;
    size_t dx = (x1 < x2)? x2 - x1 : x1 - x2;
    size_t dy = (y1 < y2)? y2 - y1 : y1 - y2;
    int rem = 0;
    if (dx >= dy)
    {
      if (x2 < x1) x = x1, y = y1, x1 = x2, y1 = y2, x2 = x, y2 = y;
      for (x = x1, y = y1; x <= x2; x++)
      {
        put_pixel(x, y, color);
        rem += dy;
        if (rem >= dx) { rem -= dx; y += (y2 >= y1)? 1 : -1; put_pixel(x, y, color); }
      }
      put_pixel(x2, y2, color);
    }
    else
    {
      if (y2 < y1) x = x1, y = y1, x1 = x2, y1 = y2, x2 = x, y2 = y;
      for (x = x1, y = y1; y <= y2; y++)
      {
        put_pixel(x, y, color);
        rem += dx;
        if (rem >= dy) { rem -= dy; x += (x2 >= x1)? 1 : -1; put_pixel(x, y, color); }
      }
      put_pixel(x2, y2, color);
    }
  }

private:
  CImg<uint8_t> bitmap_;
  CImgDisplay   display_;
};

class rasterizater
{
public:
  inline static void clear(const ColorRGB& color) { window_->clear(color); }
  inline static void set_current_context(Window* window) { window_ = window; color_buffer_ = &window->bitmap_; }
  inline static void swap_buffer() { if (!window_) return; window_->swap_buffer(); }

  inline static void input_primitives(const std::vector<Primitive>& primitives) { primitives_ = primitives; }

  inline static void set_shader(Shader* shader) { shader_ = shader; }

  inline static void draw_call()
  {
    for (Primitive& primitive : primitives_)
    {
      for (Vertex& vertex : primitive.vertices)
      {
        // 顶点着色
        shader_->vertex_shader(vertex);
        vertex.rhw = 1 / vertex.position.w;
        vertex.viewport.x = ((vertex.position.x + 1) / 2) * window_->width();
        vertex.viewport.y = ((vertex.position.y + 1) / 2) * window_->height();
      }
      auto& vertices = primitive.vertices;
      window_->draw_line(vertices[0].viewport.x, vertices[0].viewport.y, vertices[1].viewport.x, vertices[1].viewport.y, {255, 0, 255});
      window_->draw_line(vertices[0].viewport.x, vertices[0].viewport.y, vertices[2].viewport.x, vertices[2].viewport.y, {255, 0, 255});
      window_->draw_line(vertices[2].viewport.x, vertices[2].viewport.y, vertices[1].viewport.x, vertices[1].viewport.y, {255, 0, 255});
    }
  }

public:
  inline static Window*     window_;
  inline static Shader*     shader_;
  inline static ColorBuffer color_buffer_;
  inline static std::vector<Primitive> primitives_;
};
}

#endif //RASTERIZATOY_HPP