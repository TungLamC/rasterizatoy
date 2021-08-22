#ifndef RASTERIZATOY_HPP
#define RASTERIZATOY_HPP

#include "cimg.h"

#include <map>
#include <cmath>
#include <chrono>
#include <vector>
#include <string>
#include <variant>
#include <cassert>
#include <iostream>
#include <functional>
#include <type_traits>

namespace rasterizatoy
{
using namespace cimg_library;

using integer = int;
using decimal = double;

template<size_t N, typename T>
struct Vector;
template<size_t ROW, size_t COLUMN, typename T>
class Matrix;
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

// attributes 每个顶点都有attributes, 只读
// uniforms 全局shader共用 只读
// varying 顶点着色器可写 片段着色器只读
}

//-------------------------------------------------------- math --------------------------------------------------------
namespace rasterizatoy
{
namespace derived
{
template<size_t N, typename T>
struct VectorN { T components[N]; };
template<typename T>
struct Vector2 { union { struct { T x, y; }; struct { T u, v; }; T components[2]; }; };
template<typename T>
struct Vector3 { union { struct { T x, y, z; }; struct { T r, g, b; }; T components[3]; }; };
template<typename T>
struct Vector4 { union { struct { T x, y, z, w; }; struct { T r, g, b, a; }; T components[4]; }; };
}

template<size_t N, typename T>
struct Vector: public
  std::conditional_t<N == 2, derived::Vector2<T>,
  std::conditional_t<N == 3, derived::Vector3<T>,
  std::conditional_t<N == 4, derived::Vector4<T>,
  derived::VectorN<N, T>>>>
{
  inline Vector() = default;
  inline Vector(const T* pointer) { for (size_t i = 0; i < N; i++) this->components[i] = pointer[i]; }
  inline Vector(const Vector<N, T>& other) { for (size_t i = 0; i < N; i++) this->components[i] = other.components[i]; }
  inline Vector(const std::initializer_list<T>& initializer)
  {
    for (size_t i = 0; i < N; i++)
      this->components[i] = *(initializer.begin() + i);
  }

  inline T& operator[](size_t index)
  {
    assert(index < N);
    return this->components[index];
  }
  inline const T& operator[](size_t index) const
  {
    assert(index < N);
    return this->components[index];
  }

  template<typename U>
  inline Vector<N, U> as() const { return *this; }

  template<typename U>
  inline operator Vector<N, U>() const
  {
    static_assert(std::is_convertible_v<T, U>);
    Vector<N, U> result;
    for (size_t i = 0; i < N; i++)
    {
      if constexpr (std::is_floating_point_v<T> && std::is_integral_v<U>)
        result[i] = std::lround(this->components[i]);
      else
        result[i] = static_cast<U>(this->components[i]);
    }
    return result;
  }

  inline Vector<N, T> normalize() { return *this / length(); }
  inline Vector<N, T> normalize(std::in_place_t) { return *this /= length(); }

  inline T length(bool square = false) const
  {
    T result = 0;
    for (size_t i = 0; i < N; i++) result += this->components[i] * this->components[i];
    return square ? result : std::sqrt(result);
  }
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
inline Vector<N, T> operator/(T scalar, const Vector<N, T>& vector)
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
inline Vector<N, T>& operator/=(Vector<N, T>& vector, T scalar)
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
  static_assert(N == 2 || N == 3 || N == 4);
  if constexpr (N == 2)
    return lhs.x * rhs.y - lhs.y * rhs.x;
  if constexpr (N == 3)
    return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x};
  if constexpr (N == 4)
    return {lhs.y * rhs.z - lhs.z * rhs.y, lhs.z * rhs.x - lhs.x * rhs.z, lhs.x * rhs.y - lhs.y * rhs.x, lhs.w};
}

template<size_t N, typename T>
inline Vector<N, uint8_t> to_color(const Vector<N, T>& color)
{
  static_assert(N == 3 || N == 4);
  if constexpr (N == 3)
  {
    auto r =static_cast<uint8_t>(std::clamp(color.r, 0.0, 255.0));
    auto g =static_cast<uint8_t>(std::clamp(color.g, 0.0, 255.0));
    auto b =static_cast<uint8_t>(std::clamp(color.b, 0.0, 255.0));
    return {r, g, b};
  }
  if constexpr (N == 4)
  {
    auto r =static_cast<uint8_t>(std::clamp(color.r, 0.0, 255.0));
    auto g =static_cast<uint8_t>(std::clamp(color.g, 0.0, 255.0));
    auto b =static_cast<uint8_t>(std::clamp(color.b, 0.0, 255.0));
    auto a =static_cast<uint8_t>(std::clamp(color.a, 0.0, 255.0));
    return {r, g, b, a};
  }
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
  inline Vector<COLUMN, T>& operator[](size_t row)
  {
    assert(row < ROW);
    return value[row];
  }
  inline const Vector<COLUMN, T>& operator[](size_t row) const
  {
    assert(row < ROW);
    return value[row];
  }

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
inline Matrix<ROW, RHS_COLUMN, T>
operator*(const Matrix<ROW, LHS_COLUMN, T>& lhs, const Matrix<LHS_COLUMN, RHS_COLUMN, T>& rhs)
{
  Matrix<ROW, RHS_COLUMN, T> result;
  for (size_t row = 0; row < ROW; row++)
    for (size_t column = 0; column < RHS_COLUMN; column++)
      result[row][column] = dot(lhs[row], rhs.column_at(column));
  return result;
}

template<typename T>
inline T radians(T degrees)
{
  static_assert(std::numeric_limits<T>::is_iec559, "'radians' only accpet floating-point input.");
  return degrees * static_cast<T>(0.01745329251994329576923690768489);
}

inline Matrix4D look_at(const Vector3D& eye, const Vector3D& center, const Vector3D& world_up)
{
  Vector3D z_axis = (eye - center).normalize(std::in_place);
  Vector3D x_axis = cross(world_up, z_axis).normalize(std::in_place);
  Vector3D y_axis = cross(z_axis, x_axis).normalize(std::in_place);
  return {
    {x_axis.x, x_axis.y, x_axis.z, -dot(x_axis, eye)},
    {y_axis.x, y_axis.y, y_axis.z, -dot(y_axis, eye)},
    {z_axis.x, z_axis.y, z_axis.z, -dot(z_axis, eye)},
    {0,               0,        0,                 1}
  };
}

inline Matrix4D perspective(decimal fovy, decimal aspect, decimal z_near, decimal z_far)
{
  return {
    {1 / (aspect * std::tan(fovy / 2)), 0, 0, 0},
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
enum class Facing { Back, Front };

struct Vertex
{
  inline Vertex(const Vector3D& position, const Vector4D& color = {0, 0, 0, 0})
    : rhw{1}, normal{0, 0}, texcoord{0, 0}, position{position.x, position.y, position.z, 1}, color{color}, viewport{0, 0} {}

  decimal rhw;
  Vector4D color;
  Vector4D normal;
  Vector2D texcoord;
  Vector2D viewport;
  Vector4D position;
};

struct Primitive
{
  inline Primitive(const Vertex& vertex0, const Vertex& vertex1, const Vertex& vertex2)
    : vertices{vertex0, vertex1, vertex2} {}

  inline Facing facing() const
  {
    auto v0 = vertices[0].viewport;
    auto v1 = vertices[1].viewport;
    auto v2 = vertices[2].viewport;
    return cross(v1 - v0, v2 - v0) < 0 ? Facing::Back : Facing::Front;
  }

  Vertex vertices[3];
};

struct Fragment
{
  Vector3D color;
  Vector4D normal;
  Vector2D texcoord;
  Vector2I viewport;
  Vector4D position;
};

struct Shader
{
  virtual void vertex_shader(Vertex& vertex) {}
  virtual void fragment_shader(Fragment& fragment) {}
};

class Window
  {
  public:
    friend class rasterizater;

    inline Window(uint32_t width, uint32_t height, const char* title = "rasterizatoy")
    : bitmap_(width, height, 1, 3, 0), display_(bitmap_, title) { }

    void swap_buffer() { bitmap_.display(display_); }

    inline void close() { display_.close(); }
    inline uint32_t width() const { return display_.width(); }
    inline uint32_t height() const { return display_.height(); }
    inline bool should_close() const { return display_.is_closed(); }
    inline void clear(const Vector3D& color) { auto rgb = to_color(color); bitmap_.fill(rgb.r, rgb.g, rgb.b); }

    inline void put_pixel(uint32_t x, uint32_t y, const Vector3D& color)
    {
      auto rgb = to_color(color);
      y = height() - y - 1;
      bitmap_(x, y, 0) = rgb.r;
      bitmap_(x, y, 1) = rgb.g;
      bitmap_(x, y, 2) = rgb.b;
    }

    void draw_line(uint32_t x1, uint32_t y1, uint32_t x2, uint32_t y2, const Vector3D& color)
    {
      if (x1 == x2 && y1 == y2) { put_pixel(x1, y1, color); return; }
      if (x1 == x2) { for (auto y = y1; y != y2; y += (y1 <= y2 ? 1 : -1)) put_pixel(x1, y, color); return; }
      if (y1 == y2) { for (auto x = x1; x != x2; x += (x1 <= x2 ? 1 : -1)) put_pixel(x, y1, color); return; }

      uint32_t x, y;
      uint32_t rem = 0;
      uint32_t dx = (x1 < x2)? x2 - x1 : x1 - x2;
      uint32_t dy = (y1 < y2)? y2 - y1 : y1 - y2;
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
  inline static void set_shader(Shader* shader) { shader_ = shader; }

  inline static void set_current_context(Window* window) { window_ = window; }

  inline static void input_primitives(const std::vector<Primitive>& primitives) { primitives_ = primitives; }

  inline static void clear(const Vector3D& color) { window_->clear(color); }

  inline static void swap_buffer(bool fps = true) { if (fps) display_fps(); window_->swap_buffer(); }

  inline static void draw_call(bool fps = true)
  {
    for (auto primitive : primitives_)
    {
      for (auto& vertex : primitive.vertices)
      {
        // 顶点着色
        shader_->vertex_shader(vertex);
        vertex.rhw = 1 / vertex.position.w;

        // todo 裁剪

        // 透视除法
        vertex.position *= vertex.rhw;

        // ndc space -> viewport space
        vertex.viewport.x = (vertex.position.x + 1.0) * window_->width() * 0.5;
        vertex.viewport.y = (vertex.position.y + 1.0) * window_->height() * 0.5;
      }

      // 面剔除
      if (primitive.facing() == Facing::Back) continue;

      // edge equation
      if (primitive.facing() == Facing::Back) std::swap(primitive.vertices[1], primitive.vertices[2]);

      auto& vertices = primitive.vertices;

      auto [min_x, min_y, max_x, max_y] = bounding_box(primitive);

      Fragment fragment{};
      for (fragment.viewport.x = min_x; fragment.viewport.x <= max_x; fragment.viewport.x++)
      {
        for (fragment.viewport.y = min_y; fragment.viewport.y <= max_y; fragment.viewport.y++)
        {
          // 计算viewport space下的插值系数
          Vector2D point{fragment.viewport.x + 0.5, fragment.viewport.y + 0.5};
          Vector2D ap = point - vertices[0].viewport;
          Vector2D bp = point - vertices[1].viewport;
          Vector2D cp = point - vertices[2].viewport;
          decimal a = cross(bp, cp);
          decimal b = cross(cp, ap);
          decimal c = cross(ap, bp);

          if (a < 0 || b < 0 || c < 0) continue;

          decimal s = a + b + c;
          a /= s;
          b /= s;
          c /= s;

          // 计算透视校正的插值系数 https://zhuanlan.zhihu.com/p/380634903
          decimal rhw = a * vertices[0].rhw + b * vertices[1].rhw + c * vertices[2].rhw;
          decimal w = 1.0 / rhw;
          a = vertices[0].rhw * a * w;
          b = vertices[1].rhw * b * w;
          c = vertices[2].rhw * c * w;

          // 插值varying变量
//
//          Vector4D color = a * vertices[0].color + b * vertices[1].color + c * vertices[2].color;
//          fragment.color = {color.r, color.g, color.b};
//          shader_->fragment_shader(fragment);

          window_->put_pixel(fragment.viewport.x, fragment.viewport.y, fragment.color);
        }
      }
    }
  }

private:
  inline static std::tuple<integer, integer, integer, integer> bounding_box(const Primitive& primitive)
  {
    const auto& v0 = primitive.vertices[0].viewport.as<integer>();
    const auto& v1 = primitive.vertices[1].viewport.as<integer>();
    const auto& v2 = primitive.vertices[2].viewport.as<integer>();
    integer min_x = std::min({v0.x, v1.x, v2.x});
    integer min_y = std::min({v0.y, v1.y, v2.y});
    integer max_x = std::max({v0.x, v1.x, v2.x});
    integer max_y = std::max({v0.y, v1.y, v2.y});
    // todo 与屏幕宽高做对比
    return std::make_tuple(min_x, min_y, max_x, max_y);
  }

  inline static void display_fps()
  {
    static auto last = std::chrono::system_clock::now();
    auto now = std::chrono::system_clock::now();
    auto fps = 1000000 / (decimal)(std::chrono::duration_cast<std::chrono::microseconds>(now - last).count());
    const decimal foreground[] = {0 , 200, 255};
    const decimal background[] = {0, 0, 0};
    window_->bitmap_.draw_text(0, 0, std::to_string(fps).c_str(), foreground, background, 1, 24);
    last = now;
  }

private:
  inline static Window* window_;
  inline static Shader* shader_;
  inline static std::vector<Primitive> primitives_;
};
}

#endif //RASTERIZATOY_HPP
