#ifndef RASTERIZATOY_HPP
#define RASTERIZATOY_HPP

#include "cimg.h"

#include <any>
#include <map>
#include <vector>
#include <chrono>
#include <string>
#include <cassert>
#include <iostream>
#include <algorithm>
#include <functional>
#include <type_traits>

namespace rasterizatoy
{
using integer = int;
using decimal = double;

template<uint32_t N, typename T> struct Vector;
template<uint32_t ROW, uint32_t COLUMN, typename T> class Matrix;

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
using Matrix2D = Matrix<2, 2, decimal>;
using Matrix3D = Matrix<3, 3, decimal>;
using Matrix4D = Matrix<4, 4, decimal>;

using namespace cimg_library;
}

//-------------------------------------------------------- math --------------------------------------------------------
namespace rasterizatoy
{
template<typename T, typename U>
inline T numeric_cast(U source)
{
  if constexpr (std::is_floating_point_v<U> && std::is_integral_v<T>)
    return std::lround(source);
  return static_cast<T>(source);
}

template<typename T>
class RectArray
{
public:
  inline RectArray(uint32_t width, uint32_t height): width_(width), height_(height), value_(width, std::vector<T>(height)) {}

  inline uint32_t width() const { return width_; }
  inline uint32_t height() const { return height_; }

  inline std::vector<T>& operator[](uint32_t x) { return value_[x]; }
  inline const std::vector<T>& operator[](uint32_t x) const { return value_[x]; }

protected:
  uint32_t width_;
  uint32_t height_;
  std::vector<std::vector<T>> value_;
};

template<uint32_t N, typename T> struct VectorN { T components[N]; };
template<typename T> struct VectorN<2, T> { union { struct { T x, y; }; struct { T u, v; }; T components[2]; }; };
template<typename T> struct VectorN<3, T> { union { struct { T x, y, z; }; struct { T r, g, b; }; T components[3]; }; };
template<typename T> struct VectorN<4, T> { union { struct { T x, y, z, w; }; struct { T r, g, b, a; }; T components[4]; }; };

template<uint32_t N, typename T>
struct Vector: VectorN<N, T>
{
  inline Vector() = default;

  template<typename... U, typename = std::enable_if_t<sizeof...(U) == N && (std::is_convertible_v<T, U> && ...)>>
  inline Vector(U... args) { uint32_t i = 0; ((this->components[i++] = numeric_cast<T>(args)), ...); }

  template<typename U>
  inline Vector(const Vector<N, U>& other) { for (auto i = 0; i < N; i++) this->components[i] = numeric_cast<T>(other.components[i]); }

  inline T& operator[](uint32_t index) { assert(index < N); return this->components[index]; }
  inline const T& operator[](uint32_t index) const { assert(index < N); return this->components[index]; }

  inline Vector<N, T> normalize(std::in_place_t) { return *this /= length(); }

  template<typename U = T>
  inline Vector<N, U> normalize() const { return Vector<N, U>{*this} / length<U>(); }

  template<typename U = T>
  inline U length(bool square = false) const
  {
    U result{};
    for (auto i = 0; i < N; i++) result += numeric_cast<U>(std::pow(this->components[i], 2));
    return square ? result : std::sqrt(result);
  }
};

template<uint32_t N, typename T>
inline Vector<N, T> operator+(const Vector<N, T>& vector)
{
  return vector;
}

template<uint32_t N, typename T>
inline Vector<N, T> operator-(const Vector<N, T>& vector)
{
  Vector<N, T> result{};
  for (auto i = 0; i < N; i++) result[i] = -vector[i];
  return result;
}

template<uint32_t N, typename T>
inline bool operator==(const Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  for (auto i = 0; i < N; i++) if (lhs[i] != rhs[i]) return false;
  return true;
}

template<uint32_t N, typename T>
inline bool operator!=(const Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  return !(lhs == rhs);
}

template<uint32_t N, typename T>
inline Vector<N, T> operator+(const Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  Vector<N, T> result{};
  for (auto i = 0; i < N; i++) result[i] = lhs[i] + rhs[i];
  return result;
}

template<uint32_t N, typename T>
inline Vector<N, T>& operator+=(Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  for (auto i = 0; i < N; i++) lhs[i] += rhs[i];
  return lhs;
}

template<uint32_t N, typename T>
inline Vector<N, T> operator-(const Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  return lhs + (-rhs);
}

template<uint32_t N, typename T>
inline Vector<N, T>& operator-=(Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  for (auto i = 0; i < N; i++) lhs[i] -= rhs[i];
  return lhs;
}

template<uint32_t N, typename T>
inline Vector<N, T> operator*(const Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  Vector<N, T> result{};
  for (auto i = 0; i < N; i++) result[i] = lhs[i] * rhs[i];
  return result;
}

template<uint32_t N, typename T>
inline Vector<N, T>& operator*=(Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  for (auto i = 0; i < N; i++) lhs[i] *= rhs[i];
  return lhs;
}

template<uint32_t N, typename T>
inline Vector<N, T> operator*(const Vector<N, T>& vector, T scalar)
{
  Vector<N, T> result{};
  for (auto i = 0; i < N; i++) result[i] = numeric_cast<T>(vector[i] * scalar);
  return result;
}

template<uint32_t N, typename T>
inline Vector<N, T> operator*(T scalar, const Vector<N, T>& vector)
{
  return vector * scalar;
}

template<uint32_t N, typename T>
inline Vector<N, T>& operator*=(Vector<N, T>& vector, T scalar)
{
  for (auto i = 0; i < N; i++) vector[i] = numeric_cast<T>(vector[i] * scalar);
  return vector;
}

template<uint32_t N, typename T>
inline Vector<N, T> operator/(T scalar, const Vector<N, T>& vector)
{
  Vector<N, T> result{};
  for (uint32_t i = 0; i < N; i++) result[i] = scalar / vector;
  return result;
}

template<uint32_t N, typename T>
inline Vector<N, T> operator/(const Vector<N, T>& vector, T scalar)
{
  Vector<N, T> result{};
  for (auto i = 0; i < N; i++) result[i] = result[i] / scalar;
  return result;
}

template<uint32_t N, typename T>
inline Vector<N, T>& operator/=(Vector<N, T>& vector, T scalar)
{
  for (auto i = 0; i < N; i++) vector[i] = vector[i] / scalar;
  return vector;
}

template<uint32_t N, typename T>
inline T dot(const Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
  T result{};
  for (auto i = 0; i < N; i++) result += lhs[i] * rhs[i];
  return result;
}

template<uint32_t N, typename T>
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

template<uint32_t N, typename T>
inline Vector<N, uint8_t> to_color(const Vector<N, T>& color)
{
  static_assert(N == 3 || N == 4);
  if constexpr (N == 3)
  {
    auto r = numeric_cast<uint8_t>(std::clamp(color.r, 0.0, 255.0));
    auto g = numeric_cast<uint8_t>(std::clamp(color.g, 0.0, 255.0));
    auto b = numeric_cast<uint8_t>(std::clamp(color.b, 0.0, 255.0));
    return {r, g, b};
  }
  if constexpr (N == 4)
  {
    auto r = numeric_cast<uint8_t>(std::clamp(color.r, 0.0, 255.0));
    auto g = numeric_cast<uint8_t>(std::clamp(color.g, 0.0, 255.0));
    auto b = numeric_cast<uint8_t>(std::clamp(color.b, 0.0, 255.0));
    auto a = numeric_cast<uint8_t>(std::clamp(color.a, 0.0, 255.0));
    return {r, g, b, a};
  }
}

template<uint32_t N, typename T>
inline std::ostream& operator<<(std::ostream& os, const Vector<N, T>& vector)
{
  os << "[";
  for (auto i = 0; i < N; i++) os << vector[i] << (i < N - 1 ? ", " : "");
  os << "]";
  return os;
}

template<uint32_t ROW, uint32_t COLUMN, typename T>
class Matrix
{
public:
  inline Matrix() = default;

  template<typename... U>
  inline Matrix(U... args)
  {
    uint32_t cursor = 0;
    ((value_at(cursor++) = args), ...);
  }

  inline Vector<COLUMN, T>& operator[](uint32_t row) { assert(row < ROW); return value_[row]; }
  inline const Vector<COLUMN, T>& operator[](uint32_t row) const { assert(row < ROW); return value_[row]; }

  inline Vector<ROW, T> column_at(uint32_t column) const
  {
    Vector<ROW, T> result{};
    for (auto row = 0; row < ROW; row++)
      result[row] = value_[row][column];
    return result;
  }

private:
  inline T& value_at(uint32_t position)
  {
    assert(position < ROW * COLUMN);
    return value_[position / ROW][position % ROW];
  }

private:
  Vector<ROW, Vector<COLUMN, T>> value_;
};

template<uint32_t ROW, uint32_t COLUMN, typename T>
inline bool operator==(const Matrix<ROW, COLUMN, T>& lhs, const Matrix<ROW, COLUMN, T>& rhs)
{
  for (auto row = 0; row < ROW; row++)
    for (auto column = 0; column < COLUMN; column++)
      if (lhs[row][column] != rhs[row][column]) return false;
  return true;
}

template<uint32_t ROW, uint32_t COLUMN, typename T>
inline bool operator!=(const Matrix<ROW, COLUMN, T>& lhs, const Matrix<ROW, COLUMN, T>& rhs)
{
  return !(lhs == rhs);
}

template<uint32_t ROW, uint32_t COLUMN, typename T>
inline Matrix<ROW, COLUMN, T> operator+(const Matrix<ROW, COLUMN, T>& matrix)
{
  return matrix;
}

template<uint32_t ROW, uint32_t COLUMN, typename T>
inline Matrix<ROW, COLUMN, T> operator-(const Matrix<ROW, COLUMN, T>& matrix)
{
  Matrix<ROW, COLUMN, T> result{};
  for (auto row = 0; row < ROW; row++)
    for (auto column = 0; column < COLUMN; column++)
      result[row][column] = -matrix[row][column];
  return true;
}

template<uint32_t ROW, uint32_t COLUMN, typename T>
inline Matrix<ROW, COLUMN, T> operator+(const Matrix<ROW, COLUMN, T>& lhs, const Matrix<ROW, COLUMN, T>& rhs)
{
  Matrix<ROW, COLUMN, T> result{};
  for (size_t row = 0; row < ROW; row++)
    for (size_t column = 0; column < COLUMN; column++)
      result[row][column] = lhs[row][column] + rhs[row][column];
    return result;
}

template<uint32_t ROW, uint32_t COLUMN, typename T>
inline Matrix<ROW, COLUMN, T> operator-(const Matrix<ROW, COLUMN, T>& lhs, const Matrix<ROW, COLUMN, T>& rhs)
{
  return lhs + (-rhs);
}

template<uint32_t ROW, uint32_t LHS_COLUMN, uint32_t RHS_COLUMN, typename T>
inline Matrix<ROW, RHS_COLUMN, T> operator*(const Matrix<ROW, LHS_COLUMN, T>& lhs, const Matrix<LHS_COLUMN, RHS_COLUMN, T>& rhs)
{
  Matrix<ROW, RHS_COLUMN, T> result;
  for (auto row = 0; row < ROW; row++)
    for (auto column = 0; column < RHS_COLUMN; column++)
      result[row][column] = dot(lhs[row], rhs.column_at(column));
  return result;
}

template<uint32_t ROW, uint32_t COLUMN, typename T>
inline Vector<ROW, T> operator*(const Matrix<ROW, COLUMN, T>& lhs, const Vector<COLUMN, T>& rhs)
{
  Vector<ROW, T> result{};
  for (auto i = 0; i < ROW; i++) result[i] = dot(lhs[i], rhs);
  return result;
}

template<typename T>
inline T radians(T degrees)
{
  static_assert(std::numeric_limits<T>::is_iec559, "'radians' only accpet floating-point input.");
  return degrees * static_cast<T>(0.01745329251994329576923690768489);
}

template<uint32_t ROW, uint32_t COLUMN, typename T>
inline std::ostream& operator<<(std::ostream& os, const Matrix<ROW, COLUMN, T>& matrix)
{
  for (auto row = 0; row < ROW; row++) os << matrix[row] << std::endl;
  return os;
}

inline Matrix4D look_at(const Vector3D& eye, const Vector3D& center, const Vector3D& world_up)
{
  Vector3D z_axis = (eye - center).normalize(std::in_place);
  Vector3D x_axis = cross(world_up, z_axis).normalize(std::in_place);
  Vector3D y_axis = cross(z_axis, x_axis).normalize(std::in_place);
  return {
    x_axis.x, x_axis.y, x_axis.z, -dot(x_axis, eye),
    y_axis.x, y_axis.y, y_axis.z, -dot(y_axis, eye),
    z_axis.x, z_axis.y, z_axis.z, -dot(z_axis, eye),
    0,               0,        0,                 1,
  };
}

inline Matrix4D perspective(decimal fovy, decimal aspect, decimal z_near, decimal z_far)
{
  return {
    1 / (aspect * std::tan(fovy / 2)), 0, 0, 0,
    0, 1 / (std::tan(fovy / 2)), 0, 0,
    0, 0, -(z_far + z_near) / (z_far - z_near), -(2 * z_far * z_near) / (z_far - z_near),
    0, 0, -1, 0,
  };
}
}

//------------------------------------------------------- render -------------------------------------------------------
namespace rasterizatoy
{
template<typename... T> struct BufferLayout { using type = std::tuple<T...>; };

template<typename V, typename A>
struct Primitive
{
  using Varyings =   typename V::type;
  using Attributes = typename A::type;

  inline Primitive(const Attributes& attributes0, const Attributes& attributes1, const Attributes& attributes2)
    : varyings{}, attributes{attributes0, attributes1, attributes2} {}

  Varyings   varyings[3];
  Attributes attributes[3];
};

class Window
{
public:
  inline Window(uint32_t width, uint32_t height, const char* title = "rasterizatoy")
  : bitmap_(width, height, 1, 3, 0), display_(bitmap_, title) { }

  void swap_buffer(bool fps = true) { if (fps) display_fps(); bitmap_.display(display_); }

  inline void close() { display_.close(); }
  inline uint32_t width() const { return display_.width(); }
  inline uint32_t height() const { return display_.height(); }
  inline bool should_close() const { return display_.is_closed(); }

  inline void clear(const Vector3D& color)
  {
    auto rgb = to_color(color);
    cimg_forXY(bitmap_, x, y) { bitmap_(x, y, 0) = rgb.r; bitmap_(x, y, 1) = rgb.g; bitmap_(x, y, 2) = rgb.b; }
  }

  inline void display_fps()
  {
    static auto last = std::chrono::system_clock::now();
    auto now = std::chrono::system_clock::now();
    auto fps = 1000000 / (decimal)(std::chrono::duration_cast<std::chrono::microseconds>(now - last).count());
    const decimal foreground[] = {0 , 200, 255};
    const decimal background[] = {0, 0, 0};
    bitmap_.draw_text(0, 0, std::to_string(fps).c_str(), foreground, background, 1, 24);
    last = now;
  }


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

#define VARYING_LAYOUT(...) std::tuple<__VA_ARGS__>
#define LOCATION(location, varying) std::get<location>(varying)

template<typename Varying>
class Rasterizater
{
public:
  using VertexShader = std::function<Vector4D(uint32_t index, Varying& varying)>;
  using FragmentShader = std::function<Vector4D (const Varying&)>;

  inline Rasterizater(Window* window): window_(window), vertex_shader_(nullptr), fragment_shader_(nullptr) {}

  inline void set_vertex_shader(const VertexShader& shader) { vertex_shader_ = shader; }
  inline void set_fragment_shader(const FragmentShader& shader) { fragment_shader_ = shader; }

  inline void clear(const Vector3D& color) { window_->clear(color); }
  inline void clear(decimal r, decimal g, decimal b) { window_->clear({r, g, b}); }

  inline void swap_buffer(bool fps = true) { window_->swap_buffer(fps); }

  inline void draw_call(uint32_t vertices_count)
  {
    assert(vertices_count % 3 == 0);
    for (auto i =0; i < vertices_count; i += 3)
    {
      // 顶点着色 local space -> clip space
      Varying varying0{}; Vector4D position0 = vertex_shader_(i + 0, varying0);
      Varying varying1{}; Vector4D position1 = vertex_shader_(i + 1, varying1);
      Varying varying2{}; Vector4D position2 = vertex_shader_(i + 2, varying2);

      // todo 裁剪

      // 计算w的倒数
      decimal rhw0 = 1 / position0.w; decimal rhw1 = 1 / position1.w; decimal rhw2 = 1 / position2.w;

      // 透视除法 clip space -> ndc space
      position0 *= rhw0; position1 *= rhw1; position2 *= rhw2;

      // 计算视口坐标
      Vector2D viewport0 = {(position0.x + 1.0) * window_->width() * 0.5, (position0.y + 1.0) * window_->height() * 0.5};
      Vector2D viewport1 = {(position1.x + 1.0) * window_->width() * 0.5, (position1.y + 1.0) * window_->height() * 0.5};
      Vector2D viewport2 = {(position2.x + 1.0) * window_->width() * 0.5, (position2.y + 1.0) * window_->height() * 0.5};

      // todo 面剔除

      // edge equation
      auto [min_x, min_y, max_x, max_y] = bounding_box(viewport0, viewport1, viewport2);
      for (auto x = min_x; x <= max_x; x++)
      {
        for (auto y = min_y; y <= max_y; y++)
        {
          Vector2D viewport{x + 0.5, y + 0.5};
          Vector2D v_to_0 = viewport - viewport0;
          Vector2D v_to_1 = viewport - viewport1;
          Vector2D v_to_2 = viewport - viewport2;
          decimal a = cross(v_to_1, v_to_2);
          decimal b = cross(v_to_2, v_to_0);
          decimal c = cross(v_to_0, v_to_1);
          if (a < 0 || b < 0 || c < 0) continue;
          decimal total = a + b + c;
          a /= total; b /= total; c /= total;

          // 透视校正
          decimal rhw = a * rhw0 + b * rhw1 + c * rhw2;
          decimal w = 1.0 / rhw;
          a = rhw0 * a * w; b = rhw1 * b * w; c = rhw2 * c * w;

          // 插值varying
          Varying varying = interpolate_varying(a, b, c, varying0, varying1, varying2);

          // 片段着色
          Vector4D color = fragment_shader_(varying);

          window_->put_pixel(x, y, {color.r, color.g, color.b});
        }
      }
    }
  }

private:
  inline static std::tuple<integer, integer, integer, integer> bounding_box(Vector2I v0, Vector2I v1, Vector2I v2)
  {
    integer min_x = std::min({v0.x, v1.x, v2.x}); integer min_y = std::min({v0.y, v1.y, v2.y});
    integer max_x = std::max({v0.x, v1.x, v2.x}); integer max_y = std::max({v0.y, v1.y, v2.y});
    return std::make_tuple(min_x, min_y, max_x, max_y);
  }

  inline static Varying interpolate_varying(decimal a, decimal b, decimal c, Varying v0, Varying v1, Varying v2)
  {
    return interpolate_varying(a, b, c, v0, v1, v2, std::make_integer_sequence<uint32_t, std::tuple_size_v<Varying>>{});
  }

  template<uint32_t... I>
  inline static Varying interpolate_varying(decimal a, decimal b, decimal c, Varying v0, Varying v1, Varying v2, std::integer_sequence<uint32_t, I...>)
  {
    Varying varying{};
    ((std::get<I>(varying) = a * std::get<I>(v0) + b * std::get<I>(v1) + c * std::get<I>(v2)), ...);
    return varying;
  }

  Window* window_;
  VertexShader vertex_shader_;
  FragmentShader fragment_shader_;
};
}

#endif //RASTERIZATOY_HPP