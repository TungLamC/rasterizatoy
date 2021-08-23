#ifndef RASTERIZATOY_HPP
#define RASTERIZATOY_HPP

#include <cassert>
#include <iostream>
#include <algorithm>

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
using Matrix4D = Matrix<4, 4, decimal>;
}

namespace rasterizatoy
{
template<typename T, typename U>
inline T numeric_cast(U source)
{
  if constexpr (std::is_floating_point_v<U> && std::is_integral_v<T>)
    return std::lround(source);
  return static_cast<T>(source);
}

template<uint32_t N, typename T> struct VectorN { T components[N]; };
template<typename T> struct VectorN<2, T> { union { struct { T x, y; }; }; struct { T u, v; }; T components[2]; };
template<typename T> struct VectorN<3, T> { union { struct { T x, y, z; }; }; struct { T r, g, b; }; T components[3]; };
template<typename T> struct VectorN<4, T> { union { struct { T x, y, z, w; }; }; struct { T r, g, b, a; }; T components[4]; };

template<uint32_t N, typename T>
struct Vector: VectorN<N, T>
{
  inline Vector() = default;

  template<typename... U, typename = std::enable_if_t<sizeof...(U) == N && (std::is_convertible_v<T, U> && ...)>>
  inline Vector(U... args) { uint32_t i = 0; ((this->components[i++] = numeric_cast<T>(args)), ...); }

  template<typename U>
  inline Vector(const Vector<N, U>& other) { for (auto i = 0; i < N; i++) this->components = numeric_cast<T>(other.components[i]); }

  inline T& operator[](uint32_t index) { assert(index < N); return this->components[index]; }
  inline const T& operator[](uint32_t index) const { assert(index < N); return this->components[index]; }

  inline Vector<N, T> normalize(std::in_place_t) { return *this / length(); }

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

template<uint32_t ROW, uint32_t LHS_COLUMN, size_t RHS_COLUMN, typename T>
inline Matrix<ROW, RHS_COLUMN, T> operator*(const Matrix<ROW, LHS_COLUMN, T>& lhs, const Matrix<LHS_COLUMN, RHS_COLUMN, T>& rhs)
{
  Matrix<ROW, RHS_COLUMN, T> result;
  for (auto row = 0; row < ROW; row++)
    for (size_t column = 0; column < RHS_COLUMN; column++)
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

namespace rasterizatoy
{
}

#endif //RASTERIZATOY_HPP