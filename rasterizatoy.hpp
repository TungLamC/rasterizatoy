#ifndef RASTERIZATOY_HPP
#define RASTERIZATOY_HPP

#include <cassert>
#include <iostream>

namespace rasterizatoy {
using decimal = double;

template<uint32_t ROW, uint32_t COLUMN, typename T> class Matrix;

using Matrix4D = Matrix<4, 4, decimal>;
}

//------------------------------------------------------- 数学库 -------------------------------------------------------
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
  inline Vector(const Vector<N, U>& other) { for (auto i = 0; i < N; i++) this->components = numeric_cast<T>(other.components[i]); }

  inline T& operator[](uint32_t index) { assert(index < N); return this->components[index]; }
  inline const T& operator[](uint32_t index) const { assert(index < N); return this->components[index]; }

  inline Vector<N, T> normalize() const { return *this / length(); }
  inline Vector<N, T> normalize(std::in_place_t) const { return *this /= length(); }

  template<typename U = T>
  inline U length(bool square = false) const
  {
    U result{};
    for (auto i = 0; i < N; i++) result += std::pow(this->components[i], 2);
    return result;
  }
};

template<uint32_t N, typename T>
inline std::ostream& operator<<(std::ostream& os, const Vector<N, T>& vector)
{
  os << "[";
  for (auto i = 0; i < N; i++) os << vector[i] << (i < N - 1 ? ", " : "");
  os << "]";
  return os;
}
}

#endif //RASTERIZATOY_HPP
