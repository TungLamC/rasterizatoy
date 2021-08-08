#ifndef RASTERIZATOY_HPP
#define RASTERIZATOY_HPP

#include <string>
#include <cassert>
#include <type_traits>

namespace rasterizatoy
{
using integer = long;
using decimal = float;
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
}

#endif //RASTERIZATOY_HPP
