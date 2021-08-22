#pragma once

#include "cimg.h"
#include <cassert>
#include <iostream>

namespace rasterizatoy {
using integer = int;
using decimal = double;

template<uint32_t N, typename T> struct Vector;
template<uint32_t ROW, uint32_t COLUMN, typename T> class Matrix;

using Vector2I = Vector<2, integer>;
using Vector3I = Vector<3, integer>;
using Vector4I = Vector<4, integer>;
using Vector2D = Vector<2, decimal>;
using Vector3D = Vector<3, decimal>;
using Vector4D = Vector<4, decimal>;
using Matrix4D = Matrix<4, 4, decimal>;
}

namespace rasterizatoy {

template<uint32_t N, typename T> struct VectorN { T components[N]; };
template<typename T> struct VectorN<2, T> { union { struct { T x, y; }; struct { T u, v; }; T components[2]; }; };
template<typename T> struct VectorN<3, T> { union { struct { T x, y, z; }; struct { T r, g, b; }; T components[3]; }; };
template<typename T> struct VectorN<4, T> { union { struct { T x, y, z, w; }; struct { T r, g, b, a; }; T components[4]; }; };

template<typename T, typename U, typename = std::enable_if_t<std::is_arithmetic_v<T> && std::is_arithmetic_v<U>>>
inline T component_cast(U source) {
	if constexpr (std::is_floating_point_v<U> && std::is_integral_v<T>)
		return std::lround(source);
	return static_cast<T>(source);
}

///------------------------------------------------------ 通用矢量 ------------------------------------------------------
template<uint32_t N, typename T>
struct Vector: VectorN<N, T> {
	inline Vector() { for (auto i = 0; i < N; i++) this->components[i] = T{}; }
	inline Vector(const Vector<N, T>& other) { for (auto i = 0; i < N; i++) this->components[i] = other.components[i]; }

	template<typename... U, typename = std::enable_if_t<sizeof...(U) == N && (std::is_convertible_v<T, U> && ...)>>
	inline Vector(U... args) {
		uint32_t i = 0;
		((this->components[i++] = component_cast<T>(args)), ...);
	}

	inline T& operator[](uint32_t index) { assert(index < N); return this->components[index]; }
	inline const T& operator[](uint32_t index) const { assert(index < N); return this->components[index]; }

	inline Vector<N, T> normalize() { return *this / length(); }
	inline Vector<N, T> normalize(std::in_place_t) { return *this /= length(); }

	template<typename U = T>
	inline U length(bool square = false) const {
		U result = 0;
		for (auto i = 0; i < N; i++) result += this->components[i] * this->components[i];
		return square ? result : std::sqrt(result);
	}
};

///------------------------------------------------------ 矢量运算 ------------------------------------------------------

template<uint32_t N, typename T>
inline Vector<N, T> operator+(const Vector<N, T>& vector) {
	return vector;
}

template<uint32_t N, typename T, uint32_t... I>
inline Vector<N, T> operator-(const Vector<N, T>& vector) {
	Vector<N, T> result{};
	for (auto i = 0; i < N; i++) result[i] = -vector[i];
	return result;
}

template<uint32_t N, typename T>
inline bool operator==(const Vector<N, T>& lhs, const Vector<N, T>& rhs) {
	for (auto i = 0; i < N; i++) if (lhs[i] != rhs[i]) return false;
	return true;
}

template<uint32_t N, typename T>
inline bool operator!=(const Vector<N, T>& lhs, const Vector<N, T>& rhs) {
	return !(lhs == rhs);
}

template<uint32_t N, typename T>
inline Vector<N, T> operator+(const Vector<N, T>& lhs, const Vector<N, T>& rhs) {
	Vector<N, T> result{};
	for (auto i = 0; i < N; i++)
		result[i] = lhs[i] + rhs[i];
	return result;
}
template<uint32_t N, typename T>
inline Vector<N, T>& operator+=(Vector<N, T>& lhs, const Vector<N, T>& rhs) {
	for (size_t i = 0; i < N; i++) lhs[i] = lhs[i] + rhs[i];
	return lhs;
}

template<uint32_t N, typename T>
inline Vector<N, T> operator-(const Vector<N, T>& lhs, const Vector<N, T>& rhs) {
	return lhs + (-rhs);
}

template<uint32_t N, typename T>
inline Vector<N, T> operator-=(Vector<N, T>& lhs, const Vector<N, T>& rhs) {
	return lhs += (-rhs);
}

template<uint32_t N, typename T>
inline Vector<N, T> operator*(const Vector<N, T>& lhs, const Vector<N, T>& rhs)
{
	Vector<N, T> result;
	for (auto i = 0; i < N; i++) result[i] = lhs[i] * rhs[i];
	return result;
}

template<uint32_t N, typename T>
inline Vector<N, T>& operator*=(Vector<N, T>& lhs, const Vector<N, T>& rhs) {
	for (size_t i = 0; i < N; i++) lhs[i] = lhs[i] * rhs[i];
	return lhs;
}

template<uint32_t N, typename T, typename U>
inline Vector<N, T>& operator*=(Vector<N, T>& vector, U scalar) {
	for (size_t i = 0; i < N; i++) vector[i] = component_cast<T>(vector[i] * scalar);
	return vector;
}

template<uint32_t N, typename T, typename U>
inline Vector<N, T>& operator*=(U scalar, Vector<N, T>& vector) {
	return vector * scalar;
}

template<uint32_t N, typename T, typename U = T>
inline U dot(const Vector<N, T>& lhs, const Vector<N, T>& rhs) {
	
}

template<uint32_t N, typename T>
inline std::ostream& operator<<(std::ostream& os, const Vector<N, T>& vector) {
	os << "[";
	for (auto i = 0; i < N; i++) os << vector[i] << (i < N - 1 ? ", " : "");
	os << "]";
	return os;
}
}