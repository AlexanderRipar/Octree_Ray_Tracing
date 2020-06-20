#pragma once

#include <cstdint>

namespace och
{
	template<typename T>
	struct vec3
	{
		using value_t = T;

		vec3() {}

		vec3(T x, T y, T z) : x{ x }, y{ y }, z{ z } {}

		T x, y, z;

		vec3<T> operator*(const T v) const
		{
			return { x * v, y * v, z * v };
		}

		void operator*=(const T v)
		{
			x *= v;
			y *= v;
			z *= v;
		}

		vec3<T> operator+(const vec3<T>& v) const
		{
			return { x + v.x, y + v.y, z + v.z };
		}

		vec3<T> operator-(const vec3<T>& v) const
		{
			return { x - v.x, y - v.y, z - v.z };
		}

		vec3<T> operator+(T& s) const
		{
			return { x + s, y + s, z + s };
		}

		vec3<T> operator-(T& s) const
		{
			return { x - s, y - s, z - s };
		}

		void operator+=(const vec3<T>& v)
		{
			x += v.x;
			y += v.y;
			z += v.z;
		}

		void operator-=(const vec3<T>& v)
		{
			x -= v.x;
			y -= v.y;
			z -= v.z;
		}
	};

	template<typename T>
	struct vec2
	{
		using value_t = T;

		vec2() {}

		vec2(T x, T y) : x{ x }, y{ y } {}

		T x, y;

		vec2<T> operator*(const T v) const
		{
			return { x * v, y * v };
		}

		void operator*=(const T v)
		{
			x *= v;
			y *= v;
		}

		vec2<T> operator+(const vec2<T>& v) const
		{
			return { x + v.x, y + v.y };
		}

		vec2<T> operator-(const vec2<T>& v) const 
		{
			return { x - v.x, y - v.y };
		}

		void operator+=(const vec2<T>& v)
		{
			x += v.x;
			y += v.y;
		}

		void operator-=(const vec2<T>& v)
		{
			x -= v.x;
			y -= v.y;
		}
	};

	template<typename T>
	T dot(const vec3<T> v, const vec3<T> w)
	{
		return v.x * w.x + v.y * w.y + v.z * w.z;
	}

	using vec3f = vec3<float>;

	using vec3i32 = vec3<int32_t>;

	using float3 = vec3f;

	using float2 = vec2<float>;
}