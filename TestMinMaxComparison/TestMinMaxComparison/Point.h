#pragma once

#include <array>
/*
*
* Copyright (C) 2018 Carlos Carbonera - All Rights Reserved
* You may use this code with explicit consent from the author.
* It is illegal to reproduce or modify this code without the
* author's explicit consent.
*
*/

#pragma once
template<typename ScalarType, int Dimension>
using Point = std::array<ScalarType, Dimension>;

template <typename PT>
struct Dimension {
	static int constexpr value = 0;
};

template <typename PT>
int constexpr Dimension_v = Dimension<PT>::value;

template <typename PT>
struct PointTraits {
};

template <typename PT>
struct IsPoint
{
	static bool constexpr value = false;
};

template <typename PT>
int constexpr IsPoint_v = IsPoint<PT>::value;

template <typename PT>
using PointTraits_t = typename PointTraits<PT>::type;

template <typename PT>
using ScalarType_t = typename PT::value_type;

template <typename PT>
inline PointTraits_t<PT> & operator-(PT const & a)
{
	PT minus_a;
	for (auto i = 0; i < Dimension_v<PT>; ++i)
		minus_a[i] = -a[i];
	return minus_a;
}

template <typename PT>
inline PointTraits_t<PT> & operator+(PT const & a)
{
	return a;
}

template <typename PT>
inline PointTraits_t<PT> & operator+=(PT & a, PT const & other)
{
	for (auto i = 0; i < Dimension_v<PT>; ++i)
		a[i] += other[i];
	return a;
}

template <typename PT>
inline PointTraits_t<PT> & operator-=(PT & a, PT const & other)
{
	for (auto i = 0; i < Dimension_v<PT>; ++i)
		a[i] -= other[i];
	return a;
}

template <typename PT>
inline PointTraits_t<PT> & operator*=(PT & a, ScalarType_t<PT> value)
{
	for (auto i = 0; i < Dimension_v<PT>; ++i)
		a[i] *= value;
	return a;
}

template <typename PT>
inline PointTraits_t<PT> & operator/=(PT & a, ScalarType_t<PT> value)
{
	for (auto i = 0; i < Dimension_v<PT>; ++i)
		a[i] /= value;
	return a;
}

template <typename PT>
inline PointTraits_t<PT> operator+(PT a, PT const & b)
{
	return (a += b);
}

template <typename PT>
inline PointTraits_t<PT> operator-(PT a, PT const & b)
{
	return (a -= b);
}

template <typename PT>
inline PointTraits_t<PT> operator*(PT a, ScalarType_t<PT> value)
{
	return (a *= value);
}

template <typename PT>
inline PointTraits_t<PT> operator*(ScalarType_t<PT> value, PT a)
{
	return (a *= value);
}

template <typename PT>
inline PointTraits_t<PT> operator/(PT a, ScalarType_t<PT> value)
{
	return (a /= value);
}

template <typename DiagonalMatrix, typename PT>
PT diagonal_multiplication(DiagonalMatrix const & diagonal_matrix, PT const & pt)
{
	PT r;
	r[0] = pt[0] * diagonal_matrix[0];
	r[1] = pt[1] * diagonal_matrix[1];
	if (Dimension_v<PT> > 2)
	{
		r[2] = pt[2] * diagonal_matrix[2];
		for (auto i = 3; i < Dimension_v<PT>; ++i)
			r[2] *= diagonal_matrix[2];
	}
	return r;
}

template <typename PT>
ScalarType_t<PT> dot(PT const & a, PT const & b)
{
	auto result = a[0] * b[0];
	for (auto i = 1; i < Dimension_v<PT>; ++i)
		result += a[i] * b[i];

	return result;
}

template <typename ScalarType>
Point<ScalarType, 3> cross(Point<ScalarType, 3> const & a, Point<ScalarType, 3> const & b)
{
	return Point<ScalarType, 3>({ a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0] });
}

template <typename PT>
ScalarType_t<PT> length_sq(PT const & a)
{
	return dot(a, a);
}

template <typename PT>
ScalarType_t<PT> length(PT const & a)
{
	return static_cast<ScalarType_t<PT>>(std::sqrt(static_cast<double>(dot(a, a))));
}

template <typename PT>
inline PT minimum(PT const & a, PT const & b)
{
	if (Dimension_v<PT> == 3)
	{
		PT c = { std::min(a[0], b[0]), std::min(a[1], b[1]), std::min(a[2], b[2]) };
		return c;
	}
	else
	{
		PT c;
		for (auto i = 0; i < Dimension_v<PT>; ++i)
			c[i] = std::min(a[i], b[i]);
		return c;
	}
}

template <typename PT, typename ... ptArgs>
inline PT minimum(PT const & pt, ptArgs const & ... args)
{
	return minimum(pt, minimum(args ...));
}

template <typename PT>
inline PT maximum(PT const & a, PT const & b)
{
	if (Dimension_v<PT> == 3)
	{
		PT c = { std::max(a[0], b[0]), std::max(a[1], b[1]), std::max(a[2], b[2]) };
		return c;
	}
	else
	{
		PT c;
		for (auto i = 0; i < Dimension_v<PT>; ++i)
			c[i] = std::max(a[i], b[i]);
		return c;
	}
}

template <typename PT, typename ... ptArgs>
PT maximum(PT const & pt, ptArgs const & ... args)
{
	return maximum(pt, maximum(args ...));
}

template <typename PT>
PT Epsilon()
{
	PT pt;
	for (auto i = 0; i < Dimension<PT>::value; ++i)
		pt[i] = std::numeric_limits<ScalarType_t<PT>>::epsilon();

	return pt;
}
template <typename PT>
PT Max()
{
	PT pt;
	for (auto i = 0; i < Dimension<PT>::value; ++i)
		pt[i] = std::numeric_limits<ScalarType_t<PT>>::max();

	return pt;
}

template <typename PT>
PT Lowest()
{
	PT pt;
	for (auto i = 0; i < Dimension<PT>::value; ++i)
		pt[i] = std::numeric_limits<ScalarType_t<PT>>::lowest();

	return pt;
}

template <typename ElementType>
struct PointType
{

};

template <typename ElementType>
using PointType_t = typename PointType<ElementType>::type;


#define DECLARE_POINT(SCALAR, DIM)                  \
template <> struct PointTraits<Point<SCALAR, DIM>>  \
{                                                   \
    using type = Point<SCALAR, DIM>;                \
};                                                  \
template <> struct Dimension<Point<SCALAR, DIM>>    \
{                                                   \
    static int constexpr value = DIM;               \
};                                                  \
template <> struct IsPoint<Point<SCALAR, DIM>>      \
{                                                   \
    static bool constexpr value = true;           \
};

DECLARE_POINT(double, 3)
DECLARE_POINT(float, 3)
DECLARE_POINT(double, 2)

using Pointd2 = Point<double, 2>;
using Pointd3 = Point<double, 3>;
using Pointf3 = Point<float, 3>;
