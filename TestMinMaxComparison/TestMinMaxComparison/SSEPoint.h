/*
*
* Copyright (C) 2018 Carlos Carbonera - All Rights Reserved
* You may use this code with explicit consent from the author.
* It is illegal to reproduce or modify this code without the
* author's explicit consent.
*
*/

#pragma once

#include <array>
#include <functional>
#include <iostream>
#include <random>
#include <vector>
#include <ppl.h>

#include <xmmintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>

using SSEPointd3 = __m256d;

constexpr SSEPointd3 OneD = { 1, 1, 1, 1 };
constexpr SSEPointd3 ZeroD = { 0,0,0,0 };

inline const SSEPointd3 _vectorcall operator+ (SSEPointd3 &pt)
{
	return pt;
}

inline const SSEPointd3 &_vectorcall operator++ (SSEPointd3 &pt)
{
	pt = _mm256_add_pd(OneD, pt);
	return pt;
}


inline const SSEPointd3 operator+ (const SSEPointd3 &ptA, const SSEPointd3 &ptB)
{
	return _mm256_add_pd(ptA, ptB);
}

inline SSEPointd3 &_vectorcall operator+= (SSEPointd3 &ptA, const SSEPointd3 &ptB)
{
	ptA = _mm256_add_pd(ptA, ptB);

	return ptA;
}

inline const SSEPointd3 _vectorcall operator- (SSEPointd3 &pt)
{
	return _mm256_sub_pd(ZeroD, pt);

}
inline const SSEPointd3 &_vectorcall operator-- (SSEPointd3 &pt)
{
	pt = _mm256_sub_pd(OneD, pt);
	return pt;
}

inline const SSEPointd3 _vectorcall operator- (const SSEPointd3 &ptA, const SSEPointd3 &ptB)
{
	return _mm256_sub_pd(ptA, ptB);
}

inline SSEPointd3 &_vectorcall operator-= (SSEPointd3 &ptA, const SSEPointd3 &ptB)
{
	ptA = _mm256_sub_pd(ptA, ptB);
	return ptA;
}

inline const SSEPointd3 _vectorcall operator* (float a, const SSEPointd3 &ptB)
{
	return _mm256_mul_pd(_mm256_set1_pd(a), ptB);
}

inline const SSEPointd3 _vectorcall operator* (const SSEPointd3 &ptA, float b)
{
	return _mm256_mul_pd(_mm256_set1_pd(b), ptA);
}

inline SSEPointd3 &_vectorcall operator*= (SSEPointd3 &ptA, float b)
{
	ptA = _mm256_mul_pd(_mm256_set1_pd(b), ptA);

	return ptA;
}

inline const SSEPointd3 _vectorcall operator/ (const SSEPointd3 &ptA, float b)
{
	return _mm256_div_pd(ptA, _mm256_set1_pd(b));
}

inline SSEPointd3 &_vectorcall operator/= (SSEPointd3 &ptA, float b)
{
	ptA = _mm256_div_pd(ptA, _mm256_set1_pd(b));

	return ptA;
}

/// Float version
using SSEPointf3 = __m128;

inline const SSEPointf3 _vectorcall operator+ (SSEPointf3 &pt)
{
	return pt;
}

inline const SSEPointf3 &_vectorcall operator++ (SSEPointf3 &pt)
{
	constexpr SSEPointf3 one{ 1, 1, 1, 1 };
	pt = _mm_add_ps(one, pt);

	return pt;
}

inline const SSEPointf3 _vectorcall operator- (SSEPointf3 &pt)
{
	constexpr SSEPointf3 zero{ 0, 0, 0, 0 };
	return _mm_sub_ps(zero, pt);
}

inline const SSEPointf3 &_vectorcall operator-- (SSEPointf3 &pt)
{
	constexpr SSEPointf3 one{ 1.f, 1.f, 1.f, 1.f };

	pt = _mm_sub_ps(one, pt);
	return pt;
}

inline const SSEPointf3 _vectorcall operator+ (const SSEPointf3 &ptA, const SSEPointf3 &ptB)
{
	return _mm_add_ps(ptA, ptB);
}

inline SSEPointf3 &_vectorcall operator+= (SSEPointf3 &ptA, const SSEPointf3 &ptB)
{

	ptA = _mm_add_ps(ptA, ptB);

	return ptA;
}

inline const SSEPointf3 _vectorcall operator- (const SSEPointf3 &ptA, const SSEPointf3 &ptB)
{
	return _mm_sub_ps(ptA, ptB);
}

inline SSEPointf3 &_vectorcall operator-= (SSEPointf3 &ptA, const SSEPointf3 &ptB)
{
	ptA = _mm_sub_ps(ptA, ptB);

	return ptA;
}

inline const SSEPointf3 _vectorcall operator* (float a, const SSEPointf3 &ptB)
{
	return _mm_mul_ps(_mm_set_ps1(a), ptB);
}

inline const SSEPointf3 _vectorcall operator* (const SSEPointf3 &ptA, float b)
{
	return _mm_mul_ps(_mm_set_ps1(b), ptA);
}

inline SSEPointf3 &_vectorcall operator*= (SSEPointf3 &ptA, float b)
{
	ptA = _mm_mul_ps(_mm_set_ps1(b), ptA);

	return ptA;
}

inline const SSEPointf3 _vectorcall operator/ (const SSEPointf3 &ptA, float b)
{
	return _mm_div_ps(ptA, _mm_set_ps1(b));
}

inline SSEPointf3 &_vectorcall operator/= (SSEPointf3 &ptA, float b)
{
	ptA = _mm_div_ps(ptA, _mm_set_ps1(b));

	return ptA;
}

template <>
struct ScalarType<SSEPointd3>
{
	using value_type = double;
};
template <>
struct ScalarType<SSEPointf3>
{
	using value_type = float;
};

