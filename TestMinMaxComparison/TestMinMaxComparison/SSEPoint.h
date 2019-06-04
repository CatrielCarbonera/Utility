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

struct __declspec(align(16)) SSEPointd3
{
	double m_Pt[3];
};

inline SSEPointd3 make_point(double x, double y, double z)
{
	return { x, y, z };
}

inline const SSEPointd3 operator+ (SSEPointd3 &pt)
{
	return pt;
}

inline const SSEPointd3 operator++ (SSEPointd3 &pt)
{
	constexpr SSEPointd3 one{ 1, 1, 1 };
	SSEPointd3 result;
	reinterpret_cast<__m128d &> (result.m_Pt) = _mm_add_pd((__m128d&) one.m_Pt, (__m128d&) pt.m_Pt);
	result.m_Pt[2] = pt.m_Pt[2] + 1;
}

inline const SSEPointd3 operator- (SSEPointd3 &pt)
{
	SSEPointd3 result;
	result.m_Pt[0] = -pt.m_Pt[0];
	result.m_Pt[1] = -pt.m_Pt[1];
	result.m_Pt[2] = -pt.m_Pt[2];
	return result;
}

inline const SSEPointd3 operator-- (SSEPointd3 &pt)
{
	constexpr SSEPointd3 one{ 1, 1, 1 };
	SSEPointd3 result;
	reinterpret_cast<__m128d &> (result.m_Pt) = _mm_sub_pd((__m128d&) one.m_Pt, (__m128d&) pt.m_Pt);
	return result;
}

inline const SSEPointd3 operator+ (const SSEPointd3 &ptA, const SSEPointd3 &ptB)
{
	SSEPointd3 result;

	reinterpret_cast<__m128d &> (result.m_Pt) = _mm_add_pd((__m128d&) ptA.m_Pt, (__m128d&) ptB.m_Pt);
	result.m_Pt[2] = ptA.m_Pt[2] + ptB.m_Pt[2];

	return result;
}

inline SSEPointd3 &operator+= (SSEPointd3 &ptA, const SSEPointd3 &ptB)
{
	reinterpret_cast<__m128d &> (ptA.m_Pt) = _mm_add_pd((__m128d&) ptA.m_Pt, (__m128d&) ptB.m_Pt);
	ptA.m_Pt[2] = ptA.m_Pt[2] + ptB.m_Pt[2];

	return ptA;
}

inline const SSEPointd3 operator- (const SSEPointd3 &ptA, const SSEPointd3 &ptB)
{
	SSEPointd3 result;

	reinterpret_cast<__m128d &> (result.m_Pt) = _mm_sub_pd((__m128d&) ptA.m_Pt, (__m128d&) ptB.m_Pt);
	result.m_Pt[2] = ptA.m_Pt[2] - ptB.m_Pt[2];

	return result;
}

inline SSEPointd3 &operator-= (SSEPointd3 &ptA, const SSEPointd3 &ptB)
{
	reinterpret_cast<__m128d &> (ptA.m_Pt) = _mm_sub_pd((__m128d&) ptA.m_Pt, (__m128d&) ptB.m_Pt);
	ptA.m_Pt[2] = ptA.m_Pt[2] - ptB.m_Pt[2];

	return ptA;
}

/// Float version
using SSEPointf3 = __m128;

inline const SSEPointf3 operator+ (SSEPointf3 &pt)
{
	return pt;
}

inline const SSEPointf3 operator++ (SSEPointf3 &pt)
{
	constexpr SSEPointf3 one{ 1, 1, 1, 1 };
	return _mm_add_ps(one, pt);
}

inline const SSEPointf3 operator- (SSEPointf3 &pt)
{
	constexpr SSEPointf3 zero{ 0, 0, 0, 0 };
	return _mm_sub_ps(zero, pt);
}

inline const SSEPointf3 operator-- (SSEPointf3 &pt)
{
	constexpr SSEPointf3 one{ 1.f, 1.f, 1.f, 1.f };

	pt = _mm_sub_ps(one, pt);
	return pt;
}

inline const SSEPointf3 operator+ (const SSEPointf3 &ptA, const SSEPointf3 &ptB)
{
	return _mm_add_ps(ptA, ptB);
}

inline SSEPointf3 &operator+= (SSEPointf3 &ptA, const SSEPointf3 &ptB)
{

	ptA = _mm_add_ps(ptA, ptB);

	return ptA;
}

inline const SSEPointf3 operator- (const SSEPointf3 &ptA, const SSEPointf3 &ptB)
{
	return _mm_sub_ps(ptA, ptB);
}

inline SSEPointf3 &operator-= (SSEPointf3 &ptA, const SSEPointf3 &ptB)
{
	ptA = _mm_sub_ps(ptA, ptB);

	return ptA;
}

inline const SSEPointf3 operator* (float a, const SSEPointf3 &ptB)
{
	return _mm_mul_ps(_mm_set_ps1(a), ptB);
}

inline const SSEPointf3 operator* (const SSEPointf3 &ptA, float b)
{
	return _mm_mul_ps(_mm_set_ps1(b), ptA);
}

inline SSEPointf3 &operator*= (SSEPointf3 &ptA, float b)
{
	ptA = _mm_mul_ps(_mm_set_ps1(b), ptA);

	return ptA;
}

inline const SSEPointf3 operator/ (const SSEPointf3 &ptA, float b)
{
	return _mm_div_ps(ptA, _mm_set_ps1(b));
}

inline SSEPointf3 &operator/= (SSEPointf3 &ptA, float b)
{
	ptA = _mm_div_ps(ptA, _mm_set_ps1(b));

	return ptA;
}

