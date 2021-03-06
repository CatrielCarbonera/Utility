/*
*
* Copyright (C) 2018 Carlos Carbonera - All Rights Reserved
* You may use this code with explicit consent from the author.
* It is illegal to reproduce or modify this code without the
* author's explicit consent.
*
*/

#include <array>
#include <functional>
#include <iostream>
#include <random>
#include <vector>
#include <ppl.h>

#include <xmmintrin.h>
#include <emmintrin.h>
#include <smmintrin.h>

#include "Timer.h"

template <int Dim, typename Data>
using VectorBase = std::array<Data, Dim>;
using Vector3Df = VectorBase<3, float>;

void TestSSE();

struct P3Df
{
	using value_type = float;
	std::array<float, 3> pt;

	float  operator[] (int i) const { return pt[i]; }
	float &operator[] (int i) { return pt[i]; }

	constexpr std::size_t size() const { return pt.size(); }
};

template <typename PT>
struct BoundingBox
{
	PT m_Min;
	PT m_Max;
};

template <typename PT>
BoundingBox<PT> MakeBoundingBox(const PT&ptA, const PT &ptB)
{
	BoundingBox<PT> result{ ptA, ptA };
	for (auto i = 0; i < ptB.size(); ++i)
	{
		if (result.m_Min[i] > ptB[i])
			result.m_Min[i] = ptB[i];

		if (result.m_Max[i] < ptB[i])
			result.m_Max[i] = ptB[i];
	}

	return result;
};

struct UV3Df 
{
	__m128 m_Pt;
	float operator[] (int i) const
	{
		return m_Pt.m128_f32[i];
	}
	float &operator[] (int i)
	{
		return m_Pt.m128_f32[i];
	}
};

template <>
struct BoundingBox<Vector3Df>
{
	__m128 m_Min;
	__m128 m_Max;
};

using BBX3VDf = BoundingBox<Vector3Df>;

template <>
BBX3VDf MakeBoundingBox<Vector3Df>(const Vector3Df &ptA, const Vector3Df &ptB)
{
	BBX3VDf result{ };

	for (auto i = 0; i < 3; ++i)
		result.m_Max.m128_f32[i] = std::max(ptA[i], ptB[i]);
	result.m_Max.m128_f32[3] = 0.f;

	for (auto i = 0; i < 3; ++i)
		result.m_Min.m128_f32[i] = std::min(ptA[i], ptB[i]);
	result.m_Min.m128_f32[3] = 0.f;

	return result;
}

using AltBBxP3Df = BoundingBox<P3Df>;


bool Intersect(const AltBBxP3Df &bbx0, const AltBBxP3Df &bbx1)
{
	return
		std::max(bbx0.m_Min[0], bbx1.m_Min[0]) <= std::min(bbx0.m_Max[0], bbx1.m_Max[0]) &&
		std::max(bbx0.m_Min[1], bbx1.m_Min[1]) <= std::min(bbx0.m_Max[1], bbx1.m_Max[1]) &&
		std::max(bbx0.m_Min[2], bbx1.m_Min[2]) <= std::min(bbx0.m_Max[2], bbx1.m_Max[2]);
}

bool Intersect(const BBX3VDf &bbx0, const BBX3VDf &bbx1)
{
	constexpr __m128i ref_mask = {
		(short)-1, (short)-1, (short)-1, (short)-1,
		(short)-1, (short)-1, (short)-1, (short)-1,
		(short)-1, (short)-1, (short)-1, (short)-1,
		(short)-1, (short)-1, (short)-1, (short)-1,
	};
	auto const mmmin = _mm_max_ps(bbx0.m_Min, bbx1.m_Min);
	auto const mmmax = _mm_min_ps(bbx0.m_Max, bbx1.m_Max);
	auto const lte = _mm_castps_si128(_mm_cmple_ps(mmmin, mmmax));
	return _mm_testc_si128(lte, ref_mask);
}

template <typename PT>
PT GenRandomPoint()
{
	std::random_device rd;
	std::uniform_real_distribution<typename PT::value_type> dist(-1., 1.);

	PT pt;
	for (auto i = 0; i < pt.size(); ++i)
		pt[i] = dist(rd);

	return pt;
}

template <bool DoParallel, typename BBX>
int test(std::string title, std::vector<BBX> & bbxA, const std::vector<BBX> &bbxB)
{
	Timer timer(title);

	if constexpr (DoParallel)
	{
		concurrency::combinable<int> counter;
		concurrency::parallel_for(0, (int) bbxA.size(), 
			[&](int i)
		{
			if (Intersect(bbxA[i], bbxB[i]))
				counter.local() += Intersect(bbxA[i], bbxB[i]);
		});

		return counter.combine(std::plus<int>());
	}
	else
	{
		int count = 0;
		for (auto i = 0; i < bbxA.size(); ++i)
		{
			count += Intersect(bbxA[i], bbxB[i]) ? 1 : 0;
		}

		return count;
	}
}

/*
num intersects = SSE elapsed time:    0.000857855
29539
num intersects = Serial elapsed time: 0.00171457
29539
*/
int main()
{
	TestSSE();

	constexpr int NUM = 1000000;
	std::vector<BBX3VDf> bbxA;
	std::vector<BBX3VDf> bbxB;
	std::vector<AltBBxP3Df> bbxA1;
	std::vector<AltBBxP3Df> bbxB1;

	for (auto i = 0; i < NUM; ++i)
	{
		const auto ptA = GenRandomPoint<Vector3Df>();
		const auto ptB = GenRandomPoint<Vector3Df>();
		bbxA.push_back(MakeBoundingBox(ptA, ptB));

		const auto ptC = GenRandomPoint<Vector3Df>();
		const auto ptD = GenRandomPoint<Vector3Df>();
		bbxB.push_back(MakeBoundingBox(ptC, ptD));

		const P3Df ptA1{ ptA[0], ptA[1], ptA[2] };
		const P3Df ptB1{ ptB[0], ptB[1], ptB[2] };
		bbxA1.push_back(MakeBoundingBox(ptA1, ptB1));

		const P3Df ptC1{ ptC[0], ptC[1], ptC[2] };
		const P3Df ptD1{ ptD[0], ptD[1], ptD[2] };
		bbxB1.push_back(MakeBoundingBox(ptC1, ptD1));
	}


	std::cout << "num intersects = " << test<false>("Serial", bbxA1, bbxB1) << std::endl;
	std::cout << "num intersects = " << test<false>("SSE", bbxA, bbxB) << std::endl;

	return 0;
}

