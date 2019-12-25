/*
*
* Copyright (C) 2018 Carlos Carbonera - All Rights Reserved
* You may use this code with explicit consent from the author.
* It is illegal to reproduce or modify this code without the
* author's explicit consent.
*
*/

#include "Timer.h"
#include "Point.h"
#include "SSEPoint.h"
#include <ppl.h>
#include <omp.h>

template <typename PointType>
double Performance(
	std::string title,
	std::vector<PointType> &ptA,
	const std::vector<PointType> &ptB,
	std::vector<PointType> &ptC
	)
{
	Timer timer(title);
	using scalar_type = ScalarType_t<PointType>;

	auto ptrA = ptA.data();
	auto ptrB = ptB.data();
	auto ptrC = ptC.data();
#define OMP_DEF 0
#if PARALLEL // much slower!
	concurrency::parallel_for(0, static_cast<int> (ptA.size()),
		[&](int i)
	{
		ptrC[i] = saxpy(static_cast<scalar_type>(2.125), ptrA[i], ptrB[i]);
		ptrC[i] += ptrB[i];
		ptrC[i] = saxpy(static_cast<scalar_type>(2.125), (ptrA[i] + ptrB[i]), ptrC[i]);
	});
#elif OMP_DEF // much slower!
#pragma omp parallel for schedule(static, 16)
	for (auto i = 0; i < static_cast<int>(ptA.size()); ++i)
	{
		ptrC[i] = saxpy(static_cast<scalar_type>(2.125), ptrA[i], ptrB[i]);
		ptrC[i] += ptrB[i];
		ptrC[i] = saxpy(static_cast<scalar_type>(2.125), (ptrA[i] + ptrB[i]), ptrC[i]);
	}
#else
	for (auto i = 0; i <ptA.size(); ++i)
	{
#ifdef ORIG
		ptrC[i] = ptrA[i];
		ptrC[i] *= static_cast<scalar_type>(2.125);
		ptrC[i] += ptrB[i];
		ptrC[i] += ptrB[i];
		ptrC[i] = saxpy (static_cast<scalar_type>(2.125), (ptrA[i] + ptrB[i]), ptrC[i]);
#else
		ptrC[i] = saxpy(static_cast<scalar_type>(2.125), ptrA[i], ptrB[i]);
		ptrC[i] += ptrB[i];
		ptrC[i] = saxpy(static_cast<scalar_type>(2.125), (ptrA[i] + ptrB[i]), ptrC[i]);
#endif
	}
#endif
	return timer.elapsed();
}

constexpr int N = 10000000;

void TestSSEFloat()
{
	std::vector<Pointf3> ptA;
	std::vector<Pointf3> ptB;
	std::vector<Pointf3> ptC;
	std::vector<SSEPointf3> ptASSE;
	std::vector<SSEPointf3> ptBSSE;
	std::vector<SSEPointf3> ptCSSE;

	for (float i = 0; i < static_cast<float>(N); ++i)
	{
		ptA.push_back(Pointf3{ i, i + 1, i + 2 });
		ptB.push_back(Pointf3{ i + 1, i + 2, i + 3 });
		ptC.push_back(Pointf3{});
		ptASSE.push_back(SSEPointf3{ i, i + 1, i + 2, 0 });
		ptBSSE.push_back(SSEPointf3{ i + 1, i + 2, i + 3, 0 });
		ptCSSE.push_back(SSEPointf3{});
	}

	const auto timePointf3 = Performance<Pointf3>("Pointf3", ptA, ptB, ptC);
	const auto timeSSEPoint3f = Performance<SSEPointf3>("SSEPoint3f", ptASSE, ptBSSE, ptCSSE);

	std::cout << "Percentage improvement = " << static_cast<int>(100 * (timePointf3 - timeSSEPoint3f) / timePointf3) << "%" << std::endl;

	float maxError = 0.f;

	for (auto i = 0; i < ptA.size(); ++i)
	{
		for (auto j = 0; j < 3; ++j)
		{
			const auto diff = std::abs(ptC[i][j] - ptCSSE[i].m128_f32[j]);
			if (diff > maxError)
				maxError = diff;
		}
	}

	std::cout << "Maximum Error = " << maxError << std::endl;
}

void TestSSEDouble()
{
	std::vector<Pointd3> ptA;
	std::vector<Pointd3> ptB;
	std::vector<Pointd3> ptC;
	std::vector<SSEPointd3> ptASSE;
	std::vector<SSEPointd3> ptBSSE;
	std::vector<SSEPointd3> ptCSSE;

	for (float i = 0; i < static_cast<float>(N); ++i)
	{
		ptA.push_back(Pointd3{ i, i + 1, i + 2 });
		ptB.push_back(Pointd3{ i + 1, i + 2, i + 3 });
		ptC.push_back(Pointd3{});
		ptASSE.push_back(SSEPointd3{ i, i + 1, i + 2 });
		ptBSSE.push_back(SSEPointd3{ i + 1, i + 2, i + 3 });
		ptCSSE.push_back(SSEPointd3{});
	}

	const auto timePointf3 = Performance<Pointd3>("Pointd3", ptA, ptB, ptC);
	const auto timeSSEPoint3f = Performance<SSEPointd3>("SSEPointd3", ptASSE, ptBSSE, ptCSSE);

	std::cout << "Percentage improvement = " << static_cast<int>(100 * (timePointf3 - timeSSEPoint3f) / timePointf3) << "%" << std::endl;

	double maxError = 0;
	for (auto i = 0; i < ptC.size(); ++i)
	{
		for (auto j = 0; j < 3; ++j)
		{
			const auto diff = std::abs(ptC[i][j] - ptCSSE[i].m256d_f64[j]);
			if (diff > maxError)
				maxError = diff;
		}
	}

	std::cout << "Maximum Error = " << maxError << std::endl;
}
void TestSSE ()
{
	TestSSEDouble();
	TestSSEFloat();
}