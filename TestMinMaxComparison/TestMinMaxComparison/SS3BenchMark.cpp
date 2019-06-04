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

template <typename PointType>
void Performance(
	std::string title,
	std::vector<PointType> &ptA,
	const std::vector<PointType> &ptB,
	std::vector<PointType> &ptC
	)
{
	Timer timer(title);

	auto ptrA = ptA.data();
	auto ptrB = ptB.data();
	auto ptrC = ptC.data();

	for (auto i = 0; i < ptA.size(); ++i)
	{
		ptrC[i] = (2.125f * ptrA[i]) + ptrB[i];
		ptrC[i] += ptrB[i];
		ptrC[i] += (2.125f * ptrA[i]) + ptrB[i];
	}

	return;
}


void TestSSEFloat()
{
	constexpr int N = 1000000;

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

	Performance<Pointf3>("Pointf3", ptA, ptB, ptC);
	Performance<SSEPointf3>("SSEPoint3f", ptASSE, ptBSSE, ptCSSE);

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
	constexpr int N = 1000000;
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

	//Performance<Pointd3>("Pointd3", ptA, ptB, ptC);
	//Performance<SSEPointd3>("SSEPointd3", ptASSE, ptBSSE, ptCSSE);

	double maxError = 0;
	for (auto i = 0; i < ptC.size(); ++i)
	{
		for (auto j = 0; j < 3; ++j)
		{
			const auto diff = std::abs(ptC[i][j] - ptCSSE[i].m_Pt[j]);
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