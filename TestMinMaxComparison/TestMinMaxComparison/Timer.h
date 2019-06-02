
#pragma once
#include <chrono>
#include <iostream>
#include <string>

struct Timer
{
	using clock = std::chrono::high_resolution_clock;
	using time_point = std::chrono::time_point<clock>;

	std::string title;
	time_point start;
	Timer(std::string title)
		: title(title),
		start(clock::now())
	{

	}

	~Timer()
	{
		std::cout << title << " elapsed time: " << std::chrono::duration<double>(clock::now() - start).count() << std::endl;
	}
};
