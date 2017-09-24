/*
 * Copyright 2016 Grzegorz Kozikowski, University of Manchester. All rights reserved.
 */

#ifndef CPU_TIMER_WIN32_H
#define CPU_TIMER_WIN32_H

#include <windows.h>
#include <iostream>
using std::cout;

class CPUTimer
{
public:
	CPUTimer();
	~CPUTimer();

	void startTime();
	unsigned int calculateElapsedTime();
private:
	LARGE_INTEGER timerFreq_;
	LARGE_INTEGER counterAtStart_;
};

CPUTimer::CPUTimer() {};
CPUTimer::~CPUTimer() {};

void CPUTimer::startTime()
{
	QueryPerformanceFrequency(&timerFreq_);
	QueryPerformanceCounter(&counterAtStart_);
	TIMECAPS ptc;
	UINT cbtc = 8;
	MMRESULT result = timeGetDevCaps(&ptc, cbtc);
}

unsigned int CPUTimer::calculateElapsedTime()
{
	if (timerFreq_.QuadPart == 0)
	{
		return -1;
	}
	else
	{
		LARGE_INTEGER c;
		QueryPerformanceCounter(&c);
		return static_cast<unsigned int>( (c.QuadPart - counterAtStart_.QuadPart) * 1000 / timerFreq_.QuadPart );
	}
}
#endif
