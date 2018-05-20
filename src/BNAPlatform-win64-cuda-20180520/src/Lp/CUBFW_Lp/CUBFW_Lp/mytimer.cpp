#include "mytimer.h"
#include <stdio.h>
#include <time.h>

char *TimerGetLocalTime(char *fmt)
{
	time_t tp;
	struct tm *t;

	if (NULL == fmt)
	{
		return NULL;
	}
	
	time(&tp);
	t = localtime(&tp);
	sprintf(fmt, "%04d-%02d-%02d %02d:%02d:%02d", \
		1900+t->tm_year, 1+t->tm_mon, t->tm_mday, t->tm_hour, t->tm_min, t->tm_sec);

	return fmt;
}

int TimerInit(STimer *tmr)
{
	if (NULL == tmr)
	{
		return -1;
	}

#ifdef _WIN32
	tmr->start = 0;
	tmr->stop = 0;
	QueryPerformanceFrequency((LARGE_INTEGER *)&(tmr->freq));
	if (0 == tmr->freq) return -1;
#else
	tmr->start.tv_sec = 0;
	tmr->stop.tv_sec = 0;
	tmr->start.tv_usec = 0;
	tmr->stop.tv_usec = 0;
#endif

	return 0;
}

int TimerStart(STimer *tmr)
{
	if (NULL == tmr)
	{
		return -1;
	}

#ifdef _WIN32
	QueryPerformanceCounter((LARGE_INTEGER *)&(tmr->start));
#else
	gettimeofday(&(tmr->start), NULL);
#endif

	return 0;
}

int TimerStop(STimer *tmr)
{
	if (NULL == tmr)
	{
		return -1;
	}

#ifdef _WIN32
	QueryPerformanceCounter((LARGE_INTEGER *)&(tmr->stop));
#else
	gettimeofday(&(tmr->stop), NULL);
#endif

	return 0;
}

double TimerGetRuntime(STimer *tmr)
{
	if (NULL == tmr)
	{
		return -1.0;
	}

#ifdef _WIN32
	return ((double)(tmr->stop)-(double)(tmr->start))/(tmr->freq);
#else
	return (tmr->stop.tv_sec-tmr->start.tv_sec) + \
		(tmr->stop.tv_usec-tmr->start.tv_usec)/1000000.;
#endif

}
