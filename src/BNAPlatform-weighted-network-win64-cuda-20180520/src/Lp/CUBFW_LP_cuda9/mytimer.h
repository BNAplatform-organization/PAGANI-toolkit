#ifndef __TIMER_INCLUDE__
#define __TIMER_INCLUDE__

#ifdef _WIN32 /*for windows*/

#define WIN32_LEAN_AND_MEAN
#include <Windows.h>

typedef struct tagSTimer
{
	__int64 freq;
	__int64 start;
	__int64 stop;
} STimer, *LPSTimer;

#else /*for linux*/

#include <sys/time.h>
#include <time.h>

typedef struct tagSTimer
{
	struct timeval start;
	struct timeval stop;
} STimer, *LPSTimer;

#endif

#ifdef __cplusplus
extern "C" {
#endif

char	*TimerGetLocalTime(char *fmt);	/*the length of fmt must larger than 20*/
int		TimerInit(STimer *tmr);
int		TimerStart(STimer *tmr);
int		TimerStop(STimer *tmr);
double	TimerGetRuntime(STimer *tmr);	/*in second*/

#ifdef __cplusplus
}
#endif

#endif
