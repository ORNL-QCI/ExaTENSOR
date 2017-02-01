/** This file contains Linux OS wrapped calls. **/

#ifndef NO_LINUX

#include <stddef.h>
#include <stdlib.h>
#include <sys/sysinfo.h>
#include <sys/time.h>
#include <time.h>

#include "sys_service.h"

int get_memory_stat(size_t *total_ram, size_t *free_ram, size_t *used_swap){ //not very accurate
 int i;
 struct sysinfo info;
 i=sysinfo(&info);
 if(i==0){
  *total_ram=(size_t)(info.totalram*info.mem_unit);
  *free_ram=(size_t)(info.freeram*info.mem_unit);
  *used_swap=(size_t)((info.totalswap-info.freeswap)*info.mem_unit);
 }else{
  return 1;
 };
 return 0;
}

double accu_time(void){
 struct timeval timer;
 if(gettimeofday(&timer,NULL)) return -1.0;
 return (((double)timer.tv_sec)+((double)timer.tv_usec)*(1.0E-6));
}

/*
double system_clock()
{
 struct timespec tp;
 if(clock_gettime(CLOCK_MONOTONIC,&tp)) return -1.0;
 return (((double)tp.tv_sec)+((double)tp.tv_nsec)*(1.0E-9));
}
*/

#endif
