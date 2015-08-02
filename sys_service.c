/** This file contains Linux OS wrapped calls. **/

#include <stddef.h>
#include <stdlib.h>
#include <sys/sysinfo.h>
#include <sys/time.h>

#ifdef __cplusplus
extern "C"{
#endif
 int get_memory_stat(size_t *total_ram, size_t *free_ram, size_t *used_swap);
 double accu_time(void);
#ifdef __cplusplus
}
#endif

int get_memory_stat(size_t *total_ram, size_t *free_ram, size_t *used_swap){
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
 return (timer.tv_sec+timer.tv_usec*(1.0E-6));
}
