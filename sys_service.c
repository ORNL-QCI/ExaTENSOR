/** This file contains Linux OS wrapped calls. **/

#include <stdlib.h>
#include <sys/sysinfo.h>

#ifdef __cplusplus
extern "C"{
#endif
 int get_memory_stat(size_t *total_ram, size_t *free_ram, size_t *used_swap);
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
