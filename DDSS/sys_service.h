#ifndef _SYS_SERVICE_H
#define _SYS_SERVICE_H

#ifndef NO_LINUX

#ifdef __cplusplus
extern "C"{
#endif
 int get_memory_stat(size_t *total_ram, size_t *free_ram, size_t *used_swap);
 double accu_time(void);
 double system_clock(void);
#ifdef __cplusplus
}
#endif

#endif

#endif
