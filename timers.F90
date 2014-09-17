       module timers
!Timing services (OpenMP omp_get_wtime() based).
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2014/09/16
!FUNCTIONS:
! # integer timer_start(real8:time_set, integer:time_handle);
! # logical time_is_off(integer:time_handle, integer:ierr[, logical:destroy]);
! # integer timer_destroy(integer:time_handle);
! # real8 timer_tick_sec();
! # real8 thread_wtime([real8:tbase]);

!PARAMETERS:
        integer, parameter, private:: max_timers=8192
        integer, parameter, public:: timers_err_invalid_arg=1
        integer, parameter, public:: timers_err_no_timers_left=2
        integer, parameter, public:: timers_err_timer_null=3
!TYPES:
        type, private:: timer_t
         real(8), private:: beg_time      !time the timer started (sec)
         real(8), private:: time_interval !time the timer is set for (sec)
        end type timer_t
!DATA:
        integer, private:: j_
        type(timer_t), private:: timer(0:max_timers-1)=(/(timer_t(-1d0,-1d0),j_=0,max_timers-1)/)
        integer, private:: handle_stack(0:max_timers-1)=(/(j_,j_=0,max_timers-1)/)
        integer, private:: handle_sp=0
        real(8), private:: timer_tick=-1d0
        real(8), external, private:: omp_get_wtime,omp_get_wtick

       contains
!---------------------------------------------------------
        integer function timer_start(time_set,time_handle)
!This function sets up a timer limited to <time_set> seconds and returns its handle in <time_handle>.
        implicit none
        real(8), intent(in):: time_set     !requested time in microseconds
        integer, intent(out):: time_handle !timer handle
        real(8):: val
        if(time_set.ge.0d0) then
         if(handle_sp.ge.0.and.handle_sp.lt.max_timers) then
          time_handle=handle_stack(handle_sp); handle_sp=handle_sp+1
          val=omp_get_wtime(); timer(time_handle)=timer_t(val,time_set)
          timer_start=0
         else
          timer_start=timers_err_no_timers_left
         endif
        else
         timer_start=timers_err_invalid_arg
        endif
        return
        end function timer_start
!-------------------------------------------------------------
        logical function time_is_off(time_handle,ierr,destroy)
!This function tests whether a given timer has expired.
!If <destroy> is present and .true., timer handle will be destroyed if the timer has expired.
        implicit none
        integer, intent(inout):: time_handle !timer handle
        integer, intent(inout):: ierr
        logical, intent(in), optional:: destroy
        real(8):: tm
        time_is_off=.false.
        if(time_handle.ge.0.and.time_handle.lt.max_timers) then !valid range
         if(timer(time_handle)%time_interval.ge.0d0) then !valid handle
          ierr=0; tm=omp_get_wtime()
          if(tm.ge.timer(time_handle)%beg_time+timer(time_handle)%time_interval) time_is_off=.true.
          if(time_is_off.and.present(destroy)) then
           if(destroy) then
            timer(time_handle)=timer_t(-1d0,-1d0)
            handle_sp=handle_sp-1; handle_stack(handle_sp)=time_handle
           endif
          endif
         else
          ierr=timers_err_timer_null
         endif
        else
         ierr=timers_err_invalid_arg
        endif
        return
        end function time_is_off
!--------------------------------------------------
        integer function timer_destroy(time_handle)
!This function frees a time handle.
        implicit none
        integer, intent(in):: time_handle
        timer_destroy=0
        if(time_handle.ge.0.and.time_handle.lt.max_timers) then !valid range
         if(timer(time_handle)%time_interval.ge.0d0) then !valid handle
          timer(time_handle)=timer_t(-1d0,-1d0)
          handle_sp=handle_sp-1; handle_stack(handle_sp)=time_handle
         else
          timer_destroy=timers_err_timer_null
         endif
        else
         timer_destroy=timers_err_invalid_arg
        endif
        return
        end function timer_destroy
!----------------------------------------
        real(8) function timer_tick_sec()
!This function returns the wall clock tick length in seconds.
        if(timer_tick.le.0d0) timer_tick=omp_get_wtick()
        timer_tick_sec=timer_tick
        return
        end function timer_tick_sec
!-------------------------------------------
!DIR$ ATTRIBUTES OFFLOAD:mic:: thread_wtime
        real(8) function thread_wtime(tbase)
!This function returns the current wall clock time in seconds.
        implicit none
        real(8), intent(in), optional:: tbase
        real(8) tm
#ifndef NO_OMP
        thread_wtime=omp_get_wtime()
#else
#ifdef USE_GNU
        thread_wtime=real(secnds(0.),8)
#else
        call cpu_time(tm); thread_wtime=tm
#endif
#endif
        if(present(tbase)) thread_wtime=thread_wtime-tbase
        return
        end function thread_wtime

       end module timers
