       module timers
!Timing services.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2014/07/11

!PARAMETERS:
        integer, parameter, private:: max_timers=8192
!TYPES:
        type, private:: timer_t
         real(8), private:: beg_time      !time the timer started (sec)
         real(8), private:: time_interval !time the timer is set for (sec)
        end type timer_t
!DATA:
        integer, private:: j_
        type(timer_t), private:: timers(0:max_timers-1)=timer_t(-1d0,-1d0)
        integer, private:: handle_stack(0:max_timers-1)=(/(j_,j_=0,max_timers-1)/)
        integer, private:: handle_sp=0
        real(8), private:: timer_tick=0d0
        real(8), external, private:: omp_get_wtime,omp_get_wtick

       contains
!---------------------------------------------------------
        integer function timer_start(time_set,time_handle)
!This function sets up a timer limited to <time_set> seconds and returns its handle in <time_handle>.
        implicit none
        real(8), intent(in):: time_set     !requested time in microseconds
        integer, intent(out):: time_handle !timer handle
        if(time_set.ge.0d0) then
         if(handle_sp.ge.0.and.handle_sp.lt.max_timers) then
          time_handle=handle_stack(handle_sp); handle_sp=handle_sp+1
          timers(time_handle)=timer_t(omp_get_wtime(),time_set)
          timer_start=0
         else
          timer_start=1
         endif
        else
         timer_start=2
        endif
        return
        end function timer_start
!-----------------------------------------------------
        logical function time_is_off(time_handle,ierr)
!This function tests whether a given timer has expired. Note that if the timer
!has expired, its handle will be immediately destroyd here (you cannot use it after that).
        implicit none
        integer, intent(inout):: time_handle !timer handle
        integer, intent(inout):: ierr
        real(8):: tm
        time_is_off=.false.
        if(time_handle.ge.0.and.time_handle.lt.max_timers) then !valid range
         if(timers(time_handle)%time_interval.ge.0d0) then !valid handle
          ierr=0; tm=omp_get_wtime()
          if(tm.ge.timers(time_handle)%beg_time+timers(time_handle)%time_interval) then
           timers(time_handle)=timer_t(-1d0,-1d0); handle_sp=handle_sp-1; handle_stack(handle_sp)=time_handle
           time_handle=-1; time_is_off=.true.
          endif
         else
          ierr=2
         endif
        else
         ierr=1
        endif
        return
        end function time_is_off
!----------------------------------------
        real(8) function timer_tick_sec()
!This function returns the wall clock tick length in seconds.
        if(timer_tick.le.0d0) timer_tick=omp_get_wtick()
        timer_tick_sec=timer_tick
        return
        end function timer_tick_sec

       end module timers
