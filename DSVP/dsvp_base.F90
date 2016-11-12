!Domain-specific virtual processor (DSVP): Abstract base module.
!This module provides the infrastructure for the tensor algebra processor.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2016/11/12

!Copyright (C) 2014-2016 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2016 Oak Ridge National Laboratory (UT-Battelle)

!This file is part of ExaTensor.

!ExaTensor is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published
!by the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.

!ExaTensor is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!GNU Lesser General Public License for more details.

!You should have received a copy of the GNU Lesser General Public License
!along with ExaTensor. If not, see <http://www.gnu.org/licenses/>.

       module dsvp_base
        use dil_basic
        use timers
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !default output for this module
        integer(INTD), private:: DEBUG=0    !debugging mode
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
 !Error codes:
        integer(INTD), parameter, public:: DSVP_SUCCESS=SUCCESS       !success
        integer(INTD), parameter, public:: DSVP_ERROR=GENERIC_ERROR   !generic error code
        integer(INTD), parameter, public:: DSVP_ERR_INVALID_ARGS=1    !invalid procedure arguments
        integer(INTD), parameter, public:: DSVP_ERR_INVALID_REQ=2     !invalid request
        integer(INTD), parameter, public:: DSVP_ERR_MEM_ALLOC_FAIL=3  !failed memory allocation
 !DSVP status:
        integer(INTD), parameter, public:: DSVP_STAT_OFF=0    !DSVP is off (either not initialized or turned off)
        integer(INTD), parameter, public:: DSVP_STAT_ACTIVE=1 !DSVP has been initialized and is active now
        integer(INTD), parameter, public:: DSVP_STAT_ERROR=2  !DSVP encountered an error
 !DS operand:
        integer(INTD), parameter, public:: DS_OPRND_NO_COMM=0      !no pending communication on the domain-specific operand
        integer(INTD), parameter, public:: DS_OPRND_FETCHING=100   !domain-specific operand is being fetched
        integer(INTD), parameter, public:: DS_OPRND_UPLOADING=-100 !domain-specific operand is being uploaded
!DERIVED TYPES:
 !Domain-specific operand:
        type, abstract, public:: ds_oprnd_t
         logical:: active=.FALSE.                                   !TRUE is the operand is defined (active), FALSE otherwise (undefined)
         logical:: delivered=.FALSE.                                !TRUE when the domain-specific operand becomes locally available
         integer(INTD):: in_route=DS_OPRND_NO_COMM                  !communication status: {DS_OPRND_NO_COMM,DS_OPRND_FETCHING,DS_OPRND_UPLOADING}
         contains
          !public:
          procedure(ds_oprnd_self_i), deferred, public:: is_remote  !checks whether the domain-specific operand is local or remote
          procedure(ds_oprnd_self_i), deferred, public:: prefetch   !starts prefetching a remote domain-specific operand (acquires local resources!)
          procedure(ds_oprnd_self_i), deferred, public:: upload     !starts uploading the domain-specific operand to its remote location
          procedure(ds_oprnd_self_i), deferred, public:: sync       !synchronizes the currently pending communication on the domain-specific operand
          procedure(ds_oprnd_self_i), deferred, public:: release    !destroys the local copy of the domain-specific operand (releases local resources!), but the operand will still be defined
          procedure, public:: is_active=>DSOprndIsActive            !returns TRUE if the domain-specific operand is active (defined)
          procedure, public:: is_delivered=>DSOprndIsDelivered      !returns TRUE if the domain-specific operand is locally available
          !protected:
          procedure, public:: mark_active_=>DSOprndMarkActive           !marks the domain-specific operand active (defined)
          procedure, public:: mark_inactive_=>DSOprndMarkInactive       !marks the domain-specific operand inactive (undefined), local resources are released
          procedure, public:: mark_delivered_=>DSOprndMarkDelivered     !marks the domain-specific operand locally available
          procedure, public:: mark_undelivered_=>DSOprndMarkUndelivered !marks the domain-specific operand locally unavailable, local resources are released
          procedure, public:: set_comm_stat_=>DSOprndSetCommStat        !sets the communication status
        end type ds_oprnd_t
 !Domain-specific instruction:
        type, abstract, public:: ds_instr_t
         
        end type ds_instr_t
 !Domain-specific virtual processor:
        type, abstract, public:: dsvp_t
         integer(INTD), private:: stat=DSVP_STAT_OFF          !current DSVP status
         integer(INTL), private:: id=-1                       !DSVP unique ID
         integer(INTL), private:: instr_received=0_INTL       !total number of received domain-specific instructions
         integer(INTL), private:: instr_processed=0_INTL      !total number of processed (retired) instructions (both successful and failed)
         integer(INTL), private:: instr_failed=0_INTL         !total number of failed instructions
         real(8), private:: time_start                        !start time stamp (sec)
         character(:), allocatable, private:: description     !symbolic description of the DSVP
         contains
          !public:
          procedure(dsvp_init_i), deferred, public:: init                               !initializes DSVP to an active state
          procedure(dsvp_init_i), deferred, public:: shutdown                           !shutdowns DSVP
          procedure(dsvp_comm_instr_i), deferred, public:: fetch_instructions           !fetches a block of domain-specific instructions from another DSVP
          procedure(dsvp_comm_instr_i), deferred, public:: return_retired_instructions  !returns back a block of retired instructions with their statuses
          procedure(dsvp_comm_instr_i), deferred, public:: send_instructions            !sends a block of domain-specific instructions to another DSVP for execution
          procedure(dsvp_comm_instr_i), deferred, public:: receive_retired_instructions !receives back a block of retired instructions with their statuses
          procedure, public:: set_description=>DSVPSetDescription                !sets DSVP ID and symbolic description
          procedure, public:: get_description=>DSVPGetDescription                !gets DSVP ID and symbolic description
          procedure, public:: get_status=>DSVPGetStatus                          !returns the current status of the DSVP
          procedure, public:: time_active=>DSVPTimeActive                        !returns the time DSVP is active in seconds
          !protected:
          procedure, public:: start_time_=>DSVPStartTime                         !starts the time when DSVP is initialized
          procedure, public:: clean_=>DSVPClean                                  !cleans the DSVP state after the destruction
          procedure, public:: set_status_=>DSVPSetStatus                         !sets the DSVP status
          procedure, public:: incr_recv_instr_counter_=>DSVPIncrRecvInstrCounter !increments the receieved instruction counter
          procedure, public:: incr_rtrd_instr_counter_=>DSVPIncrRtrdInstrCounter !increments the processed (retired) instruction counter
          procedure, public:: incr_fail_instr_counter_=>DSVPIncrFailInstrCounter !increments the failed instruction counter
        end type dsvp_t
!INTERFACES:
 !Abstract:
        abstract interface
  !ds_oprnd_t:
   !self:
         subroutine ds_oprnd_self_i(this,ierr)
          import:: ds_oprnd_t,INTD
          implicit none
          class(ds_oprnd_t), intent(inout):: this
          integer(INTD), intent(out), optional:: ierr
         end subroutine ds_oprnd_self_i
  !ds_instr_t:

  !dsvp_t:
   !init and shutdown:
         subroutine dsvp_init_i(this,ierr)
          import:: dsvp_t,INTD
          implicit none
          class(dsvp_t), intent(inout):: this
          integer(INTD), intent(out), optional:: ierr
         end subroutine dsvp_init_i
   !instruction packet send/receive:
         subroutine dsvp_comm_instr_i(this,dsvp_id,ierr)
          import:: dsvp_t,INTD,INTL
          implicit none
          class(dsvp_t), intent(inout):: this
          integer(INTL), intent(in):: dsvp_id
          integer(INTD), intent(out), optional:: ierr
         end subroutine dsvp_comm_instr_i
        end interface
!DATA:

!VISIBILITY:
 !ds_oprnd_t:
        private DSOprndIsActive
        private DSOprndIsDelivered
        private DSOprndMarkActive
        private DSOprndMarkInactive
        private DSOprndMarkDelivered
        private DSOprndMarkUndelivered
        private DSOprndSetCommStat
 !ds_instr_t:
        
 !dsvp_t:
        private DSVPStartTime
        private DSVPClean
        private DSVPSetDescription
        private DSVPGetDescription
        private DSVPSetStatus
        private DSVPGetStatus
        private DSVPTimeActive
        private DSVPIncrRecvInstrCounter
        private DSVPIncrRtrdInstrCounter
        private DSVPIncrFailInstrCounter
!IMPLEMENTATION:
       contains
![ds_oprnd_t]==========================================
        function DSOprndIsActive(this,ierr) result(res)
!Returns TRUE if the domain-specific operand is active (defined).
         implicit none
         logical:: res                               !out: answer
         class(ds_oprnd_t), intent(in):: this        !in: domain-specific operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         res=this%active
         if(present(ierr)) ierr=errc
         return
        end function DSOprndIsActive
!---------------------------------------------------------
        function DSOprndIsDelivered(this,ierr) result(res)
!Returns TRUE if the domain-specific operand is locally available.
         implicit none
         logical:: res                               !out: answer
         class(ds_oprnd_t), intent(in):: this        !in: domain-specific operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS; res=.FALSE.
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) res=this%delivered
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end function DSOprndIsDelivered
!----------------------------------------------
        subroutine DSOprndMarkActive(this,ierr)
!Marks the domain-specific operand as active (defined).
         implicit none
         class(ds_oprnd_t), intent(inout):: this     !inout: domain-specific operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         this%active=.TRUE.
         if(present(ierr)) ierr=errc
         return
        end subroutine DSOprndMarkActive
!------------------------------------------------
        subroutine DSOprndMarkInactive(this,ierr)
!Marks the domain-specific operand as inactive (undefined).
!The local resources will be automatically released.
         implicit none
         class(ds_oprnd_t), intent(inout):: this     !inout: domain-specific operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,ier

         errc=DSVP_SUCCESS
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(this%in_route.eq.DS_OPRND_NO_COMM) then
            if(this%is_delivered(ier)) call this%mark_undelivered_(errc)
            if(errc.eq.DSVP_SUCCESS) errc=ier
            this%active=.FALSE.
           else
            errc=DSVP_ERR_INVALID_REQ
           endif
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSOprndMarkInactive
!------------------------------------------------
        subroutine DSOprndMarkDelivered(this,ierr)
!Marks the domain-specific operand as delivered (locally available).
         implicit none
         class(ds_oprnd_t), intent(inout):: this     !inout: domain-specific operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           call this%set_comm_stat_(DS_OPRND_NO_COMM,errc)
           this%delivered=.TRUE.
          endif
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSOprndMarkDelivered
!---------------------------------------------------
        subroutine DSOprndMarkUndelivered(this,ierr)
!Marks the domain-specific operand as undelivered (locally unavailable).
!The local resources will be automatically released.
         implicit none
         class(ds_oprnd_t), intent(inout):: this     !inout: domain-specific operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,ier

         errc=DSVP_SUCCESS
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(this%in_route.eq.DS_OPRND_NO_COMM) then
            if(this%is_delivered(ier)) call this%release(errc)
            if(errc.eq.DSVP_SUCCESS) errc=ier
            this%delivered=.FALSE.
           else
            errc=DSVP_ERR_INVALID_REQ
           endif
          endif
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSOprndMarkUndelivered
!---------------------------------------------------
        subroutine DSOprndSetCommStat(this,sts,ierr)
!Sets the communication status on the domain-specific operand.
         implicit none
         class(ds_oprnd_t), intent(inout):: this
         integer(INTD), intent(in):: sts
         integer(INTD), intent(out), optional:: ierr
         integer(INTD):: errc

         select case(sts)
         case(DS_OPRND_NO_COMM,DS_OPRND_FETCHING,DS_OPRND_UPLOADING)
          this%in_route=sts
         case default
          errc=DSVP_ERR_INVALID_ARGS
         end select
         if(present(ierr)) ierr=errc
         return
        end subroutine DSOprndSetCommStat
![dsvp_t]==================================
        subroutine DSVPStartTime(this,ierr)
!Starts the time after initializing the DSVP.
         implicit none
         class(dsvp_t), intent(inout):: this         !inout: DSVP
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         this%time_start=thread_wtime()
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPStartTime
!--------------------------------------
        subroutine DSVPClean(this,ierr)
!Cleans the DSVP state after desctruction.
         implicit none
         class(dsvp_t), intent(inout):: this         !inout: DSVP
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         if(this%get_status(errc).eq.DSVP_STAT_OFF) then
          this%id=-1
          this%instr_received=0
          this%instr_processed=0
          this%instr_failed=0
          if(allocated(this%description)) deallocate(this%description)
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPClean
!--------------------------------------------------------
        subroutine DSVPSetDescription(this,id,descr,ierr)
!Sets the DSVP ID and description.
         implicit none
         class(dsvp_t), intent(inout):: this         !inout: DSVP
         integer(INTL), intent(in):: id              !in: unique ID
         character(*), intent(in):: descr            !in: symbolic description
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,l
         integer:: ier

         errc=DSVP_SUCCESS
         if(this%get_status(errc).ne.DSVP_STAT_OFF) then
          if(allocated(this%description)) deallocate(this%description)
          l=len_trim(descr)
          if(l.gt.0) then
           allocate(character(len=l)::this%description,STAT=ier)
           if(ier.eq.0) then
            this%description(1:l)=descr(1:l)
           else
            if(errc.eq.DSVP_SUCCESS) errc=DSVP_ERR_MEM_ALLOC_FAIL
           endif
          endif
          this%id=id
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPSetDescription
!------------------------------------------------------------------
        subroutine DSVPGetDescription(this,id,descr,descr_len,ierr)
!Gets the DSVP ID and description.
         implicit none
         class(dsvp_t), intent(in):: this            !in: DSVP
         integer(INTL), intent(out):: id             !out: DSVP ID
         character(*), intent(inout):: descr         !out: symbolic description
         integer(INTD), intent(out):: descr_len      !out: length of the description string
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,l

         errc=DSVP_SUCCESS
         if(this%get_status(errc).ne.DSVP_STAT_OFF) then
          descr_len=0
          if(allocated(this%description)) then
           descr_len=len_trim(this%description)
           if(descr_len.gt.0) descr(1:descr_len)=this%description(1:descr_len)
          endif
          id=this%id
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPGetDescription
!----------------------------------------------
        subroutine DSVPSetStatus(this,sts,ierr)
!Sets the DSVP status.
         implicit none
         class(dsvp_t), intent(inout):: this         !inout: DSVP
         integer(INTD), intent(in):: sts             !in: status
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         this%stat=sts
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPSetStatus
!----------------------------------------------------
        function DSVPGetStatus(this,ierr) result(sts)
!Returns the current status of the DSVP.
         implicit none
         integer(INTD):: sts                         !out: DSVP current status
         class(dsvp_t), intent(in):: this            !in: DSVP
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS
         sts=this%stat
         if(present(ierr)) ierr=errc
         return
        end function DSVPGetStatus
!----------------------------------------------------
        function DSVPTimeActive(this,ierr) result(tm)
!Returns the time DSVP is active in seconds.
         implicit none
         real(8):: tm                                !out: time active in seconds
         class(dsvp_t), intent(in):: this            !in: DSVP
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=DSVP_SUCCESS; tm=-1d0
         if(this%get_status(errc).ne.DSVP_STAT_OFF) then
          if(errc.eq.DSVP_SUCCESS) tm=thread_wtime(this%time_start)
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end function DSVPTimeActive
!----------------------------------------------------------
        subroutine DSVPIncrRecvInstrCounter(this,ierr,incr)
!Increments the received instruction counter.
         implicit none
         class(dsvp_t), intent(inout):: this         !inout: DSVP
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTL), intent(in), optional:: incr  !in: specific increment
         integer:: errc

         errc=DSVP_SUCCESS
         if(this%get_status(errc).ne.DSVP_STAT_OFF) then
          if(present(incr)) then
           this%instr_received=this%instr_received+incr
          else
           this%instr_received=this%instr_received+1_INTL
          endif
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPIncrRecvInstrCounter
!----------------------------------------------------------
        subroutine DSVPIncrRtrdInstrCounter(this,ierr,incr)
!Increments the retired instruction counter.
         implicit none
         class(dsvp_t), intent(inout):: this         !inout: DSVP
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTL), intent(in), optional:: incr  !in: specific increment
         integer:: errc

         errc=DSVP_SUCCESS
         if(this%get_status(errc).ne.DSVP_STAT_OFF) then
          if(present(incr)) then
           this%instr_processed=this%instr_processed+incr
          else
           this%instr_processed=this%instr_processed+1_INTL
          endif
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPIncrRtrdInstrCounter
!----------------------------------------------------------
        subroutine DSVPIncrFailInstrCounter(this,ierr,incr)
!Increments the failed instruction counter.
         implicit none
         class(dsvp_t), intent(inout):: this         !inout: DSVP
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTL), intent(in), optional:: incr  !in: specific increment
         integer:: errc

         errc=DSVP_SUCCESS
         if(this%get_status(errc).ne.DSVP_STAT_OFF) then
          if(present(incr)) then
           this%instr_failed=this%instr_failed+incr
          else
           this%instr_failed=this%instr_failed+1_INTL
          endif
         else
          errc=DSVP_ERR_INVALID_REQ
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DSVPIncrFailInstrCounter

       end module dsvp_base
