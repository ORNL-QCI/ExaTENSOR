!Distributed data storage primitives.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2015/05/06 (started 2015/03/18)
!All rights reserved!
        module distributed
!       use, intrinsic:: ISO_C_BINDING
        use service_mpi !includes ISO_C_BINDING & MPI (must stay public)
        use:: tensor_algebra, only: NO_TYPE,R4,R8,C8 !tensor data kinds
        implicit none
        public !because of mpi.mod (or mpif.h)
!PARAMETERS:
 !Output:
        integer, private:: CONS_OUT=6     !default output for this module
        logical, private:: VERBOSE=.true. !verbosity for errors
        logical, private:: DEBUG=.true.   !debugging mode
 !Data transfers:
  !Messaging:
        integer(INT_MPI), parameter, public:: MAX_MPI_MSG_VOL=2**23 !max number of elements in a single MPI message (larger to be split)
  !Lock types:
        integer(INT_MPI), parameter, public:: NO_LOCK=0
        integer(INT_MPI), parameter, public:: SHARED_LOCK=1
        integer(INT_MPI), parameter, public:: EXCLUSIVE_LOCK=2
  !Asynchronous data requests:
        integer(INT_MPI), parameter, public:: MPI_ASYNC_NOT=0 !blocking data request (default)
        integer(INT_MPI), parameter, public:: MPI_ASYNC_NRM=1 !non-blocking data request without a request handle
        integer(INT_MPI), parameter, public:: MPI_ASYNC_REQ=2 !non-blocking data request with a request handle
  !Data request status:
        integer(INT_MPI), parameter, public:: MPI_STAT_ONESIDED_ERR=-1  !one-sided communication error occurred
        integer(INT_MPI), parameter, public:: MPI_STAT_NONE=0           !data request has not been processed by MPI yet
        integer(INT_MPI), parameter, public:: MPI_STAT_PROGRESS_NRM=1   !data request is in progress without a request handle
        integer(INT_MPI), parameter, public:: MPI_STAT_PROGRESS_REQ=2   !data request is in progress with a request handle
        integer(INT_MPI), parameter, public:: MPI_STAT_COMPLETED_ORIG=3 !data request has completed at the origin
        integer(INT_MPI), parameter, public:: MPI_STAT_COMPLETED=4      !data request has completed both at the origin and target
!TYPES:
 !MPI window:
        type, private:: WinMPI_t
         integer(INT_MPI), private:: Window           !MPI window handle
         integer(INT_MPI), private:: DataType=NO_TYPE !data type: {R4,R8,C8}
         integer(INT_MPI), private:: DispUnit=0       !displacement unit size in bytes
         integer(INT_MPI), private:: CommMPI          !MPI communicator
         logical, private:: Dynamic                   !.true. if the MPI window is dynamic, .false. otherwise
         contains
          procedure, private:: clean=>WinMPIClean     !clean MPI window info
        end type WinMPI_t
 !Local MPI window descriptor:
        type, public:: WinLoc_t
         integer(INT_ADDR), private:: WinSize=0       !current size (in bytes) of the local part of the MPI window
         type(WinMPI_t), private:: WinMPI             !MPI window info
         contains
          procedure, public:: create=>DataWinCreate   !create an MPI window
          procedure, public:: destroy=>DataWinDestroy !destroy an MPI window
          procedure, public:: attach=>DataWinAttach   !attach a data segment to the (dynamic) MPI window
          procedure, public:: detach=>DataWinDetach   !detach a data segment from the (dynamic) MPI window
        end type WinLoc_t
 !Global data location descriptor:
        type, public:: DataDescr_t
         integer(INT_MPI), private:: RankMPI=-1 !MPI rank on which the data resides
         type(WinMPI_t), private:: WinMPI       !MPI window the data is exposed with
         integer(INT_ADDR), private:: Offset    !offset in the MPI window (in displacement units)
         integer(INT_COUNT), private:: DataVol  !number of entries (data volume)
         integer(INT_MPI), private:: Locked     !lock type set on the data: {NO_LOCK, SHARED_LOCK, EXCLUSIVE_LOCK}
         integer(INT_MPI), private:: StatMPI    !status of the data request (see parameters above)
         integer(INT_MPI), private:: ReqHandle  !MPI request handle (for communications with a request handle)
         contains
          procedure, public:: clean=>DataDescrClean            !clean a data descriptor
          procedure, public:: init=>DataDescrInit              !set up a data descriptor
          procedure, public:: locked=>DataDescrLocked          !mark MPI lock set/unset for this data descriptor
          procedure, public:: lock_data=>DataDescrLockData     !lock data in an MPI wnidow
          procedure, public:: unlock_data=>DataDescrUnlockData !unlock data in an MPI window
          procedure, public:: flush_data=>DataDescrFlushData   !complete an outstanding passive target communication without unlocking
!         procedure, public:: put_data=>DataDescrPutData       !store data from a local buffer to the location specified by a data descriptor
          procedure, public:: get_data=>DataDescrGetData       !load data referred to by a data descriptor into a local buffer
          procedure, public:: acc_data=>DataDescrAccData       !accumulate data from a local buffer to the location specified by a data descriptor
        end type DataDescr_t
!GLOBAL DATA:
!       ...
!FUNCTION VISIBILITY:
 !WinMPI_t:
        private WinMPIClean
 !WinLoc_t:
        private DataWinCreate
        private DataWinDestroy
        private DataWinAttach
        private DataWinDetach
 !DataDescr_t:
        private DataDescrClean
        private DataDescrInit
        private DataDescrLocked
        private DataDescrLockData
        private DataDescrUnlockData
        private DataDescrFlushData
!       private DataDescrPutData
        private DataDescrGetData
        private DataDescrAccData

        contains
!METHODS:
!========================================
        subroutine WinMPIClean(this,ierr)
!Cleans an MPI window info.
        implicit none
        class(WinMPI_t), intent(out):: this              !out: empty MPI window info
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0
        this%DataType=NO_TYPE
        this%DispUnit=0
        if(present(ierr)) ierr=errc
        return
        end subroutine WinMPIClean
!======================================================================
        subroutine DataWinCreate(this,data_type,ierr,comm_mpi,win_size)
!This (collective) subroutine creates an MPI data window.
!The MPI window is dynamic, unless <win_size> is specified.
        implicit none
        class(DataWin_t), intent(inout):: this             !out: Data window object
        integer(C_INT), intent(in):: data_type             !in: Data type: {R4,R8,C8}
        integer(INT_MPI), intent(inout), optional:: ierr   !out: error code (0:success)
        integer(INT_MPI), intent(in), optional:: comm_mpi  !in: MPI communicator (defaults to MPI_COMM_WORLD)
        integer(INT_ADDR), intent(in), optional:: win_size !in: if present, defines the size of the static MPI window in bytes
        integer(INT_MPI):: errc

        errc=0
          if(mpi_window.ne.MPI_WIN_NULL) then
          this%PrRank=process_rank
          this%Window=mpi_window
          if(present(comm_mpi)) then; comm=comm_mpi; else; comm=MPI_COMM_WORLD; endif
          call MPI_COMM_GROUP(comm,my_group,errc)
          if(errc.eq.0) then
           call MPI_WIN_GET_ATTR(mpi_window,MPI_WIN_DISP_UNIT,attr,flag,errc)
           if(errc.eq.0.and.flag) then
            this%DispUnit=int(attr,INT_MPI)
            call MPI_WIN_GET_GROUP(mpi_window,win_group,errc)
            if(errc.eq.0) then
             if(present(check_thorough)) then; flag=check_thorough; else; flag=.false.; endif
             res=MPI_UNEQUAL
             if(flag) then
              call MPI_GROUP_COMPARE(my_group,win_group,res,errc)
             else
              if(my_group.eq.win_group) res=MPI_IDENT
             endif
        return
        end subroutine DataWinCreate
!===========================================
        subroutine DataDescrClean(this,ierr)
!Cleans a data descriptor.
        implicit none
        class(DataDescr_t), intent(out):: this           !out: empty data descriptor
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0
        this%RankMPI=-1
        call this%WinMPI%clean(errc); if(errc.ne.0) errc=1
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrClean
!-----------------------------------------------------------------------------------
        subroutine DataDescrInit(this,process_rank,win_mpi,win_offset,data_vol,ierr)
!Data descriptor constructor.
        implicit none
        class(DataDescr_t), intent(inout):: this         !out: filled data descriptor
        integer(INT_MPI), intent(in):: process_rank      !in: MPI process rank where the data resides
        type(WinMPI_t), intent(in):: win_mpi             !in: MPI window descriptor the data is exposed with
        integer(INT_ADDR), intent(in):: win_offset       !in: offset in the MPI window (in displacement units)
        integer(INT_COUNT), intent(in):: data_vol        !in: data volume (number of elements)
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)

        errc=0
        if(process_rank.ge.0) then
         if(data_vol.ge.0) then
          this%RankMPI=process_rank
          this%WinMPI=win_mpi
          this%Offset=win_offset
          this%DataVol=data_vol
          this%Locked=NO_LOCK !initial state is always NOT_LOCKED
          this%StatMPI=MPI_STAT_NONE
         else
          errc=1
         endif
        else
         errc=2
        endif
        if(errc.ne.0) call this%clean()
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrInit
!------------------------------------------------------
        subroutine DataDescrLocked(this,lock_type,ierr)
!Locally marks the distributed data as locked/unlocked when it is locked/unlocked explicitly.
        implicit none
        class(DataDescr_t), intent(inout):: this         !inout: distributed data descriptor
        integer(INT_MPI), intent(in):: lock_type         !in: current data lock status: {NO_LOCK, SHARED_LOCK, EXCLUSIVE_LOCK}
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0
        if(this%RankMPI.ge.0) then !non-empty data descriptor
         if(lock_type.eq.NO_LOCK.or.lock_type.eq.SHARED_LOCK.or.lock_type.eq.EXCLUSIVE_LOCK) then
          this%Locked=lock_type
         else
          errc=1
         endif
        else
         errc=2
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrLocked
!----------------------------------------------------------
        subroutine DataDescrLockData(this,lock_type,ierr) !`Must also update the local (win,rank) lock table
!Locks the data referred to by a data descriptor, unless
!the corresponding (window,rank) has been already locked.
        implicit none
        class(DataDescr_t), intent(inout):: this         !inout: distributed data descriptor
        integer(INT_MPI), intent(in):: lock_type         !in: lock type: {NO_LOCK, SHARED_LOCK, EXCLUSIVE_LOCK}
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc
        logical:: already_locked

        errc=0
        if(this%RankMPI.ge.0) then
         if(this%Locked.eq.NO_LOCK) then
          already_locked=.false. !`Check whether this (win,rank) is already locked
          if(lock_type.eq.SHARED_LOCK) then
           if(.not.already_locked) call MPI_WIN_LOCK(MPI_LOCK_SHARED,this%RankMPI,0,this%WinMPI%Window,errc)
           if(errc.eq.0) then; this%Locked=lock_type; else; this%StatMPI=MPI_STAT_ONESIDED_ERR; errc=1; endif
          elseif(lock_type.eq.EXCLUSIVE_LOCK) then
           if(.not.already_locked) call MPI_WIN_LOCK(MPI_LOCK_EXCLUSIVE,this%RankMPI,0,this%WinMPI%Window,errc)
           if(errc.eq.0) then; this%Locked=lock_type; else; this%StatMPI=MPI_STAT_ONESIDED_ERR; errc=2; endif
          else
           errc=3
          endif
         else
          errc=4
         endif
        else
         errc=5
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrLockData
!------------------------------------------------
        subroutine DataDescrUnlockData(this,ierr) !`Must also update the local (win,rank) lock table
!Unlocks the data referred to by a data descriptor.
        implicit none
        class(DataDescr_t), intent(inout):: this         !inout: distributed data descriptor
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc
        logical:: locked

        errc=0
        if(this%RankMPI.ge.0) then
         if(this%Locked.ne.NO_LOCK) then
          locked=.true. !`Check whether this (win,rank) is actually locked
          if(locked) call MPI_WIN_UNLOCK(this%RankMPI,this%WinMPI%Window,errc)
          if(errc.eq.0) then
           this%Locked=NO_LOCK
           this%StatMPI=MPI_STAT_COMPLETED
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR
           errc=1
          endif
         else
          errc=2
         endif
        else
         errc=3
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrUnlockData
!-----------------------------------------------------
        subroutine DataDescrFlushData(this,ierr,local)
!Flushes the MPI passive target communication referred to by a data descriptor.
        implicit none
        class(DataDescr_t), intent(inout):: this         !inout: distributed data descriptor
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        logical, intent(in), optional:: local            !in: if .true., the data is flushed only at the origin
        integer(INT_MPI):: errc
        logical:: lcl

        errc=0
        if(present(local)) then; lcl=local; else; lcl=.false.; endif !default is flushing both at the origin and target
        if(this%RankMPI.ge.0) then
         if(this%Locked.ne.NO_LOCK) then
          if(lcl) then
           call MPI_WIN_FLUSH_LOCAL(this%RankMPI,this%WinMPI%Window,errc) !complete at the origin only
           if(errc.eq.0) then; this%StatMPI=MPI_STAT_COMPLETED_ORIG; else; this%StatMPI=MPI_STAT_ONESIDED_ERR; errc=1; endif
          else
           call MPI_WIN_FLUSH(this%RankMPI,this%WinMPI%Window,errc) !complete both at the origin and target
           if(errc.eq.0) then; this%StatMPI=MPI_STAT_COMPLETED; else; this%StatMPI=MPI_STAT_ONESIDED_ERR; errc=2; endif
          endif
         else
          errc=3
         endif
        else
         errc=4
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrFlushData
!-----------------------------------------------------------
        subroutine DataDescrGetData(this,loc_ptr,ierr,async)
!Initiates a load of distributed data into a local buffer.
        implicit none
        class(DataDescr_t), intent(inout):: this         !inout: distributed data descriptor
        type(C_PTR), intent(in):: loc_ptr                !in: pointer to a local buffer
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI), intent(in), optional:: async   !in: if present, the passive target communication is not completed here
        integer(INT_MPI):: errc,asnc
        real(4), pointer:: r4_ptr(:)
        real(8), pointer:: r8_ptr(:)
        complex(8), pointer:: c8_ptr(:)

        errc=0
        if(present(async)) then; asnc=async; else; asnc=MPI_ASYNC_NOT; endif
        if(asnc.eq.MPI_ASYNC_NOT.or.asnc.eq.MPI_ASYNC_NRM.or.asnc.eq.MPI_ASYNC_REQ) then
         if(.not.c_associated(loc_ptr,C_NULL_PTR)) then
          if(this%RankMPI.ge.0) then
           if(this%StatMPI.eq.MPI_STAT_NONE.or.this%StatMPI.eq.MPI_STAT_COMPLETED) then
            if(this%Locked.eq.NO_LOCK) then
             call this%lock_data(SHARED_LOCK,errc); if(errc.ne.0) errc=1
            endif
            if(errc.eq.0) then
             if(this%DataVol.gt.0) then
              select case(this%WinMPI%DataType)
              case(R4)
               call c_f_pointer(loc_ptr,r4_ptr,(/this%DataVol/))
               
              case(R8)
               call c_f_pointer(loc_ptr,r8_ptr,(/this%DataVol/))
               call start_get_r8(r8_ptr,errc); if(errc.ne.0) errc=1
              case(C8)
               call c_f_pointer(loc_ptr,c8_ptr,(/this%DataVol/))
               
              case(NO_TYPE)
               errc=1
              case default
               errc=1
              endif
              if(errc.eq.0.and.asnc.eq.MPI_ASYNC_NOT)
               call this%unlock_data(errc); if(errc.ne.0) errc=1
              endif
             elseif(this%DataVol.lt.0) then
              errc=1
             endif
            else
             errc=1
            endif
           else
            errc=1
           endif
          else
           errc=1
          endif
         else
          errc=1
         endif
        else
         errc=1
        endif
        if(present(ierr)) ierr=errc
        return

        contains

         subroutine start_get_r8(r8_arr,jerr)
         real(8), intent(inout), asynchronous:: r8_arr(*)
         integer(INT_MPI), intent(out):: jerr
         integer(INT_COUNT):: ji,js
         integer(INT_MPI):: jv
         jerr=0; ji=int(MAX_MPI_MSG_VOL,INT_COUNT)
         if(asnc.ne.MPI_ASYNC_REQ) then !regular
          do js=1,this%DataVol,ji
           jv=int(min(this%DataVol-js+1,ji),INT_MPI)
           call MPI_GET(r8_arr,jv,MPI_REAL8,this%RankMPI,this%Offset,jv,MPI_REAL8,this%WinMPI%Window,jerr)
           if(jerr.ne.0) exit
          enddo
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_NRM
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=1
          endif
         else !request-handle
          call MPI_GET(r8_arr,jv,MPI_REAL8,this%RankMPI,this%Offset,jv,MPI_REAL8,this%WinMPI%Window,this%ReqHandle,jerr)
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_REQ
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=2
          endif
         endif
         return
         end subroutine start_get_r8

        end subroutine DataDescrGetData

        end module distributed
