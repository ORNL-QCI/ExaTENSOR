!Distributed data storage infrastructure.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2015/05/07 (started 2015/03/18)
!All rights reserved!
!CONCEPTS:
! * Each MPI process can participate in one or more distributed spaces,
!   where each distributed space is defined within a specific MPI communicator.
! * Within each distributed space, each MPI process opens a number of dynamic MPI windows.
! * Whenever a new (persistent) tensor block is allocated on the MPI process,
!   it is attached to one of the dynamic MPI windows of a specified distributed space
!   and the corresponding data descriptor will be sent to the manager for registration.
! * Upon a request from the manager, a tensor block can be detached
!   from an MPI window and destroyed.
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
 !Distributed data storage:
        integer(INT_MPI), parameter, protected:: DISTR_SPACE_NAME_LEN=128 !max length of a distributed space name
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
         integer(INT_MPI), private:: DispUnit=0       !displacement unit size in bytes
         integer(INT_MPI), private:: CommMPI          !MPI communicator the window is associated with
         logical, private:: Dynamic                   !.true. if the MPI window is dynamic, .false. otherwise
         contains
          procedure, private:: clean=>WinMPIClean     !clean MPI window info
        end type WinMPI_t
 !Local MPI window descriptor:
        type, public:: DataWin_t
         integer(INT_ADDR), private:: WinSize=-1      !current size (in bytes) of the local part of the MPI window
         type(WinMPI_t), private:: WinMPI             !MPI window info
         contains
          procedure, public:: create=>DataWinCreate   !create an MPI window (collective)
          procedure, public:: destroy=>DataWinDestroy !destroy an MPI window (collective)
          procedure, public:: attach=>DataWinAttach   !attach a data segment to the (dynamic) MPI window
          procedure, public:: detach=>DataWinDetach   !detach a data segment from the (dynamic) MPI window
          procedure, public:: sync=>DataWinSync       !synchronize the private copy of the MPI window with its public copy
        end type DataWin_t
 !Distributed space descriptor:
        type, public:: DistrSpace_t
         integer(INT_MPI), private:: NumWins=0                !number of data windows in the distributed space
         integer(INT_MPI), private:: CommMPI                  !MPI communicator
         type(DataWin_t), allocatable, private:: DataWins(:)  !local data windows
         character(DISTR_SPACE_NAME_LEN), private:: SpaceName !distributed space name
         contains
          procedure, public:: create=>DistrSpaceCreate        !create a distributed space (collective)
          procedure, public:: destroy=>DistrSpaceDestroy      !destroy a distributed space (collecive)
          procedure, public:: local_size=>DistrSpaceLocalSize !get the local size of the distributed space
        end type DistrSpace_t
 !Global data location descriptor:
        type, public:: DataDescr_t
         integer(INT_MPI), private:: RankMPI=-1 !MPI rank on which the data resides
         type(WinMPI_t), private:: WinMPI       !MPI window the data is exposed with
         integer(INT_ADDR), private:: Offset    !offset in the MPI window (in displacement units)
         integer(INT_COUNT), private:: DataVol  !number of entries (data volume)
         integer(INT_MPI), private:: DataType   !data type of each entry: {R4,R8,C8}
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
          procedure, public:: test_data=>DataDescrTestData     !test whether the data has been transferred to/from the origin (request)
          procedure, public:: wait_data=>DataDescrWaitData     !wait until the data has been transferred to/from the origin (request)
!         procedure, public:: put_data=>DataDescrPutData       !store data from a local buffer to the location specified by a data descriptor
          procedure, public:: get_data=>DataDescrGetData       !load data referred to by a data descriptor into a local buffer
          procedure, public:: acc_data=>DataDescrAccData       !accumulate data from a local buffer to the location specified by a data descriptor
        end type DataDescr_t
!GLOBAL DATA:
!       ...
!FUNCTION VISIBILITY:
 !WinMPI_t:
        private WinMPIClean
 !DataWin_t:
        private DataWinCreate
        private DataWinDestroy
        private DataWinAttach
        private DataWinDetach
        private DataWinSync
 !DistrSpace_t:
        private DistrSpaceCreate
        private DistrSpaceDestroy
        private DistrSpaceLocalSize
 !DataDescr_t:
        private DataDescrClean
        private DataDescrInit
        private DataDescrLocked
        private DataDescrLockData
        private DataDescrUnlockData
        private DataDescrFlushData
        private DataDescrTestData
        private DataDescrWaitData
!       private DataDescrPutData !`Do I really need put?
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
        this%DispUnit=0
        if(present(ierr)) ierr=errc
        return
        end subroutine WinMPIClean
!===================================================
        subroutine DataWinCreate(this,comm_mpi,ierr) !COLLECTIVE
!This (collective) subroutine creates a dynamic MPI data window.
        implicit none
        class(DataWin_t), intent(inout):: this             !out: data window
        integer(INT_MPI), intent(in):: comm_mpi            !in: MPI communicator the data window is created over
        integer(INT_MPI), intent(inout), optional:: ierr   !out: error code (0:success)
        integer(INT_MPI):: errc
        integer(INT_ADDR):: attr
        logical:: flag

        call MPI_BARRIER(comm_mpi,errc)
        if(errc.eq.0) then
         if(this%WinSize.lt.0) then
          call MPI_WIN_CREATE_DYNAMIC(MPI_INFO_NULL,comm_mpi,this%WinMPI%Window,errc)
          if(errc.eq.0) then
           call MPI_WIN_GET_ATTR(this%WinMPI%Window,MPI_WIN_DISP_UNIT,attr,flag,errc)
           if(errc.eq.0.and.flag) then
            this%WinMPI%DispUnit=int(attr,INT_MPI)
            this%WinMPI%CommMPI=comm_mpi
            this%WinMPI%Dynamic=.true.
            this%WinSize=0
           else
            errc=1
           endif
          else
           errc=2
          endif
         else
          errc=3
         endif
        else
         errc=4
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine DataWinCreate
!-------------------------------------------
        subroutine DataWinDestroy(this,ierr) !COLLECTIVE
!This (collective) subroutine destroys a dynamic MPI data window.
        implicit none
        class(DataWin_t), intent(inout):: this             !inout: data window
        integer(INT_MPI), intent(inout), optional:: ierr   !out: error code (0:success)
        integer(INT_MPI):: errc

        call MPI_BARRIER(this%WinMPI%CommMPI,errc)
        if(errc.eq.0) then
         if(this%WinSize.eq.0) then
          call MPI_WIN_FREE(this%WinMPI%Window,errc)
          if(errc.eq.0) then
           call this%WinMPI%clean(errc)
           this%WinSize=-1
          else
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
        end subroutine DataWinDestroy
!-----------------------------------------------------------
        subroutine DataWinAttach(this,loc_ptr,loc_size,ierr)
!Attaches a chunk of memory to the MPI data window.
        implicit none
        class(DataWin_t), intent(inout):: this           !inout: data window
        type(C_PTR), intent(in):: loc_ptr                !in: C pointer to the local buffer
        integer(INT_ADDR), intent(in):: loc_size         !in: size of the local buffer in bytes
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc
        integer(INT_ADDR):: vol
        real(4), pointer, contiguous:: r4_ptr(:)

        errc=0
        if(.not.c_associated(loc_ptr,C_NULL_PTR)) then
         if(loc_size.gt.0.and.mod(loc_size,4).eq.0) then
          vol=loc_size/4 !volume in 4-byte words
          if(this%WinSize.ge.0) then !data window is active
           call c_f_pointer(loc_ptr,r4_ptr,(/vol/))
           call attach_buffer(r4_ptr,loc_size,errc)
           if(errc.eq.0) then
            this%WinSize=this%WinSize+loc_size
           else
            errc=1
           endif
          else
           errc=2
          endif
         else
          errc=3
         endif
        else
         errc=4
        endif
        if(present(ierr)) ierr=errc
        return

        contains

         subroutine attach_buffer(r4_arr,bsize,jerr)
         real(4), intent(in):: r4_arr(*)
         integer(INT_ADDR), intent(in):: bsize
         integer(INT_MPI), intent(out):: jerr
         jerr=0
         call MPI_WIN_ATTACH(this%WinMPI%Window,r4_arr,bsize,jerr)
         return
         end subroutine attach_buffer

        end subroutine DataWinAttach
!-----------------------------------------------------------
        subroutine DataWinDetach(this,loc_ptr,loc_size,ierr)
!Detaches a chunk of memory from the MPI data window.
        implicit none
        class(DataWin_t), intent(inout):: this           !inout: data window
        type(C_PTR), intent(in):: loc_ptr                !in: C pointer to the local buffer
        integer(INT_ADDR), intent(in):: loc_size         !in: size of the local buffer in bytes
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc
        integer(INT_ADDR):: vol
        real(4), pointer, contiguous:: r4_ptr(:)

        errc=0
        if(.not.c_associated(loc_ptr,C_NULL_PTR)) then
         if(loc_size.gt.0.and.mod(loc_size,4).eq.0) then
          vol=loc_size/4 !volume in 4-byte words
          if(this%WinSize.ge.0) then !data window is active
           call c_f_pointer(loc_ptr,r4_ptr,(/vol/))
           call detach_buffer(r4_ptr,errc)
           if(errc.eq.0) then
            this%WinSize=this%WinSize+loc_size
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

         subroutine detach_buffer(r4_arr,jerr)
         real(4), intent(in):: r4_arr(*)
         integer(INT_MPI), intent(out):: jerr
         jerr=0
         call MPI_WIN_DETACH(this%WinMPI%Window,r4_arr,jerr)
         return
         end subroutine detach_buffer

        end subroutine DataWinDetach
!------------------------------------------------------
        subroutine DataWinSync(this,ierr)
!Synchronizes the private and public copy of the MPI data window.
        implicit none
        class(DataWin_t), intent(in):: this              !in: data window
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0
        if(this%WinSize.ge.0) then !data window is active
         call MPI_WIN_SYNC(this%WinMPI%Window,errc); if(errc.ne.0) errc=1
        else
         errc=2
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine DataWinSync
!==========================================================================
        subroutine DistrSpaceCreate(this,comm_mpi,num_wins,space_name,ierr) !COLLECTIVE
!Creates a distributed space with <num_wins> windows.
!This is a collective call. Every MPI process must receive the same <num_wins>!
        implicit none
        class(DistrSpace_t), intent(inout):: this              !out: distributed space
        integer(INT_MPI), intent(in):: comm_mpi                !in: MPI communicator the space is created over
        integer(INT_MPI), intent(in):: num_wins                !in: number of data windows in the distributed space
        character(*), intent(in):: space_name                  !in: distributed space name
        integer(INT_MPI), intent(inout), optional:: ierr       !out: error code (0:success)
        integer(INT_MPI):: i,errc

        call MPI_BARRIER(comm_mpi,errc)
        if(errc.eq.0) then
         if(num_wins.gt.0) then
          allocate(this%DataWins(1:num_wins),STAT=errc)
          if(errc.eq.0) then
           do i=1,num_wins
            call this%DataWins(i)%create(...,errc)
            if(errc.ne.0) exit
           enddo
           if(errc.eq.0) then
            this%NumWins=num_wins
            this%CommMPI=comm_mpi
            this%SpaceName=space_name(1:min(len_trim(space_name),DISTR_SPACE_NAME_LEN))
           else
            do i=num_wins,1,-1
             call this%DataWins(i)%destroy(errc)
            enddo
            errc=1
           endif
          else
           errc=2
          endif
         else
          errc=3
         endif
        else
         errc=4
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine DistrSpaceCreate
!------------------------------------------
        subroutine DistrSpaceDestroy() !COLLECTIVE
!Destroys a distributed space. This is a collective call.
        implicit none
        class(DistrSpace_t), intent(inout):: this        !inout: distributed space
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: i,errc

        errc=0
        
        if(present(ierr)) ierr=errc
        return
        end subroutine DistrSpaceDestroy
!----------------------------------------------------------------
        integer(INT_ADDR) function DistrSpaceLocalSize(this,ierr)
!Returns the total occupied size (bytes) of all local data windows
!belonging to the distributed space.
        implicit none
        class(DistrSpace_t), intent(in):: this           !in: distributed space
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: i,errc
        integer(INT_ADDR):: loc_size

        errc=0; DistrSpaceLocalSize=0
        if(this%NumWins.gt.0) then
         loc_size=0
         do i=1,this%NumWins
          if(this%DataWins(i)%WinSize.ge.0) then
           loc_size=loc_size+this%DataWins(i)%WinSize
          else
           errc=1
           exit
          endif
         enddo
         if(errc.eq.0) DistrSpaceLocalSize=loc_size
        else
         errc=2
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine DistrSpaceLocalSize
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
        real(4), pointer, contiguous:: r4_ptr(:)
        real(8), pointer, contiguous:: r8_ptr(:)
        complex(8), pointer, contiguous:: c8_ptr(:)

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
              select case(this%DataType)
              case(R4)
               call c_f_pointer(loc_ptr,r4_ptr,(/this%DataVol/))
               call start_get_r4(r4_ptr,errc); if(errc.ne.0) errc=2
              case(R8)
               call c_f_pointer(loc_ptr,r8_ptr,(/this%DataVol/))
               call start_get_r8(r8_ptr,errc); if(errc.ne.0) errc=3
              case(C8)
               call c_f_pointer(loc_ptr,c8_ptr,(/this%DataVol/))
               call start_get_c8(c8_ptr,errc); if(errc.ne.0) errc=4
              case(NO_TYPE)
               errc=5
              case default
               errc=6
              endif
              if(errc.eq.0.and.asnc.eq.MPI_ASYNC_NOT)
               call this%unlock_data(errc); if(errc.ne.0) errc=7
              endif
             elseif(this%DataVol.lt.0) then
              errc=8
             endif
            else
             errc=9
            endif
           else
            errc=10
           endif
          else
           errc=11
          endif
         else
          errc=12
         endif
        else
         errc=13
        endif
        if(present(ierr)) ierr=errc
        return

        contains

         subroutine start_get_r4(r4_arr,jerr)
         real(4), intent(inout), asynchronous:: r4_arr(*)
         integer(INT_MPI), intent(out):: jerr
         integer(INT_COUNT):: ji,js
         integer(INT_MPI):: jv
         jerr=0; ji=int(MAX_MPI_MSG_VOL,INT_COUNT)
         if(asnc.ne.MPI_ASYNC_REQ) then !regular
          do js=1,this%DataVol,ji
           jv=int(min(this%DataVol-js+1,ji),INT_MPI)
           call MPI_GET(r4_arr,jv,MPI_REAL4,this%RankMPI,this%Offset,jv,MPI_REAL4,this%WinMPI%Window,jerr)
           if(jerr.ne.0) exit
          enddo
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_NRM
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=1
          endif
         else !request-handle
          call MPI_GET(r4_arr,jv,MPI_REAL4,this%RankMPI,this%Offset,jv,MPI_REAL4,this%WinMPI%Window,this%ReqHandle,jerr)
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_REQ
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=2
          endif
         endif
         return
         end subroutine start_get_r4

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

         subroutine start_get_c8(c8_arr,jerr)
         complex(8), intent(inout), asynchronous:: c8_arr(*)
         integer(INT_MPI), intent(out):: jerr
         integer(INT_COUNT):: ji,js
         integer(INT_MPI):: jv
         jerr=0; ji=int(MAX_MPI_MSG_VOL,INT_COUNT)
         if(asnc.ne.MPI_ASYNC_REQ) then !regular
          do js=1,this%DataVol,ji
           jv=int(min(this%DataVol-js+1,ji),INT_MPI)
           call MPI_GET(c8_arr,jv,MPI_COMPLEX8,this%RankMPI,this%Offset,jv,MPI_COMPLEX8,this%WinMPI%Window,jerr)
           if(jerr.ne.0) exit
          enddo
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_NRM
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=1
          endif
         else !request-handle
          call MPI_GET(c8_arr,jv,MPI_COMPLEX8,this%RankMPI,this%Offset,jv,MPI_COMPLEX8,this%WinMPI%Window,this%ReqHandle,jerr)
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_REQ
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=2
          endif
         endif
         return
         end subroutine start_get_c8

        end subroutine DataDescrGetData

        end module distributed
