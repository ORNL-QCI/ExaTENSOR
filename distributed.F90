!Distributed data storage infrastructure.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2015/05/12 (started 2015/03/18)
!CONCEPTS:
! * Each MPI process can participate in one or more distributed memory spaces,
!   where each distributed memory space is defined within a specific MPI communicator.
!   Within each distributed memory space, each MPI process opens a number
!   of dynamic MPI windows, which are opaque to the user.
! * Whenever a new (persistent) tensor block is allocated on the MPI process,
!   it is attached to one of the dynamic MPI windows of a specified distributed memory space
!   and the corresponding data descriptor will be sent to the manager for registration.
! * Upon a request from the manager, a tensor block can be detached from
!   an MPI window of the corresponding distributed memory space and destroyed.
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
 !MPI window info:
        type, private:: WinMPI_t
         integer(INT_MPI), private:: Window            !MPI window handle
         integer(INT_MPI), private:: DispUnit=0        !displacement unit size in bytes
         integer(INT_MPI), private:: CommMPI           !MPI communicator the window is associated with
         logical, private:: Dynamic                    !.true. if the MPI window is dynamic, .false. otherwise
         contains
          procedure, private:: clean=>WinMPIClean      !clean MPI window info
        end type WinMPI_t
 !Local MPI window descriptor:
        type, private:: DataWin_t
         integer(INT_ADDR), private:: WinSize=-1       !current size (in bytes) of the local part of the MPI window
         type(WinMPI_t), private:: WinMPI              !MPI window info
         contains
          procedure, private:: clean=>DataWinClean     !clean MPI data window info (for uninitialized windows only)
          procedure, private:: create=>DataWinCreate   !create an MPI window (collective)
          procedure, private:: destroy=>DataWinDestroy !destroy an MPI window (collective)
          procedure, private:: attach=>DataWinAttach   !attach a data segment to the (dynamic) MPI window
          procedure, private:: detach=>DataWinDetach   !detach a data segment from the (dynamic) MPI window
          procedure, private:: sync=>DataWinSync       !synchronize the private copy of the MPI window with its public copy
        end type DataWin_t
 !Distributed memory space descriptor:
        type, public:: DistrSpace_t
         integer(INT_MPI), private:: NumWins=0                !number of data windows in the distributed space
         integer(INT_MPI), private:: CommMPI                  !MPI communicator
         type(DataWin_t), allocatable, private:: DataWins(:)  !local data (MPI) windows
         character(DISTR_SPACE_NAME_LEN), private:: SpaceName !distributed space name
         contains
          procedure, public:: create=>DistrSpaceCreate        !create a distributed space (collective)
          procedure, public:: destroy=>DistrSpaceDestroy      !destroy a distributed space (collective)
          procedure, public:: local_size=>DistrSpaceLocalSize !get the local size (bytes) of the distributed memory space
          procedure, public:: attach=>DistrSpaceAttach        !attach a local buffer to the distributed memory space
          procedure, public:: detach=>DistrSpaceDetach        !detach a local buffer from the distributed memory space
        end type DistrSpace_t
 !Global data location descriptor:
        type, public:: DataDescr_t
         integer(INT_MPI), private:: RankMPI=-1 !MPI rank on which the data resides
         type(WinMPI_t), private:: WinMPI       !MPI window the data is exposed with
         integer(INT_ADDR), private:: Offset    !offset in the MPI window (in displacement units)
         integer(INT_COUNT), private:: DataVol  !data volume (number of typed elements)
         integer(INT_MPI), private:: DataType   !data type of each element: {R4,R8,C8,...}
         integer(INT_MPI), private:: Locked     !lock type set on the data: {NO_LOCK, SHARED_LOCK, EXCLUSIVE_LOCK}
         integer(INT_MPI), private:: StatMPI    !status of the data request (see parameters above)
         integer(INT_MPI), private:: ReqHandle  !MPI request handle (for communications with a request handle)
         type(C_PTR), private:: LocPtr          !local C pointer to the data buffer (internal use only)
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
          procedure, public:: pack=>DataDescrPack              !pack the DataDescr_t object into a plain-byte packet
          procedure, public:: unpack=>DataDescrUnpack          !unpack the DataDescr_t object from a plain-byte packet
        end type DataDescr_t
!GLOBAL DATA:
!       ...
!FUNCTION VISIBILITY:
 !Global:
        public data_type_size
 !WinMPI_t:
        private WinMPIClean
 !DataWin_t:
        private DataWinClean
        private DataWinCreate
        private DataWinDestroy
        private DataWinAttach
        private DataWinDetach
        private DataWinSync
 !DistrSpace_t:
        private DistrSpaceCreate
        private DistrSpaceDestroy
        private DistrSpaceLocalSize
        private DistrSpaceAttach
        private DistrSpaceDetach
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
        private DataDescrPack
        private DataDescrUnpack

        contains
!METHODS:
!===============================================================
        integer(INT_MPI) function data_type_size(data_type,ierr)
        implicit none
        integer(INT_MPI), intent(in):: data_type         !in: data type handle: {R4,R8,C8,...}
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        real(4), parameter:: r4_(1)=0.0
        real(8), parameter:: r8_(1)=0d0
        complex(8), parameter:: c8_(1)=(0d0,0d0)
        integer(INT_MPI):: errc

        errc=0; data_type_size=0
        select case(data_type)
        case(R4)
         data_type_size=storage_size(r4_,kind=INT_MPI)
        case(R8)
         data_type_size=storage_size(r8_,kind=INT_MPI)
        case(C8)
         data_type_size=storage_size(c8_,kind=INT_MPI)
        case default
         if(VERBOSE) write(CONS_OUT,'("#ERROR(distributed::data_type_size): Unknown data type: ",i11)') data_type
         errc=1
        end select
        if(present(ierr)) ierr=errc
        return
        end function data_type_size
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
!=========================================
        subroutine DataWinClean(this,ierr)
!Cleans an uninitialized MPI data window. Should be used
!after an allocation of a DataWin_t object, just in case.
        implicit none
        class(DataWin_t), intent(inout):: this           !inout: MPI data window
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0
        this%WinSize=-1
        call this%WinMPI%clean(errc)
        if(present(ierr)) ierr=errc
        return
        end subroutine DataWinClean
!---------------------------------------------------
        subroutine DataWinCreate(this,comm_mpi,ierr) !COLLECTIVE
!This (collective) subroutine creates a dynamic MPI data window.
        implicit none
        class(DataWin_t), intent(inout):: this           !inout: MPI data window
        integer(INT_MPI), intent(in):: comm_mpi          !in: MPI communicator the data window to be created over
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc
        integer(INT_ADDR):: attr
        logical:: flag

        call MPI_BARRIER(comm_mpi,errc) !test the validity of the MPI communicator
        if(errc.eq.0) then
         if(this%WinSize.lt.0) then !uninitialized window
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
!The MPI data window must be empty at this point.
        implicit none
        class(DataWin_t), intent(inout):: this           !inout: data window
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        call MPI_BARRIER(this%WinMPI%CommMPI,errc) !test the validity of the MPI communicator
        if(errc.eq.0) then
         if(this%WinSize.eq.0) then !MPI window must be empty
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
!Attaches a local (contiguous) data buffer to the MPI data window.
!Allowed data types: {R4,R8,C8,...}: Must be at least 4-byte long!
        implicit none
        class(DataWin_t), intent(inout):: this           !inout: data window
        type(C_PTR), intent(in):: loc_ptr                !in: C pointer to the local data buffer
        integer(INT_ADDR), intent(in):: loc_size         !in: size of the local data buffer in bytes: Must be multiple of 4
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc
        integer(INT_ADDR):: vol
        real(4), pointer, contiguous:: r4_ptr(:)

        errc=0
        if(.not.c_associated(loc_ptr,C_NULL_PTR)) then
         if(loc_size.gt.0.and.mod(loc_size,4).eq.0) then
          vol=loc_size/4 !volume in 4-byte words
          if(this%WinSize.ge.0) then !data window has been initialized
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
!Detaches a previously attached local (contiguous) data buffer from the MPI data window.
!The size of the data buffer in bytes must be a multiple of 4.
        implicit none
        class(DataWin_t), intent(inout):: this           !inout: data window
        type(C_PTR), intent(in):: loc_ptr                !in: C pointer to the local data buffer to be detached
        integer(INT_ADDR), intent(in):: loc_size         !in: size of the local data buffer in bytes: Must be multiple of 4
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc
        integer(INT_ADDR):: vol
        real(4), pointer, contiguous:: r4_ptr(:)

        errc=0
        if(.not.c_associated(loc_ptr,C_NULL_PTR)) then
         if(loc_size.gt.0.and.mod(loc_size,4).eq.0) then
          vol=loc_size/4 !volume in 4-byte words
          if(this%WinSize.ge.loc_size) then !data window is active
           call c_f_pointer(loc_ptr,r4_ptr,(/vol/))
           call detach_buffer(r4_ptr,errc)
           if(errc.eq.0) then
            this%WinSize=this%WinSize-loc_size
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

         subroutine detach_buffer(r4_arr,jerr)
         real(4), intent(in):: r4_arr(*)
         integer(INT_MPI), intent(out):: jerr
         jerr=0
         call MPI_WIN_DETACH(this%WinMPI%Window,r4_arr,jerr)
         return
         end subroutine detach_buffer

        end subroutine DataWinDetach
!----------------------------------------
        subroutine DataWinSync(this,ierr)
!Synchronizes the private and public copy of the MPI data window.
        implicit none
        class(DataWin_t), intent(in):: this              !in: data window
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0
        if(this%WinSize.ge.0) then !data window is initialized
         call MPI_WIN_SYNC(this%WinMPI%Window,errc); if(errc.ne.0) errc=1
        else
         errc=2
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine DataWinSync
!==========================================================================
        subroutine DistrSpaceCreate(this,comm_mpi,num_wins,space_name,ierr) !COLLECTIVE
!Creates a distributed space with <num_wins> dynamic windows.
!This is a collective call and every MPI process must receive the same <num_wins>!
        implicit none
        class(DistrSpace_t), intent(inout):: this        !inout: distributed memory space
        integer(INT_MPI), intent(in):: comm_mpi          !in: MPI communicator the space to be created over
        integer(INT_MPI), intent(in):: num_wins          !in: number of data windows in the distributed memory space
        character(*), intent(in):: space_name            !in: distributed memory space name
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: i,errc

        call MPI_BARRIER(comm_mpi,errc) !test the validity of the MPI communicator
        if(errc.eq.0) then
         if(num_wins.gt.0) then
          allocate(this%DataWins(1:num_wins),STAT=errc)
          if(errc.eq.0) then
           do i=1,num_wins
            call this%DataWins(i)%clean(errc); if(errc.ne.0) exit
           enddo
           if(errc.eq.0) then
            do i=1,num_wins
             call this%DataWins(i)%create(comm_mpi,errc); if(errc.ne.0) exit
            enddo
            if(errc.eq.0) then
             this%NumWins=num_wins
             this%CommMPI=comm_mpi
             this%SpaceName=space_name(1:min(len_trim(space_name),DISTR_SPACE_NAME_LEN)) !alphabetic name for convenience
            else
             do i=num_wins,1,-1
              call this%DataWins(i)%destroy()
             enddo
             deallocate(this%DataWins)
             errc=1
            endif
           else
            deallocate(this%DataWins)
            errc=2
           endif
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
        end subroutine DistrSpaceCreate
!----------------------------------------------
        subroutine DistrSpaceDestroy(this,ierr) !COLLECTIVE
!Destroys a distributed memory space. This is a collective call.
!The distributed memory space descriptor passed here must be valid!
        implicit none
        class(DistrSpace_t), intent(inout):: this        !inout: distributed memory space
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: i,j,errc

        call MPI_BARRIER(this%CommMPI,errc) !test the validity of the MPI communicator
        if(errc.eq.0) then
         if(this%NumWins.gt.0) then !initialized distributed memory space
          if(this%local_size(errc).eq.0) then !distributed memory space must be empty
           if(errc.eq.0) then
            do i=this%NumWins,1,-1
             call this%DataWins(i)%destroy(j); if(j.ne.0) errc=1
            enddo
            if(errc.eq.0) then
             deallocate(this%DataWins,STAT=errc)
             if(errc.eq.0) then
              this%NumWins=0
              this%CommMPI=MPI_COMM_NULL
              this%SpaceName=' '
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
         else
          errc=5
         endif
        else
         errc=6
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine DistrSpaceDestroy
!----------------------------------------------------------------
        integer(INT_ADDR) function DistrSpaceLocalSize(this,ierr)
!Returns the total size of all local MPI data windows (in bytes)
!belonging to the distributed space.
        implicit none
        class(DistrSpace_t), intent(in):: this           !in: distributed memory space
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
!-----------------------------------------------------------------------------------
        subroutine DistrSpaceAttach(this,loc_ptr,data_type,data_vol,data_descr,ierr)
!Attaches a local (contiguous) buffer to the initialized distributed memory space.
!On success, returns a valid data descriptor that can be used for remote/local accesses.
        implicit none
        class(DistrSpace_t), intent(inout):: this        !inout: distributed memory space
        type(C_PTR), intent(in):: loc_ptr                !in: C pointer to the local data buffer
        integer(INT_MPI), intent(in):: data_type         !in: pre-existing data type: {R4,R8,C8,...}
        integer(INT_COUNT), intent(in):: data_vol        !in: positive data volume (number of typed elements)
        class(DataDescr_t), intent(out):: data_descr     !out: filled data descriptor
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: i,m,my_rank,errc
        integer(INT_ADDR):: min_mem,loc_size

        errc=0
        if(this%NumWins.gt.0) then !initialized distributed memory space
         m=1; min_mem=this%DataWins(m)%WinSize
         do i=2,this%NumWins !select the least occupied MPI data window
          if(this%DataWins(i)%WinSize.lt.min_mem) then
           min_mem=this%DataWins(i)%WinSize; m=i
          endif
         enddo
         loc_size=data_vol*data_type_size(data_type,errc) !data buffer size in bytes
         if(errc.eq.0.and.loc_size.gt.0) then
          call this%DataWins(m)%attach(loc_ptr,loc_size,errc)
          if(errc.eq.0) then
           call MPI_COMM_RANK(this%CommMPI,my_rank,errc)
           if(errc.eq.0) then
            call data_descr%init(my_rank,this%DataWins(m)%WinMPI,loc_ptr,data_type,data_vol,errc)
            if(errc.ne.0) then
             call this%DataWins(m)%detach(loc_ptr,loc_size)
             errc=1
            endif
           else
            call this%DataWins(m)%detach(loc_ptr,loc_size)
            errc=2
           endif
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
        end subroutine DistrSpaceAttach
!--------------------------------------------------------
        subroutine DistrSpaceDetach(this,data_descr,ierr)
!Detaches a previously attached local (contiguous) data buffer from the distributed memory space.
!The data buffer is specified via a data descriptor.
        implicit none
        class(DistrSpace_t), intent(inout):: this        !inout: distributed memory space
        class(DataDescr_t), intent(inout):: data_descr   !inout: valid data descriptor
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: i,m,my_rank,errc
        integer(INT_ADDR):: loc_size

        errc=0
        if(this%NumWins.gt.0) then !initialized distributed memory space
         call MPI_COMM_RANK(this%CommMPI,my_rank,errc)
         if(errc.eq.0) then
          if(data_descr%RankMPI.eq.my_rank) then
           loc_size=data_descr%DataVol*data_type_size(data_descr%DataType,errc) !data buffer size in bytes
           if(errc.eq.0.and.loc_size.gt.0) then
            m=0
            do i=1,this%NumWins
             if(this%DataWins(i)%WinMPI%Window.eq.data_descr%WinMPI%Window) then; m=i; exit; endif
            enddo
            if(m.gt.0) then
             call this%DataWins(m)%detach(data_descr%LocPtr,loc_size,errc)
             if(errc.eq.0) then
              call data_descr%clean()
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
         else
          errc=5
         endif
        else
         errc=6
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine DistrSpaceDetach
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
!------------------------------------------------------------------------------------------
        subroutine DataDescrInit(this,process_rank,win_mpi,loc_ptr,data_type,data_vol,ierr)
!Data descriptor constructor: Data is typed and its volume must be positive.
        implicit none
        class(DataDescr_t), intent(inout):: this         !out: filled data descriptor
        integer(INT_MPI), intent(in):: process_rank      !in: MPI process rank where the data resides
        type(WinMPI_t), intent(in):: win_mpi             !in: MPI window descriptor the data is exposed with
        type(C_PTR), intent(in):: loc_ptr                !in: C pointer to the local data buffer
        integer(INT_MPI), intent(in):: data_type         !in: data type: {R4,R8,C8,...}
        integer(INT_COUNT), intent(in):: data_vol        !in: positive data volume (number of elements)
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc,elem_size

        errc=0
        if(process_rank.ge.0) then
         if(data_vol.gt.0) then
          elem_size=data_type_size(data_type,errc)
          if(errc.eq.0.and.elem_size.gt.0) then
           this%RankMPI=process_rank
           this%WinMPI=win_mpi
           this%DataVol=data_vol
           this%DataType=data_type
           this%Locked=NO_LOCK !initial state is always NOT LOCKED
           this%StatMPI=MPI_STAT_NONE
           this%LocPtr=loc_ptr
           call MPI_Get_Displacement(loc_ptr,this%Offset,errc)
          else
           errc=1
          endif
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
!Locally marks the distributed data as locked/unlocked when it is locked/unlocked outside.
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
