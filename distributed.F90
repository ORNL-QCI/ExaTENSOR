!Distributed data storage service (DDSS).
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2015/11/16 (started 2015/03/18)
!Copyright (C) 2015 Dmitry I. Lyakh (email: quant4me@gmail.com)
!Copyright (C) 2015 Oak Ridge National Laboratory (UT-Battelle)
!LICENSE: GPLv2

!CONCEPTS:
! * Each MPI process can participate in one or more distributed memory spaces (DMS),
!   where each distributed memory space is defined within a specific MPI communicator.
!   Internally, within each distributed memory space, each MPI process opens a number
!   of dynamic MPI windows for data storage, which are opaque to the user.
! * Whenever new persistent data (e.g., a tensor block) is allocated by the MPI process,
!   it can be attached to one of the dynamic MPI windows in a specific distributed memory space
!   and the corresponding data descriptor (DD) should be sent to the manager for registration.
!   The registered data descriptor can be used by other MPI processes to remotely access the data.
! * Upon a request from the manager, data (e.g., a tensor block) can be detached from
!   the corresponding distributed memory space and subsequently destroyed (if needed).
! * Data communication is accomplished via data transfer requests (DTR) and
!   data transfer completion requests (DTCR), using data descriptors. On each MPI process,
!   all data transfer requests are enumerated sequentially in the order they were issued (starting at 1).
! * All procedures return error codes where special return statuses must be distinguished,
!   for example TRY_LATER. Normally, error codes are smaller by absolute value integers.
!   Contrary, special return statuses should be closer to the HUGE by their absolute values.
! * Currently, chunks of memory that can be attached to a distributed memory space must be
!   multiple of 4 bytes. This is because the memory chunks are mapped to 32-bit words internally.
        module distributed
!       use, intrinsic:: ISO_C_BINDING
        use service_mpi !includes ISO_C_BINDING & MPI
        use:: tensor_algebra, only: TRY_LATER,NO_TYPE,R4,R8,C8,R4_,R8_,C8_ !some basic types and statuses
        implicit none
        private
!EXPOSE some <service_mpi>:
        public INT_MPI,INT_ADDR,INT_OFFSET,INT_COUNT !MPI integer kinds
        public jo                  !process log output device
        public impis               !size of the global MPI communicator
        public impir               !global MPI rank of the process
        public max_threads         !max number of threads per process
        public gpu_start,gpu_count !GPU range assigned to the process
!PARAMETERS:
 !Output:
        integer, private:: CONS_OUT=6     !default output for this module
        logical, private:: VERBOSE=.true. !verbosity for errors
        logical, private:: DEBUG=.true.   !debugging mode
 !Distributed memory spaces:
        integer(INT_MPI), parameter, public:: DISTR_SPACE_NAME_LEN=128 !max length of a distributed space name (multiple of 8)
 !Data transfers:
  !Active (rank,window) management:
        integer(INT_MPI), parameter, private:: HASH_MOD=256  !number of hash bins for (rank,window) searches
        integer(INT_MPI), parameter, private:: READ_SIGN=+1  !incoming traffic sign (reading direction)
        integer(INT_MPI), parameter, private:: WRITE_SIGN=-1 !outgoing traffic sign (writing direction)
  !Messaging:
        integer(INT_COUNT), parameter, private:: MAX_MPI_MSG_VOL=2**23 !max number of elements in a single MPI message (larger to be split)
        integer(INT_MPI), parameter, private:: MAX_ONESIDED_REQS=4096  !max number of outstanding one-sided data transfer requests per process
  !Lock types:
        integer(INT_MPI), parameter, public:: NO_LOCK=0        !no MPI lock (must be zero)
        integer(INT_MPI), parameter, public:: SHARED_LOCK=1    !shared MPI lock (must be positive)
        integer(INT_MPI), parameter, public:: EXCLUSIVE_LOCK=2 !exclusive MPI lock (must be positive)
  !Asynchronous data transfer requests:
        integer(INT_MPI), parameter, public:: MPI_ASYNC_NOT=0  !blocking data transfer request (default)
        integer(INT_MPI), parameter, public:: MPI_ASYNC_NRM=1  !non-blocking data transfer request without a request handle
        integer(INT_MPI), parameter, public:: MPI_ASYNC_REQ=2  !non-blocking data transfer request with a request handle
  !Data request status:
        integer(INT_MPI), parameter, public:: MPI_STAT_ONESIDED_ERR=-1  !one-sided communication error occurred
        integer(INT_MPI), parameter, public:: MPI_STAT_NONE=0           !data request has not been processed by MPI yet
        integer(INT_MPI), parameter, public:: MPI_STAT_PROGRESS_NRM=1   !data request is in progress without a request handle
        integer(INT_MPI), parameter, public:: MPI_STAT_PROGRESS_REQ=2   !data request is in progress with a request handle
        integer(INT_MPI), parameter, public:: MPI_STAT_COMPLETED_ORIG=3 !data request has completed at the origin
        integer(INT_MPI), parameter, public:: MPI_STAT_COMPLETED=4      !data request has completed both at the origin and target
!TYPES:
 !One-sided data transfer bookkeeping (internal use only):
  !Rank/window descriptor:
        type, private:: RankWin_t
         integer(INT_MPI), private:: Rank     !MPI rank: >=0
         integer(INT_MPI), private:: Window   !MPI window handle
         integer(INT_MPI), private:: LockType !current lock type: {NO_LOCK,SHARED_LOCK,EXCLUSIVE_LOCK} * sign{READ_SIGN,WRITE_SIGN}
         integer(INT_MPI), private:: RefCount !data transfer reference count (number of bound data descriptors): >=0
         integer(8), private:: LastSync       !data transfer ID for which the last unlock/flush was performed (>0|0:never_sync)
         contains
          procedure, private:: init=>RankWinInit !initialize/clean a (rank,window) descriptor
        end type RankWin_t
  !Rank/window linked list (active one-sided comms at origin):
        type, private:: RankWinList_t
         integer(8), private:: TransCount=0                         !total number of posted data transfer requests
         real(8), private:: TransSize=0d0                           !total size of all posted data transfers in bytes
         integer(INT_MPI), private:: NumEntries=0                   !number of active entries in the list
         integer(INT_MPI), private:: FirstFree=1                    !first free (inactive) entry
         type(RankWin_t), private:: RankWins(1:MAX_ONESIDED_REQS)   !(rank,window) entries
         integer(INT_MPI), private:: NextEntry(1:MAX_ONESIDED_REQS) !next entry within a bin (linked list)
         integer(INT_MPI), private:: PrevEntry(1:MAX_ONESIDED_REQS) !previous entry within a bin (linked list)
         integer(INT_MPI), private:: HashBin(0:HASH_MOD-1)=0        !first element in each hash bin
         contains
          procedure, private:: clean=>RankWinListClean !clean the (rank,window) list (initialization)
          procedure, private:: test=>RankWinListTest   !test whether a given (rank,window) entry is in the list (with an optional append)
          procedure, private:: delete=>RankWinListDel  !delete a given active (rank,window) entry
          procedure, private:: new_transfer=>RankWinNewTrans !register a new data transfer (increment the global transfer ID)
        end type RankWinList_t
 !Basic MPI window info:
        type, private:: WinMPI_t
         integer(INT_MPI), private:: Window            !MPI window handle
         integer(INT_MPI), private:: DispUnit=0        !MPI window displacement unit size in bytes
         integer(INT_MPI), private:: CommMPI           !MPI communicator the MPI window is associated with
         logical, private:: Dynamic                    !.true. if the MPI window is dynamic, .false. otherwise
         contains
          procedure, private:: clean=>WinMPIClean      !clean MPI window info
        end type WinMPI_t
 !Local MPI data window descriptor:
        type, private:: DataWin_t
         integer(INT_ADDR), private:: WinSize=-1       !current size (in bytes) of the local part of the MPI window
         type(WinMPI_t), private:: WinMPI              !MPI window info
         contains
          procedure, private:: clean=>DataWinClean     !clean the MPI data window descriptor
          procedure, private:: create=>DataWinCreate   !create an MPI window (collective)
          procedure, private:: destroy=>DataWinDestroy !destroy an MPI window (collective)
          procedure, private:: attach=>DataWinAttach   !attach a data segment to the (dynamic) MPI window
          procedure, private:: detach=>DataWinDetach   !detach a data segment from the (dynamic) MPI window
          procedure, private:: sync=>DataWinSync       !synchronize the private copy of the MPI window with its public copy
        end type DataWin_t
 !Distributed memory space descriptor:
        type, public:: DistrSpace_t
         integer(INT_MPI), private:: NumWins=0                !number of MPI windows in the distributed space
         integer(INT_MPI), private:: CommMPI                  !MPI communicator the distributed memory space is created over
         type(DataWin_t), allocatable, private:: DataWins(:)  !local MPI data windows
         character(DISTR_SPACE_NAME_LEN), private:: SpaceName !distributed memory space name
         contains
          procedure, public:: create=>DistrSpaceCreate        !create a distributed memory space (collective)
          procedure, public:: destroy=>DistrSpaceDestroy      !destroy a distributed memory space (collective)
          procedure, public:: local_size=>DistrSpaceLocalSize !get the local size (bytes) of the distributed memory space
          procedure, public:: attach=>DistrSpaceAttach        !attach a local buffer to the distributed memory space
          procedure, public:: detach=>DistrSpaceDetach        !detach a local buffer from the distributed memory space
        end type DistrSpace_t
 !Global data location descriptor:
        type, public:: DataDescr_t
         type(C_PTR), private:: LocPtr          !local C pointer to the original data buffer (internal use only)
         integer(INT_MPI), private:: RankMPI=-1 !MPI rank on which the data resides
         type(WinMPI_t), private:: WinMPI       !info on the MPI window the data is exposed with
         integer(INT_ADDR), private:: Offset    !offset in the MPI window (in displacement units)
         integer(INT_COUNT), private:: DataVol  !data volume (number of typed elements)
         integer(INT_MPI), private:: DataType   !data type of each element: {R4,R8,C8,...}, see <tensor_algebra.F90>
         integer(8), private:: TransID          !data transfer request ID (sequential data transfer number within this process)
         integer(INT_MPI), private:: StatMPI    !status of the data transfer request (see MPI_STAT_XXX parameters above)
         integer(INT_MPI), private:: ReqHandle  !MPI request handle (for MPI communications with a request handle)
         contains
          procedure, private:: clean=>DataDescrClean            !clean a data descriptor
          procedure, private:: init=>DataDescrInit              !set up a data descriptor (initialization)
          procedure, public:: flush_data=>DataDescrFlushData    !complete an outstanding data transfer request
          procedure, public:: test_data=>DataDescrTestData      !test whether the data has been transferred to/from the origin (request)
          procedure, public:: wait_data=>DataDescrWaitData      !wait until the data has been transferred to/from the origin (request)
          procedure, public:: get_data=>DataDescrGetData        !load data referred to by a data descriptor into a local buffer
          procedure, public:: acc_data=>DataDescrAccData        !accumulate data from a local buffer to the location specified by a data descriptor
!          procedure, public:: pack=>DataDescrPack               !pack the DataDescr_t object into a plain-byte packet
!          procedure, public:: unpack=>DataDescrUnpack           !unpack the DataDescr_t object from a plain-byte packet
        end type DataDescr_t
!GLOBAL DATA:
 !MPI one-sided data transfer bookkeeping (master thread only):
        type(RankWinList_t), target, private:: RankWinRefs !container for active one-sided communications initiated at the local origin
!FUNCTION VISIBILITY:
 !Global:
        public data_type_size
 !RankWin_t:
        private RankWinInit
 !RankWinList_t:
        private RankWinListClean
        private RankWinListTest
        private RankWinListDel
        private RankWinNewTrans
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
        private DataDescrFlushData
        private DataDescrTestData
        private DataDescrWaitData
        private DataDescrGetData
        private DataDescrAccData
!        private DataDescrPack    !`Write
!        private DataDescrUnpack  !`Write

        contains
!METHODS:
!===============================================================
        integer(INT_MPI) function data_type_size(data_type,ierr) !Fortran2008: storage_size()
!Returns the local storage size of a given (pre-registered) data type in bytes.
!All registered data types must possess a size that is multiple of 4 (in bytes)!
        implicit none
        integer(INT_MPI), intent(in):: data_type         !in: pre-registered data type handle: {NO_TYPE,R4,R8,C8,...}
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0; data_type_size=0
        select case(data_type)
        case(R4)
         data_type_size=storage_size(R4_,kind=INT_MPI) !bits
        case(R8)
         data_type_size=storage_size(R8_,kind=INT_MPI) !bits
        case(C8)
         data_type_size=storage_size(C8_,kind=INT_MPI) !bits
        case(NO_TYPE)
         data_type_size=0
        case default
         if(VERBOSE) write(CONS_OUT,'("#ERROR(distributed::data_type_size): Unknown data type: ",i11)') data_type
         errc=1; data_type_size=-1
        end select
        if(errc.eq.0) then
         if(mod(data_type_size,8).eq.0) then
          data_type_size=data_type_size/8 !convert bits to bytes
          if(mod(data_type_size,4).ne.0) errc=2 !data type size must be a multiple of 4
         else
          if(VERBOSE) write(CONS_OUT,'("#ERROR(distributed::data_type_size): Fractional type size detected: ",i11)') data_type_size
          errc=3 !fractional data type size is not allowed
         endif
        endif
        if(present(ierr)) ierr=errc
        return
        end function data_type_size
!=================================================
        subroutine RankWinInit(this,rank,win,ierr)
!Initializes a (rank,window) descriptor (both <rank> and <win> must be present) OR
!cleans a (rank,window) descriptor (neither <rank> nor <win> should appear here).
        implicit none
        class(RankWin_t), intent(inout):: this           !inout: (rank,window) descriptor
        integer(INT_MPI), intent(in), optional:: rank    !in: MPI rank (>=0)
        integer(INT_MPI), intent(in), optional:: win     !in: MPI window handle
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0
        if(present(rank)) then
         if(rank.ge.0) then
          if(present(win)) then !active descriptor
           this%Rank=rank
           this%Window=win
           this%LockType=NO_LOCK
           this%RefCount=0
           this%LastSync=0 !0 means "never synced"
          else
           errc=1
          endif
         else
          errc=2
         endif
        else
         if(present(win)) then
          errc=3
         else !empty descriptor
          this%Rank=-1
          this%Window=0
          this%LockType=NO_LOCK
          this%RefCount=0
          this%LastSync=0
         endif
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine RankWinInit
!=============================================
        subroutine RankWinListClean(this,ierr)
!Cleans the (rank,window) list (must be called before use).
        implicit none
        class(RankWinList_t), intent(inout):: this       !inout: (rank,window) list
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: i

        this%TransCount=0; this%TransSize=0d0
        this%NumEntries=0; this%FirstFree=1; this%HashBin(:)=0
        this%PrevEntry(1)=0; do i=2,MAX_ONESIDED_REQS; this%PrevEntry(i)=i-1; enddo
        do i=1,MAX_ONESIDED_REQS-1; this%NextEntry(i)=i+1; enddo; this%NextEntry(MAX_ONESIDED_REQS)=0
        do i=1,MAX_ONESIDED_REQS; call this%RankWins(i)%init(); enddo !init all entries to null
        if(present(ierr)) ierr=0
        return
        end subroutine RankWinListClean
!---------------------------------------------------------------------------
        integer(INT_MPI) function RankWinListTest(this,rank,win,ierr,append)
!Looks up a given pair (rank,window) in the active (rank,window) list and
!returns the corresponding entry number (>0), if found. If not found, returns 0.
!If <append> is present and TRUE, a new entry will be created (only if not found)
!and its number will be returned. In case when there are no free entries in the table,
!a special return status TRY_LATER will be returned to postpone the request for later.
        implicit none
        class(RankWinList_t), intent(inout):: this       !inout: (rank,window) list
        integer(INT_MPI), intent(in):: rank              !in: MPI rank
        integer(INT_MPI), intent(in):: win               !in: MPI window handle
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        logical, intent(in), optional:: append           !in: if .true., a new entry will be appended if not found
        integer(INT_MPI):: i,j,m,n,errc
        logical:: apnd

        errc=0; RankWinListTest=0
        if(present(append)) then; apnd=append; else; apnd=.false.; endif
        if(rank.ge.0) then
         n=mod(rank,HASH_MOD); i=this%HashBin(n); m=0
         do while(i.ge.1.and.i.le.MAX_ONESIDED_REQS)
          if(this%RankWins(i)%Rank.lt.rank) then
           m=i; i=this%NextEntry(i)
          elseif(this%RankWins(i)%Rank.eq.rank) then
           if(this%RankWins(i)%Window.lt.win) then
            m=i; i=this%NextEntry(i)
           elseif(this%RankWins(i)%Window.eq.win) then !match found
            RankWinListTest=i
            exit
           else
            exit
           endif
          else
           exit
          endif
         enddo
         if(RankWinListTest.le.0.and.apnd) then !append if not found
          if(this%NumEntries.lt.MAX_ONESIDED_REQS) then
           j=this%FirstFree; this%FirstFree=this%NextEntry(j)
           this%PrevEntry(this%FirstFree)=0
           this%NumEntries=this%NumEntries+1
           if(m.gt.0) then !previous entry exists
            this%NextEntry(m)=j
           else !no previous entry
            this%HashBin(n)=j
           endif
           this%PrevEntry(j)=m; this%NextEntry(j)=i
           if(i.gt.0) this%PrevEntry(i)=j
           call this%RankWins(j)%init(rank,win,errc)
           if(errc.eq.0) then
            RankWinListTest=j
            if(DEBUG) write(jo,'("#DEBUG(distributed::RankWinList.Test)[",i7,"]: New (rank,window) registered: ",i9,1x,i13)')&
                      &impir,rank,win
           else
            call this%delete(j); errc=1
           endif
          else
           errc=TRY_LATER !no more free entries, try later (not an error)
          endif
         endif
        else
         errc=2
        endif
        if(present(ierr)) ierr=errc
        return
        end function RankWinListTest
!-----------------------------------------------------
        subroutine RankWinListDel(this,entry_num,ierr)
!Deletes entry <entry_num> from the active (rank,window) list.
        implicit none
        class(RankWinList_t), intent(inout):: this       !inout: (rank,window) list
        integer(INT_MPI), intent(in):: entry_num         !in: number of the entry to be deleted
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: i,j,m,errc

        errc=0
        if(entry_num.ge.1.and.entry_num.le.MAX_ONESIDED_REQS) then
         if(this%RankWins(entry_num)%Rank.ge.0) then !active entry
          m=mod(this%RankWins(entry_num)%Rank,HASH_MOD)
          i=this%PrevEntry(entry_num)
          j=this%NextEntry(entry_num)
          if(i.gt.0) this%NextEntry(i)=j
          if(j.gt.0) this%PrevEntry(j)=i
          if(this%HashBin(m).eq.entry_num) this%HashBin(m)=j
          this%NextEntry(entry_num)=this%FirstFree
          if(this%FirstFree.gt.0) this%PrevEntry(this%FirstFree)=entry_num
          this%FirstFree=entry_num; this%NumEntries=this%NumEntries-1
          if(DEBUG) write(jo,'("#DEBUG(distributed::RankWinList.Del)[",i7,"]: (rank,window) deleted: ",i9,1x,i13)')&
                    &impir,this%RankWins(entry_num)%Rank,this%RankWins(entry_num)%Window
          call this%RankWins(entry_num)%init() !clean entry
         else
          errc=1 !empty entries cannot be deleted
         endif
        else
         errc=2
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine RankWinListDel
!---------------------------------------------------
        subroutine RankWinNewTrans(this,dd,rwe,ierr)
!This subroutine increments the global transfer ID and communicated data size
!and assigns the current transfer ID to the given data descriptor.
        implicit none
        class(RankWinList_t), intent(inout):: this       !inout: (rank,window) list
        class(DataDescr_t), intent(inout):: dd           !inout: valid data descriptor
        integer(INT_MPI), intent(in):: rwe               !in: the corresponding (rank,window) entry
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(8), parameter:: int8=0_8
        integer(INT_MPI):: n,errc

        errc=0
        if(this%TransCount.lt.huge(int8)) then
         if(rwe.ge.1.and.rwe.le.MAX_ONESIDED_REQS) then
          this%TransCount=this%TransCount+1
          n=data_type_size(dd%DataType)
          this%TransSize=this%TransSize+real(n,8)*real(dd%DataVol,8)
          dd%TransID=this%TransCount !global transaction ID
          this%RankWins(rwe)%RefCount=this%RankWins(rwe)%RefCount+1
          if(DEBUG) write(jo,'("#DEBUG(distributed::RankWin.NewTrans)[",i7,"]: New transfer: ",i9,1x,i13,1x,i5,1x,i13)')&
                    &impir,this%RankWins(rwe)%Rank,this%RankWins(rwe)%Window,this%RankWins(rwe)%RefCount,dd%TransID
         else
          errc=1
         endif
        else
         if(VERBOSE) write(CONS_OUT,'("#FATAL(distributed::RankWin.NewTrans): Max int8 MPI message count exceeded!")')
         errc=2
         call quit(-1,'#FATAL: Unable to continue: Distributed Data Service failed!')
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine RankWinNewTrans
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
            if(DEBUG) write(jo,'("#DEBUG(distributed::DataWin.Create)[",i7,"]: MPI window created: ",i11,1x,i11,1x,i2,1x,i11)')&
                      &impir,this%WinMPI%Window,this%WinMPI%CommMPI,this%WinMPI%DispUnit,this%WinSize
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
           if(DEBUG) write(jo,'("#DEBUG(distributed::DataWin.Destroy)[",i7,"]: MPI window destroyed: ",i11,1x,i11,1x,i2,1x,i11)')&
                     &impir,this%WinMPI%Window,this%WinMPI%CommMPI,this%WinMPI%DispUnit,this%WinSize
           call this%clean(errc)
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
!The size of the data buffer in bytes must be a multiple of 4.
        implicit none
        class(DataWin_t), intent(inout):: this           !inout: data window
        type(C_PTR), intent(in):: loc_ptr                !in: C pointer to the local data buffer
        integer(INT_ADDR), intent(in):: loc_size         !in: size of the local data buffer in bytes: Must be multiple of 4
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc
        integer(INT_ADDR):: vol
        real(4), pointer, contiguous:: r4_ptr(:) !the (arbitrary) data will be mapped as a real(4) array

        errc=0
        if(.not.c_associated(loc_ptr,C_NULL_PTR)) then
         if(loc_size.gt.0.and.mod(loc_size,4).eq.0) then
          vol=loc_size/4 !volume in 4-byte words
          if(this%WinSize.ge.0) then !data window has been initialized
           call c_f_pointer(loc_ptr,r4_ptr,(/vol/))
           call attach_buffer(r4_ptr,loc_size,errc)
           if(errc.eq.0) then
            this%WinSize=this%WinSize+loc_size
            if(DEBUG) write(jo,'("#DEBUG(distributed::DataWin.Attach)[",i7,"]: Buffer attached: ",i11,1x,i11,1x,i11)')&
                      &impir,loc_size,this%WinMPI%Window,this%WinSize
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
         real(4), intent(in):: r4_arr(*) !`asynchronous attribute
         integer(INT_ADDR), intent(in):: bsize
         integer(INT_MPI), intent(out):: jerr
         jerr=0
         call MPI_WIN_ATTACH(this%WinMPI%Window,r4_arr,bsize,jerr)
         return
         end subroutine attach_buffer

        end subroutine DataWinAttach
!-----------------------------------------------------------
        subroutine DataWinDetach(this,loc_ptr,loc_size,ierr)
!Detaches a previously attached local (contiguous) data buffer from an MPI data window.
!The size of the data buffer in bytes must be a multiple of 4.
        implicit none
        class(DataWin_t), intent(inout):: this           !inout: data window
        type(C_PTR), intent(in):: loc_ptr                !in: C pointer to the local data buffer to be detached
        integer(INT_ADDR), intent(in):: loc_size         !in: size of the local data buffer in bytes: Must be multiple of 4
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc
        integer(INT_ADDR):: vol
        real(4), pointer, contiguous:: r4_ptr(:) !the (arbitrary) data is mapped as a real(4) array

        errc=0
        if(.not.c_associated(loc_ptr,C_NULL_PTR)) then
         if(loc_size.gt.0.and.mod(loc_size,4).eq.0) then
          vol=loc_size/4 !volume in 4-byte words
          if(this%WinSize.ge.loc_size) then
           call c_f_pointer(loc_ptr,r4_ptr,(/vol/))
           call detach_buffer(r4_ptr,errc)
           if(errc.eq.0) then
            this%WinSize=this%WinSize-loc_size
            if(DEBUG) write(jo,'("#DEBUG(distributed::DataWin.Detach)[",i7,"]: Buffer detached: ",i11,1x,i11,1x,i11)')&
                      &impir,loc_size,this%WinMPI%Window,this%WinSize
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
         real(4), intent(in):: r4_arr(*) !`asynchronous attribute
         integer(INT_MPI), intent(out):: jerr
         jerr=0
         call MPI_WIN_DETACH(this%WinMPI%Window,r4_arr,jerr)
         return
         end subroutine detach_buffer

        end subroutine DataWinDetach
!----------------------------------------
        subroutine DataWinSync(this,ierr)
!Synchronizes the local private and public copy of an MPI data window.
        implicit none
        class(DataWin_t), intent(in):: this              !in: data window
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0
        if(this%WinSize.ge.0) then !data window is initialized
         call MPI_WIN_SYNC(this%WinMPI%Window,errc)
         if(errc.eq.0) then
          if(DEBUG) write(jo,'("#DEBUG(distributed::DataWin.Sync)[",i7,"]: MPI window synced: ",i11,1x,i11)')&
                    &impir,this%WinMPI%Window,this%WinSize
         else
          errc=1
         endif
        else
         errc=2
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine DataWinSync
!==========================================================================
        subroutine DistrSpaceCreate(this,comm_mpi,num_wins,space_name,ierr) !COLLECTIVE
!Creates a distributed memory space with <num_wins> dynamic windows.
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
             if(DEBUG)write(jo,'("#DEBUG(distributed::DistrSpace.Create)[",i7,"]: Distributed space created: ",i11,1x,i4,1x,A16)')&
                       &impir,this%CommMPI,this%NumWins,this%SpaceName(1:min(DISTR_SPACE_NAME_LEN,16))
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
!Also the distributed memory space itself must be empty (0-byte size)!
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
              if(DEBUG)&
               &write(jo,'("#DEBUG(distributed::DistrSpace.Destroy)[",i7,"]: Distributed space destroyed: ",i11,1x,i4,1x,A16)')&
               &impir,this%CommMPI,this%NumWins,this%SpaceName(1:min(DISTR_SPACE_NAME_LEN,16))
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
!belonging to a given distributed space.
        implicit none
        class(DistrSpace_t), intent(in):: this           !in: distributed memory space
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: i,errc
        integer(INT_ADDR):: loc_size

        errc=0; DistrSpaceLocalSize=-1
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
        else
         errc=2
        endif
        if(errc.eq.0) DistrSpaceLocalSize=loc_size
        if(present(ierr)) ierr=errc
        return
        end function DistrSpaceLocalSize
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
!Detaches a previously attached local (contiguous) data buffer from a distributed memory space.
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
            do i=1,this%NumWins !find the MPI window
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
         if(.not.c_associated(loc_ptr,C_NULL_PTR)) then
          if(data_vol.gt.0) then
           elem_size=data_type_size(data_type,errc)
           if(errc.eq.0.and.elem_size.gt.0) then
            this%LocPtr=loc_ptr
            this%RankMPI=process_rank
            this%WinMPI=win_mpi
            this%DataVol=data_vol
            this%DataType=data_type
            this%TransID=0 !0 means "never assigned a Transfer ID"
            this%StatMPI=MPI_STAT_NONE
            this%ReqHandle=MPI_REQUEST_NULL
            call MPI_Get_Displacement(loc_ptr,this%Offset,errc); if(errc.ne.0) errc=1
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
        if(errc.ne.0) call this%clean()
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrInit
!-----------------------------------------------------
        subroutine DataDescrFlushData(this,ierr,local)
!Completes an MPI one-sided communication specified by a given data descriptor.
!Other outstanding MPI communications on the same (rank,window) will be completed as well.
!It is safe to flush the same (valid) data descriptor multiple times. However,
!if the data descriptor was flushed the first time at the origin only (local=.true.),
!one will no longer be able to request a full completion (both at origin and target).
!In case of a one-sided synchronization error, the data request will still be marked
!as pending, thus blocking the corresponding (rank,window) entry from release!
        implicit none
        class(DataDescr_t), intent(inout):: this         !inout: data descriptor
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        logical, intent(in), optional:: local            !in: if .true., the data is flushed only at the origin
        integer(INT_MPI):: rwe,errc
        type(RankWin_t), pointer:: rw_entry
        logical:: lcl

        errc=0
        if(present(local)) then; lcl=local; else; lcl=.false.; endif !default is flushing both at the origin and target
        if(this%RankMPI.ge.0) then
         if(this%StatMPI.eq.MPI_STAT_PROGRESS_NRM.or.this%StatMPI.eq.MPI_STAT_PROGRESS_REQ) then !first flush
          rwe=RankWinRefs%test(this%RankMPI,this%WinMPI%Window,errc) !get the (rank,window) entry
          if(rwe.gt.0.and.errc.eq.0) then
           rw_entry=>RankWinRefs%RankWins(rwe)
           if(rw_entry%RefCount.gt.1) then !not the last reference to this (rank,window)
            if(lcl) then !complete at origin only
             if(this%TransID.gt.rw_entry%LastSync) call MPI_WIN_FLUSH_LOCAL(rw_entry%Rank,rw_entry%Window,errc) !complete at the origin only
             if(errc.eq.0) then
              this%StatMPI=MPI_STAT_COMPLETED_ORIG
              rw_entry%RefCount=rw_entry%RefCount-1
             else
              this%StatMPI=MPI_STAT_ONESIDED_ERR; errc=1
             endif
            else !complete at both origin and target
             if(this%TransID.gt.rw_entry%LastSync) call MPI_WIN_FLUSH(rw_entry%Rank,rw_entry%Window,errc) !complete both at the origin and target
             if(errc.eq.0) then
              this%StatMPI=MPI_STAT_COMPLETED
              rw_entry%RefCount=rw_entry%RefCount-1
             else
              this%StatMPI=MPI_STAT_ONESIDED_ERR; errc=2
             endif
            endif
           else !the last reference to this (rank,window)
            if(this%TransID.gt.rw_entry%LastSync) call MPI_WIN_UNLOCK(rw_entry%Rank,rw_entry%Window,errc) !complete both at origin and target
            if(errc.eq.0) then
             this%StatMPI=MPI_STAT_COMPLETED
             rw_entry%RefCount=rw_entry%RefCount-1
            else
             this%StatMPI=MPI_STAT_ONESIDED_ERR; errc=3
            endif
           endif
           if(errc.eq.0) rw_entry%LastSync=RankWinRefs%TransCount
           if(rw_entry%RefCount.eq.0) then !delete the (rank,window) entry if no references are attached to it
            call RankWinRefs%delete(rwe,errc); if(errc.ne.0) errc=4
           elseif(rw_entry%RefCount.lt.0) then
            if(VERBOSE) write(CONS_OUT,'("#FATAL(distributed::DataDescr.FlushData): Negative reference count: ",i12)')&
                        &rw_entry%RefCount
            errc=5
            call quit(-1,'#FATAL: Unable to continue: Distributed Data Service failed!')
           endif
           nullify(rw_entry)
          else
           errc=6
          endif
         elseif(this%StatMPI.eq.MPI_STAT_NONE.or.this%StatMPI.eq.MPI_STAT_ONESIDED_ERR) then !invalid data descriptor
          errc=7
         endif
        else
         errc=8
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrFlushData
!----------------------------------------------------
        logical function DataDescrTestData(this,ierr)
!Tests a completion (at the origin) of a data request with a request handle.
!In case of a one-sided synchronization error, the data request will still be marked
!as pending, thus blocking the corresponding (rank,window) entry from release!
        implicit none
        class(DataDescr_t), intent(inout):: this         !inout: data descriptor
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: mpi_stat(MPI_STATUS_SIZE),rwe,errc
        type(RankWin_t), pointer:: rw_entry
        logical:: compl

        errc=0; DataDescrTestData=.false.
        if(this%RankMPI.ge.0) then
         if(this%StatMPI.eq.MPI_STAT_PROGRESS_REQ) then
          rwe=RankWinRefs%test(this%RankMPI,this%WinMPI%Window,errc) !get the (rank,window) entry
          if(rwe.gt.0.and.errc.eq.0) then
           rw_entry=>RankWinRefs%RankWins(rwe)
           if(rw_entry%RefCount.gt.0) then
            if(this%TransID.gt.rw_entry%LastSync) then
             call MPI_TEST(this%ReqHandle,compl,mpi_stat,errc)
             if(errc.eq.0) then
              if(compl) then
               this%StatMPI=MPI_STAT_COMPLETED_ORIG
               rw_entry%RefCount=rw_entry%RefCount-1
               DataDescrTestData=.true.
              endif
             else
              this%StatMPI=MPI_STAT_ONESIDED_ERR; errc=1
             endif
            else
             this%StatMPI=MPI_STAT_COMPLETED_ORIG
             rw_entry%RefCount=rw_entry%RefCount-1
             DataDescrTestData=.true.
            endif
            if(rw_entry%RefCount.eq.0) then !delete the (rank,window) entry if no references are attached to it
             call RankWinRefs%delete(rwe,errc); if(errc.ne.0) errc=2
            endif
            nullify(rw_entry)
           else
            if(VERBOSE) write(CONS_OUT,'("#FATAL(distributed::DataDescr.TestData): Invalid reference count: ",i12)')&
                        &rw_entry%RefCount
            errc=3
            call quit(-1,'#FATAL: Unable to continue: Distributed Data Service failed!')
           endif
          else
           errc=4
          endif
         elseif(this%StatMPI.eq.MPI_STAT_COMPLETED.or.this%StatMPI.eq.MPI_STAT_COMPLETED_ORIG) then
          DataDescrTestData=.true.
         else
          errc=5
         endif
        else
         errc=6
        endif
        if(present(ierr)) ierr=errc
        return
        end function DataDescrTestData
!----------------------------------------------
        subroutine DataDescrWaitData(this,ierr)
!Waits for a completion of a data request with a request handle.
!In case of a one-sided synchronization error, the data request will still be marked
!as pending, thus blocking the corresponding (rank,window) entry from release!
        implicit none
        class(DataDescr_t), intent(inout):: this         !inout: data descriptor
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: mpi_stat(MPI_STATUS_SIZE),rwe,errc
        type(RankWin_t), pointer:: rw_entry

        errc=0
        if(this%RankMPI.ge.0) then
         if(this%StatMPI.eq.MPI_STAT_PROGRESS_REQ) then
          rwe=RankWinRefs%test(this%RankMPI,this%WinMPI%Window,errc) !get the (rank,window) entry
          if(rwe.gt.0.and.errc.eq.0) then
           rw_entry=>RankWinRefs%RankWins(rwe)
           if(rw_entry%RefCount.gt.0) then
            if(this%TransID.gt.rw_entry%LastSync) then
             call MPI_WAIT(this%ReqHandle,mpi_stat,errc)
             if(errc.eq.0) then
              this%StatMPI=MPI_STAT_COMPLETED_ORIG
              rw_entry%RefCount=rw_entry%RefCount-1
             else
              this%StatMPI=MPI_STAT_ONESIDED_ERR; errc=1
             endif
            else
             this%StatMPI=MPI_STAT_COMPLETED_ORIG
             rw_entry%RefCount=rw_entry%RefCount-1
            endif
            if(rw_entry%RefCount.eq.0) then !delete the (rank,window) entry if no references are attached to it
             call RankWinRefs%delete(rwe,errc); if(errc.ne.0) errc=2
            endif
            nullify(rw_entry)
           else
            if(VERBOSE) write(CONS_OUT,'("#FATAL(distributed::DataDescr.WaitData): Invalid reference count: ",i12)')&
                        &rw_entry%RefCount
            errc=3
            call quit(-1,'#FATAL: Unable to continue: Distributed Data Service failed!')
           endif
          else
           errc=4
          endif
         elseif(this%StatMPI.ne.MPI_STAT_COMPLETED.and.this%StatMPI.ne.MPI_STAT_COMPLETED_ORIG) then
          errc=5
         endif
        else
         errc=6
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrWaitData
!-----------------------------------------------------------
        subroutine DataDescrGetData(this,loc_ptr,ierr,async)
!Initiates a (remote) data fetch into a local buffer.
!The default (synchronous) call will complete the transfer within this call.
!If <async> is present and equal to MPI_ASYNC_NRM or MPI_ASYNC_REQ, then
!another (completion) call will be required to complete the communication.
!If the internal tables for communication tracking no longer have free entries,
!the status TRY_LATER is returned, meaning that one needs to wait until later.
        implicit none
        class(DataDescr_t), intent(inout):: this         !inout: data descriptor
        type(C_PTR), intent(in):: loc_ptr                !in: pointer to a local buffer
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success, TRY_LATER:resource is currently busy)
        integer(INT_MPI), intent(in), optional:: async   !in: asynchronisity: {MPI_ASYNC_NOT,MPI_ASYNC_NRM,MPI_ASYNC_REQ}
        integer(INT_MPI), parameter:: MPI_ASSER=0        !MPI lock assertion
        integer(INT_MPI):: rwe,errc,asnc
        real(4), pointer, contiguous:: r4_ptr(:)
        real(8), pointer, contiguous:: r8_ptr(:)
        complex(8), pointer, contiguous:: c8_ptr(:)

        errc=0
        if(present(async)) then; asnc=async; else; asnc=MPI_ASYNC_NOT; endif !default is synchronous communication
        if(asnc.eq.MPI_ASYNC_NOT.or.asnc.eq.MPI_ASYNC_NRM.or.asnc.eq.MPI_ASYNC_REQ) then
         if(.not.c_associated(loc_ptr,C_NULL_PTR)) then
          if(this%RankMPI.ge.0) then
           if(this%StatMPI.eq.MPI_STAT_NONE.or.this%StatMPI.eq.MPI_STAT_COMPLETED.or.&
             &this%StatMPI.eq.MPI_STAT_COMPLETED_ORIG) then
            rwe=RankWinRefs%test(this%RankMPI,this%WinMPI%Window,errc,append=.true.) !get the (rank,window) entry
            if(errc.eq.0) then
             if(this%DataVol.gt.0) then
              if(.not.(asnc.eq.MPI_ASYNC_REQ.and.this%DataVol.gt.huge(asnc))) then
               call modify_lock(rwe,errc) !modify the (rank,window) lock status if needed
               if(errc.eq.0) then
                select case(this%DataType)
                case(R4)
                 call c_f_pointer(loc_ptr,r4_ptr,(/this%DataVol/))
                 call start_get_r4(r4_ptr,errc); if(errc.ne.0) errc=1
                case(R8)
                 call c_f_pointer(loc_ptr,r8_ptr,(/this%DataVol/))
                 call start_get_r8(r8_ptr,errc); if(errc.ne.0) errc=2
                case(C8)
                 call c_f_pointer(loc_ptr,c8_ptr,(/this%DataVol/))
                 call start_get_c8(c8_ptr,errc); if(errc.ne.0) errc=3
                case(NO_TYPE)
                 errc=4
                case default
                 errc=5
                end select
                if(errc.eq.0) then
                 call RankWinRefs%new_transfer(this,rwe) !register a new transfer (will also set this%TransID field)
                 if(asnc.eq.MPI_ASYNC_NOT) then
                  call this%flush_data(errc); if(errc.ne.0) errc=6
                 endif
                endif
               else
                errc=7
               endif
              else
               errc=8
              endif
             elseif(this%DataVol.eq.0) then
              this%StatMPI=MPI_STAT_COMPLETED
             else
              errc=9
             endif
            else
             if(errc.ne.TRY_LATER) errc=10
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
        else
         errc=14
        endif
        if(present(ierr)) ierr=errc
        return

        contains

         subroutine modify_lock(rw,jerr)
         integer(INT_MPI), intent(in):: rw
         integer(INT_MPI), intent(out):: jerr
         type(RankWin_t), pointer:: rw_entry
         jerr=0; rw_entry=>RankWinRefs%RankWins(rw)
         if(rw_entry%LockType*READ_SIGN.lt.0) then !communication direction change
          call MPI_WIN_UNLOCK(rw_entry%Rank,rw_entry%Window,jerr)
          if(jerr.eq.0) then
           rw_entry%LockType=NO_LOCK
           rw_entry%LastSync=RankWinRefs%TransCount
          else
           jerr=1
          endif
         endif
         if(jerr.eq.0.and.rw_entry%LockType.eq.NO_LOCK) then
          call MPI_WIN_LOCK(MPI_LOCK_SHARED,rw_entry%Rank,MPI_ASSER,rw_entry%Window,jerr)
          if(jerr.eq.0) then
           rw_entry%LockType=SHARED_LOCK*READ_SIGN
          else
           jerr=2
          endif
         endif
         nullify(rw_entry)
         return
         end subroutine modify_lock

         subroutine start_get_r4(r4_arr,jerr)
         real(4), intent(inout):: r4_arr(*) !asynchronous
         integer(INT_MPI), intent(out):: jerr
         integer(INT_COUNT):: ji,js
         integer(INT_MPI):: jv
         jerr=0; ji=int(MAX_MPI_MSG_VOL,INT_COUNT)
         if(asnc.ne.MPI_ASYNC_REQ) then !regular
          do js=1,this%DataVol,ji
           jv=int(min(this%DataVol-js+1,ji),INT_MPI)
           call MPI_GET(r4_arr(js:js+jv-1),jv,MPI_REAL4,this%RankMPI,this%Offset,jv,MPI_REAL4,this%WinMPI%Window,jerr)
           if(jerr.ne.0) exit
          enddo
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_NRM
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=1
          endif
         else !request-handle
          jv=this%DataVol
          call MPI_RGET(r4_arr,jv,MPI_REAL4,this%RankMPI,this%Offset,jv,MPI_REAL4,this%WinMPI%Window,this%ReqHandle,jerr)
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_REQ
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=2
          endif
         endif
         return
         end subroutine start_get_r4

         subroutine start_get_r8(r8_arr,jerr)
         real(8), intent(inout):: r8_arr(*) !asynchronous
         integer(INT_MPI), intent(out):: jerr
         integer(INT_COUNT):: ji,js
         integer(INT_MPI):: jv
         jerr=0; ji=int(MAX_MPI_MSG_VOL,INT_COUNT)
         if(asnc.ne.MPI_ASYNC_REQ) then !regular
          do js=1,this%DataVol,ji
           jv=int(min(this%DataVol-js+1,ji),INT_MPI)
           call MPI_GET(r8_arr(js:js+jv-1),jv,MPI_REAL8,this%RankMPI,this%Offset,jv,MPI_REAL8,this%WinMPI%Window,jerr)
           if(jerr.ne.0) exit
          enddo
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_NRM
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=1
          endif
         else !request-handle
          jv=this%DataVol
          call MPI_RGET(r8_arr,jv,MPI_REAL8,this%RankMPI,this%Offset,jv,MPI_REAL8,this%WinMPI%Window,this%ReqHandle,jerr)
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_REQ
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=2
          endif
         endif
         return
         end subroutine start_get_r8

         subroutine start_get_c8(c8_arr,jerr)
         complex(8), intent(inout):: c8_arr(*) !asynchronous
         integer(INT_MPI), intent(out):: jerr
         integer(INT_COUNT):: ji,js
         integer(INT_MPI):: jv
         jerr=0; ji=int(MAX_MPI_MSG_VOL,INT_COUNT)
         if(asnc.ne.MPI_ASYNC_REQ) then !regular
          do js=1,this%DataVol,ji
           jv=int(min(this%DataVol-js+1,ji),INT_MPI)
           call MPI_GET(c8_arr(js:js+jv-1),jv,MPI_COMPLEX8,this%RankMPI,this%Offset,jv,MPI_COMPLEX8,this%WinMPI%Window,jerr)
           if(jerr.ne.0) exit
          enddo
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_NRM
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=1
          endif
         else !request-handle
          jv=this%DataVol
          call MPI_RGET(c8_arr,jv,MPI_COMPLEX8,this%RankMPI,this%Offset,jv,MPI_COMPLEX8,this%WinMPI%Window,this%ReqHandle,jerr)
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_REQ
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=2
          endif
         endif
         return
         end subroutine start_get_c8

        end subroutine DataDescrGetData
!-----------------------------------------------------------
        subroutine DataDescrAccData(this,loc_ptr,ierr,async)
!Initiates a (remote) data accumulate from a local buffer.
!The default (synchronous) call will complete the transfer within this call.
!If <async> is present and equal to MPI_ASYNC_NRM or MPI_ASYNC_REQ, then
!another (completion) call will be required to complete the communication.
!If the internal tables for communication tracking no longer have free entries,
!the status TRY_LATER is returned, meaning that one needs to wait until later.
        implicit none
        class(DataDescr_t), intent(inout):: this         !inout: data descriptor
        type(C_PTR), intent(in):: loc_ptr                !in: pointer to a local buffer
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success, TRY_LATER:resource is currently busy)
        integer(INT_MPI), intent(in), optional:: async   !in: asynchronisity: {MPI_ASYNC_NOT,MPI_ASYNC_NRM,MPI_ASYNC_REQ}
        integer(INT_MPI), parameter:: MPI_ASSER=0        !MPI lock assertion
        integer(INT_MPI):: rwe,errc,asnc
        real(4), pointer, contiguous:: r4_ptr(:)
        real(8), pointer, contiguous:: r8_ptr(:)
        complex(8), pointer, contiguous:: c8_ptr(:)

        errc=0
        if(present(async)) then; asnc=async; else; asnc=MPI_ASYNC_NOT; endif !default is synchronous communication
        if(asnc.eq.MPI_ASYNC_NOT.or.asnc.eq.MPI_ASYNC_NRM.or.asnc.eq.MPI_ASYNC_REQ) then
         if(.not.c_associated(loc_ptr,C_NULL_PTR)) then
          if(this%RankMPI.ge.0) then
           if(this%StatMPI.eq.MPI_STAT_NONE.or.this%StatMPI.eq.MPI_STAT_COMPLETED.or.&
              this%StatMPI.eq.MPI_STAT_COMPLETED_ORIG) then
            rwe=RankWinRefs%test(this%RankMPI,this%WinMPI%Window,errc,append=.true.) !get the (rank,window) entry
            if(errc.eq.0) then
             if(this%DataVol.gt.0) then
              if(.not.(asnc.eq.MPI_ASYNC_REQ.and.this%DataVol.gt.huge(asnc))) then
               call modify_lock(rwe,errc) !modify the (rank,window) lock status if needed
               if(errc.eq.0) then
                select case(this%DataType)
                case(R4)
                 call c_f_pointer(loc_ptr,r4_ptr,(/this%DataVol/))
                 call start_acc_r4(r4_ptr,errc); if(errc.ne.0) errc=1
                case(R8)
                 call c_f_pointer(loc_ptr,r8_ptr,(/this%DataVol/))
                 call start_acc_r8(r8_ptr,errc); if(errc.ne.0) errc=2
                case(C8)
                 call c_f_pointer(loc_ptr,c8_ptr,(/this%DataVol/))
                 call start_acc_c8(c8_ptr,errc); if(errc.ne.0) errc=3
                case(NO_TYPE)
                 errc=4
                case default
                 errc=5
                end select
                if(errc.eq.0) then
                 call RankWinRefs%new_transfer(this,rwe) !register a new transfer (will also set this%TransID field)
                 if(asnc.eq.MPI_ASYNC_NOT) then
                  call this%flush_data(errc); if(errc.ne.0) errc=6
                 endif
                endif
               else
                errc=7
               endif
              else
               errc=8
              endif
             elseif(this%DataVol.eq.0) then
              this%StatMPI=MPI_STAT_COMPLETED
             else
              errc=9
             endif
            else
             if(errc.ne.TRY_LATER) errc=10
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
        else
         errc=14
        endif
        if(present(ierr)) ierr=errc
        return

        contains

         subroutine modify_lock(rw,jerr)
         integer(INT_MPI), intent(in):: rw
         integer(INT_MPI), intent(out):: jerr
         type(RankWin_t), pointer:: rw_entry
         jerr=0; rw_entry=>RankWinRefs%RankWins(rw)
         if(rw_entry%LockType*WRITE_SIGN.lt.0) then !communication direction change
          call MPI_WIN_UNLOCK(rw_entry%Rank,rw_entry%Window,jerr)
          if(jerr.eq.0) then
           rw_entry%LockType=NO_LOCK
           rw_entry%LastSync=RankWinRefs%TransCount
          else
           jerr=1
          endif
         endif
         if(jerr.eq.0.and.rw_entry%LockType.eq.NO_LOCK) then
          call MPI_WIN_LOCK(MPI_LOCK_SHARED,rw_entry%Rank,MPI_ASSER,rw_entry%Window,jerr)
          if(jerr.eq.0) then
           rw_entry%LockType=SHARED_LOCK*WRITE_SIGN
          else
           jerr=2
          endif
         endif
         nullify(rw_entry)
         return
         end subroutine modify_lock

         subroutine start_acc_r4(r4_arr,jerr)
         real(4), intent(inout):: r4_arr(*) !asynchronous
         integer(INT_MPI), intent(out):: jerr
         integer(INT_COUNT):: ji,js
         integer(INT_MPI):: jv
         jerr=0; ji=int(MAX_MPI_MSG_VOL,INT_COUNT)
         if(asnc.ne.MPI_ASYNC_REQ) then !regular
          do js=1,this%DataVol,ji
           jv=int(min(this%DataVol-js+1,ji),INT_MPI)
           call MPI_ACCUMULATE(r4_arr(js:js+jv-1),jv,MPI_REAL4,this%RankMPI,this%Offset,jv,MPI_REAL4,MPI_SUM,&
                              &this%WinMPI%Window,jerr)
           if(jerr.ne.0) exit
          enddo
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_NRM
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=1
          endif
         else !request-handle
          jv=this%DataVol
          call MPI_RACCUMULATE(r4_arr,jv,MPI_REAL4,this%RankMPI,this%Offset,jv,MPI_REAL4,MPI_SUM,&
                              &this%WinMPI%Window,this%ReqHandle,jerr)
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_REQ
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=2
          endif
         endif
         return
         end subroutine start_acc_r4

         subroutine start_acc_r8(r8_arr,jerr)
         real(8), intent(inout):: r8_arr(*) !asynchronous
         integer(INT_MPI), intent(out):: jerr
         integer(INT_COUNT):: ji,js
         integer(INT_MPI):: jv
         jerr=0; ji=int(MAX_MPI_MSG_VOL,INT_COUNT)
         if(asnc.ne.MPI_ASYNC_REQ) then !regular
          do js=1,this%DataVol,ji
           jv=int(min(this%DataVol-js+1,ji),INT_MPI)
           call MPI_ACCUMULATE(r8_arr(js:js+jv-1),jv,MPI_REAL8,this%RankMPI,this%Offset,jv,MPI_REAL8,MPI_SUM,&
                              &this%WinMPI%Window,jerr)
           if(jerr.ne.0) exit
          enddo
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_NRM
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=1
          endif
         else !request-handle
          jv=this%DataVol
          call MPI_RACCUMULATE(r8_arr,jv,MPI_REAL8,this%RankMPI,this%Offset,jv,MPI_REAL8,MPI_SUM,&
                              &this%WinMPI%Window,this%ReqHandle,jerr)
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_REQ
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=2
          endif
         endif
         return
         end subroutine start_acc_r8

         subroutine start_acc_c8(c8_arr,jerr)
         complex(8), intent(inout):: c8_arr(*) !asynchronous
         integer(INT_MPI), intent(out):: jerr
         integer(INT_COUNT):: ji,js
         integer(INT_MPI):: jv
         jerr=0; ji=int(MAX_MPI_MSG_VOL,INT_COUNT)
         if(asnc.ne.MPI_ASYNC_REQ) then !regular
          do js=1,this%DataVol,ji
           jv=int(min(this%DataVol-js+1,ji),INT_MPI)
           call MPI_ACCUMULATE(c8_arr(js:js+jv-1),jv,MPI_COMPLEX8,this%RankMPI,this%Offset,jv,MPI_COMPLEX8,MPI_SUM,&
                              &this%WinMPI%Window,jerr)
           if(jerr.ne.0) exit
          enddo
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_NRM
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=1
          endif
         else !request-handle
          jv=this%DataVol
          call MPI_RACCUMULATE(c8_arr,jv,MPI_COMPLEX8,this%RankMPI,this%Offset,jv,MPI_COMPLEX8,MPI_SUM,&
                              &this%WinMPI%Window,this%ReqHandle,jerr)
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_REQ
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=2
          endif
         endif
         return
         end subroutine start_acc_c8

        end subroutine DataDescrAccData

        end module distributed
