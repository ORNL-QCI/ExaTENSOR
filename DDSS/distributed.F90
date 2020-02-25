!Distributed data storage service (DDSS).
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2020/02/24 (started 2015/03/18)

!Copyright (C) 2014-2020 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2020 Oak Ridge National Laboratory (UT-Battelle)

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

!IMPORTANT NOTES:
! * Cloning DataDescr_t, either as a standalone object or as a component of a larger object,
!   must employ DataDescr_t.clone() method. Otherwise, for example when using sourced memory
!   allocation, one must subsequently call DataDescr_t.clear_lock() method on the just
!   created DataDescr_t clone!
!CONCEPTS:
! * Each MPI process can participate in one or more distributed memory spaces (DMS),
!   where each distributed memory space is defined within a specific MPI communicator.
!   Internally, within each distributed memory space, each MPI process opens a number
!   of dynamic MPI windows for data storage, which are opaque to the user.
! * Whenever new persistent data (e.g., a tensor block) is allocated by an MPI process,
!   it can be attached to one of the dynamic MPI windows in a specific distributed memory space
!   and the corresponding data descriptor (DD) should be sent to the manager for registration.
!   The registered data descriptor can be used by other MPI processes to remotely access the data.
!   The size of the data being attached to a distributed space must be aligned to DDSS_BUFFER_ALIGN.
! * A data descriptor is communicated between MPI processes as a plain integer packet,
!   with integer kind ELEM_PACK_SIZE. The first integer in the packet always contains
!   the length of the rest of the packet (number of integer elements carrying the information).
!   Besides, multiple data descriptors can be communicated simultaneously in a single super-packet,
!   which consists of the first integer specifying the number of simple packets in the super-packet
!   followed by the simple packets themselves. Note that the data packed in packets, including the
!   packet length (first element), cannot be read directly, but requires special access functions!
!   The number of packets in a super-packet can be read directly though (the absolute value of the
!   first integer of integer kind ELEM_PACK_SIZE in the data packet container). Additionally,
!   a tagged (marked) super-packet can be created (marked data packet container), in which each
!   simple packet is prefixed with an additional integer element (tag). Tagged data containers
!   are distinguished from untagged ones by the sign of the first integer in the data container,
!   the one whose absolute value shows the number of simple packets stored in the data container.
!   For optimal performance, allocate just enough memory in the data packet container.
! * Upon a request from the manager, data (e.g., a tensor block) can be detached from
!   the corresponding distributed memory space and subsequently destroyed (if needed).
! * Data communication is accomplished via data transfer requests (DTR) and
!   data transfer completion requests (DTCR), using data descriptors. On each MPI process,
!   all data transfer requests are enumerated sequentially in the order they were issued (starting from 1).
! * All procedures return error codes where special return statuses must be distinguished,
!   for example TRY_LATER or NOT_CLEAN. Normally, error codes are smaller by absolute value integers.
!   Contrary, special return statuses are made closer to the HUGE by their absolute values.
! * Currently, chunks of memory that can be attached to a distributed memory space must be
!   multiples of 4 bytes. This is because the memory chunks are mapped to 32-bit words internally.
!   Also, so far only the REAL(4), REAL(8), COMPLEX(4), and COMPLEX(8) data types are explicitly
!   supported, but this is more like an artificial restriction which can be removed relatively easy.
!TYPICAL USAGE WORKFLOW:
! 0. Initialize the DDSS (over MPI) parallel service via the procedure provided in "service_mpi.F90".
! 1. Create one or more distributed spaces <DistrSpace_t> over specific MPI communicators (collective).
! 2. Each MPI process participating in a distributed space can attach a contiguous array of data to
!    that space and obtain the associated global data descriptor <DataDescr_t>.
! 3. Data descriptors can be packed into plain integer packets and collected in a data packet container <PackCont_t>.
!    A data packet container can be communicated between MPI processes using communication handles <CommHandle_t>.
! 4. Upon receival, data descriptors can be unpacked from the data packet container back into <DataDescr_t> objects.
!    These data descriptors can be used for remotely accessing the corresponding data stored on other MPI processes
!    in an asynchronous one-sided manner. It is the user responsibility to ensure proper synchronization between
!    conflicting data accesses (avoid races)!
! 5. Once no longer needed, the data can be detached from the distributed memory space by the owning MPI process.
! 6. Once a distributed memory space is empty and no longer needed, it can be destoyed (it must be empty!).
! 7. Finalize the DDSS/MPI parallel service via the procedure provided in "service_mpi.F90".
!FOR DEVELOPERS ONLY:
! * ELEM_PACK_SIZE is chosen large enough (typically 8 bytes) to capture all basic Fortran data types:
!   INTEGER(1,2,4,8), REAL(4,8), COMPLEX(8,16) (16 is split in two), LOGICAL, CHARACTER(1).
! * Packing of smaller than ELEM_PACK_SIZE data types must be done at the starting address
!   of the ELEM_PACK_SIZE integer word to avoid problems with endianess. In particular, the
!   length of a simple data packet is encoded as INT_COUNT stored in an INTEGER(ELEM_PACK_SIZE)
!   and the number of packets in the data packet container (super-packet) is encoded as INT_MPI
!   stored as INTEGER(ELEM_PACK_SIZE).
! * Currently the entire data packet container volume is communicated. A better choice would be
!   to commmunicate only the used part of it, thus preventing unnecessary communication
!   (in case the capacity of the data container is significant but mostly unused).
       module distributed
!       use, intrinsic:: ISO_C_BINDING
        use service_mpi !includes ISO_C_BINDING & MPI & dil_basic
        use pack_prim
        use timers
!       Depends on stsubs.F90, extern_names.F90 in some procedures
        implicit none
        private
!EXPOSE some <service_mpi>:
        public MPI_COMM_NULL,MPI_REQUEST_NULL !MPI null objects
        public INT_MPI,INT_ADDR,INT_OFFSET,INT_COUNT !MPI integer kinds
        public jo                  !process log output device
        public impis               !size of the global MPI communicator
        public impir               !global MPI rank of the process
        public max_threads         !max number of threads per process
        public gpu_start,gpu_count !GPU range assigned to the process
!PARAMETERS:
 !Output:
        integer, private:: CONS_OUT=6     !default output for this module
        logical, private:: VERBOSE=.TRUE. !verbosity for errors
        integer, private:: DEBUG=0        !debugging mode
        integer, private:: LOGGING=0      !logging mode
 !MPI errors:
        logical, private:: DDSS_MPI_ERR_FATAL=.TRUE. !if .TRUE., MPI errors will cause termination, otherwise just error code
 !Packing/unpacking:
        integer(INT_MPI), parameter, public:: ELEM_PACK_SIZE=max(8,max(C_SIZE_T,max(INT_ADDR,INT_COUNT))) !packing size for integers/logicals/reals/C_pointers/C_sizes
 !Data alignment:
        integer(INT_MPI), parameter, public:: DDSS_DATA_TYPE_ALIGN=4 !data type alignment in bytes
        integer(INT_MPI), parameter, public:: DDSS_BUFFER_ALIGN=max(4,DDSS_DATA_TYPE_ALIGN) !alignment for attached memory buffers in bytes
 !Distributed memory spaces:
        integer(INT_MPI), parameter, public:: DISTR_SPACE_NAME_LEN=128 !max length of a distributed space name (multiple of 8)
 !Data transfers:
  !Active (rank,window) management:
        integer(INT_MPI), parameter, private:: HASH_MOD=256  !number of hash bins for (rank,window) searches
        integer(INT_MPI), parameter, private:: READ_SIGN=+1  !incoming traffic sign (reading direction)
        integer(INT_MPI), parameter, private:: WRITE_SIGN=-1 !outgoing traffic sign (writing direction)
  !Messaging:
        logical, parameter, private:: LAZY_LOCKING=.TRUE.                 !lazy MPI window locking
        logical, parameter, private:: TEST_AND_FLUSH=.FALSE.              !MPI_Test() will call MPI_Win_flush() on completion when entry reference count becomes 0 (not necessary)
        logical, parameter, private:: NO_FLUSH_AFTER_READ_EPOCH=.TRUE.    !if TRUE, there will be no MPI_Win_flush() after the read epoch, thus mandating external synchronization
        logical, parameter, private:: NO_FLUSH_AFTER_WRITE_EPOCH=.FALSE.  !if TRUE, there will be no MPI_Win_flush() after the write epoch, thus mandating external synchronization
        integer(INT_COUNT), parameter, private:: MAX_MPI_MSG_VOL=2**27    !max number of elements in a single MPI message (larger to be split)
        integer(INT_MPI), parameter, private:: MAX_ONESIDED_REQS=4096     !max number of outstanding one-sided data transfer requests per process
        integer(INT_MPI), parameter, public:: DEFAULT_MPI_TAG=0           !default communication tag (for P2P MPI communications)
        integer(INT_MPI), parameter, private:: MPI_ASSER=MPI_MODE_NOCHECK !MPI assertion for locking
       !integer(INT_MPI), parameter, private:: MPI_ASSER=0                !MPI assertion for locking
  !Lock types:
        integer(INT_MPI), parameter, public:: NO_LOCK=0        !no MPI lock (must be zero)
        integer(INT_MPI), parameter, public:: SHARED_LOCK=1    !shared MPI lock (must be positive)
        integer(INT_MPI), parameter, public:: EXCLUSIVE_LOCK=2 !exclusive MPI lock (must be positive)
  !Asynchronous data transfer requests:
        integer(INT_MPI), parameter, public:: MPI_ASYNC_NOT=0  !blocking data transfer request (default)
        integer(INT_MPI), parameter, public:: MPI_ASYNC_NRM=1  !non-blocking data transfer request without a request handle
        integer(INT_MPI), parameter, public:: MPI_ASYNC_REQ=2  !non-blocking data transfer request with a request handle
  !Data transfer communication status:
        integer(INT_MPI), parameter, public:: DDSS_COMM_NONE=0   !no outstanding communication
        integer(INT_MPI), parameter, public:: DDSS_COMM_READ=+1  !outstanding read communication
        integer(INT_MPI), parameter, public:: DDSS_COMM_WRITE=-1 !outstanding write communication
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
          procedure, private:: init=>RankWinInit        !initialize/clean a (rank,window) descriptor
          procedure, private:: print_it=>RankWinPrintIt !prints
        end type RankWin_t
  !Rank/window linked list (active one-sided comms at origin):
        type, private:: RankWinList_t
         integer(8), private:: TransCount=0                         !total number of posted data transfer requests
         real(8), private:: TransSize=0d0                           !total size of all posted data transfers in bytes
         integer(INT_MPI), private:: NumEntries=0                   !number of active entries in the list
         integer(INT_MPI), private:: FirstFree=-1                   !first free (inactive) entry: Must be set to -1 when not initalized
         type(RankWin_t), private:: RankWins(1:MAX_ONESIDED_REQS)   !(rank,window) entries
         integer(INT_MPI), private:: NextEntry(1:MAX_ONESIDED_REQS) !next entry within a bin (linked list)
         integer(INT_MPI), private:: PrevEntry(1:MAX_ONESIDED_REQS) !previous entry within a bin (linked list)
         integer(INT_MPI), private:: HashBin(0:HASH_MOD-1)=0        !first element in each hash bin
         contains
          procedure, private:: init=>RankWinListInit             !clean the (rank,window) list (initialization)
          procedure, private:: test=>RankWinListTest             !test whether a given (rank,window) entry is in the list (with an optional append)
          procedure, private:: new_transfer=>RankWinListNewTrans !register a new data transfer (increment the global transfer ID)
          procedure, private:: delete=>RankWinListDelete         !delete a given active (rank,window) entry
          procedure, private:: delete_all=>RankWinListDeleteAll  !delete all (rank,window) entries
          procedure, private:: flush_all=>RankWinListFlushAll    !flushes all active (rank,window) entries
          procedure, private:: print_all=>RankWinListPrintAll    !print all active communications
        end type RankWinList_t
 !Basic MPI window info:
        type, private:: WinMPI_t
         integer(INT_MPI), private:: Window            !MPI window handle
         integer(INT_MPI), private:: DispUnit=0        !MPI window displacement unit size in bytes
         integer(INT_MPI), private:: CommMPI           !MPI communicator the MPI window is associated with
         logical, private:: Dynamic                    !.TRUE. if the MPI window is dynamic, .FALSE. otherwise
         contains
          procedure, private:: clean=>WinMPIClean                 !clean MPI window info
          procedure, private:: WinMPIPackNew                      !packs the object into a packet (obj_pack_t)
          procedure, private:: WinMPIPackInt                      !packs the object into a plain integer packet
          procedure, private:: WinMPIPack                         !packs the object into SimplePack_t
          generic, private:: pack=>WinMPIPackNew,WinMPIPackInt,WinMPIPack
          procedure, private:: WinMPIUnpackNew                    !unpacks the object from a packet (obj_pack_t)
          procedure, private:: WinMPIUnpackInt                    !unpacks the object from a plain integer packet
          procedure, private:: WinMPIUnpack                       !unpacks the object from SimplePack_t
          generic, private:: unpack=>WinMPIUnpackNew,WinMPIUnpackInt,WinMPIUnpack
          procedure, private:: print_it=>WinMPIPrint              !prints the object data
        end type WinMPI_t
        integer(INT_MPI), parameter, private:: WinMPI_PACK_LEN=4  !packed length of WinMPI_t (in packing integers)
        type(WinMPI_t), parameter, public:: win_mpi_rnd_=WinMPI_t(1983,8,1979,.TRUE.) !random WinMPI_t object for internal testing only
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
         type(C_PTR), private:: LocPtr=C_NULL_PTR !local C pointer to the original data buffer
         integer(INT_MPI), private:: RankMPI=-1   !MPI rank on which the data resides
         type(WinMPI_t), private:: WinMPI         !info on the MPI window the data is exposed with
         integer(INT_ADDR), private:: Offset      !offset in the MPI window (in displacement units)
         integer(INT_COUNT), private:: DataVol    !data volume (number of typed elements)
         integer(INT_MPI), private:: DataType     !data type of each element: {R4,R8,C4,C8,...}, see <dil_basic.F90>
         integer(8), private:: TransID=0_8        !absolute value = data transfer request ID (or zero if none); sign = data transfer direction {READ_SIGN,WRITE_SIGN}
         integer(INT_MPI), private:: StatMPI=MPI_STAT_NONE !status of the data transfer request (see MPI_STAT_XXX parameters above)
         integer(INT_MPI), private:: ReqHandle=MPI_REQUEST_NULL !MPI request handle (for MPI communications with a request handle)
         real(8), private:: TimeStarted=-1d0      !time stamp of communication initiation
         real(8), private:: TimeSynced=-1d0       !time stamp of communication synchronization
         type(object_lock_t), private:: ObjLock   !object lock
         contains
          procedure, private:: clean=>DataDescrClean            !clean a data descriptor
          procedure, private:: init=>DataDescrInit              !set up a data descriptor (initialization)
          procedure, public:: clone=>DataDescrClone             !clones a data descriptor
          procedure, public:: is_set=>DataDescrIsSet            !returns TRUE if the data descriptor is set, FALSE otherwise
          procedure, public:: data_type=>DataDescrDataType      !returns the type of the data associated with the data descriptor
          procedure, public:: data_volume=>DataDescrDataVol     !returns the data volume associated with the data descriptor
          procedure, public:: data_size=>DataDescrDataSize      !returns the data size in bytes
          procedure, public:: get_data_ptr=>DataDescrGetDataPtr !returns a C pointer to the local data buffer
          procedure, public:: get_comm_stat=>DataDescrGetCommStat !returns the current communication status of the data descriptor
          procedure, public:: flush_data=>DataDescrFlushData    !completes an outstanding data transfer request
          procedure, public:: sync_data=>DataDescrSyncData      !synchronizes the private and public data views (in case the data was modified locally)
          procedure, public:: test_data=>DataDescrTestData      !tests whether the data has been transferred to/from the origin (request based)
          procedure, public:: wait_data=>DataDescrWaitData      !waits until the data has been transferred to/from the origin (request based)
          procedure, public:: get_data=>DataDescrGetData        !loads data referred to by a data descriptor into a local buffer
          procedure, public:: acc_data=>DataDescrAccData        !accumulates data from a local buffer to the location specified by a data descriptor
          procedure, private:: DataDescrPackNew                 !packs the DataDescr_t object into a packet (obj_pack_t)
          procedure, private:: DataDescrPackInt                 !packs the DataDescr_t object into a plain integer packet (simple packet)
          procedure, private:: DataDescrPack                    !packs the DataDescr_t object into SimplePack_t (simple packet)
          generic, public:: pack=>DataDescrPackNew,DataDescrPackInt,DataDescrPack
          procedure, private:: DataDescrUnpackNew               !unpacks the DataDescr_t object from a packet (obj_pack_t)
          procedure, private:: DataDescrUnpackInt               !unpacks the DataDescr_t object from a plain integer packet (simple packet)
          procedure, private:: DataDescrUnpack                  !unpacks the DataDescr_t object from SimplePack_t (simple packet)
          generic, public:: unpack=>DataDescrUnpackNew,DataDescrUnpackInt,DataDescrUnpack
          procedure, public:: print_it=>DataDescrPrint          !prints the data descriptor
          procedure, private:: lock=>DataDescrLock              !locks the data descriptor
          procedure, private:: unlock=>DataDescrUnlock          !unlocks the data descriptor
          procedure, public:: clear_lock=>DataDescrClearLock    !clears the lock after cloning DataDescr_t
          final:: DataDescrDtor                                 !dtor
        end type DataDescr_t
        integer(INT_MPI), parameter, private:: DataDescr_PACK_LEN=6+WinMPI_PACK_LEN !packed length of DataDescr_t (in packing integers)
        !type(DataDescr_t), protected:: data_descr_rnd_=&        !random DataDescr_t object for internal testing only
            !&DataDescr_t(C_NULL_PTR,13,win_mpi_rnd_,1024_INT_ADDR,256_INT_COUNT,R8,0_8,MPI_STAT_NONE,MPI_REQUEST_NULL,-1d0,-1d0,&
            !&object_lock_null)
 !Simple packet (plain integer array):
        type, public:: SimplePack_t
         integer(ELEM_PACK_SIZE), pointer, contiguous, private:: Packet(:)=>NULL() !plain integer packet (1d array)
         logical, private:: Alloc                                                  !if .TRUE., the packet was allocated (not associated)
         contains
          procedure, public:: reserve_mem=>SimplePackReserve  !reserves a memory buffer %Packet(:)
          procedure, public:: buf_len=>SimplePackBufLen       !returns the volume of the memory buffer
          procedure, public:: pack_len=>SimplePackFullLen     !returns the full length of the packet in the memory buffer
          procedure, public:: body_len=>SimplePackBodyLen     !returns the body length of the packet in the memory buffer
          procedure, public:: discard=>SimplePackDiscard      !discards an existing packet from the memory buffer
          procedure, public:: destroy=>SimplePackDestroy      !destroys the memory buffer
        end type SimplePack_t
 !Communication handle:
        type, public:: CommHandle_t
         integer(INT_MPI), private:: ReqHandle=MPI_REQUEST_NULL !current MPI request handle
         integer(INT_MPI), private:: LastReq=MPI_REQUEST_NULL   !the most recently completed MPI request handle
         integer(INT_MPI), private:: MPIRank=-1                 !MPI process rank
         integer(INT_MPI), private:: CommMPI                    !MPI communicator
         integer(INT_MPI), private:: CommTag                    !MPI P2P communication tag
         integer(INT_MPI), private:: MPIStat(MPI_STATUS_SIZE)   !MPI status
         class(PackCont_t), pointer, private:: DataCont=>NULL() !pointer to the corresponding data container CommHandle refers to
         real(8), private:: TimeInitiated                       !time stamp when the communication was initiated
         real(8), private:: TimeFoundCompleted                  !time stamp when the communication is confirmed completed
         contains
          procedure, public:: clean=>CommHandleClean !clean the commmunication handle
          procedure, public:: wait=>CommHandleWait   !wait upon completion of the communication associated with the handle
          procedure, public:: test=>CommHandleTest_  !test the completion of the communication associated with the handle
        end type CommHandle_t
 !Data packet container (collection of plain data packets):
        type, public:: PackCont_t
         integer(ELEM_PACK_SIZE), pointer, contiguous, private:: Packets(:)=>NULL() !packet container (0:max): packets(0) = NumPackets
         integer(INT_COUNT), private:: ffe=-1      !first free element in <packets(0:)> = active length of the container
         logical, private:: Alloc=.FALSE.          !TRUE: packets() array was allocated, FALSE: packets(:) array was associated to an external buffer
         integer(INT_MPI), private:: NumPackets=-1 !number of packets in the packet container
         logical, private:: Marked=.FALSE.         !tells whether individual packets are marked (tagged) or not
         contains
          procedure, public:: clean=>PackContClean          !reset the data container to an empty state
          procedure, public:: reserve_mem=>PackContReserve  !reserve memory for the packet container (either allocate or external)
          procedure, public:: num_packets=>PackContNumPacks !returns the total number of packets in the packet container
          procedure, public:: append=>PackContAppend        !add a data packet to the packet container
          procedure, public:: remove=>PackContRemove        !remove a data packet from the packet container
          procedure, public:: send=>PackContSend            !send the packet container to other MPI process(es)
          procedure, public:: receive=>PackContRecv         !receive a packet container from other MPI process
          procedure, private:: register_arrived=>PackContRegArrived !register an arrived data packet container
        end type PackCont_t
!GLOBAL DATA:
 !MPI one-sided data transfer bookkeeping (master thread only):
        type(RankWinList_t), target, private:: RankWinRefs !container for active one-sided communications initiated at the local origin
 !MPI one-sided data transfer statistics:
        real(8), private:: comm_bytes_in=0d0  !amount of data (bytes) one-sided communicated in by the process
        real(8), private:: comm_bytes_out=0d0 !amount of data (bytes) one-sided communicated out by the process
        real(8), private:: comm_time_in=0d0   !time (sec) spent in incoming one-sided communications
        real(8), private:: comm_time_out=0d0  !time (sec) spent in outgoing one-sided communications
!FUNCTION VISIBILITY:
 !Global:
        public data_type_size
        public ddss_flush_all
        public ddss_update_stat
        public ddss_print_stat
 !Auxiliary:
        private get_mpi_int_datatype
        public packet_full_len
        public num_packs_in_container
 !RankWin_t:
        private RankWinInit
        private RankWinPrintIt
 !RankWinList_t:
        private RankWinListInit
        private RankWinListTest
        private RankWinListNewTrans
        private RankWinListDelete
        private RankWinListDeleteAll
        private RankWinListFlushAll
        private RankWinListPrintAll
 !WinMPI_t:
        private WinMPIClean
        private WinMPIPackNew
        private WinMPIPackInt
        private WinMPIPack
        private WinMPIUnpackNew
        private WinMPIUnpackInt
        private WinMPIUnpack
        private WinMPIPrint
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
        private DataDescrClone
        private DataDescrIsSet
        private DataDescrDataType
        private DataDescrDataVol
        private DataDescrDataSize
        private DataDescrGetDataPtr
        private DataDescrGetCommStat
        private DataDescrFlushData
        private DataDescrSyncData
        private DataDescrTestData
        private DataDescrWaitData
        private DataDescrGetData
        private DataDescrAccData
        private DataDescrPackNew
        private DataDescrPackInt
        private DataDescrPack
        private DataDescrUnpackNew
        private DataDescrUnpackInt
        private DataDescrUnpack
        private DataDescrPrint
        private DataDescrLock
        private DataDescrUnlock
        private DataDescrClearLock
        public DataDescrDtor
 !SimplePack_t:
        private SimplePackReserve
        private SimplePackBufLen
        private SimplePackFullLen
        private SimplePackBodyLen
        private SimplePackDiscard
        private SimplePackDestroy
 !CommHandle_t:
        private CommHandleClean
        private CommHandleWait
        private CommHandleTest_
 !PackCont_t:
        private PackContClean
        private PackContReserve
        private PackContNumPacks
        private PackContAppend
        private PackContRemove
        private PackContSend
        private PackContRecv
        private PackContRegArrived

        contains
!METHODS:
!===============================================================
        integer(INT_MPI) function data_type_size(data_type,ierr) !Fortran2008: storage_size()
!Returns the local storage size of a given (pre-registered) data type in bytes.
!All registered data types must possess a size that is multiple of 4 (in bytes)!
        implicit none
        integer(INT_MPI), intent(in):: data_type         !in: pre-registered data type handle: {NO_TYPE,R4,R8,C4,C8,...}
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0; data_type_size=0
        select case(data_type)
        case(R4)
         data_type_size=storage_size(R4_,kind=INT_MPI) !bits
        case(R8)
         data_type_size=storage_size(R8_,kind=INT_MPI) !bits
        case(C4)
         data_type_size=storage_size(C4_,kind=INT_MPI) !bits
        case(C8)
         data_type_size=storage_size(C8_,kind=INT_MPI) !bits
        case(NO_TYPE)
         data_type_size=0
        case default
         if(VERBOSE) write(CONS_OUT,'("#ERROR(distributed:data_type_size): Unknown data type: ",i11)') data_type
         errc=1; data_type_size=-1
        end select
        if(errc.eq.0) then
         if(mod(data_type_size,8).eq.0) then
          data_type_size=data_type_size/8 !convert bits to bytes
          if(mod(data_type_size,DDSS_DATA_TYPE_ALIGN).ne.0) errc=2 !data type size must be a multiple of DDSS_DATA_TYPE_ALIGN
         else
          if(VERBOSE) write(CONS_OUT,'("#ERROR(distributed:data_type_size): Fractional type size detected: ",i11)') data_type_size
          errc=3 !fractional data type size is not allowed
         endif
        endif
        if(present(ierr)) ierr=errc
        return
        end function data_type_size
!--------------------------------------
        subroutine ddss_flush_all(ierr)
!Flushes all active cached communication entries.
         implicit none
         integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
         integer(INT_MPI):: errc

         call RankWinRefs%flush_all(errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine ddss_flush_all
!----------------------------------------------
        subroutine ddss_update_stat(descr,ierr)
         implicit none
         class(DataDescr_t), intent(inout):: descr        !in: active data descriptor
         integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
         integer(INT_MPI):: errc,dir,rk,cm,mr
         integer(INT_COUNT):: data_size

         if(descr%is_set(errc,rk,cm)) then
          if(errc.eq.0) then
           call MPI_Comm_rank(cm,mr,errc)
           if(errc.eq.MPI_SUCCESS) then
            if(mr.ne.rk) then !remote communication
             dir=descr%get_comm_stat(errc)
             if(errc.eq.0) then
              data_size=descr%data_size(errc)
              if(errc.eq.0) then
               if(dir.eq.DDSS_COMM_READ) then
                comm_bytes_in=comm_bytes_in+real(data_size,8)
                comm_time_in=comm_time_in+real(descr%TimeSynced-descr%TimeStarted,8)
               elseif(dir.eq.DDSS_COMM_WRITE) then
                comm_bytes_out=comm_bytes_out+real(data_size,8)
                comm_time_out=comm_time_out+real(descr%TimeSynced-descr%TimeStarted,8)
               else
                errc=6
               endif
              else
               errc=5
              endif
             else
              errc=4
             endif
            endif
           else
            errc=3
           endif
          else
           errc=2
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ddss_update_stat
!-----------------------------------------------
        subroutine ddss_print_stat(ierr,dev_out)
!Prints the current DDSS communication statistics.
        implicit none
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI), intent(in), optional:: dev_out !in: output device (defaults to <jo>)
        integer(INT_MPI):: devo,errc

        errc=0
        if(present(dev_out)) then; devo=dev_out; else; devo=jo; endif
        if(comm_time_in.gt.0d0) then
         write(devo,'("#INFO(DDSS): Incoming one-sided communication effective bandwidth (GB/s) = ",F10.3)')&
         &comm_bytes_in/(comm_time_in*(1024d0*1024d0*1024d0))
        endif
        if(comm_time_out.gt.0d0) then
         write(devo,'("#INFO(DDSS): Outgoing one-sided communication effective bandwidth (GB/s) = ",F10.3)')&
         &comm_bytes_out/(comm_time_out*(1024d0*1024d0*1024d0))
        endif
        write(devo,'("#INFO(DDSS): Incoming one-sided communication volume (GB) = ",F12.3)') comm_bytes_in/(1024d0*1024d0*1024d0)
        write(devo,'("#INFO(DDSS): Outgoing one-sided communication volume (GB) = ",F12.3)') comm_bytes_out/(1024d0*1024d0*1024d0)
        flush(devo)
        call RankWinRefs%print_all(errc,devo)
        if(present(ierr)) ierr=errc
        return
        end subroutine ddss_print_stat
!========================================================================
        function get_mpi_int_datatype(int_kind,mpi_data_typ) result(ierr)
!Given an integer kind, returns the corresponing MPI integer data type handle.
!If not found, returns a negatie error code <ierr>.
        implicit none
        integer(INT_MPI):: ierr
        integer(INT_MPI), intent(in):: int_kind
        integer(INT_MPI), intent(out):: mpi_data_typ

        ierr=0; mpi_data_typ=0
        select case(int_kind)
        case(1)
         mpi_data_typ=MPI_INTEGER1
        case(2)
         mpi_data_typ=MPI_INTEGER2
        case(4)
         mpi_data_typ=MPI_INTEGER4
        case(8)
         mpi_data_typ=MPI_INTEGER8
        case(16)
         mpi_data_typ=MPI_INTEGER16
        case default
         ierr=-1
        end select
        return
        end function get_mpi_int_datatype
!-------------------------------------------------------------
        function packet_full_len(packet,body_len) result(plen)
!Returns the full packet length (number of elements);
!each element is an INTEGER(ELEM_PACK_SIZE).
        implicit none
        integer(INT_COUNT):: plen                                 !out: data packet full length
        integer(ELEM_PACK_SIZE), intent(in), target:: packet(0:*) !in: data packet
        integer(INT_COUNT), intent(out), optional:: body_len      !out: data packet body length (useful legnth)
        integer(INT_COUNT), pointer:: len_p
        type(C_PTR):: cptr

        cptr=c_loc(packet(0)); call c_f_pointer(cptr,len_p) !body length has integer kind INT_COUNT
        plen=1+len_p !full length (header + body)
        if(present(body_len)) body_len=max(0_INT_COUNT,len_p) !body length (in packing integers)
        nullify(len_p)
        return
        end function packet_full_len
!-----------------------------------------------------------------------
        function num_packs_in_container(pack_cont,tagged) result(npacks)
!Returns the total number of packets in the data packet container (super-packet).
        implicit none
        integer(INT_MPI):: npacks                            !out: number of data packets in the data packet container
        integer(ELEM_PACK_SIZE), intent(in):: pack_cont(0:*) !in: plain data packet container (integer array)
        logical, intent(out), optional:: tagged              !out: TRUE if the packet container is tagged, FALSE otherwise

        npacks=int(pack_cont(0),INT_MPI) !properly convert the storing integer into the result
        if(npacks.lt.0) then !tagged packet container
         npacks=-npacks; if(present(tagged)) tagged=.TRUE.
        else !tag-free packet container
         if(present(tagged)) tagged=.FALSE.
        endif
        return
        end function num_packs_in_container
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
!---------------------------------------------------
        subroutine RankWinPrintIt(this,ierr,dev_out)
         implicit none
         class(RankWin_t), intent(in):: this     !inout: (rank,window) descriptor
         integer, intent(out), optional:: ierr   !out: error code
         integer, intent(in), optional:: dev_out !in: output device id (6:screen, default)
         integer:: errc,dev

         errc=0
         dev=6; if(present(dev_out)) dev=dev_out
         write(dev,'("WinRank{",i13,1x,i5,1x,"{",i2,"}: Ref = ",i4,"; Sync = ",i11,"}")')&
         &this%Window,this%Rank,this%LockType,this%RefCount,this%LastSync
         flush(dev)
         if(present(ierr)) ierr=errc
         return
        end subroutine RankWinPrintIt
!============================================
        subroutine RankWinListInit(this,ierr)
!Cleans the (rank,window) list (must be called before use).
        implicit none
        class(RankWinList_t), intent(inout):: this       !inout: (rank,window) list
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: i

        this%TransCount=0; this%TransSize=0d0
        this%NumEntries=0; this%FirstFree=1; this%HashBin(:)=0 !setting .FirstFree to 1 means initialized
        this%PrevEntry(1)=0; do i=2,MAX_ONESIDED_REQS; this%PrevEntry(i)=i-1; enddo !linked list
        do i=1,MAX_ONESIDED_REQS-1; this%NextEntry(i)=i+1; enddo; this%NextEntry(MAX_ONESIDED_REQS)=0 !linked list
        do i=1,MAX_ONESIDED_REQS; call this%RankWins(i)%init(); enddo !init all entries to null
        if(DEBUG.ge.2) then
         write(jo,'("#DEBUG(distributed:RankWinList.Init)[",i7,"]: (rank,win)-list initialized.")') impir
         flush(jo)
        endif
        if(present(ierr)) ierr=0
        return
        end subroutine RankWinListInit
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
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success) or TRY_LATER
        logical, intent(in), optional:: append           !in: if .TRUE., a new entry will be appended if not found
        integer(INT_MPI):: i,j,m,n,errc
        logical:: apnd

        errc=0; RankWinListTest=0
        if(present(append)) then; apnd=append; else; apnd=.FALSE.; endif
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
         if(DEBUG.ge.1.and.RankWinListTest.gt.0) then
          write(jo,'("#DEBUG(distributed:RankWinList.Test)[",i5,":",i3,"]: Existing (window,rank) found: ",i5,": ",i13,1x,i7)')&
          &impir,thread_id,RankWinListTest,win,rank
          flush(jo)
         endif
         if(RankWinListTest.le.0.and.apnd) then !append if not found
          if(DEBUG.ge.2) then
           write(jo,'("#DEBUG(distribiuted::RankWinList.Test)[",i7,"]: Registering new (rank,window) entry: "'//&
           &',i7,1x,i13,3x,i6,1x,i6)') impir,rank,win,this%NumEntries,this%FirstFree
           flush(jo)
          endif
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
            if(DEBUG.ge.1) then
             write(jo,'("#DEBUG(distributed:RankWinList.Test)[",i5,":",i3,"]: New (window,rank) registered: ",i5,": ",i13,1x,i7)')&
             &impir,thread_id,RankWinListTest,win,rank
             flush(jo)
            endif
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
!-----------------------------------------------------------
        subroutine RankWinListNewTrans(this,dd,rwe,dir,ierr)
!This subroutine increments the global transfer ID and communicated data size
!and assigns the current transfer ID to the given data descriptor.
        implicit none
        class(RankWinList_t), intent(inout):: this       !inout: (rank,window) list
        class(DataDescr_t), intent(inout):: dd           !inout: valid data descriptor
        integer(INT_MPI), intent(in):: rwe               !in: the corresponding (rank,window) entry
        integer(INT_MPI), intent(in):: dir               !in: data transfer direction: {READ_SIGN,WRITE_SIGN}
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(8), parameter:: int8=0_8
        integer(INT_MPI):: n,errc

        errc=0
        if(this%TransCount.lt.huge(int8)) then
         if(rwe.ge.1.and.rwe.le.MAX_ONESIDED_REQS) then
          if(abs(dir).eq.1) then
           this%TransCount=this%TransCount+1
           n=data_type_size(dd%DataType)
           this%TransSize=this%TransSize+real(n,8)*real(dd%DataVol,8)
           dd%TransID=this%TransCount*int(dir,8) !transfer id = global transaction ID * data transfer sign
           this%RankWins(rwe)%RefCount=this%RankWins(rwe)%RefCount+1
           if(DEBUG.ge.1) then
            write(jo,'("#DEBUG(distributed:RankWinList.NewTrans)[",i5,":",i3,"]: New transfer ",i13,": ")',ADVANCE='NO')&
            &impir,thread_id,dd%TransID
            call this%RankWins(rwe)%print_it(dev_out=jo)
            flush(jo)
           endif
          else
           errc=1
          endif
         else
          errc=2
         endif
        else
         if(VERBOSE) write(CONS_OUT,'("#FATAL(distributed:RankWinList.NewTrans): Max int8 MPI message count exceeded!")')
         errc=3
         call quit(-1,'#FATAL: Unable to continue: Distributed Data Service failed!')
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine RankWinListNewTrans
!--------------------------------------------------------
        subroutine RankWinListDelete(this,entry_num,ierr)
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
          this%PrevEntry(entry_num)=0; this%NextEntry(entry_num)=this%FirstFree
          if(this%FirstFree.gt.0) this%PrevEntry(this%FirstFree)=entry_num
          this%FirstFree=entry_num; this%NumEntries=this%NumEntries-1
          if(DEBUG.ge.1) then
           write(jo,'("#DEBUG(distributed:RankWinList.Delete)[",i5,":",i3,"]: (window,rank) entry ",i5," deleted: ")',ADVANCE='NO')&
           &impir,thread_id,entry_num
           call this%RankWins(entry_num)%print_it(dev_out=jo)
           flush(jo)
          endif
          call this%RankWins(entry_num)%init() !clean entry
         else
          errc=1 !empty entries cannot be deleted
         endif
        else
         errc=2
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine RankWinListDelete
!--------------------------------------------------------
        subroutine RankWinListDeleteAll(this,ierr,window)
!Deletes all (rank,window) entries with a proper synchronization when needed.
!If <window> is present, only the (rank,window) entries corresponding to the
!given <window> will be deleted with a proper synchronization when needed.
        implicit none
        class(RankWinList_t), intent(inout):: this       !inout: (rank,window) list
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI), intent(in), optional:: window  !in: specific window
        integer(INT_MPI):: errc,i,rnk,win

        errc=0
        if(present(window)) then
         do i=lbound(this%RankWins,1),ubound(this%RankWins,1)
          rnk=this%RankWins(i)%Rank; win=this%RankWins(i)%Window
          if(rnk.ge.0) then !active entry
           if(win.eq.window) then
            if(this%RankWins(i)%RefCount.gt.0.or.LAZY_LOCKING) call MPI_Win_unlock(rnk,win,errc)
            if(errc.eq.0) then
             call this%delete(i,errc); if(errc.ne.0) then; errc=1; exit; endif
            else
             errc=2; exit
            endif
           endif
          endif
         enddo
        else
         do i=lbound(this%RankWins,1),ubound(this%RankWins,1)
          rnk=this%RankWins(i)%Rank; win=this%RankWins(i)%Window
          if(rnk.ge.0) then !active entry
           if(this%RankWins(i)%RefCount.gt.0.or.LAZY_LOCKING) call MPI_Win_unlock(rnk,win,errc)
           if(errc.eq.0) then
            call this%delete(i,errc); if(errc.ne.0) then; errc=3; exit; endif
           else
            errc=4; exit
           endif
          endif
         enddo
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine RankWinListDeleteAll
!------------------------------------------------
        subroutine RankWinListFlushAll(this,ierr)
!Flushes all active (rank,window) entries.
        implicit none
        class(RankWinList_t), intent(inout):: this       !inout: (rank,window) list
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc,i,rnk,win

        errc=0
        do i=lbound(this%RankWins,1),ubound(this%RankWins,1)
         rnk=this%RankWins(i)%Rank; win=this%RankWins(i)%Window
         if(rnk.ge.0) then !active entry
          if(this%RankWins(i)%RefCount.gt.0.or.LAZY_LOCKING) call MPI_Win_flush(rnk,win,errc)
          if(errc.eq.0) then
           call this%delete(i,errc); if(errc.ne.0) then; errc=1; exit; endif
          else
           errc=2; exit
          endif
         endif
        enddo
        if(present(ierr)) ierr=errc
        return
        end subroutine RankWinListFlushAll
!--------------------------------------------------------
        subroutine RankWinListPrintAll(this,ierr,dev_out)
!This subroutine prints the current state of the RankWinList_t,
!that is, all active communications initiated at origin.
        implicit none
        class(RankWinList_t), intent(in):: this          !in: RankWinList object
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI), intent(in), optional:: dev_out !in: output device (defaults to <jo>)
        integer(INT_MPI):: i,devo,errc

        errc=0
        if(present(dev_out)) then; devo=dev_out; else; devo=jo; endif
        write(devo,'("#Printing the current Rank-Window list (active communications):")')
        write(devo,'(1x,"Entry",4x,"Rank",6x,"Window",4x,"Lock",1x,"Refs",3x,"Last Synced")')
        do i=1,this%NumEntries
         write(devo,'(1x,i4,3x,i7,1x,i13,2x,i2,2x,i4,1x,i13)') i,this%RankWins(i)%Rank,this%RankWins(i)%Window,&
         &this%RankWins(i)%LockType,this%RankWins(i)%RefCount,this%RankWins(i)%LastSync
        enddo
        write(devo,'(1x,"Total number of transactions registered  = ",i18)') this%TransCount
        write(devo,'(1x,"Total size of communicated data (Mbytes) = ",D18.6)') this%TransSize/(1024d0*1024d0)
        if(present(ierr)) ierr=errc
        return
        end subroutine RankWinListPrintAll
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
!------------------------------------------------
       subroutine WinMPIPackNew(this,packet,ierr)
!Packs the WinMPI_t object into a packet.
        implicit none
        class(WinMPI_t), intent(in):: this             !in: WinMPI_t object
        class(obj_pack_t), intent(inout):: packet      !inout: packet
        integer(INT_MPI), intent(out), optional:: ierr !out: error code
        integer(INT_MPI):: errc

        call pack_builtin(packet,this%Window,errc)
        if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%DispUnit,errc)
        if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%CommMPI,errc)
        if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%Dynamic,errc)
        if(present(ierr)) ierr=errc
        return
       end subroutine WinMPIPackNew
!----------------------------------------------------------
        subroutine WinMPIPackInt(this,packet,ierr,pack_len)
!Packs a WinMPI_t object into a plain integer packet <packet>.
!The first integer is always the useful length of the packet, that is,
!the number of the following integer elements storing the information.
!It is the user responsibility to provide a large enough packet buffer.
        implicit none
        class(WinMPI_t), intent(in):: this                           !in: WinMPI_t object
        integer(ELEM_PACK_SIZE), intent(inout), target:: packet(0:*) !out: large enough data packet (length + information)
        integer(INT_MPI), intent(inout), optional:: ierr             !out: error code (0:success)
        integer(INT_COUNT), intent(out), optional:: pack_len         !out: full data packet length (in packing integers)
        integer(INT_MPI):: errc
        integer(INT_COUNT):: pl
        type(C_PTR):: cptr
        integer(INT_COUNT), pointer:: len_p
        integer(INT_MPI), pointer:: impi_p
        logical, pointer:: log_p

        errc=0; packet(0)=0; pl=0
        pl=pl+1; packet(pl)=0; cptr=c_loc(packet(pl)); call c_f_pointer(cptr,impi_p); impi_p=this%Window
        pl=pl+1; packet(pl)=0; cptr=c_loc(packet(pl)); call c_f_pointer(cptr,impi_p); impi_p=this%DispUnit
        pl=pl+1; packet(pl)=0; cptr=c_loc(packet(pl)); call c_f_pointer(cptr,impi_p); impi_p=this%CommMPI
        pl=pl+1; packet(pl)=0; cptr=c_loc(packet(pl)); call c_f_pointer(cptr,log_p); log_p=this%Dynamic
        cptr=c_loc(packet(0)); call c_f_pointer(cptr,len_p); len_p=pl !packet body length
        len_p=>NULL(); log_p=>NULL(); impi_p=>NULL()
        if(present(pack_len)) pack_len=1+pl !header integer + packet body
        if(present(ierr)) ierr=errc
        return
        end subroutine WinMPIPackInt
!-------------------------------------------------------
        subroutine WinMPIPack(this,packet,ierr,pack_len)
!Packs a WinMPI_t object into a simple packet <packet>.
!The first integer is always the useful length of the packet, that is,
!the number of the following integer elements storing the information.
!It is the user responsibility to provide a large enough packet buffer.
        implicit none
        class(WinMPI_t), intent(in):: this                   !in: WinMPI_t object
        class(SimplePack_t), intent(inout), target:: packet  !out: data packet (length + information)
        integer(INT_MPI), intent(inout), optional:: ierr     !out: error code (0:success)
        integer(INT_COUNT), intent(out), optional:: pack_len !out: full data packet length (in packing integers)
        integer(INT_MPI):: errc

        errc=0
        if(packet%buf_len().lt.1+WinMPI_PACK_LEN) call packet%reserve_mem(int(1+WinMPI_PACK_LEN,INT_MPI),errc)
        if(errc.eq.0) then
         if(present(pack_len)) then
          call this%WinMPIPackInt(packet%Packet,errc,pack_len)
         else
          call this%WinMPIPackInt(packet%Packet,errc)
         endif
        else
         if(present(pack_len)) pack_len=0
         errc=1
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine WinMPIPack
!---------------------------------------------------
        subroutine WinMPIUnpackNew(this,packet,ierr)
!Unpacks a WinMPI_t object from a packet.
         implicit none
         class(WinMPI_t), intent(out):: this            !out: WinMPI_t object
         class(obj_pack_t), intent(inout):: packet      !inout: packet
         integer(INT_MPI), intent(out), optional:: ierr !out: error code
         integer(INT_MPI):: errc

         call unpack_builtin(packet,this%Window,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%DispUnit,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%CommMPI,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%Dynamic,errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine WinMPIUnpackNew
!------------------------------------------------------------
        subroutine WinMPIUnpackInt(this,packet,ierr,pack_len)
!Unpacks a WinMPI_t object from a plain integer packet <packet>.
        implicit none
        class(WinMPI_t), intent(inout):: this                     !out: unpacked WinMPI_t object
        integer(ELEM_PACK_SIZE), intent(in), target:: packet(0:*) !in: plain integer data packet (length + information)
        integer(INT_MPI), intent(inout), optional:: ierr          !out: error code (0:success)
        integer(INT_COUNT), intent(out), optional:: pack_len      !out: full data packet length (in packing integers)
        integer(INT_MPI):: errc
        integer(INT_COUNT):: pl
        type(C_PTR):: cptr
        integer(INT_COUNT), pointer:: len_p
        integer(INT_MPI), pointer:: impi_p
        logical, pointer:: log_p

        errc=0
!Check the length:
        cptr=c_loc(packet(0)); call c_f_pointer(cptr,len_p)
        if(present(pack_len)) pack_len=1+len_p
        if(len_p.eq.WinMPI_PACK_LEN) then
!Unpack:
         call this%clean(); pl=0
         pl=pl+1; cptr=c_loc(packet(pl)); call c_f_pointer(cptr,impi_p); this%Window=impi_p
         pl=pl+1; cptr=c_loc(packet(pl)); call c_f_pointer(cptr,impi_p); this%DispUnit=impi_p
         pl=pl+1; cptr=c_loc(packet(pl)); call c_f_pointer(cptr,impi_p); this%CommMPI=impi_p
         pl=pl+1; cptr=c_loc(packet(pl)); call c_f_pointer(cptr,log_p); this%Dynamic=log_p
         len_p=>NULL(); log_p=>NULL(); impi_p=>NULL()
        else
         errc=1
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine WinMPIUnpackInt
!---------------------------------------------------------
        subroutine WinMPIUnpack(this,packet,ierr,pack_len)
!Unpacks a WinMPI_t object from a simple packet <packet>.
        implicit none
        class(WinMPI_t), intent(inout):: this                !out: unpacked WinMPI_t object
        class(SimplePack_t), intent(in), target:: packet     !in: data packet (length + information)
        integer(INT_MPI), intent(inout), optional:: ierr     !out: error code (0:success)
        integer(INT_COUNT), intent(out), optional:: pack_len !out: full data packet length (in packing integers)
        integer(INT_MPI):: errc

        errc=0
        if(packet%body_len().ge.WinMPI_PACK_LEN) then
         if(present(pack_len)) then
          call this%WinMPIUnpackInt(packet%Packet,errc,pack_len)
         else
          call this%WinMPIUnpackInt(packet%Packet,errc)
         endif
        else
         errc=-1
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine WinMPIUnpack
!-------------------------------------------------
        subroutine WinMPIPrint(this,dev_out,space)
!Prints the object data.
        use stsubs, only: numchar
        implicit none
        class(WinMPI_t), intent(in):: this               !in: object to print
        integer(INT_MPI), intent(in), optional:: dev_out !in: output device
        integer(INT_MPI), intent(in), optional:: space   !in: left alignment (formatting)
        character(32):: sfmt
        integer(INT_MPI):: devo,sp
        integer:: fl

        devo=jo; if(present(dev_out)) devo=dev_out
        sp=0; if(present(space)) sp=space
        if(sp.gt.0) then; call numchar(int(sp),fl,sfmt); sfmt(fl+1:fl+2)='x,'; fl=fl+2; else; fl=0; endif
        write(devo,'('//sfmt(1:fl)//'"#Printing WinMPI_t object:")')
        write(devo,'('//sfmt(1:fl)//'"  MPI window       : ",i18)') this%Window
        write(devo,'('//sfmt(1:fl)//'"  Displacement unit: ",i18)') this%DispUnit
        write(devo,'('//sfmt(1:fl)//'"  MPI communicator : ",i18)') this%CommMPI
        write(devo,'('//sfmt(1:fl)//'"  Dynamic window?  : ",l18)') this%Dynamic
        return
        end subroutine WinMPIPrint
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

        call MPI_Barrier(comm_mpi,errc) !test the validity of the MPI communicator
        if(errc.eq.0) then
         if(this%WinSize.lt.0) then !uninitialized window
          call MPI_Win_create_dynamic(MPI_INFO_NULL,comm_mpi,this%WinMPI%Window,errc)
          if(errc.eq.0) then
           call MPI_Win_get_attr(this%WinMPI%Window,MPI_WIN_DISP_UNIT,attr,flag,errc)
           if(errc.eq.0.and.flag) then
            this%WinMPI%DispUnit=int(attr,INT_MPI)
            this%WinMPI%CommMPI=comm_mpi
            this%WinMPI%Dynamic=.TRUE.
            this%WinSize=0
            if(DEBUG.ge.2) then
             write(jo,'("#DEBUG(distributed:DataWin.Create)[",i7,"]: MPI window created: ",i11,1x,i11,1x,i2,1x,i11)')&
             &impir,this%WinMPI%Window,this%WinMPI%CommMPI,this%WinMPI%DispUnit,this%WinSize
             flush(jo)
            endif
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

        call MPI_Barrier(this%WinMPI%CommMPI,errc) !test the validity of the MPI communicator
        if(errc.eq.0) then
         if(this%WinSize.eq.0) then !MPI window must be empty
          if(DEBUG.ge.2) then
           write(jo,'("#DEBUG(distributed:DataWin.Destroy)[",i7,"]: MPI window to destroy: ",i11,1x,i11,1x,i2,1x,i11)')&
           &impir,this%WinMPI%Window,this%WinMPI%CommMPI,this%WinMPI%DispUnit,this%WinSize
           flush(jo)
          endif
          call MPI_Win_free(this%WinMPI%Window,errc)
          if(errc.eq.0) then
           if(DEBUG.ge.2) then
            write(jo,'("#DEBUG(distributed:DataWin.Destroy)[",i7,"]: MPI window destroyed: ",i11,1x,i11,1x,i2,1x,i11)')&
            &impir,this%WinMPI%Window,this%WinMPI%CommMPI,this%WinMPI%DispUnit,this%WinSize
            flush(jo)
           endif
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
!The size of the data buffer in bytes must be a multiple of DDSS_BUFFER_ALIGN.
        implicit none
        class(DataWin_t), intent(inout):: this           !inout: data window
        type(C_PTR), intent(in):: loc_ptr                !in: C pointer to the local data buffer
        integer(INT_ADDR), intent(in):: loc_size         !in: size of the local data buffer in bytes: Must be multiple of DDSS_BUFFER_ALIGN
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc,ier
        integer(INT_ADDR):: vol,disp
        real(4), pointer, contiguous:: r4_ptr(:) !the (arbitrary) data will be mapped as a real(4) array

        errc=0
        if(.not.c_associated(loc_ptr,C_NULL_PTR)) then
         if(loc_size.gt.0.and.mod(loc_size,int(DDSS_BUFFER_ALIGN,INT_ADDR)).eq.0) then
          vol=loc_size/4 !volume in 4-byte (32-bit) words
          if(this%WinSize.ge.0) then !data window has been initialized
           call c_f_pointer(loc_ptr,r4_ptr,(/vol/))
           call attach_buffer(r4_ptr,loc_size,errc)
           if(errc.eq.0) then
            this%WinSize=this%WinSize+loc_size
            if(DEBUG.ge.1) then
             call MPI_Get_Displacement(loc_ptr,disp,ier); if(ier.ne.0) disp=-666
             write(jo,'("#DEBUG(distributed:DataWin.Attach)[",i5,":",i3,"]: Buffer attached: ",i11,1x,i11,1x,i13,": ",i18)')&
             &impir,thread_id,loc_size,this%WinMPI%Window,this%WinSize,disp
             flush(jo)
            endif
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
         call MPI_Win_attach(this%WinMPI%Window,r4_arr,bsize,jerr)
         return
         end subroutine attach_buffer

        end subroutine DataWinAttach
!-----------------------------------------------------------
        subroutine DataWinDetach(this,loc_ptr,loc_size,ierr)
!Detaches a previously attached local (contiguous) data buffer from an MPI data window.
!The size of the data buffer in bytes must be a multiple of DDSS_BUFFER_ALIGN.
        implicit none
        class(DataWin_t), intent(inout):: this           !inout: data window
        type(C_PTR), intent(in):: loc_ptr                !in: C pointer to the local data buffer to be detached
        integer(INT_ADDR), intent(in):: loc_size         !in: size of the local data buffer in bytes: Must be multiple of DDSS_BUFFER_ALIGN
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc
        integer(INT_ADDR):: vol
        real(4), pointer, contiguous:: r4_ptr(:) !the (arbitrary) data is mapped as a real(4) array

        errc=0
        if(.not.c_associated(loc_ptr,C_NULL_PTR)) then
         if(loc_size.gt.0.and.mod(loc_size,int(DDSS_BUFFER_ALIGN,INT_ADDR)).eq.0) then
          vol=loc_size/4 !volume in 4-byte words
          if(this%WinSize.ge.loc_size) then
           call c_f_pointer(loc_ptr,r4_ptr,(/vol/))
           call detach_buffer(r4_ptr,errc)
           if(errc.eq.0) then
            this%WinSize=this%WinSize-loc_size
            if(DEBUG.ge.1) then
             write(jo,'("#DEBUG(distributed:DataWin.Detach)[",i5,":",i3,"]: Buffer detached: ",i11,1x,i11,1x,i13)')&
             &impir,thread_id,loc_size,this%WinMPI%Window,this%WinSize
             flush(jo)
            endif
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
         call MPI_Win_detach(this%WinMPI%Window,r4_arr,jerr)
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
         call MPI_Win_sync(this%WinMPI%Window,errc)
         if(errc.eq.0) then
          if(DEBUG.ge.2) then
           write(jo,'("#DEBUG(distributed:DataWin.Sync)[",i7,"]: MPI window synced: ",i11,1x,i11)')&
           &impir,this%WinMPI%Window,this%WinSize
           flush(jo)
          endif
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

        call MPI_Barrier(comm_mpi,errc) !test the validity of the MPI communicator
        if(errc.eq.0) then
         if(RankWinRefs%FirstFree.lt.0) call RankWinRefs%init() !init the (rank,win) list
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
             if(DEBUG.ge.2) then
              write(jo,'("#DEBUG(distributed:DistrSpace.Create)[",i7,"]: Distributed space created: ",i11,1x,i4,1x,A32)')&
              &impir,this%CommMPI,this%NumWins,this%SpaceName(1:min(min(len_trim(this%SpaceName),DISTR_SPACE_NAME_LEN),32))
              flush(jo)
             endif
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
         call MPI_Barrier(comm_mpi,errc)
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

        call RankWinRefs%delete_all(errc)
        if(errc.eq.0) then
         call MPI_Barrier(this%CommMPI,errc) !test the validity of the MPI communicator
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
               if(DEBUG.ge.2) then
                write(jo,'("#DEBUG(distributed:DistrSpace.Destroy)[",i7,"]: Distributed space destroyed: ",i11,1x,i4,1x,A32)')&
                &impir,this%CommMPI,this%NumWins,this%SpaceName(1:min(min(len_trim(this%SpaceName),DISTR_SPACE_NAME_LEN),32))
                flush(jo)
               endif
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
        else
         errc=7
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
        integer(INT_MPI), intent(in):: data_type         !in: pre-existing data type: {R4,R8,C4,C8,...}
        integer(INT_COUNT), intent(in):: data_vol        !in: positive data volume (number of typed elements)
        class(DataDescr_t), intent(out):: data_descr     !out: filled data descriptor
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: i,m,my_rank,errc
        integer(INT_ADDR):: min_mem,loc_size
!$OMP FLUSH
        errc=0
        if(this%NumWins.gt.0) then !initialized distributed memory space
         m=1; min_mem=this%DataWins(m)%WinSize
         do i=2,this%NumWins !select the least occupied MPI data window
          if(this%DataWins(i)%WinSize.lt.min_mem) then
           min_mem=this%DataWins(i)%WinSize; m=i
          endif
         enddo
         call RankWinRefs%delete_all(errc,this%DataWins(m)%WinMPI%Window) !sync and delete all (rank,win)-entries on this window
         if(errc.eq.0) then
          loc_size=data_vol*data_type_size(data_type,errc) !data buffer size in bytes
          if(errc.eq.0.and.loc_size.gt.0) then
           call this%DataWins(m)%attach(loc_ptr,loc_size,errc)
           if(errc.eq.0) then
            call MPI_Comm_rank(this%CommMPI,my_rank,errc)
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
        else
         errc=6
        endif
        if(present(ierr)) ierr=errc
!$OMP FLUSH
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
!$OMP FLUSH
        errc=0
        if(this%NumWins.gt.0) then !initialized distributed memory space
         call MPI_Comm_rank(this%CommMPI,my_rank,errc)
         if(errc.eq.0) then
          if(data_descr%RankMPI.eq.my_rank) then
           loc_size=data_descr%DataVol*data_type_size(data_descr%DataType,errc) !data buffer size in bytes
           if(errc.eq.0.and.loc_size.gt.0) then
            m=0
            do i=1,this%NumWins !find the MPI window
             if(this%DataWins(i)%WinMPI%Window.eq.data_descr%WinMPI%Window) then; m=i; exit; endif
            enddo
            if(m.gt.0) then
             call RankWinRefs%delete_all(errc,this%DataWins(m)%WinMPI%Window) !sync and delete all (rank,win)-entries on this window
             if(errc.eq.0) then
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
        else
         errc=7
        endif
        if(present(ierr)) ierr=errc
!$OMP FLUSH
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
        call this%lock()
        this%RankMPI=-1
        call this%WinMPI%clean(errc); if(errc.ne.0) errc=1
        call this%unlock()
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
        integer(INT_MPI), intent(in):: data_type         !in: data type: {R4,R8,C4,C8,...}
        integer(INT_COUNT), intent(in):: data_vol        !in: positive data volume (number of elements)
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc,elem_size

        errc=0
        if(process_rank.ge.0) then
         if(.not.c_associated(loc_ptr,C_NULL_PTR)) then
          if(data_vol.gt.0) then
           elem_size=data_type_size(data_type,errc)
           if(errc.eq.0.and.elem_size.gt.0) then
            call this%lock()
            this%LocPtr=loc_ptr
            this%RankMPI=process_rank
            this%WinMPI=win_mpi
            this%DataVol=data_vol
            this%DataType=data_type
            this%TransID=0 !0 means "never assigned a Transfer ID"
            this%StatMPI=MPI_STAT_NONE
            this%ReqHandle=MPI_REQUEST_NULL
            this%TimeStarted=-1d0
            this%TimeSynced=-1d0
            call MPI_Get_Displacement(loc_ptr,this%Offset,errc); if(errc.ne.0) errc=1
            call this%unlock()
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
!---------------------------------------------------
        subroutine DataDescrClone(this,another,ierr)
        implicit none
        class(DataDescr_t), intent(inout):: this               !in: source data descriptor
        class(DataDescr_t), allocatable, intent(out):: another !out: cloned data descriptor
        integer(INTD), intent(out), optional:: ierr            !out: error code
        integer(INTD):: errc

        call this%lock()
        allocate(another,SOURCE=this,STAT=errc)
        call another%clear_lock()
        call this%unlock()
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrClone
!------------------------------------------------------------------------
        function DataDescrIsSet(this,ierr,proc_rank,mpi_comm) result(ans)
!Returns TRUE if the data descriptor is set, FALSE otherwise.
         implicit none
         logical:: ans                                       !out: answer
         class(DataDescr_t), intent(inout):: this            !in: data descriptor
         integer(INT_MPI), intent(out), optional:: ierr      !out: error code (0:success)
         integer(INT_MPI), intent(out), optional:: proc_rank !out: MPI process rank
         integer(INT_MPI), intent(out), optional:: mpi_comm  !out: MPI communicator

         call this%lock()
         ans=(this%RankMPI.ge.0)
         if(ans.and.present(proc_rank)) proc_rank=this%RankMPI
         if(ans.and.present(mpi_comm)) mpi_comm=this%WinMPI%CommMPI
         call this%unlock()
         if(present(ierr)) ierr=0
         return
        end function DataDescrIsSet
!-------------------------------------------------------------
        function DataDescrDataType(this,ierr) result(data_typ)
!Returns the type of the data associated with the data descriptor.
        implicit none
        integer(INT_MPI):: data_typ                      !out: data type
        class(DataDescr_t), intent(inout):: this         !in: data descriptor
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0
        call this%lock()
        if(this%RankMPI.ge.0) then
         data_typ=this%DataType
        else
         data_typ=NO_TYPE; errc=-1
        endif
        call this%unlock()
        if(present(ierr)) ierr=errc
        return
        end function DataDescrDataType
!------------------------------------------------------------
        function DataDescrDataVol(this,ierr) result(data_vol)
!Returns the volume of the data associated with the data descriptor.
        implicit none
        integer(INT_COUNT):: data_vol                    !out: data volume (number of typed elements)
        class(DataDescr_t), intent(inout):: this         !in: data descriptor
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0
        call this%lock()
        if(this%RankMPI.ge.0) then
         data_vol=this%DataVol
        else
         data_vol=-1; errc=-1
        endif
        call this%unlock()
        if(present(ierr)) ierr=errc
        return
        end function DataDescrDataVol
!--------------------------------------------------------------
        function DataDescrDataSize(this,ierr) result(data_size)
!Returns the size (in bytes) of the data associated with the data descriptor.
        implicit none
        integer(INT_COUNT):: data_size                   !out: data size in bytes
        class(DataDescr_t), intent(inout):: this         !in: data descriptor
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0
        call this%lock()
        if(this%RankMPI.ge.0) then
         data_size=this%DataVol*int(data_type_size(this%DataType,errc),INT_COUNT)
         if(errc.ne.0) data_size=-1
        else
         data_size=-1; errc=-1
        endif
        call this%unlock()
        if(present(ierr)) ierr=errc
        return
        end function DataDescrDataSize
!---------------------------------------------------------------
        function DataDescrGetDataPtr(this,ierr) result(data_ptr)
!Returns a C pointer to the local data buffer. If the data descriptor
!corresponds to remote data, the result is undefined.
        implicit none
        type(C_PTR):: data_ptr                         !out: data pointer
        class(DataDescr_t), intent(inout):: this       !in: data descriptor
        integer(INT_MPI), intent(out), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0
        call this%lock()
        if(this%RankMPI.ge.0) then
         data_ptr=this%LocPtr
        else
         data_ptr=C_NULL_PTR; errc=-1
        endif
        call this%unlock()
        if(present(ierr)) ierr=errc
        return
        end function DataDescrGetDataPtr
!---------------------------------------------------------------
        function DataDescrGetCommStat(this,ierr,req) result(sts)
!Returns the current communication status of the data descriptor.
        implicit none
        integer(INT_MPI):: sts                         !out: communication status: {DDSS_COMM_NONE,DDSS_COMM_READ,DDSS_COMM_WRITE}
        class(DataDescr_t), intent(inout):: this       !in: data descriptor
        integer(INT_MPI), intent(out), optional:: ierr !out: error code (0:success)
        integer(INT_MPI), intent(out), optional:: req  !out: MPI communication request handle, if set
        integer(INT_MPI):: errc,creq

        errc=0; sts=DDSS_COMM_NONE; creq=MPI_REQUEST_NULL
        call this%lock()
        if(this%StatMPI.eq.MPI_STAT_PROGRESS_NRM.or.this%StatMPI.eq.MPI_STAT_PROGRESS_REQ) then
         if(this%TransID.gt.0) then
          sts=DDSS_COMM_READ
         elseif(this%TransID.lt.0) then
          sts=DDSS_COMM_WRITE
         else
          errc=-2
         endif
         if(this%StatMPI.eq.MPI_STAT_PROGRESS_REQ) creq=this%ReqHandle
        elseif(this%StatMPI.eq.MPI_STAT_ONESIDED_ERR) then
         errc=-1
        endif
        call this%unlock()
        if(present(req)) req=creq
        if(present(ierr)) ierr=errc
        return
        end function DataDescrGetCommStat
!-----------------------------------------------------
        subroutine DataDescrFlushData(this,ierr,local)
!Completes an MPI one-sided communication specified by a given data descriptor.
!Other outstanding MPI communications on the same (rank,window) will be completed as well.
!It is safe to flush the same (valid) data descriptor multiple times. However,
!if the data descriptor was flushed the first time at the origin only (local=.TRUE.),
!one will no longer be able to request a full completion (both at origin and target).
!In case of a one-sided synchronization error, the data request will still be marked
!as pending, thus blocking the corresponding (rank,window) entry from release!
        implicit none
        class(DataDescr_t), intent(inout):: this         !inout: data descriptor
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        logical, intent(in), optional:: local            !in: if .TRUE., the data is flushed only at the origin
        integer(INT_MPI):: rwe,errc
        type(RankWin_t), pointer:: rw_entry
        logical:: lcl,synced
        real(8):: tm

        errc=0; synced=.FALSE.
        if(present(local)) then; lcl=local; else; lcl=.FALSE.; endif !default is flushing both at the origin and target
        call this%lock()
        if(this%RankMPI.ge.0) then
         if(this%StatMPI.eq.MPI_STAT_PROGRESS_NRM.or.this%StatMPI.eq.MPI_STAT_PROGRESS_REQ) then !first flush
          rwe=RankWinRefs%test(this%RankMPI,this%WinMPI%Window,errc) !get the (rank,window) entry
          if(rwe.gt.0.and.errc.eq.0) then
           rw_entry=>RankWinRefs%RankWins(rwe)
           if(DEBUG.ge.1) then
            write(jo,'("#DEBUG(distributed:DataDescr.FlushData)[",i5,":",i3,"]: Syncing transfer ",i13,": ")',ADVANCE='NO')&
            &impir,thread_id,this%TransID
            call rw_entry%print_it(dev_out=jo)
            flush(jo)
           endif
           if(rw_entry%RefCount.gt.1) then !not the last reference to this (rank,window)
            if(lcl) then !complete at origin only
             if(abs(this%TransID).gt.rw_entry%LastSync) then
              if(DEBUG.ge.1) then
               write(jo,'("#DEBUG(distributed:DataDescr.FlushData)[",i5,":",i3,"]: WIN_FLUSH_LOCAL: ")',ADVANCE='NO')&
               &impir,thread_id
               call rw_entry%print_it(dev_out=jo)
               flush(jo)
              endif
              call nvtx_push('MPI_Win_flush_local'//CHAR_NULL,1)
              tm=time_sys_sec()
              call MPI_Win_flush_local(rw_entry%Rank,rw_entry%Window,errc) !complete at the origin only
              tm=time_sys_sec()-tm
              call nvtx_pop()
              if(LOGGING.gt.0) then
               write(jo,'("#MSG(DDSS::flush)[",i5,":",i3,"]: MPI_WIN_FLUSH_LOCAL(",i13,",",i5,") time (sec) = ",F8.4)')&
               &impir,thread_id,rw_entry%Window,rw_entry%Rank,tm
               flush(jo)
              endif
              if(errc.ne.0.and.DDSS_MPI_ERR_FATAL) then
               call quit(errc,'#FATAL(distributed:DataDescr.FlushData): MPI_Win_flush_local failed!')
              else
               synced=(errc.eq.0)
              endif
             endif
             if(errc.eq.0) then
              this%TimeSynced=time_sys_sec()
              call ddss_update_stat(this)
              this%StatMPI=MPI_STAT_COMPLETED_ORIG
              rw_entry%RefCount=rw_entry%RefCount-1
             else
              this%StatMPI=MPI_STAT_ONESIDED_ERR; errc=1
             endif
            else !complete at both origin and target
             if(abs(this%TransID).gt.rw_entry%LastSync) then
              if(DEBUG.ge.1) then
               write(jo,'("#DEBUG(distributed:DataDescr.FlushData)[",i5,":",i3,"]: WIN_FLUSH: ")',ADVANCE='NO') impir,thread_id
               call rw_entry%print_it(dev_out=jo)
               flush(jo)
              endif
              call nvtx_push('MPI_Win_flush'//CHAR_NULL,2)
              tm=time_sys_sec()
              call MPI_Win_flush(rw_entry%Rank,rw_entry%Window,errc) !complete both at the origin and target
              tm=time_sys_sec()-tm
              call nvtx_pop()
              if(LOGGING.gt.0) then
               write(jo,'("#MSG(DDSS::flush)[",i5,":",i3,"]: MPI_WIN_FLUSH(",i13,",",i5,") time (sec) = ",F8.4)')&
               &impir,thread_id,rw_entry%Window,rw_entry%Rank,tm
               flush(jo)
              endif
              if(errc.ne.0.and.DDSS_MPI_ERR_FATAL) then
               call quit(errc,'#FATAL(distributed:DataDescr.FlushData): MPI_Win_flush failed!')
              else
               synced=(errc.eq.0)
              endif
             endif
             if(errc.eq.0) then
              this%TimeSynced=time_sys_sec()
              call ddss_update_stat(this)
              this%StatMPI=MPI_STAT_COMPLETED
              rw_entry%RefCount=rw_entry%RefCount-1
             else
              this%StatMPI=MPI_STAT_ONESIDED_ERR; errc=2
             endif
            endif
           else !the last reference to this (rank,window): Mandatory UNLOCK
            if(DEBUG.ge.1) then
             write(jo,'("#DEBUG(distributed:DataDescr.FlushData)[",i5,":",i3,"]: WIN_UNLOCK(.flush): ")',ADVANCE='NO')&
             &impir,thread_id
             call rw_entry%print_it(dev_out=jo)
             flush(jo)
            endif
            if(LAZY_LOCKING) then
             call nvtx_push('MPI_Win_flush'//CHAR_NULL,2)
             tm=time_sys_sec()
             call MPI_Win_flush(rw_entry%Rank,rw_entry%Window,errc) !complete both at origin and target
             tm=time_sys_sec()-tm
             call nvtx_pop()
             if(LOGGING.gt.0) then
              write(jo,'("#MSG(DDSS::flush)[",i5,":",i3,"]: MPI_WIN_FLUSH(",i13,",",i5,") time (sec) = ",F8.4)')&
              &impir,thread_id,rw_entry%Window,rw_entry%Rank,tm
              flush(jo)
             endif
            else
             call nvtx_push('MPI_Win_unlock'//CHAR_NULL,3)
             tm=time_sys_sec()
             call MPI_Win_unlock(rw_entry%Rank,rw_entry%Window,errc) !complete both at origin and target
             tm=time_sys_sec()-tm
             call nvtx_pop()
             if(LOGGING.gt.0) then
              write(jo,'("#MSG(DDSS::flush)[",i5,":",i3,"]: MPI_WIN_UNLOCK(",i13,",",i5,") time (sec) = ",F8.4)')&
              &impir,thread_id,rw_entry%Window,rw_entry%Rank,tm
              flush(jo)
             endif
            endif
            if(errc.ne.0.and.DDSS_MPI_ERR_FATAL) then
             call quit(errc,'#FATAL(distributed:DataDescr.FlushData): MPI_Win_unlock/MPI_Win_flush failed!')
            else
             synced=(errc.eq.0)
            endif
            if(errc.eq.0) then
             this%TimeSynced=time_sys_sec()
             call ddss_update_stat(this)
             this%StatMPI=MPI_STAT_COMPLETED
             rw_entry%RefCount=rw_entry%RefCount-1
            else
             this%StatMPI=MPI_STAT_ONESIDED_ERR; errc=3
            endif
           endif
           if(errc.eq.0) then
            if(synced) rw_entry%LastSync=RankWinRefs%TransCount !update the last sync event for this (rank,win)
            if(rw_entry%RefCount.eq.0) then !delete the (rank,window) entry if no references are attached to it
             nullify(rw_entry)
             if(.not.LAZY_LOCKING) then
              call RankWinRefs%delete(rwe,errc); if(errc.ne.0) errc=4
             endif
            elseif(rw_entry%RefCount.lt.0) then
             if(VERBOSE) write(CONS_OUT,'("#FATAL(distributed:DataDescr.FlushData): Negative reference count: ",i12)')&
                         &rw_entry%RefCount
             nullify(rw_entry)
             errc=5
             call quit(-1,'#FATAL: Unable to continue: Distributed Data Service failed!')
            endif
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
        call this%unlock()
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrFlushData
!----------------------------------------------
        subroutine DataDescrSyncData(this,ierr)
!Synchronizes the private and public views of the data in case the data was modified locally.
        implicit none
        class(DataDescr_t), intent(inout):: this         !inout: data descriptor
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0
        call this%lock()
        if(this%RankMPI.ge.0) then
         if(DEBUG.ge.2) then
          write(jo,'("#DEBUG[",i7,"]: WIN_SYNC: ",i11)') impir,this%WinMPI%Window !debug
          flush(jo)
         endif
         call MPI_Win_sync(this%WinMPI%Window,errc)
         if(errc.ne.0) errc=1
        else
         errc=-1
        endif
        call this%unlock()
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrSyncData
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
        real(8):: tm

        errc=0; DataDescrTestData=.FALSE.
        call this%lock()
        if(this%RankMPI.ge.0) then
         if(this%StatMPI.eq.MPI_STAT_PROGRESS_REQ) then
          rwe=RankWinRefs%test(this%RankMPI,this%WinMPI%Window,errc) !get the (rank,window) entry
          if(rwe.gt.0.and.errc.eq.0) then
           rw_entry=>RankWinRefs%RankWins(rwe)
           if(rw_entry%RefCount.gt.0) then
            if(abs(this%TransID).gt.rw_entry%LastSync) then
             if(DEBUG.ge.1) then
              write(jo,'("#DEBUG(distributed:DataDescr.TestData)[",i5,":",i3,"]: MPI_TEST: ")',ADVANCE='NO') impir,thread_id
              call rw_entry%print_it(dev_out=jo)
              flush(jo)
             endif
             call MPI_Test(this%ReqHandle,compl,mpi_stat,errc)
             if(errc.eq.0) then
              if(compl) then
               this%TimeSynced=time_sys_sec()
               call ddss_update_stat(this)
               this%StatMPI=MPI_STAT_COMPLETED_ORIG
               rw_entry%RefCount=rw_entry%RefCount-1
               DataDescrTestData=.TRUE.
               if(LOGGING.gt.0) then
                write(jo,'("#MSG(DDSS::test)[",i5,":",i3,"]: MPI_TEST() completed with time (sec) = ",F8.4)')&
                &impir,thread_id,this%TimeSynced-this%TimeStarted
                flush(jo)
               endif
              endif
             else
              this%StatMPI=MPI_STAT_ONESIDED_ERR; errc=1
             endif
            else
             this%TimeSynced=time_sys_sec()
             call ddss_update_stat(this)
             this%StatMPI=MPI_STAT_COMPLETED_ORIG
             rw_entry%RefCount=rw_entry%RefCount-1
             DataDescrTestData=.TRUE.
             if(LOGGING.gt.0) then
              write(jo,'("#MSG(DDSS::test)[",i5,":",i3,"]: MPI_TEST() overcompleted with time (sec) = ",F8.4)')&
              &impir,thread_id,this%TimeSynced-this%TimeStarted
              flush(jo)
             endif
            endif
            if(errc.eq.0) then
             if(rw_entry%RefCount.eq.0) then !delete the (rank,window) entry if no references are attached to it
              if(LAZY_LOCKING) then
               if(TEST_AND_FLUSH) then
                if(DEBUG.ge.1) then
                 write(jo,'("#DEBUG(distributed:DataDescr.TestData)[",i5,":",i3,"]: WIN_FLUSH(.test): ")',ADVANCE='NO')&
                 &impir,thread_id
                 call rw_entry%print_it(dev_out=jo)
                 flush(jo)
                endif
                call nvtx_push('MPI_Win_flush'//CHAR_NULL,2)
                tm=time_sys_sec()
                call MPI_Win_flush(rw_entry%Rank,rw_entry%Window,errc)
                tm=time_sys_sec()-tm
                call nvtx_pop()
                if(LOGGING.gt.0) then
                 write(jo,'("#MSG(DDSS::test)[",i5,":",i3,"]: MPI_WIN_FLUSH(",i13,",",i5,") time (sec) = ",F8.4)')&
                 &impir,thread_id,rw_entry%Window,rw_entry%Rank,tm
                 flush(jo)
                endif
                if(errc.eq.0) then
                 this%StatMPI=MPI_STAT_COMPLETED
                 rw_entry%LastSync=RankWinRefs%TransCount
                else
                 if(DDSS_MPI_ERR_FATAL) call quit(errc,'#FATAL(distributed:DataDescr.TestData): MPI_Win_flush failed!')
                 errc=2
                endif
               endif
              else
               if(DEBUG.ge.1) then
                write(jo,'("#DEBUG(distributed:DataDescr.TestData)[",i5,":",i3,"]: WIN_UNLOCK(.test): ")',ADVANCE='NO')&
                &impir,thread_id
                call rw_entry%print_it(dev_out=jo)
                flush(jo)
               endif
               call nvtx_push('MPI_Win_unlock'//CHAR_NULL,3)
               tm=time_sys_sec()
               call MPI_Win_unlock(rw_entry%Rank,rw_entry%Window,errc)
               tm=time_sys_sec()-tm
               call nvtx_pop()
               if(LOGGING.gt.0) then
                write(jo,'("#MSG(DDSS::test)[",i5,":",i3,"]: MPI_WIN_UNLOCK(",i13,",",i5,") time (sec) = ",F8.4)')&
                &impir,thread_id,rw_entry%Window,rw_entry%Rank,tm
                flush(jo)
               endif
               if(errc.eq.0) then
                this%StatMPI=MPI_STAT_COMPLETED
                nullify(rw_entry)
                call RankWinRefs%delete(rwe,errc); if(errc.ne.0) errc=3
               else
                if(DDSS_MPI_ERR_FATAL) call quit(errc,'#FATAL(distributed:DataDescr.TestData): MPI_Win_unlock failed!')
                errc=4
               endif
              endif
             endif
            endif
           else
            if(VERBOSE) write(CONS_OUT,'("#FATAL(distributed:DataDescr.TestData): Invalid reference count: ",i12)')&
                        &rw_entry%RefCount
            nullify(rw_entry)
            errc=5
            call quit(-1,'#FATAL: Unable to continue: Distributed Data Service failed!')
           endif
           nullify(rw_entry)
          else
           errc=6
          endif
         elseif(this%StatMPI.eq.MPI_STAT_COMPLETED.or.this%StatMPI.eq.MPI_STAT_COMPLETED_ORIG) then
          DataDescrTestData=.TRUE.
         else
          errc=7
         endif
        else
         errc=8
        endif
        call this%unlock()
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
        real(8):: tm

        errc=0
        call this%lock()
        if(this%RankMPI.ge.0) then
         if(this%StatMPI.eq.MPI_STAT_PROGRESS_REQ) then
          rwe=RankWinRefs%test(this%RankMPI,this%WinMPI%Window,errc) !get the (rank,window) entry
          if(rwe.gt.0.and.errc.eq.0) then
           rw_entry=>RankWinRefs%RankWins(rwe)
           if(rw_entry%RefCount.gt.0) then
            if(abs(this%TransID).gt.rw_entry%LastSync) then
             if(DEBUG.ge.1) then
              write(jo,'("#DEBUG(distributed:DataDescr.WaitData)[",i5,":",i3,"]: MPI_WAIT: ")',ADVANCE='NO') impir,thread_id
              call rw_entry%print_it(dev_out=jo)
              flush(jo)
             endif
             call MPI_Wait(this%ReqHandle,mpi_stat,errc)
             if(errc.eq.0) then
              this%TimeSynced=time_sys_sec()
              call ddss_update_stat(this)
              this%StatMPI=MPI_STAT_COMPLETED_ORIG
              rw_entry%RefCount=rw_entry%RefCount-1
              if(LOGGING.gt.0) then
               write(jo,'("#MSG(DDSS::test)[",i5,":",i3,"]: MPI_WAIT() completed with time (sec) = ",F8.4)')&
               &impir,thread_id,this%TimeSynced-this%TimeStarted
               flush(jo)
              endif
             else
              this%StatMPI=MPI_STAT_ONESIDED_ERR; errc=1
             endif
            else
             this%TimeSynced=time_sys_sec()
             call ddss_update_stat(this)
             this%StatMPI=MPI_STAT_COMPLETED_ORIG
             rw_entry%RefCount=rw_entry%RefCount-1
             if(LOGGING.gt.0) then
              write(jo,'("#MSG(DDSS::test)[",i5,":",i3,"]: MPI_WAIT() overcompleted with time (sec) = ",F8.4)')&
              &impir,thread_id,this%TimeSynced-this%TimeStarted
              flush(jo)
             endif
            endif
            if(errc.eq.0) then
             if(rw_entry%RefCount.eq.0) then !delete the (rank,window) entry if no references are attached to it
              if(LAZY_LOCKING) then
               if(TEST_AND_FLUSH) then
                if(DEBUG.ge.1) then
                 write(jo,'("#DEBUG(distributed:DataDescr.WaitData)[",i5,":",i3,"]: WIN_FLUSH(.wait): ")',ADVANCE='NO')&
                 &impir,thread_id
                 call rw_entry%print_it(dev_out=jo)
                 flush(jo)
                endif
                call nvtx_push('MPI_Win_flush'//CHAR_NULL,2)
                tm=time_sys_sec()
                call MPI_Win_flush(rw_entry%Rank,rw_entry%Window,errc)
                tm=time_sys_sec()-tm
                call nvtx_pop()
                if(LOGGING.gt.0) then
                 write(jo,'("#MSG(DDSS::wait)[",i5,":",i3,"]: MPI_WIN_FLUSH(",i13,",",i5,") time (sec) = ",F8.4)')&
                 &impir,thread_id,rw_entry%Window,rw_entry%Rank,tm
                 flush(jo)
                endif
                if(errc.eq.0) then
                 this%StatMPI=MPI_STAT_COMPLETED
                 rw_entry%LastSync=RankWinRefs%TransCount
                else
                 if(DDSS_MPI_ERR_FATAL) call quit(errc,'#FATAL(distributed:DataDescr.WaitData): MPI_Win_flush failed!')
                 errc=2
                endif
               endif
              else
               if(DEBUG.ge.1) then
                write(jo,'("#DEBUG(distributed:DataDescr.WaitData)[",i5,":",i3,"]: WIN_UNLOCK(.wait): ")',ADVANCE='NO')&
                &impir,thread_id
                call rw_entry%print_it(dev_out=jo)
                flush(jo)
               endif
               call nvtx_push('MPI_Win_unlock'//CHAR_NULL,3)
               tm=time_sys_sec()
               call MPI_Win_unlock(rw_entry%Rank,rw_entry%Window,errc)
               tm=time_sys_sec()-tm
               call nvtx_pop()
               if(LOGGING.gt.0) then
                write(jo,'("#MSG(DDSS::wait)[",i5,":",i3,"]: MPI_WIN_UNLOCK(",i13,",",i5,") time (sec) = ",F8.4)')&
                &impir,thread_id,rw_entry%Window,rw_entry%Rank,tm
                flush(jo)
               endif
               if(errc.eq.0) then
                this%StatMPI=MPI_STAT_COMPLETED
                nullify(rw_entry)
                call RankWinRefs%delete(rwe,errc); if(errc.ne.0) errc=3
               else
                if(DDSS_MPI_ERR_FATAL) call quit(errc,'#FATAL(distributed:DataDescr.WaitData): MPI_Win_unlock failed!')
                errc=4
               endif
              endif
             endif
            endif
           else
            if(VERBOSE) write(CONS_OUT,'("#FATAL(distributed:DataDescr.WaitData): Invalid reference count: ",i12)')&
                        &rw_entry%RefCount
            errc=5
            nullify(rw_entry)
            call quit(-1,'#FATAL: Unable to continue: Distributed Data Service failed!')
           endif
           nullify(rw_entry)
          else
           errc=6
          endif
         elseif(this%StatMPI.ne.MPI_STAT_COMPLETED.and.this%StatMPI.ne.MPI_STAT_COMPLETED_ORIG) then
          errc=7
         endif
        else
         errc=8
        endif
        call this%unlock()
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
        integer(INT_MPI):: rwe,errc,asnc
        real(4), pointer, contiguous:: r4_ptr(:)
        real(8), pointer, contiguous:: r8_ptr(:)
        complex(4), pointer, contiguous:: c4_ptr(:)
        complex(8), pointer, contiguous:: c8_ptr(:)

        errc=0
        if(present(async)) then; asnc=async; else; asnc=MPI_ASYNC_NOT; endif !default is synchronous communication
        call this%lock()
        if(asnc.eq.MPI_ASYNC_NOT.or.asnc.eq.MPI_ASYNC_NRM.or.asnc.eq.MPI_ASYNC_REQ) then
         if(.not.c_associated(loc_ptr,C_NULL_PTR)) then
          if(this%RankMPI.ge.0) then
           if(this%StatMPI.eq.MPI_STAT_NONE.or.this%StatMPI.eq.MPI_STAT_COMPLETED.or.&
             &this%StatMPI.eq.MPI_STAT_COMPLETED_ORIG) then
            rwe=RankWinRefs%test(this%RankMPI,this%WinMPI%Window,errc,append=.TRUE.) !get the (rank,window) entry
            if(errc.eq.0) then
             if(this%DataVol.gt.0) then
              if(.not.(asnc.eq.MPI_ASYNC_REQ.and.this%DataVol.gt.huge(asnc))) then
               call modify_lock(rwe,errc) !modify the (rank,window) lock status if needed
               if(errc.eq.0) then
                if(DEBUG.ge.1) then
                 write(jo,'("#DEBUG(distributed:DataDescr.GetData)[",i5,":",i3,"]: Get: ",i18,1x,i9,": ")',ADVANCE='NO')&
                 &impir,thread_id,this%Offset,this%DataVol
                 call RankWinRefs%RankWins(rwe)%print_it(dev_out=jo)
                 flush(jo)
                endif
                select case(this%DataType)
                case(R4)
                 call c_f_pointer(loc_ptr,r4_ptr,(/this%DataVol/))
                 call start_get_r4(r4_ptr,errc); if(errc.ne.0) errc=1
                case(R8)
                 call c_f_pointer(loc_ptr,r8_ptr,(/this%DataVol/))
                 call start_get_r8(r8_ptr,errc); if(errc.ne.0) errc=2
                case(C4)
                 call c_f_pointer(loc_ptr,c4_ptr,(/this%DataVol/))
                 call start_get_c4(c4_ptr,errc); if(errc.ne.0) errc=3
                case(C8)
                 call c_f_pointer(loc_ptr,c8_ptr,(/this%DataVol/))
                 call start_get_c8(c8_ptr,errc); if(errc.ne.0) errc=4
                case(NO_TYPE)
                 errc=5
                case default
                 errc=6
                end select
                if(errc.eq.0) then
                 call RankWinRefs%new_transfer(this,rwe,READ_SIGN,errc) !register a new transfer (will also set this%TransID field)
                 if(errc.eq.0) then
                  if(asnc.eq.MPI_ASYNC_NOT) then
                   call this%flush_data(errc); if(errc.ne.0) errc=7
                  endif
                 else
                  errc=8
                 endif
                endif
               else
                errc=9
               endif
              else
               errc=10
              endif
             elseif(this%DataVol.eq.0) then
              this%StatMPI=MPI_STAT_COMPLETED
             else
              errc=11
             endif
            else
             if(errc.ne.TRY_LATER) errc=12
            endif
           else
            errc=13
           endif
          else
           errc=14
          endif
         else
          errc=15
         endif
        else
         errc=16
        endif
        call this%unlock()
        if(present(ierr)) ierr=errc
        return

        contains

         subroutine modify_lock(rw,jerr)
         integer(INT_MPI), intent(in):: rw
         integer(INT_MPI), intent(out):: jerr
         type(RankWin_t), pointer:: rw_entry
         real(8):: tm

         jerr=0; rw_entry=>RankWinRefs%RankWins(rw)
         if(rw_entry%LockType*READ_SIGN.lt.0) then !communication direction change
          if(LAZY_LOCKING) then
           if(DEBUG.ge.1) then
            write(jo,'("#DEBUG(distributed:DataDescr.GetData)[",i5,":",i3,"]: WIN_FLUSH(.get): ")',ADVANCE='NO') impir,thread_id
            call rw_entry%print_it(dev_out=jo)
            flush(jo)
           endif
           call nvtx_push('MPI_Win_flush'//CHAR_NULL,2)
           tm=time_sys_sec()
           if(.not.NO_FLUSH_AFTER_WRITE_EPOCH) call MPI_Win_flush(rw_entry%Rank,rw_entry%Window,jerr)
           tm=time_sys_sec()-tm
           call nvtx_pop()
           if(LOGGING.gt.0) then
            write(jo,'("#MSG(DDSS::get)[",i5,":",i3,"]: MPI_WIN_FLUSH(",i13,",",i5,") time (sec) = ",F8.4)')&
            &impir,thread_id,rw_entry%Window,rw_entry%Rank,tm
            flush(jo)
           endif
          else
           if(DEBUG.ge.1) then
            write(jo,'("#DEBUG(distributed:DataDescr.GetData)[",i5,":",i3,"]: WIN_UNLOCK(.get): ")',ADVANCE='NO') impir,thread_id
            call rw_entry%print_it(dev_out=jo)
            flush(jo)
           endif
           call nvtx_push('MPI_Win_unlock'//CHAR_NULL,3)
           tm=time_sys_sec()
           call MPI_Win_unlock(rw_entry%Rank,rw_entry%Window,jerr)
           tm=time_sys_sec()-tm
           call nvtx_pop()
           if(LOGGING.gt.0) then
            write(jo,'("#MSG(DDSS::get)[",i5,":",i3,"]: MPI_WIN_UNLOCK(",i13,",",i5,") time (sec) = ",F8.4)')&
            &impir,thread_id,rw_entry%Window,rw_entry%Rank,tm
            flush(jo)
           endif
          endif
          if(jerr.eq.0) then
           if(LAZY_LOCKING) then
            rw_entry%LockType=SHARED_LOCK*READ_SIGN
           else
            rw_entry%LockType=NO_LOCK
           endif
           rw_entry%LastSync=RankWinRefs%TransCount
          else
           if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.GetData): MPI_Win_unlock/MPI_Win_flush failed!')
           jerr=1
          endif
         endif
         if(jerr.eq.0.and.rw_entry%LockType.eq.NO_LOCK) then
          if(DEBUG.ge.1) then
           write(jo,'("#DEBUG(distributed:DataDescr.GetData)[",i5,":",i3,"]: WIN_LOCK(.get): ")',ADVANCE='NO') impir,thread_id
           call rw_entry%print_it(dev_out=jo)
           flush(jo)
          endif
          call nvtx_push('MPI_Win_lock'//CHAR_NULL,0)
          tm=time_sys_sec()
          call MPI_Win_lock(MPI_LOCK_SHARED,rw_entry%Rank,MPI_ASSER,rw_entry%Window,jerr)
          tm=time_sys_sec()-tm
          call nvtx_pop()
          if(LOGGING.gt.0) then
           write(jo,'("#MSG(DDSS::get)[",i5,":",i3,"]: MPI_WIN_LOCK(",i13,",",i5,") time (sec) = ",F8.4)')&
           &impir,thread_id,rw_entry%Window,rw_entry%Rank,tm
           flush(jo)
          endif
          if(jerr.eq.0) then
           rw_entry%LockType=SHARED_LOCK*READ_SIGN
          else
           if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.GetData): MPI_Win_lock failed!')
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
         integer(INT_MPI):: jv,jdts,jdu
         integer(INT_ADDR):: jtarg
         real(8):: tm

         jerr=0; jdts=data_type_size(R4)
         jdu=this%WinMPI%DispUnit
         ji=MAX_MPI_MSG_VOL
         if(asnc.ne.MPI_ASYNC_REQ) then !regular
          do js=1,this%DataVol,ji
           jv=int(min(this%DataVol-js+1,ji),INT_MPI)
           jtarg=this%Offset+((js-1)*jdts)/jdu
           call nvtx_push('MPI_Get'//CHAR_NULL,4)
           tm=time_sys_sec()
           call MPI_Get(r4_arr(js:js+jv-1),jv,MPI_REAL4,this%RankMPI,jtarg,jv,MPI_REAL4,this%WinMPI%Window,jerr)
           tm=time_sys_sec()-tm
           call nvtx_pop()
           if(LOGGING.gt.0) then
            write(jo,'("#MSG(DDSS::get)[",i5,":",i3,"]: MPI_GET(",i13,",",i5,") time (sec) = ",F8.4)')&
            &impir,thread_id,this%WinMPI%Window,this%RankMPI,tm
            flush(jo)
           endif
           if(jerr.ne.0) then
            if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.GetData): MPI_Get failed!')
            exit
           endif
          enddo
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_NRM
           this%TimeStarted=time_sys_sec()
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=1
          endif
         else !request-handle
          jv=this%DataVol
          call nvtx_push('MPI_Rget'//CHAR_NULL,5)
          tm=time_sys_sec()
          call MPI_Rget(r4_arr,jv,MPI_REAL4,this%RankMPI,this%Offset,jv,MPI_REAL4,this%WinMPI%Window,this%ReqHandle,jerr)
          tm=time_sys_sec()-tm
          call nvtx_pop()
          if(LOGGING.gt.0) then
           write(jo,'("#MSG(DDSS:get)[",i5,":",i3,"]: MPI_RGET(",i13,",",i5,") time (sec) = ",F8.4)')&
           &impir,thread_id,this%WinMPI%Window,this%RankMPI,tm
           flush(jo)
          endif
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_REQ
           this%TimeStarted=time_sys_sec()
          else
           if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.GetData): MPI_Rget failed!')
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=2
          endif
         endif
         return
         end subroutine start_get_r4

         subroutine start_get_r8(r8_arr,jerr)
         real(8), intent(inout):: r8_arr(*) !asynchronous
         integer(INT_MPI), intent(out):: jerr
         integer(INT_COUNT):: ji,js
         integer(INT_MPI):: jv,jdts,jdu
         integer(INT_ADDR):: jtarg
         real(8):: tm

         jerr=0; jdts=data_type_size(R8)
         jdu=this%WinMPI%DispUnit
         ji=MAX_MPI_MSG_VOL
         if(asnc.ne.MPI_ASYNC_REQ) then !regular
          do js=1,this%DataVol,ji
           jv=int(min(this%DataVol-js+1,ji),INT_MPI)
           jtarg=this%Offset+((js-1)*jdts)/jdu
           call nvtx_push('MPI_Get'//CHAR_NULL,4)
           tm=time_sys_sec()
           call MPI_Get(r8_arr(js:js+jv-1),jv,MPI_REAL8,this%RankMPI,jtarg,jv,MPI_REAL8,this%WinMPI%Window,jerr)
           tm=time_sys_sec()-tm
           call nvtx_pop()
           if(LOGGING.gt.0) then
            write(jo,'("#MSG(DDSS::get)[",i5,":",i3,"]: MPI_GET(",i13,",",i5,") time (sec) = ",F8.4)')&
            &impir,thread_id,this%WinMPI%Window,this%RankMPI,tm
            flush(jo)
           endif
           if(jerr.ne.0) then
            if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.GetData): MPI_Get failed!')
            exit
           endif
          enddo
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_NRM
           this%TimeStarted=time_sys_sec()
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=1
          endif
         else !request-handle
          jv=this%DataVol
          call nvtx_push('MPI_Rget'//CHAR_NULL,5)
          tm=time_sys_sec()
          call MPI_Rget(r8_arr,jv,MPI_REAL8,this%RankMPI,this%Offset,jv,MPI_REAL8,this%WinMPI%Window,this%ReqHandle,jerr)
          tm=time_sys_sec()-tm
          call nvtx_pop()
          if(LOGGING.gt.0) then
           write(jo,'("#MSG(DDSS::get)[",i5,":",i3,"]: MPI_RGET(",i13,",",i5,") time (sec) = ",F8.4)')&
           &impir,thread_id,this%WinMPI%Window,this%RankMPI,tm
           flush(jo)
          endif
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_REQ
           this%TimeStarted=time_sys_sec()
          else
           if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.GetData): MPI_Rget failed!')
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=2
          endif
         endif
         return
         end subroutine start_get_r8

         subroutine start_get_c4(c4_arr,jerr)
         complex(4), intent(inout):: c4_arr(*) !asynchronous
         integer(INT_MPI), intent(out):: jerr
         integer(INT_COUNT):: ji,js
         integer(INT_MPI):: jv,jdts,jdu
         integer(INT_ADDR):: jtarg
         real(8):: tm

         jerr=0; jdts=data_type_size(C4)
         jdu=this%WinMPI%DispUnit
         ji=MAX_MPI_MSG_VOL
         if(asnc.ne.MPI_ASYNC_REQ) then !regular
          do js=1,this%DataVol,ji
           jv=int(min(this%DataVol-js+1,ji),INT_MPI)
           jtarg=this%Offset+((js-1)*jdts)/jdu
           call nvtx_push('MPI_Get'//CHAR_NULL,4)
           tm=time_sys_sec()
           call MPI_Get(c4_arr(js:js+jv-1),jv,MPI_COMPLEX8,this%RankMPI,jtarg,jv,MPI_COMPLEX8,this%WinMPI%Window,jerr)
           tm=time_sys_sec()-tm
           call nvtx_pop()
           if(LOGGING.gt.0) then
            write(jo,'("#MSG(DDSS::get)[",i5,":",i3,"]: MPI_GET(",i13,",",i5,") time (sec) = ",F8.4)')&
            &impir,thread_id,this%WinMPI%Window,this%RankMPI,tm
            flush(jo)
           endif
           if(jerr.ne.0) then
            if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.GetData): MPI_Get failed!')
            exit
           endif
          enddo
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_NRM
           this%TimeStarted=time_sys_sec()
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=1
          endif
         else !request-handle
          jv=this%DataVol
          call nvtx_push('MPI_Rget'//CHAR_NULL,5)
          tm=time_sys_sec()
          call MPI_Rget(c4_arr,jv,MPI_COMPLEX8,this%RankMPI,this%Offset,jv,MPI_COMPLEX8,this%WinMPI%Window,this%ReqHandle,jerr)
          tm=time_sys_sec()-tm
          call nvtx_pop()
          if(LOGGING.gt.0) then
           write(jo,'("#MSG(DDSS::get)[",i5,":",i3,"]: MPI_RGET(",i13,",",i5,") time (sec) = ",F8.4)')&
           &impir,thread_id,this%WinMPI%Window,this%RankMPI,tm
           flush(jo)
          endif
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_REQ
           this%TimeStarted=time_sys_sec()
          else
           if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.GetData): MPI_Rget failed!')
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=2
          endif
         endif
         return
         end subroutine start_get_c4

         subroutine start_get_c8(c8_arr,jerr)
         complex(8), intent(inout):: c8_arr(*) !asynchronous
         integer(INT_MPI), intent(out):: jerr
         integer(INT_COUNT):: ji,js
         integer(INT_MPI):: jv,jdts,jdu
         integer(INT_ADDR):: jtarg
         real(8):: tm

         jerr=0; jdts=data_type_size(C8)
         jdu=this%WinMPI%DispUnit
         ji=MAX_MPI_MSG_VOL
         if(asnc.ne.MPI_ASYNC_REQ) then !regular
          do js=1,this%DataVol,ji
           jv=int(min(this%DataVol-js+1,ji),INT_MPI)
           jtarg=this%Offset+((js-1)*jdts)/jdu
           call nvtx_push('MPI_Get'//CHAR_NULL,4)
           tm=time_sys_sec()
           call MPI_Get(c8_arr(js:js+jv-1),jv,MPI_COMPLEX16,this%RankMPI,jtarg,jv,MPI_COMPLEX16,this%WinMPI%Window,jerr)
           tm=time_sys_sec()-tm
           call nvtx_pop()
           if(LOGGING.gt.0) then
            write(jo,'("#MSG(DDSS::get)[",i5,":",i3,"]: MPI_GET(",i13,",",i5,") time (sec) = ",F8.4)')&
            &impir,thread_id,this%WinMPI%Window,this%RankMPI,tm
            flush(jo)
           endif
           if(jerr.ne.0) then
            if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.GetData): MPI_Get failed!')
            exit
           endif
          enddo
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_NRM
           this%TimeStarted=time_sys_sec()
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=1
          endif
         else !request-handle
          jv=this%DataVol
          call nvtx_push('MPI_Rget'//CHAR_NULL,5)
          tm=time_sys_sec()
          call MPI_Rget(c8_arr,jv,MPI_COMPLEX16,this%RankMPI,this%Offset,jv,MPI_COMPLEX16,this%WinMPI%Window,this%ReqHandle,jerr)
          tm=time_sys_sec()-tm
          call nvtx_pop()
          if(LOGGING.gt.0) then
           write(jo,'("#MSG(DDSS::get)[",i5,":",i3,"]: MPI_RGET(",i13,",",i5,") time (sec) = ",F8.4)')&
           &impir,thread_id,this%WinMPI%Window,this%RankMPI,tm
           flush(jo)
          endif
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_REQ
           this%TimeStarted=time_sys_sec()
          else
           if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.GetData): MPI_Rget failed!')
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
        integer(INT_MPI):: rwe,errc,asnc
        real(4), pointer, contiguous:: r4_ptr(:)
        real(8), pointer, contiguous:: r8_ptr(:)
        complex(4), pointer, contiguous:: c4_ptr(:)
        complex(8), pointer, contiguous:: c8_ptr(:)

        errc=0
        if(present(async)) then; asnc=async; else; asnc=MPI_ASYNC_NOT; endif !default is synchronous communication
        call this%lock()
        if(asnc.eq.MPI_ASYNC_NOT.or.asnc.eq.MPI_ASYNC_NRM.or.asnc.eq.MPI_ASYNC_REQ) then
         if(.not.c_associated(loc_ptr,C_NULL_PTR)) then
          if(this%RankMPI.ge.0) then
           if(this%StatMPI.eq.MPI_STAT_NONE.or.this%StatMPI.eq.MPI_STAT_COMPLETED.or.&
             &this%StatMPI.eq.MPI_STAT_COMPLETED_ORIG) then
            rwe=RankWinRefs%test(this%RankMPI,this%WinMPI%Window,errc,append=.TRUE.) !get the (rank,window) entry
            if(errc.eq.0) then
             if(this%DataVol.gt.0) then
              if(.not.(asnc.eq.MPI_ASYNC_REQ.and.this%DataVol.gt.huge(asnc))) then
               call modify_lock(rwe,errc) !modify the (rank,window) lock status if needed
               if(errc.eq.0) then
                if(DEBUG.ge.1) then
                 write(jo,'("#DEBUG(distributed:DataDescr.AccData)[",i5,":",i3,"]: Accumulate: ",i18,1x,i9,": ")',ADVANCE='NO')&
                 &impir,thread_id,this%Offset,this%DataVol
                 call RankWinRefs%RankWins(rwe)%print_it(dev_out=jo)
                 flush(jo)
                endif
                select case(this%DataType)
                case(R4)
                 call c_f_pointer(loc_ptr,r4_ptr,(/this%DataVol/))
                 call start_acc_r4(r4_ptr,errc); if(errc.ne.0) errc=1
                case(R8)
                 call c_f_pointer(loc_ptr,r8_ptr,(/this%DataVol/))
                 call start_acc_r8(r8_ptr,errc); if(errc.ne.0) errc=2
                case(C4)
                 call c_f_pointer(loc_ptr,c4_ptr,(/this%DataVol/))
                 call start_acc_c4(c4_ptr,errc); if(errc.ne.0) errc=3
                case(C8)
                 call c_f_pointer(loc_ptr,c8_ptr,(/this%DataVol/))
                 call start_acc_c8(c8_ptr,errc); if(errc.ne.0) errc=4
                case(NO_TYPE)
                 errc=5
                case default
                 errc=6
                end select
                if(errc.eq.0) then
                 call RankWinRefs%new_transfer(this,rwe,WRITE_SIGN,errc) !register a new transfer (will also set this%TransID field)
                 if(errc.eq.0) then
                  if(asnc.eq.MPI_ASYNC_NOT) then
                   call this%flush_data(errc); if(errc.ne.0) errc=7
                  endif
                 else
                  errc=8
                 endif
                endif
               else
                errc=9
               endif
              else
               errc=10
              endif
             elseif(this%DataVol.eq.0) then
              this%StatMPI=MPI_STAT_COMPLETED
             else
              errc=11
             endif
            else
             if(errc.ne.TRY_LATER) errc=12
            endif
           else
            errc=13
           endif
          else
           errc=14
          endif
         else
          errc=15
         endif
        else
         errc=16
        endif
        call this%unlock()
        if(present(ierr)) ierr=errc
        return

        contains

         subroutine modify_lock(rw,jerr)
         integer(INT_MPI), intent(in):: rw
         integer(INT_MPI), intent(out):: jerr
         type(RankWin_t), pointer:: rw_entry
         real(8):: tm

         jerr=0; rw_entry=>RankWinRefs%RankWins(rw)
         if(rw_entry%LockType*WRITE_SIGN.lt.0) then !communication direction change
          if(LAZY_LOCKING) then
           if(DEBUG.ge.1) then
            write(jo,'("#DEBUG(distributed:DataDescr.AccData)[",i5,":",i3,"]: WIN_FLUSH(.acc): ")',ADVANCE='NO') impir,thread_id
            call rw_entry%print_it(dev_out=jo)
            flush(jo)
           endif
           call nvtx_push('MPI_Win_flush'//CHAR_NULL,2)
           tm=time_sys_sec()
           if(.not.NO_FLUSH_AFTER_READ_EPOCH) call MPI_Win_flush(rw_entry%Rank,rw_entry%Window,jerr)
           tm=time_sys_sec()-tm
           call nvtx_pop()
           if(LOGGING.gt.0) then
            write(jo,'("#MSG(DDSS::accumulate)[",i5,":",i3,"]: MPI_WIN_FLUSH(",i13,",",i5,") time (sec) = ",F8.4)')&
            &impir,thread_id,rw_entry%Window,rw_entry%Rank,tm
            flush(jo)
           endif
          else
           if(DEBUG.ge.1) then
            write(jo,'("#DEBUG(distributed:DataDescr.AccData)[",i5,":",i3,"]: WIN_UNLOCK(.acc): ")',ADVANCE='NO') impir,thread_id
            call rw_entry%print_it(dev_out=jo)
            flush(jo)
           endif
           call nvtx_push('MPI_Win_unlock'//CHAR_NULL,3)
           tm=time_sys_sec()
           call MPI_Win_unlock(rw_entry%Rank,rw_entry%Window,jerr)
           tm=time_sys_sec()-tm
           call nvtx_pop()
           if(LOGGING.gt.0) then
            write(jo,'("#MSG(DDSS::accumulate)[",i5,":",i3,"]: MPI_WIN_UNLOCK(",i13,",",i5,") time (sec) = ",F8.4)')&
            &impir,thread_id,rw_entry%Window,rw_entry%Rank,tm
            flush(jo)
           endif
          endif
          if(jerr.eq.0) then
           if(LAZY_LOCKING) then
            rw_entry%LockType=SHARED_LOCK*WRITE_SIGN
           else
            rw_entry%LockType=NO_LOCK
           endif
           rw_entry%LastSync=RankWinRefs%TransCount
          else
           if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.AccData): MPI_Win_unlock/MPI_Win_flush failed!')
           jerr=1
          endif
         endif
         if(jerr.eq.0.and.rw_entry%LockType.eq.NO_LOCK) then
          if(DEBUG.ge.1) then
           write(jo,'("#DEBUG(distributed:DataDescr.AccData)[",i5,":",i3,"]: WIN_LOCK(.acc): ")',ADVANCE='NO') impir,thread_id
           call rw_entry%print_it(dev_out=jo)
           flush(jo)
          endif
          call nvtx_push('MPI_Win_lock'//CHAR_NULL,0)
          tm=time_sys_sec()
          call MPI_Win_lock(MPI_LOCK_SHARED,rw_entry%Rank,MPI_ASSER,rw_entry%Window,jerr)
          tm=time_sys_sec()-tm
          call nvtx_pop()
          if(LOGGING.gt.0) then
           write(jo,'("#MSG(DDSS::accumulate)[",i5,":",i3,"]: MPI_WIN_LOCK(",i13,",",i5,") time (sec) = ",F8.4)')&
           &impir,thread_id,rw_entry%Window,rw_entry%Rank,tm
           flush(jo)
          endif
          if(jerr.eq.0) then
           rw_entry%LockType=SHARED_LOCK*WRITE_SIGN
          else
           if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.AccData): MPI_lock failed!')
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
         integer(INT_MPI):: jv,jdts,jdu
         integer(INT_ADDR):: jtarg
         real(8):: tm

         jerr=0; jdts=data_type_size(R4)
         jdu=this%WinMPI%DispUnit
         ji=MAX_MPI_MSG_VOL
         if(asnc.ne.MPI_ASYNC_REQ) then !regular
          do js=1,this%DataVol,ji
           jv=int(min(this%DataVol-js+1,ji),INT_MPI)
           jtarg=this%Offset+((js-1)*jdts)/jdu
           call nvtx_push('MPI_Accumulate'//CHAR_NULL,6)
           tm=time_sys_sec()
           call MPI_Accumulate(r4_arr(js:js+jv-1),jv,MPI_REAL4,this%RankMPI,jtarg,jv,MPI_REAL4,MPI_SUM,&
                              &this%WinMPI%Window,jerr)
           tm=time_sys_sec()-tm
           call nvtx_pop()
           if(LOGGING.gt.0) then
            write(jo,'("#MSG(DDSS::accumulate)[",i5,":",i3,"]: MPI_ACCUMULATE(",i13,",",i5,") time (sec) = ",F8.4)')&
            &impir,thread_id,this%WinMPI%Window,this%RankMPI,tm
            flush(jo)
           endif
           if(jerr.ne.0) then
            if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.AccData): MPI_Accumulate failed!')
            exit
           endif
          enddo
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_NRM
           this%TimeStarted=time_sys_sec()
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=1
          endif
         else !request-handle
          jv=this%DataVol
          call nvtx_push('MPI_Raccumulate'//CHAR_NULL,7)
          tm=time_sys_sec()
          call MPI_Raccumulate(r4_arr,jv,MPI_REAL4,this%RankMPI,this%Offset,jv,MPI_REAL4,MPI_SUM,&
                              &this%WinMPI%Window,this%ReqHandle,jerr)
          tm=time_sys_sec()-tm
          call nvtx_pop()
          if(LOGGING.gt.0) then
           write(jo,'("#MSG(DDSS::accumulate)[",i5,":",i3,"]: MPI_RACCUMULATE(",i13,",",i5,") time (sec) = ",F8.4)')&
           &impir,thread_id,this%WinMPI%Window,this%RankMPI,tm
           flush(jo)
          endif
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_REQ
           this%TimeStarted=time_sys_sec()
          else
           if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.AccData): MPI_Raccumulate failed!')
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=2
          endif
         endif
         return
         end subroutine start_acc_r4

         subroutine start_acc_r8(r8_arr,jerr)
         real(8), intent(inout):: r8_arr(*) !asynchronous
         integer(INT_MPI), intent(out):: jerr
         integer(INT_COUNT):: ji,js
         integer(INT_MPI):: jv,jdts,jdu
         integer(INT_ADDR):: jtarg
         real(8):: tm

         jerr=0; jdts=data_type_size(R8)
         jdu=this%WinMPI%DispUnit
         ji=MAX_MPI_MSG_VOL
         if(asnc.ne.MPI_ASYNC_REQ) then !regular
          do js=1,this%DataVol,ji
           jv=int(min(this%DataVol-js+1,ji),INT_MPI)
           jtarg=this%Offset+((js-1)*jdts)/jdu
           call nvtx_push('MPI_Accumulate'//CHAR_NULL,6)
           tm=time_sys_sec()
           call MPI_Accumulate(r8_arr(js:js+jv-1),jv,MPI_REAL8,this%RankMPI,jtarg,jv,MPI_REAL8,MPI_SUM,&
                              &this%WinMPI%Window,jerr)
           tm=time_sys_sec()-tm
           call nvtx_pop()
           if(LOGGING.gt.0) then
            write(jo,'("#MSG(DDSS::accumulate)[",i5,":",i3,"]: MPI_ACCUMULATE(",i13,",",i5,") time (sec) = ",F8.4)')&
            &impir,thread_id,this%WinMPI%Window,this%RankMPI,tm
            flush(jo)
           endif
           if(jerr.ne.0) then
            if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.AccData): MPI_Accumulate failed!')
            exit
           endif
          enddo
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_NRM
           this%TimeStarted=time_sys_sec()
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=1
          endif
         else !request-handle
          jv=this%DataVol
          call nvtx_push('MPI_Raccumulate'//CHAR_NULL,7)
          tm=time_sys_sec()
          call MPI_Raccumulate(r8_arr,jv,MPI_REAL8,this%RankMPI,this%Offset,jv,MPI_REAL8,MPI_SUM,&
                              &this%WinMPI%Window,this%ReqHandle,jerr)
          tm=time_sys_sec()-tm
          call nvtx_pop()
          if(LOGGING.gt.0) then
           write(jo,'("#MSG(DDSS::accumulate)[",i5,":",i3,"]: MPI_RACCUMULATE(",i13,",",i5,") time (sec) = ",F8.4)')&
           &impir,thread_id,this%WinMPI%Window,this%RankMPI,tm
           flush(jo)
          endif
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_REQ
           this%TimeStarted=time_sys_sec()
          else
           if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.AccData): MPI_Raccumulate failed!')
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=2
          endif
         endif
         return
         end subroutine start_acc_r8

         subroutine start_acc_c4(c4_arr,jerr)
         complex(4), intent(inout):: c4_arr(*) !asynchronous
         integer(INT_MPI), intent(out):: jerr
         integer(INT_COUNT):: ji,js
         integer(INT_MPI):: jv,jdts,jdu
         integer(INT_ADDR):: jtarg
         real(8):: tm

         jerr=0; jdts=data_type_size(C4)
         jdu=this%WinMPI%DispUnit
         ji=MAX_MPI_MSG_VOL
         if(asnc.ne.MPI_ASYNC_REQ) then !regular
          do js=1,this%DataVol,ji
           jv=int(min(this%DataVol-js+1,ji),INT_MPI)
           jtarg=this%Offset+((js-1)*jdts)/jdu
           call nvtx_push('MPI_Accumulate'//CHAR_NULL,6)
           tm=time_sys_sec()
           call MPI_Accumulate(c4_arr(js:js+jv-1),jv,MPI_COMPLEX8,this%RankMPI,jtarg,jv,MPI_COMPLEX8,MPI_SUM,&
                              &this%WinMPI%Window,jerr)
           tm=time_sys_sec()-tm
           call nvtx_pop()
           if(LOGGING.gt.0) then
            write(jo,'("#MSG(DDSS::accumulate)[",i5,":",i3,"]: MPI_ACCUMULATE(",i13,",",i5,") time (sec) = ",F8.4)')&
            &impir,thread_id,this%WinMPI%Window,this%RankMPI,tm
            flush(jo)
           endif
           if(jerr.ne.0) then
            if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.AccData): MPI_Accumulate failed!')
            exit
           endif
          enddo
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_NRM
           this%TimeStarted=time_sys_sec()
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=1
          endif
         else !request-handle
          jv=this%DataVol
          call nvtx_push('MPI_Raccumulate'//CHAR_NULL,7)
          tm=time_sys_sec()
          call MPI_Raccumulate(c4_arr,jv,MPI_COMPLEX8,this%RankMPI,this%Offset,jv,MPI_COMPLEX8,MPI_SUM,&
                              &this%WinMPI%Window,this%ReqHandle,jerr)
          tm=time_sys_sec()-tm
          call nvtx_pop()
          if(LOGGING.gt.0) then
           write(jo,'("#MSG(DDSS::accumulate)[",i5,":",i3,"]: MPI_RACCUMULATE(",i13,",",i5,") time (sec) = ",F8.4)')&
           &impir,thread_id,this%WinMPI%Window,this%RankMPI,tm
           flush(jo)
          endif
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_REQ
           this%TimeStarted=time_sys_sec()
          else
           if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.AccData): MPI_Raccumulate failed!')
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=2
          endif
         endif
         return
         end subroutine start_acc_c4

         subroutine start_acc_c8(c8_arr,jerr)
         complex(8), intent(inout):: c8_arr(*) !asynchronous
         integer(INT_MPI), intent(out):: jerr
         integer(INT_COUNT):: ji,js
         integer(INT_MPI):: jv,jdts,jdu
         integer(INT_ADDR):: jtarg
         real(8):: tm

         jerr=0; jdts=data_type_size(C8)
         jdu=this%WinMPI%DispUnit
         ji=MAX_MPI_MSG_VOL
         if(asnc.ne.MPI_ASYNC_REQ) then !regular
          do js=1,this%DataVol,ji
           jv=int(min(this%DataVol-js+1,ji),INT_MPI)
           jtarg=this%Offset+((js-1)*jdts)/jdu
           call nvtx_push('MPI_Accumulate'//CHAR_NULL,6)
           tm=time_sys_sec()
           call MPI_Accumulate(c8_arr(js:js+jv-1),jv,MPI_COMPLEX16,this%RankMPI,jtarg,jv,MPI_COMPLEX16,MPI_SUM,&
                              &this%WinMPI%Window,jerr)
           tm=time_sys_sec()-tm
           call nvtx_pop()
           if(LOGGING.gt.0) then
            write(jo,'("#MSG(DDSS::accumulate)[",i5,":",i3,"]: MPI_ACCUMULATE(",i13,",",i5,") time (sec) = ",F8.4)')&
            &impir,thread_id,this%WinMPI%Window,this%RankMPI,tm
            flush(jo)
           endif
           if(jerr.ne.0) then
            if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.AccData): MPI_Accumulate failed!')
            exit
           endif
          enddo
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_NRM
           this%TimeStarted=time_sys_sec()
          else
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=1
          endif
         else !request-handle
          jv=this%DataVol
          call nvtx_push('MPI_Raccumulate'//CHAR_NULL,7)
          tm=time_sys_sec()
          call MPI_Raccumulate(c8_arr,jv,MPI_COMPLEX16,this%RankMPI,this%Offset,jv,MPI_COMPLEX16,MPI_SUM,&
                              &this%WinMPI%Window,this%ReqHandle,jerr)
          tm=time_sys_sec()-tm
          call nvtx_pop()
          if(LOGGING.gt.0) then
           write(jo,'("#MSG(DDSS::accumulate)[",i5,":",i3,"]: MPI_RACCUMULATE(",i13,",",i5,") time (sec) = ",F8.4)')&
           &impir,thread_id,this%WinMPI%Window,this%RankMPI,tm
           flush(jo)
          endif
          if(jerr.eq.0) then
           this%StatMPI=MPI_STAT_PROGRESS_REQ
           this%TimeStarted=time_sys_sec()
          else
           if(DDSS_MPI_ERR_FATAL) call quit(jerr,'#FATAL(distributed:DataDescr.AccData): MPI_Raccumulate failed!')
           this%StatMPI=MPI_STAT_ONESIDED_ERR; jerr=2
          endif
         endif
         return
         end subroutine start_acc_c8

        end subroutine DataDescrAccData
!-------------------------------------------------------------
        subroutine DataDescrPackInt(this,packet,ierr,pack_len)
!Packs a data descriptor into a plain packet of integers of kind ELEM_PACK_SIZE.
!The first integer in the packet is always the useful length of the packet,
!that is, the number of following integer elements carrying the information.
!It is the user responsibility to provide a large enough packet buffer.
        use extern_names, only: c_ptr_value
        implicit none
        class(DataDescr_t), intent(inout):: this                     !in: data descriptor
        integer(ELEM_PACK_SIZE), intent(inout), contiguous, target:: packet(0:) !out: plain integer packet (length + descriptor data)
        integer(INT_MPI), intent(inout), optional:: ierr             !out: error code (0:success)
        integer(INT_COUNT), intent(out), optional:: pack_len         !out: total packet length (in packing integers)
        integer(INT_MPI):: errc
        integer(INT_COUNT):: pl,wl
        integer(C_SIZE_T), pointer:: isize_p
        integer(INT_ADDR), pointer:: iaddr_p
        integer(INT_MPI), pointer:: impi_p
        integer(INT_COUNT), pointer:: len_p
        type(C_PTR):: cptr

        errc=0; if(present(pack_len)) pack_len=0
        call this%lock()
        packet(0)=0; pl=0
        pl=pl+1; packet(pl)=0; cptr=c_loc(packet(pl)); call c_f_pointer(cptr,isize_p); isize_p=c_ptr_value(this%LocPtr)
        pl=pl+1; packet(pl)=0; cptr=c_loc(packet(pl)); call c_f_pointer(cptr,impi_p); impi_p=this%RankMPI
        pl=pl+1; call this%WinMPI%pack(packet(pl:),errc,wl)
        if(errc.eq.0) then
         pl=pl+wl; packet(pl)=0; cptr=c_loc(packet(pl)); call c_f_pointer(cptr,iaddr_p); iaddr_p=this%Offset
         pl=pl+1; packet(pl)=0; cptr=c_loc(packet(pl)); call c_f_pointer(cptr,len_p); len_p=this%DataVol
         pl=pl+1; packet(pl)=0; cptr=c_loc(packet(pl)); call c_f_pointer(cptr,impi_p); impi_p=this%DataType
         cptr=c_loc(packet(0)); call c_f_pointer(cptr,len_p); len_p=pl !packet body length
         len_p=>NULL(); impi_p=>NULL(); iaddr_p=>NULL(); isize_p=>NULL()
         if(present(pack_len)) pack_len=1+pl !header integer + packet body
        else
         errc=1
        endif
        call this%unlock()
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrPackInt
!----------------------------------------------------
        subroutine DataDescrPackNew(this,packet,ierr)
!Packs an object into a packet.
        implicit none
        class(DataDescr_t), intent(inout):: this       !in: data descriptor
        class(obj_pack_t), intent(inout):: packet      !inout: packet
        integer(INT_MPI), intent(out), optional:: ierr !out: error code
        integer(INT_MPI):: errc

        errc=0
        call this%lock()
        call pack_builtin(packet,this%RankMPI,errc)
        if(errc.eq.PACK_SUCCESS) call this%WinMPI%pack(packet,errc)
        if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%Offset,errc)
        if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%DataVol,errc)
        if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%DataType,errc)
        call this%unlock()
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrPackNew
!----------------------------------------------------------
        subroutine DataDescrPack(this,packet,ierr,pack_len)
!Packs a data descriptor into a plain packet of integers of kind ELEM_PACK_SIZE.
!The first integer in the packet is always the useful length of the packet,
!that is, the number of following integer elements carrying the information.
!It is the user responsibility to provide a large enough packet buffer.
        use extern_names, only: c_ptr_value
        implicit none
        class(DataDescr_t), intent(inout):: this             !in: data descriptor
        class(SimplePack_t), intent(inout), target:: packet  !out: plain integer packet (length + descriptor data)
        integer(INT_MPI), intent(inout), optional:: ierr     !out: error code (0:success)
        integer(INT_COUNT), intent(out), optional:: pack_len !out: total packet length (in packing integers)
        integer(INT_MPI):: errc
        integer(INT_COUNT):: pl,wl
        integer(C_SIZE_T), pointer:: isize_p
        integer(INT_ADDR), pointer:: iaddr_p
        integer(INT_MPI), pointer:: impi_p
        integer(INT_COUNT), pointer:: len_p
        type(C_PTR):: cptr

        errc=0; if(present(pack_len)) pack_len=0
        call this%lock()
        if(packet%buf_len().lt.1+DataDescr_PACK_LEN) call packet%reserve_mem(int(1+DataDescr_PACK_LEN,INT_MPI),errc)
        if(errc.eq.0) then
         packet%Packet(0)=0; pl=0
         pl=pl+1; packet%Packet(pl)=0; cptr=c_loc(packet%Packet(pl))
         call c_f_pointer(cptr,isize_p); isize_p=c_ptr_value(this%LocPtr)
         pl=pl+1; packet%Packet(pl)=0; cptr=c_loc(packet%Packet(pl))
         call c_f_pointer(cptr,impi_p); impi_p=this%RankMPI
         pl=pl+1; call this%WinMPI%pack(packet%Packet(pl:),errc,wl)
         if(errc.eq.0) then
          pl=pl+wl; packet%Packet(pl)=0; cptr=c_loc(packet%Packet(pl))
          call c_f_pointer(cptr,iaddr_p); iaddr_p=this%Offset
          pl=pl+1; packet%Packet(pl)=0; cptr=c_loc(packet%Packet(pl))
          call c_f_pointer(cptr,len_p); len_p=this%DataVol
          pl=pl+1; packet%Packet(pl)=0; cptr=c_loc(packet%Packet(pl))
          call c_f_pointer(cptr,impi_p); impi_p=this%DataType
          cptr=c_loc(packet%Packet(0)); call c_f_pointer(cptr,len_p); len_p=pl !packet body length
          len_p=>NULL(); impi_p=>NULL(); iaddr_p=>NULL(); isize_p=>NULL()
          if(present(pack_len)) pack_len=1+pl !header integer + packet body
         else
          errc=2
         endif
        else
         errc=1
        endif
        call this%unlock()
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrPack
!------------------------------------------------------
        subroutine DataDescrUnpackNew(this,packet,ierr)
!Unpacks an object from a packet.
        implicit none
        class(DataDescr_t), intent(out):: this         !out: data descriptor
        class(obj_pack_t), intent(inout):: packet      !inout: packet
        integer(INT_MPI), intent(out), optional:: ierr !out: error code
        integer(INT_MPI):: errc

        errc=0
        call this%lock()
        call unpack_builtin(packet,this%RankMPI,errc)
        if(errc.eq.PACK_SUCCESS) call this%WinMPI%unpack(packet,errc)
        if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%Offset,errc)
        if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%DataVol,errc)
        if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%DataType,errc)
        call this%unlock()
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrUnpackNew
!---------------------------------------------------------------
        subroutine DataDescrUnpackInt(this,packet,ierr,pack_len)
!Unpacks a data descriptor from a plain integer packet.
        use extern_names, only: c_ptr_set
        implicit none
        class(DataDescr_t), intent(inout):: this                  !out: unpacked data descriptor
        integer(ELEM_PACK_SIZE), intent(in), contiguous, target:: packet(0:) !in: plain integer packet containing the data descriptor information
        integer(INT_MPI), intent(inout), optional:: ierr          !out: error code (0:success)
        integer(INT_COUNT), intent(out), optional:: pack_len      !out: total packet length (in packing integers)
        integer(INT_MPI):: errc
        integer(INT_COUNT):: pl,wl
        integer(C_SIZE_T), pointer:: isize_p
        integer(INT_ADDR), pointer:: iaddr_p
        integer(INT_MPI), pointer:: impi_p
        integer(INT_COUNT), pointer:: len_p
        type(C_PTR):: cptr

        errc=0
        call this%lock()
        cptr=c_loc(packet(0)); call c_f_pointer(cptr,len_p)
        if(present(pack_len)) pack_len=1+len_p
        if(len_p.eq.DataDescr_PACK_LEN) then
         call this%clean(); pl=0
         pl=pl+1; cptr=c_loc(packet(pl)); call c_f_pointer(cptr,isize_p); call c_ptr_set(isize_p,this%LocPtr)
         pl=pl+1; cptr=c_loc(packet(pl)); call c_f_pointer(cptr,impi_p); this%RankMPI=impi_p
         pl=pl+1; call this%WinMPI%unpack(packet(pl:),errc,wl)
         if(errc.eq.0) then
          pl=pl+wl; cptr=c_loc(packet(pl)); call c_f_pointer(cptr,iaddr_p); this%Offset=iaddr_p
          pl=pl+1; cptr=c_loc(packet(pl)); call c_f_pointer(cptr,len_p); this%DataVol=len_p
          pl=pl+1; cptr=c_loc(packet(pl)); call c_f_pointer(cptr,impi_p); this%DataType=impi_p
         else
          errc=2
         endif
         len_p=>NULL(); impi_p=>NULL(); iaddr_p=>NULL(); isize_p=>NULL()
        else
         errc=1
        endif
        call this%unlock()
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrUnpackInt
!------------------------------------------------------------
        subroutine DataDescrUnpack(this,packet,ierr,pack_len)
!Unpacks a data descriptor from a plain integer packet.
        use extern_names, only: c_ptr_set
        implicit none
        class(DataDescr_t), intent(inout):: this             !out: unpacked data descriptor
        class(SimplePack_t), intent(in), target:: packet     !in: plain integer packet containing the data descriptor information
        integer(INT_MPI), intent(inout), optional:: ierr     !out: error code (0:success)
        integer(INT_COUNT), intent(out), optional:: pack_len !out: total packet length (in packing integers)
        integer(INT_MPI):: errc
        integer(INT_COUNT):: pl,wl
        integer(C_SIZE_T), pointer:: isize_p
        integer(INT_ADDR), pointer:: iaddr_p
        integer(INT_MPI), pointer:: impi_p
        integer(INT_COUNT), pointer:: len_p
        type(C_PTR):: cptr

        errc=0
        call this%lock()
        if(packet%body_len().ge.DataDescr_PACK_LEN) then
         cptr=c_loc(packet%Packet(0)); call c_f_pointer(cptr,len_p)
         if(present(pack_len)) pack_len=1+len_p
         if(len_p.eq.DataDescr_PACK_LEN) then
          call this%clean(); pl=0
          pl=pl+1; cptr=c_loc(packet%Packet(pl)); call c_f_pointer(cptr,isize_p); call c_ptr_set(isize_p,this%LocPtr)
          pl=pl+1; cptr=c_loc(packet%Packet(pl)); call c_f_pointer(cptr,impi_p); this%RankMPI=impi_p
          pl=pl+1; call this%WinMPI%unpack(packet%Packet(pl:),errc,wl)
          if(errc.eq.0) then
           pl=pl+wl; cptr=c_loc(packet%Packet(pl)); call c_f_pointer(cptr,iaddr_p); this%Offset=iaddr_p
           pl=pl+1; cptr=c_loc(packet%Packet(pl)); call c_f_pointer(cptr,len_p); this%DataVol=len_p
           pl=pl+1; cptr=c_loc(packet%Packet(pl)); call c_f_pointer(cptr,impi_p); this%DataType=impi_p
          else
           errc=2
          endif
          len_p=>NULL(); impi_p=>NULL(); iaddr_p=>NULL(); isize_p=>NULL()
         else
          errc=1
         endif
        else
         errc=-1
        endif
        call this%unlock()
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrUnpack
!----------------------------------------------------
        subroutine DataDescrPrint(this,dev_out,space)
!Prints object data.
        use stsubs, only: numchar
!       use extern_names, only: print_c_ptr
        implicit none
        class(DataDescr_t), intent(inout):: this         !in: object to print
        integer(INT_MPI), intent(in), optional:: dev_out !in: output device
        integer(INT_MPI), intent(in), optional:: space   !in: left alignment
        character(32):: sfmt
        integer(INT_MPI):: devo,sp
        integer:: fl

        devo=jo; if(present(dev_out)) devo=dev_out
        sp=0; if(present(space)) sp=space
        call this%lock()
        if(sp.gt.0) then; call numchar(int(sp),fl,sfmt); sfmt(fl+1:fl+2)='x,'; fl=fl+2; else; fl=0; endif
        write(devo,'('//sfmt(1:fl)//'"#Printing DataDescr_t object:")')
!       write(devo,'('//sfmt(1:fl)//'"  C pointer at data origin: ")',ADVANCE='NO'); write(devo,*) this%LocPtr
        write(devo,'('//sfmt(1:fl)//'"  Origin MPI rank         : ",i18)') this%RankMPI
        write(devo,'('//sfmt(1:fl)//'"  Origin data displacement: ",i18)') this%Offset
        write(devo,'('//sfmt(1:fl)//'"  Data volume (elements)  : ",i18)') this%DataVol
        write(devo,'('//sfmt(1:fl)//'"  Data type               : ",i18)') this%DataType
        call this%WinMPI%print_it(devo,sp+2)
        write(devo,'('//sfmt(1:fl)//'"  Current data transfer ID: ",i18)') this%TransID
        write(devo,'('//sfmt(1:fl)//'"  MPI status              : ",i18)') this%StatMPI
        write(devo,'('//sfmt(1:fl)//'"  MPI request handle      : ",i18)') this%ReqHandle
        call this%unlock()
        return
        end subroutine DataDescrPrint
!-------------------------------------
        subroutine DataDescrLock(this)
        implicit none
        class(DataDescr_t), intent(inout):: this

        !write(*,*, ADVANCE='NO') 'Locking data descriptor ...' !debug
        call this%ObjLock%lock()
!$OMP FLUSH(this)
        !write(*,*) 'Locked' !debug
        return
        end subroutine DataDescrLock
!---------------------------------------
        subroutine DataDescrUnlock(this)
        implicit none
        class(DataDescr_t), intent(inout):: this

        !write(*,*) 'Unlocked data descriptor' !debug
!$OMP FLUSH(this)
        call this%ObjLock%unlock()
        return
        end subroutine DataDescrUnlock
!------------------------------------------
        subroutine DataDescrClearLock(this)
        implicit none
        class(DataDescr_t), intent(inout):: this
!$OMP FLUSH(this)
        call this%ObjLock%clear()
        return
        end subroutine DataDescrClearLock
!-------------------------------------
        subroutine DataDescrDtor(this)
        implicit none
        type(DataDescr_t):: this

        return
        end subroutine DataDescrDtor
!==========================================================
        subroutine SimplePackReserve(this,vol,ierr,ext_buf)
!Reserves a new memory buffer. If the memory buffer had already been reserved,
!it will be freed and the new one will be reserved.
        implicit none
        class(SimplePack_t), intent(inout):: this        !inout: simple packet
        integer(INT_MPI), intent(in):: vol               !in: requested volume of the memory buffer (in integers of kind ELEM_PACK_SIZE)
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(ELEM_PACK_SIZE), contiguous, target, intent(in), optional:: ext_buf(:) !in: external buffer
        integer(INT_MPI):: l,errc
        integer:: i

        errc=0
        if(vol.gt.0) then
         if(associated(this%Packet)) then
          if(this%Alloc) deallocate(this%Packet)
         endif
         nullify(this%Packet); this%Alloc=.FALSE.
         if(present(ext_buf)) then
          if(size(ext_buf).ge.vol) then
           l=lbound(ext_buf,1)
           this%Packet(0:vol-1)=>ext_buf(l:l+vol-1)
           this%Alloc=.FALSE.
          else
           errc=3
          endif
         else
          allocate(this%Packet(0:vol-1),STAT=i)
          if(i.eq.0) then; this%Alloc=.TRUE.; else; nullify(this%Packet); errc=2; endif
         endif
         if(associated(this%Packet)) this%Packet(0)=0 !empty packet flag
        else
         errc=1
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine SimplePackReserve
!--------------------------------------------------------
        function SimplePackBufLen(this,ierr) result(blen)
!Returns the memory buffer length in integers of kind ELEM_PACK_SIZE.
        implicit none
        integer(INT_MPI):: blen                          !out: memory buffer length
        class(SimplePack_t), intent(in):: this           !in: simple packet
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0; blen=0
        if(associated(this%Packet)) blen=size(this%Packet)
        if(present(ierr)) ierr=errc
        return
        end function SimplePackBufLen
!---------------------------------------------------------
        function SimplePackFullLen(this,ierr) result(blen)
!Returns the full length of the packet in integers of kind ELEM_PACK_SIZE.
!0 means no packet resides in the memory buffer.
        implicit none
        integer(INT_MPI):: blen                          !out: full packet length
        class(SimplePack_t), intent(in):: this           !in: simple packet
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc
        integer(INT_COUNT):: l

        errc=0; blen=0
        if(associated(this%Packet)) then
         blen=packet_full_len(this%Packet,l)
         if(l.le.0) blen=0
        else
         errc=-1
        endif
        if(present(ierr)) ierr=errc
        return
        end function SimplePackFullLen
!--------------------------------------------------------
        function SimplePackBodyLen(this,ierr) result(blen)
!Returns the body length of an existing packet in integers of kind ELEM_PACK_SIZE.
!0 means no packet resides in the memory buffer.
        implicit none
        integer(INT_MPI):: blen                          !out: body length of the packet
        class(SimplePack_t), intent(in):: this           !in: simple packet
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc
        integer(INT_COUNT):: l

        errc=0; blen=0
        if(associated(this%Packet)) then
         blen=packet_full_len(this%Packet,l); blen=l
        else
         errc=-1
        endif
        if(present(ierr)) ierr=errc
        return
        end function SimplePackBodyLen
!----------------------------------------------
        subroutine SimplePackDiscard(this,ierr)
!Discards the packet from the memory buffer.
        implicit none
        class(SimplePack_t), intent(inout):: this        !in: simple packet
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0
        if(associated(this%Packet)) then
         this%Packet(0)=0
        else
         errc=-1
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine SimplePackDiscard
!----------------------------------------------
        subroutine SimplePackDestroy(this,ierr)
!Destroys the memory buffer.
        implicit none
        class(SimplePack_t), intent(inout):: this        !in: simple packet
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0
        if(associated(this%Packet)) then
         if(this%Alloc) deallocate(this%Packet)
         nullify(this%Packet); this%Alloc=.FALSE.
        else
         errc=-1
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine SimplePackDestroy
!==================================================
        subroutine CommHandleClean(this,ierr,force)
!Cleans a communication handle, unless it is still participating in a communication.
!In the latter case, a status TRY_LATER is returned, unless <force>=TRUE, in which
!case the handle will be cleaned anyways.
        implicit none
        class(CommHandle_t), intent(inout):: this        !inout: communication handle
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success) or TRY_LATER
        logical, intent(in), optional:: force            !in: if TRUE, cleaning will be enforced
        integer(INT_MPI):: errc
        logical:: frc

        frc=.FALSE.; if(present(force)) frc=force
        errc=this%test()
        if(errc.ne.MPI_STAT_PROGRESS_REQ.or.frc) then
         call MPI_Request_free(this%ReqHandle,errc)
         this%ReqHandle=MPI_REQUEST_NULL
         this%LastReq=MPI_REQUEST_NULL
         this%MPIRank=-1
         this%DataCont=>NULL()
         this%TimeInitiated=0d0; this%TimeFoundCompleted=-1d0
         errc=0
        else
         errc=TRY_LATER !communication handle is still participating a communication
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine CommHandleClean
!------------------------------------------------------
        subroutine CommHandleWait(this,ierr,ignore_old)
!Waits upon a completion of an active communication (at origin).
!It is safe to wait upon completion of the same message twice and more,
!unless <ignore_old>=TRUE. In the latter case, the message is assumed
!never waited upon completion before.
        implicit none
        class(CommHandle_t), intent(inout):: this        !inout: communication handle
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        logical, intent(in), optional:: ignore_old       !in: if TRUE, the completion of older messages will be ignored
        integer(INT_MPI):: errc
        logical:: igo

        errc=0
        igo=.FALSE.; if(present(ignore_old)) igo=ignore_old
        if(this%ReqHandle.ne.MPI_REQUEST_NULL) then
         call MPI_Wait(this%ReqHandle,this%MPIStat,errc)
         if(errc.eq.0) then
          this%LastReq=this%ReqHandle
          this%ReqHandle=MPI_REQUEST_NULL
          this%TimeFoundCompleted=thread_wtime()
          if(associated(this%DataCont)) then
           call this%DataCont%register_arrived(errc); if(errc.ne.0) errc=3 !corrupted data packet container
          else
           errc=2 !no associated data packet container found
          endif
         else
          errc=1 !MPI_WAIT failed
         endif
        else
         if(.not.igo.and.this%LastReq.ne.MPI_REQUEST_NULL) then
          errc=0 !redundant wait (already completed)
         else
          errc=-1 !empty communication handle
         endif
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine CommHandleWait
!---------------------------------------------------------------------
        function CommHandleTest_(this,ierr,ignore_old) result(msg_stat)
!Tests an active communication for a completion. It is safe to test for
!completion the same already completed message twice and more, unless <ignore_old>=TRUE.
!Possible successful return statuses:
! # MPI_STAT_PROGRESS_REQ: Communication is still in progress (with a request handle);
! # MPI_STAT_COMPLETED_ORIG: Communication has completed locally (MPI P2P semantics);
!A negative return status indicates an error.
        implicit none
        integer(INT_MPI):: msg_stat                      !out: message status
        class(CommHandle_t), intent(inout):: this        !inout: communication handle
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        logical, intent(in), optional:: ignore_old       !in: if TRUE, the completion of older messages will be ignored
        integer(INT_MPI):: errc
        logical:: igo,fin

        errc=0; msg_stat=MPI_STAT_PROGRESS_REQ
        igo=.FALSE.; if(present(ignore_old)) igo=ignore_old
        if(this%ReqHandle.ne.MPI_REQUEST_NULL) then
         call MPI_Test(this%ReqHandle,fin,this%MPIStat,errc)
         if(errc.eq.0) then
          if(fin) then !completed
           this%LastReq=this%ReqHandle
           this%ReqHandle=MPI_REQUEST_NULL
           this%TimeFoundCompleted=thread_wtime()
           msg_stat=MPI_STAT_COMPLETED_ORIG
           if(associated(this%DataCont)) then
            call this%DataCont%register_arrived(errc); if(errc.ne.0) errc=3 !corrupted data packet container
           else
            errc=2 !no associated data packet container found
           endif
          endif
         else
          errc=1 !MPI_TEST failed
         endif
        else
         if(.not.igo.and.this%LastReq.ne.MPI_REQUEST_NULL) then
          msg_stat=MPI_STAT_COMPLETED_ORIG
          errc=0 !redundant wait (already completed)
         else
          msg_stat=-1 !error:
          errc=-1     !empty communication handle
         endif
        endif
        if(present(ierr)) ierr=errc
        return
        end function CommHandleTest_
!======================================================
        subroutine PackContClean(this,ierr,keep_buffer)
!Cleans a data packet container. If <keep_buffer>=TRUE,
!the buffer memory will not be released, so it can be reused.
!In any case, all packets will be lost.
        implicit none
        class(PackCont_t), intent(inout):: this          !inout: data packet container
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        logical, intent(in), optional:: keep_buffer      !in: TRUE keeps the memory buffer (defaults to FALSE)
        integer(INT_MPI):: errc
        logical:: free_buf

        errc=0
        free_buf=.TRUE.; if(present(keep_buffer)) free_buf=.not.(keep_buffer.and.associated(this%Packets))
        if(free_buf) then
         if(associated(this%Packets).and.this%Alloc) deallocate(this%Packets)
         nullify(this%Packets)
         this%ffe=-1
         this%Alloc=.FALSE.
         this%NumPackets=-1 !buffer does not exist
         this%Marked=.FALSE.
        else !keep the memory buffer
         this%ffe=1 !packets(0) is reserved to contain <NumPackets>
         this%NumPackets=0  !buffer exists but empty
         this%Marked=.FALSE.
         this%Packets(0)=int(this%NumPackets,ELEM_PACK_SIZE)
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine PackContClean
!----------------------------------------------------------------------
        subroutine PackContReserve(this,buf_len,ierr,ext_buf,resize_it)
!Reserves memory for the data packet container. Returns TRY_LATER in case
!the memory allocation was unsuccessful. If the memory had already been
!previously reserved, one can resize it by supplying <resize_it>=TRUE,
!otherwise the request will fail.
        implicit none
        class(PackCont_t), intent(inout):: this                             !inout: data packet container
        integer(INT_MPI), intent(in):: buf_len                              !in: buffer length (in integers of kind ELEM_PACK_SIZE)
        integer(INT_MPI), intent(inout), optional:: ierr                    !out: error code (0:success) or TRY_LATER
        integer(ELEM_PACK_SIZE), target, contiguous, optional:: ext_buf(1:) !in: external buffer
        logical, intent(in), optional:: resize_it                           !in: if TRUE, a previously allocated/associated buffer will be released
        integer(INT_MPI):: errc
        integer:: erc
        logical:: rsz

        errc=0
        if(buf_len.gt.0) then
         if(present(resize_it)) then; rsz=resize_it; else; rsz=.FALSE.; endif
         if(this%NumPackets.le.0) then
          if(associated(this%Packets)) then
           if(rsz) then
            call this%clean(errc); if(errc.ne.0) errc=5
           else
            errc=4 !memory already reserved and resizing is not asked for
           endif
          endif
          if(errc.eq.0) then
           if(present(ext_buf)) then !associate to an exeternally provided buffer
            if(size(ext_buf).ge.buf_len) then
             this%Packets(0:buf_len-1)=>ext_buf(1:buf_len)
             this%ffe=1 !packets(0) is reserved to contain <NumPackets>
             this%Alloc=.FALSE.
             this%NumPackets=0
             this%Marked=.FALSE.
            else
             errc=3 !exetrnal buffer is not large enough
            endif
           else !allocate the buffer
            allocate(this%Packets(0:buf_len-1),STAT=erc)
            if(erc.eq.0) then
             this%ffe=1 !packets(0) is reserved to contain <NumPackets>
             this%Alloc=.TRUE.
             this%NumPackets=0
             this%Marked=.FALSE.
            else
             nullify(this%Packets)
             call this%clean()
             errc=TRY_LATER
            endif
           endif
           if(errc.eq.0) this%Packets(0)=int(this%NumPackets,ELEM_PACK_SIZE)
          endif
         else
          errc=2 !non-empty data container, clean it first
         endif
        else
         errc=1
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine PackContReserve
!----------------------------------------------------------
        function PackContNumPacks(this,ierr) result(npacks)
!Returns the total number of data packets in the data packet container.
!A negative return value means either an empty container or an error.
        implicit none
        integer(INT_MPI):: npacks                        !out: total number of packets in the packet container
        class(PackCont_t), intent(in):: this             !in: data packet container
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0; npacks=this%NumPackets
        if(associated(this%Packets)) then
         if(num_packs_in_container(this%Packets(0:)).ne.npacks) errc=1
        else
         errc=-1; npacks=-1
        endif
        if(present(ierr)) ierr=errc
        return
        end function PackContNumPacks
!------------------------------------------------------
        subroutine PackContAppend(this,packet,ierr,tag)
!Appends a packet to the data container. One cannot append a tag-marked
!packet to a non-empty unmarked packet container, and vice versa.
!That is, either all stored packets are tagged or all are not tagged.
!An unmarked data container has the following format:
! [num_packets; packet1; packet2; ...]
!A marked data packet container has the following format:
! [num_packets; tag1, packet1; tag2, packet2; ...]
!Each tag occupies a single integer element of kind ELEM_PACK_SIZE.
!Packet numeration starts from 1.
        implicit none
        class(PackCont_t), intent(inout):: this                   !inout: data packet container (either marked or unmarked)
        integer(ELEM_PACK_SIZE), intent(in), target:: packet(0:*) !in: data packet
        integer(INT_MPI), intent(inout), optional:: ierr          !out: error code (0:success)
        integer(ELEM_PACK_SIZE), intent(in), optional:: tag       !in: marking tag (to mark the packet)
        integer(INT_MPI):: errc
        integer(INT_COUNT):: i,l,n,body_len

        errc=0
        if(associated(this%Packets)) then
         if(present(tag)) then
          if(this%NumPackets.le.0) this%Marked=.TRUE.
          if(this%Marked) then; i=1; else; errc=1; endif
         else
          if(this%NumPackets.le.0) this%Marked=.FALSE.
          if(.not.this%Marked) then; i=0; else; errc=2; endif
         endif
         if(errc.eq.0) then
          n=packet_full_len(packet,body_len)
          if(body_len.gt.0) then
           if(this%ffe+i+n-1.le.ubound(this%Packets,1)) then
            if(i.gt.0) then; this%Packets(this%ffe)=tag; this%ffe=this%ffe+1; endif
            this%Packets(this%ffe:this%ffe+n-1)=packet(0:n-1); this%ffe=this%ffe+n
            this%NumPackets=this%NumPackets+1
            this%Packets(0)=int(this%NumPackets*(1-i*2),ELEM_PACK_SIZE) !inverts the sign when data container is marked (tagged)
           else
            errc=3 !insufficient buffer space: Resize the buffer
           endif
          else
           errc=4 !empty packet
          endif
         endif
        else
         errc=5 !no buffer been reserved
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine PackContAppend
!---------------------------------------------------------------
        subroutine PackContRemove(this,ierr,pack_num,packet,tag)
!Removes a data packet from the data packet container.
!If <pack_num> is absent, the last packet will be removed.
!Packet numeration starts from 1.
!If <packet> is present, it will contain the removed packet.
!For marked (tagged) packet containers, the tag can be returned via <tag>.
        implicit none
        class(PackCont_t), intent(inout):: this                        !inout: data packet container
        integer(INT_MPI), intent(inout), optional:: ierr               !out: error code (0:success)
        integer(INT_MPI), intent(in), optional:: pack_num              !in: optional packet number
        integer(ELEM_PACK_SIZE), intent(inout), optional:: packet(0:*) !out: data packet
        integer(ELEM_PACK_SIZE), intent(out), optional:: tag           !out: data packet tag (if any)
        integer(INT_MPI):: i,pn,errc
        integer(INT_COUNT):: l,k,j,m,n

        errc=0
        if(this%NumPackets.gt.0.and.associated(this%Packets)) then
         if(present(pack_num)) then; pn=pack_num; else; pn=this%NumPackets; endif
         if(pn.gt.0.and.pn.le.this%NumPackets) then
          n=ubound(this%Packets,1)
          l=1; if(this%Marked) then; k=1; else; k=0; endif
          do i=1,pn-1
           if(l+k.le.n) then
            j=packet_full_len(this%Packets(l+k:)); if(j.le.1) then; errc=1; exit; endif !invalid packet length
            j=l+k+j; if(j.gt.l) then; l=j; else; errc=2; exit; endif !corrupted container
           else
            errc=3; exit !search is out-of-bounds
           endif
          enddo
          if(errc.eq.0) then
           j=packet_full_len(this%Packets(l+k:))
           if(j.gt.1) then
            if(k.gt.0.and.present(tag)) tag=this%Packets(l)
            if(present(packet)) packet(0:j-1)=this%Packets(l+k:l+k+j-1)
            if(pn.lt.this%NumPackets) then
             do m=l+k+j,this%ffe-1
              this%Packets(l)=this%Packets(m)
              l=l+1
             enddo
            endif
            this%ffe=l; this%NumPackets=this%NumPackets-1
            this%Packets(0)=int(this%NumPackets*(1-k*2),ELEM_PACK_SIZE)
           else
            errc=4 !invalid packet length
           endif
          endif
         else
          errc=5 !invalid packet number
         endif
        else
         errc=6 !empty data packet container
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine PackContRemove
!---------------------------------------------------------------------------------
        subroutine PackContSend(this,comm_hl,ierr,recv_rank,msg_tag,comm_mpi,sync)
!Sends a data packet container to other MPI process(es).
!If <recv_rank> is not specified, a broadcast will be performed.
!Otherwise, a non-blocking P2P communication will be initiated,
!and optionally completed (<sync>=TRUE). If <msg_tag> is absent,
!the DEFAULT_MPI_TAG will be used as the P2P message tag. Note that
!in case of broadcast all other participating processes must call
!the corresponding receive method!
        implicit none
        class(PackCont_t), intent(in), target:: this       !in: non-empty data packet container
        type(CommHandle_t), intent(inout):: comm_hl        !out: communication handle
        integer(INT_MPI), intent(inout), optional:: ierr   !out: error code (0:success)
        integer(INT_MPI), intent(in), optional:: recv_rank !in: receiver MPI rank (defaults to broadcast)
        integer(INT_MPI), intent(in), optional:: msg_tag   !in: P2P communication tag (defaults to DEFAULT_MPI_TAG)
        integer(INT_MPI), intent(in), optional:: comm_mpi  !in: MPI communicator (defaults to GLOBAL_MPI_COMM)
        logical, intent(in), optional:: sync               !in: if TRUE, the communication will be origin-completed here (default:FALSE)
        integer(INT_MPI):: comm,comm_sz,proc_rank,rx_rank,ctag,errc

        errc=0
        comm=GLOBAL_MPI_COMM; if(present(comm_mpi)) comm=comm_mpi
        if(present(recv_rank)) then !P2P
         rx_rank=recv_rank
         if(comm.eq.GLOBAL_MPI_COMM) then
          comm_sz=impis
         else
          call MPI_Comm_size(comm,comm_sz,errc)
         endif
         if(errc.eq.0) then
          if(rx_rank.lt.0.or.rx_rank.ge.comm_sz) errc=1 !invalid reciever MPI rank
         else
          errc=2 !failed to determine the size of the MPI communicator
         endif
        else !broadcast
         rx_rank=-1
         if(comm.eq.GLOBAL_MPI_COMM) then
          comm_sz=impis; proc_rank=impir
         else
          call MPI_Comm_size(comm,comm_sz,errc)
          call MPI_Comm_rank(comm,proc_rank,errc)
         endif
        endif
        if(errc.eq.0) then
         ctag=DEFAULT_MPI_TAG; if(present(msg_tag)) ctag=msg_tag
         errc=comm_hl%test()
         if(errc.ne.MPI_STAT_PROGRESS_REQ) then
          errc=0
          comm_hl%TimeInitiated=thread_wtime()
          if(rx_rank.ge.0) then !initiate a P2P send
           call send_message(this%Packets,errc)
          else !(initiate) a sending broadcast
           call broadcast_message(this%Packets,errc)
           rx_rank=MPI_ANY_SOURCE
          endif
          if(errc.eq.0) then
           comm_hl%MPIRank=rx_rank !MPI_ANY_SOURCE will mean broadcast
           comm_hl%CommMPI=comm
           comm_hl%CommTag=ctag
           comm_hl%DataCont=>this
           if(present(sync)) then
            if(sync) then
             call comm_hl%wait(errc,ignore_old=.TRUE.)
             if(errc.ne.0) errc=3 !synchronization failed: Communication is lost
            endif
           endif
          else
           call MPI_Request_free(comm_hl%ReqHandle,errc)
           comm_hl%ReqHandle=MPI_REQUEST_NULL
           errc=4 !MPI sending failed
          endif
         else
          errc=5 !communication handle is associated with a message which is currently being communicated
         endif
        endif
        if(present(ierr)) ierr=errc
        return

        contains

         subroutine send_message(msg,jerr)
          integer(ELEM_PACK_SIZE), intent(in):: msg(1:*)
          integer(INT_MPI), intent(out):: jerr
          integer(INT_MPI):: data_typ
          jerr=get_mpi_int_datatype(ELEM_PACK_SIZE,data_typ)
          if(jerr.eq.0) call MPI_Isend(msg,int(this%ffe),data_typ,rx_rank,ctag,comm,comm_hl%ReqHandle,jerr)
          return
         end subroutine send_message

         subroutine broadcast_message(msg,jerr)
          integer(ELEM_PACK_SIZE), intent(in):: msg(1:*)
          integer(INT_MPI), intent(out):: jerr
          integer(INT_MPI):: data_typ
          jerr=get_mpi_int_datatype(ELEM_PACK_SIZE,data_typ)
          if(jerr.eq.0) call MPI_Ibcast(msg,int(this%ffe),data_typ,proc_rank,comm,comm_hl%ReqHandle,jerr)
          return
         end subroutine broadcast_message

        end subroutine PackContSend
!---------------------------------------------------------------------------------------
        subroutine PackContRecv(this,comm_hl,ierr,send_rank,msg_tag,comm_mpi,bcast,sync)
!Posts a receive for a data packet container. If <send_rank> is present,
!the message is expected to come from that MPI process (either P2P or broadcast),
!otherwise it can come from any MPI process (P2P only). If <msg_tag> is present,
!the message is expected to have that MPI tag, otherwise DEFAULT_MPI_TAG will be used.
!If <bcast> is TRUE, a broadcast from <send_rank> is assumed. If <sync> is TRUE,
!the receive will be completed here, otherwise one will need to test
!the communication handle <comm_hl> for completion later.
        implicit none
        class(PackCont_t), intent(inout), target:: this    !in: data packet container
        type(CommHandle_t), intent(inout):: comm_hl        !out: communication handle
        integer(INT_MPI), intent(inout), optional:: ierr   !out: error code (0:success)
        integer(INT_MPI), intent(in), optional:: send_rank !in: sender MPI rank (either P2P or broadcast)
        integer(INT_MPI), intent(in), optional:: msg_tag   !in: P2P communication tag (defaults to MPI_ANY_TAG)
        integer(INT_MPI), intent(in), optional:: comm_mpi  !in: MPI communicator (defaults to GLOBAL_MPI_COMM)
        logical, intent(in), optional:: bcast              !in: if TRUE, a broadcast from <send_rank> is asssumed (default:FALSE)
        logical, intent(in), optional:: sync               !in: if TRUE, the communication will be completed here (default:FALSE)
        integer(INT_MPI):: comm,comm_sz,sx_rank,ctag,buf_vol,errc
        logical:: bcs

        errc=0
        if(present(bcast)) then; bcs=bcast; else; bcs=.FALSE.; endif
        comm=GLOBAL_MPI_COMM; if(present(comm_mpi)) comm=comm_mpi
        if(present(send_rank)) then
         sx_rank=send_rank
         if(comm.eq.GLOBAL_MPI_COMM) then
          comm_sz=impis
         else
          call MPI_Comm_size(comm,comm_sz,errc)
         endif
         if(errc.eq.0) then
          if(sx_rank.lt.0.or.sx_rank.ge.comm_sz) errc=1 !invalid sender MPI rank
         else
          errc=2 !failed to determine the size of the MPI communicator
         endif
        else
         if(.not.bcs) then
          sx_rank=MPI_ANY_SOURCE
         else
          errc=3 !broadcast requires <send_rank>
         endif
        endif
        if(errc.eq.0) then
         ctag=DEFAULT_MPI_TAG; if(present(msg_tag)) ctag=msg_tag
         errc=comm_hl%test()
         if(errc.ne.MPI_STAT_PROGRESS_REQ) then
          buf_vol=size(this%Packets)
          if(buf_vol.gt.1) then
           errc=0
           comm_hl%TimeInitiated=thread_wtime()
           if(bcs) then !(initiate) a receiving broadcast
            call broadcast_message(this%Packets,errc)
           else !initiate a P2P receive
            call receive_message(this%Packets,errc)
           endif
           if(errc.eq.0) then
            comm_hl%MPIRank=sx_rank !sender rank (P2P or BCAST) or MPI_ANY_SOURCE (P2P only)
            comm_hl%CommMPI=comm
            comm_hl%CommTag=ctag
            comm_hl%DataCont=>this
            if(present(sync)) then
             if(sync) then
              call comm_hl%wait(errc,ignore_old=.TRUE.)
              if(errc.ne.0) errc=4 !synchronization failed: Communication is lost
             endif
            endif
           else
            call MPI_Request_free(comm_hl%ReqHandle,errc)
            comm_hl%ReqHandle=MPI_REQUEST_NULL
            errc=5 !MPI sending failed
           endif
          else
           errc=6
          endif
         else
          errc=7 !communication handle is associated with a message which is currently being communicated
         endif
        endif
        if(present(ierr)) ierr=errc
        return

        contains

         subroutine receive_message(msg,jerr)
          integer(ELEM_PACK_SIZE), intent(inout):: msg(1:*)
          integer(INT_MPI), intent(out):: jerr
          integer(INT_MPI):: data_typ
          jerr=get_mpi_int_datatype(ELEM_PACK_SIZE,data_typ)
          if(jerr.eq.0) call MPI_Irecv(msg,buf_vol,data_typ,sx_rank,ctag,comm,comm_hl%ReqHandle,jerr)
          return
         end subroutine receive_message

         subroutine broadcast_message(msg,jerr)
          integer(ELEM_PACK_SIZE), intent(inout):: msg(1:*)
          integer(INT_MPI), intent(out):: jerr
          integer(INT_MPI):: data_typ
          jerr=get_mpi_int_datatype(ELEM_PACK_SIZE,data_typ)
          if(jerr.eq.0) call MPI_Ibcast(msg,buf_vol,data_typ,sx_rank,comm,comm_hl%ReqHandle,jerr)
          return
         end subroutine broadcast_message

        end subroutine PackContRecv
!-----------------------------------------------
        subroutine PackContRegArrived(this,ierr)
!Registers an arrived data packet container (fills in the fields of PackCont_t).
        implicit none
        class(PackCont_t), intent(inout):: this          !inout: data packet container with arrived data
        integer(INT_MPI), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: i,errc
        integer(INT_COUNT):: k,l,m

        errc=0
        if(associated(this%Packets)) then
         this%NumPackets=int(this%Packets(0),INT_MPI)
         if(this%NumPackets.lt.0) then !tagged container
          this%NumPackets=-this%NumPackets
          this%Marked=.TRUE.; k=1
         else !tag-free container
          this%Marked=.FALSE.; k=0
         endif
         this%ffe=1; m=ubound(this%Packets,1)
         do i=1,this%NumPackets
          if(this%ffe+k.le.m) then
           l=packet_full_len(this%Packets(this%ffe+k:))
           if(l.gt.0) then !`this allows empty packets in the container
            this%ffe=this%ffe+k+l
           else
            errc=3; exit
           endif
          else
           errc=2; exit
          endif
         enddo
         if(errc.ne.0) then; this%NumPackets=0; this%ffe=1; endif !mark the data packet container as empty
        else
         errc=1
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine PackContRegArrived

       end module distributed
