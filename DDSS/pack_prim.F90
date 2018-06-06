!Basic object packing/unpacking primitives.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2018/06/05

!Copyright (C) 2014-2017 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2017 Oak Ridge National Laboratory (UT-Battelle)

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

!DESCRIPTION:
! This module provides basic primitives for packing/unpacking objects
! of built-in and derived types for a subsequent communication between
! MPI processes. The main purpose is to provide a necessary API for implementing
! public communication methods for derived types: Object%send() and Object%receive().
! These public communication API may require packing/unpacking of data under the hood.
! Thus, the packing/unpacking primitives provided in this module are intended
! to facilitate the implementation of the .send() and .receive() methods
! in Fortran classes containing a number of (recursive) components of smaller size.
! The size of the components being packed should be O(1). O(N) and larger
! class components should be communicated directly in separate messages.
! Thus, the main purpose of these packing/unpacking primitives is to
! aggregate a number of small objects into a larger plain message.
! Both .send() and .receive() methods are non-blocking and require
! an additional synchronization via either test() or wait() methods.
! Basic abstractions:
!  # Packet (obj_pack_t): Plain byte array which the data is packed into
!    and unpacked from. Normally a packet contains a single object.
!  # Packet envelope (pack_env_t): Container of packets with some
!    additional layout information. An individual packet space is acquired
!    from an existing packet container. Then the data can be packed into
!    the packet space and sealed. Then the packet container (envelope),
!    which contains one or more packets, can be sent to a different MPI process which
!    will be able to unpack any packet from the packet envelope back to an object.
!  # Communication handle (comm_handle_t): Returned by asynchronous .send() and
!    .receive() member procedures, and can later be used for checking
!    the completion of the communication (on both sides separately).
! Behavior:
!  # In order to pack/unpack Fortran built-in and derived types, a packet envelope
!    needs to be created first via the .reserve_mem() member procedure. Then a packet
!    needs to be allocated in the packet envelope via the .acquire_packet() member
!    procedure. The packet can subsequently be used for packing Fortran objects.
!    Once the objects of built-in or derived types have been packed into the packet,
!    the packet needs to be sealed via the .seal_packet() member procedure. If other
!    objects need to be packed, a new packet should be acquired in the packet envelope,
!    filled in, and sealed. Once all objects have been packed into the packet envelope,
!    the latter can be communicated between two MPI processes via the .send() and .receive()
!    member procedures which follow non-blocking semantics. Each of the two procedures
!    returns a communication handle that can be used for checking the completion
!    of the corresponding communication. Once received, any packet from the delivered
!    packet envelope can be unpacked back into the original objects of built-in or derived types.
!    Each packet in the packet envelope has an optional integer tag. Packet envelopes
!    participating in an on-going (non-blocking) communication must not be used for
!    local packing/unpacking until that communication is completed.
!  # It may happen that the free space provided by a packet is insufficient for packing
!    the objects of interest. In this case, the .resize() member procedure needs to be
!    invoked on the packet envelope the packet is part of. The packet envelope will
!    be extended to a user-specified capacity, in turn resulting in a larger buffer
!    in the packet one is currently filling in with data. Similarly, if one
!    exceeds the max number of packets that can be stored in the packet envelope,
!    this resource can also be extended by calling the same member procedure .resize().
!    The packet buffer space overflow occurs during packing data objects into the packet.
!    The max-packets-per-envelope overflow occurs during sealing the packet.
! #  As mentioned, the packing/unpacking primitives are intended for remote cloning
!    of Fortran objects, which is done by value. However, Fortran objects containing
!    non-allocated pointer components associated with local targets cannot be cloned
!    in this way. In order to clone such Fortran objects, the object needs to be
!    constructed in multiple stages on the remote side (builder), based on the
!    information received from the other MPI process. All dependencies have to be
!    constructed, possibly in multiple stages. Then the pointer components can
!    be associated with local targets on the remote process.

       module pack_prim
        use, intrinsic:: ISO_C_BINDING, only: C_PTR,C_INT,C_CHAR,C_NULL_PTR,c_loc,c_f_pointer
        use stsubs, only: size_of
        use dil_basic, only: INTD,INTL
#ifdef USE_MPI_MOD
#ifdef FORTRAN2008
        use mpi_f08      !MPI Fortran 2008 interface `Will not work
#else
        use mpi          !MPI Fortran interface
#endif
        implicit none
        private
#else
        implicit none
        private
        include 'mpif.h' !MPI Fortran interface
#endif
!PARAMETERS:
 !General:
        integer, private:: CONS_OUT=6     !output device
        logical, private:: VERBOSE=.TRUE. !verbosity for errors
        integer, private:: DEBUG=0        !debugging level (0:none)
 !Integers:
        integer, parameter, private:: INT_MPI=INTD !4-byte default integer
 !Error codes:
        integer(INTD), parameter, public:: PACK_SUCCESS=0       !success
        integer(INTD), parameter, public:: PACK_ERROR=-666      !generic error
        integer(INTD), parameter, public:: PACK_NULL=-1         !null object
        integer(INTD), parameter, public:: PACK_OVERFLOW=-2     !overflow
        integer(INTD), parameter, public:: PACK_INVALID_ARGS=-3 !invalid arguments
        integer(INTD), parameter, public:: PACK_ALLOC_FAILED=-4 !memory allocation failed
        integer(INTD), parameter, public:: PACK_FREE_FAILED=-5  !memory deallocation failed
        integer(INTD), parameter, public:: PACK_BUSY=-6         !object is in use by others
        integer(INTD), parameter, public:: PACK_IDLE=-7         !object is idle (not in use)
        integer(INTD), parameter, public:: PACK_MPI_ERR=-8      !MPI communication error
 !Packet:
        integer(INTL), parameter, public:: PACK_BASE=1_INTL     !first offset in a packet is 1 (not 0)
 !Packet envelope configuration:
        integer(INTD), parameter, private:: DEFAULT_MAX_PACKETS=1024       !default max number of packets per envelope
        integer(INTD), parameter, private:: DEFAULT_AVERAGE_PACK_SIZE=1024 !default average packet volume in bytes
        integer(INTL), parameter, private:: DEFAULT_ENVELOPE_CAPACITY=DEFAULT_AVERAGE_PACK_SIZE*DEFAULT_MAX_PACKETS !default envelope capacity
        integer(INTL), parameter, private:: DEFAULT_PACKET_TAG=0 !default packet tag
 !MPI:
        integer(INT_MPI), parameter, private:: DEFAULT_MPI_TAG=0
!TYPES:
 !Packet (local):
        type, public:: obj_pack_t
         integer(INTL), private:: length=0_INTL !used length of the packet buffer (bytes)
         character(C_CHAR), pointer, contiguous, private:: buffer(:)=>NULL() !packet buffer
         integer(INTL), private:: offset=0_INTL !current offset during unpacking
         contains
          procedure, private:: construct=>ObjPackConstruct     !packet constructor (internal)
          procedure, public:: clean=>ObjPackClean              !packet cleaner
          procedure, public:: is_healthy=>ObjPackIsHealthy     !returns TRUE if the packet is healthy (consistent)
          procedure, public:: get_capacity=>ObjPackGetCapacity !returns the capacity of the packet buffer in bytes
          procedure, public:: get_length=>ObjPackGetLength     !returns the current length of the packet in bytes
          procedure, public:: has_room=>ObjPackHasRoom         !TRUE means one can add data to the packet, FALSE otherwise
          procedure, public:: space_left=>ObjPackSpaceLeft     !returns the amount of free space left in the packet buffer in bytes
          procedure, public:: reset=>ObjPackReset              !resets the unpacking offset to the beginning of the packet
        end type obj_pack_t
 !Packet envelope (communicable):
        type, public:: pack_env_t
         integer(INTL), private:: length=0_INTL       !in-use length of the packet envelope (bytes)
         integer(INTD), private:: num_packets=0       !number of packets in the packet envelope
         integer(INTD), private:: tag=DEFAULT_MPI_TAG !optional tag
         class(obj_pack_t), pointer, private:: curr_packet=>NULL() !current packet (set when in-use)
         logical, private:: busy=.FALSE.              !.TRUE. when there is an active packet being filled in (in-use flag)
         integer(INTL), pointer, contiguous, private:: pack_offset(:)=>NULL() !offset of each packet in the envelope (byte)
         integer(INTL), pointer, contiguous, private:: pack_len(:)=>NULL()    !length of each packet present in the envelope (bytes)
         integer(INTL), pointer, contiguous, private:: pack_tag(:)=>NULL()    !tag for each packet present in the envelope
         character(C_CHAR), pointer, contiguous, private:: buffer(:)=>NULL()  !buffer
         contains
          procedure, private:: resize=>PackEnvResize                !resize either the packet storage buffer or the tables or both
          procedure, public:: get_capacity=>PackEnvGetCapacity      !get the current buffer capacity (bytes)
          procedure, public:: get_length=>PackEnvGetLength          !get the current used length of the buffer (bytes)
          procedure, public:: get_max_packets=>PackEnvGetMaxPackets !get the max limit on the amount of packets in the envelope
          procedure, public:: get_num_packets=>PackEnvGetNumPackets !get the current number of packets in the envelope
          procedure, public:: get_tag=>PackEnvGetTag                !get the packet envelope tag
          procedure, public:: set_tag=>PackEnvSetTag                !set the packet envelope tag
          procedure, public:: is_busy=>PackEnvIsBusy                !check whether the packet envelope is currently in use
          procedure, public:: is_healthy=>PackEnvIsHealthy          !check whether the object is healthy (consistent)
          procedure, public:: reserve_mem=>PackEnvReserveMem        !reserve memory for the data buffer and/or layout tables
          procedure, public:: clean=>PackEnvClean                   !clean the packet envelope without releasing the memory
          procedure, public:: destroy=>PackEnvDestroy               !destroy the packet envelope completely
          procedure, public:: acquire_packet=>PackEnvAcquirePacket  !acquire a packet in the packet envelope
          procedure, public:: discard_packet=>PackEnvDiscardPacket  !discard a packet (either finished or unfinished)
          procedure, public:: seal_packet=>PackEnvSealPacket        !seal a packet (finalize)
          procedure, public:: extract_packet=>PackEnvExtractPacket  !extract a packet from the packet envelope
          procedure, public:: send=>PackEnvSend                     !send the packet envelope to another process
          procedure, public:: receive=>PackEnvReceive               !receive a packet envelope from another process
          procedure, private:: decode_mpi_msg=>PackEnvDecodeMPIMsg  !decodes the incoming MPI message back into the object
        end type pack_env_t
 !Data communication handle:
        type, public:: comm_handle_t
         logical, private:: active=.FALSE.                  !switches to .TRUE. once the communication handle is associated with a message
         integer(INT_MPI), private:: req=MPI_REQUEST_NULL   !MPI request handle
         integer(INT_MPI), private:: comm=MPI_COMM_NULL     !MPI communicator
         integer(INT_MPI), private:: stat(MPI_STATUS_SIZE)  !MPI status
         class(pack_env_t), pointer:: recv_pack_env=>NULL() !receive operation pointer to the packet envelope being received
         contains
          procedure, private:: construct=>CommHandleConstruct !construct the communication handle (internal)
          procedure, public:: clean=>CommHandleClean          !clean the communication handle
          procedure, public:: is_clean=>CommHandleIsClean     !tests whether the communication handle is clean
          procedure, public:: is_active=>CommHandleIsActive   !tests whether the communication handle is associated with an active message
          procedure, public:: test=>CommHandleTest            !tests the completion of the communication
          procedure, public:: wait=>CommHandleWait            !waits upon the completion of the communication
          procedure, private:: attach_pack_env=>CommHandleAttachPackEnv !attaches the parental packet envelope (for receive operations)
          procedure, public:: print_it=>CommHandlePrintIt     !prints the communication handle info
        end type comm_handle_t
!INTERFACES:
 !Packing for built-in types:
        interface pack_builtin
         module procedure pack_integer1
         module procedure pack_integer2
         module procedure pack_integer4
         module procedure pack_integer4_arr1
         module procedure pack_integer8
         module procedure pack_logical
         module procedure pack_real4
         module procedure pack_real8
         module procedure pack_complex4
         module procedure pack_complex8
         module procedure pack_string
        end interface pack_builtin
        public pack_builtin
        public pack_string
 !Unpacking for built-in types:
        interface unpack_builtin
         module procedure unpack_integer1
         module procedure unpack_integer2
         module procedure unpack_integer4
         module procedure unpack_integer4_arr1
         module procedure unpack_integer8
         module procedure unpack_logical
         module procedure unpack_real4
         module procedure unpack_real8
         module procedure unpack_complex4
         module procedure unpack_complex8
         module procedure unpack_string
        end interface unpack_builtin
        public unpack_builtin
        public unpack_string

       contains
!DEFINITION:
![obj_pack_t]============================================
        subroutine ObjPackConstruct(this,buf,ierr,length)
!Constructs a packet (either empty or filled in).
         implicit none
         class(obj_pack_t), intent(inout):: this         !inout: packet (in: empty packet, out: allocated packet)
         character(C_CHAR), target, contiguous:: buf(1:) !in: external buffer space
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTL), intent(in), optional:: length    !if present, the packet will be given this legnth (pre-existing data buffer)
         integer(INTD):: errc
         integer(INTL):: bs

         errc=PACK_SUCCESS
         if(this%get_capacity().le.0) then !empty packet
          bs=size(buf)
          if(bs.gt.0) then
           this%buffer(1:)=>buf(:); this%length=0_INTL; this%offset=0_INTL
           if(present(length)) then !non-empty packet constructor (empty if <length> = 0)
            if(length.ge.0.and.length.le.bs) then
             this%length=length
            else
             errc=PACK_INVALID_ARGS
            endif
           endif
          else
           errc=PACK_INVALID_ARGS
          endif
         else
          errc=PACK_BUSY !non-empty packet
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ObjPackConstruct
!-----------------------------------------
        subroutine ObjPackClean(this,ierr)
!Cleans a packet (resets it to empty).
         implicit none
         class(obj_pack_t), intent(inout):: this     !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=PACK_SUCCESS
         this%length=0_INTL; this%buffer=>NULL(); this%offset=0_INTL
         if(present(ierr)) ierr=errc
         return
        end subroutine ObjPackClean
!-------------------------------------------------------
        function ObjPackIsHealthy(this,ierr) result(res)
!Returns TRUE if the packet is healthy (consistent).
!Null packets are not considered healthy: <ierr> = PACK_NULL.
         implicit none
         logical:: res                               !out: result
         class(obj_pack_t), intent(in):: this        !in: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=PACK_SUCCESS; res=.FALSE.
         if(associated(this%buffer)) then
          if(this%get_capacity().ge.this%get_length()) then
           if(this%offset.ge.0_INTL.and.this%offset.le.this%length) then
            res=.TRUE.
           else
            errc=PACK_ERROR
           endif
          else
           errc=PACK_ERROR
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end function ObjPackIsHealthy
!-----------------------------------------------------------
        function ObjPackGetCapacity(this,ierr) result(bytes)
!Returns the packet capacity (max length) in bytes.
         implicit none
         integer(INTL):: bytes                       !out: packet buffer capacity
         class(obj_pack_t), intent(in):: this        !in: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=PACK_SUCCESS; bytes=0_INTL
         if(associated(this%buffer)) then
          bytes=size(this%buffer)
         endif
         if(present(ierr)) ierr=errc
         return
        end function ObjPackGetCapacity
!---------------------------------------------------------
        function ObjPackGetLength(this,ierr) result(bytes)
!Returns the current length of the packet in bytes.
         implicit none
         integer(INTL):: bytes                       !out: current packet length
         class(obj_pack_t), intent(in):: this        !in: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=PACK_SUCCESS; bytes=0_INTL
         if(this%get_capacity().gt.0) then
          bytes=this%length
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end function ObjPackGetLength
!------------------------------------------------------
        function ObjPackHasRoom(this,ierr) result(answ)
!Returns .TRUE. if the packet can be used for packing data, .FALSE. otherwise.
!A packet can be used for packing data if it has been assigned an external
!buffer during construction.
         implicit none
         logical:: answ                              !out: answer
         class(obj_pack_t), intent(in):: this        !in: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer(INTL):: cap

         answ=.TRUE.; cap=this%get_capacity(errc)
         if(cap.le.0_INTL.or.errc.ne.PACK_SUCCESS) answ=.FALSE.
         if(this%length.ge.cap) then
          if(this%length.gt.cap) errc=PACK_ERROR !corrupted
          answ=.FALSE.
         endif
         if(present(ierr)) ierr=errc
         return
        end function ObjPackHasRoom
!---------------------------------------------------------
        function ObjPackSpaceLeft(this,ierr) result(bytes)
!Returns the number of bytes left in the packet buffer.
!A negative return means that the packet has not been constructed yet.
         implicit none
         integer(INTL):: bytes                       !out: number of bytes left in the packet buffer
         class(obj_pack_t), intent(in):: this        !in: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer(INTL):: cap

         bytes=-1_INTL; cap=this%get_capacity(errc)
         if(cap.gt.0.and.errc.eq.PACK_SUCCESS) then
          bytes=cap-this%length
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end function ObjPackSpaceLeft
!-----------------------------------------
        subroutine ObjPackReset(this,ierr)
!Resets the unpacking offset to the beginning of the packet.
         implicit none
         class(obj_pack_t), intent(inout):: this     !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code

         this%offset=0_INTL
         if(present(ierr)) ierr=PACK_SUCCESS
         return
        end subroutine ObjPackReset
![pack_env_t]===================================================
        subroutine PackEnvResize(this,ierr,buf_size,max_packets)
!Resizes either the packet buffer or the packet layout tables or both.
!The buffer is used for storing data. The layout tables contain
!bookeeping information and the size of these tables limits
!the max number of packets that can be stored in the packet envelope.
!The packet envelope passed into this procedure is allowed to be in
!the in-use state (with an active open packet). If neither <buf_size>
!nor <max_packets> is passed here, default values will be used for both.
         implicit none
         class(pack_env_t), intent(inout):: this           !inout: packet envelope
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTL), intent(in), optional:: buf_size    !in: new buffer size in bytes
         integer(INTD), intent(in), optional:: max_packets !in: new max limit on the number of packets stored
         integer:: jerr
         integer(INTD):: mpk,errc
         integer(INTL):: fl,bfs
         character(C_CHAR), pointer, contiguous:: chp(:)
         integer(INTL), pointer, contiguous:: i8p(:)

         if(this%is_healthy(errc)) then
          fl=this%length
          if(associated(this%curr_packet)) fl=fl+this%curr_packet%get_length()
          if(present(buf_size)) then
           bfs=buf_size
           if(present(max_packets)) then
            mpk=max_packets
           else
            mpk=-1
           endif
          else
           if(present(max_packets)) then
            mpk=max_packets; bfs=-1_INTL
           else
            mpk=this%get_max_packets()+DEFAULT_MAX_PACKETS
            bfs=this%get_capacity()+DEFAULT_ENVELOPE_CAPACITY
           endif
          endif
 !Buffer:
          if(bfs.ge.0) then
           if(fl.le.bfs) then
            allocate(chp(1:bfs),STAT=jerr)
            if(jerr.eq.0) then
             chp(1:fl)=this%buffer(1:fl)
             if(associated(this%curr_packet)) this%curr_packet%buffer(1:)=>chp(this%length+1_INTL:)
             deallocate(this%buffer); this%buffer=>chp; chp=>NULL()
            else
             errc=PACK_ALLOC_FAILED
            endif
           else
            errc=PACK_INVALID_ARGS
           endif
          endif
 !Layout tables:
          if(mpk.ge.0) then
           if(this%num_packets.le.mpk) then
            allocate(i8p(1:mpk),STAT=jerr)
            if(jerr.eq.0) then
             i8p(1:this%num_packets)=this%pack_offset(1:this%num_packets)
             deallocate(this%pack_offset); this%pack_offset=>i8p; i8p=>NULL()
             allocate(i8p(1:mpk),STAT=jerr)
             if(jerr.eq.0) then
              i8p(1:this%num_packets)=this%pack_len(1:this%num_packets)
              deallocate(this%pack_len); this%pack_len=>i8p; i8p=>NULL()
              allocate(i8p(1:mpk),STAT=jerr)
              if(jerr.eq.0) then
               i8p(1:this%num_packets)=this%pack_tag(1:this%num_packets)
               deallocate(this%pack_tag); this%pack_tag=>i8p; i8p=>NULL()
              else
               errc=PACK_ALLOC_FAILED
              endif
             else
              errc=PACK_ALLOC_FAILED
             endif
            else
             errc=PACK_ALLOC_FAILED
            endif
           else
            errc=PACK_INVALID_ARGS
           endif
          endif
         else
          errc=PACK_ERROR
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine PackEnvResize
!-----------------------------------------------------------
        function PackEnvGetCapacity(this,ierr) result(bytes)
!Returns the capacity of the packet envelope in bytes.
         implicit none
         integer(INTL):: bytes                       !out: packet envelope capacity in bytes
         class(pack_env_t), intent(in):: this        !in: packet envelope
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=PACK_SUCCESS; bytes=0_INTL
         if(associated(this%buffer)) bytes=size(this%buffer)
         if(present(ierr)) ierr=errc
         return
        end function PackEnvGetCapacity
!---------------------------------------------------------
        function PackEnvGetLength(this,ierr) result(bytes)
!Returns the current length of the packet envelope in bytes.
         implicit none
         integer(INTL):: bytes                       !out: envelope length in bytes
         class(pack_env_t), intent(in):: this        !in: envelope
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=PACK_SUCCESS; bytes=this%length
         if(bytes.gt.this%get_capacity()) then; errc=PACK_ERROR; bytes=-1_INTL; endif !trap
         if(present(ierr)) ierr=errc
         return
        end function PackEnvGetLength
!-----------------------------------------------------------
        function PackEnvGetMaxPackets(this,ierr) result(num)
!Returns the max number of packets assumed for the packet envelope.
         implicit none
         integer(INTD):: num                         !out: max number of packets
         class(pack_env_t), intent(in):: this        !in: packet envelope
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         logical:: as1,as2,as3

         errc=PACK_SUCCESS; num=0
         as1=associated(this%pack_offset)
         as2=associated(this%pack_len)
         as3=associated(this%pack_tag)
         if(as1.and.as2.and.as3) then
          num=size(this%pack_offset)
          if(size(this%pack_len).ne.num.or.size(this%pack_tag).ne.num) then; errc=PACK_ERROR; num=-1; endif !trap
         else
          if(as1.or.as2.or.as3) then; errc=PACK_ERROR; num=-1; endif
         endif
         if(present(ierr)) ierr=errc
         return
        end function PackEnvGetMaxPackets
!-----------------------------------------------------------
        function PackEnvGetNumPackets(this,ierr) result(num)
!Returns the current number of packets sealed in the packet envelope.
         implicit none
         integer(INTD):: num                         !out: number of sealed packets in the envelope
         class(pack_env_t), intent(in):: this        !in: packet envelope
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=PACK_SUCCESS; num=this%num_packets
         if(num.gt.this%get_max_packets(errc)) then; errc=PACK_ERROR; num=-1; endif !trap
         if(present(ierr)) ierr=errc
         return
        end function PackEnvGetNumPackets
!----------------------------------------------------
        function PackEnvGetTag(this,ierr) result(tag)
!Returns the packet envelope tag.
         implicit none
         integer(INTD):: tag                         !out: tag
         class(pack_env_t), intent(in):: this        !in: packet envelope
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=PACK_SUCCESS
         tag=this%tag
         if(present(ierr)) ierr=errc
         return
        end function PackEnvGetTag
!----------------------------------------------
        subroutine PackEnvSetTag(this,tag,ierr)
!Sets the packet envelope tag.
         implicit none
         class(pack_env_t), intent(inout):: this     !inout: packet envelope
         integer(INTD), intent(in):: tag             !in: tag
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=PACK_SUCCESS
         this%tag=tag
         if(present(ierr)) ierr=errc
         return
        end subroutine PackEnvSetTag
!-----------------------------------------------------
        function PackEnvIsBusy(this,ierr) result(answ)
!Returns .TRUE. if the envelope is blocked by an active packet.
!An active packet is blocking the packet envelope until the packet
!is sealed or discarded.
         implicit none
         logical:: answ                              !out: answer
         class(pack_env_t), intent(in):: this        !in: packet envelope
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         answ=.TRUE.
         if(this%is_healthy(errc)) then
          if(errc.eq.PACK_SUCCESS) answ=this%busy
         else
          errc=PACK_ERROR
         endif
         if(present(ierr)) ierr=errc
         return
        end function PackEnvIsBusy
!--------------------------------------------------------
        function PackEnvIsHealthy(this,ierr) result(answ)
!Returns .TRUE. if the object is healthy, .FALSE. otherwise (corrupted or null).
!Null packet envelopes are not considered healthy (<ierr> = PACK_NULL).
         implicit none
         logical:: answ                              !out: answer
         class(pack_env_t), intent(in):: this        !in: packet envelope
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         answ=.FALSE.; errc=PACK_SUCCESS
         if(this%get_capacity().gt.0) then
          if(this%get_capacity().lt.this%get_length()) errc=PACK_ERROR
          if(this%get_max_packets().lt.this%get_num_packets()) errc=PACK_ERROR
          if((this%busy.and.(.not.associated(this%curr_packet))).or.&
            &((.not.this%busy).and.associated(this%curr_packet))) errc=PACK_ERROR
          if(errc.eq.PACK_SUCCESS) answ=.TRUE.
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end function PackEnvIsHealthy
!-------------------------------------------------------------------------------
        subroutine PackEnvReserveMem(this,ierr,mem_size,max_packets,ignore_less)
!Reserves memory for a packet envelope. If the packet envelope has its memory
!already reserved, it will be re-reserved if the new amount of memory is larger
!than the original one, otherwise an error PACK_INVALID_ARGS will be returned
!without any action (supressed by <ignore_less> = .TRUE.). The same policy applies
!to <max_packets>. If <mem_size> or <max_packets> are absent, the default values
!DEFAULT_ENVELOPE_CAPACITY and DEFAULT_MAX_PACKETS will be used, respectively.
         implicit none
         class(pack_env_t), intent(inout):: this           !inout: packet envelope
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTL), intent(in), optional:: mem_size    !in: memory size in bytes to reserve for the buffer
         integer(INTD), intent(in), optional:: max_packets !in: max number of packets to assume
         logical, intent(in), optional:: ignore_less       !in: if .TRUE., no error will be raised if new < old size
         integer(INTD):: errc
         integer(INTL), pointer, contiguous:: i8p(:),j8p(:),k8p(:)
         character(C_CHAR), pointer, contiguous:: chp(:)
         integer(INTL):: buf_sz
         integer(INTD):: max_pkts
         integer:: jerr
         logical:: igless

         if(present(ignore_less)) then; igless=ignore_less; else; igless=.FALSE.; endif
         if(this%get_capacity(errc).gt.0) then
          if(errc.eq.PACK_SUCCESS) then
           if(.not.this%is_healthy(errc)) errc=PACK_ERROR
          endif
         endif
         if(errc.eq.PACK_SUCCESS) then
          if(.not.this%busy) then
           if(present(mem_size)) then
            if(mem_size.le.0) errc=PACK_INVALID_ARGS
            buf_sz=mem_size
           else
            buf_sz=DEFAULT_ENVELOPE_CAPACITY
           endif
           if(present(max_packets)) then
            if(max_packets.le.0) errc=PACK_INVALID_ARGS
            max_pkts=max_packets
           else
            max_pkts=DEFAULT_MAX_PACKETS
           endif
           if(errc.eq.PACK_SUCCESS) then
 !Buffer:
            if(buf_sz.gt.this%get_capacity()) then
             allocate(chp(1:buf_sz),STAT=jerr)
             if(jerr.eq.0) then
              if(this%get_capacity().gt.0) then
               chp(1:this%length)=this%buffer(1:this%length)
               deallocate(this%buffer)
              endif
              this%buffer=>chp; chp=>NULL()
             else
              errc=PACK_ALLOC_FAILED
             endif
            elseif((.not.igless).and.buf_sz.lt.this%get_capacity()) then
             errc=PACK_INVALID_ARGS
            endif
 !Tables:
            if(max_pkts.gt.this%get_max_packets(errc)) then
             if(errc.eq.PACK_SUCCESS) then
              allocate(i8p(1:max_pkts),STAT=jerr)
              if(jerr.eq.0) then
               allocate(j8p(1:max_pkts),STAT=jerr)
               if(jerr.eq.0) then
                allocate(k8p(1:max_pkts),STAT=jerr)
                if(jerr.eq.0) then
                 if(this%get_max_packets().gt.0) then
                  i8p(1:this%num_packets)=this%pack_offset(1:this%num_packets)
                  j8p(1:this%num_packets)=this%pack_len(1:this%num_packets)
                  k8p(1:this%num_packets)=this%pack_tag(1:this%num_packets)
                  deallocate(this%pack_tag); deallocate(this%pack_len); deallocate(this%pack_offset)
                 endif
                 this%pack_offset=>i8p; this%pack_len=>j8p; this%pack_tag=>k8p
                 i8p=>NULL(); j8p=>NULL(); k8p=>NULL()
                else
                 deallocate(j8p); deallocate(i8p); errc=PACK_ALLOC_FAILED
                endif
               else
                deallocate(i8p); errc=PACK_ALLOC_FAILED
               endif
              else
               errc=PACK_ALLOC_FAILED
              endif
             endif
            elseif((.not.igless).and.max_pkts.lt.this%get_max_packets(errc)) then
             if(errc.eq.PACK_SUCCESS) errc=PACK_INVALID_ARGS
            endif
           endif
          else
           errc=PACK_BUSY
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine PackEnvReserveMem
!-----------------------------------------
        subroutine PackEnvClean(this,ierr)
!Cleans the packet envelope without deallocating the memory.
         implicit none
         class(pack_env_t), intent(inout):: this     !inout: packet envelope
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=PACK_SUCCESS
         if(.not.this%busy) then
          this%length=0_INTL
          this%num_packets=0
          this%tag=DEFAULT_MPI_TAG
          this%busy=.FALSE.
          if(associated(this%curr_packet)) then
           call this%curr_packet%clean()
           this%curr_packet=>NULL()
          endif
         else
          errc=PACK_BUSY
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine PackEnvClean
!-------------------------------------------
        subroutine PackEnvDestroy(this,ierr)
!Destroys the packet envelope.
         implicit none
         class(pack_env_t), intent(inout):: this     !inout: packet envelope
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(.not.this%busy) then
          if(associated(this%pack_tag)) deallocate(this%pack_tag)
          if(associated(this%pack_len)) deallocate(this%pack_len)
          if(associated(this%pack_offset)) deallocate(this%pack_offset)
          if(associated(this%buffer)) deallocate(this%buffer)
          call this%clean(errc)
         else
          errc=PACK_BUSY
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine PackEnvDestroy
!--------------------------------------------------------------
        subroutine PackEnvAcquirePacket(this,pkt,ierr,preclean)
!Acquires a packet space in the packet envelope and returns
!a packet object. The packet object must be clean on entrance,
!otherwise a status PACK_BUSY will be returned. However, if
!<preclean>=TRUE, the packet will be precleaned on entrance.
         implicit none
         class(pack_env_t), intent(inout):: this        !inout: packet envelope
         class(obj_pack_t), target, intent(inout):: pkt !inout: in:clean (empty) packet, out:active empty packet
         integer(INTD), intent(out), optional:: ierr    !out: error code
         logical, intent(in), optional:: preclean       !in: if TRUE the packet will be cleaned here before use
         integer(INTD):: errc
         integer(INTL):: cap
         character(C_CHAR), pointer, contiguous:: chp(:)

         cap=this%get_capacity(errc)
         if(cap.gt.0.and.errc.eq.PACK_SUCCESS) then
          if(this%is_healthy(errc)) then
           if(errc.eq.PACK_SUCCESS) then
            if(.not.this%busy) then
             if(this%length.lt.cap) then
              if(present(preclean)) then; if(preclean) call pkt%clean(); endif !preclean the packet if asked for
              chp(1:cap-this%length)=>this%buffer(this%length+1_INTL:cap)
              call pkt%construct(chp,errc) !construct the packet
              if(errc.eq.PACK_SUCCESS) then
               this%curr_packet=>pkt !associate the packet with the packet envelope
               this%busy=.TRUE. !mark the packet envelope as busy
              else
               call pkt%clean()
              endif
              chp=>NULL()
             else
              errc=PACK_OVERFLOW
             endif
            else
             errc=PACK_BUSY
            endif
           endif
          else
           errc=PACK_ERROR
          endif
         else
          if(errc.eq.PACK_SUCCESS) errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine PackEnvAcquirePacket
!----------------------------------------------------------
        subroutine PackEnvDiscardPacket(this,ierr,pack_num)
!Discards a packet from the packet envelope. If <pack_num> is present,
!that specific (sealed) packet will be discarded. Otherwise, the currently
!active (unfinished) packet, if any, will be discarded. If <pack_num>
!is present, the packet envelope shall not be in the in-use (busy) state.
!If <pack_num> is absent, the packet envelope must be in the in-use state.
         implicit none
         class(pack_env_t), intent(inout):: this        !inout: packet envelope
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTL), intent(in), optional:: pack_num !in: packet number [1:MAX]
         integer(INTD):: i,errc
         integer(INTL):: bg,ln,l,cap

         cap=this%get_capacity(errc)
         if(cap.gt.0.and.errc.eq.PACK_SUCCESS) then
          if(this%is_healthy(errc)) then
           if(errc.eq.PACK_SUCCESS) then
            if(present(pack_num)) then !discard a sealed packet
             if(.not.this%busy) then
              if(pack_num.ge.1.and.pack_num.le.this%num_packets) then
               bg=this%pack_offset(pack_num); ln=this%pack_len(pack_num)
               do l=bg,this%length-ln; this%buffer(l)=this%buffer(l+ln); enddo
               this%length=this%length-ln
               do i=pack_num,this%num_packets-1
                this%pack_offset(i)=this%pack_offset(i+1)-ln
                this%pack_len(i)=this%pack_len(i+1)
                this%pack_tag(i)=this%pack_tag(i+1)
               enddo
               this%num_packets=this%num_packets-1
              else
               errc=PACK_INVALID_ARGS
              endif
             else
              errc=PACK_BUSY
             endif
            else !discard an unfinished (current) packet
             if(this%busy) then
              if(associated(this%curr_packet)) then
               call this%curr_packet%clean(errc); this%curr_packet=>NULL()
               this%busy=.FALSE.
              else
               errc=PACK_ERROR
              endif
             else
              errc=PACK_IDLE
             endif
            endif
           endif
          else
           errc=PACK_ERROR
          endif
         else
          if(errc.eq.PACK_SUCCESS) errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine PackEnvDiscardPacket
!--------------------------------------------------
        subroutine PackEnvSealPacket(this,ierr,tag)
!Seals a packet in the packet envelope and releases the busy flag.
!The packet object will be cleaned automatically at the end.
         implicit none
         class(pack_env_t), intent(inout):: this     !inout: packet envelope
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTL), intent(in), optional:: tag   !in: packet tag
         integer(INTD):: n,errc
         integer(INTL):: l,cap

         cap=this%get_capacity(errc)
         if(cap.gt.0.and.errc.eq.PACK_SUCCESS) then
          if(this%is_healthy(errc)) then
           if(errc.eq.PACK_SUCCESS) then
            if(this%busy) then
             if(this%length.le.cap) then
              if(associated(this%curr_packet)) then
               l=this%curr_packet%get_length(errc)
               if(l.gt.0.and.errc.eq.PACK_SUCCESS) then !packet contains information
                n=this%get_max_packets(errc)
                if(this%num_packets.lt.n.and.errc.eq.PACK_SUCCESS) then
                 this%num_packets=this%num_packets+1
                 this%pack_offset(this%num_packets)=this%length+1_INTL
                 this%pack_len(this%num_packets)=l
                 if(present(tag)) then
                  this%pack_tag(this%num_packets)=tag
                 else
                  this%pack_tag(this%num_packets)=DEFAULT_PACKET_TAG
                 endif
                 this%length=this%length+l
                 call this%curr_packet%clean(); this%curr_packet=>NULL()
                 this%busy=.FALSE.
                else
                 if(errc.eq.PACK_SUCCESS) errc=PACK_OVERFLOW
                endif
               else
                if(errc.eq.PACK_SUCCESS) errc=PACK_NULL
               endif
              else
               errc=PACK_ERROR
              endif
             else
              errc=PACK_OVERFLOW
             endif
            else
             errc=PACK_IDLE
            endif
           endif
          else
           errc=PACK_ERROR
          endif
         else
          if(errc.eq.PACK_SUCCESS) errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine PackEnvSealPacket
!--------------------------------------------------------------------------
        subroutine PackEnvExtractPacket(this,pkt_num,pkt,ierr,tag,preclean)
!Extracts a packet from a packet envelope. This means that
!an object of class(obj_pack_t) will be returned and its buffer
!will be associated with the corresponding part of the packet
!envelope buffer. That is, the packet is returned by reference
!and the packet data still resides in the packet envelope.
!The packet must be clean on entrance. Since the extraction
!is done by reference, a packet can be extracted even if the
!packet envelope is marked as in-use (an unsealed active
!packet is being currently filled in).
         implicit none
         class(pack_env_t), intent(in):: this        !in: packet envelope
         integer(INTD), intent(in):: pkt_num         !in: packet number [1..MAX]
         class(obj_pack_t), intent(inout):: pkt      !inout: packet (must be clean on entrance, unless <preclean> = TRUE)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTL), intent(out), optional:: tag  !out: packet tag
         logical, intent(in), optional:: preclean    !in: if TRUE the packet will be forcefully cleaned on entrance
         integer(INTD):: errc
         integer(INTL):: cap,bg,ln
         character(C_CHAR), pointer, contiguous:: buf_p(:)=>NULL()

         cap=this%get_capacity(errc)
         if(cap.gt.0.and.errc.eq.PACK_SUCCESS) then
          if(this%is_healthy(errc)) then
           if(errc.eq.PACK_SUCCESS) then
            if(pkt_num.ge.1.and.pkt_num.le.this%get_num_packets()) then
             if(present(preclean)) then; if(preclean) call pkt%clean(); endif
             if(present(tag)) tag=this%pack_tag(pkt_num)
             bg=this%pack_offset(pkt_num); ln=this%pack_len(pkt_num)
             buf_p=>this%buffer(bg:bg+ln-1_INTL)
             call pkt%construct(buf_p,errc,ln)
             if(errc.ne.PACK_SUCCESS) call pkt%clean()
            else
             if(VERBOSE) write(CONS_OUT,'("#ERROR(pack_env_t.extract_packet): Invalid packet number: ",i9,1x,i9)')&
                         &pkt_num,this%get_num_packets()
             errc=PACK_INVALID_ARGS
            endif
           endif
          else
           errc=PACK_ERROR
          endif
         else
          if(errc.eq.PACK_SUCCESS) errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine PackEnvExtractPacket
!-----------------------------------------------------------------------
        subroutine PackEnvSend(this,proc_rank,comm_handle,ierr,tag,comm)
!Initiates a send of the packet envelope to the specified process. The packet
!envelope must not contain an unsealed (active) packet. The MPI message
!will contain the packet envelope data buffer as a plain CHARACTER array
!succeeded by the data layout tables appended to the tail. The last
!three 8-byte fields are: Length of the data buffer, Number of packets, Tag.
!The communication handle must not be active on entrance.
!This is a NON-BLOCKING method!
         implicit none
         class(pack_env_t), intent(inout):: this           !in: packet envelope
         integer(INT_MPI), intent(in):: proc_rank          !in: process rank
         class(comm_handle_t), intent(inout):: comm_handle !out: communication handle (must not be active on entrance)
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INT_MPI), intent(in), optional:: tag      !in: MPI message tag
         integer(INT_MPI), intent(in), optional:: comm     !in: MPI communicator (defaults to MPI_COMM_WORLD)
         integer(INTD):: errc
         integer(INTL):: cap,fl
         integer(INT_MPI):: rk,tg,cm,rh

         cap=this%get_capacity(errc)
         if(cap.gt.0.and.errc.eq.PACK_SUCCESS) then
          if(this%is_healthy(errc)) then
           if(errc.eq.PACK_SUCCESS) then
            if(.not.this%busy) then
             if(this%get_length().gt.0.and.this%get_num_packets().gt.0) then
              if(proc_rank.ge.0) then
               rk=proc_rank
               if(.not.comm_handle%is_active(errc)) then
                if(errc.eq.PACK_SUCCESS) then
                 call comm_handle%clean(errc) !clean the communication handle before usage
                 if(errc.eq.PACK_SUCCESS) then
                  if(present(tag)) then; tg=tag; else; tg=DEFAULT_MPI_TAG; endif
                  if(present(comm)) then; cm=comm; else; cm=MPI_COMM_WORLD; endif
                  call pack_layout_tables(fl,errc)
                  if(errc.eq.PACK_SUCCESS) then
                   call send_mpi_message(cm,rk,fl,this%buffer,tg,rh,errc)
                   if(errc.eq.PACK_SUCCESS) call comm_handle%construct(cm,rh,errc)
                  endif
                 endif
                endif
               else
                errc=PACK_INVALID_ARGS !communication handle is active
               endif
              else
               errc=PACK_INVALID_ARGS
              endif
             else
              errc=PACK_NULL
             endif
            else
             errc=PACK_BUSY
            endif
           endif
          else
           errc=PACK_ERROR
          endif
         else
          if(errc.eq.PACK_SUCCESS) errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return

         contains

          subroutine pack_layout_tables(jl,jerr)
           implicit none
           integer(INTL), intent(out):: jl
           integer(INTD), intent(out):: jerr
           integer(INTD):: js
           type(C_PTR):: cptr
           integer(INTL), pointer:: j8p
           character(C_CHAR), pointer, contiguous:: fptr(:)

           jerr=PACK_SUCCESS; jl=0_INTL
           js=size_of(jl) !size of integer(INTL) in bytes
           jl=this%length+js*int(this%num_packets*3,INTL) !3 tables: .pack_offset, .pack_len, .pack_tag
           jl=jl+js*2 !the last 2 INTL-byte words are this%length and this%num_packets
 !Resize the buffer if needed:
           if(jl.gt.cap) then
            call this%resize(jerr,buf_size=jl)
            if(jerr.eq.PACK_SUCCESS) cap=this%get_capacity(jerr)
           endif
           if(jerr.eq.PACK_SUCCESS) then
            jl=this%length
 !Packet offsets:
            cptr=c_loc(this%pack_offset)
            call c_f_pointer(cptr,fptr,(/js*this%num_packets/))
            this%buffer(jl+1:jl+js*this%num_packets)=fptr(1:js*this%num_packets)
            jl=jl+js*this%num_packets
 !Packet lengths:
            cptr=c_loc(this%pack_len)
            call c_f_pointer(cptr,fptr,(/js*this%num_packets/))
            this%buffer(jl+1:jl+js*this%num_packets)=fptr(1:js*this%num_packets)
            jl=jl+js*this%num_packets
 !Packet tags:
            cptr=c_loc(this%pack_tag)
            call c_f_pointer(cptr,fptr,(/js*this%num_packets/))
            this%buffer(jl+1:jl+js*this%num_packets)=fptr(1:js*this%num_packets)
            jl=jl+js*this%num_packets
 !Length of the data buffer:
            fptr(1:)=>this%buffer(jl+1:); cptr=c_loc(fptr)
            call c_f_pointer(cptr,j8p); j8p=this%length; jl=jl+js
 !Number of packets in the packet envelope:
            fptr(1:)=>this%buffer(jl+1:); cptr=c_loc(fptr)
            call c_f_pointer(cptr,j8p); j8p=int(this%num_packets,INTL); jl=jl+js
 !Tag of the packet envelope:
            fptr(1:)=>this%buffer(jl+1:); cptr=c_loc(fptr)
            call c_f_pointer(cptr,j8p); j8p=int(this%tag,INTL); jl=jl+js
           else
            jl=0_INTL
           endif
           return
          end subroutine pack_layout_tables

          subroutine send_mpi_message(jc,jr,jl,jbuf,jtag,jreq,jerr)
           implicit none
           integer(INT_MPI), intent(in):: jc         !in: MPI communicator
           integer(INT_MPI), intent(in):: jr         !in: MPI rank
           integer(INTL), intent(in):: jl            !in: count
           character(C_CHAR), intent(in):: jbuf(1:*) !in: buffer (C_CHAR)
           integer(INT_MPI), intent(in):: jtag       !in: tag
           integer(INT_MPI), intent(in):: jreq       !in: request handle
           integer(INTD), intent(out):: jerr         !out: error code
           integer(INT_MPI):: jcnt,jcs,jer
           logical:: intercomm

           jerr=PACK_SUCCESS; jcnt=0
           call MPI_Comm_test_inter(jc,intercomm,jer)
           if(jer.eq.MPI_SUCCESS) then
            if(intercomm) then
             call MPI_Comm_remote_size(jc,jcs,jer)
            else
             call MPI_Comm_size(jc,jcs,jer)
            endif
            if(jer.eq.MPI_SUCCESS.and.jr.ge.0.and.jr.lt.jcs) then
             if(jl.le.int(huge(jcnt),INTL)) then
              jcnt=int(jl,INT_MPI)
              call MPI_Isend(jbuf,jcnt,MPI_CHARACTER,jr,jtag,jc,jreq,jer)
              if(jer.ne.MPI_SUCCESS) jerr=PACK_MPI_ERR
             else
              jerr=PACK_OVERFLOW
             endif
            else
             jerr=PACK_INVALID_ARGS
            endif
           else
            jerr=PACK_MPI_ERR
           endif
           return
          end subroutine send_mpi_message

        end subroutine PackEnvSend
!------------------------------------------------------------------------------------------
        function PackEnvReceive(this,comm_handle,ierr,proc_rank,tag,comm) result(delivered)
!Initiates a receive of a packet envelope from some MPI process, either specified or any,
!but only if the corresponding message is already available, otherwise simply returns FALSE.
!The packet envelope passed into this subroutine must be clean on entrance (zero length).
!The communication handle must not be active on entrance.
!This is a NON-BLOCKING method!
         implicit none
         logical:: delivered                                !TRUE if the message is available for pickup, FALSE otherwise
         class(pack_env_t), intent(inout):: this            !out: packet envelope (clean on entrance)
         class(comm_handle_t), intent(inout):: comm_handle  !out: communication handle (must not be active on entrance)
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INT_MPI), intent(in), optional:: proc_rank !in: process rank (defaults to MPI_ANY_SOURCE)
         integer(INT_MPI), intent(in), optional:: tag       !in: MPI message tag (defaults to MPI_ANY_TAG)
         integer(INT_MPI), intent(in), optional:: comm      !in: MPI communicator (defaults to MPI_COMM_WORLD)
         integer(INTD):: errc
         integer(INTL):: cap
         integer(INT_MPI):: rk,tg,cm,rh,ml,err_mpi

         delivered=.FALSE.; cap=this%get_capacity(errc)
         if(cap.gt.0.and.errc.eq.PACK_SUCCESS) then
          if(this%get_length(errc).eq.0) then
           if(this%get_num_packets(errc).eq.0) then
            if(.not.this%is_busy(errc)) then
             if(errc.eq.PACK_SUCCESS) then
              if(.not.comm_handle%is_active(errc)) then
               if(errc.eq.PACK_SUCCESS) then
                call comm_handle%clean(errc) !clean the communication handle before usage
                if(errc.eq.PACK_SUCCESS) then
                 if(present(proc_rank)) then; rk=proc_rank; else; rk=MPI_ANY_SOURCE; endif
                 if(present(tag)) then; tg=tag; else; tg=MPI_ANY_TAG; endif
                 if(present(comm)) then; cm=comm; else; cm=MPI_COMM_WORLD; endif
                 call MPI_Iprobe(rk,tg,cm,delivered,comm_handle%stat,err_mpi)
                 if(err_mpi.eq.MPI_SUCCESS.and.delivered) then
                  call MPI_Get_Count(comm_handle%stat,MPI_CHARACTER,ml,err_mpi)
                  if(err_mpi.eq.MPI_SUCCESS) then
                   if(cap.lt.int(ml,INTL)) call this%resize(errc,buf_size=int(ml,INTL))
                   if(errc.eq.PACK_SUCCESS) then
                    call receive_mpi_message(cm,rk,ml,this%buffer,tg,rh,errc)
                    if(errc.eq.PACK_SUCCESS) call comm_handle%construct(cm,rh,errc,this)
                   endif
                  else
                   errc=PACK_MPI_ERR
                  endif
                 else
                  if(err_mpi.ne.MPI_SUCCESS) errc=PACK_MPI_ERR
                 endif
                endif
               endif
              else
               if(errc.eq.PACK_SUCCESS) errc=PACK_INVALID_ARGS
              endif
             endif
            else
             if(errc.eq.PACK_SUCCESS) errc=PACK_BUSY
            endif
           else
            if(errc.eq.PACK_SUCCESS) errc=PACK_INVALID_ARGS
           endif
          else
           if(errc.eq.PACK_SUCCESS) errc=PACK_INVALID_ARGS
          endif
         else
          if(errc.eq.PACK_SUCCESS) errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return

         contains

          subroutine receive_mpi_message(jc,jr,jl,buf,jt,jreq,jerr)
           implicit none
           integer(INT_MPI), intent(in):: jc           !MPI communicator
           integer(INT_MPI), intent(in):: jr           !MPI rank
           integer(INT_MPI), intent(in):: jl           !length of the buffer
           character(C_CHAR), intent(inout):: buf(1:*) !buffer
           integer(INT_MPI), intent(in):: jt           !MPI message tag
           integer(INT_MPI), intent(inout):: jreq      !MPI request handle
           integer(INTD), intent(out):: jerr           !error code
           integer(INT_MPI):: jer

           jerr=PACK_SUCCESS
           if(jl.ge.0) then
            call MPI_Irecv(buf,jl,MPI_CHARACTER,jr,jt,jc,jreq,jer)
            if(jer.ne.MPI_SUCCESS) jerr=PACK_MPI_ERR
           else
            jerr=PACK_OVERFLOW
           endif
           return
          end subroutine receive_mpi_message

        end function PackEnvReceive
!--------------------------------------------------------
        subroutine PackEnvDecodeMPIMsg(this,msg_len,ierr)
!Decodes the received MPI message back into the <pack_env_t> object.
!Namely, it decodes the layout tables stored in the tail of the message.
         implicit none
         class(pack_env_t), intent(inout):: this     !out: decoded packet envelope (must be empty on entrance)
         integer(INTL), intent(in):: msg_len         !in: MPI message length
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: npkts,errc
         integer(INTL):: cap,sil,offs
         character(C_CHAR), pointer, contiguous:: chp(:)
         integer(INTL), pointer, contiguous:: j8p(:)
         integer(INTL), pointer:: i8p
         type(C_PTR):: cptr

         cap=this%get_capacity(errc)
         if(cap.ge.msg_len.and.errc.eq.PACK_SUCCESS) then
          if(this%get_length(errc).eq.0.and.msg_len.ge.int(1+3*1*8+3*8,INTL)) then !buffer (min 1 byte) + three tables (min 1 byte each) + three int(8) words
           if(errc.eq.PACK_SUCCESS) then
            sil=size_of(cap) !size of INTL integer
 !Tag of the packet envelope:
            chp(1:)=>this%buffer(msg_len-sil+1:msg_len)
            cptr=c_loc(chp); call c_f_pointer(cptr,i8p); this%tag=i8p
 !Number of packets in the packet envelope:
            chp(1:)=>this%buffer(msg_len-sil*2+1:msg_len-sil)
            cptr=c_loc(chp); call c_f_pointer(cptr,i8p); npkts=i8p
 !Length of the data buffer:
            chp(1:)=>this%buffer(msg_len-sil*3+1:msg_len-sil*2)
            cptr=c_loc(chp); call c_f_pointer(cptr,i8p); this%length=i8p
            i8p=>NULL()
            if(npkts.gt.this%get_max_packets()) call this%resize(errc,max_packets=npkts)
            if(errc.eq.PACK_SUCCESS) then
             this%num_packets=npkts
 !Packet tags:
             offs=msg_len-sil*3-sil*npkts+1
             chp(1:)=>this%buffer(offs:)
             cptr=c_loc(chp); call c_f_pointer(cptr,j8p,(/npkts/))
             this%pack_tag(1:npkts)=j8p(1:npkts)
 !Packet lengths:
             offs=offs-sil*npkts
             chp(1:)=>this%buffer(offs:)
             cptr=c_loc(chp); call c_f_pointer(cptr,j8p,(/npkts/))
             this%pack_len(1:npkts)=j8p(1:npkts)
 !Packet offsets:
             offs=offs-sil*npkts
             chp(1:)=>this%buffer(offs:)
             cptr=c_loc(chp); call c_f_pointer(cptr,j8p,(/npkts/))
             this%pack_offset(1:npkts)=j8p(1:npkts)
             chp=>NULL(); j8p=>NULL()
            endif
           endif
          else
           errc=PACK_INVALID_ARGS !non-empty packet envelope or invalid message length
          endif
         else
          errc=PACK_ERROR
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine PackEnvDecodeMPIMsg
![comm_handle_t]==========================================================
        subroutine CommHandleConstruct(this,comm,req_handle,ierr,pack_env)
!Constructs a communication handle. The communication handle
!must be clean on entrance. Note that upon construction the
!communication handle automatically becomes active (pending message).
!For receive operations, the <pack_env> argument is the packet envelope
!for which the (receive) communication has just been initiated.
         implicit none
         class(comm_handle_t), intent(inout):: this  !inout: communication handle (must be clean on entrance)
         integer(INT_MPI), intent(in):: comm         !in: MPI communicator
         integer(INT_MPI), intent(in):: req_handle   !in: MPI request handle
         integer(INTD), intent(out), optional:: ierr !out: error code
         class(pack_env_t), intent(in), target, optional:: pack_env !in: packet envelope (for receive operations only)
         integer(INTD):: errc

         if(.not.this%is_active(errc)) then
          if(errc.eq.PACK_SUCCESS) then
           if(this%is_clean(errc)) then
            this%req=req_handle
            this%comm=comm
            if(present(pack_env)) this%recv_pack_env=>pack_env
            this%active=.TRUE.
           else
            errc=PACK_INVALID_ARGS
           endif
          endif
         else
          errc=PACK_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine CommHandleConstruct
!--------------------------------------------------------
        subroutine CommHandleClean(this,ierr,force_clean)
!Cleans an inactive communication handle. If <force_clean> = TRUE,
!an active communication handle will be forcefully cleaned as well.
         implicit none
         class(comm_handle_t), intent(inout):: this  !inout: communication handle
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: force_clean !in: if TRUE, cleaning will be enforced
         integer(INTD):: errc
         logical:: fcl

         errc=PACK_SUCCESS
         if(present(force_clean)) then; fcl=force_clean; else; fcl=.FALSE.; endif
         if(fcl.or.(.not.this%is_active(errc))) then
          this%active=.FALSE.
          this%req=MPI_REQUEST_NULL
          this%comm=MPI_COMM_NULL
          this%recv_pack_env=>NULL()
         else
          errc=PACK_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine CommHandleClean
!---------------------------------------------------------
        function CommHandleIsClean(this,ierr) result(answ)
!Returns TRUE if the communication handle is clean, FALSE otherwise.
         implicit none
         logical:: answ                              !out: result
         class(comm_handle_t), intent(in):: this     !in: communication handle
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=PACK_SUCCESS; answ=.TRUE.
         if(this%comm.ne.MPI_COMM_NULL) then
          answ=.FALSE.
         else
          if(this%active) errc=PACK_ERROR
         endif
         if(present(ierr)) ierr=errc
         return
        end function CommHandleIsClean
!----------------------------------------------------------
        function CommHandleIsActive(this,ierr) result(answ)
!Returns TRUE if the communication handle is active, FALSE otherwise.
!An active communication handle is the one associated with an active
!non-blocking MPI message.
         implicit none
         logical:: answ                              !out: result
         class(comm_handle_t), intent(in):: this     !in: communication handle
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=PACK_SUCCESS; answ=.FALSE.
         if(this%active) then
          if(this%comm.ne.MPI_COMM_NULL) then
           answ=.TRUE.
          else
           errc=PACK_ERROR
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end function CommHandleIsActive
!------------------------------------------------------
        function CommHandleTest(this,ierr) result(answ)
!Tests the completion of an active communication. An error
!code PACK_NULL is returned if the communication handle is clean.
!If completed, the communication handle becomes inactive, but not clean.
         implicit none
         logical:: answ                              !out: result
         class(comm_handle_t), intent(inout):: this  !inout: communication handle
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer(INT_MPI):: err_mpi
         integer(INT_MPI):: ml

         errc=PACK_SUCCESS; answ=.FALSE.
         if(.not.this%is_clean(errc)) then
          if(errc.eq.PACK_SUCCESS) then
           if(this%is_active(errc)) then !an active communication handle
            if(errc.eq.PACK_SUCCESS) then
             call MPI_Test(this%req,answ,this%stat,err_mpi)
             if(err_mpi.eq.MPI_SUCCESS) then
              if(answ) then
               if(associated(this%recv_pack_env)) then !receive operation required decoding
                call MPI_Get_Count(this%stat,MPI_CHARACTER,ml,err_mpi)
                if(err_mpi.eq.MPI_SUCCESS) then
                 call this%recv_pack_env%decode_mpi_msg(int(ml,INTL),errc)
                 if(errc.eq.PACK_SUCCESS) this%recv_pack_env=>NULL()
                else
                 errc=PACK_MPI_ERR
                endif
               endif
               this%active=.FALSE.
              endif
             else
              errc=PACK_MPI_ERR
             endif
            endif
           else
            answ=.TRUE. !already completed communication handle
           endif
          endif
         else
          errc=PACK_NULL !communication handle is clean
         endif
         if(present(ierr)) ierr=errc
         return
        end function CommHandleTest
!-------------------------------------------
        subroutine CommHandleWait(this,ierr)
!Waits for the completion of an active communication. An error
!code PACK_NULL is returned if the communication handle is clean.
!When completed, the communication handle becomes inactive, but not clean.
         implicit none
         class(comm_handle_t), intent(inout):: this  !inout: communication handle
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer(INT_MPI):: err_mpi
         integer(INT_MPI):: ml

         errc=PACK_SUCCESS
         if(.not.this%is_clean(errc)) then
          if(errc.eq.PACK_SUCCESS) then
           if(this%is_active(errc)) then !an active communication handle
            if(errc.eq.PACK_SUCCESS) then
             call MPI_Wait(this%req,this%stat,err_mpi)
             if(err_mpi.eq.MPI_SUCCESS) then
              if(associated(this%recv_pack_env)) then !receive operation required decoding
               call MPI_Get_Count(this%stat,MPI_CHARACTER,ml,err_mpi)
               if(err_mpi.eq.MPI_SUCCESS) then
                call this%recv_pack_env%decode_mpi_msg(int(ml,INTL),errc)
                if(errc.eq.PACK_SUCCESS) this%recv_pack_env=>NULL()
               else
                errc=PACK_MPI_ERR
               endif
              endif
              this%active=.FALSE.
             else
              errc=PACK_MPI_ERR
             endif
            endif
           endif
          endif
         else
          errc=PACK_NULL !communication handle is clean
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine CommHandleWait
!-------------------------------------------------------------
        subroutine CommHandleAttachPackEnv(this,pack_env,ierr)
!Attaches the parental packet envelope for each the receive operation
!has just been initiated that corresponds to this communication handle.
         implicit none
         class(comm_handle_t), intent(inout):: this       !inout: communication handle
         class(pack_env_t), target, intent(in):: pack_env !in: packet envelope
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD):: errc

         errc=PACK_SUCCESS
         this%recv_pack_env=>pack_env
         if(present(ierr)) ierr=errc
         return
        end subroutine CommHandleAttachPackEnv
!-----------------------------------------
        subroutine CommHandlePrintIt(this)
         implicit none
         class(comm_handle_t), intent(in):: this

         write(CONS_OUT,'("comm_handle_t{",l1,1x,i11,1x,i11,"}")') this%active,this%req,this%comm
         return
        end subroutine CommHandlePrintIt
![non-member]====================================
!PACKING/UNPACKING for built-in types:
        subroutine pack_integer1(packet,obj,ierr)
!Packs object <obj> into packet <packet>. The length of the packet
!is increased by the storage size of the object in bytes.
         implicit none
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(1), intent(in):: obj                !in: builtin type object
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: sl
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         integer(1), pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          sl=packet%space_left(errc)
          if(errc.eq.PACK_SUCCESS) then
           if(sl.ge.int(obj_size,INTL)) then
            chp(1:)=>packet%buffer(packet%length+1_INTL:)
            cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
            fptr=obj; fptr=>NULL()
            packet%length=packet%length+obj_size
           else
            errc=PACK_OVERFLOW
           endif
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine pack_integer1
!--------------------------------------------------
        subroutine unpack_integer1(packet,obj,ierr)
!Unpacks object <obj> from packet <packet>. After unpacking, the internal
!packet offset is automatically incremented to the next field.
         implicit none
         class(obj_pack_t), intent(inout):: packet       !inout: packet
         integer(1), intent(out):: obj                   !out: builtin type object
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: ppos
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         integer(1), pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          ppos=PACK_BASE+packet%offset
          if(ppos.gt.0_INTL.and.ppos+obj_size-1.le.packet%get_length(errc)) then
           if(errc.eq.PACK_SUCCESS) then
            chp(1:)=>packet%buffer(ppos:)
            cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
            obj=fptr; fptr=>NULL()
            packet%offset=packet%offset+obj_size
           endif
          else
           errc=PACK_OVERFLOW
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine unpack_integer1
!------------------------------------------------
        subroutine pack_integer2(packet,obj,ierr)
!Packs object <obj> into packet <packet>. The length of the packet
!is increased by the storage size of the object in bytes.
         implicit none
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(2), intent(in):: obj                !in: builtin type object
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: sl
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         integer(2), pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          sl=packet%space_left(errc)
          if(errc.eq.PACK_SUCCESS) then
           if(sl.ge.int(obj_size,INTL)) then
            chp(1:)=>packet%buffer(packet%length+1_INTL:)
            cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
            fptr=obj; fptr=>NULL()
            packet%length=packet%length+obj_size
           else
            errc=PACK_OVERFLOW
           endif
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine pack_integer2
!--------------------------------------------------
        subroutine unpack_integer2(packet,obj,ierr)
!Unpacks object <obj> from packet <packet>. After unpacking, the internal
!packet offset is automatically incremented to the next field.
         implicit none
         class(obj_pack_t), intent(inout):: packet       !inout: packet
         integer(2), intent(out):: obj                   !out: builtin type object
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: ppos
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         integer(2), pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          ppos=PACK_BASE+packet%offset
          if(ppos.gt.0_INTL.and.ppos+obj_size-1.le.packet%get_length(errc)) then
           if(errc.eq.PACK_SUCCESS) then
            chp(1:)=>packet%buffer(ppos:)
            cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
            obj=fptr; fptr=>NULL()
            packet%offset=packet%offset+obj_size
           endif
          else
           errc=PACK_OVERFLOW
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine unpack_integer2
!------------------------------------------------
        subroutine pack_integer4(packet,obj,ierr)
!Packs object <obj> into packet <packet>. The length of the packet
!is increased by the storage size of the object in bytes.
         implicit none
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(4), intent(in):: obj                !in: builtin type object
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: sl
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         integer(4), pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          sl=packet%space_left(errc)
          if(errc.eq.PACK_SUCCESS) then
           if(sl.ge.int(obj_size,INTL)) then
            chp(1:)=>packet%buffer(packet%length+1_INTL:)
            cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
            fptr=obj; fptr=>NULL()
            packet%length=packet%length+obj_size
           else
            errc=PACK_OVERFLOW
           endif
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine pack_integer4
!---------------------------------------------------------------
        subroutine pack_integer4_arr1(packet,objs,num_objs,ierr)
!Packs <num_objs> objects <objs> into packet <packet>. The length
!of the packet is increased by the storage size of all objects in bytes.
         implicit none
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(4), intent(in):: objs(1:)           !in: builtin type objects
         integer(4), intent(in):: num_objs           !in: number of objects to pack
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: obj_size,errc,i
         integer(INTL):: sl
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         integer(4), pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(objs(1)) !size of the object in bytes
         if(obj_size.gt.0) then
          sl=packet%space_left(errc)
          if(errc.eq.PACK_SUCCESS) then
           if(sl.ge.int(obj_size,INTL)*int(num_objs,INTL)) then
            do i=1,num_objs
             chp(1:)=>packet%buffer(packet%length+1_INTL:)
             cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
             fptr=objs(i); fptr=>NULL()
             packet%length=packet%length+obj_size
            enddo
           else
            errc=PACK_OVERFLOW
           endif
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine pack_integer4_arr1
!--------------------------------------------------
        subroutine unpack_integer4(packet,obj,ierr)
!Unpacks object <obj> from packet <packet>. After unpacking, the internal
!packet offset is automatically incremented to the next field.
         implicit none
         class(obj_pack_t), intent(inout):: packet       !inout: packet
         integer(4), intent(out):: obj                   !out: builtin type object
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: ppos
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         integer(4), pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          ppos=PACK_BASE+packet%offset
          if(ppos.gt.0_INTL.and.ppos+obj_size-1.le.packet%get_length(errc)) then
           if(errc.eq.PACK_SUCCESS) then
            chp(1:)=>packet%buffer(ppos:)
            cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
            obj=fptr; fptr=>NULL()
            packet%offset=packet%offset+obj_size
           endif
          else
           errc=PACK_OVERFLOW
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine unpack_integer4
!-----------------------------------------------------------------
        subroutine unpack_integer4_arr1(packet,objs,num_objs,ierr)
!Unpacks <num_objs> objects <objs> from packet <packet>. After unpacking,
!the internal packet offset is automatically incremented to the next field.
         implicit none
         class(obj_pack_t), intent(inout):: packet       !inout: packet
         integer(4), intent(inout):: objs(1:)            !out: builtin type objects
         integer(4), intent(in):: num_objs               !in: number of objects to unpack
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: obj_size,errc,i
         integer(INTL):: ppos
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         integer(4), pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(objs(1)) !size of the object in bytes
         if(obj_size.gt.0) then
          ppos=PACK_BASE+packet%offset
          if(ppos.gt.0_INTL.and.ppos+(obj_size*num_objs)-1.le.packet%get_length(errc)) then
           if(errc.eq.PACK_SUCCESS) then
            do i=1,num_objs
             chp(1:)=>packet%buffer(ppos:)
             cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
             objs(i)=fptr; fptr=>NULL()
             packet%offset=packet%offset+obj_size
             ppos=PACK_BASE+packet%offset
            enddo
           endif
          else
           errc=PACK_OVERFLOW
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine unpack_integer4_arr1
!------------------------------------------------
        subroutine pack_integer8(packet,obj,ierr)
!Packs object <obj> into packet <packet>. The length of the packet
!is increased by the storage size of the object in bytes.
         implicit none
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(8), intent(in):: obj                !in: builtin type object
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: sl
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         integer(8), pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          sl=packet%space_left(errc)
          if(errc.eq.PACK_SUCCESS) then
           if(sl.ge.int(obj_size,INTL)) then
            chp(1:)=>packet%buffer(packet%length+1_INTL:)
            cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
            fptr=obj; fptr=>NULL()
            packet%length=packet%length+obj_size
           else
            errc=PACK_OVERFLOW
           endif
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine pack_integer8
!--------------------------------------------------
        subroutine unpack_integer8(packet,obj,ierr)
!Unpacks object <obj> from packet <packet>. After unpacking, the internal
!packet offset is automatically incremented to the next field.
         implicit none
         class(obj_pack_t), intent(inout):: packet       !inout: packet
         integer(8), intent(out):: obj                   !out: builtin type object
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: ppos
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         integer(8), pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          ppos=PACK_BASE+packet%offset
          if(ppos.gt.0_INTL.and.ppos+obj_size-1.le.packet%get_length(errc)) then
           if(errc.eq.PACK_SUCCESS) then
            chp(1:)=>packet%buffer(ppos:)
            cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
            obj=fptr; fptr=>NULL()
            packet%offset=packet%offset+obj_size
           endif
          else
           errc=PACK_OVERFLOW
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine unpack_integer8
!-----------------------------------------------
        subroutine pack_logical(packet,obj,ierr)
!Packs object <obj> into packet <packet>. The length of the packet
!is increased by the storage size of the object in bytes.
         implicit none
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         logical, intent(in):: obj                   !in: builtin type object
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: sl
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         logical, pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          sl=packet%space_left(errc)
          if(errc.eq.PACK_SUCCESS) then
           if(sl.ge.int(obj_size,INTL)) then
            chp(1:)=>packet%buffer(packet%length+1_INTL:)
            cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
            fptr=obj; fptr=>NULL()
            packet%length=packet%length+obj_size
           else
            errc=PACK_OVERFLOW
           endif
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine pack_logical
!-------------------------------------------------
        subroutine unpack_logical(packet,obj,ierr)
!Unpacks object <obj> from packet <packet>. After unpacking, the internal
!packet offset is automatically incremented to the next field.
         implicit none
         class(obj_pack_t), intent(inout):: packet       !inout: packet
         logical, intent(out):: obj                      !out: builtin type object
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: ppos
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         logical, pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          ppos=PACK_BASE+packet%offset
          if(ppos.gt.0_INTL.and.ppos+obj_size-1.le.packet%get_length(errc)) then
           if(errc.eq.PACK_SUCCESS) then
            chp(1:)=>packet%buffer(ppos:)
            cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
            obj=fptr; fptr=>NULL()
            packet%offset=packet%offset+obj_size
           endif
          else
           errc=PACK_OVERFLOW
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine unpack_logical
!---------------------------------------------
        subroutine pack_real4(packet,obj,ierr)
!Packs object <obj> into packet <packet>. The length of the packet
!is increased by the storage size of the object in bytes.
         implicit none
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         real(4), intent(in):: obj                   !in: builtin type object
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: sl
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         real(4), pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          sl=packet%space_left(errc)
          if(errc.eq.PACK_SUCCESS) then
           if(sl.ge.int(obj_size,INTL)) then
            chp(1:)=>packet%buffer(packet%length+1_INTL:)
            cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
            fptr=obj; fptr=>NULL()
            packet%length=packet%length+obj_size
           else
            errc=PACK_OVERFLOW
           endif
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine pack_real4
!-----------------------------------------------
        subroutine unpack_real4(packet,obj,ierr)
!Unpacks object <obj> from packet <packet>. After unpacking, the internal
!packet offset is automatically incremented to the next field.
         implicit none
         class(obj_pack_t), intent(inout):: packet       !inout: packet
         real(4), intent(out):: obj                      !out: builtin type object
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: ppos
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         real(4), pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          ppos=PACK_BASE+packet%offset
          if(ppos.gt.0_INTL.and.ppos+obj_size-1.le.packet%get_length(errc)) then
           if(errc.eq.PACK_SUCCESS) then
            chp(1:)=>packet%buffer(ppos:)
            cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
            obj=fptr; fptr=>NULL()
            packet%offset=packet%offset+obj_size
           endif
          else
           errc=PACK_OVERFLOW
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine unpack_real4
!---------------------------------------------
        subroutine pack_real8(packet,obj,ierr)
!Packs object <obj> into packet <packet>. The length of the packet
!is increased by the storage size of the object in bytes.
         implicit none
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         real(8), intent(in):: obj                   !in: builtin type object
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: sl
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         real(8), pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          sl=packet%space_left(errc)
          if(errc.eq.PACK_SUCCESS) then
           if(sl.ge.int(obj_size,INTL)) then
            chp(1:)=>packet%buffer(packet%length+1_INTL:)
            cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
            fptr=obj; fptr=>NULL()
            packet%length=packet%length+obj_size
           else
            errc=PACK_OVERFLOW
           endif
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine pack_real8
!-----------------------------------------------
        subroutine unpack_real8(packet,obj,ierr)
!Unpacks object <obj> from packet <packet>. After unpacking, the internal
!packet offset is automatically incremented to the next field.
         implicit none
         class(obj_pack_t), intent(inout):: packet       !inout: packet
         real(8), intent(out):: obj                      !out: builtin type object
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: ppos
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         real(8), pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          ppos=PACK_BASE+packet%offset
          if(ppos.gt.0_INTL.and.ppos+obj_size-1.le.packet%get_length(errc)) then
           if(errc.eq.PACK_SUCCESS) then
            chp(1:)=>packet%buffer(ppos:)
            cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
            obj=fptr; fptr=>NULL()
            packet%offset=packet%offset+obj_size
           endif
          else
           errc=PACK_OVERFLOW
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine unpack_real8
!------------------------------------------------
        subroutine pack_complex4(packet,obj,ierr)
!Packs object <obj> into packet <packet>. The length of the packet
!is increased by the storage size of the object in bytes.
         implicit none
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         complex(4), intent(in):: obj                !in: builtin type object
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: sl
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         complex(4), pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          sl=packet%space_left(errc)
          if(errc.eq.PACK_SUCCESS) then
           if(sl.ge.int(obj_size,INTL)) then
            chp(1:)=>packet%buffer(packet%length+1_INTL:)
            cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
            fptr=obj; fptr=>NULL()
            packet%length=packet%length+obj_size
           else
            errc=PACK_OVERFLOW
           endif
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine pack_complex4
!--------------------------------------------------
        subroutine unpack_complex4(packet,obj,ierr)
!Unpacks object <obj> from packet <packet>. After unpacking, the intenal
!packet offset is automatically incremented to the next field.
         implicit none
         class(obj_pack_t), intent(inout):: packet       !inout: packet
         complex(4), intent(out):: obj                   !out: builtin type object
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: ppos
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         complex(4), pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          ppos=PACK_BASE+packet%offset
          if(ppos.gt.0_INTL.and.ppos+obj_size-1.le.packet%get_length(errc)) then
           if(errc.eq.PACK_SUCCESS) then
            chp(1:)=>packet%buffer(ppos:)
            cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
            obj=fptr; fptr=>NULL()
            packet%offset=packet%offset+obj_size
           endif
          else
           errc=PACK_OVERFLOW
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine unpack_complex4
!------------------------------------------------
        subroutine pack_complex8(packet,obj,ierr)
!Packs object <obj> into packet <packet>. The length of the packet
!is increased by the storage size of the object in bytes.
         implicit none
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         complex(8), intent(in):: obj                !in: builtin type object
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: sl
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         complex(8), pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          sl=packet%space_left(errc)
          if(errc.eq.PACK_SUCCESS) then
           if(sl.ge.int(obj_size,INTL)) then
            chp(1:)=>packet%buffer(packet%length+1_INTL:)
            cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
            fptr=obj; fptr=>NULL()
            packet%length=packet%length+obj_size
           else
            errc=PACK_OVERFLOW
           endif
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine pack_complex8
!--------------------------------------------------
        subroutine unpack_complex8(packet,obj,ierr)
!Unpacks object <obj> from packet <packet>. After unpacking, the internal
!packet offset is automatically incremented to the next field.
         implicit none
         class(obj_pack_t), intent(inout):: packet       !inout: packet
         complex(8), intent(out):: obj                   !out: builtin type object
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: ppos
         type(C_PTR):: cptr
         character(C_CHAR), pointer, contiguous:: chp(:)
         complex(8), pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          ppos=PACK_BASE+packet%offset
          if(ppos.gt.0_INTL.and.ppos+obj_size-1.le.packet%get_length(errc)) then
           if(errc.eq.PACK_SUCCESS) then
            chp(1:)=>packet%buffer(ppos:)
            cptr=c_loc(chp); call c_f_pointer(cptr,fptr)
            obj=fptr; fptr=>NULL()
            packet%offset=packet%offset+obj_size
           endif
          else
           errc=PACK_OVERFLOW
          endif
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine unpack_complex8
!----------------------------------------------
        subroutine pack_string(packet,obj,ierr)
!Packs object <obj> into packet <packet>. The length of the packet
!is increased by the storage size of the object in bytes.
         implicit none
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         character(*), intent(in):: obj              !in: builtin type object
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: sl,l,i
         character(1):: ch

         errc=PACK_SUCCESS; l=len(obj)
         if(l.gt.0) then
          ch=obj(1:1); obj_size=storage_size(ch)/8 !size_of(ch) `Cray compiler bug !size of the object in bytes
          if(obj_size.eq.1) then
           sl=packet%space_left(errc)
           if(errc.eq.PACK_SUCCESS) then
            if(sl.ge.int(INTL+l,INTL)) then !8 is the leading integer(8) that contains the length of the string
             call pack_builtin(packet,l,errc)
             if(errc.eq.PACK_SUCCESS) then
              do i=1,l; packet%buffer(packet%length+i)=obj(i:i); enddo
              packet%length=packet%length+l
             endif
            else
             errc=PACK_OVERFLOW
            endif
           endif
          else
           write(*,'("#FATAL(pack_prim::pack_string): Fortran character size is not equal to 1: ",i11)') obj_size
           stop
          endif
         else !empty string
          sl=packet%space_left(errc)
          if(errc.eq.PACK_SUCCESS) then
           if(sl.ge.INTL) then
            l=0_INTL
            call pack_builtin(packet,l,errc)
           else
            errc=PACK_OVERFLOW
           endif
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine pack_string
!-----------------------------------------------------
        subroutine unpack_string(packet,obj,strl,ierr)
!Unpacks object <obj> from packet <packet>. After unpacking, the internal
!packet offset is automatically incremented to the next field.
         implicit none
         class(obj_pack_t), intent(inout):: packet       !inout: packet
         character(*), intent(inout):: obj               !out: builtin type object
         integer(INTL), intent(out):: strl               !out: unpacked string length
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: ppos,i,l,lo
         character(1):: ch

         errc=PACK_SUCCESS; strl=0_INTL; lo=len(obj)
         if(lo.gt.0) then
          ch=obj(1:1); obj_size=storage_size(ch)/8 !size_of(ch) `Cray compiler bug !size of the object in bytes
          if(obj_size.eq.1) then
           call unpack_builtin(packet,l,errc)
           if(errc.eq.PACK_SUCCESS) then
            if(l.gt.0) then
             if(l.le.lo) then
              ppos=PACK_BASE+packet%offset
              if(ppos+l-1_INTL.le.packet%get_length(errc)) then
               if(errc.eq.PACK_SUCCESS) then
                do i=1,l; obj(i:i)=packet%buffer(ppos+i-1_INTL); enddo
                packet%offset=packet%offset+l
                strl=l
               endif
              else
               if(errc.eq.PACK_SUCCESS) errc=PACK_ERROR
              endif
             else
              errc=PACK_OVERFLOW
             endif
            elseif(l.lt.0) then
             errc=PACK_ERROR
            endif
           endif
          else
           write(*,'("#FATAL(pack_prim::unpack_string): Fortran character size is not equal to 1: ",i11)') obj_size
           stop
          endif
         else
          errc=PACK_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine unpack_string

       end module pack_prim
!===================================================================================
!===================================================================================
!TESTING:
       module pack_prim_test
        use pack_prim
        use dil_basic, only: INTD,INTL
#ifdef USE_MPI_MOD
#ifdef FORTRAN2008
        use mpi_f08      !MPI Fortran 2008 interface `Will not work
#else
        use mpi          !MPI Fortran interface
#endif
        implicit none
        private
#else
        implicit none
        private
        include 'mpif.h' !MPI Fortran interface
#endif
        public test_pack_prim

       contains

        function test_pack_prim(errc) result(ierr)
         implicit none
         integer:: ierr
         integer(INTD), intent(out):: errc
         real(4), parameter:: EPS4=1d-6
         real(8), parameter:: EPS8=1d-13
         integer(1), parameter:: i1=-63_1
         integer(2), parameter:: i2=-1645_2
         integer(4), parameter:: i4=-716894563_4
         integer(8), parameter:: i8=-1143557645657_8
         logical, parameter:: ld=.TRUE.
         real(4), parameter:: r4=-13.767357
         real(8), parameter:: r8=-0.8347853456D-5
         complex(4), parameter:: c4=cmplx(r4,-r4,4)
         complex(8), parameter:: c8=cmplx(r8,-r8,8)
         character(27), parameter:: s27='You better work correctly!!'
!--------------------------------------------------------------------
         integer(1):: ii1
         integer(2):: ii2
         integer(4):: ii4
         integer(8):: ii8
         logical:: lld
         real(4):: rr4
         real(8):: rr8
         complex(4):: cc4
         complex(8):: cc8
         character(128):: str=' '
         integer(INTD):: my_rank,comm_size,i,n
         integer(INTL):: mtag,sl
         logical:: delivered
         type(pack_env_t):: envelope
         type(obj_pack_t):: packet
         type(comm_handle_t), allocatable:: comm_hl(:)

         ierr=0; errc=PACK_SUCCESS
         call MPI_Comm_Size(MPI_COMM_WORLD,comm_size,errc); if(errc.ne.PACK_SUCCESS) then; ierr=1; return; endif
         if(comm_size.lt.2) then; ierr=-1; errc=-555; return; endif
         call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,errc); if(errc.ne.PACK_SUCCESS) then; ierr=2; return; endif
         call envelope%reserve_mem(errc); if(errc.ne.PACK_SUCCESS) then; ierr=3; return; endif
         if(my_rank.eq.0) then
 !Process 0 is packing built-in types and sending the envelope to other processes:
          allocate(comm_hl(1:comm_size-1),STAT=ierr); if(ierr.ne.0) then; ierr=4; errc=-999; return; endif
  !Pack integer1 -> packet 1:
          call envelope%acquire_packet(packet,errc); if(errc.ne.PACK_SUCCESS) then; ierr=5; return; endif
          call pack_builtin(packet,i1,errc); if(errc.ne.PACK_SUCCESS) then; ierr=6; return; endif
          call envelope%seal_packet(errc,tag=1_INTL); if(errc.ne.PACK_SUCCESS) then; ierr=7; return; endif
  !Pack integer2 -> packet 2:
          call envelope%acquire_packet(packet,errc); if(errc.ne.PACK_SUCCESS) then; ierr=8; return; endif
          call pack_builtin(packet,i2,errc); if(errc.ne.PACK_SUCCESS) then; ierr=9; return; endif
          call envelope%seal_packet(errc,tag=2_INTL); if(errc.ne.PACK_SUCCESS) then; ierr=10; return; endif
  !Pack integer4 -> packet 3:
          call envelope%acquire_packet(packet,errc); if(errc.ne.PACK_SUCCESS) then; ierr=11; return; endif
          call pack_builtin(packet,i4,errc); if(errc.ne.PACK_SUCCESS) then; ierr=12; return; endif
          call envelope%seal_packet(errc,tag=3_INTL); if(errc.ne.PACK_SUCCESS) then; ierr=13; return; endif
  !Pack integer8 -> packet 4:
          call envelope%acquire_packet(packet,errc); if(errc.ne.PACK_SUCCESS) then; ierr=14; return; endif
          call pack_builtin(packet,i8,errc); if(errc.ne.PACK_SUCCESS) then; ierr=15; return; endif
          call envelope%seal_packet(errc,tag=4_INTL); if(errc.ne.PACK_SUCCESS) then; ierr=16; return; endif
  !Pack logical -> packet 5:
          call envelope%acquire_packet(packet,errc); if(errc.ne.PACK_SUCCESS) then; ierr=17; return; endif
          call pack_builtin(packet,ld,errc); if(errc.ne.PACK_SUCCESS) then; ierr=18; return; endif
          call envelope%seal_packet(errc,tag=5_INTL); if(errc.ne.PACK_SUCCESS) then; ierr=19; return; endif
  !Pack real4 -> packet 6:
          call envelope%acquire_packet(packet,errc); if(errc.ne.PACK_SUCCESS) then; ierr=20; return; endif
          call pack_builtin(packet,r4,errc); if(errc.ne.PACK_SUCCESS) then; ierr=21; return; endif
          call envelope%seal_packet(errc,tag=6_INTL); if(errc.ne.PACK_SUCCESS) then; ierr=22; return; endif
  !Pack real8 -> packet 7:
          call envelope%acquire_packet(packet,errc); if(errc.ne.PACK_SUCCESS) then; ierr=23; return; endif
          call pack_builtin(packet,r8,errc); if(errc.ne.PACK_SUCCESS) then; ierr=24; return; endif
          call envelope%seal_packet(errc,tag=7_INTL); if(errc.ne.PACK_SUCCESS) then; ierr=25; return; endif
  !Pack complex4 -> packet 8:
          call envelope%acquire_packet(packet,errc); if(errc.ne.PACK_SUCCESS) then; ierr=26; return; endif
          call pack_builtin(packet,c4,errc); if(errc.ne.PACK_SUCCESS) then; ierr=27; return; endif
          call envelope%seal_packet(errc,tag=8_INTL); if(errc.ne.PACK_SUCCESS) then; ierr=28; return; endif
  !Pack complex8 -> packet 9:
          call envelope%acquire_packet(packet,errc); if(errc.ne.PACK_SUCCESS) then; ierr=29; return; endif
          call pack_builtin(packet,c8,errc); if(errc.ne.PACK_SUCCESS) then; ierr=30; return; endif
          call envelope%seal_packet(errc,tag=9_INTL); if(errc.ne.PACK_SUCCESS) then; ierr=31; return; endif
  !Pack string -> packet 10:
          call envelope%acquire_packet(packet,errc); if(errc.ne.PACK_SUCCESS) then; ierr=32; return; endif
          call pack_builtin(packet,s27,errc); if(errc.ne.PACK_SUCCESS) then; ierr=33; return; endif
          call envelope%seal_packet(errc,tag=10_INTL); if(errc.ne.PACK_SUCCESS) then; ierr=34; return; endif
  !Send the packet envelope to other MPI processes (10 packets):
          do i=1,comm_size-1
           call envelope%send(i,comm_hl(i),errc,tag=13); if(errc.ne.PACK_SUCCESS) then; ierr=35; return; endif
          enddo
  !Test the completion of all sends:
          n=comm_size-1
          do while(n.gt.0)
           do i=1,comm_size-1
            if(comm_hl(i)%is_active(errc)) then
             if(errc.ne.PACK_SUCCESS) then; ierr=36; return; endif
             if(comm_hl(i)%test(errc)) then
              if(errc.ne.PACK_SUCCESS) then; ierr=37; return; endif
              n=n-1; call comm_hl(i)%clean(errc)
              if(errc.ne.PACK_SUCCESS) then; ierr=38; return; endif
             endif
            endif
            if(errc.ne.PACK_SUCCESS) then; ierr=39; return; endif
           enddo
          enddo
          deallocate(comm_hl)
         else
 !Other processes are receiving the packet envelope from process 0:
          allocate(comm_hl(1),STAT=ierr); if(ierr.ne.0) then; ierr=40; errc=-888; return; endif
  !Initiate a receive of the packet envelope:
          do while(.not.envelope%receive(comm_hl(1),errc,proc_rank=0)) !wait until the message has been initiated
           if(errc.ne.PACK_SUCCESS) then; ierr=41; return; endif
          enddo
          if(errc.ne.PACK_SUCCESS) then; ierr=42; return; endif
  !Wait upon the completion of the receive:
          call comm_hl(1)%wait(errc); if(errc.ne.PACK_SUCCESS) then; ierr=43; return; endif
  !Unpack packets from the envelope (10 packets):
          n=0; n=envelope%get_num_packets(errc); if(errc.ne.PACK_SUCCESS) then; ierr=44; return; endif
          if(n.ne.10) then; ierr=45; errc=-777; return; endif
   !Unpack integer1 (packet 1):
          call envelope%extract_packet(1,packet,errc,tag=mtag,preclean=.TRUE.)
          if(errc.ne.PACK_SUCCESS) then; ierr=46; return; endif
          call unpack_builtin(packet,ii1,errc); if(errc.ne.PACK_SUCCESS) then; ierr=47; return; endif
          !write(*,'("#DEBUG[",i3,"]: i1 = ",i11)') my_rank,ii1 !debug
          if(ii1.ne.i1) then; ierr=48; errc=1001; return; endif
   !Unpack integer2 (packet 2):
          call envelope%extract_packet(2,packet,errc,tag=mtag,preclean=.TRUE.)
          if(errc.ne.PACK_SUCCESS) then; ierr=49; return; endif
          call unpack_builtin(packet,ii2,errc); if(errc.ne.PACK_SUCCESS) then; ierr=50; return; endif
          if(ii2.ne.i2) then; ierr=51; errc=1001; return; endif
   !Unpack integer4 (packet 3):
          call envelope%extract_packet(3,packet,errc,tag=mtag,preclean=.TRUE.)
          if(errc.ne.PACK_SUCCESS) then; ierr=52; return; endif
          call unpack_builtin(packet,ii4,errc); if(errc.ne.PACK_SUCCESS) then; ierr=53; return; endif
          if(ii4.ne.i4) then; ierr=54; errc=1001; return; endif
   !Unpack integer8 (packet 4):
          call envelope%extract_packet(4,packet,errc,tag=mtag,preclean=.TRUE.)
          if(errc.ne.PACK_SUCCESS) then; ierr=55; return; endif
          call unpack_builtin(packet,ii8,errc); if(errc.ne.PACK_SUCCESS) then; ierr=56; return; endif
          if(ii8.ne.i8) then; ierr=57; errc=1001; return; endif
   !Unpack logical (packet 5):
          call envelope%extract_packet(5,packet,errc,tag=mtag,preclean=.TRUE.)
          if(errc.ne.PACK_SUCCESS) then; ierr=58; return; endif
          call unpack_builtin(packet,lld,errc); if(errc.ne.PACK_SUCCESS) then; ierr=59; return; endif
          if(lld.neqv.ld) then; ierr=60; errc=1001; return; endif
   !Unpack real4 (packet 6):
          call envelope%extract_packet(6,packet,errc,tag=mtag,preclean=.TRUE.)
          if(errc.ne.PACK_SUCCESS) then; ierr=61; return; endif
          call unpack_builtin(packet,rr4,errc); if(errc.ne.PACK_SUCCESS) then; ierr=62; return; endif
          if(rr4.ne.r4) then; ierr=63; errc=1001; return; endif
   !Unpack real8 (packet 7):
          call envelope%extract_packet(7,packet,errc,tag=mtag,preclean=.TRUE.)
          if(errc.ne.PACK_SUCCESS) then; ierr=64; return; endif
          call unpack_builtin(packet,rr8,errc); if(errc.ne.PACK_SUCCESS) then; ierr=65; return; endif
          if(rr8.ne.r8) then; ierr=66; errc=1001; return; endif
   !Unpack complex4 (packet 8):
          call envelope%extract_packet(8,packet,errc,tag=mtag,preclean=.TRUE.)
          if(errc.ne.PACK_SUCCESS) then; ierr=67; return; endif
          call unpack_builtin(packet,cc4,errc); if(errc.ne.PACK_SUCCESS) then; ierr=68; return; endif
          if(cc4.ne.c4) then; ierr=69; errc=1001; return; endif
   !Unpack complex8 (packet 9):
          call envelope%extract_packet(9,packet,errc,tag=mtag,preclean=.TRUE.)
          if(errc.ne.PACK_SUCCESS) then; ierr=70; return; endif
          call unpack_builtin(packet,cc8,errc); if(errc.ne.PACK_SUCCESS) then; ierr=71; return; endif
          if(cc8.ne.c8) then; ierr=72; errc=1001; return; endif
   !Unpack string (packet 10):
          call envelope%extract_packet(10,packet,errc,tag=mtag,preclean=.TRUE.)
          if(errc.ne.PACK_SUCCESS) then; ierr=73; return; endif
          call unpack_builtin(packet,str,sl,errc); if(errc.ne.PACK_SUCCESS) then; ierr=74; return; endif
          !write(*,'("#DEBUG[",i3,"]: str = ",A27)') my_rank,str(1:27) !debug
          if(sl.ne.27.or.str(1:sl).ne.s27) then; ierr=75; errc=1001; return; endif
          deallocate(comm_hl)
         endif
         call envelope%destroy(errc); if(errc.ne.PACK_SUCCESS) then; ierr=76; return; endif
         return
        end function test_pack_prim

       end module pack_prim_test
