!Basic object packing/unpacking primitives.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2016/07/21

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
! Basic abstractions:
!  # Packet (obj_pack_t): Plain array which the data is packed into and unpacked from.
!  # Packet envelope (pack_env_t): Container of packets with some
!    additional layout information. A packet space is acquired from
!    an existing packet container. Then the data can be packed into
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
!    Once the object of a built-in or derived type has been packed into the packet,
!    the packet needs to be sealed via the .seal_packet() member procedure. If another
!    object needs to be packed, a new packet should be acquired, filled in, and sealed.
!    Once all objects have been packed into the packet envelope, the latter can be
!    communicated between two MPI processes via the .send(), .receive(), and
!    .test_completion() member procedures which follow asynchronous semantics. Once
!    received by another process, any packet from the delivered packet envelope can
!    be unpacked back into the original object of a built-in or derived type. Each
!    packet in the packet envelope has an optional integer tag.
!  # It may happen that the free space provided by a packet is insufficient for packing
!    the object of interest. In this case, the .resize() member procedure needs to be
!    invoked on the packet envelope the packet is part of. The packet envelope will
!    be extended to a user-specified capacity, in turn resulting in a larger buffer
!    in the packet that one is currently filling in with data. Similarly, if one
!    exceeds the max number of packets that can be stored in the packet envelope,
!    this resource can also be extended by calling the same member procedure .resize().
!    The packet buffer space overflow occurs during packing data objects into the packet.
!    The max-packets-per-envelope overflow occurs during sealing the packet.
       module pack_prim
#if 0
        use, intrinsic:: ISO_C_BINDING, only: C_PTR,C_INT,C_CHAR,C_NULL_PTR,c_loc,c_f_pointer
        use stsubs, only: size_of
        implicit none
        private
!PARAMETERS:
 !General:
        integer, private:: CONS_OUT=6     !output device
        logical, private:: VERBOSE=.TRUE. !verbosity for errors
        integer, private:: DEBUG=0        !debugging level (0:none)
 !Integers:
        integer, parameter, private:: INTD=4
        integer, parameter, private:: INTL=8
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
 !Packet envelope configuration:
        integer(INTD), parameter, private:: DEFAULT_MAX_PACKETS=1024 !default max number of packets per envelope
        integer(INTL), parameter, private:: DEFAULT_ENVELOPE_CAPACITY=1024_INTL*DEFAULT_MAX_PACKETS !default envelope capacity
        integer(INTL), parameter, private:: DEFAULT_PACKET_TAG=0 !default packet tag
!TYPES:
 !Packet:
        type, public:: obj_pack_t
         integer(INTL), private:: length=0 !used length of the packet buffer (bytes)
         character(C_CHAR), pointer, contiguous, private:: buffer(:)=>NULL() !buffer
         contains
          procedure, private:: construct=>ObjPackConstruct     !packet constructor (internal)
          procedure, private:: clean=>ObjPackClean             !packet cleaner (internal)
          procedure, public:: get_capacity=>ObjPackGetCapacity !returns the capacity of the packet buffer in bytes
          procedure, public:: get_length=>ObjPackGetLength     !returns the current length of the packet in bytes
          procedure, public:: has_room=>ObjPackHasRoom         !.TRUE. means one can add data to the packet, .FALSE. otherwise
          procedure, public:: space_left=>ObjPackSpaceLeft     !returns the amount of free space left in the packet buffer in bytes
        end type obj_pack_t
 !Packet envelope (communicable):
        type, public:: pack_env_t
         integer(INTL), private:: length=0      !used length of the packet envelope (bytes)
         integer(INTD), private:: num_packets=0 !number of packets in the packet envelope
         class(obj_pack_t), pointer, private:: curr_packet=>NULL() !current packet (set when in-use)
         logical, private:: busy=.FALSE.        !.TRUE. when there is an active packet being filled in (in-use flag)
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
          procedure, public:: is_busy=>PackEnvIsBusy                !check whether the packet envelope is currently in use
          procedure, public:: is_healthy=>PackEnvIsHealthy          !check whether the object is healthy (consistent)
          procedure, public:: reserve_mem=>PackEnvReserveMem        !reserve memory for the data buffer and/or layout tables
          procedure, public:: clean=>PackEnvClean                   !clean the packet envelope without releasing the memory
          procedure, public:: destroy=>PackEnvDestroy               !destroy the packet envelope completely
          procedure, public:: acquire_packet=>PackEnvAcquirePacket  !acquire a packet in the packet envelope
          procedure, public:: discard_packet=>PackEnvDiscardPacket  !discard a packet (either finished or unfinished)
          procedure, public:: seal_packet=>PackEnvSealPacket        !seal a packet (finalize)
          procedure, public:: extract_packet=>PackEnvExtractPacket  !extract a packet from the packet envelope
          procedure, public:: send=>PackEnvSend                     !send a packet to another process
          procedure, public:: receive=>PackEnvReceive               !receive a packet from another process
          procedure, public:: test_completion=>PackEnvCommCompleted !test the completion of the send/receive operation
        end type pack_env_t
 !Data communication handle:
        type, public:: comm_handle_t
         logical, private:: tested=.FALSE.                      !set to .TRUE. when an active communication handle has been tested at least once
        contains
         procedure, public:: construct=>PackCommHandleConstruct !construct the communication handle (internal)
         procedure, public:: clean=>PackCommHandleClean         !clean the communication handle
         procedure, public:: test=>PackCommHandleTest           !test the completion of the communication
         procedure, public:: wait=>PackCommHandleWait           !wait upon the completion of the communication
        end type comm_handle_t
!INTERFACES:
 !Packing for built-in types:
        interface pack_builtin
         module procedure pack_integer1
         module procedure pack_integer2
         module procedure pack_integer4
         module procedure pack_integer8
         module procedure pack_logical
         module procedure pack_real4
         module procedure pack_real8
         module procedure pack_complex4
         module procedure pack_complex8
         module procedure pack_string
        end interface pack_builtin
        public pack_builtin
        public pack_universal
 !Unpacking for built-in types:
        interface unpack_builtin
         module procedure unpack_integer1
         module procedure unpack_integer2
         module procedure unpack_integer4
         module procedure unpack_integer8
         module procedure unpack_logical
         module procedure unpack_real4
         module procedure unpack_real8
         module procedure unpack_complex4
         module procedure unpack_complex8
         module procedure unpack_string
        end interface unpack_builtin
        public unpack_builtin
        public unpack_universal

       contains
!DEFINITION:
!==========================================================
!CLASS obj_pack_t:
        subroutine ObjPackConstruct(this,buf_p,ierr,length)
!Constructs a packet (either empty or filled in).
         implicit none
         class(obj_pack_t), intent(inout):: this           !inout: packet (in: empty packet, out: allocated packet)
         character(C_CHAR), pointer, contiguous:: buf_p(:) !in: external buffer space
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTL), intent(in), optional:: length      !if present, the packet will be given this legnth (pre-existing data buffer)
         integer(INTD):: errc

         errc=PACK_SUCCESS
         if(this%get_capacity().le.0) then !empty packet
          if(associated(buf_p)) then
           this%buffer(1:)=>buf_p(:); this%length=0
           if(present(length)) then !non-empty packet constructor (empty if <length> = 0)
            if(length.ge.0.and.length.le.size(this%buffer)) then
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
         this%length=0; this%buffer=>NULL()
         if(present(ierr)) ierr=errc
         return
        end subroutine ObjPackClean
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
         if(this%get_capacity().gt.0) bytes=this%length
         if(bytes.lt.0) errc=PACK_ERROR
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
!===============================================================
!CLASS pack_env_t:
        subroutine PackEnvResize(this,ierr,buf_size,max_packets)
!Resizes either the packet buffer or the packet layout tables or both.
!The buffer is used for storing data. The layout tables contain
!bookeeping information and the size of these tables limits
!the max number of packets that can be stored in the packet envelope.
!The packet envelope passed into this procedure is allowed to be in
!the in-use state (with an active open packet).
         implicit none
         class(pack_env_t), intent(inout):: this           !inout: packet envelope
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTL), intent(in), optional:: buf_size    !in: new buffer size in bytes
         integer(INTD), intent(in), optional:: max_packets !in: new max limit on the number of packets stored
         integer:: jerr
         integer(INTD):: errc
         integer(INTL):: fl
         character(C_CHAR), pointer, contiguous:: chp(:)
         integer(INTL), pointer, contiguous:: i8p(:)

         errc=PACK_SUCCESS
         if(this%is_healthy(errc)) then
          fl=this%length
          if(associated(this%curr_packet)) fl=fl+this%curr_packet%get_length()
 !Buffer:
          if(present(buf_size)) then
           if(fl.le.buf_size) then
            allocate(chp(1:buf_size),STAT=jerr)
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
 !Tables:
          if(present(max_packets)) then
           if(this%num_packets.le.max_packets) then
            allocate(i8p(1:max_packets),STAT=jerr)
            if(jerr.eq.0) then
             i8p(1:this%num_packets)=this%pack_offset(1:this%num_packets)
             deallocate(this%pack_offset); this%pack_offset=>i8p; i8p=>NULL()
             allocate(i8p(1:max_packets),STAT=jerr)
             if(jerr.eq.0) then
              i8p(1:this%num_packets)=this%pack_len(1:this%num_packets)
              deallocate(this%pack_len); this%pack_len=>i8p; i8p=>NULL()
              allocate(i8p(1:max_packets,STAT=jerr)
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
         integer(INTL), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=PACK_SUCCESS; bytes=0
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
         if(bytes.gt.this%get_capacity()) then; errc=PACK_ERROR; bytes=-1_INTL; endif
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
          if(size(this%pack_len).ne.num.or.size(this%pack_tag).ne.num) then; errc=PACK_ERROR; num=-1; endif
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
         if(num.gt.this%get_max_packets(errc)) then; errc=PACK_ERROR; num=-1; endif
         if(present(ierr)) ierr=errc
         return
        end function PackEnvGetNumPackets
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

         answ=.FALSE.
         if(this%get_capacity().gt.0) then
          if(this%get_capacity().lt.this%get_length()) errc=PACK_ERROR
          if(this%get_max_packets().lt.this%get_num_packets()) errc=PACK_ERROR
          if((this%busy.and.(.not.associated(this%curr_packet))).or.&
             ((.not.this%busy).and.associated(this%curr_packet))) errc=PACK_ERROR
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
         integer(INTL), intent(in), optional:: mem_size    !in: memory size to reserve for the buffer
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
         if(this%is_healthy(errc)) then
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
               chp(1:this%length)=this%buffer(1:this%length)
               deallocate(this%buffer); this%buffer=>chp; chp=>NULL()
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
                  i8p(1:this%num_packets)=this%pack_offset(1:this%num_packets)
                  j8p(1:this%num_packets)=this%pack_len(1:this%num_packets)
                  k8p(1:this%num_packets)=this%pack_tag(1:this%num_packets)
                  deallocate(this%pack_tag); deallocate(this%pack_len); deallocate(this%pack_offset)
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
         else
          errc=PACK_ERROR
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
          this%length=0
          this%num_packets=0
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
         class(pack_env_t), intent(inout):: this     !inout: packet envelope
         class(obj_pack_t), intent(inout):: pkt      !inout: in:clean packet, out: active empty packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: preclean    !in: if TRUE the packet will be cleaned here before use
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
         integer(INTL):: bg,ln

         cap=this%get_capacity(errc)
         if(cap.gt.0.and.errc.eq.PACK_SUCCESS) then
          if(this%is_healthy(errc)) then
           if(errc.eq.PACK_SUCCESS) then
            if(present(preclean)) then; if(preclean) call pkt%clean(); endif
            if(pkt_num.ge.1.and.pkt_num.le.this%get_num_packets()) then
             if(present(tag)) tag=this%pack_tag(pkt_num)
             bg=this%pack_offset(pkt_num); ln=this%pack_len(pkt_num)
             call pkt%construct(this%buffer(bg:bg+ln-1_INTL),errc,ln)
             if(errc.ne.PACK_SUCCESS) call pkt%clean()
            else
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
!================================================
!PACKING/UNPACKING for built-in types:
        subroutine pack_integer1(packet,obj,ierr)
         implicit none
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(1), intent(in):: obj                !in: builtin type object
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: obj_size,errc
         integer(INTL):: sl
         type(C_PTR):: cptr
         integer(1), pointer:: fptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          sl=packet%space_left(errc)
          if(errc.eq.PACK_SUCCESS) then
           if(sl.ge.int(obj_size,INTL)) then
            cptr=c_loc(packet%buffer(packet%length+1_INTL))
            call c_f_pointer(cptr,fptr)
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
#endif
       end module pack_prim
!===================================================================================
!TESTING:
       module pack_prim_test
        use pack_prim
        implicit none
        private
        public test_pack_prim

       contains

        function test_pack_prim() result(passed)
         implicit none
         logical:: passed
         real(4), parameter:: EPS4=1d-6
         real(8), parameter:: EPS8=1d-13
         integer(1), parameter:: i1=-63_1
         integer(2), parameter:: i2=-1645_2
         integer(4), parameter:: i4=-716894563_4
         integer(8), parameter:: i8=-1143557645657_8
         logical, parameter:: l=.TRUE.
         real(4), parameter:: r4=-13.767357
         real(8), parameter:: r8=-0.8347853456D-5
         complex(4), parameter:: c4=cmplx(r4,-r4,4)
         complex(8), parameter:: c8=cmplx(r8,-r8,8)
         character(27), parameter:: s27='You better work correctly!!'

         passed=.TRUE.
         
         return
        end function test_pack_prim

       end module pack_prim_test
