!Basic object packing/unpacking primitives.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2016/07/13

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
! This module provides basic primitives for packing objects of derived types
! for a subsequent communication between MPI processes. Its main purpose
! is to provide necessary private API for implementing public communication
! methods in derived types: Object%send() and Object%receive(). These public
! communication API may require packing/unpacking of data under the hood.
! Thus, the packing/unpacking primitives provided in this module are intended
! to facilitate the implementation of the .send() and .receive() methods
! in Fortran classes.
! Basic abstractions:
!  # Packet (obj_pack_t): Plain array which the data is packed into.
!  # Packet envelope (pack_env_t): Container of packets with some
!    additional layout information. A packet space is acquired from
!    an existing packet container. Then the data is packed into
!    that packet space. Then the packet container (envelope) can
!    be sent to a different MPI process that will be able to unpack
!    any packet from the packet envelope back to an object.
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
        integer(INTD), parameter, public:: PACK_BUSY=-6         !object in use by others
 !Packet envelope configuration:
        integer(INTD), parameter, private:: DEFAULT_MAX_PACKETS=1024 !default max number of packets per envelope
        integer(INTL), parameter, private:: DEFAULT_ENVELOPE_CAPACITY=1024_INTL*DEFAULT_MAX_PACKETS !default envelope capacity
!TYPES:
 !Packet:
        type, public:: obj_pack_t
         integer(INTL), private:: length=0 !used length of the packet buffer (bytes)
         character(C_CHAR), pointer, contiguous, private:: buffer(:)=>NULL() !buffer
         contains
          procedure, private:: construct=>ObjPackConstruct     !constructor
          procedure, private:: clean=>ObjPackClean             !cleaner
          procedure, public:: get_capacity=>ObjPackGetCapacity !returns the capacity of the buffer in bytes
          procedure, public:: has_room=>ObjPackHasRoom         !.TRUE. means one can add data to the packet, .FALSE. otherwise
          procedure, public:: space_left=>ObjPackSpaceLeft     !frees space left in the buffer in bytes
        end type obj_pack_t
 !Packet envelope (communicable):
        type, public:: pack_env_t
         integer(INTL), private:: length=0      !used length of the packet envelope (bytes)
         integer(INTL), private:: num_packets=0 !number of packets in the packet envelope
         class(obj_pack_t), pointer, private:: curr_packet=>NULL() !current packet
         logical, private:: busy=.FALSE.        !.TRUE. when there is an active packet being filled in (in use flag)
         integer(INTL), pointer, contiguous, private:: pack_offset(:)=>NULL() !offset of each packet in the envelope (byte)
         integer(INTL), pointer, contiguous, private:: pack_len(:)=>NULL()    !length of each packet present in the envelope (bytes)
         integer(INTL), pointer, contiguous, private:: pack_tag(:)=>NULL()    !tag for each packet present in the envelope
         character(C_CHAR), pointer, contiguous, private:: buffer(:)=>NULL()  !buffer
         contains
          procedure, public:: get_capacity=>PackEnvGetCapacity      !get the current buffer capacity (bytes)
          procedure, public:: get_length=>PackEnvGetLength          !get the current used length of the buffer (bytes)
          procedure, public:: get_max_packets=>PackEnvGetMaxPackets !get the max limit on the amount of packets in the envelope
          procedure, public:: get_num_packets=>PackEnvGetNumPackets !get the current number of packets in the envelope
          procedure, public:: is_busy=>PackEnvIsBusy                !check whether the packet envelope is currently in use
          procedure, public:: is_healthy=>PackEnvIsHealthy          !check whether the object is healthy (consistent)
          procedure, public:: reserve_mem=>PackEnvReserveMem        !reserve memory for the buffer
          procedure, public:: clean=>PackEnvClean                   !clean the packet envelope without releasing the memory
          procedure, public:: destroy=>PackEnvDestroy               !destroy the packet envelope completely
          procedure, public:: acquire_packet=>PackEnvAcquirePacket  !acquire a packet in the packet envelope
          procedure, public:: discard_packet=>PackEnvDiscardPacket  !discard a packet (either finished or unfinished)
          procedure, public:: seal_packet=>PackEnvSealPacket        !seal a packet (finalize)
          procedure, public:: extract_packet=>PackEnvExtractPacket  !extract a packet from the envelope
          procedure, public:: send=>PackEnvSend                     !send a packet to another process
          procedure, public:: receive=>PackEnvReceive               !receive a packet from another process
          procedure, public:: comm_completed=>PackEnvCommCompleted  !test the completion of the send/receive operation
        end type pack_env_t
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

       contains
!DEFINITION:
!===================================================
!CLASS obj_pack_t:
        subroutine ObjPackConstruct(this,buf_p,ierr)
!Constructs a packet.
         implicit none
         class(obj_pack_t), intent(inout):: this           !inout: packet (in: empty packet, out: allocated packet)
         character(C_CHAR), pointer, contiguous:: buf_p(:) !in: buffer space
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc

         errc=PACK_SUCCESS
         if(this%get_capacity().le.0) then !empty packet
          if(associated(buf_p)) then
           this%buffer(1:)=>buf_p(:)
          else
           errc=PACK_INVALID_ARGS
          endif
         else
          errc=PACK_BUSY
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
         integer(INTL):: bytes                       !out: buffer capacity
         class(obj_pack_t), intent(in):: this        !in: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=PACK_SUCCESS; bytes=0
         if(associated(this%buffer)) then
          bytes=size(this%buffer)
         endif
         if(present(ierr)) ierr=errc
         return
        end function ObjPackGetCapacity
!------------------------------------------------------
        function ObjPackHasRoom(this,ierr) result(answ)
!Returns .TRUE. if the packet can be used for packing data, .FALSE. otherwise.
         implicit none
         logical:: answ                              !out: answer
         class(obj_pack_t), intent(in):: this        !in: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer(INTL):: cap

         errc=PACK_SUCCESS; answ=.TRUE.
         cap=this%get_capacity(errc)
         if(cap.le.0_INTL) answ=.FALSE.
         if(this%length.ge.cap) then
          if(this%length.gt.cap) errc=PACK_ERROR !corrupted
          answ=.FALSE.
         endif
         if(present(ierr)) ierr=errc
         return
        end function ObjPackHasRoom
!---------------------------------------------------------
        function ObjPackSpaceLeft(this,ierr) result(bytes)
!Returns the number of bytes left in the packet.
!A negative return means that the packet is not reserved.
         implicit none
         integer(INTL):: bytes                       !out: number of bytes left in the packet
         class(obj_pack_t), intent(in):: this        !in: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         bytes=-1_INTL
         if(this%has_room(errc)) then
          if(errc.eq.PACK_SUCCESS) bytes=this%get_capacity(errc)-this%length
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end function ObjPackSpaceLeft
!===========================================================
!CLASS pack_env_t:
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
         implicit none
         logical:: answ                              !out: answer
         class(pack_env_t), intent(in):: this        !in: packet envelope
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         answ=.TRUE.
         if(this%is_healthy(errc)) then
          if(errc.eq.PACK_SUCCESS) then
           answ=this%busy
           if(answ.and.this%get_capacity().le.0) errc=PACK_ERROR
          endif
         else
          errc=PACK_ERROR
         endif
         if(present(ierr)) ierr=errc
         return
        end function PackEnvIsBusy
!--------------------------------------------------------
        function PackEnvIsHealthy(this,ierr) result(answ)
!Returns .TRUE. if the object is healthy, .FALSE. otherwise (corrupted).
         implicit none
         logical:: answ                              !out: answer
         class(pack_env_t), intent(in):: this        !in: packet envelope
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         answ=.FALSE.
         if(this%get_capacity().lt.this%get_length(errc)) errc=PACK_ERROR
         if(this%get_max_packets().lt.this%get_num_packets(errc)) errc=PACK_ERROR
         if(errc.eq.PACK_SUCCESS) answ=.TRUE.
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
         logical, intent(in), optional:: ignore_less       !in: if .TRUE. the resize will only occur if new > old (no error otherwise)
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
!-----------------------------------------------------
        subroutine PackEnvAcquirePacket(this,pkt,ierr)
!Acquires a packet space in the packet envelope.
         implicit none
         class(pack_env_t), intent(inout):: this     !inout: packet envelope
         class(obj_pack_t), intent(inout):: pkt      !out: clean packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer(INTL):: cap
         character(C_CHAR), pointer, contiguous:: chp(:)

         cap=this%get_capacity(errc)
         if(cap.gt.0) then
          if(this%is_healthy(errc)) then
           if(errc.eq.PACK_SUCCESS) then
            if(.not.this%busy) then
             if(this%length.lt.cap) then
              chp(1:cap-this%length)=>this%buffer(this%length+1:cap)
              call pkt%construct(chp,errc)
              if(errc.eq.PACK_SUCCESS) then
               this%current_packet=>pkt
               this%busy=.TRUE.
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
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine PackEnvAcquirePacket
!----------------------------------------------------------
        subroutine PackEnvDiscardPacket(this,ierr,pack_num)
!Discards a packet from the packet envelope. If <pack_num> is present,
!that specific (sealed) packet will be discarded. Otherwise, the currently
!active (unfinished) packet, if any, will be discarded.
         implicit none
         class(pack_env_t), intent(inout):: this        !inout: packet envelope
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTL), intent(in), optional:: pack_num !in: packet number
         integer(INTD):: errc

         

         return
        end subroutine PackEnvDiscardPacket
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
