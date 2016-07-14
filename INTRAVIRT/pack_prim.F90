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
! for a subsequent communication between MPI processes. It helps 
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
 !Packet envelope:
        integer, parameter, private:: DEFAULT_MAX_PACKETS=1024 !default max number of packets per envelope
        integer, parameter, private:: DEFAULT_ENVELOPE_CAPACITY=1048576 !default envelope capacity
!TYPES:
 !Packet:
        type, public:: obj_pack_t
         integer(INTL), private:: length=0   !used length of the packet buffer (bytes)
         character(C_CHAR), pointer, contiguous, private:: buffer(:)=>NULL() !buffer
         contains
          procedure, public:: get_capacity=>ObjPackGetCapacity !returns the capacity of the buffer in bytes
          procedure, public:: is_reserved=>ObjPackIsReserved   !.TRUE. means it can be used, .FALSE. it is empty (cannot be used)
          procedure, public:: space_left=>ObjPackSpaceLeft     !space left in the buffer in bytes
        end type obj_pack_t
 !Packet envelope (communicable):
        type, public:: pack_env_t
         integer(INTL), private:: length=0      !used length of the packet envelope (bytes)
         integer(INTL), private:: num_packets=0 !number of packets in the packet envelope
         logical, private:: busy=.FALSE.        !.TRUE. when there is an active packet being filled in
         integer(INTL), pointer, contiguous, private:: pack_len(:)=>NULL() !length of each packet present in the envelope
         integer(INTL), pointer, contiguous, private:: pack_tag(:)=>NULL() !tag for each packet present in the envelope
         character(C_CHAR), pointer, contiguous, private:: buffer(:)=>NULL() !buffer
         contains
          procedure, public:: get_capacity=>PackEnvGetCapacity      !get the current buffer capacity
          procedure, public:: get_length=>PackEnvGetLength          !get the current used length of the buffer
          procedure, public:: get_max_packets=>PackEnvGetMaxPackets !get the max amount of packets in the envelope
          procedure, public:: get_num_packets=>PackEnvGetNumPackets !get the current number of packets in the envelope
          procedure, public:: is_busy=>PackEnvIsBusy                !checks whether the packet envelope is busy with an actively filled in packet
          procedure, public:: is_healthy=>PackEnvIsHealthy          !checks whether the object is healthy
          procedure, public:: reserve_mem=>PackEnvReserveMem        !reserve memory for the buffer
          procedure, public:: clean=>PackEnvClean                   !cleans the packet envelope without releasing the memory
          procedure, public:: destroy=>PackEnvDestroy               !destroy the packet envelope completely
          procedure, public:: acquire_packet=>PackEnvAcquirePacket  !acquires a packet in the packet envelope
          procedure, public:: discard_packet=>PackEnvDiscardPacket  !discard unfinished packet
          procedure, public:: seal_packet=>PackEnvSealPacket        !seals a packet
          procedure, public:: extract_packet=>PackEnvExtractPacket  !extract a packet from the envelope
          procedure, public:: send=>PackEnvSend                     !sends a packet to another process
          procedure, public:: receive=>PackEnvReceive               !receives a packet from another process
          procedure, public:: comm_completed=>PackEnvCommCompleted  !tests the completion of the send/receive operation
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
!===========================================================
!CLASS obj_pack_t:
        function ObjPackGetCapacity(this,ierr) result(bytes)
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
        end function ObjPackGetCapacity
!---------------------------------------------------------
        function ObjPackIsReserved(this,ierr) result(answ)
!Returns .TRUE. if the packet is ready to be filled in, .FALSE. otherwise.
         implicit none
         logical:: answ                              !out: answer
         class(obj_pack_t), intent(in):: this        !in: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=PACK_SUCCESS; answ=.TRUE.
         if(this%get_capacity().le.0) answ=.FALSE.
         if(this%length.gt.this%get_capacity(errc)) then; errc=PACK_ERROR; answ=.FALSE.; endif
         if(present(ierr)) ierr=errc
         return
        end function ObjPackIsReserved
!---------------------------------------------------------
        function ObjPackSpaceLeft(this,ierr) result(bytes)
!Returns the number of bytes left free in the packet.
!A negative return means that the packet is not reserved.
         implicit none
         integer(INTL):: bytes                       !out: number of bytes left in the packet
         class(obj_pack_t), intent(in):: this        !in: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         bytes=-1
         if(this%is_reserved(errc)) then
          if(errc.eq.PACK_SUCCESS) bytes=this%get_capacity(errc)-this%length
         else
          errc=PACK_NULL
         endif
         if(present(ierr)) ierr=errc
         return
        end function ObjPackSpaceLeft
!======================================================
!CLASS pack_env_t:
        function PackEnvGetCapacity(this) result(bytes)
!Returns the capacity in bytes of the packet envelope.
         implicit none
         integer(INTL):: bytes                !out: envelope capacity in bytes
         class(pack_env_t), intent(in):: this !in: envelope
         bytes=0
         if(associated(this%buffer)) bytes=size(this%buffer)
         return
        end function PackEnvGetCapacity
!---------------------------------------------------------
        function PackEnvGetLength(this,ierr) result(bytes)
!Returns the current length in bytes of the packet envelope.
         implicit none
         integer(INTL):: bytes                       !out: envelope length in bytes
         class(pack_env_t), intent(in):: this        !in: envelope
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=PACK_SUCCESS; bytes=this%length
         if(bytes.gt.this%get_capacity()) then; errc=PACK_ERROR; bytes=-1; endif
         if(present(ierr)) ierr=errc
         return
        end function PackEnvGetLength
!-----------------------------------------------------------
        function PackEnvGetMaxPackets(this,ierr) result(num)
!Returns the max number of packets assumed for this envelope.
         implicit none
         integer(INTD):: num                         !out: max number of packets
         class(pack_env_t), intent(in):: this        !in: envelope
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         logical:: as1,as2

         errc=PACK_SUCCESS; num=0
         as1=associated(this%pack_len); as2=associated(this%pack_tag)
         if(as1.and.as2) then
          num=size(this%pack_len)
          if(size(this%pack_tag).ne.num) then; errc=PACK_ERROR; num=-1; endif
         else
          if(as1.or.as2) then; errc=PACK_ERROR; num=-1; endif
         endif
         if(present(ierr)) ierr=errc
         return
        end function PackEnvGetMaxPackets
!-----------------------------------------------------------
        function PackEnvGetNumPackets(this,ierr) result(num)
!Returns the current number of packets sealed in the envelope.
         implicit none
         integer(INTD):: num                         !out: number of sealed packets in the envelope
         class(pack_env_t), intent(in):: this        !in: envelope
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         num=this%num_packets
         if(num.gt.this%get_max_packets(errc)) then; errc=PACK_ERROR; num=-1; endif
         if(present(ierr)) ierr=errc
         return
        end function PackEnvGetNumPackets
!-----------------------------------------------------
        function PackEnvIsBusy(this,ierr) result(answ)
!Returns .TRUE. if the envelope is blocked by an active packet.
         implicit none
         logical:: answ                              !out: answer
         class(pack_env_t), intent(in):: this        !in: envelope
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
!Returns .TRUE. if the object if healthy, .FALSE. otherwise (corrupted).
         implicit none
         logical:: answ                              !out: answer
         class(pack_env_t), intent(in):: this        !in: envelope
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         answ=.TRUE.
         if(this%get_capacity().lt.this%get_length(errc)) then; errc=PACK_ERROR; answ=.FALSE.; endif
         if(this%get_max_packets().lt.this%get_num_packets(errc)) then; errc=PACK_ERROR; answ=.FALSE.; endif
         if(present(ierr)) ierr=errc
         return
        end function PackEnvIsHealthy
!-------------------------------------------------------------------
        subroutine PackEnvReserveMem(this,ierr,mem_size,max_packets)
!Reserves memory for a packet envelope. If the packet envelope has its memory
!already reserved, it will be re-reserved if the new amount of memory is larger
!than the original one, otherwise an error PACK_INVALID_ARGS will be returned
!without any action. The same policy applies to <max_packets>.
         implicit none
         class(pack_env_t), intent(inout):: this           !inout: packet envelope
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTL), intent(in), optional:: mem_size    !in: memory size to reserve
         integer(INTD), intent(in), optional:: max_packets !in: max number of packets to assume
         integer(INTD):: errc
         integer(INTL), pointer, contiguous:: i8p(:),j8p(:)
         character(C_CHAR), pointer, contiguous:: chp(:)
         integer:: jerr

         if(this%is_healthy(errc)) then
          if(errc.eq.PACK_SUCCESS) then
           if(present(mem_size)) then; if(mem_size.le.0) errc=PACK_INVALID_ARGS; endif
           if(present(max_packets)) then; if(max_packets.le.0) errc=PACK_INVALID_ARGS; endif
           if(errc.eq.PACK_SUCCESS) then
            if(present(mem_size)) then
 !Buffer:
             if(mem_size.gt.this%get_capacity()) then
              allocate(chp(1:mem_size),STAT=jerr)
              if(jerr.eq.0) then
               chp(1:this%length)=this%buffer(1:this%length)
               deallocate(this%buffer); this%buffer=>chp; chp=>NULL()
              else
               errc=PACK_ALLOC_FAILED
              endif
             elseif(mem_size.lt.this%get_capacity()) then
              errc=PACK_INVALID_ARGS
             endif
            endif
 !TABLES:
            if(present(max_packets)) then
             if(max_packets.gt.this%get_max_packets(errc)) then
              if(errc.eq.PACK_SUCCESS) then
               allocate(i8p(1:max_packets),STAT=jerr)
               if(jerr.eq.0) then
                allocate(j8p(1:max_packets),STAT=jerr)
                if(jerr.eq.0) then
                 i8p(1:this%num_packets)=this%pack_len(1:this%num_packets)
                 j8p(1:this%num_packets)=this%pack_tag(1:this%num_packets)
                 deallocate(this%pack_len); deallocate(this%pack_tag)
                 this%pack_len=>i8p; this%pack_tag=>j8p; i8p=>NULL(); j8=>NULL()
                else
                 deallocate(i8p); errc=PACK_ALLOC_FAILED
                endif
               else
                errc=PACK_ALLOC_FAILED
               endif
              endif
             elseif(max_packets.lt.this%get_max_packets(errc)) then
              if(errc.eq.PACK_SUCCESS) errc=PACK_INVALID_ARGS
             endif
            endif
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
!Cleans the envelope without deallocating the memory.
         implicit none
         class(pack_env_t), intent(inout):: this     !inout: envelope
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=PACK_SUCCESS
         this%length=0
         this%num_packets=0
         this%busy=.FALSE.
         if(present(ierr)) ierr=errc
         return
        end subroutine PackEnvClean
!-------------------------------------------
        subroutine PackEnvDestroy(this,ierr)
!Destroys the envelope.
         implicit none
         class(pack_env_t), intent(inout):: this     !inout: envelope
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(associated(this%pack_tag)) deallocate(this%pack_tag)
         if(associated(this%pack_len)) deallocate(this%pack_len)
         if(associated(this%buffer)) deallocate(this%buffer)
         call this%clean(errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine PackEnvDestroy
!================================================
!PACKING/UNPACKING for built-in types:
        subroutine pack_integer1(packet,obj,ierr)
         implicit none
         class(obj_pack_t), intent(inout):: packet
         integer(1), intent(in):: obj
         integer(INTD), intent(out), optional:: ierr
         integer(INTD):: obj_size,errc
         integer(INTL):: sl
         integer(1), pointer:: fptr
         type(C_PTR):: cptr

         errc=PACK_SUCCESS
         obj_size=size_of(obj) !size of the object in bytes
         if(obj_size.gt.0) then
          sl=packet%space_left(errc)
          if(errc.eq.PACK_SUCCESS) then
           if(sl.ge.int(obj_size,INTL)) then
            cptr=c_loc(packet%buffer(packet%length+1_INTL))
            call c_f_pointer(cptr,fptr)
            fptr=obj; packet%length=packet%length+obj_size
            fptr=>NULL()
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
!====================================================================================
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
         integer(1), parameter:: i1=-63
         integer(2), parameter:: i2=-1645
         integer(4), parameter:: i4=-716894563
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
