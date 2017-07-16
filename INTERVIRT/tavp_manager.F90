!ExaTENSOR: TAVP Manager
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/07/16

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

       module tavp_manager
        use virta
        use gfc_base
        use gfc_list
        use gfc_dictionary
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6
        integer(INTD), private:: DEBUG=0
        logical, private:: VERBOSE=.TRUE.
!TYPES:
 !Tensor owner:
        type, private:: tens_owner_t
         class(tens_header_t), pointer, private:: tens_header=>NULL() !pointer to the tensor header
         integer(INTD), private:: owner_id=-1                         !non-negative tensor owner id
         type(tens_status_t), private:: tens_status                   !tensor status
         contains
          procedure, private:: TensOwnerCtor                          !ctor
          generic, public:: tens_owner_ctor=>TensOwnerCtor
          procedure, public:: is_set=>TensOwnerIsSet                  !returns TRUE if object is set
          procedure, public:: get_tens_header=>TensOwnerGetTensHeader !returns a pointer to the tensor header
          procedure, public:: get_owner_id=>TensOwnerGetOwnerId       !returns the tensor owner id
          final:: tens_owner_dtor                                     !dtor
        end type tens_owner_t
 !Tensor cache entry:
        type, private:: tens_entry_t
         class(tens_rcrsv_t), pointer, private:: tensor=>NULL() !composite tensor (consists of subtensors specified by their headers)
         type(list_bi_t), private:: owner_list                  !list of subtensor owners (tens_owner_t)
         type(tens_status_t), private:: tens_status             !tensor status
         logical, private:: tens_alloc=.FALSE.                  !TRUE if the tensor was allocated, FALSE if associated
         contains
          procedure, private:: TensEntryCtor                        !ctor
          generic, public:: tens_entry_ctor=>TensEntryCtor
          procedure, public:: is_set=>TensEntryIsSet                !returns TRUE if the tensor entry is set
          procedure, public:: get_tensor=>TensEntryGetTensor        !returns a pointer to the tensor
          procedure, public:: get_owner_list=>TensEntryGetOwnerList !returns a pointer to the subtensor owner list
          final:: tens_entry_dtor
        end type tens_entry_t
 !Tensor instruction (realization of a tensor operation for a specific TAVP):
        type, extends(ds_instr_t), private:: tens_instr_t
        contains
         procedure, private:: TensInstrCtor                     !ctor: constructs a tensor instruction from the specification of a tensor operation
         generic, public:: tens_instr_ctor=>TensInstrCtor
         procedure, public:: decode=>TensInstrDecode            !decoding procedure: Unpacks the raw byte packet (bytecode) and constructs a TAVP instruction
         procedure, public:: encode=>TensInstrEncode            !encoding procedure: Packs the TAVP instruction into a raw byte packet (bytecode)
         final:: tens_instr_dtor                                !dtor
        end type tens_instr_t
!INTERFACES:
!VISIBILITY:
 !tens_owner_t:
        private TensOwnerCtor
        private TensOwnerIsSet
        private TensOwnerGetTensHeader
        private TensOwnerGetOwnerId
        private tens_owner_dtor
 !tens_entry_t:
        private TensEntryCtor
        private TensEntryIsSet
        private TensEntryGetTensor
        private TensEntryGetOwnerList
        private tens_entry_dtor
 !tens_instr_t:
        private TensInstrCtor
        private TensInstrDecode
        private TensInstrEncode
        private tens_instr_dtor

!IMPLEMENTATION:
       contains
![tens_owner_t]=================================================
        subroutine TensOwnerCtor(this,tens_header,owner_id,ierr)
         implicit none
         class(tens_owner_t), intent(out):: this
         class(tens_header_t), intent(in), target:: tens_header
         integer(INTD), intent(in):: owner_id
         integer(INTD), intent(out), optional:: ierr
         integer(INTD):: errc

         errc=0
         if(owner_id.ge.0) then
          this%owner_id=owner_id
          this%tens_header=>tens_header
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOwnerCtor
!-----------------------------------------------------
        function TensOwnerIsSet(this,ierr) result(ans)
         implicit none
         logical:: ans
         class(tens_owner_t), intent(in):: this
         integer(INTD), intent(out), optional:: ierr
         integer(INTD):: errc

         errc=0
         ans=(this%owner_id.ge.0)
         if(present(ierr)) ierr=errc
         return
        end function TensOwnerIsSet
!---------------------------------------------------------------------
        function TensOwnerGetTensHeader(this,ierr) result(tens_header)
         implicit none
         class(tens_header_t), pointer:: tens_header
         class(tens_owner_t), intent(in):: this
         integer(INTD), intent(out), optional:: ierr
         integer(INTD):: errc

         errc=0
         if(this%owner_id.ge.0) then
          tens_header=>this%tens_header
         else
          tens_header=>NULL(); errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensOwnerGetTensHeader
!---------------------------------------------------------------
        function TensOwnerGetOwnerId(this,ierr) result(owner_id)
         implicit none
         integer(INTD):: owner_id
         class(tens_owner_t), intent(in):: this
         integer(INTD), intent(out), optional:: ierr
         integer(INTD):: errc

         errc=0; owner_id=this%owner_id
         if(owner_id.lt.0) errc=-1
         if(present(ierr)) ierr=errc
         return
        end function TensOwnerGetOwnerId
!---------------------------------------
        subroutine tens_owner_dtor(this)
         implicit none
         type(tens_owner_t):: this

         this%tens_header=>NULL()
         this%owner_id=-1
         return
        end subroutine tens_owner_dtor
![tens_entry_t]===================================
        subroutine TensEntryCtor(this,ierr,tensor)
!Constructs a tensor cache entry. If <tensor> is present,
!it will be assumed persistent and will be pointed to.
!If <tensor> is absent, it will be allocated here.
         implicit none
         class(tens_entry_t), intent(out):: this                    !out: tensor cache entry
         integer(INTD), intent(out), optional:: ierr                !out: error code
         class(tens_rcrsv_t), intent(in), target, optional:: tensor !in: persistent tensor
         integer(INTD):: errc

         errc=0
         if(present(tensor)) then
          this%tensor=>tensor; this%tens_alloc=.FALSE.
         else
          allocate(this%tensor,STAT=errc)
          if(errc.eq.0) then; this%tens_alloc=.TRUE.; else; errc=-1; endif
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensEntryCtor
!-----------------------------------------------------
        function TensEntryIsSet(this,ierr) result(ans)
         implicit none
         logical:: ans                               !out: answer
         class(tens_entry_t), intent(in):: this      !in: tensor cache entry
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         ans=associated(this%tensor)
         if(present(ierr)) ierr=errc
         return
        end function TensEntryIsSet
!------------------------------------------------------------
        function TensEntryGetTensor(this,ierr) result(tensor)
         implicit none
         class(tens_rcrsv_t), pointer:: tensor       !out: pointer to the tensor
         class(tens_entry_t), intent(in):: this      !in: tensor cache entry
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%is_set(errc)) then
          if(errc.eq.0) then
           tensor=>this%tensor
          else
           tensor=>NULL(); errc=-2
          endif
         else
          tensor=>NULL(); errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensEntryGetTensor
!-------------------------------------------------------------------
        function TensEntryGetOwnerList(this,ierr) result(owner_list)
         implicit none
         class(list_bi_t), pointer:: owner_list         !out: pointer to the owner list
         class(tens_entry_t), intent(in), target:: this !in: tensor cache entry
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD):: errc

         owner_list=>NULL()
         if(this%is_set(errc)) then
          if(errc.eq.0) then
           owner_list=>this%owner_list
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensEntryGetOwnerList
!----------------------------------------------
        subroutine tens_entry_dtor(this)
         implicit none
         type(tens_entry_t):: this
         type(list_iter_t):: lit
         integer(INTD):: errc

         if(this%tens_alloc.and.associated(this%tensor)) deallocate(this%tensor)
         this%tensor=>NULL(); this%tens_alloc=.FALSE.
         errc=lit%init(this%owner_list); if(errc.eq.GFC_SUCCESS) errc=lit%delete_all()
         errc=lit%release()
         return
        end subroutine tens_entry_dtor
![tens_instr_t]============================================
        subroutine TensInstrCtor(this,op_code,ierr,op_spec)
!Constructs a tensor instruction from a given tensor operation.
!The tensor instruction is a realization of the given tensor operation for TAVP.
         implicit none
         class(tens_instr_t), intent(inout):: this        !out: tensor instruction (must be empty on entrance)
         integer(INTD), intent(in):: op_code              !in: instruction code (see top)
         integer(INTD), intent(out), optional:: ierr      !out: error code
         class(*), intent(in), target, optional:: op_spec !in: operation specification
         integer(INTD):: errc

         if(this%is_empty(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
!Construct instruction fields:
           select case(op_code)
           case(TAVP_INSTR_NOOP)
           case(TAVP_INSTR_STOP)

           case(TAVP_INSTR_CREATE,TAVP_INSTR_DESTROY)

           case(TAVP_INSTR_CONTRACT)

           case default
            errc=-5 !invalid operation (or not implemented)
           end select
!Activate the instruction:
           if(errc.eq.0) then
            call this%set_code(op_code,errc)
            if(errc.eq.DSVP_SUCCESS) then
             call this%set_status(DS_INSTR_NEW,errc)
             if(errc.ne.DSVP_SUCCESS) then
              call this%set_status(DS_INSTR_RETIRED)
              call tens_instr_dtor(this)
              errc=-4
             endif
            else
             call this%set_status(DS_INSTR_RETIRED)
             call tens_instr_dtor(this)
             errc=-3
            endif
           else
            call this%set_status(DS_INSTR_RETIRED)
            call tens_instr_dtor(this)
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensInstrCtor
!---------------------------------------------------------
        subroutine TensInstrDecode(this,instr_packet,ierr)
!Decodes a tensor instruction from the bytecode packet.
         implicit none
         class(tens_instr_t), intent(inout):: this       !out: tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet !in: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine TensInstrDecode
!---------------------------------------------------------
        subroutine TensInstrEncode(this,instr_packet,ierr)
!Encodes a tensor instruction into the bytecode packet.
         implicit none
         class(tens_instr_t), intent(in):: this          !in: tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet !out: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine TensInstrEncode
!---------------------------------------
        subroutine tens_instr_dtor(this)
         implicit none
         type(tens_instr_t):: this
         integer(INTD):: sts,errc

         sts=this%get_status(errc)
         if((sts.eq.DS_INSTR_EMPTY.or.sts.eq.DS_INSTR_RETIRED).and.errc.eq.DSVP_SUCCESS) then
          call this%clean(errc)
          if(errc.ne.0) call quit(errc,'#ERROR(tavp_manager:tens_instr_dtor): Tensor instruction destruction failed!')
         else
          call quit(-1,'#ERROR(tavp_manager:tens_instr_dtor): Attempt to destroy an active TAVP instruction!')
         endif
         return
        end subroutine tens_instr_dtor

       end module tavp_manager
