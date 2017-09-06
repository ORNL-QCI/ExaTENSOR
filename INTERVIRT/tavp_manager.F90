!ExaTENSOR: TAVP-Manager (TAVP-MNG) implementation
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/09/06

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
        integer(INTD), private:: CONS_OUT=6 !default output device
        integer(INTD), private:: DEBUG=0    !debugging mode
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
!TYPES:
 !Tensor argument cache entry (TAVP-specific):
        type, extends(tens_cache_entry_t), private:: tens_entry_mng_t
         integer(INTD), private:: owner_id=-1                      !tensor owner id (non-negative)
         contains
          procedure, private:: TensEntryMngCtor                    !ctor
          generic, public:: tens_entry_mng_ctor=>TensEntryMngCtor
          procedure, public:: get_owner_id=>TensEntryMngGetOwnerId !returns the owner id
          final:: tens_entry_mng_dtor
        end type tens_entry_mng_t
 !Tensor operand (encapsulated tensor data processible by a specific TAVP):
        type, extends(ds_oprnd_t), private:: tens_oprnd_t
         class(tens_rcrsv_t), pointer, private:: tensor=>NULL() !non-owning pointer to a persistent recursive tensor
         integer(INTD), private:: owner_id=-1                   !tensor owner id (non-negative)
         contains
          procedure, private:: TensOprndCtor                    !ctor
          generic, public:: tens_oprnd_ctor=>TensOprndCtor
          procedure, public:: get_tensor=>TensOprndGetTensor    !returns a pointer to the tensor
          procedure, public:: get_owner_id=>TensOprndGetOwnerId !returns the tensor owner id
          procedure, public:: set_owner_id=>TensOprndSetOwnerId !sets the tensor owner id
          procedure, public:: is_remote=>TensOprndIsRemote      !returns TRUE if the tensor operand is remote
          procedure, public:: acquire_rsc=>TensOprndAcquireRsc  !explicitly acquires local resources for the tensor operand
          procedure, public:: prefetch=>TensOprndPrefetch       !starts prefetching the remote tensor operand (acquires local resources!)
          procedure, public:: upload=>TensOprndUpload           !starts uploading the tensor operand to its remote location
          procedure, public:: sync=>TensOprndSync               !synchronizes the currently pending communication on the tensor operand
          procedure, public:: release=>TensOprndRelease         !destroys the present local copy of the tensor operand (releases local resources!), but the operand stays defined
          procedure, public:: destruct=>TensOprndDestruct       !performs complete destruction back to an empty state
          final:: tens_oprnd_dtor
        end type tens_oprnd_t
 !Tensor instruction (realization of a tensor operation for a specific TAVP):
        type, extends(ds_instr_t), private:: tens_instr_t
         integer(INTD), private:: num_out_oprnds=0                    !number of the output tensor instruction operands
         integer(INTD), private:: out_oprnds(0:MAX_TENSOR_OPERANDS-1) !tensor instruction operands which are considered output (to sync on them for completion)
         contains
          procedure, private:: TensInstrCtor               !ctor: constructs a tensor instruction from the specification of a tensor operation
          generic, public:: tens_instr_ctor=>TensInstrCtor
          procedure, public:: encode=>TensInstrEncode      !encoding procedure: Packs the TAVP instruction into a raw byte packet
          final:: tens_instr_dtor                          !dtor
        end type tens_instr_t
!MODULE DATA:
 !TAVP-MNG microcode (static) table, set by dsvp.configure():
        type(ds_microcode_t), target, private:: microcode(0:TAVP_ISA_SIZE-1)
!VISIBILITY:
 !tens_entry_mng_t:
        private TensEntryMngCtor
        private TensEntryMngGetOwnerId
        public tens_entry_mng_dtor
        private tens_entry_mng_alloc
 !tens_oprnd_t:
        private TensOprndCtor
        private TensOprndGetTensor
        private TensOprndGetOwnerId
        private TensOprndSetOwnerId
        private TensOprndIsRemote
        private TensOprndAcquireRsc
        private TensOprndPrefetch
        private TensOprndUpload
        private TensOprndSync
        private TensOprndRelease
        private TensOprndDestruct
        public tens_oprnd_dtor
 !tens_instr_t:
        private TensInstrCtor
        private TensInstrEncode
        public tens_instr_dtor
!IMPLEMENTATION:
       contains
![non-member:Microcode Implementation]==============
        subroutine acquire_resource_dummy(this,ierr)
!Dummy procedure for acquiring no resource.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine acquire_resource_dummy
!---------------------------------------------------
        subroutine acquire_resource_basic(this,ierr)
!Acquires resource for each tensor instruction operand.
!If some resources cannot be acquired now, returns TRY_LATER.
!In that case, the successfully acquired resources will be kept,
!unless an error other than TRY_LATER occurred.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code, possibly TRY_LATER
         integer(INTD):: errc,ier,n
         class(ds_oprnd_t), pointer:: oprnd

         n=this%get_num_operands(errc)
         if(errc.eq.DSVP_SUCCESS) then
          do while(n.gt.0)
           oprnd=>this%get_operand(n-1,ier)
           if(ier.eq.DSVP_SUCCESS) then
            call oprnd%acquire_rsc(ier)
            if(ier.ne.0) then
             if(ier.eq.TRY_LATER) then
              errc=ier
             else
              errc=-3; exit
             endif
            endif
           else
            errc=-2; exit
           endif
           n=n-1
          enddo
          if(errc.ne.0.and.errc.ne.TRY_LATER) call this%release_resource(ier)
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine acquire_resource_basic
!-------------------------------------------------
        subroutine prefetch_input_dummy(this,ierr)
!Dummy procedure for prefetching no input.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine prefetch_input_dummy
!-------------------------------------------------
        subroutine prefetch_input_basic(this,ierr)
!Starts prefetching input tensor operands.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,ier,n
         class(ds_oprnd_t), pointer:: oprnd

         n=this%get_num_operands(errc)
         if(errc.eq.DSVP_SUCCESS) then
          do while(n.gt.1) !operand 0 is output
           n=n-1
           oprnd=>this%get_operand(n,ier)
           if(ier.eq.DSVP_SUCCESS) then
            call oprnd%prefetch(ier)
            if(ier.ne.0.and.errc.eq.0) errc=-3
           else
            errc=-2; exit
           endif
          enddo
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine prefetch_input_basic
!---------------------------------------------------------------
        function sync_prefetch_dummy(this,ierr,wait) result(res)
!Dummy procedure for syncing no input prefetch.
         implicit none
         logical:: res                               !out: TRUE if synchronized
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: WAIT or TEST (defaults to WAIT)
         integer(INTD):: errc

         errc=0; res=.TRUE.
         if(present(ierr)) ierr=errc
         return
        end function sync_prefetch_dummy
!---------------------------------------------------------------
        function sync_prefetch_basic(this,ierr,wait) result(res)
!Synchronization on the input prefetch, either WAIT or TEST.
         implicit none
         logical:: res                               !out: TRUE if synchronized (all input operands)
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: WAIT or TEST (defaults to WAIT)
         integer(INTD):: errc,ier,n
         class(ds_oprnd_t), pointer:: oprnd
         logical:: wt

         errc=0; res=.FALSE.
         wt=.TRUE.; if(present(wait)) wt=wait
         n=this%get_num_operands(errc)
         if(errc.eq.DSVP_SUCCESS) then
          res=.TRUE.
          do while(n.gt.1)
           n=n-1
           oprnd=>this%get_operand(n,ier)
           if(ier.eq.DSVP_SUCCESS) then
            res=res.and.oprnd%sync(ier,wt)
            if(ier.ne.0.and.errc.eq.0) then; res=.FALSE.; errc=-3; endif
           else
            res=.FALSE.; if(errc.eq.0) errc=-2
           endif
          enddo
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function sync_prefetch_basic
!------------------------------------------------
        subroutine upload_output_dummy(this,ierr)
!Dummy procedure for uploading no output.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine upload_output_dummy
!------------------------------------------------
        subroutine upload_output_basic(this,ierr)
!Starts uploading the output tensor operand 0.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         class(ds_oprnd_t), pointer:: oprnd

         oprnd=>this%get_operand(0,errc) !output operand
         if(errc.eq.DSVP_SUCCESS) then
          call oprnd%upload(errc); if(errc.ne.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine upload_output_basic
!-------------------------------------------------------------
        function sync_upload_dummy(this,ierr,wait) result(res)
!Dummy procedure for syncing no output upload.
         implicit none
         logical:: res                               !out: TRUE if synchronized
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: WAIT or TEST (defaults to WAIT)
         integer(INTD):: errc

         errc=0; res=.TRUE.
         if(present(ierr)) ierr=errc
         return
        end function sync_upload_dummy
!-------------------------------------------------------------
        function sync_upload_basic(this,ierr,wait) result(res)
!Synchronization on the output upload.
         implicit none
         logical:: res                               !out: TRUE if synchronized
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: WAIT or TEST (defaults to WAIT)
         integer(INTD):: errc
         class(ds_oprnd_t), pointer:: oprnd
         logical:: wt

         errc=0; res=.FALSE.
         wt=.TRUE.; if(present(wait)) wt=wait
         oprnd=>this%get_operand(0,errc)
         if(errc.eq.DSVP_SUCCESS) then
          res=oprnd%sync(errc,wt); if(errc.ne.0) then; res=.FALSE.; errc=-2; endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function sync_upload_basic
!---------------------------------------------------
        subroutine release_resource_dummy(this,ierr)
!Dummy procedure for releasing no resource.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine release_resource_dummy
!---------------------------------------------------
        subroutine release_resource_basic(this,ierr)
!Releases resources occupied by tensor instruction operands,
!but the operands stay defined.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,ier,n
         class(ds_oprnd_t), pointer:: oprnd

         n=this%get_num_operands(errc)
         if(errc.eq.DSVP_SUCCESS) then
          do while(n.gt.0)
           n=n-1
           oprnd=>this%get_operand(n,ier)
           if(ier.eq.DSVP_SUCCESS) then
            call oprnd%release(ier)
            if(ier.ne.0.and.errc.eq.0) errc=-5
           else
            if(errc.eq.0) errc=-4
           endif
          enddo
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine release_resource_basic
!----------------------------------------------------------------
        function sync_execution_dummy(this,ierr,wait) result(res)
!Dummy procedure for syncing no execution.
         implicit none
         logical:: res                               !out: TRUE if synchronized
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: WAIT or TEST (defaults to WAIT)
         integer(INTD):: errc

         errc=0; res=.TRUE.
         if(present(ierr)) ierr=errc
         return
        end function sync_execution_dummy
!----------------------------------------------------------------
        function sync_execution_basic(this,ierr,wait) result(res)
!Synchronization on the tensor instruction execution.
         implicit none
         logical:: res                               !out: TRUE if synchronized
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: WAIT or TEST (defaults to WAIT)
         integer(INTD):: errc,sts,ans
         logical:: wt

         errc=0; res=.FALSE.
         wt=.TRUE.; if(present(wait)) wt=wait
         select type(this)
         class is(tens_instr_t)
          !`Implement: Test/wait on all parts of the destination tensor having been computed (check the task handle)
         class default
          errc=-1
         end select
         if(present(ierr)) ierr=errc
         return
        end function sync_execution_basic
!------------------------------------------
        subroutine execute_dummy(this,ierr)
!Dummy procedure for executing nothing.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine execute_dummy
!--------------------------------------------------
        subroutine execute_tensor_create(this,ierr)
!Executes tensor creation.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         !`Implement
         if(present(ierr)) ierr=errc
         return
        end subroutine execute_tensor_create
!---------------------------------------------------
        subroutine execute_tensor_destroy(this,ierr)
!Executes tensor destruction.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         !`Implement
         if(present(ierr)) ierr=errc
         return
        end subroutine execute_tensor_destroy
!----------------------------------------------------
        subroutine execute_tensor_contract(this,ierr)
!Executes tensor contraction.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         !`Implement
         if(present(ierr)) ierr=errc
         return
        end subroutine execute_tensor_contract
!--------------------------------------
        subroutine init_microcode(ierr)
!Initializes the global TAVP microcode bindings table when configuring TAVP.
!There are three kinds of microcode bindings (suffix):
! 1. DUMMY: Do nothing.
! 2. BASIC: Default microcode.
! 3. SPECIAL: Specialized microcode.
         implicit none
         integer(INTD), intent(out), optional:: ierr
         integer(INTD):: errc

         errc=0
 !TENSOR CREATE:
         microcode(TAVP_INSTR_CREATE)%acquire_resource=>acquire_resource_basic
         microcode(TAVP_INSTR_CREATE)%prefetch_input=>prefetch_input_dummy
         microcode(TAVP_INSTR_CREATE)%sync_prefetch=>sync_prefetch_dummy
         microcode(TAVP_INSTR_CREATE)%execute=>execute_tensor_create
         microcode(TAVP_INSTR_CREATE)%sync_execution=>sync_execution_basic
         microcode(TAVP_INSTR_CREATE)%upload_output=>upload_output_dummy
         microcode(TAVP_INSTR_CREATE)%sync_upload=>sync_upload_dummy
         microcode(TAVP_INSTR_CREATE)%release_resource=>release_resource_dummy
 !TENSOR DESTROY:
         microcode(TAVP_INSTR_DESTROY)%acquire_resource=>acquire_resource_dummy
         microcode(TAVP_INSTR_DESTROY)%prefetch_input=>prefetch_input_dummy
         microcode(TAVP_INSTR_DESTROY)%sync_prefetch=>sync_prefetch_dummy
         microcode(TAVP_INSTR_DESTROY)%execute=>execute_tensor_destroy
         microcode(TAVP_INSTR_DESTROY)%sync_execution=>sync_execution_basic
         microcode(TAVP_INSTR_DESTROY)%upload_output=>upload_output_dummy
         microcode(TAVP_INSTR_DESTROY)%sync_upload=>sync_upload_dummy
         microcode(TAVP_INSTR_DESTROY)%release_resource=>release_resource_basic
 !TENSOR CONTRACT:
         microcode(TAVP_INSTR_CONTRACT)%acquire_resource=>acquire_resource_basic
         microcode(TAVP_INSTR_CONTRACT)%prefetch_input=>prefetch_input_basic
         microcode(TAVP_INSTR_CONTRACT)%sync_prefetch=>sync_prefetch_basic
         microcode(TAVP_INSTR_CONTRACT)%execute=>execute_tensor_contract
         microcode(TAVP_INSTR_CONTRACT)%sync_execution=>sync_execution_basic
         microcode(TAVP_INSTR_CONTRACT)%upload_output=>upload_output_basic
         microcode(TAVP_INSTR_CONTRACT)%sync_upload=>sync_upload_basic
         microcode(TAVP_INSTR_CONTRACT)%release_resource=>release_resource_basic
 !DONE
         if(present(ierr)) ierr=errc
         return
        end subroutine init_microcode
![tens_entry_mng_t]========================================
        subroutine TensEntryMngCtor(this,tensor,owner,ierr)
!Constructs a <tens_entry_mng_t>. Note move semantics for <tensor>!
         implicit none
         class(tens_entry_mng_t), intent(out):: this          !out: specialized tensor cache entry
         class(tens_rcrsv_t), pointer, intent(inout):: tensor !inout: pointer to an allocated tensor (ownership transfer will occur here!)
         integer(INTD), intent(in):: owner                    !in: tensor owner id
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTD):: errc

         errc=0
         if(associated(tensor)) then
          if(owner.ge.0) then
           call this%set_tensor(tensor,errc)
           if(errc.eq.0) then
            tensor=>NULL() !transfer the ownership
            this%owner_id=owner
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensEntryMngCtor
!------------------------------------------------------------
        function TensEntryMngGetOwnerId(this,ierr) result(id)
         implicit none
         integer(INTD):: id                          !out: tensor owner id
         class(tens_entry_mng_t), intent(in):: this  !in: specialized tensor cache entry
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%is_set(errc)) then
          id=this%owner_id
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensEntryMngGetOwnerId
!-------------------------------------------
        subroutine tens_entry_mng_dtor(this)
         implicit none
         type(tens_entry_mng_t):: this !inout: specialized tensor cache entry

         call this%nullify_tensor(.TRUE.) !deallocate the tensor component (if it is set)
         this%owner_id=-1
         return
        end subroutine tens_entry_mng_dtor
!-------------------------------------------------------------
        function tens_entry_mng_alloc(tens_entry) result(ierr)
!Non-member allocator for tens_entry_mng_t.
         implicit none
         integer(INTD):: ierr
         class(tens_cache_entry_t), allocatable, intent(out):: tens_entry

         allocate(tens_entry_mng_t::tens_entry,STAT=ierr)
         return
        end function tens_entry_mng_alloc
![tens_oprnd_t]=========================================
        subroutine TensOprndCtor(this,tensor,ierr,owner)
!Constructs a tensor operand. The <tensor> must be set.
!The owner id is optional.
         implicit none
         class(tens_oprnd_t), intent(inout):: this        !inout: undefined tensor operand (on entrance)
         class(tens_rcrsv_t), target, intent(in):: tensor !in: defined tensor
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD), intent(in), optional:: owner      !in: tensor owner id (no restrictions)
         integer(INTD):: errc

         if(.not.this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(tensor%is_set(errc)) then
            if(errc.eq.TEREC_SUCCESS) then
             this%tensor=>tensor
             if(present(owner)) this%owner_id=owner
             call this%mark_active(errc)
             if(errc.ne.DSVP_SUCCESS) errc=-5
            else
             errc=-4
            endif
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndCtor
!------------------------------------------------------------
        function TensOprndGetTensor(this,ierr) result(tens_p)
!Returns a pointer to the tensor.
         implicit none
         class(tens_rcrsv_t), pointer:: tens_p       !out: pointer to the tensor
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; tens_p=>this%tensor
         if(.not.associated(tens_p)) errc=-1
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetTensor
!---------------------------------------------------------
        function TensOprndGetOwnerId(this,ierr) result(id)
!Returns the tensor owner id.
         implicit none
         integer(INTD):: id                          !out: tensor owner id
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(associated(this%tensor)) then
          id=this%owner_id
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetOwnerId
!------------------------------------------------------
        subroutine TensOprndSetOwnerId(this,owner,ierr)
!Sets the tensor owner id (must be non-negative).
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: active tensor operand
         integer(INTD), intent(in):: owner           !in: tensor owner id (>=0)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(owner.ge.0) then
            this%owner_id=owner
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndSetOwnerId
!--------------------------------------------------------
        function TensOprndIsRemote(this,ierr) result(res)
!Returns TRUE if the tensor operand is remote, FALSE otherwise.
         implicit none
         logical:: res                               !out: result
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,id

         res=.TRUE. !assume remote by default
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           id=this%get_owner_id(errc)
           if(errc.eq.0) then
            res=id.ne.role_rank !tensor is owned by a different TAVP-MNG processor
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensOprndIsRemote
!------------------------------------------------
        subroutine TensOprndAcquireRsc(this,ierr)
!Acquires local resources for the remote tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         !`No local resources are currently needed
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndAcquireRsc
!----------------------------------------------
        subroutine TensOprndPrefetch(this,ierr)
!Starts prefetching the (remote) tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(.not.this%is_present(errc)) then
            if(errc.eq.DSVP_SUCCESS) then
             if(this%get_comm_stat().eq.DS_OPRND_NO_COMM) then
              call this%set_comm_stat(DS_OPRND_FETCHING,errc); if(errc.ne.DSVP_SUCCESS) errc=-6
             else
              errc=-5
             endif
            else
             errc=-4
            endif
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndPrefetch
!----------------------------------------------
        subroutine TensOprndUpload(this,ierr)
!Starts uploading the (remote) tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%is_present(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(this%get_comm_stat().eq.DS_OPRND_NO_COMM) then
            call this%set_comm_stat(DS_OPRND_UPLOADING,errc); if(errc.ne.DSVP_SUCCESS) errc=-4
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndUpload
!---------------------------------------------------------
        function TensOprndSync(this,ierr,wait) result(res)
!Synchronizes a pending prefetch/upload.
         implicit none
         logical:: res                               !out: TRUE on communication completion, FALSE otherwise
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: TRUE activates WAIT instead of TEST synchronization (default)
         integer(INTD):: errc

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(this%get_comm_stat().ne.DS_OPRND_NO_COMM) then
            call this%mark_delivered(errc); if(errc.ne.DSVP_SUCCESS) errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensOprndSync
!---------------------------------------------
        subroutine TensOprndRelease(this,ierr)
!Releases local resources acquired for the remote tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         logical:: delivered

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(this%get_comm_stat().ne.DS_OPRND_NO_COMM) then
            delivered=this%sync(errc,wait=.TRUE.)
            if((.not.delivered).or.(errc.ne.0)) errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         !`No local resources are currently needed
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndRelease
!----------------------------------------------
        subroutine TensOprndDestruct(this,ierr)
!Destructs the tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,ier
         logical:: delivered

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(this%get_comm_stat().ne.DS_OPRND_NO_COMM) then
            delivered=this%sync(errc,wait=.TRUE.)
            if((.not.delivered).or.(errc.ne.0)) errc=-4
           endif
           call this%mark_empty(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-3
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndDestruct
!---------------------------------------
        subroutine tens_oprnd_dtor(this)
         implicit none
         type(tens_oprnd_t):: this
         integer(INTD):: errc

         call this%destruct(errc)
         if(errc.ne.0) call quit(errc,'#FATAL(TAVP-MNG:tens_oprnd_dtor): Destructor failed!')
        end subroutine tens_oprnd_dtor
![tens_instr_t]============================================
        subroutine TensInstrCtor(this,op_code,ierr,op_spec)
!Constructs a tensor instruction from a given tensor operation.
!The tensor instruction is a realization of a given tensor operation
!for a specific TAVP kind.
         implicit none
         class(tens_instr_t), intent(inout):: this        !out: tensor instruction (must be empty on entrance)
         integer(INTD), intent(in):: op_code              !in: instruction code (see top of this module)
         integer(INTD), intent(out), optional:: ierr      !out: error code
         class(*), intent(in), target, optional:: op_spec !in: formal operation specification
         integer(INTD):: errc

         if(this%is_empty(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
!Construct the instruction:
           select case(op_code)
           case(TAVP_INSTR_NOOP)
           case(TAVP_INSTR_STOP)
           case(TAVP_INSTR_CREATE,TAVP_INSTR_DESTROY)
            call construct_instr_create_destroy(errc); if(errc.ne.0) errc=-6
           case(TAVP_INSTR_CONTRACT)
            call construct_instr_contract(errc); if(errc.ne.0) errc=-5
           case default
            errc=-4 !invalid instruction opcode (or not implemented)
           end select
!Activate the instruction:
           if(errc.eq.0) then
            call this%activate(op_code,microcode(op_code),errc); if(errc.ne.0) errc=-3
           else
            call this%set_status(DS_INSTR_RETIRED,errc,TAVP_ERR_GEN_FAILURE)
           endif
           if(errc.ne.0) call tens_instr_dtor(this)
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return

        contains

         subroutine construct_instr_create_destroy(jerr)
          !CREATE/DESTROY a tensor:
          !op_spec={tens_rcrsv_t}
          integer(INTD), intent(out):: jerr
          class(ds_oprnd_t), pointer:: oprnd
          class(tens_oprnd_t), pointer:: tens_oprnd
          class(tens_rcrsv_t), pointer:: tensor

          jerr=0; tensor=>NULL()
          select type(op_spec); class is(tens_rcrsv_t); tensor=>op_spec; end select
          if(associated(tensor)) then
           if(tensor%is_set()) then
            call this%alloc_operands(1,jerr)
            if(jerr.eq.DSVP_SUCCESS) then
             allocate(tens_oprnd,STAT=jerr)
             if(jerr.eq.0) then
              call tens_oprnd%tens_oprnd_ctor(tensor,jerr) !`tensor owner is omitted here
              if(jerr.eq.0) then
               oprnd=>tens_oprnd
               call this%set_operand(0,oprnd,jerr); if(jerr.ne.DSVP_SUCCESS) jerr=-6
               oprnd=>NULL() !<oprnd> pointer was saved in the tensor instruction and will later be deallocated
               this%num_out_oprnds=1; this%out_oprnds(0:this%num_out_oprnds-1)=(/0/)
              else
               jerr=-5
              endif
              tens_oprnd=>NULL()
             else
              jerr=-4
             endif
            else
             jerr=-3
            endif
           else
            jerr=-2
           endif
           tensor=>NULL() !<tensor> pointed to an external object
          else
           jerr=-1
          endif
          return
         end subroutine construct_instr_create_destroy

         subroutine construct_instr_contract(jerr)
          !CONTRACT two tensors into another tensor:
          !op_spec={tens_contraction_t}
          integer(INTD), intent(out):: jerr
          integer(INTD):: jj
          class(ds_oprnd_t), pointer:: oprnd
          class(ds_instr_ctrl_t), pointer:: instr_ctrl
          class(contr_ptrn_ext_t), pointer:: contr_ptrn
          class(tens_contraction_t), pointer:: tens_contr
          class(ctrl_tens_contr_t), pointer:: tens_contr_ctrl
          class(tens_rcrsv_t), pointer:: tensor
          class(tens_oprnd_t), pointer:: tens_oprnd

          jerr=0; tens_contr=>NULL()
          select type(op_spec); class is(tens_contraction_t); tens_contr=>op_spec; end select
          if(associated(tens_contr)) then
           if(tens_contr%is_set()) then
            contr_ptrn=>tens_contr%get_ext_contr_ptrn(jerr)
            if(jerr.eq.TEREC_SUCCESS) then
             allocate(tens_contr_ctrl,STAT=jerr)
             if(jerr.eq.0) then
              call tens_contr_ctrl%ctrl_tens_contr_ctor(contr_ptrn,jerr,tens_contr%get_prefactor()) !contraction pattern is cloned by value
              if(jerr.eq.0) then
               instr_ctrl=>tens_contr_ctrl
               call this%set_control(instr_ctrl,jerr) !ownership transfer for instr_ctrl=tens_contr_ctrl
               if(jerr.eq.DSVP_SUCCESS) then
                call this%alloc_operands(3,jerr)
                if(jerr.eq.DSVP_SUCCESS) then
                 do jj=0,2
                  tensor=>tens_contr%get_argument(jj,jerr); if(jerr.ne.TEREC_SUCCESS) exit
                  allocate(tens_oprnd,STAT=jerr); if(jerr.ne.0) exit
                  call tens_oprnd%tens_oprnd_ctor(tensor,jerr); if(jerr.ne.0) exit !`tensor owner is omitted here
                  oprnd=>tens_oprnd; call this%set_operand(jj,oprnd,jerr); if(jerr.ne.DSVP_SUCCESS) exit !ownership transfer for oprnd=tens_oprnd
                  tensor=>NULL(); tens_oprnd=>NULL(); oprnd=>NULL() !<oprnd> pointer was saved in the tensor instruction and will later be deallocated
                 enddo
                 this%num_out_oprnds=1; this%out_oprnds(0:this%num_out_oprnds-1)=(/0/)
                else
                 jerr=-7
                endif
               else
                jerr=-6
               endif
              else
               jerr=-5
              endif
              instr_ctrl=>NULL(); tens_contr_ctrl=>NULL() !<tens_contr_ctrl> pointer was saved in the tensor instruction and will later be deallocated
             else
              jerr=-4
             endif
            else
             jerr=-3
            endif
            contr_ptrn=>NULL()
           else
            jerr=-2
           endif
           tens_contr=>NULL()
          else
           jerr=-1
          endif
          return
         end subroutine construct_instr_contract

        end subroutine TensInstrCtor
!---------------------------------------------------------
        subroutine TensInstrEncode(this,instr_packet,ierr)
!Encodes a tensor instruction into the bytecode packet.
         implicit none
         class(tens_instr_t), intent(in):: this          !in: defined tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet !out: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,op_code

!Pack the instruction code (op_code):
         if(.not.this%is_empty(errc)) then
          op_code=this%get_code(errc)
          if(errc.eq.DSVP_SUCCESS) then
           call pack_builtin(instr_packet,op_code,errc)
           if(errc.eq.0) then
!Pack the instruction body:
            select case(op_code)
            case(TAVP_INSTR_NOOP)
            case(TAVP_INSTR_STOP)
            case(TAVP_INSTR_CREATE,TAVP_INSTR_DESTROY)
             call encode_instr_create_destroy(errc); if(errc.ne.0) errc=-6
            case(TAVP_INSTR_CONTRACT)
             call encode_instr_contract(errc); if(errc.ne.0) errc=-5
            case default
             errc=-4
            end select
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return

        contains

         subroutine encode_instr_create_destroy(jerr)
          !CREATE/DESTROY a tensor:
          !Packet format: {op_code|tensor}
          integer(INTD), intent(out):: jerr
          class(ds_oprnd_t), pointer:: oprnd
          class(tens_rcrsv_t), pointer:: tensor

          jerr=0
          oprnd=>this%get_operand(0,jerr)
          if(jerr.eq.DSVP_SUCCESS) then
           select type(oprnd)
           class is(tens_oprnd_t)
            tensor=>oprnd%get_tensor(jerr)
            if(jerr.eq.0) then
             if(tensor%is_set()) then
              call tensor%pack(instr_packet,jerr)
              if(jerr.ne.0) jerr=-5
             else
              jerr=-4
             endif
             tensor=>NULL()
            else
             jerr=-3
            endif
           class default
            jerr=-2
           end select
           oprnd=>NULL()
          else
           jerr=-1
          endif
          return
         end subroutine encode_instr_create_destroy

         subroutine encode_instr_contract(jerr)
          !CONTRACT two tensors: tensor0+=tensor1*tensor2*scalar:
          !Packed format: {op_code|ctrl_tens_contr_t|tensor0,tensor1,tensor2}
          integer(INTD), intent(out):: jerr
          integer(INTD):: jj
          class(ds_oprnd_t), pointer:: oprnd
          class(tens_rcrsv_t), pointer:: tensor
          class(ds_instr_ctrl_t), pointer:: tens_contr_ctrl

          jerr=0
          tens_contr_ctrl=>this%get_control(jerr)
          if(jerr.eq.DSVP_SUCCESS) then
           select type(tens_contr_ctrl)
           class is(ctrl_tens_contr_t)
            call tens_contr_ctrl%pack(instr_packet,jerr); if(jerr.ne.0) jerr=-8
           class default
            jerr=-7
           end select
           if(jerr.eq.0) then
            do jj=0,2
             oprnd=>this%get_operand(jj,jerr); if(jerr.ne.DSVP_SUCCESS) then; jerr=-6; exit; endif
             select type(oprnd)
             class is(tens_oprnd_t)
              tensor=>oprnd%get_tensor(jerr); if(jerr.ne.0) then; jerr=-5; exit; endif
              if(.not.tensor%is_set()) then; jerr=-4; exit; endif !trap
              call tensor%pack(instr_packet,jerr); if(jerr.ne.0) then; jerr=-3; exit; endif
             class default
              jerr=-2; exit
             end select
            enddo
            tensor=>NULL(); oprnd=>NULL()
           endif
           tens_contr_ctrl=>NULL()
          else
           jerr=-1
          endif
          return
         end subroutine encode_instr_contract

        end subroutine TensInstrEncode
!---------------------------------------
        subroutine tens_instr_dtor(this)
         implicit none
         type(tens_instr_t):: this !inout: empty or retired tensor instruction
         integer(INTD):: sts,errc

         sts=this%get_status(errc)
         if((sts.eq.DS_INSTR_EMPTY.or.sts.eq.DS_INSTR_RETIRED).and.errc.eq.DSVP_SUCCESS) then
          call this%clean(errc)
          if(errc.ne.DSVP_SUCCESS) call quit(errc,'#FATAL(TAVP-MNG:tens_instr_dtor): Tensor instruction destruction failed!')
         else
          call quit(-1,'#FATAL(TAVP-MNG:tens_instr_dtor): Attempt to destroy an active TAVP instruction!')
         endif
         return
        end subroutine tens_instr_dtor

       end module tavp_manager
