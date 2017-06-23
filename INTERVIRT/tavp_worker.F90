!ExaTENSOR: TAVP Worker
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/06/23

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

       module tavp_worker
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
 !Elementary tensor instruction granularity classification:
        real(8), public:: EXA_FLOPS_MEDIUM=1d10 !minimal number of Flops to consider the operation as medium-cost
        real(8), public:: EXA_FLOPS_HEAVY=1d12  !minimal number of Flops to consider the operation as heavy-cost
        real(8), public:: EXA_COST_TO_SIZE=1d2  !minimal cost (Flops) to size (Words) ratio to consider the operation arithmetically intensive
!TYPES:
 !Tensor cache entry:
        type, private:: tens_entry_t
         class(tens_rcrsv_t), pointer, private:: tensor=>NULL()   !tensor (either allocated or associated)
         class(tens_resrc_t), pointer, private:: resource=>NULL() !resource (either allocated or associated)
         logical, private:: tens_alloc=.FALSE.                    !TRUE if the tensor pointer is allocated
         logical, private:: res_alloc=.FALSE.                     !TRUE if the resource pointer is allocated
         contains
          procedure, private:: TensEntryCtor                      !ctor
          generic, public:: tens_entry_ctor=>TensEntryCtor
          procedure, public:: is_set=>TensEntryIsSet              !returns TRUE if the tensor cache entry is set
          procedure, public:: get_tensor=>TensEntryGetTensor      !returns a pointer to the tensor
          procedure, public:: get_resource=>TensEntryGetResource  !returns a pointer to the resource
          final:: tens_entry_dtor                                 !dtor
        end type tens_entry_t
 !Tensor cache:
        type, private:: tens_cache_t
         type(dictionary_t), private:: map                        !cache dictionary
         contains
          procedure, public:: lookup=>TensCacheLookup             !looks up a given tensor in the cache
          procedure, public:: store=>TensCacheStore               !stores a given tensor in the cache
          procedure, public:: evict=>TensCacheEvict               !evicts a given tensor from the cache
          procedure, public:: erase=>TensCacheErase               !erases everything from the cache (regardless of MPI communications)
        end type tens_cache_t
 !Tensor instruction (realization of a tensor operation for a specific TAVP):
        type, extends(ds_instr_t), private:: tens_instr_t
         contains
          procedure, private:: TensInstrCtor                        !ctor: constructs a tensor instruction from the specification of a tensor operation
          generic, public:: tens_instr_ctor=>TensInstrCtor
          procedure, public:: decode=>TensInstrDecode               !decoding procedure: Unpacks the raw byte packet (bytecode) and constructs a TAVP instruction
          procedure, public:: encode=>TensInstrEncode               !encoding procedure: Packs the TAVP instruction into a raw byte packet (bytecode)
          procedure, private:: set_microcode=>TensInstrSetMicrocode !sets up instruction dynamic bindings to the corresponding microcode
          procedure, private:: activate=>TensInstrActivate          !activates the instruction
          final:: tens_instr_dtor                                   !dtor
        end type tens_instr_t
 !TAVP instruction microcode binding (set by the TAVP initialization):
        type, private:: microcode_bind_t
         procedure(ds_instr_self_i), nopass, pointer:: acquire_resource=>NULL() !acquires local resources for instruction operands
         procedure(ds_instr_self_i), nopass, pointer:: prefetch_input=>NULL()   !starts prefetching input operands
         procedure(ds_instr_self_i), nopass, pointer:: sync_prefetch=>NULL()    !synchronizes the input prefetch (either test or wait)
         procedure(ds_instr_self_i), nopass, pointer:: execute=>NULL()          !executes the domain-specific instruction
         procedure(ds_instr_self_i), nopass, pointer:: sync_execution=>NULL()   !synchronizes the execution (either test or wait)
         procedure(ds_instr_self_i), nopass, pointer:: upload_output=>NULL()    !starts uploading the output
         procedure(ds_instr_self_i), nopass, pointer:: sync_upload=>NULL()      !synchronizes the output upload (either test or wait)
         procedure(ds_instr_self_i), nopass, pointer:: release_resource=>NULL() !releases local resources occupied by instruction operands
        end type microcode_bind_t
!INTERFACES:
!DATA:
 !Tensor cache (both persistent and temporary tensors):
        type(tens_cache_t), private:: tens_cache
 !TAVP instruction microcode bindings (set by the TAVP initialization):
        type(microcode_bind_t), private:: microcode(0:TAVP_ISA_SIZE-1)
!VISIBILITY:
 !tens_entry_t:
        private TensEntryCtor
        private TensEntryIsSet
        private TensEntryGetTensor
        private TensEntryGetResource
        private tens_entry_dtor
 !tens_cache_t:
        private TensCacheLookup
        private TensCacheStore
        private TensCacheEvict
        private TensCacheErase
 !tens_instr_t:
        private TensInstrCtor
        private TensInstrDecode
        private TensInstrEncode
        private TensInstrSetMicrocode
        private TensInstrActivate
        private tens_instr_dtor

!IMPLEMENTATION:
       contains
![tens_entry_t]============================================
        subroutine TensEntryCtor(this,ierr,tensor,resource)
!CTOR: If either <tensor> or <resource> are not present,
!the corresponding components of <this> will be allocated,
!otherwise pointer associated.
         implicit none
         class(tens_entry_t), intent(out):: this                       !out: tensor cache entry
         integer(INTD), intent(out), optional:: ierr                   !out: error code
         class(tens_rcrsv_t), intent(in), pointer, optional:: tensor   !in: pointer to the tensor
         class(tens_resrc_t), intent(in), pointer, optional:: resource !in: pointer to the resource
         integer(INTD):: errc

         errc=0
         if(present(tensor)) then
          if(associated(tensor)) then
           this%tensor=>tensor; this%tens_alloc=.FALSE.
          else
           errc=-1
          endif
         else
          allocate(this%tensor,STAT=errc)
          if(errc.eq.0) then; this%tens_alloc=.TRUE.; else; errc=-2; endif
         endif
         if(errc.eq.0) then
          if(present(resource)) then
           if(associated(resource)) then
            this%resource=>resource; this%res_alloc=.FALSE.
           else
            errc=-3
           endif
          else
           allocate(this%resource,STAT=errc)
           if(errc.eq.0) then; this%res_alloc=.TRUE.; else; errc=-4; endif
          endif
         endif
         if(errc.ne.0) call tens_entry_dtor(this)
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
!--------------------------------------------------------------
        function TensEntryGetTensor(this,ierr) result(tensor_p)
         implicit none
         class(tens_rcrsv_t), pointer:: tensor_p     !out: pointer to the tensor
         class(tens_entry_t), intent(in):: this      !in: tensor cache entry
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         tensor_p=>this%tensor
         if(.not.associated(tensor_p)) errc=-1
         if(present(ierr)) ierr=errc
         return
        end function TensEntryGetTensor
!------------------------------------------------------------------
        function TensEntryGetResource(this,ierr) result(resource_p)
         implicit none
         class(tens_resrc_t), pointer:: resource_p   !out: pointer to the resource
         class(tens_entry_t), intent(in):: this      !in: tensor cache entry
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         resource_p=>this%resource
         if(.not.associated(resource_p)) errc=-1
         if(present(ierr)) ierr=errc
         return
        end function TensEntryGetResource
!---------------------------------------
        subroutine tens_entry_dtor(this)
         implicit none
         type(tens_entry_t):: this

         if(associated(this%tensor).and.this%tens_alloc) deallocate(this%tensor)
         if(associated(this%resource).and.this%res_alloc) deallocate(this%resource)
         this%tensor=>NULL(); this%resource=>NULL()
         this%tens_alloc=.FALSE.; this%res_alloc=.FALSE.
         return
        end subroutine tens_entry_dtor
![tens_cache_t]================================================================
        function TensCacheLookup(this,tens_signature,ierr) result(tens_entry_p)
!Looks up a given tensor signature in the tensor cache. If found, returns
!a pointer to the corresponding tensor cache entry. If not found, returns NULL.
         implicit none
         class(tens_entry_t), pointer:: tens_entry_p                  !out: pointer to the tensor cache entry or NULL
         class(tens_cache_t), intent(in):: this                       !in: tensor cache
         class(tens_signature_t), intent(in), target:: tens_signature !in: tensor signature to look up
         integer(INTD), intent(out), optional:: ierr                  !out: error code
         integer(INTD):: errc,res
         class(*), pointer:: uptr
         type(dictionary_iter_t):: dit

         tens_entry_p=>NULL()
         errc=dit%init(this%map)
         if(errc.eq.GFC_SUCCESS) then
          uptr=>NULL()
          res=dit%search(GFC_DICT_FETCH_IF_FOUND,cmp_tens_signatures,tens_signature,value_out=uptr)
          if(res.eq.GFC_FOUND) then
           select type(uptr); class is(tens_entry_t); tens_entry_p=>uptr; end select
           if(.not.associated(tens_entry_p)) errc=-4
          else
           if(res.ne.GFC_NOT_FOUND) errc=-3
          endif
          res=dit%release(); if(res.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensCacheLookup
!-------------------------------------------------------------------------------------
        function TensCacheStore(this,tensor,ierr,resource,tens_entry_p) result(stored)
!Stores a tensor (and optionally its associated resource) in the tensor cache,
!unless the tensor is already present in the tensor cache. In any case, an optional
!<tens_entry_p> will point to either existing or newly created tensor cache entry.
!If the resource is absent, and empty resource will be allocated for the newly created
!tensor cache entry. <stored> is TRUE on successful new storage, FALSE otherwise.
         implicit none
         logical:: stored                                                   !out: TRUE on successful new store, FALSE otherwise
         class(tens_cache_t), intent(inout):: this                          !inout: tensor cache
         class(tens_rcrsv_t), pointer, intent(in):: tensor                  !in: tensor
         integer(INTD), intent(out), optional:: ierr                        !out: error code
         class(tens_resrc_t), pointer, intent(in), optional:: resource      !in: resource associated with the tensor
         class(tens_entry_t), pointer, intent(out), optional:: tens_entry_p !out: tensor cache entry (new or existing)
         integer(INTD):: errc,res
         class(*), pointer:: uptr
         class(tens_header_t), pointer:: thp
         class(tens_signature_t), pointer:: tsp
         type(dictionary_iter_t):: dit
         type(tens_entry_t):: tens_entry

         stored=.FALSE.
         errc=dit%init(this%map)
         if(errc.eq.GFC_SUCCESS) then
          if(associated(tensor)) then
           if(tensor%is_set()) then
            thp=>tensor%get_header(errc)
            if(errc.eq.TEREC_SUCCESS) then
             tsp=>thp%get_signature(errc)
             if(errc.eq.TEREC_SUCCESS) then
              if(present(resource)) then
               call tens_entry%tens_entry_ctor(errc,tensor,resource)
              else
               call tens_entry%tens_entry_ctor(errc,tensor)
              endif
              if(errc.eq.0) then
               uptr=>NULL()
               res=dit%search(GFC_DICT_ADD_IF_NOT_FOUND,cmp_tens_signatures,tsp,tens_entry,GFC_BY_VAL,value_out=uptr)
               if(res.eq.GFC_NOT_FOUND) then; stored=.TRUE.; else; if(res.ne.GFC_FOUND) errc=-9; endif
               if(present(tens_entry_p).and.errc.eq.0) then
                tens_entry_p=>NULL()
                select type(uptr); class is(tens_entry_t); tens_entry_p=>uptr; end select
                if(.not.associated(tens_entry_p)) errc=-8
               endif
              else
               errc=-7
              endif
             else
              errc=-6
             endif
            else
             errc=-5
            endif
           else
            errc=-4
           endif
          else
           errc=-3
          endif
          res=dit%release(); if(res.ne.0.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensCacheStore
!----------------------------------------------------------------
        function TensCacheEvict(this,tensor,ierr) result(evicted)
!Evicts a specific tensor cache entry from the tensor cache.
         implicit none
         logical:: evicted                                !TRUE if the tensor cache entry was found and evicted, FALSE otherwise
         class(tens_cache_t), intent(inout):: this        !inout: tensor cache
         class(tens_rcrsv_t), intent(in), target:: tensor !in: tensor to find and evict from the cache
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD):: errc,res
         class(tens_header_t), pointer:: thp
         class(tens_signature_t), pointer:: tsp
         type(dictionary_iter_t):: dit

         evicted=.FALSE.
         errc=dit%init(this%map)
         if(errc.eq.GFC_SUCCESS) then
          if(tensor%is_set()) then
           thp=>tensor%get_header(errc)
           if(errc.eq.TEREC_SUCCESS) then
            tsp=>thp%get_signature(errc)
            if(errc.eq.TEREC_SUCCESS) then
             res=dit%search(GFC_DICT_DELETE_IF_FOUND,cmp_tens_signatures,tsp)
             if(res.eq.GFC_FOUND) then; evicted=.TRUE.; else; if(res.ne.GFC_NOT_FOUND) errc=-6; endif
            else
             errc=-5
            endif
           else
            errc=-4
           endif
          else
           errc=-3
          endif
          res=dit%release(); if(res.ne.0.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensCacheEvict
!-------------------------------------------
        subroutine TensCacheErase(this,ierr)
!Erases the tensor cache completely and unconditionally.
         implicit none
         class(tens_cache_t), intent(inout):: this   !inout: tensor cache
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,res
         type(dictionary_iter_t):: dit

         errc=dit%init(this%map)
         if(errc.eq.GFC_SUCCESS) then
          call dit%delete_all(errc); if(errc.ne.GFC_SUCCESS) errc=-3
          res=dit%release(); if(res.ne.0.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCacheErase
![tens_instr_t]============================================
        subroutine TensInstrCtor(this,op_code,ierr,op_spec)
!Constructs a tensor instruction from a given tensor operation.
!The tensor instruction is a realization of a given tensor operation
!for a specific TAVP kind.
         implicit none
         class(tens_instr_t), intent(inout):: this        !out: tensor instruction (must be empty on entrance)
         integer(INTD), intent(in):: op_code              !in: instruction code (see top)
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
            call construct_instr_create(errc); if(errc.ne.0) errc=-6
           case(TAVP_INSTR_CONTRACT)
            call construct_instr_contract(errc); if(errc.ne.0) errc=-5
           case default
            errc=-4 !invalid operation (or not implemented)
           end select
!Activate the instruction:
           if(errc.eq.0) then
            call this%activate(op_code,errc); if(errc.ne.0) errc=-3
           else
            call this%set_status(DS_INSTR_RETIRED,errc,-1)
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

        contains

         subroutine construct_instr_create(jerr)
          !op_spec={tens_rcrsv_t}
          integer(INTD), intent(out):: jerr
          class(tens_oprnd_t), pointer:: oprnd
          class(tens_rcrsv_t), pointer:: tensor

          jerr=0
          tensor=>NULL()
          select type(op_spec); class is(tens_rcrsv_t); tensor=>op_spec; end select
          if(associated(tensor)) then
           if(tensor%is_set()) then
            call this%alloc_operands(1,jerr)
            if(jerr.eq.DSVP_SUCCESS) then
             allocate(oprnd,STAT=jerr)
             if(jerr.eq.0) then
              call oprnd%tens_oprnd_ctor(tensor,jerr)
              if(jerr.eq.0) then
               call this%set_operand(0,oprnd,jerr)
               if(jerr.ne.DSVP_SUCCESS) jerr=-6
              else
               jerr=-5
              endif
              oprnd=>NULL() !<oprnd> pointer was saved in the tensor instruction and will later be deallocated
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
         end subroutine construct_instr_create

         subroutine construct_instr_contract(jerr)
          !op_spec={tens_contraction_t}
          integer(INTD), intent(out):: jerr
          integer(INTD):: jj
          class(tens_oprnd_t), pointer:: oprnd
          class(tens_rcrsv_t), pointer:: tensor
          class(tens_contraction_t), pointer:: tens_contr
          class(contr_ptrn_ext_t), pointer:: contr_ptrn
          class(ctrl_tens_contr_t), pointer:: tens_contr_ctrl

          jerr=0
          tens_contr=>NULL()
          select type(op_spec); class is(tens_contraction_t); tens_contr=>op_spec; end select
          if(associated(tens_contr)) then
           if(tens_contr%is_set()) then
            contr_ptrn=>tens_contr%get_ext_contr_ptrn(jerr)
            if(jerr.eq.TEREC_SUCCESS) then
             allocate(tens_contr_ctrl,STAT=jerr)
             if(jerr.eq.0) then
              call tens_contr_ctrl%ctrl_tens_contr_ctor(contr_ptrn,jerr,tens_contr%get_prefactor()) !contraction pattern is cloned by value
              if(jerr.eq.0) then
               call this%set_control(tens_contr_ctrl,jerr)
               if(jerr.eq.DSVP_SUCCESS) then
                call this%alloc_operands(3,jerr)
                if(jerr.eq.DSVP_SUCCESS) then
                 do jj=0,2
                  tensor=>tens_contr%get_argument(jj,jerr); if(jerr.ne.TEREC_SUCCESS) exit
                  allocate(oprnd,STAT=jerr); if(jerr.ne.0) exit
                  call oprnd%tens_oprnd_ctor(tensor,jerr); if(jerr.ne.0) exit
                  call this%set_operand(jj,oprnd,jerr); if(jerr.ne.DSVP_SUCCESS) exit
                  tensor=>NULL(); oprnd=>NULL() !<oprnd> pointer was saved in the tensor instruction and will later be deallocated
                 enddo
                else
                 jerr=-7
                endif
               else
                jerr=-6
               endif
              else
               jerr=-5
              endif
              tens_contr_ctrl=>NULL() !<tens_contr_ctrl> pointer was saved in the tensor instruction and will later be deallocated
             else
              jerr=-4
             endif
             contr_ptrn=>NULL()
            else
             jerr=-3
            endif
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
        subroutine TensInstrDecode(this,instr_packet,ierr)
!Decodes a tensor instruction from the bytecode packet.
         implicit none
         class(tens_instr_t), intent(inout):: this       !out: tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet !in: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,op_code

         errc=0
!Read the instruction op_code:
         call unpack_builtin(instr_packet,op_code,errc)
!Select TAVP microcode to execute:
         if(errc.eq.0) then
          select case(op_code)
          case(TAVP_INSTR_NOOP)
          case(TAVP_INSTR_STOP)
          case(TAVP_INSTR_CREATE,TAVP_INSTR_DESTROY)
          case(TAVP_INSTR_CONTRACT)
          case default
           errc=-3 !unknown instruction (or not implemented)
          end select
!Activate the instruction:
          if(errc.eq.0) then
           call this%activate(op_code,errc); if(errc.ne.0) errc=-2
          else
           call this%set_status(DS_INSTR_RETIRED,errc,-1)
           call tens_instr_dtor(this)
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensInstrDecode
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
             call encode_instr_create(errc); if(errc.ne.0) errc=-6
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

         subroutine encode_instr_create(jerr)
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
         end subroutine encode_instr_create

         subroutine encode_instr_contract(jerr)
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
            call tens_contr_ctrl%pack(instr_packet,jerr)
            if(jerr.ne.0) jerr=-8
           class default
            jerr=-7
           end select
           if(jerr.eq.0) then
            do jj=0,2
             oprnd=>this%get_operand(jj,jerr); if(jerr.ne.DSVP_SUCCESS) then; jerr=-6; exit; endif
             select type(oprnd)
             class is(tens_oprnd_t)
              tensor=>oprnd%get_tensor(jerr); if(jerr.ne.0) then; jerr=-5; exit; endif
              if(.not.tensor%is_set()) then; jerr=-4; exit; endif
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
!--------------------------------------------------
        subroutine TensInstrSetMicrocode(this,ierr)
!Sets up instruction dynamic bindings to the corresponding microcode.
         implicit none
         class(tens_instr_t), intent(inout):: this   !inout: defined tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,op_code

         op_code=this%get_code(errc)
         if(errc.eq.DSVP_SUCCESS) then
          if(op_code.ge.0.and.op_code.lt.TAVP_ISA_SIZE) then
           this%acquire_resource=>microcode(op_code)%acquire_resource
           this%prefetch_input=>microcode(op_code)%prefetch_input
           this%sync_prefetch=>microcode(op_code)%sync_prefetch
           this%execute=>microcode(op_code)%execute
           this%sync_execution=>microcode(op_code)%sync_execution
           this%upload_output=>microcode(op_code)%upload_output
           this%sync_upload=>microcode(op_code)%sync_upload
           this%release_resource=>microcode(op_code)%release_resource
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensInstrSetMicrocode
!------------------------------------------------------
        subroutine TensInstrActivate(this,op_code,ierr)
!Activates the instruction after it has been constructed or decoded.
         implicit none
         class(tens_instr_t), intent(inout):: this   !inout: defined tensor instruction
         integer(INTD), intent(in):: op_code         !in: instruction code
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         call this%set_code(op_code,errc)
         if(errc.eq.DSVP_SUCCESS) then
          call this%set_microcode(errc)
          if(errc.eq.0) then
           call this%set_status(DS_INSTR_NEW,errc,DSVP_SUCCESS)
           if(errc.ne.DSVP_SUCCESS) then
            call this%set_status(DS_INSTR_RETIRED,errc,-1)
            call tens_instr_dtor(this)
            errc=-3
           endif
          else
           call this%set_status(DS_INSTR_RETIRED,errc,-1)
           call tens_instr_dtor(this)
           errc=-2
          endif
         else
          call this%set_status(DS_INSTR_RETIRED,errc,-1)
          call tens_instr_dtor(this)
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensInstrActivate
!---------------------------------------
        subroutine tens_instr_dtor(this)
         implicit none
         type(tens_instr_t):: this !inout: empty or retired tensor instruction
         integer(INTD):: sts,errc

         sts=this%get_status(errc)
         if((sts.eq.DS_INSTR_EMPTY.or.sts.eq.DS_INSTR_RETIRED).and.errc.eq.DSVP_SUCCESS) then
          call this%clean(errc)
          if(errc.ne.0) call quit(errc,'#ERROR(tavp_worker:tens_instr_dtor): Tensor instruction destruction failed!')
         else
          call quit(-1,'#ERROR(tavp_worker:tens_instr_dtor): Attempt to destroy an active TAVP instruction!')
         endif
         return
        end subroutine tens_instr_dtor

       end module tavp_worker
