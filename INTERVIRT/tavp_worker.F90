!ExaTENSOR: TAVP Worker
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/06/22

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
          procedure, private:: TensInstrCtor                     !ctor: constructs a tensor instruction from the specification of a tensor operation
          generic, public:: tens_instr_ctor=>TensInstrCtor
          procedure, public:: decode=>TensInstrDecode            !decoding procedure: Unpacks the raw byte packet (bytecode) and constructs a TAVP instruction
          procedure, public:: encode=>TensInstrEncode            !encoding procedure: Packs the TAVP instruction into a raw byte packet (bytecode)
          final:: tens_instr_dtor                                !dtor
        end type tens_instr_t
!INTERFACES:
!DATA:
 !Tensor cache:
        type(tens_cache_t), private:: tens_cache
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
!Erases the tensor cache fully and unconditionally.
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
            select type(op_spec)
            class is(tens_rcrsv_t)
             
            class default
             errc=-1
            end select
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
         integer(INTD):: errc,op_code

         errc=0
!Read the instruction op_code:
         call unpack_builtin(instr_packet,op_code,errc)
!Select TAVP microcode to execute:
         if(errc.eq.0) then
          select case(op_code)

          case default
           errc=-1 !unknown instruction
          end select
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
         class(tens_instr_t), intent(in):: this          !in: tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet !out: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         errc=0
!Pack the instruction code (op_code):

!Pack the instruction body:

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
          if(errc.ne.0) call quit(errc,'#ERROR(tavp_worker:tens_instr_dtor): Tensor instruction destruction failed!')
         else
          call quit(-1,'#ERROR(tavp_worker:tens_instr_dtor): Attempt to destroy an active TAVP instruction!')
         endif
         return
        end subroutine tens_instr_dtor

       end module tavp_worker
