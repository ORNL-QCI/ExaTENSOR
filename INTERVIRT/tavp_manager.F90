!ExaTENSOR: TAVP-Manager (TAVP-MNG) implementation
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/11/08

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
        integer(INTD), private:: DEBUG=1    !debugging mode
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
 !Bytecode:
        integer(INTL), parameter, private:: MAX_BYTECODE_SIZE=32_INTL*(1024_INTL*1024_INTL) !max size of an incoming/outgoing bytecode envelope (bytes)
        integer(INTD), parameter, private:: MAX_BYTECODE_INSTR=65536                        !max number of tensor instructions in a bytecode envelope
        integer(INTD), private:: MAX_LOCATE_NEW_INSTR=4096           !max number of new tensor instructions per bytecode envelope in the locating cycle
        integer(INTD), private:: MAX_LOCATE_DEF_INSTR=1024           !max number of deferred tensor instructions per bytecode envelope in the locating cycle
        integer(INTD), private:: MAX_DECOMPOSE_INSTR=512             !max number of input tensor instructions in the decomposition phase
        integer(INTD), private:: MAX_ISSUE_INSTR=4096                !max number of tensor instructions per bytecode issued to the child nodes
!TYPES:
 !Tensor argument cache entry (TAVP-specific):
        type, extends(tens_cache_entry_t), private:: tens_entry_mng_t
         integer(INTD), private:: owner_id=-1                         !tensor owner id (non-negative TAVP-MNG id), negative means the tensor is remote with an unknown location
         contains
          procedure, private:: TensEntryMngCtor                       !ctor
          generic, public:: tens_entry_mng_ctor=>TensEntryMngCtor
          procedure, public:: get_owner_id=>TensEntryMngGetOwnerId    !returns the owner id
          procedure, public:: holds_remote=>TensEntryMngHoldsRemote   !returns TRUE if the stored tensor is remote, FALSE otherwise
          final:: tens_entry_mng_dtor
        end type tens_entry_mng_t
 !Tensor operand (encapsulated tensor data processible by a specific TAVP):
        type, extends(ds_oprnd_t), private:: tens_oprnd_t
         class(tens_rcrsv_t), pointer, private:: tensor=>NULL()          !non-owning pointer to a persistent recursive tensor
         class(tens_entry_mng_t), pointer, private:: cache_entry=>NULL() !non-owning pointer to a tensor cache entry (optional)
         integer(INTD), private:: owner_id=-1                            !non-negative tensor owner (TAVP-MNG) id (optional)
         contains
          procedure, private:: TensOprndCtor                    !ctor
          generic, public:: tens_oprnd_ctor=>TensOprndCtor
          procedure, public:: get_tensor=>TensOprndGetTensor    !returns a pointer to the tensor
          procedure, public:: get_cache_entry=>TensOprndGetCacheEntry !returns a pointer to the tensor cache entry (may be NULL)
          procedure, public:: get_owner_id=>TensOprndGetOwnerId !returns the tensor owner id
          procedure, public:: set_owner_id=>TensOprndSetOwnerId !sets the tensor owner id
          procedure, public:: is_located=>TensOprndIsLocated    !returns TRUE if the tensor operand has been located (its structure is known)
          procedure, public:: is_remote=>TensOprndIsRemote      !returns TRUE if the tensor operand is remote
          procedure, public:: is_valued=>TensOprndIsValued      !returns TRUE if the tensor operand is set to some value (neither undefined nor being updated)
          procedure, public:: acquire_rsc=>TensOprndAcquireRsc  !explicitly acquires local resources for the tensor operand
          procedure, public:: prefetch=>TensOprndPrefetch       !starts prefetching the remote tensor operand (acquires local resources!)
          procedure, public:: upload=>TensOprndUpload           !starts uploading the tensor operand to its remote location
          procedure, public:: sync=>TensOprndSync               !synchronizes the currently pending communication on the tensor operand
          procedure, public:: release=>TensOprndRelease         !destroys the present local copy of the tensor operand (releases local resources!), but the operand stays defined
          procedure, public:: destruct=>TensOprndDestruct       !performs complete destruction back to an empty state
          final:: tens_oprnd_dtor                               !dtor
        end type tens_oprnd_t
 !Tensor instruction (realization of a tensor operation for a specific TAVP):
        type, extends(ds_instr_t), public:: tens_instr_t
         integer(INTD), private:: num_out_oprnds=0                    !number of the output tensor instruction operands
         integer(INTD), private:: out_oprnds(0:MAX_TENSOR_OPERANDS-1) !tensor instruction operands which are considered output (to sync on them for completion)
         contains
          procedure, private:: TensInstrCtor                          !ctor: constructs a tensor instruction from the specification of a tensor operation
          generic, public:: tens_instr_ctor=>TensInstrCtor
          procedure, public:: encode=>TensInstrEncode                 !encoding procedure: Packs the TAVP instruction into a raw byte packet
          procedure, public:: fully_located=>TensInstrFullyLocated    !returns TRUE if the tensor instruction operands have been fully located, FALSE otherwise
          final:: tens_instr_dtor                                     !dtor
        end type tens_instr_t
 !TAVP-MNG decoder:
        type, extends(ds_decoder_t), private:: tavp_mng_decoder_t
         integer(INTD), public:: num_ports=1                        !number of ports: Port 0 <- self
         integer(INTD), private:: source_comm                       !bytecode source communicator
         integer(INTD), private:: source_rank=-1                    !bytecode source process rank (negative means ANY)
         type(pack_env_t), private:: bytecode                       !incoming bytecode
         contains
          procedure, public:: configure=>TAVPMNGDecoderConfigure    !configures TAVP-MNG decoder
          procedure, public:: start=>TAVPMNGDecoderStart            !starts and lives TAVP-MNG decoder
          procedure, public:: shutdown=>TAVPMNGDecoderShutdown      !shuts down TAVP-MNG decoder
          procedure, public:: decode=>TAVPMNGDecoderDecode          !decodes the DS bytecode into DS instructions
        end type tavp_mng_decoder_t
 !TAVP-MNG decoder configuration:
        type, extends(dsv_conf_t), private:: tavp_mng_decoder_conf_t
         integer(INTD), public:: source_comm                        !MPI communicator of the source process
         integer(INTD), public:: source_rank                        !source process rank from which the bytecode is coming
         class(ds_unit_t), pointer, public:: acceptor=>NULL()       !non-owning pointer to the acceptor DS unit
         integer(INTD), public:: acceptor_port_id                   !associated acceptor unit port id
        end type tavp_mng_decoder_conf_t
 !TAVP-MNG retirer:
        type, extends(ds_encoder_t), private:: tavp_mng_retirer_t
         integer(INTD), public:: num_ports=1                        !number of ports: Port 0 <- collector
         integer(INTD), private:: retire_comm                       !retired bytecode destination communicator
         integer(INTD), private:: retire_rank=-1                    !retired bytecode destination process rank
         type(pack_env_t), private:: bytecode                       !outgoing bytecode
         contains
          procedure, public:: configure=>TAVPMNGRetirerConfigure    !configures TAVP-MNG retirer
          procedure, public:: start=>TAVPMNGRetirerStart            !starts and lives TAVP-MNG retirer
          procedure, public:: shutdown=>TAVPMNGRetirerShutdown      !shuts down TAVP-MNG retirer
          procedure, public:: encode=>TAVPMNGRetirerEncode          !encodes a DS instruction into the DS bytecode
        end type tavp_mng_retirer_t
 !TAVP-MNG retirer configuration:
        type, extends(dsv_conf_t), private:: tavp_mng_retirer_conf_t
         integer(INTD), public:: retire_comm                        !MPI communicator of the retired bytecode destination process
         integer(INTD), public:: retire_rank                        !destination process rank to which the retired bytecode is going
        end type tavp_mng_retirer_conf_t
 !TAVP-MNG locator:
        type, extends(ds_encoder_t), private:: tavp_mng_locator_t
         integer(INTD), public:: num_ports=2                        !number of ports: Port 0 <- udecoder; Port 1 <- ldecoder
         integer(INTD), private:: ring_comm                         !MPI communicator of the locating ring (tree level)
         integer(INTD), private:: ring_send=-1                      !MPI process rank in the locating ring to which instructions are sent
         integer(INTD), private:: ring_recv=-1                      !MPI process rank in the locating ring from which instructions are received
         type(pack_env_t), private:: bytecode                       !circulating bytecode
         type(list_bi_t), private:: located_list                    !list of tensor instructions decoded by ldecoder
         type(list_iter_t), private:: loc_list                      !iterator for <located_list>
         type(list_bi_t), private:: deferred_list                   !list of deferred tensor instructions due to unsatisfied data dependencies
         type(list_iter_t), private:: def_list                      !iterator for <deferred_list>
         type(list_bi_t), private:: missing_list                    !list of tensor instructions with missing tensor operands after the full location cycle
         type(list_iter_t), private:: mis_list                      !iterator for <missing_list>
         contains
          procedure, public:: configure=>TAVPMNGLocatorConfigure    !configures TAVP-MNG locator
          procedure, public:: start=>TAVPMNGLocatorStart            !starts and lives TAVP-MNG locator
          procedure, public:: shutdown=>TAVPMNGLocatorShutdown      !shuts down TAVP-MNG locator
          procedure, public:: encode=>TAVPMNGLocatorEncode          !encodes a DS instruction into the DS bytecode
          procedure, public:: send=>TAVPMNGLocatorSend              !sends a packet of DS instructions to the next process in the locating ring
        end type tavp_mng_locator_t
 !TAVP-MNG locator configuration:
        type, extends(dsv_conf_t), private:: tavp_mng_locator_conf_t
         integer(INTD), public:: ring_comm                          !MPI communicator of the locating ring (tree level)
         integer(INTD), public:: ring_send                          !MPI process rank in the locating ring to which instructions are sent
         integer(INTD), public:: ring_recv                          !MPI process rank in the locating ring from which instructions are received
        end type tavp_mng_locator_conf_t
 !TAVP-MNG decomposer:
        type, extends(ds_unit_t), private:: tavp_mng_decomposer_t
         integer(INTD), public:: num_ports=1                        !number of ports: Port 0 <- locator
         contains
          procedure, public:: configure=>TAVPMNGDecomposerConfigure !configures TAVP-MNG decomposer
          procedure, public:: start=>TAVPMNGDecomposerStart         !starts and lives TAVP-MNG decomposer
          procedure, public:: shutdown=>TAVPMNGDecomposerShutdown   !shuts down TAVP-MNG decomposer
          procedure, public:: decompose=>TAVPMNGDecomposerDecompose !decomposes a tensor instruction into smaller pieces
        end type tavp_mng_decomposer_t
 !TAVP-MNG decomposer configuration:
        type, extends(dsv_conf_t), private:: tavp_mng_decomposer_conf_t
        end type tavp_mng_decomposer_conf_t
 !TAVP-MNG dispatcher:
        type, extends(ds_encoder_t), private:: tavp_mng_dispatcher_t
         integer(INTD), public:: num_ports=1                        !number of ports: Port 0 <- decomposer
         integer(INTD), private:: dispatch_comm                     !MPI communicator of the processes dispatched to
         integer(INTD), private:: num_ranks=0                       !number of MPI ranks to dispatch instructions to
         integer(INTD), allocatable, private:: dispatch_rank(:)     !MPI process ranks of the processes dispatched to
         type(pack_env_t), allocatable, private:: bytecode(:)       !outgoing bytecode for each dispatched rank
         contains
          procedure, public:: configure=>TAVPMNGDispatcherConfigure !configures TAVP-MNG dispatcher
          procedure, public:: start=>TAVPMNGDispatcherStart         !starts and lives TAVP-MNG dispatcher
          procedure, public:: shutdown=>TAVPMNGDispatcherShutdown   !shuts down TAVP-MNG dispatcher
          procedure, public:: encode=>TAVPMNGDispatcherEncode       !encodes a DS instruction into the DS bytecode
          procedure, public:: dispatch=>TAVPMNGDispatcherDispatch   !dispatches a DS instruction to a specific lower-level TAVP
          procedure, public:: issue=>TAVPMNGDispatcherIssue         !issues (sends) instruction bytecode to a lower-level TAVP
        end type tavp_mng_dispatcher_t
! TAVP-MNG dispatcher configuration:
        type, extends(dsv_conf_t), private:: tavp_mng_dispatcher_conf_t
         integer(INTD), public:: dispatch_comm                      !MPI communicator of the processes dispatched to
         integer(INTD), allocatable, public:: dispatch_rank(:)      !MPI process ranks of the processes dispatched to
        end type tavp_mng_dispatcher_conf_t
 !TAVP-MNG replicator:
        type, extends(ds_encoder_t), private:: tavp_mng_replicator_t
         integer(INTD), public:: num_ports=2                        !number of ports: Port 0 <- ???; Port 1 <- rdecoder
         integer(INTD), private:: repl_comm                         !MPI communicator for tensor replication activity
         type(pack_env_t), private:: bytecode                       !tensor instruction bytecode
         contains
          procedure, public:: configure=>TAVPMNGReplicatorConfigure !configures TAVP-MNG replicator
          procedure, public:: start=>TAVPMNGReplicatorStart         !starts and lives TAVP-MNG replicator
          procedure, public:: shutdown=>TAVPMNGReplicatorShutdown   !shuts down TAVP-MNG replicator
          procedure, public:: encode=>TAVPMNGReplicatorEncode       !encodes a DS instruction into the DS bytecode
          procedure, public:: replicate=>TAVPMNGReplicatorReplicate !replicates tensors
        end type tavp_mng_replicator_t
 !TAVP-MNG replicator configuration:
        type, extends(dsv_conf_t), private:: tavp_mng_replicator_conf_t
         integer(INTD), private:: repl_comm                         !MPI communicator for tensor replication activity
        end type tavp_mng_replicator_conf_t
 !TAVP-MNG collector:
        type, extends(ds_unit_t), private:: tavp_mng_collector_t
         integer(INTD), public:: num_ports=2                        !number of ports: Port 0 <- Decomposer; Port 1 <- cdecoder
         contains
          procedure, public:: configure=>TAVPMNGCollectorConfigure  !configures TAVP-MNG collector
          procedure, public:: start=>TAVPMNGCollectorStart          !starts and lives TAVP-MNG collector
          procedure, public:: shutdown=>TAVPMNGCollectorShutdown    !shutsdown TAVP-MNG collector
          procedure, public:: collect=>TAVPMNGCollectorCollect      !collects processed tensor instructions from lower-level TAVPs
        end type tavp_mng_collector_t
 !TAVP-MNG collector configuration:
        type, extends(dsv_conf_t), private:: tavp_mng_collector_conf_t
        end type tavp_mng_collector_conf_t
 !TAVP-MNG:
        type, extends(dsvp_t), public:: tavp_mng_t
         type(tens_cache_t), private:: tens_cache                   !tensor argument cache (SHARED RESOURCE!)
         type(tavp_mng_decoder_t), private:: udecoder               !DSVU: decodes incoming tensor instructions from the higher-level manager
         type(tavp_mng_retirer_t), private:: retirer                !DSVU: retires processed tensor instructions and sends them back to the manager
         type(tavp_mng_locator_t), private:: locator                !DSVU: locates metadata for remote tensor arguments
         type(tavp_mng_decoder_t), private:: ldecoder               !DSVU: decodes tensor instructions from the locating ring
         type(tavp_mng_decomposer_t), private:: decomposer          !DSVU: decomposes tensors and tensor instructions into smaller pieces
         type(tavp_mng_dispatcher_t), private:: dispatcher          !DSVU: dispatches ready to be executed tensor instructions to the lower-level TAVPs
         type(tavp_mng_replicator_t), private:: replicator          !DSVU: replicates tensor blocks for communicaton avoiding
         type(tavp_mng_decoder_t), private:: rdecoder               !DSVU: decodes tensor create/destroy/copy instructions from the replication workflow
         type(tavp_mng_decoder_t), private:: cdecoder               !DSVU: decodes processed bytecode from the lower-level TAVPs for the collector
         type(tavp_mng_collector_t), private:: collector            !DSVU: collects processed tensor instructions from the lower-level TAVPs and updates argument cache
         contains
          procedure, public:: configure=>TAVPMNGConfigure           !configures the TAVP-MNG DSVP
        end type tavp_mng_t
 !TAVP-MNG configuration:
        type, extends(dsv_conf_t), public:: tavp_mng_conf_t
         character(:), allocatable, public:: description       !TAVP description
         integer(INTD), public:: tavp_id                       !TAVP id
         integer(INTD), public:: source_comm                   !MPI communicator of the bytecode source (higher-level manager)
         integer(INTD), public:: source_rank                   !MPI process rank of the bytecode source (higher-level manager)
         integer(INTD), public:: retire_comm                   !MPI communicator of the retired bytecode destination (higher-level manager)
         integer(INTD), public:: retire_rank                   !MPI process rank of the retired bytecode destination (higher-level manager)
         integer(INTD), public:: ring_comm                     !MPI communicator of the locating ring (tree level)
         integer(INTD), public:: ring_send_rank                !MPI process rank in the locating ring to which instructions are sent
         integer(INTD), public:: ring_recv_rank                !MPI process rank in the locating ring from which instructions are received
         integer(INTD), public:: dispatch_comm                 !MPI communicator of the processes to which the bytecode is further dispatched (lower-level)
         integer(INTD), allocatable, public:: dispatch_rank(:) !MPI ranks of the processes to which the bytecode is further dispatched (lower-level)
         integer(INTD), public:: collect_comm                  !MPI communicator of the processes from which the retired bytecode is collected (lower-level)
        end type tavp_mng_conf_t
!VISIBILITY:
 !non-member:
        public tavp_mng_reset_output
 !tens_entry_mng_t:
        private TensEntryMngCtor
        private TensEntryMngGetOwnerId
        private TensEntryMngHoldsRemote
        public tens_entry_mng_dtor
        private tens_entry_mng_alloc
 !tens_oprnd_t:
        private TensOprndCtor
        private TensOprndGetTensor
        private TensOprndGetCacheEntry
        private TensOprndGetOwnerId
        private TensOprndSetOwnerId
        private TensOprndIsLocated
        private TensOprndIsRemote
        private TensOprndIsValued
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
        private TensInstrFullyLocated
        public tens_instr_dtor
 !tavp_mng_decoder_t:
        private TAVPMNGDecoderConfigure
        private TAVPMNGDecoderStart
        private TAVPMNGDecoderShutdown
        private TAVPMNGDecoderDecode
 !tavp_mng_retirer_t:
        private TAVPMNGRetirerConfigure
        private TAVPMNGRetirerStart
        private TAVPMNGRetirerShutdown
        private TAVPMNGRetirerEncode
 !tavp_mng_locator_t:
        private TAVPMNGLocatorConfigure
        private TAVPMNGLocatorStart
        private TAVPMNGLocatorShutdown
        private TAVPMNGLocatorEncode
        private TAVPMNGLocatorSend
 !tavp_mng_decomposer_t:
        private TAVPMNGDecomposerConfigure
        private TAVPMNGDecomposerStart
        private TAVPMNGDecomposerShutdown
        private TAVPMNGDecomposerDecompose
 !tavp_mng_dispatcher_t:
        private TAVPMNGDispatcherConfigure
        private TAVPMNGDispatcherStart
        private TAVPMNGDispatcherShutdown
        private TAVPMNGDispatcherEncode
        private TAVPMNGDispatcherDispatch
        private TAVPMNGDispatcherIssue
 !tavp_mng_replicator_t:
        private TAVPMNGReplicatorConfigure
        private TAVPMNGReplicatorStart
        private TAVPMNGReplicatorShutdown
        private TAVPMNGReplicatorEncode
        private TAVPMNGReplicatorReplicate
 !tavp_mng_collector_t:
        private TAVPMNGCollectorConfigure
        private TAVPMNGCollectorStart
        private TAVPMNGCollectorShutdown
        private TAVPMNGCollectorCollect
 !tavp_mng_t:
        private TAVPMNGConfigure
!IMPLEMENTATION:
       contains
![non-member]=================================
        subroutine tavp_mng_reset_output(devo)
         implicit none
         integer(INTD), intent(in)::devo
         CONS_OUT=devo
         return
        end subroutine tavp_mng_reset_output
![tens_entry_mng_t]========================================
        subroutine TensEntryMngCtor(this,tensor,owner,ierr)
!Constructs a <tens_entry_mng_t>. Note move semantics for <tensor>!
         implicit none
         class(tens_entry_mng_t), intent(out):: this          !out: specialized tensor cache entry
         class(tens_rcrsv_t), pointer, intent(inout):: tensor !inout: pointer to an allocated tensor (ownership transfer will occur here!)
         integer(INTD), intent(in):: owner                    !in: tensor owner id (>=0)
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
!--------------------------------------------------------------
        function TensEntryMngHoldsRemote(this,ierr) result(res)
!Returns TRUE if the stored tensor is remote, FALSE otherwise.
!This information is inferred from the .owner_id field which
!has to be set in order to provide a meaningful answer.
         implicit none
         logical:: res                               !out: answer
         class(tens_entry_mng_t), intent(in):: this  !in: specialized tensor cache entry
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         res=(this%owner_id.ne.role_rank)
         if(present(ierr)) ierr=errc
         return
        end function TensEntryMngHoldsRemote
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
![tens_oprnd_t]==========================================================
        subroutine TensOprndCtor(this,tensor,ierr,tens_cache_entry,owner)
!Constructs a tensor operand. The <tensor> must be set.
!The owner id is optional.
         implicit none
         class(tens_oprnd_t), intent(inout):: this        !inout: undefined tensor operand (on entrance)
         class(tens_rcrsv_t), target, intent(in):: tensor !in: defined tensor
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD), intent(in), optional:: owner      !in: tensor owner id (no restrictions)
         class(tens_entry_mng_t), intent(inout), target, optional:: tens_cache_entry !in: tensor cache entry owning the tensor
         integer(INTD):: errc

         if(.not.this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(tensor%is_set(errc)) then
            if(errc.eq.TEREC_SUCCESS) then
             this%tensor=>tensor
             this%owner_id=-1; if(present(owner)) this%owner_id=owner
             if(present(tens_cache_entry)) then
              this%cache_entry=>tens_cache_entry
              call this%cache_entry%incr_ref_count()
             else
              this%cache_entry=>NULL()
             endif
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
!-----------------------------------------------------------------------
        function TensOprndGetCacheEntry(this,ierr) result(cache_entry_p)
!Returns a pointer to the tensor cache entry (may be NULL).
         implicit none
         class(tens_entry_mng_t), pointer:: cache_entry_p !out: pointer to the tensor cache entry
         class(tens_oprnd_t), intent(in):: this           !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD):: errc

         errc=0
         cache_entry_p=>this%cache_entry
         if(.not.associated(this%tensor)) errc=-1
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetCacheEntry
!---------------------------------------------------------
        function TensOprndGetOwnerId(this,ierr) result(id)
!Returns the tensor owner id.
         implicit none
         integer(INTD):: id                          !out: tensor owner id (negative value means self)
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; id=-1
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
!Sets the tensor owner id.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: active tensor operand
         integer(INTD), intent(in):: owner           !in: tensor owner id (no restrictions)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           this%owner_id=owner
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndSetOwnerId
!---------------------------------------------------------
        function TensOprndIsLocated(this,ierr) result(res)
!Returns TRUE if the tensor operand has been located, FALSE otherwise.
!By being located, it means its structure is known.
         implicit none
         logical:: res                               !out: result
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         class(tens_rcrsv_t), pointer:: tensor
         class(tens_body_t), pointer:: tens_body

         res=.FALSE.
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           tensor=>this%get_tensor(errc)
           if(errc.eq.0) then
            tens_body=>tensor%get_body(errc)
            if(errc.eq.TEREC_SUCCESS) then
             res=(tens_body%get_num_subtensors(errc).gt.0) !subtensors define the structure of the tensor
             if(errc.ne.TEREC_SUCCESS) then; res=.FALSE.; errc=-5; endif
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
        end function TensOprndIsLocated
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
            res=(id.ne.role_rank)
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
!--------------------------------------------------------
        function TensOprndIsValued(this,ierr) result(res)
!Returns TRUE if the tensor operand is set to some value, FALSE otherwise.
!By being set to some value, it means it is neither undefined nor being updated.
         implicit none
         logical:: res                               !out: result
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,sts
         logical:: laid,locd

         res=.FALSE.
         if(this%is_active(errc)) then
          if(this%tensor%is_set(errc,layed=laid,located=locd)) then
           if(errc.eq.TEREC_SUCCESS) then
            res=laid.and.locd.and.(this%tensor%get_state(errc).eq.TEREC_BODY_DEF)
            if(errc.ne.TEREC_SUCCESS) then; res=.FALSE.; errc=-4; endif
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
        end function TensOprndIsValued
!------------------------------------------------
        subroutine TensOprndAcquireRsc(this,ierr)
!Acquires local resources for the remote tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         !No local resources are currently needed
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

         res=.FALSE.
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(this%get_comm_stat().ne.DS_OPRND_NO_COMM) then
            call this%mark_delivered(errc); if(errc.ne.DSVP_SUCCESS) errc=-3
           endif
           if(errc.eq.0) res=.TRUE.
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
         !No local resources are currently needed
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
           this%tensor=>NULL()
           if(associated(this%cache_entry)) then
            call this%cache_entry%decr_ref_count()
            this%cache_entry=>NULL()
           endif
           this%owner_id=-1
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
         return
        end subroutine tens_oprnd_dtor
![tens_instr_t]================================================
        subroutine TensInstrCtor(this,op_code,ierr,op_spec,iid)
!Constructs a tensor instruction from a given tensor operation.
!The tensor instruction is a realization of a given tensor operation
!for a specific TAVP kind.
         implicit none
         class(tens_instr_t), intent(inout):: this        !out: tensor instruction (must be empty on entrance)
         integer(INTD), intent(in):: op_code              !in: instruction code (see top of this module)
         integer(INTD), intent(out), optional:: ierr      !out: error code
         class(*), intent(in), target, optional:: op_spec !in: formal operation specification
         integer(INTL), intent(in), optional:: iid        !in: instruction id (>=0)
         integer(INTD):: errc,ier

         if(this%is_empty(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
!Construct the instruction:
           select case(op_code)
           case(TAVP_INSTR_NOOP)
           case(TAVP_INSTR_CTRL_STOP)
            call construct_instr_stop(errc); if(errc.ne.0) errc=-8
           case(TAVP_INSTR_TENS_CREATE,TAVP_INSTR_TENS_DESTROY)
            call construct_instr_tens_create_destroy(errc); if(errc.ne.0) errc=-7
           case(TAVP_INSTR_TENS_CONTRACT)
            call construct_instr_tens_contract(errc); if(errc.ne.0) errc=-6
           case default
            errc=-5 !invalid instruction opcode (or not implemented)
           end select
!Activate the instruction:
           if(errc.eq.0) then
            if(present(iid)) then
             call this%activate(op_code,errc,iid=iid); if(errc.ne.0) errc=-4
            else
             call this%activate(op_code,errc); if(errc.ne.0) errc=-3
            endif
           endif
           if(errc.ne.0) call this%set_status(DS_INSTR_RETIRED,ier,TAVP_ERR_GEN_FAILURE)
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return

        contains

         subroutine construct_instr_stop(jerr)
          !STOP TAVP:
          !No op_spec, no control, no operands
          integer(INTD), intent(out):: jerr

          jerr=0
          return
         end subroutine construct_instr_stop

         subroutine construct_instr_tens_create_destroy(jerr)
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
              call tens_oprnd%tens_oprnd_ctor(tensor,jerr) !`tensor owner id is omitted here
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
         end subroutine construct_instr_tens_create_destroy

         subroutine construct_instr_tens_contract(jerr)
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
                  call tens_oprnd%tens_oprnd_ctor(tensor,jerr); if(jerr.ne.0) exit
                  oprnd=>tens_oprnd; call this%set_operand(jj,oprnd,jerr); if(jerr.ne.DSVP_SUCCESS) exit !ownership transfer for oprnd=tens_oprnd
                  oprnd=>NULL(); tens_oprnd=>NULL(); tensor=>NULL() !<oprnd> pointer was saved in the tensor instruction and will later be deallocated
                 enddo
                 this%num_out_oprnds=1; this%out_oprnds(0:this%num_out_oprnds-1)=(/0/)
                else
                 jerr=-7
                endif
               else
                jerr=-6
               endif
               instr_ctrl=>NULL()
              else
               jerr=-5
              endif
              tens_contr_ctrl=>NULL() !<tens_contr_ctrl> pointer was saved in the tensor instruction and will later be deallocated
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
         end subroutine construct_instr_tens_contract

        end subroutine TensInstrCtor
!---------------------------------------------------------
        subroutine TensInstrEncode(this,instr_packet,ierr)
!Encodes a tensor instruction into the bytecode packet:
! 0. Instruction id;
! 1. Instruction code;
! 2. Instruction status;
! 3. Instruction error code;
! 4. Instruction control field (optional);
! 5. Instruction operands (optional).
         implicit none
         class(tens_instr_t), intent(in):: this          !in: defined tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet !out: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,op_code,stat,err_code
         integer(INTL):: iid

!Pack the instruction attributes (id,op_code,status,error):
         if(.not.this%is_empty(errc)) then
          iid=this%get_id(errc)
          if(errc.eq.DSVP_SUCCESS) then
           call pack_builtin(instr_packet,iid,errc)
           if(errc.eq.0) then
            op_code=this%get_code(errc)
            if(errc.eq.DSVP_SUCCESS) then
             call pack_builtin(instr_packet,op_code,errc)
             if(errc.eq.0) then
              stat=this%get_status(errc,err_code)
              if(errc.eq.DSVP_SUCCESS) then
               call pack_builtin(instr_packet,stat,errc)
               if(errc.eq.0) then
                call pack_builtin(instr_packet,err_code,errc)
                if(errc.eq.0) then
!Pack the instruction body:
                 select case(op_code)
                 case(TAVP_INSTR_NOOP)
                 case(TAVP_INSTR_CTRL_STOP)
                  call encode_instr_stop(errc); if(errc.ne.0) errc=-12
                 case(TAVP_INSTR_TENS_CREATE,TAVP_INSTR_TENS_DESTROY)
                  call encode_instr_tens_create_destroy(errc); if(errc.ne.0) errc=-11
                 case(TAVP_INSTR_TENS_CONTRACT)
                  call encode_instr_tens_contract(errc); if(errc.ne.0) errc=-10
                 case default
                  errc=-9 !invalid instruction opcode (or not implemented)
                 end select
                else
                 errc=-8
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
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return

        contains

         subroutine encode_instr_stop(jerr)
          !STOP TAVP:
          !Packet format: {id|op_code|status|error}
          integer(INTD), intent(out):: jerr

          jerr=0
          return
         end subroutine encode_instr_stop

         subroutine encode_instr_tens_create_destroy(jerr)
          !CREATE/DESTROY a tensor:
          !Packet format: {id|op_code|status|error|tensor}
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
         end subroutine encode_instr_tens_create_destroy

         subroutine encode_instr_tens_contract(jerr)
          !CONTRACT two tensors: tensor0+=tensor1*tensor2*scalar:
          !Packed format: {id|op_code|status|error|ctrl_tens_contr_t|tensor0,tensor1,tensor2}
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
            do jj=0,2 !loop over the tensor instruction operands: 0:Destination, 1:LeftInput, 2:RightInput
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
         end subroutine encode_instr_tens_contract

        end subroutine TensInstrEncode
!----------------------------------------------------------------------------
        function TensInstrFullyLocated(this,ierr,input_ready) result(located)
!Returns TRUE if the tensor instruction has been fully located, FALSE otherwise.
!Being fully located means that each tensor operand has been located.
!<input_ready> is set to TRUE when each input tensor operand is currently defined,
!that is, it is neither undefined nor being updated. Note that <input_ready> is
!generally independent of <located> (an input tensor operand can be ready,
!yet not located), although it is uncommon.
         implicit none
         logical:: located                            !out: located or not
         class(tens_instr_t), intent(in):: this       !in: tensor instruction
         integer(INTD), intent(out), optional:: ierr  !out: error code
         logical, intent(out), optional:: input_ready !out: TRUE if all input tensors are ready (their values are defined)
         integer(INTD):: errc,n,i,arg_ready(0:MAX_TENSOR_OPERANDS-1)
         class(ds_oprnd_t), pointer:: oprnd
         logical:: ready

         located=.FALSE.; ready=.FALSE.
         if(.not.this%is_empty(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           n=this%get_num_operands(errc)
           if(errc.eq.DSVP_SUCCESS) then
            arg_ready(0:n-1)=0; located=.TRUE.
            do i=0,n-1 !loop over tensor operands
             oprnd=>this%get_operand(i,errc); if(errc.ne.DSVP_SUCCESS) then; located=.FALSE.; errc=-4; exit; endif
             if(oprnd%is_valued(errc)) arg_ready(i)=1; if(errc.ne.0) then; located=.FALSE.; errc=-3; exit; endif
             located=oprnd%is_located(errc); if(errc.ne.0) then; located=.FALSE.; errc=-2; exit; endif
             if(.not.located) exit
            enddo
            if(errc.eq.0) then
             ready=.TRUE.
             arg_ready(this%out_oprnds(0:this%num_out_oprnds-1))=1 !ingore output operands
             do i=0,n-1; ready=(arg_ready(i).eq.1); if(.not.ready) exit; enddo
            endif
           endif
          endif
         else
          errc=-1
         endif
         if(present(input_ready)) input_ready=ready
         if(present(ierr)) ierr=errc
         return
        end function TensInstrFullyLocated
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
          if(DEBUG.gt.0) write(CONS_OUT,'("#FATAL(TAVP-MNG:tens_instr_dtor): TAVP instruction is still active: ",i5,1x,i5)')&
          &this%get_code(),this%get_status()
          call quit(-1,'#FATAL(TAVP-MNG:tens_instr_dtor): Attempt to destroy an active TAVP instruction!')
         endif
         return
        end subroutine tens_instr_dtor
![tavp_mng_decoder_t]=====================================
        subroutine TAVPMNGDecoderConfigure(this,conf,ierr)
!Configures this DSVU.
         implicit none
         class(tavp_mng_decoder_t), intent(inout):: this !out: configured DSVU (must not be configured on entrance)
         class(dsv_conf_t), intent(in):: conf            !in: specific DSVU configuration
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         errc=0
         select type(conf)
         type is(tavp_mng_decoder_conf_t)
          if(associated(conf%acceptor)) then
           this%source_comm=conf%source_comm
           this%source_rank=conf%source_rank
           call this%set_acceptor(conf%acceptor,conf%acceptor_port_id,errc); if(errc.ne.DSVP_SUCCESS) errc=-3
          else
           errc=-2
          endif
         class default
          errc=-1
         end select
         if(DEBUG.gt.1) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decoder configured: ",i11,1x,i6)') impir,this%source_comm,this%source_rank
          flush(CONS_OUT)
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDecoderConfigure
!------------------------------------------------
        subroutine TAVPMNGDecoderStart(this,ierr)
!Starts and lives this DSVU, calls .shutdown() at the end.
         implicit none
         class(tavp_mng_decoder_t), intent(inout):: this !inout: TAVP-MNG decoder DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier,num_packets,i,j,opcode,sts,thid,port_id
         logical:: active,stopping,new
         class(dsvp_t), pointer:: dsvp
         class(*), pointer:: uptr
         class(ds_unit_t), pointer:: acceptor
         type(list_iter_t):: iqueue
         type(obj_pack_t):: instr_packet
         type(comm_handle_t):: comm_hl
         type(tens_instr_t), pointer:: tens_instr
         type(tens_instr_t):: tens_instr_empty
         type(list_bi_t):: ctrlq
         type(list_iter_t):: lit

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decoder started as DSVU # ",i2," (thread ",i2,"): Listening to ",i11,1x,i6)')&
          &impir,this%get_id(),thid,this%source_comm,this%source_rank
          flush(CONS_OUT)
         endif
!Reserve a bytecode buffer:
         call this%bytecode%reserve_mem(ier,MAX_BYTECODE_SIZE,MAX_BYTECODE_INSTR); if(ier.ne.0.and.errc.eq.0) errc=-37
!Initialize queues and ports:
         call this%init_queue(this%num_ports,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-36
!Wait on other TAVP units:
         dsvp=>this%get_dsvp()
         call dsvp%sync_units(errc,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-35
!Work loop:
         active=((errc.eq.0).and.(this%source_comm.ne.MPI_COMM_NULL)); stopping=(.not.active)
         wloop: do while(active)
          if(.not.stopping) then
 !Receive new bytecode (if posted):
           call comm_hl%clean(ier); if(ier.ne.0.and.errc.eq.0) then; errc=-34; exit wloop; endif
           new=this%bytecode%receive(comm_hl,ier,proc_rank=this%source_rank,comm=this%source_comm)
           if(ier.ne.0.and.errc.eq.0) then; errc=-33; exit wloop; endif
           if(new) then !new bytecode is available
            call comm_hl%wait(ier); if(ier.ne.0.and.errc.eq.0) then; errc=-32; exit wloop; endif
            num_packets=this%bytecode%get_num_packets(ier); if(ier.ne.0.and.errc.eq.0) then; errc=-31; exit wloop; endif
            if(DEBUG.gt.0) then
             write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decoder unit ",i2," received ",i9," new instructions")')&
             &impir,this%get_id(),num_packets
             flush(CONS_OUT)
            endif
            if(num_packets.gt.0) then
 !Decode new bytecode:
             do i=1,num_packets
              call this%bytecode%extract_packet(i,instr_packet,ier,preclean=.TRUE.)
              if(ier.ne.0.and.errc.eq.0) then; errc=-30; exit wloop; endif
              ier=this%iqueue%append(tens_instr_empty); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-29; exit wloop; endif
              ier=this%iqueue%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-28; exit wloop; endif
              uptr=>this%iqueue%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-27; exit wloop; endif
              tens_instr=>NULL(); select type(uptr); type is(tens_instr_t); tens_instr=>uptr; end select
              if(.not.associated(tens_instr).and.errc.eq.0) then; errc=-26; exit wloop; endif
              call this%decode(tens_instr,instr_packet,ier); if(ier.ne.0.and.errc.eq.0) then; errc=-25; exit wloop; endif
              sts=tens_instr%get_status(ier,j); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-24; exit wloop; endif
              opcode=tens_instr%get_code(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-23; exit wloop; endif
              if(opcode.eq.TAVP_INSTR_CTRL_STOP) then !copy the CTRL_STOP instruction into own port 0
               if(i.ne.num_packets.and.errc.eq.0) then; errc=-22; exit wloop; endif
               ier=lit%init(ctrlq); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-21; exit wloop; endif
               ier=lit%append(tens_instr); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-20; exit wloop; endif
               ier=this%load_port(0,lit); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-19; exit wloop; endif
               ier=lit%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-18; exit wloop; endif
              endif
             enddo
 !Pass decoded instructions to the acceptor unit:
             ier=this%get_acceptor(acceptor,port_id); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-17; exit wloop; endif
             if(associated(acceptor)) then
              ier=acceptor%load_port(port_id,this%iqueue); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-16; exit wloop; endif
              if(this%iqueue%get_status().ne.GFC_IT_EMPTY.and.errc.eq.0) then; errc=-15; exit wloop; endif !trap
             else
              if(errc.eq.0) then; errc=-14; exit wloop; endif
             endif
            else !empty bytecode
             if(errc.eq.0) then; errc=-13; exit wloop; endif
            endif
            call this%bytecode%clean(ier); if(ier.ne.0.and.errc.eq.0) then; errc=-12; exit wloop; endif
           endif
          else !stopping
           ier=this%iqueue%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-11; exit wloop; endif
           if(this%iqueue%get_status().eq.GFC_IT_EMPTY) active=.FALSE.
           if(.not.this%port_empty(0,ier).and.errc.eq.0) then; errc=-10; exit wloop; endif !trap
          endif
 !Check own port for control instructions:
          ier=this%flush_port(0); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-9; exit wloop; endif
          ier=this%iqueue%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-8; exit wloop; endif
          ier=this%iqueue%get_status()
          do while(ier.eq.GFC_IT_ACTIVE)
           uptr=>NULL(); uptr=>this%iqueue%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-7; exit wloop; endif
           tens_instr=>NULL(); select type(uptr); type is(tens_instr_t); tens_instr=>uptr; end select
           if(.not.associated(tens_instr).and.errc.eq.0) then; errc=-6; exit wloop; endif
           opcode=tens_instr%get_code(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-5; exit wloop; endif
           if(opcode.eq.TAVP_INSTR_CTRL_STOP) then
            stopping=.TRUE.
            call tens_instr%set_status(DS_INSTR_RETIRED,ier,DSVP_SUCCESS)
            if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-4; exit wloop; endif
           else
            if(errc.eq.0) then; errc=-3; exit wloop; endif
           endif
           ier=this%iqueue%next(); ier=this%iqueue%get_status()
          enddo
 !Clear the instruction queue:
          if(ier.eq.GFC_IT_DONE) then
           ier=this%iqueue%delete_all(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-2; exit wloop; endif
          endif
         enddo wloop
!Record the error:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(TAVP_MNG)[",i6,"]: Decoder error ",i11," by thread ",i2)')&
         &impir,errc,thid
!Shutdown:
         call this%shutdown(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDecoderStart
!---------------------------------------------------
        subroutine TAVPMNGDecoderShutdown(this,ierr)
!Stops DSVU (returns back a clean configured state).
         implicit none
         class(tavp_mng_decoder_t), intent(inout):: this !inout: TAVP-MNG decoder DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier,thid

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decoder stopped as DSVU # ",i2," (thread ",i2,")")') impir,this%get_id(),thid
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decoder DSVU # ",i2," (thread ",i2,"): Time stats (sec): ")',ADVANCE='NO')&
          &impir,this%get_id(),thid; call this%print_timing(CONS_OUT)
          !write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decoder DSVU # ",i2,": Port empty = ",l1)')&
          !&impir,this%get_id(),this%port_empty(0) !debug
          flush(CONS_OUT)
         endif
         call this%release_queue(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-2
         call this%bytecode%destroy(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDecoderShutdown
!-----------------------------------------------------------------------
        subroutine TAVPMNGDecoderDecode(this,ds_instr,instr_packet,ierr)
!Decodes a tensor instruction from a plain bytecode.
!The bytecode format is specified in tens_instr_t.encode().
         implicit none
         class(tavp_mng_decoder_t), intent(inout):: this     !inout: TAVP-MNG decoder
         class(ds_instr_t), intent(inout), target:: ds_instr !out: decoded tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet     !in: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr         !out: error code
         integer(INTD):: errc,op_code,stat,err_code
         class(dsvp_t), pointer:: dsvp
         class(tens_cache_t), pointer:: arg_cache
         integer(INTL):: iid
         real(8):: tm

         tm=thread_wtime()
         if(ds_instr%is_empty(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
!Retrieve the TAVP argument cache:
           arg_cache=>NULL(); dsvp=>this%get_dsvp() !host DSVP
           select type(dsvp); class is(tavp_mng_t); arg_cache=>dsvp%tens_cache; end select
           if(associated(arg_cache)) then
!Extract the instruction attributes (id,opcode,status,error):
            call unpack_builtin(instr_packet,iid,errc)
            if(errc.eq.0) then
             call unpack_builtin(instr_packet,op_code,errc)
             if(errc.eq.0) then
              call unpack_builtin(instr_packet,stat,errc)
              if(errc.eq.0) then
               call unpack_builtin(instr_packet,err_code,errc)
!Extract the instruction body:
               if(errc.eq.0) then
                select case(op_code)
                case(TAVP_INSTR_NOOP)
                case(TAVP_INSTR_CTRL_STOP)
                case(TAVP_INSTR_TENS_CREATE,TAVP_INSTR_TENS_DESTROY)
                 call decode_instr_tens_create_destroy(errc); if(errc.ne.0) errc=-11
                case(TAVP_INSTR_TENS_CONTRACT)
                 call decode_instr_tens_contract(errc); if(errc.ne.0) errc=-10
                case default
                 call ds_instr%set_status(DS_INSTR_RETIRED,errc,TAVP_ERR_GEN_FAILURE)
                 errc=-9 !unknown instruction opcode (or not implemented)
                end select
!Activate the instruction:
                if(errc.eq.0) then
                 call ds_instr%activate(op_code,errc,stat,err_code,iid)
                 if(errc.ne.DSVP_SUCCESS) then
                  call ds_instr%set_status(DS_INSTR_RETIRED,errc,TAVP_ERR_GEN_FAILURE)
                  errc=-8
                 endif
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
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         tm=thread_wtime(tm); call this%update_timing(tm)
         return

         contains

          subroutine decode_instr_operands(jerr)
           !Decodes tensor instruction operands
           integer(INTD), intent(out):: jerr
           class(tens_rcrsv_t), pointer:: tensor
           class(tens_rcrsv_t), pointer:: tensor_tmp
           class(tens_oprnd_t), pointer:: tens_oprnd
           class(ds_oprnd_t), pointer:: oprnd
           class(tens_cache_entry_t), pointer:: tens_entry
           class(tens_entry_mng_t), pointer:: tens_mng_entry
           integer(INTD):: jj,jn
           logical:: stored

           jn=ds_instr%get_num_operands(jerr)
           if(jerr.eq.DSVP_SUCCESS.and.jn.ge.0) then
            tensor_tmp=>NULL()
!$OMP CRITICAL (TAVP_MNG_CACHE)
            do jj=0,jn-1 !loop over tensor operands
             allocate(tensor_tmp,STAT=jerr) !tensor will either be owned by a tensor cache entry or deallocated here
             if(jerr.ne.0) then
              tensor_tmp=>NULL()
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE); jerr=-1; exit
             endif
             call tensor_tmp%tens_rcrsv_ctor(instr_packet,jerr) !unpack tensor information into a temporary tensor
             if(jerr.ne.TEREC_SUCCESS) then
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_BTC_BAD); jerr=-2; exit
             endif
             tens_entry=>NULL(); tens_entry=>arg_cache%lookup(tensor_tmp,jerr)
             if(jerr.ne.0) then
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_CHE_FAILURE); jerr=-3; exit
             endif
             if(.not.associated(tens_entry)) then !tensor is absent in the tensor cache: Create an entry for it
              stored=arg_cache%store(tensor_tmp,tens_entry_mng_alloc,jerr,tens_entry_p=tens_entry) !tensor ownership is moved to the tensor cache entry
              if((.not.stored).or.(jerr.ne.0).or.(.not.associated(tens_entry))) then
               call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_CHE_FAILURE); jerr=-4; exit
              else
               tensor_tmp=>NULL() !tensor ownership has been transferred to the tensor cache
              endif
             else
              stored=.FALSE.
             endif
             tens_mng_entry=>NULL(); select type(tens_entry); class is(tens_entry_mng_t); tens_mng_entry=>tens_entry; end select
             if(.not.associated(tens_mng_entry)) then
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_CHE_FAILURE); jerr=-5; exit
             endif
             tensor=>tens_mng_entry%get_tensor() !use the tensor from the tensor cache
             if(.not.stored) then !the tensor was already in the tensor cache before, update it by the information from the just decoded tensor
              call tensor%update(tensor_tmp,jerr); deallocate(tensor_tmp); tensor_tmp=>NULL() !deallocate temporary tensor after importing its information
              if(jerr.ne.TEREC_SUCCESS) then
               call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE); jerr=-6; exit
              endif
             endif
             allocate(tens_oprnd,STAT=jerr) !tensor operand will be owned by the tensor instruction
             if(jerr.ne.0) then
              tens_oprnd=>NULL()
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE); jerr=-7; exit
             endif
             call tens_oprnd%tens_oprnd_ctor(tensor,jerr,tens_mng_entry) !tensor is owned by the tensor cache
             if(jerr.ne.0) then
              deallocate(tens_oprnd)
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE); jerr=-8; exit
             endif
             oprnd=>tens_oprnd; call ds_instr%set_operand(jj,oprnd,jerr) !tensor operand ownership is moved to the tensor instruction
             if(jerr.ne.DSVP_SUCCESS) then
              deallocate(tens_oprnd)
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE); jerr=-9; exit
             endif
             tens_oprnd=>NULL()
            enddo
            if(associated(tensor_tmp)) deallocate(tensor_tmp) !deallocate the temporary tensor (in case of error)
!$OMP END CRITICAL (TAVP_MNG_CACHE)
           else
            call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE); jerr=-10
           endif
           return
          end subroutine decode_instr_operands

          subroutine decode_instr_tens_create_destroy(jerr)
           !CREATE/DESTROY a tensor
           integer(INTD), intent(out):: jerr

           call ds_instr%alloc_operands(1,jerr)
           if(jerr.eq.DSVP_SUCCESS) then
            call decode_instr_operands(jerr)
            if(jerr.eq.0) then !mark output operands
             select type(ds_instr)
             class is(tens_instr_t)
              ds_instr%num_out_oprnds=1; ds_instr%out_oprnds(0:ds_instr%num_out_oprnds-1)=(/0/)
             end select
            endif
           else
            call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE); jerr=-1
           endif
           return
          end subroutine decode_instr_tens_create_destroy

          subroutine decode_instr_tens_contract(jerr)
           !CONTRACT two tensors and accumulate the result into another tensor
           integer(INTD), intent(out):: jerr
           class(ds_instr_ctrl_t), pointer:: instr_ctrl
           class(ctrl_tens_contr_t), pointer:: tens_contr_ctrl

           allocate(tens_contr_ctrl,STAT=jerr) !tensor contraction control will be owned by the tensor instruction
           if(jerr.eq.0) then
            call tens_contr_ctrl%unpack(instr_packet,jerr)
            if(jerr.eq.0) then
             instr_ctrl=>tens_contr_ctrl; call ds_instr%set_control(instr_ctrl,jerr) !tensor contraction control ownership is moved to the tensor instruction
             if(jerr.eq.DSVP_SUCCESS) then
              tens_contr_ctrl=>NULL()
              call ds_instr%alloc_operands(3,jerr)
              if(jerr.eq.DSVP_SUCCESS) then
               call decode_instr_operands(jerr)
               if(jerr.eq.0) then !mark output operands
                select type(ds_instr)
                class is(tens_instr_t)
                 ds_instr%num_out_oprnds=1; ds_instr%out_oprnds(0:ds_instr%num_out_oprnds-1)=(/0/)
                end select
               endif
              else
               call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE); jerr=-1
              endif
             else
              deallocate(tens_contr_ctrl)
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE); jerr=-2
             endif
            else
             deallocate(tens_contr_ctrl)
             call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_BTC_BAD); jerr=-3
            endif
           else
            call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE); jerr=-4
           endif
           return
          end subroutine decode_instr_tens_contract

        end subroutine TAVPMNGDecoderDecode
![tavp_mng_retirer_t]=====================================
        subroutine TAVPMNGRetirerConfigure(this,conf,ierr)
!Configures this DSVU.
         implicit none
         class(tavp_mng_retirer_t), intent(inout):: this !out: configured DSVU (must not be configured on entrance)
         class(dsv_conf_t), intent(in):: conf            !in: specific DSVU configuration
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         errc=0
         select type(conf)
         type is(tavp_mng_retirer_conf_t)
          if(conf%retire_rank.ge.0) then
           this%retire_rank=conf%retire_rank
           this%retire_comm=conf%retire_comm
          else
           errc=-2
          endif
         class default
          errc=-1
         end select
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGRetirerConfigure
!------------------------------------------------
        subroutine TAVPMNGRetirerStart(this,ierr)
!Starts and lives this DSVU, calls .shutdown() at the end.
         implicit none
         class(tavp_mng_retirer_t), intent(inout):: this !inout: TAVP-MNG retirer DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier,thid
         logical:: active,stopping
         class(dsvp_t), pointer:: dsvp

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Retirer started as DSVU # ",i2," (thread ",i2,")")') impir,this%get_id(),thid
          flush(CONS_OUT)
         endif
!Reserve a bytecode buffer:
         call this%bytecode%reserve_mem(ier,MAX_BYTECODE_SIZE,MAX_BYTECODE_INSTR); if(ier.ne.0.and.errc.eq.0) errc=-1
!Initialize queues and ports:
         call this%init_queue(this%num_ports,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
!Wait on other TAVP units:
         dsvp=>this%get_dsvp()
         call dsvp%sync_units(errc,ier); if(ier.ne.0.and.errc.eq.0) errc=-1
!Work loop:
         active=((errc.eq.0).and.(this%retire_comm.ne.MPI_COMM_NULL)); stopping=(.not.active)
         wloop: do while(active)
          !`Implement
         enddo wloop
!Record the error:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(TAVP_MNG)[",i6,"]: Retirer error ",i11," by thread ",i2)')&
         &impir,errc,thid
!Shutdown:
         call this%shutdown(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGRetirerStart
!---------------------------------------------------
        subroutine TAVPMNGRetirerShutdown(this,ierr)
!Stops DSVU (returns back a clean configured state).
         implicit none
         class(tavp_mng_retirer_t), intent(inout):: this !inout: TAVP-MNG retirer DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier,thid

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Retirer stopped as DSVU # ",i2," (thread ",i2,")")') impir,this%get_id(),thid
          flush(CONS_OUT)
         endif
         call this%release_queue(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-2
         call this%bytecode%destroy(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGRetirerShutdown
!-----------------------------------------------------------------------
        subroutine TAVPMNGRetirerEncode(this,ds_instr,instr_packet,ierr)
!Encodes a tensor instruction into a plain bytecode.
         implicit none
         class(tavp_mng_retirer_t), intent(inout):: this  !inout: retirer
         class(ds_instr_t), intent(in), target:: ds_instr !in: defined tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet  !out: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD):: errc

         if(ds_instr%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           call ds_instr%encode(instr_packet,errc); if(errc.ne.PACK_SUCCESS) errc=-3
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGRetirerEncode
![tavp_mng_locator_t]=====================================
        subroutine TAVPMNGLocatorConfigure(this,conf,ierr)
!Configures this DSVU.
         implicit none
         class(tavp_mng_locator_t), intent(inout):: this !out: configured DSVU (must not be configured on entrance)
         class(dsv_conf_t), intent(in):: conf            !in: specific DSVU configuration
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         errc=0
         select type(conf)
         type is(tavp_mng_locator_conf_t)
          if(conf%ring_send.ge.0.and.conf%ring_recv.ge.0) then
           this%ring_comm=conf%ring_comm
           this%ring_send=conf%ring_send
           this%ring_recv=conf%ring_recv
          else
           errc=-2
          endif
         class default
          errc=-1
         end select
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGLocatorConfigure
!------------------------------------------------
        subroutine TAVPMNGLocatorStart(this,ierr)
!Starts and lives this DSVU, calls .shutdown() at the end.
         implicit none
         class(tavp_mng_locator_t), intent(inout):: this !inout: TAVP-MNG locator DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier,thid,opcode,rot_num
         logical:: active,stopping,ring_exists,located,valued
         class(*), pointer:: uptr
         class(dsvp_t), pointer:: dsvp
         class(tavp_mng_t), pointer:: tavp
         class(tens_instr_t), pointer:: tens_instr
         type(list_bi_t):: control_list
         type(list_iter_t):: ctrl_list

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Locator started as DSVU # ",i2," (thread ",i2,"): Send/Recv = ",i6,"/",i6)')&
          &impir,this%get_id(),thid,this%ring_send,this%ring_recv
          flush(CONS_OUT)
         endif
!Reserve a bytecode buffer:
         call this%bytecode%reserve_mem(ier,MAX_BYTECODE_SIZE,MAX_BYTECODE_INSTR); if(ier.ne.0.and.errc.eq.0) errc=-10
!Initialize queues and ports:
         call this%init_queue(this%num_ports,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-9
!Initialize the located instruction list iterator:
         ier=this%loc_list%init(this%located_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-8
!Initialize the deferred instruction list iterator:
         ier=this%def_list%init(this%deferred_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-7
!Initialize the missing-operand instruction list iterator:
         ier=this%mis_list%init(this%missing_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-6
!Initialize the control list:
         ier=ctrl_list%init(control_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-6
!Check ring topology:
         if(this%ring_send.eq.role_rank.and.this%ring_recv.eq.role_rank) then !root TAVP-MNG (no ring)
          ring_exists=.FALSE.
         else !non-root TAVP-MNG (non-trivial ring)
          ring_exists=.TRUE.
          if(this%ring_send.eq.role_rank.or.this%ring_send.lt.0.or.&
            &this%ring_recv.eq.role_rank.or.this%ring_recv.lt.0) errc=-1 !trap
         endif
!Wait on other TAVP units:
         dsvp=>this%get_dsvp(); select type(dsvp); class is(tavp_mng_t); tavp=>dsvp; end select
         call tavp%sync_units(errc,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
!Work loop:
         active=((errc.eq.0).and.(this%ring_comm.ne.MPI_COMM_NULL)); stopping=(.not.active)
         wloop: do while(active)
 !Absorb new tensor instructions from udecoder (port 0) into the main queue:
          ier=this%flush_port(0); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-5; exit wloop; endif
 !Move the leading subset of new tensor instructions from the main queue into the locating list:
          ier=this%loc_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-3; exit wloop; endif
          ier=this%iqueue%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-3; exit wloop; endif
          ier=this%iqueue%get_status()
          if(ier.eq.GFC_IT_ACTIVE) then !queue is not empty
           ier=this%iqueue%move_list(this%loc_list,MAX_LOCATE_NEW_INSTR) !at most MAX_LOCATE_NEW_INSTR tensor instructions will be moved to the locating list
           if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-2; exit wloop; endif
          endif
 !Move a limited number of tensor instructions from the deferred list back into the locating list:
          ier=this%loc_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-3; exit wloop; endif
          ier=this%def_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-3; exit wloop; endif
          ier=this%def_list%get_status()
          if(ier.eq.GFC_IT_ACTIVE) then
           ier=this%def_list%move_list(this%loc_list,MAX_LOCATE_DEF_INSTR) !at most MAX_LOCATE_DEF_INSTR tensor instructions will be moved to the locating list from the deferred list
           if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-2; exit wloop; endif
          endif
 !Perform a full rotation of the tensor instructions being located at the given NAT level:
          if(ring_exists) then
           rot_num=0
  !Encode tensor instructions from the locating list into bytecode (also deactivate control instructions):
           ier=this%loc_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-3; exit wloop; endif
           ier=this%loc_list%get_status()
           if(ier.eq.GFC_IT_ACTIVE) then !add a dummy tensor instruction when no tensor instructions are found
            ier=GFC_SUCCESS
           else
            !`Implement
            ier=this%loc_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-3; exit wloop; endif
           endif
           eloop: do while(ier.eq.GFC_SUCCESS)
            uptr=>this%loc_list%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-3; exit wloop; endif
            tens_instr=>NULL(); select type(uptr); class is(tens_instr_t); tens_instr=>uptr; end select
            if(.not.associated(tens_instr)) then; errc=-1; exit wloop; endif
            opcode=tens_instr%get_code(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
            if(opcode.eq.TAVP_INSTR_CTRL_STOP) then
             ier=this%loc_list%move_elem(ctrl_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-3; exit wloop; endif
             cycle eloop
            endif
            ier=this%loc_list%next()
           enddo eloop
           if(ier.ne.GFC_NO_MOVE.and.errc.eq.0) then; errc=-1; exit wloop; endif
  !Send the bytecode to the next NAT node at the same tree level (ring):

  !Clean the locating list and evict the relevant remote tensors with zero reference count from the tensor cache:

  !Absorb located tensor insructions from port 1:

  !Check whether these tensor instructions belong to me (thus finish the rotation):

          endif
 !Move partially located tensor instructions from the locating list to the deferred list:
          ier=this%def_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-3; exit wloop; endif
          ier=this%loc_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-3; exit wloop; endif
          ier=this%loc_list%get_status()
          do while(ier.eq.GFC_IT_ACTIVE)
           uptr=>this%loc_list%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-3; exit wloop; endif
           tens_instr=>NULL(); select type(uptr); class is(tens_instr_t); tens_instr=>uptr; end select
           if(.not.associated(tens_instr)) then; errc=-1; exit wloop; endif
           located=tens_instr%fully_located(ier,valued); if(ier.ne.0.and.errc.eq.0) then; errc=-1; exit wloop; endif
           if(.not.(located.and.valued)) then
            ier=this%loc_list%move_elem(this%def_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-3; exit wloop; endif
           else
            ier=this%loc_list%next()
           endif
           ier=this%loc_list%get_status()
          enddo
 !Move the remaining fully located tensor instructions from the locating list to the decomposer:
          ier=this%loc_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-3; exit wloop; endif
          ier=this%loc_list%get_status()
          if(ier.eq.GFC_IT_ACTIVE) then
           ier=tavp%decomposer%load_port(0,this%loc_list); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-19; exit wloop; endif
          endif
         enddo wloop
!Release the control list:
         ier=ctrl_list%delete_all(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-1
         ier=ctrl_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-1
!Record the error:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(TAVP_MNG)[",i6,"]: Locator error ",i11," by thread ",i2)')&
         &impir,errc,thid
!Shutdown:
         call this%shutdown(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGLocatorStart
!---------------------------------------------------
        subroutine TAVPMNGLocatorShutdown(this,ierr)
!Stops DSVU (returns back a clean configured state).
         implicit none
         class(tavp_mng_locator_t), intent(inout):: this !inout: TAVP-MNG locator DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier,thid

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Locator stopped as DSVU # ",i2," (thread ",i2,")")') impir,this%get_id(),thid
          flush(CONS_OUT)
         endif
 !Deactivate the missing list:
         ier=this%mis_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%mis_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-11
           ier=this%mis_list%delete_all()
          endif
          ier=this%mis_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-10
         else
          if(errc.eq.0) errc=-9
         endif
 !Deactivate the deferred list:
         ier=this%def_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%def_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-8
           ier=this%def_list%delete_all()
          endif
          ier=this%def_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-7
         else
          if(errc.eq.0) errc=-6
         endif
 !Deactivate the locating list:
         ier=this%loc_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%loc_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-5
           ier=this%loc_list%delete_all()
          endif
          ier=this%loc_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-4
         else
          if(errc.eq.0) errc=-3
         endif
 !Deactivate the main instruction queue:
         call this%release_queue(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-2
 !Destroy the bytecode buffer:
         call this%bytecode%destroy(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
 !Record an error, if any:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGLocatorShutdown
!-----------------------------------------------------------------------
        subroutine TAVPMNGLocatorEncode(this,ds_instr,instr_packet,ierr)
!Encodes a tensor instruction into a plain bytecode.
         implicit none
         class(tavp_mng_locator_t), intent(inout):: this  !inout: TAVP-MNG locator DSVU
         class(ds_instr_t), intent(in), target:: ds_instr !in: defined tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet  !out: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD):: errc

         if(ds_instr%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           call ds_instr%encode(instr_packet,errc); if(errc.ne.0) errc=-3
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGLocatorEncode
!-----------------------------------------------
        subroutine TAVPMNGLocatorSend(this,ierr)
!Sends the current bytecode envelope to the next locating process.
         implicit none
         class(tavp_mng_locator_t), intent(inout):: this !inout: TAVP-MNG locator DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         errc=0
         !`Implement
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGLocatorSend
![tavp_mng_decomposer_t]=====================================
        subroutine TAVPMNGDecomposerConfigure(this,conf,ierr)
!Configures this DSVU.
         implicit none
         class(tavp_mng_decomposer_t), intent(inout):: this !out: configured DSVU (must not be configured on entrance)
         class(dsv_conf_t), intent(in):: conf               !in: specific DSVU configuration
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc

         errc=0
         select type(conf)
         type is(tavp_mng_decomposer_conf_t)
         class default
          errc=-1
         end select
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDecomposerConfigure
!---------------------------------------------------
        subroutine TAVPMNGDecomposerStart(this,ierr)
!Starts and lives this DSVU, calls .shutdown() at the end.
         implicit none
         class(tavp_mng_decomposer_t), intent(inout):: this !inout: TAVP-MNG decomposer DSVU
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc,ier,thid
         logical:: active,stopping
         class(dsvp_t), pointer:: dsvp

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer started as DSVU # ",i2," (thread ",i2,")")') impir,this%get_id(),thid
          flush(CONS_OUT)
         endif
!Initialize queues:
         call this%init_queue(this%num_ports,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
!Wait on other TAVP units:
         dsvp=>this%get_dsvp()
         call dsvp%sync_units(errc,ier); if(ier.ne.0.and.errc.eq.0) errc=-1
!Work loop:
         active=(errc.eq.0); stopping=(.not.active)
         wloop: do while(active)
          !`Implement
         enddo wloop
!Record the error:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(TAVP_MNG)[",i6,"]: Decomposer error ",i11," by thread ",i2)')&
         &impir,errc,thid
!Shutdown:
         call this%shutdown(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDecomposerStart
!------------------------------------------------------
        subroutine TAVPMNGDecomposerShutdown(this,ierr)
!Stops DSVU (returns back a clean configured state).
         implicit none
         class(tavp_mng_decomposer_t), intent(inout):: this !inout: TAVP-MNG decomposer DSVU
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc,ier,thid

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer stopped as DSVU # ",i2," (thread ",i2,")")')&
          &impir,this%get_id(),thid
          flush(CONS_OUT)
         endif
         call this%release_queue(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDecomposerShutdown
!-------------------------------------------------------
        subroutine TAVPMNGDecomposerDecompose(this,ierr)
!Scans DS instructions and fills in information for own tensors.
         implicit none
         class(tavp_mng_decomposer_t), intent(inout):: this !inout: TAVP-MNG decomposer DSVU
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc

         errc=0
         !`Implement
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDecomposerDecompose
![tavp_mng_dispatcher_t]=====================================
        subroutine TAVPMNGDispatcherConfigure(this,conf,ierr)
!Configures this DSVU.
         implicit none
         class(tavp_mng_dispatcher_t), intent(inout):: this !out: configured DSVU (must not be configured on entrance)
         class(dsv_conf_t), intent(in):: conf               !in: specific DSVU configuration
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc,i

         errc=0
         select type(conf)
         type is(tavp_mng_dispatcher_conf_t)
          if(allocated(conf%dispatch_rank)) then
           this%dispatch_comm=conf%dispatch_comm
           this%dispatch_rank=conf%dispatch_rank
           this%num_ranks=size(this%dispatch_rank)
           if(this%num_ranks.gt.0) then
            do i=lbound(this%dispatch_rank,1),ubound(this%dispatch_rank,1)
             if(this%dispatch_rank(i).lt.0) then; errc=-4; exit; endif
            enddo
            if(errc.ne.0) deallocate(this%dispatch_rank)
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         class default
          errc=-1
         end select
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDispatcherConfigure
!---------------------------------------------------
        subroutine TAVPMNGDispatcherStart(this,ierr)
!Starts and lives this DSVU, calls .shutdown() at the end.
         implicit none
         class(tavp_mng_dispatcher_t), intent(inout):: this !inout: TAVP-MNG dispatcher DSVU
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc,ier,thid,i
         logical:: active,stopping
         class(dsvp_t), pointer:: dsvp

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Dispatcher started as DSVU # ",i2," (thread ",i2,")")') impir,this%get_id(),thid
          flush(CONS_OUT)
         endif
!Reserve bytecode buffers:
         allocate(this%bytecode(this%num_ranks),STAT=ier); if(ier.ne.0.and.errc.eq.0) errc=-1
         if(ier.eq.0) then
          do i=1,this%num_ranks
           call this%bytecode(i)%reserve_mem(ier,MAX_BYTECODE_SIZE,MAX_BYTECODE_INSTR); if(ier.ne.0.and.errc.eq.0) errc=-1
          enddo
         endif
!Initialize queues:
         call this%init_queue(this%num_ports,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
!Wait on other TAVP units:
         dsvp=>this%get_dsvp()
         call dsvp%sync_units(errc,ier); if(ier.ne.0.and.errc.eq.0) errc=-1
!Work loop:
         active=((errc.eq.0).and.(this%dispatch_comm.ne.MPI_COMM_NULL)); stopping=(.not.active)
         wloop: do while(active)
          !`Implement
         enddo wloop
!Record the error:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(TAVP_MNG)[",i6,"]: Dispatcher error ",i11," by thread ",i2)')&
         &impir,errc,thid
!Shutdown:
         call this%shutdown(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDispatcherStart
!------------------------------------------------------
        subroutine TAVPMNGDispatcherShutdown(this,ierr)
!Stops DSVU (returns back a clean configured state).
         implicit none
         class(tavp_mng_dispatcher_t), intent(inout):: this !inout: TAVP-MNG dispacher DSVU
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc,ier,thid,i

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Dispatcher stopped as DSVU # ",i2," (thread ",i2,")")')&
          &impir,this%get_id(),thid
          flush(CONS_OUT)
         endif
         call this%release_queue(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-2
         if(allocated(this%bytecode)) then
          do i=ubound(this%bytecode,1),lbound(this%bytecode,1),-1
           call this%bytecode(i)%destroy(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
          enddo
         endif
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDispatcherShutdown
!--------------------------------------------------------------------------
        subroutine TAVPMNGDispatcherEncode(this,ds_instr,instr_packet,ierr)
!Encodes a tensor instruction into a plain bytecode.
         implicit none
         class(tavp_mng_dispatcher_t), intent(inout):: this  !inout: dispatcher
         class(ds_instr_t), target, intent(in):: ds_instr    !in: defined tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet     !out: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr         !out: error code
         integer(INTD):: errc

         if(ds_instr%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           call ds_instr%encode(instr_packet,errc); if(errc.ne.PACK_SUCCESS) errc=-3
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDispatcherEncode
!-------------------------------------------------------------------------
        subroutine TAVPMNGDispatcherDispatch(this,tens_instr,channel,ierr)
!Dispatches a tensor instruction to a specific lower-level TAVP
!and encodes it into its bytecode container.
         implicit none
         class(tavp_mng_dispatcher_t), intent(inout):: this   !inout: dispatcher
         class(tens_instr_t), target, intent(in):: tens_instr !in: defined tensor instruction
         integer(INTD), intent(in):: channel                  !in: offset in this.bytecode(:)
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTD):: errc
         type(obj_pack_t):: instr_packet

         errc=0
         if(channel.ge.lbound(this%dispatch_rank,1).and.channel.le.ubound(this%dispatch_rank,1)) then
          call this%bytecode(channel)%acquire_packet(instr_packet,errc,preclean=.TRUE.)
          if(errc.eq.PACK_SUCCESS) then
           call this%encode(tens_instr,instr_packet,errc)
           if(errc.eq.0) then
            call this%bytecode(channel)%seal_packet(errc); if(errc.ne.PACK_SUCCESS) errc=-4
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
        end subroutine TAVPMNGDispatcherDispatch
!-----------------------------------------------------------
        subroutine TAVPMNGDispatcherIssue(this,channel,ierr)
!Issues the encoded bytecode to a specific dispatch channel.
         implicit none
         class(tavp_mng_dispatcher_t), intent(inout):: this !inout: dispatcher
         integer(INTD), intent(in):: channel                !in: offset in this.bytecode(:)
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc

         errc=0
         if(channel.ge.lbound(this%dispatch_rank,1).and.channel.le.ubound(this%dispatch_rank,1)) then
          !`Implement: Send the current bytecode(channel) to the recepient
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDispatcherIssue
![tavp_mng_replicator_t]=====================================
        subroutine TAVPMNGReplicatorConfigure(this,conf,ierr)
!Configures this DSVU.
         implicit none
         class(tavp_mng_replicator_t), intent(inout):: this !out: configured DSVU (must not be configured on entrance)
         class(dsv_conf_t), intent(in):: conf               !in: specific DSVU configuration
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc

         errc=0
         select type(conf)
         type is(tavp_mng_replicator_conf_t)
          this%repl_comm=conf%repl_comm
         class default
          errc=-1
         end select
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGReplicatorConfigure
!---------------------------------------------------
        subroutine TAVPMNGReplicatorStart(this,ierr)
!Starts and lives this DSVU, calls .shutdown() at the end.
         implicit none
         class(tavp_mng_replicator_t), intent(inout):: this !inout: TAVP-MNG replicator DSVU
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc,ier,thid
         logical:: active,stopping
         class(dsvp_t), pointer:: dsvp

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Replicator started as DSVU # ",i2," (thread ",i2,")")') impir,this%get_id(),thid
          flush(CONS_OUT)
         endif
!Reserve a bytecode buffer:
         call this%bytecode%reserve_mem(ier,MAX_BYTECODE_SIZE,MAX_BYTECODE_INSTR); if(ier.ne.0.and.errc.eq.0) errc=-1
!Initialize queues:
         call this%init_queue(this%num_ports,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
!Wait on other TAVP units:
         dsvp=>this%get_dsvp()
         call dsvp%sync_units(errc,ier); if(ier.ne.0.and.errc.eq.0) errc=-1
!Work loop:
         active=((errc.eq.0).and.(this%repl_comm.ne.MPI_COMM_NULL)); stopping=(.not.active)
         wloop: do while(active)
          !`Implement
         enddo wloop
!Record the error:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(TAVP_MNG)[",i6,"]: Replicator error ",i11," by thread ",i2)')&
         &impir,errc,thid
!Shutdown:
         call this%shutdown(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGReplicatorStart
!------------------------------------------------------
        subroutine TAVPMNGReplicatorShutdown(this,ierr)
!Stops DSVU (returns back a clean configured state).
         implicit none
         class(tavp_mng_replicator_t), intent(inout):: this !inout: TAVP-MNG replicator DSVU
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc,ier,thid

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Replicator stopped as DSVU # ",i2," (thread ",i2,")")')&
          &impir,this%get_id(),thid
          flush(CONS_OUT)
         endif
         call this%release_queue(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-2
         call this%bytecode%destroy(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGReplicatorShutdown
!--------------------------------------------------------------------------
        subroutine TAVPMNGReplicatorEncode(this,ds_instr,instr_packet,ierr)
!Encodes a tensor instruction into a plain bytecode.
         implicit none
         class(tavp_mng_replicator_t), intent(inout):: this  !inout: replicator DSVU
         class(ds_instr_t), intent(in), target:: ds_instr    !in: defined tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet     !out: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr         !out: error code
         integer(INTD):: errc

         if(ds_instr%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           call ds_instr%encode(instr_packet,errc); if(errc.ne.PACK_SUCCESS) errc=-3
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGReplicatorEncode
!--------------------------------------------------------------
        subroutine TAVPMNGReplicatorReplicate(this,tensor,ierr)
!Replicates tensors.
         implicit none
         class(tavp_mng_replicator_t), intent(inout):: this !inout: replicator DSVU
         class(tens_rcrsv_t), intent(in):: tensor           !in: tensor
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc

         if(tensor%is_set(errc)) then
          !`Implement
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGReplicatorReplicate
![tavp_mng_collector_t]=====================================
        subroutine TAVPMNGCollectorConfigure(this,conf,ierr)
!Configures this DSVU.
         implicit none
         class(tavp_mng_collector_t), intent(inout):: this !out: configured DSVU (must not be configured on entrance)
         class(dsv_conf_t), intent(in):: conf              !in: specific DSVU configuration
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc

         errc=0
         select type(conf)
         type is(tavp_mng_collector_conf_t)
         class default
          errc=-1
         end select
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGCollectorConfigure
!--------------------------------------------------
        subroutine TAVPMNGCollectorStart(this,ierr)
!Starts and lives this DSVU, calls .shutdown() at the end.
         implicit none
         class(tavp_mng_collector_t), intent(inout):: this !inout: TAVP-MNG collector DSVU
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc,ier,thid
         logical:: active,stopping
         class(dsvp_t), pointer:: dsvp

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Collector started as DSVU # ",i2," (thread ",i2,")")') impir,this%get_id(),thid
          flush(CONS_OUT)
         endif
!Initialize queues:
         call this%init_queue(this%num_ports,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
!Wait on other TAVP units:
         dsvp=>this%get_dsvp()
         call dsvp%sync_units(errc,ier); if(ier.ne.0.and.errc.eq.0) errc=-1
!Work loop:
         active=(errc.eq.0); stopping=(.not.active)
         wloop: do while(active)
          !`Implement
         enddo wloop
!Record the error:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(TAVP_MNG)[",i6,"]: Collector error ",i11," by thread ",i2)')&
         &impir,errc,thid
!Shutdown:
         call this%shutdown(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGCollectorStart
!-----------------------------------------------------
        subroutine TAVPMNGCollectorShutdown(this,ierr)
!Stops DSVU (returns back a clean configured state).
         implicit none
         class(tavp_mng_collector_t), intent(inout):: this !inout: TAVP-MNG collector DSVU
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc,ier,thid

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Collector stopped as DSVU # ",i2," (thread ",i2,")")') impir,this%get_id(),thid
          flush(CONS_OUT)
         endif
         call this%release_queue(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGCollectorShutdown
!----------------------------------------------------
        subroutine TAVPMNGCollectorCollect(this,ierr)
!Collects processed tensor instructions from lower-level TAVPs.
         implicit none
         class(tavp_mng_collector_t), intent(inout):: this !inout: TAVP-MNG collector DSVU
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc

         errc=0
         !`Implement
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGCollectorCollect
![tavp_mng_t]======================================
        subroutine TAVPMNGConfigure(this,conf,ierr)
!Configures TAVP-MNG DSVP:
! * Configures static DSVU;
! * Allocates and configures dynamic DSVU;
! * Sets up global DSVU table in DSVP;
! * Sets up DSVP description and global id;
         implicit none
         class(tavp_mng_t), intent(inout), target:: this !out: configured DSVP (must not be configured on entrance)
         class(dsv_conf_t), intent(in):: conf            !in: specific DSVP configuration
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,num_units
         class(ds_unit_t), pointer:: decode_acceptor
         type(tavp_mng_decoder_conf_t):: decoder_conf
         type(tavp_mng_retirer_conf_t):: retirer_conf
         type(tavp_mng_locator_conf_t):: locator_conf
         type(tavp_mng_decomposer_conf_t):: decomposer_conf
         type(tavp_mng_dispatcher_conf_t):: dispatcher_conf
         type(tavp_mng_replicator_conf_t):: replicator_conf
         type(tavp_mng_collector_conf_t):: collector_conf

         if(.not.this%is_configured(errc)) then
          if(errc.eq.0) then
           select type(conf)
           type is(tavp_mng_conf_t)
            if(conf%tavp_id.ge.0.and.allocated(conf%description)) then
             num_units=0 !increment by one after each unit configuration
 !Configure static DSVU:
  !Up-decoder:
             decode_acceptor=>this%locator
             decoder_conf=tavp_mng_decoder_conf_t(source_comm=conf%source_comm,source_rank=conf%source_rank,&
                         &acceptor=decode_acceptor,acceptor_port_id=0)
             call this%udecoder%configure(decoder_conf,errc)
             if(errc.eq.0) then
              num_units=num_units+1
  !Retirer:
              retirer_conf=tavp_mng_retirer_conf_t(conf%retire_comm,conf%retire_rank)
              call this%retirer%configure(retirer_conf,errc)
              if(errc.eq.0) then
               num_units=num_units+1
  !Locator:
               locator_conf=tavp_mng_locator_conf_t(conf%ring_comm,conf%ring_send_rank,conf%ring_recv_rank)
               call this%locator%configure(locator_conf,errc)
               if(errc.eq.0) then
                num_units=num_units+1
  !Locating-decoder:
                decode_acceptor=>this%locator
                decoder_conf=tavp_mng_decoder_conf_t(source_comm=conf%ring_comm,source_rank=conf%ring_recv_rank,&
                            &acceptor=decode_acceptor,acceptor_port_id=1)
                call this%ldecoder%configure(decoder_conf,errc)
                if(errc.eq.0) then
                 num_units=num_units+1
  !Decomposer:
                 decomposer_conf=tavp_mng_decomposer_conf_t()
                 call this%decomposer%configure(decomposer_conf,errc)
                 if(errc.eq.0) then
                  num_units=num_units+1
  !Dispatcher:
                  dispatcher_conf=tavp_mng_dispatcher_conf_t(conf%dispatch_comm,conf%dispatch_rank)
                  call this%dispatcher%configure(dispatcher_conf,errc)
                  if(errc.eq.0) then
                   num_units=num_units+1
  !Replicator:
                   replicator_conf=tavp_mng_replicator_conf_t(role_comm)
                   call this%replicator%configure(replicator_conf,errc)
                   if(errc.eq.0) then
                    num_units=num_units+1
  !Replicating-decoder:
                    decode_acceptor=>this%replicator
                    decoder_conf=tavp_mng_decoder_conf_t(source_comm=role_comm,source_rank=MPI_ANY_SOURCE,&
                                &acceptor=decode_acceptor,acceptor_port_id=1)
                    call this%rdecoder%configure(decoder_conf,errc)
                    if(errc.eq.0) then
                     num_units=num_units+1
  !Collecting-decoder:
                     decode_acceptor=>this%collector
                     decoder_conf=tavp_mng_decoder_conf_t(source_comm=conf%collect_comm,source_rank=MPI_ANY_SOURCE,&
                                 &acceptor=decode_acceptor,acceptor_port_id=1)
                     call this%cdecoder%configure(decoder_conf,errc)
                     if(errc.eq.0) then
                      num_units=num_units+1
  !Collector:
                      collector_conf=tavp_mng_collector_conf_t()
                      call this%collector%configure(collector_conf,errc)
                      if(errc.eq.0) then
                       num_units=num_units+1
 !Set up global DSVU table (references to all DSVU):
                       call this%alloc_units(num_units,errc)
                       if(errc.eq.DSVP_SUCCESS) then
                        call this%set_unit(this%udecoder,errc)
                        if(errc.eq.DSVP_SUCCESS) then
                         call this%set_unit(this%retirer,errc)
                         if(errc.eq.DSVP_SUCCESS) then
                          call this%set_unit(this%locator,errc)
                          if(errc.eq.DSVP_SUCCESS) then
                           call this%set_unit(this%ldecoder,errc)
                           if(errc.eq.DSVP_SUCCESS) then
                            call this%set_unit(this%decomposer,errc)
                            if(errc.eq.DSVP_SUCCESS) then
                             call this%set_unit(this%dispatcher,errc)
                             if(errc.eq.DSVP_SUCCESS) then
                              call this%set_unit(this%replicator,errc)
                              if(errc.eq.DSVP_SUCCESS) then
                               call this%set_unit(this%rdecoder,errc)
                               if(errc.eq.DSVP_SUCCESS) then
                                call this%set_unit(this%cdecoder,errc)
                                if(errc.eq.DSVP_SUCCESS) then
                                 call this%set_unit(this%collector,errc)
                                 if(errc.eq.DSVP_SUCCESS) then
 !Set the DSVP id and description:
                                  call this%set_description(int(conf%tavp_id,INTL),conf%description,errc)
                                  if(errc.ne.DSVP_SUCCESS) errc=-26
                                 else
                                  errc=-25
                                 endif
                                else
                                 errc=-24
                                endif
                               else
                                errc=-23
                               endif
                              else
                               errc=-22
                              endif
                             else
                              errc=-21
                             endif
                            else
                             errc=-20
                            endif
                           else
                            errc=-19
                           endif
                          else
                           errc=-18
                          endif
                         else
                          errc=-17
                         endif
                        else
                         errc=-16
                        endif
                       else
                        errc=-15
                       endif
                      else
                       errc=-14
                      endif
                     else
                      errc=-13
                     endif
                    else
                     errc=-12
                    endif
                   else
                    errc=-11
                   endif
                  else
                   errc=-10
                  endif
                 else
                  errc=-9
                 endif
                else
                 errc=-8
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
           class default
            errc=-3
           end select
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         !write(CONS_OUT,*) '#DEBUG(TAVPMNGConfigure): Exit status ',errc !debug
         if(errc.ne.0) call this%destroy()
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGConfigure

       end module tavp_manager
