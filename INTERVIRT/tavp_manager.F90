!ExaTENSOR: TAVP-Manager (TAVP-MNG) implementation
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2018/09/13

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
!NOTES:
! # TENSOR INSTRUCTION FORMAT:
!    0. Instruction id;
!    1. Instruction code (opcode);
!    2. Instruction status;
!    3. Instruction error code;
!    4. Instruction control field (optional);
!    5. Instruction operands (optional):
!       {Owner_id,Read_count,Write_count,Tensor} for each tensor operand.
! # TAVP-MNG virtual processor processes the following classes of instructions:
!   (A) CONTROL instructions: Each instruction stalls the TAVP-MNG pipeline and
!       is executed sequentially (in-order) and individually. The STOP instruction
!       may not be followed by any other instruction. The RESUME instruction does
!       nothing and is normally used for control markup. Control instructions
!       retire locally and are not returned to the upper level of hierarchy.
!   (B) AUXILIARY instructions: Each instruction stalls the TAVP-MNG pipeline
!       and is executed sequentially in-order, but not necessarily individually.
!       That is, a contiguous batch of auxiliary instructions can be issued together.
!       Auxiliary instructions normally carry definitions and other auxiliary operations
!       for subsequent numerical tensor algebra instructions. Auxiliary instructions
!       retire locally and are not returned to the upper level of hierarchy.
!   (C) TENSOR instructions: Tensor instructions represent numerical tensor algebra
!       computations and are executed asynchronously out-of-order, only obeying
!       the data dependencies. In general, in order for a tensor instruction to get
!       issued to the lower level of hierarchy, all its input arguments must be defined,
!       that is, they must neither be undefined nor be currently updated. The only
!       exception is when multiple tensor instructions with data dependencies are issued
!       to the same TAVP such that the latter can track the data dependencies locally.
!       The output tensor operands do not have to be defined, that is, they can be either
!       undefined or currently updated, or they even may not exist. In case an output
!       operand does not exist, it will be created on-the-fly and initialized to zero.
!       In case, an output operand is undefined but is existing, it will be initialized
!       to zero before update. A tensor instruction is considered completed when all its
!       subinstructions have completed. A tensor instruction is considered completed
!       successfully when all its subinstructions have completed successfully, otherwise it
!       is considered completed with error (some subinstructions have completed with error).
! # TENSOR INSTRUCTION NUMERATION:
!   (A) The DRIVER MPI process constructs all instructions, gives them their IDs,
!       sends them to the root TAVP-MNG, and receives them back retired.
!   (B) Each TAVP-MNG accepts instructions from the upper level, decomposes them
!       into subinstructions, assigns the subinstructions their new IDs based on
!       the local numeration, and sends the subinstructions to the next level.
!   (C) Each TAVP-WRK accepts subinstructions from its corresponding TAVP-MNG
!       and executes them on all computing units available on the node.
! # TENSOR ARGUMENT CACHE:
!   (A) All tensors from tensor instructions are owned by the tensor argument cache
!       which is owned by the TAVP. Tensor argument cache is a shared resource among
!       TAVP virtual units. Consequently, all operations on the tensor cache, including
!       creation/deletion of cache entries as well as a read/update of the content of a
!       cache entry (e.g., tensor meta-data update) must be protected by locks.
!   (B) Each tensor argument cache entry has a reference count for the number
!       of active tensor instruction operands currently associated with the cache entry.
!   (C) Each tensor argument cache entry has a use count for the number of active
!       instances of returned pointers pointing to the cache entry.
!   (D) Each tensor argument cache entry has an <owner_id> field referring to a TAVP-MNG
!       which owns the tensor metadata stored in that cache entry. A negative value
!       of this field means that the owning TAVP-MNG is not (yet) known.
!   (E) Tensor argument cache entries storing local tensors (owned by the current TAVP)
!       are only deleted when the tensor is destroyed. Tensor argument cache entries
!       storing subtensors of local tensors are only deleted when the tensor is destroyed.
!       Other tensor argument cache entries are deleted once the reference/use count is zero.
!       The protection of the local tensors/subtensors is done via the persistency flag.
! # DATA/TASK DECOMPOSITION:
!   (A) Each level of the TAVP-MNG hierarchy decomposes tasks (tensor instructions) into
!       subtasks (tensor subinstructions) guided by the original data (tensor) decomposition
!       created by the execution of the TENSOR_CREATE instruction. The decomposition is
!       preceded by the intra-level metadata location cycle which obtains the information
!       on the tensor composition. At the last level of the TAVP-MNG hierarchy the just
!       decomposed tensor subinstructions undergo an additional (final) metadata location
!       cycle before being issued to the TAVP-WRK level.
! # METADATA LOCATION CYCLE:
!   (A) The tensor metadata is distributed horizontally at each level of the TAVP-MNG hierarchy.
!       During the metadata location cycle, each tensor instruction gets its tensor operand
!       metadata located. If specific tensor metadata is present in multiple locations, the
!       most complete and closest from the left instance will be used by the tensor instruction.
!ISSUES:
! # Cloning <tens_instr_t> when calling container.append(tens_instr_t): Check clonability.

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
 !Locator:
        integer(INTD), private:: MAX_LOCATE_NEW_INSTR=4096  !max number of new tensor instructions per bytecode envelope in the locating cycle
        integer(INTD), private:: MAX_LOCATE_DEF_INSTR=1024  !max number of deferred tensor instructions per bytecode envelope in the locating cycle
        integer(INTD), private:: MAX_LOCATE_SUB_INSTR=8192  !max number of subinstructions per bytecode envelope in the locating cycle (bottom TAVP-MNG NAT level only)
        integer(INTD), private:: MIN_LOCATE_INSTR=256       !min number of instructions in the locating cycle to start locating rotations
        real(8), private:: MAX_LOCATE_WAIT_TIME=0.5d-3      !max wait time (sec) until performing a locating rotation even if there is no enough instructions
 !Decomposer:
        integer(INTD), private:: MAX_DECOMPOSE_PRNT_INSTR=256   !max number of processed parent instructions in the decomposition phase
        integer(INTD), private:: MAX_DECOMPOSE_CHLD_INSTR=16384 !max number of created child instructions in the decomposition phase
        real(8), private:: MAX_DECOMPOSE_PHASE_TIME=0.5d-3      !max time (sec) before passing instructions to Dispatcher
 !Dispatcher:
        logical, private:: DISPATCH_RANDOM=.TRUE.               !if TRUE the round-robin dispatch will be replaced by the random dispatch
        integer(INTD), private:: MAX_ISSUE_INSTR=4096           !max number of tensor instructions in the bytecode issued to a child node
        integer(INTD), private:: MIN_ISSUE_INSTR=512            !min number of tensor instructions being currently processed by a child node
 !Collector:
        integer(INTD), private:: MAX_COLLECT_INSTR=8192         !max number of active tensor (sub-)instructions in the collection phase
 !Retirer:
        integer(INTD), private:: MAX_RETIRE_INSTR=4096          !max number of tensor instructions in the retirement phase
        real(8), private:: MAX_RETIRE_PHASE_TIME=0.1d-3         !max time (sec) before sending the retired instructions to the upper level
!TYPES:
 !Tensor argument cache entry (TAVP-specific):
        type, extends(tens_cache_entry_t), private:: tens_entry_mng_t
         integer(INTD), private:: owner_id=-1                         !tensor meta-data owner id (non-negative TAVP-MNG id), negative means the tensor is remote with an unknown location
         contains
          procedure, private:: TensEntryMngCtor                       !ctor
          generic, public:: tens_entry_mng_ctor=>TensEntryMngCtor
          procedure, public:: get_owner_id=>TensEntryMngGetOwnerId    !returns the tensor owner id
          procedure, public:: set_owner_id=>TensEntryMngSetOwnerId    !set the tensor owner id
          procedure, public:: holds_remote=>TensEntryMngHoldsRemote   !returns TRUE if the stored tensor is remote, FALSE otherwise
          procedure, public:: print_it=>TensEntryMngPrintIt           !prints
          final:: tens_entry_mng_dtor
        end type tens_entry_mng_t
 !Reference to the tensor argument cache entry:
        type, private:: tens_entry_mng_ref_t
         class(tens_entry_mng_t), pointer, public:: cache_entry=>NULL() !non-owning pointer to a tensor cache entry
        end type tens_entry_mng_ref_t
 !Tensor operand (encapsulated tensor data processible by a specific TAVP):
        type, extends(ds_oprnd_t), private:: tens_oprnd_t
         class(tens_entry_mng_t), pointer, private:: cache_entry=>NULL() !non-owning pointer to a tensor cache entry where the tensor is stored (optional)
         class(tens_rcrsv_t), pointer, private:: tensor=>NULL()          !non-owning pointer to a persistent recursive tensor (normally stored in the tensor cache)
         integer(INTD), private:: owner_id=-1                            !non-negative tensor meta-data owner id (TAVP-MNG id), normally a copy of the value from the tensor cache entry (optional)
         contains
          procedure, private:: TensOprndCtorTensor                       !ctor by tensor only
          procedure, private:: TensOprndCtorCache                        !ctor by cache entry only
          generic, public:: tens_oprnd_ctor=>TensOprndCtorTensor,TensOprndCtorCache
          procedure, public:: get_tensor=>TensOprndGetTensor             !returns a pointer to the tensor
          procedure, public:: get_cache_entry=>TensOprndGetCacheEntry    !returns a pointer to the tensor cache entry (may be NULL)
          procedure, public:: set_cache_entry=>TensOprndSetCacheEntry    !sets the associated tensor cache entry (may be NULL)
          procedure, public:: get_owner_id=>TensOprndGetOwnerId          !returns the tensor owner id
          procedure, public:: set_owner_id=>TensOprndSetOwnerId          !sets the tensor owner id (with or without cache update)
          procedure, public:: sync_owner_id=>TensOprndSyncOwnerId        !synchronizes the tensor owner id between public cache and private reference
          procedure, public:: reset_persistency=>TensOprndResetPersistency !resets the persistency status of the underlying tensor cache entry
          procedure, public:: get_data_descriptor=>TensOprndGetDataDescriptor !returns a pointer to the tensor data descriptor or NULL if not available yet
          procedure, public:: register_read=>TensOprndRegisterRead       !registers a new read access on the tensor operand
          procedure, public:: unregister_read=>TensOprndUnregisterRead   !unregisters a read access on the tensor operand
          procedure, public:: get_read_count=>TensOprndGetReadCount      !returns the current read access count on the tensor operand
          procedure, public:: register_write=>TensOprndRegisterWrite     !registers a new write access on the tensor operand
          procedure, public:: unregister_write=>TensOprndUnregisterWrite !unregisters a write access on the tensor operand
          procedure, public:: get_write_count=>TensOprndGetWriteCount    !returns the current read access count on the tensor operand
          procedure, public:: is_active=>TensOprndIsActive               !returns TRUE if the tensor operand is active (defined)
          procedure, public:: is_located=>TensOprndIsLocated             !returns TRUE if the tensor operand has been located (its structure is known)
          procedure, public:: get_comm_stat=>TensOprndGetCommStat        !returns the current communication status on the tensor operand data
          procedure, public:: acquire_rsc=>TensOprndAcquireRsc  !acquires local resource for the tensor operand
          procedure, public:: prefetch=>TensOprndPrefetch       !starts prefetching the remote tensor operand (may acquire local resource!)
          procedure, public:: upload=>TensOprndUpload           !starts uploading the tensor operand to its remote location
          procedure, public:: sync=>TensOprndSync               !synchronizes the currently pending communication on the tensor operand data
          procedure, public:: release_rsc=>TensOprndReleaseRsc  !releases local resource but the operand stays defined
          procedure, public:: destruct=>TensOprndDestruct       !performs complete destruction back to an empty state
          procedure, public:: lock=>TensOprndLock               !sets the lock for accessing/updating the tensor operand content
          procedure, public:: unlock=>TensOprndUnlock           !releases the access lock
          procedure, public:: print_it=>TensOprndPrintIt        !prints
          final:: tens_oprnd_dtor                               !dtor
        end type tens_oprnd_t
 !Tensor instruction (realization of a tensor operation for a specific TAVP):
        type, extends(ds_instr_t), public:: tens_instr_t
         integer(INTD), private:: num_out_oprnds=0                    !number of the output tensor instruction operands
         integer(INTD), private:: out_oprnds(0:MAX_TENSOR_OPERANDS-1) !positions of the tensor instruction operands which are considered output
         contains
          procedure, private:: TensInstrCtor                              !ctor: constructs a tensor instruction from the specification of a tensor operation
          generic, public:: tens_instr_ctor=>TensInstrCtor
          procedure, public:: encode=>TensInstrEncode                     !encoding procedure: Packs the tensor instruction into a raw byte packet
          procedure, public:: fully_located=>TensInstrFullyLocated        !returns TRUE if the tensor instruction operands have been fully located (their structure is known), FALSE otherwise
          procedure, public:: get_output_operands=>TensInstrGetOutputOperands !returns the list of the output operands by their positions
          procedure, public:: get_flops=>TensInstrGetFlops                !returns an estimate of the total number of required Flops (mul/add) and memory Words
          procedure, public:: get_operation=>TensInstrGetOperation        !returns back the encapsulated tensor operation
          procedure, public:: print_it=>TensInstrPrintIt                  !prints
          procedure, private:: extract_cache_entries=>TensInstrExtractCacheEntries !returns an array of references to tensor cache entries used by the tensor operands for subsequent eviction
          procedure, private:: remove_persistency=>TensInstrRemovePersistency !removes the persistent status from the tensor instruction operands
          final:: tens_instr_dtor                                         !dtor
        end type tens_instr_t
 !TAVP-MNG decoder:
        type, extends(ds_decoder_t), private:: tavp_mng_decoder_t
         integer(INTD), public:: num_ports=1                        !number of ports: Port 0 <- self
         integer(INTD), private:: source_comm                       !bytecode source communicator
         integer(INTD), private:: source_rank=-1                    !bytecode source process rank (negative means ANY)
         integer(INT_MPI), private:: msg_tag=TAVP_DEFAULT_TAG       !target message tag
         type(pack_env_t), private:: bytecode                       !incoming bytecode buffer
         class(tens_cache_t), pointer, private:: arg_cache=>NULL()  !non-owning pointer to the tensor argument cache
         type(list_bi_t), private:: control_list                    !list of control instructions
         type(list_iter_t), private:: ctrl_list                     !iterator for <control_list>
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
         integer(INT_MPI), public:: msg_tag                         !target message tag
         class(ds_unit_t), pointer, public:: acceptor=>NULL()       !non-owning pointer to the acceptor DS unit
         integer(INTD), public:: acceptor_port_id                   !associated acceptor unit port id
        end type tavp_mng_decoder_conf_t
 !TAVP-MNG retirer:
        type, extends(ds_encoder_t), private:: tavp_mng_retirer_t
         integer(INTD), public:: num_ports=1                        !number of ports: Port 0 <- collector
         integer(INTD), private:: retire_comm                       !retired bytecode destination communicator
         integer(INTD), private:: retire_rank=-1                    !retired bytecode destination process rank
         type(pack_env_t), private:: bytecode                       !outgoing bytecode buffer
         class(tens_cache_t), pointer, private:: arg_cache=>NULL()  !non-owning pointer to the tensor argument cache
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
         integer(INTD), public:: num_ports=3                        !number of ports: Port 0 <- udecoder; Port 1 <- ldecoder; Port 2 <- decomposer (bottom TAVP-MNG NAT level only)
         integer(INTD), private:: ring_comm                         !MPI communicator of the locating ring (tree level)
         integer(INTD), private:: ring_send=-1                      !MPI process rank in the locating ring to which instructions are sent
         integer(INTD), private:: ring_recv=-1                      !MPI process rank in the locating ring from which instructions are received
         type(pack_env_t), private:: bytecode                       !circulating bytecode buffer (location cycle)
         class(tens_cache_t), pointer, private:: arg_cache=>NULL()  !non-owning pointer to the tensor argument cache
         type(list_bi_t), private:: located_list                    !list of tensor instructions decoded by ldecoder
         type(list_iter_t), private:: loc_list                      !iterator for <located_list>
         type(list_bi_t), private:: deferred_list                   !list of deferred tensor instructions due to unsatisfied data dependencies
         type(list_iter_t), private:: def_list                      !iterator for <deferred_list>
         type(list_bi_t), private:: control_list                    !list of control instructions
         type(list_iter_t), private:: ctrl_list                     !iterator for <control_list>
         contains
          procedure, public:: configure=>TAVPMNGLocatorConfigure    !configures TAVP-MNG locator
          procedure, public:: start=>TAVPMNGLocatorStart            !starts and lives TAVP-MNG locator
          procedure, public:: shutdown=>TAVPMNGLocatorShutdown      !shuts down TAVP-MNG locator
          procedure, public:: encode=>TAVPMNGLocatorEncode          !encodes a DS instruction into the DS bytecode
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
         class(tens_cache_t), pointer, private:: arg_cache=>NULL()  !non-owning pointer to the tensor argument cache
         type(list_bi_t), private:: subinstruction_list             !list of newly created tensor subinstructions
         type(list_iter_t), private:: sub_list                      !iterator for <subinstruction_list>
         type(list_bi_t), private:: returning_list                  !list of previously decomposed tensor instructions returned by Locator (bottom TAVP-MNG only)
         type(list_iter_t), private:: ret_list                      !iterator for <returning_list>
         type(list_bi_t), private:: collecting_list                 !list of processed parent tensor instructions
         type(list_iter_t), private:: col_list                      !iterator for <collecting_list>
         type(list_bi_t), private:: auxiliary_list                  !list of auxiliary instructions
         type(list_iter_t), private:: aux_list                      !iterator for <auxiliary_list>
         type(list_bi_t), private:: control_list                    !list of control instructions
         type(list_iter_t), private:: ctrl_list                     !iterator for <control_list>
         logical, private:: tavp_is_top                             !TRUE if the unit belongs to the top (root) TAVP-MNG, FALSE otherwise
         logical, private:: tavp_is_bottom                          !TRUE if the unit belongs to a bottom TAVP-MNG, FALSE otherwise
         contains
          procedure, public:: configure=>TAVPMNGDecomposerConfigure !configures TAVP-MNG decomposer
          procedure, public:: start=>TAVPMNGDecomposerStart         !starts and lives TAVP-MNG decomposer
          procedure, public:: shutdown=>TAVPMNGDecomposerShutdown   !shuts down TAVP-MNG decomposer
          procedure, public:: decompose=>TAVPMNGDecomposerDecompose !decomposes a tensor instruction into smaller pieces (subinstructions)
        end type tavp_mng_decomposer_t
 !TAVP-MNG decomposer configuration:
        type, extends(dsv_conf_t), private:: tavp_mng_decomposer_conf_t
        end type tavp_mng_decomposer_conf_t
 !TAVP-MNG dispatcher:
        type, extends(ds_encoder_t), private:: tavp_mng_dispatcher_t
         integer(INTD), public:: num_ports=1                        !number of ports: Port 0 <- decomposer
         integer(INTD), private:: dispatch_comm                     !MPI communicator of the processes dispatched to
         integer(INTD), private:: num_ranks=0                       !number of MPI ranks to dispatch instructions to
         integer(INTD), allocatable, private:: dispatch_rank(:)     !MPI ranks of the processes dispatched to (within their communicator)
         integer(INTL), allocatable, private:: dispatch_count(:)    !current number of tensor instructions dispatched to each MPI rank
         real(8), allocatable, private:: dispatch_flops(:)          !Flop count for currently dispatched tensor instructions for each MPI rank
         type(pack_env_t), allocatable, private:: bytecode(:)       !outgoing bytecode buffer for each dispatched MPI rank
         type(comm_handle_t), allocatable, private:: comm_hl(:)     !communication handle for each dispatched MPI rank
         class(tens_cache_t), pointer, private:: arg_cache=>NULL()  !non-owning pointer to the tensor argument cache
         integer(INTD), private:: next_channel=1                    !next channel to dispatch to in the round-robin dispatch
         logical, private:: tavp_is_bottom                          !TRUE if the unit belongs to a bottom TAVP-MNG, FALSE otherwise
         contains
          procedure, public:: configure=>TAVPMNGDispatcherConfigure  !configures TAVP-MNG dispatcher
          procedure, public:: start=>TAVPMNGDispatcherStart          !starts and lives TAVP-MNG dispatcher
          procedure, public:: shutdown=>TAVPMNGDispatcherShutdown    !shuts down TAVP-MNG dispatcher
          procedure, public:: encode=>TAVPMNGDispatcherEncode        !encodes a DS instruction into the DS bytecode
          procedure, public:: map_instr=>TAVPMNGDispatcherMapInstr   !maps a DS instruction to a specific lower-level TAVP
          procedure, public:: dispatch=>TAVPMNGDispatcherDispatch    !dispatches a DS instruction to a specific lower-level TAVP bytecode buffer
          procedure, public:: issue=>TAVPMNGDispatcherIssue          !issues (sends) instructions bytecode to a lower-level TAVP (async)
          procedure, public:: sync_issue=>TAVPMNGDispatcherSyncIssue !synchronizes asynchronous instruction bytecode issue to lower-level TAVPs
          procedure, public:: update_dispatch_count=>TAVPMNGUpdateDispatchCount !updates the current instruction dispatch count
          procedure, public:: update_dispatch_flops=>TAVPMNGUpdateDispatchFlops !updates the current instruction dispatch Flop count
        end type tavp_mng_dispatcher_t
 !TAVP-MNG dispatcher configuration:
        type, extends(dsv_conf_t), private:: tavp_mng_dispatcher_conf_t
         integer(INTD), public:: dispatch_comm                      !MPI communicator of the processes dispatched to
         integer(INTD), allocatable, public:: dispatch_rank(:)      !MPI process ranks of the processes dispatched to
        end type tavp_mng_dispatcher_conf_t
 !TAVP-MNG replicator:
        type, extends(ds_encoder_t), private:: tavp_mng_replicator_t
         integer(INTD), public:: num_ports=2                        !number of ports: Port 0 <- ???; Port 1 <- rdecoder
         integer(INTD), private:: repl_comm                         !MPI communicator for tensor replication activity
         type(pack_env_t), private:: bytecode                       !tensor instruction bytecode buffer
         class(tens_cache_t), pointer, private:: arg_cache=>NULL()  !non-owning pointer to the tensor argument cache
         type(list_bi_t), private:: control_list                    !list of control instructions
         type(list_iter_t), private:: ctrl_list                     !iterator for <control_list>
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
 !TAVP-MNG collector: Entry value of the dictionary of parent instructions:
        type, private:: tavp_mng_collector_entry_t
         integer(INTL), public:: children_count=0_INTL              !number of subinstructions created by this TAVP from the parent instruction
         type(list_pos_t), public:: list_elem                       !bookmark of the parent instruction in the COLLECTOR's parent instruction list (main queue)
        end type tavp_mng_collector_entry_t
 !TAVP-MNG collector:
        type, extends(ds_unit_t), private:: tavp_mng_collector_t
         integer(INTD), public:: num_ports=2                        !number of ports: Port 0 <- Decomposer; Port 1 <- cdecoder
         class(tens_cache_t), pointer, private:: arg_cache=>NULL()  !non-owning pointer to the tensor argument cache
         type(dictionary_t), private:: parent_instr_map             !map: Parent Instruction ID --> tavp_mng_collector_entry_t
         type(dictionary_iter_t), private:: parent_instr            !iterator for <parent_instr_map>
         type(list_bi_t), private:: matching_list                   !list of parent instructions being matched against their child instructions
         type(list_iter_t), private:: mat_list                      !iterator for <matching_list>
         type(list_bi_t), private:: deferred_list                   !list of previously unmatched child instructions to be matched later
         type(list_iter_t), private:: def_list                      !iterator for <deferred_list>
         type(list_bi_t), private:: retired_list                    !list of retired parent instructions
         type(list_iter_t), private:: ret_list                      !iterator for <retired_list>
         type(list_bi_t), private:: control_list                    !list of control instructions
         type(list_iter_t), private:: ctrl_list                     !iterator for <control_list>
         contains
          procedure, public:: configure=>TAVPMNGCollectorConfigure  !configures TAVP-MNG collector
          procedure, public:: start=>TAVPMNGCollectorStart          !starts and lives TAVP-MNG collector
          procedure, public:: shutdown=>TAVPMNGCollectorShutdown    !shutsdown TAVP-MNG collector
          procedure, public:: register_instr=>TAVPMNGCollectorRegisterInstr     !registers a parent instruction in the <parent_instr_map> dictionary
          procedure, public:: unregister_instr=>TAVPMNGCollectorUnregisterInstr !unregisters a parent instruction in the <parent_instr_map> dictionary
          procedure, public:: match_subinstr=>TAVPMNGCollectorMatchSubinstr     !matches a subinstruction against its parent instruction
        end type tavp_mng_collector_t
 !TAVP-MNG collector configuration:
        type, extends(dsv_conf_t), private:: tavp_mng_collector_conf_t
        end type tavp_mng_collector_conf_t
 !TAVP-MNG:
        type, extends(dsvp_t), public:: tavp_mng_t
         type(tens_cache_t), private:: tens_cache                  !tensor argument cache (SHARED RESOURCE!)
         type(dictionary_t), private:: instr_map                   !map: Child Instruction ID --> Parent Instruction ID (SHARED RESOURCE!)
         type(tavp_mng_decoder_t), private:: udecoder              !DSVU: decodes incoming tensor instructions from the higher-level manager
         type(tavp_mng_retirer_t), private:: retirer               !DSVU: retires processed tensor instructions and sends them back to the manager
         type(tavp_mng_locator_t), private:: locator               !DSVU: locates metadata for remote tensor arguments
         type(tavp_mng_decoder_t), private:: ldecoder              !DSVU: decodes tensor instructions from the locating ring
         type(tavp_mng_decomposer_t), private:: decomposer         !DSVU: decomposes tensors and tensor instructions into smaller pieces
         type(tavp_mng_dispatcher_t), private:: dispatcher         !DSVU: dispatches ready to be executed tensor instructions to the lower-level TAVPs
#if 0
         type(tavp_mng_replicator_t), private:: replicator         !DSVU: replicates tensor blocks for communicaton avoiding
         type(tavp_mng_decoder_t), private:: rdecoder              !DSVU: decodes tensor create/destroy/copy instructions from the replication workflow
#endif
         type(tavp_mng_decoder_t), private:: cdecoder              !DSVU: decodes processed bytecode from the lower-level TAVPs for the collector
         type(tavp_mng_collector_t), private:: collector           !DSVU: collects processed tensor instructions from the lower-level TAVPs and updates argument cache
         contains
          procedure, public:: configure=>TAVPMNGConfigure              !configures the TAVP-MNG DSVP
          procedure, public:: register_instr=>TAVPMNGRegisterInstr     !registers a new child instruction (subinstruction)
          procedure, public:: unregister_instr=>TAVPMNGUnregisterInstr !unregisters a subinstruction
          procedure, public:: map_instr=>TAVPMNGMapInstr               !returns the parent instruction ID for a subinstruction
          procedure, public:: is_top=>TAVPMNGIsTop                     !returns TRUE if the TAVP-MNG is the root of the TAVP-MNG hierarchy
          procedure, public:: is_bottom=>TAVPMNGIsBottom               !returns TRUE if the TAVP-MNG is a leaf (bottom) of the TAVP-MNG hierarchy
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
        private TensEntryMngSetOwnerId
        private TensEntryMngHoldsRemote
        private TensEntryMngPrintIt
        public tens_entry_mng_dtor
        private tens_entry_mng_alloc
 !tens_oprnd_t:
        private TensOprndCtorTensor
        private TensOprndCtorCache
        private TensOprndGetTensor
        private TensOprndGetCacheEntry
        private TensOprndSetCacheEntry
        private TensOprndGetOwnerId
        private TensOprndSetOwnerId
        private TensOprndSyncOwnerId
        private TensOprndResetPersistency
        private TensOprndGetDataDescriptor
        private TensOprndRegisterRead
        private TensOprndUnregisterRead
        private TensOprndGetReadCount
        private TensOprndRegisterWrite
        private TensOprndUnregisterWrite
        private TensOprndGetWriteCount
        private TensOprndIsActive
        private TensOprndIsLocated
        private TensOprndGetCommStat
        private TensOprndAcquireRsc
        private TensOprndPrefetch
        private TensOprndUpload
        private TensOprndSync
        private TensOprndReleaseRsc
        private TensOprndDestruct
        private TensOprndLock
        private TensOprndUnlock
        private TensOprndPrintIt
        public tens_oprnd_dtor
 !tens_instr_t:
        private TensInstrCtor
        private TensInstrEncode
        private TensInstrFullyLocated
        private TensInstrGetOutputOperands
        private TensInstrGetFlops
        private TensInstrGetOperation
        private TensInstrPrintIt
        private TensInstrExtractCacheEntries
        private TensInstrRemovePersistency
        public tens_instr_dtor
        private tens_instr_locator_mark
        private tens_instr_print
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
        private TAVPMNGDispatcherMapInstr
        private TAVPMNGDispatcherDispatch
        private TAVPMNGDispatcherIssue
        private TAVPMNGDispatcherSyncIssue
        private TAVPMNGUpdateDispatchCount
        private TAVPMNGUpdateDispatchFlops
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
        private TAVPMNGCollectorRegisterInstr
        private TAVPMNGCollectorUnregisterInstr
        private TAVPMNGCollectorMatchSubinstr
 !tavp_mng_t:
        private TAVPMNGConfigure
        private TAVPMNGRegisterInstr
        private TAVPMNGUnregisterInstr
        private TAVPMNGMapInstr
        private TAVPMNGIsTop
        private TAVPMNGIsBottom
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
         integer(INTD), intent(in):: owner                    !in: tensor owner id
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTD):: errc

         errc=0
         if(associated(tensor)) then
          call this%init_lock()
          call this%lock()
          call this%set_tensor(tensor,errc)
          if(errc.eq.0) then
           tensor=>NULL() !transfer the ownership
           this%owner_id=owner
          else
           errc=-2
          endif
          call this%unlock()
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensEntryMngCtor
!------------------------------------------------------------
        function TensEntryMngGetOwnerId(this,ierr) result(id)
!Returns the tensor owner id.
         implicit none
         integer(INTD):: id                            !out: tensor owner id
         class(tens_entry_mng_t), intent(inout):: this !in: specialized tensor cache entry
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: errc

         id=-1
         call this%lock()
         if(this%is_set(errc)) then
          if(errc.eq.0) id=this%owner_id
         else
          errc=-1
         endif
         call this%unlock()
         if(present(ierr)) ierr=errc
         return
        end function TensEntryMngGetOwnerId
!------------------------------------------------------
        subroutine TensEntryMngSetOwnerId(this,id,ierr)
!Sets the tensor owner id.
         implicit none
         class(tens_entry_mng_t), intent(inout):: this !inout: specialized tensor cache entry
         integer(INTD), intent(in):: id                !in: tensor owner id
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: errc

         call this%lock()
         if(this%is_set(errc)) then
          if(errc.eq.0) this%owner_id=id
         else
          errc=-1
         endif
         call this%unlock()
         if(present(ierr)) ierr=errc
         return
        end subroutine TensEntryMngSetOwnerId
!--------------------------------------------------------------
        function TensEntryMngHoldsRemote(this,ierr) result(res)
!Returns TRUE if the stored tensor is remote, FALSE otherwise.
!This information is inferred from the .owner_id field which
!has to be set in order to provide a meaningful answer.
         implicit none
         logical:: res                                 !out: answer
         class(tens_entry_mng_t), intent(inout):: this !in: specialized tensor cache entry
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: errc,id

         errc=0
         call this%lock()
         id=this%owner_id
         call this%unlock()
         res=(id.ne.role_rank)
         if(present(ierr)) ierr=errc
         return
        end function TensEntryMngHoldsRemote
!---------------------------------------------------------------
        subroutine TensEntryMngPrintIt(this,ierr,dev_id,nspaces)
!Prints the tensor cache entry.
         implicit none
         class(tens_entry_mng_t), intent(inout):: this !in: tensor cache entry
         integer(INTD), intent(out), optional:: ierr   !out: errror code
         integer(INTD), intent(in), optional:: dev_id  !in: output device id
         integer(INTD), intent(in), optional:: nspaces !in: left alignment
         integer(INTD):: errc,devo,nsp,j,refc,usec,rwc,drwc
         class(tens_rcrsv_t), pointer:: tensor
         logical:: pers
!$OMP FLUSH
         errc=0
         devo=6; if(present(dev_id)) devo=dev_id
         nsp=0; if(present(nspaces)) nsp=nspaces
         call this%lock()
         do j=1,nsp; write(devo,'(" ")',ADVANCE='NO'); enddo
         write(devo,'("TENSOR CACHE ENTRY{")')
 !Tensor:
         tensor=>this%get_tensor(errc)
         if(errc.eq.0) call tensor%print_it(errc,devo,nsp+1)
 !Metadata owner id:
         do j=1,nsp+1; write(devo,'(" ")',ADVANCE='NO'); enddo
         write(devo,'("Metadata owner TAVP-MNG id = ",i6)') this%owner_id
 !Counters:
         pers=this%is_persistent(); refc=this%get_ref_count(); usec=this%get_use_count()
         rwc=this%get_rw_counter(); drwc=this%get_rw_counter(defer=.TRUE.)
         do j=1,nsp+1; write(devo,'(" ")',ADVANCE='NO'); enddo
         write(devo,'("Persist = ",l1,". Counters: Ref = ",i5,"; Use = ",i2,"; RW/DRW = ",i3,1x,i3)') pers,refc,usec,rwc,drwc
         do j=1,nsp; write(devo,'(" ")',ADVANCE='NO'); enddo
         write(devo,'("}")')
         call this%unlock()
         flush(devo)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensEntryMngPrintIt
!-------------------------------------------
        subroutine tens_entry_mng_dtor(this)
!Tensor cache entry dtor is only expected to be called during cache entry eviction,
!thus it is protected by the cache lock.
         implicit none
         type(tens_entry_mng_t):: this !inout: specialized tensor cache entry

         this%owner_id=-1
         call this%destroy(.TRUE.) !deallocate the tensor component (if it is set)
         return
        end subroutine tens_entry_mng_dtor
!-------------------------------------------------------------
        function tens_entry_mng_alloc(tens_entry) result(ierr)
!Non-member allocator for tens_entry_mng_t (called by the store() method of tensor cache).
         implicit none
         integer(INTD):: ierr
         class(tens_cache_entry_t), allocatable, intent(out):: tens_entry

         allocate(tens_entry_mng_t::tens_entry,STAT=ierr)
         return
        end function tens_entry_mng_alloc
![tens_oprnd_t]===============================================
        subroutine TensOprndCtorTensor(this,tensor,ierr,owner)
!Constructs a tensor operand. The <tensor> must be set (defined).
         implicit none
         class(tens_oprnd_t), intent(inout):: this        !inout: empty tensor operand (on entrance)
         class(tens_rcrsv_t), intent(in), target:: tensor !in: defined persistent tensor
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD), intent(in), optional:: owner      !in: tensor owner id (no restrictions)
         integer(INTD):: errc

         if(.not.this%is_active(errc)) then
          if(errc.eq.0) then
           if(tensor%is_set(errc)) then
            if(errc.eq.TEREC_SUCCESS) then
             this%cache_entry=>NULL()
             this%tensor=>tensor
             if(present(owner)) then
              this%owner_id=owner !metadata owner (TAVP-MNG id)
             else
              this%owner_id=-1
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
        end subroutine TensOprndCtorTensor
!----------------------------------------------------------------
        subroutine TensOprndCtorCache(this,tens_cache_entry,ierr)
!Constructs a tensor operand by importing the content of a tensor cache entry.
         implicit none
         class(tens_oprnd_t), intent(inout):: this                         !inout: empty tensor operand (on entrance)
         class(tens_entry_mng_t), intent(inout), target:: tens_cache_entry !in: tensor cache entry owning the tensor
         integer(INTD), intent(out), optional:: ierr                       !out: error code
         integer(INTD):: errc,refc

         if(.not.this%is_active(errc)) then
          if(errc.eq.0) then
           call tens_cache_entry%lock()
           call tens_cache_entry%incr_ref_count()
           if(tens_cache_entry%is_set(errc)) then
            if(errc.eq.0) then
             this%cache_entry=>tens_cache_entry
             this%tensor=>this%cache_entry%get_tensor(errc)
             if(errc.eq.0) then
              this%owner_id=this%cache_entry%get_owner_id()
              if(DEBUG.gt.1.and.errc.eq.0) then
               refc=this%cache_entry%get_ref_count()
               write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Tensor operand associated with a cache entry (",i4,"):")') impir,refc
               call this%tensor%print_it(dev_id=CONS_OUT); flush(CONS_OUT)
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
           if(errc.ne.0) call tens_cache_entry%decr_ref_count()
           call tens_cache_entry%unlock()
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndCtorCache
!------------------------------------------------------------
        function TensOprndGetTensor(this,ierr) result(tens_p)
!Returns a pointer to the tensor.
         implicit none
         class(tens_rcrsv_t), pointer:: tens_p       !out: pointer to the tensor
         class(tens_oprnd_t), intent(inout):: this   !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
!$OMP FLUSH
         errc=0; tens_p=>this%tensor
         if(.not.associated(tens_p)) errc=-1
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetTensor
!-----------------------------------------------------------------------
        function TensOprndGetCacheEntry(this,ierr) result(cache_entry_p)
!Returns a pointer to the tensor cache entry (may be NULL).
         implicit none
         class(tens_entry_mng_t), pointer:: cache_entry_p !out: pointer to the tensor cache entry (may be NULL)
         class(tens_oprnd_t), intent(in):: this           !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD):: errc
!$OMP FLUSH
         errc=0; cache_entry_p=>this%cache_entry
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetCacheEntry
!---------------------------------------------------------------
        subroutine TensOprndSetCacheEntry(this,cache_entry,ierr)
!Sets the associated tensor cache entry (may be NULL). If the cache entry
!is not NULL, it must be the one containing the tensor already set in the tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this                  !inout: active tensor operand
         class(tens_entry_mng_t), intent(in), pointer:: cache_entry !in: pointer to a tensor cache entry containing the same tensor (may be NULL)
         integer(INTD), intent(out), optional:: ierr                !out: error code
         integer(INTD):: errc
         class(tens_rcrsv_t), pointer:: tensor

         if(this%is_active(errc)) then
          if(errc.eq.0) then
           if(associated(this%cache_entry)) then
            call this%cache_entry%decr_ref_count(); this%cache_entry=>NULL()
           endif
           if(associated(cache_entry)) then
            call cache_entry%lock()
            tensor=>cache_entry%get_tensor(errc)
            if(errc.eq.0.and.associated(this%tensor,tensor)) then !cache entry corresponds to the same tensor
             call cache_entry%incr_ref_count()
             this%cache_entry=>cache_entry
             this%owner_id=this%cache_entry%get_owner_id()
            else
             errc=-3
            endif
            call cache_entry%unlock()
           else
            this%cache_entry=>NULL() !`owner_id is kept unchanged in the tensor operand
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndSetCacheEntry
!---------------------------------------------------------
        function TensOprndGetOwnerId(this,ierr) result(id)
!Returns the tensor owner id.
         implicit none
         integer(INTD):: id                          !out: tensor owner id
         class(tens_oprnd_t), intent(inout):: this   !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; id=-1
         !call this%lock() !only needed if syncing owner_id
         if(associated(this%tensor)) then
          !if(this%owner_id.lt.0.and.associated(this%cache_entry)) this%owner_id=this%cache_entry%get_owner_id(errc) !sync owner_id
          id=this%owner_id
         else
          errc=-1
         endif
         !call this%unlock() !only needed if syncing owner_id
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetOwnerId
!-------------------------------------------------------------------
        subroutine TensOprndSetOwnerId(this,owner,ierr,update_cache)
!Sets the tensor owner id, with or without tensor cache update.
         implicit none
         class(tens_oprnd_t), intent(inout):: this    !inout: active tensor operand
         integer(INTD), intent(in):: owner            !in: tensor owner id (no restrictions)
         integer(INTD), intent(out), optional:: ierr  !out: error code
         logical, intent(in), optional:: update_cache !if TRUE, the owner id in the tensor cache will be updated as well (defaults to FALSE)
         integer(INTD):: errc

         if(this%is_active(errc)) then
          if(errc.eq.0) then
           this%owner_id=owner
           if(present(update_cache)) then
            if(update_cache) then
             if(associated(this%cache_entry)) then
              call this%cache_entry%set_owner_id(owner,errc); if(errc.ne.0) errc=-4
             else
              errc=-3
             endif
            endif
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
!-------------------------------------------------
        subroutine TensOprndSyncOwnerId(this,ierr)
!Synchronizes the tensor owner id between public cache and private reference,
!that is, overwrites the private reference with the public cache value, if any.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%is_active(errc)) then
          if(errc.eq.0) then
           if(associated(this%cache_entry)) then
            this%owner_id=this%cache_entry%get_owner_id(errc); if(errc.ne.0) errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndSyncOwnerId
!------------------------------------------------------------------
        subroutine TensOprndResetPersistency(this,persistency,ierr)
!Resets the persistency status of the underlying tensor cache entry.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         logical, intent(in):: persistency           !in: new persistency status
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         call this%lock()
         if(associated(this%cache_entry)) then
          call this%cache_entry%set_persistency(persistency)
         else
          errc=-1
         endif
         call this%unlock()
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndResetPersistency
!-------------------------------------------------------------------
        function TensOprndGetDataDescriptor(this,ierr) result(descr)
!Returns a pointer to the tensor data descriptor or NULL if not available yet (not an error).
         implicit none
         class(DataDescr_t), pointer:: descr         !out: pointer to the tensor data descriptor
         class(tens_oprnd_t), intent(inout):: this   !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         class(tens_rcrsv_t), pointer:: tensor
         class(tens_layout_t), pointer:: layout
         logical:: locd

         descr=>NULL()
         call this%lock()
         tensor=>this%get_tensor(errc)
         if(errc.eq.0) then
          if(tensor%is_set(errc,located=locd)) then
           if(errc.eq.TEREC_SUCCESS) then
            if(locd) then
             layout=>tensor%get_layout(errc)
             if(errc.eq.TEREC_SUCCESS) then
              descr=>layout%get_data_descr(errc); if(errc.ne.TEREC_SUCCESS) errc=-5
             else
              errc=-4
             endif
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
         call this%unlock()
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetDataDescriptor
!--------------------------------------------------
        subroutine TensOprndRegisterRead(this,ierr)
!Registers a read access on the tensor operand. In case the tensor
!operand is not associated with a tensor cache entry, does nothing.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(associated(this%cache_entry)) call this%cache_entry%incr_read_count()
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndRegisterRead
!----------------------------------------------------
        subroutine TensOprndUnregisterRead(this,ierr)
!Unregisters a read access on the tensor operand. In case the tensor
!operand is not associated with a tensor cache entry, does nothing.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(associated(this%cache_entry)) call this%cache_entry%decr_read_count()
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndUnregisterRead
!--------------------------------------------------------------
        function TensOprndGetReadCount(this,ierr) result(count)
!Returns the current read access count on the tensor operand. In case the tensor
!operand is not associated with a tensor cache entry, returns zero.
         implicit none
         integer(INTD):: count                       !out: read count
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; count=0
         if(associated(this%cache_entry)) count=this%cache_entry%get_read_count()
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetReadCount
!---------------------------------------------------
        subroutine TensOprndRegisterWrite(this,ierr)
!Registers a write access on the tensor operand. In case the tensor
!operand is not associated with a tensor cache entry, does nothing.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(associated(this%cache_entry)) call this%cache_entry%incr_write_count()
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndRegisterWrite
!-----------------------------------------------------
        subroutine TensOprndUnregisterWrite(this,ierr)
!Unregisters a write access on the tensor operand. In case the tensor
!operand is not associated with a tensor cache entry, does nothing.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(associated(this%cache_entry)) call this%cache_entry%decr_write_count()
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndUnregisterWrite
!---------------------------------------------------------------
        function TensOprndGetWriteCount(this,ierr) result(count)
!Returns the current write access count on the tensor operand. In case the tensor
!operand is not associated with a tensor cache entry, returns zero.
         implicit none
         integer(INTD):: count                       !out: write count
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; count=0
         if(associated(this%cache_entry)) count=this%cache_entry%get_write_count()
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetWriteCount
!--------------------------------------------------------
        function TensOprndIsActive(this,ierr) result(ans)
!Returns TRUE if the tensor operand is active (defined).
         implicit none
         logical:: ans                               !out: answer
         class(tens_oprnd_t), intent(inout):: this   !in: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
!$OMP FLUSH
         ans=associated(this%tensor)
         if(present(ierr)) ierr=0
         return
        end function TensOprndIsActive
!-----------------------------------------------------------------------
        function TensOprndIsLocated(this,ierr,remote,valued) result(res)
!Returns TRUE if the tensor operand has been located, FALSE otherwise,
!plus some additional attributes. By being located, it means that the
!tensor structure (its subtensor composition) is known.
         implicit none
         logical:: res                               !out: result
         class(tens_oprnd_t), intent(inout):: this   !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(out), optional:: remote     !out: TRUE if the operand is remote
         logical, intent(out), optional:: valued     !out: TRUE if the operand is valued (neither undefined nor being updated)
         integer(INTD):: errc
         logical:: laid,locd

         res=.FALSE.
         if(this%is_active(errc)) then
          if(errc.eq.0) then
           call this%lock()
           if(associated(this%tensor)) then
            res=(this%tensor%get_num_subtensors(errc).gt.0) !subtensors define the structure of the tensor (which is being located)
            if(errc.eq.TEREC_SUCCESS) then
             if(present(valued)) then
              if(this%tensor%is_set(errc,layed=laid,located=locd)) then
               if(errc.eq.TEREC_SUCCESS) then
                valued=(this%get_write_count().eq.0) !no one is currently writing to this tensor => DEFINED
               else
                errc=-7
               endif
              else
               errc=-6
              endif
             endif
             if(present(remote).and.errc.eq.0) then
              remote=(this%get_owner_id(errc).ne.role_rank)
              if(errc.ne.0) errc=-5
             endif
            else
             res=.FALSE.; errc=-4
            endif
           else
            errc=-3
           endif
           call this%unlock()
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(errc.ne.0) then
          if(present(remote)) remote=.TRUE.
          if(present(valued)) valued=.FALSE.
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensOprndIsLocated
!----------------------------------------------------------------
        function TensOprndGetCommStat(this,ierr,req) result(stat)
!Returns the current communication status on the tensor operand data.
         implicit none
         integer(INTD):: stat                        !out: communication status: {DS_OPRND_NO_COMM,DS_OPRND_FETCHING,DS_OPRND_UPLOADING}
         class(tens_oprnd_t), intent(inout):: this   !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD), intent(out), optional:: req  !out: communication request handle
         integer(INTD):: errc

         stat=DS_OPRND_NO_COMM !TAVP-MNG does not perform tensor body data communications
         if(this%is_active(errc)) then
          if(errc.ne.0) errc=-2
         else
          errc=-1
         endif
         if(present(req)) req=MPI_REQUEST_NULL
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetCommStat
!---------------------------------------------------------
        subroutine TensOprndAcquireRsc(this,ierr,init_rsc)
!Acquires local resources for the remote tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: init_rsc    !in: if TRUE, the acquired resource will be initialized (not used)
         integer(INTD):: errc

         if(this%is_active(errc)) then
          if(errc.ne.0) errc=-2 !No local resources are currently needed
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndAcquireRsc
!----------------------------------------------
        subroutine TensOprndPrefetch(this,ierr)
!Starts prefetching the (remote) tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%is_active(errc)) then
          if(errc.ne.0) errc=-2 !No tensor body data prefetch is needed
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
         class(tens_oprnd_t), intent(inout):: this   !inout: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%is_active(errc)) then
          if(errc.ne.0) errc=-2 !No tensor body data upload is needed
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
         class(tens_oprnd_t), intent(inout):: this   !inout: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: TRUE activates WAIT instead of TEST synchronization (default)
         integer(INTD):: errc

         res=.FALSE.
         if(this%is_active(errc)) then
          if(errc.eq.0) then
           res=.TRUE.
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensOprndSync
!------------------------------------------------
        subroutine TensOprndReleaseRsc(this,ierr)
!Releases local resources acquired for the remote tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%is_active(errc)) then
          if(errc.ne.0) errc=-2 !No local resources are currently needed
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndReleaseRsc
!----------------------------------------------
        subroutine TensOprndDestruct(this,ierr)
!Destructs the tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
!$OMP FLUSH
         if(this%is_active(errc)) then
          if(errc.eq.0) then
           call this%release_rsc(errc)
           if(errc.eq.0) then
            this%owner_id=-1
            if(associated(this%cache_entry)) then
             call this%cache_entry%decr_ref_count()
             this%cache_entry=>NULL()
            endif
            this%tensor=>NULL()
           else
            errc=-2
           endif
          else
           errc=-1
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndDestruct
!-------------------------------------
        subroutine TensOprndLock(this)
!Sets the lock for accessing/updating the tensor operand content.
         implicit none
         class(tens_oprnd_t), intent(inout):: this !inout: tensor operand

         if(associated(this%cache_entry)) call this%cache_entry%lock()
         return
        end subroutine TensOprndLock
!---------------------------------------
        subroutine TensOprndUnlock(this)
!Releases the access lock.
         implicit none
         class(tens_oprnd_t), intent(inout):: this !inout: tensor operand

         if(associated(this%cache_entry)) call this%cache_entry%unlock()
         return
        end subroutine TensOprndUnlock
!------------------------------------------------------------
        subroutine TensOprndPrintIt(this,ierr,dev_id,nspaces)
         implicit none
         class(tens_oprnd_t), intent(inout):: this     !in: tensor operand
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD), intent(in), optional:: dev_id  !in: output device id
         integer(INTD), intent(in), optional:: nspaces !in: left alignment
         integer(INTD):: errc,devo,nsp,j,sts
         logical:: actv
!$OMP FLUSH
         errc=0
         devo=6; if(present(dev_id)) devo=dev_id
         nsp=0; if(present(nspaces)) nsp=nspaces
         call this%lock()
         do j=1,nsp; write(devo,'(" ")',ADVANCE='NO'); enddo
         write(devo,'("TENSOR OPERAND{")')
         actv=this%is_active(); sts=this%get_comm_stat()
         do j=1,nsp+1; write(devo,'(" ")',ADVANCE='NO'); enddo
         write(devo,'("Active = ",l1,"; Communication = ",i2,"; Metadata owner = ",i6)') actv,sts,this%owner_id
         if(associated(this%tensor)) then
          call this%tensor%print_it(errc,devo,nsp+1); if(errc.ne.TEREC_SUCCESS) errc=-1
         else
          do j=1,nsp+1; write(devo,'(" ")',ADVANCE='NO'); enddo
          write(devo,'("No Tensor!")')
         endif
         do j=1,nsp; write(devo,'(" ")',ADVANCE='NO'); enddo
         write(devo,'("}")')
         call this%unlock()
         flush(devo)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndPrintIt
!---------------------------------------
        subroutine tens_oprnd_dtor(this)
         implicit none
         type(tens_oprnd_t):: this
         integer(INTD):: errc

         call this%destruct(errc)
         if(errc.ne.0) then
          if(VERBOSE) then
           write(CONS_OUT,'("#ERROR(TAVP-MNG:tens_oprnd_dtor): Destruction error ",i11)') errc
           call this%print_it(dev_id=CONS_OUT)
           flush(CONS_OUT)
          endif
          call quit(errc,'#FATAL(TAVP-MNG:tens_oprnd_dtor): Tensor operand destructor failed!')
         endif
         return
        end subroutine tens_oprnd_dtor
![tens_instr_t]==============================================================
        subroutine TensInstrCtor(this,op_code,ierr,op_spec,iid,stat,err_code)
!Constructs a tensor instruction from a given tensor operation.
!The tensor instruction is a realization of a given tensor operation
!for a specific TAVP kind. The tensors in the tensor operation
!must be persistent (associated by aggregation)). Note that the
!tensor operands will not have any information on tensor ownership,
!it should be set separately, if needed.
         implicit none
         class(tens_instr_t), intent(inout):: this        !out: tensor instruction (must be empty on entrance)
         integer(INTD), intent(in):: op_code              !in: instruction code (see top of this module)
         integer(INTD), intent(out), optional:: ierr      !out: error code
         class(*), intent(in), target, optional:: op_spec !in: formal operation specification
         integer(INTL), intent(in), optional:: iid        !in: instruction id (>=0)
         integer(INTD), intent(in), optional:: stat       !in: instruction status to set (defaults to DS_INSTR_NEW)
         integer(INTD), intent(in), optional:: err_code   !in: instruction error code to set (defaults to DSVP_SUCCESS)
         integer(INTD):: errc,ier

         if(this%is_empty(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
!Construct the instruction:
           select case(op_code)
           case(TAVP_INSTR_NOOP)
           case(TAVP_INSTR_CTRL_RESUME,TAVP_INSTR_CTRL_STOP,TAVP_INSTR_CTRL_DUMP_CACHE)
            call construct_instr_ctrl(errc); if(errc.ne.0) errc=-11
           case(TAVP_INSTR_TENS_CREATE,TAVP_INSTR_TENS_DESTROY)
            call construct_instr_tens_create_destroy(errc); if(errc.ne.0) errc=-10
           case(TAVP_INSTR_TENS_INIT)
            call construct_instr_tens_transform(errc); if(errc.ne.0) errc=-9
           case(TAVP_INSTR_TENS_CONTRACT)
            call construct_instr_tens_contract(errc); if(errc.ne.0) errc=-8
           case default
            errc=-7 !invalid instruction opcode (or not implemented)
           end select
!Activate the instruction:
           if(errc.eq.0) then
            ier=DSVP_SUCCESS; if(present(err_code)) ier=err_code
            if(present(stat)) then
             if(present(iid)) then
              call this%activate(op_code,errc,stat=stat,err_code=ier,iid=iid); if(errc.ne.0) errc=-6
             else
              call this%activate(op_code,errc,stat=stat,err_code=ier); if(errc.ne.0) errc=-5
             endif
            else
             if(present(iid)) then
              call this%activate(op_code,errc,err_code=ier,iid=iid); if(errc.ne.0) errc=-4
             else
              call this%activate(op_code,errc,err_code=ier); if(errc.ne.0) errc=-3
             endif
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

         subroutine construct_instr_ctrl(jerr)
          !TAVP CTRL instruction:
          !No op_spec, no control, no operands
          integer(INTD), intent(out):: jerr

          jerr=0
          return
         end subroutine construct_instr_ctrl

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
              call tens_oprnd%tens_oprnd_ctor(tensor,jerr) !tensor owner id is omitted here
              if(jerr.eq.0) then
               oprnd=>tens_oprnd; call this%set_operand(0,oprnd,jerr) !ownership transfer for oprnd=tens_oprnd
               if(jerr.eq.DSVP_SUCCESS) then
                if(op_code.eq.TAVP_INSTR_TENS_CREATE) then
                 this%num_out_oprnds=1; this%out_oprnds(0:this%num_out_oprnds-1)=(/0/) !operand 0 is output
                elseif(op_code.eq.TAVP_INSTR_TENS_DESTROY) then
                 this%num_out_oprnds=0 !no output operands
                endif
               else
                deallocate(tens_oprnd); jerr=-6
               endif
               oprnd=>NULL() !<oprnd> pointer was saved in the tensor instruction and will later be deallocated
              else
               deallocate(tens_oprnd); jerr=-5
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

         subroutine construct_instr_tens_transform(jerr)
          !INIT/TRANSFORM a tensor:
          !op_spec={tens_transformation_t}
          integer(INTD), intent(out):: jerr
          class(ds_oprnd_t), pointer:: oprnd
          class(ds_instr_ctrl_t), pointer:: instr_ctrl
          class(tens_transformation_t), pointer:: tens_trans
          class(ctrl_tens_trans_t), pointer:: tens_trans_ctrl
          class(tens_rcrsv_t), pointer:: tensor
          class(tens_oprnd_t), pointer:: tens_oprnd
          character(:), allocatable:: method_name
          complex(8):: scalar
          logical:: defined

          jerr=0; tens_trans=>NULL()
          select type(op_spec); class is(tens_transformation_t); tens_trans=>op_spec; end select
          if(associated(tens_trans)) then
           if(tens_trans%is_set()) then
            call tens_trans%get_method(defined,scalar,method_name,jerr)
            if(jerr.eq.TEREC_SUCCESS) then
             allocate(tens_trans_ctrl,STAT=jerr)
             if(jerr.eq.0) then
              if(allocated(method_name)) then
               call tens_trans_ctrl%ctrl_tens_trans_ctor(jerr,scalar,method_name)
              else
               call tens_trans_ctrl%ctrl_tens_trans_ctor(jerr,scalar)
              endif
              if(jerr.eq.0) then
               instr_ctrl=>tens_trans_ctrl; call this%set_control(instr_ctrl,jerr) !ownership transfer for instr_ctrl=tens_trans_ctrl
               if(jerr.eq.DSVP_SUCCESS) then
                call this%alloc_operands(1,jerr)
                if(jerr.eq.DSVP_SUCCESS) then
                 tensor=>tens_trans%get_argument(0,jerr)
                 if(jerr.eq.TEREC_SUCCESS) then
                  allocate(tens_oprnd,STAT=jerr)
                  if(jerr.eq.0) then
                   call tens_oprnd%tens_oprnd_ctor(tensor,jerr)
                   if(jerr.eq.0) then
                    oprnd=>tens_oprnd; call this%set_operand(0,oprnd,jerr) !ownership transfer for oprnd=tens_oprnd
                    if(jerr.eq.DSVP_SUCCESS) then
                     this%num_out_oprnds=1; this%out_oprnds(0:this%num_out_oprnds-1)=(/0/) !operand #0 is output
                    else
                     deallocate(tens_oprnd); jerr=-11
                    endif
                    oprnd=>NULL()
                   else
                    deallocate(tens_oprnd); jerr=-10
                   endif
                   tens_oprnd=>NULL()
                  else
                   jerr=-9
                  endif
                 else
                  jerr=-8
                 endif
                 tensor=>NULL()
                else
                 jerr=-7
                endif
               else
                deallocate(tens_trans_ctrl); jerr=-6
               endif
               instr_ctrl=>NULL()
              else
               deallocate(tens_trans_ctrl); jerr=-5
              endif
              tens_trans_ctrl=>NULL()
             else
              jerr=-4
             endif
            else
             jerr=-3
            endif
           else
            jerr=-2
           endif
           tens_trans=>NULL()
          else
           jerr=-1
          endif
          return
         end subroutine construct_instr_tens_transform

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
               instr_ctrl=>tens_contr_ctrl; call this%set_control(instr_ctrl,jerr) !ownership transfer for instr_ctrl=tens_contr_ctrl
               if(jerr.eq.DSVP_SUCCESS) then
                call this%alloc_operands(3,jerr)
                if(jerr.eq.DSVP_SUCCESS) then
                 do jj=0,2
                  tensor=>tens_contr%get_argument(jj,jerr); if(jerr.ne.TEREC_SUCCESS) exit
                  allocate(tens_oprnd,STAT=jerr); if(jerr.ne.0) exit
                  call tens_oprnd%tens_oprnd_ctor(tensor,jerr); if(jerr.ne.0) then; deallocate(tens_oprnd); exit; endif
                  oprnd=>tens_oprnd; call this%set_operand(jj,oprnd,jerr) !ownership transfer for oprnd=tens_oprnd
                  if(jerr.ne.DSVP_SUCCESS) then; oprnd=>NULL(); deallocate(tens_oprnd); exit; endif
                  oprnd=>NULL(); tens_oprnd=>NULL(); tensor=>NULL() !<oprnd> pointer was saved in the tensor instruction and will later be deallocated
                 enddo
                 this%num_out_oprnds=1; this%out_oprnds(0:this%num_out_oprnds-1)=(/0/) !operand 0 is output
                else
                 jerr=-7
                endif
               else
                deallocate(tens_contr_ctrl); jerr=-6
               endif
               instr_ctrl=>NULL()
              else
               deallocate(tens_contr_ctrl); jerr=-5
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
!Encodes a tensor instruction into a bytecode packet:
! 0. Instruction id;
! 1. Instruction code;
! 2. Instruction status;
! 3. Instruction error code;
! 4. Instruction control field (optional);
! 5. Instruction operands (optional): {Owner_id,Read_count,Write_count,Tensor} for each tensor operand.
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
                 case(TAVP_INSTR_CTRL_RESUME,TAVP_INSTR_CTRL_STOP,TAVP_INSTR_CTRL_DUMP_CACHE)
                  call encode_instr_ctrl(errc); if(errc.ne.0) errc=-13
                 case(TAVP_INSTR_TENS_CREATE,TAVP_INSTR_TENS_DESTROY)
                  call encode_instr_tens_create_destroy(errc); if(errc.ne.0) errc=-12
                 case(TAVP_INSTR_TENS_INIT)
                  call encode_instr_tens_transform(errc); if(errc.ne.0) errc=-11
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

         subroutine encode_instr_ctrl(jerr)
          !TAVP CTRL instruction:
          !Packet format: {id|op_code|status|error}
          integer(INTD), intent(out):: jerr

          jerr=0
          return
         end subroutine encode_instr_ctrl

         subroutine encode_instr_tens_create_destroy(jerr)
          !CREATE/DESTROY a tensor:
          !Packet format: {id|op_code|status|error|tensor_operand}
          integer(INTD), intent(out):: jerr
          class(ds_oprnd_t), pointer:: oprnd
          class(tens_rcrsv_t), pointer:: tensor

          jerr=0
          oprnd=>this%get_operand(0,jerr)
          if(jerr.eq.DSVP_SUCCESS) then
           select type(oprnd)
           class is(tens_oprnd_t)
            call oprnd%lock()
            tensor=>oprnd%get_tensor(jerr)
            if(jerr.eq.0) then
             if(tensor%is_set()) then
              call pack_builtin(instr_packet,oprnd%get_owner_id(),jerr)
              if(jerr.eq.0) call pack_builtin(instr_packet,oprnd%get_read_count(),jerr)
              if(jerr.eq.0) call pack_builtin(instr_packet,oprnd%get_write_count(),jerr)
              if(jerr.eq.0) call tensor%pack(instr_packet,jerr)
              if(jerr.ne.0) jerr=-5
             else
              jerr=-4
             endif
             tensor=>NULL()
            else
             jerr=-3
            endif
            call oprnd%unlock()
           class default
            jerr=-2
           end select
           oprnd=>NULL()
          else
           jerr=-1
          endif
          return
         end subroutine encode_instr_tens_create_destroy

         subroutine encode_instr_tens_transform(jerr)
          !TRANSFORM/INIT a tensor: tensor0=value
          !Packed format: {id|op_code|status|error|ctrl_tens_trans_t|tensor_operand0}
          integer(INTD), intent(out):: jerr
          class(ds_oprnd_t), pointer:: oprnd
          class(tens_rcrsv_t), pointer:: tensor
          class(ds_instr_ctrl_t), pointer:: tens_trans_ctrl

          jerr=0
          tens_trans_ctrl=>this%get_control(jerr)
          if(jerr.eq.DSVP_SUCCESS) then
           select type(tens_trans_ctrl)
           class is(ctrl_tens_trans_t)
            call tens_trans_ctrl%pack(instr_packet,jerr); if(jerr.ne.0) jerr=-8
           class default
            jerr=-7
           end select
           if(jerr.eq.0) then
            tensor=>NULL()
            oprnd=>this%get_operand(0,jerr)
            if(jerr.eq.DSVP_SUCCESS) then
             select type(oprnd)
             class is(tens_oprnd_t)
              call oprnd%lock()
              tensor=>oprnd%get_tensor(jerr)
              if(jerr.eq.0) then
               if(tensor%is_set()) then !trap
                call pack_builtin(instr_packet,oprnd%get_owner_id(),jerr)
                if(jerr.eq.0) call pack_builtin(instr_packet,oprnd%get_read_count(),jerr)
                if(jerr.eq.0) call pack_builtin(instr_packet,oprnd%get_write_count(),jerr)
                if(jerr.eq.0) call tensor%pack(instr_packet,jerr)
                if(jerr.ne.0) jerr=-6
               else
                jerr=-5
               endif
              else
               jerr=-4
              endif
              tensor=>NULL()
              call oprnd%unlock()
             class default
              jerr=-3
             end select
            else
             jerr=-2
            endif
            oprnd=>NULL()
           endif
           tens_trans_ctrl=>NULL()
          else
           jerr=-1
          endif
          return
         end subroutine encode_instr_tens_transform

         subroutine encode_instr_tens_contract(jerr)
          !CONTRACT two tensors: tensor0+=tensor1*tensor2*scalar:
          !Packed format: {id|op_code|status|error|ctrl_tens_contr_t|tensor_operand0,tensor_operand1,tensor_operand2}
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
            tensor=>NULL()
            do jj=0,2 !loop over the tensor instruction operands: 0:Destination, 1:LeftInput, 2:RightInput
             oprnd=>this%get_operand(jj,jerr); if(jerr.ne.DSVP_SUCCESS) then; jerr=-6; exit; endif
             select type(oprnd)
             class is(tens_oprnd_t)
              call oprnd%lock()
              tensor=>oprnd%get_tensor(jerr)
              if(jerr.eq.0) then
               if(.not.tensor%is_set()) then; jerr=-5; exit; endif !trap
               call pack_builtin(instr_packet,oprnd%get_owner_id(),jerr)
               if(jerr.eq.0) call pack_builtin(instr_packet,oprnd%get_read_count(),jerr)
               if(jerr.eq.0) call pack_builtin(instr_packet,oprnd%get_write_count(),jerr)
               if(jerr.eq.0) call tensor%pack(instr_packet,jerr)
               if(jerr.ne.0) then; jerr=-4; exit; endif
               tensor=>NULL()
              else
               tensor=>NULL(); call oprnd%unlock(); jerr=-3; exit
              endif
              call oprnd%unlock()
             class default
              jerr=-2; exit
             end select
             oprnd=>NULL()
            enddo
            if(associated(tensor)) then !in case of error
             select type(oprnd); class is(tens_oprnd_t); call oprnd%unlock(); end select
             tensor=>NULL()
            endif
            oprnd=>NULL()
           endif
           tens_contr_ctrl=>NULL()
          else
           jerr=-1
          endif
          return
         end subroutine encode_instr_tens_contract

        end subroutine TensInstrEncode
!------------------------------------------------------------------------------------------
        function TensInstrFullyLocated(this,ierr,input_located,input_ready) result(located)
!Returns TRUE if the tensor instruction has been fully located, FALSE otherwise.
!Being fully located means that each tensor operand has been located (its info available).
!<input_located> is set to TRUE when all input tensor operands have been located.
!<input_ready> is set to TRUE when each input tensor operand is currently defined,
!that is, it is neither undefined nor being updated. Note that <input_ready> is
!generally independent of <located> (an input tensor operand can be ready,
!yet not located), although it is probably unlikely.
         implicit none
         logical:: located                              !out: fully located or not
         class(tens_instr_t), intent(in):: this         !in: tensor instruction
         integer(INTD), intent(out), optional:: ierr    !out: error code
         logical, intent(out), optional:: input_located !out: TRUE if all input tensors have been located
         logical, intent(out), optional:: input_ready   !out: TRUE if all input tensors are ready (their values are defined)
         integer(INTD):: errc,n,i
         integer(INTD):: arg_located(0:MAX_TENSOR_OPERANDS-1),arg_ready(0:MAX_TENSOR_OPERANDS-1)
         class(ds_oprnd_t), pointer:: oprnd
         logical:: ilocated,iready,vald

         located=.FALSE.; ilocated=.FALSE.; iready=.FALSE.
         if(.not.this%is_empty(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           n=this%get_num_operands(errc)
           if(errc.eq.DSVP_SUCCESS) then
            arg_located(0:n-1)=0; arg_ready(0:n-1)=0
            do i=0,n-1 !loop over tensor operands
             oprnd=>this%get_operand(i,errc); if(errc.ne.DSVP_SUCCESS) then; errc=-3; exit; endif
             if(oprnd%is_located(errc,valued=vald)) arg_located(i)=1; if(errc.ne.0) then; errc=-2; exit; endif
             if(vald) arg_ready(i)=1
            enddo
            if(errc.eq.0) then
             located=.TRUE.; ilocated=.TRUE.; iready=.TRUE.
             do i=0,n-1
              located=(arg_located(i).eq.1); if(.not.located) exit
             enddo
             arg_located(this%out_oprnds(0:this%num_out_oprnds-1))=1 !ingore output operands
             arg_ready(this%out_oprnds(0:this%num_out_oprnds-1))=1 !ingore output operands
             do i=0,n-1
              ilocated=ilocated.and.(arg_located(i).eq.1)
              iready=iready.and.(arg_ready(i).eq.1)
             enddo
            endif
           endif
          endif
         else
          errc=-1
         endif
         if(present(input_located)) input_located=ilocated
         if(present(input_ready)) input_ready=iready
         if(present(ierr)) ierr=errc
         return
        end function TensInstrFullyLocated
!-------------------------------------------------------------------------------
        function TensInstrGetOutputOperands(this,ierr,num_oprs) result(out_oprs)
!Returns the list of the output tensor operands by their positions.
         implicit none
         integer(INTD), pointer:: out_oprs(:)            !out: positions of the output tensor instruction operands
         class(tens_instr_t), intent(in), target:: this  !in: tensor instruction
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD), intent(out), optional:: num_oprs !out: number of the output tensor instruction operands
         integer(INTD):: errc,n

         out_oprs=>NULL(); n=0
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           n=this%num_out_oprnds
           if(n.gt.0) out_oprs(0:)=>this%out_oprnds(0:n-1)
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(num_oprs)) num_oprs=n
         if(present(ierr)) ierr=errc
         return
        end function TensInstrGetOutputOperands
!-------------------------------------------------------------------------------------
        function TensInstrGetFlops(this,ierr,arithm_intensity,tot_words) result(flops)
!Returns an estimate of the Flop count as well as arithmetic intensity and total word count.
         implicit none
         real(8):: flops                                   !out: estimate of the total number of Flops
         class(tens_instr_t), intent(in):: this            !in: active tensor instruction
         integer(INTD), intent(out), optional:: ierr       !out: error code
         real(8), intent(out), optional:: arithm_intensity !out: arithmetic intensity estimate (Flops/words ratio)
         real(8), intent(out), optional:: tot_words        !out: total words estimate
         integer(INTD):: errc,opcode,i,j,n
         integer(INTL):: dims(1:MAX_TENSOR_RANK)
         class(ds_oprnd_t), pointer:: tens_oprnd
         class(tens_rcrsv_t), pointer:: tensor
         real(8):: vol,tvol,words

         flops=0d0; words=0d0
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           opcode=this%get_code(errc)
           if(errc.eq.DSVP_SUCCESS) then
            select case(opcode)
            case(TAVP_INSTR_TENS_CONTRACT)
             tvol=1d0
             oloop: do i=0,2 !loop over tensor operands
              tens_oprnd=>this%get_operand(i,errc); if(errc.ne.DSVP_SUCCESS) then; errc=-7; exit oloop; endif
              select type(tens_oprnd)
              class is(tens_oprnd_t)
               call tens_oprnd%lock()
               tensor=>tens_oprnd%get_tensor(errc)
               if(errc.eq.0) then
                call tensor%get_dims(dims,n,errc); if(errc.ne.TEREC_SUCCESS) errc=-6
               else
                errc=-5
               endif
               tensor=>NULL()
               call tens_oprnd%unlock()
               if(errc.eq.0) then
                vol=1d0; do j=1,n; vol=vol*real(dims(j),8); enddo
                tvol=tvol*vol; words=words+vol
               else
                exit oloop
               endif
              class default
               errc=-4; exit oloop
              end select
             enddo oloop
             if(errc.eq.0) flops=dsqrt(tvol)*2d0 !factor of 2 because of additions (along with multiplications)
            case default
             !`Implement Flop counting for other relevant tensor instructions
             flops=1d0 !default (but meaningless) value
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
         if(present(arithm_intensity)) then
          if(errc.eq.0.and.words.gt.0d0) then
           arithm_intensity=flops/words
          else
           arithm_intensity=-1d0
          endif
         endif
         if(present(tot_words).and.errc.eq.0) tot_words=words
         if(present(ierr)) ierr=errc
         return
        end function TensInstrGetFlops
!-----------------------------------------------------------------
        subroutine TensInstrGetOperation(this,tens_operation,ierr)
!Given a tensor instruction, returns back the encapsulated tensor operation
!expressed in terms of the very same tensors (by pointer association).
         implicit none
         class(tens_instr_t), intent(in):: this                                     !in: tensor instruction
         class(tens_operation_t), allocatable, target, intent(out):: tens_operation !out: corresponding (encapsulated) tensor operation
         integer(INTD), intent(out), optional:: ierr                                !out: error code
         integer(INTD):: errc,opcode

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           opcode=this%get_code(errc)
           if(errc.eq.DSVP_SUCCESS) then
            select case(opcode)
            case(TAVP_INSTR_TENS_CREATE) !no associated tensor operation
            case(TAVP_INSTR_TENS_DESTROY) !no associated tensor operation
            case(TAVP_INSTR_TENS_INIT)
             call get_tens_transformation(errc); if(errc.ne.0) errc=-6
            case(TAVP_INSTR_TENS_CONTRACT)
             call get_tens_contraction(errc); if(errc.ne.0) errc=-5
            case default
             errc=-4
             call quit(errc,'#FATAL(TAVP-MNG:tens_instr_t.get_operation): Not implemented!') !`Implement for other relevant tensor instructions
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
         if(errc.ne.0.and.allocated(tens_operation)) deallocate(tens_operation)
         if(present(ierr)) ierr=errc
         return

         contains

          subroutine get_tens_transformation(jerr)
           integer(INTD), intent(out):: jerr
           integer(INTD):: sl
           character(EXA_MAX_METHOD_NAME_LEN):: method_name
           complex(8):: alpha
           class(tens_rcrsv_t), pointer:: tensor
           class(ds_oprnd_t), pointer:: oprnd
           class(ds_instr_ctrl_t), pointer:: instr_ctrl

           jerr=0
           allocate(tens_transformation_t::tens_operation)
           select type(tens_operation)
           class is(tens_transformation_t)
            oprnd=>this%get_operand(0,jerr)
            if(jerr.eq.DSVP_SUCCESS) then
             select type(oprnd)
             class is(tens_oprnd_t)
              call oprnd%lock()
              tensor=>oprnd%get_tensor(jerr)
              if(jerr.eq.0.and.associated(tensor)) then
               call tens_operation%set_argument(tensor,jerr); if(jerr.ne.TEREC_SUCCESS) jerr=-9
               tensor=>NULL()
              else
               jerr=-8
              endif
              call oprnd%unlock()
             class default
              jerr=-7
             end select
             if(jerr.eq.0) then
              instr_ctrl=>this%get_control(jerr)
              if(jerr.eq.DSVP_SUCCESS) then
               select type(instr_ctrl)
               class is(ctrl_tens_trans_t)
                call instr_ctrl%get_method(method_name,sl,jerr,scalar=alpha)
                if(jerr.eq.0) then
#if !(defined(__GNUC__) && __GNUC__ < 8)
                 call tens_operation%set_method(jerr,alpha,.FALSE.,method_name(1:sl),method_map_f)
#else
                 call tens_operation%set_method(jerr,alpha,.FALSE.,method_name(1:sl))
#endif
                 if(jerr.ne.TEREC_SUCCESS) jerr=-6
                else
                 jerr=-5
                endif
               class default
                jerr=-4
               end select
              else
               jerr=-3
              endif
              instr_ctrl=>NULL()
             endif
            else
             jerr=-2
            endif
           class default
            jerr=-1
           end select
           return
          end subroutine get_tens_transformation

          subroutine get_tens_contraction(jerr)
           integer(INTD), intent(out):: jerr
           integer(INTD):: numo,i
           complex(8):: alpha
           class(tens_rcrsv_t), pointer:: tensor
           class(ds_oprnd_t), pointer:: oprnd
           class(ds_instr_ctrl_t), pointer:: instr_ctrl
           type(contr_ptrn_ext_t), pointer:: contr_ptrn

           jerr=0
           allocate(tens_contraction_t::tens_operation)
           select type(tens_operation)
           class is(tens_contraction_t)
            numo=this%get_num_operands(jerr)
            if(jerr.eq.DSVP_SUCCESS) then
             do i=0,numo-1
              oprnd=>this%get_operand(i,jerr); if(jerr.ne.DSVP_SUCCESS) then; jerr=-10; exit; endif
              select type(oprnd)
              class is(tens_oprnd_t)
               call oprnd%lock()
               tensor=>oprnd%get_tensor(jerr)
               if(jerr.eq.0.and.associated(tensor)) then
                call tens_operation%set_argument(tensor,jerr); if(jerr.ne.TEREC_SUCCESS) jerr=-9
                tensor=>NULL()
               else
                jerr=-8
               endif
               call oprnd%unlock()
              class default
               jerr=-7
              end select
              if(jerr.ne.0) exit
             enddo
            else
             jerr=-6
            endif
            if(jerr.eq.0) then
             instr_ctrl=>this%get_control(jerr)
             if(jerr.eq.DSVP_SUCCESS) then
              select type(instr_ctrl)
              class is(ctrl_tens_contr_t)
               contr_ptrn=>instr_ctrl%get_contr_ptrn(jerr,alpha)
               if(jerr.eq.0) then
                call tens_operation%set_contr_ptrn(contr_ptrn,jerr,alpha); if(jerr.ne.TEREC_SUCCESS) jerr=-5
               else
                jerr=-4
               endif
               nullify(contr_ptrn)
              class default
               jerr=-3
              end select
             else
              jerr=-2
             endif
             nullify(instr_ctrl)
            endif
           class default
            jerr=-1
           end select
           return
          end subroutine get_tens_contraction

        end subroutine TensInstrGetOperation
!------------------------------------------------------------
        subroutine TensInstrPrintIt(this,ierr,dev_id,nspaces)
         implicit none
         class(tens_instr_t), intent(in):: this        !in: tensor instruction
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD), intent(in), optional:: dev_id  !in: output device id
         integer(INTD), intent(in), optional:: nspaces !in: left alignment
         integer(INTD):: errc,ier,devo,nsp,i,n,opcode,sts
         integer(INTL):: iid
         class(ds_instr_ctrl_t), pointer:: ctrl
         class(ds_oprnd_t), pointer:: oprnd

         errc=0
         devo=6; if(present(dev_id)) devo=dev_id
         nsp=0; if(present(nspaces)) nsp=nspaces
!$OMP CRITICAL (TAVP_PRINT)
         do i=1,nsp; write(devo,'(" ")',ADVANCE='NO'); enddo
         write(devo,'("TENSOR INSTRUCTION{")')
         iid=this%get_id(); opcode=this%get_code(); sts=this%get_status(errc,ier)
         do i=1,nsp+1; write(devo,'(" ")',ADVANCE='NO'); enddo
         write(devo,'("id = ",i11,"; opcode = ",i4,"; stat = ",i6,"; err = ",i11)') iid,opcode,sts,ier
         n=this%get_num_operands(errc)
         if(errc.eq.DSVP_SUCCESS) then
          ctrl=>this%get_control()
          if(associated(ctrl)) call ctrl%print_it(errc,devo,nsp+1)
          do i=0,n-1
           oprnd=>NULL(); oprnd=>this%get_operand(i)
           if(associated(oprnd)) call oprnd%print_it(errc,devo,nsp+1)
          enddo
         else
          do i=1,nsp+1; write(devo,'(" ")',ADVANCE='NO'); enddo
          write(devo,'("Error: Unable to determine the number of operands!")')
         endif
         do i=1,nsp; write(devo,'(" ")',ADVANCE='NO'); enddo
         write(devo,'("}")')
         flush(devo)
!$OMP END CRITICAL (TAVP_PRINT)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensInstrPrintIt
!-----------------------------------------------------------------------------------
        subroutine TensInstrExtractCacheEntries(this,cache_entries,num_entries,ierr)
!Returns an array of references to tensor cache entries associated with the tensor operands for subsequent eviction.
!Only temporary tensor cache entries with reference count of 1 are returned.
         implicit none
         class(tens_instr_t), intent(inout):: this                     !in: tensor instruction
         type(tens_entry_mng_ref_t), intent(inout):: cache_entries(1:) !out: references to the tensor cache entries to be evicted
         integer(INTD), intent(out):: num_entries                      !out: number of tensor cache entries returned
         integer(INTD), intent(out), optional:: ierr                   !out: error code
         integer(INTD):: errc,i,n
         class(ds_oprnd_t), pointer:: oprnd
         class(tens_entry_mng_t), pointer:: tens_entry

         num_entries=0
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           n=this%get_num_operands(errc)
           if(errc.eq.DSVP_SUCCESS) then
            do i=0,n-1
             oprnd=>this%get_operand(i,errc); if(errc.ne.DSVP_SUCCESS) then; errc=-6; exit; endif
             select type(oprnd)
             class is(tens_oprnd_t)
              call oprnd%lock()
              tens_entry=>oprnd%get_cache_entry(errc)
              if(errc.eq.0.and.associated(tens_entry)) then
               !if(tens_entry%get_ref_count().le.1.and.tens_entry%get_use_count().eq.0.and.&
                 !&(.not.tens_entry%is_persistent())) then
               if((tens_entry%get_ref_count().le.1).and.(.not.tens_entry%is_persistent())) then
                call tens_entry%incr_use_count() !to protect from repeated evictions (will be decremented by .evict())
                num_entries=num_entries+1; cache_entries(num_entries)%cache_entry=>tens_entry
               endif
              else
               errc=-5
              endif
              tens_entry=>NULL()
              call oprnd%unlock()
              if(errc.ne.0) exit
             class default
              errc=-4; exit
             end select
            enddo
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
        end subroutine TensInstrExtractCacheEntries
!-------------------------------------------------------
        subroutine TensInstrRemovePersistency(this,ierr)
!Removes the persistent status from the tensor instruction operands.
         implicit none
         class(tens_instr_t), intent(inout):: this   !in: active tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,n,i
         class(ds_oprnd_t), pointer:: oprnd
         class(tens_entry_mng_t), pointer:: cache_entry

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           n=this%get_num_operands(errc)
           if(errc.eq.DSVP_SUCCESS) then
            if(n.gt.0) then
             do i=0,n-1
              oprnd=>this%get_operand(i,errc)
              if(errc.eq.DSVP_SUCCESS.and.associated(oprnd)) then
               select type(oprnd)
               class is(tens_oprnd_t)
                call oprnd%lock()
                cache_entry=>oprnd%get_cache_entry(errc)
                if(errc.eq.0) then
                 call cache_entry%set_persistency(.FALSE.); cache_entry=>NULL()
                else
                 errc=-6
                endif
                call oprnd%unlock()
               class default
                errc=-5
               end select
              else
               errc=-4
              endif
              if(errc.ne.0) exit
             enddo
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
        end subroutine TensInstrRemovePersistency
!---------------------------------------
        subroutine tens_instr_dtor(this)
!DTOR: Expects tensor instruction status to be either DS_INSTR_EMPTY or DS_INSTR_RETIRED.
         implicit none
         type(tens_instr_t):: this !inout: empty or retired tensor instruction
         integer(INTD):: errc,sts,opcode

         opcode=this%get_code(); sts=this%get_status(errc)
         if((sts.eq.DS_INSTR_EMPTY.or.sts.eq.DS_INSTR_RETIRED).and.errc.eq.DSVP_SUCCESS) then
          call this%clean(errc)
          if(errc.ne.DSVP_SUCCESS) call quit(errc,'#ERROR(TAVP-MNG:tens_instr_dtor): Tensor instruction destruction failed!')
         else
          if(errc.eq.DSVP_SUCCESS.and.VERBOSE) then
           write(CONS_OUT,'("#ERROR(TAVP-MNG:tens_instr_dtor): TAVP instruction is still active: code = ",i5,", stat = ",i5)')&
           &opcode,sts
           call this%print_it(dev_id=CONS_OUT)
           flush(CONS_OUT)
          endif
          !call crash() !debug
          call quit(-1,'#FATAL(TAVP-MNG:tens_instr_dtor): Tensor instruction destructor failed!')
         endif
         return
        end subroutine tens_instr_dtor
!----------------------------------------------------------------
        function tens_instr_locator_mark(tens_instr) result(pred)
!Returns GFC_TRUE if the tensor instruction is a special markup instruction in the location cycle.
         implicit none
         integer(INTD):: pred                      !out: {GFC_FALSE,GFC_TRUE,GFC_ERROR}
         class(*), intent(in), target:: tens_instr !in: tensor instruction
         integer(INTD):: sts

         pred=GFC_FALSE
         select type(tens_instr)
         class is(tens_instr_t)
          sts=tens_instr%get_status()
          if(tens_instr%get_code().eq.TAVP_INSTR_CTRL_RESUME.and.(sts.eq.DS_INSTR_SPECIAL.or.sts.eq.DS_INSTR_TERMINAL))&
          &pred=GFC_TRUE
         class default
          pred=GFC_ERROR
         end select
         return
        end function tens_instr_locator_mark
!---------------------------------------------------------
        function tens_instr_print(tens_instr) result(ierr)
!GFC printer for tens_instr_t.
         implicit none
         integer(INTD):: ierr
         class(*), intent(inout), target:: tens_instr

         ierr=0
         select type(tens_instr)
         class is(tens_instr_t)
          call tens_instr%print_it(ierr,dev_id=CONS_OUT)
         class default
          ierr=-1
         end select
         return
        end function tens_instr_print
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
           this%msg_tag=conf%msg_tag
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
         class(tavp_mng_decoder_t), intent(inout):: this !inout: TAVP-MNG Decoder DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier,thid,num_packets,opcode,sts,port_id,i,j,uid
         logical:: active,stopping,new
         class(dsvp_t), pointer:: dsvp
         class(tavp_mng_t), pointer:: tavp
         class(*), pointer:: uptr
         class(ds_unit_t), pointer:: acceptor
         type(obj_pack_t):: instr_packet
         type(comm_handle_t):: comm_hl
         type(tens_instr_t), pointer:: tens_instr
         type(tens_instr_t):: tens_instr_empty

         errc=0; thid=omp_get_thread_num(); uid=this%get_id()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decoder started as DSVU # ",i2," (thread ",i2,"): Listening to ",i11,1x,i6)')&
          &impir,uid,thid,this%source_comm,this%source_rank
          flush(CONS_OUT)
         endif
!Reserve a bytecode buffer:
         call this%bytecode%reserve_mem(ier,MAX_BYTECODE_SIZE,MAX_BYTECODE_INSTR); if(ier.ne.0.and.errc.eq.0) errc=-41
!Initialize queues and ports:
         call this%init_queue(this%num_ports,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-40
!Initialize the control list:
         ier=this%ctrl_list%init(this%control_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-39
!Set up tensor argument cache and wait on other TAVP units:
         tavp=>NULL(); dsvp=>this%get_dsvp(); select type(dsvp); class is(tavp_mng_t); tavp=>dsvp; end select
         if(associated(tavp)) then
          this%arg_cache=>tavp%tens_cache
          call tavp%sync_units(errc,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-38
         else
          this%arg_cache=>NULL(); if(errc.eq.0) errc=-37
         endif
!Work loop:
         active=((errc.eq.0).and.(this%source_comm.ne.MPI_COMM_NULL)); stopping=(.not.active)
         wloop: do while(active)
          if(.not.stopping) then
 !Receive new bytecode (if posted):
           call comm_hl%clean(ier); if(ier.ne.0.and.errc.eq.0) then; errc=-36; exit wloop; endif
           new=this%bytecode%receive(comm_hl,ier,proc_rank=this%source_rank,tag=this%msg_tag,comm=this%source_comm)
           if(ier.ne.0.and.errc.eq.0) then; errc=-35; exit wloop; endif
           if(new) then !new bytecode is available
            call comm_hl%wait(ier); if(ier.ne.0.and.errc.eq.0) then; errc=-34; exit wloop; endif
            num_packets=this%bytecode%get_num_packets(ier); if(ier.ne.0.and.errc.eq.0) then; errc=-33; exit wloop; endif
            if(DEBUG.gt.1) then
             write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decoder unit ",i2," received ",i9," new instructions")')&
             &impir,uid,num_packets
             flush(CONS_OUT)
            endif
            if(num_packets.gt.0) then
 !Decode new bytecode:
             ier=this%ctrl_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-32; exit wloop; endif
             ier=this%iqueue%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-31; exit wloop; endif
             do i=1,num_packets
  !Extract an instruction:
              call this%bytecode%extract_packet(i,instr_packet,ier,preclean=.TRUE.)
              if(ier.ne.0.and.errc.eq.0) then
               write(CONS_OUT,'("#ERROR(TAVP-MNG:Decoder)[",i6,":",i2,"]: Instruction packet extraction error ",i11,'//&
               &'": ",i11,1x,i11)') impir,uid,ier,i,num_packets
               flush(CONS_OUT)
               errc=-30; exit wloop
              endif
  !Append an empty instruction at the end of the main queue:
              ier=this%iqueue%append(tens_instr_empty); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-29; exit wloop; endif
              ier=this%iqueue%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-28; exit wloop; endif
              uptr=>this%iqueue%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-27; exit wloop; endif
              tens_instr=>NULL(); select type(uptr); type is(tens_instr_t); tens_instr=>uptr; end select
              if(.not.associated(tens_instr).and.errc.eq.0) then; errc=-26; exit wloop; endif !trap
  !Construct the instruction by decoding the extracted instruction:
              call this%decode(tens_instr,instr_packet,ier); if(ier.ne.0.and.errc.eq.0) then; errc=-25; exit wloop; endif
              sts=tens_instr%get_status(ier,j); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-24; exit wloop; endif
              opcode=tens_instr%get_code(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-23; exit wloop; endif
              if(DEBUG.gt.0) then
               if(this%get_id().eq.0) then !uDecoder only
                write(CONS_OUT,'("#DEBUG(TAVP-MNG:uDecoder): Decoded a new tensor instruction:")')
                call tens_instr%print_it(dev_id=CONS_OUT)
                flush(CONS_OUT)
               endif
              endif
  !Clone CONTROL instructions for own port (uDecoder only):
              if(opcode.ge.TAVP_ISA_CTRL_FIRST.and.opcode.le.TAVP_ISA_CTRL_LAST) then !copy control instructions to own port
               if(opcode.eq.TAVP_INSTR_CTRL_STOP.and.i.ne.num_packets.and.errc.eq.0) then; errc=-22; exit wloop; endif !STOP must be the last instruction in the packet
               ier=this%ctrl_list%append(tens_instr); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-21; exit wloop; endif
              endif
             enddo
 !Pass cloned control instructions to own port (uDecoder only):
             ier=this%ctrl_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-20; exit wloop; endif
             if(this%ctrl_list%get_status().eq.GFC_IT_ACTIVE) then
              ier=this%load_port(0,this%ctrl_list); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-19; exit wloop; endif
             endif
 !Pass all decoded instructions to the acceptor unit (whichever it is):
             ier=this%iqueue%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-18; exit wloop; endif
             if(this%iqueue%get_status().eq.GFC_IT_ACTIVE) then
              ier=this%get_acceptor(acceptor,port_id); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-17; exit wloop; endif
              if(associated(acceptor)) then
               ier=acceptor%load_port(port_id,this%iqueue); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-16; exit wloop; endif
               if(this%iqueue%get_status().ne.GFC_IT_EMPTY.and.errc.eq.0) then; errc=-15; exit wloop; endif !trap
              else
               if(errc.eq.0) then; errc=-14; exit wloop; endif
              endif
             endif
            endif
            call this%bytecode%clean(ier); if(ier.ne.0.and.errc.eq.0) then; errc=-13; exit wloop; endif
           endif !new bytecode
          else !stopping
           ier=this%iqueue%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-12; exit wloop; endif
           if(this%iqueue%get_status().eq.GFC_IT_EMPTY) then
            active=.FALSE.
           else
            if(errc.eq.0) then; errc=-11; exit wloop; endif
           endif
           if(.not.this%port_empty(0,ier).and.errc.eq.0) then; errc=-10; exit wloop; endif !trap
          endif
 !Check own port for control instructions:
          ier=this%flush_port(0); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-9; exit wloop; endif
          ier=this%iqueue%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-8; exit wloop; endif
          ier=this%iqueue%get_status()
          do while(ier.eq.GFC_IT_ACTIVE)
           uptr=>NULL(); uptr=>this%iqueue%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-7; exit wloop; endif
           tens_instr=>NULL(); select type(uptr); type is(tens_instr_t); tens_instr=>uptr; end select
           if(.not.associated(tens_instr).and.errc.eq.0) then; errc=-6; exit wloop; endif !trap
           opcode=tens_instr%get_code(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-5; exit wloop; endif
           if(opcode.eq.TAVP_INSTR_CTRL_STOP) stopping=.TRUE. !only STOP instruction is expected
           call tens_instr%set_status(DS_INSTR_RETIRED,ier,DSVP_SUCCESS)
           if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-3; exit wloop; endif
           ier=this%iqueue%next(); ier=this%iqueue%get_status() !in general more control instructions can be expected
          enddo
 !Clear the instruction queue (control instructions only):
          if(ier.eq.GFC_IT_DONE) then !control list was not empty
           ier=this%iqueue%delete_all(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-2; exit wloop; endif
          endif
         enddo wloop
!Record the error:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(TAVP-MNG)[",i6,"]: Decoder error ",i11," by thread ",i2)')&
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
         class(tavp_mng_decoder_t), intent(inout):: this !inout: TAVP-MNG Decoder DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier,thid,uid
         logical:: empt

         errc=0; thid=omp_get_thread_num(); uid=this%get_id()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decoder stopped as DSVU # ",i2," (thread ",i2,")")') impir,uid,thid
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decoder DSVU # ",i2," (thread ",i2,"): Time stats (sec): ")',ADVANCE='NO')&
          &impir,uid,thid; call this%print_timing(CONS_OUT)
          empt=this%port_empty(0)
          !write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decoder DSVU # ",i2,": Port empty = ",l1)') impir,uid,empt !debug
          flush(CONS_OUT)
         endif
!Release the tensor argument cache pointer:
         this%arg_cache=>NULL()
!Deactivate the control list:
         ier=this%ctrl_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%ctrl_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-5
           ier=this%ctrl_list%delete_all()
          endif
          ier=this%ctrl_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-4
         else
          if(errc.eq.0) errc=-3
         endif
!Release queues:
         call this%release_queue(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-2
!Release the bytecode buffer:
         call this%bytecode%destroy(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
!Record an error, if any:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDecoderShutdown
!-----------------------------------------------------------------------
        subroutine TAVPMNGDecoderDecode(this,ds_instr,instr_packet,ierr)
!Decodes a tensor instruction from a plain bytecode.
!The bytecode format is specified in tens_instr_t.encode().
         implicit none
         class(tavp_mng_decoder_t), intent(inout):: this     !inout: TAVP-MNG Decoder
         class(ds_instr_t), intent(inout), target:: ds_instr !out: decoded tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet     !in: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr         !out: error code
         integer(INTD):: errc,op_code,stat,err_code,uid
         integer(INTL):: iid
         real(8):: tm

         tm=thread_wtime(); uid=this%get_id()
         if(ds_instr%is_empty(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(instr_packet%is_healthy(errc)) then
!Retrieve the TAVP argument cache:
            if(associated(this%arg_cache)) then
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
                 case(TAVP_INSTR_CTRL_RESUME,TAVP_INSTR_CTRL_STOP,TAVP_INSTR_CTRL_DUMP_CACHE)
                 case(TAVP_INSTR_TENS_CREATE,TAVP_INSTR_TENS_DESTROY)
                  call decode_instr_tens_create_destroy(errc); if(errc.ne.0) errc=-13
                 case(TAVP_INSTR_TENS_INIT)
                  call decode_instr_tens_transform(errc); if(errc.ne.0) errc=-12
                 case(TAVP_INSTR_TENS_CONTRACT)
                  call decode_instr_tens_contract(errc); if(errc.ne.0) errc=-11
                 case default
                  call ds_instr%set_status(DS_INSTR_RETIRED,errc,TAVP_ERR_GEN_FAILURE)
                  errc=-10 !unknown instruction opcode (or not implemented)
                 end select
!Activate the instruction:
                 if(errc.eq.0) then
                  call ds_instr%activate(op_code,errc,stat,err_code,iid)
                  if(errc.ne.DSVP_SUCCESS) then
                   call ds_instr%set_status(DS_INSTR_RETIRED,errc,TAVP_ERR_GEN_FAILURE)
                   errc=-9
                  endif
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
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(DEBUG.gt.0.and.errc.ne.0) then
          write(CONS_OUT,'("#ERROR(TAVP-MNG:Decoder.decode)[",i6,":",i3,"]: Error ",i11)')&
          &impir,omp_get_thread_num(),errc
          flush(CONS_OUT)
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
           integer(INTD):: jj,jn,jown,jread,jwrite,jpwn
           logical:: stored,updated

           jn=ds_instr%get_num_operands(jerr)
           if(jerr.eq.DSVP_SUCCESS.and.jn.ge.0) then
            tensor_tmp=>NULL()
            do jj=0,jn-1 !loop over tensor operands
             allocate(tensor_tmp,STAT=jerr) !tensor will either be owned by a tensor cache entry or deallocated here
             if(jerr.ne.0) then
              tensor_tmp=>NULL()
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE); jerr=-12; exit
             endif
             call unpack_builtin(instr_packet,jown,jerr) !tensor owner id
             if(jerr.ne.PACK_SUCCESS) then; call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_BTC_BAD); jerr=-11; exit; endif
             call unpack_builtin(instr_packet,jread,jerr) !tensor read access count
             if(jerr.ne.PACK_SUCCESS) then; call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_BTC_BAD); jerr=-10; exit; endif
             call unpack_builtin(instr_packet,jwrite,jerr) !tensor write access count
             if(jerr.ne.PACK_SUCCESS) then; call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_BTC_BAD); jerr=-9; exit; endif
             call tensor_tmp%tens_rcrsv_ctor(instr_packet,jerr) !unpack tensor information into a temporary tensor
             if(jerr.ne.TEREC_SUCCESS) then
              if(DEBUG.gt.0) then
               write(CONS_OUT,'("#ERROR(TAVP-MNG:Decoder.decode.decode_instr_operands)[",i6,":",i3,"]: TensCtor Error ",i11)')&
               &impir,omp_get_thread_num(),jerr
               flush(CONS_OUT)
              endif
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_BTC_BAD); jerr=-8; exit
             endif
             tens_entry=>NULL()
             stored=this%arg_cache%store(tensor_tmp,tens_entry_mng_alloc,jerr,tens_entry_p=tens_entry) !tensor ownership is moved to the tensor cache entry
             if(associated(tens_entry)) then
              call tens_entry%lock()
              if(stored) tensor_tmp=>NULL() !tensor ownership has been transferred to the tensor cache
             else
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_CHE_FAILURE); jerr=-7; exit
             endif
             updated=stored
             tens_mng_entry=>NULL(); select type(tens_entry); class is(tens_entry_mng_t); tens_mng_entry=>tens_entry; end select
             if(.not.associated(tens_mng_entry)) then
              if(associated(tens_entry)) then; call tens_entry%unlock(); call this%arg_cache%release_entry(tens_entry); endif
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_CHE_FAILURE); jerr=-6; exit
             endif
             tensor=>tens_mng_entry%get_tensor() !use the tensor from the tensor cache
             if(.not.stored) then !the tensor was already in the tensor cache before, update it by the information from the just decoded tensor
              if(DEBUG.gt.1) then
               write(CONS_OUT,'("#DEBUG(TAVP-MNG:Decoder): Tensor info update in the tensor cache:")')
               call tensor%print_it(dev_id=CONS_OUT)
               call tensor_tmp%print_it(dev_id=CONS_OUT)
               flush(CONS_OUT)
              endif
              call tensor%update(tensor_tmp,jerr,updated) !tensor metadata update is done inside the tensor cache
              if(DEBUG.gt.1.and.jerr.eq.TEREC_SUCCESS) then
               call tensor%print_it(dev_id=CONS_OUT)
               jpwn=tens_mng_entry%get_owner_id()
               write(CONS_OUT,'("Update status = ",l1,": Metadata exporter = ",i5,", prev = ",i5," (Decoder id = ",i2,")")')&
               &updated,jown,jpwn,uid
               flush(CONS_OUT)
              endif
              deallocate(tensor_tmp); tensor_tmp=>NULL() !deallocate temporary tensor after importing its information into the cache
              if(jerr.ne.TEREC_SUCCESS) then
               if(DEBUG.gt.0) then
                write(CONS_OUT,'("#ERROR(TAVP-MNG:Decoder.decode.decode_instr_operands)[",i6,":",i3,'//&
                &'"]: tens_rcrsv_t.update() error ",i11)') impir,omp_get_thread_num(),jerr
                flush(CONS_OUT)
               endif
               if(associated(tens_entry)) then; call tens_entry%unlock(); call this%arg_cache%release_entry(tens_entry); endif
               call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE); jerr=-5; exit
              endif
             endif
             if(updated.and.jown.ge.0) then !anonymous metadata updates (jown < 0) do not change the current owner
              call tens_mng_entry%set_owner_id(jown) !update the tensor owner id if the tensor info was updated in the cache
              call tens_mng_entry%reset_rw_counters(jread,jwrite) !update the tensor access counters if the tensor info was updated in the cache
             endif
             allocate(tens_oprnd,STAT=jerr) !tensor operand will be owned by the tensor instruction
             if(jerr.ne.0) then
              tens_oprnd=>NULL()
              if(associated(tens_entry)) then; call tens_entry%unlock(); call this%arg_cache%release_entry(tens_entry); endif
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE); jerr=-4; exit
             endif
             call tens_oprnd%tens_oprnd_ctor(tens_mng_entry,jerr) !tensor is owned by the tensor cache (owner id will be imported from the tensor cache)
             if(jerr.ne.0) then
              deallocate(tens_oprnd); tens_oprnd=>NULL()
              if(associated(tens_entry)) then; call tens_entry%unlock(); call this%arg_cache%release_entry(tens_entry); endif
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE); jerr=-3; exit
             endif
             oprnd=>tens_oprnd; call ds_instr%set_operand(jj,oprnd,jerr) !tensor operand ownership is moved to the tensor instruction
             if(jerr.ne.DSVP_SUCCESS) then
              deallocate(tens_oprnd); tens_oprnd=>NULL()
              if(associated(tens_entry)) then; call tens_entry%unlock(); call this%arg_cache%release_entry(tens_entry); endif
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE); jerr=-2; exit
             endif
             tens_oprnd=>NULL()
             if(associated(tens_entry)) then
              call tens_entry%unlock(); call this%arg_cache%release_entry(tens_entry); tens_entry=>NULL()
             endif
            enddo
            if(associated(tensor_tmp)) deallocate(tensor_tmp) !deallocate the temporary tensor (in case of error)
           else
            call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE); jerr=-1
           endif
           if(DEBUG.gt.0.and.jerr.ne.0) then
            write(CONS_OUT,'("#ERROR(TAVP-MNG:Decoder.decode.decode_instr_operands)[",i6,":",i3,"]: Error ",i11)')&
            &impir,omp_get_thread_num(),jerr
            flush(CONS_OUT)
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
              if(op_code.eq.TAVP_INSTR_TENS_CREATE) then
               ds_instr%num_out_oprnds=1; ds_instr%out_oprnds(0:ds_instr%num_out_oprnds-1)=(/0/) !the single tensor operand of TENSOR_CREATE is an output operand
              elseif(op_code.eq.TAVP_INSTR_TENS_DESTROY) then
               ds_instr%num_out_oprnds=0 !TENSOR_DESTROY does not have output operands
              endif
             end select
            else
             call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE); jerr=-2
            endif
           else
            call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE); jerr=-1
           endif
           if(DEBUG.gt.0.and.jerr.ne.0) then
            write(CONS_OUT,'("#ERROR(TAVP-MNG:Decoder.decode.decode_instr_tens_create_destroy)[",i6,":",i3,"]: Error ",i11)')&
            &impir,omp_get_thread_num(),jerr
            flush(CONS_OUT)
           endif
           return
          end subroutine decode_instr_tens_create_destroy

          subroutine decode_instr_tens_transform(jerr)
           !TRANSFORMS a tensor, either from a defined or undefined state
           integer(INTD), intent(out):: jerr
           class(ds_instr_ctrl_t), pointer:: instr_ctrl
           class(ctrl_tens_trans_t), pointer:: tens_trans_ctrl

           allocate(tens_trans_ctrl,STAT=jerr) !tensor transformation control will be owned by the tensor instruction
           if(jerr.eq.0) then
            call tens_trans_ctrl%unpack(instr_packet,jerr)
            if(jerr.eq.0) then
             instr_ctrl=>tens_trans_ctrl; call ds_instr%set_control(instr_ctrl,jerr) !tensor transformation control ownership is moved to the tensor instruction
             if(jerr.eq.DSVP_SUCCESS) then
              tens_trans_ctrl=>NULL()
              call ds_instr%alloc_operands(1,jerr)
              if(jerr.eq.DSVP_SUCCESS) then
               call decode_instr_operands(jerr)
               if(jerr.eq.0) then !mark output operands
                select type(ds_instr)
                class is(tens_instr_t)
                 ds_instr%num_out_oprnds=1; ds_instr%out_oprnds(0:ds_instr%num_out_oprnds-1)=(/0/) !tensor operand 0 is the output operand
                end select
               else
                call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_BTC_BAD); jerr=-5
               endif
              else
               call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE); jerr=-4
              endif
             else
              deallocate(tens_trans_ctrl)
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE); jerr=-3
             endif
            else
             deallocate(tens_trans_ctrl)
             call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_BTC_BAD); jerr=-2
            endif
           else
            call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE); jerr=-1
           endif
           if(jerr.ne.0.and.VERBOSE) then
            write(CONS_OUT,'("#ERROR(TAVP-MNG:Decoder.decode.decode_instr_tens_transform)[",i6,":",i3,"]: Error ",i11)')&
            &impir,omp_get_thread_num(),jerr
            flush(CONS_OUT)
           endif
           return
          end subroutine decode_instr_tens_transform

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
                 ds_instr%num_out_oprnds=1; ds_instr%out_oprnds(0:ds_instr%num_out_oprnds-1)=(/0/) !tensor operand 0 is the output operand
                end select
               else
                call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_BTC_BAD); jerr=-5
               endif
              else
               call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE); jerr=-4
              endif
             else
              deallocate(tens_contr_ctrl)
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE); jerr=-3
             endif
            else
             deallocate(tens_contr_ctrl)
             call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_BTC_BAD); jerr=-2
            endif
           else
            call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE); jerr=-1
           endif
           if(jerr.ne.0.and.VERBOSE) then
            write(CONS_OUT,'("#ERROR(TAVP-MNG:Decoder.decode.decode_instr_tens_contract)[",i6,":",i3,"]: Error ",i11)')&
            &impir,omp_get_thread_num(),jerr
            flush(CONS_OUT)
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
         class(tavp_mng_retirer_t), intent(inout):: this !inout: TAVP-MNG Retirer DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier,thid,opcode,sts,iec,n,i,uid
         logical:: active,stopping,evicted
         class(dsvp_t), pointer:: dsvp
         class(tavp_mng_t), pointer:: tavp
         class(tens_rcrsv_t), pointer:: tensor
         class(tens_instr_t), pointer:: tens_instr
         type(tens_entry_mng_ref_t):: cache_entries(1:MAX_TENSOR_OPERANDS)
         type(obj_pack_t):: instr_packet
         type(comm_handle_t):: comm_hl
         class(*), pointer:: uptr

         errc=0; thid=omp_get_thread_num(); uid=this%get_id()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Retirer started as DSVU # ",i2," (thread ",i2,"): Reporting to ",i11,1x,i6)')&
          &impir,uid,thid,this%retire_comm,this%retire_rank
          flush(CONS_OUT)
         endif
!Reserve a bytecode buffer:
         call this%bytecode%reserve_mem(ier,MAX_BYTECODE_SIZE,MAX_BYTECODE_INSTR); if(ier.ne.0.and.errc.eq.0) errc=-30
!Initialize queues and ports:
         call this%init_queue(this%num_ports,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-29
!Set up tensor argument cache and wait on other TAVP units:
         tavp=>NULL(); dsvp=>this%get_dsvp(); select type(dsvp); class is(tavp_mng_t); tavp=>dsvp; end select
         if(associated(tavp)) then
          this%arg_cache=>tavp%tens_cache
          call tavp%sync_units(errc,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-28
         else
          this%arg_cache=>NULL(); if(errc.eq.0) errc=-27
         endif
!Work loop:
         active=((errc.eq.0).and.(this%retire_comm.ne.MPI_COMM_NULL)); stopping=(.not.active)
         wloop: do while(active)
          active=.not.stopping
 !Get a new batch of retired instructions (and control instructions) from Collector (port 0):
          ier=this%iqueue%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-26; exit wloop; endif
          ier=this%flush_port(0,num_moved=n); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-25; exit wloop; endif
          if(DEBUG.gt.0.and.n.gt.0) then
           write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Retirer unit ",i2," received ",i9," instructions from Collector")')&
           &impir,uid,n
           !ier=this%iqueue%reset(); ier=this%iqueue%scanp(action_f=tens_instr_print) !print all instructions
           flush(CONS_OUT)
          endif
 !Encode the retired instructions into bytecode and send them to the upper level:
          ier=this%iqueue%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-24; exit wloop; endif
          ier=this%iqueue%get_status(); if(stopping.and.ier.ne.GFC_IT_EMPTY.and.errc.eq.0) then; errc=-23; exit wloop; endif
          rloop: do while(ier.eq.GFC_IT_ACTIVE)
  !Extract an instruction:
           uptr=>this%iqueue%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-22; exit wloop; endif
           tens_instr=>NULL(); select type(uptr); class is(tens_instr_t); tens_instr=>uptr; end select
           if((.not.associated(tens_instr)).and.errc.eq.0) then; errc=-21; exit wloop; endif !trap
           sts=tens_instr%get_status(ier,iec); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-20; exit wloop; endif
           if(sts.ne.DS_INSTR_COMPLETED.and.errc.eq.0) then; errc=-19; exit wloop; endif !trap
           call tens_instr%set_status(DS_INSTR_RETIRED,ier,iec)
           if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-18; exit wloop; endif
           opcode=tens_instr%get_code(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-17; exit wloop; endif
           if(DEBUG.gt.0) then
            write(CONS_OUT,'("#DEBUG(TAVP-MNG:Retirer): Retired instruction:")')
            call tens_instr%print_it(dev_id=CONS_OUT)
            flush(CONS_OUT)
           endif
  !Encode a retired tensor instruction into bytecode:
           n=0 !number of tensor cache entries that may need eviction
           if(opcode.ge.TAVP_ISA_TENS_FIRST.and.opcode.le.TAVP_ISA_TENS_LAST) then
            if(opcode.eq.TAVP_INSTR_TENS_DESTROY) then
             call tens_instr%remove_persistency(ier); if(ier.ne.0.and.errc.eq.0) then; errc=-16; exit wloop; endif
            endif
            call tens_instr%extract_cache_entries(cache_entries,n,ier) !eviction candidates
            if(ier.ne.0.and.errc.eq.0) then; errc=-15; exit wloop; endif
            call this%bytecode%acquire_packet(instr_packet,ier,preclean=.TRUE.)
            if(ier.ne.PACK_SUCCESS.and.errc.eq.0) then; errc=-14; exit wloop; endif
            call this%encode(tens_instr,instr_packet,ier); if(ier.ne.0.and.errc.eq.0) then; errc=-13; exit wloop; endif
            call this%bytecode%seal_packet(ier); if(ier.ne.PACK_SUCCESS.and.errc.eq.0) then; errc=-12; exit wloop; endif
  !Check for control instructions:
           elseif(opcode.ge.TAVP_ISA_CTRL_FIRST.and.opcode.le.TAVP_ISA_CTRL_LAST) then
            if(opcode.eq.TAVP_INSTR_CTRL_STOP) stopping=.TRUE.
  !Other instructions are not expected here:
           else
            if(errc.eq.0) then; errc=-11; exit wloop; endif
           endif
  !Delete the instruction from the main queue:
           ier=this%iqueue%delete(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-10; exit wloop; endif
  !Evict unneded remote tensor cache entries:
           if(n.gt.0) then
            do i=1,n
             tensor=>cache_entries(i)%cache_entry%get_tensor(ier)
             if(ier.eq.0) then
              evicted=this%arg_cache%evict(tensor,ier,decr_use=.TRUE.); if(ier.ne.0.and.errc.eq.0) errc=-9
             else
              if(errc.eq.0) errc=-8
             endif
             cache_entries(i)%cache_entry=>NULL()
            enddo
            if(errc.ne.0) exit wloop
           endif
  !Send the bytecode to the parent TAVP at the upper level:
           n=this%bytecode%get_num_packets()
           if(n.gt.0.and.(n.ge.MAX_RETIRE_INSTR.or.this%iqueue%get_status().eq.GFC_IT_EMPTY)) then
            call comm_hl%clean(ier); if(ier.ne.PACK_SUCCESS.and.errc.eq.0) then; errc=-7; exit wloop; endif
            call this%bytecode%send(this%retire_rank,comm_hl,ier,tag=TAVP_COLLECT_TAG,comm=this%retire_comm)
            if(ier.ne.PACK_SUCCESS.and.errc.eq.0) then; errc=-6; exit wloop; endif
  !Synchronize the bytecode send:
            call comm_hl%wait(ier); if(ier.ne.PACK_SUCCESS.and.errc.eq.0) then; errc=-5; exit wloop; endif
            if(DEBUG.gt.0) then
             write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Retirer unit ",i2," retired ",i9," instructions")') impir,uid,n
             flush(CONS_OUT)
            endif
            call comm_hl%clean(ier); if(ier.ne.PACK_SUCCESS.and.errc.eq.0) then; errc=-4; exit wloop; endif
            call this%bytecode%clean(ier); if(ier.ne.PACK_SUCCESS.and.errc.eq.0) then; errc=-3; exit wloop; endif
           endif
           ier=this%iqueue%get_status()
          enddo rloop
          if(ier.ne.GFC_IT_EMPTY.and.errc.eq.0) then; errc=-2; exit wloop; endif
         enddo wloop
!Record the error:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(TAVP-MNG)[",i6,"]: Retirer error ",i11," by thread ",i2)')&
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
         class(tavp_mng_retirer_t), intent(inout):: this !inout: TAVP-MNG Retirer DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier,thid,uid

         errc=0; thid=omp_get_thread_num(); uid=this%get_id()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Retirer stopped as DSVU # ",i2," (thread ",i2,")")') impir,uid,thid
          flush(CONS_OUT)
         endif
!Release the tensor argument cache pointer:
         this%arg_cache=>NULL()
!Release queues:
         call this%release_queue(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-2
!Release the bytecode buffer:
         call this%bytecode%destroy(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
!Record an error, if any:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGRetirerShutdown
!-----------------------------------------------------------------------
        subroutine TAVPMNGRetirerEncode(this,ds_instr,instr_packet,ierr)
!Encodes a tensor instruction into a plain bytecode.
         implicit none
         class(tavp_mng_retirer_t), intent(inout):: this  !inout: TAVP-MNG Retirer DSVU
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
         class(tavp_mng_locator_t), intent(inout):: this !inout: TAVP-MNG Locator DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier,thid,opcode,rot_num,num_loc_instr,num_def_instr,i,n,num_sent,num_recv,sts,term_num,uid
         integer(INTL):: bytecode_tag,num_c_entries
         integer:: loc_wait
         logical:: active,stopping,ring_exists,located,inp_located,inp_valued,evicted,stalled,last
         type(tens_entry_mng_ref_t):: cache_entries(1:MAX_TENSOR_OPERANDS)
         class(tens_rcrsv_t), pointer:: tensor
         class(*), pointer:: uptr
         class(dsvp_t), pointer:: dsvp
         class(tavp_mng_t), pointer:: tavp
         class(tens_instr_t), pointer:: tens_instr
         type(tens_instr_t):: tens_instr_dummy
         type(obj_pack_t):: instr_packet
         type(comm_handle_t):: comm_hl

         errc=0; thid=omp_get_thread_num(); uid=this%get_id()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Locator started as DSVU # ",i2," (thread ",i2,"): Recv/Send = ",i6,"/",i6)')&
          &impir,uid,thid,this%ring_recv,this%ring_send
          flush(CONS_OUT)
         endif
!Reserve a bytecode buffer:
         call this%bytecode%reserve_mem(ier,MAX_BYTECODE_SIZE,MAX_BYTECODE_INSTR); if(ier.ne.0.and.errc.eq.0) errc=-93
!Initialize queues and ports:
         call this%init_queue(this%num_ports,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-92
!Initialize the located instruction list iterator:
         ier=this%loc_list%init(this%located_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-91
!Initialize the deferred instruction list iterator:
         ier=this%def_list%init(this%deferred_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-90
!Initialize the control list:
         ier=this%ctrl_list%init(this%control_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-89
!Set up a dummy instruction (RESUME) to carry the TAVP id info:
         call tens_instr_dummy%tens_instr_ctor(TAVP_INSTR_CTRL_RESUME,ier,iid=int(role_rank,INTL),stat=DS_INSTR_SPECIAL) !special dummy instruction with id = role_rank (TAVP-MNG id) and status = DS_INSTR_SPECIAL
         if(ier.ne.0.and.errc.eq.0) errc=-88
         if(tens_instr_locator_mark(tens_instr_dummy).ne.GFC_TRUE.and.errc.eq.0) errc=-87 !trap
!Check ring topology:
         if(this%ring_send.eq.role_rank.and.this%ring_recv.eq.role_rank) then !root TAVP-MNG (no ring)
          ring_exists=.FALSE.
         else !non-root TAVP-MNG (non-trivial ring)
          ring_exists=.TRUE.
          if((this%ring_send.eq.role_rank.or.this%ring_send.lt.0.or.&
            &this%ring_recv.eq.role_rank.or.this%ring_recv.lt.0).and.errc.eq.0) errc=-86 !trap
         endif
!Set up tensor argument cache and wait on other TAVP units:
         tavp=>NULL(); dsvp=>this%get_dsvp(); select type(dsvp); class is(tavp_mng_t); tavp=>dsvp; end select
         if(associated(tavp)) then
          this%arg_cache=>tavp%tens_cache
          call tavp%sync_units(errc,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-85
         else
          this%arg_cache=>NULL(); if(errc.eq.0) errc=-84
         endif
!Work loop:
         num_loc_instr=0; num_def_instr=0
         ier=timer_start(loc_wait,MAX_LOCATE_WAIT_TIME); if(ier.ne.TIMERS_SUCCESS.and.errc.eq.0) errc=-83
         active=((errc.eq.0).and.(this%ring_comm.ne.MPI_COMM_NULL)); stopping=(.not.active); stalled=.FALSE.
         wloop: do while(active)
 !Append new instructions from uDecoder (port 0) at the end of the main queue:
          ier=this%iqueue%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-82; exit wloop; endif
          ier=this%flush_port(0,num_moved=i); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-81; exit wloop; endif
          if(DEBUG.gt.0.and.i.gt.0) then
           write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Locator unit ",i2," appended ",i9,'//&
           &'" new instructions from uDecoder into main queue")') impir,uid,i
           !ier=this%iqueue%reset(); ier=this%iqueue%scanp(action_f=tens_instr_print); ier=this%iqueue%reset() !print all instructions
           flush(CONS_OUT)
          endif
 !Move the leading subset of subinstructions from port 2 (from Decomposer) into the locating list (bottom TAVP-MNG only):
          ier=this%loc_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-80; exit wloop; endif
          ier=this%unload_port(2,this%loc_list,max_items=MAX_LOCATE_SUB_INSTR,num_moved=n)
          if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-79; exit wloop; endif
          if(DEBUG.gt.0.and.n.gt.0) then
           write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Locator unit ",i2," received ",i9," instructions back from Decomposer")')&
           &impir,uid,n
           flush(CONS_OUT)
          endif
          num_loc_instr=num_loc_instr+n; n=0
 !Move the leading subset of new tensor instructions from the main queue into the locating list:
          ier=this%loc_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-78; exit wloop; endif
          ier=this%ctrl_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-77; exit wloop; endif
          if(this%ctrl_list%get_status().ne.GFC_IT_EMPTY.and.errc.eq.0) then; errc=-76; exit wloop; endif !trap
          ier=this%iqueue%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-75; exit wloop; endif
          if(this%iqueue%get_status().eq.GFC_IT_ACTIVE) then !main queue is not empty
           n=0
           mloop: do i=1,MAX_LOCATE_NEW_INSTR
  !Extract an instruction:
            uptr=>this%iqueue%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-74; exit wloop; endif
            tens_instr=>NULL(); select type(uptr); class is(tens_instr_t); tens_instr=>uptr; end select
            if((.not.associated(tens_instr)).and.errc.eq.0) then; errc=-73; exit wloop; endif
            opcode=tens_instr%get_code(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-72; exit wloop; endif
            if(DEBUG.gt.0) then
             write(CONS_OUT,'("#DEBUG(TAVP-MNG:Locator): An instruction is received for operand location:")')
             call tens_instr%print_it(dev_id=CONS_OUT)
             flush(CONS_OUT)
            endif
  !Move the instruction into an appropriate list (locating or control):
            if(opcode.ge.TAVP_ISA_TENS_FIRST.and.opcode.le.TAVP_ISA_TENS_LAST) then !tensor instruction
             if(stalled) exit mloop
             call tens_instr%set_status(DS_INSTR_INPUT_WAIT,ier)
             if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-71; exit wloop; endif
             ier=this%iqueue%move_elem(this%loc_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-70; exit wloop; endif
             n=n+1; num_loc_instr=num_loc_instr+1 !increment the number of instructions waiting for location
            elseif(opcode.ge.TAVP_ISA_CTRL_FIRST.and.opcode.le.TAVP_ISA_CTRL_LAST) then !control instruction
             if(stalled) exit mloop
             call tens_instr%set_status(DS_INSTR_READY_TO_EXEC,ier)
             if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-69; exit wloop; endif
             if(opcode.eq.TAVP_INSTR_CTRL_STOP) then
              if(DEBUG.gt.0) then
               write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Locator unit ",i2," received STOP instruction!")') impir,uid
               flush(CONS_OUT)
              endif
              stopping=.TRUE.
              call tens_instr_dummy%set_status(DS_INSTR_TERMINAL,ier)
              if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-68; exit wloop; endif
             endif
             ier=this%iqueue%move_elem(this%ctrl_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-67; exit wloop; endif
             n=n+1; stalled=.TRUE.; exit mloop
            else !auxiliary instruction
             call tens_instr%set_status(DS_INSTR_READY_TO_EXEC,ier)
             if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-66; exit wloop; endif
             ier=this%iqueue%move_elem(this%ctrl_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-65; exit wloop; endif
             n=n+1; stalled=.TRUE.
            endif
            if(this%iqueue%get_status().ne.GFC_IT_ACTIVE) exit mloop
           enddo mloop
           if(DEBUG.gt.0.and.n.gt.0) then
            write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Locator unit ",i2," extracted ",i9," new instructions from main queue")')&
            &impir,uid,n
            flush(CONS_OUT)
           endif
           call tavp%incr_recv_instr_counter(ier,int(n,INTL))
           if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-64; exit wloop; endif
          endif !main queue is not empty
 !Move a limited number of tensor instructions from the deferred list back into the locating list:`Deferred instructions are appended at the end (out-of-order)
          ier=this%loc_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-63; exit wloop; endif
          ier=this%def_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-62; exit wloop; endif
          if(this%def_list%get_status().eq.GFC_IT_ACTIVE.and.num_def_instr.lt.MAX_LOCATE_DEF_INSTR) then
           ier=this%def_list%move_list(this%loc_list,MAX_LOCATE_DEF_INSTR,n) !at most MAX_LOCATE_DEF_INSTR tensor instructions will be moved to the locating list from the deferred list
           if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-61; exit wloop; endif
           if(DEBUG.gt.0.and.n.gt.0) then
            write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Locator unit ",i2," extracted ",i9," instructions from deferred queue")')&
            &impir,uid,n
            flush(CONS_OUT)
           endif
           num_def_instr=num_def_instr+n !increment the number of instructions taken from the deferred instruction list
           num_loc_instr=num_loc_instr+n !increment the number of instructions waiting for location
          endif
 !Check the locating list length and whether the timer has expired:
          if(num_loc_instr.le.0) then !no tensor instructions for location: Dump control/auxiliary instructions to Decomposer
           ier=this%ctrl_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-60; exit wloop; endif
           if(this%ctrl_list%get_status().eq.GFC_IT_ACTIVE) then
            ier=tavp%decomposer%load_port(0,this%ctrl_list)
            if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-59; exit wloop; endif
           endif
           stalled=.FALSE. !since all control and auxiliary instructions have been dumped further into the pipeline
          endif
          if(timer_expired(loc_wait,ier,destroy=.TRUE.)) then
           ier=timer_start(loc_wait,MAX_LOCATE_WAIT_TIME); if(ier.ne.TIMERS_SUCCESS.and.errc.eq.0) then; errc=-58; exit wloop; endif
          else !timer has not expired yet
           if(num_loc_instr.lt.MIN_LOCATE_INSTR.and.(.not.stalled)) cycle wloop !no enough tensor instructions to initiate locating rotations
          endif
 !Perform a full rotation cycle for the tensor instructions being located at a given NAT level:
          term_num=0 !number of terminated Locators
          rot_num=0 !rotation number
          if(ring_exists) then
           if(DEBUG.gt.1) then
            num_c_entries=this%arg_cache%get_num_entries()
            write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Locator unit ",i2,": Entered rotation cycle; Cache volume = ",i12)')&
            &impir,uid,num_c_entries
            flush(CONS_OUT)
           endif
  !Insert the trailing dummmy instruction (CTRL_RESUME), which will tag the bytecode, into the locating list:
           ier=this%loc_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-57; exit wloop; endif
           ier=this%loc_list%append(tens_instr_dummy); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-56; exit wloop; endif
  !Begin rotations:
           rloop: do !rotations
  !Encode tensor instructions from the locating list into bytecode:
            ier=this%loc_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-55; exit wloop; endif
            if(this%loc_list%get_status().ne.GFC_IT_ACTIVE.and.errc.eq.0) then; errc=-54; exit wloop; endif !trap
            eloop: do while(ier.eq.GFC_SUCCESS)
             uptr=>this%loc_list%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-53; exit wloop; endif
             tens_instr=>NULL(); select type(uptr); class is(tens_instr_t); tens_instr=>uptr; end select
             if(.not.associated(tens_instr).and.errc.eq.0) then; errc=-52; exit wloop; endif !trap
             call this%bytecode%acquire_packet(instr_packet,ier); if(ier.ne.0.and.errc.eq.0) then; errc=-51; exit wloop; endif
             call this%encode(tens_instr,instr_packet,ier); if(ier.ne.0.and.errc.eq.0) then; errc=-50; exit wloop; endif
             call this%bytecode%seal_packet(ier); if(ier.ne.0.and.errc.eq.0) then; errc=-49; exit wloop; endif
             ier=this%loc_list%next()
            enddo eloop
            if(ier.ne.GFC_NO_MOVE.and.errc.eq.0) then; errc=-48; exit wloop; endif
  !Send the bytecode to the next (cousin) NAT node at the same tree level (in a ring):
            call comm_hl%clean(ier); if(ier.ne.0.and.errc.eq.0) then; errc=-47; exit wloop; endif
            call this%bytecode%send(this%ring_send,comm_hl,ier,tag=TAVP_LOCATE_TAG,comm=this%ring_comm)
            if(ier.ne.0.and.errc.eq.0) then; errc=-46; exit wloop; endif
            num_sent=this%bytecode%get_num_packets()
            if(DEBUG.gt.1.and.num_sent.gt.0) then !one instruction is always the special markup instruction
             write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Locator unit ",i2,": Rotation ",i3,": Sent ",i9," instructions")')&
             &impir,uid,rot_num+1,num_sent-1
             flush(CONS_OUT)
            endif
  !Clean the locating list and evict the relevant remote tensors with zero reference count from the tensor cache:
            ier=this%loc_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-45; exit wloop; endif
            ier=this%loc_list%get_status()
            dloop: do while(ier.eq.GFC_IT_ACTIVE)
   !Delete tensor instruction:
             uptr=>this%loc_list%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-44; exit wloop; endif
             tens_instr=>NULL(); select type(uptr); class is(tens_instr_t); tens_instr=>uptr; end select
             if(.not.associated(tens_instr).and.errc.eq.0) then; errc=-43; exit wloop; endif !trap
             call tens_instr%extract_cache_entries(cache_entries,n,ier) !eviction candidates
             if(ier.ne.0.and.errc.eq.0) then; errc=-42; exit wloop; endif
             call tens_instr%set_status(DS_INSTR_RETIRED,ier,DSVP_SUCCESS)
             if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-41; exit wloop; endif
             ier=this%loc_list%delete(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-40; exit wloop; endif
   !Evict remote tensor cache entries with zero reference and use counts:
             do i=1,n
              tensor=>cache_entries(i)%cache_entry%get_tensor(ier)
              if(ier.eq.0) then
               if(DEBUG.gt.1) then
                write(CONS_OUT,'("#DEBUG(TAVP-MNG:Locator): Trying to evict the following tensor:")')
                call tensor%print_it(dev_id=CONS_OUT); flush(CONS_OUT)
               endif
               evicted=this%arg_cache%evict(tensor,ier,decr_use=.TRUE.); if(ier.ne.0.and.errc.eq.0) errc=-39
               if(DEBUG.gt.1.and.ier.eq.0) then
                write(CONS_OUT,'("Eviction status: ",l1)') evicted
                flush(CONS_OUT)
               endif
              else
               if(errc.eq.0) errc=-38
              endif
              cache_entries(i)%cache_entry=>NULL()
             enddo
             if(errc.ne.0) exit wloop
             ier=this%loc_list%get_status()
            enddo dloop
            if(ier.ne.GFC_IT_EMPTY.and.errc.eq.0) then; errc=-37; exit wloop; endif
  !Synchronize the bytecode send and clean the bytecode buffer:
            call comm_hl%wait(ier); if(ier.ne.0.and.errc.eq.0) then; errc=-36; exit wloop; endif
            call comm_hl%clean(ier); if(ier.ne.0.and.errc.eq.0) then; errc=-35; exit wloop; endif
            call this%bytecode%clean(ier); if(ier.ne.0.and.errc.eq.0) then; errc=-34; exit wloop; endif
            if(DEBUG.gt.1.and.num_sent.gt.1) then
             write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Locator unit ",i2,": Rotation ",i3,": Send synced")') impir,uid,rot_num+1
             flush(CONS_OUT)
            endif
  !Absorb located tensor insructions from port 1 (delivered by lDecoder):
   !Wait until port 1 is filled by lDecoder:
            ier=this%loc_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-33; exit wloop; endif
            ier=this%loc_list%get_status(); if(ier.ne.GFC_IT_EMPTY.and.errc.eq.0) then; errc=-32; exit wloop; endif !trap
            do while(ier.eq.GFC_IT_EMPTY) !wait until located tensor instructions are delivered by lDecoder
             ier=this%unload_port(1,this%loc_list,stop_predicate=tens_instr_locator_mark,num_moved=num_recv)
             if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-31; exit wloop; endif
             if(DEBUG.gt.1.and.num_recv.gt.0) then
              write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Locator unit ",i2,": Rotation ",i3,": Received ",i6," instructions")')&
              &impir,uid,rot_num+1,num_recv-1
              flush(CONS_OUT)
             endif
             ier=this%loc_list%get_status()
            enddo
   !Read the last (dummy) CTRL_RESUME instruction with the bytecode tag written in its id (status = DS_INSTR_SPECIAL or DS_INSTR_TERMINAL):
            ier=this%loc_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-30; exit wloop; endif
            if(this%loc_list%get_status().eq.GFC_IT_ACTIVE) then
             uptr=>this%loc_list%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-29; exit wloop; endif
             tens_instr=>NULL(); select type(uptr); class is(tens_instr_t); tens_instr=>uptr; end select
             if(.not.associated(tens_instr).and.errc.eq.0) then; errc=-28; exit wloop; endif !trap
             opcode=tens_instr%get_code(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-27; exit wloop; endif
             sts=tens_instr%get_status(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-26; exit wloop; endif
             if(opcode.eq.TAVP_INSTR_CTRL_RESUME.and.(sts.eq.DS_INSTR_SPECIAL.or.sts.eq.DS_INSTR_TERMINAL)) then !dummy instruction is always (fake) CTRL_RESUME with status DS_INSTR_SPECIAL or DS_INSTR_TERMINAL
              if(sts.eq.DS_INSTR_TERMINAL) term_num=term_num+1 !one more terminated Locator detected
              rot_num=rot_num+1 !next rotation step (send/recv) has happened
              bytecode_tag=tens_instr%get_id(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-25; exit wloop; endif
              if(DEBUG.gt.1.and.num_recv.gt.0) then
               write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Locator unit ",i2,": Rotation ",i3,": Sender TAVP-MNG id = ",i4)')&
               &impir,uid,rot_num,bytecode_tag
               flush(CONS_OUT)
              endif
              if(bytecode_tag.eq.role_rank) then !TAVP's own instructions are back => remove the header dummy instruction
               call tens_instr%set_status(DS_INSTR_RETIRED,ier,DSVP_SUCCESS)
               if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-24; exit wloop; endif
               ier=this%loc_list%delete(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-23; exit wloop; endif
               exit rloop !rotations are finished
              endif
             else
              if(errc.eq.0) then; errc=-22; exit wloop; endif
             endif
            else
             if(errc.eq.0) then; errc=-21; exit wloop; endif
            endif
           enddo rloop
           if(DEBUG.gt.1) then
            num_c_entries=this%arg_cache%get_num_entries()
            write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Locator unit ",i2,": ",i4,"-rotation cycle completed; Cache volume = ",i12)')&
            &impir,uid,rot_num,num_c_entries
            flush(CONS_OUT)
           endif
          else !ring does not exist (single TAVP-MNG level)
           if(DEBUG.gt.1) then
            num_c_entries=this%arg_cache%get_num_entries()
            write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Locator unit ",i2,": No rotations; Cache volume = ",i12)')&
            &impir,uid,num_c_entries
            flush(CONS_OUT)
           endif
          endif
 !Move partially located tensor instructions from the locating list to the deferred list:
          ier=this%def_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-20; exit wloop; endif
          ier=this%loc_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-19; exit wloop; endif
          do while(this%loc_list%get_status().eq.GFC_IT_ACTIVE)
           uptr=>this%loc_list%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-18; exit wloop; endif
           tens_instr=>NULL(); select type(uptr); class is(tens_instr_t); tens_instr=>uptr; end select
           if(.not.associated(tens_instr).and.errc.eq.0) then; errc=-17; exit wloop; endif !trap
           located=tens_instr%fully_located(ier,inp_located,inp_valued)
           if(ier.ne.0.and.errc.eq.0) then
            if(DEBUG.gt.0) then
             write(CONS_OUT,'("#ERROR(TAVP-MNG:Locator)[",i6,":",i3,"]: tens_instr.fully_located() error ",i11)')&
             &impir,omp_get_thread_num(),ier
             flush(CONS_OUT)
            endif
            errc=-16; exit wloop
           endif
           if(.not.(inp_located.and.inp_valued)) then !input tensors must have been located and they must be defined
            if(DEBUG.gt.0) then
             write(CONS_OUT,'("#DEBUG(TAVP-MNG:Locator): An instruction is deferred (operands not ready):")')
             call tens_instr%print_it(dev_id=CONS_OUT)
             flush(CONS_OUT)
            endif
            last=this%loc_list%on_last()
            ier=this%loc_list%move_elem(this%def_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-15; exit wloop; endif
            if(last) ier=this%loc_list%next() !to make iterator DONE
           else !instruction is ready to be executed (issued to the next level)
            call tens_instr%set_status(DS_INSTR_READY_TO_EXEC,ier)
            if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-14; exit wloop; endif
            ier=this%loc_list%next()
           endif
          enddo
          ier=this%loc_list%get_status()
          if((ier.ne.GFC_IT_EMPTY.and.ier.ne.GFC_IT_DONE).and.errc.eq.0) then; errc=-13; exit wloop; endif
 !Print located and deferred tensor instructions (debug):
          if(DEBUG.gt.1.and.errc.eq.0) then
           ier=this%loc_list%reset()
           if(this%loc_list%get_status().eq.GFC_IT_ACTIVE) then
            write(CONS_OUT,'("#DEBUG(TAVP-MNG:Locator): Queue of located tensor instructions:")')
            ier=this%loc_list%scanp(action_f=tens_instr_print); ier=this%loc_list%reset() !print all instructions
           endif
           ier=this%def_list%reset()
           if(this%def_list%get_status().eq.GFC_IT_ACTIVE) then
            write(CONS_OUT,'("#DEBUG(TAVP-MNG:Locator): Queue of deferred tensor instructions:")')
            ier=this%def_list%scanp(action_f=tens_instr_print); ier=this%def_list%reset() !print all instructions
           endif
           flush(CONS_OUT)
          endif
 !Append saved control instructions at the end of the locating list:
          ier=this%def_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-12; exit wloop; endif
          if(this%def_list%get_status().eq.GFC_IT_EMPTY) then !cannot proceed until all located tensor instructions have been issued
           ier=this%loc_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-11; exit wloop; endif
           ier=this%ctrl_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-10; exit wloop; endif
           if(this%ctrl_list%get_status().eq.GFC_IT_ACTIVE) then
            ier=this%ctrl_list%move_list(this%loc_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-9; exit wloop; endif
           endif
           stalled=.FALSE.
          endif
 !Move located tensor instructions plus auxiliary/control instructions from the locating list to Decomposer:
          ier=this%loc_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-8; exit wloop; endif
          if(this%loc_list%get_status().eq.GFC_IT_ACTIVE) then
           if(DEBUG.gt.0) then
            !ier=this%loc_list%scanp(action_f=tens_instr_print); ier=this%loc_list%reset() !print instructions
           endif
           ier=tavp%decomposer%load_port(0,this%loc_list); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-7; exit wloop; endif
          endif
 !Reset locating instruction counters:
          num_loc_instr=0; num_def_instr=0
          if(stopping.and.(.not.stalled).and.(term_num.eq.rot_num)) active=.FALSE.
         enddo wloop
 !Send STOP instruction to lDecoder:
         call tens_instr_dummy%set_code(TAVP_INSTR_CTRL_STOP,ier)
         if(ier.eq.DSVP_SUCCESS) then
          call tens_instr_dummy%set_status(DS_INSTR_NEW,ier,DSVP_SUCCESS)
          if(ier.eq.DSVP_SUCCESS) then
           ier=this%ctrl_list%append(tens_instr_dummy)
           if(ier.eq.GFC_SUCCESS) then
            ier=tavp%ldecoder%load_port(0,this%ctrl_list) !send STOP instruction to lDecoder
            if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-6
           else
            if(errc.eq.0) errc=-5
           endif
          else
           if(errc.eq.0) errc=-4
          endif
         else
          if(errc.eq.0) errc=-3
         endif
!Retire the special dummy instruction:
         call tens_instr_dummy%set_status(DS_INSTR_RETIRED,ier,DSVP_SUCCESS); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-2
!Record the error:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(TAVP-MNG)[",i6,"]: Locator error ",i11," by thread ",i2)')&
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
         class(tavp_mng_locator_t), intent(inout):: this !inout: TAVP-MNG Locator DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier,thid,uid

         errc=0; thid=omp_get_thread_num(); uid=this%get_id()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Locator stopped as DSVU # ",i2," (thread ",i2,")")') impir,uid,thid
          flush(CONS_OUT)
         endif
!Release the tensor argument cache pointer:
         this%arg_cache=>NULL()
!Deactivate the control list:
         ier=this%ctrl_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%ctrl_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-11
           ier=this%ctrl_list%delete_all()
          endif
          ier=this%ctrl_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-10
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
         class(tavp_mng_locator_t), intent(inout):: this  !inout: TAVP-MNG Locator DSVU
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
         class(tavp_mng_decomposer_t), intent(inout):: this !inout: TAVP-MNG Decomposer DSVU
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc,ier,thid,iec,sts,opcode,num_processed,base_created,i,n,uid
         integer:: dec_timer
         logical:: active,stopping,expired,last_instr
         class(dsvp_t), pointer:: dsvp
         class(tavp_mng_t), pointer:: tavp
         class(tens_rcrsv_t), pointer:: tensor
         class(tens_instr_t), pointer:: tens_instr
         class(*), pointer:: uptr

         errc=0; thid=omp_get_thread_num(); uid=this%get_id()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer started as DSVU # ",i2," (thread ",i2,")")') impir,uid,thid
          flush(CONS_OUT)
         endif
!Initialize queues:
         call this%init_queue(this%num_ports,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-59
!Initialize the subinstruction list iterator:
         ier=this%sub_list%init(this%subinstruction_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-58
!Initialize the returning list iterator:
         ier=this%ret_list%init(this%returning_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-57
!Initialize the collecting list iterator:
         ier=this%col_list%init(this%collecting_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-56
!Initialize the auxiliary list iterator:
         ier=this%aux_list%init(this%auxiliary_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-55
!Initialize the control list iterator:
         ier=this%ctrl_list%init(this%control_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-54
!Set up tensor argument cache and wait on other TAVP units:
         tavp=>NULL(); dsvp=>this%get_dsvp(); select type(dsvp); class is(tavp_mng_t); tavp=>dsvp; end select
         if(associated(tavp)) then
          this%arg_cache=>tavp%tens_cache
          call tavp%sync_units(errc,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-53
         else
          this%arg_cache=>NULL(); if(errc.eq.0) errc=-52
         endif
!Check whether this TAVP-MNG is at the bottom of TAVP-MNG hierarchy:
         this%tavp_is_top=tavp%is_top(ier); if(ier.ne.0.and.errc.eq.0) errc=-51
         this%tavp_is_bottom=tavp%is_bottom(ier); if(ier.ne.0.and.errc.eq.0) errc=-50
!Acquire a timer:
         ier=timer_start(dec_timer,MAX_DECOMPOSE_PHASE_TIME); if(ier.ne.TIMERS_SUCCESS.and.errc.eq.0) errc=-49
!Work loop:
         active=(errc.eq.0); stopping=(.not.active)
         wloop: do while(active)
 !Get a batch of instructions from Locator into the main queue (new + previously decomposed possibly):
          ier=this%iqueue%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-48; exit wloop; endif
          ier=this%flush_port(0,num_moved=i); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-47; exit wloop; endif
          if(DEBUG.gt.0.and.i.gt.0) then
           write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer unit ",i2," received ",i9," instructions from Locator")')&
           &impir,uid,i
           !ier=this%iqueue%reset(); ier=this%iqueue%scanp(action_f=tens_instr_print) !print all instructions
           flush(CONS_OUT)
          endif
 !Decompose parent tensor instructions, creating new child subinstructions (but filter previously decomposed subinstructions):
          num_processed=0; base_created=tavp%get_crtd_instr_counter()
          ier=this%iqueue%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-46; exit wloop; endif
          ier=this%iqueue%get_status()
          dloop: do while(ier.eq.GFC_IT_ACTIVE)
  !Extract an instruction:
           last_instr=this%iqueue%on_last(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-45; exit wloop; endif
           uptr=>this%iqueue%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-44; exit wloop; endif
           tens_instr=>NULL(); select type(uptr); class is(tens_instr_t); tens_instr=>uptr; end select
           if((.not.associated(tens_instr)).and.errc.eq.0) then; errc=-43; exit wloop; endif !trap
  !Separate new parent instructions from old child subinstructions that just passed the location phase (bottom TAVP-MNG only):
           sts=tens_instr%get_status(ier,iec); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-42; exit wloop; endif
           if(sts.ne.DS_INSTR_READY_TO_EXEC.and.errc.eq.0) then; errc=-41; exit wloop; endif !trap
           if(iec.eq.TAVP_ERR_TAG_ONE) then !special tag to distinguish previously decomposed child tensor subinstructions from new parent instructions
            if(DEBUG.gt.0) then
             write(CONS_OUT,'("#DEBUG(TAVP-MNG:Decomposer): A previously decomposed instruction is received back (bottom):")')
             call tens_instr%print_it(dev_id=CONS_OUT)
             flush(CONS_OUT)
            endif
            call tens_instr%set_status(sts,ier,DSVP_SUCCESS) !remove the tag and reset the error code back to success
            if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-40; exit wloop; endif
            ier=this%iqueue%move_elem(this%ret_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-39; exit wloop; endif
           else !new instruction
  !Process the extracted parent instruction according to its category:
            opcode=tens_instr%get_code(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-38; exit wloop; endif
   !Decompose the parent tensor instruction into child tensor instructions and append them into the subinstruction list:
            if(opcode.ge.TAVP_ISA_TENS_FIRST.and.opcode.le.TAVP_ISA_TENS_LAST) then
             call this%decompose(tens_instr,ier); if(ier.ne.0.and.errc.eq.0) then; errc=-37; exit wloop; endif
             if(DEBUG.gt.0) then
              write(CONS_OUT,'("#DEBUG(TAVP-MNG:Decomposer): A new instruction was decomposed:")')
              call tens_instr%print_it(dev_id=CONS_OUT)
              flush(CONS_OUT)
             endif
   !Move processed parent instruction to the collecting list:
             ier=this%iqueue%move_elem(this%col_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-36; exit wloop; endif
             num_processed=num_processed+1
   !Move control instructions into the control list:
            elseif(opcode.ge.TAVP_ISA_CTRL_FIRST.and.opcode.le.TAVP_ISA_CTRL_LAST) then
             if(opcode.eq.TAVP_INSTR_CTRL_STOP) stopping=.TRUE.
             ier=this%iqueue%move_elem(this%ctrl_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-35; exit wloop; endif
   !Move other instructions into the auxiliary list:
            else
             ier=this%iqueue%move_elem(this%aux_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-34; exit wloop; endif
            endif
           endif
  !Pass a batch of new child subinstructions to Dispatcher (port 0) or Locator (port 2) and a batch of processed parent instructions to Collector (port 0):
           expired=timer_expired(dec_timer,ier); if(ier.ne.TIMERS_SUCCESS.and.errc.eq.0) then; errc=-33; exit wloop; endif
           if(expired.or.last_instr.or.tavp%get_crtd_instr_counter()-base_created.ge.MAX_DECOMPOSE_CHLD_INSTR.or.&
             &num_processed.ge.MAX_DECOMPOSE_PRNT_INSTR.or.stopping) then
            if(DEBUG.gt.0) then
             write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer unit ",i2," entered distribution phase")')&
             &impir,uid
             flush(CONS_OUT)
            endif
   !Clone control instructions for Collector:
            ier=this%col_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-32; exit wloop; endif
            ier=this%ctrl_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-31; exit wloop; endif
            if(this%ctrl_list%get_status(ier).eq.GFC_IT_ACTIVE) then
             n=0
             do while(ier.eq.GFC_SUCCESS)
              uptr=>this%ctrl_list%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-30; exit wloop; endif
              tens_instr=>NULL(); select type(uptr); class is(tens_instr_t); tens_instr=>uptr; end select
              if((.not.associated(tens_instr)).and.errc.eq.0) then; errc=-29; exit wloop; endif !trap
              call tens_instr%set_status(DS_INSTR_ISSUED,ier,DSVP_SUCCESS) !all instructions passed to Collector are marked DS_INSTR_ISSUED
              if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-28; exit wloop; endif
              ier=this%col_list%append(tens_instr); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-27; exit wloop; endif
              ier=this%col_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-26; exit wloop; endif
              call tens_instr%set_status(DS_INSTR_READY_TO_EXEC,ier,DSVP_SUCCESS) !revert instruction status back to the original state
              if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-25; exit wloop; endif
              n=n+1; ier=this%ctrl_list%next()
             enddo
             if(ier.ne.GFC_NO_MOVE.and.errc.eq.0) then; errc=-24; exit wloop; endif
             if(DEBUG.gt.0.and.n.gt.0) then
              write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer unit ",i2," copied ",i9," control instructions for Collector")')&
              &impir,uid,n
              flush(CONS_OUT)
             endif
            endif
   !Pass processed parent instructions (together with cloned control instructions) to Collector port 0:
            ier=this%col_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-23; exit wloop; endif
            if(this%col_list%get_status().eq.GFC_IT_ACTIVE) then
             ier=tavp%collector%load_port(0,this%col_list,num_moved=n)
             if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-22; exit wloop; endif
             if(DEBUG.gt.0.and.n.gt.0) then
              write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer unit ",i2," passed ",i9," parent instructions to Collector")')&
              &impir,uid,n
              flush(CONS_OUT)
             endif
            endif
            if(this%tavp_is_bottom) then !bottom TAVP-MNG: Previous tensor subinstructions --> Dispatcher port 0; New tensor subinstructions --> Locator port 2
   !Pass previously decomposed child subinstructions to Dispatcher port 0 and new child subinstructions to Locator port 2 (bottom TAVP-MNG only):
    !Append auxiliary instructions to the end of the returning list:
             ier=this%ret_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-21; exit wloop; endif
             ier=this%aux_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-20; exit wloop; endif
             ier=this%aux_list%move_list(this%ret_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-19; exit wloop; endif
    !Append control instructions to the end of the returning list:
             ier=this%ret_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-18; exit wloop; endif
             ier=this%ctrl_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-17; exit wloop; endif
             ier=this%ctrl_list%move_list(this%ret_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-16; exit wloop; endif
    !Pass returning tensor subinstructions, auxiliary instructions and control instructions to Dispatcher:
             ier=this%ret_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-15; exit wloop; endif
             if(this%ret_list%get_status().eq.GFC_IT_ACTIVE) then
              ier=tavp%dispatcher%load_port(0,this%ret_list,num_moved=n)
              if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-14; exit wloop; endif
              if(DEBUG.gt.0.and.n.gt.0) then
               write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer unit ",i2," passed ",i9," subinstructions to Dispatcher")')&
               &impir,uid,n
               flush(CONS_OUT)
              endif
             endif
    !Pass new tensor subinstructions to Locator port 2 (for metadata locating):
             ier=this%sub_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-13; exit wloop; endif
             if(this%sub_list%get_status().eq.GFC_IT_ACTIVE) then
              ier=tavp%locator%load_port(2,this%sub_list,num_moved=n)
              if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-12; exit wloop; endif
              if(DEBUG.gt.0.and.n.gt.0) then
               write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer unit ",i2," passed ",i9," subinstructions back to Locator")')&
               &impir,uid,n
               flush(CONS_OUT)
              endif
             endif
            else !not bottom TAVP-MNG: New tensor subinstructions --> Dispatcher port 0
   !Pass new child subinstructions to Dispatcher port 0:
    !Append auxiliary instructions to the end of the subinstruction list:
             ier=this%sub_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-11; exit wloop; endif
             ier=this%aux_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-10; exit wloop; endif
             ier=this%aux_list%move_list(this%sub_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-9; exit wloop; endif
    !Append control instructions to the end of the subinstruction list:
             ier=this%sub_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-8; exit wloop; endif
             ier=this%ctrl_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-7; exit wloop; endif
             ier=this%ctrl_list%move_list(this%sub_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-6; exit wloop; endif
    !Pass tensor subinstructions, auxiliary instructions and control instructions to Dispatcher:
             ier=this%sub_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-5; exit wloop; endif
             if(this%sub_list%get_status().eq.GFC_IT_ACTIVE) then
              ier=tavp%dispatcher%load_port(0,this%sub_list,num_moved=n)
              if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-4; exit wloop; endif
              if(DEBUG.gt.0.and.n.gt.0) then
               write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer unit ",i2," passed ",i9," subinstructions to Dispatcher")')&
               &impir,uid,n
               flush(CONS_OUT)
              endif
             endif
            endif
  !Clear local counters:
            num_processed=0; base_created=tavp%get_crtd_instr_counter()
  !Restart the timer:
            if(expired) then
             ier=timer_reset(dec_timer,MAX_DECOMPOSE_PHASE_TIME)
             if(ier.ne.TIMERS_SUCCESS.and.errc.eq.0) then; errc=-3; exit wloop; endif
            endif
            if(DEBUG.gt.0) then
             write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer unit ",i2," exited distribution phase")') impir,uid
             flush(CONS_OUT)
            endif
           endif !distribution phase
           ier=this%iqueue%get_status()
          enddo dloop
          if(ier.ne.GFC_IT_EMPTY.and.errc.eq.0) then; errc=-2; exit wloop; endif
          active=.not.stopping
         enddo wloop
!Destroy the timer:
         ier=timer_destroy(dec_timer)
!Record the error:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(TAVP-MNG)[",i6,"]: Decomposer error ",i11," by thread ",i2)')&
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
         class(tavp_mng_decomposer_t), intent(inout):: this !inout: TAVP-MNG Decomposer DSVU
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc,ier,thid,uid

         errc=0; thid=omp_get_thread_num(); uid=this%get_id()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer stopped as DSVU # ",i2," (thread ",i2,")")')&
          &impir,uid,thid
          flush(CONS_OUT)
         endif
!Release the tensor argument cache pointer:
         this%arg_cache=>NULL()
!Deactivate the control list:
         ier=this%ctrl_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%ctrl_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-16
           ier=this%ctrl_list%delete_all()
          endif
          ier=this%ctrl_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-15
         else
          if(errc.eq.0) errc=-14
         endif
!Deactivate the auxiliary list:
         ier=this%aux_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%aux_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-13
           ier=this%aux_list%delete_all()
          endif
          ier=this%aux_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-12
         else
          if(errc.eq.0) errc=-11
         endif
!Deactivate the collecting list:
         ier=this%col_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%col_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-10
           ier=this%col_list%delete_all()
          endif
          ier=this%col_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-9
         else
          if(errc.eq.0) errc=-8
         endif
!Deactivate the returning list:
         ier=this%ret_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%ret_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-7
           ier=this%ret_list%delete_all()
          endif
          ier=this%ret_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-6
         else
          if(errc.eq.0) errc=-5
         endif
!Deactivate the subinstruction list:
         ier=this%sub_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%sub_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-4
           ier=this%sub_list%delete_all()
          endif
          ier=this%sub_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-3
         else
          if(errc.eq.0) errc=-2
         endif
!Release queues:
         call this%release_queue(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
!Record an error, if any:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDecomposerShutdown
!------------------------------------------------------------------
        subroutine TAVPMNGDecomposerDecompose(this,tens_instr,ierr)
!Decomposes a parental tensor instruction <tens_instr> into child instructions (subinstructions)
!and appends them into the subinstruction list. The input tensor instruction operands
!must already contain their subtensor composition lists. The output tensor instruction
!operand(s) may either already contain or still lack their subtensor composition lists.
!In the latter case, they will be generated here. All subinstructions generated from
!the parental tensor instruction are registered with the TAVP-MNG instruction map.
!The parental tensor instruction saves the number of its subinstructions in the
!.error_code field which will be reset back to DSVP_SUCCESS by Collector.
         implicit none
         class(tavp_mng_decomposer_t), intent(inout):: this !inout: TAVP-MNG Decomposer DSVU
         class(tens_instr_t), intent(inout):: tens_instr    !in: parental tensor instruction (.error_code field will be set to the subinstruction count)
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc,opcode,init_stat,init_tag,uid
         integer(INTL):: parent_id
         class(dsvp_t), pointer:: dsvp
         class(tavp_mng_t), pointer:: tavp

         uid=this%get_id()
         if(tens_instr%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           dsvp=>this%get_dsvp(errc)
           if(errc.eq.DSVP_SUCCESS.and.associated(dsvp)) then
            tavp=>NULL(); select type(dsvp); class is(tavp_mng_t); tavp=>dsvp; end select
            if(tavp%is_bottom()) then
             init_stat=DS_INSTR_INPUT_WAIT; init_tag=TAVP_ERR_TAG_ONE
            else
             init_stat=DS_INSTR_READY_TO_EXEC; init_tag=DSVP_SUCCESS
            endif
!Decompose structureless output tensor operands (set their composition lists), if needed:
            call decompose_output_tensors(errc)
!Decompose the tensor instruction:
            if(errc.eq.0) then
             parent_id=tens_instr%get_id(errc)
             if(errc.eq.DSVP_SUCCESS) then
              opcode=tens_instr%get_code(errc)
              if(errc.eq.DSVP_SUCCESS) then
               select case(opcode)
 !TENSOR CREATE:
               case(TAVP_INSTR_TENS_CREATE) !new subinstructions will go into the subinstruction list
                call decompose_instr_tens_create(errc)
 !TENSOR DESTROY:
               case(TAVP_INSTR_TENS_DESTROY) !new subinstructions will go into the subinstruction list
                call decompose_instr_tens_destroy(errc)
 !TENSOR INIT:
               case(TAVP_INSTR_TENS_INIT) !new subinstructions will go into the subinstruction list
                call decompose_instr_tens_transform(errc)
 !TENSOR CONTRACT:
               case(TAVP_INSTR_TENS_CONTRACT) !new subinstructions will go into the subinstruction list
                call decompose_instr_tens_contract(errc)
               case default
                !`Call decomposition for other relevant tensor instructions
                write(CONS_OUT,&
                &'("#FATAL(TAVP-MNG:Decomposer:decompose)[",i6,"]: Tensor instruction code ",i3," is not implemented!")')&
                &impir,opcode
                errc=-7
               end select
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

         subroutine decompose_output_tensors(jerr)
         !In case the output tensor(s) do not have internal structure yet,
         !this subroutine will decompose them into subtensors, based on
         !the universal tensor dimension strength assessing function.
         !The newly created subtensors will then be cached in the tensor cache.
         !`Different TAVPs may create the same output subtensors, causing duplication.
          implicit none
          integer(INTD), intent(out):: jerr !out: error code
          integer(INTD):: jj,ns,split_dims(1:MAX_TENSOR_RANK)
          real(8):: total_strength,dim_strength(1:MAX_TENSOR_RANK)
          class(ds_oprnd_t), pointer:: oprnd
          class(tens_rcrsv_t), pointer:: tensor

          jerr=0
          tloop: do jj=0,tens_instr%num_out_oprnds-1 !loop over the output tensor operands
           oprnd=>tens_instr%get_operand(tens_instr%out_oprnds(jj),jerr)
           if(jerr.ne.DSVP_SUCCESS) then; jerr=-8; exit tloop; endif
           select type(oprnd)
           class is(tens_oprnd_t)
            call oprnd%lock()
            if(this%tavp_is_top) then !root TAVP-MNG owns all incoming output tensors
             if(oprnd%get_owner_id(jerr).lt.0) then
              if(jerr.eq.0) then
               call oprnd%set_owner_id(role_rank,jerr,update_cache=.TRUE.); if(jerr.ne.0) jerr=-7
              else
               jerr=-6
              endif
             else
              if(jerr.ne.0) jerr=-5
             endif
            endif
            if(jerr.eq.0) then
             tensor=>oprnd%get_tensor(jerr)
             if(jerr.eq.0) then
              if(tensor%get_num_subtensors().le.0) then !tensor does not have an internal structure yet
               total_strength=tens_dim_strength_assess(tensor,dim_strength,jerr,tens_dim_strength_thresh,ns,split_dims)
               if(jerr.eq.0) then
                call tensor%decompose(split_dims(1:ns),jerr); if(jerr.ne.TEREC_SUCCESS) jerr=-4
               else
                jerr=-3
               endif
              endif
             else
              jerr=-2
             endif
             tensor=>NULL()
            endif
            call oprnd%unlock()
           class default
            jerr=-1
           end select
           if(jerr.ne.0) exit tloop
          enddo tloop
          return
         end subroutine decompose_output_tensors

         subroutine decompose_instr_tens_create(jerr)
         !Decomposes TENSOR_CREATE instruction into subinstructions,
         !subsequently appending them into the subinstruction list.
          implicit none
          integer(INTD), intent(out):: jerr !out: error code
          integer(INTD):: jj,num_subinstr,data_type
          integer(INTL):: iid
          class(ds_oprnd_t), pointer:: oprnd,tens_oprnd
          class(tens_header_t), pointer:: header
          class(tens_rcrsv_t), pointer:: tensor,subtensor
          class(tens_cache_entry_t), pointer:: tens_entry
          class(tens_entry_mng_t), pointer:: tens_entry_mng
          class(tens_instr_t), pointer:: subinstr
          type(tens_instr_t):: tens_instr_empty
          type(list_bi_t), pointer:: subtensors
          type(list_iter_t):: lit
          class(*), pointer:: uptr
          logical:: stored

          jerr=this%sub_list%reset_back()
          if(jerr.eq.GFC_SUCCESS) then
           oprnd=>tens_instr%get_operand(0,jerr)
           if(jerr.eq.DSVP_SUCCESS) then
            select type(oprnd)
            class is(tens_oprnd_t)
             call oprnd%lock()
             call oprnd%reset_persistency(.TRUE.) !mark the parental tensor cache entry as persistent in the tensor cache
             tensor=>oprnd%get_tensor(jerr) !parental tensor
             if(jerr.eq.0) then
              !if(DEBUG.gt.0) then
               !write(CONS_OUT,'("#DEBUG(TAVP-MNG:Decomposer.decompose)[",i6,"]: Decomposing TENS_CREATE for tensor:")') impir
               !call tensor%print_it(dev_id=CONS_OUT)
               !flush(CONS_OUT)
              !endif
              data_type=tensor%get_data_type() !parental tensor data type, if any, will propagate into subtensors
              if(tensor%get_num_subtensors().gt.0) then !tensor must have an internal structure
               subtensors=>tensor%get_subtensors(jerr) !list of subtensors in terms of tensor headers
               if(jerr.eq.TEREC_SUCCESS) then
                jerr=lit%init(subtensors)
                if(jerr.eq.GFC_SUCCESS) then
                 num_subinstr=0; tens_entry=>NULL()
                 cloop: do while(jerr.eq.GFC_SUCCESS)
                  uptr=>lit%get_value(jerr); if(jerr.ne.GFC_SUCCESS) then; jerr=-26; exit cloop; endif
                  header=>NULL(); select type(uptr); class is(tens_header_t); header=>uptr; end select
                  if(.not.associated(header)) then; jerr=-25; exit cloop; endif !trap
                  allocate(tens_rcrsv_t::subtensor,STAT=jerr); if(jerr.ne.0) then; jerr=-24; exit cloop; endif
                  call subtensor%tens_rcrsv_ctor(header,jerr)
                  if(jerr.ne.TEREC_SUCCESS) then; deallocate(subtensor); jerr=-23; exit cloop; endif
                  if(data_type.ne.EXA_DATA_KIND_NN) then
                   call subtensor%set_data_type(data_type,jerr)
                   if(jerr.ne.TEREC_SUCCESS) then; deallocate(subtensor); jerr=-22; exit cloop; endif
                  endif
                  !if(DEBUG.gt.0) then
                   !write(CONS_OUT,'("#DEBUG(TAVP-MNG:Decomposer.decompose)[",i6,"]: TENS_CREATE created a subtensor:")') impir
                   !call subtensor%print_it(dev_id=CONS_OUT)
                   !flush(CONS_OUT)
                  !endif
                  stored=this%arg_cache%store(subtensor,tens_entry_mng_alloc,jerr,tens_entry_p=tens_entry) !new subtensor: Ownership transferred to the tensor cache
                  if(.not.(jerr.eq.0.and.stored.and.associated(tens_entry))) then
                   deallocate(subtensor); jerr=-21; exit cloop
                  endif
                  call tens_entry%lock()
                  call tens_entry%set_persistency(.TRUE.) !mark subtensor as persistent in the tensor cache
                  tens_entry_mng=>NULL()
                  select type(tens_entry); class is(tens_entry_mng_t); tens_entry_mng=>tens_entry; end select
                  jerr=this%sub_list%append(tens_instr_empty); if(jerr.ne.GFC_SUCCESS) then; jerr=-20; exit cloop; endif
                  jerr=this%sub_list%reset_back(); if(jerr.ne.GFC_SUCCESS) then; jerr=-19; exit cloop; endif
                  uptr=>this%sub_list%get_value(jerr); if(jerr.ne.GFC_SUCCESS) then; jerr=-18; exit cloop; endif
                  subinstr=>NULL(); select type(uptr); class is(tens_instr_t); subinstr=>uptr; end select
                  if(.not.associated(subinstr)) then; jerr=-17; exit cloop; endif !trap
                  iid=dsvp%get_crtd_instr_counter()
                  call subinstr%tens_instr_ctor(TAVP_INSTR_TENS_CREATE,jerr,subtensor,iid)
                  if(jerr.ne.0) then; jerr=-16; exit cloop; endif
                  tens_oprnd=>subinstr%get_operand(0,jerr); if(jerr.ne.DSVP_SUCCESS) then; jerr=-15; exit cloop; endif
                  select type(tens_oprnd)
                  class is(tens_oprnd_t)
                   call tens_oprnd%set_cache_entry(tens_entry_mng,jerr); if(jerr.ne.0) then; jerr=-14; exit cloop; endif
                  class default
                   jerr=-13; exit cloop
                  end select
                  call subinstr%set_status(init_stat,jerr,init_tag) !instruction status and error code will depend on whether the TAVP-MNG is bottom or not
                  if(jerr.ne.DSVP_SUCCESS) then; jerr=-12; exit cloop; endif
                  call tavp%register_instr(iid,tens_instr%get_id(),jerr); if(jerr.ne.0) then; jerr=-11; exit cloop; endif
                  call tens_entry%unlock(); call this%arg_cache%release_entry(tens_entry); tens_entry=>NULL()
                  call dsvp%incr_crtd_instr_counter(); num_subinstr=num_subinstr+1 !new subinstruction has been created
                  jerr=lit%next()
                 enddo cloop
                 if(associated(tens_entry)) then !in case of error
                  call tens_entry%unlock(); call this%arg_cache%release_entry(tens_entry); tens_entry=>NULL()
                 endif
                 if(jerr.eq.GFC_NO_MOVE) then
                  call tens_instr%set_status(DS_INSTR_ISSUED,jerr,num_subinstr) !.error_code of the parent instruction will store the number of subinstructions (in Collector)
                  if(jerr.ne.DSVP_SUCCESS) jerr=-10
                 else
                  jerr=-9
                 endif
                 jj=lit%release(); if(jj.ne.GFC_SUCCESS.and.jerr.eq.0) jerr=-8
                else
                 jerr=-7
                endif
               else
                jerr=-6
               endif
              else
               jerr=-5
              endif
             else
              jerr=-4
             endif
             call oprnd%unlock()
            class default
             jerr=-3
            end select
           else
            jerr=-2
           endif
          else
           jerr=-1
          endif
          if(DEBUG.gt.0) then
           if(jerr.eq.0) then
            write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer unit ",i2," created ",i9," TENS_CREATE subinstructions")')&
            &impir,uid,num_subinstr
           else
            write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer unit ",i2,'//&
            &'" failed to create TENS_CREATE subinstructions: Error ",i11)') impir,uid,jerr
           endif
           flush(CONS_OUT)
          endif
          return
         end subroutine decompose_instr_tens_create

         subroutine decompose_instr_tens_destroy(jerr)
         !Decomposes TENSOR_DESTROY instruction into subinstructions,
         !subsequently appending them into the subinstruction list.
          implicit none
          integer(INTD), intent(out):: jerr !out: error code
          integer(INTD):: jj,num_subinstr
          integer(INTL):: iid
          class(ds_oprnd_t), pointer:: oprnd,tens_oprnd
          class(tens_header_t), pointer:: header
          class(tens_rcrsv_t), pointer:: tensor,subtensor
          class(tens_cache_entry_t), pointer:: tens_entry
          class(tens_entry_mng_t), pointer:: tens_entry_mng
          class(tens_instr_t), pointer:: subinstr
          type(tens_instr_t):: tens_instr_empty
          type(list_bi_t), pointer:: subtensors
          type(list_iter_t):: lit
          class(*), pointer:: uptr

          jerr=this%sub_list%reset_back()
          if(jerr.eq.GFC_SUCCESS) then
           oprnd=>tens_instr%get_operand(0,jerr)
           if(jerr.eq.DSVP_SUCCESS) then
            select type(oprnd)
            class is(tens_oprnd_t)
             call oprnd%lock()
             tensor=>oprnd%get_tensor(jerr)
             if(jerr.eq.0) then
              if(tensor%get_num_subtensors().gt.0) then !tensor must have an internal structure
               subtensors=>tensor%get_subtensors(jerr) !list of subtensors in terms of tensor headers
               if(jerr.eq.TEREC_SUCCESS) then
                jerr=lit%init(subtensors)
                if(jerr.eq.GFC_SUCCESS) then
                 num_subinstr=0; tens_entry=>NULL()
                 cloop: do while(jerr.eq.GFC_SUCCESS)
                  uptr=>lit%get_value(jerr); if(jerr.ne.GFC_SUCCESS) then; jerr=-27; exit cloop; endif
                  header=>NULL(); select type(uptr); class is(tens_header_t); header=>uptr; end select
                  if(.not.associated(header)) then; jerr=-26; exit cloop; endif !trap
                  allocate(tens_rcrsv_t::subtensor,STAT=jerr); if(jerr.ne.0) then; jerr=-25; exit cloop; endif
                  call subtensor%tens_rcrsv_ctor(header,jerr)
                  if(jerr.ne.TEREC_SUCCESS) then; deallocate(subtensor); jerr=-24; exit cloop; endif
                  tens_entry=>this%arg_cache%lookup(subtensor,jerr) !subtensor must be present in the tensor cache since its creation
                  if((jerr.ne.0).or.(.not.associated(tens_entry))) then; deallocate(subtensor); jerr=-23; exit cloop; endif !trap
                  call tens_entry%lock()
                  deallocate(subtensor); subtensor=>tens_entry%get_tensor(jerr); if(jerr.ne.0) then; jerr=-22; exit cloop; endif
                  tens_entry_mng=>NULL()
                  select type(tens_entry); class is(tens_entry_mng_t); tens_entry_mng=>tens_entry; end select
                  if(.not.associated(tens_entry_mng)) then; jerr=-21; exit cloop; endif !trap
                  jerr=this%sub_list%append(tens_instr_empty); if(jerr.ne.GFC_SUCCESS) then; jerr=-20; exit cloop; endif
                  jerr=this%sub_list%reset_back(); if(jerr.ne.GFC_SUCCESS) then; jerr=-19; exit cloop; endif
                  uptr=>this%sub_list%get_value(jerr); if(jerr.ne.GFC_SUCCESS) then; jerr=-18; exit cloop; endif
                  subinstr=>NULL(); select type(uptr); class is(tens_instr_t); subinstr=>uptr; end select
                  if(.not.associated(subinstr)) then; jerr=-17; exit cloop; endif !trap
                  iid=dsvp%get_crtd_instr_counter()
                  call subinstr%tens_instr_ctor(TAVP_INSTR_TENS_DESTROY,jerr,subtensor,iid)
                  if(jerr.ne.0) then; jerr=-16; exit cloop; endif
                  tens_oprnd=>subinstr%get_operand(0,jerr); if(jerr.ne.DSVP_SUCCESS) then; jerr=-15; exit cloop; endif
                  select type(tens_oprnd)
                  class is(tens_oprnd_t)
                   call tens_oprnd%set_cache_entry(tens_entry_mng,jerr); if(jerr.ne.0) then; jerr=-14; exit cloop; endif
                  class default
                   jerr=-13; exit cloop
                  end select
                  call subinstr%set_status(init_stat,jerr,init_tag) !instruction status and error code will depend on whether the TAVP-MNG is bottom or not
                  if(jerr.ne.DSVP_SUCCESS) then; jerr=-12; exit cloop; endif
                  call tavp%register_instr(iid,tens_instr%get_id(),jerr); if(jerr.ne.0) then; jerr=-11; exit cloop; endif
                  call tens_entry%unlock(); call this%arg_cache%release_entry(tens_entry); tens_entry=>NULL()
                  call dsvp%incr_crtd_instr_counter(); num_subinstr=num_subinstr+1 !new subinstruction has been created
                  jerr=lit%next()
                 enddo cloop
                 if(associated(tens_entry)) then !in case of error
                  call tens_entry%unlock(); call this%arg_cache%release_entry(tens_entry); tens_entry=>NULL()
                 endif
                 if(jerr.eq.GFC_NO_MOVE) then
                  call tens_instr%set_status(DS_INSTR_ISSUED,jerr,num_subinstr) !.error_code of the parent instruction will store the number of subinstructions (in Collector)
                  if(jerr.ne.DSVP_SUCCESS) jerr=-10
                 else
                  jerr=-9
                 endif
                 jj=lit%release(); if(jj.ne.GFC_SUCCESS.and.jerr.eq.0) jerr=-8
                else
                 jerr=-7
                endif
               else
                jerr=-6
               endif
              else
               jerr=-5
              endif
             else
              jerr=-4
             endif
             call oprnd%unlock()
            class default
             jerr=-3
            end select
           else
            jerr=-2
           endif
          else
           jerr=-1
          endif
          if(DEBUG.gt.0) then
           if(jerr.eq.0) then
            write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer unit ",i2," created ",i9," TENS_DESTROY subinstructions")')&
            &impir,uid,num_subinstr
           else
            write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer unit ",i2,'//&
            &'" failed to create TENS_DESTROY subinstructions: Error ",i11)') impir,uid,jerr
           endif
           flush(CONS_OUT)
          endif
          return
         end subroutine decompose_instr_tens_destroy

         subroutine decompose_instr_tens_transform(jerr)
         !Decomposes TENSOR_INIT or TENSOR_UNARY instruction into subinstructions,
         !subsequently appending them into the subinstruction list.
          implicit none
          integer(INTD), intent(out):: jerr !out: error code
          integer(INTD):: num_subinstr,num_subtrans,jer,jj,jn
          integer(INTL):: iid
          type(tens_entry_mng_ref_t):: cache_entries(0:0) !tensor transformation has 1 argument
          class(tens_rcrsv_t), pointer:: subtensor
          class(tens_operation_t), allocatable:: tens_trans
          class(tens_transformation_t), pointer:: subtrans
          class(tens_cache_entry_t), pointer:: tens_entry
          class(ds_oprnd_t), pointer:: tens_oprnd
          class(tens_instr_t), pointer:: subinstr
          type(tens_instr_t):: tens_instr_empty
          type(list_bi_t):: subtransformations
          type(list_iter_t):: subit
          class(*), pointer:: uptr
          logical:: stored

          jerr=this%sub_list%reset_back()
          if(jerr.eq.GFC_SUCCESS) then
 !Extract the underlying tensor operation (tensor transformation):
           call tens_instr%get_operation(tens_trans,jerr)
           if(jerr.eq.0) then
            select type(tens_trans)
            class is(tens_transformation_t)
 !Decompose the tensor operation into suboperations:
             num_subinstr=0; num_subtrans=0
             call tens_trans%split(subtransformations,jerr,num_subtrans) !internal decomposition induced by the tensor structure
             if(jerr.eq.TEREC_SUCCESS) then
 !Generate subinstructions from the suboperations:
              jerr=subit%init(subtransformations)
              if(jerr.eq.GFC_SUCCESS) then
               tens_entry=>NULL(); jerr=subit%get_status()
               sloop: do while(jerr.eq.GFC_IT_ACTIVE)
                uptr=>subit%get_value(jerr); if(jerr.ne.GFC_SUCCESS) then; jerr=-27; exit sloop; endif
                subtrans=>NULL(); select type(uptr); class is(tens_transformation_t); subtrans=>uptr; end select
                if(.not.associated(subtrans)) then; jerr=-26; exit sloop; endif !trap
  !Store subtensors in the tensor cache if they are not there:
                tens_entry=>NULL()
                jn=subtrans%get_num_args(jerr); if(jerr.ne.TEREC_SUCCESS) then; jerr=-25; exit sloop; endif
                do jj=0,jn-1
                 subtensor=>subtrans%get_argument(jj,jerr); if(jerr.ne.TEREC_SUCCESS) then; jerr=-24; exit sloop; endif
                 stored=this%arg_cache%store(subtensor,tens_entry_mng_alloc,jerr,tens_entry_p=tens_entry) !<stored> assumes tensor ownership be moved
                 if(jerr.ne.0) then; jerr=-23; exit sloop; endif
                 call tens_entry%lock()
                 if(stored) then !release ownership of the subtensor since it has been moved to the tensor cache
                  call subtrans%donate_argument(jj,jerr); if(jerr.ne.TEREC_SUCCESS) then; jerr=-22; exit sloop; endif
                 else !subtensor already existed in the tensor cache: Use it
                  subtensor=>tens_entry%get_tensor(jerr); if(jerr.ne.0) then; jerr=-21; exit sloop; endif
                  call subtrans%reset_argument(subtensor,jj,jerr); if(jerr.ne.TEREC_SUCCESS) then; jerr=-20; exit sloop; endif
                 endif
                 select type(tens_entry)
                 class is(tens_entry_mng_t)
                  cache_entries(jj)%cache_entry=>tens_entry
                 class default
                  jerr=-19; exit sloop
                 end select
                 call tens_entry%unlock()
                 tens_entry=>NULL()
                enddo
  !Construct a new subinstruction:
                jerr=this%sub_list%append(tens_instr_empty); if(jerr.ne.GFC_SUCCESS) then; jerr=-18; exit sloop; endif
                jerr=this%sub_list%reset_back(); if(jerr.ne.GFC_SUCCESS) then; jerr=-17; exit sloop; endif
                uptr=>this%sub_list%get_value(jerr); if(jerr.ne.GFC_SUCCESS) then; jerr=-16; exit sloop; endif
                subinstr=>NULL(); select type(uptr); class is(tens_instr_t); subinstr=>uptr; end select
                if(.not.associated(subinstr)) then; jerr=-15; exit sloop; endif !trap
                iid=dsvp%get_crtd_instr_counter()
                call subinstr%tens_instr_ctor(TAVP_INSTR_TENS_INIT,jerr,subtrans,iid)
                if(jerr.ne.0) then; jerr=-14; exit sloop; endif
  !Back-associate tensor operands with their corresponding tensor cache entries:
                do jj=0,jn-1
                 tens_oprnd=>subinstr%get_operand(jj,jerr); if(jerr.ne.DSVP_SUCCESS) then; jerr=-13; exit sloop; endif
                 select type(tens_oprnd)
                 class is(tens_oprnd_t)
                  call tens_oprnd%set_cache_entry(cache_entries(jj)%cache_entry,jerr)
                  if(jerr.ne.0) then; jerr=-12; exit sloop; endif
                  tens_entry=>cache_entries(jj)%cache_entry
                  call this%arg_cache%release_entry(tens_entry); tens_entry=>NULL(); cache_entries(jj)%cache_entry=>NULL()
                 class default
                  jerr=-11; exit sloop
                 end select
                enddo
  !Register the new subinstruction:
                call subinstr%set_status(init_stat,jerr,init_tag) !instruction status and error code will depend on whether the TAVP-MNG is bottom or not
                if(jerr.ne.DSVP_SUCCESS) then; jerr=-10; exit sloop; endif
                call tavp%register_instr(iid,tens_instr%get_id(),jerr); if(jerr.ne.0) then; jerr=-9; exit sloop; endif
                call dsvp%incr_crtd_instr_counter(); num_subinstr=num_subinstr+1 !new subinstruction has been created
  !Delete the processed suboperation:
                jerr=subit%delete(); if(jerr.ne.GFC_SUCCESS) then; jerr=-8; exit sloop; endif
                jerr=subit%get_status()
               enddo sloop
               if(jerr.ne.0) then
                if(associated(tens_entry)) then !in case of error
                 call tens_entry%unlock(); call this%arg_cache%release_entry(tens_entry); tens_entry=>NULL()
                endif
                do jj=ubound(cache_entries,1),lbound(cache_entries,1),-1 !in case of error
                 if(associated(cache_entries(jj)%cache_entry)) then
                  tens_entry=>cache_entries(jj)%cache_entry
                  call this%arg_cache%release_entry(tens_entry); cache_entries(jj)%cache_entry=>NULL(); tens_entry=>NULL()
                 endif
                enddo
               endif
               if(jerr.eq.GFC_IT_EMPTY) jerr=subit%get_status(); if(jerr.eq.GFC_IT_EMPTY) jerr=0
               if(jerr.eq.0) then
                call tens_instr%set_status(DS_INSTR_ISSUED,jerr,num_subinstr) !.error_code of the parent instruction will store the number of subinstructions (in Collector)
                if(jerr.ne.DSVP_SUCCESS) jerr=-7
               endif
               jer=subit%release(); if(jer.ne.GFC_SUCCESS.and.jerr.eq.0) jerr=-6
              else
               jerr=-5
              endif
             else
              jerr=-4
             endif
            class default
             jerr=-3
            end select
           else
            jerr=-2
           endif
          else
           jerr=-1
          endif
          if(DEBUG.gt.0) then
           if(jerr.eq.0) then
            write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer unit ",i2," created ",i9," TENS_INIT subinstructions")')&
            &impir,uid,num_subinstr
           else
            write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer unit ",i2,'//&
            &'" failed to create TENS_INIT subinstructions with error ",i11," for parent tensor instruction:")')&
            &impir,uid,jerr
            call tens_instr%print_it(dev_id=CONS_OUT)
           endif
           flush(CONS_OUT)
          endif
          return
         end subroutine decompose_instr_tens_transform

         subroutine decompose_instr_tens_contract(jerr)
         !Decomposes TENSOR_CONTRACT instruction into subinstructions,
         !subsequently appending them into the subinstruction list.
          implicit none
          integer(INTD), intent(out):: jerr !out: error code
          integer(INTD):: jj,jn,num_subinstr,num_subcontr
          integer(INTL):: iid
          type(tens_entry_mng_ref_t):: cache_entries(0:2) !tensor contraction has 3 arguments
          class(tens_rcrsv_t), pointer:: subtensor
          class(tens_operation_t), allocatable:: tens_contr
          class(tens_contraction_t), pointer:: subcontr
          class(tens_cache_entry_t), pointer:: tens_entry
          class(tens_entry_mng_t), pointer:: tens_entry_mng
          class(ds_oprnd_t), pointer:: tens_oprnd
          class(tens_instr_t), pointer:: subinstr
          type(tens_instr_t):: tens_instr_empty
          type(list_bi_t):: subcontractions
          type(list_iter_t):: subit
          class(*), pointer:: uptr
          logical:: stored

          jerr=this%sub_list%reset_back()
          if(jerr.eq.GFC_SUCCESS) then
 !Extract the underlying tensor operation (tensor contraction):
           call tens_instr%get_operation(tens_contr,jerr)
           if(jerr.eq.0) then
            select type(tens_contr)
            class is(tens_contraction_t)
 !Decompose the tensor operation into suboperations:
             num_subinstr=0; num_subcontr=0
             call tens_contr%split(subcontractions,jerr,num_subcontr) !internal decomposition induced by the tensor structure
             if(jerr.eq.TEREC_SUCCESS) then
 !Generate subinstructions from the suboperations:
              jerr=subit%init(subcontractions)
              if(jerr.eq.GFC_SUCCESS) then
               tens_entry=>NULL(); jerr=subit%get_status()
               sloop: do while(jerr.eq.GFC_IT_ACTIVE)
                uptr=>subit%get_value(jerr); if(jerr.ne.GFC_SUCCESS) then; jerr=-27; exit sloop; endif
                subcontr=>NULL(); select type(uptr); class is(tens_contraction_t); subcontr=>uptr; end select
                if(.not.associated(subcontr)) then; jerr=-26; exit sloop; endif !trap
  !Store subtensors in the tensor cache if they are not there:
                tens_entry=>NULL()
                jn=subcontr%get_num_args(jerr); if(jerr.ne.TEREC_SUCCESS) then; jerr=-25; exit sloop; endif
                do jj=0,jn-1
                 subtensor=>subcontr%get_argument(jj,jerr); if(jerr.ne.TEREC_SUCCESS) then; jerr=-24; exit sloop; endif
                 stored=this%arg_cache%store(subtensor,tens_entry_mng_alloc,jerr,tens_entry_p=tens_entry) !<stored> assumes tensor ownership be moved to the tensor cache
                 if(jerr.ne.0) then; jerr=-23; exit sloop; endif
                 call tens_entry%lock()
                 if(stored) then !release ownership of the subtensor since it has been moved to the tensor cache
                  call subcontr%donate_argument(jj,jerr); if(jerr.ne.TEREC_SUCCESS) then; jerr=-22; exit sloop; endif
                 else !subtensor already existed in the tensor cache: Use it
                  subtensor=>tens_entry%get_tensor(jerr); if(jerr.ne.0) then; jerr=-21; exit sloop; endif
                  call subcontr%reset_argument(subtensor,jj,jerr); if(jerr.ne.TEREC_SUCCESS) then; jerr=-20; exit sloop; endif
                 endif
                 select type(tens_entry)
                 class is(tens_entry_mng_t)
                  cache_entries(jj)%cache_entry=>tens_entry
                 class default
                  jerr=-19; exit sloop
                 end select
                 call tens_entry%unlock()
                 tens_entry=>NULL()
                enddo
  !Construct a new subinstruction:
                jerr=this%sub_list%append(tens_instr_empty); if(jerr.ne.GFC_SUCCESS) then; jerr=-18; exit sloop; endif
                jerr=this%sub_list%reset_back(); if(jerr.ne.GFC_SUCCESS) then; jerr=-17; exit sloop; endif
                uptr=>this%sub_list%get_value(jerr); if(jerr.ne.GFC_SUCCESS) then; jerr=-16; exit sloop; endif
                subinstr=>NULL(); select type(uptr); class is(tens_instr_t); subinstr=>uptr; end select
                if(.not.associated(subinstr)) then; jerr=-15; exit sloop; endif !trap
                iid=dsvp%get_crtd_instr_counter()
                call subinstr%tens_instr_ctor(TAVP_INSTR_TENS_CONTRACT,jerr,subcontr,iid)
                if(jerr.ne.0) then; jerr=-14; exit sloop; endif
  !Back-associate tensor operands with their corresponding tensor cache entries:
                do jj=0,jn-1
                 tens_oprnd=>subinstr%get_operand(jj,jerr); if(jerr.ne.DSVP_SUCCESS) then; jerr=-13; exit sloop; endif
                 select type(tens_oprnd)
                 class is(tens_oprnd_t)
                  call tens_oprnd%set_cache_entry(cache_entries(jj)%cache_entry,jerr)
                  if(jerr.ne.0) then; jerr=-12; exit sloop; endif
                  tens_entry=>cache_entries(jj)%cache_entry
                  call this%arg_cache%release_entry(tens_entry); tens_entry=>NULL(); cache_entries(jj)%cache_entry=>NULL()
                 class default
                  jerr=-11; exit sloop
                 end select
                enddo
  !Register the new subinstruction:
                call subinstr%set_status(init_stat,jerr,init_tag) !instruction status and error code will depend on whether the TAVP-MNG is bottom or not
                if(jerr.ne.DSVP_SUCCESS) then; jerr=-10; exit sloop; endif
                call tavp%register_instr(iid,tens_instr%get_id(),jerr); if(jerr.ne.0) then; jerr=-9; exit sloop; endif
                call dsvp%incr_crtd_instr_counter(); num_subinstr=num_subinstr+1 !new subinstruction has been created
  !Delete the processed suboperation:
                jerr=subit%delete(); if(jerr.ne.GFC_SUCCESS) then; jerr=-8; exit sloop; endif
                jerr=subit%get_status()
               enddo sloop
               if(jerr.ne.0) then
                if(associated(tens_entry)) then !in case of error
                 call tens_entry%unlock(); call this%arg_cache%release_entry(tens_entry); tens_entry=>NULL()
                endif
                do jj=ubound(cache_entries,1),lbound(cache_entries,1),-1 !in case of error
                 if(associated(cache_entries(jj)%cache_entry)) then
                  tens_entry=>cache_entries(jj)%cache_entry
                  call this%arg_cache%release_entry(tens_entry); cache_entries(jj)%cache_entry=>NULL(); tens_entry=>NULL()
                 endif
                enddo
               endif
               if(jerr.eq.GFC_IT_EMPTY) jerr=subit%get_status(); if(jerr.eq.GFC_IT_EMPTY) jerr=0
               if(jerr.eq.0) then
                call tens_instr%set_status(DS_INSTR_ISSUED,jerr,num_subinstr) !.error_code of the parent instruction will store the number of subinstructions (in Collector)
                if(jerr.ne.DSVP_SUCCESS) jerr=-7
               endif
               jj=subit%release(); if(jj.ne.GFC_SUCCESS.and.jerr.eq.0) jerr=-6
              else
               jerr=-5
              endif
             else
              jerr=-4
             endif
            class default
             jerr=-3
            end select
           else
            jerr=-2
           endif
 !Delete the tensor operation:
           if(allocated(tens_contr)) deallocate(tens_contr)
          else
           jerr=-1
          endif
          if(DEBUG.gt.0) then
           if(jerr.eq.0) then
            write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer unit ",i2," created ",i9," TENS_CONTRACT subinstructions")')&
            &impir,uid,num_subinstr
           else
            write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Decomposer unit ",i2,'//&
            &'" failed to create TENS_CONTRACT subinstructions with error ",i11," for parent tensor instruction:")')&
            &impir,uid,jerr
            call tens_instr%print_it(dev_id=CONS_OUT)
           endif
           flush(CONS_OUT)
          endif
          return
         end subroutine decompose_instr_tens_contract

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
         class(tavp_mng_dispatcher_t), intent(inout):: this !inout: TAVP-MNG Dispatcher DSVU
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc,ier,thid,i,n,opcode,sts,iec,channel,uid
         logical:: active,stopping,synced
         class(dsvp_t), pointer:: dsvp
         class(tavp_mng_t), pointer:: tavp
         class(tens_instr_t), pointer:: tens_instr
         class(*), pointer:: uptr

         errc=0; thid=omp_get_thread_num(); uid=this%get_id()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Dispatcher started as DSVU # ",i2," (thread ",i2,'//&
          &'") with ",i4," channels over communicator ",i11)') impir,uid,thid,this%num_ranks,this%dispatch_comm
          flush(CONS_OUT)
         endif
!Reserve bytecode buffers and clean communication handles:
         allocate(this%bytecode(this%num_ranks),this%comm_hl(this%num_ranks),STAT=ier)
         if(ier.eq.0) then
          do i=1,this%num_ranks
           call this%bytecode(i)%reserve_mem(ier,MAX_BYTECODE_SIZE,MAX_BYTECODE_INSTR); if(ier.ne.0.and.errc.eq.0) errc=-29
          enddo
          do i=1,this%num_ranks
           call this%comm_hl(i)%clean(ier); if(ier.ne.0.and.errc.eq.0) errc=-28
          enddo
         else
          if(errc.eq.0) errc=-27
         endif
!Allocate dispatched instruction counters:
         allocate(this%dispatch_count(this%num_ranks),this%dispatch_flops(this%num_ranks),STAT=ier)
         if(ier.eq.0) then
          this%dispatch_count(:)=0; this%dispatch_flops(:)=0d0
         else
          if(errc.eq.0) errc=-26
         endif
!Reset the next channel for the round-robin dispatch:
         this%next_channel=1 !channels: [1:this%num_ranks]
!Initialize queues:
         call this%init_queue(this%num_ports,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-25
!Set up tensor argument cache and wait on other TAVP units:
         tavp=>NULL(); dsvp=>this%get_dsvp(); select type(dsvp); class is(tavp_mng_t); tavp=>dsvp; end select
         if(associated(tavp)) then
          this%arg_cache=>tavp%tens_cache
          call tavp%sync_units(errc,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-24
         else
          this%arg_cache=>NULL(); if(errc.eq.0) errc=-23
         endif
!Check whether this TAVP-MNG is at the bottom of the TAVP-MNG hierarchy:
         this%tavp_is_bottom=tavp%is_bottom(ier); if(ier.ne.0.and.errc.eq.0) errc=-22
!Work loop:
         active=((errc.eq.0).and.(this%dispatch_comm.ne.MPI_COMM_NULL)); stopping=(.not.active)
         wloop: do while(active)
 !Receive newly created child subinstructions from Decomposer into port 0:
          ier=this%iqueue%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-21; exit wloop; endif
          ier=this%flush_port(0,num_moved=i); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-20; exit wloop; endif
          if(DEBUG.gt.0.and.i.gt.0) then
           write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Dispatcher unit ",i2," received ",i9," instructions from Decomposer")')&
           &impir,uid,i
           !ier=this%iqueue%reset(); ier=this%iqueue%scanp(action_f=tens_instr_print) !print all instructions
           flush(CONS_OUT)
          endif
 !Dispatch/encode the instructions into the bytecode buffers and issue bytecode to the child TAVPs:
          ier=this%iqueue%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-19; exit wloop; endif
          ier=this%iqueue%get_status()
          dloop: do while(ier.eq.GFC_IT_ACTIVE)
  !Extract an instruction:
           uptr=>this%iqueue%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-18; exit wloop; endif
           tens_instr=>NULL(); select type(uptr); class is(tens_instr_t); tens_instr=>uptr; end select
           if((.not.associated(tens_instr)).and.errc.eq.0) then; errc=-17; exit wloop; endif !trap
           sts=tens_instr%get_status(ier,iec); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-16; exit wloop; endif
           if((sts.ne.DS_INSTR_READY_TO_EXEC.or.iec.ne.DSVP_SUCCESS).and.errc.eq.0) then; errc=-15; exit wloop; endif !trap
           call tens_instr%set_status(DS_INSTR_NEW,ier,iec) !reset the instruction status to NEW before dispatching to the lower level
           if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-14; exit wloop; endif
           opcode=tens_instr%get_code(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-13; exit wloop; endif
           if(DEBUG.gt.0) then
            write(CONS_OUT,'("#DEBUG(TAVP-MNG:Dispatcher): A new instruction is dispatched:")')
            call tens_instr%print_it(dev_id=CONS_OUT)
            flush(CONS_OUT)
           endif
           if(opcode.ge.TAVP_ISA_TENS_FIRST.and.opcode.le.TAVP_ISA_TENS_LAST) then
  !Encode a tensor instruction and dispatch it to the appropriate channel:
            channel=this%map_instr(tens_instr,ier); if(ier.ne.0.and.errc.eq.0) then; errc=-12; exit wloop; endif
            if((channel.lt.lbound(this%dispatch_rank,1).or.channel.gt.ubound(this%dispatch_rank,1)).and.errc.eq.0) then
             errc=-11; exit wloop !trap
            endif
            call this%dispatch(tens_instr,channel,ier); if(ier.ne.0.and.errc.eq.0) then; errc=-10; exit wloop; endif
           else
  !Encode an auxiliary/control instruction and dispatch it to all channels:
            do i=lbound(this%dispatch_rank,1),ubound(this%dispatch_rank,1)
             call this%dispatch(tens_instr,i,ier); if(ier.ne.0.and.errc.eq.0) then; errc=-9; exit wloop; endif
            enddo
  !Check on control instructions:
            if(opcode.eq.TAVP_INSTR_CTRL_STOP) then
             stopping=.TRUE.
            elseif(opcode.eq.TAVP_INSTR_CTRL_DUMP_CACHE) then
             write(jo,'("#DEBUG(TAVP-MNG): TENSOR CACHE DUMP:")')
             call this%arg_cache%print_it()
             flush(jo)
            endif
           endif
  !Delete the dispatched instruction from the main queue:
           call tens_instr%set_status(DS_INSTR_RETIRED,ier,DSVP_SUCCESS)
           if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-8; exit wloop; endif
           if(DEBUG.gt.1) then
            write(CONS_OUT,'("#MSG(TAVP-MNG:Dispatcher): Deleting dispatched instruction:")')
            call tens_instr%print_it(dev_id=CONS_OUT); flush(CONS_OUT)
           endif
           ier=this%iqueue%delete(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-7; exit wloop; endif
  !Issue (send) the bytecode to the child TAVPs:
           do i=1,this%num_ranks !loop over dispatch channels
            n=this%bytecode(i)%get_num_packets(ier); if(ier.ne.PACK_SUCCESS.and.errc.eq.0) then; errc=-6; exit wloop; endif
            if(n.gt.0) then
             if(n.ge.MAX_ISSUE_INSTR.or.this%dispatch_count(i).le.MIN_ISSUE_INSTR.or.this%iqueue%get_status().eq.GFC_IT_EMPTY) then
              call this%issue(i,ier); if(ier.ne.0.and.errc.eq.0) then; errc=-5; exit wloop; endif
             endif
            endif
           enddo
  !Synchronize the bytecode sends to each TAVP:
           synced=this%sync_issue(ier); if(ier.ne.0.and.errc.eq.0) then; errc=-4; exit wloop; endif
           ier=this%iqueue%get_status()
           if(stopping.and.ier.eq.GFC_IT_ACTIVE.and.errc.eq.0) then; errc=-3; exit wloop; endif !trap
           active=.not.stopping
          enddo dloop
          if(ier.ne.GFC_IT_EMPTY.and.errc.eq.0) then; errc=-2; exit wloop; endif
         enddo wloop
!Record the error:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(TAVP-MNG)[",i6,"]: Dispatcher error ",i11," by thread ",i2)')&
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
         class(tavp_mng_dispatcher_t), intent(inout):: this !inout: TAVP-MNG Dispatcher DSVU
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc,ier,thid,i,uid

         errc=0; thid=omp_get_thread_num(); uid=this%get_id()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Dispatcher stopped as DSVU # ",i2," (thread ",i2,")")') impir,uid,thid
          flush(CONS_OUT)
         endif
!Release the tensor argument cache pointer:
         this%arg_cache=>NULL()
!Release queues:
         call this%release_queue(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-3
!Deallocate dispatched instruction counters:
         if(allocated(this%dispatch_flops)) deallocate(this%dispatch_flops)
         if(allocated(this%dispatch_count)) deallocate(this%dispatch_count)
!Release bytecode buffers and communication handles:
         if(allocated(this%comm_hl)) then
          do i=ubound(this%comm_hl,1),lbound(this%comm_hl,1),-1
           call this%comm_hl(i)%clean(ier); if(ier.ne.0.and.errc.eq.0) errc=-2
          enddo
          deallocate(this%comm_hl)
         endif
         if(allocated(this%bytecode)) then
          do i=ubound(this%bytecode,1),lbound(this%bytecode,1),-1
           call this%bytecode(i)%destroy(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
          enddo
          deallocate(this%bytecode)
         endif
!Record an error, if any:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDispatcherShutdown
!--------------------------------------------------------------------------
        subroutine TAVPMNGDispatcherEncode(this,ds_instr,instr_packet,ierr)
!Encodes a tensor instruction into a plain bytecode.
         implicit none
         class(tavp_mng_dispatcher_t), intent(inout):: this  !inout: TAVP-MNG Dispatcher DSVU
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
!-------------------------------------------------------------------------------
        function TAVPMNGDispatcherMapInstr(this,tens_instr,ierr) result(channel)
!Maps a tensor instruction to a specific lower-level (child) TAVP,
!based on the tensor argument locality and load imbalance.
         implicit none
         integer(INTD):: channel                            !out: channel: TAVP id can be retrieved from this.dispatch_rank(channel)
         class(tavp_mng_dispatcher_t), intent(inout):: this !inout: TAVP-MNG Dispatcher DSVU
         class(tens_instr_t), intent(in):: tens_instr       !in: active tensor instruction
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc,i,num_args
         class(ds_oprnd_t), pointer:: tens_oprnd
         integer(INTD):: owner_ids(0:MAX_TENSOR_OPERANDS-1)
         class(DataDescr_t), pointer:: descr

         channel=-1 !negative value means undefined
         if(tens_instr%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           num_args=tens_instr%get_num_operands(errc)
           if(errc.eq.DSVP_SUCCESS) then
 !Determine which TAVP owns each tensor argument:
            tens_oprnd=>NULL()
            do i=0,num_args-1
             tens_oprnd=>tens_instr%get_operand(i,errc)
             if(errc.eq.DSVP_SUCCESS) then
              select type(tens_oprnd)
              class is(tens_oprnd_t)
               call tens_oprnd%lock()
               if(this%tavp_is_bottom) then
                descr=>tens_oprnd%get_data_descriptor(errc)
                if(errc.eq.0) then
                 if(associated(descr)) then
                  if(descr%is_set(errc,proc_rank=owner_ids(i))) then !TAVP-WRK id
                   if(errc.ne.0) then; errc=-9; exit; endif
                  else
                   owner_ids(i)=-1 !data descriptor is not set yet
                  endif
                 else
                  owner_ids(i)=-1 !data descriptor is not set yet
                 endif
                else
                 errc=-8; exit
                endif
               else
                owner_ids(i)=tens_oprnd%get_owner_id(errc) !TAVP-MNG id
                if(errc.ne.0) then; errc=-7; exit; endif
               endif
               call tens_oprnd%unlock()
              class default
               errc=-6; exit
              end select
              tens_oprnd=>NULL()
             else
              errc=-5; exit
             endif
            enddo
            if(errc.ne.0.and.associated(tens_oprnd)) then !in case of error
             select type(tens_oprnd); class is(tens_oprnd_t); call tens_oprnd%unlock(); end select
             tens_oprnd=>NULL()
            endif
 !Decide which dispatch channel to map the instruction to:
            if(errc.eq.0.and.num_args.gt.0) then
             channel=map_by_arg_order(errc); if(errc.ne.0) errc=-4
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

         function map_by_arg_order(jerr) result(chnl)
          !Selects the dispatch channel by argument locality with argument
          !priority following the ascending argument order: 0,1,2,...
          implicit none
          integer(INTD):: chnl              !out: dispatch channel
          integer(INTD), intent(out):: jerr !out: error code
          integer(INTD):: ja,jj
          logical:: jf

          jerr=0; chnl=-1
          aloop: do ja=0,num_args-1
           do jj=lbound(this%dispatch_rank,1),ubound(this%dispatch_rank,1)
            if(owner_ids(ja).eq.this%dispatch_rank(jj)) then; chnl=jj; exit aloop; endif
           enddo
          enddo aloop
          if(chnl.lt.0) then
           if(DISPATCH_RANDOM) then
            chnl=map_by_random(jerr)
           else
            chnl=map_by_round(jerr)
           endif
          endif
          return
         end function map_by_arg_order

         function map_by_round(jerr) result(chnl)
          implicit none
          integer(INTD):: chnl              !out: dispatch channel
          integer(INTD), intent(out):: jerr !out: error code

          jerr=0; chnl=this%next_channel
          this%next_channel=1+mod(this%next_channel,this%num_ranks)
          return
         end function map_by_round

         function map_by_random(jerr) result(chnl)
          implicit none
          integer(INTD):: chnl              !out: dispatch channel
          integer(INTD), intent(out):: jerr !out: error code
          integer(INTD):: js
          real(8):: jdelta,jrnd

          jerr=0; chnl=-1
          js=size(this%dispatch_rank)
          if(js.gt.0) then
           jdelta=1d0/real(js,8)
           call random_number(jrnd)
           chnl=lbound(this%dispatch_rank,1)+min(int(jrnd/jdelta,INTD),js-1)
          else
           jerr=-1
          endif
          return
         end function map_by_random

        end function TAVPMNGDispatcherMapInstr
!-------------------------------------------------------------------------
        subroutine TAVPMNGDispatcherDispatch(this,tens_instr,channel,ierr)
!Dispatches a tensor instruction to a specific lower-level TAVP
!and encodes it into its bytecode buffer.
         implicit none
         class(tavp_mng_dispatcher_t), intent(inout):: this      !inout: TAVP-MNG Dispatcher DSVU
         class(tens_instr_t), target, intent(inout):: tens_instr !in: defined tensor instruction
         integer(INTD), intent(in):: channel                     !in: dispatch channel: offset in this.bytecode(1:max)
         integer(INTD), intent(out), optional:: ierr             !out: error code
         integer(INTD):: errc,opcode,i,n,home
         type(obj_pack_t):: instr_packet
         integer(INTD), pointer:: out_oprs(:)
         class(ds_oprnd_t), pointer:: oprnd

 !Set the home for orphaned output tensor operands:
         opcode=tens_instr%get_code(errc)
         if(errc.eq.DSVP_SUCCESS) then
          if(opcode.ge.TAVP_ISA_TENS_FIRST.and.opcode.le.TAVP_ISA_TENS_LAST) then !tensor instruction
           if(channel.ge.lbound(this%dispatch_rank,1).and.channel.le.ubound(this%dispatch_rank,1)) then
            out_oprs=>tens_instr%get_output_operands(errc,num_oprs=n)
            if(errc.eq.0) then
             if(n.gt.0) then
              if(this%tavp_is_bottom) then
               home=role_rank !this TAVP-MNG id
              else
               home=this%dispatch_rank(channel) !child TAVP-MNG id
              endif
              do i=lbound(out_oprs,1),ubound(out_oprs,1)
               oprnd=>tens_instr%get_operand(out_oprs(i),errc); if(errc.ne.DSVP_SUCCESS) then; errc=-12; exit; endif
               select type(oprnd)
               class is(tens_oprnd_t)
                if(oprnd%get_owner_id(errc).lt.0) then !no owner
                 if(errc.eq.0) then
                  call oprnd%set_owner_id(home,errc,update_cache=.TRUE.); if(errc.ne.0) then; errc=-11; exit; endif
                 else
                  errc=-10; exit
                 endif
                else
                 if(errc.ne.0) then; errc=-9; exit; endif
                endif
               class default
                errc=-8; exit
               end select
              enddo
             endif
            else
             errc=-7
            endif
           else
            errc=-6
           endif
          endif
         else
          errc=-5
         endif
         if(errc.eq.0) then
 !Encode and dispatch the tensor instruction to the corresponding bytecode channel:
          call this%bytecode(channel)%acquire_packet(instr_packet,errc,preclean=.TRUE.)
          if(errc.eq.PACK_SUCCESS) then
           call this%encode(tens_instr,instr_packet,errc)
           if(errc.eq.0) then
            call this%bytecode(channel)%seal_packet(errc)
            if(errc.eq.PACK_SUCCESS) then
 !Update dispatch stats:
             call this%update_dispatch_count(channel,1_INTL) !update current instruction count for this channel `Collector must decrement this counter
             call this%update_dispatch_flops(channel,tens_instr%get_flops(errc)) !update current Flop count for this channel `Collector must decrement this counter
             if(errc.ne.0) errc=-4
            else
             errc=-3
            endif
           else
            errc=-2
           endif
          else
           if(VERBOSE) then
            write(CONS_OUT,'("#ERROR(TAVP-MNG:Dispatcher.dispatch)[",i6,"]: Unable to acquire packet: Error ",i11)') impir,errc
            flush(CONS_OUT)
           endif
           errc=-1
          endif
         endif
         if(VERBOSE.and.errc.ne.0) then
          write(CONS_OUT,'("#ERROR(TAVP-MNG:Dispatcher.dispatch)[",i6,"]: Dispatch to channel ",i2," failed: Error ",i11)')&
          &impir,channel,errc
          flush(CONS_OUT)
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDispatcherDispatch
!-----------------------------------------------------------
        subroutine TAVPMNGDispatcherIssue(this,channel,ierr)
!Issues the encoded bytecode to a specific dispatch channel.
         implicit none
         class(tavp_mng_dispatcher_t), intent(inout):: this !inout: TAVP-MNG Dispatcher DSVU
         integer(INTD), intent(in):: channel                !in: offset in this.bytecode(1:max)
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc
         integer(INTL):: npck

         errc=0
         if(channel.ge.lbound(this%dispatch_rank,1).and.channel.le.ubound(this%dispatch_rank,1)) then
          if(this%comm_hl(channel)%is_clean(errc)) then
           if(errc.eq.0) then
            call this%bytecode(channel)%send(this%dispatch_rank(channel),this%comm_hl(channel),errc,&
                                            &tag=TAVP_DISPATCH_TAG,comm=this%dispatch_comm)
            if(errc.ne.0) errc=-4
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(DEBUG.gt.0) then
          if(errc.eq.0) then
           npck=this%bytecode(channel)%get_num_packets()
           write(CONS_OUT,'("#MSG(TAVP-MNG:Dispatcher.issue)[",i6,"]: Issued ",i6," instructions to channel ",i2,": Rank = ",i4)')&
           &impir,npck,channel,this%dispatch_rank(channel)
          else
           write(CONS_OUT,'("#ERROR(TAVP-MNG:Dispatcher.issue)[",i6,"]: Unable to issue bytecode to channel ",i2,": Error ",i11)')&
           &impir,channel,errc
           write(CONS_OUT,'("#INFO(TAVP-MNG:Dispatcher.issue)[",i6,"]: Target rank = ",i6,", Target comm = ",i11)')&
           &impir,this%dispatch_rank(channel),this%dispatch_comm
          endif
          flush(CONS_OUT)
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGDispatcherIssue
!----------------------------------------------------------------------------
        function TAVPMNGDispatcherSyncIssue(this,ierr,channel) result(synced)
!Synchronizes asynchronous instruction bytecode issue for channel <channel>.
!If <channel> is not provided, all active channels will be sychronized.
!Upon success, the corresponding bytecode buffer(s) will be cleaned.
         implicit none
         logical:: synced                                   !out: TRUE if the channel(s) has(ve) been synchronized
         class(tavp_mng_dispatcher_t), intent(inout):: this !inout: TAVP-MNG Dispatcher DSVU
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD), intent(in), optional:: channel      !in: specific dispatch channel to synchronize
         integer(INTD):: errc,i

         errc=0; synced=.FALSE.
         if(present(channel)) then !sync a specific channel
          if(channel.ge.lbound(this%dispatch_rank,1).and.channel.le.ubound(this%dispatch_rank,1)) then
           call this%comm_hl(channel)%wait(errc); if(errc.ne.0) errc=-9
           if(errc.eq.0) then; call this%comm_hl(channel)%clean(errc); if(errc.ne.0) errc=-8; endif
           if(errc.eq.0) then; call this%bytecode(channel)%clean(errc); if(errc.ne.0) errc=-7; endif
          else
           errc=-6
          endif
         else !sync all channels
          do i=lbound(this%dispatch_rank,1),ubound(this%dispatch_rank,1)
           if(this%comm_hl(i)%is_active(errc)) then
            if(errc.eq.PACK_SUCCESS) then
             call this%comm_hl(i)%wait(errc); if(errc.ne.0) errc=-5
             if(errc.eq.0) then; call this%comm_hl(i)%clean(errc); if(errc.ne.0) errc=-4; endif
             if(errc.eq.0) then; call this%bytecode(i)%clean(errc); if(errc.ne.0) errc=-3; endif
            else
             errc=-2
            endif
           else
            if(errc.ne.PACK_SUCCESS) errc=-1
           endif
           if(errc.ne.0) exit
          enddo
         endif
         if(errc.eq.0) synced=.TRUE.
         if(DEBUG.gt.0) then
          if(errc.eq.0) then
           if(present(channel)) then
            write(CONS_OUT,'("#MSG(TAVP-MNG:Dispatcher.sync_issue)[",i6,"]: Synced channel ",i2,": Rank = ",i4)')&
            &impir,channel,this%dispatch_rank(channel)
           else
            write(CONS_OUT,'("#MSG(TAVP-MNG:Dispatcher.sync_issue)[",i6,"]: Synced all channels")') impir
           endif
          else
           write(CONS_OUT,'("#ERROR(TAVP-MNG:Dispatcher.sync_issue)[",i6,"]: Sync error ",i11)') impir,errc
          endif
          flush(CONS_OUT)
         endif
         if(present(ierr)) ierr=errc
         return
        end function TAVPMNGDispatcherSyncIssue
!---------------------------------------------------------------------
        subroutine TAVPMNGUpdateDispatchCount(this,channel,delta,ierr)
!Updates the current instruction dispatch count.
         implicit none
         class(tavp_mng_dispatcher_t), intent(inout):: this !inout: TAVP-MNG dispatcher DSVU
         integer(INTD), intent(in):: channel                !in: channel to update count for
         integer(INTL), intent(in):: delta                  !in: delta to update
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc

         errc=0
         if(channel.ge.lbound(this%dispatch_rank,1).and.channel.le.ubound(this%dispatch_rank,1)) then
!!!$OMP ATOMIC UPDATE SEQ_CST
!$OMP ATOMIC UPDATE
          this%dispatch_count(channel)=this%dispatch_count(channel)+delta
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGUpdateDispatchCount
!---------------------------------------------------------------------
        subroutine TAVPMNGUpdateDispatchFlops(this,channel,delta,ierr)
!Updates the current instruction dispatch Flop count.
         implicit none
         class(tavp_mng_dispatcher_t), intent(inout):: this !inout: TAVP-MNG dispatcher DSVU
         integer(INTD), intent(in):: channel                !in: channel to update Flop count for
         real(8), intent(in):: delta                        !in: delta to update
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc

         errc=0
         if(channel.ge.lbound(this%dispatch_rank,1).and.channel.le.ubound(this%dispatch_rank,1)) then
!!!$OMP ATOMIC UPDATE SEQ_CST
!$OMP ATOMIC UPDATE
          this%dispatch_flops(channel)=this%dispatch_flops(channel)+delta
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGUpdateDispatchFlops
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
         integer(INTD):: errc,ier,thid,uid
         logical:: active,stopping
         class(dsvp_t), pointer:: dsvp
         class(tavp_mng_t), pointer:: tavp

         errc=0; thid=omp_get_thread_num(); uid=this%get_id()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Replicator started as DSVU # ",i2," (thread ",i2,")")') impir,uid,thid
          flush(CONS_OUT)
         endif
!Reserve a bytecode buffer:
         call this%bytecode%reserve_mem(ier,MAX_BYTECODE_SIZE,MAX_BYTECODE_INSTR); if(ier.ne.0.and.errc.eq.0) errc=-1
!Initialize queues:
         call this%init_queue(this%num_ports,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
!Initialize the control list:
         ier=this%ctrl_list%init(this%control_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-1
!Set up tensor argument cache and wait on other TAVP units:
         tavp=>NULL(); dsvp=>this%get_dsvp(); select type(dsvp); class is(tavp_mng_t); tavp=>dsvp; end select
         if(associated(tavp)) then
          this%arg_cache=>tavp%tens_cache
          call tavp%sync_units(errc,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
         else
          this%arg_cache=>NULL(); if(errc.eq.0) errc=-1
         endif
!Work loop:
         active=((errc.eq.0).and.(this%repl_comm.ne.MPI_COMM_NULL)); stopping=(.not.active)
         wloop: do while(active)
          !`Implement
          stopping=.TRUE.
          active=(.not.stopping)
         enddo wloop
!Record the error:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(TAVP-MNG)[",i6,"]: Replicator error ",i11," by thread ",i2)')&
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
         integer(INTD):: errc,ier,thid,uid

         errc=0; thid=omp_get_thread_num(); uid=this%get_id()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Replicator stopped as DSVU # ",i2," (thread ",i2,")")') impir,uid,thid
          flush(CONS_OUT)
         endif
!Release the tensor argument cache pointer:
         this%arg_cache=>NULL()
!Deactivate the control list:
         ier=this%ctrl_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%ctrl_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-5
           ier=this%ctrl_list%delete_all()
          endif
          ier=this%ctrl_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-4
         else
          if(errc.eq.0) errc=-3
         endif
!Release queues:
         call this%release_queue(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-2
!Release the bytecode buffer:
         call this%bytecode%destroy(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
!Record an error, if any:
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
         class(tavp_mng_collector_t), intent(inout):: this !inout: TAVP-MNG Collector DSVU
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc,ier,thid,cnt,iec,sts,opcode,i,n,uid
         integer(INTL):: pid,cid
         logical:: active,stopping,matched,evicted
         class(dsvp_t), pointer:: dsvp
         class(tavp_mng_t), pointer:: tavp
         class(tens_rcrsv_t), pointer:: tensor
         class(tens_instr_t), pointer:: tens_instr,par_instr
         type(tens_entry_mng_ref_t):: cache_entries(1:MAX_TENSOR_OPERANDS)
         class(*), pointer:: uptr
         type(list_pos_t):: list_pos

         errc=0; thid=omp_get_thread_num(); uid=this%get_id()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Collector started as DSVU # ",i2," (thread ",i2,")")') impir,uid,thid
          flush(CONS_OUT)
         endif
!Initialize queues:
         call this%init_queue(this%num_ports,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-68
!Initialize the matching list iterator:
         ier=this%mat_list%init(this%matching_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-67
!Initialize the deferred list iterator:
         ier=this%def_list%init(this%deferred_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-66
!Initialize the retired list iterator:
         ier=this%ret_list%init(this%retired_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-65
!Initialize the control list:
         ier=this%ctrl_list%init(this%control_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-64
!Initialize parent instruction map iterator:
         ier=this%parent_instr%init(this%parent_instr_map); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-63
!Set up tensor argument cache and wait on other TAVP units:
         tavp=>NULL(); dsvp=>this%get_dsvp(); select type(dsvp); class is(tavp_mng_t); tavp=>dsvp; end select
         if(associated(tavp)) then
          this%arg_cache=>tavp%tens_cache
          call tavp%sync_units(errc,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-62
         else
          this%arg_cache=>NULL(); if(errc.eq.0) errc=-61
         endif
!Work loop:
         active=(errc.eq.0); stopping=(.not.active)
         wloop: do while(active)
 !Process a new batch of parent tensor instructions (and possibly control instructions):
  !Get a new batch of parent tensor instructions from Decomposer (port 0) into the main queue:
          ier=this%iqueue%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-60; exit wloop; endif
          ier=this%flush_port(0,num_moved=i); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-59; exit wloop; endif
          if(DEBUG.gt.0.and.i.gt.0) then
           write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Collector unit ",i2," received ",i9," parent instructions from Decomposer")')&
           &impir,uid,i
           !ier=this%iqueue%reset(); ier=this%iqueue%scanp(action_f=tens_instr_print) !print all instructions
           flush(CONS_OUT)
          endif
  !Register parent tensor instructions and move them into the matching list, move other instructions elsewhere:
          ier=this%mat_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-58; exit wloop; endif
          ier=this%iqueue%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-57; exit wloop; endif
          ier=this%iqueue%get_status()
          ploop: do while(ier.eq.GFC_IT_ACTIVE)
   !Extract a parent tensor instruction (.error_code field contains the number of child subinstructions):
           uptr=>this%iqueue%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-56; exit wloop; endif
           tens_instr=>NULL(); select type(uptr); class is(tens_instr_t); tens_instr=>uptr; end select
           if((.not.associated(tens_instr)).and.errc.eq.0) then; errc=-55; exit wloop; endif !trap
           pid=tens_instr%get_id(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-54; exit wloop; endif
           sts=tens_instr%get_status(ier,cnt); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-53; exit wloop; endif
           opcode=tens_instr%get_code(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-52; exit wloop; endif
           if(DEBUG.gt.0) then
            write(CONS_OUT,'("#DEBUG(TAVP-MNG:Collector): A parent instruction was received:")')
            call tens_instr%print_it(dev_id=CONS_OUT)
            flush(CONS_OUT)
           endif
   !Register decomposable parent tensor instructions and move them into the matching list:
           if(opcode.ge.TAVP_ISA_TENS_FIRST.and.opcode.le.TAVP_ISA_TENS_LAST) then
            if(stopping.and.errc.eq.0) then; errc=-51; exit wloop; endif !trap
            call tens_instr%set_status(sts,ier,DSVP_SUCCESS) !.error_code field was used as children count, reset it back to DSVP_SUCCESS
            if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-50; exit wloop; endif
            if(cnt.le.0.and.errc.eq.0) then; errc=-49; exit wloop; endif !trap
            ier=this%iqueue%move_elem(this%mat_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-48; exit wloop; endif
            ier=this%mat_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-47; exit wloop; endif
            call this%register_instr(pid,cnt,this%mat_list%bookmark(),ier)
            if(ier.ne.0.and.errc.eq.0) then; errc=-46; exit wloop; endif
   !Move control instructions into the control list (or delete them):
           elseif(opcode.ge.TAVP_ISA_CTRL_FIRST.and.opcode.le.TAVP_ISA_CTRL_LAST) then
            if(stopping.and.errc.eq.0) then; errc=-45; exit wloop; endif !trap
            select case(opcode)
            case(TAVP_INSTR_CTRL_STOP)
             call tens_instr%set_status(DS_INSTR_COMPLETED,ier)
             if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-44; exit wloop; endif
             ier=this%ctrl_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-43; exit wloop; endif
             ier=this%ctrl_list%append(tens_instr); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-42; exit wloop; endif
             ier=tavp%cdecoder%load_port(0,this%ctrl_list); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-41; exit wloop; endif
             ier=this%iqueue%move_elem(this%ctrl_list)
             stopping=.TRUE.
            case(TAVP_INSTR_CTRL_RESUME,TAVP_INSTR_CTRL_DUMP_CACHE)
             call tens_instr%set_status(DS_INSTR_RETIRED,ier,DSVP_SUCCESS)
             if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-40; exit wloop; endif
             ier=this%iqueue%delete(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-39; exit wloop; endif
            case default !unknown control instruction
             if(errc.eq.0) then; errc=-38; exit wloop; endif
            end select
   !Other instructions are not expected here:
           else !other instructions should not appear here
            if(errc.eq.0) then; errc=-37; exit wloop; endif
           endif
           ier=this%iqueue%get_status()
          enddo ploop
          if(ier.ne.GFC_IT_EMPTY.and.errc.eq.0) then; errc=-36; exit wloop; endif
 !Process a new batch of child instructions (subinstructions):
  !Move previously received subinstructions from the deferred list back into the main queue:
          ier=this%iqueue%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-35; exit wloop; endif
          if(this%iqueue%get_status().ne.GFC_IT_EMPTY.and.errc.eq.0) then; errc=-34; exit wloop; endif !trap
          ier=this%def_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-33; exit wloop; endif
          if(this%def_list%get_status().eq.GFC_IT_ACTIVE) then
           ier=this%def_list%move_list(this%iqueue); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-32; exit wloop; endif
           ier=this%iqueue%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-31; exit wloop; endif
          endif
  !Append a new batch of subinstructions from cDecoder (port 1) into the main queue:`Do I need MAX_COLLECT_INSTR here?
          ier=this%flush_port(1,num_moved=i); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-30; exit wloop; endif
          if(DEBUG.gt.0.and.i.gt.0) then
           write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Collector unit ",i2," received ",i9," subinstructions from lower level")')&
           &impir,uid,i
           !ier=this%iqueue%reset(); ier=this%iqueue%scanp(action_f=tens_instr_print) !print all instructions
           flush(CONS_OUT)
          endif
  !Match subinstructions with their respective parent instructions and possibly retire the latter to the retired list:
          ier=this%iqueue%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-29; exit wloop; endif
          ier=this%iqueue%get_status()
          mloop: do while(ier.eq.GFC_IT_ACTIVE)
   !Extract a child tensor instruction:
           uptr=>this%iqueue%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-28; exit wloop; endif
           tens_instr=>NULL(); select type(uptr); class is(tens_instr_t); tens_instr=>uptr; end select
           if((.not.associated(tens_instr)).and.errc.eq.0) then; errc=-27; exit wloop; endif !trap
           opcode=tens_instr%get_code(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-26; exit wloop; endif
           if(opcode.eq.TAVP_INSTR_TENS_DESTROY) then !remove tensor persistency for deleted tensors
            call tens_instr%remove_persistency(ier); if(ier.ne.0.and.errc.eq.0) then; errc=-25; exit wloop; endif
           endif
           cid=tens_instr%get_id(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-24; exit wloop; endif
           sts=tens_instr%get_status(ier,iec); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-23; exit wloop; endif
   !Map the child tensor instruction to its parent instruction and decrement the reference count:
           call list_pos%clean() !list position will refer to the position in the matching list
           matched=this%match_subinstr(cid,ier,iec,list_pos) !list position is only set when the reference count is zero
           if(ier.ne.0.and.errc.eq.0) then; errc=-22; exit wloop; endif
           if(matched) then
            if(DEBUG.gt.0) then
             write(CONS_OUT,'("#DEBUG(TAVP-MNG:Collector): Matched a child instruction:")')
             call tens_instr%print_it(dev_id=CONS_OUT)
             flush(CONS_OUT)
            endif
   !Move parent tensor instruction to the retired list when all child instructions have been matched:
            if(list_pos%is_set()) then !list position set: all child instructions have been matched
             ier=this%mat_list%jump(list_pos); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-21; exit wloop; endif
             uptr=>this%mat_list%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-20; exit wloop; endif
             par_instr=>NULL(); select type(uptr); class is(tens_instr_t); par_instr=>uptr; end select
             if((.not.associated(par_instr)).and.errc.eq.0) then; errc=-19; exit wloop; endif !trap
             call par_instr%set_status(DS_INSTR_COMPLETED,ier)
             if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-18; exit wloop; endif
             ier=this%mat_list%move_elem(this%ret_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-17; exit wloop; endif
             par_instr=>NULL()
            endif
   !Unregister the matched child instruction within TAVP:
            call tavp%unregister_instr(cid,ier); if(ier.ne.0.and.errc.eq.0) then; errc=-16; exit wloop; endif
   !Delete the matched child instruction and evict unneeded remote tensor cache entries (if any):
    !Delete the matched child instruction:
            call tens_instr%extract_cache_entries(cache_entries,n,ier)
            if(ier.ne.0.and.errc.eq.0) then; errc=-15; exit wloop; endif
            call tens_instr%set_status(DS_INSTR_RETIRED,ier)
            if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-14; exit wloop; endif
            ier=this%iqueue%delete(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-13; exit wloop; endif
    !Evict remote tensor cache entries with zero reference count (if any):
            do i=1,n
             tensor=>cache_entries(i)%cache_entry%get_tensor(ier)
             if(ier.eq.0) then
              evicted=this%arg_cache%evict(tensor,ier,decr_use=.TRUE.); if(ier.ne.0.and.errc.eq.0) errc=-12
             else
              if(errc.eq.0) errc=-11
             endif
             cache_entries(i)%cache_entry=>NULL()
            enddo
            if(errc.ne.0) exit wloop
           else !unmatched child instruction
   !Move unmatched child instruction to the deferred list:
            ier=this%iqueue%move_elem(this%def_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-10; exit wloop; endif
           endif
           ier=this%iqueue%get_status()
          enddo mloop
          if(ier.ne.GFC_IT_EMPTY.and.errc.eq.0) then; errc=-9; exit wloop; endif
 !Check whether the Collector should stop:
          ier=this%mat_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-8; exit wloop; endif
          ier=this%def_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-7; exit wloop; endif
          if(this%mat_list%get_status().eq.GFC_IT_EMPTY.and.this%def_list%get_status().eq.GFC_IT_EMPTY) then
           if(stopping) then
            ier=this%ret_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-6; exit wloop; endif
            ier=this%ctrl_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-5; exit wloop; endif
            if(this%ctrl_list%get_status().eq.GFC_IT_ACTIVE) then
             ier=this%ctrl_list%move_list(this%ret_list) !append deferred control instructions to the end
             if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-4; exit wloop; endif
            endif
            active=.FALSE.
           endif
          endif
 !Pass retired parent instructions from the retired list to Retirer (port 0):
          ier=this%ret_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-3; exit wloop; endif
          if(this%ret_list%get_status().eq.GFC_IT_ACTIVE) then
           ier=tavp%retirer%load_port(0,this%ret_list); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-2; exit wloop; endif
          endif
         enddo wloop
!Record the error:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(TAVP-MNG)[",i6,"]: Collector error ",i11," by thread ",i2)')&
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
         class(tavp_mng_collector_t), intent(inout):: this !inout: TAVP-MNG Collector DSVU
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc,ier,thid,uid

         errc=0; thid=omp_get_thread_num(); uid=this%get_id()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-MNG)[",i6,"]: Collector stopped as DSVU # ",i2," (thread ",i2,")")') impir,uid,thid
          flush(CONS_OUT)
         endif
!Release the tensor argument cache pointer:
         this%arg_cache=>NULL()
!Release parent instruction map iterator:
         ier=this%parent_instr%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%parent_instr%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-16
           ier=this%parent_instr%delete_all()
          endif
          ier=this%parent_instr%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-15
         else
          if(errc.eq.0) errc=-14
         endif
!Deactivate the control list:
         ier=this%ctrl_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%ctrl_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-13
           ier=this%ctrl_list%delete_all()
          endif
          ier=this%ctrl_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-12
         else
          if(errc.eq.0) errc=-11
         endif
!Deactivate the retired list:
         ier=this%ret_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%ret_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-10
           ier=this%ret_list%delete_all()
          endif
          ier=this%ret_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-9
         else
          if(errc.eq.0) errc=-8
         endif
!Deactivate the deferred list:
         ier=this%def_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%def_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-7
           ier=this%def_list%delete_all()
          endif
          ier=this%def_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-6
         else
          if(errc.eq.0) errc=-5
         endif
!Deactivate the matching list:
         ier=this%mat_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%mat_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-4
           ier=this%mat_list%delete_all()
          endif
          ier=this%mat_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-3
         else
          if(errc.eq.0) errc=-2
         endif
!Release queues:
         call this%release_queue(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
!Record an error, if any:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGCollectorShutdown
!----------------------------------------------------------------------------------------
        subroutine TAVPMNGCollectorRegisterInstr(this,instr_id,child_count,list_pos,ierr)
!Registers a parent instruction with the Collector.
         implicit none
         class(tavp_mng_collector_t), intent(inout):: this !inout: TAVP-MNG Collector
         integer(INTL), intent(in):: instr_id              !in: parent instruction id
         integer(INTD), intent(in):: child_count           !in: number of subinstructions spawned by this parent instruction (>=0)
         type(list_pos_t), intent(in):: list_pos           !in: bookmarked position of the parent instruction in Collector's matching list
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc

         errc=0
         if(child_count.ge.0) then
          errc=this%parent_instr%search(GFC_DICT_ADD_IF_NOT_FOUND,cmp_integers,instr_id,&
                                       &tavp_mng_collector_entry_t(child_count,list_pos))
          if(errc.eq.GFC_NOT_FOUND) then; errc=0; else; errc=-2; endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGCollectorRegisterInstr
!---------------------------------------------------------------------
        subroutine TAVPMNGCollectorUnregisterInstr(this,instr_id,ierr)
!Explicitly unregisters a parent instruction with the Collector.
         implicit none
         class(tavp_mng_collector_t), intent(inout):: this !inout: TAVP-MNG Collector
         integer(INTL), intent(in):: instr_id              !in: parent instruction id
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc

         errc=this%parent_instr%search(GFC_DICT_DELETE_IF_FOUND,cmp_integers,instr_id)
         if(errc.eq.GFC_FOUND) then; errc=0; else; errc=-1; endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGCollectorUnregisterInstr
!----------------------------------------------------------------------------------------------------------
        function TAVPMNGCollectorMatchSubinstr(this,subinstr_id,ierr,subinstr_err,list_pos) result(matched)
!Matches a subinstruction against its parent instruction in Collector's dictionary and
!decrements the number of active subinstructions associated with the parent instruction.
!If that reference count becomes zero and an optional <list_pos> argument is present,
!the entry will be destroyed and the bookmark to the parent instruction position in
!Collector's matching list will be returned in <list_pos>. If <list_pos> is absent,
!the entry will not be destroyed, making it necessary to destroy it explicitly later.
!If <list_pos> is present but the reference count is still positive, it will return nothing.
!An optional subinstruction error code will propagate to the parent instruction.
         implicit none
         logical:: matched                                  !out: matched or not
         class(tavp_mng_collector_t), intent(inout):: this  !inout: TAVP-MNG Collector
         integer(INTL), intent(in):: subinstr_id            !in: subinstruction id
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD), intent(in), optional:: subinstr_err !in: subinstruction error code
         type(list_pos_t), intent(out), optional:: list_pos !out: list position bookmark (parent instruction position in Collector's matching list)
         integer(INTD):: errc,ier,sts
         integer(INTL):: parent_instr_id
         class(dsvp_t), pointer:: dsvp
         class(tavp_mng_collector_entry_t), pointer:: col_entry
         class(*), pointer:: uptr
         class(tens_instr_t), pointer:: tens_instr
         type(list_iter_t):: lit

         matched=.FALSE.; dsvp=>this%get_dsvp(errc)
         if(errc.eq.DSVP_SUCCESS) then
!Look up the parent instruction ID:
          select type(dsvp)
          class is(tavp_mng_t)
           parent_instr_id=dsvp%map_instr(subinstr_id,errc); if(errc.ne.0) errc=-15
          class default
           errc=-14
          end select
!Find the corresponding entry in Collector's dictionary:
          if(errc.eq.0.and.parent_instr_id.ge.0) then
           errc=this%parent_instr%search(GFC_DICT_FETCH_IF_FOUND,cmp_integers,parent_instr_id,value_out=uptr)
           if(errc.eq.GFC_FOUND) then
            errc=0; matched=.TRUE.
            select type(uptr)
            class is(tavp_mng_collector_entry_t)
             col_entry=>uptr
            class default
             errc=-13
            end select
!Decrement the reference count:
            if(errc.eq.0) then
             if(col_entry%children_count.gt.0) then
              col_entry%children_count=col_entry%children_count-1
!Propagate the error code to the parent instruction:
              if(present(subinstr_err)) then
               errc=lit%init(this%matching_list)
               if(errc.eq.GFC_SUCCESS) then
                errc=lit%jump(col_entry%list_elem)
                if(errc.eq.GFC_SUCCESS) then
                 uptr=>lit%get_value(errc)
                 if(errc.eq.GFC_SUCCESS) then
                  tens_instr=>NULL(); select type(uptr); class is(tens_instr_t); tens_instr=>uptr; end select
                  if(associated(tens_instr)) then
                   sts=tens_instr%get_status(errc,ier)
                   if(errc.eq.DSVP_SUCCESS) then
                    if(ier.eq.DSVP_SUCCESS) then
                     call tens_instr%set_status(sts,errc,subinstr_err); if(errc.ne.DSVP_SUCCESS) errc=-12
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
                ier=lit%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=-7
               else
                errc=-6
               endif
              endif
!Destroy the entry if the reference count is zero (only when <list_pos> is present!):
              if(present(list_pos).and.errc.eq.0) then
               if(col_entry%children_count.eq.0) then
                list_pos=col_entry%list_elem
                errc=this%parent_instr%search(GFC_DICT_DELETE_IF_FOUND,cmp_integers,parent_instr_id)
                if(errc.eq.GFC_FOUND) then; errc=0; else; errc=-5; endif
               endif
              endif
             else
              errc=-4
             endif
            endif
           else
            if(errc.eq.GFC_NOT_FOUND) then; errc=-3; else; errc=-2; endif
           endif
          endif
         else
          errc=-1
         endif
         if(errc.ne.0) write(CONS_OUT,'("#ERROR(TAVP-MNG:Collector.match_subinstr): Error ",i11)') errc !debug
         if(present(ierr)) ierr=errc
         return
        end function TAVPMNGCollectorMatchSubinstr
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
                         &msg_tag=TAVP_DISPATCH_TAG,acceptor=decode_acceptor,acceptor_port_id=0)
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
                            &msg_tag=TAVP_LOCATE_TAG,acceptor=decode_acceptor,acceptor_port_id=1)
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
#if 0
  !Replicator:
                   replicator_conf=tavp_mng_replicator_conf_t(role_comm)
                   call this%replicator%configure(replicator_conf,errc)
                   if(errc.eq.0) then
                    num_units=num_units+1
  !Replicating-decoder:
                    decode_acceptor=>this%replicator
                    decoder_conf=tavp_mng_decoder_conf_t(source_comm=role_comm,source_rank=MPI_ANY_SOURCE,&
                                &msg_tag=TAVP_REPLICA_TAG,acceptor=decode_acceptor,acceptor_port_id=1)
                    call this%rdecoder%configure(decoder_conf,errc)
                    if(errc.eq.0) then
                     num_units=num_units+1
#endif
  !Collecting-decoder:
                     decode_acceptor=>this%collector
                     decoder_conf=tavp_mng_decoder_conf_t(source_comm=conf%collect_comm,source_rank=MPI_ANY_SOURCE,&
                                 &msg_tag=TAVP_COLLECT_TAG,acceptor=decode_acceptor,acceptor_port_id=1)
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
#if 0
                              call this%set_unit(this%replicator,errc)
                              if(errc.eq.DSVP_SUCCESS) then
                               call this%set_unit(this%rdecoder,errc)
                               if(errc.eq.DSVP_SUCCESS) then
#endif
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
#if 0
                               else
                                errc=-23
                               endif
                              else
                               errc=-22
                              endif
#endif
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
#if 0
                    else
                     errc=-12
                    endif
                   else
                    errc=-11
                   endif
#endif
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
!--------------------------------------------------------------------
        subroutine TAVPMNGRegisterInstr(this,child_id,parent_id,ierr)
!Registers a new child instruction (subinstruction).
         implicit none
         class(tavp_mng_t), intent(inout):: this     !inout: TAVP-MNG
         integer(INTL), intent(in):: child_id        !in: child instruction id (key)
         integer(INTL), intent(in):: parent_id       !in: associated parent instruction id (value)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,ier
         type(dictionary_iter_t):: dit

         if(this%get_status(errc).ne.DSVP_STAT_OFF) then
          if(errc.eq.DSVP_SUCCESS) then
!$OMP CRITICAL (TAVP_MNG_INSTR)
           errc=dit%init(this%instr_map)
           if(errc.eq.GFC_SUCCESS) then
            errc=dit%search(GFC_DICT_ADD_IF_NOT_FOUND,cmp_integers,child_id,parent_id)
            if(errc.eq.GFC_NOT_FOUND) then; errc=0; else; errc=-3; endif
            ier=dit%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
           endif
!$OMP END CRITICAL (TAVP_MNG_INSTR)
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGRegisterInstr
!------------------------------------------------------------
        subroutine TAVPMNGUnregisterInstr(this,child_id,ierr)
!Unregisters a child instruction (subinstruction).
         implicit none
         class(tavp_mng_t), intent(inout):: this     !inout: TAVP-MNG
         integer(INTL), intent(in):: child_id        !in: child instruction id (key)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,ier
         type(dictionary_iter_t):: dit

         if(this%get_status(errc).ne.DSVP_STAT_OFF) then
          if(errc.eq.DSVP_SUCCESS) then
!$OMP CRITICAL (TAVP_MNG_INSTR)
           errc=dit%init(this%instr_map)
           if(errc.eq.GFC_SUCCESS) then
            errc=dit%search(GFC_DICT_DELETE_IF_FOUND,cmp_integers,child_id)
            if(errc.eq.GFC_FOUND) then; errc=0; else; errc=-3; endif
            ier=dit%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
           endif
!$OMP END CRITICAL (TAVP_MNG_INSTR)
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPMNGUnregisterInstr
!---------------------------------------------------------------
        function TAVPMNGMapInstr(this,child_id,ierr) result(pid)
!Returns the parent instruction ID (>=0) for a subinstruction.
!If the parent instruction is not found, returns -1 (but no error).
         implicit none
         integer(INTL):: pid                         !out: parent instruction id (>=0 or -1 if not found)
         class(tavp_mng_t), intent(in):: this        !in: TAVP-MNG
         integer(INTL), intent(in):: child_id        !in: child instruction id (key)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,ier
         type(dictionary_iter_t):: dit
         class(*), pointer:: uptr

         pid=-1_INTL
         if(this%get_status(errc).ne.DSVP_STAT_OFF) then
          if(errc.eq.DSVP_SUCCESS) then
!$OMP CRITICAL (TAVP_MNG_INSTR)
           errc=dit%init(this%instr_map)
           if(errc.eq.GFC_SUCCESS) then
            uptr=>NULL()
            errc=dit%search(GFC_DICT_FETCH_IF_FOUND,cmp_integers,child_id,value_out=uptr)
            if(errc.eq.GFC_FOUND) then
             if(associated(uptr)) then
              errc=0
              select type(uptr)
              type is(integer(INTL))
               pid=uptr
              class default
               errc=-6
              end select
             else
              errc=-5
             endif
            else
             if(errc.eq.GFC_NOT_FOUND) then; errc=0; else; errc=-4; endif
            endif
            ier=dit%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-3
           endif
!$OMP END CRITICAL (TAVP_MNG_INSTR)
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(errc.ne.0) write(CONS_OUT,'("#ERROR(TAVP-MNG.map_instr): Error ",i11)') errc !debug
         if(present(ierr)) ierr=errc
         return
        end function TAVPMNGMapInstr
!---------------------------------------------------
        function TAVPMNGIsTop(this,ierr) result(res)
!Returns TRUE if the TAVP-MNG is the root of the TAVP-MNG tree.
         implicit none
         logical:: res                               !out: result
         class(tavp_mng_t), intent(in):: this        !in: active TAVP-MNG
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         res=.FALSE.
         if(this%get_status(errc).eq.DSVP_STAT_ON) then
          if(errc.eq.DSVP_SUCCESS) then
           res=(this%retirer%retire_comm.ne.role_comm) !top TAVP-MNG retires instructions outside (to Driver)
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TAVPMNGIsTop
!------------------------------------------------------
        function TAVPMNGIsBottom(this,ierr) result(res)
!Returns TRUE if the TAVP-MNG is a leaf of the TAVP-MNG tree.
         implicit none
         logical:: res                               !out: result
         class(tavp_mng_t), intent(in):: this        !in: active TAVP-MNG
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         res=.FALSE.
         if(this%get_status(errc).eq.DSVP_STAT_ON) then
          if(errc.eq.DSVP_SUCCESS) then
           res=(this%dispatcher%dispatch_comm.ne.role_comm) !leaf (bottom) TAVP-MNG dispatches instructions outside (to Workers)
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TAVPMNGIsBottom

       end module tavp_manager
