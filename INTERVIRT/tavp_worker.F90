!ExaTENSOR: TAVP-Worker (TAVP-WRK) implementation
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2018/04/12

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
!NOTES:
! # TENSOR INSTRUCTION FORMAT:
!    0. Instruction id;
!    1. Instruction code (opcode);
!    2. Instruction status;
!    3. Instruction error code;
!    4. Instruction control field (optional);
!    5. Instruction operands (optional):
!       {Owner_id,Read_count,Write_count,Tensor} for each tensor operand.
! # TENSOR INSTRUCTION ISSUE:
!   (1) Check that the input tensor operands are not currently updated/undefined locally;
!   (2) Check that the output tensor operand(s) is(are) not currently in use locally;
!   (3) Upon first appearance of each distinct output tensor operand (ref_count=1),
!       create and initialize to zero a local accumulator tensor of the same shape/layout.
!   (4) Upon any appearance of each distinct output tensor operand, create and initialize
!       to zero a new local temporary tensor of the same shape/layout. Then replace the
!       original output tensor operand with the just created local temporary tensor.
!       Subsequently, inject a TENSOR_ADD instruction in order to accumulate the local
!       temporary tensor into the local accumulator tensor.
!       Name mangling rules for the output tensor substitution:
!        Accumulator: "TENSOR3" --> "TENSOR3#0"
!        Temporary: "TENSOR3" --> "TENSOR3#1", "TENSOR3#2", "TENSOR3#3", etc.
!   (5) Upon completion of the last local tensor instruction writing to a specific output
!       tensor, inject a TENSOR_ADD instruction in order to upload the local accumulator
!       tensor into the corresponding persistent (local or remote) output tensor.
!    Example:
!     D3 += L1 * R1 -> D3#0 = 0, D3#1 = 0, D3#1 += L1 * R1, D3#0 += D3#1, ~D3#1;
!     D3 += L2 * R2 ->           D3#2 = 0, D3#2 += L2 * R2, D3#0 += D3#2, ~D3#2, D3 += D3#0, ~D3#0;
! # TENSOR INSTRUCTION DEPENDENCY:
!   (1) A tensor instruction has a data dependency if any of the following applies:
!       (a) Its input operand has currently a non-zero WRITE count;
!       (b) Its output operand has currently a non-zero READ count;
!       Such a tensor instruction will be issued in the deferred queue until the data dependency
!       is gone. Yet, it will affect the READ/WRITE counters of its tensor operands.
!   (2) A tensor instruction has a blocking data dependency if at least one of its
!       operands has both READ and WRITE counters positive. Such a tensor instruction
!       will not be issued and will not affect the READ/WRITE counters.

        use virta
        use gfc_base
        use gfc_list
        use gfc_vector
        use gfc_dictionary
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !default console output
        integer(INTD), private:: DEBUG=1    !debugging mode
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
 !Distributed memory space:
        integer(INTD), parameter, private:: TAVP_WRK_NUM_WINS=1 !number of MPI windows in the DDSS distributed space
 !On-node pinned Host memory buffer:
        integer(INTL), parameter, private:: TAVP_WRK_HOST_BUF_SIZE=1_INTL*(1024_INTL*1024_INTL*1024_INTL) !default Host buffer size in bytes
 !Elementary tensor instruction granularity classification:
        real(8), protected:: TAVP_WRK_FLOPS_HEAVY=1d12  !minimal number of Flops to consider the operation as heavy-cost
        real(8), protected:: TAVP_WRK_FLOPS_MEDIUM=1d10 !minimal number of Flops to consider the operation as medium-cost
        real(8), protected:: TAVP_WRK_COST_TO_SIZE=1d2  !minimal cost (Flops) to size (Words) ratio to consider the operation arithmetically intensive
 !Bytecode:
        integer(INTL), parameter, private:: MAX_BYTECODE_SIZE=32_INTL*(1024_INTL*1024_INTL) !max size of an incoming/outgoing bytecode envelope (bytes)
        integer(INTD), parameter, private:: MAX_BYTECODE_INSTR=65536                        !max number of tensor instructions in a bytecode envelope
 !Resourcer:
        integer(INTD), private:: MAX_RESOURCER_INSTR=64  !max number of instructions during a single new resource allocation phase
        real(8), private:: MAX_RESOURCER_PHASE_TIME=1d-3 !max time spent in a single new resource allocation phase
 !Communicator:
        logical, private:: COMMUNICATOR_BLOCKING=.TRUE.         !switches between blocking and non-blocking communications
        integer(INTD), private:: MAX_COMMUNICATOR_PREFETCHES=16 !max number of outstanding prefetches issued by Communicator
        integer(INTD), private:: MAX_COMMUNICATOR_UPLOADS=8     !max number of outstanding uploads issued by Communicator
        real(8), private:: MAX_COMMUNICATOR_PHASE_TIME=1d-3     !max time spent by Communicator in each subphase
 !Dispatcher:
        integer(INTD), private:: MAX_DISPATCHER_INTAKE=64       !max number of instructions taken from the port at a time
!TYPES:
 !Tensor resource (local resource):
        type, extends(ds_resrc_t), private:: tens_resrc_t
         type(C_PTR), private:: base_addr=C_NULL_PTR   !C pointer to a local buffer for tensor body storage
         integer(C_SIZE_T), private:: bytes=0_C_SIZE_T !size of the tensor body storage buffer in bytes
         logical, private:: pinned=.FALSE.             !whether or not the buffer is pinned
         integer(C_INT), private:: dev_id=DEV_NULL     !flat device id where the buffer resides
         integer(C_INT), private:: ref_count=0         !reference count (how many tensor operands are associated with this resource)
         contains
          procedure, public:: is_empty=>TensResrcIsEmpty               !returns TRUE if the tensor resource is empty (unallocated)
          procedure, public:: allocate_buffer=>TensResrcAllocateBuffer !allocates a local buffer for tensor body storage
          procedure, public:: free_buffer=>TensResrcFreeBuffer         !frees the local buffer (at most one tensor operand can be associated with this resource at this time)
          procedure, public:: get_mem_ptr=>TensResrcGetMemPtr          !returns a C pointer to the local memory buffer
          procedure, public:: get_mem_size=>TensResrcGetMemSize        !returns the size of the memory buffer in bytes
          procedure, private:: incr_ref_count=>TensResrcIncrRefCount   !increments the reference count (number of tensor operands associated with the resource)
          procedure, private:: decr_ref_count=>TensResrcDecrRefCount   !decrements the reference count (number of tensor operands associated with the resource)
          procedure, private:: get_ref_count=>TensResrcGetRefCount     !returns the current reference count
#if !(defined(__GNUC__) && __GNUC__ < 8)
          final:: tens_resrc_dtor
#endif
        end type tens_resrc_t
 !Tensor argument cache entry (TAVP-specific):
        type, extends(tens_cache_entry_t), private:: tens_entry_wrk_t
         type(tens_resrc_t), private:: resource                                !tensor resource
         contains
          procedure, private:: TensEntryWrkCtor                                !ctor
          procedure, public:: tens_entry_wrk_ctor=>TensEntryWrkCtor
          procedure, public:: get_resource=>TensEntryWrkGetResource            !returns a non-owning pointer to the resource
          procedure, public:: set_tensor_layout=>TensEntryWrkSetTensorLayout   !sets the tensor layout, if not already set
          procedure, public:: acquire_resource=>TensEntryWrkAcquireResource    !acquires resource for the tensor cache entry
          procedure, public:: release_resource=>TensEntryWrkReleaseResource    !releases resource for the tensor cache entry
          final:: tens_entry_wrk_dtor                                          !dtor
        end type tens_entry_wrk_t
 !Reference to the tensor argument cache entry:
        type, private:: tens_entry_wrk_ref_t
         class(tens_entry_wrk_t), pointer, public:: cache_entry=>NULL() !non-owning pointer to a tensor cache entry
        end type tens_entry_wrk_ref_t
 !Tensor operand (encapsulated tensor data processible by a specific TAVP):
        type, extends(ds_oprnd_t), private:: tens_oprnd_t
         class(tens_entry_wrk_t), pointer, private:: cache_entry=>NULL() !non-owning pointer to a tensor cache entry where the tensor and tensor resource are stored (optional)
         class(tens_rcrsv_t), pointer, private:: tensor=>NULL()   !non-owning pointer to a persistent recursive tensor (normally stored in the tensor cache)
         class(tens_resrc_t), pointer, private:: resource=>NULL() !non-owning pointer to a persistent local tensor resource (normally stored in the tensor cache)
         type(talsh_tens_t), private:: talsh_tens                 !TAL-SH tensor object (for performing actual computations)
         contains
          procedure, private:: TensOprndCtorTensor                       !ctor by tensor only
          procedure, private:: TensOprndCtorCache                        !ctor by cache entry only
          generic, public:: tens_oprnd_ctor=>TensOprndCtorTensor,TensOprndCtorCache
          procedure, public:: get_tensor=>TensOprndGetTensor             !returns a non-owning pointer to the tensor
          procedure, public:: set_cache_entry=>TensOprndSetCacheEntry    !sets the associated tensor cache entry (may be NULL)
          procedure, public:: get_cache_entry=>TensOprndGetCacheEntry    !returns a non-owning pointer to the tensor cache entry (may be NULL)
          procedure, public:: set_resource=>TensOprndSetResource         !sets the resource component if it has not been set via the constructor
          procedure, public:: get_resource=>TensOprndGetResource         !returns a non-owning pointer to the tensor resource
          procedure, public:: set_talsh_tens=>TensOprndSetTalshTens      !sets up the TAL-SH tensor object for further processing with TAL-SH
          procedure, public:: set_tensor_layout=>TensOprndSetTensorLayout!sets the tensor layout, if not already set
          procedure, public:: register_read=>TensOprndRegisterRead       !registers a new read access on the tensor operand
          procedure, public:: unregister_read=>TensOprndUnregisterRead   !unregisters a read access on the tensor operand
          procedure, public:: get_read_count=>TensOprndGetReadCount      !returns the current read access count on the tensor operand
          procedure, public:: register_write=>TensOprndRegisterWrite     !registers a new write access on the tensor operand
          procedure, public:: unregister_write=>TensOprndUnregisterWrite !unregisters a write access on the tensor operand
          procedure, public:: get_write_count=>TensOprndGetWriteCount    !returns the current read access count on the tensor operand
          procedure, public:: is_located=>TensOprndIsLocated             !returns TRUE if the tensor operand has been located (its physical layout/location is known)
          procedure, public:: is_remote=>TensOprndIsRemote               !returns TRUE if the tensor operand is remote
          procedure, public:: is_valued=>TensOprndIsValued               !returns TRUE if the tensor operand is set to some value (neither undefined nor being updated)
          procedure, public:: has_resource=>TensOprndHasResource         !returns TRUE if the tensor operand has been allocated an actual local resource
          procedure, public:: acquire_rsc=>TensOprndAcquireRsc      !explicitly acquires local resources for the tensor operand
          procedure, public:: prefetch=>TensOprndPrefetch           !starts prefetching the remote tensor operand (acquires local resources!)
          procedure, public:: upload=>TensOprndUpload               !starts uploading the tensor operand to its remote location
          procedure, public:: sync=>TensOprndSync                   !synchronizes the currently pending communication on the tensor operand
          procedure, public:: release=>TensOprndRelease             !destroys the present local copy of the tensor operand (releases local resources!), but the operand stays defined
          procedure, public:: destruct=>TensOprndDestruct           !performs complete destruction back to an empty state
          procedure, public:: upload_local=>TensOprndUploadLocal    !performs a local accumulation of the tensor into another cached tensor
          procedure, public:: lock=>TensOprndLock                   !sets the lock for accessing/updating the tensor operand content
          procedure, public:: unlock=>TensOprndUnlock               !releases the access lock
          procedure, public:: print_it=>TensOprndPrintIt            !prints
          procedure, private:: reset_tmp_tensor=>TensOprndResetTmpTensor !resets the tensor in a tensor operand by providing another tensor cache entry with a temporary tensor (internal use only)
          final:: tens_oprnd_dtor                                   !dtor
        end type tens_oprnd_t
 !Tensor instruction (realization of a tensor operation for a specific TAVP):
        type, extends(ds_instr_t), public:: tens_instr_t
         integer(INTD), private:: num_out_oprnds=0                    !number of the output tensor instruction operands
         integer(INTD), private:: out_oprnds(0:MAX_TENSOR_OPERANDS-1) !positions of the tensor instruction operands which are considered output
         type(talsh_task_t), private:: talsh_task                     !TAL-SH task
         contains
          procedure, private:: TensInstrCtor                                 !ctor: constructs a tensor instruction from the specification of a tensor operation
          generic, public:: tens_instr_ctor=>TensInstrCtor
          procedure, public:: encode=>TensInstrEncode                        !encoding procedure: Packs the TAVP instruction into a raw byte packet (bytecode)
          procedure, public:: get_num_out_operands=>TensInstrGetNumOutOperands!returns the total number of output tensor operands in the tensor instruction (normally 1 or 0)
          procedure, public:: operand_is_output=>TensInstrOperandIsOutput    !returns TRUE if the specific tensor instruction operand is output, FALSE otherwise
          procedure, public:: get_output_operands=>TensInstrGetOutputOperands!returns the list of the output operands by their positions
          procedure, public:: lay_output_operands=>TensInstrLayOutputOperands!sets up storage layout for non-existing output tensor operands
          procedure, public:: get_cache_entries=>TensInstrGetCacheEntries    !returns an array of references to tensor cache entries used by the tensor operands
          procedure, public:: get_flops=>TensInstrGetFlops                   !returns an estimate of the total number of required Flops (mul/add) and memory Words
          procedure, public:: get_operation=>TensInstrGetOperation           !returns back the encapsulated tensor operation
          procedure, public:: mark_issue=>TensInstrMarkIssue                 !updates the tensor access counters for all tensor instruction operands due to the instruction issue
          procedure, public:: mark_completion=>TensInstrMarkCompletion       !updates the tensor access counters for all tensor instruction operands due to the instruction completion
          procedure, public:: dependency_free=>TensInstrDependencyFree       !returns TRUE if the tensor instruction is dependency-free
          procedure, public:: is_substitutable=>TensInstrIsSubstitutable     !returns TRUE if the tensor instruction enables output substitution (rename)
          procedure, public:: output_substituted=>TensInstrOutputSubstituted !returns TRUE if the output tensor(s) is(are) substituted with a temporary one(s)
          procedure, public:: print_it=>TensInstrPrintIt                     !prints
          final:: tens_instr_dtor                                            !dtor
        end type tens_instr_t
 !TAVP-WRK decoder:
        type, extends(ds_decoder_t), private:: tavp_wrk_decoder_t
         integer(INTD), public:: num_ports=1                        !number of ports: Port 0 <- Self
         integer(INTD), private:: source_comm                       !bytecode source communicator
         integer(INTD), private:: source_rank=-1                    !bytecode source process rank
         type(pack_env_t), private:: bytecode                       !incoming bytecode buffer
         class(tens_cache_t), pointer, private:: arg_cache=>NULL()  !non-owning pointer to the tensor argument cache
         type(list_bi_t), private:: control_list                    !list of control instructions
         type(list_iter_t), private:: ctrl_list                     !iterator for <control_list>
         contains
          procedure, public:: configure=>TAVPWRKDecoderConfigure    !configures TAVP-WRK decoder
          procedure, public:: start=>TAVPWRKDecoderStart            !starts TAVP-WRK decoder
          procedure, public:: shutdown=>TAVPWRKDecoderShutdown      !shuts down TAVP-WRK decoder
          procedure, public:: decode=>TAVPWRKDecoderDecode          !decodes the DS bytecode into DS instructions
        end type tavp_wrk_decoder_t
 !TAVP-WRK decoder configuration:
        type, extends(dsv_conf_t), private:: tavp_wrk_decoder_conf_t
         integer(INTD), public:: source_comm                  !MPI communicator of the source process
         integer(INTD), public:: source_rank                  !source process rank from which the bytecode is coming
         class(ds_unit_t), pointer, public:: acceptor=>NULL() !non-owning pointer to the acceptor DS unit for which the decoding is done
         integer(INTD), public:: acceptor_port_id             !associated acceptor port id
        end type tavp_wrk_decoder_conf_t
 !TAVP-WRK retirer:
        type, extends(ds_encoder_t), private:: tavp_wrk_retirer_t
         integer(INTD), public:: num_ports=1                        !number of ports: Port 0 <- Communicator (Tens) + Dispatcher (Ctrl,Aux)
         integer(INTD), private:: retire_comm                       !retired bytecode destination communicator
         integer(INTD), private:: retire_rank=-1                    !retired bytecode destination process rank
         type(pack_env_t), private:: bytecode                       !outgoing bytecode
         class(tens_cache_t), pointer, private:: arg_cache=>NULL()  !non-owning pointer to the tensor argument cache
         contains
          procedure, public:: configure=>TAVPWRKRetirerConfigure    !configures TAVP-WRK retirer
          procedure, public:: start=>TAVPWRKRetirerStart            !starts TAVP-WRK retirer
          procedure, public:: shutdown=>TAVPWRKRetirerShutdown      !shuts down TAVP-WRK retirer
          procedure, public:: encode=>TAVPWRKRetirerEncode          !encodes a DS instruction into the DS bytecode
        end type tavp_wrk_retirer_t
 !TAVP-WRK retirer configuration:
        type, extends(dsv_conf_t), private:: tavp_wrk_retirer_conf_t
         integer(INTD), public:: retire_comm                        !MPI communicator of the retired bytecode destination process
         integer(INTD), public:: retire_rank                        !destination process rank to which the retired bytecode is going
        end type tavp_wrk_retirer_conf_t
 !TAVP-WRK resourcer:
        type, extends(ds_unit_t), private:: tavp_wrk_resourcer_t
         integer(INTD), public:: num_ports=2                        !number of ports: Port 0 <- Decoder (Tens,Ctrl,Aux), Port 1 <- Communicator (Tens)
         integer(INTL), private:: host_ram_size=0_INTL              !size of the usable Host RAM memory in bytes
         integer(INTL), private:: nvram_size=0_INTL                 !size of the usable NVRAM memory (if any) in bytes
         class(tens_cache_t), pointer, private:: arg_cache=>NULL()  !non-owning pointer to the tensor argument cache
         type(list_bi_t), private:: staged_list                     !list of staged instructions ready for subsequent processing
         type(list_iter_t), private:: stg_list                      !iterator for <staged_list>
         type(list_bi_t), private:: deferred_list                   !list of deferred instructions
         type(list_iter_t), private:: def_list                      !iterator for <deferred_list>
         type(list_bi_t), private:: release_list                    !list of completed tensor instructions expecting resource release
         type(list_iter_t), private:: rls_list                      !iterator for <release_list>
         contains
          procedure, public:: configure=>TAVPWRKResourcerConfigure                !configures TAVP-WRK resourcer
          procedure, public:: start=>TAVPWRKResourcerStart                        !starts TAVP-WRK resourcer
          procedure, public:: shutdown=>TAVPWRKResourcerShutdown                  !shuts down TAVP-WRK resourcer
          procedure, public:: acquire_resources=>TAVPWRKResourcerAcquireResources !acquires local resources for a tensor instruction
          procedure, public:: release_resources=>TAVPWRKResourcerReleaseResources !releases local resources from a tensor instruction
          procedure, public:: substitute_output=>TAVPWRKResourcerSubstituteOutput !substitutes the persistent output tensor with a temporary one
          procedure, public:: restore_output=>TAVPWRKResourcerRestoreOutput       !restores back the original (persistent) output tensor
        end type tavp_wrk_resourcer_t
 !TAVP-WRK resourcer configuration:
        type, extends(dsv_conf_t), private:: tavp_wrk_resourcer_conf_t
         integer(INTL), public:: host_ram_size !size of the usable Host RAM memory in bytes
         integer(INTL), public:: nvram_size    !size of the usable NVRAM memory (if any) in bytes
        end type tavp_wrk_resourcer_conf_t
 !TAVP-WRK communicator:
        type, extends(ds_unit_t), private:: tavp_wrk_communicator_t
         integer(INTD), public:: num_ports=2                        !number of ports: Port 0 <- Resourcer (Tens,Ctrl,Aux), Port 1 <- Dispatcher (Tens)
         integer(INTD), private:: num_mpi_windows=TAVP_WRK_NUM_WINS !number of dynamic MPI windows per global addressing space
         class(DistrSpace_t), pointer, private:: addr_space=>NULL() !non-owning pointer to the global address space
         class(tens_cache_t), pointer, private:: arg_cache=>NULL()  !non-owning pointer to the tensor argument cache
         type(list_bi_t), private:: prefetch_list                   !list of tensor instructions undergoing input prefetch
         type(list_iter_t), private:: fet_list                      !iterator for <prefetch_list>
         type(list_bi_t), private:: upload_list                     !list of tensor instructions undergoing output upload
         type(list_iter_t), private:: upl_list                      !iterator for <upload_list>
         type(list_bi_t), private:: dispatch_list                   !list of tensor instructions ready to be dispatched
         type(list_iter_t), private:: dsp_list                      !iterator for <dispatch_list>
         type(list_bi_t), private:: retire_list                     !list of tensor instructions ready to be retired
         type(list_iter_t), private:: ret_list                      !iterator for <retire_list>
         contains
          procedure, public:: configure=>TAVPWRKCommunicatorConfigure          !configures TAVP-WRK communicator
          procedure, public:: start=>TAVPWRKCommunicatorStart                  !starts TAVP-WRK communicator
          procedure, public:: shutdown=>TAVPWRKCommunicatorShutdown            !shuts down TAVP-WRK communicator
          procedure, public:: prefetch_input=>TAVPWRKCommunicatorPrefetchInput !starts prefetching input arguments
          procedure, public:: sync_prefetch=>TAVPWRKCommunicatorSyncPrefetch   !synchronizes on the input prefetch
          procedure, public:: upload_output=>TAVPWRKCommunicatorUploadOutput   !starts uploading the output argument
          procedure, public:: sync_upload=>TAVPWRKCommunicatorSyncUpload       !synchronizes on the output upload
        end type tavp_wrk_communicator_t
 !TAVP-WRK communicator configuration:
        type, extends(dsv_conf_t), private:: tavp_wrk_communicator_conf_t
         integer(INTD), public:: num_mpi_windows                               !number of dynamic MPI windows per global addressing space
        end type tavp_wrk_communicator_conf_t
 !TAVP-WRK dispatcher instruction execution procedure:
        type, private:: tavp_wrk_dispatch_proc_t
         procedure(tavp_wrk_dispatch_proc_i), nopass, pointer, public:: instr_proc=>NULL() !procedure pointer to the instruction execution procedure
        end type tavp_wrk_dispatch_proc_t
 !TAVP-WRK dispatcher:
        type, extends(ds_unit_t), private:: tavp_wrk_dispatcher_t
         integer(INTD), public:: num_ports=1                                    !number of ports: Port 0 <- Communicator (Tens,Ctrl,Aux)
         integer(INTL), private:: host_buf_size=TAVP_WRK_HOST_BUF_SIZE          !size of the pinned Host argument buffer
         integer(INTD), private:: host_arg_max=0                                !max number of tensor arguments in the pinned Host argument buffer
         integer(INTD), allocatable, private:: gpu_list(:)                      !list of available NVIDIA GPU
         integer(INTD), allocatable, private:: amd_list(:)                      !list of available AMD GPU
         integer(INTD), allocatable, private:: mic_list(:)                      !list of available INTEL MIC
         type(tavp_wrk_dispatch_proc_t), private:: microcode(0:TAVP_ISA_SIZE-1) !instruction execution microcode bindings
         class(tens_cache_t), pointer, private:: arg_cache=>NULL()              !non-owning pointer to the tensor argument cache
         type(list_bi_t), private:: issued_list                                 !list of issued tensor instructions
         type(list_iter_t), private:: iss_list                                  !iterator for <issued_list>
         type(list_bi_t), private:: completed_list                              !list of (locally) completed tensor instructions
         type(list_iter_t), private:: cml_list                                  !iterator for <completed_list>
         contains
          procedure, public:: configure=>TAVPWRKDispatcherConfigure     !configures TAVP-WRK dispatcher
          procedure, public:: start=>TAVPWRKDispatcherStart             !starts TAVP-WRK dispatcher
          procedure, public:: shutdown=>TAVPWRKDispatcherShutdown       !shuts down TAVP-WRK dispatcher
          procedure, public:: issue_instr=>TAVPWRKDispatcherIssueInstr  !issues a tensor instruction to a specific computing device
          procedure, public:: sync_instr=>TAVPWRKDispatcherSyncInstr    !synchronizes on tensor instruction execution, either TEST or WAIT
        end type tavp_wrk_dispatcher_t
 !TAVP-WRK dispatcher configuration:
        type, extends(dsv_conf_t), private:: tavp_wrk_dispatcher_conf_t
         integer(INTL), public:: host_buf_size                          !size of the pinned Host argument buffer
         integer(INTD), allocatable, public:: gpu_list(:)               !list of available NVIDIA GPU
         integer(INTD), allocatable, public:: amd_list(:)               !list of available AMD GPU
         integer(INTD), allocatable, public:: mic_list(:)               !list of available INTEL MIC
        end type tavp_wrk_dispatcher_conf_t
 !TAVP-WRK:
        type, extends(dsvp_t), public:: tavp_wrk_t
         type(tens_cache_t), private:: tens_cache                 !tensor argument cache (SHARED RESOURCE!)
         type(DistrSpace_t), private:: addr_space                 !global (distributed) address space
         type(tavp_wrk_decoder_t), private:: decoder              !DSVU: decodes incoming tensor instructions from the manager
         type(tavp_wrk_retirer_t), private:: retirer              !DSVU: retires processed tensor instructions and sends them back to the manager
         type(tavp_wrk_resourcer_t), private:: resourcer          !DSVU: allocates local resources for tensor instructions (only for defined tensor operands)
         type(tavp_wrk_communicator_t), private:: communicator    !DSVU: fetches/uploads remote tensor operands
         type(tavp_wrk_dispatcher_t), private:: dispatcher        !DSVU: dispatches and executes ready to be executed tensor instructions to compute devices (runs the microcode)
         contains
          procedure, public:: configure=>TAVPWRKConfigure         !configures the TAVP-WRK DSVP
        end type tavp_wrk_t
 !TAVP-WRK configuration:
        type, extends(dsv_conf_t), public:: tavp_wrk_conf_t
         character(:), allocatable, public:: description    !TAVP description
         integer(INTD), public:: tavp_id                    !TAVP id
         integer(INTD), public:: source_comm                !MPI communicator of the bytecode source
         integer(INTD), public:: source_rank                !MPI process rank of the bytecode source
         integer(INTD), public:: retire_comm                !MPI communicator of the retired bytecode destination
         integer(INTD), public:: retire_rank                !MPI process rank of the retired bytecode destination
         integer(INTL), public:: host_ram_size              !size of the usable Host RAM memory in bytes
         integer(INTL), public:: nvram_size                 !size of the usable NVRAM memory (if any) in bytes
         integer(INTD), public:: num_mpi_windows            !number of dynamic MPI windows per global addressing space
         integer(INTL), public:: host_buf_size              !pinned Host argument buffer size in bytes
         integer(INTD), allocatable, public:: gpu_list(:)   !list of the accesible NVIDIA GPU devices
         integer(INTD), allocatable, public:: amd_list(:)   !list of the accesible AMD GPU devices
         integer(INTD), allocatable, public:: mic_list(:)   !list of the accesible Intel MIC devices
        end type tavp_wrk_conf_t
!INTERFACES:
        abstract interface
 !tavp_wrk_dispatcher_t:
         subroutine tavp_wrk_dispatch_proc_i(this,tens_instr,ierr,dev_id)
          import:: tavp_wrk_dispatcher_t,tens_instr_t,INTD
          class(tavp_wrk_dispatcher_t), intent(inout):: this !inout: TAVP-WRK dispatcher DSVU
          class(tens_instr_t), intent(inout):: tens_instr    !inout: tensor instruction
          integer(INTD), intent(out), optional:: ierr        !out: error code
          integer(INTD), intent(in), optional:: dev_id       !in: flat device id
         end subroutine tavp_wrk_dispatch_proc_i
        end interface
!VISIBILITY:
 !non-member test/debug:
        private test_carma
        public tavp_wrk_reset_output
 !tens_resrc_t:
        private TensResrcIsEmpty
        private TensResrcAllocateBuffer
        private TensResrcFreeBuffer
        private TensResrcGetMemPtr
        private TensResrcGetMemSize
        private TensResrcIncrRefCount
        private TensResrcDecrRefCount
        private TensResrcGetRefCount
        public tens_resrc_dtor
 !tens_entry_wrk_t:
        private TensEntryWrkCtor
        private TensEntryWrkGetResource
        private TensEntryWrkSetTensorLayout
        private TensEntryWrkAcquireResource
        private TensEntryWrkReleaseResource
        public tens_entry_wrk_dtor
        public tens_entry_wrk_alloc
 !tens_oprnd_t:
        private TensOprndCtorTensor
        private TensOprndCtorCache
        private TensOprndGetTensor
        private TensOprndSetCacheEntry
        private TensOprndGetCacheEntry
        private TensOprndSetResource
        private TensOprndGetResource
        private TensOprndSetTalshTens
        private TensOprndSetTensorLayout
        private TensOprndRegisterRead
        private TensOprndUnregisterRead
        private TensOprndGetReadCount
        private TensOprndRegisterWrite
        private TensOprndUnregisterWrite
        private TensOprndGetWriteCount
        private TensOprndIsLocated
        private TensOprndIsRemote
        private TensOprndIsValued
        private TensOprndHasResource
        private TensOprndAcquireRsc
        private TensOprndPrefetch
        private TensOprndUpload
        private TensOprndSync
        private TensOprndRelease
        private TensOprndDestruct
        private TensOprndUploadLocal
        private TensOprndLock
        private TensOprndUnlock
        private TensOprndPrintIt
        private TensOprndResetTmpTensor
        public tens_oprnd_dtor
 !tens_instr_t:
        private TensInstrCtor
        private TensInstrEncode
        private TensInstrGetNumOutOperands
        private TensInstrOperandIsOutput
        private TensInstrGetOutputOperands
        private TensInstrLayOutputOperands
        private TensInstrGetCacheEntries
        private TensInstrGetFlops
        private TensInstrGetOperation
        private TensInstrMarkIssue
        private TensInstrMarkCompletion
        private TensInstrDependencyFree
        private TensInstrIsSubstitutable
        private TensInstrOutputSubstituted
        private TensInstrPrintIt
        public tens_instr_dtor
        private tens_instr_print
 !tavp_wrk_decoder_t:
        private TAVPWRKDecoderConfigure
        private TAVPWRKDecoderStart
        private TAVPWRKDecoderShutdown
        private TAVPWRKDecoderDecode
 !tavp_wrk_retirer_t:
        private TAVPWRKRetirerConfigure
        private TAVPWRKRetirerStart
        private TAVPWRKRetirerShutdown
        private TAVPWRKRetirerEncode
 !tavp_wrk_resourcer_t:
        private TAVPWRKResourcerConfigure
        private TAVPWRKResourcerStart
        private TAVPWRKResourcerShutdown
        private TAVPWRKResourcerAcquireResources
        private TAVPWRKResourcerReleaseResources
        private TAVPWRKResourcerSubstituteOutput
        private TAVPWRKResourcerRestoreOutput
 !tavp_wrk_communicator_t:
        private TAVPWRKCommunicatorConfigure
        private TAVPWRKCommunicatorStart
        private TAVPWRKCommunicatorShutdown
        private TAVPWRKCommunicatorPrefetchInput
        private TAVPWRKCommunicatorSyncPrefetch
        private TAVPWRKCommunicatorUploadOutput
        private TAVPWRKCommunicatorSyncUpload
 !tavp_wrk_dispatcher_t:
        private TAVPWRKDispatcherConfigure
        private TAVPWRKDispatcherStart
        private TAVPWRKDispatcherShutdown
        private TAVPWRKDispatcherIssueInstr
        private TAVPWRKDispatcherSyncInstr
        private TAVPWRKExecTensorCreate
        private TAVPWRKExecTensorDestroy
        private TAVPWRKExecTensorInit
        private TAVPWRKExecTensorContract
        private tavp_wrk_dispatch_proc_i
 !tavp_wrk_t:
        private TAVPWRKConfigure
!IMPLEMENTATION:
       contains
![non-member:Test/Debug]===========
        subroutine test_carma(ierr)
!DEBUG: Brute-force implementation of a single distributed tensor contraction.
         implicit none
         integer(INTD), intent(out):: ierr
         integer(INTD), parameter:: TENS_RANK=4
         integer(INTD), parameter:: DIM_SEG_SIZE=8
         integer(INTD), parameter:: DIM_NUM_LEVELS=4
         integer(INTD), parameter:: NUM_DIM_SEGS=2**(DIM_NUM_LEVELS-1)
         integer(INTD), parameter:: SPLIT_BASE=2**TENS_RANK
         integer(INTD), parameter:: TOTAL_BLOCKS=SPLIT_BASE**(DIM_NUM_LEVELS-1)
         integer(INTD), parameter:: BLOCK_VOL=DIM_SEG_SIZE**TENS_RANK
         integer(INTD):: num_procs,num_blocks,i,j,k,l,n,lid,rid
         integer(INTD):: dsg(1:TENS_RANK),lsg(1:TENS_RANK),rsg(1:TENS_RANK)
         integer(INTL):: tg
         type(DistrSpace_t):: tavp_addr_space
         type(DataDescr_t):: dd,ld,rd
         type(DataDescr_t), allocatable:: ddes(:),ldes(:),rdes(:)
         type(DataDescr_t), allocatable:: ddesa(:),ldesa(:),rdesa(:) !root only
         real(8), allocatable, target:: dtens(:),ltens(:),rtens(:)
         type(pack_env_t):: packenv
         type(obj_pack_t):: packet
         type(comm_handle_t):: chl
         real(8), pointer, contiguous:: block_p(:)
         type(C_PTR):: mem_p,dmem_p,lmem_p,rmem_p
         type(talsh_tens_t):: dtns,ltns,rtns
         type(talsh_task_t):: tsk0
         real(8):: tms,tm,tcs,tc,tts,tt,tas,ta

         ierr=0
         num_procs=role_size
         if(ierr.ne.0) call quit(-1,'Bad CARMA!')
         if(mod(TOTAL_BLOCKS,num_procs).ne.0) call quit(-2,'Bad CARMA!')
         num_blocks=TOTAL_BLOCKS/num_procs !number of tensor blocks per process
         call tavp_addr_space%create(role_comm,TAVP_WRK_NUM_WINS,'WorkAddressSpace',ierr)
         if(ierr.ne.0) call quit(-3,'Bad CARMA!')
!Create <num_blocks> tensor blocks on each process:
 !Destination tensor:
         allocate(dtens(0:BLOCK_VOL*num_blocks-1)); allocate(ddes(0:num_blocks-1))
         dtens(:)=0d0
         do i=0,num_blocks-1
          block_p(0:)=>dtens(BLOCK_VOL*i:BLOCK_VOL*(i+1)-1); mem_p=c_loc(block_p)
          call tavp_addr_space%attach(mem_p,R8,int(BLOCK_VOL,8),ddes(i),ierr)
          if(ierr.ne.0) call quit(-3,'Bad CARMA!')
         enddo
 !Left tensor:
         allocate(ltens(0:BLOCK_VOL*num_blocks-1)); allocate(ldes(0:num_blocks-1))
         ltens(:)=1d-3
         do i=0,num_blocks-1
          block_p(0:)=>ltens(BLOCK_VOL*i:BLOCK_VOL*(i+1)-1); mem_p=c_loc(block_p)
          call tavp_addr_space%attach(mem_p,R8,int(BLOCK_VOL,8),ldes(i),ierr)
          if(ierr.ne.0) call quit(-4,'Bad CARMA!')
         enddo
 !Right tensor:
         allocate(rtens(0:BLOCK_VOL*num_blocks-1)); allocate(rdes(0:num_blocks-1))
         rtens(:)=1d-4
         do i=0,num_blocks-1
          block_p(0:)=>rtens(BLOCK_VOL*i:BLOCK_VOL*(i+1)-1); mem_p=c_loc(block_p)
          call tavp_addr_space%attach(mem_p,R8,int(BLOCK_VOL,8),rdes(i),ierr)
          if(ierr.ne.0) call quit(-5,'Bad CARMA!')
         enddo
!Root collects data descriptors from all other processes:
 !Destination tensor:
         if(role_rank.eq.0) then
          allocate(ddesa(0:TOTAL_BLOCKS-1))
          do i=1,num_procs-1 !remote descriptors
           call packenv%reserve_mem(ierr); if(ierr.ne.0) call quit(-6,'Bad CARMA!')
           do while(.not.packenv%receive(chl,ierr,comm=role_comm))
            if(ierr.ne.0) call quit(-7,'Bad CARMA!')
           enddo
           if(ierr.ne.0) call quit(-8,'Bad CARMA!')
           call chl%wait(ierr); if(ierr.ne.0) call quit(-9,'Bad CARMA!')
           n=packenv%get_num_packets(ierr); if(ierr.ne.0) call quit(-10,'Bad CARMA!')
           do j=1,n
            call packenv%extract_packet(j,packet,ierr,tg,preclean=.TRUE.); if(ierr.ne.0) call quit(-11,'Bad CARMA!')
            call ddesa(tg)%unpack(packet,ierr); if(ierr.ne.0) call quit(-12,'Bad CARMA!')
            call packet%clean(ierr); if(ierr.ne.0) call quit(-13,'Bad CARMA!')
           enddo
           call chl%clean(ierr); if(ierr.ne.0) call quit(-14,'Bad CARMA!')
           call packenv%destroy(ierr); if(ierr.ne.0) call quit(-15,'Bad CARMA!')
          enddo
          do i=0,num_blocks-1 !local descriptors
           ddesa(i)=ddes(i)
          enddo
         else
          call packenv%reserve_mem(ierr); if(ierr.ne.0) call quit(-16,'Bad CARMA!')
          do i=0,num_blocks-1
           call packenv%acquire_packet(packet,ierr); if(ierr.ne.0) call quit(-17,'Bad CARMA!')
           call ddes(i)%pack(packet,ierr); if(ierr.ne.0) call quit(-18,'Bad CARMA!')
           call packenv%seal_packet(ierr,tag=int(role_rank*num_blocks+i,INTL)); if(ierr.ne.0) call quit(-19,'Bad CARMA!')
          enddo
          call packenv%send(0,chl,ierr,comm=role_comm); if(ierr.ne.0) call quit(-20,'Bad CARMA!')
          call chl%wait(ierr); if(ierr.ne.0) call quit(-21,'Bad CARMA!')
          call chl%clean(ierr); if(ierr.ne.0) call quit(-22,'Bad CARMA!')
          call packenv%destroy(ierr); if(ierr.ne.0) call quit(-23,'Bad CARMA!')
         endif
         call role_barrier()
 !Left tensor:
         if(role_rank.eq.0) then
          allocate(ldesa(0:TOTAL_BLOCKS-1))
          do i=1,num_procs-1 !remote descriptors
           call packenv%reserve_mem(ierr); if(ierr.ne.0) call quit(-24,'Bad CARMA!')
           do while(.not.packenv%receive(chl,ierr,comm=role_comm))
            if(ierr.ne.0) call quit(-25,'Bad CARMA!')
           enddo
           if(ierr.ne.0) call quit(-26,'Bad CARMA!')
           call chl%wait(ierr); if(ierr.ne.0) call quit(-27,'Bad CARMA!')
           n=packenv%get_num_packets(ierr); if(ierr.ne.0) call quit(-28,'Bad CARMA!')
           do j=1,n
            call packenv%extract_packet(j,packet,ierr,tg,preclean=.TRUE.); if(ierr.ne.0) call quit(-29,'Bad CARMA!')
            call ldesa(tg)%unpack(packet,ierr); if(ierr.ne.0) call quit(-30,'Bad CARMA!')
            call packet%clean(ierr); if(ierr.ne.0) call quit(-31,'Bad CARMA!')
           enddo
           call chl%clean(ierr); if(ierr.ne.0) call quit(-32,'Bad CARMA!')
           call packenv%destroy(ierr); if(ierr.ne.0) call quit(-33,'Bad CARMA!')
          enddo
          do i=0,num_blocks-1 !local descriptors
           ldesa(i)=ldes(i)
          enddo
         else
          call packenv%reserve_mem(ierr); if(ierr.ne.0) call quit(-34,'Bad CARMA!')
          do i=0,num_blocks-1
           call packenv%acquire_packet(packet,ierr); if(ierr.ne.0) call quit(-35,'Bad CARMA!')
           call ldes(i)%pack(packet,ierr); if(ierr.ne.0) call quit(-36,'Bad CARMA!')
           call packenv%seal_packet(ierr,tag=int(role_rank*num_blocks+i,INTL)); if(ierr.ne.0) call quit(-37,'Bad CARMA!')
          enddo
          call packenv%send(0,chl,ierr,comm=role_comm); if(ierr.ne.0) call quit(-38,'Bad CARMA!')
          call chl%wait(ierr); if(ierr.ne.0) call quit(-39,'Bad CARMA!')
          call chl%clean(ierr); if(ierr.ne.0) call quit(-40,'Bad CARMA!')
          call packenv%destroy(ierr); if(ierr.ne.0) call quit(-41,'Bad CARMA!')
         endif
         call role_barrier()
 !Right tensor:
         if(role_rank.eq.0) then
          allocate(rdesa(0:TOTAL_BLOCKS-1))
          do i=1,num_procs-1 !remote descriptors
           call packenv%reserve_mem(ierr); if(ierr.ne.0) call quit(-42,'Bad CARMA!')
           do while(.not.packenv%receive(chl,ierr,comm=role_comm))
            if(ierr.ne.0) call quit(-43,'Bad CARMA!')
           enddo
           if(ierr.ne.0) call quit(-44,'Bad CARMA!')
           call chl%wait(ierr); if(ierr.ne.0) call quit(-45,'Bad CARMA!')
           n=packenv%get_num_packets(ierr); if(ierr.ne.0) call quit(-46,'Bad CARMA!')
           do j=1,n
            call packenv%extract_packet(j,packet,ierr,tg,preclean=.TRUE.); if(ierr.ne.0) call quit(-47,'Bad CARMA!')
            call rdesa(tg)%unpack(packet,ierr); if(ierr.ne.0) call quit(-48,'Bad CARMA!')
            call packet%clean(ierr); if(ierr.ne.0) call quit(-49,'Bad CARMA!')
           enddo
           call chl%clean(ierr); if(ierr.ne.0) call quit(-50,'Bad CARMA!')
           call packenv%destroy(ierr); if(ierr.ne.0) call quit(-51,'Bad CARMA!')
          enddo
          do i=0,num_blocks-1 !local descriptors
           rdesa(i)=rdes(i)
          enddo
         else
          call packenv%reserve_mem(ierr); if(ierr.ne.0) call quit(-52,'Bad CARMA!')
          do i=0,num_blocks-1
           call packenv%acquire_packet(packet,ierr); if(ierr.ne.0) call quit(-53,'Bad CARMA!')
           call rdes(i)%pack(packet,ierr); if(ierr.ne.0) call quit(-54,'Bad CARMA!')
           call packenv%seal_packet(ierr,tag=int(role_rank*num_blocks+i,INTL)); if(ierr.ne.0) call quit(-55,'Bad CARMA!')
          enddo
          call packenv%send(0,chl,ierr,comm=role_comm); if(ierr.ne.0) call quit(-56,'Bad CARMA!')
          call chl%wait(ierr); if(ierr.ne.0) call quit(-57,'Bad CARMA!')
          call chl%clean(ierr); if(ierr.ne.0) call quit(-58,'Bad CARMA!')
          call packenv%destroy(ierr); if(ierr.ne.0) call quit(-59,'Bad CARMA!')
         endif
         call role_barrier()
!Root creates and sends tasks to all processes:
         if(role_rank.eq.0) then
          do i=1,num_procs-1 !process rank
           call packenv%reserve_mem(ierr,mem_size=16*1048576_INTL,max_packets=65536); if(ierr.ne.0) call quit(-60,'Bad CARMA!')
           do j=i*num_blocks,(i+1)*num_blocks-1 !destination block flat id
            do k=0,NUM_DIM_SEGS-1
             do l=0,NUM_DIM_SEGS-1
              call flat2signa(j,dsg)
              lsg(1:TENS_RANK)=(/dsg(4),k,dsg(2),l/); lid=signa2flat(lsg)
              rsg(1:TENS_RANK)=(/dsg(3),l,dsg(1),k/); rid=signa2flat(rsg)
              call packenv%acquire_packet(packet,ierr); if(ierr.ne.0) call quit(-61,'Bad CARMA!')
              call ddesa(j)%pack(packet,ierr); if(ierr.ne.0) call quit(-62,'Bad CARMA!')
              call ldesa(lid)%pack(packet,ierr); if(ierr.ne.0) call quit(-63,'Bad CARMA!')
              call rdesa(rid)%pack(packet,ierr); if(ierr.ne.0) call quit(-64,'Bad CARMA!')
              call packenv%seal_packet(ierr); if(ierr.ne.0) then; print *,ierr; call quit(-65,'Bad CARMA!'); endif
             enddo
            enddo
           enddo
           call packenv%send(i,chl,ierr,comm=role_comm); if(ierr.ne.0) call quit(-66,'Bad CARMA!')
           call chl%wait(ierr); if(ierr.ne.0) call quit(-67,'Bad CARMA!')
           call chl%clean(ierr); if(ierr.ne.0) call quit(-68,'Bad CARMA!')
           call packenv%destroy(ierr); if(ierr.ne.0) call quit(-69,'Bad CARMA!')
          enddo
          call packenv%reserve_mem(ierr,mem_size=16*1048576_INTL,max_packets=65536); if(ierr.ne.0) call quit(-70,'Bad CARMA!')
          do j=0,num_blocks-1
           do k=0,NUM_DIM_SEGS-1
            do l=0,NUM_DIM_SEGS-1
             call flat2signa(j,dsg)
             lsg(1:TENS_RANK)=(/dsg(4),k,dsg(2),l/); lid=signa2flat(lsg)
             rsg(1:TENS_RANK)=(/dsg(3),l,dsg(1),k/); rid=signa2flat(rsg)
             call packenv%acquire_packet(packet,ierr); if(ierr.ne.0) call quit(-71,'Bad CARMA!')
             call ddesa(j)%pack(packet,ierr); if(ierr.ne.0) call quit(-72,'Bad CARMA!')
             call ldesa(lid)%pack(packet,ierr); if(ierr.ne.0) call quit(-73,'Bad CARMA!')
             call rdesa(rid)%pack(packet,ierr); if(ierr.ne.0) call quit(-74,'Bad CARMA!')
             call packenv%seal_packet(ierr); if(ierr.ne.0) then; print *,ierr; call quit(-75,'Bad CARMA!'); endif
            enddo
           enddo
          enddo
         else
          call packenv%reserve_mem(ierr,mem_size=16*1048576_INTL,max_packets=65536); if(ierr.ne.0) call quit(-76,'Bad CARMA!')
          do while(.not.packenv%receive(chl,ierr,comm=role_comm))
           if(ierr.ne.0) call quit(-77,'Bad CARMA!')
          enddo
          if(ierr.ne.0) call quit(-78,'Bad CARMA!')
          call chl%wait(ierr); if(ierr.ne.0) call quit(-79,'Bad CARMA!')
          call chl%clean(ierr); if(ierr.ne.0) call quit(-80,'Bad CARMA!')
         endif
         call role_barrier()
!All processes synchronize and execute their tasks:
         n=packenv%get_num_packets(ierr); if(ierr.ne.0) call quit(-81,'Bad CARMA!'); print *,role_rank,n
         ierr=mem_allocate(talsh_flat_dev_id(DEV_HOST,0),int(BLOCK_VOL*8,C_SIZE_T),YEP,lmem_p)
         if(ierr.ne.0) call quit(-82,'Bad CARMA!')
         ierr=mem_allocate(talsh_flat_dev_id(DEV_HOST,0),int(BLOCK_VOL*8,C_SIZE_T),YEP,rmem_p)
         if(ierr.ne.0) call quit(-83,'Bad CARMA!')
         call role_barrier()
         tc=0d0; tt=0d0; ta=0d0; tms=thread_wtime()
         do i=1,n
          j=mod((i-1)+role_rank,n)+1
          call packenv%extract_packet(j,packet,ierr,preclean=.TRUE.); if(ierr.ne.0) call quit(-84,'Bad CARMA!')
          call dd%unpack(packet,ierr); if(ierr.ne.0) call quit(-85,'Bad CARMA!')
          call ld%unpack(packet,ierr); if(ierr.ne.0) call quit(-86,'Bad CARMA!')
          call rd%unpack(packet,ierr); if(ierr.ne.0) call quit(-87,'Bad CARMA!')
          tts=thread_wtime()
          call ld%get_data(lmem_p,ierr); if(ierr.ne.0) call quit(-88,'Bad CARMA!')
          call rd%get_data(rmem_p,ierr); if(ierr.ne.0) call quit(-89,'Bad CARMA!')
          tt=tt+thread_wtime(tts)
          tas=thread_wtime()
          dmem_p=dd%get_data_ptr(ierr); if(ierr.ne.0) call quit(-90,'Bad CARMA!')
          ierr=talsh_tensor_construct(dtns,R8,(/DIM_SEG_SIZE,DIM_SEG_SIZE,DIM_SEG_SIZE,DIM_SEG_SIZE/),&
                                     &talsh_flat_dev_id(DEV_HOST,0),dmem_p)
          ierr=talsh_tensor_construct(ltns,R8,(/DIM_SEG_SIZE,DIM_SEG_SIZE,DIM_SEG_SIZE,DIM_SEG_SIZE/),&
                                     &talsh_flat_dev_id(DEV_HOST,0),lmem_p)
          if(ierr.ne.0) call quit(-91,'Bad CARMA!')
          ierr=talsh_tensor_construct(rtns,R8,(/DIM_SEG_SIZE,DIM_SEG_SIZE,DIM_SEG_SIZE,DIM_SEG_SIZE/),&
                                     &talsh_flat_dev_id(DEV_HOST,0),rmem_p)
          if(ierr.ne.0) call quit(-92,'Bad CARMA!')
          ta=ta+thread_wtime(tas)
!          tts=thread_wtime()
!          call ld%flush_data(ierr,local=.FALSE.); if(ierr.ne.0) call quit(-93,'Bad CARMA!')
!          call rd%flush_data(ierr,local=.FALSE.); if(ierr.ne.0) call quit(-94,'Bad CARMA!')
!          tt=tt+thread_wtime(tts)
          tcs=thread_wtime()
          ierr=talsh_tensor_contract('D(a,b,c,d)+=L(d,i,b,j)*R(c,j,a,i)',dtns,ltns,rtns,dev_id=0,dev_kind=DEV_HOST,&
                                    &copy_ctrl=COPY_TTT)
          if(ierr.ne.0) call quit(-95,'Bad CARMA!')
          tc=tc+thread_wtime(tcs)
          tas=thread_wtime()
          ierr=talsh_tensor_destruct(rtns); if(ierr.ne.0) call quit(-96,'Bad CARMA!')
          ierr=talsh_tensor_destruct(ltns); if(ierr.ne.0) call quit(-97,'Bad CARMA!')
          ierr=talsh_tensor_destruct(dtns); if(ierr.ne.0) call quit(-98,'Bad CARMA!')
          ta=ta+thread_wtime(tas)
         enddo
         tm=thread_wtime(tms); print *,'Rank ',role_rank,': Total ',tm,' s: Contract ',tc,' s: Comm ',tt,' s: Alloc ',ta,' s'
         call role_barrier()
         ierr=mem_free(talsh_flat_dev_id(DEV_HOST,0),rmem_p); if(ierr.ne.0) call quit(-99,'Bad CARMA!')
         ierr=mem_free(talsh_flat_dev_id(DEV_HOST,0),lmem_p); if(ierr.ne.0) call quit(-100,'Bad CARMA!')
         call packenv%destroy(ierr); if(ierr.ne.0) call quit(-101,'Bad CARMA!')
!Detach tensor blocks:
 !Right tensor:
         do i=0,num_blocks-1
          call tavp_addr_space%detach(rdes(i),ierr)
          if(ierr.ne.0) call quit(-102,'Bad CARMA!')
         enddo
         if(allocated(rdesa)) deallocate(rdesa)
         deallocate(rdes); deallocate(rtens)
 !Left tensor:
         do i=0,num_blocks-1
          call tavp_addr_space%detach(ldes(i),ierr)
          if(ierr.ne.0) call quit(-103,'Bad CARMA!')
         enddo
         if(allocated(ldesa)) deallocate(ldesa)
         deallocate(ldes); deallocate(ltens)
 !Destination tensor:
         do i=0,num_blocks-1
          call tavp_addr_space%detach(ddes(i),ierr)
          if(ierr.ne.0) call quit(-104,'Bad CARMA!')
         enddo
         if(allocated(ddesa)) deallocate(ddesa)
         deallocate(ddes); deallocate(dtens)
         call tavp_addr_space%destroy(ierr)
         if(ierr.ne.0) call quit(-105,'Bad CARMA!')
         return

        contains

         subroutine flat2tuple(flat,tuple)
          implicit none
          integer(INTD), intent(in):: flat
          integer(INTD), intent(inout):: tuple(1:*)
          integer(INTD):: i,n

          n=flat
          do i=DIM_NUM_LEVELS-1,1,-1
           tuple(i)=mod(n,SPLIT_BASE)
           n=n/SPLIT_BASE
          enddo
          return
         end subroutine flat2tuple

         function tuple2flat(tuple) result(flat)
          implicit none
          integer(INTD):: flat
          integer(INTD), intent(in):: tuple(1:*)
          integer(INTD):: i

          flat=0
          do i=1,DIM_NUM_LEVELS-1
           flat=flat*SPLIT_BASE+tuple(i)
          enddo
          return
         end function tuple2flat

         function tuple2segment(tuple,dimsn) result(segment)
          implicit none
          integer(INTD):: segment
          integer(INTD), intent(in):: tuple(1:*)
          integer(INTD), intent(in):: dimsn
          integer(INTD):: i,j,m

          segment=0
          do i=1,DIM_NUM_LEVELS-1
           m=tuple(i)
           do j=1,dimsn-1; m=m/2; enddo
           segment=segment*2+mod(m,2)
          enddo
          return
         end function tuple2segment

         function flat2segment(flat,dimsn) result(segment)
          implicit none
          integer(INTD):: segment
          integer(INTD), intent(in):: flat
          integer(INTD), intent(in):: dimsn
          integer(INTD):: tuple(1:DIM_NUM_LEVELS-1)

          call flat2tuple(flat,tuple)
          segment=tuple2segment(tuple,dimsn)
          return
         end function flat2segment

         subroutine tuple2signa(tuple,signa)
          implicit none
          integer(INTD), intent(in):: tuple(1:*)
          integer(INTD), intent(inout):: signa(1:*)
          integer(INTD):: i

          do i=1,TENS_RANK
           signa(i)=tuple2segment(tuple,i)
          enddo
          return
         end subroutine tuple2signa

         subroutine signa2tuple(signa,tuple)
          implicit none
          integer(INTD), intent(in):: signa(1:*)
          integer(INTD), intent(inout):: tuple(1:*)
          integer(INTD):: i,j,m

          tuple(1:DIM_NUM_LEVELS-1)=0
          do i=TENS_RANK,1,-1
           m=signa(i)
           do j=DIM_NUM_LEVELS-1,1,-1
            tuple(j)=tuple(j)*2+mod(m,2)
            m=m/2
           enddo
          enddo
          return
         end subroutine signa2tuple

         subroutine flat2signa(flat,signa)
          implicit none
          integer(INTD), intent(in):: flat
          integer(INTD), intent(inout):: signa(1:*)
          integer(INTD):: tuple(1:DIM_NUM_LEVELS-1)

          call flat2tuple(flat,tuple)
          call tuple2signa(tuple,signa)
          return
         end subroutine flat2signa

         function signa2flat(signa) result(flat)
          implicit none
          integer(INTD):: flat
          integer(INTD), intent(in):: signa(1:*)
          integer(INTD):: tuple(1:DIM_NUM_LEVELS-1)

          call signa2tuple(signa,tuple)
          flat=tuple2flat(tuple)
          return
         end function signa2flat

        end subroutine test_carma
!---------------------------------------------
        subroutine tavp_wrk_reset_output(devo)
         implicit none
         integer(INTD), intent(in):: devo
         CONS_OUT=devo
         return
        end subroutine tavp_wrk_reset_output
!tens_resrc_t]==========================================
        function TensResrcIsEmpty(this,ierr) result(ans)
!Returns TRUE if the tensor resource is empty (unacquired).
         implicit none
         logical:: ans                               !out: answer
         class(tens_resrc_t), intent(in):: this      !in: tensor resource
         integer(INTD), intent(out), optional:: ierr !out: error code

         ans=(this%bytes.le.0_C_SIZE_T)
         if(present(ierr)) ierr=0
         return
        end function TensResrcIsEmpty
!---------------------------------------------------------------------------
        subroutine TensResrcAllocateBuffer(this,bytes,ierr,in_buffer,dev_id)
!Allocates local memory either from a system or from a custom buffer.
         implicit none
         class(tens_resrc_t), intent(inout):: this    !inout: tensor resource
         integer(INTL), intent(in):: bytes            !in: size in bytes
         integer(INTD), intent(out), optional:: ierr  !out: error code
         logical, intent(in), optional:: in_buffer    !in: if TRUE the memory will be allocated from a custom buffer, FALSE from the system
         integer(INTD), intent(in), optional:: dev_id !in: flat device id (defaults to Host)
         integer(INTD):: errc
         integer(C_INT):: in_buf,dev
         type(C_PTR):: addr

         if(this%is_empty(errc)) then
          if(bytes.gt.0_INTL) then
           in_buf=NOPE; if(present(in_buffer)) then; if(in_buffer) in_buf=YEP; endif
           dev=talsh_flat_dev_id(DEV_HOST,0); if(present(dev_id)) dev=dev_id
           errc=mem_allocate(dev,int(bytes,C_SIZE_T),in_buf,addr)
           if(errc.eq.0) then
            this%base_addr=addr
            this%bytes=bytes
            this%pinned=(in_buf.ne.NOPE)
            this%dev_id=dev
           else
            errc=-1 !`Here TRY_LATER should be distinguished from fatal errors
           endif
          else
           errc=-2
          endif
         else
          errc=-3
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensResrcAllocateBuffer
!------------------------------------------------
        subroutine TensResrcFreeBuffer(this,ierr)
!Frees the tensor resource buffer.
         implicit none
         class(tens_resrc_t), intent(inout):: this    !inout: tensor resource
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         if(.not.this%is_empty(errc)) then !free only allocated resources
          if(this%ref_count.le.1) then !at most one tensor operand can still be associated with this resource
           errc=mem_free(this%dev_id,this%base_addr)
           if(errc.eq.0) then
            this%base_addr=C_NULL_PTR
            this%bytes=0_C_SIZE_T
            this%pinned=.FALSE.
            this%dev_id=DEV_NULL
           else
            errc=-2
           endif
          else
           errc=-1
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensResrcFreeBuffer
!-----------------------------------------------------------------
        function TensResrcGetMemPtr(this,ierr,bytes) result(mem_p)
!Returns a C pointer to the local memory buffer used by the resource.
         implicit none
         type(C_PTR):: mem_p                          !out: C pointer
         class(tens_resrc_t), intent(in):: this       !in: tensor resource
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTL), intent(out), optional:: bytes !out: number of bytes
         integer(INTD):: errc

         mem_p=C_NULL_PTR
         if(.not.this%is_empty(errc)) then
          mem_p=this%base_addr
          if(present(bytes)) bytes=this%bytes
         else
          if(present(bytes)) bytes=0_INTL
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensResrcGetMemPtr
!------------------------------------------------------------
        function TensResrcGetMemSize(this,ierr) result(bytes)
!Returns the size of the resource memory buffer in bytes.
         implicit none
         integer(INTL):: bytes                       !out: size in bytes
         class(tens_resrc_t), intent(in):: this      !in: tensor resource
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(.not.this%is_empty(errc)) then
          bytes=this%bytes
         else
          bytes=0_INTL; errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensResrcGetMemSize
!---------------------------------------------
        subroutine TensResrcIncrRefCount(this)
!Increments the reference count.
         implicit none
         class(tens_resrc_t), intent(inout):: this !inout: tensor resource

!$OMP ATOMIC UPDATE
         this%ref_count=this%ref_count+1
         return
        end subroutine TensResrcIncrRefCount
!---------------------------------------------
        subroutine TensResrcDecrRefCount(this)
!Decrements the reference count.
         implicit none
         class(tens_resrc_t), intent(inout):: this !inout: tensor resource

!$OMP ATOMIC UPDATE
         this%ref_count=this%ref_count-1
         return
        end subroutine TensResrcDecrRefCount
!------------------------------------------------------
        function TensResrcGetRefCount(this) result(cnt)
!Returns the current reference count.
         implicit none
         integer(INTD):: cnt                    !out: reference count
         class(tens_resrc_t), intent(in):: this !in: tensor resource

!$OMP ATOMIC READ
         cnt=this%ref_count
         return
        end function TensResrcGetRefCount
!---------------------------------------
        subroutine tens_resrc_dtor(this)
         implicit none
         type(tens_resrc_t):: this
         integer(INTD):: errc

         call this%free_buffer(errc)
         if(errc.ne.0)&
         &call quit(errc,'#FATAL(TAVP-WRK:tens_resrc_dtor): Attempt to free a resource associated with multiple tensor operands!')
         return
        end subroutine tens_resrc_dtor
![tens_entry_wrk_t]==================================
        subroutine TensEntryWrkCtor(this,tensor,ierr)
!Constructs a <tens_entry_wrk_t>. Note move semantics for <tensor>
!due to the transfer of the tensor ownership to the tensor cache entry!
         implicit none
         class(tens_entry_wrk_t), intent(out):: this          !out: specialized tensor cache entry
         class(tens_rcrsv_t), pointer, intent(inout):: tensor !inout: pointer to an allocated tensor (ownership transfer will occur here!)
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTD):: errc

         errc=0
         if(associated(tensor)) then
!$OMP CRITICAL (TAVP_WRK_CACHE)
          call this%incr_use_count()
          call this%set_tensor(tensor,errc); if(errc.ne.0) errc=-2
          if(errc.eq.0) tensor=>NULL() !transfer the ownership
          call this%decr_use_count()
!$OMP END CRITICAL (TAVP_WRK_CACHE)
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensEntryWrkCtor
!---------------------------------------------------------------------
        function TensEntryWrkGetResource(this,ierr) result(resource_p)
!Returns a pointer to the tensor cache entry resource.
         implicit none
         class(tens_resrc_t), pointer:: resource_p          !out: pointer to the tensor resource
         class(tens_entry_wrk_t), intent(in), target:: this !in: specialized tensor cache entry
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc

         resource_p=>NULL()
         if(this%is_set(errc)) then
          if(errc.eq.0) then
           resource_p=>this%resource
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensEntryWrkGetResource
!---------------------------------------------------------------
        subroutine TensEntryWrkSetTensorLayout(this,ierr,tensor)
!Sets the tensor layout, either the default one or imported from another tensor.
!If the tensor stored in the tensor cache entry already has a layout, nothing will be done.
         implicit none
         class(tens_entry_wrk_t), intent(inout):: this      !inout: tensor cache entry
         integer(INTD), intent(out), optional:: ierr        !out: error code
         class(tens_rcrsv_t), intent(in), optional:: tensor !in: prototype tensor whose layout to be imported
         integer(INTD):: errc
         class(tens_rcrsv_t), pointer:: tens
         class(tens_header_t), pointer:: header
         logical:: laid

!$OMP CRITICAL (TAVP_WRK_CACHE)
         tens=>this%get_tensor(errc)
         if(errc.eq.0.and.associated(tens)) then
          if(tens%is_set(errc,layed=laid)) then
           if(errc.eq.TEREC_SUCCESS) then
            if(.not.laid) then
             if(present(tensor)) then !import layout from an existing tensor
              call quit(-1,'#FATAL(TAVP-WRK:tens_entry_wrk_t.set_tensor_layout): Tensor layout import is not implemented!') !`Implement layout import
             else !set the default layout
 !Set tensor composition, if not set:
              if(.not.tens%has_structure(errc)) then
               if(errc.eq.TEREC_SUCCESS) then
                header=>tens%get_header(errc)
                if(errc.eq.TEREC_SUCCESS) then
                 call tens%add_subtensor(header,errc); if(errc.ne.TEREC_SUCCESS) errc=-8
                endif
               else
                errc=-7
               endif
              else
               if(errc.ne.TEREC_SUCCESS) errc=-6
              endif
 !Resolve tensor dimensions, if not resolved:
              if(errc.eq.0) then
               errc=tens_dim_resolve(tens)
 !Set physical layout:
               if(errc.eq.TEREC_SUCCESS) then
                call tens%set_layout(TEREC_LAY_FDIMS,errc)
                if(errc.eq.TEREC_SUCCESS) then
                 if(DEBUG.gt.1) then
                  write(CONS_OUT,'("#DEBUG(TAVP-WRK:tens_entry_wrk_t.set_tensor_layout)[",i6,'//&
                  &'"]: Tensor storage layout successfully created for tensor:")') impir
                  call tens%print_it(dev_id=CONS_OUT)
                  flush(CONS_OUT)
                 endif
                else
                 errc=-5
                endif
               else
                errc=-4
               endif
              endif
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
!$OMP END CRITICAL (TAVP_WRK_CACHE)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensEntryWrkSetTensorLayout
!--------------------------------------------------------
        subroutine TensEntryWrkAcquireResource(this,ierr)
!Acquires resource for the tensor cache entry if it has not been acquired yet.
         implicit none
         class(tens_entry_wrk_t), intent(inout):: this !inout: tensor cache entry
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: errc
         integer(INTL):: buf_size
         class(tens_rcrsv_t), pointer:: tensor
         class(tens_body_t), pointer:: body
         class(tens_layout_t), pointer:: layout

!$OMP CRITICAL (TAVP_WRK_CACHE)
         if(this%resource%is_empty(errc)) then
          if(errc.eq.0) then
           tensor=>this%get_tensor(errc)
           if(errc.eq.0.and.associated(tensor)) then
            body=>tensor%get_body(errc)
            if(errc.eq.TEREC_SUCCESS.and.associated(body)) then
             layout=>body%get_layout(errc)
             if(errc.eq.TEREC_SUCCESS.and.associated(layout)) then
              buf_size=0; buf_size=layout%get_body_size(errc)
              if(errc.eq.TEREC_SUCCESS.and.buf_size.gt.0) then
               call this%resource%allocate_buffer(buf_size,errc); if(errc.ne.0) errc=-7
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
          if(errc.ne.0) errc=-1
         endif
!$OMP END CRITICAL (TAVP_WRK_CACHE)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensEntryWrkAcquireResource
!--------------------------------------------------------
        subroutine TensEntryWrkReleaseResource(this,ierr)
!Releases resource for the tensor cache entry if it has not been released yet.
         implicit none
         class(tens_entry_wrk_t), intent(inout):: this !inout: tensor cache entry
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: errc

!$OMP CRITICAL (TAVP_WRK_CACHE)
         if(.not.this%resource%is_empty(errc)) then
          if(errc.eq.0) then
           if((.not.this%is_persistent()).and.(this%get_ref_count().le.1)) then !at most one tensor operand can still be associated with this temporary cache entry
            call this%resource%free_buffer(errc); if(errc.ne.0) errc=-3
           endif
          else
           errc=-2
          endif
         else
          if(errc.ne.0) errc=-1
         endif
!$OMP END CRITICAL (TAVP_WRK_CACHE)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensEntryWrkReleaseResource
!-------------------------------------------
        subroutine tens_entry_wrk_dtor(this)
         implicit none
         type(tens_entry_wrk_t):: this
         integer(INTD):: errc

         call this%release_resource(errc) !release tensor cache entry resource
         if(errc.ne.0) then
          call quit(errc,'#FATAL(tens_entry_wrk_t.dtor): Unable to release tensor cache entry resource!')
         endif
         call this%destroy(.TRUE.) !deallocate the tensor component (if it is set)
         return
        end subroutine tens_entry_wrk_dtor
!-------------------------------------------------------------
        function tens_entry_wrk_alloc(tens_entry) result(ierr)
!Non-member allocator for tens_entry_wrk_t.
         implicit none
         integer(INTD):: ierr
         class(tens_cache_entry_t), allocatable, intent(out):: tens_entry

         allocate(tens_entry_wrk_t::tens_entry,STAT=ierr)
         return
        end function tens_entry_wrk_alloc
![tens_oprnd_t]=======================================================
        subroutine TensOprndCtorTensor(this,tensor,ierr,tens_resource)
!Constructs a tensor operand. The <tensor> must be set (defined).
!The associated tensor resource is optional and may still be empty.
         implicit none
         class(tens_oprnd_t), intent(inout):: this                            !inout: empty tensor operand (on entrance)
         class(tens_rcrsv_t), intent(in), target:: tensor                     !in: defined tensor
         integer(INTD), intent(out), optional:: ierr                          !out: error code
         class(tens_resrc_t), intent(inout), target, optional:: tens_resource !in: local tensor resource (may still be empty)
         integer(INTD):: errc

         if(.not.this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(tensor%is_set(errc)) then
            if(errc.eq.TEREC_SUCCESS) then
             this%cache_entry=>NULL()
             this%tensor=>tensor
             if(present(tens_resource)) then
              call tens_resource%incr_ref_count()
              this%resource=>tens_resource
             else
              this%resource=>NULL()
             endif
             call this%mark_active(errc); if(errc.ne.DSVP_SUCCESS) errc=-5
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
         class(tens_entry_wrk_t), intent(inout), target:: tens_cache_entry !in: tensor cache entry owning the tensor (and its resource)
         integer(INTD), intent(out), optional:: ierr                       !out: error code
         integer(INTD):: errc

         if(.not.this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(tens_cache_entry%is_set(errc)) then
            if(errc.eq.0) then
             call tens_cache_entry%incr_ref_count()
             this%cache_entry=>tens_cache_entry
             this%tensor=>this%cache_entry%get_tensor(errc)
             if(errc.eq.0) then
              this%resource=>this%cache_entry%get_resource() !may be empty resource
              if(associated(this%resource)) call this%resource%incr_ref_count()
              call this%mark_active(errc); if(errc.ne.DSVP_SUCCESS) errc=-6
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
        end subroutine TensOprndCtorCache
!------------------------------------------------------------
        function TensOprndGetTensor(this,ierr) result(tens_p)
!Returns a pointer to the tensor.
         implicit none
         class(tens_rcrsv_t), pointer:: tens_p       !out: pointer to the tensor
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; tens_p=>this%tensor
         if(.not.associated(this%tensor)) errc=-1
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetTensor
!---------------------------------------------------------------
        subroutine TensOprndSetCacheEntry(this,cache_entry,ierr)
!Sets the associated tensor cache entry (may be NULL). If the cache entry
!is not NULL, it must be the one containing the tensor set in the tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this                  !inout: active tensor operand
         class(tens_entry_wrk_t), intent(in), pointer:: cache_entry !in: pointer to a tensor cache entry containing the same tensor (may be NULL)
         integer(INTD), intent(out), optional:: ierr                !out: error code
         integer(INTD):: errc
         class(tens_rcrsv_t), pointer:: tensor

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(associated(this%cache_entry)) then
            call this%cache_entry%decr_ref_count(); this%cache_entry=>NULL()
           endif
           if(associated(cache_entry)) then
            tensor=>cache_entry%get_tensor(errc)
            if(errc.eq.0.and.associated(this%tensor,tensor)) then !cache entry must correspond to the same tensor
             call cache_entry%incr_ref_count()
            else
             errc=-3
            endif
           endif
           if(errc.eq.0) this%cache_entry=>cache_entry
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndSetCacheEntry
!-----------------------------------------------------------------------
        function TensOprndGetCacheEntry(this,ierr) result(cache_entry_p)
!Returns a pointer to the tensor cache entry (may be NULL).
         implicit none
         class(tens_entry_wrk_t), pointer:: cache_entry_p !out: pointer to the tensor cache entry
         class(tens_oprnd_t), intent(in):: this           !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD):: errc

         errc=0; cache_entry_p=>this%cache_entry
         if(.not.associated(this%tensor)) errc=-1 !trap
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetCacheEntry
!---------------------------------------------------------------
        subroutine TensOprndSetResource(this,tens_resource,ierr)
!Sets the resource component if it has not been set via the constructor.
         implicit none
         class(tens_oprnd_t), intent(inout):: this                  !inout: active tensor operand
         class(tens_resrc_t), intent(inout), target:: tens_resource !inout: local tensor resource (may still be empty)
         integer(INTD), intent(out), optional:: ierr                !out: error code
         integer(INTD):: errc

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(.not.associated(this%resource)) then
            call tens_resource%incr_ref_count()
            this%resource=>tens_resource
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
        end subroutine TensOprndSetResource
!------------------------------------------------------------------
        function TensOprndGetResource(this,ierr) result(resource_p)
!Returns a pointer to the tensor resource (may be NULL).
         implicit none
         class(tens_resrc_t), pointer:: resource_p   !out: pointer to the tensor resource
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; resource_p=>this%resource
         if(.not.associated(this%tensor)) errc=-1 !trap
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetResource
!--------------------------------------------------
        subroutine TensOprndSetTalshTens(this,ierr)
!Sets up the TAL-SH tensor object for further processing with TAL-SH.
!The local tensor body location is imported from the tensor resource.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,nd,unres,lay,data_kind,dims(1:MAX_TENSOR_RANK)
         integer(INTL):: dim_ext(1:MAX_TENSOR_RANK)
         logical:: shpd,layd,locd
         type(C_PTR):: mem_p
         class(tens_rcrsv_t), pointer:: tens_p
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p
         class(tens_header_t), pointer:: header_p

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(talsh_tensor_is_empty(this%talsh_tens)) then
            tens_p=>this%tensor
            if(tens_p%is_set(errc,num_dims=nd,shaped=shpd,unresolved=unres,layed=layd,located=locd)) then
             if((errc.eq.0).and.(unres.eq.0).and.shpd.and.layd) then
              body_p=>tens_p%get_body(errc)
              if(errc.eq.TEREC_SUCCESS) then
               layout_p=>body_p%get_layout(errc)
               if(errc.eq.TEREC_SUCCESS) then
                lay=layout_p%get_layout_kind(errc)
                if(errc.eq.TEREC_SUCCESS.and.lay.eq.TEREC_LAY_FDIMS) then
                 data_kind=layout_p%get_data_type(errc)
                 if(errc.eq.TEREC_SUCCESS) then
                  if(.not.this%resource%is_empty()) then
                   mem_p=this%resource%get_mem_ptr(errc)
                   if(errc.eq.0) then
                    header_p=>tens_p%get_header(errc)
                    if(errc.eq.TEREC_SUCCESS) then
                     call header_p%get_dims(dim_ext,nd,errc)
                     if(errc.eq.TEREC_SUCCESS) then
                      dims(1:nd)=dim_ext(1:nd)
                      errc=talsh_tensor_construct(this%talsh_tens,data_kind,dims(1:nd),ext_mem=mem_p)
                      if(errc.ne.TALSH_SUCCESS) errc=-14
                     else
                      errc=-13
                     endif
                    else
                     errc=-12
                    endif
                    header_p=>NULL()
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
               layout_p=>NULL()
              else
               errc=-6
              endif
              body_p=>NULL()
             else
              errc=-5
             endif
            else
             errc=-4
            endif
            tens_p=>NULL()
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
        end subroutine TensOprndSetTalshTens
!------------------------------------------------------------
        subroutine TensOprndSetTensorLayout(this,ierr,tensor)
!Sets the tensor layout, if not already set, either the default one
!or by importing it from an existing tensor.
         implicit none
         class(tens_oprnd_t), intent(inout):: this          !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr        !out: error code
         class(tens_rcrsv_t), intent(in), optional:: tensor !in: tensor whose layout to be imported
         integer(INTD):: errc

         errc=0
         if(associated(this%cache_entry)) then
          if(present(tensor)) then
           call this%cache_entry%set_tensor_layout(errc,tensor)
          else
           call this%cache_entry%set_tensor_layout(errc)
          endif
         else
          call quit(-1,'#FATAL(TAVP-WRK:tens_oprnd_t.set_tensor_layout): Direct layout setup is not implemented!') !`Implement direct tensor layout setup
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndSetTensorLayout
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
!---------------------------------------------------------
        function TensOprndIsLocated(this,ierr) result(res)
!Returns TRUE if the tensor operand has been located, FALSE otherwise.
!By being located, it means its physical location is known.
         implicit none
         logical:: res                               !out: result
         class(tens_oprnd_t), intent(inout):: this   !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         logical:: laid,locd

         res=.FALSE.
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
!$OMP CRITICAL (TAVP_WRK_CACHE)
           if(associated(this%cache_entry)) call this%cache_entry%incr_use_count()
           if(this%tensor%is_set(errc,layed=laid,located=locd)) then
            if(errc.eq.TEREC_SUCCESS) then
             res=(laid.and.locd)
            else
             errc=-4
            endif
           else
            errc=-3
           endif
           if(associated(this%cache_entry)) call this%cache_entry%decr_use_count()
!$OMP END CRITICAL (TAVP_WRK_CACHE)
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
         class(tens_oprnd_t), intent(inout):: this   !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer(INT_MPI):: host_proc_rank,mpi_comm,my_rank
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p
         class(DataDescr_t), pointer:: descr_p

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
!$OMP CRITICAL (TAVP_WRK_CACHE)
           if(associated(this%cache_entry)) call this%cache_entry%incr_use_count()
           body_p=>this%tensor%get_body(errc)
           if(errc.eq.TEREC_SUCCESS.and.associated(body_p)) then
            layout_p=>body_p%get_layout(errc)
            if(errc.eq.TEREC_SUCCESS.and.associated(layout_p)) then
             descr_p=>layout_p%get_data_descr(errc)
             if(errc.eq.TEREC_SUCCESS.and.associated(descr_p)) then
              if(descr_p%is_set(errc,host_proc_rank,mpi_comm)) then
               if(errc.eq.0) then
                call MPI_Comm_Rank(mpi_comm,my_rank,errc)
                if(errc.eq.0) then
                 res=.not.(host_proc_rank.eq.my_rank)
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
           if(associated(this%cache_entry)) call this%cache_entry%decr_use_count()
!$OMP END CRITICAL (TAVP_WRK_CACHE)
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
!In order to valued for TAVP-WRK, the tensor has to be laid out and located.
         implicit none
         logical:: res                               !out: result
         class(tens_oprnd_t), intent(inout):: this   !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         logical:: laid,locd

         res=.FALSE.
         if(this%is_active(errc)) then
!$OMP CRITICAL (TAVP_WRK_CACHE)
          if(associated(this%cache_entry)) call this%cache_entry%incr_use_count()
          if(this%tensor%is_set(errc,layed=laid,located=locd)) then
           if(errc.eq.TEREC_SUCCESS) then
            res=(laid.and.locd.and.(this%get_write_count().eq.0)) !tensor is laid out, located, and no one currently writes into it => DEFINED
           else
            errc=-3
           endif
          else
           errc=-2
          endif
          if(associated(this%cache_entry)) call this%cache_entry%decr_use_count()
!$OMP END CRITICAL (TAVP_WRK_CACHE)
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensOprndIsValued
!-----------------------------------------------------------
        function TensOprndHasResource(this,ierr) result(res)
!Returns TRUE if the tensor operand has been allocated an actual local resource.
         implicit none
         logical:: res                               !out: result
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         res=.FALSE.
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(associated(this%resource)) then
            res=(.not.this%resource%is_empty(errc)); if(errc.ne.0) errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensOprndHasResource
!------------------------------------------------
        subroutine TensOprndAcquireRsc(this,ierr)
!Acquires local resources for a tensor operand.
!If the resources have already been allocated, does nothing.
!If the resource component is not set, an error will be returned.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand (with an associated resource component)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer(INTL):: buf_size
         integer(INT_MPI):: host_proc_rank
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p

         if(this%is_active(errc)) then
          if(associated(this%resource)) then
           if(this%resource%is_empty()) then
            body_p=>this%tensor%get_body(errc)
            if(errc.eq.TEREC_SUCCESS.and.associated(body_p)) then
             layout_p=>body_p%get_layout(errc)
             if(errc.eq.TEREC_SUCCESS.and.associated(layout_p)) then
              if(layout_p%is_set(errc)) then
               if(errc.eq.TEREC_SUCCESS) then
                buf_size=layout_p%get_body_size(errc)
                if(errc.eq.TEREC_SUCCESS.and.buf_size.gt.0_INTL) then
                 call this%resource%allocate_buffer(buf_size,errc); if(errc.ne.0) errc=-1
                 if(errc.eq.0.and.DEBUG.gt.0) then
                  write(CONS_OUT,'("#DEBUG(TAVP-WRK:tens_oprnd_t:acquire_rsc)[",i6,"]: Memory acquired: Size (B) = ",i13)')&
                  &impir,buf_size
                  flush(CONS_OUT)
                 endif
                else
                 errc=-2
                endif
               else
                errc=-3
               endif
              else
               errc=-4
              endif
             else
              errc=-5
             endif
            else
             errc=-6
            endif
           endif
          else
           errc=-7
          endif
         else
          errc=-8
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndAcquireRsc
!----------------------------------------------
        subroutine TensOprndPrefetch(this,ierr)
!Starts prefetching the (remote) tensor operand using the local tensor resource.
!If the resource component has not been set, an error will be returned.
!If the local resource has not been allocated, it will be allocated here.
!If the tensor operand has been delivered before, does nothing.
!If there is a pending communication on the tensor operand, returns an error.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand (with an associated resource component)
         integer(INTD), intent(out), optional:: ierr !out: error code, may return TRY_LATER
         integer(INTD):: errc
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p
         class(DataDescr_t), pointer:: descr_p
         type(C_PTR):: cptr

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(.not.this%is_present(errc)) then
            if(errc.eq.DSVP_SUCCESS) then
             if(associated(this%resource)) then
              body_p=>this%tensor%get_body(errc)
              if(errc.eq.TEREC_SUCCESS.and.associated(body_p)) then
               layout_p=>body_p%get_layout(errc)
               if(errc.eq.TEREC_SUCCESS.and.associated(layout_p)) then
                descr_p=>layout_p%get_data_descr(errc)
                if(errc.eq.TEREC_SUCCESS.and.associated(descr_p)) then
                 if(descr_p%is_set(errc)) then
                  if(errc.eq.0) then
                   if(this%resource%is_empty()) call this%acquire_rsc(errc)
                   if(errc.eq.0) then
                    if(this%get_comm_stat().eq.DS_OPRND_NO_COMM) then
                     cptr=this%resource%get_mem_ptr(errc)
                     if(errc.eq.0) then
                      call descr_p%get_data(cptr,errc,MPI_ASYNC_REQ)
                      if(errc.ne.0.and.errc.ne.TRY_LATER) errc=-1
                      if(errc.eq.0) then
                       call this%set_comm_stat(DS_OPRND_FETCHING,errc); if(errc.ne.DSVP_SUCCESS) errc=-2
                      endif
                     else
                      errc=-3
                     endif
                    else
                     errc=-4
                    endif
                   else
                    errc=-5
                   endif
                  else
                   errc=-6
                  endif
                 else
                  errc=-7
                 endif
                else
                 errc=-8
                endif
               else
                errc=-9
               endif
              else
               errc=-10
              endif
             else
              errc=-11
             endif
            else
             errc=-12
            endif
           else
            if(errc.ne.DSVP_SUCCESS) errc=-13
           endif
          else
           errc=-14
          endif
         else
          errc=-15
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndPrefetch
!--------------------------------------------
        subroutine TensOprndUpload(this,ierr)
!Starts uploading the (remote) tensor operand from the local tensor resource.
!The tensor operand must be marked as delivered (present), even if it is local.
!If there is a pending communication on the tensor operand, returns an error.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p
         class(DataDescr_t), pointer:: descr_p
         type(C_PTR):: cptr

         if(this%is_present(errc)) then !assumes that the local resource is allocated
          if(errc.eq.DSVP_SUCCESS) then
           body_p=>this%tensor%get_body(errc)
           if(errc.eq.TEREC_SUCCESS.and.associated(body_p)) then
            layout_p=>body_p%get_layout(errc)
            if(errc.eq.TEREC_SUCCESS.and.associated(layout_p)) then
             descr_p=>layout_p%get_data_descr(errc)
             if(errc.eq.TEREC_SUCCESS.and.associated(descr_p)) then
              if(descr_p%is_set(errc)) then
               if(errc.eq.0) then
                if(.not.this%resource%is_empty()) then !trap
                 if(this%get_comm_stat().eq.DS_OPRND_NO_COMM) then
                  cptr=this%resource%get_mem_ptr(errc)
                  if(errc.eq.0) then
                   call descr_p%acc_data(cptr,errc,MPI_ASYNC_REQ)
                   if(errc.ne.0.and.errc.ne.TRY_LATER) errc=-1
                   if(errc.eq.0) then
                    call this%set_comm_stat(DS_OPRND_UPLOADING,errc); if(errc.ne.DSVP_SUCCESS) errc=-2
                   endif
                  else
                   errc=-3
                  endif
                 else
                  errc=-4
                 endif
                else
                 errc=-5
                endif
               else
                errc=-6
               endif
              else
               errc=-7
              endif
             else
              errc=-8
             endif
            else
             errc=-9
            endif
           else
            errc=-10
           endif
          else
           errc=-11
          endif
         else
          errc=-12
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndUpload
!---------------------------------------------------------
        function TensOprndSync(this,ierr,wait) result(res)
!Synchronizes a pending prefetch/upload, either TEST or WAIT.
!A successful synchronization on prefetch will mark the tensor operand
!as delivered (present). A successful synchronization on upload will
!not change the status of the tensor operand (which is present).
!`An attempt to synchronize a non-existing communication will cause an error.
         implicit none
         logical:: res                               !out: TRUE on communication completion, FALSE otherwise
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: TRUE activates WAIT instead of TEST synchronization (default)
         integer(INTD):: errc,sts
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p
         class(DataDescr_t), pointer:: descr_p
         logical:: tw

         res=.FALSE.; tw=.FALSE.; if(present(wait)) tw=wait
         if(this%is_active(errc)) then
          body_p=>this%tensor%get_body(errc)
          if(errc.eq.TEREC_SUCCESS.and.associated(body_p)) then
           layout_p=>body_p%get_layout(errc)
           if(errc.eq.TEREC_SUCCESS.and.associated(layout_p)) then
            descr_p=>layout_p%get_data_descr(errc)
            if(errc.eq.TEREC_SUCCESS.and.associated(descr_p)) then
             if(descr_p%is_set(errc)) then
              if(errc.eq.0) then
               sts=this%get_comm_stat()
               if(sts.ne.DS_OPRND_NO_COMM) then
                if(tw) then
                 call descr_p%wait_data(errc); if(errc.eq.0) then; res=.TRUE.; else; errc=-1; endif
                else
                 res=descr_p%test_data(errc); if(errc.ne.0) errc=-2
                endif
                if(sts.eq.DS_OPRND_FETCHING.and.res) then
                 call this%mark_delivered(errc); if(errc.ne.DSVP_SUCCESS) errc=-3
                endif
               else
                errc=-4 !`if there is no pending communication, it is an error?
               endif
              else
               errc=-5
              endif
             else
              errc=-6
             endif
            else
             errc=-7
            endif
           else
            errc=-8
           endif
          else
           errc=-9
          endif
         else
          errc=-10
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensOprndSync
!---------------------------------------------
        subroutine TensOprndRelease(this,ierr)
!Releases local tensor resources occupied by the tensor operand,
!unless there are other active tensor operands sharing the same resource.
!In the latter case, nothing will be done and no error raised.
!Also, if the resource component is not set, nothing will be done either.
!Note that the resource object itself (either allocated or unallocated) is
!still associated with this tensor operand, until it is destructed.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         logical:: delivered,persistent

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(associated(this%resource)) then
            if(this%get_comm_stat().ne.DS_OPRND_NO_COMM) then
             delivered=this%sync(errc,wait=.TRUE.)
             if((.not.delivered).or.(errc.ne.0)) errc=-1
            endif
            persistent=.FALSE.; if(associated(this%cache_entry)) persistent=this%cache_entry%is_persistent()
            if((.not.persistent).and.(this%resource%get_ref_count().eq.1)) then !only one (last) tensor operand is associated with this resource
             call this%resource%free_buffer(errc); if(errc.ne.0) errc=-2 !free the resource memory buffer
            endif
           endif
          else
           errc=-3
          endif
         else
          if(errc.ne.DSVP_SUCCESS) errc=-4
         endif
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
            if(.not.delivered.or.errc.ne.DSVP_SUCCESS) errc=-1
           endif
           call this%mark_empty(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-2
           if(.not.talsh_tensor_is_empty(this%talsh_tens)) then
            ier=talsh_tensor_destruct(this%talsh_tens)
            if(ier.ne.TALSH_SUCCESS.and.errc.eq.0) errc=-3
           endif
           if(associated(this%resource)) then
            call this%resource%decr_ref_count()
            this%resource=>NULL()
           endif
           if(associated(this%cache_entry)) then
            call this%cache_entry%decr_ref_count()
            this%cache_entry=>NULL()
           endif
           this%tensor=>NULL()
          else
           errc=-4
          endif
         else
          if(errc.ne.DSVP_SUCCESS) errc=-5
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndDestruct
!-------------------------------------------------------------
        subroutine TensOprndUploadLocal(this,cache_entry,ierr)
!Performs a local accumulation of the tensor into another cached tensor.
         implicit none
         class(tens_oprnd_t), intent(in):: this               !in: tensor operand being accumulated
         class(tens_entry_wrk_t), intent(inout):: cache_entry !inout: tensor cache entry containing the tensor being accumulated into
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTD):: errc,dtk,dts
         integer(INTL):: src_size,dst_size,vol,i
         class(tens_resrc_t), pointer:: dst_resource
         type(C_PTR):: src,dst
         real(4), pointer:: src_r4(:),dst_r4(:)
         real(8), pointer:: src_r8(:),dst_r8(:)
         complex(4), pointer:: src_c4(:),dst_c4(:)
         complex(8), pointer:: src_c8(:),dst_c8(:)

         if(cache_entry%is_set(errc)) then
          if(errc.eq.0) then
           call cache_entry%incr_use_count()
           src=this%resource%get_mem_ptr(errc,src_size)
           if(errc.eq.0) then
            dst_resource=>cache_entry%get_resource(errc)
            if(errc.eq.0.and.associated(dst_resource)) then
             dst=dst_resource%get_mem_ptr(errc,dst_size)
             if(errc.eq.0) then
              if(dst_size.eq.src_size) then
               dtk=this%tensor%get_data_type(errc)
               if(errc.eq.TEREC_SUCCESS) then
                if(talsh_valid_data_kind(dtk,dts).eq.YEP) then
                 if(mod(src_size,int(dts,INTL)).eq.0) then
                  vol=src_size/int(dts,INTL)
                  select case(dtk)
                  case(R4)
                   call c_f_pointer(src,src_r4,(/vol/))
                   call c_f_pointer(dst,dst_r4,(/vol/))
                   do i=1,vol
                    dst_r4(i)=dst_r4(i)+src_r4(i)
                   enddo
                  case(R8)
                   call c_f_pointer(src,src_r8,(/vol/))
                   call c_f_pointer(dst,dst_r8,(/vol/))
                   do i=1,vol
                    dst_r8(i)=dst_r8(i)+src_r8(i)
                   enddo
                  case(C4)
                   call c_f_pointer(src,src_c4,(/vol/))
                   call c_f_pointer(dst,dst_c4,(/vol/))
                   do i=1,vol
                    dst_c4(i)=dst_c4(i)+src_c4(i)
                   enddo
                  case(C8)
                   call c_f_pointer(src,src_c8,(/vol/))
                   call c_f_pointer(dst,dst_c8,(/vol/))
                   do i=1,vol
                    dst_c8(i)=dst_c8(i)+src_c8(i)
                   enddo
                  case default
                   errc=-10
                  end select
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
           else
            errc=-3
           endif
           call cache_entry%decr_use_count()
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndUploadLocal
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
         class(tens_oprnd_t), intent(in):: this
         integer(INTD), intent(out), optional:: ierr
         integer(INTD), intent(in), optional:: dev_id
         integer(INTD), intent(in), optional:: nspaces
         integer(INTD):: errc,devo,nsp,j

         errc=0
         devo=6; if(present(dev_id)) devo=dev_id
         nsp=0; if(present(nspaces)) nsp=nspaces
         if(associated(this%tensor)) then
          call this%tensor%print_it(errc,devo,nsp); if(errc.ne.TEREC_SUCCESS) errc=-1
         else
          do j=1,nsp; write(devo,'(" ")',ADVANCE='NO'); enddo
          write(devo,'("No Tensor!")')
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndPrintIt
!------------------------------------------------------------------------
        subroutine TensOprndResetTmpTensor(this,cache_entry,forward,ierr)
!Resets the tensor in a tensor operand by providing another tensor cache entry
!with a temporary tensor (<forward> = TRUE) or vice versa (<forward> = FALSE).
!The reference count of the persistent tensor cache entry is kept unchanged.
         implicit none
         class(tens_oprnd_t), intent(inout):: this                  !inout: active tensor operand
         class(tens_entry_wrk_t), pointer, intent(in):: cache_entry !in: tensor cache entry containing the new tensor
         logical, intent(in):: forward                              !in: TRUE: persistent --> temporary; FALSE: temporary --> persistent
         integer(INTD), intent(out), optional:: ierr                !out: error code
         integer(INTD):: errc
         class(tens_entry_wrk_t), pointer:: tens_entry
         class(tens_rcrsv_t), pointer:: tensor_new

         errc=0
         if(associated(cache_entry)) then
          tens_entry=>this%get_cache_entry(errc)
          if(errc.eq.0.and.associated(tens_entry)) then
           if(forward) then
            call tens_entry%incr_ref_count() !set_cache_entry() will then decrement this reference count
           else
            call cache_entry%decr_ref_count() !set_cache_entry() will then increment this reference count
           endif
           tensor_new=>cache_entry%get_tensor(errc)
           if(errc.eq.0.and.associated(tensor_new)) then
            this%tensor=>tensor_new
            call this%set_cache_entry(cache_entry,errc); if(errc.ne.0) errc=-4
           else
            if(forward) then
             call tens_entry%decr_ref_count()
            else
             call cache_entry%incr_ref_count()
            endif
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
        end subroutine TensOprndResetTmpTensor
!---------------------------------------
        subroutine tens_oprnd_dtor(this)
         implicit none
         type(tens_oprnd_t):: this
         integer(INTD):: errc

         call this%destruct(errc)
         if(errc.ne.0) call quit(errc,'#FATAL(TAVP-WRK:tens_oprnd_dtor): Destructor failed!')
         return
        end subroutine tens_oprnd_dtor
![tens_instr_t]==============================================================
        subroutine TensInstrCtor(this,op_code,ierr,op_spec,iid,stat,err_code)
!Constructs a tensor instruction from a given tensor operation.
!The tensor instruction is a realization of a given tensor operation
!for a specific TAVP kind.
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
           case(TAVP_INSTR_CTRL_RESUME,TAVP_INSTR_CTRL_STOP)
            call construct_instr_ctrl(errc); if(errc.ne.0) errc=-10
           case(TAVP_INSTR_TENS_CREATE,TAVP_INSTR_TENS_DESTROY)
            call construct_instr_tens_create_destroy(errc); if(errc.ne.0) errc=-9
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
              call tens_oprnd%tens_oprnd_ctor(tensor,jerr)
              if(jerr.eq.0) then
               oprnd=>tens_oprnd
               call this%set_operand(0,oprnd,jerr)
               if(jerr.eq.DSVP_SUCCESS) then
                if(op_code.eq.TAVP_INSTR_TENS_CREATE) then
                 this%num_out_oprnds=1; this%out_oprnds(0:this%num_out_oprnds-1)=(/0/)
                elseif(op_code.eq.TAVP_INSTR_TENS_DESTROY) then
                 this%num_out_oprnds=0
                endif
               else
                jerr=-6
               endif
               oprnd=>NULL() !<oprnd> pointer was saved in the tensor instruction and will later be deallocated
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
! 5. Instruction operands (optional): {Owner_id,Read_count,Write_count,Tensor} for each tensor operand.
!NOTE: Owner_ID for each tensor operand originates from the TAVP-MNG tensor instruction format
!      since they need to comply. It is irrelevant for TAVP-WRK, resulting in a default value (-1).
         implicit none
         class(tens_instr_t), intent(in):: this          !in: defined tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet !out: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,op_code,stat,err_code
         integer(INTL):: iid

!Pack the instruction attributes (id,opcode,status,error):
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
                 case(TAVP_INSTR_CTRL_RESUME,TAVP_INSTR_CTRL_STOP)
                  call encode_instr_ctrl(errc); if(errc.ne.0) errc=-12
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

         subroutine encode_instr_ctrl(jerr)
          !TAVP CTRL instruction:
          !Packet format: {id|op_code|status|error}
          integer(INTD), intent(out):: jerr

          jerr=0
          return
         end subroutine encode_instr_ctrl

         subroutine encode_instr_tens_create_destroy(jerr)
          !CREATE/DESTROY a tensor:
          !Packet format: {id|opcode|status|error|tensor_operand0}
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
              call pack_builtin(instr_packet,-1,jerr) !owner id
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
          !Packed format: {id|opcode|status|error|ctrl_tens_contr_t|tensor_operand0,tensor_operand1,tensor_operand2}
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
              call pack_builtin(instr_packet,-1,jerr) !owner id
              if(jerr.eq.0) call pack_builtin(instr_packet,oprnd%get_read_count(),jerr)
              if(jerr.eq.0) call pack_builtin(instr_packet,oprnd%get_write_count(),jerr)
              if(jerr.eq.0) call tensor%pack(instr_packet,jerr)
              if(jerr.ne.0) then; jerr=-3; exit; endif
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
!------------------------------------------------------------------------
        function TensInstrGetNumOutOperands(this,ierr) result(num_oprnds)
!Returns the total number of output tensor operands in the tensor instruction.
         implicit none
         integer(INTD):: num_oprnds                  !out: total number of output tensor operands
         class(tens_instr_t), intent(in):: this      !in: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code

         num_oprnds=this%num_out_oprnds
         if(present(ierr)) ierr=0
         return
        end function TensInstrGetNumOutOperands
!----------------------------------------------------------------------
        function TensInstrOperandIsOutput(this,op_num,ierr) result(ans)
!Returns TRUE if the specific tensor instruction operand is output, FALSE otherwise.
         implicit none
         logical:: ans                               !out: answer
         class(tens_instr_t), intent(in):: this      !in: active tensor instruction
         integer(INTD), intent(in):: op_num          !in: operand number: [0..max]
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,n,i

         ans=.FALSE.
         n=this%get_num_operands(errc)
         if(errc.eq.DSVP_SUCCESS) then
          if(op_num.ge.0.and.op_num.lt.n) then
           do i=0,this%num_out_oprnds-1
            if(this%out_oprnds(i).eq.op_num) then; ans=.TRUE.; exit; endif
           enddo
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensInstrOperandIsOutput
!-------------------------------------------------------------------------------
        function TensInstrGetOutputOperands(this,ierr,num_oprs) result(out_oprs)
!Returns the list of the output tensor operands by their positions.
         implicit none
         integer(INTD), pointer:: out_oprs(:)            !out: positions of the output tensor instruction operands
         class(tens_instr_t), intent(in), target:: this  !in: active tensor instruction
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD), intent(out), optional:: num_oprs !out: number of the output tensor instruction operands
         integer(INTD):: errc,n

         out_oprs=>NULL(); n=0
         if(.not.this%is_empty(errc)) then
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
!-------------------------------------------------------
        subroutine TensInstrLayOutputOperands(this,ierr)
!Sets up the storage layout for non-existing output tensor operands.
         implicit none
         class(tens_instr_t), intent(inout):: this   !inout: active tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,i,j,opcode
         class(ds_oprnd_t), pointer:: oprnd

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           opcode=this%get_code(errc)
           if(errc.eq.DSVP_SUCCESS) then
            if(opcode.ge.TAVP_ISA_TENS_FIRST.and.opcode.le.TAVP_ISA_TENS_LAST) then
             select case(opcode)
             case(TAVP_INSTR_TENS_CREATE)
              do i=0,this%num_out_oprnds-1
               j=this%out_oprnds(i) !position of the output operand #i
               oprnd=>this%get_operand(j,errc); if(errc.ne.DSVP_SUCCESS) then; errc=-7; exit; endif
               select type(oprnd)
               class is(tens_oprnd_t)
                call oprnd%set_tensor_layout(errc); if(errc.ne.0) then; errc=-6; exit; endif
               class default
                errc=-5; exit
               end select
              enddo
             case default
              call quit(-1,'#FATAL(TAVP-WRK:tens_instr_t.lay_output_operands): Layout inferrence not implemented yet!') !`Implement output tensor layout inferrence for other tensor instructions
             end select
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
        end subroutine TensInstrLayOutputOperands
!-------------------------------------------------------------------------------
        subroutine TensInstrGetCacheEntries(this,cache_entries,num_entries,ierr)
!Returns an array of references to tensor cache entries associated with the tensor operands.
         implicit none
         class(tens_instr_t), intent(in):: this                        !in: tensor instruction
         type(tens_entry_wrk_ref_t), intent(inout):: cache_entries(1:) !out: references to tensor cache entries
         integer(INTD), intent(out):: num_entries                      !out: number of tensor cache entries returned
         integer(INTD), intent(out), optional:: ierr                   !out: error code
         integer(INTD):: errc,i
         class(ds_oprnd_t), pointer:: oprnd

         num_entries=0
         if(.not.this%is_empty(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           num_entries=this%get_num_operands(errc)
           if(errc.eq.DSVP_SUCCESS) then
            do i=0,num_entries-1
             oprnd=>this%get_operand(i,errc); if(errc.ne.DSVP_SUCCESS) then; errc=-4; exit; endif
             select type(oprnd)
             class is(tens_oprnd_t)
              cache_entries(1+i)%cache_entry=>oprnd%get_cache_entry(errc)
              if(errc.ne.0) then; errc=-3; exit; endif
             class default
              errc=-2; exit
             end select
            enddo
           endif
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensInstrGetCacheEntries
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
               tensor=>tens_oprnd%get_tensor(errc); if(errc.ne.0) then; errc=-6; exit oloop; endif
               call tensor%get_dims(dims,n,errc); if(errc.ne.TEREC_SUCCESS) then; errc=-5; exit oloop; endif
               vol=1d0; do j=1,n; vol=vol*real(dims(j),8); enddo
               tvol=tvol*vol; words=words+vol
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
         integer(INTD):: errc,opcode,numo,i
         complex(8):: alpha
         class(tens_rcrsv_t), pointer:: tensor
         class(ds_oprnd_t), pointer:: oprnd
         class(ds_instr_ctrl_t), pointer:: instr_ctrl
         type(contr_ptrn_ext_t), pointer:: contr_ptrn

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           numo=this%get_num_operands(errc)
           if(errc.eq.DSVP_SUCCESS) then
            opcode=this%get_code(errc)
            if(errc.eq.DSVP_SUCCESS) then
             select case(opcode)
             case(TAVP_INSTR_TENS_CREATE) !no associated tensor operation
             case(TAVP_INSTR_TENS_DESTROY) !no associated tensor operation
             case(TAVP_INSTR_TENS_CONTRACT)
              allocate(tens_contraction_t::tens_operation)
              select type(tens_operation)
              class is(tens_contraction_t)
               do i=0,numo-1
                oprnd=>this%get_operand(i,errc); if(errc.ne.DSVP_SUCCESS) exit
                tensor=>NULL(); select type(oprnd); class is(tens_oprnd_t); tensor=>oprnd%get_tensor(errc); end select
                if(.not.associated(tensor)) errc=-10; if(errc.ne.0) exit
                call tens_operation%set_argument(tensor,errc); if(errc.ne.TEREC_SUCCESS) exit
               enddo
               if(errc.eq.0) then
                instr_ctrl=>this%get_control(errc)
                if(errc.eq.DSVP_SUCCESS) then
                 select type(instr_ctrl)
                 class is(ctrl_tens_contr_t)
                  contr_ptrn=>instr_ctrl%get_contr_ptrn(errc,alpha)
                  if(errc.eq.0) call tens_operation%set_contr_ptrn(contr_ptrn,errc,alpha)
                  nullify(contr_ptrn)
                 class default
                  errc=-9
                 end select
                else
                 errc=-8
                endif
                nullify(instr_ctrl)
               else
                errc=-7
               endif
              class default
               errc=-6
              end select
             case default
              !`Implement for other relevant tensor instructions
              errc=-5
             end select
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
        end subroutine TensInstrGetOperation
!------------------------------------------------------------
        function TensInstrMarkIssue(this,ierr) result(passed)
!Updates the tensor access counters for each tensor operand when the tensor instruction is issued.
         implicit none
         logical:: passed                            !out: TRUE if there were no data dependency problem, FALSE otherwise
         class(tens_instr_t), intent(inout):: this   !inout: active tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,i,n,rd,wr
         class(ds_oprnd_t), pointer:: oprnd

         passed=.TRUE.
         if(.not.this%is_empty(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           n=this%get_num_operands(errc)
           if(errc.eq.DSVP_SUCCESS) then
            if(n.gt.0) then
!$OMP CRITICAL (TAVP_WRK_CACHE)
             do i=0,n-1
              oprnd=>this%get_operand(i,errc)
              if(errc.eq.DSVP_SUCCESS.and.associated(oprnd)) then
               select type(oprnd)
               class is(tens_oprnd_t)
                rd=oprnd%get_read_count(); wr=oprnd%get_write_count()
                if(this%operand_is_output(i,errc)) then !output tensor operand
                 passed=(passed.and.(rd.eq.0))
                 call oprnd%register_write()
                else !input tensor operand
                 passed=(passed.and.(wr.eq.0))
                 call oprnd%register_read()
                endif
                if(errc.ne.0) then; errc=-6; exit; endif
               class default
                errc=-5; exit
               end select
              else
               errc=-4; exit
              endif
             enddo
!$OMP END CRITICAL (TAVP_WRK_CACHE)
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
        end function TensInstrMarkIssue
!----------------------------------------------------
        subroutine TensInstrMarkCompletion(this,ierr)
!Updates the tensor access counters for each tensor operand when the tensor instruction is completed.
         implicit none
         class(tens_instr_t), intent(inout):: this   !inout: active tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,i,n
         class(ds_oprnd_t), pointer:: oprnd

         if(.not.this%is_empty(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           n=this%get_num_operands(errc)
           if(errc.eq.DSVP_SUCCESS) then
            if(n.gt.0) then
!$OMP CRITICAL (TAVP_WRK_CACHE)
             do i=0,n-1
              oprnd=>this%get_operand(i,errc)
              if(errc.eq.DSVP_SUCCESS.and.associated(oprnd)) then
               select type(oprnd)
               class is(tens_oprnd_t)
                if(this%operand_is_output(i,errc)) then !output tensor operand
                 call oprnd%unregister_write()
                else !input tensor operand
                 call oprnd%unregister_read()
                endif
                if(errc.ne.0) then; errc=-6; exit; endif
               class default
                errc=-5; exit
               end select
              else
               errc=-4; exit
              endif
             enddo
!$OMP END CRITICAL (TAVP_WRK_CACHE)
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
        end subroutine TensInstrMarkCompletion
!----------------------------------------------------------------------
        function TensInstrDependencyFree(this,ierr,blocked) result(res)
!Returns TRUE if the tensor instruction is dependency-free (locally).
!Argument <blocked> is set to TRUE if at least one of the tensor operands
!is READ/WRITE blocked, that is, there is at least one pending READ and
!at least one pending WRITE registered on it simultaneously.
         implicit none
         logical:: res                               !out: answer
         class(tens_instr_t), intent(in):: this      !in: active tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(out), optional:: blocked    !out: TRUE if at least one of the tensor operands is READ/WRITE blocked
         integer(INTD):: errc,n,i,rd,wr
         class(ds_oprnd_t), pointer:: oprnd
         logical:: blk

         res=.TRUE.; blk=.FALSE.
         if(.not.this%is_empty(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           n=this%get_num_operands(errc)
           if(errc.eq.DSVP_SUCCESS) then
            if(n.gt.0) then
             do i=0,n-1
              oprnd=>this%get_operand(i,errc)
              if(errc.eq.DSVP_SUCCESS) then
               select type(oprnd)
               class is(tens_oprnd_t)
                rd=oprnd%get_read_count(); wr=oprnd%get_write_count()
                blk=(blk.or.(rd.gt.0.and.wr.gt.0))
                if(this%operand_is_output(i,errc)) then
                 res=(res.and.(.not.(rd.gt.0))) !output tensor operand must not be in-use
                else
                 res=(res.and.(.not.(wr.gt.0))) !input tensor operand must not be currently updated
                endif
                if(errc.ne.0) then; errc=-6; exit; endif
               class default
                errc=-5; exit
               end select
              else
               errc=-4; exit
              endif
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
         if(present(blocked)) blocked=blk
         if(present(ierr)) ierr=errc
         return
        end function TensInstrDependencyFree
!---------------------------------------------------------------
        function TensInstrIsSubstitutable(this,ierr) result(res)
!Returns TRUE if the tensor instruction enables output substitution (rename).
         implicit none
         logical:: res                               !out: result
         class(tens_instr_t), intent(in):: this      !active tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,opcode

         res=.FALSE.
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           opcode=this%get_code(errc)
           if(errc.eq.DSVP_SUCCESS) then
            if(opcode.ge.TAVP_ISA_TENS_FIRST.and.opcode.le.TAVP_ISA_TENS_LAST) then !tensor instruction
             if(this%num_out_oprnds.gt.0) then !output operand(s) exist
              if(opcode.ne.TAVP_INSTR_TENS_CREATE) res=.TRUE. !TENS_CREATE does not require output upload since it creates its output locally
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
         if(present(ierr)) ierr=errc
         return
        end function TensInstrIsSubstitutable
!-----------------------------------------------------------------
        function TensInstrOutputSubstituted(this,ierr) result(res)
!Returns TRUE if the (persistent) output tensor operand is substituted with a temporary tensor.
         implicit none
         logical:: res                               !out: result
         class(tens_instr_t), intent(in):: this      !in: active tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,n,l,i,j,k
         class(ds_oprnd_t), pointer:: oprnd
         class(tens_rcrsv_t), pointer:: tensor
         character(TEREC_MAX_TENS_NAME_LEN):: tname

         res=.FALSE.
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           n=this%get_num_operands(errc)
           if(errc.eq.DSVP_SUCCESS.and.n.gt.0) then
            if(this%num_out_oprnds.eq.1) then
             oprnd=>this%get_operand(this%out_oprnds(0),errc) !`checks only the first output operand
             if(errc.eq.DSVP_SUCCESS) then
              select type(oprnd)
              class is(tens_oprnd_t)
               tensor=>oprnd%get_tensor(errc)
               if(errc.eq.0) then
                call tensor%get_name(tname,l,errc)
                if(errc.eq.TEREC_SUCCESS) then
                 if(l.gt.0) then
                  i=l; do while(i.gt.0); if(tname(i:i).eq.'#') exit; i=i-1; enddo
                  if(i.gt.0) then
                   if(i.lt.l) then
                    k=l-i; j=icharnum(k,tname(i+1:l))
                    if(k.eq.l-i) then
                     if(j.ge.0) then; res=.TRUE.; else; errc=-12; endif
                    else
                     errc=-11
                    endif
                   else
                    errc=-10
                   endif
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
              class default
               errc=-6
              end select
             else
              errc=-5
             endif
            else
             if(this%num_out_oprnds.gt.1) errc=-4 !`cannot deal with more than one output tensor operand
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
        end function TensInstrOutputSubstituted
!------------------------------------------------------------
        subroutine TensInstrPrintIt(this,ierr,dev_id,nspaces)
         implicit none
         class(tens_instr_t), intent(in):: this        !in: tensor instruction
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD), intent(in), optional:: dev_id  !in: output device id
         integer(INTD), intent(in), optional:: nspaces !in: left alignment
         integer(INTD):: errc,devo,nsp,i,n
         class(ds_instr_ctrl_t), pointer:: ctrl
         class(ds_oprnd_t), pointer:: oprnd

         errc=0
         devo=6; if(present(dev_id)) devo=dev_id
         nsp=0; if(present(nspaces)) nsp=nspaces
         do i=1,nsp; write(devo,'(" ")',ADVANCE='NO'); enddo
         write(devo,'("TENSOR INSTRUCTION{")')
         do i=1,nsp+1; write(devo,'(" ")',ADVANCE='NO'); enddo
         write(devo,'("id = ",i11,"; opcode = ",i4,"; stat = ",i6,"; err = ",i11)')&
         &this%get_id(),this%get_code(),this%get_status(errc,i),i
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
         if(present(ierr)) ierr=errc
         return
        end subroutine TensInstrPrintIt
!---------------------------------------
        subroutine tens_instr_dtor(this)
         implicit none
         type(tens_instr_t):: this !inout: empty or retired tensor instruction
         integer(INTD):: sts,errc

         sts=this%get_status(errc)
         if((sts.eq.DS_INSTR_EMPTY.or.sts.eq.DS_INSTR_RETIRED).and.errc.eq.DSVP_SUCCESS) then
          call this%clean(errc)
          if(errc.ne.DSVP_SUCCESS) call quit(errc,'#FATAL(TAVP-WRK:tens_instr_dtor): Tensor instruction destruction failed!')
         else
          if(DEBUG.gt.0)&
          &write(CONS_OUT,'("#FATAL(TAVP-WRK:tens_instr_dtor): TAVP instruction is still active: code = ",i5,", stat = ",i5)')&
          &this%get_code(),this%get_status()
          call quit(-1,'#FATAL(TAVP-WRK:tens_instr_dtor): Attempt to destroy an active TAVP instruction!')
         endif
         return
        end subroutine tens_instr_dtor
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
![tavp_wrk_decoder_t]=====================================
        subroutine TAVPWRKDecoderConfigure(this,conf,ierr)
!Configures this DSVU.
         implicit none
         class(tavp_wrk_decoder_t), intent(inout):: this !out: configured DSVU (must not be configured on entrance)
         class(dsv_conf_t), intent(in):: conf            !in: specific DSVU configuration
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         errc=0
         select type(conf)
         type is(tavp_wrk_decoder_conf_t)
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
          write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Decoder configured: ",i11,1x,i6)') impir,this%source_comm,this%source_rank
          flush(CONS_OUT)
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKDecoderConfigure
!------------------------------------------------
        subroutine TAVPWRKDecoderStart(this,ierr)
!Starts and lives this DSVU, calls .shutdown() at the end.
         implicit none
         class(tavp_wrk_decoder_t), intent(inout):: this !inout: TAVP-WRK decoder DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier,thid,num_packets,opcode,sts,port_id,i,j
         logical:: active,stopping,new
         class(dsvp_t), pointer:: dsvp
         class(tavp_wrk_t), pointer:: tavp
         class(*), pointer:: uptr
         class(ds_unit_t), pointer:: acceptor
         type(obj_pack_t):: instr_packet
         type(comm_handle_t):: comm_hl
         type(tens_instr_t), pointer:: tens_instr
         type(tens_instr_t):: tens_instr_empty

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Decoder started as DSVU # ",i2," (thread ",i2,"): Listening to ",i11,1x,i6)')&
          &impir,this%get_id(),thid,this%source_comm,this%source_rank
          flush(CONS_OUT)
         endif
!Reserve a bytecode buffer:
         call this%bytecode%reserve_mem(ier,MAX_BYTECODE_SIZE,MAX_BYTECODE_INSTR); if(ier.ne.0.and.errc.eq.0) errc=-41
!Initialize queues and ports:
         call this%init_queue(this%num_ports,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-40
!Initialize the control list:
         ier=this%ctrl_list%init(this%control_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-39
!Set up tensor argument cache and wait on other TAVP units:
         tavp=>NULL(); dsvp=>this%get_dsvp(); select type(dsvp); class is(tavp_wrk_t); tavp=>dsvp; end select
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
           new=this%bytecode%receive(comm_hl,ier,proc_rank=this%source_rank,tag=TAVP_DISPATCH_TAG,comm=this%source_comm)
           if(ier.ne.0.and.errc.eq.0) then; errc=-35; exit wloop; endif
           if(new) then !new bytecode is available
            call comm_hl%wait(ier); if(ier.ne.0.and.errc.eq.0) then; errc=-34; exit wloop; endif
            num_packets=this%bytecode%get_num_packets(ier); if(ier.ne.0.and.errc.eq.0) then; errc=-33; exit wloop; endif
            if(DEBUG.gt.0) then
             write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Decoder unit ",i2," received ",i9," new instructions")')&
             &impir,this%get_id(),num_packets
             flush(CONS_OUT)
            endif
            if(num_packets.gt.0) then
 !Decode new bytecode:
             ier=this%ctrl_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-32; exit wloop; endif
             ier=this%iqueue%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-31; exit wloop; endif
             do i=1,num_packets
  !Extract an instruction:
              call this%bytecode%extract_packet(i,instr_packet,ier,preclean=.TRUE.)
              if(ier.ne.0.and.errc.eq.0) then; errc=-30; exit wloop; endif
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
  !Clone CONTROL instructions for own port:
              if(opcode.ge.TAVP_ISA_CTRL_FIRST.and.opcode.le.TAVP_ISA_CTRL_LAST) then !copy control instructions to own port
               if(opcode.eq.TAVP_INSTR_CTRL_STOP.and.i.ne.num_packets.and.errc.eq.0) then; errc=-22; exit wloop; endif !STOP must be the last instruction in the packet
               ier=this%ctrl_list%append(tens_instr); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-21; exit wloop; endif
              endif
             enddo
 !Pass cloned control instructions to own port:
             ier=this%ctrl_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-20; exit wloop; endif
             if(this%ctrl_list%get_status().eq.GFC_IT_ACTIVE) then
              ier=this%load_port(0,this%ctrl_list); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-19; exit wloop; endif
             endif
 !Pass all decoded instructions to the acceptor unit (Resourcer in this case):
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
           if(opcode.eq.TAVP_INSTR_CTRL_STOP) then !only STOP instruction is expected
            stopping=.TRUE.
           else
            if(opcode.ne.TAVP_INSTR_CTRL_RESUME) then !`RESUME currently does nothing
             if(errc.eq.0) then; errc=-4; exit wloop; endif
            endif
           endif
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
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(TAVP_WRK)[",i6,"]: Decoder error ",i11," by thread ",i2)')&
         &impir,errc,thid
!Shutdown:
         call this%shutdown(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKDecoderStart
!---------------------------------------------------
        subroutine TAVPWRKDecoderShutdown(this,ierr)
!Stops DSVU (returns back a clean configured state).
         implicit none
         class(tavp_wrk_decoder_t), intent(inout):: this !inout: TAVP-WRK decoder DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier,thid

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Decoder stopped as DSVU # ",i2," (thread ",i2,")")') impir,this%get_id(),thid
          write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Decoder DSVU # ",i2," (thread ",i2,"): Time stats (sec): ")',ADVANCE='NO')&
          &impir,this%get_id(),thid; call this%print_timing(CONS_OUT)
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
!Release queueus:
         call this%release_queue(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-2
!Release the bytecode buffer:
         call this%bytecode%destroy(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
!Record an error, if any:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKDecoderShutdown
!-----------------------------------------------------------------------
        subroutine TAVPWRKDecoderDecode(this,ds_instr,instr_packet,ierr)
!Decodes a tensor instruction from a plain bytecode.
!The bytecode format is specified in tens_instr_t.encode().
         implicit none
         class(tavp_wrk_decoder_t), intent(inout):: this     !inout: TAVP-WRK decoder
         class(ds_instr_t), intent(inout), target:: ds_instr !out: decoded tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet     !in: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr         !out: error code
         integer(INTD):: errc,op_code,stat,err_code
         integer(INTL):: iid
         real(8):: tm

         tm=thread_wtime()
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
                 case(TAVP_INSTR_CTRL_RESUME,TAVP_INSTR_CTRL_STOP)
                 case(TAVP_INSTR_TENS_CREATE,TAVP_INSTR_TENS_DESTROY)
                  call decode_instr_tens_create_destroy(errc); if(errc.ne.0) errc=-12
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
          write(CONS_OUT,'("#ERROR(TAVP-WRK:Decoder.decode)[",i6,":",i3,"]: Error ",i11)')&
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
           class(tens_entry_wrk_t), pointer:: tens_wrk_entry
           integer(INTD):: jj,jn,jown,jread,jwrite
           logical:: stored,updated

           jn=ds_instr%get_num_operands(jerr)
           if(jerr.eq.DSVP_SUCCESS.and.jn.ge.0) then
            tensor_tmp=>NULL()
            do jj=0,jn-1 !loop over tensor operands
             allocate(tensor_tmp,STAT=jerr) !tensor will either be owned by a tensor cache entry or deallocated here
             if(jerr.ne.0) then
              tensor_tmp=>NULL()
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE); jerr=-11; exit
             endif
             call unpack_builtin(instr_packet,jown,jerr) !tensor owner id (not used in TAVP-WRK)
             if(jerr.ne.PACK_SUCCESS) then; call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_BTC_BAD); jerr=-10; exit; endif
             call unpack_builtin(instr_packet,jread,jerr) !tensor read access count (ignored in TAVP-WRK)
             if(jerr.ne.PACK_SUCCESS) then; call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_BTC_BAD); jerr=-10; exit; endif
             call unpack_builtin(instr_packet,jwrite,jerr) !tensor write access count (ignored in TAVP-WRK)
             if(jerr.ne.PACK_SUCCESS) then; call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_BTC_BAD); jerr=-10; exit; endif
             call tensor_tmp%tens_rcrsv_ctor(instr_packet,jerr) !unpack tensor information into a temporary tensor
             if(jerr.ne.TEREC_SUCCESS) then; call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_BTC_BAD); jerr=-9; exit; endif
             tens_entry=>NULL(); tens_entry=>this%arg_cache%lookup(tensor_tmp,jerr)
             if(jerr.ne.0) then
              if(associated(tens_entry)) call this%arg_cache%release_entry(tens_entry)
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_CHE_FAILURE); jerr=-8; exit
             endif
             if(.not.associated(tens_entry)) then !tensor is absent in the tensor cache: Create an entry for it
              stored=this%arg_cache%store(tensor_tmp,tens_entry_wrk_alloc,jerr,tens_entry_p=tens_entry) !tensor ownership is moved to the tensor cache entry
              if((.not.stored).or.(jerr.ne.0).or.(.not.associated(tens_entry))) then
               if(associated(tens_entry)) call this%arg_cache%release_entry(tens_entry)
               call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_CHE_FAILURE); jerr=-7; exit
              else
               tensor_tmp=>NULL() !tensor ownership has been transferred to the tensor cache
              endif
             else
              stored=.FALSE.
             endif
             updated=stored
             tens_wrk_entry=>NULL(); select type(tens_entry); class is(tens_entry_wrk_t); tens_wrk_entry=>tens_entry; end select
             if(.not.associated(tens_wrk_entry)) then !trap
              if(associated(tens_entry)) call this%arg_cache%release_entry(tens_entry)
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_CHE_FAILURE); jerr=-6; exit
             endif
             tensor=>tens_wrk_entry%get_tensor() !use the tensor from the tensor cache
             if(.not.stored) then !the tensor was already in the tensor cache before, update it by the information from the just decoded tensor
!$OMP CRITICAL (TAVP_WRK_CACHE)
              call tensor%update(tensor_tmp,jerr,updated) !tensor metadata update is done inside the tensor cache
!$OMP END CRITICAL (TAVP_WRK_CACHE)
              deallocate(tensor_tmp); tensor_tmp=>NULL() !deallocate temporary tensor after importing its information into the cache
              if(jerr.ne.TEREC_SUCCESS) then
               if(associated(tens_entry)) call this%arg_cache%release_entry(tens_entry)
               call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE); jerr=-5; exit
              endif
             endif
             allocate(tens_oprnd,STAT=jerr) !tensor operand will be owned by the tensor instruction
             if(jerr.ne.0) then
              tens_oprnd=>NULL()
              if(associated(tens_entry)) call this%arg_cache%release_entry(tens_entry)
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_RSC_UNAVAILABLE); jerr=-4; exit
             endif
             call tens_oprnd%tens_oprnd_ctor(tens_wrk_entry,jerr) !tensor and tensor resource are owned by the tensor cache
             if(jerr.ne.0) then
              deallocate(tens_oprnd); tens_oprnd=>NULL()
              if(associated(tens_entry)) call this%arg_cache%release_entry(tens_entry)
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE); jerr=-3; exit
             endif
             oprnd=>tens_oprnd; call ds_instr%set_operand(jj,oprnd,jerr) !tensor operand ownership is moved to the tensor instruction
             if(jerr.ne.DSVP_SUCCESS) then
              deallocate(tens_oprnd); tens_oprnd=>NULL()
              if(associated(tens_entry)) call this%arg_cache%release_entry(tens_entry)
              call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE); jerr=-2; exit
             endif
             tens_oprnd=>NULL()
             if(associated(tens_entry)) call this%arg_cache%release_entry(tens_entry); tens_entry=>NULL()
            enddo !loop over tensor operands
            if(associated(tensor_tmp)) deallocate(tensor_tmp) !deallocate the temporary tensor (in case of error)
           else
            call ds_instr%set_status(DS_INSTR_RETIRED,jerr,TAVP_ERR_GEN_FAILURE); jerr=-1
           endif
           if(DEBUG.gt.0.and.jerr.ne.0) then
            write(CONS_OUT,'("#ERROR(TAVP-WRK:Decoder.decode.decode_instr_operands)[",i6,":",i3,"]: Error ",i11)')&
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
            if(jerr.eq.0) then
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
            write(CONS_OUT,'("#ERROR(TAVP-WRK:Decoder.decode.decode_instr_tens_create_destroy)[",i6,":",i3,"]: Error ",i11)')&
            &impir,omp_get_thread_num(),jerr
            flush(CONS_OUT)
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
               if(jerr.eq.0) then
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
           if(DEBUG.gt.0.and.jerr.ne.0) then
            write(CONS_OUT,'("#ERROR(TAVP-WRK:Decoder.decode.decode_instr_tens_contract)[",i6,":",i3,"]: Error ",i11)')&
            &impir,omp_get_thread_num(),jerr
            flush(CONS_OUT)
           endif
           return
          end subroutine decode_instr_tens_contract

        end subroutine TAVPWRKDecoderDecode
![tavp_wrk_retirer_t]=====================================
        subroutine TAVPWRKRetirerConfigure(this,conf,ierr)
!Configures this DSVU.
         implicit none
         class(tavp_wrk_retirer_t), intent(inout):: this !out: configured DSVU (must not be configured on entrance)
         class(dsv_conf_t), intent(in):: conf            !in: specific DSVU configuration
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,i,n

         errc=0
         select type(conf)
         type is(tavp_wrk_retirer_conf_t)
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
        end subroutine TAVPWRKRetirerConfigure
!------------------------------------------------
        subroutine TAVPWRKRetirerStart(this,ierr)
!Starts and lives this DSVU, calls .shutdown() at the end.
         implicit none
         class(tavp_wrk_retirer_t), intent(inout):: this !inout: TAVP-WRK retirer DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier,thid
         logical:: active,stopping
         class(dsvp_t), pointer:: dsvp
         class(tavp_wrk_t), pointer:: tavp

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Retirer started as DSVU # ",i2," (thread ",i2,"): Reporting to ",i11,1x,i6)')&
          &impir,this%get_id(),thid,this%retire_comm,this%retire_rank
          flush(CONS_OUT)
         endif
!Reserve a bytecode buffer:
         call this%bytecode%reserve_mem(ier,MAX_BYTECODE_SIZE,MAX_BYTECODE_INSTR); if(ier.ne.0.and.errc.eq.0) errc=-1
!Initialize queues and ports:
         call this%init_queue(this%num_ports,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
!Set up tensor argument cache and wait on other TAVP units:
         tavp=>NULL(); dsvp=>this%get_dsvp(); select type(dsvp); class is(tavp_wrk_t); tavp=>dsvp; end select
         if(associated(tavp)) then
          this%arg_cache=>tavp%tens_cache
          call tavp%sync_units(errc,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
         else
          this%arg_cache=>NULL(); if(errc.eq.0) errc=-1
         endif
!Work loop:
         active=((errc.eq.0).and.(this%retire_comm.ne.MPI_COMM_NULL)); stopping=(.not.active)
         wloop: do while(active)
          exit wloop
         enddo wloop
!Record the error:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(TAVP-WRK)[",i6,"]: Retirer error ",i11," by thread ",i2)')&
         &impir,errc,thid
!Shutdown:
         call this%shutdown(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKRetirerStart
!---------------------------------------------------
        subroutine TAVPWRKRetirerShutdown(this,ierr)
!Stops DSVU (returns back a clean configured state).
         implicit none
         class(tavp_wrk_retirer_t), intent(inout):: this !inout: TAVP-WRK retirer DSVU
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier,thid

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Retirer stopped as DSVU # ",i2," (thread ",i2,")")') impir,this%get_id(),thid
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
        end subroutine TAVPWRKRetirerShutdown
!-----------------------------------------------------------------------
        subroutine TAVPWRKRetirerEncode(this,ds_instr,instr_packet,ierr)
!Encodes a tensor instruction into a plain bytecode.
         implicit none
         class(tavp_wrk_retirer_t), intent(inout):: this  !inout: retirer
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
        end subroutine TAVPWRKRetirerEncode
![tavp_wrk_resourcer_t]=====================================
        subroutine TAVPWRKResourcerConfigure(this,conf,ierr)
!Configures this DSVU.
         implicit none
         class(tavp_wrk_resourcer_t), intent(inout):: this !out: configured DSVU (must not be configured on entrance)
         class(dsv_conf_t), intent(in):: conf              !in: specific DSVU configuration
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc

         errc=0
         select type(conf)
         type is(tavp_wrk_resourcer_conf_t)
          if(conf%host_ram_size.ge.0.and.conf%nvram_size.ge.0) then
           this%host_ram_size=conf%host_ram_size
           this%nvram_size=conf%nvram_size
          else
           errc=-2
          endif
         class default
          errc=-1
         end select
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKResourcerConfigure
!--------------------------------------------------
        subroutine TAVPWRKResourcerStart(this,ierr)
!Starts and lives this DSVU, calls .shutdown() at the end.
         implicit none
         class(tavp_wrk_resourcer_t), intent(inout):: this !inout: TAVP-WRK resourcer DSVU
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc,ier,thid,n,num_staged,opcode,sts,errcode
         integer:: rsc_timer
         logical:: active,stopping,auxiliary,deferd,dependent,blocked,passed,expired
         type(tens_instr_t):: instr_fence
         class(tens_instr_t), pointer:: instr
         class(dsvp_t), pointer:: dsvp
         class(tavp_wrk_t), pointer:: tavp
         class(*), pointer:: uptr

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Resourcer started as DSVU # ",i2," (thread ",i2,"): Max memory (B) = ",i15)')&
          impir,this%get_id(),thid,this%host_ram_size
          flush(CONS_OUT)
         endif
!Initialize queues and ports:
         call this%init_queue(this%num_ports,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-80
!Initialize the staged list:
         ier=this%stg_list%init(this%staged_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-79
!Initialize the deferred list:
         ier=this%def_list%init(this%deferred_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-78
!Initialize the release list:
         ier=this%rls_list%init(this%release_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-77
!Create a special FENCE instruction:
         call instr_fence%tens_instr_ctor(TAVP_INSTR_CTRL_RESUME,ier,iid=0_INTL,stat=DS_INSTR_SPECIAL)
         if(ier.ne.0.and.errc.eq.0) errc=-76
!Set up tensor argument cache and wait on other TAVP units:
         tavp=>NULL(); dsvp=>this%get_dsvp(); select type(dsvp); class is(tavp_wrk_t); tavp=>dsvp; end select
         if(associated(tavp)) then
          this%arg_cache=>tavp%tens_cache
          call tavp%sync_units(errc,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-75
         else
          this%arg_cache=>NULL(); if(errc.eq.0) errc=-74
         endif
!Work loop:
         ier=timer_start(rsc_timer,MAX_RESOURCER_PHASE_TIME); if(ier.ne.TIMERS_SUCCESS.and.errc.eq.0) errc=-73
         active=(errc.eq.0); stopping=(.not.active); num_staged=0
         wloop: do while(active)
 !Process the deferred queue (check data dependencies and try acquiring resources for input tensor operands again):
          ier=this%def_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-72; exit wloop; endif
          deferd=.FALSE.; ier=this%def_list%get_status()
          dloop: do while(ier.eq.GFC_IT_ACTIVE)
           deferd=.TRUE.
  !Extract a deferred tensor instruction:
           uptr=>this%def_list%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-71; exit wloop; endif
           instr=>NULL(); select type(uptr); class is(tens_instr_t); instr=>uptr; end select
           if((.not.associated(instr)).and.errc.eq.0) then; errc=-70; exit wloop; endif !trap
           opcode=instr%get_code(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-69; exit wloop; endif
           sts=instr%get_status(ier,errcode); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-68; exit wloop; endif
           if(sts.ne.DS_INSTR_RSC_WAIT.and.errc.eq.0) then; errc=-67; exit wloop; endif !trap
  !Process the deferred tensor instruction:
           if(opcode.ge.TAVP_ISA_TENS_FIRST.and.opcode.le.TAVP_ISA_TENS_LAST) then !tensor instruction
   !Check tensor instruction data dependencies:
            dependent=.not.instr%dependency_free(ier,blocked); if(ier.ne.0.and.errc.eq.0) then; errc=-66; exit wloop; endif
            if(.not.dependent) then
   !Acquire resources for tensor operands (if available):
             call this%acquire_resources(instr,ier,omit_output=.FALSE.)
             if(ier.eq.0) then !resources have been acquired: issue into the staged list
              call instr%set_status(DS_INSTR_INPUT_WAIT,ier)
              if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-65; exit wloop; endif
              ier=this%def_list%move_elem(this%stg_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-64; exit wloop; endif
              num_staged=num_staged+1
             elseif(ier.eq.TRY_LATER) then !required resources are still unavailable
              ier=this%def_list%next()
             else
              if(errc.eq.0) then; errc=-63; exit wloop; endif
             endif
            else
             ier=this%def_list%next()
            endif
           else
            if(errc.eq.0) then; errc=-62; exit wloop; endif
           endif
           ier=this%def_list%get_status()
          enddo dloop
          if(ier.ne.GFC_IT_EMPTY.and.ier.ne.GFC_IT_DONE.and.errc.eq.0) then; errc=-61; exit wloop; endif !trap
 !Get new instructions from Decoder (port 0) into the main queue:
          ier=this%iqueue%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-60; exit wloop; endif
          ier=this%flush_port(0,num_moved=n); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-59; exit wloop; endif
          if(DEBUG.gt.0.and.n.gt.0) then
           write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Resourcer unit ",i2," received ",i9," instructions from Decoder")')&
           &impir,this%get_id(),n
           !ier=this%iqueue%reset(); ier=this%iqueue%scanp(action_f=tens_instr_print); ier=this%iqueue%reset_back() !print all instructions
           flush(CONS_OUT)
          endif
 !Process the main queue (rename output tensor operands, check data dependencies, and acquire resources for input tensor operands):
          ier=this%iqueue%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-58; exit wloop; endif
          auxiliary=.FALSE.; ier=this%iqueue%get_status()
          if(stopping.and.ier.ne.GFC_IT_EMPTY.and.errc.eq.0) then; errc=-57; exit wloop; endif !if stopping, no further instructions are expected
          mloop: do while(ier.eq.GFC_IT_ACTIVE)
  !Extract an instruction:
           uptr=>this%iqueue%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-56; exit wloop; endif
           instr=>NULL(); select type(uptr); class is(tens_instr_t); instr=>uptr; end select
           if((.not.associated(instr)).and.errc.eq.0) then; errc=-55; exit wloop; endif !trap
           opcode=instr%get_code(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-54; exit wloop; endif
           sts=instr%get_status(ier,errcode); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-53; exit wloop; endif
           if(sts.ne.DS_INSTR_NEW.and.errc.eq.0) then; errc=-52; exit wloop; endif !trap
           call instr%set_status(DS_INSTR_RSC_WAIT,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-51; exit wloop; endif
  !Process the instruction according to its category:
           if(opcode.ge.TAVP_ISA_TENS_FIRST.and.opcode.le.TAVP_ISA_TENS_LAST) then !tensor instruction
            if(.not.stopping) then
             if(auxiliary) then !auxiliary instructions stall the pipeline
              call instr%set_status(DS_INSTR_NEW,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-50; exit wloop; endif
              ier=this%iqueue%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-49; exit wloop; endif
              ier=this%iqueue%next(); ier=0 !to stall the pipeline
             else !pipeline is free
   !Substitute (rename) the output tensor operand with a temporary tensor for numerical tensor operations (for concurrency):
              if(instr%is_substitutable(ier)) then
               if(ier.eq.0) then
                if(.not.instr%output_substituted(ier)) then
                 if(ier.eq.0) then
                  call this%substitute_output(instr,ier); if(ier.ne.0.and.errc.eq.0) then; errc=-48; exit wloop; endif
                 else
                  if(errc.eq.0) then; errc=-47; exit wloop; endif
                 endif
                else
                 if(ier.ne.0.and.errc.eq.0) then; errc=-46; exit wloop; endif
                endif
               else
                if(errc.eq.0) then; errc=-45; exit wloop; endif
               endif
              else
               if(ier.ne.0.and.errc.eq.0) then; errc=-44; exit wloop; endif
              endif
   !Set up the storage layout for non-existing output tensor operands (if any):
              call instr%lay_output_operands(ier); if(ier.ne.0.and.errc.eq.0) then; errc=-43; exit wloop; endif
   !Check tensor instruction data dependencies:
              dependent=.not.instr%dependency_free(ier,blocked); if(ier.ne.0.and.errc.eq.0) then; errc=-42; exit wloop; endif
   !Update data dependencies:
              if(dependent) then !at least one tensor operand has a simple data dependency (read_count > 0 for WRITE or write_count > 0 for READ)
               if(blocked) then !at least one tensor operand has blocking data dependency (read_count > 0 and write_count > 0)
                call instr%set_status(DS_INSTR_NEW,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-41; exit wloop; endif
                ier=this%iqueue%next()
                if(ier.ne.GFC_SUCCESS.and.ier.ne.GFC_NO_MOVE.and.errc.eq.0) then; errc=-40; exit wloop; endif
               else !no blocking data dependencies: issue into the deferred instruction list
                passed=instr%mark_issue(ier); if((ier.ne.0.or.passed).and.errc.eq.0) then; errc=-39; exit wloop; endif
                ier=this%iqueue%move_elem(this%def_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-38; exit wloop; endif
               endif
              else !no data dependencies
               passed=instr%mark_issue(ier); if((.not.(ier.eq.0.and.passed)).and.errc.eq.0) then; errc=-37; exit wloop; endif
   !Acquire resources for tensor operands (if available):
               call this%acquire_resources(instr,ier,omit_output=.FALSE.)
               if(ier.eq.0) then !resources have been acquired: issue into the staged list
                call instr%set_status(DS_INSTR_INPUT_WAIT,ier)
                if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-36; exit wloop; endif
                ier=this%iqueue%move_elem(this%stg_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-35; exit wloop; endif
                num_staged=num_staged+1
               elseif(ier.eq.TRY_LATER) then !required resources are not currently available: issue into the deferred list
                ier=this%iqueue%move_elem(this%def_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-34; exit wloop; endif
                ier=0
               else
                if(errc.eq.0) then; errc=-33; exit wloop; endif
               endif
              endif
             endif
            else
             if(errc.eq.0) then; errc=-32; exit wloop; endif
            endif
           elseif(opcode.ge.TAVP_ISA_CTRL_FIRST.and.opcode.le.TAVP_ISA_CTRL_LAST) then !control instruction: no dependencies
            if(.not.stopping) then
             auxiliary=.TRUE.
             if(opcode.eq.TAVP_INSTR_CTRL_STOP) stopping=.TRUE.
             call instr%set_status(DS_INSTR_INPUT_WAIT,ier)
             if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-31; exit wloop; endif
             ier=this%iqueue%move_elem(this%stg_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-30; exit wloop; endif
             num_staged=num_staged+1
            else
             if(errc.eq.0) then; errc=-29; exit wloop; endif
            endif
           else !auixiliary instruction: no dependencies
            if(.not.stopping) then
             auxiliary=.TRUE.
             call instr%set_status(DS_INSTR_INPUT_WAIT,ier)
             if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-28; exit wloop; endif
             ier=this%iqueue%move_elem(this%stg_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-27; exit wloop; endif
             num_staged=num_staged+1
            else
             if(errc.eq.0) then; errc=-26; exit wloop; endif
            endif
           endif
  !Pass periodically staged instructions to Communicator Port 0:
           expired=timer_expired(rsc_timer,ier); if(ier.ne.TIMERS_SUCCESS.and.errc.eq.0) then; errc=-25; exit wloop; endif
           if(expired.or.auxiliary.or.num_staged.gt.MAX_RESOURCER_INSTR) then
            ier=this%stg_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-24; exit wloop; endif
            if(this%stg_list%get_status().eq.GFC_IT_ACTIVE) then
             ier=tavp%communicator%load_port(0,this%stg_list,num_moved=n)
             if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-23; exit wloop; endif
             if(DEBUG.gt.0.and.n.gt.0) then
              write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Resourcer unit ",i2," passed ",i9," instructions to Communicator")')&
              &impir,this%get_id(),n
              flush(CONS_OUT)
             endif
             num_staged=0
            endif
            if(expired) ier=timer_reset(rsc_timer,MAX_RESOURCER_PHASE_TIME)
            auxiliary=.FALSE.
           endif
           ier=this%iqueue%get_status()
          enddo mloop
          if(ier.ne.GFC_IT_EMPTY.and.ier.ne.GFC_IT_DONE.and.errc.eq.0) then; errc=-22; exit wloop; endif !trap
 !Pass the remaining staged instructions to Communicator Port 0:
          ier=this%stg_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-21; exit wloop; endif
          if(this%stg_list%get_status().eq.GFC_IT_ACTIVE) then
           ier=tavp%communicator%load_port(0,this%stg_list,num_moved=n)
           if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-20; exit wloop; endif
           if(DEBUG.gt.0.and.n.gt.0) then
            write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Resourcer unit ",i2," passed ",i9," instructions to Communicator")')&
            &impir,this%get_id(),n
            flush(CONS_OUT)
           endif
           num_staged=0
          endif
 !Process retired instructions from Communicator (release resources):
  !Get instructions for resource release from Communicator (port 1) into the release queue:
          ier=this%rls_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-19; exit wloop; endif
          ier=this%unload_port(1,this%rls_list,num_moved=n); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-18; exit wloop; endif
          if(DEBUG.gt.0.and.n.gt.0) then
           write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Resourcer unit ",i2," received ",i9," instructions from Communicator")')&
           &impir,this%get_id(),n
           !ier=this%rls_list%reset(); ier=this%rls_list%scanp(action_f=tens_instr_print); ier=this%rls_list%reset_back() !print all instructions
           flush(CONS_OUT)
          endif
  !Release resources and rename output tensor operand for retired instructions:
          ier=this%rls_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-17; exit wloop; endif
          if(this%rls_list%get_status().eq.GFC_IT_ACTIVE) then
           ier=GFC_SUCCESS
           rloop: do while(ier.eq.GFC_SUCCESS)
            uptr=>this%rls_list%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-16; exit wloop; endif
            instr=>NULL(); select type(uptr); class is(tens_instr_t); instr=>uptr; end select
            if((.not.associated(instr)).and.errc.eq.0) then; errc=-15; exit wloop; endif !trap
            opcode=instr%get_code(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-14; exit wloop; endif
            sts=instr%get_status(ier,errcode); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-13; exit wloop; endif
            if(sts.ne.DS_INSTR_RETIRED.and.errc.eq.0) then; errc=-12; exit wloop; endif !trap
            if(opcode.ge.TAVP_ISA_TENS_FIRST.and.opcode.le.TAVP_ISA_TENS_LAST) then !tensor instruction
             call this%release_resources(instr,ier); if(ier.ne.0.and.errc.eq.0) then; errc=-11; exit wloop; endif
             if(instr%output_substituted(ier)) then
              if(ier.eq.0) then
               call this%restore_output(instr,ier); if(ier.ne.0.and.errc.eq.0) then; errc=-10; exit wloop; endif
              else
               if(errc.eq.0) then; errc=-9; exit wloop; endif
              endif
             else
              if(ier.ne.0.and.errc.eq.0) then; errc=-8; exit wloop; endif
             endif
            else
             if(errc.eq.0) then; errc=-7; exit wloop; endif
            endif
            ier=this%rls_list%next()
           enddo rloop
  !Pass tensor instructions after resource release to Retirer (port 0):
           if(ier.eq.GFC_NO_MOVE) then
            ier=this%rls_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-6; exit wloop; endif
            ier=tavp%retirer%load_port(0,this%rls_list); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-5; exit wloop; endif
           else
            if(errc.eq.0) then; errc=-4; exit wloop; endif
           endif
          endif
          active=((.not.stopping).or.deferd)
         enddo wloop
!Destroy the timer:
         ier=timer_destroy(rsc_timer); if(ier.ne.TIMERS_SUCCESS.and.errc.eq.0) errc=-3
!Retire the special FENCE instruction:
         call instr_fence%set_status(DS_INSTR_RETIRED,ier,DSVP_SUCCESS); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-2
!Record the error:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(TAVP-WRK)[",i6,"]: Resourcer error ",i11," by thread ",i2)')&
         &impir,errc,thid
!Shutdown:
         call this%shutdown(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKResourcerStart
!-----------------------------------------------------
        subroutine TAVPWRKResourcerShutdown(this,ierr)
!Stops DSVU (returns back a clean configured state).
         implicit none
         class(tavp_wrk_resourcer_t), intent(inout):: this !inout: TAVP-WRK resourcer DSVU
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc,ier,thid

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Resourcer stopped as DSVU # ",i2," (thread ",i2,")")') impir,this%get_id(),thid
          flush(CONS_OUT)
         endif
!Release the tensor argument cache pointer:
         this%arg_cache=>NULL()
!Deactivate the release list:
         ier=this%rls_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%rls_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-10
           ier=this%rls_list%delete_all()
          endif
          ier=this%rls_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-9
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
!Deactivate the staged list:
         ier=this%stg_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%stg_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-4
           ier=this%stg_list%delete_all()
          endif
          ier=this%stg_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-3
         else
          if(errc.eq.0) errc=-2
         endif
!Release queues:
         call this%release_queue(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
!Record an error, if any:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKResourcerShutdown
!------------------------------------------------------------------------------------
        subroutine TAVPWRKResourcerAcquireResources(this,tens_instr,ierr,omit_output)
!Acquires local resource for each tensor instruction operand.
!If some resources cannot be acquired now, returns TRY_LATER.
!In that case, the successfully acquired resources will be kept,
!unless an error other than TRY_LATER occurred. If an operand
!already has its resource previously acquired, it will be kept so.
         implicit none
         class(tavp_wrk_resourcer_t), intent(inout):: this !inout: TAVP-WRK resourcer DSVU
         class(tens_instr_t), intent(inout):: tens_instr   !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr       !out: error code
         logical, intent(in), optional:: omit_output       !in: if TRUE, the output operand(s) will be ommitted (defaults to FALSE)
         integer(INTD):: errc,ier,n,l,k
         class(ds_oprnd_t), pointer:: oprnd
         class(tens_rcrsv_t), pointer:: tensor
         class(tens_header_t), pointer:: header
         class(tens_cache_entry_t), pointer:: tens_entry
         type(tens_rcrsv_t):: tens
         logical:: no_output,op_output,subst
         character(TEREC_MAX_TENS_NAME_LEN+2):: tname,oname,aname

         no_output=.FALSE.; if(present(omit_output)) no_output=omit_output
         n=tens_instr%get_num_operands(errc)
         if(errc.eq.DSVP_SUCCESS) then
          aloop: do while(n.gt.0)
           n=n-1; op_output=tens_instr%operand_is_output(n)
           if(op_output.and.no_output) cycle
           oprnd=>tens_instr%get_operand(n,ier)
           if(ier.eq.DSVP_SUCCESS) then
            call oprnd%acquire_rsc(ier)
            if(ier.eq.0) then
             if(op_output) then !output operand may need resource for an accumulator tensor as well
              subst=tens_instr%output_substituted(ier); if(ier.ne.0) then; errc=-17; exit aloop; endif !`assumes only one output tensor operand
              if(subst) then !if output operand is substituted, its accumulator tensor needs resources
               select type(oprnd)
               class is(tens_oprnd_t)
                tensor=>oprnd%get_tensor(ier); if(ier.ne.0) then; errc=-16; exit aloop; endif
                header=>tensor%get_header(ier); if(ier.ne.TEREC_SUCCESS) then; errc=-15; exit aloop; endif
                call tens%tens_rcrsv_ctor(header,ier); if(ier.ne.TEREC_SUCCESS) then; errc=-14; exit aloop; endif
                call tens%get_name(tname,l,ier); if(ier.ne.TEREC_SUCCESS) then; errc=-13; exit aloop; endif
                call tensor_name_unmangle_temporary(tname(1:l),oname,k,ier); if(ier.ne.0) then; errc=-12; exit aloop; endif
                call tensor_name_mangle_temporary(oname(1:k),aname,l,ier,0); if(ier.ne.0) then; errc=-11; exit aloop; endif
                call tens%rename(aname(1:l),ier); if(ier.ne.TEREC_SUCCESS) then; errc=-10; exit aloop; endif
                tens_entry=>this%arg_cache%lookup(tens,ier); if(ier.ne.0) then; errc=-9; exit aloop; endif
                if(associated(tens_entry)) then
                 select type(tens_entry)
                 class is(tens_entry_wrk_t)
                  call tens_entry%acquire_resource(ier)
                  if(ier.ne.0) then
                   if(ier.eq.TRY_LATER) then
                    errc=ier
                   else
                    write(CONS_OUT,'("#ERROR(TAVP-WRK:Resourcer.acquire_resources)[",i6,"]: Severe failure for operand # ",i2,'//&
                    &'": Error ",i11)') impir,n,ier
                    errc=-8; exit aloop
                   endif
                  endif
                 class default
                  errc=-7; exit aloop
                 end select
                 call this%arg_cache%release_entry(tens_entry,ier); if(ier.ne.0) then; errc=-6; exit aloop; endif
                else
                 errc=-5; exit aloop
                endif
               class default
                errc=-4; exit aloop
               end select
              endif
             endif
            else
             if(ier.eq.TRY_LATER) then
              errc=ier
             else
              write(CONS_OUT,'("#ERROR(TAVP-WRK:Resourcer.acquire_resources)[",i6,"]: Severe failure for operand # ",i2,'//&
              &'": Error ",i11)') impir,n,ier
              flush(CONS_OUT)
              errc=-3; exit aloop
             endif
            endif
           else
            errc=-2; exit aloop
           endif
          enddo aloop
          if(errc.ne.0.and.errc.ne.TRY_LATER) call this%release_resources(tens_instr,ier)
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKResourcerAcquireResources
!------------------------------------------------------------------------
        subroutine TAVPWRKResourcerReleaseResources(this,tens_instr,ierr)
!Releases local resources occupied by the tensor instruction operands,
!but the operands stay defined.
         implicit none
         class(tavp_wrk_resourcer_t), intent(inout):: this !inout: TAVP-WRK resourcer DSVU
         class(tens_instr_t), intent(inout):: tens_instr   !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc,ier,n,l,k
         type(tens_rcrsv_t):: tens
         class(ds_oprnd_t), pointer:: oprnd
         class(tens_rcrsv_t), pointer:: tensor
         class(tens_header_t), pointer:: header
         class(tens_cache_entry_t), pointer:: tens_entry
         character(TEREC_MAX_TENS_NAME_LEN+2):: tname,oname,aname

         n=tens_instr%get_num_operands(errc)
         if(errc.eq.DSVP_SUCCESS) then
          rloop: do while(n.gt.0)
           n=n-1; oprnd=>tens_instr%get_operand(n,ier)
           if(ier.eq.DSVP_SUCCESS) then
            call oprnd%release(ier)
            if(ier.eq.0) then
             if(tens_instr%operand_is_output(n)) then !the associated accumulator tensor may need to be released as well
              if(tens_instr%output_substituted(ier)) then
               if(ier.eq.0) then
                select type(oprnd)
                class is(tens_oprnd_t)
                 tensor=>oprnd%get_tensor(ier); if(ier.ne.0) then; errc=-18; exit rloop; endif
                 header=>tensor%get_header(ier); if(ier.ne.TEREC_SUCCESS) then; errc=-17; exit rloop; endif
                 call tens%tens_rcrsv_ctor(header,ier); if(ier.ne.TEREC_SUCCESS) then; errc=-16; exit rloop; endif
                 call tens%get_name(tname,l,ier); if(ier.ne.TEREC_SUCCESS) then; errc=-15; exit rloop; endif
                 call tensor_name_unmangle_temporary(tname(1:l),oname,k,ier); if(ier.ne.0) then; errc=-14; exit rloop; endif
                 call tensor_name_mangle_temporary(oname(1:k),aname,l,ier,0); if(ier.ne.0) then; errc=-13; exit rloop; endif
                 call tens%rename(aname(1:l),ier); if(ier.ne.TEREC_SUCCESS) then; errc=-12; exit rloop; endif
                 tens_entry=>this%arg_cache%lookup(tens,ier); if(ier.ne.0) then; errc=-11; exit rloop; endif
                 if(associated(tens_entry)) then
                  select type(tens_entry)
                  class is(tens_entry_wrk_t)
                   call tens_entry%release_resource(ier)
                   if(ier.ne.0) then
                    write(CONS_OUT,'("#ERROR(TAVP-WRK:Resourcer.release_resources)[",i6,'//&
                    &'"]: Resource release failure for operand-accumulator # ",i2)') impir,n
                    flush(CONS_OUT)
                    if(errc.eq.0) errc=-10
                   endif
                  class default
                   errc=-9; exit rloop
                  end select
                 else
                  errc=-8; exit rloop
                 endif
                class default
                 errc=-7; exit rloop
                end select
               else
                errc=-6; exit rloop
               endif
              else
               if(ier.ne.0) then; errc=-5; exit rloop; endif
              endif
             endif
            else
             write(CONS_OUT,'("#ERROR(TAVP-WRK:Resourcer.release_resources)[",i6,'//&
             &'"]: Resource release failure for operand # ",i2)') impir,n
             flush(CONS_OUT)
             if(errc.eq.0) errc=-4
            endif
           else
            errc=-3; exit rloop
           endif
          enddo rloop
          ier=talsh_task_destruct(tens_instr%talsh_task) !TAL-SH task has an implicit initialization (no explicit ctor)
          if(ier.ne.TALSH_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKResourcerReleaseResources
!------------------------------------------------------------------------
        subroutine TAVPWRKResourcerSubstituteOutput(this,tens_instr,ierr)
!Substitutes the persistent output tensor in a tensor instruction with a temporary one
!(assumes a single output tensor for now). If this is the first local substitution
!of the given persistent output tensor, an accumulator tensor will be created in the
!tensor cache in addition to the temporary tensor. The persistent tensor still stays
!in the tensor cache as its reference count is not decremented by this procedure.
         implicit none
         class(tavp_wrk_resourcer_t), intent(inout):: this !inout: TAVP-WRK Resourcer
         class(tens_instr_t), intent(inout):: tens_instr   !inout: active tensor instruction
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc,n,nou,tc
         integer(INTD), pointer:: out_oprs(:)
         class(ds_oprnd_t), pointer:: oprnd
         class(tens_entry_wrk_t), pointer:: cache_entry
         class(tens_cache_entry_t), pointer:: new_cache_entry
         class(tens_rcrsv_t), pointer:: tensor
         class(tens_header_t), pointer:: header

         if(tens_instr%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           n=tens_instr%get_num_operands(errc)
           if(errc.eq.DSVP_SUCCESS.and.n.gt.0) then
            out_oprs=>tens_instr%get_output_operands(errc,nou)
            if(errc.eq.0.and.nou.eq.1.and.lbound(out_oprs,1).eq.0) then
             oprnd=>tens_instr%get_operand(out_oprs(0),errc) !output tensor operand
             if(errc.eq.DSVP_SUCCESS) then
              select type(oprnd)
              class is(tens_oprnd_t)
               cache_entry=>oprnd%get_cache_entry(errc)
               if(errc.eq.0.and.associated(cache_entry)) then
                tensor=>oprnd%get_tensor(errc)
                if(errc.eq.0.and.associated(tensor)) then
                 header=>tensor%get_header(errc)
                 if(errc.eq.TEREC_SUCCESS.and.associated(header)) then
 !Mark the persistent output tensor as being written to (before substitution):
                  call cache_entry%incr_write_count()
 !Register the accumulator tensor, if needed:
                  call cache_entry%incr_temp_count() !new temporary tensor to be created
                  tc=cache_entry%get_temp_count()
                  if(tc.eq.1) then !first tensor instruction with this output tensor operand: Register accumulator tensor
                   call register_temp_tensor(0,new_cache_entry,errc) !register accumulator tensor in the cache
                   if(errc.eq.0) then
                    call oprnd%set_tensor_layout(errc); if(errc.ne.0) errc=-14 !set accumulator tensor layout
                    call this%arg_cache%release_entry(new_cache_entry); new_cache_entry=>NULL()
                   else
                    errc=-13
                   endif
                  endif
 !Register the temporary tensor:
                  if(errc.eq.0) then
                   call register_temp_tensor(tc,new_cache_entry,errc) !register a temporary tensor (its layout will be set later)
 !Substitute the persistent output tensor with the temporary tensor (persistent tensor cache entry reference count unchanged):
                   if(errc.eq.0) then
                    cache_entry=>NULL()
                    select type(new_cache_entry); class is(tens_entry_wrk_t); cache_entry=>new_cache_entry; end select
                    if(associated(cache_entry)) then
                     call oprnd%reset_tmp_tensor(cache_entry,.TRUE.,errc); if(errc.ne.0) errc=-12
                     cache_entry=>NULL()
                    else
                     errc=-11
                    endif
                    call this%arg_cache%release_entry(new_cache_entry); new_cache_entry=>NULL()
                   else
                    errc=-10
                   endif
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
              class default
               errc=-6
              end select
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

          subroutine register_temp_tensor(copy_num,tens_entry,jerr)
           implicit none
           integer(INTD), intent(in):: copy_num
           class(tens_cache_entry_t), pointer, intent(out):: tens_entry
           integer(INTD), intent(out):: jerr
           character(TEREC_MAX_TENS_NAME_LEN+16):: tname
           class(tens_rcrsv_t), pointer:: temptens
           integer(INTD):: jl,jn
           logical:: stored

           tens_entry=>NULL()
           allocate(tens_rcrsv_t::temptens,STAT=jerr)
           if(jerr.eq.0) then
            call temptens%tens_rcrsv_ctor(header,jerr)
            if(jerr.eq.TEREC_SUCCESS) then
             call temptens%get_name(tname,jl,jerr)
             if(jerr.eq.TEREC_SUCCESS.and.jl.le.TEREC_MAX_TENS_NAME_LEN) then
              tname(jl+1:jl+1)='#'; jl=jl+1; call numchar(copy_num,jn,tname(jl+1:)); jl=jl+jn !name mangling
              call temptens%rename(tname(1:jl),jerr)
              if(jerr.eq.TEREC_SUCCESS) then
               stored=this%arg_cache%store(temptens,tens_entry_wrk_alloc,jerr,tens_entry_p=tens_entry) !tensor ownership is moved to the tensor cache entry
               if(.not.(jerr.eq.0.and.stored.and.associated(tens_entry))) then
                if(associated(tens_entry)) call this%arg_cache%release_entry(tens_entry); tens_entry=>NULL()
                jerr=-5
               endif
               if(DEBUG.gt.0.and.jerr.eq.0) then
                write(CONS_OUT,'("#MSG(TAVP-WRK:Resourcer)[",i6,"]: Output tensor renamed to")',ADVANCE='NO') impir
                write(CONS_OUT,*) tname(1:jl)
                flush(CONS_OUT)
               endif
              else
               jerr=-4
              endif
             else
              jerr=-3
             endif
            else
             jerr=-2
            endif
           else
            jerr=-1
           endif
           return
          end subroutine register_temp_tensor

        end subroutine TAVPWRKResourcerSubstituteOutput
!---------------------------------------------------------------------
        subroutine TAVPWRKResourcerRestoreOutput(this,tens_instr,ierr)
!Restores the original (persistent) output tensor in a tensor instruction
!and destroys the temporary tensor (assumes a single output tensor for now).
         implicit none
         class(tavp_wrk_resourcer_t), intent(inout):: this !inout: TAVP-WRK Resourcer
         class(tens_instr_t), intent(inout):: tens_instr   !in: active tensor instruction
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc,n,nou,nl,onl
         integer(INTD), pointer:: out_oprs(:)
         class(ds_oprnd_t), pointer:: oprnd
         class(tens_cache_entry_t), pointer:: tens_entry
         class(tens_entry_wrk_t), pointer:: cache_entry_tmp
         class(tens_rcrsv_t), pointer:: tensor
         class(tens_header_t), pointer:: header
         type(tens_rcrsv_t):: tensor_old
         character(TEREC_MAX_TENS_NAME_LEN+16):: tname
         logical:: evicted

         if(tens_instr%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           n=tens_instr%get_num_operands(errc)
           if(errc.eq.DSVP_SUCCESS.and.n.gt.0) then
            out_oprs=>tens_instr%get_output_operands(errc,nou)
            if(errc.eq.0.and.nou.eq.1.and.lbound(out_oprs,1).eq.0) then
             oprnd=>tens_instr%get_operand(out_oprs(0),errc) !output tensor operand
             if(errc.eq.DSVP_SUCCESS) then
              select type(oprnd)
              class is(tens_oprnd_t)
               cache_entry_tmp=>oprnd%get_cache_entry(errc)
               if(errc.eq.0.and.associated(cache_entry_tmp)) then
                tensor=>oprnd%get_tensor(errc)
                if(errc.eq.0.and.associated(tensor)) then
                 header=>tensor%get_header(errc)
                 if(errc.eq.TEREC_SUCCESS.and.associated(header)) then
                  call header%get_name(tname,nl,errc); onl=nl
                  if(errc.eq.TEREC_SUCCESS.and.nl.gt.2) then !mangled names have at least three characters
                   do while(nl.gt.0); if(tname(nl:nl).eq.'#') exit; nl=nl-1; enddo
                   if(nl.gt.1) then !the original tensor name has at least one character
                    nl=nl-1 !tname(1:nl) is now the unmangled name of the persistent tensor
                    call tensor_old%tens_rcrsv_ctor(header,errc)
                    if(errc.eq.TEREC_SUCCESS) then
                     call tensor_old%rename(tname(1:nl),errc) !original name of the persistent tensor
                     if(errc.eq.TEREC_SUCCESS) then
                      tens_entry=>NULL(); tens_entry=>this%arg_cache%lookup(tensor_old,errc) !lookup original tensor
                      if(errc.eq.0.and.associated(tens_entry)) then
                       select type(tens_entry)
                       class is(tens_entry_wrk_t)
                        call oprnd%reset_tmp_tensor(tens_entry,.FALSE.,errc); if(errc.ne.0) errc=-18
                       class default
                        errc=-17
                       end select
                       call this%arg_cache%release_entry(tens_entry); tens_entry=>NULL()
                       if(errc.eq.0) then
                        if(cache_entry_tmp%get_ref_count().eq.0.and.cache_entry_tmp%get_use_count().eq.0) then
                         evicted=this%arg_cache%evict(tensor,errc)
                         if(errc.eq.0) then; if(.not.evicted) errc=-16; else; errc=-15; endif
                        endif
                        if(DEBUG.gt.0) then
                         write(CONS_OUT,'("#MSG(TAVP-WRK:Resourcer)[",i6,"]: Output tensor renamed back from")',ADVANCE='NO') impir
                         write(CONS_OUT,*) tname(1:onl)
                         flush(CONS_OUT)
                        endif
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
              class default
               errc=-6
              end select
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
        end subroutine TAVPWRKResourcerRestoreOutput
![tavp_wrk_communicator_t]=====================================
        subroutine TAVPWRKCommunicatorConfigure(this,conf,ierr)
!Configures this DSVU.
         implicit none
         class(tavp_wrk_communicator_t), intent(inout):: this !out: configured DSVU (must not be configured on entrance)
         class(dsv_conf_t), intent(in):: conf                 !in: specific DSVU configuration
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTD):: errc

         errc=0
         select type(conf)
         type is(tavp_wrk_communicator_conf_t)
          if(conf%num_mpi_windows.gt.0) then
           this%num_mpi_windows=conf%num_mpi_windows
          else
           errc=-2
          endif
         class default
          errc=-1
         end select
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKCommunicatorConfigure
!-----------------------------------------------------
        subroutine TAVPWRKCommunicatorStart(this,ierr)
!Starts and lives this DSVU, calls .shutdown() at the end.
         implicit none
         class(tavp_wrk_communicator_t), intent(inout):: this !inout: TAVP-WRK communicator DSVU
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTD):: errc,ier,thid,n,num_fetch,num_upload,opcode,sts,errcode
         integer:: com_timer
         logical:: active,stopping,really_stopping,delivered
         class(dsvp_t), pointer:: dsvp
         class(tavp_wrk_t), pointer:: tavp
         class(tens_instr_t), pointer:: tens_instr
         class(*), pointer:: uptr

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Communicator started as DSVU # ",i2," (thread ",i2,'//&
          '"): Number of MPI windows = ",i6)') impir,this%get_id(),thid,this%num_mpi_windows
          flush(CONS_OUT)
         endif
!Initialize queues and ports:
         call this%init_queue(this%num_ports,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
!Initialize the prefetch queue:
         ier=this%fet_list%init(this%prefetch_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-1
!Initialize the upload queue:
         ier=this%upl_list%init(this%upload_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-1
!Initialize the dispatch queue:
         ier=this%dsp_list%init(this%dispatch_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-1
!Initialize the retire queue:
         ier=this%ret_list%init(this%retire_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-1
!Initialize the global addressing space and set up tensor argument cache:
         tavp=>NULL(); dsvp=>this%get_dsvp(); select type(dsvp); class is(tavp_wrk_t); tavp=>dsvp; end select
         if(associated(tavp)) then
          call tavp%addr_space%create(role_comm,this%num_mpi_windows,'TAVPWRKAddressSpace',ier)
          if(ier.eq.0) then
           this%addr_space=>tavp%addr_space
           this%arg_cache=>tavp%tens_cache
           write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Communicator unit ",i2," created a global addressing space successfully!")')&
           &impir,this%get_id()
          else
           if(errc.eq.0) errc=-1
          endif
!Sync with other TAVP units:
          call tavp%sync_units(errc,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
         else
          this%arg_cache=>NULL(); if(errc.eq.0) errc=-1
         endif
!Work loop:
         ier=timer_start(com_timer,MAX_COMMUNICATOR_PHASE_TIME); if(ier.ne.TIMERS_SUCCESS.and.errc.eq.0) errc=-1
         active=(errc.eq.0); stopping=(.not.active); really_stopping=.FALSE.; num_fetch=0; num_upload=0 !number of outstanding prefetches and uploads
         wloop: do while(active)
 !Get new instructions from Resourcer (port 0) into the prefetch queue:
          ier=this%fet_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          ier=this%unload_port(0,this%fet_list,num_moved=n); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-11; exit wloop; endif
          if(DEBUG.gt.0.and.n.gt.0) then
           write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Communicator unit ",i2," received ",i9," instructions from Resourcer")')&
           &impir,this%get_id(),n
           !ier=this%fet_list%reset(); ier=this%fet_list%scanp(action_f=tens_instr_print); ier=this%fet_list%reset_back() !print all instructions
           flush(CONS_OUT)
          endif
 !Initiate input prefetch:
          ier=this%iqueue%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          ier=this%dsp_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          ier=this%fet_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          do while(this%fet_list%get_status().eq.GFC_IT_ACTIVE.and.num_fetch.lt.MAX_COMMUNICATOR_PREFETCHES)
           if(stopping.and.errc.eq.0) then; errc=-1; exit wloop; endif !no other instruction can follow STOP
           uptr=>this%fet_list%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
           tens_instr=>NULL(); select type(uptr); class is(tens_instr_t); tens_instr=>uptr; end select
           if((.not.associated(tens_instr)).and.errc.eq.0) then; errc=-1; exit wloop; endif !trap
           sts=tens_instr%get_status(ier,errcode); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
           if(sts.ne.DS_INSTR_INPUT_WAIT.and.errc.eq.0) then; errc=-1; exit wloop; endif !trap
           opcode=tens_instr%get_code(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
           if(opcode.ge.TAVP_ISA_TENS_FIRST.and.opcode.le.TAVP_ISA_TENS_LAST) then !tensor instruction
            if(opcode.ne.TAVP_INSTR_TENS_CREATE.and.opcode.ne.TAVP_INSTR_TENS_DESTROY) then
             call this%prefetch_input(tens_instr,ier); if(ier.ne.0.and.errc.eq.0) then; errc=-1; exit wloop; endif
             ier=this%fet_list%move_elem(this%iqueue); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
             num_fetch=num_fetch+1
            else
             call tens_instr%set_status(DS_INSTR_READY_TO_EXEC,ier,errcode)
             if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
             ier=this%fet_list%move_elem(this%dsp_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
            endif
           else
            if(opcode.eq.TAVP_INSTR_CTRL_STOP) stopping=.TRUE.
            call tens_instr%set_status(DS_INSTR_READY_TO_EXEC,ier,errcode)
            if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
            ier=this%fet_list%move_elem(this%dsp_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
           endif
          enddo
 !Get completed instructions from Dispatcher (port 1) into the upload queue:
          ier=this%upl_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          ier=this%unload_port(1,this%upl_list,num_moved=n); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-11; exit wloop; endif
          if(DEBUG.gt.0.and.n.gt.0) then
           write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Communicator unit ",i2," received ",i9," instructions from Dispatcher")')&
           &impir,this%get_id(),n
           !ier=this%upl_list%reset(); ier=this%upl_list%scanp(action_f=tens_instr_print); ier=this%upl_list%reset_back() !print all instructions
           flush(CONS_OUT)
          endif
 !Initiate output upload:
          ier=this%iqueue%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          ier=this%ret_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          ier=this%upl_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          do while(this%upl_list%get_status().eq.GFC_IT_ACTIVE.and.num_upload.lt.MAX_COMMUNICATOR_UPLOADS)
           if(really_stopping.and.errc.eq.0) then; errc=-1; exit wloop; endif !trap
           uptr=>this%upl_list%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
           tens_instr=>NULL(); select type(uptr); class is(tens_instr_t); tens_instr=>uptr; end select
           if((.not.associated(tens_instr)).and.errc.eq.0) then; errc=-1; exit wloop; endif !trap
           sts=tens_instr%get_status(ier,errcode); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
           if(sts.ne.DS_INSTR_COMPLETED.and.errc.eq.0) then; errc=-1; exit wloop; endif !trap
           opcode=tens_instr%get_code(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
           if(opcode.ge.TAVP_ISA_TENS_FIRST.and.opcode.le.TAVP_ISA_TENS_LAST) then !tensor instruction
            if(tens_instr%get_num_out_operands().gt.0.and.opcode.ne.TAVP_INSTR_TENS_CREATE) then
             call this%upload_output(tens_instr,ier); if(ier.ne.0.and.errc.eq.0) then; errc=-1; exit wloop; endif
             ier=this%upl_list%move_elem(this%iqueue); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
             num_upload=num_upload+1
            else
             call tens_instr%set_status(DS_INSTR_RETIRED,ier,errcode)
             if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
             ier=this%upl_list%move_elem(this%ret_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
            endif
           else
            if(opcode.eq.TAVP_INSTR_CTRL_STOP) really_stopping=.TRUE.
            call tens_instr%set_status(DS_INSTR_RETIRED,ier,errcode)
            if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
            ier=this%upl_list%move_elem(this%ret_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
           endif
          enddo
 !Test outstanding communication completion (both fetch and upload):
          ier=this%dsp_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          ier=this%ret_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          ier=this%iqueue%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          do while(this%iqueue%get_status().eq.GFC_IT_ACTIVE)
           uptr=>this%iqueue%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
           tens_instr=>NULL(); select type(uptr); class is(tens_instr_t); tens_instr=>uptr; end select
           if((.not.associated(tens_instr)).and.errc.eq.0) then; errc=-1; exit wloop; endif !trap
           opcode=tens_instr%get_code(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
           sts=tens_instr%get_status(ier,errcode); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
           if(sts.eq.DS_INSTR_INPUT_WAIT) then !fetched
            delivered=this%sync_prefetch(tens_instr,ier,wait=COMMUNICATOR_BLOCKING)
            if(ier.ne.0.and.errc.eq.0) then; errc=-1; exit wloop; endif
            if(delivered) then
             num_fetch=num_fetch-1
             call tens_instr%set_status(DS_INSTR_READY_TO_EXEC,ier,errcode)
             if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
             ier=this%iqueue%move_elem(this%dsp_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
            else
             ier=this%iqueue%next()
            endif
           elseif(sts.eq.DS_INSTR_OUTPUT_WAIT) then !uploaded
            delivered=this%sync_upload(tens_instr,ier,wait=COMMUNICATOR_BLOCKING)
            if(ier.ne.0.and.errc.eq.0) then; errc=-1; exit wloop; endif
            if(delivered) then
             num_upload=num_upload-1
             call tens_instr%set_status(DS_INSTR_RETIRED,ier,errcode)
             if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
             ier=this%iqueue%move_elem(this%ret_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
            else
             ier=this%iqueue%next()
            endif
           else
            if(errc.eq.0) then; errc=-1; exit wloop; endif
           endif
          enddo
 !Pass ready instructions to Dispatcher (port 0) for execution:
          ier=this%dsp_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          if(this%dsp_list%get_status().eq.GFC_IT_ACTIVE) then
           ier=tavp%dispatcher%load_port(0,this%dsp_list); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          endif
 !Pass completed instructions to Resourcer (port 1) for resource release:
          ier=this%ret_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          if(this%ret_list%get_status().eq.GFC_IT_ACTIVE) then
           ier=tavp%resourcer%load_port(1,this%ret_list); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          endif
 !Exit condition:
          active=.not.(stopping.and.really_stopping.and.num_fetch.eq.0.and.num_upload.eq.0)
         enddo wloop
!Destroy the timer:
         ier=timer_destroy(com_timer); if(ier.ne.TIMERS_SUCCESS.and.errc.eq.0) errc=-3
!Record the error:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(TAVP-WRK)[",i6,"]: Communicator error ",i11," by thread ",i2)')&
         &impir,errc,thid
!Shutdown:
         call this%shutdown(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKCommunicatorStart
!--------------------------------------------------------
        subroutine TAVPWRKCommunicatorShutdown(this,ierr)
!Stops DSVU (returns back a clean configured state).
         implicit none
         class(tavp_wrk_communicator_t), intent(inout):: this !inout: TAVP-WRK communicator DSVU
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTD):: errc,ier,thid

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Communicator stopped as DSVU # ",i2," (thread ",i2,")")')&
          &impir,this%get_id(),thid
          flush(CONS_OUT)
         endif
!Release the tensor argument cache pointer:
         this%arg_cache=>NULL()
!Destroy the global addressing space:
         call this%addr_space%destroy(ier); if(ier.ne.0.and.errc.eq.0) errc=-14
         this%addr_space=>NULL()
         write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Communicator unit ",i2," destroyed global addressing space: Status ",i11)')&
         &impir,this%get_id(),ier
!Deactivate the retire list:
         ier=this%ret_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%ret_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-13
           ier=this%ret_list%delete_all()
          endif
          ier=this%ret_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-12
         else
          if(errc.eq.0) errc=-11
         endif
!Deactivate the dispatch list:
         ier=this%dsp_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%dsp_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-10
           ier=this%dsp_list%delete_all()
          endif
          ier=this%dsp_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-9
         else
          if(errc.eq.0) errc=-8
         endif
!Deactivate the upload list:
         ier=this%upl_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%upl_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-7
           ier=this%upl_list%delete_all()
          endif
          ier=this%upl_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-6
         else
          if(errc.eq.0) errc=-5
         endif
!Deactivate the prefetch list:
         ier=this%fet_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%fet_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-4
           ier=this%fet_list%delete_all()
          endif
          ier=this%fet_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-3
         else
          if(errc.eq.0) errc=-2
         endif
!Release queues:
         call this%release_queue(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
!Record an error, if any:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKCommunicatorShutdown
!------------------------------------------------------------------------
        subroutine TAVPWRKCommunicatorPrefetchInput(this,tens_instr,ierr)
!Starts prefetching input arguments for a given tensor instruction.
         implicit none
         class(tavp_wrk_communicator_t), intent(inout):: this !inout: TAVP-WRK communicator DSVU
         class(tens_instr_t), intent(inout):: tens_instr      !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTD):: errc,ier,n
         class(ds_oprnd_t), pointer:: oprnd

         n=tens_instr%get_num_operands(errc)
         if(errc.eq.DSVP_SUCCESS.and.n.gt.0) then
          do while(n.gt.0)
           n=n-1; if(tens_instr%operand_is_output(n)) cycle
           oprnd=>tens_instr%get_operand(n,ier)
           if(ier.eq.DSVP_SUCCESS) then
            call oprnd%prefetch(ier); if(ier.ne.0.and.errc.eq.0) errc=-3
           else
            errc=-2; exit
           endif
          enddo
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKCommunicatorPrefetchInput
!--------------------------------------------------------------------------------------
        function TAVPWRKCommunicatorSyncPrefetch(this,tens_instr,ierr,wait) result(res)
!Synchronizes on the input arguments prefetch, either WAIT or TEST.
         implicit none
         logical:: res                                        !out: TRUE if synchronized (all input operands)
         class(tavp_wrk_communicator_t), intent(inout):: this !inout: TAVP-WRK communicator DSVU
         class(tens_instr_t), intent(inout):: tens_instr      !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr          !out: error code
         logical, intent(in), optional:: wait                 !in: WAIT or TEST (defaults to WAIT)
         integer(INTD):: errc,ier,n
         class(ds_oprnd_t), pointer:: oprnd
         logical:: wt

         errc=0; res=.FALSE.
         wt=.TRUE.; if(present(wait)) wt=wait
         n=tens_instr%get_num_operands(errc)
         if(errc.eq.DSVP_SUCCESS.and.n.gt.0) then
          res=.TRUE.
          do while(n.gt.0)
           n=n-1; if(tens_instr%operand_is_output(n)) cycle
           oprnd=>tens_instr%get_operand(n,ier)
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
        end function TAVPWRKCommunicatorSyncPrefetch
!-----------------------------------------------------------------------
        subroutine TAVPWRKCommunicatorUploadOutput(this,tens_instr,ierr)
!Starts uploading the output argument for a given tensor instruction.
         implicit none
         class(tavp_wrk_communicator_t), intent(inout):: this !inout: TAVP-WRK communicator DSVU
         class(tens_instr_t), intent(inout):: tens_instr      !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTD):: errc,ier,i,n
         integer(INTD), pointer:: out_oprs(:)
         class(ds_oprnd_t), pointer:: oprnd
         class(tens_entry_wrk_t), pointer:: cache_entry

         out_oprs(0:)=>tens_instr%get_output_operands(errc,n)
         if(errc.eq.0) then
          if(tens_instr%output_substituted(errc)) then
           if(errc.eq.0) then
            do i=0,n-1
             oprnd=>tens_instr%get_operand(out_oprs(i),ier)
             if(ier.eq.DSVP_SUCCESS) then
              select type(oprnd)
              class is(tens_oprnd_t)
               !`Find the accumulator cache entry in the tensor cache
               call oprnd%upload_local(cache_entry,ier); if(ier.ne.0.and.errc.eq.0) errc=-8
               !`Still need to initiate an upload of the accumulator tensor into the persistent tensor
              class default
               errc=-7; exit
              end select
             else
              errc=-6; exit
             endif
            enddo
           else
            errc=-5
           endif
          else
           if(errc.eq.0) then
            do i=0,n-1
             oprnd=>tens_instr%get_operand(out_oprs(i),ier)
             if(ier.eq.DSVP_SUCCESS) then
              call oprnd%upload(ier); if(ier.ne.0.and.errc.eq.0) errc=-4
             else
              errc=-3; exit
             endif
            enddo
           else
            errc=-2
           endif
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKCommunicatorUploadOutput
!------------------------------------------------------------------------------------
        function TAVPWRKCommunicatorSyncUpload(this,tens_instr,ierr,wait) result(res)
!Synchronizes on the output upload, either TEST or WAIT.
         logical:: res                                        !out: TRUE if synchronized
         class(tavp_wrk_communicator_t), intent(inout):: this !inout: TAVP-WRK communicator DSVU
         class(tens_instr_t), intent(inout):: tens_instr      !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr          !out: error code
         logical, intent(in), optional:: wait                 !in: WAIT or TEST (defaults to WAIT)
         integer(INTD):: errc,ier,i,n
         integer(INTD), pointer:: out_oprs(:)
         class(ds_oprnd_t), pointer:: oprnd
         logical:: wt

         res=.FALSE.
         wt=.TRUE.; if(present(wait)) wt=wait
         if(tens_instr%output_substituted(errc)) then
          if(errc.eq.0) then
           !`Still need to sync the upload of the accumulator tensor into the persistent tensor
           res=.TRUE.
          else
           errc=-5
          endif
         else
          if(errc.eq.0) then
           out_oprs(0:)=>tens_instr%get_output_operands(errc,n)
           if(errc.eq.0) then
            do i=0,n-1
             oprnd=>tens_instr%get_operand(out_oprs(i),ier)
             if(ier.eq.DSVP_SUCCESS) then
              res=res.and.oprnd%sync(ier,wt)
              if(ier.ne.0.and.errc.eq.0) then; res=.FALSE.; errc=-4; endif
             else
              res=.FALSE.; if(errc.eq.0) errc=-3
             endif
            enddo
           else
            errc=-2
           endif
          else
           errc=-1
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end function TAVPWRKCommunicatorSyncUpload
![tavp_wrk_dispatcher_t]=====================================
        subroutine TAVPWRKDispatcherConfigure(this,conf,ierr)
!Configures this DSVU:
! (a) Imports enabled device ids:
! (b) Sets up tensor Microcode;
         implicit none
         class(tavp_wrk_dispatcher_t), intent(inout):: this !out: configured DSVU (must not be configured on entrance)
         class(dsv_conf_t), intent(in):: conf               !in: specific DSVU configuration
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc

         errc=0
         select type(conf)
         type is(tavp_wrk_dispatcher_conf_t)
          if(conf%host_buf_size.ge.0) then
           this%host_buf_size=conf%host_buf_size
           if(allocated(this%gpu_list)) deallocate(this%gpu_list)
           if(allocated(conf%gpu_list)) this%gpu_list=conf%gpu_list
           if(allocated(this%amd_list)) deallocate(this%amd_list)
           if(allocated(conf%amd_list)) this%amd_list=conf%amd_list
           if(allocated(this%mic_list)) deallocate(this%mic_list)
           if(allocated(conf%mic_list)) this%mic_list=conf%mic_list
           call set_microcode()
          else
           errc=-2
          endif
         class default
          errc=-1
         end select
         if(present(ierr)) ierr=errc
         return

         contains

          subroutine set_microcode() !implementation of each tensor operation
           this%microcode(TAVP_INSTR_TENS_CREATE)%instr_proc=>TAVPWRKExecTensorCreate
           this%microcode(TAVP_INSTR_TENS_DESTROY)%instr_proc=>TAVPWRKExecTensorDestroy
           this%microcode(TAVP_INSTR_TENS_INIT)%instr_proc=>TAVPWRKExecTensorInit
           this%microcode(TAVP_INSTR_TENS_CONTRACT)%instr_proc=>TAVPWRKExecTensorContract
           return
          end subroutine set_microcode

        end subroutine TAVPWRKDispatcherConfigure
!---------------------------------------------------
        subroutine TAVPWRKDispatcherStart(this,ierr)
!Starts and lives this DSVU, calls .shutdown() at the end.
         implicit none
         class(tavp_wrk_dispatcher_t), intent(inout):: this !inout: TAVP-WRK dispatcher DSVU
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc,ier,thid,n,sts,opcode,errcode,num_outstanding
         logical:: active,stopping,completed
         class(dsvp_t), pointer:: dsvp
         class(tavp_wrk_t), pointer:: tavp
         class(tens_instr_t), pointer:: tens_instr
         class(*), pointer:: uptr

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Dispatcher started as DSVU # ",i2,'//&
          &'" (thread ",i2,"): Host buffer size (B) = ",i15)') impir,this%get_id(),thid,this%host_buf_size
          flush(CONS_OUT)
         endif
!Initialize queues and ports:
         call this%init_queue(this%num_ports,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
!Initialize the issued instruction queue:
         ier=this%iss_list%init(this%issued_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-1
!Initialize the completed instruction queue:
         ier=this%cml_list%init(this%completed_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-1
!Initialize the numerical computing runtime (TAL-SH) and set up tensor argument cache:
         tavp=>NULL(); dsvp=>this%get_dsvp(); select type(dsvp); class is(tavp_wrk_t); tavp=>dsvp; end select
         if(associated(tavp)) then
          ier=talsh_init(this%host_buf_size,this%host_arg_max,this%gpu_list,this%mic_list,this%amd_list)
          if(ier.eq.TALSH_SUCCESS) then
           this%arg_cache=>tavp%tens_cache
          else
           if(errc.eq.0) errc=-1
          endif
!Sync with other TAVP units:
          call tavp%sync_units(errc,ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
         else
          this%arg_cache=>NULL(); if(errc.eq.0) errc=-1
         endif
!Work loop:
         active=(errc.eq.0); stopping=(.not.active); num_outstanding=0
         wloop: do while(active)
 !Get new instructions from Communicator (port 0) into the main queue:
          ier=this%iqueue%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          ier=this%flush_port(0,max_items=MAX_DISPATCHER_INTAKE,num_moved=n)
          if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          if(DEBUG.gt.0.and.n.gt.0) then
           write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Dispatcher unit ",i2," received ",i9," instructions from Communicator")')&
           &impir,this%get_id(),n
           !ier=this%iqueue%reset(); ier=this%iqueue%scanp(action_f=tens_instr_print); ier=this%iqueue%reset_back() !print all instructions
           flush(CONS_OUT)
          endif
 !Issue instructions:
          ier=this%iss_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          ier=this%cml_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          ier=this%iqueue%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          do while(this%iqueue%get_status().eq.GFC_IT_ACTIVE)
           if(stopping.and.errc.eq.0) then; errc=-1; exit wloop; endif !no other instruction can follow STOP
           uptr=>this%iqueue%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
           tens_instr=>NULL(); select type(uptr); class is(tens_instr_t); tens_instr=>uptr; end select
           if((.not.associated(tens_instr)).and.errc.eq.0) then; errc=-1; exit wloop; endif !trap
           sts=tens_instr%get_status(ier,errcode); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
           if(sts.ne.DS_INSTR_READY_TO_EXEC.and.errc.eq.0) then; errc=-1; exit wloop; endif !trap
           opcode=tens_instr%get_code(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
           if(opcode.ge.TAVP_ISA_TENS_FIRST.and.opcode.le.TAVP_ISA_TENS_LAST) then !tensor instruction
            call tens_instr%set_status(DS_INSTR_ISSUED,ier,errcode)
            if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
            select case(opcode)
            case(TAVP_INSTR_TENS_CREATE,TAVP_INSTR_TENS_DESTROY)
             call this%issue_instr(tens_instr,ier); if(ier.ne.0.and.errc.eq.0) then; errc=-1; exit wloop; endif
             num_outstanding=num_outstanding+1
            case default
             call this%issue_instr(tens_instr,ier); if(ier.ne.0.and.errc.eq.0) then; errc=-1; exit wloop; endif
             num_outstanding=num_outstanding+1
            end select
            ier=this%iqueue%move_elem(this%iss_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
           else !auxiliary or control instruction
            if(opcode.eq.TAVP_INSTR_CTRL_STOP) stopping=.TRUE.
            !`Process ctrl/aux instructions here
            call tens_instr%set_status(DS_INSTR_COMPLETED,ier,errcode)
            if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
            ier=this%iqueue%move_elem(this%cml_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
           endif
          enddo
 !Test/wait for completion of the issued instructions:
          ier=this%iss_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          ier=this%cml_list%reset_back(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          do while(this%iss_list%get_status().eq.GFC_IT_ACTIVE)
           uptr=>this%iss_list%get_value(ier); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
           tens_instr=>NULL(); select type(uptr); class is(tens_instr_t); tens_instr=>uptr; end select
           if((.not.associated(tens_instr)).and.errc.eq.0) then; errc=-1; exit wloop; endif !trap
           sts=tens_instr%get_status(ier,errcode); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
           completed=this%sync_instr(tens_instr,ier); if(ier.ne.0.and.errc.eq.0) then; errc=-1; exit wloop; endif
           if(completed) then
            num_outstanding=num_outstanding-1
            call tens_instr%set_status(DS_INSTR_COMPLETED,ier,errcode)
            if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
            ier=this%iss_list%move_elem(this%cml_list); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
           else
            ier=this%iss_list%next()
           endif
          enddo
 !Pass completed instructions back to Communicator (port 1):
          ier=this%cml_list%reset(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) then; errc=-1; exit wloop; endif
          if(this%cml_list%get_status().eq.GFC_IT_ACTIVE) then
           !`Pass
          endif
 !Exit condition:
          active=.not.(stopping.and.num_outstanding.eq.0)
         enddo wloop
!Record the error:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(errc.ne.0.and.VERBOSE) write(CONS_OUT,'("#ERROR(TAVP-WRK)[",i6,"]: Dispatcher error ",i11," by thread ",i2)')&
         &impir,errc,thid
!Shutdown:
         call this%shutdown(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKDispatcherStart
!------------------------------------------------------
        subroutine TAVPWRKDispatcherShutdown(this,ierr)
!Stops DSVU (returns back a clean configured state).
         implicit none
         class(tavp_wrk_dispatcher_t), intent(inout):: this !inout: TAVP-WRK dispatcher DSVU
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc,ier,thid

         errc=0; thid=omp_get_thread_num()
         if(DEBUG.gt.0) then
          write(CONS_OUT,'("#MSG(TAVP-WRK)[",i6,"]: Dispatcher stopped as DSVU # ",i2," (thread ",i2,")")')&
          &impir,this%get_id(),thid
          flush(CONS_OUT)
         endif
!Release the tensor argument cache pointer:
         this%arg_cache=>NULL()
!Shutdown the numerical computing runtime (TAL-SH):
         ier=talsh_shutdown(); if(ier.ne.TALSH_SUCCESS.and.errc.eq.0) errc=-8
!Deactivate the completed list:
         ier=this%cml_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%cml_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-7
           ier=this%cml_list%delete_all()
          endif
          ier=this%cml_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-6
         else
          if(errc.eq.0) errc=-5
         endif
!Deactivate the issued list:
         ier=this%iss_list%reset()
         if(ier.eq.GFC_SUCCESS) then
          ier=this%iss_list%get_status()
          if(ier.ne.GFC_IT_EMPTY) then
           if(errc.eq.0) errc=-4
           ier=this%iss_list%delete_all()
          endif
          ier=this%iss_list%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-3
         else
          if(errc.eq.0) errc=-2
         endif
!Release queues:
         call this%release_queue(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-1
!Record an error, if any:
         ier=this%get_error(); if(ier.eq.DSVP_SUCCESS) call this%set_error(errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKDispatcherShutdown
!--------------------------------------------------------------------------
        subroutine TAVPWRKDispatcherIssueInstr(this,tens_instr,ierr,dev_id)
!Issues the given tensor instruction to a specific compute device (or default).
         implicit none
         class(tavp_wrk_dispatcher_t), intent(inout):: this !inout: TAVP-WRK dispatcher DSVU
         class(tens_instr_t), intent(inout):: tens_instr    !inout: defined tensor instruction ready to be issued
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD), intent(in), optional:: dev_id       !in: flat device id to issue the tensor instruction to
         integer(INTD):: errc,opcode
         procedure(tavp_wrk_dispatch_proc_i), pointer:: iproc

         if(tens_instr%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           opcode=tens_instr%get_code(errc)
           if(errc.eq.DSVP_SUCCESS) then
            iproc=>this%microcode(opcode)%instr_proc
            if(associated(iproc)) then
             if(present(dev_id)) then
              call iproc(this,tens_instr,errc,dev_id)
             else
              call iproc(this,tens_instr,errc)
             endif
             if(errc.ne.0) errc=-5
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
        end subroutine TAVPWRKDispatcherIssueInstr
!---------------------------------------------------------------------------------
        function TAVPWRKDispatcherSyncInstr(this,tens_instr,ierr,wait) result(res)
!Synchronization on tensor instruction execution, either TEST or WAIT.
         implicit none
         logical:: res                                      !out: TRUE if tensor instruction has completed
         class(tavp_wrk_dispatcher_t), intent(inout):: this !inout: TAVP-WRK dispatcher DSVU
         class(tens_instr_t), intent(inout):: tens_instr    !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr        !out: error code
         logical, intent(in), optional:: wait               !in: WAIT or TEST (defaults to WAIT)
         integer(INTD):: errc,sts,ans
         logical:: wt

         errc=0; res=.FALSE.
         if(talsh_task_status(tens_instr%talsh_task).eq.TALSH_TASK_EMPTY) then
          res=.TRUE.
         else
          wt=.TRUE.; if(present(wait)) wt=wait
          if(wt) then !WAIT
           errc=talsh_task_wait(tens_instr%talsh_task,sts)
           res=(errc.eq.TALSH_SUCCESS.and.(sts.eq.TALSH_TASK_COMPLETED.or.sts.eq.TALSH_TASK_ERROR))
           if(errc.ne.TALSH_SUCCESS) errc=-3
          else !TEST
           ans=talsh_task_complete(tens_instr%talsh_task,sts,errc)
           res=(ans.eq.YEP.and.errc.eq.TALSH_SUCCESS)
           if(errc.ne.TALSH_SUCCESS) errc=-2
          endif
          if(sts.eq.TALSH_TASK_ERROR.and.errc.eq.0) errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TAVPWRKDispatcherSyncInstr
!----------------------------------------------------------------------
        subroutine TAVPWRKExecTensorCreate(this,tens_instr,ierr,dev_id)
!Executes tensor creation. The tensor layout is assumed already defined.
!The tensor body location is set here via the newly created DDSS descriptor.
!The tensor body location comes from the local resource associated with the tensor.
         implicit none
         class(tavp_wrk_dispatcher_t), intent(inout):: this !inout: TAVP-WRK dispatcher DSVU
         class(tens_instr_t), intent(inout):: tens_instr    !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD), intent(in), optional:: dev_id       !in: flat device id
         integer(INTD):: errc,dtk
         integer(INTL):: bytes,vol
         class(dsvp_t), pointer:: dsvp
         class(tavp_wrk_t), pointer:: tavp
         class(ds_oprnd_t), pointer:: oprnd
         class(tens_rcrsv_t), pointer:: tensor
         class(tens_resrc_t), pointer:: resource
         class(tens_body_t), pointer:: tens_body
         class(tens_layout_t), pointer:: tens_layout
         type(DataDescr_t):: descr
         type(C_PTR):: mem_p

!$OMP FLUSH
         oprnd=>tens_instr%get_operand(0,errc)
         if(errc.eq.DSVP_SUCCESS) then
          select type(oprnd)
          class is(tens_oprnd_t)
           resource=>oprnd%get_resource(errc)
           if(errc.eq.0) then
            if(.not.resource%is_empty()) then !resource is supposed to be preallocated by Resourcer
             tensor=>oprnd%get_tensor(errc)
             if(errc.eq.0) then
              tens_body=>tensor%get_body(errc)
              if(errc.eq.TEREC_SUCCESS) then
               tens_layout=>tens_body%get_layout(errc)
               if(errc.eq.TEREC_SUCCESS) then
                vol=tens_layout%get_volume()
                dtk=tens_layout%get_data_type(errc)
                if(errc.eq.TEREC_SUCCESS) then
                 mem_p=resource%get_mem_ptr(errc,bytes) !bytes = vol * sizeof(data_kind)
                 if(errc.eq.0) then
                  dsvp=>this%get_dsvp(errc)
                  if(errc.eq.DSVP_SUCCESS.and.associated(dsvp)) then
                   tavp=>NULL(); select type(dsvp); class is(tavp_wrk_t); tavp=>dsvp; end select
                   if(associated(tavp)) then
                    call tavp%addr_space%attach(mem_p,dtk,vol,descr,errc)
                    if(errc.eq.0) then
!$OMP CRITICAL (TAVP_WRK_CACHE)
                     call tensor%set_location(descr,errc) !tensor has been located
!$OMP END CRITICAL (TAVP_WRK_CACHE)
                     if(errc.eq.TEREC_SUCCESS) then
                      if(DEBUG.gt.0) then
                       write(CONS_OUT,'("#DEBUG(TAVP-WRK:TENSOR_CREATE)[",i6,"]: Tensor created: Size (bytes) = ",i13,":")')&
                       &impir,bytes
                       call tensor%print_it(dev_id=CONS_OUT)
                       flush(CONS_OUT)
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
           else
            errc=-3
           endif
          class default
           errc=-2
          end select
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKExecTensorCreate
!-----------------------------------------------------------------------
        subroutine TAVPWRKExecTensorDestroy(this,tens_instr,ierr,dev_id)
!Executes tensor destruction.
         implicit none
         class(tavp_wrk_dispatcher_t), intent(inout):: this !inout: TAVP-WRK dispatcher DSVU
         class(tens_instr_t), intent(inout):: tens_instr    !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD), intent(in), optional:: dev_id       !in: flat device id
         integer(INTD):: errc

         errc=0
         !`Implement
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKExecTensorDestroy
!--------------------------------------------------------------------
        subroutine TAVPWRKExecTensorInit(this,tens_instr,ierr,dev_id)
!Executes tensor initialization.
         implicit none
         class(tavp_wrk_dispatcher_t), intent(inout):: this !inout: TAVP-WRK dispatcher DSVU
         class(tens_instr_t), intent(inout):: tens_instr    !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD), intent(in), optional:: dev_id       !in: flat device id
         integer(INTD):: errc

         errc=0
         !`Implement
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKExecTensorInit
!------------------------------------------------------------------------
        subroutine TAVPWRKExecTensorContract(this,tens_instr,ierr,dev_id)
!Executes tensor contraction.
         implicit none
         class(tavp_wrk_dispatcher_t), intent(inout):: this !inout: TAVP-WRK dispatcher DSVU
         class(tens_instr_t), intent(inout):: tens_instr    !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD), intent(in), optional:: dev_id       !in: flat device id
         integer(INTD):: errc

         errc=0
         !`Implement
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKExecTensorContract
![tavp_wrk_t]======================================
        subroutine TAVPWRKConfigure(this,conf,ierr)
!Configures TAVP-WRK DSVP:
! * Configures static DSVU;
! * Allocates and configures dynamic DSVU;
! * Sets up global DSVU table in DSVP;
! * Sets up DSVP description and global id;
         implicit none
         class(tavp_wrk_t), intent(inout), target:: this !out: configured DSVP (must not be configured on entrance)
         class(dsv_conf_t), intent(in):: conf            !in: specific DSVP configuration
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,num_units
         class(ds_unit_t), pointer:: decode_acceptor
         type(tavp_wrk_decoder_conf_t):: decoder_conf
         type(tavp_wrk_retirer_conf_t):: retirer_conf
         type(tavp_wrk_resourcer_conf_t):: resourcer_conf
         type(tavp_wrk_communicator_conf_t):: communicator_conf
         type(tavp_wrk_dispatcher_conf_t):: dispatcher_conf

         if(.not.this%is_configured(errc)) then
          if(errc.eq.0) then
           select type(conf)
           type is(tavp_wrk_conf_t)
            if(conf%tavp_id.ge.0.and.allocated(conf%description)) then
             num_units=0 !increment by one after each unit configuration
 !Configure static DSVU:
  !Decoder:
             decode_acceptor=>this%resourcer
             decoder_conf=tavp_wrk_decoder_conf_t(source_comm=conf%source_comm,source_rank=conf%source_rank,&
                         &acceptor=decode_acceptor,acceptor_port_id=0)
             call this%decoder%configure(decoder_conf,errc)
             if(errc.eq.0) then
              num_units=num_units+1
  !Retirer:
              retirer_conf=tavp_wrk_retirer_conf_t(conf%retire_comm,conf%retire_rank)
              call this%retirer%configure(retirer_conf,errc)
              if(errc.eq.0) then
               num_units=num_units+1
  !Resourcer:
               resourcer_conf=tavp_wrk_resourcer_conf_t(conf%host_ram_size,conf%nvram_size)
               call this%resourcer%configure(resourcer_conf,errc)
               if(errc.eq.0) then
                num_units=num_units+1
  !Communicator:
                communicator_conf=tavp_wrk_communicator_conf_t(conf%num_mpi_windows)
                call this%communicator%configure(communicator_conf,errc)
                if(errc.eq.0) then
                 num_units=num_units+1
  !Dispatcher:
                 dispatcher_conf=tavp_wrk_dispatcher_conf_t(conf%host_buf_size,conf%gpu_list,conf%amd_list,conf%mic_list)
                 call this%dispatcher%configure(dispatcher_conf,errc)
                 if(errc.eq.0) then
                  num_units=num_units+1
 !Set up global DSVU table (references to all DSVU):
                  call this%alloc_units(num_units,errc)
                  if(errc.eq.DSVP_SUCCESS) then
                   call this%set_unit(this%decoder,errc)
                   if(errc.eq.DSVP_SUCCESS) then
                    call this%set_unit(this%retirer,errc)
                    if(errc.eq.DSVP_SUCCESS) then
                     call this%set_unit(this%resourcer,errc)
                     if(errc.eq.DSVP_SUCCESS) then
                      call this%set_unit(this%communicator,errc)
                      if(errc.eq.DSVP_SUCCESS) then
                       call this%set_unit(this%dispatcher,errc)
                       if(errc.eq.DSVP_SUCCESS) then
 !Set the DSVP id and description:
                        call this%set_description(int(conf%tavp_id,INTL),conf%description,errc)
                        if(errc.ne.DSVP_SUCCESS) errc=-16
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
         !write(CONS_OUT,*) '#DEBUG(TAVPWRKConfigure): Exit status ',errc !debug
         if(errc.ne.0) call this%destroy()
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWRKConfigure

       end module tavp_worker
