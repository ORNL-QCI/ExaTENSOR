!ExaTENSOR: Parallel Virtual Processing for Scale-Adaptive Hierarchical Tensor Algebra:
!Derives from abstract domain-specific virtual processing (DSVP).
!This module provides basic infrastructure for ExaTENSOR, tensor algebra virtual processor (TAVP).
!Different specializations (roles) of tensor algebra virtual processors derive from this module.
!Note that, in general, different specializations (roles) of the domain-specific virtual processor
!may have differing instruction sets (non-overlapping, overlapping, or identical). The instruction
!codes provided in this module are common for all specializations of the tensor algebra virtual processors.
!However, different specializations always have different microcodes, even for the same instruction codes.

!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2018/08/17

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

!NOTES:
! # Tensor Cache:
!   (a) Tensor cache operations, namely, lookup(), store(), and evict() are serialized
!       by a cache-wide lock.
!   (b) Returning a pointer to a tensor cache entry in lookup() and store() results in
!       an increment of that entry's USE_COUNT whereas releasing the obtained pointer
!       via the member procedure .release_entry() decrements that entry's USE_COUNT.
!   (c) Associating a tensor operand with a tensor cache entry increments that entry's
!       REF_COUNT whereas destroying that tensor operand decrements the REF_COUNT.
!   (d) Tensor cache entries created by the TENS_CREATE instruction are marked persistent.
!       Persistent tensor cache entries as well as tensor cache entries with a non-zero
!       USE_COUNT or REF_COUNT cannot be evicted/destroyed.
!   (e) Tensor cache entries created by the TENS_CREATE instruction are marked up-to-date
!       upon creation. Temporary tensor cache entries are marked up-to-date only upon the
!       completion of the first prefetch of their data.
!   (f) Accessing/updating the content of a tensor cache entry via the tensor cache entry
!       itself must be protected by the lock provided by the tensor cache entry:
!       1. Get cache entry pointer;
!       2. Lock cache entry (thread-exclusive access guaranteed);
!       3. Access/update the content of the cache entry by the thread;
!       4. Unlock cache entry;
!       5. Release the cache entry pointer.
!   (g) Accessing/updating the content of a tensor cache entry indirectly via
!       accessing/updating the pointer components of the tensor operand associated
!       with the tensor cache entry must be protected by the lock provided by the tensor
!       operand which delegates locking/unlocking to the associated tensor cache entry:
!       1. Lock tensor operand (locks the associated tensor cache entry for thread-exclusive access);
!       2. Access/update the content of the tensor operand;
!       3. Unlock tensor operand (unlocks the associated tensor cache entry).
!   (h) An issue of a tensor instruction increments READ/WRITE counters of the tensor
!       cache entries associated with the INPUT/OUTPUT tensor operands, respectively.
!       The completion of the tensor instruction decrements those counters.
!       A deferrence of a tensor instruction increments deferred READ/WRITE counters
!       of the tensor cache entries associated with the INPUT/OUTPUT tensor operands,
!       respectively. The subsequent actual issue of the tensor instruction decrements
!       those counters.

       module virta !VIRtual Tensor Algebra
        use tensor_algebra     !basic constants
        use talsh              !on-node heterogeneous numeric tensor algebra
        use pack_prim          !object packing/unpacking primitives
        use distributed        !distributed one-sided communication layer
        use service_mpi        !basic MPI service
        use hardware           !hardware abstraction and aggregation
        use subspaces          !hierarchical vector space representation
        use tensor_recursive   !recursive (hierarchical) tensors
        use dsvp_base          !abstract domain-specific virtual processor (DSVP)
        use gfc_base           !GFC base
        use gfc_dictionary     !GFC dictionary
        use gfc_vector         !GFC vector
        implicit none
        public
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !default output for this module
        integer(INTD), private:: DEBUG=1    !debugging mode
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
 !Runtime errors (ExaTENSOR aliases of basic DSVP errors):
        integer(INTD), parameter, public:: EXA_SUCCESS=DSVP_SUCCESS                         !success
        integer(INTD), parameter, public:: EXA_ERROR=DSVP_ERROR                             !generic error
        integer(INTD), parameter, public:: EXA_ERR_INVALID_ARGS=DSVP_ERR_INVALID_ARGS       !invalid arguments passed
        integer(INTD), parameter, public:: EXA_ERR_INVALID_REQ=DSVP_ERR_INVALID_REQ         !invalid request
        integer(INTD), parameter, public:: EXA_ERR_MEM_ALLOC_FAIL=DSVP_ERR_MEM_ALLOC_FAIL   !memory allocation failed
        integer(INTD), parameter, public:: EXA_ERR_MEM_FREE_FAIL=DSVP_ERR_MEM_FREE_FAIL     !memory deallocation failed
        integer(INTD), parameter, public:: EXA_ERR_BROKEN_OBJ=DSVP_ERR_BROKEN_OBJ           !broken object
        integer(INTD), parameter, public:: EXA_ERR_UNABLE_COMPLETE=DSVP_ERR_UNABLE_COMPLETE !unable to complete
        integer(INTD), parameter, public:: EXA_ERR_RSC_EXCEEDED=DSVP_ERR_RSC_EXCEEDED       !resource exceeded
 !Tensor-algebra virtual processor (TAVP) kinds (roles):
        integer(INTD), parameter, public:: EXA_NO_ROLE=DSVP_NO_KIND !undefined role
        integer(INTD), parameter, public:: EXA_DRIVER=0             !domain-specific driver process (not a TAVP)
        integer(INTD), parameter, public:: EXA_MANAGER=1            !manager (logical) process (TAVP)
        integer(INTD), parameter, public:: EXA_WORKER=2             !worker (numerical) process (TAVP)
        integer(INTD), parameter, public:: EXA_HELPER=3             !helper (auxiliary) process (TAVP)
 !Tensor data kinds:
        integer(INTD), parameter, public:: EXA_DATA_KIND_NN=NO_TYPE !none
        integer(INTD), parameter, public:: EXA_DATA_KIND_R4=R4      !single precision real
        integer(INTD), parameter, public:: EXA_DATA_KIND_R8=R8      !double precision real
        integer(INTD), parameter, public:: EXA_DATA_KIND_C4=C4      !single precision complex
        integer(INTD), parameter, public:: EXA_DATA_KIND_C8=C8      !double precision complex
 !External methods:
        integer(INTD), parameter, public:: EXA_MAX_METHOD_NAME_LEN=64 !max length of an external method name
 !TAVP hierarchy configuration:
        integer(INTD), public:: EXA_MAX_WORK_GROUP_SIZE=2 !maximal size of a work group (max number of workers per manager)
        integer(INTD), public:: EXA_MANAGER_BRANCH_FACT=2 !branching factor for the managing hierarchy
 !TAVP identification:
        integer(INTD), parameter, public:: TAVP_ANY_ID=-1         !any TAVP
 !TAVP MPI message tags:
        integer(INT_MPI), parameter, public:: TAVP_DEFAULT_TAG=0  !default MPI message tag
        integer(INT_MPI), parameter, public:: TAVP_DISPATCH_TAG=1 !MPI message containing dispatched instructions
        integer(INT_MPI), parameter, public:: TAVP_COLLECT_TAG=2  !MPI message containing instructions being collected after execution
        integer(INT_MPI), parameter, public:: TAVP_LOCATE_TAG=3   !MPI message containing instructions undergoing meta-data location
        integer(INT_MPI), parameter, public:: TAVP_REPLICA_TAG=4  !MPI message containing instructions for data replication
 !TAVP instruction error codes [-1:-100]:
        integer(INTD), parameter, public:: TAVP_ERR_GEN_FAILURE=-1     !unspecified generic failure
        integer(INTD), parameter, public:: TAVP_ERR_BTC_BAD=-2         !bad instruction bytecode
        integer(INTD), parameter, public:: TAVP_ERR_CHE_FAILURE=-3     !argument cache failure
        integer(INTD), parameter, public:: TAVP_ERR_ARG_UNDEFINED=-4   !instruction argument is undefined in the argument cache
        integer(INTD), parameter, public:: TAVP_ERR_ARG_DEFINED=-5     !instrcution argument is already defined in the argument cache
        integer(INTD), parameter, public:: TAVP_ERR_RSC_UNAVAILABLE=-6 !instruction is unable to obtain needed resources
        integer(INTD), parameter, public:: TAVP_ERR_RSC_FAILURE=-7     !instruction is unable to release resources cleanly
        integer(INTD), parameter, public:: TAVP_ERR_COM_FAILURE=-8     !instruction communication failed
        integer(INTD), parameter, public:: TAVP_ERR_EXC_FAILURE=-9     !instruction computation (execution) failed
 !TAVP special instruction markup (via error code):
        integer(INTD), parameter, public:: TAVP_ERR_TAG_ONE=-101       !special tag 1 (for instruction markup, not an error)
        integer(INTD), parameter, public:: TAVP_ERR_TAG_TWO=-102       !special tag 2 (for instruction markup, not an error)
 !Tensor algebra virtual processor (TAVP) ISA:
  !General:
        integer(INTD), parameter, public:: TAVP_ISA_SIZE=256       !max number of TAVP instruction codes [0:TAVP_ISA_SIZE-1]
        integer(INTD), parameter, public:: TAVP_ISA_CTRL_FIRST=0   !first TAVP_INSTR_CTRL_XXX code
        integer(INTD), parameter, public:: TAVP_ISA_CTRL_LAST=15   !last TAVP_INSTR_CTRL_XXX code
        integer(INTD), parameter, public:: TAVP_ISA_SPACE_FIRST=16 !first TAVP_INSTR_SPACE_XXX code
        integer(INTD), parameter, public:: TAVP_ISA_SPACE_LAST=63  !last TAVP_INSTR_SPACE_XXX code
        integer(INTD), parameter, public:: TAVP_ISA_TENS_FIRST=64  !first TAVP_INSTR_TENS_XXX code
        integer(INTD), parameter, public:: TAVP_ISA_TENS_LAST=255  !last TAVP_INSTR_TENS_XXX code
  !TAVP instruction code (opcode), must be non-negative (consult TAProL spec), limited by TAVP_ISA_SIZE:
   !NOOP:
        integer(INTD), parameter, public:: TAVP_INSTR_NOOP=DS_INSTR_NOOP !no operation (empty instruction): Negative opcode
   !General control [0-15]:
        integer(INTD), parameter, public:: TAVP_INSTR_CTRL_RESUME=0     !resume TAVP execution (used for special purposes or has no effect)
        integer(INTD), parameter, public:: TAVP_INSTR_CTRL_STOP=1       !stop TAVP (finishes current instructions and shuts down TAVP)
        integer(INTD), parameter, public:: TAVP_INSTR_CTRL_DUMP_CACHE=2 !dumps the cache of each TAVP into the log file
   !Auxiliary space definitions [16-63]:
        integer(INTD), parameter, public:: TAVP_INSTR_SPACE_CREATE=16   !create a (hierarchical) vector space
        integer(INTD), parameter, public:: TAVP_INSTR_SPACE_DESTROY=17  !destroy a (hierarchical) vector space
   !Tensor operations (decomposable) [64-255]:
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_CREATE=64    !tensor creation
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_DESTROY=65   !tensor destruction
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_LOAD=66      !tensor loading from persistent storage
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_SAVE=67      !tensor saving in persistent storage
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_INIT=68      !tensor initialization (assignment of a value)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_NORM1=69     !tensor 1-norm (sum of absolute values)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_NORM2=70     !tensor 2-norm (square root of the sum of squares of absolute values)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_MIN=71       !tensor min element
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_MAX=72       !tensor max element
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_FOLD=73      !tensor dimension folding (flattening)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_UNFOLD=74    !tensor dimension unfolding (tensorization)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_SLICE=75     !tensor slicing (extracting a slice of a tensor with an optional index permutation)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_INSERT=76    !tensor insertion (inserting a tensor slice in a tensor with an optional index permutation)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_COPY=77      !tensor copy (copying a tensor into another tensor with an optional index permutation)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_PERMUTE=78   !tensor dimension permutation (in-place)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_SCALE=79     !tensor scaling (multiplication by a real/complex number)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_ADD=80       !tensor addition (with an optional index permutation)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_TRACE=81     !tensor trace (produces a lower-rank tensor by tracing over some/all tensor indices)
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_CONTRACT=82  !tensor contraction (also includes tensor product, tensor addition, and tensor scaling as special cases)
   !Special tensor operations (TAVP-WRK only):
        integer(INTD), parameter, public:: TAVP_INSTR_TENS_ACCUMULATE=255 !tensor accumulation (accumulates one tensor into another one with no index permutation)
!TYPES:
 !Tensor status:
        type, public:: tens_status_t
         logical, public:: created=.FALSE. !TRUE if the tensor has been created (allocated physical memory), FALSE otherwise
         logical, public:: defined=.FALSE. !TRUE if the tensor value has been defined, FALSE otherwise
         logical, public:: replica=.FALSE. !TRUE if the tensor is a replica (another defined instance of this same tensor exists), FALSE otherwise
         logical, public:: updated=.FALSE. !TRUE if the tensor value is currently being updated, FALSE otherwise
         integer(INTD), public:: is_used=0 !number of read-only references which are currently using the tensor
        end type tens_status_t
        type(tens_status_t), parameter:: TENS_STATUS_NONE=tens_status_t() !tensor status null
        public TENS_STATUS_NONE
 !Tensor instruction control fields:
  !Tensor transformation (initialization, scaling, etc.) control field:
        type, extends(ds_instr_ctrl_t), public:: ctrl_tens_trans_t
         class(tens_method_uni_t), pointer, private:: definer=>NULL() !non-owning pointer to the defining unary tensor method (initialization/transformation)
         character(EXA_MAX_METHOD_NAME_LEN), private:: method_name    !name of the defining method
         complex(8), private:: alpha=(0d0,0d0)                        !either an initialization scalar or an alpha prefactor
         contains
          procedure, private:: CtrlTensTransCtor                    !ctor
          generic, public:: ctrl_tens_trans_ctor=>CtrlTensTransCtor
          procedure, public:: map_method=>CtrlTensTransMapMethod    !maps the defining method name to the actual definer object
          procedure, public:: pack=>CtrlTensTransPack               !packs the instruction control field into a plain byte packet
          procedure, public:: unpack=>CtrlTensTransUnpack           !unpacks the instruction control field from a plain byte packet
          procedure, public:: print_it=>CtrlTensTransPrintIt        !prints
          final:: ctrl_tens_trans_dtor                              !dtor
        end type ctrl_tens_trans_t
  !Tensor addition/copy/slice/insert control field:
        type, extends(ds_instr_ctrl_t), public:: ctrl_tens_add_t
         type(permutation_t), private:: permutation                 !tensor dimension permutation
         complex(8), private:: alpha=(1d0,0d0)                      !alpha prefactor
         contains
          procedure, private:: CtrlTensAddCtor                           !ctor
          generic, public:: ctrl_tens_add_ctor=>CtrlTensAddCtor
          procedure, public:: get_permutation=>CtrlTensAddGetPermutation !returns a non-owning pointer to the encapsulated permutation
          procedure, public:: pack=>CtrlTensAddPack                      !packs the instruction control field into a plain byte packet
          procedure, public:: unpack=>CtrlTensAddUnpack                  !unpacks the instruction control field from a plain byte packet
          procedure, public:: print_it=>CtrlTensAddPrintIt               !prints
          final:: ctrl_tens_add_dtor                                     !dtor
        end type ctrl_tens_add_t
  !Tensor contraction control field:
        type, extends(ds_instr_ctrl_t), public:: ctrl_tens_contr_t
         type(contr_ptrn_ext_t), private:: contr_ptrn           !extended tensor contraction pattern
         complex(8), private:: alpha=(1d0,0d0)                  !alpha prefactor
         contains
          procedure, private:: CtrlTensContrCtor                        !ctor
          generic, public:: ctrl_tens_contr_ctor=>CtrlTensContrCtor
          procedure, public:: get_contr_ptrn=>CtrlTensContrGetContrPtrn !returns a non-owning pointer to the extended tensor contraction pattern and, optionally, the prefactor and conjugation flags
          procedure, public:: pack=>CtrlTensContrPack                   !packs the instruction control field into a plain byte packet
          procedure, public:: unpack=>CtrlTensContrUnpack               !unpacks the instruction control field from a plain byte packet
          procedure, public:: print_it=>CtrlTensContrPrintIt            !prints
          final:: ctrl_tens_contr_dtor                                  !dtor
        end type ctrl_tens_contr_t
 !Tensor argument cache entry:
        type, abstract, public:: tens_cache_entry_t
         class(tens_rcrsv_t), pointer, private:: tensor=>NULL() !either owning or non-owning pointer to a tensor (ctor/dtor of an extended type will decide)
         integer(INTD), private:: ref_count=0                   !reference count: Number of existing tensor operands associated with this tensor cache entry
         integer(INTD), private:: use_count=0                   !use count: Number of active pointers to this tensor cache entry explicitly returned by the tensor cache
         integer(INTD), private:: read_write_count=0            !read/write count: Number of issued tensor instructions which refer to this tensor cache entry as input (positive) or output (negative)
         integer(INTD), private:: read_write_def_count=0        !deferred read/write count: Number of deferred tensor instructions which refer to this tensor cache entry as input (positive) or output (negative)
         integer(INTD), private:: temp_count=0                  !temporary count: Number of temporary tensors stemmed from this tensor cache entry (used for output rename)
         logical, private:: up_to_date=.FALSE.                  !up-to-date flag (TRUE means the tensor is defined, that is, neither undefined nor being updated)
         logical, private:: persistent=.FALSE.                  !persistency flag (persistent cache entries can only be evicted via an explicit TENS_DESTROY)
#ifndef NO_OMP
         integer(omp_nest_lock_kind), private:: entry_lock=-1   !tensor cache entry lock
         logical, private:: lock_initialized=.FALSE.            !lock initialization status
#endif
         contains
          procedure(tens_cache_entry_print_i), deferred:: print_it        !prints
          procedure, public:: is_set=>TensCacheEntryIsSet                 !returns TRUE if the tensor cache entry is set (constructed)
          procedure, public:: set_tensor=>TensCacheEntrySetTensor         !sets the pointer to a tensor, either owning or non-owning (basic ctor)
          procedure, public:: get_tensor=>TensCacheEntryGetTensor         !returns a non-owning pointer to the tensor
          procedure, public:: incr_ref_count=>TensCacheEntryIncrRefCount  !increments the reference count
          procedure, public:: decr_ref_count=>TensCacheEntryDecrRefCount  !decrements the reference count
          procedure, public:: get_ref_count=>TensCacheEntryGetRefCount    !returns the current reference count: Number of existing tensor operands associated with the tensor cache entry
          procedure, public:: incr_use_count=>TensCacheEntryIncrUseCount  !increments the use count
          procedure, public:: decr_use_count=>TensCacheEntryDecrUseCount  !decrements the use count
          procedure, public:: get_use_count=>TensCacheEntryGetUseCount    !returns the current use count: Number of active pointers to the tensor cache entry explicitly returned by the tensor cache
          procedure, public:: incr_read_count=>TensCacheEntryIncrReadCount    !increments the read count
          procedure, public:: decr_read_count=>TensCacheEntryDecrReadCount    !decrements the read count
          procedure, public:: get_read_count=>TensCacheEntryGetReadCount      !returns the current read count: Number of issued tensor instructions referring to this tensor cache entry as input
          procedure, public:: incr_write_count=>TensCacheEntryIncrWriteCount  !increments the write count
          procedure, public:: decr_write_count=>TensCacheEntryDecrWriteCount  !decrements the write count
          procedure, public:: get_write_count=>TensCacheEntryGetWriteCount    !returns the current write count: Number of issued tensor instructions referring to this tensor cache entry as output
          procedure, public:: get_rw_counter=>TensCacheEntryGetRwCounter      !returns the current value of the read/write access counter
          procedure, public:: reset_rw_counter=>TensCacheEntryResetRwCounter  !resets the read/write access counter
          procedure, public:: reset_rw_counters=>TensCacheEntryResetRwCounters!resets the read/write access counter by providing both new read and write counters separately
          procedure, public:: incr_temp_count=>TensCacheEntryIncrTempCount    !increments the temporary count (cannot be decremented)
          procedure, public:: get_temp_count=>TensCacheEntryGetTempCount      !returns the current temporary count
          procedure, public:: set_up_to_date=>TensCacheEntrySetUpToDate       !sets/resets the up-to-date status
          procedure, public:: is_up_to_date=>TensCacheEntryIsUpToDate         !returns whether or not the cache entry data is up-to-date
          procedure, public:: set_persistency=>TensCacheEntrySetPersistency   !sets/resets the persistency status
          procedure, public:: is_persistent=>TensCacheEntryIsPersistent       !returns TRUE if the tensor cache entry is persistent, FALSE otherwise
          procedure, public:: reset_all_counters=>TensCacheEntryResetAllCounters !resets all counters
          procedure, public:: destroy=>TensCacheEntryDestroy                  !destroys the tensor cache entry
          procedure, public:: init_lock=>TensCacheEntryInitLock               !initializes the tensor cache entry access lock
          procedure, public:: destroy_lock=>TensCacheEntryDestroyLock         !destroys the tensor cache entry access lock
          procedure, public:: lock=>TensCacheEntryLock                        !locks the tensor cache entry for an exclusive use/update
          procedure, public:: unlock=>TensCacheEntryUnlock                    !unlocks the tensor cache entry
          procedure, public:: is_lockable=>TensCacheEntryIsLockable           !returns TRUE if the tensor cache entry can be locked/unlocked
        end type tens_cache_entry_t
 !Tensor argument cache:
        type, public:: tens_cache_t
         type(dictionary_t), private:: map                         !cache dictionary: <tens_descr_t-->tens_cache_entry_t>
         integer(INTL), private:: num_entries=0                    !number of active entries in the tensor cache
#ifndef NO_OMP
         integer(omp_lock_kind), allocatable, private:: cache_lock !cache lock
#endif
         contains
          procedure, private:: init_lock=>TensCacheInitLock        !initializes the tensor cache lock for multithreading
          procedure, private:: lock=>TensCacheLock                 !locks the tensor cache
          procedure, private:: unlock=>TensCacheUnlock             !unlocks the tensor cache
          procedure, public:: lookup=>TensCacheLookup              !looks up a given tensor in the cache
          procedure, public:: store=>TensCacheStore                !stores a given tensor in the cache
          procedure, public:: evict=>TensCacheEvict                !evicts a given tensor from the cache
          procedure, public:: release_entry=>TensCacheReleaseEntry !releases a previously obtained pointer to a tensor cache entry
          procedure, public:: erase=>TensCacheErase                !erases everything from the cache (regardless of pending MPI communications!)
          procedure, public:: get_num_entries=>TensCacheGetNumEntries !returns the current number of active entries in the tensor cache
          procedure, public:: print_it=>TensCachePrintIt           !prints the content of the tensor cache
          final:: tens_cache_dtor                                  !dtor
        end type tens_cache_t
 !External data register:
        type, public:: data_register_t
         type(dictionary_t), private:: ext_data                           !string --> talsh_tens_data_t
         contains
          procedure, public:: register_data=>DataRegisterRegisterData     !registers external local data
          procedure, public:: unregister_data=>DataRegisterUnregisterData !unregisters previously registered data
          procedure, public:: retrieve_data=>DataRegisterRetrieveData     !retrieves previously registered data
          final:: data_register_dtor                                      !dtor
        end type data_register_t
 !External method register:
        type, public:: method_register_t
         type(dictionary_t), private:: ext_methods                              !string --> tens_method_uni_t
         contains
          procedure, public:: register_method=>MethodRegisterRegisterMethod     !registers an external tensor defining procedure
          procedure, public:: unregister_method=>MethodRegisterUnregisterMethod !unregisters a previously registered tensor defining procedure
          procedure, public:: retrieve_method=>MethodRegisterRetrieveMethod     !retrieves a previously registered tensor defining procedure
          final:: method_register_dtor                                          !dtor
        end type method_register_t
!INTERFACES:
        abstract interface
 !tens_cache_entry_t:
  !Polymorphic allocator:
         function tens_cache_entry_alloc_i(tens_cache_entry) result(ierr)
          import:: tens_cache_entry_t,INTD
          integer(INTD):: ierr
          class(tens_cache_entry_t), allocatable, intent(out):: tens_cache_entry
         end function tens_cache_entry_alloc_i
  !Printing:
         subroutine tens_cache_entry_print_i(this,ierr,dev_id,nspaces)
          import:: tens_cache_entry_t,INTD
          class(tens_cache_entry_t), intent(inout):: this !in: tensor cache entry
          integer(INTD), intent(out), optional:: ierr     !out: error code
          integer(INTD), intent(in), optional:: dev_id    !in: output device id
          integer(INTD), intent(in), optional:: nspaces   !in: left alignment
         end subroutine tens_cache_entry_print_i
 !External (user-defined) unary method (tensor initialization/transformation):
         function exatns_method_uni_i(tensor,scalar) result(ierr)
          import:: INTD,tens_rcrsv_t
          implicit none
          integer(INTD):: ierr                         !out: error code
          class(tens_rcrsv_t), intent(inout):: tensor  !inout: tensor argument
          complex(8), intent(inout), optional:: scalar !inout: scalar argument
         end function exatns_method_uni_i
 !External n-ary method (user-defined tensor operation):
         function exatns_method_i(tens_args,scal_args) result(ierr)
          import:: INTD,tens_rcrsv_t
          implicit none
          integer(INTD):: ierr                                !out: error code
          class(tens_rcrsv_t), intent(inout):: tens_args(0:)  !inout: tensor arguments
          complex(8), intent(inout), optional:: scal_args(0:) !inout: scalar arguments
         end function exatns_method_i
        end interface
        public tens_cache_entry_alloc_i
        public tens_cache_entry_print_i
        public exatns_method_uni_i
        public exatns_method_i
!GLOBAL DATA:
 !MPI process specialization (TAVP role, set by exatns_start):
        integer(INT_MPI), public:: process_role=EXA_NO_ROLE   !MPI process role (see above)
        integer(INT_MPI), public:: role_comm=MPI_COMM_NULL    !role-specific MPI communicator
        integer(INT_MPI), public:: role_size=0                !size of the role-specific MPI communicator
        integer(INT_MPI), public:: role_rank=-1               !process rank within the role-specific MPI communicator
        integer(INT_MPI), public:: driver_gl_rank=-1          !driver process rank in the global MPI communicator
        integer(INT_MPI), public:: top_manager_gl_rank=-1     !root manager process rank in the global MPI communicator
        integer(INT_MPI), public:: drv_mng_comm=MPI_COMM_NULL !MPI intercommunicator for the driver and managers
        integer(INT_MPI), public:: mng_wrk_comm=MPI_COMM_NULL !MPI intercommunicator for managers and workers
 !External universal tensor dimension strength assessing and shape resolution:
        procedure(tens_rcrsv_dim_resolve_i), pointer, public:: tens_dim_extent_resolve=>NULL() !resolves tensor dimension extents (determines actual shape of tensor blocks)
        procedure(tens_rcrsv_dim_strength_i), pointer, public:: tens_dim_strength_assess=>NULL() !assesses the strength of tensor dimensions
        real(8), public:: tens_dim_strength_thresh=0d0 !tensor dimension strength threshold above which the dimension will split (fine tensor decomposition granularity control)
 !External data register:
        type(data_register_t), public:: data_register     !string --> talsh_tens_data_t
 !External method register:
        type(method_register_t), public:: method_register !string --> tens_method_uni_t
!VISIBILITY:
 !non-member:
        public role_barrier
        public tensor_name_is_temporary
        public tensor_name_is_replica
        public tensor_name_mangle_temporary
        public tensor_name_unmangle_temporary
        public tensor_name_mangle_replica
        public tensor_name_unmangle_replica
        public opcode_control
        public opcode_auxiliary
        public opcode_tensor
 !ctrl_tens_trans_t:
        private CtrlTensTransCtor
        private CtrlTensTransMapMethod
        private CtrlTensTransPack
        private CtrlTensTransUnpack
        private CtrlTensTransPrintIt
        public ctrl_tens_trans_dtor
 !ctrl_tens_add_t:
        private CtrlTensAddCtor
        private CtrlTensAddGetPermutation
        private CtrlTensAddPack
        private CtrlTensAddUnpack
        private CtrlTensAddPrintIt
        public ctrl_tens_add_dtor
 !ctrl_tens_contr_t:
        private CtrlTensContrCtor
        private CtrlTensContrGetContrPtrn
        private CtrlTensContrPack
        private CtrlTensContrUnpack
        private CtrlTensContrPrintIt
        public ctrl_tens_contr_dtor
 !tens_cache_entry_t:
        private TensCacheEntryIsSet
        private TensCacheEntrySetTensor
        private TensCacheEntryGetTensor
        private TensCacheEntryIncrRefCount
        private TensCacheEntryDecrRefCount
        private TensCacheEntryGetRefCount
        private TensCacheEntryIncrUseCount
        private TensCacheEntryDecrUseCount
        private TensCacheEntryGetUseCount
        private TensCacheEntryIncrReadCount
        private TensCacheEntryDecrReadCount
        private TensCacheEntryGetReadCount
        private TensCacheEntryIncrWriteCount
        private TensCacheEntryDecrWriteCount
        private TensCacheEntryGetWriteCount
        private TensCacheEntryGetRwCounter
        private TensCacheEntryResetRwCounter
        private TensCacheEntryResetRwCounters
        private TensCacheEntryIncrTempCount
        private TensCacheEntryGetTempCount
        private TensCacheEntrySetUpToDate
        private TensCacheEntryIsUpToDate
        private TensCacheEntrySetPersistency
        private TensCacheEntryIsPersistent
        private TensCacheEntryResetAllCounters
        private TensCacheEntryDestroy
        private TensCacheEntryInitLock
        private TensCacheEntryDestroyLock
        private TensCacheEntryLock
        private TensCacheEntryUnlock
        private TensCacheEntryIsLockable
        public tens_cache_entry_print_f
 !tens_cache_t:
        private TensCacheInitLock
        private TensCacheLock
        private TensCacheUnlock
        private TensCacheLookup
        private TensCacheStore
        private TensCacheEvict
        private TensCacheReleaseEntry
        private TensCacheErase
        private TensCacheGetNumEntries
        private TensCachePrintIt
        public tens_cache_dtor
 !data_register_t:
        private DataRegisterRegisterData
        private DataRegisterUnregisterData
        private DataRegisterRetrieveData
        public data_register_dtor
 !method_register_t:
        private MethodRegisterRegisterMethod
        private MethodRegisterUnregisterMethod
        private MethodRegisterRetrieveMethod
        public method_register_dtor
        public method_map_f

       contains
!IMPLEMENTATION:
![non-member]========================
        subroutine role_barrier(ierr)
!MPI barrier for the role communicator.
         implicit none
         integer(INTD), intent(out), optional:: ierr
         integer(INT_MPI):: errc

         call MPI_Barrier(role_comm,errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine role_barrier
!-------------------------------------------------------------------------
        function tensor_name_is_temporary(tensor_name,ierr,id) result(res)
!Returns TRUE if the given tensor name is a name of a temporary tensor.
         implicit none
         logical:: res                               !out: answer
         character(*), intent(in):: tensor_name      !in: tensor name
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD), intent(out), optional:: id   !out: temporary tensor id
         integer(INTD):: errc,i,j

         errc=0
         i=index(tensor_name,'$',BACK=.TRUE.)
         j=index(tensor_name,'#',BACK=.TRUE.)
         if(j.gt.0) then
          if(j.gt.i) then
           res=.TRUE.
           if(present(id)) then
            i=len(tensor_name)-j
            if(i.gt.0) then
             id=icharnum(i,tensor_name(j+1:)); if(i.le.0) errc=-3
            else
             errc=-2
            endif
           endif
          else
           res=.FALSE.; errc=-1
          endif
         else
          res=.FALSE.
         endif
         if(present(ierr)) ierr=errc
         return
        end function tensor_name_is_temporary
!--------------------------------------------------------------------
        function tensor_name_is_replica(tensor_name,ierr) result(res)
!Returns TRUE if the given tensor name is a name of a replica tensor.
         implicit none
         logical:: res                               !out: answer
         character(*), intent(in):: tensor_name      !in: tensor name
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,i,j

         errc=0
         i=index(tensor_name,'#',BACK=.TRUE.)
         j=index(tensor_name,'$',BACK=.TRUE.)
         if(j.gt.0) then
          if(i.le.0.or.j.lt.i) then
           res=.TRUE.
          else
           res=.FALSE.; errc=-1
          endif
         else
          res=.FALSE.
         endif
         if(present(ierr)) ierr=errc
         return
        end function tensor_name_is_replica
!---------------------------------------------------------------------------------------------
        subroutine tensor_name_mangle_temporary(orig_name,new_name,new_name_len,ierr,instance)
!Mangles the tensor name to make it a name of a temporary tensor, if <instance> is present.
!If <instance> is absent, removes the temporary suffix from the tensor name.
!Trying to remove a non-existing temporary suffix or append more than one temporary suffix results in error.
!Mangling examples: Tens23 <--> Tens23#3, Tamp12$5 <--> Tamp12$5#2
         implicit none
         character(*), intent(in):: orig_name           !in: original tensor name
         character(*), intent(inout):: new_name         !out: new tensor name
         integer(INTD), intent(out):: new_name_len      !out: new tensor name length (must be long enough)
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD), intent(in), optional:: instance !in: temporary instance number (0:accumulator tensor; 1,2,3,..:temporary tensor)
         integer(INTD):: errc,l,i

         errc=0; new_name_len=0; l=len_trim(orig_name)
         if(l.gt.0) then
          if(present(instance)) then !persistent --> temporary name mangling
           if(instance.ge.0) then
            i=index(orig_name(1:l),'#')
            if(i.le.0.or.i.gt.l) then
             new_name(1:l+1)=orig_name(1:l)//'#'
             call numchar(instance,i,new_name(l+2:)); new_name_len=l+1+i
            else
             errc=-4
            endif
           else
            errc=-3
           endif
          else !temporary --> persistent name mangling
           i=l; do while(orig_name(i:i).ne.'#'); i=i-1; if(i.eq.0) exit; enddo
           if(i.gt.1) then
            new_name_len=i-1; new_name(1:new_name_len)=orig_name(1:new_name_len)
           else
            errc=-2
           endif
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine tensor_name_mangle_temporary
!--------------------------------------------------------------------------------------
        subroutine tensor_name_unmangle_temporary(orig_name,new_name,new_name_len,ierr)
         implicit none
         character(*), intent(in):: orig_name           !in: original tensor name
         character(*), intent(inout):: new_name         !out: new tensor name
         integer(INTD), intent(out):: new_name_len      !out: new tensor name length (must be long enough)
         integer(INTD), intent(out), optional:: ierr    !out: error code

         if(present(ierr)) then
          call tensor_name_mangle_temporary(orig_name,new_name,new_name_len,ierr)
         else
          call tensor_name_mangle_temporary(orig_name,new_name,new_name_len)
         endif
         return
        end subroutine tensor_name_unmangle_temporary
!-------------------------------------------------------------------------------------------
        subroutine tensor_name_mangle_replica(orig_name,new_name,new_name_len,ierr,instance)
!Mangles the tensor name to make it a name of a tensor replica, if <instance> is present.
!If <instance> is absent, removes the replica suffix from the tensor name.
!Trying to remove a non-existing replica suffix or append more than one replica suffix results in error.
!Mangling examples: Tens49 <--> Tens49$5, Tens49#1 <--> Tens49$5#1
         implicit none
         character(*), intent(in):: orig_name           !in: original tensor name
         character(*), intent(inout):: new_name         !out: new tensor name
         integer(INTD), intent(out):: new_name_len      !out: new tensor name length (must be long enough)
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD), intent(in), optional:: instance !in: replica instance number (>0)
         integer(INTD):: errc,l,i,k

         errc=0; new_name_len=0; l=len_trim(orig_name)
         if(l.gt.0) then
          if(present(instance)) then !non-replica --> replica name mangling
           if(instance.ge.0) then
            i=index(orig_name(1:l),'$')
            if(i.le.0.or.i.gt.l) then
             i=index(orig_name(1:l),'#')
             if(i.ge.1.and.i.le.l) then !has the temporary suffix
              if(i.gt.1.and.i.lt.l) then
               new_name(1:i)=orig_name(1:i-1)//'$'
               call numchar(instance,k,new_name(i+1:)); k=i+k
               new_name(k+1:k+(l-i+1))=orig_name(i:l); new_name_len=k+(l-i+1)
              else
               errc=-6
              endif
             else !does not have the temporary suffix
              new_name(1:l+1)=orig_name(1:l)//'$'
              call numchar(instance,i,new_name(l+2:)); new_name_len=l+1+i
             endif
            else
             errc=-5
            endif
           else
            errc=-4
           endif
          else !replica --> non-replica name mangling
           i=l; do while(orig_name(i:i).ne.'$'); i=i-1; if(i.eq.0) exit; enddo
           if(i.gt.1) then
            new_name_len=i-1; new_name(1:new_name_len)=orig_name(1:new_name_len)
            i=i+1; do while(i.le.l); if(is_it_number(orig_name(i:i)).lt.0) exit; i=i+1; enddo
            if(i.le.l) then !we still need to keep the temporary part
             if(orig_name(i:i).eq.'#'.and.i.lt.l) then
              new_name(new_name_len+1:new_name_len+(l-i+1))=orig_name(i:l); new_name_len=new_name_len+(l-i+1)
             else
              errc=-3
             endif
            endif
           else
            errc=-2
           endif
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine tensor_name_mangle_replica
!------------------------------------------------------------------------------------
        subroutine tensor_name_unmangle_replica(orig_name,new_name,new_name_len,ierr)
         implicit none
         character(*), intent(in):: orig_name           !in: original tensor name
         character(*), intent(inout):: new_name         !out: new tensor name
         integer(INTD), intent(out):: new_name_len      !out: new tensor name length (must be long enough)
         integer(INTD), intent(out), optional:: ierr    !out: error code

         if(present(ierr)) then
          call tensor_name_mangle_replica(orig_name,new_name,new_name_len,ierr)
         else
          call tensor_name_mangle_replica(orig_name,new_name,new_name_len)
         endif
         return
        end subroutine tensor_name_unmangle_replica
!-------------------------------------------------------
        function opcode_control(opcode,ierr) result(ans)
         implicit none
         logical:: ans                               !out: answer
         integer(INTD), intent(in):: opcode          !in: opcode
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; ans=.FALSE.
         if(opcode.ge.0.and.opcode.lt.TAVP_ISA_SIZE) then
          ans=(opcode.ge.TAVP_ISA_CTRL_FIRST.and.opcode.le.TAVP_ISA_CTRL_LAST)
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function opcode_control
!---------------------------------------------------------
        function opcode_auxiliary(opcode,ierr) result(ans)
         implicit none
         logical:: ans                               !out: answer
         integer(INTD), intent(in):: opcode          !in: opcode
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; ans=.FALSE.
         if(opcode.ge.0.and.opcode.lt.TAVP_ISA_SIZE) then
          ans=(opcode.ge.TAVP_ISA_SPACE_FIRST.and.opcode.le.TAVP_ISA_SPACE_LAST)
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function opcode_auxiliary
!------------------------------------------------------
        function opcode_tensor(opcode,ierr) result(ans)
         implicit none
         logical:: ans                               !out: answer
         integer(INTD), intent(in):: opcode          !in: opcode
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; ans=.FALSE.
         if(opcode.ge.0.and.opcode.lt.TAVP_ISA_SIZE) then
          ans=(opcode.ge.TAVP_ISA_TENS_FIRST.and.opcode.le.TAVP_ISA_TENS_LAST)
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function opcode_tensor
![ctrl_tens_trans_t]=============================================
        subroutine CtrlTensTransCtor(this,ierr,alpha,method_name)
!CTOR: If neither <alpha> nor <method_name> is present, initialization
!to zero is assumed. If either <alpha> or <method_name> is present,
!then the initialization to a scalar or a user-defined transformation
!is assumed, resepectively.
         implicit none
         class(ctrl_tens_trans_t), intent(out):: this     !out: tensor transformation/initialization control field
         integer(INTD), intent(out), optional:: ierr      !out: error code
         complex(8), intent(in), optional:: alpha         !in: scalar
         character(*), intent(in), optional:: method_name !in: external (registered) method name
         integer(INTD):: errc

         errc=0; this%method_name=' '
         if(present(alpha)) this%alpha=alpha
         if(present(method_name)) then
          if(len_trim(method_name).le.EXA_MAX_METHOD_NAME_LEN) then
           this%method_name=method_name(1:len_trim(method_name))
          else
           errc=-1
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine CtrlTensTransCtor
!---------------------------------------------------
        subroutine CtrlTensTransMapMethod(this,ierr)
!Maps the method name to the actual definer object.
         implicit none
         class(ctrl_tens_trans_t), intent(inout):: this !inout: tensor transformation/initialization control field
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD):: errc,l

         errc=0; l=len_trim(this%method_name)
         if(l.gt.0) then
          this%definer=>method_register%retrieve_method(this%method_name(1:l),errc); if(errc.ne.0) errc=-1
         else
          errc=-2
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine CtrlTensTransMapMethod
!-----------------------------------------------------
        subroutine CtrlTensTransPack(this,packet,ierr)
!Packer.
         implicit none
         class(ctrl_tens_trans_t), intent(in):: this !in: tensor transformation/initialization control field
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         call pack_builtin(packet,this%method_name(1:len_trim(this%method_name)),errc)
         if(errc.eq.PACK_SUCCESS) then
          call pack_builtin(packet,this%alpha,errc); if(errc.ne.PACK_SUCCESS) errc=-1
         else
          errc=-2
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine CtrlTensTransPack
!-------------------------------------------------------
        subroutine CtrlTensTransUnpack(this,packet,ierr)
!Unpacker.
         implicit none
         class(ctrl_tens_trans_t), intent(inout):: this !out: tensor transformation/initialization control field
         class(obj_pack_t), intent(inout):: packet      !inout: packet
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD):: errc
         integer(INTL):: sl

         this%method_name=' '
         call unpack_builtin(packet,this%method_name,sl,errc)
         if(errc.eq.PACK_SUCCESS) then
          if(sl.le.EXA_MAX_METHOD_NAME_LEN) then
           call unpack_builtin(packet,this%alpha,errc); if(errc.ne.PACK_SUCCESS) errc=-1
          else
           errc=-2
          endif
         else
          errc=-3
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine CtrlTensTransUnpack
!-----------------------------------------------------------------
        subroutine CtrlTensTransPrintIt(this,ierr,dev_id,nspaces)
         implicit none
         class(ctrl_tens_trans_t), intent(in):: this
         integer(INTD), intent(out), optional:: ierr
         integer(INTD), intent(in), optional:: dev_id
         integer(INTD), intent(in), optional:: nspaces
         integer(INTD):: errc,devo,nsp,l
!$OMP FLUSH
         errc=0
         devo=6; if(present(dev_id)) devo=dev_id
         nsp=0; if(present(nspaces)) nsp=nspaces
         l=len_trim(this%method_name)
         if(l.gt.0) call printl(devo,'Method: '//this%method_name(1:l))
         write(devo,'("Scalar: ",D24.14,1x,D24.14)') this%alpha
         flush(devo)
         if(present(ierr)) ierr=errc
         return
        end subroutine CtrlTensTransPrintIt
!--------------------------------------------
        subroutine ctrl_tens_trans_dtor(this)
!DTOR.
         implicit none
         type(ctrl_tens_trans_t):: this

         this%definer=>NULL()
         this%method_name=' '
         return
        end subroutine ctrl_tens_trans_dtor
![ctrl_tens_add_t]============================================
        subroutine CtrlTensAddCtor(this,permutation,ierr,alpha)
!CTOR.
         implicit none
         class(ctrl_tens_add_t), intent(out):: this    !out: tensor addition/copy control field
         type(permutation_t), intent(in):: permutation !in: tensor dimension permutation
         integer(INTD), intent(out), optional:: ierr   !out: error code
         complex(8), intent(in), optional:: alpha      !in: scalar prefactor
         integer(INTD):: errc

         errc=0
         this%permutation=permutation
         if(present(alpha)) this%alpha=alpha
         if(present(ierr)) ierr=errc
         return
        end subroutine CtrlTensAddCtor
!-------------------------------------------------------------------------
        function CtrlTensAddGetPermutation(this,ierr,alpha) result(permut)
!Returns a non-owning pointer to the encapsulated permutation.
         implicit none
         type(permutation_t), pointer:: permut             !out: pointer to the permutation
         class(ctrl_tens_add_t), intent(in), target:: this !in: tensor addition/copy control field
         integer(INTD), intent(out), optional:: ierr       !out: error code
         complex(8), intent(out), optional:: alpha         !out: scalar prefactor
         integer(INTD):: errc

         errc=0
         permut=>this%permutation
         if(present(alpha)) alpha=this%alpha
         if(present(ierr)) ierr=errc
         return
        end function CtrlTensAddGetPermutation
!---------------------------------------------------
        subroutine CtrlTensAddPack(this,packet,ierr)
!Packer.
         implicit none
         class(ctrl_tens_add_t), intent(in):: this   !in: tensor addition/copy control field
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         call this%permutation%pack(packet,errc)
         if(errc.eq.TEREC_SUCCESS) then
          call pack_builtin(packet,this%alpha,errc); if(errc.ne.PACK_SUCCESS) errc=-1
         else
          errc=-2
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine CtrlTensAddPack
!-----------------------------------------------------
        subroutine CtrlTensAddUnpack(this,packet,ierr)
!Unpacker.
         implicit none
         class(ctrl_tens_add_t), intent(inout):: this !in: tensor addition/copy control field
         class(obj_pack_t), intent(inout):: packet    !inout: packet
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         call this%permutation%unpack(packet,errc)
         if(errc.eq.TEREC_SUCCESS) then
          call unpack_builtin(packet,this%alpha,errc); if(errc.ne.PACK_SUCCESS) errc=-1
         else
          errc=-2
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine CtrlTensAddUnpack
!--------------------------------------------------------------
        subroutine CtrlTensAddPrintIt(this,ierr,dev_id,nspaces)
         implicit none
         class(ctrl_tens_add_t), intent(in):: this
         integer(INTD), intent(out), optional:: ierr
         integer(INTD), intent(in), optional:: dev_id
         integer(INTD), intent(in), optional:: nspaces
         integer(INTD):: errc,devo,nsp,i
!$OMP FLUSH
         errc=0
         devo=6; if(present(dev_id)) devo=dev_id
         nsp=0; if(present(nspaces)) nsp=nspaces
         do i=1,nsp; write(devo,'(" ")',ADVANCE='NO'); enddo
         write(devo,'("Permutation: ")',ADVANCE='NO')
         call this%permutation%print_it(errc,devo,0)
         if(errc.eq.TEREC_SUCCESS) then
          do i=1,nsp; write(devo,'(" ")',ADVANCE='NO'); enddo
          write(devo,'("Scalar: ",(D24.14,1x,D24.14))') this%alpha
         else
          errc=-1
         endif
         flush(devo)
         if(present(ierr)) ierr=errc
         return
        end subroutine CtrlTensAddPrintIt
!------------------------------------------
        subroutine ctrl_tens_add_dtor(this)
         implicit none
         type(ctrl_tens_add_t):: this

         return
        end subroutine ctrl_tens_add_dtor
![ctrl_tens_contr_t]============================================
        subroutine CtrlTensContrCtor(this,contr_ptrn,ierr,alpha)
!CTOR.
         implicit none
         class(ctrl_tens_contr_t), intent(out):: this    !out: tensor contraction control field
         type(contr_ptrn_ext_t), intent(in):: contr_ptrn !in: extended tensor contraction pattern
         integer(INTD), intent(out), optional:: ierr     !out: error code
         complex(8), intent(in), optional:: alpha        !in: scalar prefactor
         integer(INTD):: errc

         if(contr_ptrn%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) then
           this%contr_ptrn=contr_ptrn
           if(present(alpha)) this%alpha=alpha
          else
           errc=-1
          endif
         else
          errc=-2
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine CtrlTensContrCtor
!------------------------------------------------------------------------------
        function CtrlTensContrGetContrPtrn(this,ierr,alpha,conjug) result(ptrn)
!Returns a pointer to the extended tensor contraction pattern. Optionally,
!the prefactor and tensor argument conjugation flags can be returned as well.
         implicit none
         type(contr_ptrn_ext_t), pointer:: ptrn              !out: extended tensor contraction pattern
         class(ctrl_tens_contr_t), target, intent(in):: this !in: tensor contraction control field
         integer(INTD), intent(out), optional:: ierr         !out: error code
         complex(8), intent(out), optional:: alpha           !out: tensor contraction prefactor
         integer(INTD), intent(out), optional:: conjug       !out: tensor argument conjugation bits: {0:D,1:L,2:R}
         integer(INTD):: errc

         ptrn=>NULL()
         if(this%contr_ptrn%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) then
           ptrn=>this%contr_ptrn
           if(present(conjug)) then
            conjug=this%contr_ptrn%get_conjugation(errc); if(errc.ne.TEREC_SUCCESS) errc=-3
           endif
           if(present(alpha)) alpha=this%alpha
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function CtrlTensContrGetContrPtrn
!-----------------------------------------------------
        subroutine CtrlTensContrPack(this,packet,ierr)
!Packer.
         implicit none
         class(ctrl_tens_contr_t), intent(in):: this !in: tensor contraction control field
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%contr_ptrn%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) then
           call this%contr_ptrn%pack(packet,errc)
           if(errc.eq.TEREC_SUCCESS) then
            call pack_builtin(packet,this%alpha,errc); if(errc.ne.0) errc=-1
           else
            errc=-2
           endif
          else
           errc=-3
          endif
         else
          errc=-4
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine CtrlTensContrPack
!-------------------------------------------------------
        subroutine CtrlTensContrUnpack(this,packet,ierr)
!Unpacker.
         implicit none
         class(ctrl_tens_contr_t), intent(inout):: this !out: tensor contraction control field
         class(obj_pack_t), intent(inout):: packet      !inout: packet
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD):: errc

         call this%contr_ptrn%unpack(packet,errc)
         if(errc.eq.TEREC_SUCCESS) then
          call unpack_builtin(packet,this%alpha,errc); if(errc.ne.0) errc=-1
         else
          errc=-2
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine CtrlTensContrUnpack
!----------------------------------------------------------------
        subroutine CtrlTensContrPrintIt(this,ierr,dev_id,nspaces)
         implicit none
         class(ctrl_tens_contr_t), intent(in):: this
         integer(INTD), intent(out), optional:: ierr
         integer(INTD), intent(in), optional:: dev_id
         integer(INTD), intent(in), optional:: nspaces
         integer(INTD):: errc,devo,nsp,i
!$OMP FLUSH
         errc=0
         devo=6; if(present(dev_id)) devo=dev_id
         nsp=0; if(present(nspaces)) nsp=nspaces
         do i=1,nsp; write(devo,'(" ")',ADVANCE='NO'); enddo
         write(devo,'("Pattern: ")',ADVANCE='NO')
         call this%contr_ptrn%print_it(errc,devo,0)
         if(errc.eq.TEREC_SUCCESS) then
          do i=1,nsp; write(devo,'(" ")',ADVANCE='NO'); enddo
          write(devo,'("Scalar: ",D24.14,1x,D24.14)') this%alpha
         else
          errc=-1
         endif
         flush(devo)
         if(present(ierr)) ierr=errc
         return
        end subroutine CtrlTensContrPrintIt
!--------------------------------------------
        subroutine ctrl_tens_contr_dtor(this)
         implicit none
         type(ctrl_tens_contr_t):: this

         return
        end subroutine ctrl_tens_contr_dtor
![tens_cache_entry_t]======================================
        function TensCacheEntryIsSet(this,ierr) result(ans)
!Returns TRUE if the tensor is set, FALSE otherwise.
         implicit none
         logical:: ans                                   !out: answer
         class(tens_cache_entry_t), intent(inout):: this !in: tensor cache entry
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         errc=0
         call this%lock()
         ans=associated(this%tensor)
         call this%unlock()
         if(present(ierr)) ierr=errc
         return
        end function TensCacheEntryIsSet
!-----------------------------------------------------------
        subroutine TensCacheEntrySetTensor(this,tensor,ierr)
!Sets a pointer to the tensor, either owning or non-owning. The ownership
!assumption is ultimately determined by dtor's <dealloc> argument.
!If the <tensor> pointer is owning and the ownership transfer is implied,
!dtor should have its argument <dealloc>=TRUE at destruction.
         implicit none
         class(tens_cache_entry_t), intent(inout):: this   !inout: empty tensor cache entry
         class(tens_rcrsv_t), intent(in), pointer:: tensor !in: associated pointer to a tensor, either owning or non-owning
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc

         errc=0
         call this%lock()
         if((.not.associated(this%tensor)).and.associated(tensor)) then
          this%tensor=>tensor
         else
          errc=-1
         endif
         call this%unlock()
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCacheEntrySetTensor
!-------------------------------------------------------------------
        function TensCacheEntryGetTensor(this,ierr) result(tensor_p)
         implicit none
         class(tens_rcrsv_t), pointer:: tensor_p         !out: pointer to the tensor
         class(tens_cache_entry_t), intent(inout):: this !in: tensor cache entry
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         errc=0
         !call this%lock() !not needed because the lock will be set at the upper level
         tensor_p=>this%tensor
         !call this%unlock() !not needed because the lock will be unset at the upper level
         if(.not.associated(tensor_p)) errc=-1
         if(present(ierr)) ierr=errc
         return
        end function TensCacheEntryGetTensor
!--------------------------------------------------
        subroutine TensCacheEntryIncrRefCount(this)
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry

!$OMP ATOMIC UPDATE SEQ_CST
         this%ref_count=this%ref_count+1
         return
        end subroutine TensCacheEntryIncrRefCount
!--------------------------------------------------
        subroutine TensCacheEntryDecrRefCount(this)
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry

!$OMP ATOMIC UPDATE SEQ_CST
         this%ref_count=this%ref_count-1
         return
        end subroutine TensCacheEntryDecrRefCount
!----------------------------------------------------------------
        function TensCacheEntryGetRefCount(this) result(refcount)
         implicit none
         integer(INTD):: refcount                     !out: reference count (number of active tensor operands associated with this tensor cache entry)
         class(tens_cache_entry_t), intent(in):: this !in: defined tensor cache entry

!$OMP ATOMIC READ SEQ_CST
         refcount=this%ref_count
         return
        end function TensCacheEntryGetRefCount
!--------------------------------------------------
        subroutine TensCacheEntryIncrUseCount(this)
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry

!$OMP ATOMIC UPDATE SEQ_CST
         this%use_count=this%use_count+1
         return
        end subroutine TensCacheEntryIncrUseCount
!--------------------------------------------------
        subroutine TensCacheEntryDecrUseCount(this)
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry

!$OMP ATOMIC UPDATE SEQ_CST
         this%use_count=this%use_count-1
         return
        end subroutine TensCacheEntryDecrUseCount
!----------------------------------------------------------------
        function TensCacheEntryGetUseCount(this) result(usecount)
         implicit none
         integer(INTD):: usecount                     !out: use count (number of active pointers explicitly returned by cache associated with this tensor cache entry)
         class(tens_cache_entry_t), intent(in):: this !in: defined tensor cache entry

!$OMP ATOMIC READ SEQ_CST
         usecount=this%use_count
         return
        end function TensCacheEntryGetUseCount
!--------------------------------------------------------------
        subroutine TensCacheEntryIncrReadCount(this,ierr,defer)
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry
         integer(INTD), intent(out), optional:: ierr     !out: error code
         logical, intent(in), optional:: defer           !in: deferred or actual read count
         integer(INTD):: errc
         logical:: actual

         actual=.TRUE.; if(present(defer)) actual=(.not.defer)
         if(actual) then
!$OMP ATOMIC CAPTURE SEQ_CST
          errc=this%read_write_count
          this%read_write_count=this%read_write_count+1
!$OMP END ATOMIC
         else
!$OMP ATOMIC CAPTURE SEQ_CST
          errc=this%read_write_def_count
          this%read_write_def_count=this%read_write_def_count+1
!$OMP END ATOMIC
         endif
         if(errc.ge.0) then; errc=0; else; errc=-1; endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCacheEntryIncrReadCount
!--------------------------------------------------------------
        subroutine TensCacheEntryDecrReadCount(this,ierr,defer)
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry
         integer(INTD), intent(out), optional:: ierr     !out: error code
         logical, intent(in), optional:: defer           !in: deferred or actual read count
         integer(INTD):: errc
         logical:: actual

         actual=.TRUE.; if(present(defer)) actual=(.not.defer)
         if(actual) then
!$OMP ATOMIC CAPTURE SEQ_CST
          errc=this%read_write_count
          this%read_write_count=this%read_write_count-1
!$OMP END ATOMIC
         else
!$OMP ATOMIC CAPTURE SEQ_CST
          errc=this%read_write_def_count
          this%read_write_def_count=this%read_write_def_count-1
!$OMP END ATOMIC
         endif
         if(errc.gt.0) then; errc=0; else; errc=-1; endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCacheEntryDecrReadCount
!--------------------------------------------------------------------
        function TensCacheEntryGetReadCount(this,defer) result(count)
         implicit none
         integer(INTD):: count                        !out: read count
         class(tens_cache_entry_t), intent(in):: this !in: defined tensor cache entry
         logical, intent(in), optional:: defer        !in: deferred or actual read count
         logical:: actual

         actual=.TRUE.; if(present(defer)) actual=(.not.defer)
         if(actual) then
!$OMP ATOMIC READ SEQ_CST
          count=this%read_write_count
         else
!$OMP ATOMIC READ SEQ_CST
          count=this%read_write_def_count
         endif
         count=max(count,0)
         return
        end function TensCacheEntryGetReadCount
!---------------------------------------------------------------
        subroutine TensCacheEntryIncrWriteCount(this,ierr,defer)
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry
         integer(INTD), intent(out), optional:: ierr     !out: error code
         logical, intent(in), optional:: defer           !in: deferred or actual write count
         integer(INTD):: errc
         logical:: actual

         actual=.TRUE.; if(present(defer)) actual=(.not.defer)
         if(actual) then
!$OMP ATOMIC CAPTURE SEQ_CST
          errc=this%read_write_count
          this%read_write_count=this%read_write_count-1
!$OMP END ATOMIC
         else
!$OMP ATOMIC CAPTURE SEQ_CST
          errc=this%read_write_def_count
          this%read_write_def_count=this%read_write_def_count-1
!$OMP END ATOMIC
         endif
         if(errc.le.0) then; errc=0; else; errc=-1; endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCacheEntryIncrWriteCount
!---------------------------------------------------------------
        subroutine TensCacheEntryDecrWriteCount(this,ierr,defer)
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry
         integer(INTD), intent(out), optional:: ierr     !out: error code
         logical, intent(in), optional:: defer           !in: deferred or actual write count
         integer(INTD):: errc
         logical:: actual

         actual=.TRUE.; if(present(defer)) actual=(.not.defer)
         if(actual) then
!$OMP ATOMIC CAPTURE SEQ_CST
          errc=this%read_write_count
          this%read_write_count=this%read_write_count+1
!$OMP END ATOMIC
         else
!$OMP ATOMIC CAPTURE SEQ_CST
          errc=this%read_write_def_count
          this%read_write_def_count=this%read_write_def_count+1
!$OMP END ATOMIC
         endif
         if(errc.lt.0) then; errc=0; else; errc=-1; endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCacheEntryDecrWriteCount
!---------------------------------------------------------------------
        function TensCacheEntryGetWriteCount(this,defer) result(count)
         implicit none
         integer(INTD):: count                        !out: write count
         class(tens_cache_entry_t), intent(in):: this !in: defined tensor cache entry
         logical, intent(in), optional:: defer        !in: deferred or actual write count
         logical:: actual

         actual=.TRUE.; if(present(defer)) actual=(.not.defer)
         if(actual) then
!$OMP ATOMIC READ SEQ_CST
          count=this%read_write_count
         else
!$OMP ATOMIC READ SEQ_CST
          count=this%read_write_def_count
         endif
         count=max(-count,0)
         return
        end function TensCacheEntryGetWriteCount
!-------------------------------------------------------------------------------
        function TensCacheEntryGetRwCounter(this,defer) result(read_write_count)
         implicit none
         integer(INTD):: read_write_count              !out: current read/write access count
         class(tens_cache_entry_t), intent(in):: this  !in: defined tensor cache entry
         logical, intent(in), optional:: defer         !in: deferred or actual read/write count
         logical:: actual

         actual=.TRUE.; if(present(defer)) actual=(.not.defer)
         if(actual) then
!$OMP ATOMIC READ SEQ_CST
          read_write_count=this%read_write_count
         else
!$OMP ATOMIC READ SEQ_CST
          read_write_count=this%read_write_def_count
         endif
         return
        end function TensCacheEntryGetRwCounter
!---------------------------------------------------------------------------
        subroutine TensCacheEntryResetRwCounter(this,read_write_count,defer)
         implicit none
         class(tens_cache_entry_t), intent(inout):: this        !inout: defined tensor cache entry
         integer(INTD), intent(in), optional:: read_write_count !in: new read/write access count (or none)
         logical, intent(in), optional:: defer                  !in: deferred or actual read/write count
         logical:: actual

         actual=.TRUE.; if(present(defer)) actual=(.not.defer)
         if(actual) then
          if(present(read_write_count)) then
!$OMP ATOMIC WRITE SEQ_CST
           this%read_write_count=read_write_count
          else
!$OMP ATOMIC WRITE SEQ_CST
           this%read_write_count=0
          endif
         else
          if(present(read_write_count)) then
!$OMP ATOMIC WRITE SEQ_CST
           this%read_write_def_count=read_write_count
          else
!$OMP ATOMIC WRITE SEQ_CST
           this%read_write_def_count=0
          endif
         endif
         return
        end subroutine TensCacheEntryResetRwCounter
!---------------------------------------------------------------------------------------
        subroutine TensCacheEntryResetRwCounters(this,read_count,write_count,ierr,defer)
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry
         integer(INTD), intent(in):: read_count          !in: new read count
         integer(INTD), intent(in):: write_count         !in: new write count
         integer(INTD), intent(out), optional:: ierr     !out: error code
         logical, intent(in), optional:: defer           !in: deferred or actual read/write count
         integer(INTD):: errc
         logical:: actual

         if(.not.(read_count.gt.0.and.write_count.gt.0)) then !race
          actual=.TRUE.; if(present(defer)) actual=(.not.defer)
          if(actual) then
!$OMP ATOMIC WRITE SEQ_CST
           this%read_write_count=read_count-write_count
          else
!$OMP ATOMIC WRITE SEQ_CST
           this%read_write_def_count=read_count-write_count
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCacheEntryResetRwCounters
!---------------------------------------------------
        subroutine TensCacheEntryIncrTempCount(this)
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry

!$OMP ATOMIC UPDATE SEQ_CST
         this%temp_count=this%temp_count+1
         return
        end subroutine TensCacheEntryIncrTempCount
!--------------------------------------------------------------
        function TensCacheEntryGetTempCount(this) result(count)
         implicit none
         integer(INTD):: count                        !out: temporary count
         class(tens_cache_entry_t), intent(in):: this !in: defined tensor cache entry

!$OMP ATOMIC READ SEQ_CST
         count=this%temp_count
         return
        end function TensCacheEntryGetTempCount
!-----------------------------------------------------------------
        subroutine TensCacheEntrySetUpToDate(this,up_to_date_stat)
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry
         logical, intent(in):: up_to_date_stat           !in: new up-to-date status

!$OMP ATOMIC WRITE SEQ_CST
         this%up_to_date=up_to_date_stat
         return
        end subroutine TensCacheEntrySetUpToDate
!----------------------------------------------------------
        function TensCacheEntryIsUpToDate(this) result(res)
         implicit none
         logical:: res                                !out: result
         class(tens_cache_entry_t), intent(in):: this !in: defined tensor cache entry

!$OMP ATOMIC READ SEQ_CST
         res=this%up_to_date
         return
        end function TensCacheEntryIsUpToDate
!----------------------------------------------------------------
        subroutine TensCacheEntrySetPersistency(this,persistency)
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: defined tensor cache entry
         logical, intent(in):: persistency               !in: new persistency status

!$OMP ATOMIC WRITE SEQ_CST
         this%persistent=persistency
         return
        end subroutine TensCacheEntrySetPersistency
!------------------------------------------------------------
        function TensCacheEntryIsPersistent(this) result(res)
         implicit none
         logical:: res                                !out: result
         class(tens_cache_entry_t), intent(in):: this !in: defined tensor cache entry

!$OMP ATOMIC READ SEQ_CST
         res=this%persistent
         return
        end function TensCacheEntryIsPersistent
!------------------------------------------------------
        subroutine TensCacheEntryResetAllCounters(this)
!Resets all counters.
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: tensor cache entry

         this%ref_count=0
         this%use_count=0
         this%read_write_count=0
         this%read_write_def_count=0
         this%temp_count=0
         this%up_to_date=.FALSE.
         this%persistent=.FALSE.
         return
        end subroutine TensCacheEntryResetAllCounters
!----------------------------------------------------------
        subroutine TensCacheEntryDestroy(this,dealloc,ierr)
!Destroys the tensor cache entry, but only if the reference count and use count
!are both zero and the tensor cache entry is not persistent, otherwise raises error.
!Since destruction is only performed during eviction of the tensor cache entry,
!it is already race protected by the tensor cache lock.
         implicit none
         class(tens_cache_entry_t), intent(inout):: this !inout: temporary tensor cache entry with no references to it
         logical, intent(in):: dealloc                   !in: if TRUE, the .tensor field will be deallocated (assumes ownership)
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc
!$OMP FLUSH
         errc=0
         if(this%ref_count.eq.0.and.this%use_count.eq.0.and.(.not.this%persistent)) then
          if(associated(this%tensor).and.dealloc) deallocate(this%tensor) !<dealloc>=TRUE assumes tensor ownership
          this%tensor=>NULL()
          call this%reset_all_counters()
          call this%destroy_lock()
         else
          errc=-1
          write(jo,'("#ERROR(TensorCache:tens_cache_entry_t.destroy): Attempt to destroy an active tensor cache entry!")')
          write(jo,'("RefCount = ",i11,", UseCount = ",i11,", Persistency = ",l1)') this%ref_count,this%use_count,this%persistent
          call this%tensor%print_it(dev_id=jo); flush(jo)
          !call crash() !debug
          call quit(errc,'#ERROR(TensorCache:tens_cache_entry_t.destroy): Attempt to destroy an active tensor cache entry!')
         endif
!$OMP FLUSH
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCacheEntryDestroy
!----------------------------------------------
        subroutine TensCacheEntryInitLock(this)
         implicit none
         class(tens_cache_entry_t), intent(inout):: this
#ifndef NO_OMP
!$OMP FLUSH(this)
         if(.not.this%lock_initialized) then
          if(DEBUG.gt.1) then
           write(jo,'("New cache entry lock: ",i18," --> ")',ADVANCE='NO') this%entry_lock; flush(jo)
          endif
          call omp_init_nest_lock(this%entry_lock)
!$OMP ATOMIC WRITE SEQ_CST
          this%lock_initialized=.TRUE.
          if(DEBUG.gt.1) then; write(jo,'(i18)') this%entry_lock; flush(jo); endif
         endif
#endif
         return
        end subroutine TensCacheEntryInitLock
!-------------------------------------------------
        subroutine TensCacheEntryDestroyLock(this)
         implicit none
         class(tens_cache_entry_t), intent(inout):: this
#ifndef NO_OMP
!$OMP FLUSH(this)
         if(this%lock_initialized) then
          if(DEBUG.gt.1) then
           write(jo,'("Destroying cache entry lock ",i18," ... ")',ADVANCE='NO') this%entry_lock; flush(jo)
          endif
          call omp_destroy_nest_lock(this%entry_lock)
!$OMP ATOMIC WRITE SEQ_CST
          this%lock_initialized=.FALSE.
          if(DEBUG.gt.1) then; write(jo,'("Done")'); flush(jo); endif
         endif
#endif
         return
        end subroutine TensCacheEntryDestroyLock
!------------------------------------------
        subroutine TensCacheEntryLock(this)
         implicit none
         class(tens_cache_entry_t), intent(inout):: this
#ifndef NO_OMP
         if(this%lock_initialized) then
          call omp_set_nest_lock(this%entry_lock)
!$OMP FLUSH
         else
          write(jo,'("#FATAL(VIRTA:tens_cache_entry_t.lock): Attempt to set an uninitialized lock for tensor:")')
          if(associated(this%tensor)) call this%tensor%print_it(dev_id=jo); flush(jo)
          call quit(-1,'#FATAL(VIRTA:tens_cache_entry_t.lock): Attempt to set an uninitialized lock!')
         endif
#endif
         return
        end subroutine TensCacheEntryLock
!--------------------------------------------
        subroutine TensCacheEntryUnlock(this)
         implicit none
         class(tens_cache_entry_t), intent(inout):: this
#ifndef NO_OMP
         if(this%lock_initialized) then
!$OMP FLUSH
          call omp_unset_nest_lock(this%entry_lock)
         else
          write(jo,'("#FATAL(VIRTA:tens_cache_entry_t.unlock): Attempt to unset an uninitialized lock for tensor:")')
          if(associated(this%tensor)) call this%tensor%print_it(dev_id=jo); flush(jo)
          call quit(-1,'#FATAL(VIRTA:tens_cache_entry_t.lock): Attempt to unset an uninitialized lock!')
         endif
#endif
         return
        end subroutine TensCacheEntryUnlock
!---------------------------------------------------------------
        function TensCacheEntryIsLockable(this) result(lockable)
         implicit none
         logical:: lockable
         class(tens_cache_entry_t), intent(in):: this
#ifndef NO_OMP
!$OMP ATOMIC READ SEQ_CST
         lockable=this%lock_initialized
#else
         lockable=.FALSE.
#endif
         return
        end function TensCacheEntryIsLockable
!----------------------------------------------------------
        function tens_cache_entry_print_f(obj) result(ierr)
!Non-member printing for tens_cache_entry_t compatible with GFC.
         implicit none
         integer(INTD):: ierr                  !out: error code
         class(*), intent(inout), target:: obj !in: tensor cache entry

         ierr=0
         select type(obj)
         class is(tens_cache_entry_t)
          call obj%print_it(ierr,dev_id=jo)
         end select
         return
        end function tens_cache_entry_print_f
![tens_cache_t]===========================
        subroutine TensCacheInitLock(this)
!Initializes the tensor cache lock for multithreading.
         implicit none
         class(tens_cache_t), intent(inout):: this !inout: tensor cache
#ifndef NO_OMP
!$OMP CRITICAL (TAVP_CACHE)
         if(.not.allocated(this%cache_lock)) then
          allocate(this%cache_lock)
          call omp_init_lock(this%cache_lock)
         endif
!$OMP END CRITICAL (TAVP_CACHE)
#endif
         return
        end subroutine TensCacheInitLock
!-------------------------------------
        subroutine TensCacheLock(this)
!Locks the tensor cache.
         implicit none
         class(tens_cache_t), intent(inout):: this !inout: tensor cache
#ifndef NO_OMP
         call omp_set_lock(this%cache_lock)
!$OMP FLUSH
#endif
         return
        end subroutine TensCacheLock
!---------------------------------------
        subroutine TensCacheUnlock(this)
!Unlocks the tensor cache.
         implicit none
         class(tens_cache_t), intent(inout):: this !inout: tensor cache
#ifndef NO_OMP
!$OMP FLUSH
          call omp_unset_lock(this%cache_lock)
#endif
         return
        end subroutine TensCacheUnlock
!----------------------------------------------------------------------
        function TensCacheLookup(this,tensor,ierr) result(tens_entry_p)
!Looks up a given tensor in the tensor cache. If found, returns a pointer
!to the corresponding tensor cache entry. If not found, returns NULL.
         implicit none
         class(tens_cache_entry_t), pointer:: tens_entry_p !out: pointer to the found tensor cache entry or NULL
         class(tens_cache_t), intent(inout):: this         !in: tensor cache
         class(tens_rcrsv_t), intent(in):: tensor          !in: tensor to look up (via its descriptor as the key)
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc,res
         class(*), pointer:: uptr
         type(dictionary_iter_t):: dit
         type(tens_descr_t), target:: tens_descr

         tens_entry_p=>NULL()
         call this%init_lock()
         tens_descr=tensor%get_descriptor(errc,only_signature=.TRUE.) !compute tensor descriptor by tensor signature only
         if(errc.eq.TEREC_SUCCESS) then
          call this%lock()
          errc=dit%init(this%map)
          if(errc.eq.GFC_SUCCESS) then
           uptr=>NULL()
           res=dit%search(GFC_DICT_FETCH_IF_FOUND,cmp_tens_descriptors,tens_descr,value_out=uptr)
           if(res.eq.GFC_FOUND) then
            select type(uptr)
            class is(tens_cache_entry_t)
             !call uptr%lock() !can cause a deadlock
             call uptr%incr_use_count()
             tens_entry_p=>uptr
            end select
            uptr=>NULL()
           else
            if(res.ne.GFC_NOT_FOUND) errc=-4
           endif
           res=dit%release(); if(res.ne.GFC_SUCCESS.and.errc.eq.0) errc=-3
          else
           errc=-2
          endif
          call this%unlock()
         else
          if(DEBUG.gt.0) then
           write(jo,'("#ERROR(VIRTA:tens_cache_t.lookup)[",i6,":",i3,"]: Unable to get tensor descriptor: Error ",i11)')&
           &impir,omp_get_thread_num(),errc
           flush(jo)
          endif
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensCacheLookup
!-----------------------------------------------------------------------------------------------------
        function TensCacheStore(this,tensor,tens_cache_entry_alloc_f,ierr,tens_entry_p) result(stored)
!Given a tensor, checks whether it is present in the tensor cache. If yes, returns
!a pointer to the corresponding tensor cache entry. If no, allocates a new extended
!tensor cache entry, stores it in the tensor cache and returns a pointer to it.
!The allocation is done via a user-provided non-member allocator of an extension of
!tens_cache_entry_t supplied with a proper dtor. Note that the newly created extended
!tensor cache entry may need further construction. Here only the <tensor> component of
!an extended tens_cache_entry_t is set.
         implicit none
         logical:: stored                                                !out: TRUE on successful new store, FALSE otherwise
         class(tens_cache_t), intent(inout):: this                       !inout: tensor cache
         class(tens_rcrsv_t), pointer, intent(in):: tensor               !in: pointer to a tensor to be stored
         procedure(tens_cache_entry_alloc_i):: tens_cache_entry_alloc_f  !in: non-member allocator of an extended(tens_cache_entry_t) class with a proper dtor
         integer(INTD), intent(out), optional:: ierr                     !out: error code
         class(tens_cache_entry_t), pointer, intent(out), optional:: tens_entry_p !out: tensor cache entry (newly created or existing)
         integer(INTD):: errc,res
         class(*), pointer:: uptr
         type(dictionary_iter_t):: dit
         type(tens_descr_t), target:: tens_descr
         class(tens_cache_entry_t), allocatable, target:: tce
         class(tens_cache_entry_t), pointer:: tcep

         stored=.FALSE.
         call this%init_lock()
         if(associated(tensor)) then
          if(DEBUG.gt.1) then
           write(jo,'("#MSG(TensorCache)[",i6,"]: Creating cache entry for tensor")') impir; flush(jo)
           call tensor%print_it(dev_id=jo); flush(jo)
          endif
          tens_descr=tensor%get_descriptor(errc,only_signature=.TRUE.) !compute tensor descriptor by tensor signature only
          if(errc.eq.TEREC_SUCCESS) then
           call this%lock()
           errc=dit%init(this%map)
           if(errc.eq.GFC_SUCCESS) then
            errc=tens_cache_entry_alloc_f(tce) !allocates an empty instance of an extended(tens_cache_entry_t) with a proper dtor
            if(errc.eq.0) then
             uptr=>NULL()
             res=dit%search(GFC_DICT_ADD_IF_NOT_FOUND,cmp_tens_descriptors,tens_descr,tce,GFC_BY_VAL,value_out=uptr)
             tcep=>NULL(); select type(uptr); class is(tens_cache_entry_t); tcep=>uptr; end select; uptr=>NULL()
             if(associated(tcep)) then
              if(res.eq.GFC_NOT_FOUND) then
               call tcep%init_lock()
               call tcep%set_tensor(tensor,errc); if(errc.ne.0) errc=-8
               if(DEBUG.gt.1.and.errc.eq.0) then; write(jo,'("Tensor cache entry created")'); flush(jo); endif
               stored=.TRUE.; this%num_entries=this%num_entries+1
              else
               if(res.eq.GFC_FOUND) then
                if(DEBUG.gt.1) then; write(jo,'("Tensor cache entry exists")'); flush(jo); endif
               else
                errc=-7
               endif
              endif
              if(errc.eq.0.and.present(tens_entry_p)) then
               !call tcep%lock() !can cause a deadlock
               call tcep%incr_use_count()
               tens_entry_p=>tcep
              endif
              tcep=>NULL()
             else
              errc=-6
             endif
             deallocate(tce) !a clone has been stored in this%map
            else
             errc=-5
            endif
            res=dit%release(); if(res.ne.GFC_SUCCESS.and.errc.eq.0) errc=-4
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
         if(errc.ne.0.and.present(tens_entry_p)) tens_entry_p=>NULL()
         if(present(ierr)) ierr=errc
         return
        end function TensCacheStore
!-------------------------------------------------------------------------
        function TensCacheEvict(this,tensor,ierr,decr_use) result(evicted)
!Evicts a specific tensor cache entry from the tensor cache.
!If the corresponding tensor cache entry is not found or still in use,
!no error is risen, but <evicted>=FALSE.
         implicit none
         logical:: evicted                                    !out: TRUE if the tensor cache entry has been found and evicted, FALSE otherwise
         class(tens_cache_t), intent(inout):: this            !inout: tensor cache
         class(tens_rcrsv_t), pointer, intent(inout):: tensor !inout: tensor to find and evict from the cache (will become NULL upon eviction)
         integer(INTD), intent(out), optional:: ierr          !out: error code
         logical, intent(in), optional:: decr_use             !in: if TRUE, the tensor cache entry USE counter will be decremented before eviction
         integer(INTD):: errc,res
         type(dictionary_iter_t):: dit
         type(tens_descr_t), target:: tens_descr
         class(*), pointer:: uptr
         logical:: decr,ready_to_die

         evicted=.FALSE.; decr=.FALSE.; if(present(decr_use)) decr=decr_use
         call this%init_lock()
         if(DEBUG.gt.1) then
          write(jo,'("#MSG(TensorCache)[",i6,"]: Evicting cache entry for tensor")') impir; flush(jo)
          call tensor%print_it(dev_id=jo); flush(jo)
         endif
         tens_descr=tensor%get_descriptor(errc,only_signature=.TRUE.) !compute tensor descriptor by tensor signature only
         if(errc.eq.TEREC_SUCCESS) then
          call this%lock()
          errc=dit%init(this%map)
          if(errc.eq.GFC_SUCCESS) then
 !Look up the tensor cache entry first:
           res=dit%search(GFC_DICT_FETCH_IF_FOUND,cmp_tens_descriptors,tens_descr,value_out=uptr)
           if(res.eq.GFC_FOUND) then
            select type(uptr)
            class is(tens_cache_entry_t)
             if(decr) call uptr%decr_use_count()
             ready_to_die=(uptr%get_ref_count().eq.0.and.uptr%get_use_count().eq.0.and.(.not.uptr%is_persistent()))
            class default
             errc=-6
            end select
            uptr=>NULL()
 !Delete the tensor cache entry if exists and not in use:
            if(errc.eq.0) then
             if(ready_to_die) then
              res=dit%search(GFC_DICT_DELETE_IF_FOUND,cmp_tens_descriptors,tens_descr)
              if(res.eq.GFC_FOUND) then
               if(DEBUG.gt.1) then; write(jo,'("Tensor cache entry evicted")'); flush(jo); endif
               tensor=>NULL(); evicted=.TRUE.; this%num_entries=this%num_entries-1
              else
               errc=-5
              endif
             else
              if(DEBUG.gt.1) then; write(jo,'("Tensor cache entry is still in use")'); flush(jo); endif
             endif
            endif
           else
            if(res.eq.GFC_NOT_FOUND) then
             if(DEBUG.gt.1) then; write(jo,'("Tensor cache entry does not exist")'); flush(jo); endif
            else
             errc=-4
            endif
           endif
           res=dit%release(); if(res.ne.GFC_SUCCESS.and.errc.eq.0) errc=-3
          else
           errc=-2
          endif
          call this%unlock()
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensCacheEvict
!---------------------------------------------------------------
        subroutine TensCacheReleaseEntry(this,tens_entry_p,ierr)
!Releases an explicitly previously obtained pointer to a tensor cache entry.
         implicit none
         class(tens_cache_t), intent(in):: this                           !in: tensor cache
         class(tens_cache_entry_t), pointer, intent(inout):: tens_entry_p !inout: pointer to a tensor cache entry (will sink here)
         integer(INTD), intent(out), optional:: ierr                      !out: error code
         integer(INTD):: errc

         errc=0
         if(associated(tens_entry_p)) then
          call tens_entry_p%decr_use_count()
          !call tens_entry_p%unlock() !can cause a deadlock
          tens_entry_p=>NULL()
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCacheReleaseEntry
!-------------------------------------------
        subroutine TensCacheErase(this,ierr)
!Erases the tensor cache completely.
         implicit none
         class(tens_cache_t), intent(inout):: this   !inout: tensor cache
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,res
         type(dictionary_iter_t):: dit

         call this%init_lock()
         call this%lock()
         errc=dit%init(this%map)
         if(errc.eq.GFC_SUCCESS) then
          errc=dit%delete_all(); if(errc.ne.GFC_SUCCESS) errc=-3
          res=dit%release(); if(res.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
          this%num_entries=0
         else
          errc=-1
         endif
         call this%unlock()
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCacheErase
!--------------------------------------------------------
        function TensCacheGetNumEntries(this) result(num)
!Returns the current number of active entries in the tensor cache.
         implicit none
         integer(INTL):: num                       !out: number of active cache entries
         class(tens_cache_t), intent(inout):: this !in: tensor cache

         call this%init_lock()
         call this%lock()
         num=this%num_entries
         call this%unlock()
         return
        end function TensCacheGetNumEntries
!---------------------------------------------
        subroutine TensCachePrintIt(this,ierr)
!Prints the content of the tensor cache.
         implicit none
         class(tens_cache_t), intent(inout):: this   !in: tensor cache
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         type(dictionary_iter_t):: dit

         call this%init_lock()
         call this%lock()
         errc=dit%init(this%map)
         if(errc.eq.GFC_SUCCESS) then
          errc=dit%scanp(action_f=tens_cache_entry_print_f)
          errc=dit%release()
         endif
         call this%unlock()
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCachePrintIt
!---------------------------------------
        subroutine tens_cache_dtor(this)
         implicit none
         type(tens_cache_t):: this !inout: tensor cache

         call this%erase()
#ifndef NO_OMP
!$OMP CRITICAL (TAVP_CACHE)
         if(allocated(this%cache_lock)) then
          call omp_destroy_lock(this%cache_lock)
          deallocate(this%cache_lock)
         endif
!$OMP END CRITICAL (TAVP_CACHE)
#endif
         return
        end subroutine tens_cache_dtor
![data_register_t]---------------------------------------------------------
        subroutine DataRegisterRegisterData(this,data_name,extrn_data,ierr)
         implicit none
         class(data_register_t), intent(inout):: this     !inout: data register
         character(*), intent(in):: data_name             !in: data name
         type(talsh_tens_data_t), intent(in):: extrn_data !in: external data
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD):: errc,ier
         type(dictionary_iter_t):: dit

!$OMP CRITICAL (TAVP_DATA_REG)
         errc=dit%init(this%ext_data)
         if(errc.eq.GFC_SUCCESS) then
          errc=dit%search(GFC_DICT_ADD_IF_NOT_FOUND,cmp_strings,data_name,extrn_data)
          if(errc.eq.GFC_NOT_FOUND) then
           errc=0
          else
           if(errc.eq.GFC_FOUND) then; errc=-3; else; errc=-4; endif
          endif
          ier=dit%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
!$OMP END CRITICAL (TAVP_DATA_REG)
         if(present(ierr)) ierr=errc
         return
        end subroutine DataRegisterRegisterData
!-----------------------------------------------------------------
        subroutine DataRegisterUnregisterData(this,data_name,ierr)
         implicit none
         class(data_register_t), intent(inout):: this !inout: data register
         character(*), intent(in):: data_name         !in: data name
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc,ier
         type(dictionary_iter_t):: dit

!$OMP CRITICAL (TAVP_DATA_REG)
         errc=dit%init(this%ext_data)
         if(errc.eq.GFC_SUCCESS) then
          errc=dit%search(GFC_DICT_DELETE_IF_FOUND,cmp_strings,data_name)
          if(errc.eq.GFC_FOUND) then
           errc=0
          else
           if(errc.eq.GFC_NOT_FOUND) then; errc=-3; else; errc=-4; endif
          endif
          ier=dit%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
!$OMP END CRITICAL (TAVP_DATA_REG)
         if(present(ierr)) ierr=errc
         return
        end subroutine DataRegisterUnregisterData
!--------------------------------------------------------------------------------
        function DataRegisterRetrieveData(this,data_name,ierr) result(extrn_data)
         implicit none
         type(talsh_tens_data_t), pointer:: extrn_data !out: pointer to the data
         class(data_register_t), intent(in):: this     !inout: data register
         character(*), intent(in):: data_name          !in: data name
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: errc,ier
         type(dictionary_iter_t):: dit
         class(*), pointer:: uptr

         extrn_data=>NULL()
!$OMP CRITICAL (TAVP_DATA_REG)
         errc=dit%init(this%ext_data)
         if(errc.eq.GFC_SUCCESS) then
          uptr=>NULL()
          errc=dit%search(GFC_DICT_FETCH_IF_FOUND,cmp_strings,data_name,value_out=uptr)
          if(errc.eq.GFC_FOUND) then
           if(associated(uptr)) then
            select type(uptr); type is(talsh_tens_data_t); extrn_data=>uptr; end select
            if(.not.associated(extrn_data)) errc=-6
           else
            errc=-5
           endif
          else
           if(errc.eq.GFC_NOT_FOUND) then; errc=-3; else; errc=-4; endif
          endif
          ier=dit%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
!$OMP END CRITICAL (TAVP_DATA_REG)
         if(present(ierr)) ierr=errc
         return
        end function DataRegisterRetrieveData
!------------------------------------------
        subroutine data_register_dtor(this)
         implicit none
         type(data_register_t):: this  !inout: data register
         integer(INTD):: errc
         type(dictionary_iter_t):: dit

!$OMP CRITICAL (TAVP_DATA_REG)
         errc=dit%init(this%ext_data)
         if(errc.eq.GFC_SUCCESS) then
          errc=dit%delete_all()
          errc=dit%release()
         endif
!$OMP END CRITICAL (TAVP_DATA_REG)
         return
        end subroutine data_register_dtor
![method_register_t]---------------------------------------------------------------
        subroutine MethodRegisterRegisterMethod(this,method_name,extrn_method,ierr)
         implicit none
         class(method_register_t), intent(inout):: this      !inout: method register
         character(*), intent(in):: method_name              !in: method name
         class(tens_method_uni_t), intent(in):: extrn_method !in: external unary tensor method (initialization/transformation)
         integer(INTD), intent(out), optional:: ierr         !out: error code
         integer(INTD):: errc,ier
         type(dictionary_iter_t):: dit

!$OMP CRITICAL (TAVP_METHOD_REG)
         errc=dit%init(this%ext_methods)
         if(errc.eq.GFC_SUCCESS) then
          errc=dit%search(GFC_DICT_ADD_IF_NOT_FOUND,cmp_strings,method_name,extrn_method)
          if(errc.eq.GFC_NOT_FOUND) then
           errc=0
          else
           if(errc.eq.GFC_FOUND) then; errc=-3; else; errc=-4; endif
          endif
          ier=dit%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
!$OMP END CRITICAL (TAVP_METHOD_REG)
         if(present(ierr)) ierr=errc
         return
        end subroutine MethodRegisterRegisterMethod
!-----------------------------------------------------------------------
        subroutine MethodRegisterUnregisterMethod(this,method_name,ierr)
         implicit none
         class(method_register_t), intent(inout):: this !inout: method register
         character(*), intent(in):: method_name         !in: method name
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD):: errc,ier
         type(dictionary_iter_t):: dit

!$OMP CRITICAL (TAVP_METHOD_REG)
         errc=dit%init(this%ext_methods)
         if(errc.eq.GFC_SUCCESS) then
          errc=dit%search(GFC_DICT_DELETE_IF_FOUND,cmp_strings,method_name)
          if(errc.eq.GFC_FOUND) then
           errc=0
          else
           if(errc.eq.GFC_NOT_FOUND) then; errc=-3; else; errc=-4; endif
          endif
          ier=dit%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
!$OMP END CRITICAL (TAVP_METHOD_REG)
         if(present(ierr)) ierr=errc
         return
        end subroutine MethodRegisterUnregisterMethod
!----------------------------------------------------------------------------------------
        function MethodRegisterRetrieveMethod(this,method_name,ierr) result(extrn_method)
         implicit none
         class(tens_method_uni_t), pointer:: extrn_method !out: pointer to the method
         class(method_register_t), intent(in):: this      !inout: method register
         character(*), intent(in):: method_name           !in: method name
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD):: errc,ier
         type(dictionary_iter_t):: dit
         class(*), pointer:: uptr

         extrn_method=>NULL()
!$OMP CRITICAL (TAVP_METHOD_REG)
         errc=dit%init(this%ext_methods)
         if(errc.eq.GFC_SUCCESS) then
          uptr=>NULL()
          errc=dit%search(GFC_DICT_FETCH_IF_FOUND,cmp_strings,method_name,value_out=uptr)
          if(errc.eq.GFC_FOUND) then
           if(associated(uptr)) then
            select type(uptr); class is(tens_method_uni_t); extrn_method=>uptr; end select
            if(.not.associated(extrn_method)) errc=-6
           else
            errc=-5
           endif
          else
           if(errc.eq.GFC_NOT_FOUND) then; errc=-3; else; errc=-4; endif
          endif
          ier=dit%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
!$OMP END CRITICAL (TAVP_METHOD_REG)
         if(present(ierr)) ierr=errc
         return
        end function MethodRegisterRetrieveMethod
!--------------------------------------------
        subroutine method_register_dtor(this)
         implicit none
         type(method_register_t):: this  !inout: method register
         integer(INTD):: errc
         type(dictionary_iter_t):: dit

!$OMP CRITICAL (TAVP_METHOD_REG)
         errc=dit%init(this%ext_methods)
         if(errc.eq.GFC_SUCCESS) then
          errc=dit%delete_all()
          errc=dit%release()
         endif
!$OMP END CRITICAL (TAVP_METHOD_REG)
         return
        end subroutine method_register_dtor
!-------------------------------------------------------------
        function method_map_f(method_name,ierr) result(method)
!Non-member function mapping method names to method objects.
         implicit none
         class(tens_method_uni_t), pointer:: method  !out: pointer to a registered external unary tensor method
         character(*), intent(in):: method_name      !in: method name
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,l

         method=>NULL(); l=len_trim(method_name)
         if(l.gt.0) then
          method=>method_register%retrieve_method(method_name(1:l),errc); if(errc.ne.0) errc=-1
         else
          errc=-2
         endif
         if(present(ierr)) ierr=errc
         return
        end function method_map_f

       end module virta
