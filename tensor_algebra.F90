!ExaTensor::TAL-SH: Basic parameters and types:
!Keep consistent with "tensor_algebra.h"!
!REVISION: 2015/11/11
        module tensor_algebra
        use dil_kinds !contains ISO_C_BINDING
        implicit none
        public
!TENSOR ALGEBRA LIMITS (keep consistent with tensor_algebra.h):
        integer(C_INT), parameter, public:: MAX_TENSOR_RANK=32    !max allowed tensor rank (max number of indices in a tensor)
        integer(C_INT), parameter, public:: MAX_TENSOR_OPERANDS=4 !max number of tensor operands in a tensor operation
#ifndef NO_PHI
!DIR$ ATTRIBUTES OFFLOAD:mic:: MAX_TENSOR_RANK,MAX_TENSOR_OPERANDS
!DIR$ ATTRIBUTES ALIGN:128:: MAX_TENSOR_RANK,MAX_TENSOR_OPERANDS
#endif

!DEVICE KINDS (keep consistent with tensor_algebra.h):
        integer(C_INT), parameter, public:: MAX_GPUS_PER_NODE=8   !max number of NVidia GPUs on a node
        integer(C_INT), parameter, public:: MAX_MICS_PER_NODE=8   !max number of Intel MICs on a node
        integer(C_INT), parameter, public:: MAX_AMDS_PER_NODE=8   !max number of AMD APUs on a node
        integer(C_INT), parameter, public:: DEV_NULL=-1           !abstract null device
        integer(C_INT), parameter, public:: DEV_HOST=0            !multicore CPU Host (includes all self-hosted systems)
        integer(C_INT), parameter, public:: DEV_NVIDIA_GPU=1      !NVidia GPU
        integer(C_INT), parameter, public:: DEV_INTEL_MIC=2       !Intel Xeon Phi
        integer(C_INT), parameter, public:: DEV_AMD_GPU=3         !AMD APU
        integer(C_INT), parameter, public:: DEV_MAX=1+MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE+MAX_AMDS_PER_NODE

!TENSOR DATA KINDS (keep consistent with tensor_algebra.h):
        integer(C_INT), parameter, public:: NO_TYPE=0 !no type/kind
        integer(C_INT), parameter, public:: R4=4      !float data kind
        integer(C_INT), parameter, public:: R8=8      !double data kind
        integer(C_INT), parameter, public:: C8=16     !double complex data kind
        real(4), parameter, public:: R4_=0.0
        real(8), parameter, public:: R8_=0d0
        complex(8), parameter, public:: C8_=(0d0,0d0)
#ifndef NO_PHI
!DIR$ ATTRIBUTES OFFLOAD:mic:: NO_TYPE,R4,R8,C8,R4_,R8_,C8_
!DIR$ ATTRIBUTES ALIGN:128:: NO_TYPE,R4,R8,C8,R4_,R8_,C8_
#endif

!ALIASES (keep consistent with tensor_algebra.h):
        integer(C_INT), parameter, public:: TALSH_SUCCESS=0             !success
        integer(C_INT), parameter, public:: TALSH_FAILURE=-666          !failure
        integer(C_INT), parameter, public:: BLAS_ON=0                   !enables BLAS
        integer(C_INT), parameter, public:: BLAS_OFF=1                  !disables BLAS
        integer(C_INT), parameter, public:: EFF_TRN_OFF=0               !disables efficient tensor transpose algorithm
        integer(C_INT), parameter, public:: EFF_TRN_ON=1                !enables efficient tensor transpose algorithm
        integer(C_INT), parameter, public:: TRY_LATER=-918273645        !try the action later (resources are currently busy): KEEP THIS UNIQUE!
        integer(C_INT), parameter, public:: DEVICE_UNABLE=-546372819    !device is unsuitable for the given task: KEEP THIS UNIQUE!
        integer(C_INT), parameter, public:: NOT_CLEAN=-192837465        !something, like resource release, did not go right, but you can continue: KEEP THIS UNIQUE!
        integer(C_INT), parameter, public:: NOPE=0                      !"NO" answer
        integer(C_INT), parameter, public:: YEP=1                       !"YES" answer
        integer(C_INT), parameter, public:: EVERYTHING=0                !everything (source, destination, temporary)
        integer(C_INT), parameter, public:: SOURCE=1                    !source
        integer(C_INT), parameter, public:: DESTINATION=2               !destination
        integer(C_INT), parameter, public:: TEMPORARY=3                 !temporary
        integer(C_INT), parameter, public:: DEV_OFF=0                   !device status "Disabled"
        integer(C_INT), parameter, public:: DEV_ON=1                    !device status "Enabled"
        integer(C_INT), parameter, public:: DEV_ON_BLAS=2               !device status "Enabled with vendor provided BLAS"
        integer(C_INT), parameter, public:: NO_COPY_BACK=0              !keeps the tensor-result on Accelerator without updating Host
        integer(C_INT), parameter, public:: COPY_BACK=1                 !tensor-result will be copied back from Accelerator to Host (default)

        integer(C_INT), parameter, public:: COPY_D=0
        integer(C_INT), parameter, public:: COPY_M=1
        integer(C_INT), parameter, public:: COPY_T=2
        integer(C_INT), parameter, public:: COPY_K=3
        integer(C_INT), parameter, public:: COPY_DD=0
        integer(C_INT), parameter, public:: COPY_DM=1
        integer(C_INT), parameter, public:: COPY_DT=2
        integer(C_INT), parameter, public:: COPY_DK=3
        integer(C_INT), parameter, public:: COPY_MD=4
        integer(C_INT), parameter, public:: COPY_MM=5
        integer(C_INT), parameter, public:: COPY_MT=6
        integer(C_INT), parameter, public:: COPY_MK=7
        integer(C_INT), parameter, public:: COPY_TD=8
        integer(C_INT), parameter, public:: COPY_TM=9
        integer(C_INT), parameter, public:: COPY_TT=10
        integer(C_INT), parameter, public:: COPY_TK=11
        integer(C_INT), parameter, public:: COPY_KD=12
        integer(C_INT), parameter, public:: COPY_KM=13
        integer(C_INT), parameter, public:: COPY_KT=14
        integer(C_INT), parameter, public:: COPY_KK=15
        integer(C_INT), parameter, public:: COPY_DDD=0
        integer(C_INT), parameter, public:: COPY_DDM=1
        integer(C_INT), parameter, public:: COPY_DDT=2
        integer(C_INT), parameter, public:: COPY_DDK=3
        integer(C_INT), parameter, public:: COPY_DMD=4
        integer(C_INT), parameter, public:: COPY_DMM=5
        integer(C_INT), parameter, public:: COPY_DMT=6
        integer(C_INT), parameter, public:: COPY_DMK=7
        integer(C_INT), parameter, public:: COPY_DTD=8
        integer(C_INT), parameter, public:: COPY_DTM=9
        integer(C_INT), parameter, public:: COPY_DTT=10
        integer(C_INT), parameter, public:: COPY_DTK=11
        integer(C_INT), parameter, public:: COPY_DKD=12
        integer(C_INT), parameter, public:: COPY_DKM=13
        integer(C_INT), parameter, public:: COPY_DKT=14
        integer(C_INT), parameter, public:: COPY_DKK=15
        integer(C_INT), parameter, public:: COPY_MDD=16
        integer(C_INT), parameter, public:: COPY_MDM=17
        integer(C_INT), parameter, public:: COPY_MDT=18
        integer(C_INT), parameter, public:: COPY_MDK=19
        integer(C_INT), parameter, public:: COPY_MMD=20
        integer(C_INT), parameter, public:: COPY_MMM=21
        integer(C_INT), parameter, public:: COPY_MMT=22
        integer(C_INT), parameter, public:: COPY_MMK=23
        integer(C_INT), parameter, public:: COPY_MTD=24
        integer(C_INT), parameter, public:: COPY_MTM=25
        integer(C_INT), parameter, public:: COPY_MTT=26
        integer(C_INT), parameter, public:: COPY_MTK=27
        integer(C_INT), parameter, public:: COPY_MKD=28
        integer(C_INT), parameter, public:: COPY_MKM=29
        integer(C_INT), parameter, public:: COPY_MKT=30
        integer(C_INT), parameter, public:: COPY_MKK=31
        integer(C_INT), parameter, public:: COPY_TDD=32
        integer(C_INT), parameter, public:: COPY_TDM=33
        integer(C_INT), parameter, public:: COPY_TDT=34
        integer(C_INT), parameter, public:: COPY_TDK=35
        integer(C_INT), parameter, public:: COPY_TMD=36
        integer(C_INT), parameter, public:: COPY_TMM=37
        integer(C_INT), parameter, public:: COPY_TMT=38
        integer(C_INT), parameter, public:: COPY_TMK=39
        integer(C_INT), parameter, public:: COPY_TTD=40
        integer(C_INT), parameter, public:: COPY_TTM=41
        integer(C_INT), parameter, public:: COPY_TTT=42
        integer(C_INT), parameter, public:: COPY_TTK=43
        integer(C_INT), parameter, public:: COPY_TKD=44
        integer(C_INT), parameter, public:: COPY_TKM=45
        integer(C_INT), parameter, public:: COPY_TKT=46
        integer(C_INT), parameter, public:: COPY_TKK=47
        integer(C_INT), parameter, public:: COPY_KDD=48
        integer(C_INT), parameter, public:: COPY_KDM=49
        integer(C_INT), parameter, public:: COPY_KDT=50
        integer(C_INT), parameter, public:: COPY_KDK=51
        integer(C_INT), parameter, public:: COPY_KMD=52
        integer(C_INT), parameter, public:: COPY_KMM=53
        integer(C_INT), parameter, public:: COPY_KMT=54
        integer(C_INT), parameter, public:: COPY_KMK=55
        integer(C_INT), parameter, public:: COPY_KTD=56
        integer(C_INT), parameter, public:: COPY_KTM=57
        integer(C_INT), parameter, public:: COPY_KTT=58
        integer(C_INT), parameter, public:: COPY_KTK=59
        integer(C_INT), parameter, public:: COPY_KKD=60
        integer(C_INT), parameter, public:: COPY_KKM=61
        integer(C_INT), parameter, public:: COPY_KKT=62
        integer(C_INT), parameter, public:: COPY_KKK=63
#ifndef NO_PHI
!DIR$ ATTRIBUTES OFFLOAD:mic:: TALSH_SUCCESS,TALSH_FAILURE,BLAS_ON,BLAS_OFF,EFF_TRN_OFF,EFF_TRN_ON,TRY_LATER,DEVICE_UNABLE
!DIR$ ATTRIBUTES OFFLOAD:mic:: NOPE,YEP,EVERYTHING,SOURCE,DESTINATION,TEMPORARY,DEV_OFF,DEV_ON,DEV_ON_BLAS,NO_COPY_BACK,COPY_BACK
!DIR$ ATTRIBUTES ALIGN:128:: TALSH_SUCCESS,TALSH_FAILURE,BLAS_ON,BLAS_OFF,EFF_TRN_OFF,EFF_TRN_ON,TRY_LATER,DEVICE_UNABLE
!DIR$ ATTRIBUTES ALIGN:128:: NOPE,YEP,EVERYTHING,SOURCE,DESTINATION,TEMPORARY,DEV_OFF,DEV_ON,DEV_ON_BLAS,NO_COPY_BACK,COPY_BACK
#endif

!CUDA TASK STATUS (keep consistent with tensor_algebra.h):
        integer(C_INT), parameter, public:: CUDA_TASK_ERROR=-1
        integer(C_INT), parameter, public:: CUDA_TASK_EMPTY=0
        integer(C_INT), parameter, public:: CUDA_TASK_SCHEDULED=1
        integer(C_INT), parameter, public:: CUDA_TASK_STARTED=2
        integer(C_INT), parameter, public:: CUDA_TASK_INPUT_THERE=3
        integer(C_INT), parameter, public:: CUDA_TASK_OUTPUT_THERE=4
        integer(C_INT), parameter, public:: CUDA_TASK_COMPLETED=5

!TENSOR BLOCK STORAGE LAYOUT:
        integer(C_INT), parameter, public:: NOT_ALLOCATED=0   !tensor block has not been allocated/initialized
        integer(C_INT), parameter, public:: SCALAR_TENSOR=1   !scalar (rank-0 tensor)
        integer(C_INT), parameter, public:: DIMENSION_LED=2   !dense tensor block (column-major storage by default): no symmetry restrictions
        integer(C_INT), parameter, public:: BRICKED_DENSE=3   !dense tensor block (bricked storage): no symmetry restrictions
        integer(C_INT), parameter, public:: BRICKED_ORDERED=4 !symmetrically packed tensor block (bricked storage): symmetry restrictions apply
        integer(C_INT), parameter, public:: SPARSE_LIST=5     !sparse tensor block: symmetry restrictions do not apply!
        integer(C_INT), parameter, public:: COMPRESSED=6      !compressed tensor block: symmetry restrictions do not apply!
        logical, parameter, public:: FORTRAN_LIKE=.true.
        logical, parameter, public:: C_LIKE=.false.
#ifndef NO_PHI
!DIR$ ATTRIBUTES OFFLOAD:mic:: NOT_ALLOCATED,SCALAR_TENSOR,DIMENSION_LED,BRICKED_DENSE,BRICKED_ORDERED,SPARSE_LIST,COMPRESSED
!DIR$ ATTRIBUTES OFFLOAD:mic:: FORTRAN_LIKE,C_LIKE
!DIR$ ATTRIBUTES ALIGN:128:: NOT_ALLOCATED,SCALAR_TENSOR,DIMENSION_LED,BRICKED_DENSE,BRICKED_ORDERED,SPARSE_LIST,COMPRESSED
!DIR$ ATTRIBUTES ALIGN:128:: FORTRAN_LIKE,C_LIKE
#endif

!INTEROPERABLE TYPES (keep consistent with tensor_algebra.h):
 !TAL-SH tensor block:
        type, bind(C):: talsh_tens_t
         integer(C_INT):: ndev       !number of devices the tensor block resides on
         integer(C_INT):: last_write !flat device id where the last write happened, -1 means coherence on all devices where the tensor block resides
         type(C_PTR):: dev_rsc       !list of device resources occupied by the tensor block on each device
         type(C_PTR):: tensF         !pointer to Fortran <tensor_block_t>
         type(C_PTR):: tensC         !pointer to C/C++ <tensBlck_t>
        end type talsh_tens_t
 !TAL-SH task handle:
        type, bind(C):: talsh_task_t
         integer(C_INT):: dev_kind   !device kind
         type(C_PTR):: task_p        !pointer to the corresponding task object
        end type talsh_task_t

!EXTERNAL INTERFACES (keep consistent with tensor_algebra.h):
 !User-defined tensor block initialization:
        abstract interface
         subroutine talsh_tens_init_i(tens_ptr,data_type,tens_rank,tens_dims,ierr)
          import
          type(C_PTR), value:: tens_ptr                !in: pointer to the tensor elements storage
          integer(C_INT), value:: data_type            !in: data type: {R4,R8,C8}
          integer(C_INT), value:: tens_rank            !in: tensor block rank
          integer(C_INT), intent(in):: tens_dims(1:*)  !in: tensor block dimension extents
          integer(C_INT), intent(out):: ierr           !out: error code (0:success)
         end subroutine talsh_tens_init_i
        end interface
!C FUNCTION INTERFACES (for Fortran):
        interface
 !Argument buffer memory management:
  !Initialize all argument buffers on Host and all Devices (GPU constant+global, MICs, etc):
         integer(C_INT) function arg_buf_allocate(host_mem,arg_max,gpu_beg,gpu_end) bind(c,name='arg_buf_allocate')
          import
          implicit none
          integer(C_SIZE_T), intent(inout):: host_mem
          integer(C_INT), intent(out):: arg_max
          integer(C_INT), value, intent(in):: gpu_beg
          integer(C_INT), value, intent(in):: gpu_end
         end function arg_buf_allocate
  !Deallocate all argument buffers on Host and all Devices:
         integer(C_INT) function arg_buf_deallocate(gpu_beg,gpu_end) bind(c,name='arg_buf_deallocate')
          import
          implicit none
          integer(C_INT), value, intent(in):: gpu_beg
          integer(C_INT), value, intent(in):: gpu_end
         end function arg_buf_deallocate
  !Check whether the Host argument buffer is clean:
         integer(C_INT) function arg_buf_clean_host() bind(c,name='arg_buf_clean_host')
          import
          implicit none
         end function arg_buf_clean_host
#ifndef NO_GPU
  !Check whether a GPU argument buffer is clean:
         integer(C_INT) function arg_buf_clean_gpu(gpu_num) bind(c,name='arg_buf_clean_gpu')
          import
          implicit none
          integer(C_INT), value, intent(in):: gpu_num
         end function arg_buf_clean_gpu
#endif
  !Get the buffer block sizes for each level of the Host argument buffer:
         integer(C_INT) function get_blck_buf_sizes_host(blck_sizes) bind(c,name='get_blck_buf_sizes_host')
          import
          implicit none
          integer(C_SIZE_T), intent(out):: blck_sizes(*)
         end function get_blck_buf_sizes_host
#ifndef NO_GPU
  !Get the buffer block sizes for each level of the GPU argument buffer:
         integer(C_INT) function get_blck_buf_sizes_gpu(gpu_num,blck_sizes) bind(c,name='get_blck_buf_sizes_gpu')
          import
          implicit none
          integer(C_INT), value, intent(in):: gpu_num
          integer(C_SIZE_T), intent(out):: blck_sizes(*)
         end function get_blck_buf_sizes_gpu
#endif
  !Get a free argument entry in the Host argument buffer:
         integer(C_INT) function get_buf_entry_host(bsize,entry_ptr,entry_num) bind(c,name='get_buf_entry_host')
          import
          implicit none
          integer(C_SIZE_T), value, intent(in):: bsize
          type(C_PTR), intent(out):: entry_ptr
          integer(C_INT), intent(out):: entry_num
         end function get_buf_entry_host
  !Free an argument entry in the Host argument buffer:
         integer(C_INT) function free_buf_entry_host(entry_num) bind(c,name='free_buf_entry_host')
          import
          implicit none
          integer(C_INT), value, intent(in):: entry_num
         end function free_buf_entry_host
#ifndef NO_GPU
  !Get a free argument entry in the GPU argument buffer:
         integer(C_INT) function get_buf_entry_gpu(gpu_num,bsize,entry_ptr,entry_num) bind(c,name='get_buf_entry_gpu')
          import
          implicit none
          integer(C_INT), value, intent(in):: gpu_num
          integer(C_SIZE_T), value, intent(in):: bsize
          type(C_PTR), intent(out):: entry_ptr
          integer(C_INT), intent(out):: entry_num
         end function get_buf_entry_gpu
  !Free an argument entry in the GPU argument buffer:
         integer(C_INT) function free_buf_entry_gpu(gpu_num,entry_num) bind(c,name='free_buf_entry_gpu')
          import
          implicit none
          integer(C_INT), value, intent(in):: gpu_num
          integer(C_INT), value, intent(in):: entry_num
         end function free_buf_entry_gpu
  !Get a free entry from the GPU constant memory argument bank:
         integer(C_INT) function const_args_entry_get(gpu_num,entry_num) bind(c,name='const_args_entry_get')
          import
          implicit none
          integer(C_INT), value, intent(in):: gpu_num
          integer(C_INT), intent(out):: entry_num
         end function const_args_entry_get
  !Return back an entry to the GPU constant memory argument bank:
         integer(C_INT) function const_args_entry_free(gpu_num,entry_num) bind(c,name='const_args_entry_free')
          import
          implicit none
          integer(C_INT), value, intent(in):: gpu_num
          integer(C_INT), value, intent(in):: entry_num
         end function const_args_entry_free
#endif
  !Query the free buffer space in bytes on a given device:
         integer(C_INT) function mem_free_left(dev_id,free_mem) bind(c,name='mem_free_left')
          import
          implicit none
          integer(C_INT), value, intent(in):: dev_id
          integer(C_SIZE_T), intent(out):: free_mem
         end function mem_free_left
  !Print the current status of the argument buffer on a given device:
         integer(C_INT) function mem_print_stats(dev_id) bind(c,name='mem_print_stats')
          import
          implicit none
          integer(C_INT), value, intent(in):: dev_id
         end function mem_print_stats
#ifndef NO_GPU
  !Allocate pinned memory on Host:
         integer(C_INT) function host_mem_alloc_pin(cptr,bsize) bind(c,name='host_mem_alloc_pin')
          import
          implicit none
          type(C_PTR), intent(out):: cptr
          integer(C_SIZE_T), value, intent(in):: bsize !bytes
         end function host_mem_alloc_pin
  !Free pinned memory on Host:
         integer(C_INT) function host_mem_free_pin(cptr) bind(c,name='host_mem_free_pin')
          import
          implicit none
          type(C_PTR), value:: cptr
         end function host_mem_free_pin
  !Register (pin) Host memory:
         integer(C_INT) function host_mem_register(cptr,bsize) bind(c,name='host_mem_register')
          import
          implicit none
          type(C_PTR), value:: cptr
          integer(C_SIZE_T), value, intent(in):: bsize !bytes
         end function host_mem_register
  !Unregister (unpin) Host memory:
         integer(C_INT) function host_mem_unregister(cptr) bind(c,name='host_mem_unregister')
          import
          implicit none
          type(C_PTR), value:: cptr
         end function host_mem_unregister
  !Allocate memory on current GPU:
         integer(C_INT) function gpu_mem_alloc(cptr,bsize) bind(c,name='gpu_mem_alloc')
          import
          implicit none
          type(C_PTR), intent(out):: cptr
          integer(C_SIZE_T), value, intent(in):: bsize !bytes
         end function gpu_mem_alloc
  !Free memory on current GPU:
         integer(C_INT) function gpu_mem_free(cptr) bind(c,name='gpu_mem_free')
          import
          implicit none
          type(C_PTR), value:: cptr
         end function gpu_mem_free
 !NV-TAL auxiliary (`TAL-SH level):
  !Convert a kind specific device id into the flat device id:
         integer(C_INT) function encode_device_id(dev_kind,dev_num) bind(c,name='encode_device_id')
          import
          implicit none
          integer(C_INT), value, intent(in):: dev_kind
          integer(C_INT), value, intent(in):: dev_num
         end function encode_device_id
  !Convert a flat device id into the kind specific device id:
         integer(C_INT) function decode_device_id(dev_id,dev_kind) bind(c,name='decode_device_id')
          import
          implicit none
          integer(C_INT), value, intent(in):: dev_id
          integer(C_INT), intent(out):: dev_kind
         end function decode_device_id
 !NV-TAL debugging:
  !Get the current GPU error count:
         integer(C_INT) function gpu_get_error_count() bind(c,name='gpu_get_error_count')
          import
          implicit none
         end function gpu_get_error_count
  !Get the current GPU debug dump:
         integer(C_INT) function gpu_get_debug_dump(dump) bind(c,name='gpu_get_debug_dump')
          import
          implicit none
          integer(C_INT), intent(out):: dump(*)
         end function gpu_get_debug_dump
 !NV-TAL internal control:
  !Set the width of the NVidia GPU shared memory bank:
         integer(C_INT) function gpu_set_shmem_width(width) bind(c,name='gpu_set_shmem_width')
          import
          implicit none
          integer(C_INT), value, intent(in):: width
         end function gpu_set_shmem_width
  !Set the tensor transpose algorithm:
         subroutine gpu_set_transpose_algorithm(alg) bind(c,name='gpu_set_transpose_algorithm')
          import
          implicit none
          integer(C_INT), value:: alg
         end subroutine gpu_set_transpose_algorithm
  !Set the matrix multiplication algorithm:
         subroutine gpu_set_matmult_algorithm(alg) bind(c,name='gpu_set_matmult_algorithm')
          import
          implicit none
          integer(C_INT), value:: alg
         end subroutine gpu_set_matmult_algorithm
 !NV-TAL tensor block API:
  !Create an instance of tensBlck_t:
         integer(C_INT) function tensBlck_create(ctens) bind(c,name='tensBlck_create')
          import
          implicit none
          type(C_PTR), intent(out):: ctens
         end function tensBlck_create
  !Destroy an instance of tensBlck_t:
         integer(C_INT) function tensBlck_destroy(ctens) bind(c,name='tensBlck_destroy')
          import
          implicit none
          type(C_PTR), value:: ctens
         end function tensBlck_destroy
  !Construct tensBlck_t using externally provided memory pointers (custom allocation):
         integer(C_INT) function tensBlck_construct(ctens,dev_kind,dev_num,data_kind,trank,&
                                  &addr_dims,addr_divs,addr_grps,addr_prmn,addr_host,addr_gpu,&
                                  &entry_host,entry_gpu,entry_const) bind(c,name='tensBlck_construct')
          import
          implicit none
          type(C_PTR), value:: ctens
          integer(C_INT), value, intent(in):: dev_kind
          integer(C_INT), value, intent(in):: dev_num
          integer(C_INT), value, intent(in):: data_kind
          integer(C_INT), value, intent(in):: trank
          type(C_PTR), value, intent(in):: addr_dims
          type(C_PTR), value, intent(in):: addr_divs
          type(C_PTR), value, intent(in):: addr_grps
          type(C_PTR), value, intent(in):: addr_prmn
          type(C_PTR), value, intent(in):: addr_host
          type(C_PTR), value, intent(in):: addr_gpu
          integer(C_INT), value, intent(in):: entry_host
          integer(C_INT), value, intent(in):: entry_gpu
          integer(C_INT), value, intent(in):: entry_const
         end function tensBlck_construct
  !Allocate a space for tensBlck_t data in Host and GPU memory:
         integer(C_INT) function tensBlck_alloc(ctens,dev_num,data_kind,trank,dims) bind(c,name='tensBlck_alloc')
          import
          implicit none
          type(C_PTR), value:: ctens
          integer(C_INT), value, intent(in):: dev_num
          integer(C_INT), value, intent(in):: data_kind
          integer(C_INT), value, intent(in):: trank
          type(C_PTR), value, intent(in):: dims
         end function tensBlck_alloc
  !Free space occupied by the tensBlck_t data:
         integer(C_INT) function tensBlck_free(ctens) bind(c,name='tensBlck_free')
          import
          implicit none
          type(C_PTR), value:: ctens
         end function tensBlck_free
  !Get the accelerator ID and other information from tensBlck_t:
         integer(C_INT) function tensBlck_acc_id(ctens,dev_kind,entry_gpu,entry_const,data_kind,there)&
                                  &bind(c,name='tensBlck_acc_id')
          import
          implicit none
          type(C_PTR), value, intent(in):: ctens
          integer(C_INT), intent(out):: dev_kind
          integer(C_INT), intent(out):: entry_gpu
          integer(C_INT), intent(out):: entry_const
          integer(C_INT), intent(out):: data_kind
          integer(C_INT), intent(out):: there
         end function tensBlck_acc_id
  !Mark tensBlck_t as residing on GPU:
         integer(C_INT) function tensBlck_set_presence(ctens) bind(c,name='tensBlck_set_presence')
          import
          implicit none
          type(C_PTR), value:: ctens
         end function tensBlck_set_presence
  !Unmark tensBlck_t as residing on GPU:
         integer(C_INT) function tensBlck_set_absence(ctens) bind(c,name='tensBlck_set_absence')
          import
          implicit none
          type(C_PTR), value:: ctens
         end function tensBlck_set_absence
  !Check presence (residence) of tensBlck_t data on GPU:
         integer(C_INT) function tensBlck_present(ctens) bind(c,name='tensBlck_present')
          import
          implicit none
          type(C_PTR), value, intent(in):: ctens
         end function tensBlck_present
  !Free the HAB entry occupied by tensBlck_t:
         integer(C_INT) function tensBlck_hab_free(ctens) bind(c,name='tensBlck_hab_free')
          import
          implicit none
          type(C_PTR), value:: ctens
         end function tensBlck_hab_free
  !Get the number of tensor elements (volume) in tensBlck_t:
         integer(C_SIZE_T) function tensBlck_volume(ctens) bind(c,name='tensBlck_volume')
          import
          implicit none
          type(C_PTR), value, intent(in):: ctens
         end function tensBlck_volume
 !NV-TAL CUDA task API:
  !Create a CUDA task:
         integer(C_INT) function cuda_task_create(cuda_task) bind(c,name='cuda_task_create')
          import
          implicit none
          type(C_PTR), intent(out):: cuda_task
         end function cuda_task_create
  !Destroy a CUDA task:
         integer(C_INT) function cuda_task_destroy(cuda_task) bind(c,name='cuda_task_destroy')
          import
          implicit none
          type(C_PTR), value:: cuda_task
         end function cuda_task_destroy
  !Clean a CUDA task for reuse:
         integer(C_INT) function cuda_task_clean(cuda_task) bind(c,name='cuda_task_clean')
          import
          implicit none
          type(C_PTR), value:: cuda_task
         end function cuda_task_clean
  !Get the CUDA GPU ID from cudaTask_t:
         integer(C_INT) function cuda_task_gpu_id(cuda_task) bind(c,name='cuda_task_gpu_id')
          import
          implicit none
          type(C_PTR), value, intent(in):: cuda_task
         end function cuda_task_gpu_id
  !Get the CUDA task status:
         integer(C_INT) function cuda_task_status(cuda_task) bind(c,name='cuda_task_status')
          import
          implicit none
          type(C_PTR), value:: cuda_task
         end function cuda_task_status
  !Query CUDA task completion:
         integer(C_INT) function cuda_task_complete(cuda_task) bind(c,name='cuda_task_complete')
          import
          implicit none
          type(C_PTR), value:: cuda_task
         end function cuda_task_complete
  !Wait on completion of a CUDA task:
         integer(C_INT) function cuda_task_wait(cuda_task) bind(c,name='cuda_task_wait')
          import
          implicit none
          type(C_PTR), value:: cuda_task
         end function cuda_task_wait
  !Get the task timing in seconds:
         real(C_FLOAT) function cuda_task_time(cuda_task,in_copy,out_copy,comp) bind(c,name='cuda_task_time')
          import
          implicit none
          type(C_PTR), value, intent(in):: cuda_task
          real(C_FLOAT), intent(out):: in_copy
          real(C_FLOAT), intent(out):: out_copy
          real(C_FLOAT), intent(out):: comp
         end function cuda_task_time
 !NV-TAL query/action API:
  !Check whether GPU belongs to the current MPI process:
         integer(C_INT) function gpu_is_mine(gpu_num) bind(c,name='gpu_is_mine')
          import
          implicit none
          integer(C_INT), value:: gpu_num
         end function gpu_is_mine
  !Returns the ID of the least busy NVidia GPU (which belongs to the MPI process):
         integer(C_INT) function gpu_busy_least() bind(c,name='gpu_busy_least')
          import
          implicit none
         end function gpu_busy_least
  !Activate a specific GPU (only if it belongs to the MPI process):
         integer(C_INT) function gpu_activate(gpu_num) bind(c,name='gpu_activate')
          import
          implicit none
          integer(C_INT), value:: gpu_num
         end function gpu_activate
 !NV-TAL tensor operations:
  !Copy tensor block data from the Host argument buffer to a GPU argument buffer (blocking):
         integer(C_INT) function gpu_put_arg(ctens) bind(c,name='gpu_put_arg')
          import
          implicit none
          type(C_PTR), value:: ctens
         end function gpu_put_arg
  !Copy tensor block data from a GPU argument buffer to the Host argument buffer (blocking):
         integer(C_INT) function gpu_get_arg(ctens) bind(c,name='gpu_get_arg')
          import
          implicit none
          type(C_PTR), value:: ctens
         end function gpu_get_arg
  !Copy tensor block data from the Host argument buffer to a GPU argument buffer (non-blocking):
         integer(C_INT) function gpu_put_arg_(ctens,cuda_task) bind(c,name='gpu_put_arg_')
          import
          implicit none
          type(C_PTR), value:: ctens
          type(C_PTR), value:: cuda_task
         end function gpu_put_arg_
  !Copy tensor block data from a GPU argument buffer to the Host argument buffer (non-blocking):
         integer(C_INT) function gpu_get_arg_(ctens,cuda_task) bind(c,name='gpu_get_arg_')
          import
          implicit none
          type(C_PTR), value:: ctens
          type(C_PTR), value:: cuda_task
         end function gpu_get_arg_
  !Array 2-norm squared (R4):
         integer(C_INT) function gpu_array_2norm2_r4(asize,arr,norm2) bind(c,name='gpu_array_2norm2_r4')
          import
          implicit none
          integer(C_SIZE_T), value:: asize
          real(C_FLOAT), intent(in):: arr(*)
          real(C_FLOAT), intent(out):: norm2
         end function gpu_array_2norm2_r4
  !Array 2-norm squared (R8):
         integer(C_INT) function gpu_array_2norm2_r8(asize,arr,norm2) bind(c,name='gpu_array_2norm2_r8')
          import
          implicit none
          integer(C_SIZE_T), value:: asize
          real(C_DOUBLE), intent(in):: arr(*)
          real(C_DOUBLE), intent(out):: norm2
         end function gpu_array_2norm2_r8
  !Blocking matrix multiplication (TN variant, R4):
         integer(C_INT) function gpu_matrix_multiply_tn_r4(ll,lr,lc,lmat,rmat,dmat)&
                                  &bind(c,name='gpu_matrix_multiply_tn_r4')
          import
          implicit none
          integer(C_SIZE_T), value, intent(in):: ll
          integer(C_SIZE_T), value, intent(in):: lr
          integer(C_SIZE_T), value, intent(in):: lc
          real(C_FLOAT), intent(in):: lmat(*)
          real(C_FLOAT), intent(in):: rmat(*)
          real(C_FLOAT), intent(inout):: dmat(*)
         end function gpu_matrix_multiply_tn_r4
  !Blocking matrix multiplication (TN variant, R8):
         integer(C_INT) function gpu_matrix_multiply_tn_r8(ll,lr,lc,lmat,rmat,dmat)&
                                  &bind(c,name='gpu_matrix_multiply_tn_r8')
          import
          implicit none
          integer(C_SIZE_T), value, intent(in):: ll
          integer(C_SIZE_T), value, intent(in):: lr
          integer(C_SIZE_T), value, intent(in):: lc
          real(C_DOUBLE), intent(in):: lmat(*)
          real(C_DOUBLE), intent(in):: rmat(*)
          real(C_DOUBLE), intent(inout):: dmat(*)
         end function gpu_matrix_multiply_tn_r8
  !Tensor block initialization:
         integer(C_INT) function gpu_tensor_block_init_(ctens,val,copy_back,cuda_task)&
                                  &bind(c,name='gpu_tensor_block_init_')
          import
          implicit none
          type(C_PTR), value:: ctens
          real(C_DOUBLE), value:: val
          integer(C_INT), value:: copy_back
          type(C_PTR), value:: cuda_task
         end function gpu_tensor_block_init_
  !Tensor block rescaling:
         integer(C_INT) function gpu_tensor_block_scale_(ctens,val,copy_back,cuda_task)&
                                  &bind(c,name='gpu_tensor_block_scale_')
          import
          implicit none
          type(C_PTR), value:: ctens
          real(C_DOUBLE), value:: val
          integer(C_INT), value:: copy_back
          type(C_PTR), value:: cuda_task
         end function gpu_tensor_block_scale_
  !Tensor block addition:
         integer(C_INT) function gpu_tensor_block_add_dlf_(ctens0,ctens1,val,copy_back,cuda_task)&
                                  &bind(c,name='gpu_tensor_block_add_dlf_')
          import
          implicit none
          type(C_PTR), value:: ctens0
          type(C_PTR), value, intent(in):: ctens1
          real(C_DOUBLE), value:: val
          integer(C_INT), value:: copy_back
          type(C_PTR), value:: cuda_task
         end function gpu_tensor_block_add_dlf_
  !Blocking tensor transpose:
         integer(C_INT) function gpu_tensor_block_copy_dlf(dim_trn,tens_in,tens_out)&
                                  &bind(c,name='gpu_tensor_block_copy_dlf')
          import
          implicit none
          integer(C_INT), intent(in):: dim_trn(*)
          type(C_PTR), value:: tens_in
          type(C_PTR), value:: tens_out
         end function gpu_tensor_block_copy_dlf
  !Non-blocking tensor transpose:
         integer(C_INT) function gpu_tensor_block_copy_dlf_(dim_trn,tens_in,tens_out,copy_back,cuda_task)&
                                  &bind(c,name='gpu_tensor_block_copy_dlf_')
          import
          implicit none
          integer(C_INT), intent(in):: dim_trn(*)
          type(C_PTR), value:: tens_in
          type(C_PTR), value:: tens_out
          integer(C_INT), value:: copy_back
          type(C_PTR), value:: cuda_task
         end function gpu_tensor_block_copy_dlf_
  !Non-blocking tensor contraction:
         integer(C_INT) function gpu_tensor_block_contract_dlf_(cptrn,ltens,rtens,dtens,copy_back,cuda_task)&
                                  &bind(c,name='gpu_tensor_block_contract_dlf_')
          import
          implicit none
          integer(C_INT), intent(in):: cptrn(*)
          type(C_PTR), value:: ltens
          type(C_PTR), value:: rtens
          type(C_PTR), value:: dtens
          integer(C_INT), value:: copy_back
          type(C_PTR), value:: cuda_task
         end function gpu_tensor_block_contract_dlf_
#endif
        end interface

        end module tensor_algebra
