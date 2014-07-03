!This module contains interfaces for QFORCE C/C++ functions callable from Fortran.
       module extern_names

        interface
#ifndef NO_GPU
!C/C++ wrappers for CUDA Run-Time functions:
 !Get the total number of available GPU(Nvidia) devices on a node:
	 subroutine cudaGetDeviceCount(count,err_code) bind(c,name='cudagetdevicecount')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), intent(out):: count
	  integer(C_INT):: err_code
	 end subroutine cudaGetDeviceCount
 !Set active GPU(Nvidia) device:
	 subroutine cudaSetDevice(device,err_code) bind(c,name='cudasetdevice')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), value, intent(in):: device
	  integer(C_INT):: err_code
	 end subroutine cudaSetDevice
 !Get GPU(Nvidia) device properties:
	 subroutine cudaGetDeviceProperties(device,totalGlobalMem_,sharedMemPerBlock_,regsPerBlock_,warpSize_, &
	                                    maxThreadsPerBlock_,maxThreadsDim_,maxGridSize_,clockRate_,totalConstMem_, &
	                                    major_,minor_,deviceOverlap_,multiProcessorCount_,concurrentKernels_, &
	                                    ECCEnabled_,asyncEngineCount_,memoryClockRate_,memoryBusWidth_, &
	                                    maxThreadsPerMultiProcessor_,err_code) bind(c,name='cudagetdeviceproperties')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), value, intent(in):: device
	  integer(C_SIZE_T), intent(out):: totalGlobalMem_
	  integer(C_SIZE_T), intent(out):: sharedMemPerBlock_
	  integer(C_INT), intent(out):: regsPerBlock_
	  integer(C_INT), intent(out):: warpSize_
	  integer(C_INT), intent(out):: maxThreadsPerBlock_
	  integer(C_INT), intent(out):: maxThreadsDim_(3)
	  integer(C_INT), intent(out):: maxGridSize_(3)
	  integer(C_INT), intent(out):: clockRate_
	  integer(C_SIZE_T), intent(out):: totalConstMem_
	  integer(C_INT), intent(out):: major_
	  integer(C_INT), intent(out):: minor_
	  integer(C_INT), intent(out):: deviceOverlap_
	  integer(C_INT), intent(out):: multiProcessorCount_
	  integer(C_INT), intent(out):: concurrentKernels_
	  integer(C_INT), intent(out):: ECCEnabled_
	  integer(C_INT), intent(out):: asyncEngineCount_
	  integer(C_INT), intent(out):: memoryClockRate_
	  integer(C_INT), intent(out):: memoryBusWidth_
	  integer(C_INT), intent(out):: maxThreadsPerMultiProcessor_
	  integer(C_INT):: err_code
	 end subroutine cudaGetDeviceProperties
 !Device synchronization:
	 subroutine cudaDeviceSynchronize(err_code) bind(c,name='cudadevicesynchronize')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT):: err_code
	 end subroutine cudaDeviceSynchronize
#endif

#ifndef NO_GPU
!C/C++ wrappers for GPU(Nvidia) tensor algebra:
 !Array 2-norm squared (R4):
	 integer(C_INT) function gpu_array_2norm2_r4(asize,arr,norm2) bind(c,name='gpu_array_2norm2_r4')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_SIZE_T), value:: asize
	  real(C_FLOAT), intent(in):: arr(*)
	  real(C_FLOAT), intent(out):: norm2
	 end function gpu_array_2norm2_r4
 !Array 2-norm squared (R8):
	 integer(C_INT) function gpu_array_2norm2_r8(asize,arr,norm2) bind(c,name='gpu_array_2norm2_r8')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_SIZE_T), value:: asize
	  real(C_DOUBLE), intent(in):: arr(*)
	  real(C_DOUBLE), intent(out):: norm2
	 end function gpu_array_2norm2_r8
 !Blocking matrix multiplication (TN variant, R4):
	 integer(C_INT) function gpu_matrix_multiply_tn_r4(ll,lr,lc,lmat,rmat,dmat) &
	                                                   bind(c,name='gpu_matrix_multiply_tn_r4')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_SIZE_T), value, intent(in):: ll
	  integer(C_SIZE_T), value, intent(in):: lr
	  integer(C_SIZE_T), value, intent(in):: lc
	  real(C_FLOAT), intent(in):: lmat(*)
	  real(C_FLOAT), intent(in):: rmat(*)
	  real(C_FLOAT), intent(inout):: dmat(*)
	 end function gpu_matrix_multiply_tn_r4
 !Tensor block initialization:
	 integer(C_INT) function gpu_tensor_block_init_(ctens,val,copy_back,cuda_task) &
	                                               bind(c,name='gpu_tensor_block_init_')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), value:: ctens
	  real(C_DOUBLE), value:: val
	  integer(C_INT), value:: copy_back
	  type(C_PTR), value:: cuda_task
	 end function gpu_tensor_block_init_
 !Tensor block rescaling:
	 integer(C_INT) function gpu_tensor_block_scale_(ctens,val,copy_back,cuda_task) &
	                                                 bind(c,name='gpu_tensor_block_scale_')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), value:: ctens
	  real(C_DOUBLE), value:: val
	  integer(C_INT), value:: copy_back
	  type(C_PTR), value:: cuda_task
	 end function gpu_tensor_block_scale_
 !Tensor block addition:
	 integer(C_INT) function gpu_tensor_block_add_dlf_(ctens0,ctens1,val,copy_back,cuda_task) &
	                                                   bind(c,name='gpu_tensor_block_add_dlf_')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), value:: ctens0
	  type(C_PTR), value, intent(in):: ctens1
	  real(C_DOUBLE), value:: val
	  integer(C_INT), value:: copy_back
	  type(C_PTR), value:: cuda_task
	 end function gpu_tensor_block_add_dlf_
 !Blocking tensor transpose:
	 integer(C_INT) function gpu_tensor_block_copy_dlf(dim_trn,tens_in,tens_out) &
                                              bind(c,name='gpu_tensor_block_copy_dlf')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), intent(in):: dim_trn(*)
	  type(C_PTR), value:: tens_in
	  type(C_PTR), value:: tens_out
	 end function gpu_tensor_block_copy_dlf
 !Non-blocking tensor transpose:
	 integer(C_INT) function gpu_tensor_block_copy_dlf_(dim_trn,tens_in,tens_out,copy_back,cuda_task) &
                                                            bind(c,name='gpu_tensor_block_copy_dlf_')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), intent(in):: dim_trn(*)
	  type(C_PTR), value:: tens_in
	  type(C_PTR), value:: tens_out
	  integer(C_INT), value:: copy_back
	  type(C_PTR), value:: cuda_task
	 end function gpu_tensor_block_copy_dlf_
 !Non-blocking tensor contraction:
	 integer(C_INT) function gpu_tensor_block_contract_dlf_(cptrn,ltens,rtens,dtens,copy_back,cuda_task) &
                                                                bind(c,name='gpu_tensor_block_contract_dlf_')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), intent(in):: cptrn(*)
	  type(C_PTR), value:: ltens
	  type(C_PTR), value:: rtens
	  type(C_PTR), value:: dtens
	  integer(C_INT), value:: copy_back
	  type(C_PTR), value:: cuda_task
	 end function gpu_tensor_block_contract_dlf_
#endif

!Other C/C++ functions:
 !Obtain an offset byte pointer:
	 type(C_PTR) function ptr_offset(byte_ptr,byte_offset) bind(c,name='ptr_offset')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), value, intent(in):: byte_ptr
	  integer(C_SIZE_T), value, intent(in):: byte_offset
	 end function ptr_offset
 !Print C pointer from Fortran:
	 subroutine print_c_ptr(c_pointer) bind(c,name='print_c_ptr')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), value:: c_pointer
	 end subroutine print_c_ptr
 !Initialize all argument buffers on Host and Devices (GPU constant+global, MICs, etc):
	 integer(C_INT) function arg_buf_allocate(host_mem,arg_max,gpu_beg,gpu_end) bind(c,name='arg_buf_allocate')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_SIZE_T), intent(inout):: host_mem
	  integer(C_INT), intent(out):: arg_max
	  integer(C_INT), value, intent(in):: gpu_beg
	  integer(C_INT), value, intent(in):: gpu_end
	 end function arg_buf_allocate
 !Deallocate all argument buffers on Host and Devices:
	 integer(C_INT) function arg_buf_deallocate(gpu_beg,gpu_end) bind(c,name='arg_buf_deallocate')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), value, intent(in):: gpu_beg
	  integer(C_INT), value, intent(in):: gpu_end
	 end function arg_buf_deallocate
 !Check whether the Host argument buffer is clean:
	 integer(C_INT) function arg_buf_clean_host() bind(c,name='arg_buf_clean_host')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	 end function arg_buf_clean_host
#ifndef NO_GPU
 !Check whether a GPU argument buffer is clean:
	 integer(C_INT) function arg_buf_clean_gpu(gpu_num) bind(c,name='arg_buf_clean_gpu')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), value, intent(in):: gpu_num
	 end function arg_buf_clean_gpu
 !Get the current GPU error count:
	 integer(C_INT) function gpu_get_error_count() bind(c,name='gpu_get_error_count')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	 end function gpu_get_error_count
 !Get the current GPU debug dump:
	 integer(C_INT) function gpu_get_debug_dump(dump) bind(c,name='gpu_get_debug_dump')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), intent(out):: dump(*)
	 end function gpu_get_debug_dump
#endif
 !Get buffered block sizes for each level of the Host argument buffer:
	 integer(C_INT) function get_blck_buf_sizes_host(blck_sizes) bind(c,name='get_blck_buf_sizes_host')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_SIZE_T), intent(out):: blck_sizes(*)
	 end function get_blck_buf_sizes_host
#ifndef NO_GPU
 !Get buffered block sizes for each level of the GPU argument buffer:
	 integer(C_INT) function get_blck_buf_sizes_gpu(gpu_num,blck_sizes) bind(c,name='get_blck_buf_sizes_gpu')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), value, intent(in):: gpu_num
	  integer(C_SIZE_T), intent(out):: blck_sizes(*)
	 end function get_blck_buf_sizes_gpu
#endif
 !Get a free argument entry in the Host argument buffer:
	 integer(C_INT) function get_buf_entry_host(bsize,entry_ptr,entry_num) bind(c,name='get_buf_entry_host')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_SIZE_T), value, intent(in):: bsize
	  type(C_PTR), intent(out):: entry_ptr
	  integer(C_INT), intent(out):: entry_num
	 end function get_buf_entry_host
 !Free an argument entry in the Host argument buffer:
	 integer(C_INT) function free_buf_entry_host(entry_num) bind(c,name='free_buf_entry_host')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), value, intent(in):: entry_num
	 end function free_buf_entry_host
#ifndef NO_GPU
 !Get a free argument entry in the GPU argument buffer:
	 integer(C_INT) function get_buf_entry_gpu(gpu_num,bsize,entry_ptr,entry_num) bind(c,name='get_buf_entry_gpu')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), value, intent(in):: gpu_num
	  integer(C_SIZE_T), value, intent(in):: bsize
	  type(C_PTR), intent(out):: entry_ptr
	  integer(C_INT), intent(out):: entry_num
	 end function get_buf_entry_gpu
 !Free an argument entry in the GPU argument buffer:
	 integer(C_INT) function free_buf_entry_gpu(gpu_num,entry_num) bind(c,name='free_buf_entry_gpu')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), value, intent(in):: gpu_num
	  integer(C_INT), value, intent(in):: entry_num
	 end function free_buf_entry_gpu
 !Get a free entry from the GPU constant memory argument bank:
	 integer(C_INT) function const_args_entry_get(gpu_num,entry_num) bind(c,name='const_args_entry_get')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), value, intent(in):: gpu_num
	  integer(C_INT), intent(out):: entry_num
	 end function const_args_entry_get
 !Return back an entry to the GPU constant memory argument bank:
	 integer(C_INT) function const_args_entry_free(gpu_num,entry_num) bind(c,name='const_args_entry_free')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), value, intent(in):: gpu_num
	  integer(C_INT), value, intent(in):: entry_num
	 end function const_args_entry_free
 !Turn on/off CUDA timing events:
	 subroutine gpu_set_event_policy(alg) bind(c,name='gpu_set_event_policy')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), value:: alg
	 end subroutine gpu_set_event_policy
 !Set tensor transpose algorithm:
	 subroutine gpu_set_transpose_algorithm(alg) bind(c,name='gpu_set_transpose_algorithm')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), value:: alg
	 end subroutine gpu_set_transpose_algorithm
 !Set matrix multiplication algorithm:
	 subroutine gpu_set_matmult_algorithm(alg) bind(c,name='gpu_set_matmult_algorithm')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), value:: alg
	 end subroutine gpu_set_matmult_algorithm
 !Create an instance of tensBlck_t:
	 integer(C_INT) function tensBlck_create(ctens) bind(c,name='tensBlck_create')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), intent(out):: ctens
	 end function tensBlck_create
 !Destroy an instance of tensBlck_t:
	 integer(C_INT) function tensBlck_destroy(ctens) bind(c,name='tensBlck_destroy')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), value:: ctens
	 end function tensBlck_destroy
 !Construct tensBlck_t:
	 integer(C_INT) function tensBlck_construct(ctens,dev_kind,dev_num,data_kind,trank, &
                                  addr_dims,addr_divs,addr_grps,addr_prmn,addr_host,addr_gpu, &
                                  entry_gpu,entry_const) bind(c,name='tensBlck_construct')
	  use, intrinsic:: ISO_C_BINDING
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
	  integer(C_INT), value, intent(in):: entry_gpu
	  integer(C_INT), value, intent(in):: entry_const
	 end function tensBlck_construct
 !Get the Accelerator ID and other information from tensBlck_t:
	 integer(C_INT) function tensBlck_acc_id(ctens,dev_kind,entry_gpu,entry_const,data_kind,there) &
                                                 bind(c,name='tensBlck_acc_id')
	  use, intrinsic:: ISO_C_BINDING
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
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), value:: ctens
	 end function tensBlck_set_presence
 !Unmark tensBlck_t as residing on GPU:
	 integer(C_INT) function tensBlck_set_absence(ctens) bind(c,name='tensBlck_set_absence')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), value:: ctens
	 end function tensBlck_set_absence
 !Nullify the Host memory pointer in tensBlck_t:
         integer(C_INT) function tensBlck_hab_null(ctens) bind(c,name='tenBlck_hab_null')
          use, intrinsic:: ISO_C_BINDING
          implicit none
          type(C_PTR), value:: ctens
         end function tensBlck_hab_null
 !Number of tensor elements (volume) in tensBlck_t:
	 integer(C_SIZE_T) function tensBlck_volume(ctens) bind(c,name='tensBlck_volume')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), value, intent(in):: ctens
	 end function tensBlck_volume
 !Create a CUDA task:
	 integer(C_INT) function cuda_task_create(cuda_task) bind(c,name='cuda_task_create')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), intent(out):: cuda_task
	 end function cuda_task_create
 !Destroy a CUDA task:
	 integer(C_INT) function cuda_task_destroy(cuda_task) bind(c,name='cuda_task_destroy')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), value:: cuda_task
	 end function cuda_task_destroy
 !Get CUDA device ID from cudaTask_t:
	 integer(C_INT) function cuda_task_gpu_id(cuda_task) bind(c,name='cuda_task_gpu_id')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), value, intent(in):: cuda_task
	 end function cuda_task_gpu_id
 !Get CUDA task status:
	 integer(C_INT) function cuda_task_status(cuda_task) bind(c,name='cuda_task_status')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), value:: cuda_task
	 end function cuda_task_status
 !Query CUDA task completion:
	 integer(C_INT) function cuda_task_complete(cuda_task) bind(c,name='cuda_task_complete')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), value:: cuda_task
	 end function cuda_task_complete
 !Wait on accomplishment of a CUDA task:
	 integer(C_INT) function cuda_task_wait(cuda_task) bind(c,name='cuda_task_wait')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), value:: cuda_task
	 end function cuda_task_wait
 !Get task time in seconds:
	 real(C_FLOAT) function cuda_task_time(cuda_task,in_copy,out_copy,comp) bind(c,name='cuda_task_time')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), value, intent(in):: cuda_task
	  real(C_FLOAT), intent(out):: in_copy
	  real(C_FLOAT), intent(out):: out_copy
	  real(C_FLOAT), intent(out):: comp
	 end function cuda_task_time
 !Check whether GPU belongs to the current MPI process:
	 integer(C_INT) function gpu_is_mine(gpu_num) bind(c,name='gpu_is_mine')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), value:: gpu_num
	 end function gpu_is_mine
 !Returns the ID of the least busy NVidia GPU:
	 integer(C_INT) function gpu_busy_least() bind(c,name='gpu_busy_least')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	 end function gpu_busy_least
 !Activates GPU only if it belongs to the MPI process:
	 integer(C_INT) function gpu_activate(gpu_num) bind(c,name='gpu_activate')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  integer(C_INT), value:: gpu_num
	 end function gpu_activate
 !Copy tensor block data from the Host argument buffer to a GPU argument buffer (blocking):
	 integer(C_INT) function gpu_put_arg(ctens) bind(c,name='gpu_put_arg')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), value:: ctens
	 end function gpu_put_arg
 !Copy tensor block data from a GPU argument buffer to the Host argument buffer (blocking):
	 integer(C_INT) function gpu_get_arg(ctens) bind(c,name='gpu_get_arg')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), value:: ctens
	 end function gpu_get_arg
 !Copy tensor block data from the Host argument buffer to a GPU argument buffer (non-blocking):
	 integer(C_INT) function gpu_put_arg_(ctens,cuda_task) bind(c,name='gpu_put_arg_')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), value:: ctens
	  type(C_PTR), value:: cuda_task
	 end function gpu_put_arg_
 !Copy tensor block data from a GPU argument buffer to the Host argument buffer (non-blocking):
	 integer(C_INT) function gpu_get_arg_(ctens,cuda_task) bind(c,name='gpu_get_arg_')
	  use, intrinsic:: ISO_C_BINDING
	  implicit none
	  type(C_PTR), value:: ctens
	  type(C_PTR), value:: cuda_task
	 end function gpu_get_arg_
#endif
        end interface

       end module extern_names
