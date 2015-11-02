/** Tensor Algebra Library for NVidia GPU: NV-TAL (CUDA based).
REVISION: 2015/11/02
Copyright (C) 2015 Dmitry I. Lyakh (email: quant4me@gmail.com)
Copyright (C) 2015 Oak Ridge National Laboratory (UT-Battelle)

This source file is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
-------------------------------------------------------------------------------
OPTIONS:
 # -D CUDA_ARCH=350: target GPU compute capability (default is 130);
 # -D NO_GPU: disables GPU usage;
 # -D NO_BLAS: disables cuBLAS calls, they will be replaced by in-house routines (slower);
 # -D DEBUG_GPU: collection of debugging information will be activated;
NOTES:
 # Minimal required compute capability is 1.1 (1.3 for double precision);
 # cuBLAS.v2 is required when BLAS is enabled;
 # Functions without underscores at the end of their names are blocking functions;
   Functions with one underscore at the end of their names are non-blocking functions;
   Functions with two underscores at the end of their names are (non-blocking) CUDA kernels.
 # Non-blocking tensor algebra functions carry an additional input argument <coherence_ctrl>
   which controls the tensor argument synchronization accross different devices.
 # Non-blocking tensor algebra functions carry an additional output argument <cuda_task> (task handle).
 # Seems like cudaEventRecord() issued in different streams can serialize the stream
   execution for some older compute capabilities. EVENT_RECORD=0 will disable event recording.
   If GPU timing is needed, event recording has to be enabled (EVENT_RECORD=1).
FOR DEVELOPERS ONLY:
 # Current device resources:
    - Global memory pointer (any device);
    - Argument buffer entry handle (any device);
    - Multi-index-entry pointer (any device, entry length = MAX_TENSOR_RANK);
    - Constant-memory-entry handle (Nvidia GPU);
    - CUDA stream handle (Nvidia GPU);
    - CUDA event handle (Nvidia GPU).
 # A life cycle of a C object (for example, tensBlck_t):
    a) Allocate memory for the object, if needed: Suffix _alloc or _create (includes cleaning);
    b) Clean (initialize to null) an allocated (empty) object: Suffix _clean (normally included in _create);
    c) Construct (define or redefine) an existing object (resources will be acquired/released): Suffix _construct;
    d) Destruct a defined object (resources will be released, the object will be initialized to null): Suffix _destruct;
    e) Free the memory occupied by a clean (initialized to null or destructed) object: Suffix _free or _destroy (may include _destruct).
   Thus, as a rule, the device resource acquisition/release occurs solely in _construct and _destruct functions.
 # A state of a C object:
    a) Undefined: After the memory allocation (either dynamic or static);
    b) Defined nullified: After cleaning or destruction;
    c) Defined to a value: After construction;
    d) Dead: After memory deallocation (if it was allocated dynamically).
 # Resource acquisition/release:
    - Tensor block constructor/destructor acquires/releases global memory resources, including
      both pointers and buffer handles, as well as multi-index bank entries (pinned Host memory).
    - Tensor operation scheduling functions acquire GPU global memory resources,
      GPU constant memory resources, multi-index entries, CUDA stream and event handles (Nvidia GPU).
    - CUDA task completion/error check functions release GPU global memory resources,
      GPU constant memory resources, multi-index entries, CUDA stream and event handles (Nvidia GPU).
 # Functions which construct tensor blocks or perform asynchronous operations on them
   allocate resources (global/constant memory, etc). In case the corresponding resource
   allocator returns TRY_LATER or DEVICE_UNABLE (or an error), the corresponding function
   must clean the partially created tensor block or the CUDA task before returning:
   Thus, the corresponding object will be kept in its initial state if no SUCCESS.
TO BE FIXED:
 # The pinned multi-index slab is only initialized when NVidia GPU is enabled (in <init_gpus>).
   Consequently, I should probably move this initialization into <talsh_init> instead.
**/

#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#ifndef NO_GPU

#include <cuda.h>
#include <cuda_runtime.h>

#ifndef NO_BLAS
#include <cublas_v2.h>
#endif

#endif

#include "tensor_algebra.h"

#ifndef NO_GPU
//PARAMETERS:
#define GPU_DEBUG_DUMP_SIZE 128 //size of the GPU debug dump (int array)
#endif
//----------------------------------------------------------------------
//FUNCTION PROTOTYPES:
// IMPORTED:
#ifdef __cplusplus
extern "C" {
#endif
 void get_contr_permutations(int lrank, int rrank, const int *cptrn, int *dprm, int *lprm, int *rprm,
                             int *ncd, int *nlu, int *nru, int *ierr);
#ifdef __cplusplus
}
#endif
// LOCAL (PRIVATE):
static int tens_valid_data_kind(int datk);
static int prmn_convert(int n, const int *o2n, int *n2o);
static int non_trivial_prmn(int n, const int *prm);
static int tensDevRsc_create(talsh_dev_rsc_t **drsc);
static int tensDevRsc_clean(talsh_dev_rsc_t * drsc);
static int tensDevRsc_empty(talsh_dev_rsc_t * drsc);
static int tensDevRsc_attach_mem(talsh_dev_rsc_t * drsc, int dev_id, void * mem_p, int buf_entry);
static int tensDevRsc_detach_mem(talsh_dev_rsc_t * drsc);
static int tensDevRsc_allocate_mem(talsh_dev_rsc_t * drsc, int dev_id, size_t mem_size, int in_arg_buf);
static int tensDevRsc_free_mem(talsh_dev_rsc_t * drsc);
#ifndef NO_GPU
static int tensDevRsc_get_const_entry(talsh_dev_rsc_t * drsc, int dev_id);
static int tensDevRsc_release_const_entry(talsh_dev_rsc_t * drsc);
#endif
static int tensDevRsc_release_all(talsh_dev_rsc_t * drsc);
static int tensDevRsc_destroy(talsh_dev_rsc_t * drsc);
static int mi_entry_get(int ** mi_entry);
static int mi_entry_release(int * mi_entry);
static int mi_entry_pinned(int * mi_entry);
#ifndef NO_GPU
static int cuda_stream_get(int gpu_num, int * cuda_stream_handle);
static int cuda_stream_release(int gpu_num, int cuda_stream_handle);
static int cuda_event_get(int gpu_num, int * cuda_event_handle);
static int cuda_event_release(int gpu_num, int cuda_event_handle);
static void limit_cuda_blocks2d(int max_blocks, int *bx, int *by);
static int cuda_task_finalize(cudaTask_t *cuda_task, int err_code, int gpu_num);
static int cuda_task_record(cudaTask_t *cuda_task, int err_code, int gpu_num, cudaStream_t cuda_stream,
            cudaEvent_t cuda_start, cudaEvent_t cuda_comput, cudaEvent_t cuda_output, cudaEvent_t cuda_finish,
            int scr_entry_cnt, int *scr_entries);
// CUDA KERNELS:
__global__ void gpu_array_2norm2_r4__(size_t arr_size, const float *arr, float *bnorm2);
__global__ void gpu_array_2norm2_r8__(size_t arr_size, const double *arr, double *bnorm2);
__global__ void gpu_array_init_r4__(size_t tsize, float *arr, float val);
__global__ void gpu_array_init_r8__(size_t tsize, double *arr, double val);
__global__ void gpu_array_scale_r4__(size_t tsize, float *arr, float val);
__global__ void gpu_array_scale_r8__(size_t tsize, double *arr, double val);
__global__ void gpu_array_add_r4__(size_t tsize, float* __restrict__ arr0, const float* __restrict__ arr1, float val);
__global__ void gpu_array_add_r8__(size_t tsize, double* __restrict__ arr0, const double* __restrict__ arr1, double val);
__global__ void gpu_array_dot_product_r4__(size_t tsize, const float *arr1, const float *arr2, volatile float *dprod);
__global__ void gpu_array_dot_product_r8__(size_t tsize, const double *arr1, const double *arr2, volatile double *dprod);
__global__ void gpu_array_product_r4__(size_t tsize1, const float* __restrict__ arr1, size_t tsize2,
                                       const float* __restrict__ arr2, float* __restrict__ arr0);
__global__ void gpu_array_product_r8__(size_t tsize1, const double* __restrict__ arr1, size_t tsize2,
                                       const double* __restrict__ arr2, double* __restrict__ arr0);
__global__ void gpu_tensor_block_copy_dlf_r4__(int dmo, int drc, int dim_num, int const_args_pos,
                                               const float* __restrict__ tens_in, float* __restrict__ tens_out);
__global__ void gpu_tensor_block_copy_dlf_r8__(int dmo, int drc, int dim_num, int const_args_pos,
                                               const double* __restrict__ tens_in, double* __restrict__ tens_out);
__global__ void gpu_tensor_block_copy_scatter_dlf_r4__(int dmo, int drc, int dim_num, int const_args_pos,
                                                       const float* __restrict__ tens_in, float* __restrict__ tens_out);
__global__ void gpu_tensor_block_copy_scatter_dlf_r8__(int dmo, int drc, int dim_num, int const_args_pos,
                                                       const double* __restrict__ tens_in, double* __restrict__ tens_out);
__global__ void gpu_matrix_multiply_tn_r4__(size_t ll, size_t lr, size_t lc, const float* __restrict__ arg1,
                                            const float* __restrict__ arg2, float* __restrict__ arg0);
__global__ void gpu_matrix_multiply_tn_r8__(size_t ll, size_t lr, size_t lc, const double* __restrict__ arg1,
                                            const double* __restrict__ arg2, double* __restrict__ arg0);
#endif
//------------------------------------------------------------------------------------------------------
//PARAMETERS:
static int VERBOSE=1; //verbosity for error messages
#ifndef NO_GPU
//GLOBAL DATA:
// GPU control on the current MPI process:
static int gpu_up[MAX_GPUS_PER_NODE]; //GPU_OFF(0): GPU is disabled; GPU_MINE(1): GPU is enabled; GPU_MINE_CUBLAS(2): GPU is BLAS enabled
static cudaDeviceProp gpu_prop[MAX_GPUS_PER_NODE]; //properties of all GPUs present on the node
static talsh_stats_t gpu_stats[MAX_GPUS_PER_NODE]; //runtime statistics for all GPUs present on the node
#ifndef NO_BLAS
// Infrastructure for CUBLAS:
static cublasHandle_t cublas_handle[MAX_GPUS_PER_NODE]; //each GPU present on a node obtains its own cuBLAS context handle
#endif
#endif
// Slab for the multi-index storage (will be pinned):
static int miBank[MAX_GPU_ARGS*MAX_MLNDS_PER_TENS][MAX_TENSOR_RANK]; //All active .dims[], .divs[], .grps[] will be stored here
static int miFreeHandle[MAX_GPU_ARGS*MAX_MLNDS_PER_TENS]; //free entries for storing multi-indices
static int miFFE; //number of free handles left in miBank

#ifndef NO_GPU
// Slabs for the GPU asynchronous resources:
//  CUDA stream handles:
static cudaStream_t CUDAStreamBank[MAX_GPUS_PER_NODE][MAX_CUDA_TASKS]; //pre-allocated CUDA stream handles (for each CUDA device)
static int CUDAStreamFreeHandle[MAX_GPUS_PER_NODE][MAX_CUDA_TASKS]; //free CUDA stream handles
static int CUDAStreamFFE[MAX_GPUS_PER_NODE]; //number of free handles left in CUDAStreamFreeHandle
//  CUDA event handles:
static cudaEvent_t CUDAEventBank[MAX_GPUS_PER_NODE][MAX_CUDA_EVENTS]; //pre-allocated CUDA event handles (for each CUDA device)
static int CUDAEventFreeHandle[MAX_GPUS_PER_NODE][MAX_CUDA_EVENTS]; //free CUDA event handles
static int CUDAEventFFE[MAX_GPUS_PER_NODE]; //number of free handles left in CUDAEventFreeHandle

// Slab of GPU constant memory arguments for each GPU (managed by "mem_manager.cu"):
__device__ __constant__ int const_args_dims[MAX_GPU_ARGS][MAX_TENSOR_RANK]; //storage for device constant memory arguments: dimension extents
__device__ __constant__ int const_args_prmn[MAX_GPU_ARGS][MAX_TENSOR_RANK]; //storage for device constant memory arguments: permutation

// GPU error control and debugging for each GPU:
__device__ int gpu_error_count=0; //total number of CUDA errors registered on device till the current moment
__device__ int gpu_debug_dump[GPU_DEBUG_DUMP_SIZE]; //debug dump

// Global CUDA event recording policy:
static int EVENT_RECORD=1; //non-zero value enables CUDA event recording
static int PRINT_TIMING=1; //non-zero value enables time printing statements
// Infrastructure for function <gpu_tensor_block_copy_dlf> (blocking and non-blocking):
static int TRANS_SHMEM=1; //switch between shared-memory tensor transpose (1) and scatter tensor transpose (0)
// Infrastructure for <gpu_tensor_block_contract_dlf_> (non-blocking):
#ifndef NO_BLAS
static int DISABLE_BLAS=0; //non-zero value will disable cuBLAS usage (if it had been cuBLAS compiled/linked)
#else
static int DISABLE_BLAS=1; //non-zero value will disable cuBLAS usage (if it had been cuBLAS compiled/linked)
#endif
static cudaTask_t * LastTask[MAX_GPUS_PER_NODE]; //last CUDA task successfully scheduled on each GPU
__device__ __constant__ float sgemm_alpha=1.0f; //default alpha constant for SGEMM
__device__ __constant__ float sgemm_beta=1.0f;  //default beta constant SGEMM
__device__ __constant__ double dgemm_alpha=1.0; //default alpha constant for DGEMM
__device__ __constant__ double dgemm_beta=1.0;  //default beta constant DGEMM
// Infrastructure for functions <gpu_array_norm2_XX>:
__device__ int norm2_wr_lock=0; //write lock (shared by all <gpu_array_norm2_XX> running on device)
// Infrastructure for kernels <gpu_array_dot_product_XX__>:
__device__ int dot_product_wr_lock=0; //write lock (shared by all <gpu_array_dot_product_XX__> running on device)

// Infrastructure for functions <gpu_array_2norm2_XX> (blocking)`Obsolete:
__device__ float gpu_blck_norms2_r4[MAX_CUDA_BLOCKS]; //`Obsolete
__device__ double gpu_blck_norms2_r8[MAX_CUDA_BLOCKS]; //`Obsolete
static float blck_norms2_r4[MAX_CUDA_BLOCKS];  //`Obsolete `Not multi-GPU safe
static double blck_norms2_r8[MAX_CUDA_BLOCKS]; //`Obsolete `Not multi-GPU safe
#endif
//-------------------------------------------------------------------------------------
//DEVICE ID CONVERSION:
int encode_device_id(int dev_kind, int dev_num)
/** Given a device ID <dev_num> within its kind <dev_kind>, returns the flat device ID.
    DEV_MAX value on return means that the arguments were invalid. **/
{
 int dev_id=DEV_MAX; //Return of this value (= outside devices range) will mean that the arguments were invalid
 switch(dev_kind){
  case DEV_HOST: if(dev_num == 0) dev_id=0; break;
  case DEV_NVIDIA_GPU: if(dev_num >= 0 && dev_num < MAX_GPUS_PER_NODE) dev_id=1+dev_num; break;
  case DEV_INTEL_MIC: if(dev_num >= 0 && dev_num < MAX_MICS_PER_NODE) dev_id=1+MAX_GPUS_PER_NODE+dev_num; break;
  case DEV_AMD_GPU: if(dev_num >= 0 && dev_num < MAX_AMDS_PER_NODE) dev_id=1+MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE+dev_num; break;
  default: dev_id=DEV_MAX; //unknown device kind
 }
 return dev_id;
}

int decode_device_id(int dev_id, int *dev_kind)
/** Given a flat device ID <dev_id>, returns the device kind <dev_kind>
    and its kind-specific ID (>=0) as the return value.
    A negative return status means an invalid <dev_id> was passed. **/
{
 int dvn=-1; //Negative return value will correspond to an invalid <dev_id>
 int dvid=abs(dev_id); //flat device id is defined up to a sign
 if(dvid == 0){ //Host
  *dev_kind=DEV_HOST; dvn=0;
 }else if(dvid >= 1 && dvid <= MAX_GPUS_PER_NODE){ //Nvidia GPU
  *dev_kind=DEV_NVIDIA_GPU; dvn=dvid-1;
 }else if(dvid >= 1+MAX_GPUS_PER_NODE && dvid <= MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE){ //Intel MIC
  *dev_kind=DEV_INTEL_MIC; dvn=dvid-1-MAX_GPUS_PER_NODE;
 }else if(dvid >= 1+MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE && dvid <= MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE+MAX_AMDS_PER_NODE){ //AMD GPU
  *dev_kind=DEV_AMD_GPU; dvn=dvid-1-MAX_GPUS_PER_NODE-MAX_MICS_PER_NODE;
 }
 return dvn; //ID of the device within its kind
}

#ifndef NO_GPU
//GPU DEBUG FUNCTIONS:
__host__ int gpu_get_error_count()
/** Returns the total number of CUDA errors occured on current GPU.
    A negative return status means an error occurred. **/
{
 int i;
 cudaError_t err=cudaMemcpyFromSymbol((void*)&i,gpu_error_count,sizeof(gpu_error_count),0,cudaMemcpyDeviceToHost);
 if(err == cudaSuccess){return i;}else{return -1;}
}

__host__ int gpu_get_debug_dump(int *dump)
/** Returns the debug dump (int array) from current GPU.
    A positive return status is the length of the debug dump.
    A negative return status means an error occurred. **/
{
 cudaError_t err=cudaMemcpyFromSymbol((void*)dump,gpu_debug_dump,sizeof(int)*GPU_DEBUG_DUMP_SIZE,0,cudaMemcpyDeviceToHost);
 if(err == cudaSuccess){return GPU_DEBUG_DUMP_SIZE;}else{return -1;}
}
#endif

//AUXILIARY FUNCTIONS:
static int tens_valid_data_kind(int datk)
/** Returns YEP if the data kind <datk> is valid in TAL-SH, NOPE otherwise. **/
{
 if(datk == R4 || datk == R8 || datk == C8 || datk == NO_TYPE) return YEP;
 return NOPE;
}

static int prmn_convert(int n, const int *o2n, int *n2o)
/** Converts an O2N permutation into N2O (length = n). Both permutations
    are sign-free and the numeration starts from 1. **/
{
 int i,j;
 if(n >= 0){
  for(i=0;i<n;i++){j=o2n[i]-1; if(j >= 0 && j < n){n2o[j]=i+1;}else{return 1;}}
 }else{
  return 2;
 }
 return 0;
}

static int non_trivial_prmn(int n, const int *prm)
/** Returns 0 if the permutation prm[0:n-1] is trivial, 1 otherwise.
    The permutation is sign-free and the numeration starts from 1. No error check. **/
{
 int i,f=0;
 for(i=0;i<n;i++){if(prm[i] != i+1){f=1; break;}}
 return f;
}

//DEVICE RESOURCE MANAGEMENT:
static int tensDevRsc_create(talsh_dev_rsc_t **drsc)
/** Creates a new device resource descriptor and inits it to null. **/
{
 int errc;
 *drsc=(*talsh_dev_rsc_t)malloc(sizeof(talsh_dev_rsc_t)); if(*drsc == NULL) return 1;
 errc=tensDevRsc_clean(*drsc);
 return 0;
}

static int tensDevRsc_clean(talsh_dev_rsc_t * drsc)
/** Cleans (initializes to null) a device resource for tensBlck_t. **/
{
 if(drsc != NULL){
  drsc->dev_id=DEV_NULL;    //flat device id
  drsc->gmem_p=NULL;        //device global memory pointer (any device)
  drsc->buf_entry=-1;       //device argument buffer entry (any device)
  drsc->const_mem_entry=-1; //constant memory bank entry (Nvidia GPU only)
 }else{
  return -1;
 }
 return 0;
}

static int tensDevRsc_empty(talsh_dev_rsc_t * drsc)
/** Returns YEP if the device resource descriptor is empty, NOPE otherwise.
    Negative return status means an error. **/
{
 if(drsc == NULL) return -1;
 if(drsc->gmem_p != NULL || drsc->buf_entry >= 0 || drsc->const_mem_entry >= 0) return NOPE;
 drsc->dev_id=DEV_NULL;
 return YEP;
}

static int tensDevRsc_attach_mem(talsh_dev_rsc_t * drsc, int dev_id, void * mem_p, int buf_entry = -1)
/** Attaches a chunk of existing global memory to a device resource descriptor.
    If <buf_entry> >= 0, that means that the global memory is in the argument buffer.
    If the resource descriptor had already been assigned a device, the <dev_id>
    argument must match that one. **/
{
 if(drsc == NULL) return -1;
 if(dev_id < 0) return -2;
 if(mem_p == NULL) return -3;
 if(drsc->dev_id >= 0 && drsc->dev_id != dev_id) return 1; //resource was assigned to a different device
 if(drsc->gmem_p != NULL) return 2; //resource already has global memory attached
 drsc->dev_id=dev_id; drsc->gmem_p=mem_p; drsc->buf_entry=buf_entry;
 return 0;
}

static int tensDevRsc_detach_mem(talsh_dev_rsc_t * drsc)
/** Detaches a chunk of external memory from a device resource descriptor.
    Regardless of the origin of that memory, it is not released. **/
{
 int errc;

 if(drsc == NULL) return -1;
 if(drsc->dev_id < 0) return -2; //empty resource descriptor
 if(drsc->gmem_p == NULL) return 1; //no global memory attached
 drsc->gmem_p=NULL; drsc->buf_entry=-1;
 errc=tensDevRsc_empty(drsc);
 return 0;
}

static int tensDevRsc_allocate_mem(talsh_dev_rsc_t * drsc, int dev_id, size_t mem_size, int in_arg_buf = NOPE)
/** Allocates global memory on device <dev_id> and attaches it to a device resource descriptor.
    If <in_arg_buf> = YEP, the memory will be allocated via that device's argument buffer.
    A return status TRY_LATER or DEVICE_UNABLE indicates the resource shortage and is not an error. **/
{
 int i,devk,devn,errc;
 char *byte_ptr;

 if(drsc == NULL) return -1;
 if(dev_id < 0) return -2;
 if(mem_size <= 0) return -3;
 devn=decode_device_id(dev_id,&devk);
 if(devn < 0) return -4; //invalid flat device id
 if(drsc->dev_id >= 0 && drsc->dev_id != dev_id) return 1; //resource was assigned to a different device
 if(drsc->gmem_p != NULL) return 2; //resource already has global memory attached
 switch(devk){
  case DEV_HOST:
   if(in_arg_buf == NOPE){
    errc=host_mem_alloc_pin(&(drsc->gmem_p),mem_size); if(errc != 0){drsc->gmem_p = NULL; return 3;}
   }else{
    errc=get_buf_entry_host(mem_size,&byte_ptr,&i);
    if(errc != 0){if(errc == TRY_LATER || errc == DEVICE_UNABLE){return errc;}else{return 4;}}
    drsc->gmem_p=(void*)byte_ptr; drsc->buf_entry=i;
   }
   break;
  case DEV_NVIDIA_GPU:
#ifndef NO_GPU
   if(in_arg_buf == NOPE){
    errc=gpu_mem_alloc(&(drsc->gmem_p),mem_size,devn); if(errc != 0){drsc->gmem_p = NULL; return 5;}
   }else{
    errc=get_buf_entry_gpu(devn,mem_size,&byte_ptr,&i);
    if(errc != 0){if(errc == TRY_LATER || errc == DEVICE_UNABLE){return errc;}else{return 6;}}
    drsc->gmem_p=(void*)byte_ptr; drsc->buf_entry=i;
   }
   break;
#else
   return -5;
#endif
  case DEV_INTEL_MIC:
#ifndef NO_MIC
   //`Future
   break;
#else
   return -6;
#endif
  case DEV_AMD_GPU:
#ifndef NO_AMD
   //`Future
   break;
#else
   return -7;
#endif
  default:
   return -8; //unknown device kind
 }
 drsc->dev_id=dev_id;
 return 0;
}

static int tensDevRsc_free_mem(talsh_dev_rsc_t * drsc)
/** Releases global memory referred to by a device resource descriptor.
    An unsuccessful release of the global memory is marked with
    an error status NOT_CLEAN, but the corresponding components of
    the resource descriptor are cleared anyway. **/
{
 int n,devn,devk,errc;

 n=0;
 if(drsc == NULL) return -1;
 if(drsc->dev_id < 0) return -2;
 if(drsc->gmem_p == NULL) return -3;
 devn=decode_device_id(drsc->dev_id,&devk);
 if(devn < 0) return -4; //invalid flat device id
 switch(devk){
  case DEV_HOST:
   if(drsc->buf_entry >= 0){
    errc=free_buf_entry_host(drsc->buf_entry); if(errc != 0) n=NOT_CLEAN;
    drsc->gmem_p=NULL; drsc->buf_entry=-1;
   }else{
    if(drsc->gmem_p != NULL){
     errc=host_mem_free_pin(drsc->gmem_p); if(errc != 0) n=NOT_CLEAN;
     drsc->gmem_p=NULL;
    }
   }
   break;
  case DEV_NVIDIA_GPU:
#ifndef NO_GPU
   if(drsc->buf_entry >= 0){
    errc=free_buf_entry_gpu(devn,drsc->buf_entry); if(errc != 0) n=NOT_CLEAN;
    drsc->gmem_p=NULL; drsc->buf_entry=-1;
   }else{
    if(drsc->gmem_p != NULL){
     errc=gpu_mem_free(drsc->gmem_p,devn); if(errc != 0) n=NOT_CLEAN;
     drsc->gmem_p=NULL;
    }
   }
   break;
#else
   return -5;
#endif
  case DEV_INTEL_MIC:
#ifndef NO_MIC
   //`Future
   break;
#else
   return -6;
#endif
  case DEV_AMD_GPU:
#ifndef NO_AMD
   //`Future
   break;
#else
   return -7;
#endif
  default:
   return -8; //invalid device kind
 }
 errc=tensDevRsc_empty(drsc);
 return n;
}

#ifndef NO_GPU
static int tensDevRsc_get_const_entry(talsh_dev_rsc_t * drsc, int dev_id)
/** Acquires a constant memory entry on an Nvidia GPU. A return status TRY_LATER or
    DEVICE_UNABLE indicate a temporary or permanent shortage of the corresponding resource. **/
{
 int i,devn,devk,errc;

 if(drsc == NULL) return -1;
 if(dev_id < 0) return -2;
 devn=decode_device_id(dev_id,&devk);
 if(devn < 0) return -3; //invalid flat device id
 if(drsc->dev_id >= 0 && drsc->dev_id != dev_id) return 1; //resource was assigned to a different device
 if(drsc->const_mem_entry >= 0) return 2; //resource already has constant memory entry attached
 if(devk == DEV_NVIDIA_GPU){
  errc=const_args_entry_get(devn,&i);
  if(errc != 0){if(errc == TRY_LATER || errc == DEVICE_UNABLE){return errc;}else{return 3;}}
  drsc->const_mem_entry=i;
 }else{
  return 4;
 }
 drsc->dev_id=dev_id;
 return 0;
}

static int tensDevRsc_release_const_entry(talsh_dev_rsc_t * drsc)
/** Releases a GPU constant memory entry from a given device resource descriptor.
    An unsuccessful release of the GPU constant memory entry is marked with
    an error status NOT_CLEAN, but the corresponding component of the resource
    descriptor is cleared anyway. **/
{
 int n,devn,devk,errc;

 n=0;
 if(drsc == NULL) return -1;
 if(drsc->dev_id < 0) return -2;
 devn=decode_device_id(drsc->dev_id,&devk);
 if(devn < 0) return -3; //invalid flat device id
 if(drsc->const_mem_entry < 0) return 1;
 if(devk == DEV_NVIDIA_GPU){
  errc=const_args_entry_free(devn,drsc->const_mem_entry); if(errc != 0) n=NOT_CLEAN;
  drsc->const_mem_entry=-1;
 }else{
  return 2;
 }
 errc=tensDevRsc_empty(drsc);
 return n;
}
#endif

static int tensDevRsc_release_all(talsh_dev_rsc_t * drsc)
/** Releases all device resources in <drsc>. An unsuccessful release
    of one or more resources is marked with an error status NOT_CLEAN,
    but the corresponding components of the device resource descriptor
    are cleared anyway. **/
{
 int n,devk,devn,errc;

 n=0;
 if(drsc == NULL) return -1;
 if(drsc->dev_id >= 0){ //resource handle is not empty
  devn=decode_device_id(drsc->dev_id,&devk);
  if(devn < 0) return -2; //invalid flat device id
  switch(devk){
   case DEV_HOST:
    if(drsc->buf_entry >= 0){
     errc=free_buf_entry_host(drsc->buf_entry); if(errc != 0) n=NOT_CLEAN;
     drsc->gmem_p=NULL; drsc->buf_entry=-1;
    }else{
     if(drsc->gmem_p != NULL){
      errc=host_mem_free_pin(drsc->gmem_p); i(errc != 0) n=NOT_CLEAN;
      drsc->gmem_p=NULL;
     }
    }
    break;
   case DEV_NVIDIA_GPU:
#ifndef NO_GPU
    if(drsc->buf_entry >= 0){
     errc=free_buf_entry_gpu(devn,drsc->buf_entry); if(errc != 0) n=NOT_CLEAN;
     drsc->gmem_p=NULL; drsc->buf_entry=-1;
    }else{
     if(drsc->gmem_p != NULL){
      errc=gpu_mem_free(drsc->gmem_p,devn); if(errc != 0) n=NOT_CLEAN;
      drsc->gmem_p=NULL;
     }
    }
    if(drsc->const_mem_entry >= 0){
     errc=const_args_entry_free(devn,drsc->const_mem_entry); if(errc != 0) n=NOT_CLEAN;
     drsc->const_mem_entry=-1;
    }
    break;
#else
    return -3;
#endif
   case DEV_INTEL_MIC:
#ifndef NO_MIC
    //`Future
    break;
#else
    return -4;
#endif
   case DEV_AMD_GPU:
#ifndef NO_AMD
    //`Future
    break;
#else
    return -5;
#endif
   default:
    return -6; //invalid device kind
  }
  errc=tensDevRsc_clean(drsc);
 }
 return n;
}

static int tensDevRsc_destroy(talsh_dev_rsc_t * drsc)
/** Completely destroys a device resource descriptor. A return status NOT_CLEAN
    means that certain resources have not been released cleanly,
    but it is not a critical error in general (however, a leak can occur). **/
{
 int n,errc;

 n=0;
 if(drsc == NULL) return -1;
 errc=tensDevRsc_release_all(drsc); if(errc != 0) n=NOT_CLEAN;
 free(drsc);
 return n;
}

static int mi_entry_get(int ** mi_entry)
/** Obtains a pointer to an entry in the multi-index storage slab.
    The entry can fit an <int> multi-index up to MAX_TENSOR_RANK length.
    Returns TRY_LATER if no free handles are currently available. **/
{
 int m;
 *mi_entry=NULL;
 if(miFFE > 0){ //number of free handles left
  m=miFreeHandle[--miFFE];
  *mi_entry=&miBank[m][0];
 }else{
  return TRY_LATER; //currently no free handles left
 }
 return 0;
}

static int mi_entry_release(int * mi_entry)
/** Releases an entry back to the multi-index storage slab. **/
{
 int m;
 if(mi_entry != NULL){
  if(miFFE >= 0){
   m=(int)(mi_entry-&miBank[0][0]);
   if(m%MAX_TENSOR_RANK == 0){
    m/=MAX_TENSOR_RANK;
    miFreeHandle[miFFE++]=m;
   }else{
    return 1;
   }
  }else{
   return 2;
  }
 }else{
  return 3;
 }
 return 0;
}

static int mi_entry_pinned(int * mi_entry)
/** Returns YEP if the multi-index is in the multi-index bank,
    NOPE othewise. **/
{
 int n;

 n=NOPE;
 if(mi_entry != NULL){
  if((unsigned int)(mi_entry-miBank[0][0]) < MAX_GPU_ARGS*MAX_MLNDS_PER_TENS*MAX_TENSOR_RANK) n=YEP;
 }
 return n;
}

#ifndef NO_GPU
__host__ static int cuda_stream_get(int gpu_num, int * cuda_stream_handle)
/** For GPU#gpu_num, returns a usable CUDA stream handle <cuda_stream_handle>.
Non-zero return status means an error, except the return status TRY_LATER means
no free resources are currently available (not an error). **/
{
 *cuda_stream_handle=-1;
 if(gpu_num >= 0 && gpu_num < MAX_GPUS_PER_NODE){
  if(gpu_up[gpu_num] >= GPU_MINE){
   if(CUDAStreamFFE[gpu_num] > 0){ //number of free handles left on GPU#gpu_num
    *cuda_stream_handle=CUDAStreamFreeHandle[gpu_num][--CUDAStreamFFE[gpu_num]];
    if(*cuda_stream_handle < 0 || *cuda_stream_handle >= MAX_CUDA_TASKS){
     *cuda_stream_handle=-1; return 3; //invalid handle: corruption
    }
   }else{
    return TRY_LATER; //all handles are currently busy
   }
  }else{
   return 2;
  }
 }else{
  return 1;
 }
 return 0;
}

__host__ static int cuda_stream_release(int gpu_num, int cuda_stream_handle)
/** For GPU#gpu_num, releases a CUDA stream handle <cuda_stream_handle>.
Non-zero return status means an error. **/
{
 if(gpu_num >= 0 && gpu_num < MAX_GPUS_PER_NODE){
  if(gpu_up[gpu_num] >= GPU_MINE){
   if(cuda_stream_handle >= 0 && cuda_stream_handle < MAX_CUDA_TASKS){
    if(CUDAStreamFFE[gpu_num] < 0 || CUDAStreamFFE[gpu_num] > MAX_CUDA_TASKS) return 5; //corrupted
    if(CUDAStreamFFE[gpu_num] < MAX_CUDA_TASKS){
     CUDAStreamFreeHandle[gpu_num][CUDAStreamFFE[gpu_num]++]=cuda_stream_handle;
    }else{
     return 4; //an attempt to release a non-existing handle
    }
   }else{
    return 3;
   }
  }else{
   return 2;
  }
 }else{
  return 1;
 }
 return 0;
}

__host__ static int cuda_event_get(int gpu_num, int * cuda_event_handle)
/** For GPU#gpu_num, returns a usable CUDA event handle <cuda_event_handle>.
Non-zero return status means an error, except the return status TRY_LATER means
no free resources are currently available (not an error). **/
{
 *cuda_event_handle=-1;
 if(gpu_num >= 0 && gpu_num < MAX_GPUS_PER_NODE){
  if(gpu_up[gpu_num] >= GPU_MINE){
   if(CUDAEventFFE[gpu_num] > 0){ //number of free handles left on GPU#gpu_num
    *cuda_event_handle=CUDAEventFreeHandle[gpu_num][--CUDAEventFFE[gpu_num]];
    if(*cuda_event_handle < 0 || *cuda_event_handle >= MAX_CUDA_EVENTS){
     *cuda_event_handle=-1; return 3; //invalid handle: corruption
    }
   }else{
    return TRY_LATER; //all handles are currently busy
   }
  }else{
   return 2;
  }
 }else{
  return 1;
 }
 return 0;
}

__host__ static int cuda_event_release(int gpu_num, int cuda_event_handle)
/** For GPU#gpu_num, releases a CUDA event handle <cuda_event_handle>.
Non-zero return status means an error. **/
{
 if(gpu_num >= 0 && gpu_num < MAX_GPUS_PER_NODE){
  if(gpu_up[gpu_num] >= GPU_MINE){
   if(cuda_event_handle >= 0 && cuda_event_handle < MAX_CUDA_EVENTS){
    if(CUDAEventFFE[gpu_num] < 0 || CUDAEventFFE[gpu_num] > MAX_CUDA_EVENTS) return 5; //corrupted
    if(CUDAEventFFE[gpu_num] < MAX_CUDA_EVENTS){
     CUDAEventFreeHandle[gpu_num][CUDAEventFFE[gpu_num]++]=cuda_event_handle;
    }else{
     return 4; //an attempt to release a non-existing handle
    }
   }else{
    return 3;
   }
  }else{
   return 2;
  }
 }else{
  return 1;
 }
 return 0;
}

__host__ static void limit_cuda_blocks2d(int max_blocks, int *bx, int *by)
/** Limits the number of CUDA blocks in a 2d grid to <max_blocks>.
    No argument validity check! **/
{
 if(max_blocks > 1){
  double rdc = ((double)max_blocks)/(((double)(*bx))*((double)(*by)));
  if(rdc < 1.0){
   rdc=sqrt(rdc);
   if(*bx > *by){
    *by=(int)(rdc*((double)(*by))); if(*by < 1){*by=1; *bx=max_blocks; return;}
    *bx=(int)(rdc*((double)(*bx)));
   }else{
    *bx=(int)(rdc*((double)(*bx))); if(*bx < 1){*bx=1; *by=max_blocks; return;}
    *by=(int)(rdc*((double)(*by)));
   }
   if((*bx)*(*by) > max_blocks){
    if(*bx > *by){(*bx)--;}else{(*by)--;}
   }
  }
 }else{
  *bx=1; *by=1;
 }
 return;
}

//NV-TAL INITIALIZATION/SHUTDOWN (internal use only):
__host__ int init_gpus(int gpu_beg, int gpu_end)
/** Initializes all GPU contexts for the current MPI process. Returned positive value is
the number of initialized GPUs. A negative return status means an error occured.
Each enabled GPU from the range [gpu_beg:gpu_end] will obtain its own cublasHandle as well.
The first GPU from the given range will be left active at the end. If <gpu_beg> > <gpu_end>,
no GPU will be initialized. **/
{
 size_t m;
 int i,j,n,errc;
 cudaError_t err;
#ifndef NO_BLAS
 cublasStatus_t err_cublas;
#endif
 n=0; for(i=0;i<MAX_GPUS_PER_NODE;i++) gpu_up[i]=GPU_OFF; //initial GPU status
//Multi-index bank:
 miFFE=MAX_GPU_ARGS*MAX_MLNDS_PER_TENS;
 for(j=0;j<miFFE;j++) miFreeHandle[j]=j;
 m=(size_t)(miFFE*MAX_TENSOR_RANK*sizeof(int));
 errc=host_mem_register(&miBank[0][0],m);
 if(errc != 0){
  if(VERBOSE) printf("#ERROR(tensor_algebra_gpu_nvidia:init_gpus): Unable to register the multi-index bank: Error %d\n",errc);
  return -3;
 }
 if(gpu_beg >= 0 && gpu_end >= gpu_beg){
  err=cudaGetDeviceCount(&i); if(err != cudaSuccess) return -1;
  if(gpu_end >= MAX_GPUS_PER_NODE || gpu_end >= i) return -2;
  for(i=gpu_end;i>=gpu_beg;i--){
   err=cudaSetDevice(i);
   if(err == cudaSuccess){
    gpu_up[i]=GPU_MINE; err=cudaGetDeviceProperties(&(gpu_prop[i]),i); if(err != cudaSuccess) gpu_up[i]=GPU_OFF;
    if(gpu_up[i] > GPU_OFF){
//SHMEM width:
     errc=gpu_set_shmem_width(GPU_SHMEM_WIDTH);
     if(errc != 0 && VERBOSE) printf("#WARNING(tensor_algebra_gpu_nvidia:init_gpus): Unable to set GPU SHMEM width %d: Error %d \n",GPU_SHMEM_WIDTH,errc);
//cuBLAS.v2 context:
#ifndef NO_BLAS
     err_cublas=cublasCreate(&(cublas_handle[i]));
     if(err_cublas == CUBLAS_STATUS_SUCCESS){
      gpu_up[i]=GPU_MINE_CUBLAS;
      err_cublas=cublasSetPointerMode(cublas_handle[i],CUBLAS_POINTER_MODE_DEVICE);
      if(err_cublas != CUBLAS_STATUS_SUCCESS) gpu_up[i]=GPU_MINE;
     }
#endif
    }
//CUDA stream bank:
    if(gpu_up[i] > GPU_OFF){
     for(j=0;j<MAX_CUDA_TASKS;j++) CUDAStreamFreeHandle[i][j]=j; CUDAStreamFFE[i]=MAX_CUDA_TASKS;
     for(j=0;j<MAX_CUDA_TASKS;j++){
      err=cudaStreamCreate(&(CUDAStreamBank[i][j])); if(err != cudaSuccess){gpu_up[i]=GPU_OFF; break;};
     }
    }
//CUDA event bank:
    if(gpu_up[i] > GPU_OFF){
     for(j=0;j<MAX_CUDA_EVENTS;j++) CUDAEventFreeHandle[i][j]=j; CUDAEventFFE[i]=MAX_CUDA_EVENTS;
     for(j=0;j<MAX_CUDA_EVENTS;j++){
      err=cudaEventCreate(&(CUDAEventBank[i][j])); if(err != cudaSuccess){gpu_up[i]=GPU_OFF; break;};
     }
    }
//Last task:
    LastTask[i]=NULL;
//Clear GPU statistics:
    gpu_stats[i].tasks_submitted=0;
    gpu_stats[i].tasks_completed=0;
    gpu_stats[i].tasks_deferred=0;
    gpu_stats[i].tasks_failed=0;
    gpu_stats[i].flops=0.0;
    gpu_stats[i].traffic_in=0.0;
    gpu_stats[i].traffic_out=0.0;
    gpu_stats[i].time_active=0.0;
    gpu_stats[i].time_start=clock();
//Accept GPU as ready (active):
    if(gpu_up[i] > GPU_OFF) n++;
   }
  }
 }
 return n; //number of initialized GPU's
}

__host__ int free_gpus(int gpu_beg, int gpu_end)
/** Destroys all GPU/CUBLAS contexts on all GPU devices belonging to the MPI process.
A positive value returned is the number of failed GPUs; a negative one is an error.
If <gpu_beg> > <gpu_end>, nothing wil be done. **/
{
 int i,j,n,failure;
 cudaError_t err;
#ifndef NO_BLAS
 cublasStatus_t err_cublas;
#endif
 failure=0; n=0;
//Multi-index bank:
 miFFE=MAX_GPU_ARGS*MAX_MLNDS_PER_TENS;
 for(j=0;j<miFFE;j++) miFreeHandle[j]=j;
 j=host_mem_unregister(&miBank[0][0]); //`This is probably not needed
 if(j != 0){
  failure++;
  if(VERBOSE) printf("#WARNING(tensor_algebra_gpu_nvidia:free_gpus): Unable to unregister the multi-index bank: Error %d\n",j);
 }
 if(gpu_beg >= 0 && gpu_end >= gpu_beg){
  err=cudaGetDeviceCount(&i); if(err != cudaSuccess) return -1;
  if(gpu_end >= MAX_GPUS_PER_NODE || gpu_end >= i) return -2;
  for(i=gpu_beg;i<=gpu_end;i++){
   if(gpu_up[i] > GPU_OFF){
    n++; err=cudaSetDevice(i);
    if(err == cudaSuccess){
#ifndef NO_BLAS
     if(gpu_up[i] >= GPU_MINE_CUBLAS){err_cublas=cublasDestroy(cublas_handle[i]); if(err_cublas == CUBLAS_STATUS_SUCCESS) gpu_up[i]=GPU_MINE;}
#endif
//CUDA stream bank:
     if(gpu_up[i] > GPU_OFF){
      for(j=0;j<MAX_CUDA_TASKS;j++) CUDAStreamFreeHandle[i][j]=j; CUDAStreamFFE[i]=MAX_CUDA_TASKS;
      for(j=0;j<MAX_CUDA_TASKS;j++){err=cudaStreamDestroy(CUDAStreamBank[i][j]); if(err != cudaSuccess) failure++;}
     }
//CUDA event bank:
     if(gpu_up[i] > GPU_OFF){
      for(j=0;j<MAX_CUDA_EVENTS;j++) CUDAEventFreeHandle[i][j]=j; CUDAEventFFE[i]=MAX_CUDA_EVENTS;
      for(j=0;j<MAX_CUDA_EVENTS;j++){err=cudaEventDestroy(CUDAEventBank[i][j]); if(err != cudaSuccess) failure++;}
     }
//Last task:
     LastTask[i]=NULL;
     n--; err=cudaDeviceReset();
    }
    gpu_up[i]=GPU_OFF; //GPU is taken out of use regardless of its status!
   }
  }
 }
 if(failure && VERBOSE) printf("#WARNING(tensor_algebra_gpu_nvidia:free_gpus): Resource deallocation was not fully successful!");
 return n;
}

__host__ int gpu_is_mine(int gpu_num)
/** Positive return: GPU is mine; 0: GPU is not mine; -1: invalid <gpu_num>. **/
{if(gpu_num >= 0 && gpu_num < MAX_GPUS_PER_NODE){return gpu_up[gpu_num];}else{return -1;}}

__host__ int gpu_busy_least()
/** Returns the ID of the least busy GPU (non-negative) or -1 (no GPU found). **/
{
 int i,j,m,n;
 m=-1; n=-1;
 for(i=0;i<MAX_GPUS_PER_NODE;i++){
  if(gpu_up[i] != GPU_OFF){
   j=gpu_stats[i].tasks_submitted-(gpu_stats[i].tasks_completed+gpu_stats[i].tasks_deferred+gpu_stats[i].tasks_failed);
   if(m >= 0){
    if(j < m){m=j; n=i;};
   }else{
    m=j; n=i;
   }
  }
 }
 return n;
}

__host__ int gpu_activate(int gpu_num)
/** If GPU is enabled (mine), does cudaSetDevice; returns non-zero otherwise (error). **/
{
 cudaError_t err;
 if(gpu_num >= 0 && gpu_num < MAX_GPUS_PER_NODE){
  if(gpu_up[gpu_num] > GPU_OFF){err=cudaSetDevice(gpu_num); if(err != cudaSuccess) return 3;}else{return 2;}
 }else{
  return 1; //invalid <gpu_num>
 }
 return 0;
}

//NV-TAL INTERNAL CONTROL:
__host__ int gpu_set_shmem_width(int width){
/** Sets the GPU shared memory bank width:
    <width> = R4: 4 bytes;
    <width> = R8: 8 bytes. **/
 cudaError_t cerr;
 if(width == R8){
  cerr=cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeEightByte);
 }else if(width == R4){
  cerr=cudaDeviceSetSharedMemConfig(cudaSharedMemBankSizeFourByte);
 }else{
  return 1; //invalid <width> passed
 }
 if(cerr != cudaSuccess) return 2;
 return 0;
}

__host__ void gpu_set_event_policy(int alg)
/** Turns on/off timing CUDA events (1/0). **/
{if(alg == EVENTS_OFF){EVENT_RECORD=EVENTS_OFF;}else{EVENT_RECORD=EVENTS_ON;}; return;}

__host__ void gpu_set_transpose_algorithm(int alg)
/** Activates either the scatter or the shared-memory based tensor transpose algorithm. **/
{if(alg == EFF_TRN_OFF){TRANS_SHMEM=EFF_TRN_OFF;}else{TRANS_SHMEM=EFF_TRN_ON;}; return;}

__host__ void gpu_set_matmult_algorithm(int alg){
/** Activates either cuBLAS (fast) or my own (slow) BLAS CUDA kernels. **/
#ifndef NO_BLAS
 if(alg == BLAS_ON){DISABLE_BLAS=BLAS_ON;}else{DISABLE_BLAS=BLAS_OFF;};
#endif
 return;
}

__host__ int gpu_print_stats(int gpu_num = -1)
/** Prints GPU statistics for GPU#<gpu_num>. If <gpu_num>=-1,
    prints GPU statistics for all active GPUs.
    A negative return status means invalid <gpu_num>. **/
{
 int i,b,f;
 clock_t ctm;

 if(gpu_num >= 0 && gpu_num < MAX_GPUS_PER_NODE){
  b=gpu_num; f=gpu_num; //select a specific GPU
 }else if(gpu_num == -1){
  b=0; f=MAX_GPUS_PER_NODE-1; //select all GPUs
 }else{
  return -1; //invalid GPU number
 }
 for(i=b;i<=f;i++){
  if(gpu_is_mine(i) != GPU_OFF){
   ctm=clock();
   gpu_stats[i].time_active=((double)(ctm-gpu_stats[i].time_start))/CLOCKS_PER_SEC;
   printf("\n#MSG(TAL-SH::NV-TAL): Statistics on GPU #%d:\n",i);
   printf(" Number of tasks submitted: %llu\n",gpu_stats[i].tasks_submitted);
   printf(" Number of tasks completed: %llu\n",gpu_stats[i].tasks_completed);
   printf(" Number of tasks deferred : %llu\n",gpu_stats[i].tasks_deferred);
   printf(" Number of tasks failed   : %llu\n",gpu_stats[i].tasks_failed);
   printf(" Number of Flops processed: %G\n",gpu_stats[i].flops);
   printf(" Number of Bytes received : %G\n",gpu_stats[i].traffic_in);
   printf(" Number of Bytes sent     : %G\n",gpu_stats[i].traffic_out);
   printf(" Time active (sec)        : %f\n",gpu_stats[i].time_active);
  }else{
   printf("\n#MSG(TAL-SH::NV-TAL): Statistics on GPU #%d: GPU is OFF\n",i);
  }
 }
 return 0;
}
#endif

//TENSOR BLOCK API:
int tensShape_clean(talsh_tens_shape_t * tshape)
/** Cleans a tensor shape. A clean (initialized to null) tensor shape has .num_dim=-1.
    A further defined tensor shape has .num_dim >= 0. **/
{
 if(tshape != NULL){
  tshape->num_dim=-1; //tensor rank
  tshape->dims=NULL;  //tensor dimension extents
  tshape->divs=NULL;  //tensor dimension dividers (segment sizes)
  tshape->grps=NULL;  //tensor dimension groups
 }else{
  return -1;
 }
 return 0;
}

int tensShape_construct(talsh_tens_shape_t * tshape, int pinned, int rank, const int * dims = NULL,
                                                  const int * divs = NULL, const int * grps = NULL)
/** (Re-)defines a tensor shape. It is errorneous to pass an uninitialized tensor shape here,
    that is, the tensor shape *(tshape) must be either clean or previously defined. If <rank> > 0,
    <dims[rank]> must be supplied, whereas <divs[rank]> and <grps[rank]> are always optional.
    If <pinned> = YEP, then the multi-indices will be allocated via the multi-index bank (pinned),
    otherwise a regular malloc will be called. TRY_LATER or DEVICE_UNABLE return statuses are not
    errors and in this case the input tensor shape will stay unchanged. A return status NOT_CLEAN
    indicates an unsuccessful resource release that can be tolerated in general
    (the construction will still occur). **/
{
 int i,errc;
 int *mi_dims,*mi_divs,*mi_grps;

 errc=0;
//Check arguments:
 if(tshape == NULL) return -1;
 if(rank < 0) return -2;
 if(dims != NULL){for(i=0;i<rank;i++){if(dims[i] < 0) return -3;}}
 if(divs != NULL){for(i=0;i<rank;i++){if(divs[i] < 0) return -4;}}
 if(grps != NULL){for(i=0;i<rank;i++){if(grps[i] < 0) return -5;}}
 if(rank > 0 && dims == NULL) return -6; //dimension extents must be present for rank>0
//Acquire/release resources if needed:
 mi_dims=NULL; mi_divs=NULL; mi_grps=NULL;
 if(rank > 0 && tshape->num_dim <= 0){ //acquire multi-index resources
  if(pinned == NOPE){
   mi_dims=(int*)malloc(3*MAX_TENSOR_RANK*sizeof(int));
   if(mi_dims == NULL) return TRY_LATER;
   mi_divs=mi_dims+MAX_TENSOR_RANK;
   mi_grps=mi_divs+MAX_TENSOR_RANK;
  }else{
 //Multi-index "Dimension extents":
   errc=mi_entry_get(&mi_dims); //acquire a mi resource
   if(errc != 0){
    if(errc == TRY_LATER || errc == DEVICE_UNABLE){return errc;}else{return 1;}
   }
 //Multi-index "Dimension dividers":
   errc=mi_entry_get(&mi_divs); //acquire a mi resource
   if(errc != 0){
    i=mi_entry_release(mi_dims);
    if(errc == TRY_LATER || errc == DEVICE_UNABLE){return errc;}else{return 2;}
   }
 //Multi-index "Dimension groups":
   errc=mi_entry_get(&mi_grps); //acquire a mi resource
   if(errc != 0){
    i=mi_entry_release(mi_divs); i=mi_entry_release(mi_dims);
    if(errc == TRY_LATER || errc == DEVICE_UNABLE){return errc;}else{return 3;}
   }
  }
  tshape->dims=mi_dims; tshape->divs=mi_divs; tshape->grps=mi_grps;
  errc=0;
 }else if(rank == 0 && tshape->num_dim > 0){ //release multi-index resources
  errc=tensShape_destruct(tshape); if(errc != 0 && errc != NOT_CLEAN) return 4;
 }
//Define the new tensor shape:
 tshape->num_dim=rank;
 if(dims != NULL){
  for(i=0;i<rank;i++) tshape->dims[i]=dims[i];
 }
 if(divs != NULL){
  for(i=0;i<rank;i++) tshape->divs[i]=divs[i];
 }else{
  for(i=0;i<rank;i++) tshape->divs[i]=tshape->dims[i]; //default dividers (one segment per dimension)
 }
 if(grps != NULL){
  for(i=0;i<rank;i++) tshape->grps[i]=grps[i];
 }else{
  for(i=0;i<rank;i++) tshape->grps[i]=0; //default groups (all indices belong to the unrestricted group)
 }
 return errc; //either 0 or NOT_CLEAN
}

int tensShape_destruct(talsh_tens_shape_t * tshape)
/** Destructs a defined tensor shape (releases resources and cleans it).
    If the input tensor shape is initialized to null, nothing happens.
    In case of an unsuccessful resource release, a return status NOT_CLEAN
    will be returned, which can be considered as a tolerable error since
    the tensor shape will be cleaned anyway (although a leak can occur). **/
{
 int n,pinned,errc;

 n=0; //will be incremented upon an unsuccessful resource release
 if(tshape == NULL) return -1;
 if(tshape->num_dim > 0){ //need to release resources
  if(tshape->dims != NULL){
   pinned=mi_entry_pinned(tshape->dims);
   if(pinned == NOPE){
    free(tshape->dims); //will free all {dims,divs,grps}
    tshape->dims=NULL; tshape->divs=NULL; tshape->grps=NULL;
   }else{
    if(tshape->dims != NULL){errc=mi_entry_release(tshape->dims); if(errc != 0) n++; tshape->dims=NULL;} //release a mi resource
    if(tshape->divs != NULL){errc=mi_entry_release(tshape->divs); if(errc != 0) n++; tshape->divs=NULL;} //release a mi resource
    if(tshape->grps != NULL){errc=mi_entry_release(tshape->grps); if(errc != 0) n++; tshape->grps=NULL;} //release a mi resource
   }
  }else{
   return -2;
  }
 }
 if(n != 0) n=NOT_CLEAN;
 errc=tensShape_clean(tshape);
 return n; //either 0 or NOT_CLEAN
}

static size_t tensShape_volume(const talsh_tens_shape_t * tshape)
/** Returns the volume of a defined tensor shape, or 0 otherwise. **/
{
 int i;
 size_t vol;

 vol=0;
 if(tshape->num_dim >= 0 && tshape->num_dim <= MAX_TENSOR_RANK){
  vol=1;
  for(i=0;i<tshape->num_dim;i++){
   if(tshape->dims[i] > 0){
    vol*=tshape->dims[i];
   }else{
    return 0;
   }
  }
 }
 return vol;
}

int tensBlck_create(tensBlck_t **ctens)
/** Creates an empty instance of tensBlck_t and initializes it to null (on Host). **/
{
 int errc;

 *ctens=(tensBlck_t*)malloc(sizeof(tensBlck_t)); if(*ctens == NULL) return TRY_LATER;
 (*ctens)->data_kind=NO_TYPE;
 errc=tensShape_clean(&((*ctens)->shape)); if(errc != 0) return 1;
 (*ctens)->src_rsc=NULL; //source memory resource (where the tensor body is before the operation)
 (*ctens)->dst_rsc=NULL; //destination memory resource (where the tensor body will be after the operation)
 (*ctens)->tmp_rsc=NULL; //temporary memory resource (where the tensor body can be during the operation)
 (*ctens)->prmn_h=NULL;  //dimension permutation (internal use only)
 return 0;
}

int tensBlck_destroy(tensBlck_t *ctens)
/** Destroys a defined instance of tensBlck_t (either nullified or shape-defined).
    A return status NOT_CLEAN indicates an unsuccessful resource release, which
    can be considered as a tolerable error (the object will still be destroyed). **/
{
 int errc;

 errc=0;
 if(ctens == NULL) return -1;
 errc=tensBlck_destruct(ctens); if(errc != 0 && errc != NOT_CLEAN) errc=1;
 free(ctens);
 return errc;
}

int tensBlck_construct(tensBlck_t *ctens,      //pointer to defined tensor block (either nullified or defined to a value)
                       int pinned,             //YEP: tensor shape multi-indices will be pinned (for GPU), NOPE: regular malloc (not pinned)
                       int trank,              //tensor rank
                       const int *dims = NULL, //tensor dimension extents (when trank > 0)
                       const int *divs = NULL, //tensor dimension dividers (when trank > 0, optional)
                       const int *grps = NULL) //tensor dimension groups (when trank > 0, optional)
/** Constructs (defines/redefines) a tensor block without attaching its body (only the shape).
    If the tensor block is to be used on Nvidia GPUs or other asynchronous devices,
    argument <pinned> must be set to YEP (NOPE will not use pinned memory).
    A return status NOT_CLEAN indicates an unsuccessful resource release, which,
    can be considered as a tolerable error (the object will still be constructed). **/
{
 int n,errc;

 n=0;
 if(ctens == NULL) return -1;
 if(trank < 0 || trank > MAX_TENSOR_RANK) return -2; //invalid tensor rank
 if(trank > 0 && dims == NULL) return -3; //dimension extents must be present for rank>0 tensors
 errc=tensBlck_destruct(ctens); if(errc != 0){if(errc == NOT_CLEAN){n=errc;}else{return 1;}}
 errc=tensShape_construct(&(ctens->shape),pinned,trank,dims,divs,grps);
 if(errc != 0){if(errc == TRY_LATER || errc == DEVICE_UNABLE){return errc;}else{return 2;}}
 return n; //either 0 or NOT_CLEAN
}

int tensBlck_attach_body(tensBlck_t *ctens,     //pointer to a shape-defined (constructed) tensor block
                         int data_kind,         //data kind (R4,R8,C8)
                         int dev_id,            //flat device id where the body resides (or should reside)
                         void *body_ptr = NULL, //pointer to the tensor body (global memory of device <dev_id>)
                         int buf_entry = -1)    //argument buffer entry handle corresponding to the <body_ptr> (optional)
/** Attaches a body to a shape-defined tensor block. If both <body_ptr> and <buf_entry> are absent,
    a resource will be allocated on device <dev_id> in the device argument buffer (if available).
    If <buf_entry> is absent, a defined <body_ptr> points to an external memory (either pinned or not).
    If both <body_ptr> and <buf_entry> are defined, the external memory is assumed to be within that
    argument buffer entry. In all cases, the memory resource will be associated with the .src_rsc component
    of tensBlck_t. It is forbidden to attempt allocating/attaching a memory resource when an existing memory
    resource is still in use (this will result in an error). A return status of TRY_LATER or DEVICE_UNABLE
    indicates the current or permanent shortage in the necessary resources and is not an error. **/
{
 int errc;
 size_t vol,body_size;

 if(ctens == NULL) return -1;
 if(tens_valid_data_kind(data_kind) != YEP || data_kind == NO_TYPE) return -2;
 if(ctens->shape.num_dim < 0 || ctens->shape.num_dim > MAX_TENSOR_RANK) return -3; //tensor block must be shape-defined
 if(body_ptr == NULL && buf_entry >= 0) return -4; //a defined argument buffer entry must be supplied with the corresponding pointer
 vol=tensShape_volume(&(ctens->shape)); //tensor body volume (number of elements)
 body_size=vol*data_kind; //tensor body size in bytes
 if(body_ptr == NULL){ //allocate memory in the argument buffer
  errc=tensDevRsc_allocate_mem(ctens->src_rsc,dev_id,body_size,YEP);
  if(errc != 0){if(errc == TRY_LATER || errc == DEVICE_UNABLE){return errc;}else{return 1;}}
 }else{ //associate memory
  errc=tensDevRsc_attach_mem(ctens->src_rsc,dev_id,body_ptr,buf_entry);
  if(errc != 0){if(errc == TRY_LATER || errc == DEVICE_UNABLE){return errc;}else{return 2;}}
 }
 return 0;
}

int tensBlck_destruct(tensBlck_t *ctens, int release_body = YEP, int which_body = EVERYTHING)
/** Destructs a defined tensor block (releases all resources and initializes the tensor block to null).
    If <release_body> == YEP/NOPE, the global memory resources will be released/kept. Argument <which_body>
    can further regulate which tensor body to be released/kept (SOURCE, DESTINATION, TEMPORARY, EVERYTHING).
    A return status NOT_CLEAN indicates an unsuccessful resource release that may be
    considered as a tolerable error since the tensor block will be nullified anyway. **/
{
 int n,errc;

 n=0;
 if(ctens == NULL) return -1;
 if(ctens->shape.num_dim >= 0){ //shape-defined tensor block
  if(ctens->shape.num_dim > MAX_TENSOR_RANK) return -2;
//Dimension permutation (temporary resource used by Nvidia GPU):
  if(ctens->prmn_h != NULL){errc=mi_entry_release(ctens->prmn_h); if(errc != 0) n=NOT_CLEAN; ctens->prmn_h=NULL;} //release a mi resource
//Release the TEMPORARY resource:
  if(ctens->tmp_rsc != NULL &&
     ((release_body == YEP && (which_body == EVERYTHING || which_body == TEMPORARY)) ||
      (release_body == NOPE && (which_body != EVERYTHING && which_body != TEMPORARY)))){
   errc=tensDevRsc_release(ctens->tmp_rsc); if(errc != 0) n=NOT_CLEAN;
  }
  ctens->tmp_rsc=NULL;
//Release the DESTINATION resource:
  if(ctens->dst_rsc != NULL &&
     ((release_body == YEP && (which_body == EVERYTHING || which_body == DESTINATION)) ||
      (release_body == NOPE && (which_body != EVERYTHING && which_body != DESTINATION)))){
   errc=tensDevRsc_release(ctens->dst_rsc); if(errc != 0) n=NOT_CLEAN;
  }
  ctens->dst_rsc=NULL;
//Release the SOURCE resource:
  if(ctens->src_rsc != NULL &&
     ((release_body == YEP && (which_body == EVERYTHING || which_body == SOURCE)) ||
      (release_body == NOPE && (which_body != EVERYTHING && which_body != SOURCE)))){
   errc=tensDevRsc_release(ctens->src_rsc); if(errc != 0) n=NOT_CLEAN;
  }
  ctens->src_rsc=NULL;
 }else{ //nullified tensor block: All resources must have been released already
  if(ctens->src_rsc != NULL){ctens->src_rsc=NULL; n=NOT_CLEAN;}
  if(ctens->dst_rsc != NULL){ctens->dst_rsc=NULL; n=NOT_CLEAN;}
  if(ctens->tmp_rsc != NULL){ctens->tmp_rsc=NULL; n=NOT_CLEAN;}
  if(ctens->prmn_h != NULL){ctens->prmn_h=NULL; n=NOT_CLEAN;}
 }
 if(tens_valid_data_kind(ctens->data_kind) != YEP) n=NOT_CLEAN; ctens->data_kind=NO_TYPE;
 errc=tensShape_destruct(&(ctens->shape)); if(errc != 0){if(errc == NOT_CLEAN){n=NOT_CLEAN;}else{return 1;}}
 return n;
}

int tensBlck_acc_id(const tensBlck_t *ctens, int *dev_kind, int *entry_gpu, int *entry_const, int *data_kind,
                             int *there)
/** Returns Accelerator ID on which the tensor block data resides or will reside (negative return means Host residence).
INPUT:
 # ctens - pointer to an instance of tensBlck_t;
OUTPUT:
 # tensBlck_acc_id - ACC ID (non-negative) OR -1 (Host residence of the tensor block);
 # dev_kind - device (accelerator) kind (GPU, MIC, etc.);
 # entry_gpu - device global memory argument-buffer entry number;
 # entry_const - GPU constant memory entry number;
 # data_kind - 4:float, 8:double, 16:double_complex;
 # there - 0 means that tensor block elements are not in the ACC memory yet, 1 means opposite.
**/
{
 int dev_num;
 if(ctens != NULL){
  *data_kind=ctens->data_kind; *there=tensBlck_present(ctens);
  *dev_kind=DEV_HOST; dev_num=decode_device_id(ctens->device_id,dev_kind);
  if(*dev_kind == DEV_NVIDIA_GPU){ //GPU residence
   *entry_gpu=ctens->buf_entry_gpu; *entry_const=ctens->const_args_entry;
  }else{ //Not GPU residence
   *entry_gpu=-1; *entry_const=-1;
  }
 }else{
  dev_num=-1; *dev_kind=-1; *entry_gpu=-1; *entry_const=-1; *data_kind=0; *there=0;
 }
 return dev_num;
}

int tensBlck_set_presence(tensBlck_t *ctens) //Marks tensor block data as residing on GPU
{if(ctens != NULL){if(ctens->device_id < 0) ctens->device_id=-(ctens->device_id); return 0;}else{return 1;}}

int tensBlck_set_absence(tensBlck_t *ctens) //Unmarks tensor block data as residing on GPU
{if(ctens != NULL){if(ctens->device_id > 0) ctens->device_id=-(ctens->device_id); return 0;}else{return 1;}}

int tensBlck_present(const tensBlck_t *ctens) //Checks presence of the block data on Device (or Host)
{if(ctens != NULL){if(ctens->device_id >= 0){return 1;}else{return 0;}}else{return -1;}}

int tensBlck_hab_free(tensBlck_t *ctens){
/** For tensor blocks simultaneously residing on Host and GPU, frees the Host copy.
The data does not have to be present on GPU, in which case the tensor block
becomes uninitialized but still usable on the GPU. If the Host-residing tensor body
had been allocated in HAB, it frees that HAB buffer entry, otherwise deallocates
the corresponding pinned memory explicitly via <host_mem_free_pin>. **/
 int i,dev_kind,dev_num,errc;

 errc=0;
 if(ctens != NULL){
  dev_num=decode_device_id(ctens->device_id,&dev_kind);
  if(dev_num >= 0){
   if(dev_kind == DEV_NVIDIA_GPU){
    if(ctens->elems_h != NULL){
     if(ctens->buf_entry_host >= 0){
      i=free_buf_entry_host(ctens->buf_entry_host); errc+=i;
      if(i == 0) ctens->buf_entry_host=-1;
     }else{
      i=host_mem_free_pin(ctens->elems_h); errc+=i;
     }
     ctens->buf_entry_host=-1; ctens->elems_h=NULL;
    }else{
     errc+=333;
    }
   }else{
    errc+=555;
   }
  }else{
   errc+=777;
  }
 }else{
  errc+=999;
 }
 return errc;
}

size_t tensBlck_volume(const tensBlck_t *ctens) //Number of elements in a tensor block (volume)
{size_t tvol=1; for(int i=0;i<ctens->rank;i++){tvol*=(ctens->dims_h[i]);}; return tvol;}

#ifndef NO_GPU
//CUDA TASK API:
__host__ int cuda_task_create(cudaTask_t **cuda_task)
/** Creates an instance of cudaTask_t **/
{
//printf("\n#DEBUG(tensor_algebra_gpu_nvidia:cuda_task_create): sizeof(cudaTask_t) = %d",sizeof(cudaTask_t)); //debug
 *cuda_task=(cudaTask_t*)malloc(sizeof(cudaTask_t)); if(*cuda_task == NULL) return 1;
 (*cuda_task)->task_error=-1; (*cuda_task)->gpu_id=-1; (*cuda_task)->scr_entry_count=0;
 return 0;
}

__host__ int cuda_task_clean(cudaTask_t *cuda_task)
/** Cleans an existing cudaTask_t for reuse. **/
{
 int i,j,cur_gpu,err_code;
 cudaError_t err;

 err_code=0;
 if(cuda_task != NULL){
  i=cuda_task_complete(cuda_task);
  if(i == CUDA_TASK_COMPLETED){ //task has completed (successfully or not)
   if(cuda_task->gpu_id >= 0 && cuda_task->gpu_id < MAX_GPUS_PER_NODE){
    cur_gpu=-1; err=cudaGetDevice(&cur_gpu); if(err != cudaSuccess){cur_gpu=-1; err_code+=1;}
    err=cudaSuccess; if(cur_gpu != cuda_task->gpu_id) err=cudaSetDevice(cuda_task->gpu_id);
    if(err == cudaSuccess){
     err=cudaEventDestroy(cuda_task->task_finish);
     err=cudaEventDestroy(cuda_task->task_output);
     err=cudaEventDestroy(cuda_task->task_comput);
     err=cudaEventDestroy(cuda_task->task_start);
     err=cudaStreamDestroy(cuda_task->task_stream); if(err != cudaSuccess) err_code+=5;
    }else{
     err_code+=20;
    }
    if(cuda_task->scr_entry_count <= MAX_SCR_ENTRY_COUNT){
     for(i=cuda_task->scr_entry_count-1;i>=0;i--){
      j=free_buf_entry_gpu(cuda_task->gpu_id,cuda_task->scr_entry[i]); err_code+=j*100;
     }
    }else{
     err_code+=10000;
    }
    if(cur_gpu >= 0 && cur_gpu != cuda_task->gpu_id){err=cudaSetDevice(cur_gpu); if(err != cudaSuccess) err_code+=30000;}
   }else{
    err_code+=50000;
   }
   cuda_task->task_error=-1; cuda_task->gpu_id=-1; cuda_task->scr_entry_count=0;
  }else if(i == CUDA_TASK_EMPTY || i == CUDA_TASK_ERROR){ //empty task or cuda_task_complete() failed
   cuda_task->task_error=-1; cuda_task->gpu_id=-1; cuda_task->scr_entry_count=0;
  }else{ //task has not completed yet, thus cannot be destroyed
   err_code+=100000;
  }
 }else{
  err_code+=300000;
 }
 return err_code;
}

__host__ int cuda_task_destroy(cudaTask_t *cuda_task)
/** Destroys an instance of cudaTask_t if the CUDA task has completed or empty. **/
{
 int i,j,cur_gpu,err_code;
 cudaError_t err;

 err_code=0;
 if(cuda_task != NULL){
  if(EVENT_RECORD != 0){ //`This should be done atomically for thread safety
   if(cuda_task->gpu_id >= 0 && cuda_task->gpu_id < MAX_GPUS_PER_NODE){
    if(LastTask[cuda_task->gpu_id] == cuda_task) LastTask[cuda_task->gpu_id]=NULL;
   }
  }
  i=cuda_task_complete(cuda_task);
  if(i == CUDA_TASK_COMPLETED){ //task has completed (successfully or not)
   if(cuda_task->gpu_id >= 0 && cuda_task->gpu_id < MAX_GPUS_PER_NODE){
    cur_gpu=-1; err=cudaGetDevice(&cur_gpu); if(err != cudaSuccess){cur_gpu=-1; err_code+=1;}
    err=cudaSuccess; if(cur_gpu != cuda_task->gpu_id) err=cudaSetDevice(cuda_task->gpu_id);
    if(err == cudaSuccess){
     err=cudaEventDestroy(cuda_task->task_finish);
     err=cudaEventDestroy(cuda_task->task_output);
     err=cudaEventDestroy(cuda_task->task_comput);
     err=cudaEventDestroy(cuda_task->task_start);
     err=cudaStreamDestroy(cuda_task->task_stream); if(err != cudaSuccess) err_code+=5;
    }else{
     err_code+=20;
    }
    if(cuda_task->scr_entry_count <= MAX_SCR_ENTRY_COUNT){
     for(i=cuda_task->scr_entry_count-1;i>=0;i--){
      j=free_buf_entry_gpu(cuda_task->gpu_id,cuda_task->scr_entry[i]); err_code+=j*100;
     }
    }else{
     err_code+=10000;
    }
    if(cur_gpu >= 0 && cur_gpu != cuda_task->gpu_id){err=cudaSetDevice(cur_gpu); if(err != cudaSuccess) err_code+=30000;}
   }else{
    err_code+=50000;
   }
   free(cuda_task); cuda_task=NULL;
  }else if(i == CUDA_TASK_EMPTY || i == CUDA_TASK_ERROR){ //empty task or cuda_task_complete() failed
   free(cuda_task); cuda_task=NULL;
  }else{ //task has not completed yet, thus cannot be destroyed
   err_code+=100000;
  }
 }else{
  err_code+=300000;
 }
 return err_code;
}

__host__ static int cuda_task_finalize(cudaTask_t *cuda_task, int err_code, int gpu_num=-1)
/** Finalizes a CUDA task: gpu_num=-1: on Host; gpu_num>=0: on GPU#gpu_num. **/
{if(cuda_task != NULL){cuda_task->task_error=err_code; cuda_task->gpu_id=gpu_num; return 0;}else{return 1;}}

__host__ static int cuda_task_record(cudaTask_t *cuda_task, int err_code, int gpu_num, cudaStream_t cuda_stream,
                     cudaEvent_t cuda_start, cudaEvent_t cuda_comput, cudaEvent_t cuda_output, cudaEvent_t cuda_finish,
                     int scr_entry_cnt, int *scr_entries)
/** Registers a CUDA task: Launch-error-free tasks are recorded with .task_error=-1 (in normal progress). **/
{
 int i;
 if(cuda_task != NULL){
  if(err_code == 0){ //No error occured during the task scheduling
   cuda_task->task_error=-1;            //error code (<0: In progress; 0: Success; >0: Launch error (may be in progress))
  }else{
   cuda_task->task_error=err_code;      //error code (<0: In progress; 0: Success; >0: Launch error (may be in progress))
  }
  if(gpu_num >= 0 && gpu_num < MAX_GPUS_PER_NODE){
   cuda_task->gpu_id=gpu_num;            //GPU number on which the task was scheduled
   cuda_task->task_stream=cuda_stream;   //CUDA stream assinged to the task
   cuda_task->task_start=cuda_start;     //CUDA event recorded at the beginning of the task
   cuda_task->task_comput=cuda_comput;   //CUDA event recorded when the computing kernel starts (all the input data is on device)
   cuda_task->task_output=cuda_output;   //CUDA event recorded when the computing kernel finishes (before output is copied back)
   cuda_task->task_finish=cuda_finish;   //CUDA event recorded at the end of the task
   if(scr_entry_cnt >= 0 && scr_entry_cnt <= MAX_SCR_ENTRY_COUNT){
    cuda_task->scr_entry_count=scr_entry_cnt; //number of additional GPU argument buffer entries allocated by the task
    for(i=0;i<scr_entry_cnt;i++) cuda_task->scr_entry[i]=scr_entries[i]; //additional GPU argument buffer entries allocated by the task
   }else{
    return 3;
   }
  }else{
   return 2;
  }
 }else{
  return 1;
 }
 return 0;
}

__host__ int cuda_task_gpu_id(const cudaTask_t *cuda_task){return cuda_task->gpu_id;}

__host__ int cuda_task_status(cudaTask_t *cuda_task)
/** Checks the status of a CUDA task. Possible status values are listed in
tensor_algebra.h and tensor_algebra.inc: Keep them consistent! **/
{
 int task_stat,cur_gpu;
 cudaError_t err;
 if(cuda_task != NULL){
  if(cuda_task->task_error == 0) return CUDA_TASK_COMPLETED; //task has successfully completed
  cur_gpu=-1; err=cudaGetDevice(&cur_gpu); if(err != cudaSuccess) return CUDA_TASK_ERROR;
  if(cur_gpu != cuda_task->gpu_id){
   err=cudaSetDevice(cuda_task->gpu_id); if(err != cudaSuccess){err=cudaSetDevice(cur_gpu); return CUDA_TASK_ERROR;}
  }
  err=cudaEventQuery(cuda_task->task_finish);
  if(err == cudaSuccess){
   if(cuda_task->task_error < 0) cuda_task->task_error=0; //task completed successfully
   task_stat=CUDA_TASK_COMPLETED; //task completed
  }else{
   err=cudaEventQuery(cuda_task->task_output);
   if(err == cudaSuccess){
    task_stat=CUDA_TASK_OUTPUT_THERE; //computing kernel has finished
   }else{
    err=cudaEventQuery(cuda_task->task_comput);
    if(err == cudaSuccess){
     task_stat=CUDA_TASK_INPUT_THERE; //computation started, input data is on device (can be reused later)
    }else{
     err=cudaEventQuery(cuda_task->task_start);
     if(err == cudaSuccess){
      task_stat=CUDA_TASK_STARTED; //task started
     }else{
      task_stat=CUDA_TASK_SCHEDULED; //task has not started yet
     }
    }
   }
  }
  if(cur_gpu != cuda_task->gpu_id){err=cudaSetDevice(cur_gpu); if(err != cudaSuccess) return CUDA_TASK_ERROR;}
 }else{
  task_stat=CUDA_TASK_EMPTY;
 }
 return task_stat;
}

__host__ int cuda_task_complete(cudaTask_t *cuda_task)
/** Returns CUDA_TASK_COMPLETED if an existing CUDA task <cuda_task> has completed.
Note that having cuda_task->task_error=0 suggests completion without further querying!
Other possible outputs: CUDA_TASK_EMPTY, CUDA_TASK_SCHEDULED, CUDA_TASK_COMPLETED, CUDA_TASK_ERROR.
**/
{
 int err_code,cur_gpu;
 cudaError_t err;
 err_code=CUDA_TASK_EMPTY;
 if(cuda_task != NULL){
  if(cuda_task->task_error != 0){ //Negative: Task in progress or empty; Positive: Task scheduling error occured
   cur_gpu=-1; err=cudaGetDevice(&cur_gpu); if(err != cudaSuccess) return CUDA_TASK_ERROR;
   if(cur_gpu != cuda_task->gpu_id){
    err=cudaSetDevice(cuda_task->gpu_id); if(err != cudaSuccess){err=cudaSetDevice(cur_gpu); return CUDA_TASK_ERROR;}
   }
   err=cudaStreamQuery(cuda_task->task_stream);
   if(err != cudaSuccess && err != cudaErrorInvalidResourceHandle){ //task is still in progress
    err_code=CUDA_TASK_SCHEDULED;
   }else{ //task completed successfully or has never been scheduled
    if(err == cudaErrorInvalidResourceHandle){ //stream does not exist
     err_code=CUDA_TASK_EMPTY;
    }else{
     err_code=CUDA_TASK_COMPLETED; if(cuda_task->task_error < 0) cuda_task->task_error=0;
    }
   }
   if(cur_gpu != cuda_task->gpu_id){err=cudaSetDevice(cur_gpu); if(err != cudaSuccess) return CUDA_TASK_ERROR;}
  }else{
   err_code=CUDA_TASK_COMPLETED;
  }
 }
 return err_code;
}

__host__ int cuda_task_wait(cudaTask_t *cuda_task)
/** Waits on accomplishment of a CUDA task: Returns the output of cuda_task_complete().
Possible values are CUDA_TASK_COMPLETED, CUDA_TASK_ERROR, CUDA_TASK_EMPTY. **/
{
 int i,j;
 i=CUDA_TASK_SCHEDULED; j=1;
 while(j>0){
  i=cuda_task_complete(cuda_task);
  if(i == CUDA_TASK_COMPLETED || i == CUDA_TASK_ERROR || i == CUDA_TASK_EMPTY) j--;
 }
 return i;
}

__host__ int cuda_tasks_wait(int num_tasks, cudaTask_t **cuda_tasks, int* task_stats)
/** Waits upon completion of a series of CUDA tasks. Returns 0 on success. **/
{
 int i,j,n;
 if(num_tasks >= 0){
  if(num_tasks > 0){
   if(cuda_tasks != NULL && task_stats != NULL){
    for(i=0;i<num_tasks;i++){task_stats[i]=CUDA_TASK_SCHEDULED;}
    n=num_tasks;
    while(n>0){
     for(i=0;i<num_tasks;i++){
      j=task_stats[i];
      if(j != CUDA_TASK_COMPLETED && j != CUDA_TASK_ERROR && j != CUDA_TASK_EMPTY){
       if(cuda_tasks[i] != NULL){
        j=cuda_task_complete(cuda_tasks[i]); task_stats[i]=j;
        if(j == CUDA_TASK_COMPLETED || j == CUDA_TASK_ERROR || j == CUDA_TASK_EMPTY) n--;
       }else{
        return 1;
       }
      }
     }
    }
   }else{
    return 2;
   }
  }
 }else{
  return 3;
 }
 return 0;
}

__host__ float cuda_task_time(const cudaTask_t *cuda_task, float *in_copy, float *out_copy, float *comp)
/** Returns the time (in seconds) the CUDA task took to complete (only when EVENT_RECORD != 0).
Also, in_copy is input copying time, out_copy is output copying time, and comp is computing time in sec. **/
{
 int cur_gpu;
 float time_ms;
 cudaError_t err;
 if(cuda_task != NULL){
  if(EVENT_RECORD != 0){
   cur_gpu=-1; err=cudaGetDevice(&cur_gpu); if(err != cudaSuccess) return -2.0f;
   if(cur_gpu != cuda_task->gpu_id){err=cudaSetDevice(cuda_task->gpu_id); if(err != cudaSuccess) return -3.0f;}
   err=cudaEventElapsedTime(&time_ms,cuda_task->task_start,cuda_task->task_comput); //time in miliseconds
   if(err == cudaSuccess){*in_copy=time_ms/1000.0f;}else{*in_copy=-1.0f;}
   err=cudaEventElapsedTime(&time_ms,cuda_task->task_comput,cuda_task->task_output); //time in miliseconds
   if(err == cudaSuccess){*comp=time_ms/1000.0f;}else{*comp=-1.0f;}
   err=cudaEventElapsedTime(&time_ms,cuda_task->task_output,cuda_task->task_finish); //time in miliseconds
   if(err == cudaSuccess){*out_copy=time_ms/1000.0f;}else{*out_copy=-1.0f;}
   err=cudaEventElapsedTime(&time_ms,cuda_task->task_start,cuda_task->task_finish); //time in miliseconds
   if(err == cudaSuccess){time_ms/=1000.0f;}else{time_ms=-1.0f;} //time in seconds
   if(cur_gpu != cuda_task->gpu_id){err=cudaSetDevice(cur_gpu); if(err != cudaSuccess) return -4.0f;}
   return time_ms;
  }else{ //timing events are disabled
   return -5.0f;
  }
 }else{
  return -13.0f; //empty task
 }
}

//----------------------------------------------------------------------------------------
//EXPORTED FUNCTIONS (callable from Fortran):
//----------------------------------------------------
// CPU->GPU TENSOR BLOCK COPY (blocking):
__host__ int gpu_put_arg(tensBlck_t *ctens) //Blocking
/** This function copies a tensor block from the Host argument buffer into a GPU argument buffer **/
{
 int i,dev_kind,dev_num,gpu_num;
 size_t tsize;
 cudaError_t err;
 if(ctens != NULL){
  if(ctens->elems_h != NULL && ctens->elems_d != NULL){
   if(ctens->rank >= 0 && ctens->data_kind > 0){ //data_kind = {4|8|16}
    err=cudaGetLastError(); err=cudaSuccess;
    dev_kind=DEV_HOST; dev_num=decode_device_id(ctens->device_id,&dev_kind);
    if(dev_kind == DEV_NVIDIA_GPU && dev_num >= 0){
     err=cudaGetDevice(&gpu_num); if(err != cudaSuccess) return 1;
     if(dev_num != gpu_num){err=cudaSetDevice(dev_num); if(err != cudaSuccess){err=cudaSetDevice(gpu_num); return 2;}}
     tsize=tensBlck_volume(ctens); //tensor block size (elements)
     if(tsize > 0){
      err=cudaMemcpy(ctens->elems_d,ctens->elems_h,tsize*(ctens->data_kind),cudaMemcpyHostToDevice);
      if(err != cudaSuccess){err=cudaSetDevice(gpu_num); return 3;}
      i=tensBlck_set_presence(ctens); if(i != 0){err=cudaSetDevice(gpu_num); return 4;}
      if(dev_num != gpu_num){err=cudaSetDevice(gpu_num); if(err != cudaSuccess) return 5;}
     }else{
      err=cudaSetDevice(gpu_num); return 6;
     }
    }else{
     return 7;
    }
   }else{
    return 8;
   }
  }else{
   return 9;
  }
 }else{
  return 10;
 }
 return 0;
}
//----------------------------------------------------
// GPU->CPU TENSOR BLOCK COPY (blocking):
__host__ int gpu_get_arg(tensBlck_t *ctens) //Blocking
/** This function copies a tensor block from a GPU argument buffer into the Host argument buffer **/
{
 int dev_kind,dev_num,gpu_num;
 size_t tsize;
 cudaError_t err;
 if(ctens != NULL){
  if(ctens->elems_h != NULL && ctens->elems_d != NULL){
   if(ctens->rank >= 0 && ctens->data_kind > 0){ //data_kind = {4|8|16}
    err=cudaGetLastError(); err=cudaSuccess;
    dev_kind=DEV_HOST; dev_num=decode_device_id(ctens->device_id,&dev_kind);
    if(dev_kind == DEV_NVIDIA_GPU && dev_num >= 0){
     if(ctens->device_id > 0){ //tensor block must be present on GPU
      err=cudaGetDevice(&gpu_num); if(err != cudaSuccess) return 1;
      if(dev_num != gpu_num){err=cudaSetDevice(dev_num); if(err != cudaSuccess){err=cudaSetDevice(gpu_num); return 2;}}
      tsize=tensBlck_volume(ctens); //tensor block size (elements)
      if(tsize > 0){
       err=cudaMemcpy(ctens->elems_h,ctens->elems_d,tsize*(ctens->data_kind),cudaMemcpyDeviceToHost);
       if(err != cudaSuccess){err=cudaSetDevice(gpu_num); return 3;}
       if(dev_num != gpu_num){err=cudaSetDevice(gpu_num); if(err != cudaSuccess) return 4;}
      }else{
       err=cudaSetDevice(gpu_num); return 5;
      }
     }else{
      return 6;
     }
    }else{
     return 7;
    }
   }else{
    return 8;
   }
  }else{
   return 9;
  }
 }else{
  return 10;
 }
 return 0;
}
//-----------------------------------------------------------------
// CPU->GPU TENSOR BLOCK COPY (non-blocking):
__host__ int gpu_put_arg_(tensBlck_t *ctens, cudaTask_t *cuda_task) //Non-blocking
/** This function copies a tensor block from the Host argument buffer into a GPU argument buffer **/
{
 int i,dev_kind,dev_num,gpu_num;
 size_t tsize;
 cudaStream_t cuda_stream;
 cudaEvent_t cuda_start,cuda_finish;
 cudaError_t err;
 if(ctens != NULL && cuda_task != NULL){
  if(ctens->elems_h != NULL && ctens->elems_d != NULL){
   if(ctens->rank >= 0 && ctens->data_kind > 0){ //data_kind = {4|8|16}
    err=cudaGetLastError(); err=cudaSuccess;
    dev_kind=DEV_HOST; dev_num=decode_device_id(ctens->device_id,&dev_kind);
    if(dev_kind == DEV_NVIDIA_GPU && dev_num >= 0){
     err=cudaGetDevice(&gpu_num); if(err != cudaSuccess){i=cuda_task_finalize(cuda_task,1); return 1;}
     if(dev_num != gpu_num){ //set new GPU
      err=cudaSetDevice(dev_num); if(err != cudaSuccess){i=cuda_task_finalize(cuda_task,2); err=cudaSetDevice(gpu_num); return 2;}
     }
     i=0; err=cudaStreamCreate(&cuda_stream); if(err != cudaSuccess) i++;
     err=cudaEventCreate(&cuda_start); if(err != cudaSuccess) i++;
     err=cudaEventCreate(&cuda_finish); if(err != cudaSuccess) i++;
     if(i != 0){
      i=cuda_task_record(cuda_task,3,dev_num,cuda_stream,cuda_start,cuda_finish,cuda_finish,cuda_finish,0,NULL);
      err=cudaSetDevice(gpu_num); return 3;
     }
     tsize=tensBlck_volume(ctens); //tensor block size (elements)
     if(tsize > 0){
      if(EVENT_RECORD != 0){
       err=cudaEventRecord(cuda_start,cuda_stream); if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,4,dev_num,cuda_stream,cuda_start,cuda_finish,cuda_finish,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 4;
       }
      }
      err=cudaMemcpyAsync(ctens->elems_d,ctens->elems_h,tsize*(ctens->data_kind),cudaMemcpyHostToDevice,cuda_stream);
      if(err != cudaSuccess){
       i=cuda_task_record(cuda_task,5,dev_num,cuda_stream,cuda_start,cuda_finish,cuda_finish,cuda_finish,0,NULL);
       err=cudaSetDevice(gpu_num); return 5;
      }
      if(EVENT_RECORD != 0){
       err=cudaEventRecord(cuda_finish,cuda_stream); if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,6,dev_num,cuda_stream,cuda_start,cuda_finish,cuda_finish,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 6;
       }
      }
      i=cuda_task_record(cuda_task,0,dev_num,cuda_stream,cuda_start,cuda_finish,cuda_finish,cuda_finish,0,NULL);
      if(i!=0){i=cuda_task_finalize(cuda_task,7,dev_num); err=cudaSetDevice(gpu_num); return 7;}
      if(dev_num != gpu_num){err=cudaSetDevice(gpu_num); if(err != cudaSuccess){i=cuda_task_finalize(cuda_task,8,dev_num); return 8;}}
     }else{
      i=cuda_task_record(cuda_task,9,dev_num,cuda_stream,cuda_start,cuda_finish,cuda_finish,cuda_finish,0,NULL);
      err=cudaSetDevice(gpu_num); return 9;
     }
    }else{
     i=cuda_task_finalize(cuda_task,10); return 10;
    }
   }else{
    i=cuda_task_finalize(cuda_task,11); return 11;
   }
  }else{
   i=cuda_task_finalize(cuda_task,12); return 12;
  }
 }else{
  i=cuda_task_finalize(cuda_task,13); return 13;
 }
 return 0;
}
//-----------------------------------------------------------------
// GPU->CPU TENSOR BLOCK COPY (non-blocking):
__host__ int gpu_get_arg_(tensBlck_t *ctens, cudaTask_t *cuda_task) //Non-blocking
/** This function copies a tensor block from a GPU argument buffer into the Host argument buffer **/
{
 int i,dev_kind,dev_num,gpu_num;
 size_t tsize;
 cudaStream_t cuda_stream;
 cudaEvent_t cuda_start,cuda_finish;
 cudaError_t err;
 if(ctens != NULL && cuda_task != NULL){
  if(ctens->elems_h != NULL && ctens->elems_d != NULL){
   if(ctens->rank >= 0 && ctens->data_kind > 0){ //data_kind = {4|8|16}
    err=cudaGetLastError(); err=cudaSuccess;
    dev_kind=DEV_HOST; dev_num=decode_device_id(ctens->device_id,&dev_kind);
    if(dev_kind == DEV_NVIDIA_GPU && dev_num >= 0){
     if(ctens->device_id > 0){ //tensor block must be present on GPU
      err=cudaGetDevice(&gpu_num); if(err != cudaSuccess){i=cuda_task_finalize(cuda_task,1); return 1;}
      if(dev_num != gpu_num){
       err=cudaSetDevice(dev_num); if(err != cudaSuccess){i=cuda_task_finalize(cuda_task,2); err=cudaSetDevice(gpu_num); return 2;}
      }
      i=0; err=cudaStreamCreate(&cuda_stream); if(err != cudaSuccess) i++;
      err=cudaEventCreate(&cuda_start); if(err != cudaSuccess) i++;
      err=cudaEventCreate(&cuda_finish); if(err != cudaSuccess) i++;
      if(i != 0){
       i=cuda_task_record(cuda_task,3,dev_num,cuda_stream,cuda_start,cuda_start,cuda_start,cuda_finish,0,NULL);
       err=cudaSetDevice(gpu_num); return 3;
      }
      tsize=tensBlck_volume(ctens); //tensor block size (elements)
      if(tsize > 0){
       if(EVENT_RECORD != 0){
        err=cudaEventRecord(cuda_start,cuda_stream); if(err != cudaSuccess){
         i=cuda_task_record(cuda_task,4,dev_num,cuda_stream,cuda_start,cuda_start,cuda_start,cuda_finish,0,NULL);
         err=cudaSetDevice(gpu_num); return 4;
        }
       }
       err=cudaMemcpyAsync(ctens->elems_h,ctens->elems_d,tsize*(ctens->data_kind),cudaMemcpyDeviceToHost,cuda_stream);
       if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,5,dev_num,cuda_stream,cuda_start,cuda_start,cuda_start,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 5;
       }
       if(EVENT_RECORD != 0){
        err=cudaEventRecord(cuda_finish,cuda_stream); if(err != cudaSuccess){
         i=cuda_task_record(cuda_task,6,dev_num,cuda_stream,cuda_start,cuda_start,cuda_start,cuda_finish,0,NULL);
         err=cudaSetDevice(gpu_num); return 6;
        }
       }
       i=cuda_task_record(cuda_task,0,dev_num,cuda_stream,cuda_start,cuda_start,cuda_start,cuda_finish,0,NULL);
       if(i!=0){i=cuda_task_finalize(cuda_task,7,dev_num); err=cudaSetDevice(gpu_num); return 7;}
       if(dev_num != gpu_num){err=cudaSetDevice(gpu_num); if(err != cudaSuccess){i=cuda_task_finalize(cuda_task,8,dev_num); return 8;}}
      }else{
       i=cuda_task_record(cuda_task,9,dev_num,cuda_stream,cuda_start,cuda_start,cuda_start,cuda_finish,0,NULL);
       err=cudaSetDevice(gpu_num); return 9;
      }
     }else{
      i=cuda_task_finalize(cuda_task,10); return 10;
     }
    }else{
     i=cuda_task_finalize(cuda_task,11); return 11;
    }
   }else{
    i=cuda_task_finalize(cuda_task,12); return 12;
   }
  }else{
   i=cuda_task_finalize(cuda_task,13); return 13;
  }
 }else{
  i=cuda_task_finalize(cuda_task,14); return 14;
 }
 return 0;
}
//-------------------------------------------------------------------------------
// SQUARED 2-NORM OF AN ARRAY (R4) RESIDING ON GPU (blocking):
__host__ int gpu_array_2norm2_r4(size_t arr_size, const float *arr, float *norm2)
/** This function computes the sum of squared elements of a <float>
array arr(0:arr_size-1) which already resides on GPU.
Executed on the currently set GPU device. **/
{
 int i,bx;
 size_t l,mthr,sz;
 float *bnorm2;
 const char *err_msg;
 cudaError_t err;
 if(arr != NULL && norm2 != NULL && arr_size > 0){
  err=cudaGetLastError(); err=cudaSuccess;
  mthr=MAX_CUDA_BLOCKS*THRDS_ARRAY_NORM2; *norm2=0.0f;
  err=cudaGetSymbolAddress((void**)&bnorm2,gpu_blck_norms2_r4); if(err != cudaSuccess) return 1;
  for(l=0;l<arr_size;l+=mthr){
   if(l+mthr > arr_size){sz=arr_size-l;}else{sz=mthr;}; bx=1+(int)((sz-1)/THRDS_ARRAY_NORM2);
   gpu_array_2norm2_r4__<<<bx,THRDS_ARRAY_NORM2,THRDS_ARRAY_NORM2*sizeof(float)>>>(sz,&arr[l],bnorm2);
   err=cudaDeviceSynchronize();
   if(err != cudaSuccess){
    err_msg=cudaGetErrorString(err);
    if(VERBOSE) printf("\n#ERROR(tensor_algebra_gpu_nvidia:gpu_array_2norm2_r4): Kernel sync error: %s\n",err_msg);
    return 2;
   }
   err=cudaGetLastError();
   if(err != cudaSuccess){
    err_msg=cudaGetErrorString(err);
    if(VERBOSE) printf("\n#ERROR(tensor_algebra_gpu_nvidia:gpu_array_2norm2_r4): Kernel error: %s\n",err_msg);
    return 3;
   }
   err=cudaMemcpyFromSymbol((void*)blck_norms2_r4,gpu_blck_norms2_r4,bx*sizeof(float),0,cudaMemcpyDeviceToHost);
   if(err != cudaSuccess){
    err_msg=cudaGetErrorString(err);
    if(VERBOSE) printf("\n#ERROR(tensor_algebra_gpu_nvidia:gpu_array_2norm2_r4): Copy error: %s\n",err_msg);
    return 4;
   }
   for(i=0;i<bx;i++) *norm2+=blck_norms2_r4[i];
  }
 }else{
  return 5;
 }
 return 0;
}
//---------------------------------------------------------------------------------
// SQUARED 2-NORM OF AN ARRAY (R8) RESIDING ON GPU (blocking):
__host__ int gpu_array_2norm2_r8(size_t arr_size, const double *arr, double *norm2)
/** This function computes the sum of squared elements of a <double>
 array arr(0:arr_size-1) which already resides on GPU.
Executed on the currently set GPU device. **/
{
 int i,bx;
 size_t l,mthr,sz;
 double *bnorm2;
 const char *err_msg;
 cudaError_t err;
 if(arr != NULL && norm2 != NULL && arr_size > 0){
  err=cudaGetLastError(); err=cudaSuccess;
  mthr=MAX_CUDA_BLOCKS*THRDS_ARRAY_NORM2; *norm2=0.0;
  err=cudaGetSymbolAddress((void**)&bnorm2,gpu_blck_norms2_r8); if(err != cudaSuccess) return 1;
  for(l=0;l<arr_size;l+=mthr){
   if(l+mthr > arr_size){sz=arr_size-l;}else{sz=mthr;}; bx=1+(int)((sz-1)/THRDS_ARRAY_NORM2);
   gpu_array_2norm2_r8__<<<bx,THRDS_ARRAY_NORM2,THRDS_ARRAY_NORM2*sizeof(double)>>>(sz,&arr[l],bnorm2);
   err=cudaDeviceSynchronize();
   if(err != cudaSuccess){
    err_msg=cudaGetErrorString(err);
    if(VERBOSE) printf("\n#ERROR(tensor_algebra_gpu_nvidia:gpu_array_2norm2_r8): Kernel sync error: %s\n",err_msg);
    return 2;
   }
   err=cudaGetLastError();
   if(err != cudaSuccess){
    err_msg=cudaGetErrorString(err);
    if(VERBOSE) printf("\n#ERROR(tensor_algebra_gpu_nvidia:gpu_array_2norm2_r8): Kernel error: %s\n",err_msg);
    return 3;
   }
   err=cudaMemcpyFromSymbol((void*)blck_norms2_r8,gpu_blck_norms2_r8,bx*sizeof(double),0,cudaMemcpyDeviceToHost);
   if(err != cudaSuccess){
    err_msg=cudaGetErrorString(err);
    if(VERBOSE) printf("\n#ERROR(tensor_algebra_gpu_nvidia:gpu_array_2norm2_r8): Copy error: %s\n",err_msg);
    return 4;
   }
   for(i=0;i<bx;i++) *norm2+=blck_norms2_r8[i];
  }
 }else{
  return 5;
 }
 return 0;
}
//---------------------------------------------------------------------------------------
// MATRIX MULTIPLICATION 'TN' (R4) (blocking):
__host__ int gpu_matrix_multiply_tn_r4(size_t ll, size_t lr, size_t lc,
                                       const float *lmat, const float *rmat, float *dmat)
/** dmat(0:ll-1,0:lr-1)+=lmat(0:lc-1,0:ll-1)*rmat(0:lc-1,0:lr-1)
All matrices are in Host memory. Executed on the currently set GPU device. **/
{
 size_t dsize,lsize,rsize;
 float *dptr,*lptr,*rptr;
 int bx,by,err_code;
 const char *err_msg;
 cudaError_t err;
 if(lc > 0 && ll > 0 && lr > 0 && lmat != NULL && rmat != NULL && dmat != NULL){
  err=cudaGetLastError(); err=cudaSuccess;
  dsize=ll*lr*sizeof(float); lsize=lc*ll*sizeof(float); rsize=lc*lr*sizeof(float);
  err_code=gpu_mem_alloc((void**)&dptr,dsize); if(err_code != 0) return 1;
  err_code=gpu_mem_alloc((void**)&lptr,lsize); if(err_code != 0) return 2;
  err_code=gpu_mem_alloc((void**)&rptr,rsize); if(err_code != 0) return 3;
  err=cudaMemcpy((void*)dptr,(void*)dmat,dsize,cudaMemcpyHostToDevice); if(err != cudaSuccess) return 4;
  err=cudaMemcpy((void*)lptr,(void*)lmat,lsize,cudaMemcpyHostToDevice); if(err != cudaSuccess) return 5;
  err=cudaMemcpy((void*)rptr,(void*)rmat,rsize,cudaMemcpyHostToDevice); if(err != cudaSuccess) return 6;
  err_code=gpu_get_error_count();
  bx=1+(ll-1)/MAT_MULT_TILE_DIMX; by=1+(lr-1)/MAT_MULT_TILE_DIMY; limit_cuda_blocks2d(MAX_CUDA_BLOCKS,&bx,&by);
  dim3 blcks(bx,by); dim3 thrds(MAT_MULT_TILE_DIMX,MAT_MULT_TILE_DIMY);
//printf("\n#DEBUG(tensor_algebra_gpu_nvidia:gpu_matrix_multiply_tn_r4): Running GPU kernel ..."); //debug
  gpu_matrix_multiply_tn_r4__<<<blcks,thrds>>>(ll,lr,lc,lptr,rptr,dptr);
  err=cudaDeviceSynchronize(); if(err != cudaSuccess) return 7;
  err=cudaGetLastError();
  if(err!=cudaSuccess){
   err_msg=cudaGetErrorString(err);
   if(VERBOSE) printf("\n#ERROR(tensor_algebra_gpu_nvidia:gpu_matrix_multiply_tn_r4): Kernel error: %s\n",err_msg);
   return 8;
  }
  if(gpu_get_error_count() > err_code) return 9;
//printf("Done: %d",err); //debug
  err=cudaMemcpy((void*)dmat,(void*)dptr,dsize,cudaMemcpyDeviceToHost); if(err != cudaSuccess) return 10;
  err=cudaDeviceSynchronize(); if(err != cudaSuccess) return 11;
  err_code=gpu_mem_free((void*)rptr); if(err_code != 0) return 12;
  err_code=gpu_mem_free((void*)lptr); if(err_code != 0) return 13;
  err_code=gpu_mem_free((void*)dptr); if(err_code != 0) return 14;
  err=cudaDeviceSynchronize(); if(err != cudaSuccess) return 15;
 }else{
  return 16;
 }
 return 0;
}
//------------------------------------------------------------------------------------------
// MATRIX MULTIPLICATION 'TN' (R8) (blocking):
__host__ int gpu_matrix_multiply_tn_r8(size_t ll, size_t lr, size_t lc,
                                       const double *lmat, const double *rmat, double *dmat)
/** dmat(0:ll-1,0:lr-1)+=lmat(0:lc-1,0:ll-1)*rmat(0:lc-1,0:lr-1)
All matrices are in Host memory. Executed on the currently set GPU device. **/
{
 size_t dsize,lsize,rsize;
 double *dptr,*lptr,*rptr;
 int bx,by,err_code;
 const char *err_msg;
 cudaError_t err;
 if(lc > 0 && ll > 0 && lr > 0 && lmat != NULL && rmat != NULL && dmat != NULL){
  err=cudaGetLastError(); err=cudaSuccess;
  dsize=ll*lr*sizeof(double); lsize=lc*ll*sizeof(double); rsize=lc*lr*sizeof(double);
  err_code=gpu_mem_alloc((void**)&dptr,dsize); if(err_code != 0) return 1;
  err_code=gpu_mem_alloc((void**)&lptr,lsize); if(err_code != 0) return 2;
  err_code=gpu_mem_alloc((void**)&rptr,rsize); if(err_code != 0) return 3;
  err=cudaMemcpy((void*)dptr,(void*)dmat,dsize,cudaMemcpyHostToDevice); if(err != cudaSuccess) return 4;
  err=cudaMemcpy((void*)lptr,(void*)lmat,lsize,cudaMemcpyHostToDevice); if(err != cudaSuccess) return 5;
  err=cudaMemcpy((void*)rptr,(void*)rmat,rsize,cudaMemcpyHostToDevice); if(err != cudaSuccess) return 6;
  err_code=gpu_get_error_count();
  bx=1+(ll-1)/MAT_MULT_TILE_DIMX; by=1+(lr-1)/MAT_MULT_TILE_DIMY; limit_cuda_blocks2d(MAX_CUDA_BLOCKS,&bx,&by);
  dim3 blcks(bx,by); dim3 thrds(MAT_MULT_TILE_DIMX,MAT_MULT_TILE_DIMY);
//printf("\n#DEBUG(tensor_algebra_gpu_nvidia:gpu_matrix_multiply_tn_r8): Running GPU kernel ..."); //debug
  gpu_matrix_multiply_tn_r8__<<<blcks,thrds>>>(ll,lr,lc,lptr,rptr,dptr);
  err=cudaDeviceSynchronize(); if(err != cudaSuccess) return 7;
  err=cudaGetLastError();
  if(err!=cudaSuccess){
   err_msg=cudaGetErrorString(err);
   if(VERBOSE) printf("\n#ERROR(tensor_algebra_gpu_nvidia:gpu_matrix_multiply_tn_r8): Kernel error: %s\n",err_msg);
   return 8;
  }
  if(gpu_get_error_count() > err_code) return 9;
//printf("Done: %d",err); //debug
  err=cudaMemcpy((void*)dmat,(void*)dptr,dsize,cudaMemcpyDeviceToHost); if(err != cudaSuccess) return 10;
  err=cudaDeviceSynchronize(); if(err != cudaSuccess) return 11;
  err_code=gpu_mem_free((void*)rptr); if(err_code != 0) return 12;
  err_code=gpu_mem_free((void*)lptr); if(err_code != 0) return 13;
  err_code=gpu_mem_free((void*)dptr); if(err_code != 0) return 14;
  err=cudaDeviceSynchronize(); if(err != cudaSuccess) return 15;
 }else{
  return 16;
 }
 return 0;
}
//------------------------------------------------------------------------------------------------------
// TENSOR BLOCK INITIALIZATION (non-blocking):
__host__ int gpu_tensor_block_init_(tensBlck_t *ctens, double val, int copy_back, cudaTask_t *cuda_task)
/** ctens=val: The GPU part of <ctens> will be initialized to value <val>.
If <copy_back> = 0, no copy back to Host: One must use gpu_get_arg() explicitly: Careful! **/
{
 size_t tsize;
 int i,bx,dev_num,dev_kind,gpu_num;
 cudaStream_t cuda_stream;
 cudaEvent_t cuda_start,cuda_output,cuda_finish;
 cudaError_t err;
 if(ctens != NULL && cuda_task != NULL){
  if(ctens->elems_h != NULL && ctens->elems_d != NULL){
   if(ctens->rank >= 0 && ctens->data_kind > 0){
    err=cudaGetLastError(); err=cudaSuccess;
    dev_num=decode_device_id(ctens->device_id,&dev_kind);
    if(dev_kind == DEV_NVIDIA_GPU && dev_num >= 0){
     err=cudaGetDevice(&gpu_num); if(err != cudaSuccess){i=cuda_task_finalize(cuda_task,1); return 1;}
     if(dev_num != gpu_num){err=cudaSetDevice(dev_num); if(err != cudaSuccess){i=cuda_task_finalize(cuda_task,2); err=cudaSetDevice(gpu_num); return 2;}}
     i=0; err=cudaStreamCreate(&cuda_stream); if(err != cudaSuccess) i++;
     err=cudaEventCreate(&cuda_start); if(err != cudaSuccess) i++;
     err=cudaEventCreate(&cuda_output); if(err != cudaSuccess) i++;
     err=cudaEventCreate(&cuda_finish); if(err != cudaSuccess) i++;
     if(i != 0){
      i=cuda_task_record(cuda_task,3,dev_num,cuda_stream,cuda_start,cuda_start,cuda_output,cuda_finish,0,NULL);
      err=cudaSetDevice(gpu_num); return 3;
     }
     tsize=tensBlck_volume(ctens); //tensor block size (elements)
     if(tsize > 0){
      if(EVENT_RECORD != 0){
       err=cudaEventRecord(cuda_start,cuda_stream); if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,4,dev_num,cuda_stream,cuda_start,cuda_start,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 4;
       }
      }
      bx=1+(tsize-1)/THRDS_ARRAY_INIT; if(bx > MAX_CUDA_BLOCKS) bx=MAX_CUDA_BLOCKS;
      switch(ctens->data_kind){
       case R4:
        gpu_array_init_r4__<<<bx,THRDS_ARRAY_INIT,0,cuda_stream>>>(tsize,(float*)(ctens->elems_d),(float)val);
        break;
       case R8:
        gpu_array_init_r8__<<<bx,THRDS_ARRAY_INIT,0,cuda_stream>>>(tsize,(double*)(ctens->elems_d),val);
        break;
       default:
        i=cuda_task_record(cuda_task,5,dev_num,cuda_stream,cuda_start,cuda_start,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 5;
      }
      if(EVENT_RECORD != 0){
       err=cudaEventRecord(cuda_output,cuda_stream); if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,6,dev_num,cuda_stream,cuda_start,cuda_start,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 6;
       }
      }
      if(copy_back != NO_COPY_BACK){
       err=cudaMemcpyAsync(ctens->elems_h,ctens->elems_d,tsize*(ctens->data_kind),cudaMemcpyDeviceToHost,cuda_stream);
       if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,7,dev_num,cuda_stream,cuda_start,cuda_start,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 7;
       }
      }
      if(EVENT_RECORD != 0){
       err=cudaEventRecord(cuda_finish,cuda_stream); if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,8,dev_num,cuda_stream,cuda_start,cuda_start,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 8;
       }
      }
      i=cuda_task_record(cuda_task,0,dev_num,cuda_stream,cuda_start,cuda_start,cuda_output,cuda_finish,0,NULL);
      if(i!=0){i=cuda_task_finalize(cuda_task,9,dev_num); err=cudaSetDevice(gpu_num); return 9;}
      if(dev_num != gpu_num){err=cudaSetDevice(gpu_num); if(err != cudaSuccess){i=cuda_task_finalize(cuda_task,10,dev_num); return 10;}}
     }else{
      i=cuda_task_record(cuda_task,11,dev_num,cuda_stream,cuda_start,cuda_start,cuda_output,cuda_finish,0,NULL);
      err=cudaSetDevice(gpu_num); return 11;
     }
    }else{
     i=cuda_task_finalize(cuda_task,12); return 12;
    }
   }else{
    i=cuda_task_finalize(cuda_task,13); return 13;
   }
  }else{
   i=cuda_task_finalize(cuda_task,14); return 14;
  }
 }else{
  i=cuda_task_finalize(cuda_task,15); return 15;
 }
 return 0;
}
//-------------------------------------------------------------------------------------------------------
// TENSOR BLOCK RESCALING (non-blocking):
__host__ int gpu_tensor_block_scale_(tensBlck_t *ctens, double val, int copy_back, cudaTask_t *cuda_task)
/** ctens*=val: The GPU part of <ctens> will be scaled by the scalar <val>.
If <copy_back> = 0, no copy back to Host: One must use gpu_get_arg() explicitly: Careful! **/
{
 size_t tsize;
 int i,bx,dev_num,dev_kind,gpu_num;
 cudaStream_t cuda_stream;
 cudaEvent_t cuda_start,cuda_comput,cuda_output,cuda_finish;
 cudaError_t err;
 if(ctens != NULL && cuda_task != NULL){
  if(ctens->elems_h != NULL && ctens->elems_d != NULL){
   if(ctens->rank >= 0 && ctens->data_kind > 0){
    err=cudaGetLastError(); err=cudaSuccess;
    dev_num=decode_device_id(ctens->device_id,&dev_kind);
    if(dev_kind == DEV_NVIDIA_GPU && dev_num >= 0){
     err=cudaGetDevice(&gpu_num); if(err != cudaSuccess){i=cuda_task_finalize(cuda_task,1); return 1;}
     if(dev_num != gpu_num){err=cudaSetDevice(dev_num); if(err != cudaSuccess){i=cuda_task_finalize(cuda_task,2); err=cudaSetDevice(gpu_num); return 2;}}
     i=0; err=cudaStreamCreate(&cuda_stream); if(err != cudaSuccess) i++;
     err=cudaEventCreate(&cuda_start); if(err != cudaSuccess) i++;
     err=cudaEventCreate(&cuda_comput); if(err != cudaSuccess) i++;
     err=cudaEventCreate(&cuda_output); if(err != cudaSuccess) i++;
     err=cudaEventCreate(&cuda_finish); if(err != cudaSuccess) i++;
     if(i != 0){
      i=cuda_task_record(cuda_task,3,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
      err=cudaSetDevice(gpu_num); return 3;
     }
     tsize=tensBlck_volume(ctens); //tensor block size (elements)
     if(tsize > 0){
      if(EVENT_RECORD != 0){
       err=cudaEventRecord(cuda_start,cuda_stream); if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,4,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 4;
       }
      }
      if(ctens->device_id < 0){
       err=cudaMemcpyAsync(ctens->elems_d,ctens->elems_h,tsize*(ctens->data_kind),cudaMemcpyHostToDevice,cuda_stream);
       if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,5,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 5;
       }
      }
      if(EVENT_RECORD != 0){
       err=cudaEventRecord(cuda_comput,cuda_stream); if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,6,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 6;
       }
      }
      bx=1+(tsize-1)/THRDS_ARRAY_SCALE; if(bx > MAX_CUDA_BLOCKS) bx=MAX_CUDA_BLOCKS;
      switch(ctens->data_kind){
       case R4:
        gpu_array_scale_r4__<<<bx,THRDS_ARRAY_SCALE,0,cuda_stream>>>(tsize,(float*)(ctens->elems_d),(float)val);
        break;
       case R8:
        gpu_array_scale_r8__<<<bx,THRDS_ARRAY_SCALE,0,cuda_stream>>>(tsize,(double*)(ctens->elems_d),val);
        break;
       default:
        i=cuda_task_record(cuda_task,7,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 7;
      }
      if(EVENT_RECORD != 0){
       err=cudaEventRecord(cuda_output,cuda_stream); if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,8,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 8;
       }
      }
      if(copy_back != NO_COPY_BACK){
       err=cudaMemcpyAsync(ctens->elems_h,ctens->elems_d,tsize*(ctens->data_kind),cudaMemcpyDeviceToHost,cuda_stream);
       if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,9,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 9;
       }
      }
      if(EVENT_RECORD != 0){
       err=cudaEventRecord(cuda_finish,cuda_stream); if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,10,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 10;
       }
      }
      i=cuda_task_record(cuda_task,0,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
      if(i!=0){i=cuda_task_finalize(cuda_task,11,dev_num); err=cudaSetDevice(gpu_num); return 11;}
      if(dev_num != gpu_num){err=cudaSetDevice(gpu_num); if(err != cudaSuccess){i=cuda_task_finalize(cuda_task,12,dev_num); return 12;}}
     }else{
      i=cuda_task_record(cuda_task,13,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
      err=cudaSetDevice(gpu_num); return 13;
     }
    }else{
     i=cuda_task_finalize(cuda_task,14); return 14;
    }
   }else{
    i=cuda_task_finalize(cuda_task,15); return 15;
   }
  }else{
   i=cuda_task_finalize(cuda_task,16); return 16;
  }
 }else{
  i=cuda_task_finalize(cuda_task,17); return 17;
 }
 return 0;
}
//----------------------------------------------------------------------------------------
// TENSOR ADDITION [DLF] (non-blocking):
__host__ int gpu_tensor_block_add_dlf_(tensBlck_t *ctens0, tensBlck_t *ctens1, double val,
                                       int copy_back, cudaTask_t *cuda_task)
/** ctens0+=ctens1*val:
If <copy_back> = 0, no copy back to Host: One must use gpu_get_arg() explicitly: Careful! **/
{
 size_t tsize;
 int i,bx,dev_num,dev_kind,gpu_num;
 cudaStream_t cuda_stream;
 cudaEvent_t cuda_start,cuda_comput,cuda_output,cuda_finish;
 cudaError_t err;
 if(ctens0 != NULL && ctens1 != NULL && cuda_task != NULL){
  if(ctens0->elems_h != NULL && ctens0->elems_d != NULL && ctens1->elems_h != NULL && ctens1->elems_d != NULL){
   if(ctens0->rank >= 0 && ctens1->rank == ctens0->rank && ctens0->data_kind > 0 && ctens1->data_kind == ctens0->data_kind){
    for(i=0;i<ctens0->rank;i++){if(ctens0->dims_h[i] != ctens1->dims_h[i]){bx=cuda_task_finalize(cuda_task,1); return 1;}}
    err=cudaGetLastError(); err=cudaSuccess;
    dev_num=decode_device_id(ctens0->device_id,&dev_kind); i=decode_device_id(ctens1->device_id,&bx);
    if(dev_kind == DEV_NVIDIA_GPU && dev_num >= 0 && bx == dev_kind && i == dev_num){
     err=cudaGetDevice(&gpu_num); if(err != cudaSuccess){i=cuda_task_finalize(cuda_task,2); return 2;}
     if(dev_num != gpu_num){err=cudaSetDevice(dev_num); if(err != cudaSuccess){i=cuda_task_finalize(cuda_task,3); err=cudaSetDevice(gpu_num); return 3;}}
     i=0; err=cudaStreamCreate(&cuda_stream); if(err != cudaSuccess) i++;
     err=cudaEventCreate(&cuda_start); if(err != cudaSuccess) i++;
     err=cudaEventCreate(&cuda_comput); if(err != cudaSuccess) i++;
     err=cudaEventCreate(&cuda_output); if(err != cudaSuccess) i++;
     err=cudaEventCreate(&cuda_finish); if(err != cudaSuccess) i++;
     if(i != 0){
      i=cuda_task_record(cuda_task,4,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
      err=cudaSetDevice(gpu_num); return 4;
     }
     tsize=tensBlck_volume(ctens0); //tensor block size (elements)
     if(tsize > 0 && tsize == tensBlck_volume(ctens1)){
      if(EVENT_RECORD != 0){
       err=cudaEventRecord(cuda_start,cuda_stream); if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,5,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 5;
       }
      }
      if(ctens0->device_id < 0){
       err=cudaMemcpyAsync(ctens0->elems_d,ctens0->elems_h,tsize*(ctens0->data_kind),cudaMemcpyHostToDevice,cuda_stream);
       if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,6,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 6;
       }
      }
      if(ctens1->device_id < 0){
       err=cudaMemcpyAsync(ctens1->elems_d,ctens1->elems_h,tsize*(ctens1->data_kind),cudaMemcpyHostToDevice,cuda_stream);
       if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,7,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 7;
       }
      }
      if(EVENT_RECORD != 0){
       err=cudaEventRecord(cuda_comput,cuda_stream); if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,8,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 8;
       }
      }
      bx=1+(tsize-1)/THRDS_ARRAY_ADD; if(bx > MAX_CUDA_BLOCKS) bx=MAX_CUDA_BLOCKS;
      switch(ctens0->data_kind){
       case R4:
        gpu_array_add_r4__<<<bx,THRDS_ARRAY_ADD,0,cuda_stream>>>(tsize,(float*)(ctens0->elems_d),
                                                                       (float*)(ctens1->elems_d),(float)val);
        break;
       case R8:
        gpu_array_add_r8__<<<bx,THRDS_ARRAY_ADD,0,cuda_stream>>>(tsize,(double*)(ctens0->elems_d),
                                                                       (double*)(ctens1->elems_d),val);
        break;
       default:
        i=cuda_task_record(cuda_task,9,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 9;
      }
      if(EVENT_RECORD != 0){
       err=cudaEventRecord(cuda_output,cuda_stream); if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,10,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 10;
       }
      }
      if(copy_back != NO_COPY_BACK){
       err=cudaMemcpyAsync(ctens0->elems_h,ctens0->elems_d,tsize*(ctens0->data_kind),cudaMemcpyDeviceToHost,cuda_stream);
       if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,11,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 11;
       }
      }
      if(EVENT_RECORD != 0){
       err=cudaEventRecord(cuda_finish,cuda_stream); if(err != cudaSuccess){
        i=cuda_task_record(cuda_task,12,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
        err=cudaSetDevice(gpu_num); return 12;
       }
      }
      i=cuda_task_record(cuda_task,0,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
      if(i!=0){i=cuda_task_finalize(cuda_task,13,dev_num); err=cudaSetDevice(gpu_num); return 13;}
      if(dev_num != gpu_num){err=cudaSetDevice(gpu_num); if(err != cudaSuccess){i=cuda_task_finalize(cuda_task,14,dev_num); return 14;}}
     }else{
      i=cuda_task_record(cuda_task,15,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,0,NULL);
      err=cudaSetDevice(gpu_num); return 15;
     }
    }else{
     i=cuda_task_finalize(cuda_task,16); return 16;
    }
   }else{
    i=cuda_task_finalize(cuda_task,17); return 17;
   }
  }else{
   i=cuda_task_finalize(cuda_task,18); return 18;
  }
 }else{
  i=cuda_task_finalize(cuda_task,19); return 19;
 }
 return 0;
}
//---------------------------------------------------------------------------------------------------
// TENSOR TRANSPOSE [DLF] (blocking):
__host__ int gpu_tensor_block_copy_dlf(const int *dim_trn, tensBlck_t *tens_in, tensBlck_t *tens_out)
/** tens_out=TRN(tens_in):
INPUT:
 # dim_trn[0,1..dim_num] - index permutation (O2N, sign-containing, numeration starts from 1);
 # tens_in - input tensor block;
OUTPUT:
 # tens_out - output tensor block;
NOTES:
 # If both arguments are scalars or tensors of volume 1, GPU will NOT be involved!
 # In the case of error, CUDA timing events will not be destroyed (when EVENT_RECORD!=0).
**/
{
 int i,j,n,dev_in,dev_out,gpu_num,cae,bx,ibus[MAX_TENSOR_RANK];
 size_t tsize;
 cudaError_t err;
 const char *err_msg;
 cudaEvent_t time_beg,time_end;
 float time_ms;

 if(tens_in == NULL || tens_out == NULL || dim_trn == NULL) return 1;
 if(tens_in->elems_h == NULL || tens_in->elems_d == NULL || tens_out->elems_h == NULL || tens_out->elems_d == NULL) return 2;
 if(tens_in->rank != tens_out->rank) return 3;
 if(tens_in->data_kind <= 0 || tens_in->data_kind != tens_out->data_kind) return 4;
 err=cudaGetLastError(); err=cudaSuccess;
 if(tens_in->rank == 0){
//0-rank tensors (scalars):
  switch(tens_in->data_kind){
   case R4:
    ((float*)(tens_out->elems_h))[0]=((float*)(tens_in->elems_h))[0];
    break;
   case R8:
    ((double*)(tens_out->elems_h))[0]=((double*)(tens_in->elems_h))[0];
    break;
   default:
    return 5;
  }
  i=tensBlck_set_absence(tens_out); //invalidate GPU copy of the output scalar
 }else if(tens_in->rank > 0){
//Non-trivial tensors (rank>0):
  n=tens_in->rank;
//DEBUG begin:
//printf("\n#DEBUG(tensor_algebra_gpu_nvidia:gpu_tensor_block_copy_dlf): rank %d, data_kind %d\n",n,tens_in->data_kind);
//for(i=0;i<n;i++) printf(" %d",tens_in->dims_h[i]); printf("\n");
//for(i=0;i<n;i++) printf(" %d",tens_out->dims_h[i]); printf("\n");
//for(i=1;i<=n;i++) printf(" %d",dim_trn[i]); printf("\n");
//DEBUG end.
// Argument check:
  for(i=0;i<n;i++) ibus[i]=0;
  for(i=1;i<=n;i++){j=dim_trn[i]; if(j>=1&&j<=n){if(ibus[j-1]==0){ibus[j-1]=i;}else{return 6;}}else{return 7;}}
  for(i=0;i<n;i++){if(tens_in->dims_h[i] != tens_out->dims_h[dim_trn[1+i]-1]) return 8;}
// Get the tensor block size:
  tsize=tensBlck_volume(tens_in); //tensor block size (elements)
  if(tsize == 1){ //tensor of volume 1
   switch(tens_in->data_kind){
    case R4:
     ((float*)(tens_out->elems_h))[0]=((float*)(tens_in->elems_h))[0];
     break;
    case R8:
     ((double*)(tens_out->elems_h))[0]=((double*)(tens_in->elems_h))[0];
     break;
    default:
     return 9;
   }
   i=tensBlck_set_absence(tens_out); //invalidate GPU copy of the output tensor (volume of 1)
  }else if(tsize > 1){ //tensor of volume > 1
   i=decode_device_id(tens_in->device_id,&dev_in); j=decode_device_id(tens_out->device_id,&dev_out);
   if(dev_in == DEV_NVIDIA_GPU && dev_out == dev_in && i >= 0 && j == i){
    err=cudaGetDevice(&gpu_num); if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); return 10;}
    err=cudaSetDevice(i); if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 11;}
    if(EVENT_RECORD != 0){
     err=cudaEventCreate(&time_beg); if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 12;}
     err=cudaEventCreate(&time_end); if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 13;}
    }
// Set up constant memory arguments (tensor block dimension extents, permutation):
    cae=tens_in->const_args_entry; if(cae < 0 || cae >= MAX_GPU_ARGS){err=cudaSetDevice(gpu_num); return 14;}
    err=cudaMemcpyToSymbol(const_args_dims,(void*)((*tens_in).dims_h),sizeof(int)*n,sizeof(int)*MAX_TENSOR_RANK*cae,cudaMemcpyHostToDevice);
    if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 15;}
    err=cudaMemcpyToSymbol(const_args_prmn,(void*)(&dim_trn[1]),sizeof(int)*n,sizeof(int)*MAX_TENSOR_RANK*cae,cudaMemcpyHostToDevice);
    if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 16;}
//  printf("\n#DEBUG(tensor_algebra_gpu_nvidia:gpu_tensor_block_copy_dlf): Constant argument entries: %d %d\n",cae,tens_out->const_args_entry); //debug
    if(tens_in->device_id < 0){ //check whether the input tensor argument is already in GPU memory
// Copy the input tensor block into GPU global memory:
//   printf("\n#DEBUG(tensor_algebra_gpu_nvidia:gpu_tensor_block_copy_dlf): HostToDevice copy: %p %p %d\n",tens_in->elems_h,tens_in->elems_d,tsize); //debug
     switch(tens_in->data_kind){
      case R4:
       err=cudaMemcpy(tens_in->elems_d,tens_in->elems_h,tsize*sizeof(float),cudaMemcpyHostToDevice);
       break;
      case R8:
       err=cudaMemcpy(tens_in->elems_d,tens_in->elems_h,tsize*sizeof(double),cudaMemcpyHostToDevice);
       break;
      default:
       err=cudaSetDevice(gpu_num); return 17;
     }
     if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 18;}
     i=tensBlck_set_presence(tens_in);
    }
// Transpose:
    j=gpu_get_error_count(); time_ms=0.0f;
    if(EVENT_RECORD != 0){err=cudaEventRecord(time_beg); if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 19;}}
    if(TRANS_SHMEM != 0){
     bx=1+(tsize-1)/THRDS_TENSOR_COPY; if(bx > MAX_CUDA_BLOCKS) bx=MAX_CUDA_BLOCKS;
     switch(tens_in->data_kind){
      case R4:
       gpu_tensor_block_copy_dlf_r4__<<<bx,THRDS_TENSOR_COPY>>>(0,0,n,cae,(float*)(tens_in->elems_d),(float*)(tens_out->elems_d)); //shared-memory tensor transpose
       break;
      case R8:
       gpu_tensor_block_copy_dlf_r8__<<<bx,THRDS_TENSOR_COPY>>>(0,0,n,cae,(double*)(tens_in->elems_d),(double*)(tens_out->elems_d)); //shared-memory tensor transpose
       break;
      default:
       err=cudaSetDevice(gpu_num); return 20;
     }
    }else{
     bx=1+(tsize-1)/THRDS_TENSOR_COPY_SCAT; if(bx > MAX_CUDA_BLOCKS) bx=MAX_CUDA_BLOCKS;
     switch(tens_in->data_kind){
      case R4:
       gpu_tensor_block_copy_scatter_dlf_r4__<<<bx,THRDS_TENSOR_COPY_SCAT>>>(0,0,n,cae,(float*)(tens_in->elems_d),(float*)(tens_out->elems_d)); //scattering tensor transpose
       break;
      case R8:
       gpu_tensor_block_copy_scatter_dlf_r8__<<<bx,THRDS_TENSOR_COPY_SCAT>>>(0,0,n,cae,(double*)(tens_in->elems_d),(double*)(tens_out->elems_d)); //scattering tensor transpose
       break;
      default:
       err=cudaSetDevice(gpu_num); return 21;
     }
    }
    if(EVENT_RECORD != 0){err=cudaEventRecord(time_end); if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 22;}}
    err=cudaDeviceSynchronize();
    if(err!=cudaSuccess){
     err_msg=cudaGetErrorString(err);
     if(VERBOSE) printf("\n#ERROR(tensor_algebra_gpu_nvidia:gpu_tensor_block_copy_dlf): Sync error: %s\n",err_msg);
     err=cudaSetDevice(gpu_num);
     return 23;
    }
    if(EVENT_RECORD != 0){err=cudaEventElapsedTime(&time_ms,time_beg,time_end); if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 24;}}
    err=cudaGetLastError();
    if(err!=cudaSuccess){
     err_msg=cudaGetErrorString(err);
     if(VERBOSE) printf("\n#ERROR(tensor_algebra_gpu_nvidia:gpu_tensor_block_copy_dlf): Kernel error: %s\n",err_msg);
     err=cudaSetDevice(gpu_num);
     return 25;
    }
    if(gpu_get_error_count() > j){err=cudaSetDevice(gpu_num); return 26;}
    if(PRINT_TIMING) printf("#DEBUG(tensor_algebra_gpu_nvidia:gpu_tensor_block_copy_dlf): Kernel (%d): Time %f: KT/s=%f \n",TRANS_SHMEM,time_ms/1000.0f,((float)(tsize*2))/time_ms);
    i=tensBlck_set_presence(tens_out);
// Copy the output tensor block back into the Host argument buffer:
//  printf("\n#DEBUG(tensor_algebra_gpu_nvidia:gpu_tensor_block_copy_dlf): DeviceToHost copy: %p %p %d\n",tens_out->elems_h,tens_out->elems_d,tsize); //debug
    switch(tens_out->data_kind){
     case R4:
      err=cudaMemcpy(tens_out->elems_h,tens_out->elems_d,tsize*sizeof(float),cudaMemcpyDeviceToHost);
      break;
     case R8:
      err=cudaMemcpy(tens_out->elems_h,tens_out->elems_d,tsize*sizeof(double),cudaMemcpyDeviceToHost);
      break;
     default:
      err=cudaSetDevice(gpu_num); return 27;
    }
    if(err!=cudaSuccess){
     err_msg=cudaGetErrorString(err);
     if(VERBOSE) printf("\n#ERROR(tensor_algebra_gpu_nvidia:gpu_tensor_block_copy_dlf): Copy back: %s\n",err_msg);
     err=cudaSetDevice(gpu_num);
     return 28;
    }
    err=cudaDeviceSynchronize(); if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 29;}
    if(EVENT_RECORD != 0){err=cudaEventDestroy(time_beg); err=cudaEventDestroy(time_end);} //destroy CUDA events
    err=cudaSetDevice(gpu_num); if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); return 30;} //restore old GPU
   }else{
    return 31;
   }
  }else{
   return 32;
  }
 }else{
  return 33;
 }
 return 0;
}
//----------------------------------------------------------------------------------------------------
// TENSOR TRANSPOSE [DLF] (non-blocking):
__host__ int gpu_tensor_block_copy_dlf_(const int *dim_trn, tensBlck_t *tens_in, tensBlck_t *tens_out,
                                        int copy_back, cudaTask_t *cuda_task)
/** tens_out=TRN(tens_in):
INPUT:
 # dim_trn[0,1..dim_num] - index permutation (O2N, sign-containing, numeration starts from 1);
 # tens_in - input tensor block;
 # copy_back - 0: Output will not be copied back to Host (careful!); 1: It will.
OUTPUT:
 # tens_out - output tensor block;
NOTES:
 # For scalar tensors and tensors of volume 1, <copy_back> will always be TRUE.
**/
{
 int i,j,n,dev_in,dev_out,gpu_num,cae,bx,ibus[MAX_TENSOR_RANK];
 size_t tsize;
 cudaStream_t cuda_stream;
 cudaEvent_t time_beg,time_comput,time_output,time_end;
 cudaError_t err;
 const char *err_msg;

 if(tens_in == NULL || tens_out == NULL || cuda_task == NULL || dim_trn == NULL){bx=cuda_task_finalize(cuda_task,1); return 1;}
 if(tens_in->elems_h == NULL || tens_in->elems_d == NULL || tens_out->elems_h == NULL || tens_out->elems_d == NULL){bx=cuda_task_finalize(cuda_task,2); return 2;}
 if(tens_in->rank != tens_out->rank){bx=cuda_task_finalize(cuda_task,3); return 3;}
 if(tens_in->data_kind <= 0 || tens_in->data_kind != tens_out->data_kind){bx=cuda_task_finalize(cuda_task,4); return 4;}
 err=cudaGetLastError(); err=cudaSuccess;
//Trivial 0-rank tensors (scalars):
 if(tens_in->rank == 0){
  i=decode_device_id(tens_in->device_id,&dev_in); j=decode_device_id(tens_out->device_id,&dev_out);
  if(dev_in == DEV_NVIDIA_GPU && dev_out == dev_in && i >= 0 && j == i){
   dev_in=i; err=cudaGetDevice(&gpu_num); if(err != cudaSuccess){bx=cuda_task_finalize(cuda_task,5); return 5;}
   if(dev_in != gpu_num){err=cudaSetDevice(dev_in); if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); bx=cuda_task_finalize(cuda_task,6); err=cudaSetDevice(gpu_num); return 6;}}
   i=0; err=cudaStreamCreate(&cuda_stream); if(err != cudaSuccess) i++;
   err=cudaEventCreate(&time_beg); if(err != cudaSuccess) i++;
   err=cudaEventCreate(&time_end); if(err != cudaSuccess) i++;
   if(i != 0){
    i=cuda_task_record(cuda_task,7,dev_in,cuda_stream,time_beg,time_end,time_end,time_end,0,NULL);
    err=cudaSetDevice(gpu_num); return 7;
   }
   if(EVENT_RECORD != 0){
    err=cudaEventRecord(time_beg,cuda_stream); if(err != cudaSuccess){
     i=cuda_task_record(cuda_task,8,dev_in,cuda_stream,time_beg,time_end,time_end,time_end,0,NULL);
     err=cudaSetDevice(gpu_num); return 8;
    }
   }
   switch(tens_in->data_kind){
    case R4:
     ((float*)(tens_out->elems_h))[0]=((float*)(tens_in->elems_h))[0];
     break;
    case R8:
     ((double*)(tens_out->elems_h))[0]=((double*)(tens_in->elems_h))[0];
     break;
    default:
     i=cuda_task_record(cuda_task,9,dev_in,cuda_stream,time_beg,time_end,time_end,time_end,0,NULL);
     err=cudaSetDevice(gpu_num); return 9;
   }
   err=cudaMemcpyAsync(tens_out->elems_d,tens_out->elems_h,tens_out->data_kind,cudaMemcpyHostToDevice,cuda_stream);
   if(err != cudaSuccess){
    i=cuda_task_record(cuda_task,10,dev_in,cuda_stream,time_beg,time_end,time_end,time_end,0,NULL);
    err=cudaSetDevice(gpu_num); return 10;
   }
   if(EVENT_RECORD != 0){
    err=cudaEventRecord(time_end,cuda_stream); if(err != cudaSuccess){
     i=cuda_task_record(cuda_task,11,dev_in,cuda_stream,time_beg,time_end,time_end,time_end,0,NULL);
     err=cudaSetDevice(gpu_num); return 11;
    }
   }
   i=cuda_task_record(cuda_task,0,dev_in,cuda_stream,time_beg,time_end,time_end,time_end,0,NULL);
   if(i!=0){bx=cuda_task_finalize(cuda_task,12,dev_in); err=cudaSetDevice(gpu_num); return 12;}
   if(dev_in != gpu_num){err=cudaSetDevice(gpu_num); if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); bx=cuda_task_finalize(cuda_task,13,dev_in); return 13;}}
  }else{
   bx=cuda_task_finalize(cuda_task,14); return 14;
  }
//Non-trivial tensors (rank>0):
 }else if(tens_in->rank > 0){
  n=tens_in->rank;
//DEBUG begin:
//printf("\n#DEBUG(tensor_algebra_gpu_nvidia:gpu_tensor_block_copy_dlf_): rank %d, data_kind %d\n",n,tens_in->data_kind);
//for(i=0;i<n;i++) printf(" %d",tens_in->dims_h[i]); printf("\n");
//for(i=0;i<n;i++) printf(" %d",tens_out->dims_h[i]); printf("\n");
//for(i=1;i<=n;i++) printf(" %d",dim_trn[i]); printf("\n");
//DEBUG end.
// Argument check:
  for(i=0;i<n;i++) ibus[i]=0;
  for(i=1;i<=n;i++){
   j=dim_trn[i];
   if(j>=1&&j<=n){
    if(ibus[j-1]==0){ibus[j-1]=i;}else{bx=cuda_task_finalize(cuda_task,15); return 15;}
   }else{
    bx=cuda_task_finalize(cuda_task,16); return 16;
   }
  }
  for(i=0;i<n;i++){
   if(tens_in->dims_h[i] != tens_out->dims_h[dim_trn[1+i]-1]){bx=cuda_task_finalize(cuda_task,17); return 17;}
  }
// Get the tensor block size:
  tsize=tensBlck_volume(tens_in); //tensor block size (elements)
  if(tsize == 1){ //tensor of volume 1
   i=decode_device_id(tens_in->device_id,&dev_in); j=decode_device_id(tens_out->device_id,&dev_out);
   if(dev_in == DEV_NVIDIA_GPU && dev_out == dev_in && i >= 0 && j == i){
    dev_in=i; err=cudaGetDevice(&gpu_num); if(err != cudaSuccess){bx=cuda_task_finalize(cuda_task,18); return 18;}
    if(dev_in != gpu_num){err=cudaSetDevice(dev_in); if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); bx=cuda_task_finalize(cuda_task,19); err=cudaSetDevice(gpu_num); return 19;}}
    i=0; err=cudaStreamCreate(&cuda_stream); if(err != cudaSuccess) i++;
    err=cudaEventCreate(&time_beg); if(err != cudaSuccess) i++;
    err=cudaEventCreate(&time_end); if(err != cudaSuccess) i++;
    if(i != 0){
     i=cuda_task_record(cuda_task,20,dev_in,cuda_stream,time_beg,time_end,time_end,time_end,0,NULL);
     err=cudaSetDevice(gpu_num); return 20;
    }
    if(EVENT_RECORD != 0){
     err=cudaEventRecord(time_beg,cuda_stream); if(err != cudaSuccess){
      i=cuda_task_record(cuda_task,21,dev_in,cuda_stream,time_beg,time_end,time_end,time_end,0,NULL);
      err=cudaSetDevice(gpu_num); return 21;
     }
    }
    switch(tens_in->data_kind){
     case R4:
      ((float*)(tens_out->elems_h))[0]=((float*)(tens_in->elems_h))[0];
      break;
     case R8:
      ((double*)(tens_out->elems_h))[0]=((double*)(tens_in->elems_h))[0];
      break;
     default:
      i=cuda_task_record(cuda_task,22,dev_in,cuda_stream,time_beg,time_end,time_end,time_end,0,NULL);
      err=cudaSetDevice(gpu_num); return 22;
    }
    err=cudaMemcpyAsync(tens_out->elems_d,tens_out->elems_h,tens_out->data_kind,cudaMemcpyHostToDevice,cuda_stream);
    if(err != cudaSuccess){
     i=cuda_task_record(cuda_task,23,dev_in,cuda_stream,time_beg,time_end,time_end,time_end,0,NULL);
     err=cudaSetDevice(gpu_num); return 23;
    }
    if(EVENT_RECORD != 0){
     err=cudaEventRecord(time_end,cuda_stream); if(err != cudaSuccess){
      i=cuda_task_record(cuda_task,24,dev_in,cuda_stream,time_beg,time_end,time_end,time_end,0,NULL);
      err=cudaSetDevice(gpu_num); return 24;
     }
    }
    i=cuda_task_record(cuda_task,0,dev_in,cuda_stream,time_beg,time_end,time_end,time_end,0,NULL);
    if(i!=0){bx=cuda_task_finalize(cuda_task,25,dev_in); err=cudaSetDevice(gpu_num); return 25;}
    if(dev_in != gpu_num){err=cudaSetDevice(gpu_num); if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); bx=cuda_task_finalize(cuda_task,26,dev_in); return 26;}}
   }else{
    bx=cuda_task_finalize(cuda_task,27); return 27;
   }
  }else if(tsize > 1){ //tensor of volume > 1
   i=decode_device_id(tens_in->device_id,&dev_in); j=decode_device_id(tens_out->device_id,&dev_out);
   if(dev_in == DEV_NVIDIA_GPU && dev_out == dev_in && i >= 0 && j == i){
    dev_in=i; //GPU ID on which this tensor transpose will be executed (where data resides or will reside)
    err=cudaGetDevice(&gpu_num); if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); bx=cuda_task_finalize(cuda_task,28); return 28;}
    if(dev_in != gpu_num){err=cudaSetDevice(dev_in); if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); bx=cuda_task_finalize(cuda_task,29); err=cudaSetDevice(gpu_num); return 29;}}
    i=0; err=cudaStreamCreate(&cuda_stream); if(err != cudaSuccess) i++;
    err=cudaEventCreate(&time_beg); if(err != cudaSuccess) i++;
    err=cudaEventCreate(&time_comput); if(err != cudaSuccess) i++;
    err=cudaEventCreate(&time_output); if(err != cudaSuccess) i++;
    err=cudaEventCreate(&time_end); if(err != cudaSuccess) i++;
    if(i != 0){
     i=cuda_task_record(cuda_task,30,dev_in,cuda_stream,time_beg,time_comput,time_output,time_end,0,NULL);
     err=cudaSetDevice(gpu_num); return 30;
    }
    if(EVENT_RECORD != 0){
     err=cudaEventRecord(time_beg,cuda_stream); if(err != cudaSuccess){
      i=cuda_task_record(cuda_task,31,dev_in,cuda_stream,time_beg,time_comput,time_output,time_end,0,NULL);
      err=cudaSetDevice(gpu_num); return 31;
     }
    }
// Set up constant memory arguments (tensor block dimension extents, permutation):
    cae=tens_in->const_args_entry;
    if(cae < 0 || cae >= MAX_GPU_ARGS){
     i=cuda_task_record(cuda_task,32,dev_in,cuda_stream,time_beg,time_comput,time_output,time_end,0,NULL);
     err=cudaSetDevice(gpu_num); return 32;
    }
    err=cudaMemcpyToSymbolAsync(const_args_dims,(void*)((*tens_in).dims_h),sizeof(int)*n,sizeof(int)*MAX_TENSOR_RANK*cae,cudaMemcpyHostToDevice,cuda_stream);
    if(err != cudaSuccess){
     err_msg=cudaGetErrorString(err); i=cuda_task_record(cuda_task,33,dev_in,cuda_stream,time_beg,time_comput,time_output,time_end,0,NULL);
     err=cudaSetDevice(gpu_num); return 33;
    }
    err=cudaMemcpyToSymbolAsync(const_args_prmn,(void*)(&dim_trn[1]),sizeof(int)*n,sizeof(int)*MAX_TENSOR_RANK*cae,cudaMemcpyHostToDevice,cuda_stream);
    if(err != cudaSuccess){
     err_msg=cudaGetErrorString(err); i=cuda_task_record(cuda_task,34,dev_in,cuda_stream,time_beg,time_comput,time_output,time_end,0,NULL);
     err=cudaSetDevice(gpu_num); return 34;
    }
//  printf("\n#DEBUG(tensor_algebra_gpu_nvidia:gpu_tensor_block_copy_dlf_): Constant argument entries: %d %d\n",cae,tens_out->const_args_entry); //debug
    if(tens_in->device_id < 0){ //check whether the input tensor argument is already in GPU memory
// Copy the input tensor block into GPU global memory:
//   printf("\n#DEBUG(tensor_algebra_gpu_nvidia:gpu_tensor_block_copy_dlf_): HostToDevice copy: %p %p %d\n",tens_in->elems_h,tens_in->elems_d,tsize); //debug
     switch(tens_in->data_kind){
      case R4:
       err=cudaMemcpyAsync(tens_in->elems_d,tens_in->elems_h,tsize*sizeof(float),cudaMemcpyHostToDevice,cuda_stream);
       break;
      case R8:
       err=cudaMemcpyAsync(tens_in->elems_d,tens_in->elems_h,tsize*sizeof(double),cudaMemcpyHostToDevice,cuda_stream);
       break;
      default:
       i=cuda_task_record(cuda_task,35,dev_in,cuda_stream,time_beg,time_comput,time_output,time_end,0,NULL);
       err=cudaSetDevice(gpu_num); return 35;
     }
     if(err != cudaSuccess){
      err_msg=cudaGetErrorString(err); i=cuda_task_record(cuda_task,36,dev_in,cuda_stream,time_beg,time_comput,time_output,time_end,0,NULL);
      err=cudaSetDevice(gpu_num); return 36;
     }
    }
    if(EVENT_RECORD != 0){
     err=cudaEventRecord(time_comput,cuda_stream); if(err != cudaSuccess){
      i=cuda_task_record(cuda_task,37,dev_in,cuda_stream,time_beg,time_comput,time_output,time_end,0,NULL);
      err=cudaSetDevice(gpu_num); return 37;
     }
    }
// Transpose:
    j=gpu_get_error_count();
    if(TRANS_SHMEM != 0){
     bx=1+(tsize-1)/THRDS_TENSOR_COPY; if(bx > MAX_CUDA_BLOCKS) bx=MAX_CUDA_BLOCKS;
     switch(tens_in->data_kind){
      case R4:
       gpu_tensor_block_copy_dlf_r4__<<<bx,THRDS_TENSOR_COPY,0,cuda_stream>>>(0,0,n,cae,(float*)(tens_in->elems_d),(float*)(tens_out->elems_d)); //shared-memory tensor transpose
       break;
      case R8:
       gpu_tensor_block_copy_dlf_r8__<<<bx,THRDS_TENSOR_COPY,0,cuda_stream>>>(0,0,n,cae,(double*)(tens_in->elems_d),(double*)(tens_out->elems_d)); //shared-memory tensor transpose
       break;
      default:
       i=cuda_task_record(cuda_task,38,dev_in,cuda_stream,time_beg,time_comput,time_output,time_end,0,NULL);
       err=cudaSetDevice(gpu_num); return 38;
     }
    }else{
     bx=1+(tsize-1)/THRDS_TENSOR_COPY_SCAT; if(bx > MAX_CUDA_BLOCKS) bx=MAX_CUDA_BLOCKS;
     switch(tens_in->data_kind){
      case R4:
       gpu_tensor_block_copy_scatter_dlf_r4__<<<bx,THRDS_TENSOR_COPY_SCAT,0,cuda_stream>>>(0,0,n,cae,(float*)(tens_in->elems_d),(float*)(tens_out->elems_d)); //scattering tensor transpose
       break;
      case R8:
       gpu_tensor_block_copy_scatter_dlf_r8__<<<bx,THRDS_TENSOR_COPY_SCAT,0,cuda_stream>>>(0,0,n,cae,(double*)(tens_in->elems_d),(double*)(tens_out->elems_d)); //scattering tensor transpose
       break;
      default:
       i=cuda_task_record(cuda_task,39,dev_in,cuda_stream,time_beg,time_comput,time_output,time_end,0,NULL);
       err=cudaSetDevice(gpu_num); return 39;
     }
    }
    if(gpu_get_error_count() > j){
     i=cuda_task_record(cuda_task,40,dev_in,cuda_stream,time_beg,time_comput,time_output,time_end,0,NULL);
     err=cudaSetDevice(gpu_num); return 40;
    }
// Copy the output tensor block back into the Host argument buffer:
    if(EVENT_RECORD != 0){
     err=cudaEventRecord(time_output,cuda_stream); if(err != cudaSuccess){
      i=cuda_task_record(cuda_task,41,dev_in,cuda_stream,time_beg,time_comput,time_output,time_end,0,NULL);
      err=cudaSetDevice(gpu_num); return 41;
     }
    }
    if(copy_back != NO_COPY_BACK){
//   printf("\n#DEBUG(tensor_algebra_gpu_nvidia:gpu_tensor_block_copy_dlf_): DeviceToHost copy: %p %p %d\n",tens_out->elems_h,tens_out->elems_d,tsize); //debug
     switch(tens_out->data_kind){
      case R4:
       err=cudaMemcpyAsync(tens_out->elems_h,tens_out->elems_d,tsize*sizeof(float),cudaMemcpyDeviceToHost,cuda_stream);
       break;
      case R8:
       err=cudaMemcpyAsync(tens_out->elems_h,tens_out->elems_d,tsize*sizeof(double),cudaMemcpyDeviceToHost,cuda_stream);
       break;
      default:
       i=cuda_task_record(cuda_task,42,dev_in,cuda_stream,time_beg,time_comput,time_output,time_end,0,NULL);
       err=cudaSetDevice(gpu_num); return 42;
     }
     if(err != cudaSuccess){
      err_msg=cudaGetErrorString(err);
      if(VERBOSE) printf("\n#ERROR(tensor_algebra_gpu_nvidia:gpu_tensor_block_copy_dlf_): Copy back: %s\n",err_msg);
      i=cuda_task_record(cuda_task,43,dev_in,cuda_stream,time_beg,time_comput,time_output,time_end,0,NULL);
      err=cudaSetDevice(gpu_num); return 43;
     }
    }
    if(EVENT_RECORD != 0){
     err=cudaEventRecord(time_end,cuda_stream); if(err != cudaSuccess){
      i=cuda_task_record(cuda_task,44,dev_in,cuda_stream,time_beg,time_comput,time_output,time_end,0,NULL);
      err=cudaSetDevice(gpu_num); return 44;
     }
    }
    i=cuda_task_record(cuda_task,0,dev_in,cuda_stream,time_beg,time_comput,time_output,time_end,0,NULL);
    if(i!=0){bx=cuda_task_finalize(cuda_task,45,dev_in); err=cudaSetDevice(gpu_num); return 45;}
    if(dev_in != gpu_num){err=cudaSetDevice(gpu_num); if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); bx=cuda_task_finalize(cuda_task,46,dev_in); return 46;}}
   }else{
    bx=cuda_task_finalize(cuda_task,47); return 47;
   }
  }else{
   bx=cuda_task_finalize(cuda_task,48); return 48;
  }
 }else{
  bx=cuda_task_finalize(cuda_task,49); return 49;
 }
 return 0;
}
//-------------------------------------------------------------------------------------------------------------
// TENSOR CONTRACTION (non-blocking):
__host__ int gpu_tensor_block_contract_dlf_(const int *cptrn, const tensBlck_t *ltens, const tensBlck_t *rtens,
                                            tensBlck_t *dtens, int copy_back, cudaTask_t *cuda_task)
/**
dtens(:)+=ltens(:)*rtens(:)
INPUT:
 # cptrn(1:lrank+rrank) - contraction pattern: Position correspondence:
                          Uncontracted indices are positive, contracted are negative;
 # ltens - left tensor argument;
 # rtens - right tensor argument;
 # dtens - destination tensor (initialized!);
 # copy_back - 0: Output will not be copied back to Host (careful!); 1: It will.
OUTPUT:
 # dtens - modified destination tensor;
 # cuda_task - returns the relevant information on the cuda task launched.
NOTES:
 # For all scalar tensors or all tensors of volume 1, <copy_back> will always be TRUE.
**/
{
 int i,j,ncd,nlu,nru,cae,non_triv,gpu_num,dev_num,dev_kind,bx,by;
 int dprm[1+MAX_TENSOR_RANK],lprm[1+MAX_TENSOR_RANK],rprm[1+MAX_TENSOR_RANK]; //the 1st element is the sign of the permutation
 size_t dsize,lsize,rsize,lc,ll,lr;
 int scr_entry_cnt,scr_entries[MAX_SCR_ENTRY_COUNT]; //additional GPU argument buffer entries (three at most)
 void *darg,*larg,*rarg,*alpha,*beta;
 cudaStream_t cuda_stream;
 cudaEvent_t cuda_start,cuda_comput,cuda_output,cuda_finish;
 cudaError_t err;
 const char *err_msg;
#ifndef NO_BLAS
 cublasStatus_t err_cublas;
#endif

 err=cudaGetLastError(); err=cudaSuccess; scr_entry_cnt=0;
//printf("\n#DEBUG(tensor_algebra_gpu_nvidia:gpu_tensor_block_contract_dlf_): GPU Tensor Contraction:\n"); //debug
//Check arguments:
 if(cptrn == NULL || dtens == NULL || ltens == NULL || rtens == NULL || cuda_task == NULL) return 1;
 if(ltens->elems_h == NULL || ltens->elems_d == NULL || rtens->elems_h == NULL || rtens->elems_d == NULL ||
    dtens->elems_h == NULL || dtens->elems_d == NULL){i=cuda_task_finalize(cuda_task,2); return 2;}
 if((*dtens).rank < 0 || (*dtens).rank > MAX_TENSOR_RANK ||
    (*ltens).rank < 0 || (*ltens).rank > MAX_TENSOR_RANK ||
    (*rtens).rank < 0 || (*rtens).rank > MAX_TENSOR_RANK){i=cuda_task_finalize(cuda_task,3); return 3;}
 if(!(dtens->data_kind > 0 && ltens->data_kind == dtens->data_kind && rtens->data_kind == dtens->data_kind)){i=cuda_task_finalize(cuda_task,4); return 4;}
//Check contraction pattern and dimension extent correspondence:
 for(i=0;i<dtens->rank;i++) dprm[i]=0; for(i=0;i<ltens->rank;i++) lprm[i]=0; for(i=0;i<rtens->rank;i++) rprm[i]=0;
 for(i=0;i<ltens->rank;i++){
  j=cptrn[i];
  if(j > 0){ //position in dtens
   if(j > dtens->rank){bx=cuda_task_finalize(cuda_task,5); return 5;}
   if(dtens->dims_h[j-1] != ltens->dims_h[i]){bx=cuda_task_finalize(cuda_task,6); return 6;}
   if(dprm[j-1] == 0){dprm[j-1]=1;}else{bx=cuda_task_finalize(cuda_task,7); return 7;}
  }else if(j < 0){ //position in rtens
   if(-j > rtens->rank){bx=cuda_task_finalize(cuda_task,8); return 8;}
   if(rtens->dims_h[-j-1] != ltens->dims_h[i]){bx=cuda_task_finalize(cuda_task,9); return 9;}
   if(cptrn[ltens->rank+(-j-1)] != -(i+1)){bx=cuda_task_finalize(cuda_task,10); return 10;}
   if(rprm[-j-1] == 0){rprm[-j-1]=1;}else{bx=cuda_task_finalize(cuda_task,11); return 11;}
  }else{
   bx=cuda_task_finalize(cuda_task,12); return 12;
  }
 }
 for(i=0;i<rtens->rank;i++){
  j=cptrn[ltens->rank+i];
  if(j > 0){ //position in dtens
   if(j > dtens->rank){bx=cuda_task_finalize(cuda_task,13); return 13;}
   if(dtens->dims_h[j-1] != rtens->dims_h[i]){bx=cuda_task_finalize(cuda_task,14); return 14;}
   if(dprm[j-1] == 0){dprm[j-1]=1;}else{bx=cuda_task_finalize(cuda_task,15); return 15;}
  }else if(j < 0){ //position in ltens
   if(-j > ltens->rank){bx=cuda_task_finalize(cuda_task,16); return 16;}
   if(ltens->dims_h[-j-1] != rtens->dims_h[i]){bx=cuda_task_finalize(cuda_task,17); return 17;}
   if(cptrn[-j-1] != -(i+1)){bx=cuda_task_finalize(cuda_task,18); return 18;}
   if(lprm[-j-1] == 0){lprm[-j-1]=1;}else{bx=cuda_task_finalize(cuda_task,19); return 19;}
  }else{
   bx=cuda_task_finalize(cuda_task,20); return 20;
  }
 }
 for(i=0;i<dtens->rank;i++) if(dprm[i] != 1){bx=cuda_task_finalize(cuda_task,21); return 21;}
//Switch GPUs, if needed:
 dev_num=decode_device_id(dtens->device_id,&dev_kind);
 i=decode_device_id(ltens->device_id,&bx); j=decode_device_id(rtens->device_id,&by);
 if(!(dev_kind == DEV_NVIDIA_GPU && bx == dev_kind && by == dev_kind && dev_num >= 0 && i == dev_num && j == dev_num)){i=cuda_task_finalize(cuda_task,22); return 22;}
 err=cudaGetDevice(&gpu_num); if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); i=cuda_task_finalize(cuda_task,23); return 23;}
 if(dev_num != gpu_num){err=cudaSetDevice(dev_num); if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); i=cuda_task_finalize(cuda_task,24); err=cudaSetDevice(gpu_num); return 24;}}
//Create a CUDA stream and events:
 i=0; err=cudaStreamCreate(&cuda_stream); if(err != cudaSuccess) i++;
 err=cudaEventCreate(&cuda_start); if(err != cudaSuccess) i++;
 err=cudaEventCreate(&cuda_comput); if(err != cudaSuccess) i++;
 err=cudaEventCreate(&cuda_output); if(err != cudaSuccess) i++;
 err=cudaEventCreate(&cuda_finish); if(err != cudaSuccess) i++;
 if(i != 0){
  i=cuda_task_record(cuda_task,25,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
  err=cudaSetDevice(gpu_num); return 25;
 }
//Contraction case: Multiplication of scalars:
 if((*dtens).rank == 0 && (*ltens).rank == 0 && (*rtens).rank == 0){
  if(EVENT_RECORD != 0){
   err=cudaEventRecord(cuda_start,cuda_stream); if(err != cudaSuccess){
    i=cuda_task_record(cuda_task,26,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
    err=cudaSetDevice(gpu_num); return 26;
   }
   cuda_comput=cuda_start; cuda_output=cuda_start;
  }
  switch(dtens->data_kind){
   case R4: ((float*)(dtens->elems_h))[0]+=(((float*)(ltens->elems_h))[0])*(((float*)(rtens->elems_h))[0]); break;
   case R8: ((double*)(dtens->elems_h))[0]+=(((double*)(ltens->elems_h))[0])*(((double*)(rtens->elems_h))[0]); break;
   default:
    i=cuda_task_record(cuda_task,27,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
    err=cudaSetDevice(gpu_num); return 27;
  }
  err=cudaMemcpyAsync(dtens->elems_d,dtens->elems_h,dtens->data_kind,cudaMemcpyHostToDevice,cuda_stream);
  if(err != cudaSuccess){
   i=cuda_task_record(cuda_task,28,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
   err=cudaSetDevice(gpu_num); return 28;
  }
  if(EVENT_RECORD != 0){
   err=cudaEventRecord(cuda_finish,cuda_stream); if(err != cudaSuccess){
    i=cuda_task_record(cuda_task,29,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
    err=cudaSetDevice(gpu_num); return 29;
   }
  }
//Other contraction cases (tensor rescaling, full contraction, tensor product, partial contraction):
 }else{
// Get contraction permutations for all arguments:
  get_contr_permutations((*ltens).rank,(*rtens).rank,cptrn,dprm,lprm,rprm,&ncd,&nlu,&nru,&i);
  if(i != 0){
   i=cuda_task_record(cuda_task,30,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
   err=cudaSetDevice(gpu_num); return 30;
  }
  for(i=0;i<dtens->rank;i++) dtens->prmn_h[i]=dprm[1+i]; //ignore permutaion sign
  for(i=0;i<ltens->rank;i++) ltens->prmn_h[i]=lprm[1+i]; //ignore permutaion sign
  for(i=0;i<rtens->rank;i++) rtens->prmn_h[i]=rprm[1+i]; //ignore permutaion sign
  dsize=tensBlck_volume(dtens); lsize=tensBlck_volume(ltens); rsize=tensBlck_volume(rtens);
  lc=1; ll=1; for(i=0;i<(*ltens).rank;i++){if(ltens->prmn_h[i] <= ncd){lc*=((*ltens).dims_h[i]);}else{ll*=((*ltens).dims_h[i]);}}
  if(lsize > 0 && rsize > 0 && dsize > 0 && dsize%ll == 0){
   lr=dsize/ll;
   if(rsize%lr == 0 && rsize/lr == lc){
/*DEBUG begin:
    printf(" Const args (d,l,r) : %d %d %d\n",dtens->const_args_entry,ltens->const_args_entry,rtens->const_args_entry); //debug
    printf(" Block sizes (d,l,r): %d %d %d\n",dsize,lsize,rsize); //debug
    printf(" Block ranks (d,l,r): %d %d %d\n",dtens->rank,ltens->rank,rtens->rank); //debug
    printf(" Contraction pattern:"); for(i=0;i<ltens->rank+rtens->rank;i++) printf(" %d",cptrn[i]); //debug
    printf("\n Contr/uncontr(l,r) : %d %d %d: %d %d %d\n",ncd,nlu,nru,lc,ll,lr); //debug
    printf(" D-permutation      :"); for(i=0;i<dtens->rank;i++) printf(" %d",dtens->prmn_h[i]); //debug
    printf("\n L-permutation      :"); for(i=0;i<ltens->rank;i++) printf(" %d",ltens->prmn_h[i]); //debug
    printf("\n R-permutation      :"); for(i=0;i<rtens->rank;i++) printf(" %d",rtens->prmn_h[i]); //debug
//DEBUG end.*/
// Record the start event:
    if(EVENT_RECORD != 0){
     err=cudaEventRecord(cuda_start,cuda_stream); if(err != cudaSuccess){
      i=cuda_task_record(cuda_task,31,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
      err_msg=cudaGetErrorString(err);
      if(VERBOSE) printf("\n#ERROR(tensor_algebra_gpu_nvidia:gpu_tensor_block_contract_dlf_): Unable to record the start event: %s\n",err_msg);
      err=cudaSetDevice(gpu_num); return 31;
     }
     if(LastTask[dev_num] != NULL){ //`This should be done atomically for thread safety
      err=cudaStreamWaitEvent(cuda_stream,(LastTask[dev_num])->task_comput,0); //input transfers should only begin after the previous task input transfers have completed
      if(err != cudaSuccess){
       i=cuda_task_record(cuda_task,76,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
       err_msg=cudaGetErrorString(err);
       if(VERBOSE) printf("\n#ERROR(tensor_algebra_gpu_nvidia:gpu_tensor_block_contract_dlf_): Unable to create a task dependency: %s\n",err_msg);
       err=cudaSetDevice(gpu_num); return 76;
      }
     }
    }
// Copy the arguments into GPU memory, unless they are already there:
//  Left tensor argument:
    if(ltens->rank > 0){ //ignore input scalars
     cae=(*ltens).const_args_entry; if(cae < 0 || cae >= MAX_GPU_ARGS){
      i=cuda_task_record(cuda_task,32,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
      err=cudaSetDevice(gpu_num); return 32;
     }
     err=cudaMemcpyToSymbolAsync(const_args_dims,(void*)(ltens->dims_h),sizeof(int)*((*ltens).rank),sizeof(int)*MAX_TENSOR_RANK*cae,cudaMemcpyHostToDevice,cuda_stream);
     if(err != cudaSuccess){
      i=cuda_task_record(cuda_task,33,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
      err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 33;
     }
     err=cudaMemcpyToSymbolAsync(const_args_prmn,(void*)(ltens->prmn_h),sizeof(int)*((*ltens).rank),sizeof(int)*MAX_TENSOR_RANK*cae,cudaMemcpyHostToDevice,cuda_stream);
     if(err != cudaSuccess){
      i=cuda_task_record(cuda_task,34,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
      err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 34;
     }
     if((*ltens).device_id < 0){ //tensor argument is not in device memory: abs(.device_id) = (1 + active GPU number)
      err=cudaMemcpyAsync(ltens->elems_d,ltens->elems_h,lsize*(ltens->data_kind),cudaMemcpyHostToDevice,cuda_stream);
      if(err != cudaSuccess){
       i=cuda_task_record(cuda_task,35,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
       err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 35;
      }
     }
    }
//  Right tensor argument:
    if(rtens->rank > 0){ //ignore input scalars
     cae=(*rtens).const_args_entry; if(cae < 0 || cae >= MAX_GPU_ARGS){
      i=cuda_task_record(cuda_task,36,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
      err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 36;
     }
     err=cudaMemcpyToSymbolAsync(const_args_dims,(void*)(rtens->dims_h),sizeof(int)*((*rtens).rank),sizeof(int)*MAX_TENSOR_RANK*cae,cudaMemcpyHostToDevice,cuda_stream);
     if(err != cudaSuccess){
      i=cuda_task_record(cuda_task,37,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
      err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 37;
     }
     err=cudaMemcpyToSymbolAsync(const_args_prmn,(void*)(rtens->prmn_h),sizeof(int)*((*rtens).rank),sizeof(int)*MAX_TENSOR_RANK*cae,cudaMemcpyHostToDevice,cuda_stream);
     if(err != cudaSuccess){
      i=cuda_task_record(cuda_task,38,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
      err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 38;
     }
     if((*rtens).device_id < 0){ //tensor argument is not in device memory: abs(.device_id) = (1 + active GPU number)
      err=cudaMemcpyAsync(rtens->elems_d,rtens->elems_h,rsize*(rtens->data_kind),cudaMemcpyHostToDevice,cuda_stream);
      if(err != cudaSuccess){
       i=cuda_task_record(cuda_task,39,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
       err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 39;
      }
     }
    }
//  Destination tensor argument (copied in in all cases):
    cae=(*dtens).const_args_entry; if(cae < 0 || cae >= MAX_GPU_ARGS){
     i=cuda_task_record(cuda_task,40,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
     err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 40;
    }
    err=cudaMemcpyToSymbolAsync(const_args_dims,(void*)(dtens->dims_h),sizeof(int)*((*dtens).rank),sizeof(int)*MAX_TENSOR_RANK*cae,cudaMemcpyHostToDevice,cuda_stream);
    if(err != cudaSuccess){
     i=cuda_task_record(cuda_task,41,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
     err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 41;
    }
    err=cudaMemcpyToSymbolAsync(const_args_prmn,(void*)(dtens->prmn_h),sizeof(int)*((*dtens).rank),sizeof(int)*MAX_TENSOR_RANK*cae,cudaMemcpyHostToDevice,cuda_stream);
    if(err != cudaSuccess){
     i=cuda_task_record(cuda_task,42,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
     err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 42;
    }
    if((*dtens).device_id < 0){ //tensor argument is not in device memory: abs(.device_id) = (1 + active GPU number)
     err=cudaMemcpyAsync(dtens->elems_d,dtens->elems_h,dsize*(dtens->data_kind),cudaMemcpyHostToDevice,cuda_stream);
     if(err != cudaSuccess){
      i=cuda_task_record(cuda_task,43,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
      err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 43;
     }
    }
// Record the event that all the input data has been copied to device:
    if(EVENT_RECORD != 0){
     err=cudaEventRecord(cuda_comput,cuda_stream); if(err != cudaSuccess){
      i=cuda_task_record(cuda_task,44,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
      err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 44;
     }
    }
// Transpose tensor arguments, if needed:
//  Left tensor argument:
    if(non_trivial_prmn((*ltens).rank,ltens->prmn_h) != 0){
     cae=ltens->const_args_entry;
     i=get_buf_entry_gpu(dev_num,lsize*(ltens->data_kind),(char**)&larg,&j);
     if(i != 0){
      j=cuda_task_record(cuda_task,45,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
      err=cudaSetDevice(gpu_num); if(i == TRY_LATER || i == DEVICE_UNABLE){return i;}else{return 45;}
     }
     if(scr_entry_cnt < MAX_SCR_ENTRY_COUNT){
      scr_entries[scr_entry_cnt++]=j;
     }else{
      i=cuda_task_record(cuda_task,46,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
      err=cudaSetDevice(gpu_num); return 46;
     }
     if(TRANS_SHMEM != 0){
      bx=1+(lsize-1)/THRDS_TENSOR_COPY; if(bx > MAX_CUDA_BLOCKS) bx=MAX_CUDA_BLOCKS;
      switch(ltens->data_kind){
       case R4:
        gpu_tensor_block_copy_dlf_r4__<<<bx,THRDS_TENSOR_COPY,0,cuda_stream>>>(0,0,(*ltens).rank,cae,(float*)(ltens->elems_d),(float*)larg);
        break;
       case R8:
        gpu_tensor_block_copy_dlf_r8__<<<bx,THRDS_TENSOR_COPY,0,cuda_stream>>>(0,0,(*ltens).rank,cae,(double*)(ltens->elems_d),(double*)larg);
        break;
       default:
        i=cuda_task_record(cuda_task,47,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
        err=cudaSetDevice(gpu_num); return 47;
      }
     }else{
      bx=1+(lsize-1)/THRDS_TENSOR_COPY_SCAT; if(bx > MAX_CUDA_BLOCKS) bx=MAX_CUDA_BLOCKS;
      switch(ltens->data_kind){
       case R4:
        gpu_tensor_block_copy_scatter_dlf_r4__<<<bx,THRDS_TENSOR_COPY_SCAT,0,cuda_stream>>>(0,0,(*ltens).rank,cae,(float*)(ltens->elems_d),(float*)larg);
        break;
       case R8:
        gpu_tensor_block_copy_scatter_dlf_r8__<<<bx,THRDS_TENSOR_COPY_SCAT,0,cuda_stream>>>(0,0,(*ltens).rank,cae,(double*)(ltens->elems_d),(double*)larg);
        break;
       default:
        i=cuda_task_record(cuda_task,48,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
        err=cudaSetDevice(gpu_num); return 48;
      }
     }
    }else{
     larg=ltens->elems_d;
    }
//  Right tensor argument:
    if(non_trivial_prmn((*rtens).rank,rtens->prmn_h) != 0){
     cae=rtens->const_args_entry;
     i=get_buf_entry_gpu(dev_num,rsize*(rtens->data_kind),(char**)&rarg,&j);
     if(i != 0){
      j=cuda_task_record(cuda_task,49,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
      err=cudaSetDevice(gpu_num); if(i == TRY_LATER || i == DEVICE_UNABLE){return i;}else{return 49;}
     }
     if(scr_entry_cnt < MAX_SCR_ENTRY_COUNT){
      scr_entries[scr_entry_cnt++]=j;
     }else{
      i=cuda_task_record(cuda_task,50,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
      err=cudaSetDevice(gpu_num); return 50;
     }
     if(TRANS_SHMEM != 0){
      bx=1+(rsize-1)/THRDS_TENSOR_COPY; if(bx > MAX_CUDA_BLOCKS) bx=MAX_CUDA_BLOCKS;
      switch(rtens->data_kind){
       case R4:
        gpu_tensor_block_copy_dlf_r4__<<<bx,THRDS_TENSOR_COPY,0,cuda_stream>>>(0,0,(*rtens).rank,cae,(float*)(rtens->elems_d),(float*)rarg);
        break;
       case R8:
        gpu_tensor_block_copy_dlf_r8__<<<bx,THRDS_TENSOR_COPY,0,cuda_stream>>>(0,0,(*rtens).rank,cae,(double*)(rtens->elems_d),(double*)rarg);
        break;
       default:
        i=cuda_task_record(cuda_task,51,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
        err=cudaSetDevice(gpu_num); return 51;
      }
     }else{
      bx=1+(rsize-1)/THRDS_TENSOR_COPY_SCAT; if(bx > MAX_CUDA_BLOCKS) bx=MAX_CUDA_BLOCKS;
      switch(rtens->data_kind){
       case R4:
        gpu_tensor_block_copy_scatter_dlf_r4__<<<bx,THRDS_TENSOR_COPY_SCAT,0,cuda_stream>>>(0,0,(*rtens).rank,cae,(float*)(rtens->elems_d),(float*)rarg);
        break;
       case R8:
        gpu_tensor_block_copy_scatter_dlf_r8__<<<bx,THRDS_TENSOR_COPY_SCAT,0,cuda_stream>>>(0,0,(*rtens).rank,cae,(double*)(rtens->elems_d),(double*)rarg);
        break;
       default:
        i=cuda_task_record(cuda_task,52,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
        err=cudaSetDevice(gpu_num); return 52;
      }
     }
    }else{
     rarg=rtens->elems_d;
    }
//  Destination tensor argument:
    non_triv=non_trivial_prmn((*dtens).rank,dtens->prmn_h);
    if(non_triv != 0){
     cae=dtens->const_args_entry;
     i=get_buf_entry_gpu(dev_num,dsize*(dtens->data_kind),(char**)&darg,&j);
     if(i != 0){
      j=cuda_task_record(cuda_task,53,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
      err=cudaSetDevice(gpu_num); if(i == TRY_LATER || i == DEVICE_UNABLE){return i;}else{return 53;}
     }
     if(scr_entry_cnt < MAX_SCR_ENTRY_COUNT){
      scr_entries[scr_entry_cnt++]=j;
     }else{
      i=cuda_task_record(cuda_task,54,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
      err=cudaSetDevice(gpu_num); return 54;
     }
     if(TRANS_SHMEM != 0){
      bx=1+(dsize-1)/THRDS_TENSOR_COPY; if(bx > MAX_CUDA_BLOCKS) bx=MAX_CUDA_BLOCKS;
      switch(dtens->data_kind){
       case R4:
        gpu_tensor_block_copy_dlf_r4__<<<bx,THRDS_TENSOR_COPY,0,cuda_stream>>>(0,1,(*dtens).rank,cae,(float*)(dtens->elems_d),(float*)darg);
        break;
       case R8:
        gpu_tensor_block_copy_dlf_r8__<<<bx,THRDS_TENSOR_COPY,0,cuda_stream>>>(0,1,(*dtens).rank,cae,(double*)(dtens->elems_d),(double*)darg);
        break;
       default:
        i=cuda_task_record(cuda_task,55,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
        err=cudaSetDevice(gpu_num); return 55;
      }
     }else{
      bx=1+(dsize-1)/THRDS_TENSOR_COPY_SCAT; if(bx > MAX_CUDA_BLOCKS) bx=MAX_CUDA_BLOCKS;
      switch(dtens->data_kind){
       case R4:
        gpu_tensor_block_copy_scatter_dlf_r4__<<<bx,THRDS_TENSOR_COPY_SCAT,0,cuda_stream>>>(0,1,(*dtens).rank,cae,(float*)(dtens->elems_d),(float*)darg);
        break;
       case R8:
        gpu_tensor_block_copy_scatter_dlf_r8__<<<bx,THRDS_TENSOR_COPY_SCAT,0,cuda_stream>>>(0,1,(*dtens).rank,cae,(double*)(dtens->elems_d),(double*)darg);
        break;
       default:
        i=cuda_task_record(cuda_task,56,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
        err=cudaSetDevice(gpu_num); return 56;
      }
     }
    }else{
     darg=dtens->elems_d;
    }
// Invoke a CUDA kernel corresponding to the current contraction case:
//  Right tensor rescaling:
    if(ltens->rank == 0){ //rtens->rank > 0
     bx=1+(rsize-1)/THRDS_ARRAY_SCALE; if(bx > MAX_CUDA_BLOCKS) bx=MAX_CUDA_BLOCKS;
     switch(rtens->data_kind){
      case R4:
       gpu_array_scale_r4__<<<bx,THRDS_ARRAY_SCALE,0,cuda_stream>>>(rsize,(float*)(rarg),((float*)(ltens->elems_h))[0]);
       break;
      case R8:
       gpu_array_scale_r8__<<<bx,THRDS_ARRAY_SCALE,0,cuda_stream>>>(rsize,(double*)(rarg),((double*)(ltens->elems_h))[0]);
       break;
      default:
       i=cuda_task_record(cuda_task,57,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
       err=cudaSetDevice(gpu_num); return 57;
     }
//  Left tensor rescaling:
    }else if(rtens->rank == 0){ //ltens->rank > 0
     bx=1+(lsize-1)/THRDS_ARRAY_SCALE; if(bx > MAX_CUDA_BLOCKS) bx=MAX_CUDA_BLOCKS;
     switch(ltens->data_kind){
      case R4:
       gpu_array_scale_r4__<<<bx,THRDS_ARRAY_SCALE,0,cuda_stream>>>(lsize,(float*)(larg),((float*)(rtens->elems_h))[0]);
       break;
      case R8:
       gpu_array_scale_r8__<<<bx,THRDS_ARRAY_SCALE,0,cuda_stream>>>(lsize,(double*)(larg),((double*)(rtens->elems_h))[0]);
       break;
      default:
       i=cuda_task_record(cuda_task,58,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
       err=cudaSetDevice(gpu_num); return 58;
     }
//  Full tensor contraction (via vector dot-product):
    }else if(dtens->rank == 0){ //ltens->rank > 0 && rtens->rank = ltens->rank && rsize = lsize
     bx=1+(lsize-1)/THRDS_ARRAY_SCALE; if(bx > MAX_CUDA_BLOCKS) bx=MAX_CUDA_BLOCKS;
     switch(ltens->data_kind){
      case R4:
       gpu_array_dot_product_r4__<<<bx,THRDS_ARRAY_SCALE,THRDS_ARRAY_SCALE*sizeof(float),cuda_stream>>>
                                 (lsize,(float*)larg,(float*)rarg,(float*)darg);
       break;
      case R8:
       gpu_array_dot_product_r8__<<<bx,THRDS_ARRAY_SCALE,THRDS_ARRAY_SCALE*sizeof(double),cuda_stream>>>
                                 (lsize,(double*)larg,(double*)rarg,(double*)darg);
       break;
      default:
       i=cuda_task_record(cuda_task,59,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
       err=cudaSetDevice(gpu_num); return 59;
     }
//  Tensor product (no contracted indices):
    }else if(dtens->rank == ltens->rank + rtens->rank){
     bx=1+(lsize-1)/THRDS_ARRAY_PRODUCT; by=1+(rsize-1)/THRDS_ARRAY_PRODUCT;
     limit_cuda_blocks2d(MAX_CUDA_BLOCKS,&bx,&by); dim3 blcks(bx,by);
     switch(dtens->data_kind){
      case R4:
       gpu_array_product_r4__<<<blcks,THRDS_ARRAY_PRODUCT,0,cuda_stream>>>
                             (lsize,(float*)larg,rsize,(float*)rarg,(float*)darg);
       break;
      case R8:
       gpu_array_product_r8__<<<blcks,THRDS_ARRAY_PRODUCT,0,cuda_stream>>>
                             (lsize,(double*)larg,rsize,(double*)rarg,(double*)darg);
       break;
      default:
       i=cuda_task_record(cuda_task,60,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
       err=cudaSetDevice(gpu_num); return 60;
     }
//  Partial tensor contraction (via TN matrix multiplication):
    }else{
#ifndef NO_BLAS
     if(DISABLE_BLAS == 0){
      i=0;
      switch(dtens->data_kind){
       case R4:
        err=cudaGetSymbolAddress(&alpha,sgemm_alpha); if(err != cudaSuccess) i++;
        err=cudaGetSymbolAddress(&beta,sgemm_beta); if(err != cudaSuccess) i++;
        break;
       case R8:
        err=cudaGetSymbolAddress(&alpha,dgemm_alpha); if(err != cudaSuccess) i++;
        err=cudaGetSymbolAddress(&beta,dgemm_beta); if(err != cudaSuccess) i++;
        break;
       default:
        i++;
      }
      if(i != 0){
       i=cuda_task_record(cuda_task,61,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
       err=cudaSetDevice(gpu_num); return 61;
      }
      err_cublas=cublasSetStream(cublas_handle[dev_num],cuda_stream);
      if(err_cublas != CUBLAS_STATUS_SUCCESS){
       i=cuda_task_record(cuda_task,62,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
       err=cudaSetDevice(gpu_num); return 62;
      }
      switch(dtens->data_kind){
       case R4:
        err_cublas=cublasSgemm(cublas_handle[dev_num],CUBLAS_OP_T,CUBLAS_OP_N,(int)ll,(int)lr,(int)lc,
                    (float*)alpha,(float*)larg,(int)lc,(float*)rarg,(int)lc,(float*)beta,(float*)darg,(int)ll);
        break;
       case R8:
        err_cublas=cublasDgemm(cublas_handle[dev_num],CUBLAS_OP_T,CUBLAS_OP_N,(int)ll,(int)lr,(int)lc,
                    (double*)alpha,(double*)larg,(int)lc,(double*)rarg,(int)lc,(double*)beta,(double*)darg,(int)ll);
        break;
       default:
        i=cuda_task_record(cuda_task,63,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
        err=cudaSetDevice(gpu_num); return 63;
      }
      if(err_cublas != CUBLAS_STATUS_SUCCESS){
       i=cuda_task_record(cuda_task,64,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
       err=cudaSetDevice(gpu_num); return 64;
      }
     }else{
      bx=1+(ll-1)/MAT_MULT_TILE_DIMX; by=1+(lr-1)/MAT_MULT_TILE_DIMY; limit_cuda_blocks2d(MAX_CUDA_BLOCKS,&bx,&by);
//    printf("\n#DEBUG(): CUDA exec conf: %d %d %d %d\n",bx,by,MAT_MULT_TILE_DIMX,MAT_MULT_TILE_DIMY); //debug
      dim3 blcks(bx,by); dim3 thrds(MAT_MULT_TILE_DIMX,MAT_MULT_TILE_DIMY);
      switch(dtens->data_kind){
       case R4:
        gpu_matrix_multiply_tn_r4__<<<blcks,thrds,0,cuda_stream>>>(ll,lr,lc,(float*)larg,(float*)rarg,(float*)darg);
        break;
       case R8:
        gpu_matrix_multiply_tn_r8__<<<blcks,thrds,0,cuda_stream>>>(ll,lr,lc,(double*)larg,(double*)rarg,(double*)darg);
        break;
       default:
        i=cuda_task_record(cuda_task,65,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
        err=cudaSetDevice(gpu_num); return 65;
      }
     }
#else
     bx=1+(ll-1)/MAT_MULT_TILE_DIMX; by=1+(lr-1)/MAT_MULT_TILE_DIMY; limit_cuda_blocks2d(MAX_CUDA_BLOCKS,&bx,&by);
//   printf("\n#DEBUG(): CUDA exec conf: %d %d %d %d\n",bx,by,MAT_MULT_TILE_DIMX,MAT_MULT_TILE_DIMY); //debug
     dim3 blcks(bx,by); dim3 thrds(MAT_MULT_TILE_DIMX,MAT_MULT_TILE_DIMY);
     switch(dtens->data_kind){
      case R4:
       gpu_matrix_multiply_tn_r4__<<<blcks,thrds,0,cuda_stream>>>(ll,lr,lc,(float*)larg,(float*)rarg,(float*)darg);
       break;
      case R8:
       gpu_matrix_multiply_tn_r8__<<<blcks,thrds,0,cuda_stream>>>(ll,lr,lc,(double*)larg,(double*)rarg,(double*)darg);
       break;
      default:
       i=cuda_task_record(cuda_task,66,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
       err=cudaSetDevice(gpu_num); return 66;
     }
#endif
    }
// Transpose back the destination tensor argument:
    if(non_triv != 0){
     cae=dtens->const_args_entry;
     if(TRANS_SHMEM != 0){
      bx=1+(dsize-1)/THRDS_TENSOR_COPY; if(bx > MAX_CUDA_BLOCKS) bx=MAX_CUDA_BLOCKS;
      switch(dtens->data_kind){
       case R4:
        gpu_tensor_block_copy_dlf_r4__<<<bx,THRDS_TENSOR_COPY,0,cuda_stream>>>(1,0,(*dtens).rank,cae,(float*)darg,(float*)(dtens->elems_d));
        break;
       case R8:
        gpu_tensor_block_copy_dlf_r8__<<<bx,THRDS_TENSOR_COPY,0,cuda_stream>>>(1,0,(*dtens).rank,cae,(double*)darg,(double*)(dtens->elems_d));
        break;
       default:
        i=cuda_task_record(cuda_task,67,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
        err=cudaSetDevice(gpu_num); return 67;
      }
     }else{
      bx=1+(dsize-1)/THRDS_TENSOR_COPY_SCAT; if(bx > MAX_CUDA_BLOCKS) bx=MAX_CUDA_BLOCKS;
      switch(dtens->data_kind){
       case R4:
        gpu_tensor_block_copy_scatter_dlf_r4__<<<bx,THRDS_TENSOR_COPY_SCAT,0,cuda_stream>>>(1,0,(*dtens).rank,cae,(float*)darg,(float*)(dtens->elems_d));
        break;
       case R8:
        gpu_tensor_block_copy_scatter_dlf_r8__<<<bx,THRDS_TENSOR_COPY_SCAT,0,cuda_stream>>>(1,0,(*dtens).rank,cae,(double*)darg,(double*)(dtens->elems_d));
        break;
       default:
        i=cuda_task_record(cuda_task,68,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
        err=cudaSetDevice(gpu_num); return 68;
      }
     }
    }
// Record the event that the output is ready on device:
    if(EVENT_RECORD != 0){
     err=cudaEventRecord(cuda_output,cuda_stream); if(err != cudaSuccess){
      i=cuda_task_record(cuda_task,69,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
      err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 69;
     }
    }
// Copy back the destination tensor argument, if needed:
    if(copy_back != NO_COPY_BACK){
     err=cudaMemcpyAsync(dtens->elems_h,dtens->elems_d,dsize*(dtens->data_kind),cudaMemcpyDeviceToHost,cuda_stream);
     if(err != cudaSuccess){
      i=cuda_task_record(cuda_task,70,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
      err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 70;
     }
    }
    if(EVENT_RECORD != 0){
     err=cudaEventRecord(cuda_finish,cuda_stream); if(err != cudaSuccess){
      i=cuda_task_record(cuda_task,71,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
      err_msg=cudaGetErrorString(err); err=cudaSetDevice(gpu_num); return 71;
     }
    }
   }else{
    i=cuda_task_record(cuda_task,72,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
    err=cudaSetDevice(gpu_num); return 72;
   }
  }else{
   i=cuda_task_record(cuda_task,73,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
   err=cudaSetDevice(gpu_num); return 73;
  }
 }
//Register the CUDA task:
 i=cuda_task_record(cuda_task,0,dev_num,cuda_stream,cuda_start,cuda_comput,cuda_output,cuda_finish,scr_entry_cnt,scr_entries);
 if(i!=0){i=cuda_task_finalize(cuda_task,74,dev_num); err=cudaSetDevice(gpu_num); return 74;}
 if(EVENT_RECORD != 0) LastTask[dev_num]=cuda_task; //record the last task for each GPU
 if(dev_num != gpu_num){err=cudaSetDevice(gpu_num); if(err!=cudaSuccess){err_msg=cudaGetErrorString(err); i=cuda_task_finalize(cuda_task,75,dev_num); return 75;}}
//printf("\n#DEBUG(tensor_algebra_gpu_nvidia:gpu_tensor_block_contract_dlf_): Scheduled Successfully.\n"); //debug
 return 0;
}
//-------------------------------------------------------------------------------------
//CUDA KERNELS:
// SQUARED 2-NORM OF AN ARRAY (R4):`Obsolete
__global__ void gpu_array_2norm2_r4__(size_t arr_size, const float *arr, float *bnorm2)
/** Computes the squared Euclidean (Frobenius) norm of an array arr(0:arr_size-1)
INPUT:
 # arr_size - size of the array;
 # arr(0:arr_size-1) - array;
OUTPUT:
 # bnorm2[0:gridDim.x-1] - squared 2-norm of a sub-array computed by each CUDA thread block;
**/
{
 size_t i,n;
 float _thread_norm2;
 extern __shared__ float thread_norms2_r4[];

 n=gridDim.x*blockDim.x; _thread_norm2=0.0f;
 for(i=blockIdx.x*blockDim.x+threadIdx.x;i<arr_size;i+=n){_thread_norm2+=arr[i]*arr[i];}
 thread_norms2_r4[threadIdx.x]=_thread_norm2;
 __syncthreads();
 if(threadIdx.x == 0){
  bnorm2[blockIdx.x]=thread_norms2_r4[0]; for(i=1;i<blockDim.x;i++){bnorm2[blockIdx.x]+=thread_norms2_r4[i];}
 }
 __syncthreads();
 return;
}
//---------------------------------------------------------------------------------------
// SQUARED 2-NORM OF AN ARRAY (R8):`Obsolete
__global__ void gpu_array_2norm2_r8__(size_t arr_size, const double *arr, double *bnorm2)
/** Computes the squared Euclidean (Frobenius) norm of an array arr(0:arr_size-1)
INPUT:
 # arr_size - size of the array;
 # arr(0:arr_size-1) - array;
OUTPUT:
 # bnorm2[0:gridDim.x-1] - squared 2-norm of a sub-array computed by each CUDA thread block;
**/
{
 size_t i,n;
 double _thread_norm2;
 extern __shared__ double thread_norms2_r8[];

 n=gridDim.x*blockDim.x; _thread_norm2=0.0;
 for(i=blockIdx.x*blockDim.x+threadIdx.x;i<arr_size;i+=n){_thread_norm2+=arr[i]*arr[i];}
 thread_norms2_r8[threadIdx.x]=_thread_norm2;
 __syncthreads();
 if(threadIdx.x == 0){
  bnorm2[blockIdx.x]=thread_norms2_r8[0]; for(i=1;i<blockDim.x;i++){bnorm2[blockIdx.x]+=thread_norms2_r8[i];}
 }
 __syncthreads();
 return;
}
//------------------------------------------------------------------------------------
// SUM OF THE SQUARES OF ALL ARRAY ELEMENTS (R4):
__global__ void gpu_array_norm2_r8__(size_t arr_size, const float *arr, float *bnorm2)
/** Computes the squared 2-norm of array arr(0:arr_size-1)
INPUT:
 # arr_size - size of the array;
 # arr(0:arr_size-1) - array;
OUTPUT:
 # bnorm2 - squared 2-norm of the array (resides on device as well);
**/
{
 size_t i,n;
 float _thread_norm2;
 extern __shared__ float thread_norms2_r4[];

 n=gridDim.x*blockDim.x; _thread_norm2=0.0f;
 for(i=blockIdx.x*blockDim.x+threadIdx.x;i<arr_size;i+=n){_thread_norm2+=arr[i]*arr[i];}
 thread_norms2_r4[threadIdx.x]=_thread_norm2;
 __syncthreads();
 if(threadIdx.x == 0){
  _thread_norm2=thread_norms2_r4[0]; for(i=1;i<blockDim.x;i++){_thread_norm2+=thread_norms2_r4[i];}
  i=1; while(i == 1){i=atomicMax(&norm2_wr_lock,1);} //waiting for a lock to unlock, then lock
  *bnorm2+=_thread_norm2;
  __threadfence();
  i=atomicExch(&norm2_wr_lock,0); //unlock
 }
 __syncthreads();
 return;
}
//--------------------------------------------------------------------------------------
// SUM OF THE SQUARES OF ALL ARRAY ELEMENTS (R8):
__global__ void gpu_array_norm2_r8__(size_t arr_size, const double *arr, double *bnorm2)
/** Computes the squared 2-norm of array arr(0:arr_size-1)
INPUT:
 # arr_size - size of the array;
 # arr(0:arr_size-1) - array;
OUTPUT:
 # bnorm2 - squared 2-norm of the array (resides on device as well);
**/
{
 size_t i,n;
 double _thread_norm2;
 extern __shared__ double thread_norms2_r8[];

 n=gridDim.x*blockDim.x; _thread_norm2=0.0;
 for(i=blockIdx.x*blockDim.x+threadIdx.x;i<arr_size;i+=n){_thread_norm2+=arr[i]*arr[i];}
 thread_norms2_r8[threadIdx.x]=_thread_norm2;
 __syncthreads();
 if(threadIdx.x == 0){
  _thread_norm2=thread_norms2_r8[0]; for(i=1;i<blockDim.x;i++){_thread_norm2+=thread_norms2_r8[i];}
  i=1; while(i == 1){i=atomicMax(&norm2_wr_lock,1);} //waiting for a lock to unlock, then lock
  *bnorm2+=_thread_norm2;
  __threadfence();
  i=atomicExch(&norm2_wr_lock,0); //unlock
 }
 __syncthreads();
 return;
}
//----------------------------------------------------------------------
// ARRAY INITIALIZATION (R4):
__global__ void gpu_array_init_r4__(size_t tsize, float *arr, float val)
/** arr(:)=val **/
{
 size_t _ti = blockIdx.x*blockDim.x + threadIdx.x;
 size_t _gd = gridDim.x*blockDim.x;
 for(size_t l=_ti;l<tsize;l+=_gd){arr[l]=val;}
 return;
}
//------------------------------------------------------------------------
// ARRAY INITIALIZATION (R8):
__global__ void gpu_array_init_r8__(size_t tsize, double *arr, double val)
/** arr(:)=val **/
{
 size_t _ti = blockIdx.x*blockDim.x + threadIdx.x;
 size_t _gd = gridDim.x*blockDim.x;
 for(size_t l=_ti;l<tsize;l+=_gd){arr[l]=val;}
 return;
}
//-----------------------------------------------------------------------
// ARRAY RESCALING (R4):
__global__ void gpu_array_scale_r4__(size_t tsize, float *arr, float val)
/** arr(:)*=val **/
{
 size_t _ti = blockIdx.x*blockDim.x + threadIdx.x;
 size_t _gd = gridDim.x*blockDim.x;
 for(size_t l=_ti;l<tsize;l+=_gd){arr[l]*=val;}
 return;
}
//-------------------------------------------------------------------------
// ARRAY RESCALING (R8):
__global__ void gpu_array_scale_r8__(size_t tsize, double *arr, double val)
/** arr(:)*=val **/
{
 size_t _ti = blockIdx.x*blockDim.x + threadIdx.x;
 size_t _gd = gridDim.x*blockDim.x;
 for(size_t l=_ti;l<tsize;l+=_gd){arr[l]*=val;}
 return;
}
//-------------------------------------------------------------------------------------------------------------------
// ARRAY ADDITION (R4):
__global__ void gpu_array_add_r4__(size_t tsize, float* __restrict__ arr0, const float* __restrict__ arr1, float val)
/** arr0(:)+=arr1(:)*val **/
{
 size_t _ti = blockIdx.x*blockDim.x + threadIdx.x;
 size_t _gd = gridDim.x*blockDim.x;
 for(size_t l=_ti;l<tsize;l+=_gd){arr0[l]+=(arr1[l]*val);}
 return;
}
//----------------------------------------------------------------------------------------------------------------------
// ARRAY ADDITION (R8):
__global__ void gpu_array_add_r8__(size_t tsize, double* __restrict__ arr0, const double* __restrict__ arr1, double val)
/** arr0(:)+=arr1(:)*val **/
{
 size_t _ti = blockIdx.x*blockDim.x + threadIdx.x;
 size_t _gd = gridDim.x*blockDim.x;
 for(size_t l=_ti;l<tsize;l+=_gd){arr0[l]+=(arr1[l]*val);}
 return;
}
//-------------------------------------------------------------------------------------------------------------------
// ARRAY DOT-PRODUCT (R4):
__global__ void gpu_array_dot_product_r4__(size_t tsize, const float *arr1, const float *arr2, volatile float *dprod)
{
 extern __shared__ float dprs_r4[]; //volume = blockDim.x
 size_t l;
 int i,j;
 float dpr;
 dpr=0.0f; for(l=blockIdx.x*blockDim.x+threadIdx.x;l<tsize;l+=gridDim.x*blockDim.x){dpr+=arr1[l]*arr2[l];}
 dprs_r4[threadIdx.x]=dpr;
 __syncthreads();
 i=1; while(i < blockDim.x){j=threadIdx.x*(i*2); if(j+i < blockDim.x) dprs_r4[j]+=dprs_r4[j+i]; i*=2;}
 __syncthreads();
 if(threadIdx.x == 0){
  i=1; while(i == 1){i=atomicMax(&dot_product_wr_lock,1);} //waiting for a lock to unlock, then lock
  *dprod+=dprs_r4[0];
  __threadfence();
  i=atomicExch(&dot_product_wr_lock,0); //unlock
 }
 __syncthreads();
 return;
}
//----------------------------------------------------------------------------------------------------------------------
// ARRAY DOT-PRODUCT (R8):
__global__ void gpu_array_dot_product_r8__(size_t tsize, const double *arr1, const double *arr2, volatile double *dprod)
{
 extern __shared__ double dprs_r8[]; //volume = blockDim.x
 size_t l;
 int i,j;
 double dpr;
 dpr=0.0; for(l=blockIdx.x*blockDim.x+threadIdx.x;l<tsize;l+=gridDim.x*blockDim.x){dpr+=arr1[l]*arr2[l];}
 dprs_r8[threadIdx.x]=dpr;
 __syncthreads();
 i=1; while(i < blockDim.x){j=threadIdx.x*(i*2); if(j+i < blockDim.x) dprs_r8[j]+=dprs_r8[j+i]; i*=2;}
 __syncthreads();
 if(threadIdx.x == 0){
  i=1; while(i == 1){i=atomicMax(&dot_product_wr_lock,1);} //waiting for a lock to unlock, then lock
  *dprod+=dprs_r8[0];
  __threadfence();
  i=atomicExch(&dot_product_wr_lock,0); //unlock
 }
 __syncthreads();
 return;
}
//-------------------------------------------------------------------------------------------------------------
// ARRAY PRODUCT (R4):
__global__ void gpu_array_product_r4__(size_t tsize1, const float* __restrict__ arr1, size_t tsize2,
                                                      const float* __restrict__ arr2, float* __restrict__ arr0)
/** arr0[0:tsize2-1][0:tsize1-1]+=arr1[0:tsize1-1]*arr2[0:tsize2-1] **/
{
 __shared__ float lbuf[THRDS_ARRAY_PRODUCT+1],rbuf[THRDS_ARRAY_PRODUCT];
 size_t _ib,_in,_jb,_jn,_tx,_jc;
 _tx=(size_t)threadIdx.x;
// if(tsize1 >= THRDS_ARRAY_PRODUCT){ //large or medium size L
  for(_jb=blockIdx.y*THRDS_ARRAY_PRODUCT;_jb<tsize2;_jb+=gridDim.y*THRDS_ARRAY_PRODUCT){
   if(_jb+THRDS_ARRAY_PRODUCT > tsize2){_jn=tsize2-_jb;}else{_jn=THRDS_ARRAY_PRODUCT;}
   if(_tx < _jn) rbuf[_tx]=arr2[_jb+_tx];
   for(_ib=blockIdx.x*THRDS_ARRAY_PRODUCT;_ib<tsize1;_ib+=gridDim.x*THRDS_ARRAY_PRODUCT){
    if(_ib+THRDS_ARRAY_PRODUCT > tsize1){_in=tsize1-_ib;}else{_in=THRDS_ARRAY_PRODUCT;}
    if(_tx < _in) lbuf[_tx]=arr1[_ib+_tx];
    __syncthreads();
    for(_jc=0;_jc<_jn;_jc++){if(_tx < _in) arr0[(_jb+_jc)*tsize1+_ib+_tx]+=lbuf[_tx]*rbuf[_jc];}
    __syncthreads();
   }
  }
// }else{ //small size L
  //`Write
// }
 return;
}
//---------------------------------------------------------------------------------------------------------------
// ARRAY PRODUCT (R8):
__global__ void gpu_array_product_r8__(size_t tsize1, const double* __restrict__ arr1, size_t tsize2,
                                                      const double* __restrict__ arr2, double* __restrict__ arr0)
/** arr0[0:tsize2-1][0:tsize1-1]+=arr1[0:tsize1-1]*arr2[0:tsize2-1] **/
{
 __shared__ double lbuf[THRDS_ARRAY_PRODUCT+1],rbuf[THRDS_ARRAY_PRODUCT];
 size_t _ib,_in,_jb,_jn,_tx,_jc;
 _tx=(size_t)threadIdx.x;
// if(tsize1 >= THRDS_ARRAY_PRODUCT){ //large or medium size L
  for(_jb=blockIdx.y*THRDS_ARRAY_PRODUCT;_jb<tsize2;_jb+=gridDim.y*THRDS_ARRAY_PRODUCT){
   if(_jb+THRDS_ARRAY_PRODUCT > tsize2){_jn=tsize2-_jb;}else{_jn=THRDS_ARRAY_PRODUCT;}
   if(_tx < _jn) rbuf[_tx]=arr2[_jb+_tx];
   for(_ib=blockIdx.x*THRDS_ARRAY_PRODUCT;_ib<tsize1;_ib+=gridDim.x*THRDS_ARRAY_PRODUCT){
    if(_ib+THRDS_ARRAY_PRODUCT > tsize1){_in=tsize1-_ib;}else{_in=THRDS_ARRAY_PRODUCT;}
    if(_tx < _in) lbuf[_tx]=arr1[_ib+_tx];
    __syncthreads();
    for(_jc=0;_jc<_jn;_jc++){if(_tx < _in) arr0[(_jb+_jc)*tsize1+_ib+_tx]+=lbuf[_tx]*rbuf[_jc];}
    __syncthreads();
   }
  }
// }else{ //small size L
  //`Write
// }
 return;
}
//-------------------------------------------------------------------------------------------------------------
// TENSOR TRANSPOSE (R4) (shared-memory version):
__global__ void gpu_tensor_block_copy_dlf_r4__(int dmo, int drc, int dim_num, int const_args_pos,
                                               const float* __restrict__ tens_in, float* __restrict__ tens_out)
/**
Shared-memory version of tensor transpose: tens_out=TRN(tens_in):
INPUT:
 # dmo - dimension extents order (0: normal, as it is in <const_args>; not 0: permuted dimension order will be imposed);
 # drc - index permutation direction (0: normal, as it is in <const_args>; not 0: inversed permutation will be used);
 # dim_num - tensor block rank;
 # const_args_pos - entry in the __constant__ memory bank where tensor block dimension extents (const_args_dims)
                    and index permutation (const_args_prmn) are stored;
 # tens_in[0:] - input tensor;
OUTPUT:
 # tens_out[0:] - output (transposed) tensor;
NOTES:
 # Minimal CUDA execution configuration is <<<1,warpSize>>>
 # Number of threads per block must be multiple of the warpSize!
**/
{
 __shared__ float buf0[TENS_TRANSP_BUF_SIZE],val;
 __shared__ size_t base_in[MAX_TENSOR_RANK],base_out[MAX_TENSOR_RANK];
 __shared__ size_t ftb[TENS_TRANSP_TAB_SIZE],gtb[TENS_TRANSP_TAB_SIZE];
 __shared__ int htb[TENS_TRANSP_TAB_SIZE],stb[TENS_TRANSP_TAB_SIZE];
 __shared__ int dim_in[MAX_TENSOR_RANK],dim_out[MAX_TENSOR_RANK],o2n[MAX_TENSOR_RANK],n2o[MAX_TENSOR_RANK];
 __shared__ int pri[MAX_TENSOR_RANK],tmp0[MAX_TENSOR_RANK];
 __shared__ int err_code,minor,minor_in,minor_out,s1_ind,s1_ond,s1_step,s1_dim,s2_ind,s2_ond,s2_step,s2_dim,ns1,ns2;
 __shared__ size_t vol,vol_ext;
 size_t _vol,_addr_in,_addr_out,_addr,_work_piece;
 int i,j,k,l,m,n,_vol_minor,_vol_in,_vol_out,_s1,_s2;

//Determine the minor index set (only the master thread in each thread block):
 if(threadIdx.x == 0){
  err_code=0;
  if(dim_num >= 0 && dim_num <= MAX_TENSOR_RANK && blockDim.x >= warpSize && blockDim.x%warpSize == 0){
   s1_ind=dim_num+1; s2_ind=dim_num-1;
   _vol=1; for(i=0;i<dim_num;i++){
    _vol*=const_args_dims[const_args_pos][i]; if(const_args_prmn[const_args_pos][i] != i+1) s1_ind=0;
   }; vol=_vol; //total volume (number of tensor elements)
   if(s1_ind == 0){ //non-trivial permutation
// Set input/output permutations and dimension extents:
    if(drc == 0){ //normal index permutation
     for(i=0;i<dim_num;i++) o2n[i]=const_args_prmn[const_args_pos][i]-1; for(i=0;i<dim_num;i++) n2o[o2n[i]]=i;
    }else{ //inversed index permutation
     for(i=0;i<dim_num;i++) n2o[i]=const_args_prmn[const_args_pos][i]-1; for(i=0;i<dim_num;i++) o2n[n2o[i]]=i;
    }
    if(dmo == 0){ //normal dimension order
     for(i=0;i<dim_num;i++) dim_in[i]=const_args_dims[const_args_pos][i];
     for(i=0;i<dim_num;i++) dim_out[o2n[i]]=dim_in[i];
    }else{ //inversed dimension order
     for(i=0;i<dim_num;i++) dim_out[i]=const_args_dims[const_args_pos][i];
     for(i=0;i<dim_num;i++) dim_in[n2o[i]]=dim_out[i];
    }
    s1_step=dim_in[s1_ind]; s2_step=dim_in[s2_ind];
    if(_vol > TENS_TRANSP_BUF_SIZE){ //tensor block does not fit into the shared memory buffer
// Determine the input/output minor index sets and the combined minor index set:
     l=(int)(sqrt((float)TENS_TRANSP_BUF_SIZE));
     minor_in=0; _vol_in=1; for(i=0;i<dim_num;i++){j=_vol_in*dim_in[i]; if(j>l) break; minor_in++; _vol_in=j;}
     minor_out=0; _vol_out=1; for(i=0;i<dim_num;i++){j=_vol_out*dim_out[i]; if(j>l) break; minor_out++; _vol_out=j;}
     minor=minor_in; _vol_minor=_vol_in; for(i=0;i<minor_out;i++){if(n2o[i]>=minor_in){minor++; _vol_minor*=dim_out[i];}}
     m=1; _s1=0; _s2=0;
     while(_vol_minor < TENS_TRANSP_BUF_SIZE && m != 0){
      m=0;
      if(_s1 == 0){for(i=minor_in;i<dim_num;i++){if(o2n[i]<minor_out){minor_in++; _vol_in*=dim_in[i];}else{break;}}}
      if(_s2 == 0){for(i=minor_out;i<dim_num;i++){if(n2o[i]<minor_in){minor_out++; _vol_out*=dim_out[i];}else{break;}}}
      j=dim_in[minor_in]; l=dim_out[minor_out];
      if(minor_in == n2o[minor_out] && _s1+_s2 == 0){ //same candidate index to both the input and output index sets
       if(j > 1 && TENS_TRANSP_BUF_SIZE < _vol_minor*2) break;
       if(_vol_minor*j > TENS_TRANSP_BUF_SIZE){s1_ind=minor_in; s1_step=TENS_TRANSP_BUF_SIZE/_vol_minor; _s1++; _s2++;}
       minor_in++; _vol_in*=j; minor_out++; _vol_out*=j; minor++; _vol_minor*=j; m++;
      }else{ //the input and output index sets consider two different candidates
       if(_vol_minor*j*l <= TENS_TRANSP_BUF_SIZE && _s1+_s2 == 0){ //accept both, no splitting
        minor_in++; _vol_in*=j; minor_out++; _vol_out*=l; minor+=2; _vol_minor*=(j*l); m++;
       }else{ //try to accept either one of the two OR both with splitting
        if(j == 1 || l == 1){
         if(j == 1 && _s1 == 0){minor_in++; minor++; m++;}
         if(l == 1 && _s2 == 0){minor_out++; minor++; m++;}
        }else{
         if(_vol_minor*j <= TENS_TRANSP_BUF_SIZE && _vol_minor*l > TENS_TRANSP_BUF_SIZE &&
            _vol_out >= warpSize && _s1 == 0){ //accept the input index, no splitting
          minor_in++; _vol_in*=j; minor++; _vol_minor*=j; m++;
         }else if(_vol_minor*j > TENS_TRANSP_BUF_SIZE && _vol_minor*l <= TENS_TRANSP_BUF_SIZE &&
                  _vol_in >= warpSize && _s2 == 0){ //accept the output index, no splitting
          minor_out++; _vol_out*=l; minor++; _vol_minor*=l; m++;
         }else{ //splitting is unavoidable (both OR one OR none)
          if(TENS_TRANSP_BUF_SIZE >= _vol_minor*2){
           if(j >= 4 && l >= 4){ //dimension extents are large enough to be split
            if(_vol_minor*4 > TENS_TRANSP_BUF_SIZE){ //impossible to split both indices
             if(_vol_in <= _vol_out && _s1 == 0){ //split the input candidate index
              s1_ind=minor_in; s1_step=TENS_TRANSP_BUF_SIZE/_vol_minor;
              minor_in++; _vol_in*=j; minor++; _vol_minor*=j; _s1++; m++;
             }else{ //split the output candidate index
              if(_s2 == 0){
               s1_ind=n2o[minor_out]; s1_step=TENS_TRANSP_BUF_SIZE/_vol_minor;
               minor_out++; _vol_out*=l; minor++; _vol_minor*=l; _s2++; m++;
              }
             }
            }else{ //possible to split both indices
             i=(int)sqrt(((float)TENS_TRANSP_BUF_SIZE)/(float)_vol_minor); if(i < 2) i=2; //uniform splitting
             s1_step=i; s2_step=i; val=(float)_vol_out/(float)_vol_in;
             if(val < 1.0f){ //scale the initial uniform splitting to reflect the disbalance between _vol_in and _vol_out
              if(val*(float)i < 1.0f) val=1.0f/(float)i; if(val*(float)l < (float)i) val=(float)i/(float)l;
             }else{
              if(val*(float)i > (float)j) val=(float)j/(float)i; if(val > float(i)) val=(float)i;
             }
             s1_step=(int)(((float)i)*val); s2_step=(int)(((float)i)/val);
             if(s1_step >= 2 && _s1 == 0){ //&& s1_step <= dim_in[minor_in]
              s1_ind=minor_in; minor_in++; _vol_in*=j; minor++; _vol_minor*=j; _s1++; m++;
             }else{
              s1_step=dim_in[s1_ind];
             }
             if(s2_step >= 2 && _s2 == 0){ //&& s2_step <= dim_out[minor_out]
              s2_ind=n2o[minor_out]; minor_out++; _vol_out*=l; minor++; _vol_minor*=l; _s2++; m++;
             }else{
              s2_step=dim_in[s2_ind];
             }
            }
           }else if(j >= 4 && l < 4 && _s1 == 0){ //split the input candidate index
            s1_ind=minor_in; s1_step=TENS_TRANSP_BUF_SIZE/_vol_minor;
            minor_in++; _vol_in*=j; minor++; _vol_minor*=j; _s1++; m++;
           }else if(j < 4 && l >= 4 && _s2 == 0){ //split the output candidate index
            s1_ind=n2o[minor_out]; s1_step=TENS_TRANSP_BUF_SIZE/_vol_minor;
            minor_out++; _vol_out*=l; minor++; _vol_minor*=l; _s2++; m++;
           }else{ //both candidate indices have too small extent to be split: try to add one of them fully
            if(_vol_minor*j <= TENS_TRANSP_BUF_SIZE && _s1 == 0){
             minor_in++; _vol_in*=j; minor++; _vol_minor*=j; m++;
            }else if(_vol_minor*l <= TENS_TRANSP_BUF_SIZE && _s2 == 0){
             minor_out++; _vol_out*=l; minor++; _vol_minor*=l; m++;
            }
           }
          }else{ //unable to add more indices in the minor set
           break;
          }
         }
        }
       }
      }
     }
     if(s1_ind == dim_num-1 && s2_ind == dim_num-1){s2_ind=0; s2_step=dim_in[0];} //s1_ind was set while s2_ind was not
    }else{ //tensor block fits into the shared memory buffer from the beginning
     minor=dim_num; minor_in=dim_num; minor_out=dim_num; _vol_minor=_vol; _vol_in=_vol; _vol_out=_vol;
    }
// Share the tensor transpose configuration with other threads in each block:
    vol_ext=_vol/_vol_minor; s1_dim=dim_in[s1_ind]; s2_dim=dim_in[s2_ind];
// Set indexing bases (OUT:{out,in_c,ext_in}_new; IN:{in,out_c,ext_in}_old):
//  OUTPUT indexing (dim_out[], base_out[]: prioritized new numeration):
    for(i=0;i<dim_num;i++){tmp0[i]=dim_out[i];} //save output dimension extents (new numeration)
    j=0; for(i=0;i<minor_out;i++){pri[j++]=i;} //output minor index set (new numeration))
    for(i=0;i<dim_num;i++){if(o2n[i]>=minor_out) pri[j++]=o2n[i];} //{compl.input minor + external} index set (new numeration)
    j=1; for(i=0;i<dim_num;i++){dim_out[i]=j; j*=tmp0[i];} //output bases (new numeration)
    for(i=0;i<dim_num;i++){base_out[i]=dim_out[pri[i]];} //output bases (prioritized new numeration)
    for(i=0;i<dim_num;i++){dim_out[i]=tmp0[pri[i]];} //output extents (prioritized new numeration)
    for(i=0;i<dim_num;i++){if(n2o[pri[i]]==s1_ind){s1_ond=i;}else if(n2o[pri[i]]==s2_ind){s2_ond=i;}} //split indices (prioritized new numeration)
//  INPUT indexing (dim_in[], base_in[]: prioritized old numeration):
    for(i=0;i<dim_num;i++){tmp0[i]=dim_in[i];} //save input dimension extents (old numeration)
    j=0; for(i=0;i<minor_in;i++){pri[j++]=i;} //input minor index set (old numeration)
    for(i=0;i<minor_out;i++){if(n2o[i]>=minor_in) pri[j++]=n2o[i];} //compl.output minor idex set (old numeration)
    for(i=j;i<dim_num;i++){pri[i]=n2o[pri[i]];} //external index set (just convert new numbers to old ones for consistency)
    j=1; for(i=0;i<dim_num;i++){dim_in[i]=j; j*=tmp0[i];} //input bases (old numeration)
    for(i=0;i<dim_num;i++){base_in[i]=dim_in[pri[i]];} //input bases (prioritized old numeration)
    for(i=0;i<dim_num;i++){dim_in[i]=tmp0[pri[i]];} //input extents (prioritized old numeration)
    for(i=0;i<dim_num;i++){if(pri[i]==s1_ind){_s1=i;}else if(pri[i]==s2_ind){_s2=i;}} //split indices (prioritized old numeration)
    s1_ind=_s1; s2_ind=_s2;
    ns1=1+(s1_dim-1)/s1_step; //number of segments from the 1st split minor index
    ns2=1+(s2_dim-1)/s2_step; //number of segments from the 2nd split minor index
//  Index position correspondence for the minor index set (pri-new --> pri-old):
    j=0; for(i=0;i<minor_out;i++){if(n2o[i]<minor_in){pri[i]=n2o[i];}else{pri[i]=(minor_in+j); j++;}}
    j=0; for(i=0;i<minor_in;i++){if(o2n[i]<minor_out){pri[o2n[i]]=i;}else{pri[minor_out+j]=i; j++;}}
// Check tensor transpose configuration parameters:
    if(minor <= 0 || minor_in <= 0 || minor_out <= 0 || _vol <= 0 || _vol_minor <= 0) err_code+=5000; //trap
    if(s1_ind >= dim_num || s2_ind >= dim_num || s1_ond >= dim_num || s2_ond >= dim_num ||
       s1_ind == s2_ind || s1_ond == s2_ond || s1_step <= 0 || s2_step <= 0) err_code+=1000; //trap
    if((s1_step != dim_in[s1_ind] && s1_ind != minor_in-1 && s1_ond != minor_out-1) ||
       (s2_step != dim_in[s2_ind] && s2_ind != minor_in-1 && s2_ond != minor_out-1)) err_code+=500; //trap
    if((_vol_minor*s1_step*s2_step)/(s1_dim*s2_dim) > TENS_TRANSP_BUF_SIZE) err_code+=100; //trap
   } //endif: non-trivial permutation
  }else{
   err_code=1+2*blockDim.x%warpSize;
  }
 } //endif: Master thread.
#ifdef DEBUG_GPU
//DEBUG RECORD begin:
 if(blockIdx.x == 0 && threadIdx.x == 0){
  j=0; gpu_debug_dump[j++]=dim_num;
  for(i=0;i<dim_num;i++) gpu_debug_dump[j++]=const_args_dims[const_args_pos][i];
  for(i=0;i<dim_num;i++) gpu_debug_dump[j++]=const_args_prmn[const_args_pos][i];
  for(i=0;i<dim_num;i++) gpu_debug_dump[j++]=base_in[i];
  for(i=0;i<dim_num;i++) gpu_debug_dump[j++]=base_out[i];
  gpu_debug_dump[j++]=vol; gpu_debug_dump[j++]=vol_ext; gpu_debug_dump[j++]=vol/vol_ext;
  gpu_debug_dump[j++]=minor; gpu_debug_dump[j++]=minor_in; gpu_debug_dump[j++]=minor_out;
  gpu_debug_dump[j++]=s1_ind; gpu_debug_dump[j++]=s1_ond; gpu_debug_dump[j++]=s1_step; gpu_debug_dump[j++]=s1_dim;
  gpu_debug_dump[j++]=s2_ind; gpu_debug_dump[j++]=s2_ond; gpu_debug_dump[j++]=s2_step; gpu_debug_dump[j++]=s2_dim;
  for(i=0;i<dim_num;i++) gpu_debug_dump[j++]=pri[i];
  gpu_debug_dump[j++]=err_code; gpu_debug_dump[j++]=-1;
 }
//DEBUG RECORD end.
#endif
 __syncthreads();

//Proceed:
 if(err_code == 0){
  if(s1_ind > dim_num){ //tag of a trivial permutation
// Direct copy:
   _vol=vol; j=gridDim.x*blockDim.x; i=blockIdx.x*blockDim.x+threadIdx.x; _addr_in=_vol-_vol%j;
   for(_addr=0;_addr<_addr_in;_addr+=j){
    _addr_out=_addr+i; tens_out[_addr_out]=tens_in[_addr_out];
   }
   _addr_out=_addr_in+i; if(_addr_out<_vol) tens_out[_addr_out]=tens_in[_addr_out];
  }else{ //non-trivial permutation
   l=threadIdx.x/warpSize; //l: warp number
// Distribute work accross CUDA blocks (external multi-index + splitting):
   for(_work_piece=blockIdx.x;_work_piece<vol_ext*ns1*ns2;_work_piece+=gridDim.x){ //(ns1*ns2*vol_ext) is the total number of independent tasks
    _addr=_work_piece; _addr/=vol_ext; _vol=_work_piece-_addr*vol_ext; _s2=(int)(_addr/ns1); _s1=(int)(_addr-_s2*ns1); //{_addr_ext,_s1,_s2} --> tensor subblock (CUDA block)
//  Modify dimension extents due to possible dimension splitting:
    if(threadIdx.x == 0){
     if(_s1+1 == ns1){ //last segment of the 1st split index
      j=s1_dim-_s1*s1_step; dim_in[s1_ind]=j; dim_out[s1_ond]=j;
     }else{ //internal segment of the 1st split index
      dim_in[s1_ind]=s1_step; dim_out[s1_ond]=s1_step;
     }
     if(_s2+1 == ns2){ //last segment of the 2nd split index
      j=s2_dim-_s2*s2_step; dim_in[s2_ind]=j; dim_out[s2_ond]=j;
     }else{ //internal segment of the 2nd split index
      dim_in[s2_ind]=s2_step; dim_out[s2_ond]=s2_step;
     }
     j=1; for(i=0;i<minor;i++){tmp0[i]=j; j*=dim_in[i];} //minor buffer bases (pri-old)
     for(i=0;i<minor;i++) n2o[i]=tmp0[pri[i]]; //look up table to accelerate further accesses to tmp0[]
    }
    __syncthreads();
//  Mount input/output volumes and bases:
    _vol_in=dim_in[0]; for(i=1;i<minor_in;i++){_vol_in*=dim_in[i];}
    _vol_out=dim_out[0]; for(i=1;i<minor_out;i++){_vol_out*=dim_out[i];}
    _vol_minor=_vol_out; for(i=minor_out;i<minor;i++){_vol_minor*=dim_out[i];}
    _addr_in=(_s1*s1_step)*base_in[s1_ind]+(_s2*s2_step)*base_in[s2_ind]; _addr_out=_vol;
    for(i=minor;i<dim_num;i++){_addr=_vol/dim_in[i]; _addr_in+=(_vol-_addr*dim_in[i])*base_in[i]; _vol=_addr;}
    _vol=_addr_out; _addr_out=(_s1*s1_step)*base_out[s1_ond]+(_s2*s2_step)*base_out[s2_ond];
    for(i=minor;i<dim_num;i++){_addr=_vol/dim_out[i]; _addr_out+=(_vol-_addr*dim_out[i])*base_out[i]; _vol=_addr;}
    if(_vol_out > TENS_TRANSP_TAB_SIZE || _vol_minor > _vol_in*TENS_TRANSP_TAB_SIZE ||
       _vol_minor > _vol_out*TENS_TRANSP_TAB_SIZE){
//  Algorithm 0 (slower):
//   Read the minor volume into the buffer from the input tensor block:
     _vol_minor/=_vol_in; //vol_in_c
     _s1=1+(_vol_in-1)/warpSize; //number of warps (lines) which fully cover the input volume
     _s2=blockDim.x/warpSize; //number of whole warps in a thread block (each warp treats one line)
     for(j=l;j<_s1*_vol_minor;j+=_s2){ //j: Line number
      m=j/_s1; _addr=_addr_in; n=m; //n: Input column number (in_c)
      for(i=minor_in;i<minor;i++){k=m/dim_in[i]; _addr+=(m-k*dim_in[i])*base_in[i]; m=k;}
//    m=(j%_s1)*warpSize+threadIdx.x%warpSize; //elemental offset in the input volume
      m=threadIdx.x+(j-n*_s1-l)*warpSize; //elemental offset in the input volume (alternative)
      if(m < _vol_in){buf0[n*_vol_in+m]=tens_in[_addr+m];}
     }
     __syncthreads();
//   Write the minor volume from the buffer into the output tensor block:
     _vol_minor=(_vol_minor*_vol_in)/_vol_out; //vol_out_c
     _s1=1+(_vol_out-1)/warpSize; //number of warps (lines) which fully cover the output volume
     for(j=l;j<_s1*_vol_minor;j+=_s2){ //j: Line number
      n=j/_s1; _addr=_addr_out; _vol=n; _vol_in=0; //_vol: Output column number (out_c)
//    for(i=minor_out;i<minor;i++){m=n%dim_out[i]; n/=dim_out[i]; _addr+=m*base_out[i]; _vol_in+=m*tmp0[pri[i]];}
      for(i=minor_out;i<minor;i++){k=n/dim_out[i]; m=n-k*dim_out[i]; n=k; _addr+=m*base_out[i]; _vol_in+=m*n2o[i];}
//    m=(j%_s1)*warpSize+threadIdx.x%warpSize; //elemental offset in the output volume
      m=threadIdx.x+(j-(int)_vol*_s1-l)*warpSize; //elemental offset in the output volume (alternative)
      if(m < _vol_out){
       _addr+=m;
//     for(i=0;i<minor_out;i++){_vol_in+=(m%dim_out[i])*tmp0[pri[i]]; m/=dim_out[i];}
       for(i=0;i<minor_out;i++){k=m/dim_out[i]; _vol_in+=(m-k*dim_out[i])*n2o[i]; m=k;}
       tens_out[_addr]=buf0[_vol_in];
      }
     }
     __syncthreads();
    }else{
//  Algorithm 1 (presumably faster):
//   Create per-block look-up tables:
     m=_vol_minor/_vol_in; //vol_in_c
     for(j=threadIdx.x;j<m;j+=blockDim.x){ //column number (input)
      _addr=0; _s1=j;
//    for(i=minor_in;i<minor;i++){_addr+=(_s1%dim_in[i])*base_in[i]; _s1/=dim_in[i];}
      for(i=minor_in;i<minor;i++){_s2=_s1/dim_in[i]; _addr+=(_s1-_s2*dim_in[i])*base_in[i]; _s1=_s2;}
      ftb[j]=_addr;
     }
     m=_vol_minor/_vol_out; //vol_out_c
     for(j=threadIdx.x;j<m;j+=blockDim.x){ //column number (output)
      _addr=0; _s1=j;
//    for(i=minor_out;i<minor;i++){_addr+=(_s1%dim_out[i])*base_out[i]; _s1/=dim_out[i];}
      for(i=minor_out;i<minor;i++){_s2=_s1/dim_out[i]; _addr+=(_s1-_s2*dim_out[i])*base_out[i]; _s1=_s2;}
      gtb[j]=_addr;
     }
     for(j=threadIdx.x;j<m;j+=blockDim.x){ //column number (output)
      n=0; _s1=j;
//    for(i=minor_out;i<minor;i++){n+=(_s1%dim_out[i])*n2o[i]; _s1/=dim_out[i];}
      for(i=minor_out;i<minor;i++){_s2=_s1/dim_out[i]; n+=(_s1-_s2*dim_out[i])*n2o[i]; _s1=_s2;}
      htb[j]=n;
     }
     for(j=threadIdx.x;j<_vol_out;j+=blockDim.x){
      n=0; _s1=j;
//    for(i=0;i<minor_out;i++){n+=(_s1%dim_out[i])*n2o[i]; _s1/=dim_out[i];}
      for(i=0;i<minor_out;i++){_s2=_s1/dim_out[i]; n+=(_s1-_s2*dim_out[i])*n2o[i]; _s1=_s2;}
      stb[j]=n;
     }
     __syncthreads();
//   Read the minor volume into the buffer from the input tensor block:
     _vol_minor/=_vol_in; //vol_in_c
     _s1=1+(_vol_in-1)/warpSize; //number of warps (lines) which fully cover the input volume
     _s2=blockDim.x/warpSize; //number of whole warps in a thread block (each warp treats one line)
     for(j=l;j<_s1*_vol_minor;j+=_s2){ //j: Line number
      m=j/_s1; n=threadIdx.x+(j-m*_s1-l)*warpSize; //m: Input column number (in_c); n: Offset in the column
      if(n < _vol_in){_addr=_addr_in+ftb[m]+n; buf0[m*_vol_in+n]=tens_in[_addr];}
     }
     __syncthreads();
//   Write the minor volume from the buffer into the output tensor block:
     _vol_minor=(_vol_minor*_vol_in)/_vol_out; //vol_out_c
     _s1=1+(_vol_out-1)/warpSize; //number of warps (lines) which fully cover the output volume
     for(j=l;j<_s1*_vol_minor;j+=_s2){ //j: Line number
      m=j/_s1; n=threadIdx.x+(j-m*_s1-l)*warpSize; //m: Output column number (out_c); n: Offset in the column
      if(n < _vol_out){_addr=_addr_out+gtb[m]+n; _vol_in=htb[m]+stb[n]; tens_out[_addr]=buf0[_vol_in];}
     }
     __syncthreads();
    }
   } //enddo _work_piece: independent work distribution among thread blocks
  }
 }

//Record errors if occured (for each block):
 if(threadIdx.x == 0){if(err_code != 0) i=atomicAdd(&gpu_error_count,1);}
 return;
}
//---------------------------------------------------------------------------------------------------------------
// TENSOR TRANSPOSE (R8) (shared-memory version):
__global__ void gpu_tensor_block_copy_dlf_r8__(int dmo, int drc, int dim_num, int const_args_pos,
                                               const double* __restrict__ tens_in, double* __restrict__ tens_out)
/**
Shared-memory version of tensor transpose: tens_out=TRN(tens_in):
INPUT:
 # dmo - dimension extents order (0: normal, as it is in <const_args>; not 0: permuted dimension order will be imposed);
 # drc - index permutation direction (0: normal, as it is in <const_args>; not 0: inversed permutation will be used);
 # dim_num - tensor block rank;
 # const_args_pos - entry in the __constant__ memory bank where tensor block dimension extents (const_args_dims)
                    and index permutation (const_args_prmn) are stored;
 # tens_in[0:] - input tensor;
OUTPUT:
 # tens_out[0:] - output (transposed) tensor;
NOTES:
 # Minimal CUDA execution configuration is <<<1,warpSize>>>
 # Number of threads per block must be multiple of the warpSize!
**/
{
 __shared__ double buf0[TENS_TRANSP_BUF_SIZE];
 __shared__ float val;
 __shared__ size_t base_in[MAX_TENSOR_RANK],base_out[MAX_TENSOR_RANK];
 __shared__ size_t ftb[TENS_TRANSP_TAB_SIZE],gtb[TENS_TRANSP_TAB_SIZE];
 __shared__ int htb[TENS_TRANSP_TAB_SIZE],stb[TENS_TRANSP_TAB_SIZE];
 __shared__ int dim_in[MAX_TENSOR_RANK],dim_out[MAX_TENSOR_RANK],o2n[MAX_TENSOR_RANK],n2o[MAX_TENSOR_RANK];
 __shared__ int pri[MAX_TENSOR_RANK],tmp0[MAX_TENSOR_RANK];
 __shared__ int err_code,minor,minor_in,minor_out,s1_ind,s1_ond,s1_step,s1_dim,s2_ind,s2_ond,s2_step,s2_dim,ns1,ns2;
 __shared__ size_t vol,vol_ext;
 size_t _vol,_addr_in,_addr_out,_addr,_work_piece;
 int i,j,k,l,m,n,_vol_minor,_vol_in,_vol_out,_s1,_s2;
/*
SHARED MEMORY USE =
 + TENS_TRANSP_BUF_SIZE*8
 + MAX_TENSOR_RANK*(8+8+4+4+4+4+4+4)
 + TENS_TRANSP_TAB_SIZE*(8+8+4+4)
 + 4*15 + 8*2
REGISTER USE =
 + 4*4 + 4*11 + 8*5 = 100 Bytes
*/

//Determine the minor index set (only the master thread in each thread block):
 if(threadIdx.x == 0){
  err_code=0;
  if(dim_num >= 0 && dim_num <= MAX_TENSOR_RANK && blockDim.x >= warpSize && blockDim.x%warpSize == 0){
   s1_ind=dim_num+1; s2_ind=dim_num-1;
   _vol=1; for(i=0;i<dim_num;i++){
    _vol*=const_args_dims[const_args_pos][i]; if(const_args_prmn[const_args_pos][i] != i+1) s1_ind=0;
   }; vol=_vol; //total volume (number of tensor elements)
   if(s1_ind == 0){ //non-trivial permutation
// Set input/output permutations and dimension extents:
    if(drc == 0){ //normal index permutation
     for(i=0;i<dim_num;i++) o2n[i]=const_args_prmn[const_args_pos][i]-1; for(i=0;i<dim_num;i++) n2o[o2n[i]]=i;
    }else{ //inversed index permutation
     for(i=0;i<dim_num;i++) n2o[i]=const_args_prmn[const_args_pos][i]-1; for(i=0;i<dim_num;i++) o2n[n2o[i]]=i;
    }
    if(dmo == 0){ //normal dimension order
     for(i=0;i<dim_num;i++) dim_in[i]=const_args_dims[const_args_pos][i];
     for(i=0;i<dim_num;i++) dim_out[o2n[i]]=dim_in[i];
    }else{ //inversed dimension order
     for(i=0;i<dim_num;i++) dim_out[i]=const_args_dims[const_args_pos][i];
     for(i=0;i<dim_num;i++) dim_in[n2o[i]]=dim_out[i];
    }
    s1_step=dim_in[s1_ind]; s2_step=dim_in[s2_ind];
    if(_vol > TENS_TRANSP_BUF_SIZE){ //tensor block does not fit into the shared memory buffer
// Determine the input/output minor index sets and the combined minor index set:
     l=(int)(sqrt((float)TENS_TRANSP_BUF_SIZE));
     minor_in=0; _vol_in=1; for(i=0;i<dim_num;i++){j=_vol_in*dim_in[i]; if(j>l) break; minor_in++; _vol_in=j;}
     minor_out=0; _vol_out=1; for(i=0;i<dim_num;i++){j=_vol_out*dim_out[i]; if(j>l) break; minor_out++; _vol_out=j;}
     minor=minor_in; _vol_minor=_vol_in; for(i=0;i<minor_out;i++){if(n2o[i]>=minor_in){minor++; _vol_minor*=dim_out[i];}}
     m=1; _s1=0; _s2=0;
     while(_vol_minor < TENS_TRANSP_BUF_SIZE && m != 0){
      m=0;
      if(_s1 == 0){for(i=minor_in;i<dim_num;i++){if(o2n[i]<minor_out){minor_in++; _vol_in*=dim_in[i];}else{break;}}}
      if(_s2 == 0){for(i=minor_out;i<dim_num;i++){if(n2o[i]<minor_in){minor_out++; _vol_out*=dim_out[i];}else{break;}}}
      j=dim_in[minor_in]; l=dim_out[minor_out];
      if(minor_in == n2o[minor_out] && _s1+_s2 == 0){ //same candidate index to both the input and output index sets
       if(j > 1 && TENS_TRANSP_BUF_SIZE < _vol_minor*2) break;
       if(_vol_minor*j > TENS_TRANSP_BUF_SIZE){s1_ind=minor_in; s1_step=TENS_TRANSP_BUF_SIZE/_vol_minor; _s1++; _s2++;}
       minor_in++; _vol_in*=j; minor_out++; _vol_out*=j; minor++; _vol_minor*=j; m++;
      }else{ //the input and output index sets consider two different candidates
       if(_vol_minor*j*l <= TENS_TRANSP_BUF_SIZE && _s1+_s2 == 0){ //accept both, no splitting
        minor_in++; _vol_in*=j; minor_out++; _vol_out*=l; minor+=2; _vol_minor*=(j*l); m++;
       }else{ //try to accept either one of the two OR both with splitting
        if(j == 1 || l == 1){
         if(j == 1 && _s1 == 0){minor_in++; minor++; m++;}
         if(l == 1 && _s2 == 0){minor_out++; minor++; m++;}
        }else{
         if(_vol_minor*j <= TENS_TRANSP_BUF_SIZE && _vol_minor*l > TENS_TRANSP_BUF_SIZE &&
            _vol_out >= warpSize && _s1 == 0){ //accept the input index, no splitting
          minor_in++; _vol_in*=j; minor++; _vol_minor*=j; m++;
         }else if(_vol_minor*j > TENS_TRANSP_BUF_SIZE && _vol_minor*l <= TENS_TRANSP_BUF_SIZE &&
                  _vol_in >= warpSize && _s2 == 0){ //accept the output index, no splitting
          minor_out++; _vol_out*=l; minor++; _vol_minor*=l; m++;
         }else{ //splitting is unavoidable (both OR one OR none)
          if(TENS_TRANSP_BUF_SIZE >= _vol_minor*2){
           if(j >= 4 && l >= 4){ //dimension extents are large enough to be split
            if(_vol_minor*4 > TENS_TRANSP_BUF_SIZE){ //impossible to split both indices
             if(_vol_in <= _vol_out && _s1 == 0){ //split the input candidate index
              s1_ind=minor_in; s1_step=TENS_TRANSP_BUF_SIZE/_vol_minor;
              minor_in++; _vol_in*=j; minor++; _vol_minor*=j; _s1++; m++;
             }else{ //split the output candidate index
              if(_s2 == 0){
               s1_ind=n2o[minor_out]; s1_step=TENS_TRANSP_BUF_SIZE/_vol_minor;
               minor_out++; _vol_out*=l; minor++; _vol_minor*=l; _s2++; m++;
              }
             }
            }else{ //possible to split both indices
             i=(int)sqrt(((float)TENS_TRANSP_BUF_SIZE)/(float)_vol_minor); if(i < 2) i=2; //uniform splitting
             s1_step=i; s2_step=i; val=(float)_vol_out/(float)_vol_in;
             if(val < 1.0f){ //scale the initial uniform splitting to reflect the disbalance between _vol_in and _vol_out
              if(val*(float)i < 1.0f) val=1.0f/(float)i; if(val*(float)l < (float)i) val=(float)i/(float)l;
             }else{
              if(val*(float)i > (float)j) val=(float)j/(float)i; if(val > float(i)) val=(float)i;
             }
             s1_step=(int)(((float)i)*val); s2_step=(int)(((float)i)/val);
             if(s1_step >= 2 && _s1 == 0){ //&& s1_step <= dim_in[minor_in]
              s1_ind=minor_in; minor_in++; _vol_in*=j; minor++; _vol_minor*=j; _s1++; m++;
             }else{
              s1_step=dim_in[s1_ind];
             }
             if(s2_step >= 2 && _s2 == 0){ //&& s2_step <= dim_out[minor_out]
              s2_ind=n2o[minor_out]; minor_out++; _vol_out*=l; minor++; _vol_minor*=l; _s2++; m++;
             }else{
              s2_step=dim_in[s2_ind];
             }
            }
           }else if(j >= 4 && l < 4 && _s1 == 0){ //split the input candidate index
            s1_ind=minor_in; s1_step=TENS_TRANSP_BUF_SIZE/_vol_minor;
            minor_in++; _vol_in*=j; minor++; _vol_minor*=j; _s1++; m++;
           }else if(j < 4 && l >= 4 && _s2 == 0){ //split the output candidate index
            s1_ind=n2o[minor_out]; s1_step=TENS_TRANSP_BUF_SIZE/_vol_minor;
            minor_out++; _vol_out*=l; minor++; _vol_minor*=l; _s2++; m++;
           }else{ //both candidate indices have too small extent to be split: try to add one of them fully
            if(_vol_minor*j <= TENS_TRANSP_BUF_SIZE && _s1 == 0){
             minor_in++; _vol_in*=j; minor++; _vol_minor*=j; m++;
            }else if(_vol_minor*l <= TENS_TRANSP_BUF_SIZE && _s2 == 0){
             minor_out++; _vol_out*=l; minor++; _vol_minor*=l; m++;
            }
           }
          }else{ //unable to add more indices in the minor set
           break;
          }
         }
        }
       }
      }
     }
     if(s1_ind == dim_num-1 && s2_ind == dim_num-1){s2_ind=0; s2_step=dim_in[0];} //s1_ind was set while s2_ind was not
    }else{ //tensor block fits into the shared memory buffer from the beginning
     minor=dim_num; minor_in=dim_num; minor_out=dim_num; _vol_minor=_vol; _vol_in=_vol; _vol_out=_vol;
    }
// Share the tensor transpose configuration with other threads in each block:
    vol_ext=_vol/_vol_minor; s1_dim=dim_in[s1_ind]; s2_dim=dim_in[s2_ind];
// Set indexing bases (OUT:{out,in_c,ext_in}_new; IN:{in,out_c,ext_in}_old):
//  OUTPUT indexing (dim_out[], base_out[]: prioritized new numeration):
    for(i=0;i<dim_num;i++){tmp0[i]=dim_out[i];} //save output dimension extents (new numeration)
    j=0; for(i=0;i<minor_out;i++){pri[j++]=i;} //output minor index set (new numeration))
    for(i=0;i<dim_num;i++){if(o2n[i]>=minor_out) pri[j++]=o2n[i];} //{compl.input minor + external} index set (new numeration)
    j=1; for(i=0;i<dim_num;i++){dim_out[i]=j; j*=tmp0[i];} //output bases (new numeration)
    for(i=0;i<dim_num;i++){base_out[i]=dim_out[pri[i]];} //output bases (prioritized new numeration)
    for(i=0;i<dim_num;i++){dim_out[i]=tmp0[pri[i]];} //output extents (prioritized new numeration)
    for(i=0;i<dim_num;i++){if(n2o[pri[i]]==s1_ind){s1_ond=i;}else if(n2o[pri[i]]==s2_ind){s2_ond=i;}} //split indices (prioritized new numeration)
//  INPUT indexing (dim_in[], base_in[]: prioritized old numeration):
    for(i=0;i<dim_num;i++){tmp0[i]=dim_in[i];} //save input dimension extents (old numeration)
    j=0; for(i=0;i<minor_in;i++){pri[j++]=i;} //input minor index set (old numeration)
    for(i=0;i<minor_out;i++){if(n2o[i]>=minor_in) pri[j++]=n2o[i];} //compl.output minor idex set (old numeration)
    for(i=j;i<dim_num;i++){pri[i]=n2o[pri[i]];} //external index set (just convert new numbers to old ones for consistency)
    j=1; for(i=0;i<dim_num;i++){dim_in[i]=j; j*=tmp0[i];} //input bases (old numeration)
    for(i=0;i<dim_num;i++){base_in[i]=dim_in[pri[i]];} //input bases (prioritized old numeration)
    for(i=0;i<dim_num;i++){dim_in[i]=tmp0[pri[i]];} //input extents (prioritized old numeration)
    for(i=0;i<dim_num;i++){if(pri[i]==s1_ind){_s1=i;}else if(pri[i]==s2_ind){_s2=i;}} //split indices (prioritized old numeration)
    s1_ind=_s1; s2_ind=_s2;
    ns1=1+(s1_dim-1)/s1_step; //number of segments from the 1st split minor index
    ns2=1+(s2_dim-1)/s2_step; //number of segments from the 2nd split minor index
//  Index position correspondence for the minor index set (pri-new --> pri-old):
    j=0; for(i=0;i<minor_out;i++){if(n2o[i]<minor_in){pri[i]=n2o[i];}else{pri[i]=(minor_in+j); j++;}}
    j=0; for(i=0;i<minor_in;i++){if(o2n[i]<minor_out){pri[o2n[i]]=i;}else{pri[minor_out+j]=i; j++;}}
// Check tensor transpose configuration parameters:
    if(minor <= 0 || minor_in <= 0 || minor_out <= 0 || _vol <= 0 || _vol_minor <= 0) err_code+=5000; //trap
    if(s1_ind >= dim_num || s2_ind >= dim_num || s1_ond >= dim_num || s2_ond >= dim_num ||
       s1_ind == s2_ind || s1_ond == s2_ond || s1_step <= 0 || s2_step <= 0) err_code+=1000; //trap
    if((s1_step != dim_in[s1_ind] && s1_ind != minor_in-1 && s1_ond != minor_out-1) ||
       (s2_step != dim_in[s2_ind] && s2_ind != minor_in-1 && s2_ond != minor_out-1)) err_code+=500; //trap
    if((_vol_minor*s1_step*s2_step)/(s1_dim*s2_dim) > TENS_TRANSP_BUF_SIZE) err_code+=100; //trap
   } //endif: non-trivial permutation
  }else{
   err_code=1+2*blockDim.x%warpSize;
  }
 } //endif: Master thread.
#ifdef DEBUG_GPU
//DEBUG RECORD begin:
 if(blockIdx.x == 0 && threadIdx.x == 0){
  j=0; gpu_debug_dump[j++]=dim_num;
  for(i=0;i<dim_num;i++) gpu_debug_dump[j++]=const_args_dims[const_args_pos][i];
  for(i=0;i<dim_num;i++) gpu_debug_dump[j++]=const_args_prmn[const_args_pos][i];
  for(i=0;i<dim_num;i++) gpu_debug_dump[j++]=base_in[i];
  for(i=0;i<dim_num;i++) gpu_debug_dump[j++]=base_out[i];
  gpu_debug_dump[j++]=vol; gpu_debug_dump[j++]=vol_ext; gpu_debug_dump[j++]=vol/vol_ext;
  gpu_debug_dump[j++]=minor; gpu_debug_dump[j++]=minor_in; gpu_debug_dump[j++]=minor_out;
  gpu_debug_dump[j++]=s1_ind; gpu_debug_dump[j++]=s1_ond; gpu_debug_dump[j++]=s1_step; gpu_debug_dump[j++]=s1_dim;
  gpu_debug_dump[j++]=s2_ind; gpu_debug_dump[j++]=s2_ond; gpu_debug_dump[j++]=s2_step; gpu_debug_dump[j++]=s2_dim;
  for(i=0;i<dim_num;i++) gpu_debug_dump[j++]=pri[i];
  gpu_debug_dump[j++]=err_code; gpu_debug_dump[j++]=-1;
 }
//DEBUG RECORD end.
#endif
 __syncthreads();

//Proceed:
 if(err_code == 0){
  if(s1_ind > dim_num){ //tag of a trivial permutation
// Direct copy:
   _vol=vol; j=gridDim.x*blockDim.x; i=blockIdx.x*blockDim.x+threadIdx.x; _addr_in=_vol-_vol%j;
   for(_addr=0;_addr<_addr_in;_addr+=j){
    _addr_out=_addr+i; tens_out[_addr_out]=tens_in[_addr_out];
   }
   _addr_out=_addr_in+i; if(_addr_out<_vol) tens_out[_addr_out]=tens_in[_addr_out];
  }else{ //non-trivial permutation
   l=threadIdx.x/warpSize; //l: warp number
// Distribute work accross CUDA blocks (external multi-index + splitting):
   for(_work_piece=blockIdx.x;_work_piece<vol_ext*ns1*ns2;_work_piece+=gridDim.x){ //(ns1*ns2*vol_ext) is the total number of independent tasks
    _addr=_work_piece; _addr/=vol_ext; _vol=_work_piece-_addr*vol_ext; _s2=(int)(_addr/ns1); _s1=(int)(_addr-_s2*ns1); //{_addr_ext,_s1,_s2} --> tensor subblock (CUDA block)
//  Modify dimension extents due to possible dimension splitting:
    if(threadIdx.x == 0){
     if(_s1+1 == ns1){ //last segment of the 1st split index
      j=s1_dim-_s1*s1_step; dim_in[s1_ind]=j; dim_out[s1_ond]=j;
     }else{ //internal segment of the 1st split index
      dim_in[s1_ind]=s1_step; dim_out[s1_ond]=s1_step;
     }
     if(_s2+1 == ns2){ //last segment of the 2nd split index
      j=s2_dim-_s2*s2_step; dim_in[s2_ind]=j; dim_out[s2_ond]=j;
     }else{ //internal segment of the 2nd split index
      dim_in[s2_ind]=s2_step; dim_out[s2_ond]=s2_step;
     }
     j=1; for(i=0;i<minor;i++){tmp0[i]=j; j*=dim_in[i];} //minor buffer bases (pri-old)
     for(i=0;i<minor;i++) n2o[i]=tmp0[pri[i]]; //look up table to accelerate further accesses to tmp0[]
    }
    __syncthreads();
//  Mount input/output volumes and bases:
    _vol_in=dim_in[0]; for(i=1;i<minor_in;i++){_vol_in*=dim_in[i];}
    _vol_out=dim_out[0]; for(i=1;i<minor_out;i++){_vol_out*=dim_out[i];}
    _vol_minor=_vol_out; for(i=minor_out;i<minor;i++){_vol_minor*=dim_out[i];}
    _addr_in=(_s1*s1_step)*base_in[s1_ind]+(_s2*s2_step)*base_in[s2_ind]; _addr_out=_vol;
    for(i=minor;i<dim_num;i++){_addr=_vol/dim_in[i]; _addr_in+=(_vol-_addr*dim_in[i])*base_in[i]; _vol=_addr;}
    _vol=_addr_out; _addr_out=(_s1*s1_step)*base_out[s1_ond]+(_s2*s2_step)*base_out[s2_ond];
    for(i=minor;i<dim_num;i++){_addr=_vol/dim_out[i]; _addr_out+=(_vol-_addr*dim_out[i])*base_out[i]; _vol=_addr;}
    if(_vol_out > TENS_TRANSP_TAB_SIZE || _vol_minor > _vol_in*TENS_TRANSP_TAB_SIZE ||
       _vol_minor > _vol_out*TENS_TRANSP_TAB_SIZE){
//  Algorithm 0 (slower):
//   Read the minor volume into the buffer from the input tensor block:
     _vol_minor/=_vol_in; //vol_in_c
     _s1=1+(_vol_in-1)/warpSize; //number of warps (lines) which fully cover the input volume
     _s2=blockDim.x/warpSize; //number of whole warps in a thread block (each warp treats one line)
     for(j=l;j<_s1*_vol_minor;j+=_s2){ //j: Line number
      m=j/_s1; _addr=_addr_in; n=m; //n: Input column number (in_c)
      for(i=minor_in;i<minor;i++){k=m/dim_in[i]; _addr+=(m-k*dim_in[i])*base_in[i]; m=k;}
//    m=(j%_s1)*warpSize+threadIdx.x%warpSize; //elemental offset in the input volume
      m=threadIdx.x+(j-n*_s1-l)*warpSize; //elemental offset in the input volume (alternative)
      if(m < _vol_in){buf0[n*_vol_in+m]=tens_in[_addr+m];}
     }
     __syncthreads();
//   Write the minor volume from the buffer into the output tensor block:
     _vol_minor=(_vol_minor*_vol_in)/_vol_out; //vol_out_c
     _s1=1+(_vol_out-1)/warpSize; //number of warps (lines) which fully cover the output volume
     for(j=l;j<_s1*_vol_minor;j+=_s2){ //j: Line number
      n=j/_s1; _addr=_addr_out; _vol=n; _vol_in=0; //_vol: Output column number (out_c)
//    for(i=minor_out;i<minor;i++){m=n%dim_out[i]; n/=dim_out[i]; _addr+=m*base_out[i]; _vol_in+=m*tmp0[pri[i]];}
      for(i=minor_out;i<minor;i++){k=n/dim_out[i]; m=n-k*dim_out[i]; n=k; _addr+=m*base_out[i]; _vol_in+=m*n2o[i];}
//    m=(j%_s1)*warpSize+threadIdx.x%warpSize; //elemental offset in the output volume
      m=threadIdx.x+(j-(int)_vol*_s1-l)*warpSize; //elemental offset in the output volume (alternative)
      if(m < _vol_out){
       _addr+=m;
//     for(i=0;i<minor_out;i++){_vol_in+=(m%dim_out[i])*tmp0[pri[i]]; m/=dim_out[i];}
       for(i=0;i<minor_out;i++){k=m/dim_out[i]; _vol_in+=(m-k*dim_out[i])*n2o[i]; m=k;}
       tens_out[_addr]=buf0[_vol_in];
      }
     }
     __syncthreads();
    }else{
//  Algorithm 1 (presumably faster):
//   Create per-block look-up tables:
     m=_vol_minor/_vol_in; //vol_in_c
     for(j=threadIdx.x;j<m;j+=blockDim.x){ //column number (input)
      _addr=0; _s1=j;
//    for(i=minor_in;i<minor;i++){_addr+=(_s1%dim_in[i])*base_in[i]; _s1/=dim_in[i];}
      for(i=minor_in;i<minor;i++){_s2=_s1/dim_in[i]; _addr+=(_s1-_s2*dim_in[i])*base_in[i]; _s1=_s2;}
      ftb[j]=_addr;
     }
     m=_vol_minor/_vol_out; //vol_out_c
     for(j=threadIdx.x;j<m;j+=blockDim.x){ //column number (output)
      _addr=0; _s1=j;
//    for(i=minor_out;i<minor;i++){_addr+=(_s1%dim_out[i])*base_out[i]; _s1/=dim_out[i];}
      for(i=minor_out;i<minor;i++){_s2=_s1/dim_out[i]; _addr+=(_s1-_s2*dim_out[i])*base_out[i]; _s1=_s2;}
      gtb[j]=_addr;
     }
     for(j=threadIdx.x;j<m;j+=blockDim.x){ //column number (output)
      n=0; _s1=j;
//    for(i=minor_out;i<minor;i++){n+=(_s1%dim_out[i])*n2o[i]; _s1/=dim_out[i];}
      for(i=minor_out;i<minor;i++){_s2=_s1/dim_out[i]; n+=(_s1-_s2*dim_out[i])*n2o[i]; _s1=_s2;}
      htb[j]=n;
     }
     for(j=threadIdx.x;j<_vol_out;j+=blockDim.x){
      n=0; _s1=j;
//    for(i=0;i<minor_out;i++){n+=(_s1%dim_out[i])*n2o[i]; _s1/=dim_out[i];}
      for(i=0;i<minor_out;i++){_s2=_s1/dim_out[i]; n+=(_s1-_s2*dim_out[i])*n2o[i]; _s1=_s2;}
      stb[j]=n;
     }
     __syncthreads();
//   Read the minor volume into the buffer from the input tensor block:
     _vol_minor/=_vol_in; //vol_in_c
     _s1=1+(_vol_in-1)/warpSize; //number of warps (lines) which fully cover the input volume
     _s2=blockDim.x/warpSize; //number of whole warps in a thread block (each warp treats one line)
     for(j=l;j<_s1*_vol_minor;j+=_s2){ //j: Line number
      m=j/_s1; n=threadIdx.x+(j-m*_s1-l)*warpSize; //m: Input column number (in_c); n: Offset in the column
      if(n < _vol_in){_addr=_addr_in+ftb[m]+n; buf0[m*_vol_in+n]=tens_in[_addr];}
     }
     __syncthreads();
//   Write the minor volume from the buffer into the output tensor block:
     _vol_minor=(_vol_minor*_vol_in)/_vol_out; //vol_out_c
     _s1=1+(_vol_out-1)/warpSize; //number of warps (lines) which fully cover the output volume
     for(j=l;j<_s1*_vol_minor;j+=_s2){ //j: Line number
      m=j/_s1; n=threadIdx.x+(j-m*_s1-l)*warpSize; //m: Output column number (out_c); n: Offset in the column
      if(n < _vol_out){_addr=_addr_out+gtb[m]+n; _vol_in=htb[m]+stb[n]; tens_out[_addr]=buf0[_vol_in];}
     }
     __syncthreads();
    }
   } //enddo _work_piece: independent work distribution among thread blocks
  }
 }

//Record errors if occured (for each block):
 if(threadIdx.x == 0){if(err_code != 0) i=atomicAdd(&gpu_error_count,1);}
 return;
}
//---------------------------------------------------------------------------------------------------------------------
// TENSOR TRANSPOSE (R4) (scatter version):
__global__ void gpu_tensor_block_copy_scatter_dlf_r4__(int dmo, int drc, int dim_num, int const_args_pos,
                                                       const float* __restrict__ tens_in, float* __restrict__ tens_out)
/**
Scattering version of tensor transpose: tens_out=TRN(tens_in):
INPUT:
 # dmo - dimension extents order (0: normal, as it is in <const_args>; not 0: permuted dimension order will be imposed);
 # drc - index permutation direction (0: normal, as it is in <const_args>; not 0: inversed permutation will be used);
 # dim_num - tensor block rank;
 # const_args_pos - entry in the __constant__ memory bank where tensor block dimension extents (const_args_dims)
                    and index permutation (const_args_prmn) are stored;
 # tens_in[0:] - input tensor;
OUTPUT:
 # tens_out[0:] - output (transposed) tensor;
**/
{
 __shared__ int n2o[MAX_TENSOR_RANK];
 __shared__ size_t vol,base_in[MAX_TENSOR_RANK],base_out[MAX_TENSOR_RANK];
 int i,j,k;
 size_t _vol,_addr_in,_addr_out,_si;

 if(dim_num == 0){
  if(blockIdx.x == 0 && threadIdx.x == 0) tens_out[0]=tens_in[0];
 }else if(dim_num == 1){
  _vol=const_args_dims[const_args_pos][0];
  j=blockIdx.x*blockDim.x+threadIdx.x;
  for(_addr_in=j;_addr_in<_vol;_addr_in+=gridDim.x*blockDim.x){tens_out[_addr_in]=tens_in[_addr_in];}
 }else if(dim_num > 1){
  if(threadIdx.x == 0){
   k=0; for(i=0;i<dim_num;i++){j=const_args_prmn[const_args_pos][i]-1; n2o[j]=i; if(j!=i) k=1;}
   if(k == 0){ //trivial permutation
    n2o[0]=dim_num; //trivial permutation flag
    _vol=1; for(i=0;i<dim_num;i++){_vol*=const_args_dims[const_args_pos][i];}; vol=_vol;
   }else{ //non-trivial permutation
    if(dmo == 0){ //normal dimension order
     _vol=1; for(i=0;i<dim_num;i++){base_in[i]=_vol; _vol*=const_args_dims[const_args_pos][i];}; vol=_vol;
     if(drc == 0){ //normal index permutation
      _vol=1; for(i=0;i<dim_num;i++){k=n2o[i]; base_out[k]=_vol; _vol*=const_args_dims[const_args_pos][k];}
     }else{ //inversed index permutation
      _vol=1; for(i=0;i<dim_num;i++){
       k=const_args_prmn[const_args_pos][i]-1; base_out[k]=_vol; _vol*=const_args_dims[const_args_pos][k];
      }
     }
    }else{ //inversed dimension order
     if(drc == 0){ //normal index permutation
      _vol=1; for(i=0;i<dim_num;i++){
       k=const_args_prmn[const_args_pos][i]-1; base_in[i]=_vol; _vol*=const_args_dims[const_args_pos][k];
      }; vol=_vol;
      _vol=1; for(i=0;i<dim_num;i++){k=n2o[i]; base_out[k]=_vol; _vol*=const_args_dims[const_args_pos][i];}
     }else{ //inversed index permutation
      _vol=1; for(i=0;i<dim_num;i++){
       k=n2o[i]; base_in[i]=_vol; _vol*=const_args_dims[const_args_pos][k];
      }; vol=_vol;
      _vol=1; for(i=0;i<dim_num;i++){
       k=const_args_prmn[const_args_pos][i]-1; base_out[k]=_vol; _vol*=const_args_dims[const_args_pos][i];
      }
     }
    }
   }
  }
#ifdef DEBUG_GPU
//DEBUG RECORD begin:
  if(blockIdx.x == 0 && threadIdx.x == 0){
   j=0; gpu_debug_dump[j++]=dim_num;
   for(i=0;i<dim_num;i++) gpu_debug_dump[j++]=const_args_dims[const_args_pos][i];
   for(i=0;i<dim_num;i++) gpu_debug_dump[j++]=const_args_prmn[const_args_pos][i];
   for(i=0;i<dim_num;i++) gpu_debug_dump[j++]=base_in[i];
   for(i=0;i<dim_num;i++) gpu_debug_dump[j++]=base_out[i];
   gpu_debug_dump[j++]=vol; gpu_debug_dump[j++]=-1;
  }
//DEBUG RECORD end.
#endif
  __syncthreads();
  _vol=vol;
  if(n2o[0] >= dim_num){ //trivial permutation
   k=gridDim.x*blockDim.x; j=blockIdx.x*blockDim.x+threadIdx.x;
   for(_addr_in=j;_addr_in<_vol;_addr_in+=k){tens_out[_addr_in]=tens_in[_addr_in];}
  }else{ //non-trivial permutation
   j=blockIdx.x*blockDim.x+threadIdx.x;
   for(_addr_in=j;_addr_in<_vol;_addr_in+=gridDim.x*blockDim.x){
    _addr_out=0; _si=_addr_in; for(i=dim_num-1;i>=0;i--){_addr_out+=(_si/base_in[i])*base_out[i]; _si%=base_in[i];}
    tens_out[_addr_out]=tens_in[_addr_in];
   }
  }
 }else{ //dim_num < 0
  if(threadIdx.x == 0) i=atomicAdd(&gpu_error_count,1); //record an error (for each thread block)
 }
 return;
}
//-----------------------------------------------------------------------------------------------------------------------
// TENSOR TRANSPOSE (R8) (scatter version):
__global__ void gpu_tensor_block_copy_scatter_dlf_r8__(int dmo, int drc, int dim_num, int const_args_pos,
                                                       const double* __restrict__ tens_in, double* __restrict__ tens_out)
/**
Scattering version of tensor transpose: tens_out=TRN(tens_in):
INPUT:
 # dmo - dimension extents order (0: normal, as it is in <const_args>; not 0: permuted dimension order will be imposed);
 # drc - index permutation direction (0: normal, as it is in <const_args>; not 0: inversed permutation will be used);
 # dim_num - tensor block rank;
 # const_args_pos - entry in the __constant__ memory bank where tensor block dimension extents (const_args_dims)
                    and index permutation (const_args_prmn) are stored;
 # tens_in[0:] - input tensor;
OUTPUT:
 # tens_out[0:] - output (transposed) tensor;
**/
{
 __shared__ int n2o[MAX_TENSOR_RANK];
 __shared__ size_t vol,base_in[MAX_TENSOR_RANK],base_out[MAX_TENSOR_RANK];
 int i,j,k;
 size_t _vol,_addr_in,_addr_out,_si;

 if(dim_num == 0){
  if(blockIdx.x == 0 && threadIdx.x == 0) tens_out[0]=tens_in[0];
 }else if(dim_num == 1){
  _vol=const_args_dims[const_args_pos][0];
  j=blockIdx.x*blockDim.x+threadIdx.x;
  for(_addr_in=j;_addr_in<_vol;_addr_in+=gridDim.x*blockDim.x){tens_out[_addr_in]=tens_in[_addr_in];}
 }else if(dim_num > 1){
  if(threadIdx.x == 0){
   k=0; for(i=0;i<dim_num;i++){j=const_args_prmn[const_args_pos][i]-1; n2o[j]=i; if(j!=i) k=1;}
   if(k == 0){ //trivial permutation
    n2o[0]=dim_num; //trivial permutation flag
    _vol=1; for(i=0;i<dim_num;i++){_vol*=const_args_dims[const_args_pos][i];}; vol=_vol;
   }else{ //non-trivial permutation
    if(dmo == 0){ //normal dimension order
     _vol=1; for(i=0;i<dim_num;i++){base_in[i]=_vol; _vol*=const_args_dims[const_args_pos][i];}; vol=_vol;
     if(drc == 0){ //normal index permutation
      _vol=1; for(i=0;i<dim_num;i++){k=n2o[i]; base_out[k]=_vol; _vol*=const_args_dims[const_args_pos][k];}
     }else{ //inversed index permutation
      _vol=1; for(i=0;i<dim_num;i++){
       k=const_args_prmn[const_args_pos][i]-1; base_out[k]=_vol; _vol*=const_args_dims[const_args_pos][k];
      }
     }
    }else{ //inversed dimension order
     if(drc == 0){ //normal index permutation
      _vol=1; for(i=0;i<dim_num;i++){
       k=const_args_prmn[const_args_pos][i]-1; base_in[i]=_vol; _vol*=const_args_dims[const_args_pos][k];
      }; vol=_vol;
      _vol=1; for(i=0;i<dim_num;i++){k=n2o[i]; base_out[k]=_vol; _vol*=const_args_dims[const_args_pos][i];}
     }else{ //inversed index permutation
      _vol=1; for(i=0;i<dim_num;i++){
       k=n2o[i]; base_in[i]=_vol; _vol*=const_args_dims[const_args_pos][k];
      }; vol=_vol;
      _vol=1; for(i=0;i<dim_num;i++){
       k=const_args_prmn[const_args_pos][i]-1; base_out[k]=_vol; _vol*=const_args_dims[const_args_pos][i];
      }
     }
    }
   }
  }
#ifdef DEBUG_GPU
//DEBUG RECORD begin:
  if(blockIdx.x == 0 && threadIdx.x == 0){
   j=0; gpu_debug_dump[j++]=dim_num;
   for(i=0;i<dim_num;i++) gpu_debug_dump[j++]=const_args_dims[const_args_pos][i];
   for(i=0;i<dim_num;i++) gpu_debug_dump[j++]=const_args_prmn[const_args_pos][i];
   for(i=0;i<dim_num;i++) gpu_debug_dump[j++]=base_in[i];
   for(i=0;i<dim_num;i++) gpu_debug_dump[j++]=base_out[i];
   gpu_debug_dump[j++]=vol; gpu_debug_dump[j++]=-1;
  }
//DEBUG RECORD end.
#endif
  __syncthreads();
  _vol=vol;
  if(n2o[0] >= dim_num){ //trivial permutation
   k=gridDim.x*blockDim.x; j=blockIdx.x*blockDim.x+threadIdx.x;
   for(_addr_in=j;_addr_in<_vol;_addr_in+=k){tens_out[_addr_in]=tens_in[_addr_in];}
  }else{ //non-trivial permutation
   j=blockIdx.x*blockDim.x+threadIdx.x;
   for(_addr_in=j;_addr_in<_vol;_addr_in+=gridDim.x*blockDim.x){
    _addr_out=0; _si=_addr_in; for(i=dim_num-1;i>=0;i--){_addr_out+=(_si/base_in[i])*base_out[i]; _si%=base_in[i];}
    tens_out[_addr_out]=tens_in[_addr_in];
   }
  }
 }else{ //dim_num < 0
  if(threadIdx.x == 0) i=atomicAdd(&gpu_error_count,1); //record an error (for each thread block)
 }
 return;
}
//----------------------------------------------------------------------------------------------------------
// MATRIX MULTIPLICATION (R4) (shared-memory version):
__global__ void gpu_matrix_multiply_tn_r4__(size_t ll, size_t lr, size_t lc, const float* __restrict__ arg1,
                                            const float* __restrict__ arg2, float* __restrict__ arg0)
/** arg0(0:ll-1,0:lr-1)+=arg1(0:lc-1,0:ll-1)*arg2(0:lc-1,0:lr-1)
NOTES:
 # Thread block dimensions (.x and .y) must be equal to MAT_MULT_TILE_DIM(X,Y), respectively.
**/
{
 __shared__ float buf1[MAT_MULT_TILE_DIMY+1][MAT_MULT_TILE_DIMY+1],buf2[MAT_MULT_TILE_DIMY+1][MAT_MULT_TILE_DIMY+1];
 size_t k,_col,_row,_col_base,_row_base;
 int i,j,l,m;
 float _val;

 if(lc > 0 && ll > 0 && lr > 0 && blockDim.x == MAT_MULT_TILE_DIMY && blockDim.y == MAT_MULT_TILE_DIMY){
  _val=0.0f; j=threadIdx.y; i=threadIdx.x;
  _col_base=blockIdx.y*MAT_MULT_TILE_DIMY;
  while(_col_base < lr){
   _row_base=blockIdx.x*MAT_MULT_TILE_DIMY;
   while(_row_base < ll){
    for(k=0;k<lc;k+=MAT_MULT_TILE_DIMY){
     _col=_col_base+j; _row=_row_base+j;
// Load two blocks into shared memory:
     if(k+MAT_MULT_TILE_DIMY > lc){
      m=lc-k;
      if(i < m){ //(k+i)<lc
       if(_row < ll){buf1[j][i]=arg1[_row*lc+(k+i)];} // Load a block of the 1st argument into the shared memory
       if(_col < lr){buf2[j][i]=arg2[_col*lc+(k+i)];} // Load a block of the 2nd argument into the shared memory
      }
     }else{
      m=MAT_MULT_TILE_DIMY;
      if(_row < ll){buf1[j][i]=arg1[_row*lc+(k+i)];} // Load a block of the 1st argument into the shared memory
      if(_col < lr){buf2[j][i]=arg2[_col*lc+(k+i)];} // Load a block of the 2nd argument into the shared memory
     }
     __syncthreads();
// Multiply the two blocks:
     _row=_row_base+i;
     if(_col < lr){
      if(_row < ll){
       _col=_col*ll+_row;
       for(l=0;l<m;l++){_val+=buf1[i][l]*buf2[j][l];}
       arg0[_col]+=_val; _val=0.0f;
      }
     }
     __syncthreads();
    }
    _row_base+=gridDim.x*MAT_MULT_TILE_DIMY;
   }
   _col_base+=gridDim.y*MAT_MULT_TILE_DIMY;
  }
 }else{
  if(threadIdx.x == 0 && threadIdx.y == 0) i=atomicAdd(&gpu_error_count,1); //record an error (for each thread block)
 }
 return;
}
//-----------------------------------------------------------------------------------------------------------
// MATRIX MULTIPLICATION (R8) (shared-memory version):
__global__ void gpu_matrix_multiply_tn_r8__(size_t ll, size_t lr, size_t lc, const double* __restrict__ arg1,
                                            const double* __restrict__ arg2, double* __restrict__ arg0)
/** arg0(0:ll-1,0:lr-1)+=arg1(0:lc-1,0:ll-1)*arg2(0:lc-1,0:lr-1)
NOTES:
 # Thread block dimensions (.x and .y) must be equal to MAT_MULT_TILE_DIM(X,Y), respectively.
**/
{
 __shared__ double buf1[MAT_MULT_TILE_DIMX+1][MAT_MULT_TILE_DIMX+1],buf2[MAT_MULT_TILE_DIMY+1][MAT_MULT_TILE_DIMX+1];
 size_t k,_col,_row,_col_base,_row_base;
 int i,j,l,m;
 double _val;

 if(lc > 0 && ll > 0 && lr > 0 && blockDim.x == MAT_MULT_TILE_DIMX && blockDim.y == MAT_MULT_TILE_DIMY){
  _val=0.0; j=threadIdx.y; i=threadIdx.x;
  _col_base=blockIdx.y*MAT_MULT_TILE_DIMY;
  while(_col_base < lr){
   _row_base=blockIdx.x*MAT_MULT_TILE_DIMX;
   while(_row_base < ll){
    for(k=0;k<lc;k+=MAT_MULT_TILE_DIMX){
     _col=_col_base+j; _row=_row_base+j;
// Load two blocks into shared memory:
     if(k+MAT_MULT_TILE_DIMX > lc){m=lc-k;}else{m=MAT_MULT_TILE_DIMX;}
     if(i < m){ //(k+i)<lc
      for(l=0;l<MAT_MULT_TILE_DIMX;l+=MAT_MULT_TILE_DIMY){
       if(_row < ll){buf1[l+j][i]=arg1[_row*lc+(k+i)];} // Load a block of the 1st argument into the shared memory
       _row+=MAT_MULT_TILE_DIMY;
      }
      if(_col < lr){buf2[j][i]=arg2[_col*lc+(k+i)];} // Load a block of the 2nd argument into the shared memory
     }
     __syncthreads();
// Multiply the two blocks:
     _row=_row_base+i;
     if(_col < lr){
      if(_row < ll){
       _col=_col*ll+_row;
       for(l=0;l<m;l++){_val+=buf1[i][l]*buf2[j][l];}
       arg0[_col]+=_val; _val=0.0;
      }
     }
     __syncthreads();
    }
    _row_base+=gridDim.x*MAT_MULT_TILE_DIMX;
   }
   _col_base+=gridDim.y*MAT_MULT_TILE_DIMY;
  }
 }else{
  if(threadIdx.x == 0 && threadIdx.y == 0) i=atomicAdd(&gpu_error_count,1); //record an error (for each thread block)
 }
 return;
}
//-D NO_GPU
#endif
