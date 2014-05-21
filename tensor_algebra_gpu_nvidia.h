/** Parameters and derived types used in tensor_algebra_gpu_nvidia.cu and related modules. **/

/** NOTES:
 # GPU_ID is the unique CUDA GPU ID given to a specific Nvidia GPU present on the Host node:
    0<=GPU_ID<MAX_GPUS_PER_NODE; GPU_ID=-1 will refer to the (multi-)CPU Host.
    All MPI processes running on the same node will have a flat numeration of all present GPUs,
    each process having its own subrange of present GPUs.
 # MIC_ID is defined completely analogously to GPU_ID.
 # In general, ACC_ID (Accelerator ID) is its ID within its class.
 # DEVICE_ID can refer either to the (multi-)CPU Host (==0) OR
    to a specific Nvidia GPU (!=0): gpu_id=abs(DEVICE_ID)-1 OR
    to a specific Intel MIC (!=0): mic_id=abs(DEVICE_ID)-1-MAX_GPUS_PER_NODE OR
    etc. Device numeration:
    {0}(Host);{1..MAX_GPUS_PER_NODE}(Nvidia GPUs);{MAX_GPUS_PER_NODE+1:MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE}(Intel MICs);etc.
    DEVICE_ID is used in tensBlck_t: If tensor elements are already on the Device it is positive, otherwise negative.
 # MAX_SCR_ENTRY_COUNT regulates the maximal amount of additional device argument-buffer entries
    allocated per tensor operation (it is no more than 3 for tensor contractions).
 # MAX_GPU_ARGS regulates the maximal amount of occupied device argument-buffer entries.
 # CUDA_TASK is considered completed successfully if the value of the .task_error field equals zero.
    Negative .task_error means that either the CUDA task is empty or it is in progress.
    In the former case, .gpu_id=-1 and .task_stream is undefined.
    Positive .task_error means that an error occured during the task scheduling process.
**/

//GLOBAL PARAMETERS:
#define MAX_TENSOR_RANK 32         //max allowed tensor rank
#define MAX_GPUS_PER_NODE 16       //max number of Nvidia GPUs on a node
#define MAX_GPU_ARGS 128           //max total number of tensor arguments simultaneously residing on a GPU device
#define MAX_SCR_ENTRY_COUNT 3      //max allowed number of additional GPU argument entries allocated per tensor operation
#define MAX_MICS_PER_NODE 8        //max number of Intel MICs on a node
#define MAX_AMDS_PER_NODE 8        //max number of AMD GPUs on a node

//KERNEL PARAMETERS:
#define GPU_CACHE_LINE_LEN 128     //cache line length in bytes
#define MAX_CUDA_BLOCKS 1024       //max number of CUDA thread blocks per kernel
#define TENS_TRANSP_BUF_SIZE 1024  //buffer size (elements) for <gpu_tensor_block_copy_dlf_XX__>
#define TENS_TRANSP_TAB_SIZE 69    //look up table size (integers) for <gpu_tensor_block_copy_dlf_XX__>
#define MAT_MULT_TILE_DIM 16       //tile dimension size for <gpu_matrix_multiply_tn_XX__>
#define THRDS_ARRAY_PRODUCT 256    //threads per block for <gpu_array_product_XX__>
#define THRDS_ARRAY_NORM2 256      //threads per block for <gpu_array_2norm2_XX__>
#define THRDS_ARRAY_INIT 256       //threads per block for <gpu_array_init_XX__>
#define THRDS_ARRAY_SCALE 256      //threads per block for <gpu_array_scale_XX__> and <gpu_array_dot_product_XX__>
#define THRDS_ARRAY_ADD 256        //threads per block for <gpu_array_add_XX__>
#define THRDS_TENSOR_COPY 192      //threads per block for <gpu_tensor_block_copy_dlf_XX__>
#define THRDS_TENSOR_COPY_SCAT 256 //threads per block for <gpu_tensor_block_copy_scatter_dlf_XX__>

//DATA KINDS:
#define R4 4     //float data kind (keep consistent with c_process.f90:tens_blck_pack/unpack)
#define R8 8     //double data kind (keep consistent with c_process.f90:tens_blck_pack/unpack)
#define C8 16    //double complex data kind (keep consistent with c_process.f90:tens_blck_pack/unpack)

//DEVICE KINDS:
#define DEV_HOST 0
#define DEV_NVIDIA_GPU 1
#define DEV_INTEL_MIC 2
#define DEV_AMD_GPU 3
#define DEV_MAX 1+MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE+MAX_AMDS_PER_NODE

//CUDA TASK STATUS (must be consistent with qforce.f90!):
#define CUDA_TASK_ERROR -1
#define CUDA_TASK_EMPTY 0
#define CUDA_TASK_SCHEDULED 1
#define CUDA_TASK_STARTED 2
#define CUDA_TASK_INPUT_THERE 3
#define CUDA_TASK_OUTPUT_THERE 4
#define CUDA_TASK_COMPLETED 5

//ALIASES:
#define NO_COPY_BACK 0
#define COPY_BACK 1
#define NOT_REALLY 0
#define GPU_MINE 1
#define GPU_MINE_CUBLAS 2

//MACRO FUNCTIONS:
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))

//DERIVED TYPES:
// Tensor block:
typedef struct{
 int device_id; //GPU on which the block already resides (+) or will reside (-): GPU# = |device_id|-1 (device_id=0 means Host)
 int data_kind; //Element byte size: float (4), double (8), or double complex (16)
 int rank; //tensor block rank (>=0)
 int *dims_h; //tensor block dimension extents (0th index is the most minor one): HOST memory
 int *divs_h; //tensor block dimension dividers: HOST memory
 int *grps_h; //tensor block dimension groups: HOST memory
 int *prmn_h; //tensor block dimension permutation (not to be set by user!): HOST memory
 void *elems_h; //tensor block elements (dlf): HOST memory (only one element for scalars)
 void *elems_d; //tensor block elements (dlf): GPU global memory (GPU#=|device_id|-1) (only one element for scalars)
 int buf_entry_gpu; //GPU argument buffer entry pointed to by *elems_d: GPU global memory
 int const_args_entry; //entry number in const_args[]: GPU constant memory (dims[] and prmn[] arrays)
} tensBlck_t;
#ifndef NO_GPU
// CUDA task information (returned by non-blocking CUDA calling functions):
typedef struct{
 int task_error;                     //error code (<0: Task is either empty or in progress; 0: Success; >0: Launch error code)
 int gpu_id;                         //Nvidia GPU ID on which the task was scheduled (-1 means CPU Host)
 cudaStream_t task_stream;           //CUDA stream the task went into
 cudaEvent_t task_start;             //CUDA event recorded at the beginning of the task
 cudaEvent_t task_comput;            //CUDA event recorded when the computing kernels start (all input data is on Device)
 cudaEvent_t task_output;            //CUDA event recorded when the computing kernels finish (before output to the Host)
 cudaEvent_t task_finish;            //CUDA event recorded at the end of the task (full completion)
 int scr_entry_count;                //number of additional GPU argument-buffer entries allocated by the task
 int scr_entry[MAX_SCR_ENTRY_COUNT]; //additional GPU argument-buffer entries allocated by the task
} cudaTask_t;
#endif
