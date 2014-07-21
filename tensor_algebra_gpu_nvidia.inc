!This is a Fortran alias of tensor_algebra_gpu_nvidia.h (KEEP CONSISTENT!):

!GLOBAL (keep consistent with tensor_algebra_gpu_nvidia.h):
        integer(C_INT), parameter:: max_tensor_rank=32    !max allowed tensor rank (max number of indices in a tensor)
        integer(C_INT), parameter:: max_tensor_operands=4 !max number of tensor operands in a tensor operation
        integer(C_INT), parameter:: MAX_GPUS_PER_NODE=16  !max number of Nvidia GPUs on a node
        integer(C_INT), parameter:: MAX_MICS_PER_NODE=8   !max number of Intel MICs on a node
        integer(C_INT), parameter:: MAX_AMDS_PER_NODE=8   !max number of AMD GPUs on a node

!DATA KINDS (keep consistent with tensor_algebra_gpu_nvidia.h):
        integer(C_INT), parameter:: R4=4  !float data kind
        integer(C_INT), parameter:: R8=8  !double data kind
        integer(C_INT), parameter:: C8=16 !double complex data kind

!DEVICE KINDS (keep consistent with tensor_algebra_gpu_nvidia.h):
        integer(C_INT), parameter:: DEV_HOST=0
        integer(C_INT), parameter:: DEV_NVIDIA_GPU=1
        integer(C_INT), parameter:: DEV_INTEL_MIC=2
        integer(C_INT), parameter:: DEV_AMD_GPU=3
        integer(C_INT), parameter:: DEV_MAX=1+MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE+MAX_AMDS_PER_NODE

!CUDA task status (keep consistent with tensor_algebra_gpu_nvidia.h):
        integer(C_INT), parameter:: cuda_task_error=-1
        integer(C_INT), parameter:: cuda_task_empty=0
        integer(C_INT), parameter:: cuda_task_scheduled=1
        integer(C_INT), parameter:: cuda_task_started=2
        integer(C_INT), parameter:: cuda_task_input_there=3
        integer(C_INT), parameter:: cuda_task_output_there=4
        integer(C_INT), parameter:: cuda_task_completed=5

!ALIASES (keep consistent with tensor_algebra_gpu_nvidia.h):
        integer(C_INT), parameter:: NOT_REALLY=0                  !"NO" answer
        integer(C_INT), parameter:: NO_COPY_BACK=0                !keeps the tensor-result on a GPU without updating HAB
        integer(C_INT), parameter:: COPY_BACK=1                   !tensor-result will be copied back from a GPU to HAB
        integer(C_INT), parameter:: EVENTS_OFF=0                  !disables CUDA event recording
        integer(C_INT), parameter:: EVENTS_ON=1                   !enables CUDA event recording
        integer(C_INT), parameter:: BLAS_ON=0                     !enables BLAS
        integer(C_INT), parameter:: BLAS_OFF=1                    !disables BLAS
        integer(C_INT), parameter:: EFF_TRN_OFF=0                 !disables efficient tensor transpose algorithm
        integer(C_INT), parameter:: EFF_TRN_ON=1                  !enables efficient tensor transpose algorithm
        integer(C_INT), parameter:: KEEP_NO_DATA=0                !ETI task cleanup flag: delete all data for the tensor argument
        integer(C_INT), parameter:: KEEP_TBB_DATA=1               !ETI task cleanup flag: keep the TBB Fortran pointer (if any)
        integer(C_INT), parameter:: KEEP_HAB_DATA=2               !ETI task cleanup flag: keep the HAB data (if any)
        integer(C_INT), parameter:: KEEP_GPU_DATA=4               !ETI task cleanup flag: keep the GPU data (if any)
        integer(C_INT), parameter:: KEEP_ALL_DATA=7               !ETI task cleanup flag: keep everything
        integer(C_INT), parameter:: ERR_ETIQ_INVALID_SUPERPACK=1  !ETIQ error: an invalid ETI superpacket received from LR
        integer(C_INT), parameter:: ERR_ETIQ_IS_FULL=2            !ETIQ error: ETIQ is full
        integer(C_INT), parameter:: ERR_ETIQ_UNPACK_FAILED=3      !ETIQ error: ETI unpacking failed
        integer(C_INT), parameter:: ERR_ETIQ_ENTRY_UNAVAIL=4      !ETIQ error: unable to obtain a free ETIQ entry
        integer(C_INT), parameter:: ERR_ETIQ_CHANNEL_FAILED=5     !ETIQ error: unable to register an ETIQ entry in an ETIQ channel

!Tensor block storage layout:
        integer(C_INT), parameter:: not_allocated=0   !tensor block has not been initialized
        integer(C_INT), parameter:: scalar_tensor=1   !scalar (rank-0 tensor)
        integer(C_INT), parameter:: dimension_led=2   !dense tensor block (column-major storage by default): no symmetry restrictions
        integer(C_INT), parameter:: bricked_dense=3   !dense tensor block (bricked storage): no symmetry restrictions
        integer(C_INT), parameter:: bricked_ordered=4 !symmetrically packed tensor block (bricked storage): symmetry restrictions apply
        integer(C_INT), parameter:: sparse_list=5     !sparse tensor block: symmetry restrictions do not apply!
        integer(C_INT), parameter:: compressed=6      !compressed tensor block: symmetry restrictions do not apply!
        logical, parameter:: fortran_like=.true.
        logical, parameter:: c_like=.false.
