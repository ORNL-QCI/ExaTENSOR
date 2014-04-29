!This is a Fortran alias of tensor_algebra_gpu_nvidia.h (KEEP CONSISTENT!):

!GLOBAL (keep consistent with tensor_algebra_gpu_nvidia.h):
        integer(C_INT), parameter:: max_tensor_rank=32   !max allowed tensor rank (max number of indices in a tensor)
        integer(C_INT), parameter:: MAX_GPUS_PER_NODE=16 !max number of Nvidia GPUs on a node
        integer(C_INT), parameter:: MAX_MICS_PER_NODE=4  !max number of Intel MICs on a node

!DATA KINDS (keep consistent with tensor_algebra_gpu_nvidia.h):
        integer(C_INT), parameter:: R4=4  !float data kind
        integer(C_INT), parameter:: R8=8  !double data kind
        integer(C_INT), parameter:: C8=16 !double complex data kind

!Tensor block storage layout:
        integer(C_INT), parameter:: not_allocated=0   !tensor block has not been initialized
        integer(C_INT), parameter:: scalar_tensor=1   !scalar (rank-0 tensor)
        integer(C_INT), parameter:: dimension_led=2   !dense tensor block (column-major storage by default): no symmetry re
        integer(C_INT), parameter:: bricked_dense=3   !dense tensor block (bricked storage): no symmetry restrictions, %dim
        integer(C_INT), parameter:: bricked_ordered=4 !symmetrically packed tensor block (bricked storage): symmetry restri
        integer(C_INT), parameter:: sparse_list=5     !sparse tensor block: symmetry restrictions do not apply!
        integer(C_INT), parameter:: compressed=6      !compressed tensor block: symmetry restrictions do not apply!
        logical, parameter:: fortran_like=.true.
        logical, parameter:: c_like=.false.

!DEVICE KINDS (keep consistent with tensor_algebra_gpu_nvidia.h):
        integer(C_INT), parameter:: DEV_HOST=0
        integer(C_INT), parameter:: DEV_NVIDIA_GPU=1
        integer(C_INT), parameter:: DEV_INTEL_MIC=2
        integer(C_INT), parameter:: DEV_MAX=1+MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE

!CUDA task status (keep consistent with tensor_algebra_gpu_nvidia.h):
        integer(C_INT), parameter:: cuda_task_error=-1
        integer(C_INT), parameter:: cuda_task_empty=0
        integer(C_INT), parameter:: cuda_task_scheduled=1
        integer(C_INT), parameter:: cuda_task_started=2
        integer(C_INT), parameter:: cuda_task_input_there=3
        integer(C_INT), parameter:: cuda_task_output_there=4
        integer(C_INT), parameter:: cuda_task_completed=5

!ALIASES (keep consistent with tensor_algebra_gpu_nvidia.h):
        integer(C_INT), parameter:: NOT_REALLY=0
        integer(C_INT), parameter:: NO_COPY_BACK=0
        integer(C_INT), parameter:: COPY_BACK=1
        integer(C_INT), parameter:: EVENTS_OFF=0
        integer(C_INT), parameter:: EVENTS_ON=1
        integer(C_INT), parameter:: BLAS_ON=0
        integer(C_INT), parameter:: BLAS_OFF=1
        integer(C_INT), parameter:: EFF_TRN_OFF=0
        integer(C_INT), parameter:: EFF_TRN_ON=1
