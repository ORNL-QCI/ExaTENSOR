!This is a Fortran alias of tensor_algebra_gpu_nvidia.h:

!GLOBAL (keep consistent with tensor_algebra_gpu_nvidia.h):
        integer(C_INT), parameter:: max_tensor_rank=32    !max allowed tensor rank (max number of indices in a tensor)
        integer(C_INT), parameter:: max_tensor_operands=4 !max number of tensor operands in a tensor operation
!DIR$ ATTRIBUTES OFFLOAD:mic:: max_tensor_rank,max_tensor_operands
!DIR$ ATTRIBUTES ALIGN:128:: max_tensor_rank,max_tensor_operands

!DEVICE KINDS (keep consistent with tensor_algebra_gpu_nvidia.h):
        integer(C_INT), parameter:: MAX_GPUS_PER_NODE=16  !max number of Nvidia GPUs on a node
        integer(C_INT), parameter:: MAX_MICS_PER_NODE=8   !max number of Intel MICs on a node
        integer(C_INT), parameter:: MAX_AMDS_PER_NODE=8   !max number of AMD GPUs on a node
        integer(C_INT), parameter:: DEV_HOST=0
        integer(C_INT), parameter:: DEV_NVIDIA_GPU=1
        integer(C_INT), parameter:: DEV_INTEL_MIC=2
        integer(C_INT), parameter:: DEV_AMD_GPU=3
        integer(C_INT), parameter:: DEV_MAX=1+MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE+MAX_AMDS_PER_NODE

!DATA KINDS (keep consistent with tensor_algebra_gpu_nvidia.h):
        integer(C_INT), parameter:: R4=4  !float data kind
        integer(C_INT), parameter:: R8=8  !double data kind
        integer(C_INT), parameter:: C8=16 !double complex data kind
!DIR$ ATTRIBUTES OFFLOAD:mic:: R4,R8,C8
!DIR$ ATTRIBUTES ALIGN:128:: R4,R8,C8

!CUDA TASK STATUS (keep consistent with tensor_algebra_gpu_nvidia.h):
        integer(C_INT), parameter:: cuda_task_error=-1
        integer(C_INT), parameter:: cuda_task_empty=0
        integer(C_INT), parameter:: cuda_task_scheduled=1
        integer(C_INT), parameter:: cuda_task_started=2
        integer(C_INT), parameter:: cuda_task_input_there=3
        integer(C_INT), parameter:: cuda_task_output_there=4
        integer(C_INT), parameter:: cuda_task_completed=5

!ALIASES (keep consistent with tensor_algebra_gpu_nvidia.h):
        integer(C_INT), parameter:: NOPE=0                     !"NO" answer
        integer(C_INT), parameter:: YEP=1                      !"YES" answer
        integer(C_INT), parameter:: NO_COPY_BACK=0             !keeps the tensor-result on GPU without updating HAB
        integer(C_INT), parameter:: COPY_BACK=1                !tensor-result will be copied back from GPU to HAB (default)
        integer(C_INT), parameter:: EVENTS_OFF=0               !disables CUDA event recording
        integer(C_INT), parameter:: EVENTS_ON=1                !enables CUDA event recording
        integer(C_INT), parameter:: BLAS_ON=0                  !enables BLAS
        integer(C_INT), parameter:: BLAS_OFF=1                 !disables BLAS
        integer(C_INT), parameter:: EFF_TRN_OFF=0              !disables efficient tensor transpose algorithm
        integer(C_INT), parameter:: EFF_TRN_ON=1               !enables efficient tensor transpose algorithm
!DIR$ ATTRIBUTES OFFLOAD:mic:: NOPE,YEP,NO_COPY_BACK,COPY_BACK,BLAS_ON,BLAS_OFF,EFF_TRN_OFF,EFF_TRN_ON
!DIR$ ATTRIBUTES ALIGN:128:: NOPE,YEP,NO_COPY_BACK,COPY_BACK,BLAS_ON,BLAS_OFF,EFF_TRN_OFF,EFF_TRN_ON

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
!DIR$ ATTRIBUTES OFFLOAD:mic:: not_allocated,scalar_tensor,dimension_led,bricked_dense,bricked_ordered,sparse_list,compressed
!DIR$ ATTRIBUTES OFFLOAD:mic:: fortran_like,c_like
!DIR$ ATTRIBUTES ALIGN:128:: not_allocated,scalar_tensor,dimension_led,bricked_dense,bricked_ordered,sparse_list,compressed
!DIR$ ATTRIBUTES ALIGN:128:: fortran_like,c_like
