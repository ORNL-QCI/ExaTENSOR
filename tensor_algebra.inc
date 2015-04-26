!ExaTensor Public Parameters:

!BASIC TYPE SIZES:
        integer(C_INT), parameter:: INTD=4            !default integer size
        integer(C_INT), parameter:: INTL=8            !long integer size
        integer(C_INT), parameter:: REALD=8           !default real size
        integer(C_INT), parameter:: INT_MPI=4         !default MPI integer size
        integer(C_INT), parameter:: INT_ADDR=C_SIZE_T !default addressing integer
!DIR$ ATTRIBUTES OFFLOAD:mic:: INTD,INTL,REALD,INT_MPI,INT_ADDR
!DIR$ ATTRIBUTES ALIGN:128:: INTD,INTL,REALD,INT_MPI,INT_ADDR

!TENSOR ALGEBRA LIMITS (keep consistent with tensor_algebra.h):
        integer(C_INT), parameter:: MAX_TENSOR_RANK=32    !max allowed tensor rank (max number of indices in a tensor)
        integer(C_INT), parameter:: MAX_TENSOR_OPERANDS=4 !max number of tensor operands in a tensor operation
!DIR$ ATTRIBUTES OFFLOAD:mic:: MAX_TENSOR_RANK,MAX_TENSOR_OPERANDS
!DIR$ ATTRIBUTES ALIGN:128:: MAX_TENSOR_RANK,MAX_TENSOR_OPERANDS

!DEVICE KINDS (keep consistent with tensor_algebra.h):
        integer(C_INT), parameter:: MAX_GPUS_PER_NODE=8   !max number of NVidia GPUs on a node
        integer(C_INT), parameter:: MAX_MICS_PER_NODE=8   !max number of Intel MICs on a node
        integer(C_INT), parameter:: MAX_AMDS_PER_NODE=8   !max number of AMD GPUs on a node
        integer(C_INT), parameter:: DEV_HOST=0
        integer(C_INT), parameter:: DEV_NVIDIA_GPU=1
        integer(C_INT), parameter:: DEV_INTEL_MIC=2
        integer(C_INT), parameter:: DEV_AMD_GPU=3
        integer(C_INT), parameter:: DEV_MAX=1+MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE+MAX_AMDS_PER_NODE

!TENSOR DATA KINDS (keep consistent with tensor_algebra.h):
        integer(C_INT), parameter:: R4=4  !float data kind
        integer(C_INT), parameter:: R8=8  !double data kind
        integer(C_INT), parameter:: C8=16 !double complex data kind
!DIR$ ATTRIBUTES OFFLOAD:mic:: R4,R8,C8
!DIR$ ATTRIBUTES ALIGN:128:: R4,R8,C8

!ALIASES (keep consistent with tensor_algebra.h):
        integer(C_INT), parameter:: NOPE=0                     !"NO" answer
        integer(C_INT), parameter:: YEP=1                      !"YES" answer
        integer(C_INT), parameter:: NO_COPY_BACK=0             !keeps the tensor-result on Accelerator without updating Host
        integer(C_INT), parameter:: COPY_BACK=1                !tensor-result will be copied back from Accelerator to Host (default)
        integer(C_INT), parameter:: BLAS_ON=0                  !enables BLAS
        integer(C_INT), parameter:: BLAS_OFF=1                 !disables BLAS
        integer(C_INT), parameter:: EFF_TRN_OFF=0              !disables efficient tensor transpose algorithm
        integer(C_INT), parameter:: EFF_TRN_ON=1               !enables efficient tensor transpose algorithm
        integer(C_INT), parameter:: EVENTS_OFF=0               !disables CUDA event recording
        integer(C_INT), parameter:: EVENTS_ON=1                !enables CUDA event recording (default)
!DIR$ ATTRIBUTES OFFLOAD:mic:: NOPE,YEP,NO_COPY_BACK,COPY_BACK,BLAS_ON,BLAS_OFF,EFF_TRN_OFF,EFF_TRN_ON
!DIR$ ATTRIBUTES ALIGN:128:: NOPE,YEP,NO_COPY_BACK,COPY_BACK,BLAS_ON,BLAS_OFF,EFF_TRN_OFF,EFF_TRN_ON

!CUDA TASK STATUS (keep consistent with tensor_algebra.h):
        integer(C_INT), parameter:: CUDA_TASK_ERROR=-1
        integer(C_INT), parameter:: CUDA_TASK_EMPTY=0
        integer(C_INT), parameter:: CUDA_TASK_SCHEDULED=1
        integer(C_INT), parameter:: CUDA_TASK_STARTED=2
        integer(C_INT), parameter:: CUDA_TASK_INPUT_THERE=3
        integer(C_INT), parameter:: CUDA_TASK_OUTPUT_THERE=4
        integer(C_INT), parameter:: CUDA_TASK_COMPLETED=5

!TENSOR BLOCK STORATE LAYOUT:
        integer(C_INT), parameter:: NOT_ALLOCATED=0   !tensor block has not been allocated/initialized
        integer(C_INT), parameter:: SCALAR_TENSOR=1   !scalar (rank-0 tensor)
        integer(C_INT), parameter:: DIMENSION_LED=2   !dense tensor block (column-major storage by default): no symmetry restrictions
        integer(C_INT), parameter:: BRICKED_DENSE=3   !dense tensor block (bricked storage): no symmetry restrictions
        integer(C_INT), parameter:: BRICKED_ORDERED=4 !symmetrically packed tensor block (bricked storage): symmetry restrictions apply
        integer(C_INT), parameter:: SPARSE_LIST=5     !sparse tensor block: symmetry restrictions do not apply!
        integer(C_INT), parameter:: COMPRESSED=6      !compressed tensor block: symmetry restrictions do not apply!
        logical, parameter:: FORTRAN_LIKE=.true.
        logical, parameter:: C_LIKE=.false.
!DIR$ ATTRIBUTES OFFLOAD:mic:: NOT_ALLOCATED,SCALAR_TENSOR,DIMENSION_LED,BRICKED_DENSE,BRICKED_ORDERED,SPARSE_LIST,COMPRESSED
!DIR$ ATTRIBUTES OFFLOAD:mic:: FORTRAN_LIKE,C_LIKE
!DIR$ ATTRIBUTES ALIGN:128:: NOT_ALLOCATED,SCALAR_TENSOR,DIMENSION_LED,BRICKED_DENSE,BRICKED_ORDERED,SPARSE_LIST,COMPRESSED
!DIR$ ATTRIBUTES ALIGN:128:: FORTRAN_LIKE,C_LIKE
