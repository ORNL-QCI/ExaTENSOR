!ExaTensor: Basic parameters and definitions:
        module tensor_algebra
        use, intrinsic:: ISO_C_BINDING
        implicit none
!BASIC TYPE SIZES:
        integer(C_INT), parameter, public:: INTD=4  !default integer size
        integer(C_INT), parameter, public:: INTL=8  !long integer size
        integer(C_INT), parameter, public:: REALD=8 !default real size
#ifndef NO_PHI
!DIR$ ATTRIBUTES OFFLOAD:mic:: INTD,INTL,REALD
!DIR$ ATTRIBUTES ALIGN:128:: INTD,INTL,REALD
#endif

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
        integer(C_INT), parameter, public:: MAX_AMDS_PER_NODE=8   !max number of AMD GPUs on a node
        integer(C_INT), parameter, public:: DEV_HOST=0
        integer(C_INT), parameter, public:: DEV_NVIDIA_GPU=1
        integer(C_INT), parameter, public:: DEV_INTEL_MIC=2
        integer(C_INT), parameter, public:: DEV_AMD_GPU=3
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
        integer(C_INT), parameter, public:: NOPE=0                     !"NO" answer
        integer(C_INT), parameter, public:: YEP=1                      !"YES" answer
        integer(C_INT), parameter, public:: NO_COPY_BACK=0             !keeps the tensor-result on Accelerator without updating Host
        integer(C_INT), parameter, public:: COPY_BACK=1                !tensor-result will be copied back from Accelerator to Host (default)
        integer(C_INT), parameter, public:: BLAS_ON=0                  !enables BLAS
        integer(C_INT), parameter, public:: BLAS_OFF=1                 !disables BLAS
        integer(C_INT), parameter, public:: EFF_TRN_OFF=0              !disables efficient tensor transpose algorithm
        integer(C_INT), parameter, public:: EFF_TRN_ON=1               !enables efficient tensor transpose algorithm
        integer(C_INT), parameter, public:: EVENTS_OFF=0               !disables CUDA event recording
        integer(C_INT), parameter, public:: EVENTS_ON=1                !enables CUDA event recording (default)
#ifndef NO_PHI
!DIR$ ATTRIBUTES OFFLOAD:mic:: NOPE,YEP,NO_COPY_BACK,COPY_BACK,BLAS_ON,BLAS_OFF,EFF_TRN_OFF,EFF_TRN_ON
!DIR$ ATTRIBUTES ALIGN:128:: NOPE,YEP,NO_COPY_BACK,COPY_BACK,BLAS_ON,BLAS_OFF,EFF_TRN_OFF,EFF_TRN_ON
#endif

!CUDA TASK STATUS (keep consistent with tensor_algebra.h):
        integer(C_INT), parameter, public:: CUDA_TASK_ERROR=-1
        integer(C_INT), parameter, public:: CUDA_TASK_EMPTY=0
        integer(C_INT), parameter, public:: CUDA_TASK_SCHEDULED=1
        integer(C_INT), parameter, public:: CUDA_TASK_STARTED=2
        integer(C_INT), parameter, public:: CUDA_TASK_INPUT_THERE=3
        integer(C_INT), parameter, public:: CUDA_TASK_OUTPUT_THERE=4
        integer(C_INT), parameter, public:: CUDA_TASK_COMPLETED=5

!TENSOR BLOCK STORATE LAYOUT:
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
        end module tensor_algebra
