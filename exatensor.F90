!Parallel Virtual Processing for Scale-Adaptive Tensor Algebra:
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2015/06/15
        module exatensor
        use tensor_algebra_cpu_phi
        use distributed
        use subspaces
        use lists
        use dictionary
        use multords
        use extern_names
        implicit none
        public
!PARAMETERS:
 !Output:
        integer, private:: CONS_OUT=6     !default output for this module
        logical, private:: VERBOSE=.true. !verbosity for errors
        logical, private:: DEBUG=.true.   !debugging mode
 !Numeric:
        real(4), parameter, public:: EPS4=epsilon(1.0) !single precision epsilon
        real(8), parameter, public:: EPS8=epsilon(1d0) !double precision epsilon
        real(8), parameter, public:: ZERO_THRESH=1d-11 !numerical comparison threshold: should account for possible round-off errors
 !Kinds of MPI processes (process roles):
        integer(INTD), parameter, public:: GLOBAL_ROOT=1             !global root
        integer(INTD), parameter, public:: LOCAL_ROOT=2              !local root
        integer(INTD), parameter, public:: C_PROCESS_PRIVATE=3       !computing process private to his root
        integer(INTD), parameter, public:: C_PROCESS_SHARED=4        !computing process shared by multiple roots
        integer(INTD), parameter, public:: D_PROCESS=5               !I/O operating process
        integer(INTD), parameter, public:: C_PROCS_PER_LOCAL_ROOT=32 !number of C-processes per local root
 !Elementary tensor instruction (ETI) granularity:
        real(8), public:: FLOPS_MEDIUM=1d9 !minimal number of Flops to consider the operation as medium-cost
        real(8), public:: FLOPS_LARGE=1d11 !minimal number of Flops to consider the operation as large-cost
        real(8), public:: COST_TO_SIZE=1d1 !minimal cost to size ratio to consider the operation cost-efficient
 !Tensor algebra virtual processor:
  !Tensor naming:
        integer(INTD), parameter, public:: TENSOR_NAME_LEN=32    !max number of characters used for tensor names
  !Tensor instruction status (any negative value will correspond to a failure, designating the error code):
        integer(INTD), parameter, public:: INSTR_NULL=0          !uninitialized instruction (empty)
        integer(INTD), parameter, public:: INSTR_FRESH=1         !newly arrived instruction
        integer(INTD), parameter, public:: INSTR_DATA_WAIT=2     !instruction is waiting for the input data to arrive
        integer(INTD), parameter, public:: INSTR_READY_TO_EXEC=3 !instruction is ready to be executed (input data has arrived)
        integer(INTD), parameter, public:: INSTR_SCHEDULED=4     !instruction has been placed into an execution queue on a specific CU
        integer(INTD), parameter, public:: INSTR_ISSUED=5        !instruction has been issued for execution to a computing unit (CU)
        integer(INTD), parameter, public:: INSTR_COMPLETED=6     !instruction has completed (the result may still need a remote upload)
        integer(INTD), parameter, public:: INSTR_RETIRED=7       !instruction can safely be removed from the queue
        integer(INTD), parameter, public:: INSTR_STATUSES=8      !total number of instruction statuses
  !Tensor instruction code:
        integer(INTD), parameter, public:: INSTR_TENSOR_INIT=1
        integer(INTD), parameter, public:: INSTR_TENSOR_NORM1=2
        integer(INTD), parameter, public:: INSTR_TENSOR_NORM2=3
        integer(INTD), parameter, public:: INSTR_TENSOR_MIN=4
        integer(INTD), parameter, public:: INSTR_TENSOR_MAX=5
        integer(INTD), parameter, public:: INSTR_TENSOR_SCALE=6
        integer(INTD), parameter, public:: INSTR_TENSOR_SLICE=7
        integer(INTD), parameter, public:: INSTR_TENSOR_INSERT=8
        integer(INTD), parameter, public:: INSTR_TENSOR_TRACE=9
        integer(INTD), parameter, public:: INSTR_TENSOR_COPY=10
        integer(INTD), parameter, public:: INSTR_TENSOR_ADD=11
        integer(INTD), parameter, public:: INSTR_TENSOR_CMP=12
        integer(INTD), parameter, public:: INSTR_TENSOR_CONTRACT=13
!TYPES:

!DATA:

!FUNCTION VISIBILITY:

!METHODS:

        end module exatensor
