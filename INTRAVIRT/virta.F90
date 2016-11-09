!ExaTensor: Parallel Virtual Processing for Scale-Adaptive Tensor Algebra
!This module provides the infrastructure for the tensor algebra processor.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2016/11/09

!Copyright (C) 2014-2016 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2016 Oak Ridge National Laboratory (UT-Battelle)

!This file is part of ExaTensor.

!ExaTensor is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published
!by the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.

!ExaTensor is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!GNU Lesser General Public License for more details.

!You should have received a copy of the GNU Lesser General Public License
!along with ExaTensor. If not, see <http://www.gnu.org/licenses/>.

       module virta
        use dil_basic
        use pack_prim
        use talsh
        use distributed
#ifndef NO_LINUX
        use service_mpi, only: get_memory_status
#endif
        use dsvp_base
        implicit none
        public
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !default output for this module
        integer(INTD), private:: DEBUG=0    !debugging mode
        logical, private:: VERBOSE=.true.   !verbosity for errors
 !Errors:
        integer(INTD), parameter, public:: EXA_SUCCESS=0              !success
        integer(INTD), parameter, public:: EXA_ERROR=-666             !generic error (or trap error)
        integer(INTD), parameter, public:: EXA_INVALID_ARGS=-1        !invalid arguments passed to a procedure
 !Kinds of MPI processes (process roles):
        integer(INTD), parameter, public:: EXA_NO_ROLE=0              !undefined role
        integer(INTD), parameter, public:: EXA_MANAGER=1              !manager process (global root is a manager as well)
        integer(INTD), parameter, public:: EXA_WORKER=2               !worker process (aka C-process)
        integer(INTD), parameter, public:: EXA_HELPER=3               !helper process
        integer(INTD), parameter, public:: EXA_MAX_WORK_GROUP_SIZE=64 !maximal size of a work group (max number of workers per manager)
 !Elementary tensor instruction (ETI) granularity:
        real(8), public:: EXA_FLOPS_MEDIUM=1d9 !minimal number of Flops to consider the operation as medium-cost
        real(8), public:: EXA_FLOPS_HEAVY=1d11 !minimal number of Flops to consider the operation as heavy-cost
        real(8), public:: EXA_COST_TO_SIZE=1d1 !minimal cost (Flops) to size (Words) ratio to consider the operation compute intensive
 !Tensor algebra virtual processor:
  !Tensor naming:
        integer(INTD), parameter, public:: TENSOR_NAME_LEN=32    !max number of characters used for tensor names
  !Tensor instruction status (any negative value will correspond to a failure, designating the error code):
        integer(INTD), parameter, public:: INSTR_NULL=0          !uninitialized instruction (empty)
        integer(INTD), parameter, public:: INSTR_FRESH=1         !newly arrived instruction
        integer(INTD), parameter, public:: INSTR_DATA_WAIT=2     !instruction is waiting for the remote data to arrive
        integer(INTD), parameter, public:: INSTR_READY_TO_EXEC=3 !instruction is ready to be executed (input data has arrived)
        integer(INTD), parameter, public:: INSTR_SCHEDULED=4     !instruction has been dispatched to the execution queue on a specific CU
        integer(INTD), parameter, public:: INSTR_ISSUED=5        !instruction has been issued for execution to a computing unit (CU)
        integer(INTD), parameter, public:: INSTR_COMPLETED=6     !instruction has completed computation (the result may still need a remote upload)
        integer(INTD), parameter, public:: INSTR_RETIRED=7       !instruction can safely be retired
        integer(INTD), parameter, public:: INSTR_STATUSES=8      !total number of instruction statuses
  !Tensor instruction code (opcode):
        integer(INTD), parameter, public:: INSTR_TENSOR_INIT=1      !tensor initialization (definition)
        integer(INTD), parameter, public:: INSTR_TENSOR_NORM1=2     !tensor 1-norm
        integer(INTD), parameter, public:: INSTR_TENSOR_NORM2=3     !tensor 2-norm
        integer(INTD), parameter, public:: INSTR_TENSOR_MIN=4       !tensor min element
        integer(INTD), parameter, public:: INSTR_TENSOR_MAX=5       !tensor max element
        integer(INTD), parameter, public:: INSTR_TENSOR_SCALE=6     !tensor scaling (multiplication by a number)
        integer(INTD), parameter, public:: INSTR_TENSOR_SLICE=7     !tensor slicing (taking a slice of a tensor)
        integer(INTD), parameter, public:: INSTR_TENSOR_INSERT=8    !tensor insertion (inserting a slice in a tensor)
        integer(INTD), parameter, public:: INSTR_TENSOR_TRACE=9     !tensor trace (tracing over some/all tensor indices)
        integer(INTD), parameter, public:: INSTR_TENSOR_COPY=10     !tensor copy (copying a tensor into another tensor with an optional index permutation)
        integer(INTD), parameter, public:: INSTR_TENSOR_ADD=11      !tensor addition
        integer(INTD), parameter, public:: INSTR_TENSOR_CMP=12      !tensor comparison
        integer(INTD), parameter, public:: INSTR_TENSOR_CONTRACT=13 !tensot contraction (tensor product is considered a tensor contraction as well)
!TYPES:

!DATA:
 !Current role of the process:
        integer(INTD), protected:: my_role=EXA_NO_ROLE !role of this MPI process (set at run-time)
        integer(INTD), protected:: my_group=-1         !computing group the process belongs to (set at run-time): [0..MAX]
        integer(INTD), protected:: my_group_size=0     !size of the computing group the process belongs to (set at runtime): [1..EXA_MAX_WORK_GROUP_SIZE]
        integer(INTD), protected:: my_group_index=-1   !process ID within its computing group (set at run-time): [0..my_group_size-1]
!VISIBILITY:
        public exa_set_process_role

       contains
!IMPLEMENTATION:
!-------------------------------------------------
        subroutine exa_set_process_role(role,ierr)
!Sets process's role.
         implicit none
         integer(INTD), intent(in):: role  !in: role of the current MPI process
         integer(INTD), intent(out):: ierr !out: error code

         ierr=EXA_SUCCESS
         select case(role)
         case(EXA_MANAGER,EXA_WORKER,EXA_HELPER)
          my_role=role
         case(EXA_NO_ROLE)
          ierr=EXA_ERROR
         case default
          ierr=EXA_INVALID_ARGS
         end select
         return
        end subroutine exa_set_process_role

       end module virta
