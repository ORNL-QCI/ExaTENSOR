!ExaTensor: Parallel Virtual Processing for Scale-Adaptive Tensor Algebra
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2016/04/17

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
        use talsh
        use distributed
        use service_mpi, only: get_memory_status
        implicit none
        public
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !default output for this module
        integer(INTD), private:: DEBUG=0    !debugging mode
        logical, private:: VERBOSE=.true.   !verbosity for errors
 !Kinds of MPI processes (process roles):
        integer(INTD), parameter, public:: EXA_NO_ROLE=0              !undefined role
        integer(INTD), parameter, public:: EXA_MANAGER=1              !manager process (global root is a manager as well)
        integer(INTD), parameter, public:: EXA_WORKER=2               !worker process (aka C-process)
        integer(INTD), parameter, public:: EXA_HELPER=3               !helper process
        integer(INTD), parameter, public:: EXA_MAX_WORK_GROUP_SIZE=64 !maximal size of a work group (max number of workers per manager)
 !Elementary tensor instruction (ETI) granularity:
        real(8), public:: EXA_FLOPS_MEDIUM=1d9 !minimal number of Flops to consider the operation as medium-cost
        real(8), public:: EXA_FLOPS_HEAVY=1d11 !minimal number of Flops to consider the operation as heavy-cost
        real(8), public:: EXA_COST_TO_SIZE=1d1 !minimal cost to size ratio to consider the operation compute intensive
 !Tensor algebra virtual processor:
  !Tensor naming:
        integer(INTD), parameter, public:: TENSOR_NAME_LEN=32    !max number of characters used for tensor names
  !Tensor instruction status (any negative value will correspond to a failure, designating the error code):
        integer(INTD), parameter, public:: INSTR_NULL=0          !uninitialized instruction (empty)
        integer(INTD), parameter, public:: INSTR_FRESH=1         !newly arrived instruction
        integer(INTD), parameter, public:: INSTR_DATA_WAIT=2     !instruction is waiting for the remote data to arrive
        integer(INTD), parameter, public:: INSTR_READY_TO_EXEC=3 !instruction is ready to be executed (input data has arrived)
        integer(INTD), parameter, public:: INSTR_SCHEDULED=4     !instruction has been placed into the execution queue on a specific CU
        integer(INTD), parameter, public:: INSTR_ISSUED=5        !instruction has been issued for execution to a computing unit (CU)
        integer(INTD), parameter, public:: INSTR_COMPLETED=6     !instruction has completed (the result may still need a remote upload)
        integer(INTD), parameter, public:: INSTR_RETIRED=7       !instruction can safely be removed from the queue
        integer(INTD), parameter, public:: INSTR_STATUSES=8      !total number of instruction statuses
  !Tensor instruction code (opcode):
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
 !Current process role:
        integer(INTD), public:: my_role=EXA_NO_ROLE !role of this MPI process (set at run-time)
        integer(INTD), public:: my_group=-1         !computing group the process belongs to (set at run-time): [0..MAX]
        integer(INTD), public:: my_group_size=0     !size of the computing group the process belongs to (set at runtime): [1..EXA_MAX_WORK_GROUP_SIZE]
        integer(INTD), public:: my_group_index=-1   !process ID within its computing group (set at run-time): [0..my_group_size-1]
!VISIBILITY:

!IMPLEMENTATION:

       end module virta
