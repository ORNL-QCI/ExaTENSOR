!ExaTensor: Recursive data/task decomposition for tensors and tensor operations.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2016/11/10

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

       module tensor_recursive
!Acronyms:
! # NAT: Node Aggregation Tree.
! # SAT: Subspace Aggregation Tree.
! # MUD: Maximally Uniform Distribution.
        use virta
        use hardware
        use subspaces
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6
        integer(INTD), private:: DEBUG=0
        logical, private:: VERBOSE=.true.
!TYPES:
 !Tensor specification:
        type, public:: tensor_spec_t
         integer(INTD), private:: num_dim=-1               !number of tensor dimensions
         integer(INTD), allocatable, private:: subspace(:) !subspace ID for each dimension
         integer(INTD), allocatable, private:: group(:)    !restriction group (>=0) each dimension belongs to (0 is the trivial group)
         integer(INTD), allocatable, private:: label(:)    !index label for each dimension (defined in tensor operations)
        end type tensor_spec_t
 !Tensor operand:
        type, public:: tensor_operand_t
         type(tensor_spec_t), private:: spec !tensor specification
        contains
         procedure, public:: split=>TensOperSplit !splits the tensor operand into smaller tensor operands according to SAT
        end type tensor_operand_t
 !Tensor operation:
        type, public:: tensor_operation_t
         integer(INTD), private:: op_code=INSTR_NULL               !operation code
         integer(INTD), private:: num_operands=0                   !number of operands
         type(tensor_operand_t), allocatable, private:: operand(:) !operands
        end type tensor_operation_t
!INTERFACES:

!VISIBILITY:
        private TensOperSplit
!DATA:

       contains
!IMPLEMENTATION:
!--------------------------------------------------------
        subroutine TensOperSplit(this)
!Splits a tensor into smaller tensors according to SAT.
         implicit none
         class(tensor_operand_t), intent(in):: this

        end subroutine TensOperSplit

       end module tensor_recursive
