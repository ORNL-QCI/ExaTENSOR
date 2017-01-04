!ExaTENSOR: Recursive tensors
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/01/04

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
! # SAT: Subspace Aggregation Tree.
! # MUD: Maximally Uniform Distribution.
        use dil_basic
        use gfc_base
        use gfc_tree
        use subspaces
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !output device
        integer(INTD), private:: DEBUG=0    !debugging level
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
!TYPES:
 !Tensor specification:
        type, public:: tensor_spec_t
         integer(INTD), private:: num_dim=-1               !number of tensor dimensions
         integer(INTD), allocatable, private:: subspace(:) !subspace ID for each dimension
         integer(INTD), allocatable, private:: group(:)    !restriction group (>=0) each dimension belongs to (0 is the unrestricted group)
         integer(INTD), allocatable, private:: label(:)    !index label for each dimension (defined in tensor operations)
        end type tensor_spec_t
!INTERFACES:

!VISIBILITY:

!DATA:

       contains
!IMPLEMENTATION:
!--------------------------------------------------------

       end module tensor_recursive
