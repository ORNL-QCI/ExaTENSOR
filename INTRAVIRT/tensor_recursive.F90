!ExaTENSOR: Recursive tensors
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/01/20

!Copyright (C) 2014-2017 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2017 Oak Ridge National Laboratory (UT-Battelle)

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
        use tensor_algebra !includes dil_basic
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
 !Tensor naming:
        integer(INTD), parameter, public:: TEREC_MAX_TENSOR_NAME_LEN=32 !max length of the alphanumeric_ tensor name
 !Index restriction kinds:
        integer(INTD), parameter, public:: TEREC_IND_RESTR_NONE=0 !no restrictions
        integer(INTD), parameter, public:: TEREC_IND_RESTR_LT=1   !indices within the group are < ordered: i1 < i2 < i3
        integer(INTD), parameter, public:: TEREC_IND_RESTR_GT=2   !indices within the group are > ordered: i1 > i2 > i3
        integer(INTD), parameter, public:: TEREC_IND_RESTR_LE=3   !indices within the group are <= ordered: i1 <= i2 <= i3
        integer(INTD), parameter, public:: TEREC_IND_RESTR_GE=4   !indices within the group are >= ordered i1 >= i2 >= i3
!TYPES:
 !Tensor signature:
        type, public:: tens_signature_t
         character(:), allocatable, private:: char_name    !character tensor name (alphanumeric_)
         integer(INTD), private:: num_dim=-1               !number of tensor dimensions (tensor order or tensor rank)
         integer(INTL), allocatable, private:: subspace(:) !subspace id for each tensor dimension
        end type tens_signature_t
 !Tensor shape:
        type, public:: tens_shape_t
         integer(INTD), private:: num_dim=-1                 !number of tensor dimensions (tensor order or tensor rank)
         integer(INTD), private:: num_groups=0               !number of defined index restriction groups
         integer(INTL), allocatable, private:: dim_extent(:) !tensor dimension extents (resolution)
         integer(INTD), allocatable, private:: dim_group(:)  !tensor dimension groups (group 0 is default, meaning no restrictions)
         integer(INTD), allocatable, private:: group_spec(:) !specification of the restriction kind for each index restriction group
        end type tens_shape_t
 !Tensor body:
        type, public:: tens_body_t
         
        end type tens_body_t
 !Recursive tensor:
        type, public:: tens_rcrsv_t
         type(tens_signature_t), private:: signature
         type(tens_shape_t), private:: shape
         type(tens_body_t), private:: body
        end type tens_rcrsv_t
!INTERFACES:

!VISIBILITY:

!DATA:

       contains
!IMPLEMENTATION:
!--------------------------------------------------------

       end module tensor_recursive
