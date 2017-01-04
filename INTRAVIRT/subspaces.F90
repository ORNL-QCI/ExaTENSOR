!Infrastructure for a recursive adaptive vector space decomposition.
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

       module subspaces
        use dil_basic
        use gfc_base
        use gfc_tree
        implicit none
        private
!PARAMETERS:
 !Output:
        integer, private:: CONS_OUT=6     !default output for this module
        integer, private:: DEBUG=0        !debugging mode level (0:none)
        logical, private:: VERBOSE=.TRUE. !verbosity for errors
 !Basis set kind:
        integer(INTD), parameter, public:: BASIS_NONE=0    !No basis set
        integer(INTD), parameter, public:: BASIS_GAUSS=1   !Gaussian basis set
        integer(INTD), parameter, public:: BASIS_WAVELET=2 !Wavelet basis set
!TYPES:
 !Real space vector:
        type, public:: real_vec_t
         integer(INTD):: num_dim          !number of dimensions
         real(8), allocatable:: coords(:) !components
        end type real_vec_t
 !1D extent (segment[min_loc:max_loc]):
        type, public:: extent1d_t
         real(8):: min_loc !minimum
         real(8):: max_loc !maximum
        end type extent1d_t
 !Abstract basis function:
        type, abstract, public:: basis_func_abs_t
         type(real_vec_t), private:: center                 !coordinates of the center of the effective function support in the real space
         type(extent1d_t), allocatable, private:: extent(:) !effective extents of the basis function support in each dimension
        end type basis_func_abs_t
 !Gaussian basis function:
        type, extends(basis_func_abs_t), public:: basis_func_gauss_t
         integer(INTD), private:: num_prims=0         !number of contracted primitives
         integer(INTD), private:: orb_moment=-1       !orbital momentum (0,1,2,3,...)
         real(8), allocatable, private:: exponents(:) !primitive exponents
         complex(8), allocatable, private:: coefs(:)  !primitive contraction coefficients
        end type basis_func_gauss_t
 !Typed basis function:
        type, public:: basis_func_t
         integer(INTD), private:: basis_kind=BASIS_NONE                 !basis kind
         class(basis_func_abs_t), pointer, private:: basis_func=>NULL() !pointer to a basis function of this kind
        end type basis_func_t
 !Subspace basis:
        type, public:: subspace_basis_t
         integer(INTD), private:: supp_dim=0                      !dimensionality of the real space support on which the basis functions reside
         integer(INTL), private:: space_dim=0                     !number of basis functions
         type(basis_func_t), allocatable, private:: basis_func(:) !basis functions [1..space_dim]
        end type subspace_basis_t
 !Subspace:
        type, public:: subspace_t
         integer(INTL), private:: subspace_id=-1                 !subspace ID (registered ID): must be non-negative, -1 means undefined
         integer(INTL), private:: max_resolution=0               !max resolution level (max dimension): 0 means undefined
         type(real_vec_t), private:: center                      !coordinates of the center of the effective subspace support in real space
         type(extend1d_t, allocatable, private:: extent(:)       !effective extents of the subspace support in each real space dimension
         type(subspace_basis_t), allocatable, private:: basis(:) !pointer to basis sets for each resolution level [1..MaxResLevel]: Defined only for the terminal subspaces
        end type subspace_t
 !Hierarchical composite index:
        type, public:: h_index_t
         integer(INTL), private:: subspace_id=-1 !subspace ID (registered ID): -1 means undefined
         integer(INTL), private:: resolution=0   !subspace resolution level 1<=ResLevel<=MaxResLevel: 0 means undefined
         integer(INTL), private:: component=0    !subspace component number at the given level of resolution: [1..ResLevel], 0 means undefined
        end type h_index_t
!DATA:

!VISIBILITY:

        contains
!IMPLEMENTATION:

       end module subspaces
