!Infrastructure for a recursive adaptive Hilbert space decomposition.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2016/03/16

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
        use dictionary
        implicit none
        private
!PARAMETERS:
 !Output:
        integer, private:: CONS_OUT=6     !default output for this module
        logical, private:: VERBOSE=.true. !verbosity for errors
        integer, private:: DEBUG=0        !debugging mode level (0:none)
 !Basis set kind:
        integer(INTD), parameter, public:: BASIS_NONE=0    !no basis set
        integer(INTD), parameter, public:: BASIS_GAUSS=1   !Gaussian basis set
        integer(INTD), parameter, public:: BASIS_WAVELET=2 !wavelet basis set
!TYPES:
 !Abstract basis function:
        type, abstract, public:: AbsBasisFunc_t
         integer(INTD), private:: RealSpaceDim       !dimensionality of the real space where the basis functions reside
         real(8), allocatable, private:: Center(:)   !coordinates of the effective function support center in the real space
         real(8), allocatable, private:: Extent(:,:) !effective extents of the basis function support in each direction {positive/negative}x{dimensions}
        end type AbsBasisFunc_t
 !Gaussian basis function:
        type, extends(AbsBasisFunc_t), public:: BasisFuncGauss_t
         integer(INTD), private:: NumPrims=0              !number of contracted primitives
         real(8), allocatable, private:: Exponents(:)     !primitive exponents
         complex(8), allocatable, private:: ContrCoefs(:) !primitive contraction coefficients
        end type BasisFuncGauss_t
 !Typed basis function:
        type, public:: BasisFunc_t
         integer(INTD), private:: BasisKind=BASIS_NONE               !basis kind
         class(AbsBasisFunc_t), pointer, private:: BasisFunc=>NULL() !pointer to a basis function of this kind
        end type BasisFunc_t
 !Subspace basis:
        type, public:: SubspaceBasis_t
         integer(INTL), private:: BasisDim=0                !number of basis functions
         type(BasisFunc_t), allocatable, private:: Basis(:) !basis functions [1..BasisDim]
        end type SubspaceBasis_t
 !Subspace:
        type, public:: Subspace_t
         integer(INTL), private:: SubspaceID=-1      !subspace ID (registered ID): -1 means undefined
         integer(INTL), private:: MaxResLevel=0      !max resolution level (max dimension): 0 means undefined
         real(8), allocatable, private:: Center(:)   !coordinates of the effective subspace support center in the real space
         real(8), allocatable, private:: Extent(:,:) !effective extents of the subspace support in each direction {positive/negative}x{dimensions}
         type(SubspaceBasis_t), allocatable, private:: Basis(:) !pointer to basis sets for each resolution level [1..MaxResLevel]: Defined only for the terminal subspaces (leaves)
        end type Subspace_t
 !Hierarchical component index:
        type, public:: HierIndex_t
         integer(INTL), private:: SubspaceID=-1 !terminal subspace ID (registered ID): -1 means undefined
         integer(INTL), private:: ResLevel=0    !subspace resolution level 1<=ResLevel<=MaxResLevel: 0 means undefined
         integer(INTL), private:: Component=0   !subspace component number at given level of resolution: [1..ResLevel], 0 means undefined
        end type HierIndex_t
!DATA:
 !Bank of basis sets:
        type(dict_t), private:: BasisBank !KEY = Subspace ID, VALUE = type(SubspaceBasis_t): Specific subspace basis (multiresolution in general)
!VISIBILITY:

        contains
!METHODS:

        end module subspaces
