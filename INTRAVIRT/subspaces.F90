!Infrastructure for a recursive adaptive Hilbert space decomposition.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2015/09/24

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
        use dil_kinds
        implicit none
        private
!PARAMETERS:
 !Output:
        integer, private:: CONS_OUT=6     !default output for this module
        logical, private:: VERBOSE=.true. !verbosity for errors
        logical, private:: DEBUG=.true.   !debugging mode
!TYPES:
 !Hierarchical index:
        type, public:: HierIndex_t
         integer(INTL), private:: SubspaceID !subspace ID (registered ID)
         integer(INTL), private:: ResLevel   !subspace resolution level (number of basis components)
         integer(INTL), private:: Component  !subspace component at a given level of resolution: [0..ResLevel-1]
        end type HierIndex_t
 !Hilbert space with adaptive subspaces:
        type, public:: HilbertSpaceAd_t
         
        end type HilbertSpaceAd_t
!DATA:

!VISIBILITY:

        contains
!METHODS:

        end module subspaces
