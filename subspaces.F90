!Infrastructure for adaptive Hilbert space decompositions.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2015/04/27 (started 2015/03/18)
!All rights reserved!
        module subspaces
        use tensor_algebra_cpu_phi
        use distributed
        use lists
        use dictionary
        implicit none
!PARAMETERS:
 !Output:
        integer(INTD), private:: CONS_OUT=6 !default output for this module
        logical, private:: VERBOSE=.true.   !verbosity for errors
        logical, private:: DEBUG=.true.     !debugging mode
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

!FUNCTION VISIBILITY:

        contains
!METHODS:

        end module subspaces
