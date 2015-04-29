!Parallel Virtual Processing for Scale-Adaptive Tensor Algebra:
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2015/04/28
        module exatensor
        use tensor_algebra_cpu_phi
        use distributed
        use subspaces
        use lists
        use dictionary
        use extern_names
        implicit none
!PARAMETERS:
 !Output:
        integer, private:: CONS_OUT=6     !default output for this module
        logical, private:: VERBOSE=.true. !verbosity for errors
        logical, private:: DEBUG=.true.   !debugging mode
!TYPES:

!DATA:

!FUNCTION VISIBILITY:

!METHODS:

        end module exatensor
