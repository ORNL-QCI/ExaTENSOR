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
        public
!PARAMETERS:
 !Output:
        integer, private:: CONS_OUT=6     !default output for this module
        logical, private:: VERBOSE=.true. !verbosity for errors
        logical, private:: DEBUG=.true.   !debugging mode
 !Numeric:
        real(4), parameter, public:: eps4=epsilon(1.0) !single precision epsilon
        real(8), parameter, public:: eps8=epsilon(1d0) !double precision epsilon
        real(8), parameter, public:: zero_thresh=1d-11 !numerical comparison threshold: should account for possible round-off errors
 !Kinds of MPI processes (process roles):
        integer(INTD), parameter, public:: global_root=1             !global root
        integer(INTD), parameter, public:: local_root=2              !local root
        integer(INTD), parameter, public:: c_process_private=3       !computing process private to a unit cell
        integer(INTD), parameter, public:: c_process_shared=4        !computing process shared by multiple unit cells
        integer(INTD), parameter, public:: d_process=5               !I/O operating process
        integer(INTD), parameter, public:: a_process=6               !accelerating process
        integer(INTD), parameter, public:: c_procs_per_local_root=32 !number of C-processes per local root
!TYPES:

!DATA:

!FUNCTION VISIBILITY:

!METHODS:

        end module exatensor
