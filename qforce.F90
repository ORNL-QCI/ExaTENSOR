       module qforce !PROGRAM SPECIFIC MODULE
        use, intrinsic:: ISO_C_BINDING
        use STSUBS
        use combinatoric
        use service
        use extern_names
        use c_process
!PARAMETERS:
 !Numeric:
        real(4), parameter:: eps4=epsilon(1.0)     !single precision epsilon
        real(8), parameter:: eps8=epsilon(1d0)     !double precision epsilon
        real(8), parameter:: zero_thresh=1d-11     !numerical comparison threshold: should account for possible round-off errors
 !Kinds of MPI processes (process roles):
        integer, parameter:: global_root=1             !global root
        integer, parameter:: local_root=2              !local root
        integer, parameter:: c_process_private=3       !computing process private to a unit cell
        integer, parameter:: c_process_shared=4        !computing process shared by multiple unit cells
        integer, parameter:: d_process=5               !I/O operating process
        integer, parameter:: a_process=6               !accelerating process
        integer, parameter:: c_procs_per_local_root=32 !number of C-processes per local root
!TYPES:
 !Nuclei vibration descriptor:
        type vibration_t
         real(8):: frequency
         integer:: symmetry
        end type vibration_t
 !Electronic state descriptor:
        type el_state_t
         integer:: nelectrons                  !number of electrons in this state
         real(8):: net_charge                  !net charge of the molecule in this state
         real(8):: spin                        !spin: <S^2>=spin*(spin+1)
         real(8):: sz                          !Z-projection of the spin
         real(8):: reference_energy            !reference energy
         real(8):: total_energy                !total energy
         real(8), allocatable:: gradients(:)   !gradients of the energy
         real(8), allocatable:: hessian(:,:)   !hessian of the energy
         type(vibration_t), allocatable:: vibrations(:) !vibrational structure information
        end type el_state_t
 !Molecular descriptor:
        type molecule_t
         integer:: natoms                             !number of atoms in the molecular system
         integer:: number_of_states                   !number of electronic states computed (Hamiltonian roots)
         type(el_state_t), allocatable:: el_states(:) !description of the electronic states
         integer:: basis_set_type=0                   !0 - Gaussian AO basis set; [1..] - alternative basis sets
         real(8), allocatable:: atoms(:,:)            !0 - atomic number, 1..3 - X,Y,Z coordinates (a.u.) for each atom [1..natoms]
         integer, allocatable:: basis(:)              !basis set ID for each atom
         integer, allocatable:: tags(:)               !atomic tags
         character(128):: mol_data_filename           !name of the file containing the rest of the molecular data
        end type molecule_t
!DATA:

!INTERFACES:

       end module qforce
