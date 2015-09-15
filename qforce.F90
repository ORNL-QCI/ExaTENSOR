!Quantum Chemistry specific module:
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2015/04/28
       module qforce
        use exatensor
        implicit none
        public
!PARAMETERS:

!TYPES:
 !Nuclei vibration descriptor:
        type vibration_t
         real(8):: frequency
         integer(INTD):: symmetry
        end type vibration_t
 !Electronic state descriptor:
        type el_state_t
         integer(INTL):: nelectrons            !number of electrons in this state
         integer(INTL):: degrees               !number of degrees of freedom
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
         integer(INTL):: natoms                       !number of atoms in the molecular system
         integer(INTD):: nstates                      !number of electronic states computed (Hamiltonian roots)
         type(el_state_t), allocatable:: el_states(:) !description of the electronic states
         integer(INTD):: basis_set_type=0             !0 - Gaussian AO basis set; [1..] - alternative basis sets
         real(8), allocatable:: atoms(:,:)            !0 - atomic number, 1..3 - X,Y,Z coordinates (a.u.) for each atom [1..natoms]
         integer(INTD), allocatable:: basis(:)        !basis set ID for each atom
         integer(INTD), allocatable:: tags(:)         !atomic tags
         character(1024):: mol_data_filename          !name of the file containing the rest of the molecular data
        end type molecule_t
!DATA:

!FUNCTION VISIBILITY:

       end module qforce
