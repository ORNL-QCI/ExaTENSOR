!Quantum Chemistry specific module
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2015/09/24
       module qforce
        use dil_kinds
        use exatensor
        implicit none
        private
!PARAMETERS:

!TYPES:
 !Nuclear vibration descriptor:
        type, public:: vibration_t
         real(8), private:: frequency
         integer(INTD), private:: symmetry
        end type vibration_t
 !Electronic state descriptor:
        type, public:: el_state_t
         integer(INTL), private:: nelectrons            !number of electrons in this state
         integer(INTL), private:: degrees               !number of degrees of freedom
         real(8), private:: net_charge                  !net charge of the molecule in this state
         real(8), private:: spin                        !spin: <S^2>=spin*(spin+1)
         real(8), private:: sz                          !Z-projection of the spin
         real(8), private:: reference_energy            !reference energy
         real(8), private:: total_energy                !total energy
         real(8), allocatable, private:: gradients(:)   !gradients of the energy
         real(8), allocatable, private:: hessian(:,:)   !hessian of the energy
         type(vibration_t), allocatable, private:: vibrations(:) !vibrational structure information
        end type el_state_t
 !Molecular descriptor:
        type, public:: molecule_t
         integer(INTL), private:: natoms                       !number of atoms in the molecular system
         integer(INTD), private:: nstates                      !number of electronic states computed (Hamiltonian roots)
         type(el_state_t), allocatable, private:: el_states(:) !description of the electronic states
         integer(INTD), private:: basis_set_type=0             !0 - Gaussian AO basis set; [1..] - alternative basis sets
         real(8), allocatable, private:: atoms(:,:)            !0 - atomic number, 1..3 - X,Y,Z coordinates (a.u.) for each atom [1..natoms]
         integer(INTD), allocatable, private:: basis(:)        !basis set ID for each atom
         integer(INTD), allocatable, private:: tags(:)         !atomic tags
         character(1024), private:: mol_data_filename          !name of the file containing the rest of the molecular data
        end type molecule_t
!DATA:

!VISIBILITY:

       end module qforce
