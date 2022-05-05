!Quantum Chemistry specific module
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2016/03/24

!Copyright (C) 2014-2022 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2022 Oak Ridge National Laboratory (UT-Battelle)

!LICENSE: BSD 3-Clause

       module qforce
        use dil_basic
        use exatensor
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6
        integer(INTD), private:: DEBUG=0
        logical, private:: VERBOSE=.true.
!TYPES:
 !Nuclear vibration descriptor:
        type, public:: vibration_t
         real(8), private:: frequency
         real(8), private:: strength
         integer(INTD), private:: symmetry
        end type vibration_t
 !Electronic state descriptor:
        type, public:: el_state_t
         integer(INTL), private:: nelectrons          !number of electrons in this state
         integer(INTL), private:: degrees             !number of degrees of freedom
         real(8), private:: net_charge                !net charge of the molecule in this state
         real(8), private:: spin                      !spin: <S^2>=spin*(spin+1)
         real(8), private:: sz                        !Z-projection of the spin
         real(8), private:: reference_energy          !reference energy
         real(8), private:: total_energy              !total energy
         real(8), allocatable, private:: gradients(:) !gradients of the energy
         real(8), allocatable, private:: hessian(:,:) !hessian of the energy
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
!INTERFACES:

!DATA:

!VISIBILITY:

!IMPLEMENTATION:

       end module qforce
