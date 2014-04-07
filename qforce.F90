       module qforce !PROGRAM SPECIFIC MODULE
        use service
        use STSUBS
        use combinatoric
        use tensor_algebra
        use extern_names
!PARAMETERS:
 !Numeric:
        real(4), parameter:: eps4=epsilon(1.0)     !single precision epsilon
        real(8), parameter:: eps8=epsilon(1d0)     !double precision epsilon
        real(8), parameter:: zero_thresh=1d-11     !numerical comparison threshold: should account for possible round-off errors
 !Kinds of MPI processes (process role):
        integer, parameter:: global_root=1         !global root
        integer, parameter:: local_root=2          !local root
        integer, parameter:: c_process_private=3   !computing process private to a unit cell
        integer, parameter:: c_process_shared=4    !computing process shared by multiple unit cells
        integer, parameter:: d_process=5           !I/O operating process
        integer, parameter:: a_process=6           !accelerating process
 !Parallelization configuration:
        integer, parameter:: c_procs_per_local_root=32 !number of C-processes per local root
 !Tensor instruction status (any negative value will correspond to failure, designating the error code):
        integer, parameter:: instr_null=0          !uninitialized instruction
        integer, parameter:: instr_data_wait=1     !instruction is waiting for data to arrive
        integer, parameter:: instr_ready_to_exec=2 !instruction is ready to be executed (data has arrived)
        integer, parameter:: instr_scheduled=3     !instruction is being executed
        integer, parameter:: instr_complete=4      !instruction has completed
 !Tensor instruction code:
        integer, parameter:: instr_tensor_init=1
        integer, parameter:: instr_tensor_norm2=2
        integer, parameter:: instr_tensor_min=3
        integer, parameter:: instr_tensor_max=4
        integer, parameter:: instr_tensor_scale=5
        integer, parameter:: instr_tensor_slice=6
        integer, parameter:: instr_tensor_insert=7
        integer, parameter:: instr_tensor_trace=8
        integer, parameter:: instr_tensor_copy=9
        integer, parameter:: instr_tensor_add=10
        integer, parameter:: instr_tensor_cmp=11
        integer, parameter:: instr_tensor_contract=12
!TYPES:
 !Locally present tensor argument (negative buf_entry_xxx means that no argument buffer space is used by the argument):
        type tens_arg_t
         type(tensor_block_t):: tens_blck_f !tensor_block (QFORCE Fortran)
         integer(C_INT):: buf_entry_host    !host argument buffer entry number where the tensor block resides as a packet (negative: tensor block is not in any of argument buffers)
         type(C_PTR):: tens_blck_c          !C pointer to tensBlck_t (QFORCE C), see "tensor_algebra_gpu_nvidia.h"
        end type tens_arg_t
 !Dispatched tensor instruction:
        type tensor_instruction_t
         integer:: instr_status                     !instruction status (see above)
         integer:: instr_priority                   !instruction priority
         integer:: instr_code                       !instruction code (see above)
         real(8):: instr_cost                       !approx. instruction computational cost (FLOPs)
         real(8):: instr_size                       !approx. instruction memory demands (Bytes)
         real(4):: instr_time_beg                   !time when the instruction was scheduled
         real(4):: instr_time_end                   !time when the instruction completed
         integer:: args_ready                       !a bit of this word is set to 1 when the corresponding argument is available in local memory
         type(tens_arg_t), pointer:: tens_arg0      !pointer to the tensor block argument #0
         type(tens_arg_t), pointer:: tens_arg1      !pointer to the tensor block argument #1
         type(tens_arg_t), pointer:: tens_arg2      !pointer to the tensor block argument #2
         type(tens_arg_t), pointer:: tens_arg3      !pointer to the tensor block argument #3
         type(C_PTR):: cuda_task                    !pointer to CUDA task associated with this tensor instruction (if any)
        end type tensor_instruction_t
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

!-----------------------------------------------------------------------------------
!INTERFACES:
!        interface
         
!        end interface
!-----------------------------------------------------------------------------------
!       contains
!SUBROUTINES/FUNCTIONS:

       end module qforce
