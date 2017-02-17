!Hardware abstraction module
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/02/16

!Copyright (C) 2014-2017 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2017 Oak Ridge National Laboratory (UT-Battelle)

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

       module hardware
        use dil_basic
        use stsubs
        use parse_prim
        use gfc_base
        use gfc_tree
        use subspaces
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !output device
        integer(INTD), private:: DEBUG=0    !debugging level
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
!TYPES:
 !Node architecture abstraction:
  !Abstract memory unit:
        type, public:: memory_t
         integer(INTL), private:: id=-1  !memory unit id (non-negative)
         real(8), private:: capacity=0d0 !memory unit capacity in bytes (approx.)
        end type memory_t
  !Abstract data transfer link:
        type, public:: memory_channel_t
         integer(INTL), private:: id=-1   !data link id (non-negative)
         real(8), private:: latency=0d0   !average latency of the data link is seconds
         real(8), private:: bandwidth=0d0 !average bandwidth of the data link in bytes/sec
        end type memory_channel_t
  !Attached memory unit:
        type, public:: memory_access_t
         type(memory_t), pointer, private:: memory=>NULL() !pointer to the memory unit description
         type(memory_channel_t), private:: channel         !memory channel description
         logical, private:: remote                         !TRUE if the memory is not directly accessible, FALSE otherwise
        end type memory_access_t
  !Abstract compute unit:
        type, public:: compute_unit_t
         integer(INTL), private:: id=-1                              !compute unit id (non-negative)
         integer(INTD), private:: dev_kind=DEV_NULL                  !compute unit kind: {DEV_HOST, DEV_NVIDIA_GPU, DEV_INTEL_MIC, etc.}
         real(8), private:: max_flops=0d0                            !max flop/s count
         integer(INTD), private:: num_mem_units=0                    !number of attached memory units
         type(memory_access_t), allocatable, private:: mem_access(:) !description of each attached memory unit/channel: [1..num_mem_units]
        end type compute_unit_t
  !Abstract network interface:
        type, public:: nic_t
         integer(INTL), private:: id=-1   !NIC id (non-negative)
         real(8), private:: latency=0d0   !average latency of the data transfer is seconds
         real(8), private:: bandwidth=0d0 !average bandwidth of the data transfer in bytes/sec
        end type nic_t
  !Abstract compute node:
        type, public:: compute_node_t
         integer(INTL), private:: id=-1                              !compute node id (non-negative)
         integer(INTD), private:: num_mem_units=0                    !number of present memory units
         integer(INTD), private:: num_comp_units=0                   !number of present compute units
         integer(INTD), private:: num_nic_units=0                    !number of presernt NIC units
         type(memory_t), allocatable, private:: mem_unit(:)          !memory units
         type(compute_unit_t), allocatable, private:: comp_unit(:)   !compute units
         type(nic_t), allocatable, private:: nic_unit(:)             !NIC units
        end type compute_node_t
 !Hierarchical computing system representation:
        type, public:: comp_system_t
         integer(INTL), private:: num_phys_nodes=0                 !number of physical nodes in the system
         integer(INTL), private:: num_virt_nodes=0                 !number of virtual (aggregated) nodes in the system
         type(compute_node_t), allocatable, private:: virt_node(:) !virtual nodes: first <num_phys_nodes> are physical, rest are their aggregates (virtual)
         type(tree_t), private:: aggr_tree                         !node aggregation tree (NAT)
         contains
          procedure, private:: CompSystemCtorSimple
          generic, public:: comp_system_ctor=>CompSystemCtorSimple
          final:: comp_system_dtor
        end type comp_system_t
!GLOBAL:
 !Computing system:
        type(comp_system_t), public:: comp_system
!VISIBILITY:
 !comp_system_t:
        private CompSystemCtorSimple
        public comp_system_dtor

       contains
!IMPLEMENTATION:
![comp_system_t]================================================
        subroutine CompSystemCtorSimple(this,hardware_spec,ierr)
!Constructs a hierarchical (virtual) representation of a computing system
!by reading its configuration from a specification file. Simple dichotomy:
!Finds how many nodes the HPC system consists of and creates the NAT.
         implicit none
         class(comp_system_t), intent(out):: this    !out: hierarchical virtual representation of the computing system
         character(*), intent(in):: hardware_spec    !in: computing system specification file
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,l,m
         character(1024):: str

         errc=0
         open(10,file=hardware_spec(1:len_trim(hardware_spec)),form='FORMATTED',status='OLD',err=2000)
          str=' '
          do
           read(10,'(A1024)',end=100) str; l=len_trim(str)
           if(begins_with(str(1:l),'@node architecture',m,errc)) then
            if(errc.ne.0) exit

           elseif(begins_with(str(1:l),'@system architecture',m,errc)) then
            if(errc.ne.0) exit

           endif
           str(1:l)=' '
          enddo
100      close(10)
         if(present(ierr)) ierr=errc
         return
!--------------
2000     write(CONS_OUT,*)'#ERROR(ExaTENSOR::hardware): Unable to open the hardware specification file: '//&
         &hardware_spec(1:len_trim(hardware_spec)); errc=-1
         if(present(ierr)) ierr=errc
         return
        end subroutine CompSystemCtorSimple
!----------------------------------------
        subroutine comp_system_dtor(this)
!Destructor for <comp_system_t>.
         implicit none
         type(comp_system_t):: this !inout: computing system representation
         type(tree_iter_t):: tree_it
         integer(INTD):: ierr

         ierr=tree_it%init(this%aggr_tree)
         if(ierr.eq.GFC_SUCCESS) ierr=tree_it%delete_subtree()
         ierr=tree_it%release()
         if(allocated(this%virt_node)) deallocate(this%virt_node)
         this%num_virt_nodes=0
         this%num_phys_nodes=0
         return
        end subroutine comp_system_dtor

       end module hardware
