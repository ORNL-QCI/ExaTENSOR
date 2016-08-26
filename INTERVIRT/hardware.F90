!Hardware abstraction module.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2016/08/26

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

       module hardware
        use dil_basic
        use gfc_base
        use tree
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6
        integer(INTD), private:: DEBUG=0
        logical, private:: VERBOSE=.true.
 !String lengths:
        integer(INTD), parameter, public:: EXA_MAX_DEVICE_DESCR_LEN=256
 !Device role:
        integer(INTD), parameter, public:: DEVICE_NO_ROLE=0
        integer(INTD), parameter, public:: DEVICE_PROCESSOR=1
        integer(INTD), parameter, public:: DEVICE_MEMORY=2
        integer(INTD), parameter, public:: DEVICE_NIC=3
        integer(INTD), parameter, public:: DEVICE_NVRAM=4
!TYPES:
 !Hardware device:
        type, public:: hardware_device_t
         integer(INTL), private:: id=-1    !global device id
         integer(INTD), private:: role     !device role (see parameters above)
         integer(INTD), private:: dev_kind !device kind (e.g., DEV_HOST, DEV_NVIDIA_GPU, etc. when .role=DEVICE_PROCESSOR)
         integer(INTD), private:: dev_num  !kind-specific device number [0..max]
         character(EXA_MAX_DEVICE_DESCR_LEN), private:: descr !readable device description
        end type hardware_device_t
 !Abstract data transfer link:
        type, public:: data_link_t
         integer(INTL), private:: id=-1   !data link id
         real(8), private:: latency=0d0   !average latency of the data link is seconds
         real(8), private:: bandwidth=0d0 !average bandwidth of the data link in bytes/sec
        end type data_link_t
 !Hardware node (linked set of devices):
        type, public:: hardware_node_t
         integer(INTD), private:: num_devices=0 !total number of devices on the node
         type(hardware_device_t), allocatable, private:: device(:) !description of each device
         type(data_link_t), allocatable, private:: link(:,:) !description of links between all devices (those which make sense)
        end type hardware_node_t
 !Abstract memory:
        type, public:: memory_t
         integer(INTL), private:: id=-1  !memory unit id
         real(8), private:: capacity=0d0 !memory unit capacity in bytes
        end type memory_t
 !Attached memory unit:
        type, public:: memory_unit_t
         type(memory_t), private:: memory     !memory description
         type(data_link_t), private:: channel !memory channel description
         logical, private:: remote            !TRUE if the memory unit is remote w.r.t. to the compute unit, FALSE otherwise
        end type memory_unit_t
 !Compute unit:
        type, public:: compute_unit_t
         integer(INTL), private:: id=-1                          !compute unit id
         real(8), private:: max_flops=0d0                        !max flop/s count
         integer(INTD), private:: num_mem_units=0                !number of attached memory units (own + remote)
         type(memory_unit_t), allocatable, private:: mem_unit(:) !description of each attached memory unit: [1..num_mem_units]
        end type compute_unit_t
!INTERFACES:
!VISIBILITY:
        public build_nat
!       public destroy_nat
!DATA:
 !Node Aggregation Tree (NAT):
        type(tree_t), protected:: exaNAT

       contains
!IMPLEMENTATION:
!------------------------------------------------
        subroutine build_nat(hardware_fname,ierr)
!Builds the Node Aggregation Tree <exaNAT> for a specific computer cluster
!specified by the hardware specification file named <hardware_fname>.
         implicit none
         character(*), intent(in):: hardware_fname !in: name of the detailed hardware specification file
         integer(INTD), intent(out):: ierr         !out: error code (0:success)
         type(tree_iter_t):: it_nat

         ierr=it_nat%init(exaNAT); if(ierr.ne.0) then; ierr=1; return; endif
         ierr=it_nat%get_status(); if(ierr.ne.GFC_IT_EMPTY) then; ierr=2; return; endif !NAT cannot be defined twice (unless destroyed)
         
         ierr=it_nat%release()
         return
        end subroutine build_nat

       end module hardware
