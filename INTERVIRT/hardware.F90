!Hardware representation module.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2016/03/24

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
        use tree
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6
        integer(INTD), private:: DEBUG=0
        logical, private:: VERBOSE=.true.
!TYPES:
 !Memory unit:
        type, public:: memory_unit_t
         integer(INTD), private:: id=-1   !memory unit id
         real(8), private:: capacity=0d0  !memory unit capacity
         real(8), private:: latency=0d0   !(average) access latency (when attached to a compute unit)
         real(8), private:: bandwidth=0d0 !(average) access bandwidth (when attached to a compute unit)
        end type memory_unit_t
 !Compute unit:
        type, public:: compute_unit_t
         integer(INTD), private:: id=-1                          !compute unit id
         real(8), private:: max_flops=0d0                        !max flop/s count
         integer(INTD), private:: num_mem_units=0                !number of attached memory units
         type(memory_unit_t), allocatable, private:: mem_unit(:) !description of each memory unit attached
        end type compute_unit_t
 !Compute node:
        type, public:: compute_node_t
         integer(INTD), private:: id=-1                              !node id
         integer(INTD), private:: num_mem_units=0                    !number of memory units on node
         type(memory_unit_t), allocatable, private:: mem_unit(:)     !description of each memory unit
         integer(INTD), private:: num_comp_units=0                   !number of compute units on node
         type(compute_unit_t), allocatable, private:: comp_unit(:)   !description of each compute unit with its connections to memory units
         integer(INTD), private:: num_connections=0                  !number of nodes directly connected to this one
         integer(INTD), allocatable, private:: connection(:)         !id's of direclty connected nodes
        end type compute_node_t
 !Cluster of nodes:
        type, public:: node_cluster_t
         integer(INTD), private:: num_nodes=0                 !number of nodes in the cluster
         type(compute_node_t), allocatable, private:: node(:) !description of each node
        end type node_cluster_t
!INTERFACES:

!VISIBILITY:

!DATA:
 !Computer cluster (or supercomputer):
        type(node_cluster_t), protected:: node_cluster
 !Node Aggregation Tree (NAT):
        type(tree_t), protected:: nat_global

       contains
!IMPLEMENTATION:
!----------------------------------------

       end module hardware
