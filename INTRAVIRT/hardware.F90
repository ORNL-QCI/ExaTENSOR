!ExaTENSOR hardware abstraction module
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2019/01/22

!Copyright (C) 2014-2019 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2019 Oak Ridge National Laboratory (UT-Battelle)

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
!This module provides the HPC platform abstraction layer.
!An HPC platform is considered comprised of the so-called
!physical nodes which are either full hardware nodes or isolated
!parts of a full hardware node, depending on the number of
!MPI processes per hardware node, such that each physical node is
!occupied by a single MPI process. Physical nodes are numerated
!from 1: [1:N]. Physical nodes are subsequently aggregated into the
!Node Aggregation Tree (NAT). Each node of the NAT is a virtual node,
!that is, a virtual node is either a single physical node or an aggregate
!of those. Virtual nodes are numerated from 0. The first N virtual nodes [0:N-1]
!are associated with the original physical nodes [1:N] in order,
!the subsequent virtual nodes are their aggregates (inner tree nodes).
!The root of the NAT (full platform) is the virtual node N. Virtual nodes
!with numbers larger than N are the inner NAT nodes representing smaller
!aggregates of physical nodes.
        use dil_basic
        use stsubs
        use parse_prim
        use gfc_base
        use gfc_vec_tree
        use subspaces
        implicit none
        private
        public LEFT_SIBLING,RIGHT_SIBLING
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
         integer(INTL), private:: num_phys_nodes=0  !number of physical nodes in the system, N: [1:N]
         integer(INTL), private:: num_virt_nodes=0  !number of virtual (simple + aggregated) nodes in the system, M: [0:M-1] = [0:N-1] + [N:M-1]
         integer(INTD), private:: num_aggr_levels=0 !number of levels in the virtual node hierarchy (>=1), excludes leaves
         type(vec_tree_t), private:: virt_nodes     !virtual nodes: first <num_phys_nodes> are physical, rest are their aggregates (virtual)
         contains
          procedure, private:: CompSystemCtorSimple                         !simple dichotomy ctor
          generic, public:: comp_system_ctor=>CompSystemCtorSimple          !ctors
          procedure, public:: get_num_phys_nodes=>CompSystemGetNumPhysNodes !returns the total number of physical nodes in the system
          procedure, public:: get_num_virt_nodes=>CompSystemGetNumVirtNodes !returns the total number of virtual nodes in the system
          procedure, public:: get_num_aggr_nodes=>CompSystemGetNumAggrNodes !returns the total number of aggregate virtual nodes (inner nodes of the tree)
          procedure, public:: get_num_aggr_levels=>CompSystemGetNumAggrLevels!returns the number of levels in the inner virtual hierarchy (excludes leaves)
          procedure, public:: get_root_id=>CompSystemGetRootId              !returns the id of the computer system root
          procedure, public:: get_ancestor_id=>CompSystemGetAncestorId      !returns the id of an ancestor of a specific node
          procedure, public:: get_sibling_id=>CompSystemGetSiblingId        !returns the id of a left/right sibling of a specific node
          procedure, public:: get_cousin_id=>CompSystemGetCousinId          !returns the id of a left/right cousin of a specific node
          procedure, public:: get_num_children=>CompSystemGetNumChildren    !returns the number of children nodes for a specific node
          procedure, public:: get_children_ids=>CompSystemGetChildrenIds    !returns the ids of all children of a specific node in order
          procedure, public:: get_hier_level=>CompSystemGetHierLevel        !returns the tree level of a specific node (distance from the root in hops)
          procedure, public:: get_node_range=>CompSystemGetNodeRange        !returns the range of physical nodes associated with a given virtual node
          procedure, public:: print_it=>CompSystemPrintIt                   !prints the node aggregation tree
          final:: comp_system_dtor                                          !dtor
        end type comp_system_t
!GLOBAL:
 !Computing system:
        type(comp_system_t), public:: comp_system
!VISIBILITY:
 !comp_system_t:
        private CompSystemCtorSimple
        private CompSystemGetNumPhysNodes
        private CompSystemGetNumVirtNodes
        private CompSystemGetNumAggrNodes
        private CompSystemGetNumAggrLevels
        private CompSystemGetRootId
        private CompSystemGetAncestorId
        private CompSystemGetSiblingId
        private CompSystemGetCousinId
        private CompSystemGetNumChildren
        private CompSystemGetChildrenIds
        private CompSystemGetHierLevel
        private CompSystemGetNodeRange
        private CompSystemPrintIt
        public comp_system_dtor

       contains
!IMPLEMENTATION:
![comp_system_t]===================================================================================
        subroutine CompSystemCtorSimple(this,hardware_spec,num_procs,ierr,branch_fac,min_aggr_size)
!Constructs a hierarchical (virtual) representation of a computing system
!by reading its configuration from the provided specification file. Simple dichotomy:
!Finds how many nodes the HPC system consists of and creates the node aggregation tree (NAT)
!by recursively splitting the node range into two (or more) parts.
         implicit none
         class(comp_system_t), intent(out):: this            !out: hierarchical virtual representation of the computing system
         character(*), intent(in):: hardware_spec            !in: computing system specification file
         integer(INTD), intent(in):: num_procs               !in: number of MPI processes running on the computing system (must be multiple of the number of compute nodes)
         integer(INTD), intent(out), optional:: ierr         !out: error code
         integer(INTD), intent(in), optional:: branch_fac    !in: tree branching factor (>=2)
         integer(INTD), intent(in), optional:: min_aggr_size !in: min node aggregate size after which the aggregate will be fully split into individual nodes (default to 1)
         integer(INTD):: errc,l,m,npr,brf,mas,offs(1:32),lens(1:32)
         character(1024):: str,nodarch,sysarch
         logical:: match,nodarch_found,sysarch_found
         type(vec_tree_iter_t):: vt_it
         type(seg_int_t), allocatable:: segs(:)
         class(seg_int_t), pointer:: rp
         class(*), pointer:: up

         errc=0
#if 0
!Read the computing system hardware specification:
         open(10,file=hardware_spec(1:len_trim(hardware_spec)),form='FORMATTED',status='OLD',err=2000)
         str=' '; nodarch=' '; sysarch=' '
         nodarch_found=.FALSE.; sysarch_found=.FALSE.
         do
          read(10,'(A1024)',IOSTAT=errc) str; if(errc.ne.0) then; errc=max(0,errc); exit; endif
          l=len_trim(str)
          if(l.gt.0) then
           match=match_symb_pattern(str(1:l),'@node architecture [`]',npr,offs,lens,errc); if(errc.ne.0) exit
           if(match.and.npr.eq.1) then; nodarch=str(offs(1):offs(1)+lens(1)-1); nodarch_found=.TRUE.; cycle; endif
           match=match_symb_pattern(str(1:l),'@system architecture [`]',npr,offs,lens,errc); if(errc.ne.0) exit
           if(match.and.npr.eq.1) then; sysarch=str(offs(1):offs(1)+lens(1)-1); sysarch_found=.TRUE.; cycle; endif
           if(nodarch_found.and.sysarch_found) then
            match=match_symb_pattern(str(1:l),'$'//nodarch(1:len_trim(nodarch))//'*` `',npr,offs,lens,errc); if(errc.ne.0) exit
            if(match) then
             if(is_this_integer(str(offs(1):offs(1)+lens(1)-1),no_sign=.TRUE.)) then
              this%num_phys_nodes=icharnum(lens(1),str(offs(1):offs(1)+lens(1)-1))
              if(this%num_phys_nodes.le.0) errc=-21
              exit
             else
              errc=-20; exit
             endif
            endif
           endif
           str(1:l)=' '
          endif
         enddo
         close(10)
#endif
         this%num_phys_nodes=num_procs; if(this%num_phys_nodes.le.0) errc=-19
         !write(CONS_OUT,'(i8," physical nodes -> ")',ADVANCE='NO') this%num_phys_nodes
!Build the hierarchical virtual computing system representation for <num_phys_nodes> nodes:
         this%num_virt_nodes=0_INTL; this%num_aggr_levels=0
         if(errc.eq.0) then
          if(present(min_aggr_size)) then; mas=min_aggr_size; else; mas=1; endif
          if(present(branch_fac)) then; brf=branch_fac; else; brf=2; endif
          if(mas.ge.1.and.brf.ge.2) then
           allocate(segs(1:max(brf,mas)),STAT=errc)
           if(errc.eq.0) then
 !Register physical nodes first (virt node # = (phys node # - 1) in [1..num_phys_nodes]):
            errc=vt_it%init(this%virt_nodes)
            if(errc.eq.GFC_SUCCESS) then
             do while(this%num_virt_nodes.lt.this%num_phys_nodes)
              call segs(1)%set(this%num_virt_nodes,this%num_virt_nodes+1_INTL,errc); if(errc.ne.0) exit
              errc=vt_it%append(segs(1)); if(errc.ne.GFC_SUCCESS) exit
              this%num_virt_nodes=this%num_virt_nodes+1_INTL
             enddo
 !Register node aggregates (new virtual nodes):
             if(errc.eq.0) then
              if(this%num_phys_nodes.gt.1) then
               call segs(1)%set(0_INTL,this%num_phys_nodes,errc) !full range of physical nodes: (0:N] = [1:N]
               !write(CONS_OUT,*)'initial node range: ',segs(1)%lower_bound(),segs(1)%upper_bound() !debug
               if(errc.eq.0) then
                errc=vt_it%append(segs(1))
                if(errc.eq.GFC_SUCCESS) then
                 errc=vt_it%add_leaf(this%num_virt_nodes) !root (full range of nodes)
                 if(errc.eq.GFC_SUCCESS) then
                  this%num_virt_nodes=this%num_virt_nodes+1_INTL
 !Recursive decomposition (building the virtual node tree):
                  match=.TRUE.
                  tloop: do while(match)
                   match=.FALSE.
                   this%num_aggr_levels=this%num_aggr_levels+1
                   do while(errc.eq.GFC_SUCCESS)
  !Process the current tree vertex;
                    up=>vt_it%get_value(errc); if(errc.ne.GFC_SUCCESS) then; errc=-18; exit tloop; endif
                    rp=>NULL(); select type(up); class is(seg_int_t); rp=>up; end select
                    if(.not.associated(rp)) then; errc=-17; exit tloop; endif
                    m=min(int((rp%length()-1_INTL)/int(mas,INTL),INTD)+1,brf); if(m.le.1) m=int(rp%length(),INTD)
                    if(m.gt.1) then
                     call rp%split(m,segs,errc); if(errc.ne.0) then; errc=-16; exit tloop; endif
                     do l=1,m
                      !write(CONS_OUT,*)'adding new subrange: ',segs(l)%lower_bound(),segs(l)%upper_bound() !debug
                      if(segs(l)%length().gt.1_INTL) then !has to be an aggregate to be added as a new virtual node
                       errc=vt_it%append(segs(l)); if(errc.ne.GFC_SUCCESS) then; errc=-15; exit tloop; endif
                       errc=vt_it%add_leaf(this%num_virt_nodes,no_move=.TRUE.)
                       if(errc.ne.GFC_SUCCESS) then; errc=-14; exit tloop; endif
                       this%num_virt_nodes=this%num_virt_nodes+1_INTL !each node aggregate is added as a virtual node
                       match=.TRUE.
                      else !individual node
                       errc=vt_it%add_leaf(segs(l)%lower_bound(),no_move=.TRUE.)
                       if(errc.ne.GFC_SUCCESS) then; errc=-13; exit tloop; endif
                      endif
                     enddo
                    endif
  !Move to the right cousin (within the current tree level):
                    errc=vt_it%move_to_cousin()
                   enddo
                   if(errc.eq.GFC_NO_MOVE) then; errc=GFC_SUCCESS; else; errc=-12; exit tloop; endif
  !Move to the children level:
                   do while(errc.eq.GFC_SUCCESS); errc=vt_it%move_to_cousin(to_previous=.TRUE.); enddo
                   if(errc.eq.GFC_NO_MOVE) then; errc=GFC_SUCCESS; else; errc=-11; exit tloop; endif
                   cloop: do
                    errc=vt_it%move_to_child(); if(errc.eq.GFC_SUCCESS) exit cloop
                    if(errc.eq.GFC_NO_MOVE) then
                     errc=vt_it%move_to_cousin(); if(errc.eq.GFC_SUCCESS) cycle cloop
                     if(errc.eq.GFC_NO_MOVE.and.(.not.match)) errc=GFC_SUCCESS
                    endif
                    exit tloop
                   enddo cloop
                  enddo tloop
                 else
                  errc=-10
                 endif
                else
                 errc=-9
                endif
               else
                errc=-8
               endif
              else !only one physical node in the computing system (already registered as virtual node 0)
               errc=vt_it%add_leaf(0_INTL); if(errc.ne.GFC_SUCCESS) errc=-7
              endif
             else
              errc=-6
             endif
             m=vt_it%release(); if(errc.eq.0.and.m.ne.0) errc=-5
            else
             errc=-4
            endif
            if(allocated(segs)) deallocate(segs)
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         endif
         if(errc.eq.0) then
          !write(CONS_OUT,'(i8," virtual nodes. ")',ADVANCE='NO') this%num_virt_nodes
          !call this%print_it() !debug
         else
          !write(CONS_OUT,'(" Error ",i4,1x)',ADVANCE='NO') errc
          call comp_system_dtor(this)
         endif
         if(present(ierr)) ierr=errc
         return
!--------------
2000     write(CONS_OUT,*)'#ERROR(ExaTENSOR::hardware): Unable to open the hardware specification file: '//&
         &hardware_spec(1:len_trim(hardware_spec)); errc=-1
         if(present(ierr)) ierr=errc
         return
        end subroutine CompSystemCtorSimple
!----------------------------------------------------------------
        function CompSystemGetNumPhysNodes(this,ierr) result(num)
         implicit none
         integer(INTL):: num                         !out: number of physical nodes in the system
         class(comp_system_t), intent(in):: this     !in: hierarchical computer system
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         num=this%num_phys_nodes
         if(present(ierr)) ierr=errc
         return
        end function CompSystemGetNumPhysNodes
!----------------------------------------------------------------
        function CompSystemGetNumVirtNodes(this,ierr) result(num)
         implicit none
         integer(INTL):: num                         !out: number of virtual nodes in the system
         class(comp_system_t), intent(in):: this     !in: hierarchical computer system
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         num=this%num_virt_nodes
         if(present(ierr)) ierr=errc
         return
        end function CompSystemGetNumVirtNodes
!----------------------------------------------------------------
        function CompSystemGetNumAggrNodes(this,ierr) result(num)
         implicit none
         integer(INTL):: num                         !out: number of virtual aggregated nodes in the system
         class(comp_system_t), intent(in):: this     !in: hierarchical computer system
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         num=this%num_virt_nodes-this%num_phys_nodes
         if(present(ierr)) ierr=errc
         return
        end function CompSystemGetNumAggrNodes
!-----------------------------------------------------------------
        function CompSystemGetNumAggrLevels(this,ierr) result(num)
         implicit none
         integer(INTD):: num                         !out: number of aggregated levels in the inner virtual hierarchy
         class(comp_system_t), intent(in):: this     !in: hierarchical computer system
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         num=this%num_aggr_levels
         if(present(ierr)) ierr=errc
         return
        end function CompSystemGetNumAggrLevels
!---------------------------------------------------------
        function CompSystemGetRootId(this,ierr) result(id)
         implicit none
         integer(INTL):: id                          !out: id of the system root (tree root)
         class(comp_system_t), intent(in):: this     !in: hierarchical computer system
         integer(INTD), intent(out), optional:: ierr
         integer(INTD):: errc,i
         type(vec_tree_iter_t):: vtit

         id=-1_INTL; errc=vtit%init(this%virt_nodes)
         if(errc.eq.GFC_SUCCESS) then
          id=vtit%get_root_id(errc)
          i=vtit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=i
         endif
         if(present(ierr)) ierr=errc
         return
        end function CompSystemGetRootId
!---------------------------------------------------------------------------------------
        function CompSystemGetAncestorId(this,node_id,ancestor_distance,ierr) result(id)
         implicit none
         integer(INTL):: id                            !out: id of an ancestor of a specific node
         class(comp_system_t), intent(in):: this       !in: hierarchical computer system
         integer(INTL), intent(in):: node_id           !in: node id
         integer(INTD), intent(in):: ancestor_distance !in: distance to the requested ancestor (1:parent,2,3,...)
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: errc,i
         type(vec_tree_iter_t):: vtit

         id=-1_INTL
         if(ancestor_distance.gt.0) then
          errc=vtit%init(this%virt_nodes)
          if(errc.eq.GFC_SUCCESS) then
           errc=vtit%move_to(node_id)
           if(errc.eq.GFC_SUCCESS) then
            errc=vtit%move_up(ancestor_distance)
            if(errc.eq.GFC_SUCCESS) then
             id=vtit%get_offset(errc)
            else
             if(errc.eq.GFC_NO_MOVE) errc=GFC_SUCCESS
            endif
           endif
           i=vtit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=i
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function CompSystemGetAncestorId
!-------------------------------------------------------------------------
        function CompSystemGetSiblingId(this,node_id,side,ierr) result(id)
         implicit none
         integer(INTL):: id                          !out: id of the left/right sibling of a specific node
         class(comp_system_t), intent(in):: this     !in: hierarchical computer system
         integer(INTL), intent(in):: node_id         !in: node id
         logical, intent(in):: side                  !in: left/right side selector: {LEFT_SIBLING,RIGHT_SIBLING}
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,i
         type(vec_tree_iter_t):: vtit

         id=-1; errc=vtit%init(this%virt_nodes)
         if(errc.eq.GFC_SUCCESS) then
          errc=vtit%move_to(node_id)
          if(errc.eq.GFC_SUCCESS) then
           errc=vtit%move_to_sibling(to_previous=(.not.side))
           if(errc.eq.GFC_SUCCESS) then
            id=vtit%get_offset(errc)
           else
            if(errc.eq.GFC_NO_MOVE) errc=GFC_SUCCESS
           endif
          endif
          i=vtit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=i
         endif
         if(present(ierr)) ierr=errc
         return
        end function CompSystemGetSiblingId
!-----------------------------------------------------------------------------
        function CompSystemGetCousinId(this,node_id,side,ierr,ring) result(id)
!Returns the id of the left/right cousin of a specific tree node.
!The tree nodes equally distant from the tree root are cousins.
!In case the tree node <node_id> is the only one at its tree level, <id>=-1.
         implicit none
         integer(INTL):: id                          !out: id of the left/right cousin of a specific node
         class(comp_system_t), intent(in):: this     !in: hierarchical computer system
         integer(INTL), intent(in):: node_id         !in: node id
         logical, intent(in):: side                  !in: left/right side selector: {LEFT_SIBLING,RIGHT_SIBLING}
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: ring        !in: if TRUE, the ring topology will be assumed (defaults to FALSE)
         integer(INTD):: errc,i
         type(vec_tree_iter_t):: vtit

         id=-1; errc=vtit%init(this%virt_nodes)
         if(errc.eq.GFC_SUCCESS) then
          errc=vtit%move_to(node_id)
          if(errc.eq.GFC_SUCCESS) then
           errc=vtit%move_to_cousin(to_previous=(.not.side))
           if(errc.eq.GFC_SUCCESS) then
            id=vtit%get_offset(errc)
           else
            if(errc.eq.GFC_NO_MOVE) then
             errc=GFC_SUCCESS
             if(present(ring)) then
              if(ring) then
               do while(errc.eq.GFC_SUCCESS)
                errc=vtit%move_to_cousin(to_previous=side)
               enddo
               if(errc.eq.GFC_NO_MOVE) then
                id=vtit%get_offset(errc)
                if(errc.ne.GFC_SUCCESS.or.id.eq.node_id) id=-1
               endif
              endif
             endif
            endif
           endif
          endif
          i=vtit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=i
         endif
         if(present(ierr)) ierr=errc
         return
        end function CompSystemGetCousinId
!--------------------------------------------------------------------------------
        function CompSystemGetNumChildren(this,node_id,ierr) result(num_children)
         implicit none
         integer(INTL):: num_children                !out: number of children
         class(comp_system_t), intent(in):: this     !in: hierarchical computer system
         integer(INTL), intent(in):: node_id         !in: node id
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,i
         type(vec_tree_iter_t):: vtit

         num_children=0
         errc=vtit%init(this%virt_nodes)
         if(errc.eq.GFC_SUCCESS) then
          errc=vtit%move_to(node_id)
          if(errc.eq.GFC_SUCCESS) num_children=vtit%get_num_children(errc)
          i=vtit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=i
         endif
         if(present(ierr)) ierr=errc
         return
        end function CompSystemGetNumChildren
!-----------------------------------------------------------------------------------------
        function CompSystemGetChildrenIds(this,node_id,children,ierr) result(num_children)
         implicit none
         integer(INTL):: num_children                !out: number of children
         class(comp_system_t), intent(in):: this     !in: hierarchical computer system
         integer(INTL), intent(in):: node_id         !in: node id
         integer(INTL), intent(inout):: children(1:) !out: children ids
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,i
         type(vec_tree_iter_t):: vtit

         num_children=0
         errc=vtit%init(this%virt_nodes)
         if(errc.eq.GFC_SUCCESS) then
          errc=vtit%move_to(node_id)
          if(errc.eq.GFC_SUCCESS) then
           num_children=vtit%get_num_children(errc)
           if(errc.eq.GFC_SUCCESS) then
            if(size(children).ge.num_children) then
             i=0; errc=vtit%move_to_child()
             do while(errc.eq.GFC_SUCCESS)
              i=i+1; children(i)=vtit%get_offset(errc); if(errc.ne.GFC_SUCCESS) exit
              errc=vtit%move_to_sibling()
             enddo
             if(errc.eq.GFC_NO_MOVE) errc=GFC_SUCCESS
            else
             errc=GFC_OVERFLOW
            endif
           endif
          endif
          i=vtit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=i
         endif
         if(present(ierr)) ierr=errc
         return
        end function CompSystemGetChildrenIds
!-----------------------------------------------------------------------
        function CompSystemGetHierLevel(this,node_id,ierr) result(level)
         implicit none
         integer(INTD):: level                       !out: distance from the root (in hops)
         class(comp_system_t), intent(in):: this     !in: hierarchical computer system
         integer(INTL), intent(in):: node_id         !in: node id
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,i
         type(vec_tree_iter_t):: vtit

         level=-1; errc=vtit%init(this%virt_nodes)
         if(errc.eq.GFC_SUCCESS) then
          errc=vtit%move_to(node_id)
          if(errc.eq.GFC_SUCCESS) level=vtit%get_level(errc)
          i=vtit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=i
         endif
         if(present(ierr)) ierr=errc
         return
        end function CompSystemGetHierLevel
!-----------------------------------------------------------------------
        function CompSystemGetNodeRange(this,node_id,ierr) result(range)
!Returns the range of physical nodes associated with a given virtual node.
!A physical node here is essentially a computing MPI process. A virtual
!node is a node of the Node Aggregation Tree (NAT). Note that physical
!nodes are numbered from 1, but virtual nodes are numbered from 0.
         implicit none
         integer(INTD):: range(2)                    !out: integer range: [range(1):range(2)]
         class(comp_system_t), intent(in):: this     !in: hierarchical computing system
         integer(INTL), intent(in):: node_id         !in: virtual node id (>=0)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,i
         type(vec_tree_iter_t):: vtit
         class(*), pointer:: uptr
         class(seg_int_t), pointer:: seg

         errc=vtit%init(this%virt_nodes)
         if(errc.eq.GFC_SUCCESS) then
          uptr=>vtit%element_value(node_id,errc)
          if(errc.eq.GFC_SUCCESS) then
           seg=>NULL(); select type(uptr); class is(seg_int_t); seg=>uptr; end select
           if(associated(seg)) then
            range(1:2)=(/int(seg%lower_bound(),INTD)+1,int(seg%upper_bound(),INTD)/)
           else
            errc=-1
           endif
          endif
          i=vtit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=i
         endif
         if(present(ierr)) ierr=errc
         return
        end function CompSystemGetNodeRange
!------------------------------------------------------
        subroutine CompSystemPrintIt(this,ierr,dev_out)
!Prints the node aggregation tree imposed on the computing system.
         implicit none
         class(comp_system_t), intent(in):: this       !in: hierarchical computing system representation
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD), intent(in), optional:: dev_out !in: output device id (defaults to 6)
         integer(INTD):: errc,devo,i,level,range(1:2)
         integer(INTL):: node_id
         type(vec_tree_iter_t):: vt_it

         errc=0
         devo=6; if(present(dev_out)) devo=dev_out
         write(devo,'("Printing the HPC platform Node Aggregation Tree (NAT): Phys/Virt nodes = ",i6,1x,i8)')&
         &this%num_phys_nodes,this%num_virt_nodes
         errc=vt_it%init(this%virt_nodes)
         if(errc.eq.GFC_SUCCESS) then
          level=-1
          tloop: do
           level=level+1
           errc=vt_it%find_first_of_level(level); if(errc.ne.GFC_SUCCESS) exit tloop
           write(devo,'(" NAT level ",i4,":")') level
           do while(errc.eq.GFC_SUCCESS)
            node_id=vt_it%get_offset(errc); if(errc.ne.GFC_SUCCESS) exit tloop
            range(1:2)=this%get_node_range(node_id,errc); if(errc.ne.0) then; errc=-3; exit tloop; endif
            write(devo,'("  Virt Node ",i8,": Phys node range [",i6,":",i6,"]")') node_id,range(1:2)
            errc=vt_it%move_to_cousin()
           enddo
           if(errc.eq.GFC_NO_MOVE) errc=GFC_SUCCESS
          enddo tloop
          if(errc.eq.GFC_NO_MOVE) errc=GFC_SUCCESS
          i=vt_it%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=-2
         else
          errc=-1
         endif
         write(devo,'("Done: Status ",i9)') errc
         flush(devo)
         if(present(ierr)) ierr=errc
         return
        end subroutine CompSystemPrintIt
!----------------------------------------
        subroutine comp_system_dtor(this)
!Destructor for <comp_system_t>.
         implicit none
         type(comp_system_t):: this !inout: computing system representation
         type(vec_tree_iter_t):: vt_it
         integer(INTD):: ierr

         ierr=vt_it%init(this%virt_nodes)
         if(ierr.eq.GFC_SUCCESS) ierr=vt_it%delete_all()
         ierr=vt_it%release()
         this%num_virt_nodes=0
         this%num_phys_nodes=0
         return
        end subroutine comp_system_dtor

       end module hardware
