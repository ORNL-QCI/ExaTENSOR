!Hardware abstraction module
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/03/03

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
        use gfc_vector
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
         integer(INTL), private:: num_phys_nodes=0 !number of physical nodes in the system
         integer(INTL), private:: num_virt_nodes=0 !number of virtual (simple + aggregated) nodes in the system
         type(vector_t), private:: virt_nodes      !virtual nodes: first <num_phys_nodes> are physical, rest are their aggregates (virtual)
         type(tree_t), private:: aggr_tree         !node aggregation tree (NAT)
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
![comp_system_t]=========================================================================
        subroutine CompSystemCtorSimple(this,hardware_spec,ierr,branch_fac,max_aggr_size)
!Constructs a hierarchical (virtual) representation of a computing system
!by reading its configuration from a specification file. Simple dichotomy:
!Finds how many nodes the HPC system consists of and creates the NAT by
!recursively splitting the node range into two (or more) parts.
         implicit none
         class(comp_system_t), intent(out):: this            !out: hierarchical virtual representation of the computing system
         character(*), intent(in):: hardware_spec            !in: computing system specification file
         integer(INTD), intent(out), optional:: ierr         !out: error code
         integer(INTD), intent(in), optional:: branch_fac    !in: tree branching factor (>=2)
         integer(INTD), intent(in), optional:: max_aggr_size !in: max allowed number of physical nodes that does not cause splitting (defaults to 1)
         integer(INTD):: errc,l,m,npr,brf,mas,offs(1:32),lens(1:32)
         character(1024):: str,nodarch,sysarch
         logical:: match,nodarch_found,sysarch_found
         type(tree_iter_t):: nit
         type(vector_iter_t):: vit
         type(seg_int_t), allocatable:: segs(:)
         class(seg_int_t), pointer:: rp
         class(*), pointer:: up

         errc=0
!Read the HPC system hardware specification:
         open(10,file=hardware_spec(1:len_trim(hardware_spec)),form='FORMATTED',status='OLD',err=2000)
         str=' '; nodarch=' '; sysarch=' '
         nodarch_found=.FALSE.; sysarch_found=.FALSE.
         do
          read(10,'(A1024)',end=100) str; l=len_trim(str)
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
              write(*,'(i8," physical nodes -> ")',ADVANCE='NO') this%num_phys_nodes
              exit
             else
              errc=-12; exit
             endif
            endif
           endif
           str(1:l)=' '
          endif
         enddo
100      close(10)
!Build the hierarchical virtual HPC system representation:
         this%num_virt_nodes=0_INTL
         if(errc.eq.0) then
          if(present(max_aggr_size)) then; mas=max_aggr_size; else; mas=1; endif !node aggregate splitting stops at <mas> (defaults to 1)
          if(present(branch_fac)) then; brf=branch_fac; else; brf=2; endif !tree branching factor (defaults to 2)
          if(mas.ge.1.and.brf.ge.2) then
           allocate(segs(1:brf),STAT=errc)
           if(errc.eq.0) then
 !Register physical nodes first (virt node # = phys node # in [1..max]):
            errc=vit%init(this%virt_nodes)
            if(errc.eq.GFC_SUCCESS) then
             do while(this%num_virt_nodes.lt.this%num_phys_nodes)
              call segs(1)%set(this%num_virt_nodes,this%num_virt_nodes+1_INTL,errc); if(errc.ne.0) exit
              errc=vit%append(segs(1)); if(errc.ne.GFC_SUCCESS) exit
              errc=vit%reset_back(); if(errc.ne.GFC_SUCCESS) exit
              this%num_virt_nodes=this%num_virt_nodes+1_INTL
             enddo
 !Register node aggregates (new virt nodes):
             if(errc.eq.0) then
              errc=nit%init(this%aggr_tree)
              if(errc.eq.GFC_SUCCESS) then
               if(this%num_phys_nodes.gt.1) then
                call segs(1)%set(0_INTL,this%num_phys_nodes,errc) !full range of physical nodes
                !write(*,*)'initial node range: ',segs(1)%lower_bound(),segs(1)%upper_bound() !debug
                if(errc.eq.0) then
                 errc=vit%append(segs(1))
                 if(errc.eq.GFC_SUCCESS) then
                  this%num_virt_nodes=this%num_virt_nodes+1_INTL
                  errc=vit%reset_back(); up=>vit%get_value(errc)
                  if(errc.eq.GFC_SUCCESS) then
                   errc=nit%add_leaf(up,assoc_only=.TRUE.) !root (full range)
                   if(errc.eq.GFC_SUCCESS) then
 !Recursive splitting (building the virtual node tree):
                    match=.TRUE.
                    tloop: do while(match)
                     match=.FALSE.
                     do while(errc.eq.GFC_SUCCESS)
  !Process the current tree vertex;
                      up=>nit%get_value(errc); if(errc.ne.GFC_SUCCESS) exit tloop
                      rp=>NULL(); select type(up); class is(seg_int_t); rp=>up; end select
                      if(.not.associated(rp)) then; errc=-11; exit tloop; endif
                      m=int(min(rp%length(),int(brf,INTL)),INTD)
                      if(rp%length().gt.int(mas,INTL).and.m.gt.1) then
                       call rp%split(m,segs,errc); if(errc.ne.0) exit tloop
                       do l=1,m
                        !write(*,*)'adding new subrange: ',segs(l)%lower_bound(),segs(l)%upper_bound() !debug
                        if(segs(l)%length().gt.1_INTL) then !has to be an aggregate to be added as a new virtual node
                         if(segs(l)%length().gt.mas) match=.TRUE.
                         errc=vit%append(segs(l)); if(errc.ne.GFC_SUCCESS) exit tloop
                         this%num_virt_nodes=this%num_virt_nodes+1_INTL !each node aggregate is added as a virtual node
                         errc=vit%reset_back(); up=>vit%get_value(errc); if(errc.ne.GFC_SUCCESS) exit tloop
                         errc=nit%add_leaf(up,assoc_only=.TRUE.,no_move=.TRUE.); if(errc.ne.GFC_SUCCESS) exit tloop
                        endif
                       enddo
                      endif
  !Move to the right cousin (within the current tree level):
                      errc=nit%move_to_cousin()
                     enddo
                     if(errc.eq.GFC_NO_MOVE) then; errc=GFC_SUCCESS; else; exit tloop; endif
  !Move to the children level:
                     do while(errc.eq.GFC_SUCCESS); errc=nit%move_to_cousin(to_previous=.TRUE.); enddo
                     if(errc.eq.GFC_NO_MOVE) then; errc=GFC_SUCCESS; else; exit tloop; endif
                     cloop: do
                      errc=nit%move_to_child(); if(errc.eq.GFC_SUCCESS) exit cloop
                      if(errc.eq.GFC_NO_MOVE) then
                       errc=nit%move_to_cousin(); if(errc.eq.GFC_SUCCESS) cycle cloop
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
                else
                 errc=-7
                endif
               endif
               if(errc.eq.0) then
                write(*,'(i8," virtual nodes. ")',ADVANCE='NO') this%num_virt_nodes
               else
                write(*,'(" FAILED. ")',ADVANCE='NO')
               endif
              else
               errc=-6
              endif
             else
              errc=-5
             endif
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
         type(vector_iter_t):: vec_it
         integer(INTD):: ierr

         ierr=tree_it%init(this%aggr_tree)
         if(ierr.eq.GFC_SUCCESS) ierr=tree_it%delete_subtree()
         ierr=tree_it%release()
         ierr=vec_it%init(this%virt_nodes)
         if(ierr.eq.GFC_SUCCESS) ierr=vec_it%delete_all()
         ierr=vec_it%release()
         this%num_virt_nodes=0
         this%num_phys_nodes=0
         return
        end subroutine comp_system_dtor

       end module hardware
