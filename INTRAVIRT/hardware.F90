!Hardware abstraction module
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/03/16

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
        use gfc_vec_tree
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
         integer(INTL), private:: num_phys_nodes=0 !number of physical nodes in the system, N: [1:N]
         integer(INTL), private:: num_virt_nodes=0 !number of virtual (simple + aggregated) nodes in the system, M: [0:M-1] = [0:N-1] + [N:M-1]
         type(vec_tree_t), private:: virt_nodes    !virtual nodes: first <num_phys_nodes> are physical, rest are their aggregates (virtual)
         contains
          procedure, private:: CompSystemCtorSimple                 !simple dichotomy ctor
          generic, public:: comp_system_ctor=>CompSystemCtorSimple  !ctors
          procedure, public:: print_it=>CompSystemPrintIt           !prints the node aggregation tree
          final:: comp_system_dtor                                  !dtor
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
              write(*,'(i8," physical nodes -> ")',ADVANCE='NO') this%num_phys_nodes
              if(this%num_phys_nodes.le.0) errc=-13
              exit
             else
              errc=-12; exit
             endif
            endif
           endif
           str(1:l)=' '
          endif
         enddo
         close(10)
!Build the hierarchical virtual computing system representation for <num_phys_nodes> nodes:
         this%num_virt_nodes=0_INTL
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
               !write(*,*)'initial node range: ',segs(1)%lower_bound(),segs(1)%upper_bound() !debug
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
                   do while(errc.eq.GFC_SUCCESS)
  !Process the current tree vertex;
                    up=>vt_it%get_value(errc); if(errc.ne.GFC_SUCCESS) exit tloop
                    rp=>NULL(); select type(up); class is(seg_int_t); rp=>up; end select
                    if(.not.associated(rp)) then; errc=-11; exit tloop; endif
                    if(rp%length().lt.int(mas*brf,INTL)) then
                     m=int(rp%length(),INTD)
                    else
                     m=brf
                    endif
                    if(m.gt.1) then
                     call rp%split(m,segs,errc); if(errc.ne.0) exit tloop
                     do l=1,m
                      !write(*,*)'adding new subrange: ',segs(l)%lower_bound(),segs(l)%upper_bound() !debug
                      if(segs(l)%length().gt.1_INTL) then !has to be an aggregate to be added as a new virtual node
                       errc=vt_it%append(segs(l)); if(errc.ne.GFC_SUCCESS) exit tloop
                       errc=vt_it%add_leaf(this%num_virt_nodes,no_move=.TRUE.); if(errc.ne.GFC_SUCCESS) exit tloop
                       this%num_virt_nodes=this%num_virt_nodes+1_INTL !each node aggregate is added as a virtual node
                       match=.TRUE.
                      else !individual node
                       errc=vt_it%add_leaf(segs(l)%lower_bound(),no_move=.TRUE.); if(errc.ne.GFC_SUCCESS) exit tloop
                      endif
                     enddo
                    endif
  !Move to the right cousin (within the current tree level):
                    errc=vt_it%move_to_cousin()
                   enddo
                   if(errc.eq.GFC_NO_MOVE) then; errc=GFC_SUCCESS; else; exit tloop; endif
  !Move to the children level:
                   do while(errc.eq.GFC_SUCCESS); errc=vt_it%move_to_cousin(to_previous=.TRUE.); enddo
                   if(errc.eq.GFC_NO_MOVE) then; errc=GFC_SUCCESS; else; exit tloop; endif
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
          write(*,'(i8," virtual nodes. ")',ADVANCE='NO') this%num_virt_nodes
          !call this%print_it() !debug
         else
          write(*,'(" Error ",i4,1x)',ADVANCE='NO') errc
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
!-----------------------------------------
        subroutine CompSystemPrintIt(this)
!Prints the node aggregation tree imposed on the computing system.
         implicit none
         class(comp_system_t), intent(in):: this       !in: hierarchical computing system representation
         integer(INTD):: errc
         type(vec_tree_iter_t):: vt_it

         errc=vt_it%init(this%virt_nodes)
         if(errc.eq.GFC_SUCCESS) errc=vt_it%scanp(action_f=seg_int_print_range)
         errc=vt_it%release()
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
