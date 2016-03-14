!ExaTensor::Massively Parallel Virtual Processor for Scale-Adaptive Tensor Algebra
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2015/09/29

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
!-------------------------------------------------------------------------
      module exatensor
       use virta
       use c_process
!      use m_process
!      use g_process
       implicit none
       private
!PARAMETERS:
 !General:
       integer, private:: CONS_OUT=6     !output device
       logical, private:: VERBOSE=.true. !verbosity
       integer, private:: DEBUG=0        !debugging level
!TYPES:

!DATA:

!VISIBILITY:
       public exa_tensor !entrance to ExaTensor

      contains
!PROCEDURES:
!------------------------------------------
       subroutine exa_tensor(ierr,ext_comm)
!ExaTensor entry gate: Starts a process, assigns a role to it, lives it, ends the process.
!If an existing MPI communicator <ext_comm> is passed here, it will be used. Otherwise,
!the MPI_COMM_WORLD will be initialized.
!Roles of MPI processes:
! * Role "Managing Process": Managing Processes (MP) form a tree with a single Level-0 Manager
!   as the Global Root (GR). Other Managers are indexed separately on each tree level.
!   Managers at the last tree level are called Local Roots (LR). Each Local Root
!   is assigned a group of Computing Processes (CP) and one or more Data Processes (DP).
! * Role "Computing Process": Computing Processes (CP) run virtual machines that perform
!   the actual computations and report to their respective Local Roots. Computing Processes
!   are grouped into Computing Groups (CG), each Computing Group being assigned a unique Local Root.
! * Role "Data Process": Data Processes (DP) facilitate data traffic and disk operations.
!This subroutine also creates local MPI communicators for computing groups, managers, etc.
       use service_mpi, only: dil_process_start,dil_process_finish,dil_global_comm_size,dil_process_global_id,&
                             &dil_global_comm_barrier,INT_MPI,jo
       implicit none
       integer, intent(out):: ierr                       !out: error code (0:success)
       integer(INT_MPI), intent(in), optional:: ext_comm !in: existing MPI communicator (defaults to MPI_COMM_WORLD)
       integer(INT_MPI):: errc,i,my_rank

!Start the (MPI) process and init its services:
       ierr=0; errc=0
       if(present(ext_comm)) then
        call dil_process_start(errc,ext_comm)
       else
        call dil_process_start(errc)
       endif
       if(errc.ne.0) then; ierr=1; return; endif
!Sync everyone:
       my_rank=dil_process_global_id()
       write(jo,'("###EXATENSOR LAUNCHED PROCESS ",i9,"/",i9,": Syncing ... ")',ADVANCE='NO')&
            &my_rank,dil_global_comm_size()
       call dil_global_comm_barrier(errc)
       if(errc.eq.0) then; write(jo,'("Ok")'); else; write(jo,'("Failed")'); endif
!Assign a role to the (MPI) process and live it:
       if(errc.eq.0) then
 !DEBUG begin:
        my_role=C_PROCESS_PRIVATE
        call c_proc_life(ierr)
        return
 !DEBUG end.
        if(my_rank.eq.0) then !MPI master --> Global Root
         my_role=GLOBAL_ROOT; my_group=-1
!        call global_root_life(ierr)
        else !MPI slave --> Local Root OR Computing Process (C-process)
         i=my_rank-1; my_group=i/(1+C_PROCS_PER_LOCAL_ROOT); my_group_index=i-my_group*(1+C_PROCS_PER_LOCAL_ROOT)
         if(my_group_index.eq.0) then
          my_role=LOCAL_ROOT
!         call local_root_life(ierr)
         else
          my_role=C_PROCESS_PRIVATE
          call c_proc_life(ierr)
         endif
        endif
       else
        ierr=-1
       endif
!Sync everyone:
       write(jo,'("###EXATENSOR FINISHED PROCESS ",i9,"/",i9,": Status = ",i4,": Syncing ... ")',ADVANCE='NO')&
            &my_rank,dil_global_comm_size(),ierr
       call dil_global_comm_barrier(errc)
       if(errc.eq.0) then; write(jo,'("Ok")'); else; write(jo,'("Failed")'); endif
!Finish the (MPI) process:
       call dil_process_finish(errc); if(errc.ne.0) ierr=2
       return
       end subroutine exa_tensor

      end module exatensor
