!ExaTensor: Massively Parallel Virtual Processor for Scale-Adaptive Tensor Algebra
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

      module exatensor
       use virta
!      use c_process
       use m_process
       use g_process
       implicit none
       private
!PARAMETERS:
 !Basic:
       integer(INTD), private:: CONS_OUT=6 !output device
       integer(INTD), private:: DEBUG=0    !debugging level
       logical, private:: VERBOSE=.true.   !verbosity
!TYPES:

!INTERFACES:

!DATA:

!VISIBILITY:
       public exa_tensor !entry point to ExaTensor

      contains
!IMPLEMENTATION:
!------------------------------------------
       subroutine exa_tensor(ierr,ext_comm)
!ExaTensor entry point: Starts a process, assigns a role to it, lives it, ends the process.
!If an existing MPI communicator <ext_comm> is passed here, it will be used. Otherwise,
!the MPI_COMM_WORLD will be initialized, subsequently used, and finalized at the end.
       use service_mpi, only: dil_process_start,dil_process_finish,&
                             &dil_global_comm_size,dil_global_process_id,&
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
       my_rank=dil_global_process_id()
       write(jo,'("###EXATENSOR LAUNCHED PROCESS ",i9,"/",i9,": Syncing ... ")',ADVANCE='NO')&
            &my_rank,dil_global_comm_size()
       call dil_global_comm_barrier(errc)
       if(errc.eq.0) then; write(jo,'("Ok")'); else; write(jo,'("Failed")'); endif
!Assign a role to the (MPI) process and live it:
       if(errc.eq.0) then
        my_role=EXA_WORKER !debug
!       call c_proc_life(ierr)
       else
        ierr=-1
       endif
!Sync everyone:
       write(jo,'()')
       write(jo,'("###EXATENSOR FINISHED PROCESS ",i9,"/",i9,": Status = ",i4,": Syncing ... ")',ADVANCE='NO')&
            &my_rank,dil_global_comm_size(),ierr
       call dil_global_comm_barrier(errc)
       if(errc.eq.0) then; write(jo,'("Ok")'); else; write(jo,'("Failed")'); endif
!Finish the (MPI) process:
       call dil_process_finish(errc); if(errc.ne.0) ierr=2
       return
       end subroutine exa_tensor

      end module exatensor
