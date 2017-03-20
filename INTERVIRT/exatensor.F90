!ExaTENSOR: Massively Parallel Virtual Processor for Scale-Adaptive Tensor Algebra
!This is the top level API module of ExaTENSOR
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/03/20

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
       logical, private:: VERBOSE=.TRUE.   !verbosity for errors
!TYPES:

!INTERFACES:

!DATA:

!VISIBILITY:
 !Control API:
       public exa_tensor                    !entry point into ExaTensor (debug)
#if 0
       public exatns_start                  !starts the ExaTENSOR DSVP
       public exatns_stop                   !stops the ExaTENSOR DSVP
 !Hierarchical vector space:
       public exatns_space_register         !registers a vector space
       public exatns_space_unregister       !unregisters a vector space
       public exatns_subspace_register      !registers a subspace in a vector space
       public exatns_subspace_unregister    !unregisters a subspace in a vector space
       public exatns_index_register         !assigns an index label to a specific space/subspace
       public exatns_index_unregister       !unregisters an index label
#endif
      contains
!IMPLEMENTATION:
!------------------------------------------
       subroutine exa_tensor(ierr,ext_comm) !debug
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
       integer(INTD):: error_code
       integer(INTL), parameter:: TEST_SPACE_DIM=20    !debug
       type(spher_symmetry_t):: symm(1:TEST_SPACE_DIM) !debug
       type(subspace_basis_t):: full_basis             !debug
       type(h_space_t):: vec_space                     !debug

!Start the (MPI) process and init its services:
       ierr=0; errc=0
       if(present(ext_comm)) then
        call dil_process_start(errc,ext_comm)
       else
        call dil_process_start(errc)
       endif
       if(errc.ne.0) then; ierr=-1; return; endif !failed to start an MPI process
!Sync everyone:
       my_rank=dil_global_process_id()
       write(jo,'("###EXATENSOR LAUNCHED PROCESS ",i9,"/",i9,": Syncing ... ")',ADVANCE='NO')&
            &my_rank,dil_global_comm_size()
       call dil_global_comm_barrier(errc)
!Root builds NAT, SAT, and assigns roles to all processes:
       if(errc.eq.0) then
        write(jo,'("Ok")')
        if(my_rank.eq.0) then
 !Build NAT:
         write(jo,'("#MSG(exatensor): Building the Node Aggregation Tree (NAT) ... ")',ADVANCE='NO')
         call comp_system%comp_system_ctor('hardware.exaconf',dil_global_comm_size(),error_code)
 !Build SAT (debug):
         if(error_code.eq.0) then
          write(jo,'("Ok")')
          write(jo,'("#MSG(exatensor): Building the Subspace Aggregation Tree (SAT) ... ")',ADVANCE='NO')
          call register_test_space(error_code)
          if(error_code.eq.0) then; write(jo,'("Ok")'); else; write(jo,'("Failed")'); ierr=-2; endif
         else
          write(jo,'("Failed")')
          ierr=-3
         endif
        endif
       else
        write(jo,'("Failed")')
        ierr=-4
       endif
       call dil_global_comm_barrier(errc)
       if(ierr.eq.0.and.errc.eq.0) then
 !Assign roles:
        call tavp_establish_role(EXA_WORKER,errc) !debug
 !Begin life:
        !...(ierr)
       endif
!Sync everyone:
       write(jo,'()')
       write(jo,'("###EXATENSOR FINISHED PROCESS ",i9,"/",i9,": Status = ",i4,": Syncing ... ")',ADVANCE='NO')&
            &my_rank,dil_global_comm_size(),ierr
       call dil_global_comm_barrier(errc)
       if(errc.eq.0) then; write(jo,'("Ok")'); else; write(jo,'("Failed")'); endif
!Finish the (MPI) process:
       call dil_process_finish(errc); if(errc.ne.0) ierr=-5
       return

       contains

        subroutine register_test_space(jerr)
         implicit none
         integer(INTD), intent(out):: jerr
         integer(INTL):: jj

         jerr=0
         call full_basis%subspace_basis_ctor(TEST_SPACE_DIM,jerr); if(jerr.ne.0) return
         do jj=1,TEST_SPACE_DIM
          call symm(jj)%spher_symmetry_ctor(int((jj-1)/5,INTD),0,jerr); if(jerr.ne.0) return
          call full_basis%set_basis_func(jj,BASIS_ABSTRACT,jerr,symm=symm(jj)); if(jerr.ne.0) return
         enddo
         call full_basis%finalize(jerr); if(jerr.ne.0) return
         call vec_space%h_space_ctor(full_basis,jerr); if(jerr.ne.0) return
         write(*,'(i9,"-dimensional full space -> ",i9," subspaces ")',ADVANCE='NO')&
              &vec_space%get_space_dim(),vec_space%get_num_subspaces() !debug
         !call vec_space%print_it() !debug
         return
        end subroutine register_test_space

       end subroutine exa_tensor

      end module exatensor
