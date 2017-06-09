!ExaTENSOR: Massively Parallel Virtual Processor for Scale-Adaptive Tensor Algebra
!This is the top level API module of ExaTENSOR (user space API)
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/06/09

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
 !Error codes:
       public EXA_SUCCESS,EXA_ERROR,EXA_ERR_INVALID_ARGS,EXA_ERR_INVALID_REQ
       public EXA_ERR_MEM_ALLOC_FAIL,EXA_ERR_MEM_FREE_FAIL,EXA_ERR_BROKEN_OBJ
 !ExaTENSOR status:
       integer(INTD), parameter, public:: EXATNS_STAT_OFF=0  !ExaTENSOR status: OFF
       integer(INTD), parameter, public:: EXATNS_STAT_ON=1   !ExaTENSOR status: ON
       integer(INTD), parameter, public:: EXATNS_STAT_ERR=-1 !ExaTENSOR status: ERROR
 !Basic:
       integer(INTD), private:: CONS_OUT=6 !output device
       integer(INTD), private:: DEBUG=0    !debugging level
       logical, private:: VERBOSE=.TRUE.   !verbosity for errors
!TYPES:
 !ExaTENSOR status:
       type, public:: exatns_status_t
        integer(INTD):: state=EXATNS_STAT_OFF
        integer(INTD):: error_code=-1
        integer(INTD):: num_procs=0
       end type exatns_status_t
!INTERFACES:
      abstract interface
 !External method (tensor operation):
       function ext_method_i(tens_args,scal_args) result(ierr)
        import:: INTD,tens_rcrsv_t
        implicit none
        integer(INTD):: ierr                                !out: error code
        class(tens_rcrsv_t), intent(inout):: tens_args(0:)  !inout: tensor arguments
        complex(8), intent(inout), optional:: scal_args(0:) !inout: scalar arguments
       end function ext_method_i
      end interface
!DATA:

!VISIBILITY:
       public exa_tensor                  !entry point into ExaTensor (debug)
 !Control:
       public exatns_start                !starts the ExaTENSOR DSVP
       public exatns_stop                 !stops the ExaTENSOR DSVP
       public exatns_status               !returns the status of the ExaTENSOR DSVP (plus statistics, if needed)
 !Symbols:
       public exatns_symbol_exists        !checks whether a specific identifier is registered (if yes, returns its attributes)
 !External methods/data:
       public exatns_method_register      !registers an external method (it has to adhere to a predefined interface)
       public exatns_method_unregister    !unregisters an external method
       public exatns_data_register        !registers external data (for future references)
       public exatns_data_unregister      !unregisters external data
#if 0
 !Hierarchical vector space:
       public exatns_space_register       !registers a vector space
       public exatns_space_unregister     !unregisters a vector space
       public exatns_space_build_tree     !builds the subspace aggregation tree (SAT) in a vector space
       public exatns_space_status         !returns the status of the vector space
       public exatns_subspace_register    !registers a subspace in a vector space
       public exatns_subspace_unregister  !unregisters a subspace in a vector space
       public exatns_index_register       !associates an index label with a specific space/subspace
       public exatns_index_unregister     !unregisters an index label
 !Tensor:
       public exatns_tensor_create        !creates an empty tensor with an optional deferred initialization method
       public exatns_tensor_destroy       !destroys a tensor
       public exatns_tensor_load          !loads a tensor from persistent storage (create + populate)
       public exatns_tensor_save          !saves a tensor to persistent storage
       public exatns_tensor_status        !returns the status of the tensor (e.g., empty, initialized, being updated, etc.)
 !Tensor operations:
       public exatns_tensor_init          !initializes the tensor to a real/complex value or invokes an initialization method
       public exatns_tensor_copy          !copies a tensor into another tensor, possibly with permutation, slicing, or insertion
       public exatns_tensor_fold          !produces a new tensor by folding multiple tensor dimensions into a single one
       public exatns_tensor_unfold        !produces a new tensor by unfolding a tensor dimension into multiple dimensions
       public exatns_tensor_scale         !multiplies all tensor elements by a real/complex number
       public exatns_tensor_add           !adds a tensor to another tensor, possibly with permutation, slicing, or insertion
       public exatns_tensor_contract      !contracts two tensors to produce another tensor (tensor product is included as a special case)
       public exatns_tensor_operation     !custom tensor operation via user-defined (external) methods
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
       integer(INT_MPI):: errc,my_rank
       integer(INTD):: error_code

       ierr=0; errc=0
!Start the (MPI) process and init its services:
       if(present(ext_comm)) then
        call dil_process_start(errc,ext_comm)
       else
        call dil_process_start(errc)
       endif
       if(errc.ne.0) then !failed to start an MPI process
        call dil_process_finish(errc)
        ierr=-1; return
       endif
!Sync everyone:
       my_rank=dil_global_process_id()
       write(jo,'("###EXATENSOR LAUNCHED PROCESS ",i9,"/",i9,": Syncing ... ")',ADVANCE='NO')&
            &my_rank,dil_global_comm_size()
       call dil_global_comm_barrier(errc)
       if(errc.ne.0) then
        write(jo,'("Failed")')
        call dil_process_finish(errc)
        ierr=-2; return
       endif
!Root builds NAT and assigns roles to all processes:
       write(jo,'("Ok")')
       if(my_rank.eq.0) then
 !Build NAT:
        write(jo,'("#MSG(exatensor): Building the Node Aggregation Tree (NAT) ... ")',ADVANCE='NO')
        call comp_system%comp_system_ctor('hardware.exaconf',dil_global_comm_size(),error_code)
        if(error_code.ne.0) then
         write(jo,'("Failed")')
         call dil_process_finish(errc)
         ierr=-3; return
        endif
       endif
!Sync all processes:
       call dil_global_comm_barrier(errc)
       if(errc.ne.0) then
        call dil_process_finish(errc)
        ierr=-4; return
       endif
 !Assign roles:
       call tavp_establish_role(EXA_WORKER,errc) !debug
 !Begin life:
       !...(ierr)
!Sync everyone:
       write(jo,'()')
       write(jo,'("###EXATENSOR FINISHED PROCESS ",i9,"/",i9,": Status = ",i4,": Syncing ... ")',ADVANCE='NO')&
            &my_rank,dil_global_comm_size(),ierr
       call dil_global_comm_barrier(errc)
       if(errc.eq.0) then; write(jo,'("Ok")'); else; write(jo,'("Failed")'); ierr=-5; endif
!Finish the (MPI) process:
       call dil_process_finish(errc); if(errc.ne.0.and.ierr.eq.0) ierr=-6
       return
       end subroutine exa_tensor
![ExaTENSOR API]-------------------------------------------
       function exatns_start(mpi_communicator) result(ierr)
!Starts the ExaTENSOR runtime within the given MPI communicator.
!This function must be called by every MPI process from <mpi_communicator>,
!but only the master process 0 will return.
        implicit none
        integer(INTD):: ierr                                      !out: error code
        integer(INT_MPI), intent(in), optional:: mpi_communicator !in: MPI communicator (defaults to MPI_COMM_WORLD)

        ierr=EXA_SUCCESS
        return
       end function exatns_start
!-----------------------------------------
       function exatns_stop() result(ierr)
!Stops the ExaTENSOR runtime.
        implicit none
        integer(INTD):: ierr !out: error code

        ierr=EXA_SUCCESS
        return
       end function exatns_stop
!----------------------------------------------
       function exatns_status(sts) result(ierr)
!Returns the current status of the ExaTENSOR runtime.
        implicit none
        integer(INTD):: ierr                     !out: error code
        type(exatns_status_t), intent(out):: sts !out: status of the ExaTENSOR runtime

        ierr=EXA_SUCCESS
        return
       end function exatns_status
!-------------------------------------------------
       function exatns_symbol_exists() result(ans)
!Checks whether a specific identifier is registered with ExaTENSOR.
        implicit none
        logical:: ans !out: answer {TRUE|FALSE}

        ans=.FALSE.
        return
       end function exatns_symbol_exists
!---------------------------------------------------------------------------------
       function exatns_method_register(method,method_name,method_tag) result(ierr)
!Registers an external method (tensor operation) with ExaTENSOR.
        implicit none
        integer(INTD):: ierr                    !out: error code
        procedure(ext_method_i):: method        !in: external method (tensor operation)
        character(*), intent(in):: method_name  !in: symbolic method name
        integer(INTD), intent(out):: method_tag !out: method tag

        ierr=EXA_SUCCESS; method_tag=-1
        return
       end function exatns_method_register
!-----------------------------------------------------------------
       function exatns_method_unregister(method_name) result(ierr)
!Unregisters a registered external method (tensor operation).
        implicit none
        integer(INTD):: ierr                   !out: error code
        character(*), intent(in):: method_name !in: method name

        ierr=EXA_SUCCESS
        return
       end function exatns_method_unregister
!-----------------------------------------------------------------------------
       function exatns_data_register(data_ptr,data_name,data_tag) result(ierr)
!Registers an external data with ExaTENSOR.
        implicit none
        integer(INTD):: ierr                  !out: error code
        type(C_PTR), intent(in):: data_ptr    !in: pointer to the external data
        character(*), intent(in):: data_name  !in: symbolic data name
        integer(INTD), intent(out):: data_tag !out: data tag

        ierr=EXA_SUCCESS; data_tag=-1
        return
       end function exatns_data_register
!-------------------------------------------------------------
       function exatns_data_unregister(data_name) result(ierr)
!Unregisters a registered external data.
        implicit none
        integer(INTD):: ierr                 !out: error code
        character(*), intent(in):: data_name !in: data name

        ierr=EXA_SUCCESS
        return
       end function exatns_data_unregister

      end module exatensor
