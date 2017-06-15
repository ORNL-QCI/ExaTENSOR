!ExaTENSOR: Massively Parallel Virtual Processor for Scale-Adaptive Tensor Algebra
!This is the top level API module of ExaTENSOR (user-level API)
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2017/06/15

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
 !Error codes:
       public EXA_SUCCESS,&
             &EXA_ERROR,&
             &EXA_ERR_INVALID_ARGS,&
             &EXA_ERR_INVALID_REQ,&
             &EXA_ERR_MEM_ALLOC_FAIL,&
             &EXA_ERR_MEM_FREE_FAIL,&
             &EXA_ERR_BROKEN_OBJ,&
             &EXA_ERR_UNABLE_COMPLETE
!TYPES:
 !ExaTENSOR runtime status:
       type, public:: exatns_rt_status_t
        integer(INTD), public:: state=DSVP_STAT_OFF !state (see parameters above)
        integer(INTD), public:: error_code=-1       !specific error code (0 means no error)
        integer(INTD), public:: num_procs=0         !number of MPI processes involved
       end type exatns_rt_status_t
 !Vector space status:
       type, public:: exa_space_status_t
        logical, public:: built=.FALSE.               !whether or not the hierarchical (tree) structure has been imposed
        integer(INTL), public:: full_dimension=0_INTL !full dimension of the vector space (number of basis vectors)
        integer(INTD), public:: num_subspaces=0       !number of registered subspaces
       end type exa_space_status_t
 !Tensor status:
       type, public:: exa_tensor_status_t
        logical, public:: created=.FALSE. !TRUE if the tensor has been created, FALSE otherwise
        logical, public:: defined=.FALSE. !TRUE if the tensor value has been defined, FALSE otherwise
        logical, public:: in_use=.FALSE.  !TRUE if the tensor is currently participating in a computation, FALSE otherwise
       end type exa_tensor_status_t
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
 !Overloads:
      interface exatns_tensor_init
       module procedure exatns_tensor_init_scalar
       module procedure exatns_tensor_init_method
      end interface exatns_tensor_init
!DATA:
 !ExaTENSOR state:
      type(exatns_rt_status_t), protected:: exatns_rt_status
 !TAVP composition:
      integer(INTD), protected:: exa_num_workers=0  !number of worker processes
      integer(INTD), protected:: exa_num_managers=0 !number of manager processes
      integer(INTD), protected:: exa_num_helpers=0  !number of helper processes
!VISIBILITY:
 !Control:
       public exatns_start                !starts the ExaTENSOR DSVP
       public exatns_stop                 !stops the ExaTENSOR DSVP
       public exatns_status               !returns the status of the ExaTENSOR DSVP (plus statistics, if needed)
 !Parser/interpreter:
       public exatns_interpret            !interprets TAProL code (string of TAProL statements)
       public exatns_symbol_exists        !checks whether a specific identifier is registered (if yes, returns its attributes)
 !External methods/data:
       public exatns_method_register      !registers an external method (it has to adhere to a predefined interface)
       public exatns_method_unregister    !unregisters an external method
       public exatns_data_register        !registers external data (for future references)
       public exatns_data_unregister      !unregisters external data
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
       public exatns_tensor_copy          !copies the content of one tensor into another tensor, allowing for permutation, slicing, or insertion
       public exatns_tensor_fold          !produces a new tensor by folding multiple tensor dimensions into a single one
       public exatns_tensor_unfold        !produces a new tensor by unfolding a tensor dimension into multiple dimensions
       public exatns_tensor_scale         !multiplies all tensor elements by a real/complex number
       public exatns_tensor_add           !adds a tensor to another tensor, possibly with permutation, slicing, or insertion
       public exatns_tensor_contract      !contracts two tensors to produce another tensor (tensor product is included as a special case)
       public exatns_tensor_unary_op      !custom unary tensor operation via a user-defined (external) method: tensor0 = Func(tensor0,scalar)
       public exatns_tensor_binary_op     !custom binary tensor operation via a user-defined (external) method: tensor0 = Func(tensor0,tensor1,scalar)
       public exatns_tensor_ternary_op    !custom ternary tensor operation via a user-defined (external) method: tensor0 = Func(tensor0,tensor1,tensor2,scalar)
       private exatns_tensor_init_scalar  !tensor initialization by a scalar
       private exatns_tensor_init_method  !tensor initialization by a user-defined method

      contains
!IMPLEMENTATION:
![ExaTENSOR Control API]-----------------------------------
       function exatns_start(mpi_communicator) result(ierr) !called by all MPI processes
!Starts the ExaTENSOR runtime within the given MPI communicator.
!This function must be called by every MPI process from <mpi_communicator>,
!however only the Master process 0 (interpreter) will return. The other MPI
!processes will receive their active role as workers, managers, helpers, etc.
!Actions performed:
! # Initializes MPI services;
! # Determines how many TAVPs of different kinds to spawn (log-TAVP, num-TAVP, dat-TAVP, etc.);
! # Establishes a hierarchy for num-TAVPs and maps the corresponding computing domains to log-TAVPs;
! # Creates dedicated MPI communicators for log-TAVPs and for num-TAVPs;
! # Each MPI process is assigned a single TAVP of a specific kind, starts its TAVP BIOS code (init);
! # The Master MPI process returns while all other MPI processes begin their active life cycle,
!   until terminated by the Master MPI process (root-TAVP).
!This is the only ExaTENSOR API function called by all MPI processes,
!the rest are to be called by the Master MPI process only.
        implicit none
        integer(INTD):: ierr                                      !out: error code
        integer(INT_MPI), intent(in), optional:: mpi_communicator !in: MPI communicator (defaults to MPI_COMM_WORLD)
        integer(INT_MPI):: errc,num_procs,my_rank
        integer(INTD):: error_code

        ierr=EXA_SUCCESS
        if(exatns_rt_status%state.ne.DSVP_STAT_OFF) then; ierr=-1; return; endif
!Start the (MPI) process and init its distributed services:
        if(present(mpi_communicator)) then
         call dil_process_start(errc,mpi_communicator)
        else
         call dil_process_start(errc)
        endif
        if(errc.ne.0) then !failed to start an MPI process
         call dil_process_finish(errc)
         ierr=-2; return
        endif
!Sync everyone:
        num_procs=dil_global_comm_size(); my_rank=dil_global_process_id()
        if(num_procs.lt.3) then
         write(jo,'("#FATAL(exatensor): ExaTENSOR requires at least three MPI processes!")')
         call dil_process_finish(errc)
         ierr=-3; return
        endif
        write(jo,'("###EXATENSOR LAUNCHED PROCESS ",i9,"/",i9,": Syncing ... ")',ADVANCE='NO') my_rank,num_procs
        call dil_global_comm_barrier(errc)
        if(errc.eq.0) then
         write(jo,'("Ok")')
        else
         write(jo,'("Failed")')
         call dil_process_finish(errc)
         ierr=-4; return
        endif
!Determine the roles of MPI processes:
        write(jo,'("#MSG(exatensor): Determining the role of the process ... ")',ADVANCE='NO')
        call determine_process_role()
        write(jo,'(" Done: Role = ",i4,": Group rank/size = ",i9,"/",i9)') process_role,group_rank,group_size
        write(jo,'("#MSG(exatensor): Creating role specific MPI communicators ... ")',ADVANCE='NO')
        call MPI_Comm_split(GLOBAL_MPI_COMM,process_role,group_rank,group_comm,errc)
        if(errc.eq.0) then
         write(jo,'("Done")')
        else
         write(jo,'("Failed")')
         call dil_process_finish(errc)
         ierr=-5; return
        endif
!Build the Node Aggregation Tree (NAT) for compute nodes:
        if(process_role.eq.EXA_MANAGER) then
         write(jo,'("#MSG(exatensor): Building the Node Aggregation Tree (NAT) ... ")',ADVANCE='NO')
         call comp_system%comp_system_ctor('hardware.exaconf',exa_num_workers,error_code)
         if(error_code.eq.0) then
          write(jo,'("Done")')
         else
          write(jo,'("Failed")')
          call dil_process_finish(errc)
          ierr=-6; return
         endif
        endif
!Sync all MPI processes:
        call dil_global_comm_barrier(errc)
        if(errc.ne.0) then
         call dil_process_finish(errc)
         ierr=-7; return
        endif
!Activate:
        exatns_rt_status=exatns_rt_status_t(DSVP_STAT_ON,EXA_SUCCESS,num_procs)
!Life:
        ierr=EXA_SUCCESS
!Deactivate:
        exatns_rt_status=exatns_rt_status_t(DSVP_STAT_OFF,ierr,0)
!Sync everyone:
        write(jo,'()')
        write(jo,'("###EXATENSOR FINISHED PROCESS ",i9,"/",i9,": Status = ",i4,": Syncing ... ")',ADVANCE='NO')&
             &dil_global_process_id(),dil_global_comm_size(),ierr
        call dil_global_comm_barrier(errc)
        if(errc.eq.0) then; write(jo,'("Ok")'); else; write(jo,'("Failed")'); ierr=-8; endif
!Free role specific communicators:
        call MPI_Comm_free(group_comm,errc); if(errc.ne.0.and.ierr.eq.0) ierr=-9
!Finish the (MPI) process:
        call dil_process_finish(errc); if(errc.ne.0.and.ierr.eq.0) ierr=-10
        return

        contains

         subroutine determine_process_role()
          exa_num_workers=num_procs-2
          exa_num_managers=1
          if(my_rank.eq.0) then
           process_role=EXA_PARSER; group_size=1; group_rank=0
          elseif(my_rank.eq.1) then
           process_role=EXA_MANAGER; group_size=exa_num_managers; group_rank=0
          else
           process_role=EXA_WORKER; group_size=exa_num_workers; group_rank=my_rank-2
          endif
          return
         end subroutine determine_process_role

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
        integer(INTD):: ierr                        !out: error code
        type(exatns_rt_status_t), intent(out):: sts !out: current status of the ExaTENSOR runtime

        ierr=EXA_SUCCESS
        return
       end function exatns_status
![ExaTENSOR Parser/Interpreter API]------------------
       function exatns_interpret(taprol) result(ierr)
!Interprets TAProL code.
        implicit none
        integer(INTD):: ierr              !out: error code
        character(*), intent(in):: taprol !in: TAProL code (string of TAProL statements)

        ierr=EXA_SUCCESS
        return
       end function exatns_interpret
!-------------------------------------------------------
       function exatns_symbol_exists(symbol) result(ans)
!Checks whether a specific identifier is registered with ExaTENSOR.
        implicit none
        logical:: ans                     !out: answer {TRUE|FALSE}
        character(*), intent(in):: symbol !in: specific symbolic identifier

        ans=.FALSE.
        return
       end function exatns_symbol_exists
![ExaTENSOR External Method/Data API]---------------------------------------------
       function exatns_method_register(method,method_name,method_tag) result(ierr)
!Registers an external method (tensor operation) with ExaTENSOR.
        implicit none
        integer(INTD):: ierr                    !out: error code
        procedure(ext_method_i):: method        !in: external method (tensor operation)
        character(*), intent(in):: method_name  !in: symbolic method name
        integer(INTD), intent(out):: method_tag !out: method tag (non-negative on success)

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
        type(C_PTR), intent(in):: data_ptr    !in: pointer to the external data (local)
        character(*), intent(in):: data_name  !in: symbolic data name
        integer(INTD), intent(out):: data_tag !out: data tag (non-negative on success)

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
![ExaTENSOR Hierarchical Vector Space API]-----------------------------------------
       function exatns_space_register(space_name,space_basis,space_id) result(ierr)
!Registers a vector space.
        implicit none
        integer(INTD):: ierr                              !out: error code
        character(*), intent(in):: space_name             !in: vector space symbolic name
        class(subspace_basis_t), intent(in):: space_basis !in: vector space basis (fully defined)
        integer(INTD), intent(out):: space_id             !out: vector space id (non-negative)

        ierr=EXA_SUCCESS; space_id=-1
        return
       end function exatns_space_register
!---------------------------------------------------------------
       function exatns_space_unregister(space_name) result(ierr)
!Unregisters a registered vector space.
        implicit none
        integer(INTD):: ierr                  !out: error code
        character(*), intent(in):: space_name !in: vector space symbolic name

        ierr=EXA_SUCCESS
        return
       end function exatns_space_unregister
!---------------------------------------------------------------
       function exatns_space_build_tree(space_name) result(ierr)
!Builds a hierarchical representation on a registered vector space.
        implicit none
        integer(INTD):: ierr                  !out: error code
        character(*), intent(in):: space_name !in: vector space symbolic name

        ierr=EXA_SUCCESS
        return
       end function exatns_space_build_tree
!---------------------------------------------------------------
       function exatns_space_status(space_name,sts) result(ierr)
!Returns the status of a registered vector space.
        implicit none
        integer(INTD):: ierr                        !out: error code
        character(*), intent(in):: space_name       !in: vector space symbolic name
        type(exa_space_status_t), intent(out):: sts !out: vector space status

        ierr=EXA_SUCCESS
        return
       end function exatns_space_status
!------------------------------------------------------------------------------------------------------------------
       function exatns_subspace_register(space_name,subspace_name,basis_subrange,subspace_id,space_id) result(ierr)
!Registers a subspace within a registered vector space.
        implicit none
        integer(INTD):: ierr                            !out: error code
        character(*), intent(in):: space_name           !in: parental space symbolic name
        character(*), intent(in):: subspace_name        !in: subspace symbolic name
        class(seg_int_t), intent(in):: basis_subrange   !in: defining subrange of basis vectors
        integer(INTL), intent(out):: subspace_id        !out: subspace id within the parental vector space
        integer(INTD), intent(out), optional:: space_id !out: vector space id
        integer(INTD):: hspid

        ierr=EXA_SUCCESS; subspace_id=-1_INTL; hspid=-1
        if(present(space_id)) space_id=hspid
        return
       end function exatns_subspace_register
!---------------------------------------------------------------------
       function exatns_subspace_unregister(subspace_name) result(ierr)
!Unregisters a registered subspace.
        implicit none
        integer(INTD):: ierr                     !out: error code
        character(*), intent(in):: subspace_name !in: subspace symbolic name

        ierr=EXA_SUCCESS
        return
       end function exatns_subspace_unregister
!------------------------------------------------------------------------
       function exatns_index_register(space_name,index_name) result(ierr)
!Registers an index by associating it with a specific space/subspace.
        implicit none
        integer(INTD):: ierr                  !out: error code
        character(*), intent(in):: space_name !in: space/subspace symbolic name
        character(*), intent(in):: index_name !in: index symbolic name

        ierr=EXA_SUCCESS
        return
       end function exatns_index_register
!---------------------------------------------------------------
       function exatns_index_unregister(index_name) result(ierr)
!Unregisters a registers index.
        implicit none
        integer(INTD):: ierr                  !out: error code
        character(*), intent(in):: index_name !in: index symbolic name

        ierr=EXA_SUCCESS
        return
       end function exatns_index_unregister
![ExaTENSOR Tensor API]------------------------------------------------------------------------------------------------------
       function exatns_tensor_create(tensor,data_kind,tens_name,hspace,subspace,dim_extent,dim_group,group_spec) result(ierr)
!Creates a tensor.
        implicit none
        integer(INTD):: ierr                                 !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor           !out: tensor
        integer(INTD), intent(in):: data_kind                !in: data kind of tensor elements: {R4,R8,C4,C8}
        character(*), intent(in):: tens_name                 !in: symbolic tensor name
        integer(INTD), intent(in):: hspace(1:)               !in: hierarchical vector space registered id for each tensor dimension
        integer(INTL), intent(in):: subspace(1:)             !in: defining subspace registered id for each tensor dimension
        integer(INTD), intent(in), optional:: dim_extent(1:) !in: dimension extent for each tensor dimension (0 means deferred)
        integer(INTD), intent(in), optional:: dim_group(1:)  !in: symmetric group (>=0) for each tensor dimension (0 means default)
        integer(INTD), intent(in), optional:: group_spec(1:) !in: symmetric group specification for non-trivial symmetric groups (see tensor_recursive.F90)
        integer(INTD):: trank

        ierr=EXA_SUCCESS
        trank=size(hspace)
        if(size(subspace).eq.trank) then
         
        else
         ierr=EXA_ERR_INVALID_ARGS
        endif
        return
       end function exatns_tensor_create
!---------------------------------------------------------
       function exatns_tensor_destroy(tensor) result(ierr)
!Destroys a tensor.
        implicit none
        integer(INTD):: ierr                       !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor !inout: tensor

        ierr=EXA_SUCCESS
        return
       end function exatns_tensor_destroy
!---------------------------------------------------------------
       function exatns_tensor_load(tensor,filename) result(ierr)
!Loads a tensor from an external storage.
        implicit none
        integer(INTD):: ierr                       !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor !inout: tensor
        character(*), intent(in):: filename        !in: file name

        ierr=EXA_SUCCESS
        return
       end function exatns_tensor_load
!---------------------------------------------------------------
       function exatns_tensor_save(tensor,filename) result(ierr)
!Saves a tensor in an external storage.
        implicit none
        integer(INTD):: ierr                       !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor !in: tensor
        character(*), intent(in):: filename        !in: file name

        ierr=EXA_SUCCESS
        return
       end function exatns_tensor_save
!------------------------------------------------------------
       function exatns_tensor_status(tensor,sts) result(ierr)
!Returns the current status of a tensor.
        implicit none
        integer(INTD):: ierr                         !out: error code
        type(tens_rcrsv_t), intent(in):: tensor      !in: tensor
        type(exa_tensor_status_t), intent(out):: sts !out: current tensor status

        ierr=EXA_SUCCESS
        return
       end function exatns_tensor_status
![ExaTENSOR Tensor Operation API]------------------------------------
       function exatns_tensor_init_scalar(tensor,scalar) result(ierr)
!Initializes all tensor elements to a given scalar value. If <scalar> is
!absent, tensor elements will be initialized to random values.
        implicit none
        integer(INTD):: ierr                       !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor !inout: tensor
        complex(8), intent(in), optional:: scalar  !in: scalar

        ierr=EXA_SUCCESS
        return
       end function exatns_tensor_init_scalar
!--------------------------------------------------------------------
       function exatns_tensor_init_method(tensor,method) result(ierr)
!Initializes the tensor via a user-defined (registered) method.
        implicit none
        integer(INTD):: ierr                       !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor !inout: tensor
        character(*), intent(in):: method          !in: registered method name

        ierr=EXA_SUCCESS
        return
       end function exatns_tensor_init_method
!--------------------------------------------------------------------
       function exatns_tensor_copy(tensor_in,tensor_out) result(ierr)
!Copies the content of one tensor into another tensor with an option of permutation,
!and/or slicing or insertion.
        implicit none
        integer(INTD):: ierr                           !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor_in  !in: input tensor
        type(tens_rcrsv_t), intent(inout):: tensor_out !inout: output tensor

        ierr=EXA_SUCCESS
        return
       end function exatns_tensor_copy
!--------------------------------------------------------------------
       function exatns_tensor_fold(tensor_in,tensor_out) result(ierr)
!Folds two or more dimensions of the input tensor, producing an output tensor
!of a lower rank (lower order). In other words, flattening.
        implicit none
        integer(INTD):: ierr                           !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor_in  !in: input tensor
        type(tens_rcrsv_t), intent(inout):: tensor_out !inout: output tensor

        ierr=EXA_SUCCESS
        return
       end function exatns_tensor_fold
!----------------------------------------------------------------------
       function exatns_tensor_unfold(tensor_in,tensor_out) result(ierr)
!Unfolds dimensions of the input tensor, producing an output tensor
!of a higher rank (higher order).
        implicit none
        integer(INTD):: ierr                           !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor_in  !in: input tensor
        type(tens_rcrsv_t), intent(inout):: tensor_out !inout: output tensor

        ierr=EXA_SUCCESS
        return
       end function exatns_tensor_unfold
!--------------------------------------------------------------
       function exatns_tensor_scale(tensor,factor) result(ierr)
!Operation: tensor *= factor
        implicit none
        integer(INTD):: ierr                       !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor !inout: tensor
        complex(8), intent(in):: factor            !in: multiplication factor

        ierr=EXA_SUCCESS
        return
       end function exatns_tensor_scale
!---------------------------------------------------------------------
       function exatns_tensor_add(tensor0,tensor1,factor) result(ierr)
!Operation: tensor0 += tensor1 * factor
        implicit none
        integer(INTD):: ierr                        !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor0 !inout: tensor
        type(tens_rcrsv_t), intent(inout):: tensor1 !in: tensor
        complex(8), intent(in), optional:: factor   !in: scalar factor

        ierr=EXA_SUCCESS
        return
       end function exatns_tensor_add
!-------------------------------------------------------------------------------------------------------------
       function exatns_tensor_contract(tensor0,tensor1,tensor2,contr_pattern,factor,restrictions) result(ierr)
!Operation: tensor0 += tensor1 * tensor2 * factor
        implicit none
        integer(INTD):: ierr
        type(tens_rcrsv_t), intent(inout):: tensor0       !inout: tensor
        type(tens_rcrsv_t), intent(inout):: tensor1       !in: tensor
        type(tens_rcrsv_t), intent(inout):: tensor2       !in: tensor
        character(*), intent(in):: contr_pattern          !in: symbolic contraction pattern
        complex(8), intent(in), optional:: factor         !in: scalar factor
        character(*), intent(in), optional:: restrictions !in: symbolic index restrictions

        ierr=EXA_SUCCESS
        return
       end function exatns_tensor_contract
!---------------------------------------------------------------------------------
       function exatns_tensor_unary_op(tensor0,method,scalar,control) result(ierr)
!Custom operation: tensor0 = Func(tensor0,scalar)
        implicit none
        integer(INTD):: ierr                        !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor0 !inout: tensor
        character(*), intent(in):: method           !in: symbolic method name
        complex(8), intent(in), optional:: scalar   !in: scalar factor
        class(*), intent(in), optional:: control    !in: control field

        ierr=EXA_SUCCESS
        return
       end function exatns_tensor_unary_op
!------------------------------------------------------------------------------------------
       function exatns_tensor_binary_op(tensor0,tensor1,method,scalar,control) result(ierr)
!Custom operation: tensor0 = Func(tensor0,tensor1,scalar)
        implicit none
        integer(INTD):: ierr                        !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor0 !inout: tensor
        type(tens_rcrsv_t), intent(inout):: tensor1 !inout: tensor
        character(*), intent(in):: method           !in: symbolic method name
        complex(8), intent(in), optional:: scalar   !in: scalar factor
        class(*), intent(in), optional:: control    !in: control field

        ierr=EXA_SUCCESS
        return
       end function exatns_tensor_binary_op
!---------------------------------------------------------------------------------------------------
       function exatns_tensor_ternary_op(tensor0,tensor1,tensor2,method,scalar,control) result(ierr)
!Custom operation: tensor0 = Func(tensor0,tensor1,tensor2,scalar)
        implicit none
        integer(INTD):: ierr                        !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor0 !inout: tensor
        type(tens_rcrsv_t), intent(inout):: tensor1 !inout: tensor
        type(tens_rcrsv_t), intent(inout):: tensor2 !inout: tensor
        character(*), intent(in):: method           !in: symbolic method name
        complex(8), intent(in), optional:: scalar   !in: scalar factor
        class(*), intent(in), optional:: control    !in: control field

        ierr=EXA_SUCCESS
        return
       end function exatns_tensor_ternary_op

      end module exatensor
