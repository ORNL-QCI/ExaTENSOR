!ExaTENSOR: Massively Parallel Virtual Processor for Scale-Adaptive Hierarchical Tensor Algebra
!This is the top level API module of ExaTENSOR (user-level API)
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2018/09/21

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
!ExaTENSOR implements a domain-specific virtual processor specialized
!for numerical tensor algebra workloads. The entire HPC system is
!recursively virtualized as a hierarchical virtual tensor processor.
!This hierarchical virtual processor is then mapped onto the original
!set of MPI processes where each MPI process receieves a specific role:
! MPI processes [0..W-1]: TAVP-WRK/TAVP-HLP: Numerical workload;
! MPI processes [W..W+M-1]: TAVP-MNG: Task creation and scheduling;
! MPI process (W+M) (last MPI process): Driving process (interpreter).
!Thus, there are W data processors, M metadata processors, and 1 driver.
       use tavp_manager, tens_instr_mng_t=>tens_instr_t
       use tavp_worker, tens_instr_wrk_t=>tens_instr_t
       use virta
       implicit none
       private
       public EXA_NO_ROLE,EXA_DRIVER,EXA_MANAGER,EXA_WORKER,EXA_HELPER !process roles
       public EXA_DATA_KIND_NN,EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8 !tensor data kinds
       public MAX_TENSOR_RANK !symbols from tensor_algebra
       public INTD,INTL,INT_MPI,jo,impir,impis !symbols from service_mpi
!PARAMETERS:
 !Basic:
       integer(INTD), private:: CONS_OUT=6 !output device
       integer(INTD), private:: DEBUG=1    !debugging level
       logical, private:: VERBOSE=.TRUE.   !verbosity for errors
 !Error codes (user-level):
       public EXA_SUCCESS,&
             &EXA_ERROR,&
             &EXA_ERR_INVALID_ARGS,&
             &EXA_ERR_INVALID_REQ,&
             &EXA_ERR_MEM_ALLOC_FAIL,&
             &EXA_ERR_MEM_FREE_FAIL,&
             &EXA_ERR_BROKEN_OBJ,&
             &EXA_ERR_UNABLE_COMPLETE
 !Subspaces:
       public seg_int_t,orthotope_t,symmetry_t,color_symmetry_t,spher_symmetry_t,&
             &basis_func_supp_t,basis_func_gauss_t,basis_func_t,&
             &subspace_basis_t,subspace_t,h_space_t,h_index_t,&
             &BASIS_ABSTRACT
 !Tensors:
       public hspace_reg_t,tens_signature_t,tens_shape_t,tens_header_t,&
             &tens_simple_part_t,tens_layout_t,tens_layout_fdims_t,&
             &tens_body_t,tens_rcrsv_t,permutation_t,contr_ptrn_ext_t,&
             &tens_status_t,tens_method_uni_t
       public talsh_tens_signature_t,talsh_tens_shape_t,talsh_tens_data_t
!TYPES:
 !ExaTENSOR runtime status:
       type, public:: exatns_rt_status_t
        integer(INTD), public:: state=DSVP_STAT_OFF !state (see dsvp.F90)
        integer(INTD), public:: error_code=0        !specific error code (0 means no error)
        integer(INTD), public:: num_procs=0         !number of MPI processes involved
        integer(INTL), public:: num_instr=0_INTL    !number of instructions accepted from the Driver process
       end type exatns_rt_status_t
 !Vector space status:
       type, public:: exatns_space_status_t
        logical, public:: hierarchical=.FALSE.        !whether or not the hierarchical (tree) structure has been imposed
        integer(INTL), public:: full_dimension=0_INTL !full dimension of the vector space (number of basis vectors)
        integer(INTL), public:: num_subspaces=0_INTL  !number of registered subspaces
       end type exatns_space_status_t
!INTERFACES:
 !Overloads:
  !Create tensor:
       interface exatns_tensor_create
        module procedure exatns_tensor_create_scalar
        module procedure exatns_tensor_create_tensor
       end interface exatns_tensor_create
  !Initialize tensor:
       interface exatns_tensor_init
        module procedure exatns_tensor_init_scalar
        module procedure exatns_tensor_init_method
       end interface exatns_tensor_init
!GLOBAL DATA:
 !ExaTENSOR runtime status:
       type(exatns_rt_status_t), private:: exatns_rt_status
 !TAVP global composition (set by exatns_start):
       integer(INTD), protected:: exa_num_drivers=0  !number of the driver processes (1)
       integer(INTD), protected:: exa_num_managers=0 !number of the manager processes
       integer(INTD), protected:: exa_num_workers=0  !number of the worker processes
       integer(INTD), protected:: exa_num_helpers=0  !number of the helper processes
 !TAVP instance (allocated one per MPI process):
       class(dsvp_t), allocatable, private:: tavp
 !Tensor instruction log (recorded instructions by Driver):
       type(vector_t), private:: instructions
       type(vector_iter_t), private:: instr_log
       integer(INTL), protected:: num_tens_instr_issued=0 !number of tensor instructions issued by the Driver (excludes control and auxiliary instructions)
       integer(INTL), protected:: num_tens_instr_synced=0 !number of tensor instructions synchronized on completion (excludes control and auxiliary instructions)
 !Reusable bytecode buffers (for Driver):
       type(pack_env_t), private:: bytecode_out !outgoing bytecode buffer
       type(pack_env_t), private:: bytecode_in  !incoming bytecode buffer
!VISIBILITY:
 !External methods/data (called by All before exatns_start()):
       public exatns_dim_resolution_setup    !sets up the universal tensor dimension extent resolution function which determines the actual shape of tensor blocks
       public exatns_dim_strength_setup      !sets up the universal tensor dimension strength assessing function and threshold (guides recursive tensor dimension splitting)
       public exatns_dim_strength_thresh_set !sets the tensor dimension strength threshold above which the dimension will split (guides recursive tensor dimension splitting)
       public exatns_method_register         !registers an external method (it has to adhere to a predefined interface)
       public exatns_method_unregister       !unregisters an external method
       public exatns_data_register           !registers external (on-node) data (for future references)
       public exatns_data_unregister         !unregisters external data
 !Control:
       public exatns_start                !starts the ExaTENSOR DSVP (called by All)
       public exatns_stop                 !stops the ExaTENSOR DSVP (Driver only)
       public exatns_sync                 !synchronizes the ExaTENSOR DSVP such that all previously issued tensor instructions will be completed (Driver only)
       public exatns_process_role         !returns the role of the current MPI process (called by Any)
       public exatns_status               !returns the status of the ExaTENSOR runtime plus statistics, if needed (Driver only)
       public exatns_dump_cache           !dumps the tensor cache content for each TAVP into its log file (debug only)
 !Parser/interpreter (Driver only):
       public exatns_interpret            !interprets TAProL code (string of TAProL statements)
       public exatns_symbol_exists        !checks whether a specific identifier is registered (if yes, returns its attributes)
 !Hierarchical vector space (either called by All before start or then by Driver only):
       public exatns_space_register       !registers a vector space
       public exatns_space_unregister     !unregisters a vector space
       public exatns_space_status         !returns the status of the vector space
       public exatns_subspace_register    !registers a subspace in a vector space
       public exatns_subspace_unregister  !unregisters a subspace in a vector space
       public exatns_index_register       !associates an index label with a specific space/subspace
       public exatns_index_unregister     !unregisters an index label
 !Tensor (Driver only):
       public exatns_tensor_create        !creates an empty tensor with an optional deferred initialization method
       public exatns_tensor_destroy       !destroys a tensor
       public exatns_tensor_get           !returns a locally storable slice of a tensor
       public exatns_tensor_load          !loads a tensor from persistent storage (create + populate)
       public exatns_tensor_save          !saves a tensor to persistent storage
       public exatns_tensor_status        !returns the status of the tensor (e.g., empty, initialized, being updated, etc.)
       private exatns_tensor_create_scalar!creates an empty scalar (order-0 tensor)
       private exatns_tensor_create_tensor!creates an empty tensor (order-1 and higher tensors)
 !Tensor operations (Driver only):
       public exatns_tensor_init          !initializes a tensor to a real/complex value or invokes an external initialization method
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
 !Internal:
       private tavp_role_rank
       private add_new_instruction
       private issue_new_instruction

      contains
!IMPLEMENTATION:
![ExaTENSOR External Method/Data API]-------------------------------------
       function exatns_dim_resolution_setup(dim_resolution_f) result(ierr) !called by all MPI processes
!Sets up the universal tensor dimension resolution function which determines the actual shape of tensor blocks.
        implicit none
        integer(INTD):: ierr                                   !out: error code
        procedure(tens_rcrsv_dim_resolve_i):: dim_resolution_f !in: tensor dimension extent resolution function

        ierr=EXA_SUCCESS
        tens_dim_extent_resolve=>dim_resolution_f
        return
       end function exatns_dim_resolution_setup
!-------------------------------------------------------------------------------------
       function exatns_dim_strength_setup(dim_strength_f,strength_thresh) result(ierr) !called by all MPI processes
!Sets up the universal tensor dimension strength assessing function.
        implicit none
        integer(INTD):: ierr                                  !out: error code
        procedure(tens_rcrsv_dim_strength_i):: dim_strength_f !in: external universal tensor dimension strength assessing function (see interface in tensor_recursive.F90)
        real(8), intent(in), optional:: strength_thresh       !in: tensor dimension strength threshold above which the dimension will split

        ierr=EXA_SUCCESS
        tens_dim_strength_assess=>dim_strength_f
        if(present(strength_thresh)) tens_dim_strength_thresh=strength_thresh
        return
       end function exatns_dim_strength_setup
!---------------------------------------------------------------------------
       function exatns_dim_strength_thresh_set(strength_thresh) result(ierr) !called by all MPI processes
!Sets the tensor dimension strength threshold above which the dimension will split,
!thus allowing for a fine tensor decomposition granularity control.
        implicit none
        integer(INTD):: ierr                  !out: error code
        real(8), intent(in):: strength_thresh !in: tensor dimension strength threshold above which the dimension will split

        ierr=EXA_SUCCESS
        tens_dim_strength_thresh=strength_thresh
        return
       end function exatns_dim_strength_thresh_set
!---------------------------------------------------------------------------------
       function exatns_method_register(method_name,method,method_tag) result(ierr) !called by all MPI processes
!Registers an external tensor body initialization/update method with ExaTENSOR.
        implicit none
        integer(INTD):: ierr                              !out: error code
        character(*), intent(in):: method_name            !in: symbolic method name
        class(tens_method_uni_t), intent(in):: method     !in: external tensor body initialization/transformation method
        integer(INTD), intent(out), optional:: method_tag !out: method tag
        integer(INTD):: tag

        tag=-1
        call method_register%register_method(method_name,method,ierr)
        if(present(method_tag)) method_tag=tag
        return
       end function exatns_method_register
!-----------------------------------------------------------------
       function exatns_method_unregister(method_name) result(ierr) !called by all MPI processes
!Unregisters a previously registered external tensor body initialization/update method.
        implicit none
        integer(INTD):: ierr                   !out: error code
        character(*), intent(in):: method_name !in: method name

        call method_register%unregister_method(method_name,ierr)
        return
       end function exatns_method_unregister
!------------------------------------------------------------------------------
       function exatns_data_register(data_name,tens_data,data_tag) result(ierr) !called by all MPI processes
!Registers an external local data with ExaTENSOR.
        implicit none
        integer(INTD):: ierr                            !out: error code
        character(*), intent(in):: data_name            !in: symbolic data name
        type(talsh_tens_data_t), intent(in):: tens_data !in: external tensor data (loca)
        integer(INTD), intent(out), optional:: data_tag !out: data tag
        integer(INTD):: tag

        tag=-1
        call data_register%register_data(data_name,tens_data,ierr)
        if(present(data_tag)) data_tag=tag
        return
       end function exatns_data_register
!-------------------------------------------------------------
       function exatns_data_unregister(data_name) result(ierr) !called by all MPI processes
!Unregisters previously registered external data.
        implicit none
        integer(INTD):: ierr                 !out: error code
        character(*), intent(in):: data_name !in: data name

        call data_register%unregister_data(data_name,ierr)
        return
       end function exatns_data_unregister
![ExaTENSOR Control API]-----------------------------------
       function exatns_start(mpi_communicator) result(ierr) !called by all MPI processes
!Starts the ExaTENSOR runtime within the given MPI communicator.
!This function must be called by every MPI process from <mpi_communicator>,
!however only one MPI process (Driver) will return. The other MPI processes
!will receive their active TAVP roles as managers, workers, helpers, etc.
!Actions performed:
! # Initializes MPI services;
! # Determines how many TAVPs of different kinds to spawn (TAVP-MNG, TAVP-WRK, TAVP-HLP, etc.);
! # Establishes a global TAVP hierarchy and maps compute node domains to TAVP-MNG recursively;
! # Creates dedicated MPI communicators for TAVP-MNG and for TAVP-WRK as well as an intercommunicator;
! # Each MPI process is assigned a single TAVP of a specific kind which is allocated and launched;
! # One MPI process (Driver) returns while all other MPI processes begin their active
!   life cycle as TAVPs until terminated by the Driver MPI process via a call to exatns_stop().
! # The Driver MPI process schedules tensor operations and finally calls exatns_stop().
        implicit none
        integer(INTD):: ierr                                      !out: error code
        integer(INT_MPI), intent(in), optional:: mpi_communicator !in: MPI communicator (defaults to MPI_COMM_WORLD)
        integer(INT_MPI):: errc,num_procs,my_rank

        ierr=EXA_SUCCESS
!Check whether the ExaTENSOR runtime is currently OFF:
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
!Build the MPI process hierarchy and determine the roles of MPI processes:
        write(jo,'("#MSG(exatensor): Building the process hierarchy and determining roles:")')
        call determine_process_role(errc) !builds NAT and determines process roles
        if(errc.eq.0) then
         write(jo,'("#MSG(exatensor): Info: Role = ",i2,": Role rank/size = ",i9,"/",i9)') process_role,role_rank,role_size
         write(jo,'("#MSG(exatensor): Info: Global rank of the driver process = ",i9)') driver_gl_rank !this MPI process will return
         write(jo,'("#MSG(exatensor): Info: Global rank of the top manager process = ",i9)') top_manager_gl_rank
         write(jo,'("#MSG(exatensor): Creating role specific MPI communicators ... ")',ADVANCE='NO')
         call MPI_Comm_split(GLOBAL_MPI_COMM,process_role,role_rank,role_comm,errc)
         if(errc.eq.0) then
          write(jo,'("Done: Comm = ",i11)') role_comm
          write(jo,'("#MSG(exatensor): Creating MPI intercommunicators ... ")',ADVANCE='NO')
          select case(process_role)
          case(EXA_DRIVER)
           call MPI_Intercomm_create(role_comm,0,GLOBAL_MPI_COMM,exa_num_workers,12,drv_mng_comm,errc)
          case(EXA_MANAGER)
           call MPI_Intercomm_create(role_comm,0,GLOBAL_MPI_COMM,num_procs-1,12,drv_mng_comm,errc)
           if(errc.eq.0) call MPI_Intercomm_create(role_comm,0,GLOBAL_MPI_COMM,0,23,mng_wrk_comm,errc)
          case(EXA_WORKER)
           call MPI_Intercomm_create(role_comm,0,GLOBAL_MPI_COMM,exa_num_workers,23,mng_wrk_comm,errc)
          case default
           call quit(-1,'#FATAL(exatensor): Unexpected process role during intercommunicator creation!')
          end select
          if(errc.eq.0) then
           write(jo,'("Done: (Drv,Mng) = ",i11,", (Mng,Wrk) = ",i11)') drv_mng_comm,mng_wrk_comm
          else
           write(jo,'("Failed: Error ",i11)') errc
           call dil_process_finish(errc)
           ierr=-5; return
          endif
         else
          write(jo,'("Failed: Error ",i11)') errc
          call dil_process_finish(errc)
          ierr=-6; return
         endif
        else
         write(jo,'("Failed: Error ",i11)') errc
         call dil_process_finish(errc)
         ierr=-7; return
        endif
!Set the default universal tensor dimension strength assessing and shape resolution functions, if none preset by a user earlier:
        if(.not.associated(tens_dim_extent_resolve)) errc=exatns_dim_resolution_setup(tens_rcrsv_dim_resolve_default)
        if(.not.associated(tens_dim_strength_assess)) errc=exatns_dim_strength_setup(tens_rcrsv_dim_strength_default,0d0)
!Sync all MPI processes before configuring and launching TAVPs:
        call dil_global_comm_barrier(errc); if(errc.ne.0) then; call dil_process_finish(errc); ierr=-8; return; endif
!Mark the ExaTENSOR runtime active:
        exatns_rt_status=exatns_rt_status_t(DSVP_STAT_ON,EXA_SUCCESS,num_procs,0_INTL)
!Live TAVP life:
        ierr=EXA_SUCCESS
        if(process_role.eq.EXA_DRIVER) then
         ierr=instr_log%init(instructions); if(ierr.ne.GFC_SUCCESS) ierr=-9
         call bytecode_out%reserve_mem(ierr); if(ierr.ne.0) ierr=-10
         call bytecode_in%reserve_mem(ierr); if(ierr.ne.0) ierr=-11
         return !Driver process returns immediately, it will later call exatns_stop()
        elseif(process_role.eq.EXA_MANAGER) then
         call tavp_mng_reset_output(jo)
         write(jo,'("#MSG(exatensor): Preparing the TAVP-MNG virtual processor ... ")',ADVANCE='NO')
         call prepare_tavp_mng(ierr)
         if(ierr.eq.0) then; write(jo,'("Done")'); else; write(jo,'("Failed: Error ",i11)') ierr; endif
        elseif(process_role.eq.EXA_WORKER) then
         call tavp_wrk_reset_output(jo)
         write(jo,'("#MSG(exatensor): Preparing the TAVP-WRK virtual processor ... ")',ADVANCE='NO')
         call prepare_tavp_wrk(ierr)
         if(ierr.eq.0) then; write(jo,'("Done")'); else; write(jo,'("Failed: Error ",i11)') ierr; endif
        elseif(process_role.eq.EXA_HELPER) then
         !call tavp_hlp_reset_output(jo)
         write(jo,'("#MSG(exatensor): Preparing the TAVP-HLP virtual processor ... ")',ADVANCE='NO')
         !call prepare_tavp_hlp(ierr)
         call quit(-1,'#FATAL(exatensor): TAVP-HLP is not implemented yet!')
         if(ierr.eq.0) then; write(jo,'("Done")'); else; write(jo,'("Failed: Error ",i11)') ierr; endif
        endif
        if(ierr.eq.0) call tavp%start(ierr) !will later call .shutdown()
!Mark the ExaTENSOR runtime is off:
        exatns_rt_status=exatns_rt_status_t(DSVP_STAT_OFF,ierr,0,exatns_rt_status%num_instr)
!Destroy TAVP:
        if(allocated(tavp)) then
         call tavp%destroy(errc); deallocate(tavp)
        endif
!Sync everyone:
        write(jo,'()')
        write(jo,'("###EXATENSOR FINISHED PROCESS ",i9,"/",i9,": Status = ",i11,": Syncing ... ")',ADVANCE='NO')&
             &dil_global_process_id(),dil_global_comm_size(),ierr
        call dil_global_comm_barrier(errc)
        if(errc.eq.0) then; write(jo,'("Ok")'); else; write(jo,'("Failed")'); ierr=-12; endif
!Free role specific MPI communicators:
        if(drv_mng_comm.ne.MPI_COMM_NULL) then
         call MPI_Comm_free(drv_mng_comm,errc); if(errc.ne.0.and.ierr.eq.0) ierr=-13
        endif
        if(mng_wrk_comm.ne.MPI_COMM_NULL) then
         call MPI_Comm_free(mng_wrk_comm,errc); if(errc.ne.0.and.ierr.eq.0) ierr=-14
        endif
        if(role_comm.ne.MPI_COMM_NULL) then
         call MPI_Comm_free(role_comm,errc); if(errc.ne.0.and.ierr.eq.0) ierr=-15
        endif
!Finish the MPI process:
        call dil_process_finish(errc); if(errc.ne.0.and.ierr.eq.0) ierr=-16
        return

        contains

         subroutine prepare_tavp_wrk(jerr)
          implicit none
          integer(INTD), intent(out):: jerr
          type(tavp_wrk_conf_t):: tavp_wrk_conf
          integer(INTD):: ji
          integer(INTL):: aid
          character(128):: tavpname

          allocate(tavp_wrk_t::tavp,STAT=jerr)
          if(jerr.eq.0) then
           tavpname='TAVP-WRK#'; call numchar(role_rank,ji,tavpname(len_trim(tavpname)+1:))
           allocate(tavp_wrk_conf%description,SOURCE=tavpname(1:len_trim(tavpname)),STAT=jerr)
           if(jerr.eq.0) then
            tavp_wrk_conf%tavp_id=my_rank !global MPI rank
            aid=comp_system%get_ancestor_id(int(my_rank,INTL),1,jerr) !parent TAVP
            if(jerr.eq.0.and.aid.ge.0) then
             tavp_wrk_conf%source_comm=mng_wrk_comm
             tavp_wrk_conf%source_rank=tavp_role_rank(int(aid,INTD))
             tavp_wrk_conf%retire_comm=tavp_wrk_conf%source_comm
             tavp_wrk_conf%retire_rank=tavp_wrk_conf%source_rank
             tavp_wrk_conf%host_ram_size=1_INTL*(1024_INTL*1024_INTL*1024_INTL) !`Make configurable
             tavp_wrk_conf%nvram_size=0_INTL !`Make configurable
             tavp_wrk_conf%num_mpi_windows=1 !`Make configurable
             tavp_wrk_conf%host_buf_size=tavp_wrk_conf%host_ram_size !`Make configurable
             if(gpu_count.gt.0) then
              allocate(tavp_wrk_conf%gpu_list(gpu_count))
              tavp_wrk_conf%gpu_list(1:gpu_count)=(/(ji,ji=gpu_start,gpu_start+gpu_count-1)/)
             endif
             call tavp%configure(tavp_wrk_conf,jerr); if(jerr.ne.0) jerr=-4
            else
             jerr=-3
            endif
           else
            jerr=-2
           endif
          else
           jerr=-1
          endif
          return
         end subroutine prepare_tavp_wrk

         subroutine prepare_tavp_mng(jerr)
          implicit none
          integer(INTD), intent(out):: jerr
          type(tavp_mng_conf_t):: tavp_mng_conf
          integer(INTD):: ji,jrl,jl,jr
          integer(INTL):: aid,lid,rid,nch
          integer(INTL), allocatable:: chid(:)
          character(128):: tavpname

          allocate(tavp_mng_t::tavp,STAT=jerr)
          if(jerr.eq.0) then
           tavpname='TAVP-MNG#'; call numchar(role_rank,ji,tavpname(len_trim(tavpname)+1:))
           allocate(tavp_mng_conf%description,SOURCE=tavpname(1:len_trim(tavpname)),STAT=jerr)
           if(jerr.eq.0) then
            tavp_mng_conf%tavp_id=my_rank !global MPI rank
            aid=comp_system%get_ancestor_id(int(my_rank,INTL),1,jerr) !parent TAVP
            if(jerr.eq.0) then
             if(aid.lt.0) then !root manager (no parent)
              tavp_mng_conf%source_comm=drv_mng_comm
              tavp_mng_conf%source_rank=0 !assumes a single Driver process
             else !intermediate manager
              tavp_mng_conf%source_comm=role_comm
              tavp_mng_conf%source_rank=tavp_role_rank(int(aid,INTD))
             endif
             tavp_mng_conf%retire_comm=tavp_mng_conf%source_comm
             tavp_mng_conf%retire_rank=tavp_mng_conf%source_rank
             lid=comp_system%get_cousin_id(int(my_rank,INTL),LEFT_SIBLING,jerr,ring=.TRUE.)
             if(jerr.eq.0) then
              if(lid.lt.0) lid=int(my_rank,INTL) !self-reference (root node)
              rid=comp_system%get_cousin_id(int(my_rank,INTL),RIGHT_SIBLING,jerr,ring=.TRUE.)
              if(jerr.eq.0) then
               if(rid.lt.0) rid=int(my_rank,INTL) !self-reference (root node)
               tavp_mng_conf%ring_comm=role_comm
               tavp_mng_conf%ring_send_rank=tavp_role_rank(int(rid,INTD),jr) !self-reference for the root manager
               tavp_mng_conf%ring_recv_rank=tavp_role_rank(int(lid,INTD),jl) !self-reference for the root manager
               if(jl.eq.process_role.and.jr.eq.process_role) then !tree nodes on the same level must be of the same kind
                nch=comp_system%get_num_children(int(my_rank,INTL),jerr)
                if(jerr.eq.0.and.nch.gt.0) then
                 if(allocated(tavp_mng_conf%dispatch_rank)) deallocate(tavp_mng_conf%dispatch_rank)
                 allocate(tavp_mng_conf%dispatch_rank(1:nch)); allocate(chid(1:nch))
                 nch=comp_system%get_children_ids(int(my_rank,INTL),chid,jerr)
                 if(jerr.eq.0) then
                  ji=tavp_role_rank(int(chid(1),INTD),jrl)
                  if(jrl.eq.EXA_MANAGER) then
                   tavp_mng_conf%dispatch_comm=role_comm
                  elseif(jrl.eq.EXA_WORKER) then
                   tavp_mng_conf%dispatch_comm=mng_wrk_comm
                  else
                   jerr=-10
                  endif
                  if(jerr.eq.0) then
                   do ji=1,int(nch,INTD)
                    tavp_mng_conf%dispatch_rank(ji)=tavp_role_rank(int(chid(ji),INTD))
                   enddo
                   tavp_mng_conf%collect_comm=tavp_mng_conf%dispatch_comm
                   call tavp%configure(tavp_mng_conf,jerr); if(jerr.ne.0) jerr=-9
                  endif
                 else
                  jerr=-8
                 endif
                 deallocate(chid)
                else
                 jerr=-7
                endif
               else
                write(jo,'("#FATAL(exatensor): Unbalanced Node Aggregation Trees are not supported yet! ")',ADVANCE='NO')
                write(jo,'("All nodes at the same tree level must be of the same kind! Adjust the number of MPI processes!")')
                jerr=-6
               endif
              else
               jerr=-5
              endif
             else
              jerr=-4
             endif
            else
             jerr=-3
            endif
           else
            jerr=-2
           endif
          else
           jerr=-1
          endif
          return
         end subroutine prepare_tavp_mng

         subroutine determine_process_role(jerr) !all MPI processes
          implicit none
          integer(INTD), intent(out):: jerr
          integer(INTD), parameter:: MAX_NAT_CONF_CYCLES=256
          integer(INTD):: num_wrk_hlp,jc

          jerr=0
          write(jo,'("#MSG(exatensor): Building the Node Aggregation Tree (NAT) ... ")',ADVANCE='NO')
 !Set the initial composition:
          exa_num_drivers=1 !the Driver process is the last MPI process
          exa_num_helpers=0
          exa_num_managers=1 !initial number of managers
          exa_num_workers=num_procs-(exa_num_drivers+exa_num_managers+exa_num_helpers)
          num_wrk_hlp=num_procs-(exa_num_drivers+exa_num_managers) !workers + helpers
 !Build NAT:
          jc=MAX_NAT_CONF_CYCLES
          hloop: do
           if(jc.le.0) then
            if(VERBOSE) then
             write(jo,'("#ERROR(exatns_start)[",i5,"]: Unable to build a proper Node Aggregation Tree!")') my_rank
             write(jo,'("#ADVICE: Try changing the number of MPI processes.")')
            endif
            jerr=-3; exit hloop
           endif
           !write(jo,'("#DEBUG(exatns_start:determine_process_role)[",i5,"]: Trial conf:",3(1x,i6))')&
                !&my_rank,exa_num_drivers,exa_num_managers,exa_num_workers !debug
  !Build the process hierarchy:
           call comp_system%comp_system_ctor('hardware.exaconf',num_wrk_hlp,jerr,&
                                            &EXA_MANAGER_BRANCH_FACT,EXA_MAX_WORK_GROUP_SIZE)
           if(jerr.ne.0) then; jerr=-2; exit hloop; endif
  !Check the total number of managers:
           if(comp_system%get_num_aggr_nodes().eq.exa_num_managers) exit hloop !match
           exa_num_managers=exa_num_managers+1
           exa_num_workers=num_procs-(exa_num_drivers+exa_num_managers+exa_num_helpers)
           num_wrk_hlp=num_procs-(exa_num_drivers+exa_num_managers) !workers + helpers
           jc=jc-1
          enddo hloop
          if(jerr.eq.0) then
           role_rank=tavp_role_rank(int(my_rank,INTD),process_role)
           driver_gl_rank=num_procs-1
           top_manager_gl_rank=exa_num_workers
           select case(process_role)
           case(EXA_DRIVER)
            role_size=exa_num_drivers
           case(EXA_MANAGER)
            role_size=exa_num_managers
           case(EXA_WORKER)
            role_size=exa_num_workers
            !`Some Workers can now be converted into Helpers
           case default
            jerr=-1
           end select
           if(jerr.eq.0) write(jo,'("Done")')
          endif
          !write(*,*) 'Process ',my_rank,': ',exa_num_drivers,exa_num_managers,exa_num_workers !debug
          return
         end subroutine determine_process_role

       end function exatns_start
!-----------------------------------------
       function exatns_stop() result(ierr)
!Stops the ExaTENSOR runtime.
        implicit none
        integer(INTD):: ierr !out: error code
        integer(INT_MPI):: errc
        integer(INTL):: ip
        class(tens_instr_mng_t), pointer:: tens_instr

        ierr=EXA_SUCCESS
        write(jo,'("#MSG(exatensor): New Instruction: STOP ExaTENSOR: IP = ")',ADVANCE='NO'); flush(jo)
!Send the STOP instruction to the root TAVP-MNG and wait for completion:
        tens_instr=>add_new_instruction(ip,ierr)
        if(ierr.eq.0) then
         write(jo,'(i11)') ip; flush(jo) !new instruction id number
         call tens_instr%tens_instr_ctor(TAVP_INSTR_CTRL_STOP,ierr,iid=ip)
         if(ierr.eq.0) then
          call issue_new_instruction(tens_instr,ierr); if(ierr.ne.0) ierr=-14
          call tens_instr%set_status(DS_INSTR_RETIRED,errc); if(errc.ne.DSVP_SUCCESS.and.ierr.eq.0) ierr=-13
          errc=exatns_sync(); if(errc.ne.EXA_SUCCESS.and.ierr.eq.0) ierr=-12
          errc=instr_log%delete_all(); if(errc.ne.GFC_SUCCESS.and.ierr.eq.0) ierr=-11
          errc=instr_log%release(); if(errc.ne.GFC_SUCCESS.and.ierr.eq.0) ierr=-10
          call bytecode_in%destroy(errc); if(errc.ne.0.and.ierr.eq.0) ierr=-9
          call bytecode_out%destroy(errc); if(errc.ne.0.and.ierr.eq.0) ierr=-8
!Mark the ExaTENSOR runtime is off:
          exatns_rt_status=exatns_rt_status_t(DSVP_STAT_OFF,ierr,0,ip+1_INTL)
!Sync with others globally:
          write(jo,'()')
          write(jo,'("###EXATENSOR FINISHED PROCESS ",i9,"/",i9,": Status = ",i11,": Syncing ... ")',ADVANCE='NO')&
               &dil_global_process_id(),dil_global_comm_size(),ierr
          flush(jo)
          call dil_global_comm_barrier(errc)
          if(errc.eq.0) then; write(jo,'("Ok")'); else; write(jo,'("Failed")'); if(ierr.eq.0) ierr=-7; endif; flush(jo)
!Free the role-specific MPI communicators:
          if(drv_mng_comm.ne.MPI_COMM_NULL) then
           call MPI_Comm_free(drv_mng_comm,errc); if(errc.ne.0.and.ierr.eq.0) ierr=-6
          endif
          if(mng_wrk_comm.ne.MPI_COMM_NULL) then
           call MPI_Comm_free(mng_wrk_comm,errc); if(errc.ne.0.and.ierr.eq.0) ierr=-5
          endif
          if(role_comm.ne.MPI_COMM_NULL) then
           call MPI_Comm_free(role_comm,errc); if(errc.ne.0.and.ierr.eq.0) ierr=-4
          endif
!Finish the MPI process:
          call dil_process_finish(errc); if(errc.ne.0.and.ierr.eq.0) ierr=-3
         else
          ierr=-2
         endif
         tens_instr=>NULL()
        else
         ierr=-1
        endif
        if(ierr.ne.0) write(jo,'(" Failed!")')
        return
       end function exatns_stop
!-----------------------------------------
       function exatns_sync() result(ierr) !Driver only
!Synchronizes the ExaTENSOR DSVP such that all previously issued tensor instructions will be completed.
        implicit none
        integer(INTD):: ierr !out: error code
        type(comm_handle_t):: comm_hl
        type(obj_pack_t):: instr_packet
        integer(INTD):: n,i,sts,err_code
        integer(INTL):: iid
        class(*), pointer:: instr
        logical:: new

        ierr=EXA_SUCCESS; call comm_hl%clean(ierr)
        if(ierr.eq.0) then
         wloop: do while(num_tens_instr_synced.lt.num_tens_instr_issued)
          new=bytecode_in%receive(comm_hl,ierr,0,TAVP_COLLECT_TAG,drv_mng_comm) !receive bytecode from the root TAVP-MNG
          if(new) then
           call comm_hl%wait(ierr); if(ierr.ne.0) then; ierr=EXA_ERR_UNABLE_COMPLETE; exit wloop; endif
           n=bytecode_in%get_num_packets()
           num_tens_instr_synced=num_tens_instr_synced+n
           do i=1,n
            call bytecode_in%extract_packet(i,instr_packet,ierr,preclean=.TRUE.)
            if(ierr.eq.PACK_SUCCESS) then
             call unpack_builtin(instr_packet,iid,ierr)
             if(ierr.eq.PACK_SUCCESS) then
              instr=>instr_log%element_value(iid,ierr)
              if(ierr.eq.GFC_SUCCESS.and.associated(instr)) then
               select type(instr)
               class is(tens_instr_mng_t)
                sts=instr%get_status(ierr,err_code)
                if(ierr.eq.DSVP_SUCCESS) then
                 call instr%set_status(DS_INSTR_RETIRED,ierr,err_code)
                 if(ierr.ne.DSVP_SUCCESS) then; ierr=EXA_ERR_UNABLE_COMPLETE; exit wloop; endif
                else
                 ierr=EXA_ERR_UNABLE_COMPLETE; exit wloop
                endif
               class default
                ierr=EXA_ERR_UNABLE_COMPLETE; exit wloop
               end select
              else
               ierr=EXA_ERR_UNABLE_COMPLETE; exit wloop
              endif
             else
              ierr=EXA_ERR_UNABLE_COMPLETE; exit wloop
             endif
            else
             ierr=EXA_ERR_UNABLE_COMPLETE; exit wloop
            endif
           enddo
          endif
          call comm_hl%clean(ierr); if(ierr.ne.0) then; ierr=EXA_ERR_UNABLE_COMPLETE; exit wloop; endif
         enddo wloop
         call bytecode_in%clean(ierr); if(ierr.ne.0) ierr=EXA_ERR_MEM_FREE_FAIL
        else
         ierr=EXA_ERR_UNABLE_COMPLETE
        endif
        return
       end function exatns_sync
!----------------------------------------------------------------
       function exatns_process_role(role,role_total) result(ierr)
!Returns the role of the current MPI process.
        implicit none
        integer(INTD):: ierr                              !out: error code
        integer(INTD), intent(out):: role                 !out: process role
        integer(INTD), intent(out), optional:: role_total !out: total number of MPI processes of this role

        ierr=EXA_SUCCESS; role=process_role
        if(present(role_total)) then
         if(role.ne.EXA_NO_ROLE) then
          role_total=role_size
         else
          role_total=0
         endif
        endif
        return
       end function exatns_process_role
!----------------------------------------------
       function exatns_status(sts) result(ierr)
!Returns the current status of the ExaTENSOR runtime.
        implicit none
        integer(INTD):: ierr                        !out: error code
        type(exatns_rt_status_t), intent(out):: sts !out: current status of the ExaTENSOR runtime

        ierr=EXA_SUCCESS; sts=exatns_rt_status
        return
       end function exatns_status
!-----------------------------------------------
       function exatns_dump_cache() result(ierr) !DEBUG only
!Dumps the tensor cache content for each TAVP into its log file.
        implicit none
        integer(INTD):: ierr !out: error code
        integer(INTD):: errc
        integer(INTL):: ip
        class(tens_instr_mng_t), pointer:: tens_instr

        ierr=EXA_SUCCESS
        write(jo,'("#MSG(exatensor): New Instruction: DUMP TENSOR CACHE: IP = ")',ADVANCE='NO'); flush(jo)
        tens_instr=>add_new_instruction(ip,ierr)
        if(ierr.eq.0) then
         write(jo,'(i11)') ip; flush(jo) !new instruction id number
         call tens_instr%tens_instr_ctor(TAVP_INSTR_CTRL_DUMP_CACHE,ierr,iid=ip)
         if(ierr.eq.0) then
          call issue_new_instruction(tens_instr,ierr); if(ierr.ne.0) ierr=-4
          call tens_instr%set_status(DS_INSTR_RETIRED,errc); if(errc.ne.DSVP_SUCCESS.and.ierr.eq.0) ierr=-3
         else
          ierr=-2
         endif
         tens_instr=>NULL()
        else
         ierr=-1
        endif
        return
       end function exatns_dump_cache
![ExaTENSOR Parser/Interpreter API]------------------
       function exatns_interpret(taprol) result(ierr)
!Interprets symbolic TAProL code.
        implicit none
        integer(INTD):: ierr              !out: error code
        character(*), intent(in):: taprol !in: TAProL code (string of valid TAProL statements)

        ierr=EXA_SUCCESS
        call quit(-1,'FATAL(exatensor:interpret): Not implemented yet!') !`Implement
        return
       end function exatns_interpret
!-------------------------------------------------------
       function exatns_symbol_exists(symbol) result(ans)
!Checks whether a specific symbolic identifier is registered with ExaTENSOR.
        implicit none
        logical:: ans                     !out: answer {TRUE|FALSE}
        character(*), intent(in):: symbol !in: specific symbolic identifier

        ans=.FALSE.
        call quit(-1,'FATAL(exatensor:symbol_exists): Not implemented yet!') !`Implement
        return
       end function exatns_symbol_exists
![ExaTENSOR Hierarchical Vector Space API]------------------------------------------------
       function exatns_space_register(space_name,space_basis,space_id,hspace) result(ierr)
!Registers a vector space based on the provided space basis.
        implicit none
        integer(INTD):: ierr                                      !out: error code
        character(*), intent(in):: space_name                     !in: vector space symbolic name
        class(subspace_basis_t), intent(in), target:: space_basis !in: vector space basis (fully defined)
        integer(INTD), intent(out):: space_id                     !out: vector space id (non-negative)
        class(h_space_t), pointer, intent(out):: hspace           !out: pointer to the registered vector space

        ierr=EXA_SUCCESS; space_id=-1
        if(space_basis%dimsn().gt.0.and.len(space_name).gt.0) then
         space_id=hspace_register%register_space(space_name,ierr,hspace)
         if(ierr.eq.TEREC_SUCCESS.and.associated(hspace)) then
          call hspace%h_space_ctor(space_basis,ierr)
          if(ierr.eq.0) then
           !if(DEBUG.gt.0) then
            !write(jo,'("#MSG(exatensor): Registered new vector space {id = ",i4,"; dim = ",i9,"}:")',ADVANCE='NO')&
                 !&space_id,hspace%get_space_dim()
            !write(jo,*) space_name
           !endif
          else
           space_id=-1; ierr=EXA_ERR_UNABLE_COMPLETE
          endif
         else
          space_id=-1; ierr=EXA_ERR_UNABLE_COMPLETE
         endif
        else
         ierr=EXA_ERR_INVALID_ARGS
        endif
        return
       end function exatns_space_register
!---------------------------------------------------------------
       function exatns_space_unregister(space_name) result(ierr)
!Unregisters a registered vector space.
        implicit none
        integer(INTD):: ierr                  !out: error code
        character(*), intent(in):: space_name !in: vector space symbolic name

        ierr=EXA_SUCCESS
        call quit(-1,'FATAL(exatensor:space_unregister): Not implemented yet!') !`Implement: Requires .unregister() method in tensor_recursive.F90
        return
       end function exatns_space_unregister
!---------------------------------------------------------------
       function exatns_space_status(space_name,sts) result(ierr)
!Returns the status of a registered vector space.
        implicit none
        integer(INTD):: ierr                           !out: error code
        character(*), intent(in):: space_name          !in: vector space symbolic name
        type(exatns_space_status_t), intent(out):: sts !out: vector space status

        ierr=EXA_SUCCESS
        call quit(-1,'FATAL(exatensor:space_status): Not implemented yet!') !`Implement
        return
       end function exatns_space_status
!------------------------------------------------------------------------------------------------------------------
       function exatns_subspace_register(subspace_name,space_name,basis_subrange,subspace_id,space_id) result(ierr)
!Registers a subspace over a contiguous subrange of basis vectors within a registered vector space.
!If the parental vector space is hierarchical, the newly created subspace will be embedded into the
!closest aligned subspace containing it (so-called embedding subspace). The corresponding subtensors
!will be defined over the embedding subspaces, using zero padding to emulate the requested subspaces.
        implicit none
        integer(INTD):: ierr                            !out: error code
        character(*), intent(in):: subspace_name        !in: subspace symbolic name
        character(*), intent(in):: space_name           !in: parental space symbolic name
        class(seg_int_t), intent(in):: basis_subrange   !in: defining subrange of basis vectors in the parental space
        integer(INTL), intent(out):: subspace_id        !out: subspace id within the parental vector space
        integer(INTD), intent(out), optional:: space_id !out: parental vector space id
        integer(INTD):: hspid

        ierr=EXA_SUCCESS; subspace_id=-1_INTL; hspid=-1
        call quit(-1,'FATAL(exatensor:subspace_register): Not implemented yet!') !`Implement
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
        call quit(-1,'FATAL(exatensor:subspace_unregister): Not implemented yet!') !`Implement
        return
       end function exatns_subspace_unregister
!------------------------------------------------------------------------
       function exatns_index_register(index_name,space_name) result(ierr)
!Registers an index by associating it with a specific space/subspace.
        implicit none
        integer(INTD):: ierr                  !out: error code
        character(*), intent(in):: index_name !in: index symbolic name
        character(*), intent(in):: space_name !in: space/subspace symbolic name

        ierr=EXA_SUCCESS
        call quit(-1,'FATAL(exatensor:index_register): Not implemented yet!') !`Implement
        return
       end function exatns_index_register
!---------------------------------------------------------------
       function exatns_index_unregister(index_name) result(ierr)
!Unregisters a registered index.
        implicit none
        integer(INTD):: ierr                  !out: error code
        character(*), intent(in):: index_name !in: index symbolic name

        ierr=EXA_SUCCESS
        call quit(-1,'FATAL(exatensor:index_unregister): Not implemented yet!') !`Implement
        return
       end function exatns_index_unregister
![ExaTENSOR Tensor API]-------------------------------------------------------------
       function exatns_tensor_create_scalar(tensor,tens_name,data_kind) result(ierr)
!Creates a scalar (order-0 tensor).
        implicit none
        integer(INTD):: ierr                       !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor !out: tensor (must be empty on entrance)
        character(*), intent(in):: tens_name       !in: symbolic tensor name
        integer(INTD), intent(in):: data_kind      !in: data kind of tensor elements: {EXA_DATA_KIND_XX: XX = R4,R8,C4,C8}
        integer(INTD):: hspaces(1)
        integer(INTL):: subspaces(1)

        ierr=exatns_tensor_create_tensor(tensor,tens_name,hspaces(1:0),subspaces(1:0),data_kind)
        return
       end function exatns_tensor_create_scalar
!-------------------------------------------------------------------------------------------------------------------------
       function exatns_tensor_create_tensor(tensor,tens_name,hspaces,subspaces,data_kind,dim_extent,dim_group,group_spec)&
                                           &result(ierr)
!Creates a tensor.
        implicit none
        integer(INTD):: ierr                                 !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor           !out: tensor (must be empty on entrance)
        character(*), intent(in):: tens_name                 !in: symbolic tensor name
        integer(INTD), intent(in):: hspaces(1:)              !in: registered id of the defining vector space for each tensor dimension
        integer(INTL), intent(in):: subspaces(1:)            !in: id of the defining subspace for each tensor dimension
        integer(INTD), intent(in):: data_kind                !in: data kind of tensor elements: {EXA_DATA_KIND_XX: XX = R4,R8,C4,C8}
        integer(INTL), intent(in), optional:: dim_extent(1:) !in: dimension extent for each tensor dimension (0 means deferred, to be set later)
        integer(INTD), intent(in), optional:: dim_group(1:)  !in: symmetric group (>=0) for each tensor dimension (0 means default)
        integer(INTD), intent(in), optional:: group_spec(1:) !in: symmetric group specification for non-trivial symmetric groups (see tensor_recursive.F90)
        class(tens_instr_mng_t), pointer:: tens_instr
        integer(INTD):: trank
        integer(INTL):: ip

        ierr=EXA_SUCCESS
        write(jo,'("#MSG(exatensor): New Instruction: CREATE TENSOR: IP = ")',ADVANCE='NO'); flush(jo)
        if(.not.tensor%is_set(ierr)) then
         if(ierr.eq.TEREC_SUCCESS) then
!Construct the tensor object:
          trank=size(hspaces) !tensor rank (number of tensor dimensions)
          if(size(subspaces).eq.trank) then
           if(present(dim_extent)) then
            if(size(dim_extent).eq.trank) then
             if(present(dim_group)) then
              if(size(dim_group).eq.trank) then
               if(present(group_spec)) then
                call tensor%tens_rcrsv_ctor(tens_name,subspaces,hspaces,ierr,dim_extent,dim_group,group_spec)
                if(ierr.ne.TEREC_SUCCESS) ierr=-18
               else
                ierr=-17
               endif
              else
               ierr=-16
              endif
             else
              if(.not.present(group_spec)) then
               call tensor%tens_rcrsv_ctor(tens_name,subspaces,hspaces,ierr,dim_extent)
               if(ierr.ne.TEREC_SUCCESS) ierr=-15
              else
               ierr=-14
              endif
             endif
            else
             ierr=-13
            endif
           else
            if(present(dim_group)) then
             if(size(dim_group).eq.trank) then
              if(present(group_spec)) then
               call tensor%tens_rcrsv_ctor(tens_name,subspaces,hspaces,ierr,dim_group=dim_group,group_spec=group_spec)
               if(ierr.ne.TEREC_SUCCESS) ierr=-12
              else
               ierr=-11
              endif
             else
              ierr=-10
             endif
            else
             if(.not.present(group_spec)) then
              call tensor%tens_rcrsv_ctor(tens_name,subspaces,hspaces,ierr)
              if(ierr.ne.TEREC_SUCCESS) ierr=-9
             else
              ierr=-8
             endif
            endif
           endif
!Set numeric data kind:
           if(ierr.eq.EXA_SUCCESS) then
            call tensor%set_data_type(data_kind,ierr); if(ierr.ne.TEREC_SUCCESS) ierr=-7
           endif
          else
           ierr=-6
          endif
!Construct the tensor instruction:
          if(ierr.eq.0) then
           tens_instr=>add_new_instruction(ip,ierr)
           if(ierr.eq.0) then
            write(jo,'(i11)') ip; flush(jo) !new instruction id number
            call tens_instr%tens_instr_ctor(TAVP_INSTR_TENS_CREATE,ierr,tensor,iid=ip)
            if(ierr.eq.0) then
!Issue the tensor instruction to TAVP:
             call issue_new_instruction(tens_instr,ierr); if(ierr.ne.0) ierr=-5
            else
             ierr=-4
            endif
           else
            ierr=-3
           endif
           tens_instr=>NULL()
          endif
         else
          ierr=-2
         endif
        else
         ierr=-1
        endif
        return
       end function exatns_tensor_create_tensor
!---------------------------------------------------------
       function exatns_tensor_destroy(tensor) result(ierr)
!Destroys a tensor.
        implicit none
        integer(INTD):: ierr                       !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor !inout: tensor
        class(tens_instr_mng_t), pointer:: tens_instr
        integer(INTL):: ip

        ierr=EXA_SUCCESS
        write(jo,'("#MSG(exatensor): New Instruction: DESTROY TENSOR: IP = ")',ADVANCE='NO'); flush(jo)
        if(tensor%is_set(ierr)) then
         if(ierr.eq.TEREC_SUCCESS) then
!Construct the tensor instruction:
          tens_instr=>add_new_instruction(ip,ierr)
          if(ierr.eq.0) then
           write(jo,'(i11)') ip; flush(jo) !new instruction id number
           call tens_instr%tens_instr_ctor(TAVP_INSTR_TENS_DESTROY,ierr,tensor,iid=ip)
           if(ierr.eq.0) then
!Issue the tensor instruction to TAVP:
            call issue_new_instruction(tens_instr,ierr)
            if(ierr.eq.0) then
             call tens_rcrsv_dtor(tensor) !local tensor is no longer needed
            else
             ierr=-5
            endif
           else
            ierr=-4
           endif
          else
           ierr=-3
          endif
          tens_instr=>NULL()
         else
          ierr=-2
         endif
        else
         ierr=-1
        endif
        return
       end function exatns_tensor_destroy
!---------------------------------------------------------------------------------
       function exatns_tensor_get(tensor,subspace_mlndx,tensor_slice) result(ierr)
!Returns a locally storable slice of a tensor.
        implicit none
        integer(INTD):: ierr                             !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor       !in: (distributed) tensor
        integer(INTL), intent(in):: subspace_mlndx(1:)   !in: subspace multi-index identifying the requested tensor slice
        type(tens_rcrsv_t), intent(inout):: tensor_slice !out: requested tensor slice stored locally
        integer(INTD):: tens_rank

        ierr=EXA_SUCCESS
        tens_rank=tensor%get_rank(ierr)
        if(ierr.eq.TEREC_SUCCESS) then
         if(size(subspace_mlndx).eq.tens_rank) then
          call quit(-1,'FATAL(exatensor:tensor_get): Not implemented yet!') !`Implement
         else
          ierr=EXA_ERR_INVALID_ARGS
         endif
        else
         ierr=EXA_ERR_INVALID_ARGS
        endif
        return
       end function exatns_tensor_get
!---------------------------------------------------------------
       function exatns_tensor_load(tensor,filename) result(ierr)
!Loads a tensor from an external storage.
        implicit none
        integer(INTD):: ierr                       !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor !inout: tensor
        character(*), intent(in):: filename        !in: file name

        ierr=EXA_SUCCESS
        call quit(-1,'FATAL(exatensor:tensor_load): Not implemented yet!') !`Implement
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
        call quit(-1,'FATAL(exatensor:tensor_save): Not implemented yet!') !`Implement
        return
       end function exatns_tensor_save
!------------------------------------------------------------
       function exatns_tensor_status(tensor,sts) result(ierr)
!Returns the current status of a tensor.
        implicit none
        integer(INTD):: ierr                    !out: error code
        type(tens_rcrsv_t), intent(in):: tensor !in: tensor
        type(tens_status_t), intent(out):: sts  !out: current tensor status

        ierr=EXA_SUCCESS
        call quit(-1,'FATAL(exatensor:tensor_status): Not implemented yet!') !`Implement
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
        class(tens_instr_mng_t), pointer:: tens_instr
        type(tens_transformation_t):: tens_init
        integer(INTL):: ip
        complex(8):: scal

        ierr=EXA_SUCCESS
        write(jo,'("#MSG(exatensor): New Instruction: INIT TENSOR (by scalar): IP = ")',ADVANCE='NO'); flush(jo)
        if(tensor%is_set(ierr)) then
         if(ierr.eq.TEREC_SUCCESS) then
!Construct the tensor transformation (initialization) object:
          call tens_init%set_argument(tensor,ierr)
          if(ierr.eq.TEREC_SUCCESS) then
           scal=(0d0,0d0); if(present(scalar)) scal=scalar
           call tens_init%set_method(ierr,scalar_value=scal,defined=.FALSE.) !simple initialization to a value
           if(ierr.eq.TEREC_SUCCESS) then
!Construct the tensor instruction:
            tens_instr=>add_new_instruction(ip,ierr)
            if(ierr.eq.0) then
             write(jo,'(i11)') ip; flush(jo) !new instruction id number
             call tens_instr%tens_instr_ctor(TAVP_INSTR_TENS_INIT,ierr,tens_init,iid=ip)
             if(ierr.eq.0) then
!Issue the tensor instruction to TAVP:
              call issue_new_instruction(tens_instr,ierr); if(ierr.ne.0) ierr=-7
             else
              ierr=-6
             endif
            else
             ierr=-5
            endif
            tens_instr=>NULL()
           else
            ierr=-4
           endif
          else
           ierr=-3
          endif
          call tens_transformation_dtor(tens_init)
         else
          ierr=-2
         endif
        else
         ierr=-1
        endif
        return
       end function exatns_tensor_init_scalar
!--------------------------------------------------------------------
       function exatns_tensor_init_method(tensor,method) result(ierr)
!Initializes the tensor via a user-defined (registered) method.
        implicit none
        integer(INTD):: ierr                       !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor !inout: tensor
        character(*), intent(in):: method          !in: registered method name
        class(tens_instr_mng_t), pointer:: tens_instr
        type(tens_transformation_t):: tens_init
        integer(INTL):: ip

        ierr=EXA_SUCCESS
        write(jo,'("#MSG(exatensor): New Instruction: INIT TENSOR (by method): IP = ")',ADVANCE='NO'); flush(jo)
        if(tensor%is_set(ierr)) then
         if(ierr.eq.TEREC_SUCCESS) then
!Construct the tensor transformation (initialization) object:
          call tens_init%set_argument(tensor,ierr)
          if(ierr.eq.TEREC_SUCCESS) then
#if !(defined(__GNUC__) && __GNUC__ < 8)
           call tens_init%set_method(ierr,defined=.FALSE.,method_name=method,method_map=method_map_f)
#endif
           if(ierr.eq.TEREC_SUCCESS) then
!Construct the tensor instruction:
            tens_instr=>add_new_instruction(ip,ierr)
            if(ierr.eq.0) then
             write(jo,'(i11)') ip; flush(jo) !new instruction id number
             call tens_instr%tens_instr_ctor(TAVP_INSTR_TENS_INIT,ierr,tens_init,iid=ip)
             if(ierr.eq.0) then
!Issue the tensor instruction to TAVP:
              call issue_new_instruction(tens_instr,ierr); if(ierr.ne.0) ierr=-7
             else
              ierr=-6
             endif
            else
             ierr=-5
            endif
            tens_instr=>NULL()
           else
            ierr=-4
           endif
          else
           ierr=-3
          endif
          call tens_transformation_dtor(tens_init)
         else
          ierr=-2
         endif
        else
         ierr=-1
        endif
        return
       end function exatns_tensor_init_method
!----------------------------------------------------------------------------
       function exatns_tensor_copy(tensor_out,tensor_in,pattern) result(ierr)
!Copies the content of one tensor into another tensor with an option of permutation,
!and/or slicing or insertion.
        implicit none
        integer(INTD):: ierr                           !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor_out !inout: output tensor
        type(tens_rcrsv_t), intent(inout):: tensor_in  !in: input tensor
        character(*), intent(in), optional:: pattern   !in: symbolic permutation pattern

        ierr=EXA_SUCCESS
        call quit(-1,'FATAL(exatensor:tensor_copy): Not implemented yet!') !`Implement
        return
       end function exatns_tensor_copy
!----------------------------------------------------------------------------
       function exatns_tensor_fold(tensor_out,tensor_in,pattern) result(ierr)
!Folds two or more dimensions of the input tensor, producing an output tensor
!of a lower rank (lower order). In other words, flattening.
        implicit none
        integer(INTD):: ierr                           !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor_out !inout: output tensor
        type(tens_rcrsv_t), intent(inout):: tensor_in  !in: input tensor
        character(*), intent(in):: pattern             !in: symbolic folding pattern

        ierr=EXA_SUCCESS
        call quit(-1,'FATAL(exatensor:tensor_fold): Not implemented yet!') !`Implement
        return
       end function exatns_tensor_fold
!------------------------------------------------------------------------------
       function exatns_tensor_unfold(tensor_out,tensor_in,pattern) result(ierr)
!Unfolds dimensions of the input tensor, producing an output tensor
!of a higher rank (higher order).
        implicit none
        integer(INTD):: ierr                           !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor_out !inout: output tensor
        type(tens_rcrsv_t), intent(inout):: tensor_in  !in: input tensor
        character(*), intent(in):: pattern             !in: symbolic unfolding pattern

        ierr=EXA_SUCCESS
        call quit(-1,'FATAL(exatensor:tensor_unfold): Not implemented yet!') !`Implement
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
!-----------------------------------------------------------------------------
       function exatns_tensor_add(tensor0,tensor1,pattern,factor) result(ierr)
!Operation: tensor0 += tensor1 * factor
        implicit none
        integer(INTD):: ierr                         !out: error code
        type(tens_rcrsv_t), intent(inout):: tensor0  !inout: tensor
        type(tens_rcrsv_t), intent(inout):: tensor1  !in: tensor
        character(*), intent(in), optional:: pattern !in: symbolic permutation pattern
        complex(8), intent(in), optional:: factor    !in: scalar factor

        ierr=EXA_SUCCESS
        call quit(-1,'FATAL(exatensor:tensor_add): Not implemented yet!') !`Implement
        return
       end function exatns_tensor_add
!----------------------------------------------------------------------------------------------------------
       function exatns_tensor_contract(tensor0,tensor1,tensor2,pattern,prefactor,restrictions) result(ierr)
!Operation: Tensor contraction: tensor0 += tensor1 * tensor2 * prefactor
        implicit none
        integer(INTD):: ierr
        type(tens_rcrsv_t), intent(inout):: tensor0       !inout: tensor
        type(tens_rcrsv_t), intent(inout):: tensor1       !in: tensor
        type(tens_rcrsv_t), intent(inout):: tensor2       !in: tensor
        character(*), intent(in):: pattern                !in: symbolic tensor contraction pattern
        complex(8), intent(in), optional:: prefactor      !in: scalar prefactor
        character(*), intent(in), optional:: restrictions !in: symbolic specification of operational index restrictions
        class(tens_instr_mng_t), pointer:: tens_instr
        type(tens_contraction_t):: tens_contr
        integer(INTD):: errc,cpl,conj_bits,contr_ptrn(1:MAX_TENSOR_RANK*2)
        integer(INTL):: ip
        logical:: check

        ierr=EXA_SUCCESS
        write(jo,'("#MSG(exatensor): New Instruction: CONTRACT TENSORS: IP = ")',ADVANCE='NO'); flush(jo)
        check=tensor0%is_set(errc); if(errc.ne.TEREC_SUCCESS.and.ierr.eq.EXA_SUCCESS) ierr=-12
        check=tensor1%is_set(errc).and.check; if(errc.ne.TEREC_SUCCESS.and.ierr.eq.EXA_SUCCESS) ierr=-11
        check=tensor2%is_set(errc).and.check; if(errc.ne.TEREC_SUCCESS.and.ierr.eq.EXA_SUCCESS) ierr=-10
        if(check.and.ierr.eq.EXA_SUCCESS) then
!Convert the symbolic tensor contraction pattern into a digital one used by TAL-SH:
         call get_contr_pattern(pattern,contr_ptrn,cpl,ierr,conj_bits) !conj_bits: tensor conjugation bits {0:D,1:L,2:R}
         if(ierr.eq.0) then
!Construct the tensor operation object:
          call tens_contr%set_argument(tensor0,ierr)
          if(ierr.eq.TEREC_SUCCESS) then
           call tens_contr%set_argument(tensor1,ierr)
           if(ierr.eq.TEREC_SUCCESS) then
            call tens_contr%set_argument(tensor2,ierr)
            if(ierr.eq.TEREC_SUCCESS) then
             if(present(prefactor)) then
              call tens_contr%set_contr_ptrn(contr_ptrn(1:cpl),ierr,prefactor,conjug=conj_bits)
             else
              call tens_contr%set_contr_ptrn(contr_ptrn(1:cpl),ierr,conjug=conj_bits)
             endif
             if(ierr.eq.TEREC_SUCCESS) then
              if(present(restrictions)) then
               !`Decode and set the operational index restrictions
              endif
             else
              if(VERBOSE) write(jo,'("#ERROR(exatns_tensor_contract): Contraction pattern setting error ",i11,'//&
                          &'" for pattern:",32(1x,i2))') ierr,contr_ptrn(1:cpl)
              ierr=-9
             endif
            else
             ierr=-8
            endif
           else
            ierr=-7
           endif
          else
           ierr=-6
          endif
         else
          ierr=-5
         endif
!Construct the tensor instruction:
         if(ierr.eq.EXA_SUCCESS) then
          tens_instr=>add_new_instruction(ip,ierr)
          if(ierr.eq.0) then
           write(jo,'(i11)') ip; flush(jo) !new instruction id number
           call tens_instr%tens_instr_ctor(TAVP_INSTR_TENS_CONTRACT,ierr,tens_contr,iid=ip)
           if(ierr.eq.0) then
!Issue the tensor instruction to TAVP:
            call issue_new_instruction(tens_instr,ierr); if(ierr.ne.0) ierr=-4
           else
            ierr=-3
           endif
          else
           ierr=-2
          endif
          tens_instr=>NULL()
         endif
!Destroy the tensor operation object:
         call tens_contraction_dtor(tens_contr)
        else
         if(ierr.eq.EXA_SUCCESS) ierr=-1
        endif
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
        call quit(-1,'FATAL(exatensor:tensor_unary_op): Not implemented yet!') !`Implement
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
        call quit(-1,'FATAL(exatensor:tensor_binary_op): Not implemented yet!') !`Implement
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
        call quit(-1,'FATAL(exatensor:tensor_ternary_op): Not implemented yet!') !`Implement
        return
       end function exatns_tensor_ternary_op
!----------------------------------------------------------
       function tavp_role_rank(global_rank,role) result(id)
!MPI process role numeration:
! EXA_WORKER: [0:exa_num_workers-1];
! EXA_MANAGER: [exa_num_workers:exa_num_workers+exa_num_managers-1];
! EXA_DRIVER: [exa_num_workers+exa_num_managers:total_number_of_mpi_processes-1].
        implicit none
        integer(INTD):: id                          !out: TAVP role-specific rank (within role intracommunicator)
        integer(INTD), intent(in):: global_rank     !in: global MPI rank of the TAVP
        integer(INTD), intent(out), optional:: role !out: TAVP role

        if(global_rank.ge.exa_num_workers+exa_num_managers) then
         if(present(role)) role=EXA_DRIVER
         id=global_rank-(exa_num_workers+exa_num_managers)
        else
         if(global_rank.ge.exa_num_workers) then
          if(present(role)) role=EXA_MANAGER
          id=global_rank-exa_num_workers
         else
          if(present(role)) role=EXA_WORKER
          id=global_rank
         endif
        endif
        return
       end function tavp_role_rank
!-------------------------------------------------------------
       function add_new_instruction(ip,ierr) result(new_instr)
!Adds a new (empty) instruction into the instruction log.
        implicit none
        class(tens_instr_mng_t), pointer:: new_instr !out: pointer to the newly added (empty) instruction
        integer(INTL), intent(out):: ip              !out: instruction id number
        integer(INTD), intent(out), optional:: ierr  !out: error code
        integer(INTD):: errc
        class(*), pointer:: uptr
        type(tens_instr_mng_t), target:: tens_instr_empty

        ip=-1_INTL; new_instr=>NULL()
        errc=instr_log%append(tens_instr_empty)
        if(errc.eq.GFC_SUCCESS) then
         errc=instr_log%reset_back()
         if(errc.eq.GFC_SUCCESS) then
          ip=instr_log%get_offset(errc)
          if(errc.eq.GFC_SUCCESS.and.ip.ge.0) then
           uptr=>instr_log%element_value(ip,errc)
           if(errc.eq.GFC_SUCCESS) then
            select type(uptr); class is(tens_instr_mng_t); new_instr=>uptr; end select
            if(.not.associated(new_instr)) errc=-5 !trap
            uptr=>NULL()
           else
            errc=-4
           endif
          else
           errc=-3
          endif
         else
          errc=-2
         endif
        else
         errc=-1
        endif
        if(present(ierr)) ierr=errc
        return
       end function add_new_instruction
!------------------------------------------------------
       subroutine issue_new_instruction(new_instr,ierr)
!Issues a new (defined) instruction to the root TAVP-MNG.
        implicit none
        class(tens_instr_mng_t), intent(in):: new_instr !in: new instruction
        integer(INTD), intent(out), optional:: ierr     !out: error code
        integer(INTD):: errc,ier,opcode
        type(obj_pack_t):: instr_packet
        type(comm_handle_t):: comm_hl

        call bytecode_out%acquire_packet(instr_packet,errc)
        if(errc.eq.0) then
         opcode=new_instr%get_code()
         call new_instr%encode(instr_packet,errc)
         if(errc.eq.0) then
          call bytecode_out%seal_packet(errc)
          if(errc.eq.0) then
           call bytecode_out%send(0,comm_hl,errc,tag=TAVP_DISPATCH_TAG,comm=drv_mng_comm) !send to root TAVP-MNG
           if(errc.eq.0) then
            call comm_hl%wait(errc)
            if(errc.eq.0) then
             if(opcode.ge.TAVP_ISA_TENS_FIRST.and.opcode.le.TAVP_ISA_TENS_LAST) num_tens_instr_issued=num_tens_instr_issued+1
             call comm_hl%clean(errc); if(errc.ne.0) errc=-7
            else
             errc=-6
            endif
           else
            errc=-5
           endif
          else
           errc=-4
          endif
         else
          errc=-3
         endif
        else
         errc=-2
        endif
        call bytecode_out%clean(ier); if(ier.ne.0.and.errc.eq.0) errc=-1
        if(present(ierr)) ierr=errc
        return
       end subroutine issue_new_instruction

      end module exatensor
