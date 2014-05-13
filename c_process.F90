!This module provides functionality for a Computing Process (C-PROCESS, CP).
!In essence, this is a single-node elementary tensor instruction scheduler (SETIS).
!AUTHOR: Dmitry I. Lyakh (Dmytro I. Liakh): quant4me@gmail.com
!REVISION: 2014/05/13
!CONCEPTS (CP workflow):
! - Each CP stores its own tensor blocks in TBB, with a possibility of disk dump.
! - LR sends a batch of ETI to be executed on this CP unit (CP MPI Process).
!   The location of each tensor-block operand in ETI is already given there,
!   as well as some other characteristics (approx. operation cost, memory requirements, etc.).
! - CP places the ETI received into 4 queues, based on the number of remote operands.
! - Each non-trivial ETI operand (tensor block) is assigned an AAR entry that points to
!   either a <tensor_block_t> entity (to be used on CPU or Intel Xeon Phi)
!   or a <tensBlck_t> entity (to be used on Nvidia GPU). Once the corresponding
!   data (tensor block) is in local memory of CP (either in TBB, or HAB, or DEB),
!   a new AAR entry is defined and all ETI operands, referring to it, are ready.
!   Data (tensor blocks) registered in AAR can be reused while locally present.
!   AAR is the only way to access a tensor block for ETI (the tensor block itself
!   can reside in either TBB, or HAB, or DEB).
! - ETI with all operands ready are issued for execution on an appropriate device:
!    - Issuing an ETI is asynchronous;
!    - Different ETI (with ready input arguments) can be prioritized;
!    - Each particular type of instruction has its own issuing workflow for each device kind;
!      The device is chosen based on the ETI cost, ETI cost/size ratio, and device availability;
!    - Once issued successfully, the ETI obtains a query handle that can be used for completion checks.
!    - AAR entries (together with associated data) used by active ETI must not be freed.
!    - If the ETI result is remote, its destination must be updated before
!      reporting to LR that the ETI has been completed.
!NOTES:
! - Data synchronization in an instance of <tensor_block_t> (Fortran)
!   associated with a Host Argument Buffer entry can allocate only regular CPU memory
!   (the one outside the pinned Host Argument buffer). Hence that newly allocated
!   memory cannot be used with Nvidia GPU.
!Acronyms:
! - CP - Computing MPI Process;
! - MT - Master Thread;
! - STCU - Slave Threads Computing Unit;
! - NVCU - Nvidia GPU Computing Unit;
! - XPCU - Intel Xeon Phi Computing Unit;
! - AMCU - AMD GPU Computing Unit;
! - ETI - Elementary Tensor Instruction;
! - ETIS - Elementary Tensor Instruction Scheduler (SETIS, LETIS, GETIS);
! - ETIQ - Elementary Tensor Instruction Queue;
! - HAB - Host Argument Buffer;
! - GAB - GPU Argument Buffer;
! - DEB - Data Exchange Buffer (additional buffer for temporary present tensor blocks);
! - TBB - Tensor Block Bank (storage);
! - AAR - Active Argument Register;
! - TAL - Tensor Algebra Library;
       module c_process
        use, intrinsic:: ISO_C_BINDING
        use tensor_algebra
        use dictionary
        use service
        use extern_names !`dependency to be removed
        implicit none
#ifndef NO_OMP
        integer, external, private:: omp_get_max_threads,omp_get_num_threads,omp_get_thread_num
#endif
!PARAMETERS:
 !Output:
        integer, private:: jo_cp=6 !default output
        logical, private:: verbose=.true.
 !General:
        integer, parameter, private:: CZ=C_SIZE_T
 !Elementary Tensor Instruction Scheduler (ETIS):
  !Host Argument Buffer (HAB):
        integer(C_SIZE_T), parameter, private:: max_hab_size=1024_CZ*1024_CZ*1024_CZ !max size in bytes of the HAB
        integer(C_SIZE_T), parameter, private:: min_hab_size=64_CZ*1024_CZ*1024_CZ   !min size in bytes of the HAB
  !Elementary Tensor Instruction Queue (ETIQ):
        integer(C_INT), parameter, private:: etiq_max_depth=65536 !max number of simultaneously scheduled ETI at this CP
        integer(C_INT), parameter, private:: etiq_loc_levels=4    !number of locality levels in ETIQ (senior)
        integer(C_INT), parameter, private:: etiq_cost_levels=5   !number of cost levels in ETIQ (minor)
        integer(C_INT), parameter, private:: etiq_levels=etiq_cost_levels*etiq_loc_levels !total number of levels in ETIQ
        integer(C_INT), parameter, private:: etiq_stcu_max_depth=1024 !max number of simultaneously scheduled ETI on STCU
        integer(C_INT), parameter, private:: etiq_nvcu_max_depth=16 !max number of simultaneously scheduled ETI on NVCU
        integer(C_INT), parameter, private:: etiq_xpcu_max_depth=32 !max number of simultaneously scheduled ETI on XPCU
        integer(C_INT), parameter, private:: stcu_max_units=64 !max number of STCU units
  !Tensor naming:
        integer(C_INT), parameter:: tensor_name_len=32 !max number of characters used for tensor names
  !Tensor instruction status (any negative value will correspond to a failure, designating the error code):
        integer, parameter:: instr_null=0              !uninitialized instruction
        integer, parameter:: instr_data_wait=1         !instruction is waiting for input data to arrive
        integer, parameter:: instr_ready_to_exec=2     !instruction is ready to be executed (input data has arrived)
        integer, parameter:: instr_issued=3            !instruction is issued for execution on some computing unit
        integer, parameter:: instr_completed=4         !instruction has completed (but the result may need to be remotely uploaded)
        integer, parameter:: instr_dead=5              !instruction can be safely removed from the queue
  !Tensor instruction code:
        integer, parameter:: instr_tensor_init=1
        integer, parameter:: instr_tensor_norm1=2
        integer, parameter:: instr_tensor_norm2=3
        integer, parameter:: instr_tensor_min=4
        integer, parameter:: instr_tensor_max=5
        integer, parameter:: instr_tensor_scale=6
        integer, parameter:: instr_tensor_slice=7
        integer, parameter:: instr_tensor_insert=8
        integer, parameter:: instr_tensor_trace=9
        integer, parameter:: instr_tensor_copy=10
        integer, parameter:: instr_tensor_add=11
        integer, parameter:: instr_tensor_cmp=12
        integer, parameter:: instr_tensor_contract=13
!TYPES:
 !Computing unit (CU) identifier:
        type cu_t
         integer(C_INT):: device_type  !device type: see tensor_algebra_gpu_nvidia.inc
         integer(C_INT):: unit_number  !logical number of the device of the above type (0..max)
        end type cu_t
 !Tensor block identifier (key):
        type tens_blck_id_t
         character(LEN=tensor_name_len):: tens_name  !tensor name
         integer(C_INT), allocatable:: tens_mlndx(:) !tensor block multiindex (bases)
        end type tens_blck_id_t
 !Tensor Block Bank (TBB):
  !Entry of TBB (named tensor block):
        type, private:: tbb_entry_t
         type(tensor_block_t), private:: tens_blck !fortran tensor block
         integer, private:: file_handle=0          !file handle (0: stored in RAM)
         integer(8), private:: file_offset         !file offset where the corresponding packet is stored (if on disk)
        end type tbb_entry_t
 !Elementary Tensor Instruction Scheduler (ETIS):
  !Locally present tensor argument:
        type, private:: tens_arg_t
         type(tensor_block_t), pointer, private:: tens_blck_f=>NULL() !pointer to a tensor_block in TBB, or HAB, or DEB
         type(C_PTR), private:: tens_blck_c=C_NULL_PTR !C pointer to tensBlck_t in HAB (see "tensor_algebra_gpu_nvidia.h")
         integer(C_INT), private:: buf_entry_host !HAB entry number where the tensor block resides as a packet (-1: undefined)
         logical, private:: in_use                !.true. if this AAR entry is currently in use by some issued ETI
         integer, private:: times_used            !total number of times this AAR entry have been used since creation
        end type tens_arg_t
  !Tensor operand type (component of ETI):
        type tens_operand_t
         type(tens_blck_id_t):: tens_blck_id !tensor block identifier (key): set by LR
         integer(C_INT):: op_host            !MPI process rank where the tensor operand resides: set by LR
         integer(C_SIZE_T):: op_pack_size    !packed size of the tensor operand in bytes: computed by LR
         integer(C_INT):: op_tag             !MPI message tag by which the tensor operand is to be delivered (-1: local): set by LR
         integer(C_INT):: op_price           !current price of the tensor operand (tensor block): set by LR
         type(tens_arg_t), pointer, private:: op_aar_entry=>NULL() !AAR entry assigned to the tensor operand: set by CP
        end type tens_operand_t
  !Dispatched elementary tensor instruction (ETI):
        type tens_instr_t
         integer:: instr_code            !tensor instruction code (see above): set by GR
         integer:: instr_status          !tensor instruction status (see above): set by CP
         integer:: instr_priority        !tensor instruction priority: set by LR
         real(8):: instr_cost            !approx. instruction computational cost (FLOPs): set by LR
         real(8):: instr_size            !approx. instruction memory demands (Bytes): set by LR
         real(4):: instr_time_beg        !time the instruction was scheduled: set by CP
         real(4):: instr_time_end        !time the instruction was completed: set by CP
         integer, private:: instr_handle !instruction handle to query its status: set by CP
         integer, private:: args_ready   !each bit is set to 1 when the corresponding operand is in AAR: set by CP
         type(tens_operand_t):: tens_op0 !tensor-block operand #0
         type(tens_operand_t):: tens_op1 !tensor-block operand #1
         type(tens_operand_t):: tens_op2 !tensor-block operand #2
         type(tens_operand_t):: tens_op3 !tensor-block opearnd #3
        end type tens_instr_t
  !Elementary tensor instruction queue (ETIQ):
        type, private:: etiq_t
         integer(C_INT), private:: depth=0                   !total number of ETIQ entries
         integer(C_INT), private:: scheduled=0               !total number of active ETIQ entries
         integer(C_INT), private:: ffe_sp=0                  !ETIQ FFE stack pointer
         integer(C_INT), private:: last(0:etiq_levels-1)=0   !last scheduled instruction for each ETIQ level
         integer(C_INT), private:: ip(0:etiq_levels-1)=0     !instruction pointers for each ETIQ level
         integer(C_INT), private:: ic(0:etiq_levels-1)=0     !instruction counters for each ETIQ level
         integer(C_INT), allocatable, private:: ffe_stack(:) !ETIQ FFE stack
         integer(C_INT), allocatable, private:: next(:)      !ETIQ next linking
         type(tens_instr_t), allocatable, private:: eti(:)   !elementary tensor instructions
        end type etiq_t
  !In-order elementary tensor instruction queue for specific computing units:
        type, private:: etiq_cu_t
         integer(C_INT), private:: depth=0                    !total number of entries of the queue
         integer(C_INT), private:: scheduled=0                !total number of active entries (occupied)
         integer(C_INT), private:: ip=0                       !instruction pointer
         integer(C_INT), allocatable, private:: etiq_entry(:) !number of the ETIQ entry where the ETI is located
         integer(C_INT), allocatable, private:: cu_id(:)      !computing unit # (within its type) which ETI is issued to
        end type etiq_cu_t
 !STCU unit:
        type, private:: stcu_unit_t
         integer(C_INT), private:: num_omp_thrds  !number of OMP threads assigned to the STCU unit
         integer(C_INT), private:: etiq_entry     !ETIQ entry number assigned to the STCU unit
        end type stcu_unit_t
!DATA:
 !Tensor Block Bank (TBB):
        type(dict_t), private:: tbb
 !Elementary Tensor Instruction Scheduler (ETIS):
  !Active Argument Register (AAR):
        type(dict_t), private:: aar
  !Elementary Tensor Instruction Queue (ETIQ):
        type(etiq_t), private:: etiq
  !STCU ETI queue (ETIQ_STCU):
        type(etiq_cu_t), private:: etiq_stcu
  !NVCU ETI queue (ETIQ_NVCU:
        type(etiq_cu_t), private:: etiq_nvcu
  !XPCU ETI queue (ETIQ_XPCU):
        type(etiq_cu_t), private:: etiq_xpcu
  !Host Argument Buffer (HAB):
        integer(C_SIZE_T), private:: hab_size=0  !actual size in bytes of the Host argument buffer (HAB)
        integer(C_INT), private:: max_hab_args=0 !max number of arguments (of lowest-size level) that can fit in HAB
 !CP runtime:
        integer, private:: mt_error              !master thread error (if != 0)
        integer, private:: stcu_error            !slave threads computing unit error (if != 0)
        integer, private:: stcu_num_units        !current number of STCU units
        type(stcu_unit_t), private:: stcu_units(0:stcu_max_units-1) !configuration of STCU units
!------------------------------------------------------------------
       contains
!CODE:
        subroutine c_proc_life(ierr)
!This subroutine implements C-PROCESS life cycle.
!INPUT (external):
! - jo - global output device (service.mod);
! - impir - MPI rank of the process (service.mod);
! - impis - global MPI comminicator size (service.mod);
! - max_threads - max number of threads available to the current MPI process (service.mod);
! - {gpu_start:gpu_start+gpu_count-1} - range of Nvidia GPUs assigned to the current MPI process (service.mod);
!OUTPUT:
! - ierr - error code (0:success);
! - Executed elementary tensor instructions.
        implicit none
        integer, intent(inout):: ierr
!------------------------------------
        integer(C_INT), parameter:: max_arg_buf_levels=256 !max number of argument buffer levels (do not exceed C values)
!------------------------------------------
        integer(C_INT) i,j,k,l,m,n,err_code
        integer(C_SIZE_T) blck_sizes(0:max_arg_buf_levels-1)
        integer thread_num,stcu_num_units
        real(8) tm

        ierr=0; jo_cp=jo
        write(jo_cp,'("#MSG(c_process::c_proc_life): I am a C-process (Computing MPI Process): MPI rank = ",i7)') impir
!Initialization:
!       write(jo_cp,'("#MSG(c_process::c_proc_life): Initialization:")')
 !Init TAL infrastructure (TAL buffers, cuBLAS, etc.):
        write(jo_cp,'("#MSG(c_process::c_proc_life): Allocating argument buffers ... ")',advance='no')
        tm=thread_wtime()
        hab_size=max_hab_size !desired (max) HAB size
        i=arg_buf_allocate(hab_size,max_hab_args,gpu_start,gpu_start+gpu_count-1); if(i.ne.0) then; write(jo_cp,'("Failed!")'); call c_proc_quit(1); return; endif
        tm=thread_wtime()-tm; write(jo_cp,'("Ok(",F4.1," sec): Host argument buffer size (B) = ",i11)') tm,hab_size
        if(hab_size.lt.min_hab_size.or.max_hab_args.lt.6) then
         write(jo_cp,'("#FATAL(c_process::c_proc_life): Host argument buffer size is lower than minimally allowed: ",i11,1x,i11,1x,i11)') min_hab_size,hab_size,max_hab_args
         call c_proc_quit(2); return
        endif
  !Check Host argument buffer (HAB) levels:
        i=get_blck_buf_sizes_host(blck_sizes)
        if(i.le.0.or.i.gt.max_arg_buf_levels) then
         write(jo_cp,'("#ERROR(c_process::c_proc_life): Invalid number of Host argument buffer levels: ",i11,1x,i11)') max_arg_buf_levels,i
         call c_proc_quit(3); return
        else
         write(jo_cp,'("#MSG(c_process::c_proc_life): Number of Host argument buffer levels = ",i4,":")') i
         do j=0,i-1
          write(jo_cp,'("#MSG(c_process::c_proc_life): Level ",i4,": Size (B) = ",i11)') j,blck_sizes(j)
         enddo
        endif
        write(jo_cp,'("#MSG(c_process::c_proc_life): Max number of arguments in the Host argument buffer = ",i6)') max_hab_args
#ifndef NO_GPU
  !Check GPU argument buffer (GAB) levels:
        do j=gpu_start,gpu_start+gpu_count-1
         if(gpu_is_mine(j).ne.NOT_REALLY) then
          i=get_blck_buf_sizes_gpu(j,blck_sizes)
          if(i.le.0.or.i.gt.max_arg_buf_levels) then
           write(jo_cp,'("#ERROR(c_process::c_proc_life): Invalid number of GPU argument buffer levels: ",i11,1x,i11)') max_arg_buf_levels,i
           call c_proc_quit(4); return
          else
           write(jo_cp,'("#MSG(c_process::c_proc_life): Number of GPU#",i2," argument buffer levels = ",i4,":")') j,i
           do k=0,i-1
            write(jo_cp,'("#MSG(c_process::c_proc_life): Level ",i4,": Size (B) = ",i11)') k,blck_sizes(k)
           enddo
          endif
         endif
        enddo
#endif
 !Init ETIQ:
        write(jo_cp,'("#MSG(c_process::c_proc_life): Allocating Tensor Instruction Queue (ETIQ) ... ")',advance='no')
        tm=thread_wtime()
        allocate(etiq%eti(1:etiq_max_depth),etiq%next(1:etiq_max_depth),etiq%ffe_stack(1:etiq_max_depth),STAT=ierr)
        if(ierr.ne.0) then; write(jo_cp,'("#ERROR(c_process::c_proc_life): ETIQ allocation failed!")'); call c_proc_quit(5); return; endif
        etiq%depth=etiq_max_depth; etiq%scheduled=0; etiq%last(:)=0; etiq%ip(:)=0; etiq%ic(:)=0
        do i=1,etiq_max_depth; etiq%ffe_stack(i)=i; enddo; etiq%ffe_sp=1
        tm=thread_wtime()-tm; write(jo_cp,'("Ok(",F4.1," sec): ETIQ total depth = ",i7)') tm,etiq_max_depth
        write(jo_cp,'("#MSG(c_process::c_proc_life): Number of ETIQ channels = ",i3)') etiq_levels
!Test C-process functionality (debug):
!        call c_proc_test(ierr); if(ierr.ne.0) then; write(jo_cp,'("#ERROR(c_process::c_proc_life): C-process functionality test failed: ",i7)') ierr; call c_proc_quit(6); return; endif
!        call run_benchmarks(ierr); if(ierr.ne.0) then; write(jo_cp,'("#ERROR(c_process::c_proc_life): C-process benchmarking failed: ",i7)') ierr; call c_proc_quit(7); return; endif
!------------------------------
!LIFE:
        ierr=0
#ifndef NO_OMP
        call omp_set_dynamic(.false.); call omp_set_nested(.true.)
        n=omp_get_max_threads()
        if(n.ge.2.and.n.eq.max_threads) then
!Init STCU ETIQ:
         write(jo_cp,'("#MSG(c_process::c_proc_life): Allocating STCU ETIQ ... ")',advance='no')
         tm=thread_wtime()
         allocate(etiq_stcu%etiq_entry(1:etiq_stcu_max_depth),etiq_stcu%cu_id(1:etiq_stcu_max_depth),STAT=ierr)
         if(ierr.ne.0) then; write(jo_cp,'("#ERROR(c_process::c_proc_life): ETIQ_STCU allocation failed!")'); call c_proc_quit(8); return; endif
         do i=1,etiq_stcu_max_depth; etiq_stcu%etiq_entry(i)=0; enddo
         do i=1,etiq_stcu_max_depth; etiq_stcu%cu_id(i)=0; enddo
         etiq_stcu%depth=etiq_stcu_max_depth; etiq_stcu%scheduled=0; etiq_stcu%ip=0; stcu_num_units=0
         tm=thread_wtime()-tm; write(jo_cp,'("Ok(",F4.1," sec): STCU ETIQ depth = ",i6)') tm,etiq_stcu%depth
         write(jo_cp,'("#MSG(c_process::c_proc_life): Max number of STCU MIMD units = ",i5)') min(max_threads-1,stcu_max_units)
!Init NVCU ETIQ:
#ifndef NO_GPU
         write(jo_cp,'("#MSG(c_process::c_proc_life): Allocating NVCU ETIQ ... ")',advance='no')
         tm=thread_wtime()
         allocate(etiq_nvcu%etiq_entry(1:etiq_nvcu_max_depth),etiq_nvcu%cu_id(1:etiq_nvcu_max_depth),STAT=ierr)
         if(ierr.ne.0) then; write(jo_cp,'("#ERROR(c_process::c_proc_life): ETIQ_NVCU allocation failed!")'); call c_proc_quit(9); return; endif
         do i=1,etiq_nvcu_max_depth; etiq_nvcu%etiq_entry(i)=0; enddo
         do i=1,etiq_nvcu_max_depth; etiq_nvcu%cu_id(i)=0; enddo
         etiq_nvcu%depth=etiq_nvcu_max_depth; etiq_nvcu%scheduled=0; etiq_nvcu%ip=0
         tm=thread_wtime()-tm; write(jo_cp,'("Ok(",F4.1," sec): NVCU ETIQ depth = ",i6)') tm,etiq_nvcu%depth
#else
         etiq_nvcu%depth=0; etiq_nvcu%scheduled=0; etiq_nvcu%ip=0
#endif
!Init XPCU ETIQ:
#ifndef NO_PHI
         write(jo_cp,'("#MSG(c_process::c_proc_life): Allocating XPCU ETIQ ... ")',advance='no')
         tm=thread_wtime()
         allocate(etiq_xpcu%etiq_entry(1:etiq_xpcu_max_depth),etiq_xpcu%cu_id(1:etiq_xpcu_max_depth),STAT=ierr)
         if(ierr.ne.0) then; write(jo_cp,'("#ERROR(c_process::c_proc_life): ETIQ_XPCU allocation failed!")'); call c_proc_quit(10); return; endif
         do i=1,etiq_xpcu_max_depth; etiq_xpcu%etiq_entry(i)=0; enddo
         do i=1,etiq_xpcu_max_depth; etiq_xpcu%cu_id(i)=0; enddo
         etiq_xpcu%depth=etiq_xpcu_max_depth; etiq_xpcu%scheduled=0; etiq_xpcu%ip=0
         tm=thread_wtime()-tm; write(jo_cp,'("Ok(",F4.1," sec): XPCU ETIQ depth = ",i6)') tm,etiq_xpcu%depth
#else
         etiq_xpcu%depth=0; etiq_xpcu%scheduled=0; etiq_xpcu%ip=0
#endif
!Begin active life:
!$OMP PARALLEL DEFAULT(SHARED) NUM_THREADS(2)
         n=omp_get_num_threads()
         if(n.eq.2) then
          thread_num=omp_get_thread_num()
!Master thread:
          if(thread_num.eq.0) then
           mt_error=0           

!Slave threads:
          else !thread_num=1
           stcu_error=0
           life_slaves: do
            
            exit life_slaves
           enddo life_slaves
          endif
         else
          write(jo_cp,'("#ERROR(c_process::c_proc_life): invalid initial number of threads: ",i6)') n
          ierr=-2
         endif
!$OMP END PARALLEL
        else
         write(jo_cp,'("#ERROR(c_process::c_proc_life): invalid max number of threads: ",i6,1x,i6)') n,max_threads
         ierr=-1
        endif
        if(ierr.ne.0) then; call c_proc_quit(998); return; endif
#else
        write(jo_cp,'("#FATAL(c_process::c_proc_life): cannot function without OpenMP!")')
        call c_proc_quit(999); return
#endif
!-------------------
        write(jo_cp,'("#MSG(c_process::c_proc_life): Cleaning ... ")',advance='no')
        call c_proc_quit(0); if(ierr.eq.0) write(jo_cp,'("Ok")')
        return

        contains

         subroutine c_proc_quit(errc)
         integer, intent(in):: errc
         integer(C_INT) j0,j1
         ierr=errc
!HAB clean up:
         j0=arg_buf_clean_host(); if(j0.ne.0) then; write(jo_cp,'("#WARNING(c_process::c_proc_life:c_proc_quit): Host Argument buffer is not clean!")'); ierr=ierr+100; endif
#ifndef NO_GPU
!GAB clean up:
         do j1=gpu_start,gpu_start+gpu_count-1
          if(gpu_is_mine(j1).ne.NOT_REALLY) then
           j0=arg_buf_clean_gpu(j1); if(j0.ne.0) then; write(jo_cp,'("#WARNING(c_process::c_proc_life:c_proc_quit): GPU#",i2," Argument buffer is not clean!")') j1; ierr=ierr+100*(2+j1); endif
          endif
         enddo
#endif
!HAB/GAB deallocation:
         j0=arg_buf_deallocate(gpu_start,gpu_start+gpu_count-1); hab_size=0_CZ; max_hab_args=0
         if(j0.ne.0) then; write(jo_cp,'("#ERROR(c_process::c_proc_life:c_proc_quit): Deallocation of argument buffers failed!")'); ierr=ierr+7000; endif
!ETIQ:
         j1=0
         if(allocated(etiq%ffe_stack)) deallocate(etiq%ffe_stack,STAT=j0); if(j0.ne.0) j1=j1+1
         if(allocated(etiq%next)) deallocate(etiq%next,STAT=j0); if(j0.ne.0) j1=j1+1
         if(allocated(etiq%eti)) deallocate(etiq%eti,STAT=j0); if(j0.ne.0) j1=j1+1 !`relies on inheritant deallocation (Fortran 2003 automatic deallocation of all allocated components)
         if(j1.ne.0) then;  write(jo_cp,'("#ERROR(c_process::c_proc_life:c_proc_quit): ETIQ deallocation failed!")'); ierr=ierr+20000; endif
         etiq%depth=0; etiq%scheduled=0; etiq%ffe_sp=0; etiq%last(:)=0; etiq%ip(:)=0; etiq%ic(:)=0
 !ETIQ_STCU:
         j1=0
         if(allocated(etiq_stcu%etiq_entry)) deallocate(etiq_stcu%etiq_entry,STAT=j0); if(j0.ne.0) j1=j1+1
         if(allocated(etiq_stcu%cu_id)) deallocate(etiq_stcu%cu_id,STAT=j0); if(j0.ne.0) j1=j1+1
         if(j1.ne.0) then;  write(jo_cp,'("#ERROR(c_process::c_proc_life:c_proc_quit): ETIQ_STCU deallocation failed!")'); ierr=ierr+50000; endif
         etiq_stcu%depth=0; etiq_stcu%scheduled=0; etiq_stcu%ip=0; stcu_num_units=0
 !ETIQ_NVCU:
#ifndef NO_GPU
         j1=0
         if(allocated(etiq_nvcu%etiq_entry)) deallocate(etiq_nvcu%etiq_entry,STAT=j0); if(j0.ne.0) j1=j1+1
         if(allocated(etiq_nvcu%cu_id)) deallocate(etiq_nvcu%cu_id,STAT=j0); if(j0.ne.0) j1=j1+1
         if(j1.ne.0) then;  write(jo_cp,'("#ERROR(c_process::c_proc_life:c_proc_quit): ETIQ_NVCU deallocation failed!")'); ierr=ierr+100000; endif
#endif
         etiq_nvcu%depth=0; etiq_nvcu%scheduled=0; etiq_nvcu%ip=0
 !ETIQ_XPCU:
#ifndef NO_PHI
         j1=0
         if(allocated(etiq_xpcu%etiq_entry)) deallocate(etiq_xpcu%etiq_entry,STAT=j0); if(j0.ne.0) j1=j1+1
         if(allocated(etiq_xpcu%cu_id)) deallocate(etiq_xpcu%cu_id,STAT=j0); if(j0.ne.0) j1=j1+1
         if(j1.ne.0) then;  write(jo_cp,'("#ERROR(c_process::c_proc_life:c_proc_quit): ETIQ_XPCU deallocation failed!")'); ierr=ierr+300000; endif
#endif
         etiq_xpcu%depth=0; etiq_xpcu%scheduled=0; etiq_xpcu%ip=0
!AAR:
         j0=aar%destroy(destruct_key_func=destructor)
         if(j0.ne.0) then; write(jo_cp,'("#ERROR(c_process::c_proc_life:c_proc_quit): AAR destruction failed!")'); ierr=ierr+500000; endif
!TBB:
         j0=tbb%destroy(destruct_key_func=destructor,destruct_val_func=destructor)
         if(j0.ne.0) then; write(jo_cp,'("#ERROR(c_process::c_proc_life:c_proc_quit): TBB destruction failed!")'); ierr=ierr+1000000; endif
         return
         end subroutine c_proc_quit

        end subroutine c_proc_life
!----------------------------------------
        integer function destructor(item) !universal destructor
        implicit none
        class(*):: item
        integer i
        destructor=0
        select type (item)
        type is (tens_blck_id_t)
         if(allocated(item%tens_mlndx)) then
          deallocate(item%tens_mlndx,STAT=i); if(i.ne.0) then; destructor=1; return; endif
         endif
        type is (tbb_entry_t)
         call tensor_block_destroy(item%tens_blck,i); if(i.ne.0) then; destructor=2; return; endif
        end select
        return
        end function destructor
!--------------------------------------------------------------------------
        subroutine tens_blck_pack(tens,dtk,packet_size,pptr,entry_num,ierr)
!This subroutine packs a tensor block <tens> (tensor_block_t) into a linear packet
!and places it in the Host argument buffer, if there is enough free space there:
! TENSOR_BLOCK_T --> PACKET (entry in the Host Argument Buffer)
!INPUT:
! - tens - tensor block (tensor_block_t);
! - dtk - data kind ('r4','r8','c8'), saying elements of which data kind to place in the argument buffer;
!OUTPUT:
! - packet_size - size of the tensor block packet in bytes;
! - pptr - C pointer to the corresponding Host argument buffer entry;
! - entry_num - Host argument buffer entry number (allocated here);
! - packet in the argument buffer space at [pptr..pptr+packet_size-1];
! - ierr - error code (0:success).
!NOTES:
! - Tensor block packet structure:
!   C_SIZE_T: tensor packet size (bytes);
!   C_INT: data kind (4:float; 8:double; 16:double_complex);
!   C_SIZE_T: tensor block size (number of tensor elements stored);
!   C_INT: tensor block rank (number of dimensions);
!   C_INT(0:rank-1): dimension extents;
!   C_INT(0:rank-1): dimension dividers;
!   C_INT(0:rank-1): dimension groups;
!   C_INT(0:rank-1): dimension bases (absolute offsets for each dimension);
!   C_INT(0:rank-1): dimension permutation (for internal use only!);
!   C_DOUBLE(0:1): scalar value (c8);
!   {C_FLOAT|C_DOUBLE|C_DOUBLE(2)}: tensor block elements (count = tensor block size).
! - For scalars tensors, not only %scalar_value will be saved,
!   but also a single element of the corresponding data type.
        implicit none
        type(tensor_block_t), intent(inout):: tens
        character(2), intent(in):: dtk
        integer(C_SIZE_T), intent(out):: packet_size
        type(C_PTR), intent(out):: pptr
        integer(C_INT), intent(out):: entry_num
        integer, intent(inout):: ierr
        integer i
        integer(C_SIZE_T) s0
        integer(C_INT) err_code
        type(C_PTR) c_addr
        integer(C_SIZE_T) off_packet_size; integer(C_SIZE_T), pointer:: ptr_packet_size=>NULL()
        integer(C_INT) data_kind; integer(C_SIZE_T) off_data_kind; integer(C_INT), pointer:: ptr_data_kind=>NULL()
        integer(C_SIZE_T) elems_count; integer(C_SIZE_T) off_elems_count; integer(C_SIZE_T), pointer:: ptr_elems_count=>NULL()
        integer(C_INT) trank; integer(C_SIZE_T) off_trank; integer(C_INT), pointer:: ptr_trank=>NULL()
        integer(C_INT) dims(1:max_tensor_rank); integer(C_SIZE_T) off_dims; integer(C_INT), pointer:: ptr_dims(:)=>NULL()
        integer(C_INT) divs(1:max_tensor_rank); integer(C_SIZE_T) off_divs; integer(C_INT), pointer:: ptr_divs(:)=>NULL()
        integer(C_INT) grps(1:max_tensor_rank); integer(C_SIZE_T) off_grps; integer(C_INT), pointer:: ptr_grps(:)=>NULL()
        integer(C_INT) base(1:max_tensor_rank); integer(C_SIZE_T) off_base; integer(C_INT), pointer:: ptr_base(:)=>NULL()
        integer(C_SIZE_T) off_prmn;
        complex(8) sclr; integer(C_SIZE_T) off_sclr; complex(8), pointer:: ptr_sclr=>NULL()
        real(4) elems_r4; integer(C_SIZE_T) off_elems_r4; real(4), pointer:: ptr_elems_r4(:)=>NULL()
        real(8) elems_r8; integer(C_SIZE_T) off_elems_r8; real(8), pointer:: ptr_elems_r8(:)=>NULL()
        complex(8) elems_c8; integer(C_SIZE_T) off_elems_c8; complex(8), pointer:: ptr_elems_c8(:)=>NULL()

        ierr=0
!	write(jo_cp,'("#DEBUG(c_process::tens_blck_pack): Creating a tensor block packet:")') !debug
!Compute the size (in bytes) of the tensor block packet (linearized tensor_block_t):
        off_packet_size=0; packet_size=0
        off_data_kind=off_packet_size+sizeof(packet_size)
        select case(dtk)
        case('r4'); data_kind=4;
        case('r8'); data_kind=8;
        case('c8'); data_kind=16;
        case default; ierr=1; return; !invalid data kind
        end select
        off_elems_count=off_data_kind+sizeof(data_kind)
        elems_count=int(tensor_shape_size(tens,i),C_SIZE_T); if(i.ne.0) then; ierr=2; return; endif !unable to get tensor block size (elements)
        if(elems_count.ne.tens%tensor_block_size.or.elems_count.lt.0) then; ierr=3; return; endif
        off_trank=off_elems_count+sizeof(elems_count)
        trank=tens%tensor_shape%num_dim; if(trank.lt.0.or.(trank.eq.0.and.elems_count.ne.1).or.(trank.gt.0.and.elems_count.lt.1)) then; ierr=4; return; endif
        off_dims=off_trank+sizeof(trank)
        if(trank.gt.0) dims(1:trank)=tens%tensor_shape%dim_extent(1:trank)
        off_divs=off_dims+sizeof(dims(1))*trank
        if(trank.gt.0) divs(1:trank)=tens%tensor_shape%dim_divider(1:trank)
        off_grps=off_divs+sizeof(divs(1))*trank
        if(trank.gt.0) grps(1:trank)=tens%tensor_shape%dim_group(1:trank)
        off_base=off_grps+sizeof(grps(1))*trank
        if(trank.gt.0) base(1:trank)=0 !`Set proper bases
        off_prmn=off_base+sizeof(base(1))*trank
        off_sclr=off_prmn+sizeof(dims(1))*trank
        sclr=tens%scalar_value
        select case(dtk)
        case('r4'); off_elems_r4=off_sclr+sizeof(sclr); packet_size=off_elems_r4+sizeof(elems_r4)*elems_count
        case('r8'); off_elems_r8=off_sclr+sizeof(sclr); packet_size=off_elems_r8+sizeof(elems_r8)*elems_count
        case('c8'); off_elems_c8=off_sclr+sizeof(sclr); packet_size=off_elems_c8+sizeof(elems_c8)*elems_count
        end select
!DEBUG begin:
!	write(jo_cp,'(" Data kind    : ",i9,1x,i9)') off_data_kind,data_kind
!	write(jo_cp,'(" Element count: ",i9,1x,i9)') off_elems_count,elems_count
!	write(jo_cp,'(" Tensor rank  : ",i9,1x,i9)') off_trank,trank
!	write(jo_cp,'(" Dims         :",i9,3x,32(1x,i4))') off_dims,dims(1:trank)
!	write(jo_cp,'(" Divs         :",i9,3x,32(1x,i4))') off_divs,divs(1:trank)
!	write(jo_cp,'(" Grps         :",i9,3x,32(1x,i4))') off_grps,grps(1:trank)
!	write(jo_cp,'(" Base         :",i9,3x,32(1x,i4))') off_base,base(1:trank)
!	write(jo_cp,'(" Scalar       : ",i9,3x,D25.14,1x,D25.14)') off_sclr,sclr
!	select case(dtk)
!	case('r4'); write(jo_cp,'(" Elements(r4) : ",i9,1x,i9)') off_elems_r4
!	case('r8'); write(jo_cp,'(" Elements(r8) : ",i9,1x,i9)') off_elems_r8
!	case('c8'); write(jo_cp,'(" Elements(c8) : ",i9,1x,i9)') off_elems_c8
!	end select
!	write(jo_cp,'(" Packet size  : ",i9,1x,i9)') off_packet_size,packet_size
!DEBUG end.
!Allocate an argument buffer space on Host:
        err_code=get_buf_entry_host(packet_size,pptr,entry_num); if(err_code.ne.0) then; ierr=5; return; endif
!	write(jo_cp,'("#DEBUG(c_process::tens_blck_pack): Host argument buffer entry obtained: ",i7)') entry_num !debug
        if(entry_num.lt.0) then; ierr=6; return; endif
!Transfer the data into the Host argument buffer:
        c_addr=ptr_offset(pptr,off_packet_size)
!	call print_c_ptr(c_addr) !debug
        call c_f_pointer(c_addr,ptr_packet_size); ptr_packet_size=packet_size; nullify(ptr_packet_size)
        c_addr=ptr_offset(pptr,off_data_kind)
!	call print_c_ptr(c_addr) !debug
        call c_f_pointer(c_addr,ptr_data_kind); ptr_data_kind=data_kind; nullify(ptr_data_kind)
        c_addr=ptr_offset(pptr,off_elems_count)
!	call print_c_ptr(c_addr) !debug
        call c_f_pointer(c_addr,ptr_elems_count); ptr_elems_count=elems_count; nullify(ptr_elems_count)
        c_addr=ptr_offset(pptr,off_trank)
!	call print_c_ptr(c_addr) !debug
        call c_f_pointer(c_addr,ptr_trank); ptr_trank=trank; nullify(ptr_trank)
        if(trank.gt.0) then
         c_addr=ptr_offset(pptr,off_dims)
!	 call print_c_ptr(c_addr) !debug
         call c_f_pointer(c_addr,ptr_dims,shape=[trank]); ptr_dims(1:trank)=dims(1:trank); nullify(ptr_dims)
         c_addr=ptr_offset(pptr,off_divs)
!	 call print_c_ptr(c_addr) !debug
         call c_f_pointer(c_addr,ptr_divs,shape=[trank]); ptr_divs(1:trank)=divs(1:trank); nullify(ptr_divs)
         c_addr=ptr_offset(pptr,off_grps)
!	 call print_c_ptr(c_addr) !debug
         call c_f_pointer(c_addr,ptr_grps,shape=[trank]); ptr_grps(1:trank)=grps(1:trank); nullify(ptr_grps)
         c_addr=ptr_offset(pptr,off_base)
!	 call print_c_ptr(c_addr) !debug
         call c_f_pointer(c_addr,ptr_base,shape=[trank]); ptr_base(1:trank)=base(1:trank); nullify(ptr_base)
        endif
        c_addr=ptr_offset(pptr,off_sclr)
!	call print_c_ptr(c_addr) !debug
        call c_f_pointer(c_addr,ptr_sclr); ptr_sclr=sclr; nullify(ptr_sclr)
        if(trank.gt.0.and.elems_count.gt.0) then
         select case(dtk)
         case('r4')
          c_addr=ptr_offset(pptr,off_elems_r4)
!	  call print_c_ptr(c_addr) !debug
          call c_f_pointer(c_addr,ptr_elems_r4,shape=[elems_count])
          do s0=0,elems_count-1; ptr_elems_r4(1+s0)=tens%data_real4(s0); enddo; nullify(ptr_elems_r4)
         case('r8')
          c_addr=ptr_offset(pptr,off_elems_r8)
!	  call print_c_ptr(c_addr) !debug
          call c_f_pointer(c_addr,ptr_elems_r8,shape=[elems_count])
          do s0=0,elems_count-1; ptr_elems_r8(1+s0)=tens%data_real8(s0); enddo; nullify(ptr_elems_r8)
         case('c8')
          c_addr=ptr_offset(pptr,off_elems_c8)
!	  call print_c_ptr(c_addr) !debug
          call c_f_pointer(c_addr,ptr_elems_c8,shape=[elems_count])
          do s0=0,elems_count-1; ptr_elems_c8(1+s0)=tens%data_cmplx8(s0); enddo; nullify(ptr_elems_c8)
         end select
        elseif(trank.eq.0.and.elems_count.eq.1) then !scalar
         select case(dtk)
         case('r4')
          c_addr=ptr_offset(pptr,off_elems_r4); call c_f_pointer(c_addr,ptr_elems_r4,shape=[1])
          ptr_elems_r4(1)=real(cmplx8_to_real8(sclr),4); nullify(ptr_elems_r4)
         case('r8')
          c_addr=ptr_offset(pptr,off_elems_r8); call c_f_pointer(c_addr,ptr_elems_r8,shape=[1])
          ptr_elems_r8(1)=cmplx8_to_real8(sclr); nullify(ptr_elems_r8)
         case('c8')
          c_addr=ptr_offset(pptr,off_elems_c8); call c_f_pointer(c_addr,ptr_elems_c8,shape=[1])
          ptr_elems_c8(1)=sclr; nullify(ptr_elems_c8)
         end select
        else
         ierr=7; return
        endif
!	write(jo_cp,'("#DEBUG(c_process::tens_blck_pack): packet created (size/entry): ",i12,1x,i7)') packet_size,entry_num !debug
        return
        end subroutine tens_blck_pack
!-----------------------------------------------------------------------------
        subroutine tens_blck_unpack(tens,pptr,ierr)
!This subroutine creates an instance of tensor_block_t (F) <tens> by unpacking
!a tensor block packet pointed to by a C pointer <pptr>:
! PACKET (Host Argument Buffer entry) --> TENSOR_BLOCK_T
!Note that the packet will still reside in the Host argument buffer after returning from this subroutine.
!INPUT:
! - tens - an allocated (uninitialized) instance of tensor_block_t;
! - pptr - C pointer to the tensor block packet located in the Host argument buffer;
!OUTPUT:
! - tens - a filled instance of tensor_block_t;
! - ierr - error code (0:success).
!NOTES:
! - This subroutine does NOT free the corresponding Host argument buffer entry!
! - Packet structure is specified in <tens_blck_pack>.
        implicit none
        type(tensor_block_t), intent(inout):: tens
        type(C_PTR), intent(in):: pptr
        integer, intent(inout):: ierr
        integer i,j,base(1:max_tensor_rank)
        integer(C_SIZE_T) s0
        integer(C_INT) err_code
        type(C_PTR) c_addr
        integer(C_SIZE_T) off_packet_size; integer(C_SIZE_T), pointer:: ptr_packet_size=>NULL(); integer(C_SIZE_T) packet_size
        integer(C_SIZE_T) off_data_kind; integer(C_INT), pointer:: ptr_data_kind=>NULL(); integer(C_INT) data_kind
        integer(C_SIZE_T) off_elems_count; integer(C_SIZE_T), pointer:: ptr_elems_count=>NULL(); integer(C_SIZE_T) elems_count
        integer(C_SIZE_T) off_trank; integer(C_INT), pointer:: ptr_trank=>NULL(); integer(C_INT) trank
        integer(C_SIZE_T) off_dims; integer(C_INT), pointer:: ptr_dims(:)=>NULL()
        integer(C_SIZE_T) off_divs; integer(C_INT), pointer:: ptr_divs(:)=>NULL()
        integer(C_SIZE_T) off_grps; integer(C_INT), pointer:: ptr_grps(:)=>NULL()
        integer(C_SIZE_T) off_base; integer(C_INT), pointer:: ptr_base(:)=>NULL()
        integer(C_SIZE_T) off_prmn
        integer(C_SIZE_T) off_sclr; complex(8), pointer:: ptr_sclr=>NULL()
        integer(C_SIZE_T) off_elems_r4; real(4), pointer:: ptr_elems_r4(:)=>NULL(); real(4) elems_r4
        integer(C_SIZE_T) off_elems_r8; real(8), pointer:: ptr_elems_r8(:)=>NULL(); real(8) elems_r8
        integer(C_SIZE_T) off_elems_c8; complex(8), pointer:: ptr_elems_c8(:)=>NULL(); complex(8) elems_c8
        logical res

        ierr=0; err_code=0
!	write(jo_cp,'("#DEBUG(c_process::tens_blck_upack): Unpacking a tensor block packet:")') !debug
!Read the tensor block packet and form an instance of tensor_block_t:
        off_packet_size=0; c_addr=ptr_offset(pptr,off_packet_size)
        call c_f_pointer(c_addr,ptr_packet_size); packet_size=ptr_packet_size; nullify(ptr_packet_size)
        off_data_kind=off_packet_size+sizeof(packet_size); c_addr=ptr_offset(pptr,off_data_kind)
        call c_f_pointer(c_addr,ptr_data_kind); data_kind=ptr_data_kind; nullify(ptr_data_kind)
        off_elems_count=off_data_kind+sizeof(data_kind); c_addr=ptr_offset(pptr,off_elems_count)
        call c_f_pointer(c_addr,ptr_elems_count); elems_count=ptr_elems_count; nullify(ptr_elems_count)
        off_trank=off_elems_count+sizeof(elems_count); c_addr=ptr_offset(pptr,off_trank)
        call c_f_pointer(c_addr,ptr_trank); trank=ptr_trank; nullify(ptr_trank)
        tens%tensor_shape%num_dim=trank; if(trank.lt.0) then; ierr=1; return; endif
        tens%tensor_block_size=elems_count; if(elems_count.le.0.or.(trank.eq.0.and.elems_count.ne.1)) then; ierr=2; return; endif
        if(tensor_block_alloc(tens,'sp',ierr)) then
         if(ierr.ne.0) then; ierr=3; return; endif
         deallocate(tens%tensor_shape%dim_extent,STAT=ierr); if(ierr.ne.0) then; ierr=4; return; endif
         deallocate(tens%tensor_shape%dim_divider,STAT=ierr); if(ierr.ne.0) then; ierr=5; return; endif
         deallocate(tens%tensor_shape%dim_group,STAT=ierr); if(ierr.ne.0) then; ierr=6; return; endif
         res=tensor_block_alloc(tens,'sp',ierr,.false.); if(ierr.ne.0) then; ierr=7; return; endif
        else
         if(ierr.ne.0) then; ierr=8; return; endif
         nullify(tens%tensor_shape%dim_extent)
         nullify(tens%tensor_shape%dim_divider)
         nullify(tens%tensor_shape%dim_group)
        endif        
        if(trank.gt.0) then
         allocate(tens%tensor_shape%dim_extent(1:trank),STAT=ierr); if(ierr.ne.0) then; ierr=9; return; endif
         allocate(tens%tensor_shape%dim_divider(1:trank),STAT=ierr); if(ierr.ne.0) then; ierr=10; return; endif
         allocate(tens%tensor_shape%dim_group(1:trank),STAT=ierr); if(ierr.ne.0) then; ierr=11; return; endif
         res=tensor_block_alloc(tens,'sp',j,.true.); if(j.ne.0) then; ierr=12; return; endif
        endif
        select case(data_kind)
        case(R4)
         if(associated(tens%data_real4)) then
          if(tensor_block_alloc(tens,'r4',j)) then
           if(j.ne.0) then; ierr=13; return; endif
           deallocate(tens%data_real4,STAT=j); if(j.ne.0) then; ierr=14; return; endif
          else
           if(j.ne.0) then; ierr=15; return; endif
           nullify(tens%data_real4)
          endif
          res=tensor_block_alloc(tens,'r4',j,.false.); if(j.ne.0) then; ierr=16; return; endif
         endif
         if(trank.gt.0) then
          allocate(tens%data_real4(0:elems_count-1),STAT=ierr); if(ierr.ne.0) then; ierr=17; return; endif
          res=tensor_block_alloc(tens,'r4',j,.true.); if(j.ne.0) then; ierr=18; return; endif
         endif
        case(R8)
         if(associated(tens%data_real8)) then
          if(tensor_block_alloc(tens,'r8',j)) then
           if(j.ne.0) then; ierr=19; return; endif
           deallocate(tens%data_real8,STAT=j); if(j.ne.0) then; ierr=20; return; endif
          else
           if(j.ne.0) then; ierr=21; return; endif
           nullify(tens%data_real8)
          endif
          res=tensor_block_alloc(tens,'r8',j,.false.); if(j.ne.0) then; ierr=22; return; endif
         endif
         if(trank.gt.0) then
          allocate(tens%data_real8(0:elems_count-1),STAT=ierr); if(ierr.ne.0) then; ierr=23; return; endif
          res=tensor_block_alloc(tens,'r8',j,.true.); if(j.ne.0) then; ierr=24; return; endif
         endif
        case(C8)
         if(associated(tens%data_cmplx8)) then
          if(tensor_block_alloc(tens,'c8',j)) then
           if(j.ne.0) then; ierr=25; return; endif
           deallocate(tens%data_cmplx8,STAT=j); if(j.ne.0) then; ierr=26; return; endif
          else
           if(j.ne.0) then; ierr=27; return; endif
           nullify(tens%data_cmplx8)
          endif
          res=tensor_block_alloc(tens,'c8',j,.false.); if(j.ne.0) then; ierr=28; return; endif
         endif
         if(trank.gt.0) then
          allocate(tens%data_cmplx8(0:elems_count-1),STAT=ierr); if(ierr.ne.0) then; ierr=29; return; endif
          res=tensor_block_alloc(tens,'c8',j,.true.); if(j.ne.0) then; ierr=30; return; endif
         endif
        case default
         ierr=31; return
        end select
        off_dims=off_trank+sizeof(trank)
        if(trank.gt.0) then
         c_addr=ptr_offset(pptr,off_dims); call c_f_pointer(c_addr,ptr_dims,shape=[trank])
         tens%tensor_shape%dim_extent(1:trank)=ptr_dims(1:trank); nullify(ptr_dims)
        endif
        off_divs=off_dims+sizeof(err_code)*trank
        if(trank.gt.0) then
         c_addr=ptr_offset(pptr,off_divs); call c_f_pointer(c_addr,ptr_divs,shape=[trank])
         tens%tensor_shape%dim_divider(1:trank)=ptr_divs(1:trank); nullify(ptr_divs)
        endif
        off_grps=off_divs+sizeof(err_code)*trank
        if(trank.gt.0) then
         c_addr=ptr_offset(pptr,off_grps); call c_f_pointer(c_addr,ptr_grps,shape=[trank])
         tens%tensor_shape%dim_group(1:trank)=ptr_grps(1:trank); nullify(ptr_grps)
        endif
        off_base=off_grps+sizeof(err_code)*trank
        if(trank.gt.0) then
         c_addr=ptr_offset(pptr,off_base); call c_f_pointer(c_addr,ptr_base,shape=[trank])
         base(1:trank)=ptr_base(1:trank); nullify(ptr_base)
        endif
        off_prmn=off_base+sizeof(err_code)*trank
        off_sclr=off_prmn+sizeof(err_code)*trank; c_addr=ptr_offset(pptr,off_sclr)
        call c_f_pointer(c_addr,ptr_sclr); tens%scalar_value=ptr_sclr; nullify(ptr_sclr)
        select case(data_kind)
        case(R4)
         off_elems_r4=off_sclr+sizeof(elems_c8)
         if(trank.gt.0) then
          c_addr=ptr_offset(pptr,off_elems_r4)
          if(off_elems_r4+elems_count*sizeof(elems_r4).ne.packet_size) then; ierr=32; return; endif
          call c_f_pointer(c_addr,ptr_elems_r4,shape=[elems_count])
          do s0=0,elems_count-1; tens%data_real4(s0)=ptr_elems_r4(1+s0); enddo
          nullify(ptr_elems_r4)
         endif
        case(R8)
         off_elems_r8=off_sclr+sizeof(elems_c8)
         if(trank.gt.0) then
          c_addr=ptr_offset(pptr,off_elems_r8)
          if(off_elems_r8+elems_count*sizeof(elems_r8).ne.packet_size) then; ierr=33; return; endif
          call c_f_pointer(c_addr,ptr_elems_r8,shape=[elems_count])
          do s0=0,elems_count-1; tens%data_real8(s0)=ptr_elems_r8(1+s0); enddo
          nullify(ptr_elems_r8)
         endif
        case(C8)
         off_elems_c8=off_sclr+sizeof(elems_c8)
         if(trank.gt.0) then
          c_addr=ptr_offset(pptr,off_elems_c8)
          if(off_elems_c8+elems_count*sizeof(elems_c8).ne.packet_size) then; ierr=34; return; endif
          call c_f_pointer(c_addr,ptr_elems_c8,shape=[elems_count])
          do s0=0,elems_count-1; tens%data_cmplx8(s0)=ptr_elems_c8(1+s0); enddo
          nullify(ptr_elems_c8)
         endif
        end select
!DEBUG begin:
!	write(jo_cp,'(" Packet size  : ",i9,1x,i9)') off_packet_size,packet_size
!	write(jo_cp,'(" Data kind    : ",i9,1x,i9)') off_data_kind,data_kind
!	write(jo_cp,'(" Element count: ",i9,1x,i9)') off_elems_count,elems_count
!	write(jo_cp,'(" Tensor rank  : ",i9,1x,i9)') off_trank,trank
!	write(jo_cp,'(" Dims         :",i9,3x,32(1x,i4))') off_dims,tens%tensor_shape%dim_extent(1:trank)
!	write(jo_cp,'(" Divs         :",i9,3x,32(1x,i4))') off_divs,tens%tensor_shape%dim_divider(1:trank)
!	write(jo_cp,'(" Grps         :",i9,3x,32(1x,i4))') off_grps,tens%tensor_shape%dim_group(1:trank)
!	write(jo_cp,'(" Base         :",i9,3x,32(1x,i4))') off_base,base(1:trank)
!	write(jo_cp,'(" Scalar       : ",i9,3x,D25.14,1x,D25.14)') off_sclr,tens%scalar_value
!	select case(data_kind)
!	case(R4); write(jo_cp,'(" Elements(r4) : ",i9,1x,i9)') off_elems_r4
!	case(R8); write(jo_cp,'(" Elements(r8) : ",i9,1x,i9)') off_elems_r8
!	case(C8); write(jo_cp,'(" Elements(c8) : ",i9,1x,i9)') off_elems_c8
!	end select
!	write(jo_cp,'("#DEBUG(c_process::tens_blck_unpack): packet unpacked (size): ",i12)') packet_size !debug
!DEBUG end.
        return
        end subroutine tens_blck_unpack
!---------------------------------------------------------------
        subroutine tens_blck_assoc(pptr,ierr,tens,ctens,gpu_num)
!Based on the packet located at <pptr>, this subroutine fills in an instance of
!either tensBlck_t (C/C++) <ctens> or tensor_block_t (Fortran) <tens>:
! {tensor_block_t|tensBlck_t} => PACKET (Host Argument Buffer entry)
!Note that the pointer fields of {tensor_block_t|tensBlck_t} will simply point
!to the corresponding locations in the Host Argument Buffer. Hence the corresponding entry
!of the Host Argument Buffer shall not be freed during the lifetime of {tensor_block_t|tensBlck_t}.
!INPUT:
! - pptr - C pointer to a tensor block packet located in the Host argument buffer;
! - tens - an instance of tensor_block_t (F);
! - ctens - an instance of tensBlck_t (C);
! - gpu_num - GPU# on which the tensor block will be used (-1 means Host);
!OUTPUT:
! - ctens - a filled instance of tensBlck_t;
! - tens - a filled instance of tensor_block_t;
! - ierr - error code (0:success).
!NOTES:
! - Either <tens> or a pair {<ctens>,<gpu_num>} must be supplied, and they are mutually exclusive.
! - This subroutine sets all the fields of tensBlck_t, like elems_h, elems_d, const_args_entry, etc.
!   Apparently it also sets all the fields of tensor_block_t.
! - The opposite function <tens_blck_dissoc> can be used only with tensBlck_t (C).
!   For tensor_block_t (Fortran) one should invoke <tensor_block_destroy>.
        implicit none
        type(C_PTR), intent(in):: pptr
        integer, intent(inout):: ierr
        type(tensor_block_t), optional, intent(inout):: tens
        type(C_PTR), optional, intent(inout):: ctens
        integer(C_INT), optional, intent(in):: gpu_num
        integer i
        integer(C_SIZE_T) s0
        integer(C_INT) entry_gpu,entry_const,err_code !err_code must be C_INT=integer(4)
        type(C_PTR) c_addr,addr_dims,addr_divs,addr_grps,addr_base,addr_prmn,addr_host,addr_gpu
        integer(C_SIZE_T) off_packet_size; integer(C_SIZE_T), pointer:: ptr_packet_size=>NULL(); integer(C_SIZE_T) packet_size
        integer(C_SIZE_T) off_data_kind; integer(C_INT), pointer:: ptr_data_kind=>NULL(); integer(C_INT) data_kind
        integer(C_SIZE_T) off_elems_count; integer(C_SIZE_T), pointer:: ptr_elems_count=>NULL(); integer(C_SIZE_T) elems_count
        integer(C_SIZE_T) off_trank; integer(C_INT), pointer:: ptr_trank=>NULL(); integer(C_INT) trank
        integer(C_SIZE_T) off_dims
        integer(C_SIZE_T) off_divs
        integer(C_SIZE_T) off_grps
        integer(C_SIZE_T) off_base
        integer(C_SIZE_T) off_prmn
        integer(C_SIZE_T) off_sclr; complex(8), pointer:: ptr_sclr=>NULL()
        integer(C_SIZE_T) off_elems_r4; real(4), pointer:: ptr_elems_r4(:)=>NULL(); real(4) elems_r4
        integer(C_SIZE_T) off_elems_r8; real(8), pointer:: ptr_elems_r8(:)=>NULL(); real(8) elems_r8
        integer(C_SIZE_T) off_elems_c8; complex(8), pointer:: ptr_elems_c8(:)=>NULL(); complex(8) elems_c8

        ierr=0
!Extract all components from the tensor block packet:
        off_packet_size=0; c_addr=ptr_offset(pptr,off_packet_size)
        call c_f_pointer(c_addr,ptr_packet_size); packet_size=ptr_packet_size; nullify(ptr_packet_size)
        off_data_kind=off_packet_size+sizeof(packet_size); c_addr=ptr_offset(pptr,off_data_kind)
        call c_f_pointer(c_addr,ptr_data_kind); data_kind=ptr_data_kind; nullify(ptr_data_kind)
        off_elems_count=off_data_kind+sizeof(data_kind); c_addr=ptr_offset(pptr,off_elems_count)
        call c_f_pointer(c_addr,ptr_elems_count); elems_count=ptr_elems_count; nullify(ptr_elems_count)
        off_trank=off_elems_count+sizeof(elems_count); c_addr=ptr_offset(pptr,off_trank)
        call c_f_pointer(c_addr,ptr_trank); trank=ptr_trank; nullify(ptr_trank)
        if(trank.lt.0.or.(trank.eq.0.and.elems_count.ne.1).or.(trank.gt.0.and.elems_count.lt.1)) then; ierr=1; return; endif
        off_dims=off_trank+sizeof(trank); addr_dims=ptr_offset(pptr,off_dims)
        off_divs=off_dims+sizeof(err_code)*trank; addr_divs=ptr_offset(pptr,off_divs)
        off_grps=off_divs+sizeof(err_code)*trank; addr_grps=ptr_offset(pptr,off_grps)
        off_base=off_grps+sizeof(err_code)*trank; addr_base=ptr_offset(pptr,off_base)
        off_prmn=off_base+sizeof(err_code)*trank; addr_prmn=ptr_offset(pptr,off_prmn)
        off_sclr=off_prmn+sizeof(err_code)*trank; c_addr=ptr_offset(pptr,off_sclr)
        call c_f_pointer(c_addr,ptr_sclr); elems_c8=ptr_sclr; nullify(ptr_sclr)
        select case(data_kind)
        case(R4)
         off_elems_r4=off_sclr+sizeof(elems_c8); addr_host=ptr_offset(pptr,off_elems_r4); s0=sizeof(elems_r4)*elems_count
        case(R8)
         off_elems_r8=off_sclr+sizeof(elems_c8); addr_host=ptr_offset(pptr,off_elems_r8); s0=sizeof(elems_r8)*elems_count
        case(C8)
         off_elems_c8=off_sclr+sizeof(elems_c8); addr_host=ptr_offset(pptr,off_elems_c8); s0=sizeof(elems_c8)*elems_count
        case default
         ierr=2; return; !invalid data kind
        end select
!Associate:
        if(present(ctens).and.(.not.present(tens)).and.present(gpu_num)) then
#ifndef NO_GPU
 !Create an instance of tensBlck_t (C):
         err_code=get_buf_entry_gpu(gpu_num,s0,addr_gpu,entry_gpu); if(err_code.ne.0.or.entry_gpu.lt.0) then; ierr=3; return; endif
         err_code=const_args_entry_get(gpu_num,entry_const); if(err_code.ne.0.or.entry_const.lt.0) then; ierr=4; return; endif
         err_code=tensBlck_create(ctens); if(err_code.ne.0) then; ctens=C_NULL_PTR; ierr=5; return; endif
         err_code=tensBlck_construct(ctens,DEV_NVIDIA_GPU,gpu_num,data_kind,trank,addr_dims,addr_divs,addr_grps,addr_prmn,addr_host,addr_gpu,entry_gpu,entry_const)
         if(err_code.ne.0) then; ierr=6; return; endif
#else
         write(jo_cp,'("#FATAL(c_process::tens_blck_assoc): attempt to initialize a GPU-resident tensor block in GPU-free code compilation!")') !trap
         ierr=-1; return !attempt to initialize tensBlck_t in a GPU-free code compilation
#endif
        elseif(present(tens).and.(.not.present(ctens)).and.(.not.present(gpu_num))) then
 !Create an instance of tensor_block_t (Fortran):
         call tensor_block_destroy(tens,ierr); if(ierr.ne.0) then; ierr=7; return; endif
         tens%tensor_block_size=elems_count
         tens%tensor_shape%num_dim=trank
         if(trank.gt.0) then
          call c_f_pointer(addr_dims,tens%tensor_shape%dim_extent,shape=[trank])
          call c_f_pointer(addr_divs,tens%tensor_shape%dim_divider,shape=[trank])
          call c_f_pointer(addr_grps,tens%tensor_shape%dim_group,shape=[trank])
         endif
         tens%scalar_value=elems_c8
         select case(data_kind)
         case(R4)
          call c_f_pointer(addr_host,ptr_elems_r4,shape=[elems_count])
          tens%data_real4(0:)=>ptr_elems_r4; nullify(ptr_elems_r4) !Hope this is portable`
         case(R8)
          call c_f_pointer(addr_host,ptr_elems_r8,shape=[elems_count])
          tens%data_real8(0:)=>ptr_elems_r8; nullify(ptr_elems_r8) !Hope this is portable`
         case(C8)
          call c_f_pointer(addr_host,ptr_elems_c8,shape=[elems_count])
          tens%data_cmplx8(0:)=>ptr_elems_c8; nullify(ptr_elems_c8) !Hope this is portable`
         end select
        else
         ierr=8 !both <ctens> and <tens> cannot be absent/present simultaneously
        endif
        return
        end subroutine tens_blck_assoc
!--------------------------------------------------------------------------
        subroutine tens_blck_dissoc(ctens,ierr)
!This subroutine dissociates an object of type tensBlck_t (C tensor block) from
!GPU argument buffers (Global & Constant memory), frees the corresponding
!GPU argument buffer entries, and destroys the tensBlck_t object.
!Note that the corresponding Host argument buffer entry (with its content) is not freed!
!INPUT:
! - ctens - tensBlck_t;
!OUTPUT:
! - ierr - error code (0:success).
        implicit none
        type(C_PTR), intent(inout):: ctens
        integer, intent(inout):: ierr
        integer(C_INT) gn,dev_kind,entry_gpu,entry_const,data_kind,there,err_code

        ierr=0
#ifndef NO_GPU
        gn=tensBlck_acc_id(ctens,dev_kind,entry_gpu,entry_const,data_kind,there)
        if(dev_kind.eq.DEV_NVIDIA_GPU.and.gn.ge.0.and.entry_gpu.ge.0.and.entry_const.ge.0) then
         err_code=free_buf_entry_gpu(gn,entry_gpu); if(err_code.ne.0) ierr=ierr+10
         err_code=const_args_entry_free(gn,entry_const); if(err_code.ne.0) ierr=ierr+100
        else
         ierr=ierr+1
        endif
        err_code=tensBlck_destroy(ctens); if(err_code.ne.0) ierr=ierr+1000
#else
        write(jo_cp,'("#FATAL(c_process::tens_blck_dissoc): attempt to dissociate a GPU-resident tensor block in GPU-free code compilation!")') !trap
        ierr=-1; return !attempt to initialize tensBlck_t in a GPU-free code compilation
#endif
        return
        end subroutine tens_blck_dissoc
!------------------------------------------------------------------------------------------------------------
        subroutine c_proc_test(ierr)
!This subroutine tests the basic computing functionality of a C-process by running some tensor algebra tests.
        implicit none
        integer, intent(inout):: ierr
        integer, parameter:: test_args_lim=15
        integer i,j,k,l,m,n,ctrl
        integer(8) diffc
        integer(C_INT) i1,err_code,ext_beg(1:max_tensor_rank),cptrn(1:max_tensor_rank*2),gpu_id
        integer(C_INT) o2n(0:max_tensor_rank),n2o(0:max_tensor_rank),ngt(0:max_tensor_rank)
        real(8) norm2
        real(C_FLOAT) tm0,tm1,tm2
        logical cmp
        integer(C_INT):: entry_num(0:test_args_lim)=-1
        type(C_PTR):: entry_ptr(0:test_args_lim)=C_NULL_PTR
        integer(C_SIZE_T):: s0,pack_size(0:test_args_lim)=0
        type(C_PTR):: cuda_task(1:3)=C_NULL_PTR
        type(C_PTR):: ctens(0:test_args_lim)=C_NULL_PTR
        type(tensor_block_t) ftens(0:test_args_lim)
        character(256) shape0,shape1,shape2

        ierr=0; write(jo_cp,'("#MSG(c_process::c_proc_test): Testing C-process functionality ... ")',advance='no')
!TEST 0: Random contraction patterns:
!        write(jo_cp,*)''
!        do i=1,16
!         call contr_pattern_rnd(8,100000000,shape0,shape1,shape2,cptrn,ierr)
!         if(ierr.ne.0) then; write(jo_cp,*)'ERROR(c_process::c_proc_test): contr_pattern_rnd error ',ierr; return; endif
!         l=index(shape0,')'); k=index(shape1,')'); j=index(shape2,')')
!         call printl(jo_cp,shape0(1:l)//'+='//shape1(1:k)//'*'//shape2(1:j)) !debug
!         m=tensor_shape_rank(shape1(1:k),ierr)+tensor_shape_rank(shape2(1:j),ierr)
!         write(jo_cp,'("CPTRN:",64(1x,i2))') cptrn(1:m) !debug
!        enddo
!        return
!TEST 1 (CPU: tensor_block_contract):
        write(jo_cp,'("1 ")',advance='no')
        call tensor_block_create('(20,30,20,30)','r8',ftens(0),ierr,val_r8=0d0); if(ierr.ne.0) then; ierr=1; goto 999; endif
        call tensor_block_create('(15,20,25,30)','r8',ftens(1),ierr,val_r8=1d0); if(ierr.ne.0) then; ierr=2; goto 999; endif
        call tensor_block_create('(30,20,15,25)','r8',ftens(2),ierr,val_r8=1d0); if(ierr.ne.0) then; ierr=3; goto 999; endif
        call tensor_block_contract((/-3,3,-4,2,4,1,-1,-3/),ftens(1),ftens(2),ftens(0),ierr,'r8'); if(ierr.ne.0) then; ierr=4; goto 999; endif
        norm2=tensor_block_norm2(ftens(0),ierr,'r8'); if(ierr.ne.0) then; ierr=5; goto 999; endif
        if(int(norm2,8).ne.ftens(0)%tensor_block_size*(15_8*25_8)**2) then; ierr=6; goto 999; endif
        call tensor_block_destroy(ftens(0),ierr); if(ierr.ne.0) then; ierr=7; goto 999; endif
        call tensor_block_destroy(ftens(1),ierr); if(ierr.ne.0) then; ierr=8; goto 999; endif
        call tensor_block_destroy(ftens(2),ierr); if(ierr.ne.0) then; ierr=9; goto 999; endif
#ifndef NO_GPU
        gpu_id=gpu_busy_least(); call cudaSetDevice(gpu_id,err_code); if(err_code.ne.0) then; ierr=10; goto 999; endif
        call gpu_set_event_policy(EVENTS_ON)
!TEST 2 (GPU: matrix multiplication):
        write(jo_cp,'("2 ")',advance='no')
        n=2; ext_beg(1:n)=1
        call tensor_block_create('(345,697)','r8',ftens(0),ierr,val_r8=0d0); if(ierr.ne.0) then; ierr=11; goto 999; endif
        call tensor_block_create('(37,345)','r8',ftens(1),ierr); if(ierr.ne.0) then; ierr=12; goto 999; endif
        call tensor_block_create('(37,697)','r8',ftens(2),ierr); if(ierr.ne.0) then; ierr=13; goto 999; endif
        call tensor_block_create('(345,697)','r4',ftens(3),ierr,val_r4=0.0); if(ierr.ne.0) then; ierr=14; goto 999; endif
        call tensor_block_contract((/-1,1,-1,2/),ftens(1),ftens(2),ftens(0),ierr,'r8'); if(ierr.ne.0) then; ierr=15; goto 999; endif
        call tensor_block_sync(ftens(0),'r8',ierr,'r4'); if(ierr.ne.0) then; ierr=16; goto 999; endif
        call tensor_block_sync(ftens(1),'r8',ierr,'r4'); if(ierr.ne.0) then; ierr=17; goto 999; endif
        call tensor_block_sync(ftens(2),'r8',ierr,'r4'); if(ierr.ne.0) then; ierr=18; goto 999; endif
!        call tensor_block_print(jo_cp,'Matrix0:',ext_beg,ftens(0),ierr,'r4') !debug
        err_code=gpu_matrix_multiply_tn_r4(345_CZ,697_CZ,37_CZ,ftens(1)%data_real4,ftens(2)%data_real4,ftens(3)%data_real4); if(err_code.ne.0) then; write(jo_cp,*)'GPU ERROR ',err_code; ierr=19; goto 999; endif
        i1=gpu_get_error_count(); if(i1.ne.0) then; write(jo_cp,*)'GPU error count = ',i1; ierr=20; goto 999; endif
!        call tensor_block_print(jo_cp,'Matrix3:',ext_beg,ftens(3),ierr,'r4') !debug
        cmp=tensor_block_cmp(ftens(0),ftens(3),ierr,'r4',cmp_thresh=1d-5); if(ierr.ne.0) then; ierr=21; goto 999; endif
        if(.not.cmp) then; ierr=22; goto 999; endif
        call tensor_block_destroy(ftens(3),ierr); if(ierr.ne.0) then; ierr=23; goto 999; endif
        call tensor_block_destroy(ftens(2),ierr); if(ierr.ne.0) then; ierr=24; goto 999; endif
        call tensor_block_destroy(ftens(1),ierr); if(ierr.ne.0) then; ierr=25; goto 999; endif
        call tensor_block_destroy(ftens(0),ierr); if(ierr.ne.0) then; ierr=26; goto 999; endif
!TEST 3 (GPU: tensor packing/unpacking/association, gpu_tensor_transpose):
        write(jo_cp,'("3 ")',advance='no')
        n=4; ext_beg(1:n)=1; o2n(0:n)=(/+1,(j,j=1,n)/); ctrl=1
        do i=1,ifcl(n)
         call trng(ctrl,n,o2n,ngt); !write(jo_cp,'(32(1x,i2))') o2n(0:n)
         call permutation_converter(.false.,n,n2o,o2n); !write(jo_cp,'(32(1x,i2))') n2o(0:n)
         call tensor_block_create('(17,27,33,44)','r8',ftens(0),ierr); if(ierr.ne.0) then; ierr=27; goto 999; endif
         call tensor_block_create('(44,33,27,17)','r8',ftens(1),ierr,val_r8=0d0); if(ierr.ne.0) then; ierr=28; goto 999; endif
         tm0=thread_wtime()
         call tensor_block_copy(ftens(0),ftens(1),ierr,o2n); if(ierr.ne.0) then; ierr=29; goto 999; endif
!         write(jo_cp,'("#DEBUG(c_process::c_proc_test): CPU tensor transpose time = ",F10.4)') thread_wtime()-tm0 !debug
         call tensor_block_sync(ftens(0),'r8',ierr,'r4'); if(ierr.ne.0) then; ierr=30; goto 999; endif
         call tensor_block_sync(ftens(1),'r8',ierr,'r4'); if(ierr.ne.0) then; ierr=31; goto 999; endif
         call tens_blck_pack(ftens(1),'r4',pack_size(1),entry_ptr(1),entry_num(1),ierr); if(ierr.ne.0) then; ierr=32; goto 999; endif
         call tens_blck_unpack(ftens(2),entry_ptr(1),ierr); if(ierr.ne.0) then; ierr=33; goto 999; endif
         cmp=tensor_block_cmp(ftens(1),ftens(2),ierr,'r4',cmp_thresh=1d-6); if(ierr.ne.0) then; ierr=34; goto 999; endif
         if(.not.cmp) then; ierr=35; goto 999; endif
         call tensor_block_destroy(ftens(2),ierr); if(ierr.ne.0) then; ierr=36; goto 999; endif
         call tensor_block_create('(17,27,33,44)','r4',ftens(2),ierr,val_r4=0.0); if(ierr.ne.0) then; ierr=37; goto 999; endif
         call tens_blck_pack(ftens(2),'r4',pack_size(2),entry_ptr(2),entry_num(2),ierr); if(ierr.ne.0) then; ierr=38; goto 999; endif
         call tens_blck_assoc(entry_ptr(1),ierr,ctens=ctens(1),gpu_num=gpu_id); if(ierr.ne.0) then; ierr=39; goto 999; endif
         call tens_blck_assoc(entry_ptr(2),ierr,ctens=ctens(2),gpu_num=gpu_id); if(ierr.ne.0) then; ierr=40; goto 999; endif
         tm0=thread_wtime()
         err_code=gpu_tensor_block_copy_dlf(n2o,ctens(1),ctens(2)); if(err_code.ne.0) then; write(jo_cp,*)'GPU ERROR ',err_code; ierr=41; goto 999; endif
!         write(jo_cp,'("#DEBUG(c_process::c_proc_test): GPU tensor transpose time = ",F10.4)') thread_wtime()-tm0 !debug
         i1=gpu_get_error_count(); if(i1.ne.0) then; write(jo_cp,*)'GPU error count = ',i1; ierr=42; goto 999; endif
!         call print_gpu_debug_dump(jo_cp) !debug
         if(tensor_block_norm2(ftens(2),ierr,'r4').ne.0d0) then; ierr=43; goto 999; endif
         call tens_blck_unpack(ftens(2),entry_ptr(2),ierr); if(ierr.ne.0) then; ierr=44; goto 999; endif
!         write(jo_cp,*) tensor_block_norm2(ftens(0),ierr,'r4'),tensor_block_norm2(ftens(2),ierr,'r4') !debug
!         call tensor_block_print(jo_cp,'Tensor0:',ext_beg,ftens(0),ierr,'r4') !debug
!         call tensor_block_print(jo_cp,'Tensor1:',ext_beg,ftens(1),ierr,'r4') !debug
!         call tensor_block_print(jo_cp,'Tensor2:',ext_beg,ftens(2),ierr,'r4') !debug
         cmp=tensor_block_cmp(ftens(0),ftens(2),ierr,'r4',cmp_thresh=1d-5,diff_count=diffc); if(ierr.ne.0) then; ierr=45; goto 999; endif
         if(.not.cmp) then; write(jo_cp,*)'DIFF COUNT = ',diffc; ierr=46; goto 999; endif
         call tens_blck_dissoc(ctens(2),ierr); if(ierr.ne.0) then; ierr=47; goto 999; endif
         call tens_blck_dissoc(ctens(1),ierr); if(ierr.ne.0) then; ierr=48; goto 999; endif
         err_code=free_buf_entry_host(entry_num(2)); if(err_code.ne.0) then; ierr=49; goto 999; endif
         err_code=free_buf_entry_host(entry_num(1)); if(err_code.ne.0) then; ierr=50; goto 999; endif
         call tensor_block_destroy(ftens(2),ierr); if(ierr.ne.0) then; ierr=51; goto 999; endif
         call tensor_block_destroy(ftens(1),ierr); if(ierr.ne.0) then; ierr=52; goto 999; endif
         call tensor_block_destroy(ftens(0),ierr); if(ierr.ne.0) then; ierr=53; goto 999; endif
        enddo
!TEST 4 (GPU: gpu_tensor_contraction):
        write(jo_cp,'("4 ")',advance='no')
        call tensor_block_create('(3,4,31,2,37,8)','r4',ftens(0),ierr,val_r4=0.0); if(ierr.ne.0) then; ierr=54; goto 999; endif
        call tensor_block_create('(2,31,7,8,11,29)','r4',ftens(1),ierr,val_r4=1.0); if(ierr.ne.0) then; ierr=55; goto 999; endif
        call tensor_block_create('(11,37,4,7,29,3)','r4',ftens(2),ierr,val_r4=1.0); if(ierr.ne.0) then; ierr=56; goto 999; endif
        call tensor_block_create('(3,4,31,2,37,8)','r4',ftens(3),ierr,val_r4=0.0); if(ierr.ne.0) then; ierr=57; goto 999; endif
        call tensor_block_create('(2,31,7,8,11,29)','r4',ftens(4),ierr,val_r4=0.0); if(ierr.ne.0) then; ierr=58; goto 999; endif
        call tensor_block_create('(11,37,4,7,29,3)','r4',ftens(5),ierr,val_r4=0.0); if(ierr.ne.0) then; ierr=59; goto 999; endif
        call tensor_block_create('(3,4,31,2,37,8)','r4',ftens(6),ierr,val_r4=0.0); if(ierr.ne.0) then; ierr=60; goto 999; endif
        cptrn(1:12)=(/4,3,-4,6,-1,-5,-5,5,2,-3,-6,1/)
        tm0=thread_wtime()
        call tensor_block_contract(cptrn,ftens(1),ftens(2),ftens(3),ierr,'r4'); if(ierr.ne.0) then; ierr=61; goto 999; endif
!       write(jo_cp,'("#DEBUG(c_process::c_proc_test): CPU tensor contraction time = ",F10.4)') thread_wtime()-tm0 !debug
        call tens_blck_pack(ftens(0),'r4',pack_size(0),entry_ptr(0),entry_num(0),ierr); if(ierr.ne.0) then; ierr=62; goto 999; endif
        call tens_blck_pack(ftens(1),'r4',pack_size(1),entry_ptr(1),entry_num(1),ierr); if(ierr.ne.0) then; ierr=63; goto 999; endif
        call tens_blck_pack(ftens(2),'r4',pack_size(2),entry_ptr(2),entry_num(2),ierr); if(ierr.ne.0) then; ierr=64; goto 999; endif
        call tens_blck_pack(ftens(4),'r4',pack_size(4),entry_ptr(4),entry_num(4),ierr); if(ierr.ne.0) then; ierr=65; goto 999; endif
        call tens_blck_pack(ftens(5),'r4',pack_size(5),entry_ptr(5),entry_num(5),ierr); if(ierr.ne.0) then; ierr=66; goto 999; endif
        call tens_blck_pack(ftens(6),'r4',pack_size(6),entry_ptr(6),entry_num(6),ierr); if(ierr.ne.0) then; ierr=67; goto 999; endif
        call tens_blck_assoc(entry_ptr(0),ierr,ctens=ctens(0),gpu_num=gpu_id); if(ierr.ne.0) then; ierr=68; goto 999; endif
        call tens_blck_assoc(entry_ptr(1),ierr,ctens=ctens(1),gpu_num=gpu_id); if(ierr.ne.0) then; ierr=69; goto 999; endif
        call tens_blck_assoc(entry_ptr(2),ierr,ctens=ctens(2),gpu_num=gpu_id); if(ierr.ne.0) then; ierr=70; goto 999; endif
        call tens_blck_assoc(entry_ptr(4),ierr,ctens=ctens(4),gpu_num=gpu_id); if(ierr.ne.0) then; ierr=71; goto 999; endif
        call tens_blck_assoc(entry_ptr(5),ierr,ctens=ctens(5),gpu_num=gpu_id); if(ierr.ne.0) then; ierr=72; goto 999; endif
        call tens_blck_assoc(entry_ptr(6),ierr,ctens=ctens(6),gpu_num=gpu_id); if(ierr.ne.0) then; ierr=73; goto 999; endif
        o2n(0:6)=(/+1,1,2,3,4,5,6/);
        err_code=gpu_tensor_block_copy_dlf(o2n,ctens(1),ctens(4)); if(err_code.ne.0) then; ierr=74; goto 999; endif
        err_code=gpu_tensor_block_copy_dlf(o2n,ctens(2),ctens(5)); if(err_code.ne.0) then; ierr=75; goto 999; endif
        err_code=cuda_task_create(cuda_task(1)); if(err_code.ne.0) then; ierr=76; goto 999; endif
        err_code=cuda_task_create(cuda_task(2)); if(err_code.ne.0) then; ierr=77; goto 999; endif
        tm0=thread_wtime()
        call gpu_set_matmult_algorithm(BLAS_OFF)
        err_code=gpu_tensor_block_contract_dlf_(cptrn,ctens(1),ctens(2),ctens(0),COPY_BACK,cuda_task(1)); if(err_code.ne.0) then; write(jo_cp,*)'GPU ERROR = ',err_code; ierr=78; goto 999; endif
        call gpu_set_matmult_algorithm(BLAS_ON)
        err_code=gpu_tensor_block_contract_dlf_(cptrn,ctens(4),ctens(5),ctens(6),COPY_BACK,cuda_task(2)); if(err_code.ne.0) then; write(jo_cp,*)'GPU ERROR = ',err_code; ierr=79; goto 999; endif
        i1=0
        do while(i1.ne.3)
         if(cuda_task_complete(cuda_task(1)).ne.cuda_task_scheduled.and.(i1.eq.0.or.i1.eq.2)) then; tm1=thread_wtime()-tm0; i1=i1+1; endif
         if(cuda_task_complete(cuda_task(2)).ne.cuda_task_scheduled.and.(i1.eq.0.or.i1.eq.1)) then; tm2=thread_wtime()-tm0; i1=i1+2; endif
        enddo
!       write(jo_cp,*)'GPU Times/Status = ',tm1,tm2,cuda_task_complete(cuda_task(1)),cuda_task_complete(cuda_task(2)) !debug
        err_code=cuda_task_destroy(cuda_task(2)); if(err_code.ne.0) then; ierr=80; goto 999; endif
        err_code=cuda_task_destroy(cuda_task(1)); if(err_code.ne.0) then; ierr=81; goto 999; endif
        call tens_blck_unpack(ftens(0),entry_ptr(0),ierr); if(ierr.ne.0) then; ierr=82; goto 999; endif
        call tens_blck_unpack(ftens(6),entry_ptr(6),ierr); if(ierr.ne.0) then; ierr=83; goto 999; endif
!        write(jo_cp,*)'1-NORMS:',tensor_block_norm1(ftens(0),ierr,'r4'),tensor_block_norm1(ftens(3),ierr,'r4'),tensor_block_norm1(ftens(6),ierr,'r4') !debug
!        write(jo_cp,*) ftens(0)%data_real4(0:5); write(jo_cp,*) ftens(3)%data_real4(0:5); write(jo_cp,*) ftens(6)%data_real4(0:5) !debug
        cmp=tensor_block_cmp(ftens(0),ftens(3),ierr,'r4',.true.,1d-3,diffc); if(ierr.ne.0) then; ierr=84; goto 999; endif
        if(.not.cmp) then; write(jo_cp,*)'DIFF(1) COUNT = ',diffc; ierr=85; goto 999; endif
        cmp=tensor_block_cmp(ftens(6),ftens(3),ierr,'r4',.true.,1d-3,diffc); if(ierr.ne.0) then; ierr=86; goto 999; endif
        if(.not.cmp) then; write(jo_cp,*)'DIFF(2) COUNT = ',diffc; ierr=87; goto 999; endif
        call tens_blck_dissoc(ctens(6),ierr); if(ierr.ne.0) then; ierr=88; goto 999; endif
        call tens_blck_dissoc(ctens(5),ierr); if(ierr.ne.0) then; ierr=89; goto 999; endif
        call tens_blck_dissoc(ctens(4),ierr); if(ierr.ne.0) then; ierr=90; goto 999; endif
        call tens_blck_dissoc(ctens(2),ierr); if(ierr.ne.0) then; ierr=91; goto 999; endif
        call tens_blck_dissoc(ctens(1),ierr); if(ierr.ne.0) then; ierr=92; goto 999; endif
        call tens_blck_dissoc(ctens(0),ierr); if(ierr.ne.0) then; ierr=93; goto 999; endif
        err_code=free_buf_entry_host(entry_num(6)); if(err_code.ne.0) then; ierr=94; goto 999; endif
        err_code=free_buf_entry_host(entry_num(5)); if(err_code.ne.0) then; ierr=95; goto 999; endif
        err_code=free_buf_entry_host(entry_num(4)); if(err_code.ne.0) then; ierr=96; goto 999; endif
        err_code=free_buf_entry_host(entry_num(2)); if(err_code.ne.0) then; ierr=97; goto 999; endif
        err_code=free_buf_entry_host(entry_num(1)); if(err_code.ne.0) then; ierr=98; goto 999; endif
        err_code=free_buf_entry_host(entry_num(0)); if(err_code.ne.0) then; ierr=99; goto 999; endif
        call tensor_block_destroy(ftens(6),ierr); if(ierr.ne.0) then; ierr=100; goto 999; endif
        call tensor_block_destroy(ftens(5),ierr); if(ierr.ne.0) then; ierr=101; goto 999; endif
        call tensor_block_destroy(ftens(4),ierr); if(ierr.ne.0) then; ierr=102; goto 999; endif
        call tensor_block_destroy(ftens(3),ierr); if(ierr.ne.0) then; ierr=103; goto 999; endif
        call tensor_block_destroy(ftens(2),ierr); if(ierr.ne.0) then; ierr=104; goto 999; endif
        call tensor_block_destroy(ftens(1),ierr); if(ierr.ne.0) then; ierr=105; goto 999; endif
        call tensor_block_destroy(ftens(0),ierr); if(ierr.ne.0) then; ierr=106; goto 999; endif
#endif
!Clean:
999     do i=0,test_args_lim; err_code=free_buf_entry_host(entry_num(i)); enddo
        do i=0,test_args_lim; call tensor_block_destroy(ftens(i),j); enddo
        if(ierr.eq.0) then; write(jo_cp,'("Ok")'); else; write(jo_cp,'("Failed: ERROR #",i4)') ierr; endif
        return
        end subroutine c_proc_test
!--------------------------------------------------------------------------------------------------------
        subroutine run_benchmarks(ierr)
!This subroutine runs computationally intensive (single-process) benchmarks of tensor algebra on CPU/GPU.
        implicit none
        integer, intent(inout):: ierr
        integer, parameter:: num_tens_sizes=1,num_tens_ranks=8,num_dim_spreads=3
        integer(8), parameter:: tens_sizes(1:num_tens_sizes)=(/999999/)
        integer, parameter:: tens_ranks(1:num_tens_ranks)=(/2,3,4,5,6,7,8,15/)
        integer, parameter:: dim_spreads(1:num_dim_spreads)=(/1,5,15/)
        integer, parameter:: test_args_lim=15
        integer i,j,k,l,m,n
        integer tsl,tens_rank,dim_spread,o2n(0:max_tensor_rank),n2o(0:max_tensor_rank),nfail
        integer(8) tens_size,diffc
        character(256) tshape
        real(8) tm
        logical cmp
        integer(C_INT) i0,i1,gpu_id,err_code
        integer(C_INT):: entry_num(0:test_args_lim)=-1
        type(C_PTR):: entry_ptr(0:test_args_lim)=C_NULL_PTR
        integer(C_SIZE_T):: s0,pack_size(0:test_args_lim)=0
        type(C_PTR):: cuda_task(1:3)=C_NULL_PTR
        type(C_PTR):: ctens(0:test_args_lim)=C_NULL_PTR
        type(tensor_block_t) ftens(0:test_args_lim)

        ierr=0; write(jo_cp,'("#MSG(c_process::c_proc_test): Tensor algebra benchmarks:")')
#ifndef NO_GPU
        gpu_id=gpu_busy_least(); call cudaSetDevice(gpu_id,err_code); if(err_code.ne.0) then; ierr=666; goto 999; endif
        call gpu_set_event_policy(EVENTS_ON)
#endif
!TENSOR TRANSPOSE:
        write(jo_cp,'(" TENSOR TRANSPOSE:")')
        nfail=0 !will be the total number of failed transposes
        do m=1,num_tens_sizes
         tens_size=tens_sizes(m)
         do n=1,num_tens_ranks
          tens_rank=tens_ranks(n)
          do k=1,num_dim_spreads
           dim_spread=dim_spreads(k)
           call tensor_shape_rnd(tshape,tsl,ierr,tens_size,tens_rank,dim_spread); if(ierr.ne.0) then; ierr=1; goto 999; endif
!           tens_rank=3; tsl=11; tshape(1:tsl)='(61,60,61)' !debug
           call tensor_block_create(tshape(1:tsl),'r4',ftens(1),ierr); if(ierr.ne.0) then; ierr=2; goto 999; endif
           call printl(jo_cp,'  New Tensor Shape: '//tshape(1:tsl)//': ',.false.)
           write(jo_cp,'(i10,1x,F16.4)') ftens(1)%tensor_block_size,tensor_block_norm1(ftens(1),ierr,'r4')
           write(jo_cp,'(3x)',advance='no')
           tm=thread_wtime()
           call tensor_block_copy(ftens(1),ftens(0),ierr); if(ierr.ne.0) then; ierr=3; goto 999; endif
           write(jo_cp,'("#DEBUG(tensor_algebra:tensor_block_copy_dlf): Direct time ",F10.6)') thread_wtime()-tm !debug
           call random_permutation(tens_rank,o2n)
!           o2n(1:tens_rank)=(/3,2,1/) !debug
           call permutation_converter(.false.,tens_rank,n2o,o2n)
           write(jo_cp,'(3x,"Permutation:",32(1x,i2))') o2n(1:tens_rank)
           write(jo_cp,'(3x,"Permutation:",32(1x,i2))') n2o(1:tens_rank)
           call set_transpose_algorithm(EFF_TRN_OFF) !scatter
           write(jo_cp,'(3x)',advance='no')
           tm=thread_wtime()
           call tensor_block_copy(ftens(1),ftens(0),ierr,o2n); if(ierr.ne.0) then; ierr=4; goto 999; endif
           write(jo_cp,'("#DEBUG(tensor_algebra:tensor_block_copy_scatter_dlf): Time ",F10.6)') thread_wtime()-tm !debug
           write(jo_cp,'(3x)',advance='no')
           tm=thread_wtime()
           call tensor_block_copy(ftens(0),ftens(2),ierr,n2o); if(ierr.ne.0) then; ierr=5; goto 999; endif
           write(jo_cp,'("#DEBUG(tensor_algebra:tensor_block_copy_scatter_dlf): Time ",F10.6)') thread_wtime()-tm !debug
           cmp=tensor_block_cmp(ftens(1),ftens(2),ierr,'r4',.true.,1d-4,diffc); if(ierr.ne.0) then; ierr=6; goto 999; endif
           write(jo_cp,'(3x,l1,1x,i9,1x,F16.4)') cmp,diffc,tensor_block_norm1(ftens(2),ierr,'r4')
           if(.not.cmp) then; nfail=nfail+1; write(jo_cp,'(3x,"Comparison Failed!")'); endif
           call set_transpose_algorithm(EFF_TRN_ON) !cache-efficient
           write(jo_cp,'(3x)',advance='no')
           tm=thread_wtime()
           call tensor_block_copy(ftens(1),ftens(0),ierr,o2n); if(ierr.ne.0) then; ierr=7; goto 999; endif
           write(jo_cp,'("#DEBUG(tensor_algebra:tensor_block_copy_dlf): Time ",F10.6)') thread_wtime()-tm !debug
           write(jo_cp,'(3x)',advance='no')
           tm=thread_wtime()
           call tensor_block_copy(ftens(0),ftens(2),ierr,n2o); if(ierr.ne.0) then; ierr=8; goto 999; endif
           write(jo_cp,'("#DEBUG(tensor_algebra:tensor_block_copy_dlf): Time ",F10.6)') thread_wtime()-tm !debug
           cmp=tensor_block_cmp(ftens(1),ftens(2),ierr,'r4',.true.,1d-4,diffc); if(ierr.ne.0) then; ierr=9; goto 999; endif
           write(jo_cp,'(3x,l1,1x,i9,1x,F16.4)') cmp,diffc,tensor_block_norm1(ftens(2),ierr,'r4')
           if(.not.cmp) then; nfail=nfail+1; write(jo_cp,'(3x,"Comparison Failed!")'); endif
#ifndef NO_GPU
           call tensor_block_init('r4',ftens(0),ierr,val_r4=0.0); if(ierr.ne.0) then; ierr=10; goto 999; endif
           call tensor_block_init('r4',ftens(2),ierr,val_r4=0.0); if(ierr.ne.0) then; ierr=11; goto 999; endif
           call tens_blck_pack(ftens(0),'r4',pack_size(0),entry_ptr(0),entry_num(0),ierr); if(ierr.ne.0) then; ierr=12; goto 999; endif
           call tens_blck_pack(ftens(1),'r4',pack_size(1),entry_ptr(1),entry_num(1),ierr); if(ierr.ne.0) then; ierr=13; goto 999; endif
           call tens_blck_pack(ftens(2),'r4',pack_size(2),entry_ptr(2),entry_num(2),ierr); if(ierr.ne.0) then; ierr=14; goto 999; endif
           call tens_blck_assoc(entry_ptr(0),ierr,ctens=ctens(0),gpu_num=gpu_id); if(ierr.ne.0) then; ierr=15; goto 999; endif
           call tens_blck_assoc(entry_ptr(1),ierr,ctens=ctens(1),gpu_num=gpu_id); if(ierr.ne.0) then; ierr=16; goto 999; endif
           call tens_blck_assoc(entry_ptr(2),ierr,ctens=ctens(2),gpu_num=gpu_id); if(ierr.ne.0) then; ierr=17; goto 999; endif
           call gpu_set_transpose_algorithm(EFF_TRN_OFF) !scatter on GPU
           write(jo_cp,'(3x)',advance='no')
           err_code=gpu_tensor_block_copy_dlf(o2n,ctens(1),ctens(0)); if(err_code.ne.0) then; write(jo_cp,*)'GPU error ',err_code; call print_gpu_debug_dump(jo_cp); ierr=18; goto 999; endif
           write(jo_cp,'(3x)',advance='no')
           err_code=gpu_tensor_block_copy_dlf(n2o,ctens(0),ctens(2)); if(err_code.ne.0) then; write(jo_cp,*)'GPU error ',err_code; call print_gpu_debug_dump(jo_cp); ierr=19; goto 999; endif
           call tens_blck_unpack(ftens(2),entry_ptr(2),ierr); if(ierr.ne.0) then; ierr=20; goto 999; endif
           cmp=tensor_block_cmp(ftens(1),ftens(2),ierr,'r4',.true.,1d-4,diffc); if(ierr.ne.0) then; ierr=21; goto 999; endif
           write(jo_cp,'(3x,l1,1x,i9,1x,F16.4)') cmp,diffc,tensor_block_norm1(ftens(2),ierr,'r4')
           if(.not.cmp) then; nfail=nfail+1; write(jo_cp,'(3x,"Comparison Failed!")'); endif
           call gpu_set_transpose_algorithm(EFF_TRN_ON) !shared-memory on GPU
           err_code=cuda_task_create(cuda_task(1)); if(err_code.ne.0) then; ierr=22; goto 999; endif
           write(jo_cp,'(3x)',advance='no')
           tm=thread_wtime()
           err_code=gpu_tensor_block_copy_dlf_(o2n,ctens(1),ctens(0),COPY_BACK,cuda_task(1)); if(err_code.ne.0) then; write(jo_cp,*)'GPU error ',err_code; call print_gpu_debug_dump(jo_cp); ierr=23; goto 999; endif
           err_code=cuda_task_wait(cuda_task(1)); if(err_code.ne.cuda_task_completed) then; ierr=24; goto 999; endif
           write(jo_cp,'("#DEBUG(tensor_algebra_gpu_nvidia:gpu_tensor_block_copy_dlf_): Time ",F10.6)') thread_wtime()-tm
           write(jo_cp,'(3x)',advance='no')
           err_code=gpu_tensor_block_copy_dlf(n2o,ctens(0),ctens(2)); if(err_code.ne.0) then; write(jo_cp,*)'GPU error ',err_code; call print_gpu_debug_dump(jo_cp); ierr=25; goto 999; endif
           err_code=cuda_task_destroy(cuda_task(1)); if(err_code.ne.0) then; ierr=26; goto 999; endif
           call tens_blck_unpack(ftens(2),entry_ptr(2),ierr); if(ierr.ne.0) then; ierr=27; goto 999; endif
           cmp=tensor_block_cmp(ftens(1),ftens(2),ierr,'r4',.true.,1d-4,diffc); if(ierr.ne.0) then; ierr=28; goto 999; endif
           write(jo_cp,'(3x,l1,1x,i9,1x,F16.4)') cmp,diffc,tensor_block_norm1(ftens(2),ierr,'r4')
           if(.not.cmp) then; nfail=nfail+1; write(jo_cp,'(3x,"Comparison Failed!")'); endif
           write(jo_cp,'(3x)',advance='no')
           n2o(1:tens_rank)=(/(i,i=1,tens_rank)/)
           err_code=gpu_tensor_block_copy_dlf(n2o,ctens(1),ctens(2)); if(err_code.ne.0) then; ierr=29; goto 999; endif
           call tens_blck_dissoc(ctens(2),ierr); if(ierr.ne.0) then; ierr=30; goto 999; endif
           call tens_blck_dissoc(ctens(1),ierr); if(ierr.ne.0) then; ierr=31; goto 999; endif
           call tens_blck_dissoc(ctens(0),ierr); if(ierr.ne.0) then; ierr=32; goto 999; endif
           err_code=free_buf_entry_host(entry_num(2)); if(err_code.ne.0) then; ierr=33; goto 999; endif
           err_code=free_buf_entry_host(entry_num(1)); if(err_code.ne.0) then; ierr=34; goto 999; endif
           err_code=free_buf_entry_host(entry_num(0)); if(err_code.ne.0) then; ierr=35; goto 999; endif
#endif
           call tensor_block_destroy(ftens(2),ierr)
           call tensor_block_destroy(ftens(1),ierr)
           call tensor_block_destroy(ftens(0),ierr)
          enddo !dim_spread
         enddo !tens_rank
        enddo !tens_size

999     if(ierr.eq.0) then
         write(jo_cp,'("Done: ",i6," failed transposes.")') nfail
        else
         write(jo_cp,'("Failed: Error #",i7)') ierr
         call tensor_block_destroy(ftens(0),i)
         call tensor_block_destroy(ftens(1),i)
         call tensor_block_destroy(ftens(2),i)
        endif
        return
        end subroutine run_benchmarks
!---------------------------------------------------------------------------------------------
#ifndef NO_GPU
        subroutine print_gpu_debug_dump(ifh) !prints the current GPU debug dump to output #ifh
        implicit none
        integer, intent(in):: ifh
        integer(C_INT), parameter:: max_dump_size=1024
        integer(C_INT) i,dump_size,dump(1:max_dump_size)
        write(ifh,*)' '; write(ifh,'("#DEBUG(c_process::print_gpu_debug_dump): GPU DEBUG DUMP:")')
        dump_size=gpu_get_debug_dump(dump)
        do i=1,dump_size,10; write(ifh,'(10(1x,i9))') dump(i:min(dump_size,i+9)); enddo
        return
        end subroutine print_gpu_debug_dump
#endif

       end module c_process
