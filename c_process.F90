!This module provides functionality for a Computing Process (C-PROCESS, CP).
!In essence, this is a single-node elementary tensor instruction scheduler (SETIS).
!AUTHOR: Dmitry I. Lyakh (Dmytro I. Liakh): quant4me@gmail.com
!REVISION: 2014/07/10
!CONCEPTS (CP workflow):
! - Each CP stores its own tensor blocks in TBB, with a possibility of disk dump.
! - LR sends a batch of ETI to be executed on this CP unit (CP MPI Process).
!   The location of each tensor-block operand in ETI is already given there,
!   as well as some other characteristics (approx. operation cost, memory requirements, etc.).
! - CP places the ETI received into multiple queues, based on the number of remote operands and ETI cost.
! - Each non-trivial ETI operand (tensor block) is assigned an AAR entry that points to
!   either a <tensor_block_t> entity (to be used on CPU or Intel Xeon Phi)
!   or a <tensBlck_t> entity (to be used on Nvidia GPU). Once the corresponding
!   data (tensor block) is in local memory of CP (either in TBB, or HAB, or DEB),
!   a new AAR entry is defined and all ETI operands, referring to it, become ready.
!   Data (tensor blocks) registered in AAR can be reused while locally present.
!   AAR is the only way to access a tensor block for ETI (the tensor block itself
!   can physically reside in either TBB, or HAB, or DEB).
! - ETI with all operands ready can be scheduled for execution on an appropriate device:
!    - On each distinct device, ETI are issued in-order but executed asynchronously;
!    - ETI with ready input arguments can be further prioritized;
!    - Each particular type of instruction has its own issuing workflow for each device kind.
!      The device is chosen based on the ETI cost, ETI cost/size ratio, and device availability.
!    - Once issued successfully, the ETI obtains a query handle that can be used for completion checks.
!    - AAR entries and associated data used by active ETI must not be freed before completion of the ETI.
!    - If the ETI result is remote, its destination must be updated before
!      reporting to LR that the ETI has been completed.
!NOTES:
! - Data synchronization in an instance of <tensor_block_t> (Fortran)
!   associated with a Host Argument Buffer entry can allocate only regular CPU memory
!   (the one outside the pinned Host Argument buffer). Hence the newly allocated
!   memory is not in HAB, thus cannot be used with Nvidia GPUs.
!Acronyms:
! - CP - Computing (MPI) Process;
! - MT - Master Thread;
! - LST - Leading Slave Thread;
! - STCU - Slave Threads Computing Unit;
! - NVCU - Nvidia GPU Computing Unit;
! - XPCU - Intel Xeon Phi Computing Unit;
! - AMCU - AMD GPU Computing Unit;
! - ETI - Elementary Tensor Instruction;
! - ETIS - Elementary Tensor Instruction Scheduler (SETIS, LETIS, GETIS);
! - ETIQ - Elementary Tensor Instruction Queue (out-of-order issuance);
! - HAB - Host Argument Buffer;
! - GAB - GPU Argument Buffer;
! - DEB - Data Exchange Buffer (additional buffer for temporary present tensor blocks);
! - TBB - Tensor Block Bank (storage);
! - AAR - Active Argument Register;
! - TAL - Tensor Algebra Library (CPU, GPU, MIC, etc.);
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
 !Debugging:
        integer, private:: debug_level=0          !debug level (0: normal execution)
        integer, private:: debug_step_pause=1000  !pause (msec) before advancing MT to the next instruction while debugging
 !General:
        integer, parameter, private:: CZ=C_SIZE_T
 !Elementary Tensor Instruction Scheduler (ETIS):
  !Host Argument Buffer (HAB):
        integer(C_SIZE_T), parameter, private:: max_hab_size=4096_CZ*(1024_CZ*1024_CZ) !max size in bytes of the HAB
        integer(C_SIZE_T), parameter, private:: min_hab_size=64_CZ*(1024_CZ*1024_CZ)   !min size in bytes of the HAB
  !Elementary Tensor Instruction Queue (ETIQ):
        integer(C_INT), parameter, private:: etiq_max_depth=65536 !max number of simultaneously scheduled ETI at this CP
        integer(C_INT), parameter, private:: etiq_loc_levels=4    !number of locality levels in ETIQ (senior)
        integer(C_INT), parameter, private:: etiq_cost_levels=5   !number of cost levels in ETIQ (minor)
        integer(C_INT), parameter, private:: etiq_channels=etiq_cost_levels*etiq_loc_levels !total number of levels in ETIQ
        integer(C_INT), parameter, private:: etiq_stcu_max_depth=4096 !max number of simultaneously scheduled ETI on STCU
        integer(C_INT), parameter, private:: etiq_nvcu_max_depth=256  !max number of simultaneously scheduled ETI on NVCU
        integer(C_INT), parameter, private:: etiq_xpcu_max_depth=512  !max number of simultaneously scheduled ETI on XPCU
        integer(C_INT), parameter, public:: eti_buf_size=32768 !ETI buffer size (for communications with LR)
        integer(C_INT), parameter, private:: stcu_max_units=64 !max number of STCU units
  !Tensor naming:
        integer(C_INT), parameter:: tensor_name_len=32 !max number of characters used for tensor names
  !Tensor instruction status (any negative value will correspond to a failure, designating the error code):
        integer, parameter:: instr_null=0              !uninitialized instruction (empty)
        integer, parameter:: instr_fresh=1             !instruction has not been seen before (new instruction)
        integer, parameter:: instr_data_wait=2         !instruction is waiting for the input data to arrive
        integer, parameter:: instr_ready_to_exec=3     !instruction is ready to be executed (input data has arrived)
        integer, parameter:: instr_scheduled=4         !instruction is placed into an execution queue on a specific CU
        integer, parameter:: instr_issued=5            !instruction is issued for execution on some computing unit (CU)
        integer, parameter:: instr_completed=6         !instruction has completed (but the result may still need to be remotely uploaded)
        integer, parameter:: instr_dead=7              !instruction can safely be removed from the queue
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
         integer(C_INT):: device_type  !device type: see tensor_algebra_gpu_nvidia.inc/tensor_algebra_gpu_nvidia.h
         integer(C_INT):: unit_number  !logical number of the device within its type (0..max)
        end type cu_t
 !Task execution configuration:
        type te_conf_t
         type(cu_t):: cu_id            !computing unit
         integer(C_INT):: num_gangs    !coarse grain parallelism (e.g., CUDA blocks)
         integer(C_INT):: num_workers  !fine grain parallelism (e.g., OMP threads, CUDA warps)
         integer(C_INT):: vec_width    !vector width (e.g., SSE, AVX, CUDA threads in a warp)
        end type te_conf_t
 !Tensor block identifier (key):
        type tens_blck_id_t
         character(LEN=tensor_name_len):: tens_name  !tensor name
         integer(C_INT), allocatable:: tens_mlndx(:) !tensor block multi-index (0th element is the length of the multi-index)
         contains
          procedure, public:: create=>tens_blck_id_create
          procedure, public:: destroy=>tens_blck_id_destroy
        end type tens_blck_id_t
 !Tensor Block Bank (TBB):
  !Entry of TBB (named tensor block):
        type, private:: tbb_entry_t
         type(tensor_block_t), private:: tens_blck !fortran tensor block
         integer, private:: file_handle=0          !file handle (0: stored in RAM)
         integer(8), private:: file_offset=0_8     !file offset where the corresponding packet is stored (if on disk)
         integer(8), private:: stored_size=0_8     !stored (disk) size of the tensor block (packet) in bytes
        end type tbb_entry_t
 !Elementary Tensor Instruction Scheduler (ETIS):
  !Locally present tensor argument:
        type, private:: tens_arg_t
         type(tensor_block_t), pointer, private:: tens_blck_f=>NULL() !pointer to a tensor_block in TBB, or HAB, or DEB
         type(C_PTR), private:: tens_blck_c=C_NULL_PTR !C pointer to tensBlck_t in HAB (see "tensor_algebra_gpu_nvidia.h")
         integer(C_INT), private:: buf_entry_host=-1 !HAB entry number where the tensor block resides as a packet (-1: undefined)
         integer, private:: mpi_tag                  !>0: MPI_TAG for a pending MPI_IRECV; =0: MPI delivered; <0: local
         integer, private:: mpi_process              !>=0: mpi rank of the host process; <0: local
         integer, private:: times_needed             !number of times the argument has been referred to in ETIQ
         integer, private:: times_used               !number of times this AAR entry (argument) has been used since creation
        end type tens_arg_t
  !Tensor operand type (component of ETI):
        type tens_operand_t
         type(tens_blck_id_t):: tens_blck_id !tensor block identifier (key): set by LR
         integer(C_INT):: op_host            !MPI process rank where the tensor operand resides: set by LR
         integer(8):: op_pack_size           !packed size of the tensor operand in bytes: computed by LR
         integer(C_INT):: op_tag             !MPI message tag by which the tensor operand is to be delivered (-1: local): set by LR
         integer(C_INT):: op_price           !current price of the tensor operand (tensor block): set by LR
         type(tens_arg_t), pointer, private:: op_aar_entry=>NULL() !AAR entry assigned to tensor operand: set by CP:MT
         contains
          procedure, public:: pack=>tens_operand_pack
          procedure, public:: unpack=>tens_operand_unpack
        end type tens_operand_t
  !Dispatched elementary tensor instruction (ETI):
        type tens_instr_t
         integer(2):: instr_code             !tensor instruction code (see above): set by GR
         character(2):: data_kind            !data kind {R4|R8|C8}: set by {GR|LR}
         integer, allocatable:: instr_aux(:) !auxiliary instruction info (contraction pattern, permutation, etc.): set by GR
         integer:: instr_priority            !tensor instruction priority: set by LR
         real(4):: instr_cost                !approx. instruction computational cost (FLOPs): set by LR
         real(4):: instr_size                !approx. instruction memory demands (Bytes): set by LR
         integer, private:: instr_status     !tensor instruction status (see above): set by CP:{MT|LST}
         type(cu_t), private:: instr_cu      !computing unit which the instruction is issued to: set by CP:MT
         integer, private:: instr_handle     !instruction handle to query its status: set by CP:MT (accelerators only)
         real(4), private:: time_touched     !time the instruction was touched for the first time: set by MT
         real(4), private:: time_data_ready  !time the instruction has all input data arrived: set by MT
         real(4), private:: time_issued      !time the instruction was issued for execution: set by CP:{MT:ST}
         real(4), private:: time_completed   !time the instruction was completed: set by CP:{MT:ST}
         real(4), private:: time_uploaded    !time the remote result upload completed: set by MT
         integer, private:: args_ready       !each bit is set to 1 when the corresponding operand is in AAR: set by CP:MT
         type(tens_operand_t):: tens_op(0:max_tensor_operands-1) !tensor-block operands
         contains
          procedure, public:: pack=>eti_pack
          procedure, public:: unpack=>eti_unpack
        end type tens_instr_t
  !Elementary tensor instruction queue (ETIQ), out-of-order (linked list):
        type, private:: etiq_t
         integer(C_INT), private:: depth=0                    !total number of ETIQ entries
         integer(C_INT), private:: scheduled=0                !total number of active ETIQ entries
         integer(C_INT), private:: ffe_sp=0                   !ETIQ FFE stack pointer
         integer(C_INT), private:: last(0:etiq_channels-1)=0  !last scheduled instruction for each ETIQ channel
         integer(C_INT), private:: ip(0:etiq_channels-1)=0    !instruction pointers for each ETIQ channel
         integer(C_INT), private:: ic(0:etiq_channels-1)=0    !instruction counters for each ETIQ channel
         integer(C_INT), allocatable, private:: ffe_stack(:)  !ETIQ FFE stack
         integer(C_INT), allocatable, private:: next(:)       !ETIQ next linking
         type(tens_instr_t), allocatable, private:: eti(:)    !elementary tensor instructions
         contains
          procedure, private:: init=>etiq_init
        end type etiq_t
  !In-order elementary tensor instruction queue for a specific computing unit:
        type, private:: etiq_cu_t
         integer(C_INT), private:: depth=0                    !total number of entries in the queue
         integer(C_INT), private:: scheduled=0                !total number of active entries in the queue
         integer(C_INT), private:: bp=0                       !base instruction pointer (first issued unfinished instruction): async units only
         integer(C_INT), private:: ip=0                       !current instruction pointer (current instruction)
         integer(C_INT), allocatable, private:: etiq_entry(:) !number of the ETIQ entry where the ETI is located
         type(te_conf_t), allocatable, private:: te_conf(:)   !computing unit configuration which ETI will be executed on
         contains
          procedure, private:: init=>etiq_cu_init
        end type etiq_cu_t
 !NVCU task:
        type, private:: nvcu_task_t
         type(C_PTR), private:: cuda_task_handle !CUDA task handle that can be used for querying the CUDA task status
        end type nvcu_task_t
 !XPCU task:
        type, private:: xpcu_task_t !MIC task handle that can be used for querying the MIC task status
         integer, private:: mic_id
         integer, private:: mic_task_handle
        end type xpcu_task_t
!PROCEDURE VISIBILITY:
        private tens_blck_id_create
        private tens_blck_id_destroy
        private etiq_init
        private etiq_cu_init
!DATA:
 !Tensor Block Bank (TBB):
        type(dict_t), private:: tbb !permanent storage of local tensor blocks
 !Elementary Tensor Instruction Scheduler (ETIS):
  !Active Argument Register (AAR):
        type(dict_t), private:: aar !register for locally present (TBB, HAB, DEB) tensor blocks participating in current ETIs
  !Elementary Tensor Instruction Queue (ETIQ):
        type(etiq_t), target, private:: etiq     !multi-channel out-of-order ETI queue (linear numeration starts from 1)
  !STCU ETI queue (ETIQ_STCU):
        type(etiq_cu_t), private:: etiq_stcu     !circular in-order ETI queue for STCU (numeration starts from 0)
  !NVCU ETI queue (ETIQ_NVCU:
        type(etiq_cu_t), private:: etiq_nvcu     !circular in-order ETI queue for NVCU (numeration starts from 0)
  !XPCU ETI queue (ETIQ_XPCU):
        type(etiq_cu_t), private:: etiq_xpcu     !circular in-order ETI queue for XPCU (numeration starts from 0)
  !Host Argument Buffer (HAB):
        integer(C_SIZE_T), private:: hab_size=0  !actual size in bytes of the Host argument buffer (HAB)
        integer(C_INT), private:: max_hab_args=0 !max number of arguments (of lowest-size level) that can fit in HAB
  !ETI communication buffers:
        integer(4), private:: eti_buf_in(0:eti_buf_size-1)  !ETI input buffer (incoming ETI from LR)
        integer(4), private:: eti_buf_out(0:eti_buf_size-1) !ETI output buffer (outgoing done ETI to LR)
 !CP runtime:
        integer, private:: mt_error              !master thread error (if != 0)
        integer, private:: stcu_error            !slave threads computing unit (STCU) error (if != 0)
        integer, private:: stcu_num_units        !current number of STCU units
!------------------------------------------------
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
!---------------------------------------------------------
        integer(C_INT), parameter:: max_arg_buf_levels=256 !max number of argument buffer levels (do not exceed C values)
!---------------------------------------------------------
        integer(C_INT) i,j,k,l,m,n,ks,kf,err_code,thread_num !general purpose: thread private
        type(nvcu_task_t), allocatable:: nvcu_tasks(:) !parallel to etiq_nvcu
        type(xpcu_task_t), allocatable:: xpcu_tasks(:) !parallel to etiq_xpcu
        integer:: stcu_base_ip,stcu_my_ip,stcu_my_eti !STCU specific (slaves): thread private
        integer(C_SIZE_T):: blck_sizes(0:max_arg_buf_levels-1)
        integer:: tree_height,stcu_mimd_my_pass,stcu_mimd_max_pass
        integer(8):: tree_volume
        type(tens_blck_id_t):: key(0:15)
        type(tensor_block_t):: tens
        type(tensor_block_t), pointer:: tens_p
        type(tbb_entry_t):: tbb_entry
        type(tbb_entry_t), pointer:: tbb_entry_p
        type(tens_arg_t), pointer:: targ_p
        class(*), pointer:: uptr
        real(8) tm,tm0,tm1

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
        ierr=etiq%init(etiq_max_depth)
        if(ierr.ne.0) then; write(jo_cp,'("#ERROR(c_process::c_proc_life): ETIQ allocation failed!")'); call c_proc_quit(5); return; endif
        tm=thread_wtime()-tm; write(jo_cp,'("Ok(",F4.1," sec): ETIQ total depth = ",i7)') tm,etiq%depth
        write(jo_cp,'("#MSG(c_process::c_proc_life): Number of ETIQ channels = ",i3)') etiq_channels
!Test C-process functionality (debug):
!        call c_proc_test(ierr); if(ierr.ne.0) then; write(jo_cp,'("#ERROR(c_process::c_proc_life): C-process functionality test failed: ",i7)') ierr; call c_proc_quit(6); return; endif
!        call run_benchmarks(ierr); if(ierr.ne.0) then; write(jo_cp,'("#ERROR(c_process::c_proc_life): C-process benchmarking failed: ",i7)') ierr; call c_proc_quit(7); return; endif
!---------------------------------
!LIFE:
        ierr=0
#ifndef NO_OMP
        call omp_set_dynamic(.false.); call omp_set_nested(.true.)
!Init ETI queues for all computing units present on the node:
        n=omp_get_max_threads()
        if(n.ge.2.and.n.eq.max_threads) then !At least 1 master + 1 slave threads
 !Init STCU ETIQ (slave CPU threads):
         write(jo_cp,'("#MSG(c_process::c_proc_life): Allocating STCU ETIQ ... ")',advance='no')
         tm=thread_wtime()
         ierr=etiq_stcu%init(etiq_stcu_max_depth)
         if(ierr.ne.0) then; write(jo_cp,'("#ERROR(c_process::c_proc_life): ETIQ_STCU allocation failed!")'); call c_proc_quit(8); return; endif
         stcu_num_units=-1
         tm=thread_wtime()-tm; write(jo_cp,'("Ok(",F4.1," sec): STCU ETIQ depth = ",i6)') tm,etiq_stcu%depth
         write(jo_cp,'("#MSG(c_process::c_proc_life): Max number of STCU MIMD units = ",i5)') min(max_threads-1,stcu_max_units)
 !Init NVCU ETIQ (Nvidia GPU's):
#ifndef NO_GPU
         write(jo_cp,'("#MSG(c_process::c_proc_life): Allocating NVCU ETIQ ... ")',advance='no')
         tm=thread_wtime()
         ierr=etiq_nvcu%init(etiq_nvcu_max_depth)
         if(ierr.ne.0) then; write(jo_cp,'("#ERROR(c_process::c_proc_life): ETIQ_NVCU allocation failed!")'); call c_proc_quit(9); return; endif
         allocate(nvcu_tasks(0:etiq_nvcu_max_depth-1),STAT=ierr) !table for active CUDA tasks
         if(ierr.ne.0) then; write(jo_cp,'("#ERROR(c_process::c_proc_life): ETIQ_NVCU task table allocation failed!")'); call c_proc_quit(10); return; endif
         tm=thread_wtime()-tm; write(jo_cp,'("Ok(",F4.1," sec): NVCU ETIQ depth = ",i6)') tm,etiq_nvcu%depth
#endif
 !Init XPCU ETIQ (Intel Xeon Phi's):
#ifndef NO_PHI
         write(jo_cp,'("#MSG(c_process::c_proc_life): Allocating XPCU ETIQ ... ")',advance='no')
         tm=thread_wtime()
         ierr=etiq_xpcu%init(etiq_xpcu_max_depth)
         if(ierr.ne.0) then; write(jo_cp,'("#ERROR(c_process::c_proc_life): ETIQ_XPCU allocation failed!")'); call c_proc_quit(11); return; endif
         allocate(xpcu_tasks(0:etiq_xpcu_max_depth-1),STAT=ierr)
         if(ierr.ne.0) then; write(jo_cp,'("#ERROR(c_process::c_proc_life): ETIQ_XPCU task table allocation failed!")'); call c_proc_quit(12); return; endif
         tm=thread_wtime()-tm; write(jo_cp,'("Ok(",F4.1," sec): XPCU ETIQ depth = ",i6)') tm,etiq_xpcu%depth
#endif
!Begin active life (master thread + leading slave thread):
!$OMP PARALLEL NUM_THREADS(2) DEFAULT(SHARED) &
!$OMP          PRIVATE(thread_num,stcu_base_ip,stcu_my_ip,stcu_my_eti,stcu_mimd_my_pass,i,j,k,l,m,n,err_code)
         n=omp_get_num_threads()
         if(n.eq.2) then
          thread_num=omp_get_thread_num()
!Master thread:
          if(thread_num.eq.0) then
           mt_error=0 !set/used only by MT
           master_life: do !MT life cycle
            !`Write
!DEBUG begin:
 !Create keys (tensor block id):
            err_code=key( 0)%create('O',(/0/)); write(jo_cp,*) 'Key creation: ',err_code
            err_code=key( 1)%create('I',(/0/)); write(jo_cp,*) 'Key creation: ',err_code
            err_code=key( 2)%create('A',(/4,0,0,0,0/)); write(jo_cp,*) 'Key creation: ',err_code
            err_code=key( 3)%create('B',(/4,0,0,0,0/)); write(jo_cp,*) 'Key creation: ',err_code
            err_code=key( 4)%create('C',(/4,0,0,0,0/)); write(jo_cp,*) 'Key creation: ',err_code
            err_code=key( 5)%create('D',(/5,0,0,0,0,0/)); write(jo_cp,*) 'Key creation: ',err_code
            err_code=key( 6)%create('E',(/5,0,0,0,0,0/)); write(jo_cp,*) 'Key creation: ',err_code
            err_code=key( 7)%create('F',(/5,0,0,0,0,0/)); write(jo_cp,*) 'Key creation: ',err_code
            err_code=key( 8)%create('G',(/5,0,0,0,0,0/)); write(jo_cp,*) 'Key creation: ',err_code
            err_code=key( 9)%create('H',(/5,0,0,0,0,0/)); write(jo_cp,*) 'Key creation: ',err_code
            err_code=key(10)%create('J',(/5,0,0,0,0,0/)); write(jo_cp,*) 'Key creation: ',err_code
 !Register tensor blocks in TBB:
  !Zero scalar:
            call tensor_block_create('()','r8',tbb_entry%tens_blck,i,val_r8=0d0); write(jo_cp,*) 'Tensor block creation: ',i
            err_code=tbb%search(dict_add_if_not_found,tens_key_cmp,key(0),tbb_entry,value_out=uptr); write(jo_cp,*) 'TBB entry search: ',err_code
            select type (uptr)
            type is (tbb_entry_t)
             err_code=aar_register(key(0),targ_p); write(jo_cp,*) 'AAR entry created: ',err_code
             targ_p%tens_blck_f=>uptr%tens_blck
            end select
  !Unity scalar:
            call tensor_block_init('r8',tbb_entry%tens_blck,i,val_r8=1d-3); write(jo_cp,*) 'Tensor block unity init: ',i
            err_code=tbb%search(dict_add_if_not_found,tens_key_cmp,key(1),tbb_entry,value_out=uptr); write(jo_cp,*) 'TBB entry search: ',err_code
            select type (uptr)
            type is (tbb_entry_t)
             err_code=aar_register(key(1),targ_p); write(jo_cp,*) 'AAR entry created: ',err_code
             targ_p%tens_blck_f=>uptr%tens_blck
            end select
  !Non-trivial tensors:
            do j=2,10
             call tensor_block_destroy(tbb_entry%tens_blck,i); write(jo_cp,*) 'Tensor block destruction: ',i
             err_code=tbb%search(dict_add_if_not_found,tens_key_cmp,key(j),tbb_entry,value_out=uptr); write(jo_cp,*) 'TBB entry search: ',err_code
             select type (uptr)
             type is (tbb_entry_t)
              err_code=aar_register(key(j),targ_p); write(jo_cp,*) 'AAR entry created: ',err_code
              targ_p%tens_blck_f=>uptr%tens_blck
             end select
            enddo
 !Fill in tensor instructions in ETIQ:
            etiq%scheduled=0
  !ETIQ(1):
            k=1
            etiq%eti(k)%instr_code=instr_tensor_init; etiq%eti(k)%data_kind='r8'
            allocate(etiq%eti(k)%instr_aux(0:1+3*4))
            etiq%eti(k)%instr_aux(0:1+3*4)=(/14,  12, 15,55,15,55, 15,55,15,55, 0,0,0,0/)
            etiq%eti(k)%instr_cu=cu_t(DEV_HOST,1)
            etiq%eti(k)%args_ready=B'1111111111111111111111111111111'
            etiq%eti(k)%tens_op(0)%tens_blck_id=key(2)
            err_code=aar_register(etiq%eti(k)%tens_op(0)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(0)%op_aar_entry=>targ_p
            etiq%eti(k)%tens_op(3)%tens_blck_id=key(0)
            err_code=aar_register(etiq%eti(k)%tens_op(3)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(3)%op_aar_entry=>targ_p
            etiq%eti(k)%instr_handle=-1
            etiq%eti(k)%instr_status=instr_ready_to_exec
            etiq%scheduled=etiq%scheduled+1
  !ETIQ(2):
            k=2
            etiq%eti(k)%instr_code=instr_tensor_init; etiq%eti(k)%data_kind='r8'
            allocate(etiq%eti(k)%instr_aux(0:1+3*4))
            etiq%eti(k)%instr_aux(0:1+3*4)=(/14,  12, 15,55,15,55, 15,55,15,55, 0,0,0,0/)
            etiq%eti(k)%instr_cu=cu_t(DEV_HOST,0)
            etiq%eti(k)%args_ready=B'1111111111111111111111111111111'
            etiq%eti(k)%tens_op(0)%tens_blck_id=key(3)
            err_code=aar_register(etiq%eti(k)%tens_op(0)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(0)%op_aar_entry=>targ_p
            etiq%eti(k)%tens_op(3)%tens_blck_id=key(0)
            err_code=aar_register(etiq%eti(k)%tens_op(3)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(3)%op_aar_entry=>targ_p
            etiq%eti(k)%instr_handle=-1
            etiq%eti(k)%instr_status=instr_ready_to_exec
            etiq%scheduled=etiq%scheduled+1
  !ETIQ(3):
            k=3
            etiq%eti(k)%instr_code=instr_tensor_init; etiq%eti(k)%data_kind='r8'
            allocate(etiq%eti(k)%instr_aux(0:1+3*4))
            etiq%eti(k)%instr_aux(0:1+3*4)=(/14,  12, 15,55,15,55, 15,55,15,55, 0,0,0,0/)
            etiq%eti(k)%instr_cu=cu_t(DEV_HOST,1)
            etiq%eti(k)%args_ready=B'1111111111111111111111111111111'
            etiq%eti(k)%tens_op(0)%tens_blck_id=key(4)
            err_code=aar_register(etiq%eti(k)%tens_op(0)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(0)%op_aar_entry=>targ_p
            etiq%eti(k)%tens_op(3)%tens_blck_id=key(0)
            err_code=aar_register(etiq%eti(k)%tens_op(3)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(3)%op_aar_entry=>targ_p
            etiq%eti(k)%instr_handle=-1
            etiq%eti(k)%instr_status=instr_ready_to_exec
            etiq%scheduled=etiq%scheduled+1
  !ETIQ(4):
            k=4
            etiq%eti(k)%instr_code=instr_tensor_init; etiq%eti(k)%data_kind='r8'
            allocate(etiq%eti(k)%instr_aux(0:1+3*5))
            etiq%eti(k)%instr_aux(0:1+3*5)=(/17,  15, 15,10,55,45,25, 15,10,55,45,25, 0,0,0,0,0/)
            etiq%eti(k)%instr_cu=cu_t(DEV_HOST,0)
            etiq%eti(k)%args_ready=B'1111111111111111111111111111111'
            etiq%eti(k)%tens_op(0)%tens_blck_id=key(5)
            err_code=aar_register(etiq%eti(k)%tens_op(0)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(0)%op_aar_entry=>targ_p
            etiq%eti(k)%tens_op(3)%tens_blck_id=key(1)
            err_code=aar_register(etiq%eti(k)%tens_op(3)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(3)%op_aar_entry=>targ_p
            etiq%eti(k)%instr_handle=-1
            etiq%eti(k)%instr_status=instr_ready_to_exec
            etiq%scheduled=etiq%scheduled+1
  !ETIQ(5):
            k=5
            etiq%eti(k)%instr_code=instr_tensor_init; etiq%eti(k)%data_kind='r8'
            allocate(etiq%eti(k)%instr_aux(0:1+3*5))
            etiq%eti(k)%instr_aux(0:1+3*5)=(/17,  15, 15,10,55,45,25, 15,10,55,45,25, 0,0,0,0,0/)
            etiq%eti(k)%instr_cu=cu_t(DEV_HOST,0)
            etiq%eti(k)%args_ready=B'1111111111111111111111111111111'
            etiq%eti(k)%tens_op(0)%tens_blck_id=key(6)
            err_code=aar_register(etiq%eti(k)%tens_op(0)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(0)%op_aar_entry=>targ_p
            etiq%eti(k)%tens_op(3)%tens_blck_id=key(1)
            err_code=aar_register(etiq%eti(k)%tens_op(3)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(3)%op_aar_entry=>targ_p
            etiq%eti(k)%instr_handle=-1
            etiq%eti(k)%instr_status=instr_ready_to_exec
            etiq%scheduled=etiq%scheduled+1
  !ETIQ(6):
            k=6
            etiq%eti(k)%instr_code=instr_tensor_copy; etiq%eti(k)%data_kind='r8'
            allocate(etiq%eti(k)%instr_aux(0:1+5))
            etiq%eti(k)%instr_aux(0:1+5)=(/7,  5, 5,4,3,2,1/)
            etiq%eti(k)%instr_cu=cu_t(DEV_HOST,1)
            etiq%eti(k)%args_ready=B'1111111111111111111111111111111'
            etiq%eti(k)%tens_op(0)%tens_blck_id=key(7)
            err_code=aar_register(etiq%eti(k)%tens_op(0)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(0)%op_aar_entry=>targ_p
            etiq%eti(k)%tens_op(1)%tens_blck_id=key(5)
            err_code=aar_register(etiq%eti(k)%tens_op(1)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(1)%op_aar_entry=>targ_p
            etiq%eti(k)%instr_handle=-1
            etiq%eti(k)%instr_status=instr_ready_to_exec
            etiq%scheduled=etiq%scheduled+1
  !ETIQ(7):
            k=7
            etiq%eti(k)%instr_code=instr_tensor_copy; etiq%eti(k)%data_kind='r8'
            allocate(etiq%eti(k)%instr_aux(0:1+5))
            etiq%eti(k)%instr_aux(0:1+5)=(/7,  5, 3,1,5,2,4/)
            etiq%eti(k)%instr_cu=cu_t(DEV_HOST,0)
            etiq%eti(k)%args_ready=B'1111111111111111111111111111111'
            etiq%eti(k)%tens_op(0)%tens_blck_id=key(8)
            err_code=aar_register(etiq%eti(k)%tens_op(0)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(0)%op_aar_entry=>targ_p
            etiq%eti(k)%tens_op(1)%tens_blck_id=key(6)
            err_code=aar_register(etiq%eti(k)%tens_op(1)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(1)%op_aar_entry=>targ_p
            etiq%eti(k)%instr_handle=-1
            etiq%eti(k)%instr_status=instr_ready_to_exec
            etiq%scheduled=etiq%scheduled+1
  !ETIQ(8):
            k=8
            etiq%eti(k)%instr_code=instr_tensor_copy; etiq%eti(k)%data_kind='r8'
            allocate(etiq%eti(k)%instr_aux(0:1+5))
            etiq%eti(k)%instr_aux(0:1+5)=(/7,  5, 5,4,3,2,1/)
            etiq%eti(k)%instr_cu=cu_t(DEV_HOST,1)
            etiq%eti(k)%args_ready=B'1111111111111111111111111111111'
            etiq%eti(k)%tens_op(0)%tens_blck_id=key(9)
            err_code=aar_register(etiq%eti(k)%tens_op(0)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(0)%op_aar_entry=>targ_p
            etiq%eti(k)%tens_op(1)%tens_blck_id=key(5)
            err_code=aar_register(etiq%eti(k)%tens_op(1)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(1)%op_aar_entry=>targ_p
            etiq%eti(k)%instr_handle=-1
            etiq%eti(k)%instr_status=instr_ready_to_exec
            etiq%scheduled=etiq%scheduled+1
  !ETIQ(9):
            k=9
            etiq%eti(k)%instr_code=instr_tensor_copy; etiq%eti(k)%data_kind='r8'
            allocate(etiq%eti(k)%instr_aux(0:1+5))
            etiq%eti(k)%instr_aux(0:1+5)=(/7,  5, 3,1,5,2,4/)
            etiq%eti(k)%instr_cu=cu_t(DEV_HOST,0)
            etiq%eti(k)%args_ready=B'1111111111111111111111111111111'
            etiq%eti(k)%tens_op(0)%tens_blck_id=key(10)
            err_code=aar_register(etiq%eti(k)%tens_op(0)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(0)%op_aar_entry=>targ_p
            etiq%eti(k)%tens_op(1)%tens_blck_id=key(6)
            err_code=aar_register(etiq%eti(k)%tens_op(1)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(1)%op_aar_entry=>targ_p
            etiq%eti(k)%instr_handle=-1
            etiq%eti(k)%instr_status=instr_ready_to_exec
            etiq%scheduled=etiq%scheduled+1
  !ETIQ(10):
            k=10
            etiq%eti(k)%instr_code=instr_tensor_contract; etiq%eti(k)%data_kind='r8'
            allocate(etiq%eti(k)%instr_aux(0:1+10))
            etiq%eti(k)%instr_aux(0:1+10)=(/12,  10, 3,-2,2,-4,-5, 1,-2,4,-4,-5/)
            etiq%eti(k)%instr_cu=cu_t(DEV_HOST,0)
            etiq%eti(k)%args_ready=B'1111111111111111111111111111111'
            etiq%eti(k)%tens_op(0)%tens_blck_id=key(2)
            err_code=aar_register(etiq%eti(k)%tens_op(0)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(0)%op_aar_entry=>targ_p
            etiq%eti(k)%tens_op(1)%tens_blck_id=key(5)
            err_code=aar_register(etiq%eti(k)%tens_op(1)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(1)%op_aar_entry=>targ_p
            etiq%eti(k)%tens_op(2)%tens_blck_id=key(6)
            err_code=aar_register(etiq%eti(k)%tens_op(2)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(2)%op_aar_entry=>targ_p
            etiq%eti(k)%instr_handle=-1
            etiq%eti(k)%instr_status=instr_ready_to_exec
            etiq%scheduled=etiq%scheduled+1
  !ETIQ(11):
            k=11
            etiq%eti(k)%instr_code=instr_tensor_contract; etiq%eti(k)%data_kind='r8'
            allocate(etiq%eti(k)%instr_aux(0:1+10+1+1))
            etiq%eti(k)%instr_aux(0:1+10+1+1)=(/14,  10, -4,-2,4,-1,1,-4,-2,3,-1,2,  1, COPY_BACK/)
            etiq%eti(k)%instr_cu=cu_t(DEV_NVIDIA_GPU,0)
            etiq%eti(k)%args_ready=B'1111111111111111111111111111111'
            etiq%eti(k)%tens_op(0)%tens_blck_id=key(3)
            err_code=aar_register(etiq%eti(k)%tens_op(0)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(0)%op_aar_entry=>targ_p
            etiq%eti(k)%tens_op(1)%tens_blck_id=key(7)
            err_code=aar_register(etiq%eti(k)%tens_op(1)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(1)%op_aar_entry=>targ_p
            etiq%eti(k)%tens_op(2)%tens_blck_id=key(8)
            err_code=aar_register(etiq%eti(k)%tens_op(2)%tens_blck_id,targ_p); write(jo_cp,*) 'AAR entry search: ',err_code
            etiq%eti(k)%tens_op(2)%op_aar_entry=>targ_p
            etiq%eti(k)%instr_handle=-1
            etiq%eti(k)%instr_status=instr_ready_to_exec
            etiq%scheduled=etiq%scheduled+1
!$OMP FLUSH

            call set_transpose_algorithm(EFF_TRN_ON)
            call set_matmult_algorithm(BLAS_OFF)
            call gpu_set_matmult_algorithm(BLAS_OFF)
!$OMP FLUSH

 !Enqueue tensor instructions to NVCU:
            etiq_nvcu%etiq_entry(0)=11
            etiq_nvcu%te_conf(0)=te_conf_t(etiq%eti(etiq_nvcu%etiq_entry(0))%instr_cu,16,16,1)
            etiq_nvcu%scheduled=etiq_nvcu%scheduled+1
            etiq%eti(etiq_nvcu%etiq_entry(0))%instr_status=instr_scheduled
 !Enqueue tensor instructions to STCU:
            etiq_stcu%etiq_entry(9)=10
            etiq_stcu%te_conf(9)=te_conf_t(etiq%eti(etiq_stcu%etiq_entry(9))%instr_cu,1,8,1)
            etiq_stcu%scheduled=etiq_stcu%scheduled+1
            etiq%eti(etiq_stcu%etiq_entry(9))%instr_status=instr_scheduled
            do k=8,0,-1
             etiq_stcu%etiq_entry(k)=1+k
             etiq_stcu%te_conf(k)=te_conf_t(etiq%eti(etiq_stcu%etiq_entry(k))%instr_cu,1,4,1)
             etiq_stcu%scheduled=etiq_stcu%scheduled+1
             etiq%eti(etiq_stcu%etiq_entry(k))%instr_status=instr_scheduled
            enddo

!$OMP FLUSH
            write(jo_cp,*)'Instruction(s) scheduled!'
 !Wait for completion:
            i=etiq%scheduled; kf=0
            do while(i.gt.0)
             do k=1,etiq%scheduled
              if(kf.ge.3) j=nvcu_task_status(0)
!$OMP ATOMIC READ
              j=etiq%eti(k)%instr_status
              if(j.eq.instr_completed) then
!$OMP ATOMIC READ
               tm0=etiq%eti(k)%time_issued
!$OMP ATOMIC READ
               tm1=etiq%eti(k)%time_completed
               write(jo_cp,'("Instruction ",i2," completed succefully: ",D25.15)') k,tm1-tm0
!$OMP ATOMIC WRITE
               etiq%eti(k)%instr_status=instr_dead
               i=i-1
               if(k.eq.6.or.k.eq.7) kf=kf+1
               if(kf.eq.2) then
                j=nvcu_execute_eti(0); kf=kf+1; write(jo_cp,'("Instr#11 went to GPU: ",i12)') j
               endif
              elseif(j.le.0) then
               write(jo_cp,'("Instruction ",i2," failed: ",i11)') k,j
!$OMP ATOMIC WRITE
               etiq%eti(k)%instr_status=instr_dead
               i=i-1
              endif
             enddo
            enddo
            write(jo_cp,*)'Instruction(s) completed: ',i
 !Clean up GPU stuff:
            j=nvcu_task_cleanup(0); write(jo_cp,'("Instr#11 cleanup status: ",i12)') j
 !Destroy keys:
            do j=0,10; err_code=key(j)%destroy(); enddo
 !Destroy AAR:
            err_code=aar%reset()
            err_code=aar%traverse_subtree(tree_volume,tree_height); write(jo_cp,*) tree_volume,tree_height
            err_code=aar%destroy(destruct_key_func=cp_destructor)
            print *,'AAR destroyed: ',err_code
 !Destroy TBB:
            err_code=tbb%reset()
            err_code=tbb%traverse_subtree(tree_volume,tree_height); write(jo_cp,*) tree_volume,tree_height
            err_code=tbb%destroy(destruct_key_func=cp_destructor,destruct_val_func=cp_destructor)
            print *,'TBB destroyed: ',err_code
!DEBUG end.
            
            exit master_life
           enddo master_life
 !Exit MT life:
!$OMP ATOMIC WRITE
           etiq_stcu%scheduled=-1 !this will make the slave threads exit
!$OMP FLUSH
!Slave threads:
          else !thread_num=1: leading slave thread (LST)
!$OMP ATOMIC WRITE
           stcu_error=0 !set by LST, can be used by MT
!$OMP FLUSH(stcu_error)
           slave_life: do !LST life cycle
!$OMP FLUSH
            if(stcu_error.eq.0) then
!$OMP ATOMIC READ
             n=etiq_stcu%scheduled
             if(n.lt.0) exit slave_life !negative etiq_stcu%scheduled causes STCU life termination
 !Retrieve the next ETI from ETIQ_STCU (if present):
             stcu_base_ip=etiq_stcu%ip !current instruction (if any)
             l=etiq_stcu%etiq_entry(stcu_base_ip) !ETIQ entry (if >0)
             if(l.gt.0.and.l.le.etiq%depth) then !instruction is there
              if(etiq%eti(l)%instr_status.eq.instr_scheduled) then !instruction has not been processed before
               stcu_num_units=etiq_stcu%te_conf(stcu_base_ip)%cu_id%unit_number+1 !LST configures STCU
!$OMP FLUSH(stcu_num_units)
               err_code=0
               if(stcu_num_units.eq.1) then !only one STCU: SIMD execution
                stcu_my_ip=stcu_base_ip; stcu_my_eti=l
                if(verbose) write(jo_cp,'("#DEBUG(c_process::c_proc_life): STCU 0",": IP ",i5,": ETI #",i7,": thread_count=",i3)') stcu_my_ip,stcu_my_eti,etiq_stcu%te_conf(stcu_my_ip)%num_workers !debug
                call omp_set_num_threads(etiq_stcu%te_conf(stcu_my_ip)%num_workers)
                err_code=stcu_execute_eti(stcu_my_eti)
                if(err_code.eq.0) then
                 etiq_stcu%ip=etiq_stcu%ip+1
                 if(etiq_stcu%ip.ge.etiq_stcu%depth) etiq_stcu%ip=etiq_stcu%ip-etiq_stcu%depth
                else
                 write(jo_cp,'("#ERROR(c_process::c_proc_life): STCU 0/0: execution error ",i11," for ETI #",i7)') err_code,stcu_my_eti
!$OMP ATOMIC WRITE
                 stcu_error=1 !error during ETI execution in one or more STCU
!$OMP FLUSH(stcu_error)
                endif
               elseif(stcu_num_units.gt.1.and.stcu_num_units.le.stcu_max_units) then !multiple STCU: MIMD execution
                stcu_mimd_my_pass=0; stcu_mimd_max_pass=1
!$OMP PARALLEL NUM_THREADS(stcu_num_units) FIRSTPRIVATE(stcu_base_ip,stcu_mimd_my_pass,err_code) &
!$OMP          PRIVATE(thread_num,stcu_my_ip,stcu_my_eti,i,j,k) DEFAULT(SHARED)
                thread_num=omp_get_thread_num()
                if(thread_num.eq.0) then
                 if(verbose) write(jo_cp,'("#DEBUG(c_process::c_proc_life): STCU: ",i3," units created.")') stcu_num_units !debug
                endif
!$OMP BARRIER
                if(omp_get_num_threads().eq.stcu_num_units) then
                 mimd_loop: do
                  stcu_mimd_my_pass=stcu_mimd_my_pass+1
                  stcu_my_ip=stcu_base_ip+thread_num; if(stcu_my_ip.ge.etiq_stcu%depth) stcu_my_ip=stcu_my_ip-etiq_stcu%depth
                  stcu_my_eti=etiq_stcu%etiq_entry(stcu_my_ip)
                  if(verbose) write(jo_cp,'("#DEBUG(c_process::c_proc_life): STCU ",i2,"/",i2,": Pass ",i5,": IP ",i5,": ETI #",i7,": thread_count=",i3)') &
                   thread_num,stcu_num_units,stcu_mimd_my_pass,stcu_my_ip,stcu_my_eti,etiq_stcu%te_conf(stcu_my_ip)%num_workers !debug
                  call omp_set_num_threads(etiq_stcu%te_conf(stcu_my_ip)%num_workers)
                  err_code=stcu_execute_eti(stcu_my_eti)
                  if(err_code.ne.0) then
                   write(jo_cp,'("#ERROR(c_process::c_proc_life): STCU ",i2,": execution error ",i11," for ETI #",i7)') thread_num,err_code,stcu_my_eti
!$OMP ATOMIC WRITE
                   stcu_error=2 !error during ETI execution in one or more STCU
                  endif
!$OMP FLUSH(stcu_error)
                  if(stcu_error.eq.0) then !check whether the next ETIs are MIMD of the same width
                   stcu_base_ip=stcu_base_ip+stcu_num_units
                   if(stcu_base_ip.ge.etiq_stcu%depth) stcu_base_ip=stcu_base_ip-etiq_stcu%depth
                   j=0
                   wait_next_mimd: do while(j.ge.0) !check if the next ETI(s) are configured for the same MIMD execution
!$OMP FLUSH
!$OMP CRITICAL(test_next_mimd)
                    if(stcu_mimd_my_pass.lt.iabs(stcu_mimd_max_pass)) then
                     stcu_mimd_my_pass=-stcu_mimd_my_pass
                    else
                     if(stcu_mimd_max_pass.gt.0) then !I am the first here: check next MIMD batch
                      k=etiq_stcu%etiq_entry(stcu_base_ip)
                      if(k.gt.0.and.k.le.etiq%depth) then !valid ETIQ entry
                       if(etiq%eti(k)%instr_status.eq.instr_scheduled) then !new instruction scheduled
                        if(etiq_stcu%te_conf(stcu_base_ip)%cu_id%unit_number+1.eq.stcu_num_units) then
                         stcu_mimd_max_pass=stcu_mimd_max_pass+1
                         stcu_mimd_my_pass=-stcu_mimd_my_pass
                        else
                         stcu_mimd_max_pass=-stcu_mimd_max_pass; j=-2 !STCU conf. changed
                        endif
                       endif
                      endif
                     else !I am not the first here and MIMD batches are over
                      j=-3
                     endif
                    endif
!$OMP END CRITICAL(test_next_mimd)
                    if(j.ge.0) then
                     if(stcu_mimd_my_pass.lt.0) then
                      stcu_mimd_my_pass=-stcu_mimd_my_pass
                      exit wait_next_mimd
                     else
!$OMP ATOMIC READ
                      j=etiq_stcu%scheduled
                     endif
                    endif
                   enddo wait_next_mimd
                   if(j.lt.0) exit mimd_loop
                  else
                   exit mimd_loop
                  endif
                 enddo mimd_loop
                else
                 write(jo_cp,'("#ERROR(c_process::c_proc_life): STCU: Invalid number of units created: ",i11,1x,i11)') omp_get_num_threads(),stcu_num_units
!$OMP ATOMIC WRITE
                 stcu_error=3 !invalid number of STCU units created
                endif
!$OMP END PARALLEL
                if(stcu_error.eq.0) then !LST only
!$ATOMIC UPDATE
                 etiq_stcu%ip=mod(etiq_stcu%ip+stcu_num_units*iabs(stcu_mimd_max_pass),etiq_stcu%depth)
                endif
               else
                write(jo_cp,'("#ERROR(c_process::c_proc_life): STCU: Invalid number of units requested: ",i11)') stcu_num_units
!$OMP ATOMIC WRITE
                stcu_error=4 !invalid number of STCU units requested
!$OMP FLUSH(stcu_error)
               endif
              endif
             endif
            else
             write(jo_cp,'("#ERROR(c_process::c_proc_life): LST life: STCU error ",i11)') stcu_error
             !`do not exit: handle the STCU error by the LST
!DEBUG begin:
!$OMP ATOMIC WRITE
             ierr=1 !STCU error occured
!$OMP FLUSH(ierr)
             exit slave_life
!DEBUG end.
            endif
           enddo slave_life
           if(verbose) write(jo_cp,'("#DEBUG(c_process::c_proc_life): STCU life is over with status: ",i11)') stcu_error !debug
          endif !MT/LST
         else
!$OMP MASTER
          write(jo_cp,'("#ERROR(c_process::c_proc_life): invalid initial number of threads: ",i6)') n
          ierr=-2
!$OMP END MASTER
         endif
!$OMP END PARALLEL
        else
         write(jo_cp,'("#ERROR(c_process::c_proc_life): invalid max number of threads: ",i6,1x,i6)') n,max_threads
         ierr=-1
        endif
        if(ierr.ne.0) then; call c_proc_quit(998); return; endif
#else
        write(jo_cp,'("#FATAL(c_process::c_proc_life): CP cannot function without OpenMP!")')
        call c_proc_quit(999); return
#endif
!-------------------
        write(jo_cp,'("#MSG(c_process::c_proc_life): Cleaning ... ")',advance='no')
        call c_proc_quit(0); if(ierr.eq.0) write(jo_cp,'("Ok")')
        return

        contains

         integer function stcu_execute_eti(eti_loc) !STCU ETI execution workflow: ST only (thread private)
!STCU execution is blocking, but multiple ETI can be executed in parallel (MIMD);
!Once each issued ETI has completed (or failed), its status (time) is updated in ETIQ.
!Format of the %instr_aux() array from <tens_instr_t> is given below:
! * Total number of elements in %instr_aux();
!  * Number of elements in the 1st field of %instr_aux();
!   * 1st field of %instr_aux();
!  * Number of elements in the 2nd field of %instr_aux();
!   * 2nd field of %instr_aux();
!  etc.
         implicit none
         integer, intent(in):: eti_loc !entry number in <etiq> (ETI location in ETIQ)
         type(tens_instr_t), pointer:: my_eti
         type(tensor_block_t), pointer:: dtens_,ltens_,rtens_,stens_
         character(max_shape_str_len):: jtsss
         integer:: j0,j1,jl,ju,ier
         integer(LONGINT):: jdiff
         real(8):: jval
         logical:: jres
         stcu_execute_eti=0
         if(eti_loc.gt.0.and.eti_loc.le.etiq%depth) then
          my_eti=>etiq%eti(eti_loc)
          if(my_eti%instr_status.eq.instr_scheduled) then
!$OMP ATOMIC WRITE
           my_eti%instr_status=instr_issued
!$OMP ATOMIC WRITE
           my_eti%time_issued=thread_wtime()
 !Associate arguments:
           if(associated(my_eti%tens_op(0)%op_aar_entry)) then
            if(associated(my_eti%tens_op(0)%op_aar_entry%tens_blck_f)) then
             dtens_=>my_eti%tens_op(0)%op_aar_entry%tens_blck_f
            else
             dtens_=>NULL()
            endif
           else
            dtens_=>NULL()
           endif
           if(associated(my_eti%tens_op(1)%op_aar_entry)) then
            if(associated(my_eti%tens_op(1)%op_aar_entry%tens_blck_f)) then
             ltens_=>my_eti%tens_op(1)%op_aar_entry%tens_blck_f
            else
             ltens_=>NULL()
            endif
           else
            ltens_=>NULL()
           endif
           if(associated(my_eti%tens_op(2)%op_aar_entry)) then
            if(associated(my_eti%tens_op(2)%op_aar_entry%tens_blck_f)) then
             rtens_=>my_eti%tens_op(2)%op_aar_entry%tens_blck_f
            else
             rtens_=>NULL()
            endif
           else
            rtens_=>NULL()
           endif
           if(associated(my_eti%tens_op(3)%op_aar_entry)) then
            if(associated(my_eti%tens_op(3)%op_aar_entry%tens_blck_f)) then
             stens_=>my_eti%tens_op(3)%op_aar_entry%tens_blck_f
            else
             stens_=>NULL()
            endif
           else
            stens_=>NULL()
           endif
 !Execute:
           select case(my_eti%instr_code)
           case(instr_tensor_init)
  !Create/initialize a tensor block:
            if(associated(dtens_).and..not.(associated(ltens_).or.associated(rtens_))) then
             if(tensor_block_layout(dtens_,ier).eq.not_allocated) then !create & init
              if(allocated(my_eti%instr_aux)) then
               jl=lbound(my_eti%instr_aux,1); ju=ubound(my_eti%instr_aux,1)
               if(ju.ge.jl+1.and.ju-jl+1.eq.my_eti%instr_aux(jl)) then !total number of elements in <instr_aux>
                jl=jl+1; j0=my_eti%instr_aux(jl) !Field#1: {dims[1:rank],divs[1:rank],grps[1:rank]}
                if(j0.ge.0.and.mod(j0,3).eq.0.and.jl+j0.le.ju) then !j0 equals tensor_rank*3
                 if(j0.gt.0) then !true tensor
                  j0=j0/3 !now j0 = tensor rank
                  call tensor_shape_str_create(j0,my_eti%instr_aux(jl+1:jl+j0),jtsss,j1,ier, &
                        divs=my_eti%instr_aux(jl+j0+1:jl+j0*2),grps=my_eti%instr_aux(jl+j0*2+1:jl+j0*3))
                 else !scalar
                  call tensor_shape_str_create(j0,my_eti%instr_aux(jl:),jtsss,j1,ier)
                 endif
                 if(ier.eq.0) then
                  if(associated(stens_)) then
                   if(stens_%tensor_shape%num_dim.eq.0) then !must be a scalar
                    call tensor_block_create(jtsss(1:j1),my_eti%data_kind,dtens_,ier,val_c8=stens_%scalar_value)
                    if(ier.ne.0) stcu_execute_eti=1
                   else
                    stcu_execute_eti=2
                   endif
                  else
                   call tensor_block_create(jtsss(1:j1),my_eti%data_kind,dtens_,ier)
                   if(ier.ne.0) stcu_execute_eti=3
                  endif
                 else
                  stcu_execute_eti=4
                 endif
                else
                 stcu_execute_eti=5
                endif
               else
                stcu_execute_eti=6
               endif
              else
               stcu_execute_eti=7
              endif
             else !only init
              if(associated(stens_)) then
               if(stens_%tensor_shape%num_dim.eq.0) then !must be a scalar
                call tensor_block_init(my_eti%data_kind,dtens_,ier,val_c8=stens_%scalar_value)
                if(ier.ne.0) stcu_execute_eti=8
               else
                stcu_execute_eti=9
               endif
              else
               call tensor_block_init(my_eti%data_kind,dtens_,ier)
               if(ier.ne.0) stcu_execute_eti=10
              endif
             endif
            else
             stcu_execute_eti=11
            endif
           case(instr_tensor_norm1)
  !Compute the 1-norm of a tensor block:
            if(associated(dtens_).and.associated(ltens_).and.(.not.associated(rtens_))) then
             if(tensor_block_layout(ltens_,ier).ne.not_allocated) then
              if(dtens_%tensor_shape%num_dim.eq.0) then !must be a scalar
               jval=tensor_block_norm1(ltens_,ier,my_eti%data_kind)
               if(ier.eq.0) then
                dtens_%scalar_value=cmplx(jval,0d0,8)
               else
                stcu_execute_eti=12
               endif
              else
               stcu_execute_eti=13
              endif
             else
              stcu_execute_eti=14
             endif
            else
             stcu_execute_eti=15
            endif
           case(instr_tensor_norm2)
  !Compute the 2-norm of a tensor block:
            if(associated(dtens_).and.associated(ltens_).and.(.not.associated(rtens_))) then
             if(tensor_block_layout(ltens_,ier).ne.not_allocated) then
              if(dtens_%tensor_shape%num_dim.eq.0) then !must be a scalar
               jval=tensor_block_norm2(ltens_,ier,my_eti%data_kind) !squared 2-norm
               if(ier.eq.0) then
                dtens_%scalar_value=cmplx(dsqrt(jval),0d0,8)
               else
                stcu_execute_eti=16
               endif
              else
               stcu_execute_eti=17
              endif
             else
              stcu_execute_eti=18
             endif
            else
             stcu_execute_eti=19
            endif
           case(instr_tensor_min)
  !Return the min by absolute value element of a tensor block:
            if(associated(dtens_).and.associated(ltens_).and.(.not.associated(rtens_))) then
             if(tensor_block_layout(ltens_,ier).ne.not_allocated) then
              if(dtens_%tensor_shape%num_dim.eq.0) then !must be a scalar
               jval=tensor_block_min(ltens_,ier,my_eti%data_kind)
               if(ier.eq.0) then
                dtens_%scalar_value=cmplx(jval,0d0,8)
               else
                stcu_execute_eti=20
               endif
              else
               stcu_execute_eti=21
              endif
             else
              stcu_execute_eti=22
             endif
            else
             stcu_execute_eti=23
            endif
           case(instr_tensor_max)
  !Return the max by absolute value element of a tensor block:
            if(associated(dtens_).and.associated(ltens_).and.(.not.associated(rtens_))) then
             if(tensor_block_layout(ltens_,ier).ne.not_allocated) then
              if(dtens_%tensor_shape%num_dim.eq.0) then !must be a scalar
               jval=tensor_block_max(ltens_,ier,my_eti%data_kind)
               if(ier.eq.0) then
                dtens_%scalar_value=cmplx(jval,0d0,8)
               else
                stcu_execute_eti=24
               endif
              else
               stcu_execute_eti=25
              endif
             else
              stcu_execute_eti=26
             endif
            else
             stcu_execute_eti=27
            endif
           case(instr_tensor_scale)
  !Multiply a tensor block by a scalar:
            if(associated(dtens_).and.(.not.associated(ltens_)).and.(.not.associated(rtens_)).and.associated(stens_)) then
             if(tensor_block_layout(dtens_,ier).ne.not_allocated) then
              if(dtens_%tensor_shape%num_dim.ge.0) then
               call tensor_block_scale(dtens_,stens_%scalar_value,ier)
               if(ier.ne.0) stcu_execute_eti=28
              else
               stcu_execute_eti=29
              endif
             else
              stcu_execute_eti=30
             endif
            else
             stcu_execute_eti=31
            endif
           case(instr_tensor_slice)
  !Return a slice of a tensor block:
           case(instr_tensor_insert)
  !Insert a slice in a tensor block:
           case(instr_tensor_trace)
  !Partial/full trace in a tensor block:
           case(instr_tensor_copy)
  !Copy/permute a tensor block into another tensor block:
            if(associated(dtens_).and.associated(ltens_).and.(.not.associated(rtens_))) then
             if(tensor_block_layout(ltens_,ier).ne.not_allocated) then
              if(allocated(my_eti%instr_aux)) then
               jl=lbound(my_eti%instr_aux,1); ju=ubound(my_eti%instr_aux,1)
               if(ju.ge.jl+1.and.ju-jl+1.eq.my_eti%instr_aux(jl)) then
                jl=jl+1; j0=my_eti%instr_aux(jl) !tensor block rank: Field#1: trn[1:rank] is the required permutation
                if(j0.ge.0.and.jl+j0.le.ju) then
                 call tensor_block_copy(ltens_,dtens_,ier,my_eti%instr_aux(jl:)) !instr_aux(jl) will be ignored
                 if(ier.ne.0) stcu_execute_eti=32
                else
                 stcu_execute_eti=33
                endif
               else
                stcu_execute_eti=34
               endif
              else
               call tensor_block_copy(ltens_,dtens_,ier)
               if(ier.ne.0) stcu_execute_eti=35
              endif
             else
              stcu_execute_eti=36
             endif
            else
             stcu_execute_eti=37
            endif
           case(instr_tensor_add)
  !Add a tensor block to another tensor block:
            if(associated(dtens_).and.associated(ltens_).and.(.not.associated(rtens_))) then
             if(tensor_block_layout(dtens_,ier).ne.not_allocated.and.tensor_block_layout(ltens_,ier).ne.not_allocated) then
              if(associated(stens_)) then !scaling factor
               call tensor_block_add(dtens_,ltens_,ier,stens_%scalar_value,my_eti%data_kind)
               if(ier.ne.0) stcu_execute_eti=38
              else
               call tensor_block_add(dtens_,ltens_,ier,data_kind=my_eti%data_kind)
               if(ier.ne.0) stcu_execute_eti=39
              endif
             else
              stcu_execute_eti=40
             endif
            else
             stcu_execute_eti=41
            endif
           case(instr_tensor_cmp)
  !Compare two tensor blocks (scalar dtens_ will contain the number of elements that differ):
            if(associated(dtens_).and.associated(ltens_).and.associated(rtens_)) then
             if(tensor_block_layout(dtens_,ier).ne.not_allocated.and. &
                tensor_block_layout(ltens_,ier).ne.not_allocated.and.tensor_block_layout(rtens_,ier).ne.not_allocated) then
              if(dtens_%tensor_shape%num_dim.eq.0) then
               jdiff=0
               if(associated(stens_)) then !real part is the comparison threshold
                jres=tensor_block_cmp(ltens_,rtens_,ier,my_eti%data_kind,.true.,real(stens_%scalar_value,8),jdiff)
                if(ier.ne.0) stcu_execute_eti=42
               else
                jres=tensor_block_cmp(ltens_,rtens_,ier,my_eti%data_kind,rel=.true.,diff_count=jdiff)
                if(ier.ne.0) stcu_execute_eti=43
               endif
               dtens_%scalar_value=cmplx(real(jdiff,8),0d0,8)
              else
               stcu_execute_eti=44
              endif
             else
              stcu_execute_eti=45
             endif
            else
             stcu_execute_eti=46
            endif
           case(instr_tensor_contract)
  !Contract two tensor blocks and accumulate the result into another tensor block:
            if(associated(dtens_).and.associated(ltens_).and.associated(rtens_)) then
             if(tensor_block_layout(dtens_,ier).ne.not_allocated.and. &
                tensor_block_layout(ltens_,ier).ne.not_allocated.and.tensor_block_layout(rtens_,ier).ne.not_allocated) then
              if(allocated(my_eti%instr_aux)) then
               jl=lbound(my_eti%instr_aux,1); ju=ubound(my_eti%instr_aux,1)
               if(ju.ge.jl+1.and.ju-jl+1.eq.my_eti%instr_aux(jl)) then
                jl=jl+1; j0=my_eti%instr_aux(jl) !contraction pattern length: Field#1: cptrn[1:length] is the contraction pattern
                if(j0.ge.0.and.jl+j0.le.ju) then
                 if(j0.gt.0) then !at least one true tensor present
                  call tensor_block_contract(my_eti%instr_aux(jl+1:jl+j0),ltens_,rtens_,dtens_,ier,my_eti%data_kind)
                 else !all three tensors are scalars
                  call tensor_block_contract(my_eti%instr_aux(jl:),ltens_,rtens_,dtens_,ier)
                 endif
                 if(ier.ne.0) stcu_execute_eti=47
                else
                 stcu_execute_eti=48
                endif
               else
                stcu_execute_eti=49
               endif
              else
               stcu_execute_eti=50
              endif
             else
              stcu_execute_eti=51
             endif
            else
             stcu_execute_eti=52
            endif
           case default
            stcu_execute_eti=997
           end select
!$OMP ATOMIC WRITE
           my_eti%time_completed=thread_wtime()
          else
           stcu_execute_eti=998
          endif
          if(stcu_execute_eti.eq.0) then
           j0=eti_mark_aar_used(eti_loc) !mark all tensor arguments as used
!$OMP ATOMIC WRITE
           my_eti%instr_status=instr_completed
          else
!$OMP ATOMIC WRITE
           my_eti%instr_status=-(iabs(stcu_execute_eti))
          endif
          my_eti=>NULL()
         else
          stcu_execute_eti=999
         endif
         return
         end function stcu_execute_eti

#ifndef NO_GPU
         integer(C_INT) function nvcu_execute_eti(etiq_nvcu_loc) !NVCU ETI execution workflow: MT only
!NVCU is non-blocking and multiple ETI can be issued one right after another;
!If the ETI issue was unsuccessful, its status in ETIQ is either kept <scheduled> or changed to <error>;
!The only case when the status is kept <scheduled> is when the HAB packing fails (unable to pack into HAB).
         implicit none
         integer(C_INT), intent(in):: etiq_nvcu_loc !entry number in <etiq_nvcu>
         type(tens_instr_t), pointer:: my_eti
         type(C_PTR):: dtens_,ltens_,rtens_,pptr
         integer(C_INT):: d_hab_entry,l_hab_entry,r_hab_entry
         type(tensor_block_t), pointer:: stens_
         integer(C_SIZE_T):: pack_size
         integer(C_INT) j0,j1,j2,jl,ju,ier
         nvcu_execute_eti=0
         if(etiq_nvcu_loc.ge.0.and.etiq_nvcu_loc.lt.etiq_nvcu%depth) then
          j0=etiq_nvcu%etiq_entry(etiq_nvcu_loc)
          if(j0.gt.0.and.j0.le.etiq%depth) then
           my_eti=>etiq%eti(j0)
           if(my_eti%instr_status.eq.instr_scheduled) then
!$OMP ATOMIC WRITE
            my_eti%instr_status=instr_issued
!$OMP ATOMIC WRITE
            my_eti%time_issued=thread_wtime()
            if(associated(my_eti%tens_op(0)%op_aar_entry)) then
             dtens_=my_eti%tens_op(0)%op_aar_entry%tens_blck_c
             d_hab_entry=my_eti%tens_op(0)%op_aar_entry%buf_entry_host
             if(d_hab_entry.lt.0) then !tensor block is not in HAB yet
              if(associated(my_eti%tens_op(0)%op_aar_entry%tens_blck_f)) then
               call tens_blck_pack(my_eti%tens_op(0)%op_aar_entry%tens_blck_f,my_eti%data_kind,pack_size,pptr,d_hab_entry,ier)
               if(ier.eq.0) then
                call tens_blck_assoc(pptr,ier,ctens=dtens_,gpu_num=etiq_nvcu%te_conf(etiq_nvcu_loc)%cu_id%unit_number)
                if(ier.eq.0) then
                 my_eti%tens_op(0)%op_aar_entry%tens_blck_c=dtens_
                 my_eti%tens_op(0)%op_aar_entry%buf_entry_host=d_hab_entry
                else
                 nvcu_execute_eti=-1
                endif
               else
                nvcu_execute_eti=1 !unable to pack this tensor block: not necessarily an error, maybe no free HAB entries
               endif
              else
               nvcu_execute_eti=-2
              endif
             endif
            else
             dtens_=C_NULL_PTR; d_hab_entry=-1
            endif
            if(nvcu_execute_eti.eq.0.and.associated(my_eti%tens_op(1)%op_aar_entry)) then
             ltens_=my_eti%tens_op(1)%op_aar_entry%tens_blck_c
             l_hab_entry=my_eti%tens_op(1)%op_aar_entry%buf_entry_host
             if(l_hab_entry.lt.0) then !tensor block is not in HAB yet
              if(associated(my_eti%tens_op(1)%op_aar_entry%tens_blck_f)) then
               call tens_blck_pack(my_eti%tens_op(1)%op_aar_entry%tens_blck_f,my_eti%data_kind,pack_size,pptr,l_hab_entry,ier)
               if(ier.eq.0) then
                call tens_blck_assoc(pptr,ier,ctens=ltens_,gpu_num=etiq_nvcu%te_conf(etiq_nvcu_loc)%cu_id%unit_number)
                if(ier.eq.0) then
                 my_eti%tens_op(1)%op_aar_entry%tens_blck_c=ltens_
                 my_eti%tens_op(1)%op_aar_entry%buf_entry_host=l_hab_entry
                else
                 nvcu_execute_eti=-3
                endif
               else
                nvcu_execute_eti=2 !unable to pack this tensor block: not necessarily an error, maybe no free HAB entries
               endif
              else
               nvcu_execute_eti=-4
              endif
             endif
            else
             ltens_=C_NULL_PTR; l_hab_entry=-1
            endif
            if(nvcu_execute_eti.eq.0.and.associated(my_eti%tens_op(2)%op_aar_entry)) then
             rtens_=my_eti%tens_op(2)%op_aar_entry%tens_blck_c
             r_hab_entry=my_eti%tens_op(2)%op_aar_entry%buf_entry_host
             if(r_hab_entry.lt.0) then !tensor block is not in HAB yet
              if(associated(my_eti%tens_op(2)%op_aar_entry%tens_blck_f)) then
               call tens_blck_pack(my_eti%tens_op(2)%op_aar_entry%tens_blck_f,my_eti%data_kind,pack_size,pptr,r_hab_entry,ier)
               if(ier.eq.0) then
                call tens_blck_assoc(pptr,ier,ctens=rtens_,gpu_num=etiq_nvcu%te_conf(etiq_nvcu_loc)%cu_id%unit_number)
                if(ier.eq.0) then
                 my_eti%tens_op(2)%op_aar_entry%tens_blck_c=rtens_
                 my_eti%tens_op(2)%op_aar_entry%buf_entry_host=r_hab_entry
                else
                 nvcu_execute_eti=-5
                endif
               else
                nvcu_execute_eti=3 !unable to pack this tensor block: not necessarily an error, maybe no free HAB entries
               endif
              else
               nvcu_execute_eti=-6
              endif
             endif
            else
             rtens_=C_NULL_PTR; r_hab_entry=-1
            endif
            if(nvcu_execute_eti.eq.0) then
             if(associated(my_eti%tens_op(3)%op_aar_entry)) then
              if(associated(my_eti%tens_op(3)%op_aar_entry%tens_blck_f)) then
               stens_=>my_eti%tens_op(3)%op_aar_entry%tens_blck_f
              else
               stens_=>NULL()
              endif
             else
              stens_=>NULL()
             endif
             select case(my_eti%instr_code)
             case(instr_tensor_contract)
 !Tensor contraction:
              if(d_hab_entry.ge.0.and.l_hab_entry.ge.0.and.r_hab_entry.ge.0) then
               if(allocated(my_eti%instr_aux)) then
                jl=lbound(my_eti%instr_aux,1); ju=ubound(my_eti%instr_aux,1)
                if(ju.ge.jl+1.and.ju-jl+1.eq.my_eti%instr_aux(jl)) then
                 jl=jl+1; j0=my_eti%instr_aux(jl) !Field#1: contraction pattern (mandatory)
                 if(j0.ge.0.and.jl+j0.le.ju) then
                  ier=cuda_task_create(nvcu_tasks(etiq_nvcu_loc)%cuda_task_handle)
                  if(ier.eq.0) then
                   j2=jl+j0+1
                   if(j2.le.ju) then !Field#2 present: reuse flags
                    if(my_eti%instr_aux(j2).eq.1.and.j2.lt.ju) then
                     j1=my_eti%instr_aux(j2+1)
                    else
                     nvcu_execute_eti=-7
                    endif
                   else
                    j1=COPY_BACK
                   endif
                   if(nvcu_execute_eti.eq.0) then
                    if(j0.gt.0) then !at least one tensor is not a scalar
                     ier=gpu_tensor_block_contract_dlf_(my_eti%instr_aux(jl+1:), &
                                                        my_eti%tens_op(1)%op_aar_entry%tens_blck_c, &
                                                        my_eti%tens_op(2)%op_aar_entry%tens_blck_c, &
                                                        my_eti%tens_op(0)%op_aar_entry%tens_blck_c, &
                                                        j1,nvcu_tasks(etiq_nvcu_loc)%cuda_task_handle)
                    else !all three tensor are scalars
                     ier=gpu_tensor_block_contract_dlf_(my_eti%instr_aux(jl:), & !%instr_aux will be ignored
                                                        my_eti%tens_op(1)%op_aar_entry%tens_blck_c, &
                                                        my_eti%tens_op(2)%op_aar_entry%tens_blck_c, &
                                                        my_eti%tens_op(0)%op_aar_entry%tens_blck_c, &
                                                        j1,nvcu_tasks(etiq_nvcu_loc)%cuda_task_handle)
                    endif
                    if(ier.ne.0) then
                     ier=cuda_task_destroy(nvcu_tasks(etiq_nvcu_loc)%cuda_task_handle)
                     nvcu_execute_eti=-7
                    endif
                   endif
                  else
                   nvcu_execute_eti=-8
                  endif
                 else
                  nvcu_execute_eti=-9
                 endif
                else
                 nvcu_execute_eti=-10
                endif
               else
                nvcu_execute_eti=-11
               endif
              else
               nvcu_execute_eti=-12
              endif
             case default
              nvcu_execute_eti=-13
             end select
            endif
           else
            nvcu_execute_eti=-997
           endif
           if(nvcu_execute_eti.lt.0) then !error
!$OMP ATOMIC WRITE
            my_eti%instr_status=nvcu_execute_eti !error code (negative)
           elseif(nvcu_execute_eti.gt.0) then !packing into HAB was impossible at this time
!$OMP ATOMIC WRITE
            my_eti%instr_status=instr_scheduled  !insufficient resources (not an error): revert to status <scheduled>
           else !ETI has been issued successfully
!$OMP ATOMIC WRITE
            my_eti%instr_handle=etiq_nvcu_loc    !CUDA task handle can be retrieved from nvcu_tasks(my_eti%instr_handle)
           endif
           my_eti=>NULL()
          else
           nvcu_execute_eti=-998
          endif
         else
          nvcu_execute_eti=-999
         endif
         return
         end function nvcu_execute_eti

         integer(C_INT) function nvcu_task_status(nvcu_task_num) !queries the status of a CUDA task: MT only
         implicit none
         integer(C_INT), intent(in):: nvcu_task_num !entry number in <nvcu_tasks>||<etiq_nvcu>
         real(C_FLOAT):: in_copy,out_copy,comp_time
         real(4):: ttm,tti
         type(tens_instr_t), pointer:: my_eti
         integer j0,j1,j2
         if(nvcu_task_num.ge.0.and.nvcu_task_num.lt.etiq_nvcu%depth) then
          j0=etiq_nvcu%etiq_entry(nvcu_task_num) !ETIQ entry number
          if(j0.gt.0.and.j0.le.etiq%depth) then
           my_eti=>etiq%eti(j0)
!$OMP ATOMIC READ
           nvcu_task_status=my_eti%instr_status
           if(nvcu_task_status.eq.instr_issued) then
            j1=cuda_task_complete(nvcu_tasks(nvcu_task_num)%cuda_task_handle)
            if(j1.eq.cuda_task_completed) then
             do j1=0,max_tensor_operands-1 !mark all tensor arguments as present on the GPU because the ETI has completed
              if(associated(my_eti%tens_op(j1)%op_aar_entry)) then
               if(c_associated(my_eti%tens_op(j1)%op_aar_entry%tens_blck_c)) then
                j2=tensBlck_set_presence(my_eti%tens_op(j1)%op_aar_entry%tens_blck_c) !mark tensBlck_t as present on the GPU
                if(j2.ne.0) write(jo_cp,'("ERROR(c_process::c_proc_life:nvcu_task_status): tensBlck_set_presence failed: ",i10)') j2
               endif
              endif
             enddo
             j1=eti_mark_aar_used(j0); nvcu_task_status=instr_completed
            elseif(j1.eq.cuda_task_error.or.j1.eq.cuda_task_empty) then
             nvcu_task_status=-1 !error
            endif
!$OMP ATOMIC WRITE
            my_eti%instr_status=nvcu_task_status
            ttm=cuda_task_time(nvcu_tasks(nvcu_task_num)%cuda_task_handle,in_copy,out_copy,comp_time)
!$OMP ATOMIC READ
            tti=my_eti%time_issued
!$OMP ATOMIC WRITE
            my_eti%time_completed=tti+ttm
           endif
           my_eti=>NULL()
          else
           if(verbose) write(jo_cp,'("#ERROR(c_process::c_proc_life:nvcu_task_status): invalid ETIQ entry number: ",i12)') j0
           nvcu_task_status=-998
          endif
         else
          if(verbose) write(jo_cp,'("#ERROR(c_process::c_proc_life:nvcu_task_status): invalid NVCU task number: ",i12)') nvcu_task_num
          nvcu_task_status=-999
         endif
         return
         end function nvcu_task_status

         integer(C_INT) function nvcu_task_cleanup(nvcu_task_num,free_flags) !cleans up after a CUDA task: MT only
!This function destroyes the CUDA task handle, frees GPU buffers, frees HAB buffers, and destroys all tensBlck_t.
!The optional argument <free_flags> can be used to prevent freeing GPU/HAB buffers for some or all ETI operands,
!such that the corresponding tensor blocks will still occupy those buffers. Only the CUDA task handle is destroyed mandatory.
!INPUT:
! # nvcu_task_num - NVCU task number (entry number in <etiq_nvcu>||<nvcu_tasks>);
! # free_flags - (optional) flags: 3 bits per tensor operand starting from 0 (see subroutine <set_cleanup_flags>).
         implicit none
         integer, intent(in):: nvcu_task_num !NVCU task number (entry # in <nvcu_tasks>||<etiq_nvcu>)
         integer, intent(in), optional:: free_flags !controls for partial cleanup (bit flags)
         integer:: j0,je
         nvcu_task_cleanup=0
         if(nvcu_task_num.ge.0.and.nvcu_task_num.lt.etiq_nvcu%depth) then
          j0=etiq_nvcu%etiq_entry(nvcu_task_num)
          if(j0.gt.0.and.j0.le.etiq%depth) then
           je=cuda_task_destroy(nvcu_tasks(nvcu_task_num)%cuda_task_handle); if(je.ne.0) nvcu_task_cleanup=nvcu_task_cleanup+1
           if(present(free_flags)) then
            je=eti_task_cleanup(j0,free_flags)
           else
            je=eti_task_cleanup(j0)
           endif
           if(je.ne.0) then
            nvcu_task_cleanup=nvcu_task_cleanup+2
            if(verbose) write(jo_cp,'("ERROR(c_process::c_proc_life:nvcu_task_cleanup): eti_task_cleanup error: ",i6)') je
           endif
          else
           nvcu_task_cleanup=-998
          endif
         else
          nvcu_task_cleanup=-999
         endif
         return
         end function nvcu_task_cleanup
#endif

#ifndef NO_PHI
         integer(C_INT) function xpcu_execute_eti(eti_loc) !XPCU ETI execution workflow: MT only
         implicit none
         integer(C_INT), intent(in):: eti_loc
         type(tens_instr_t), pointer:: my_eti
         xpcu_execute_eti=0
         if(eti_loc.gt.0.and.eti_loc.le.etiq%depth) then
          my_eti=>etiq%eti(eti_loc)
          if(my_eti%instr_status.eq.instr_scheduled) then
!$OMP ATOMIC WRITE
           my_eti%instr_status=instr_issued
!$OMP ATOMIC WRITE
           my_eti%time_issued=thread_wtime()
           
          else
           xpcu_execute_eti=-998
          endif
          if(xpcu_execute_eti.lt.0) then
!$OMP ATOMIC WRITE
           my_eti%instr_status=xpcu_execute_eti !error
          elseif(xpcu_execute_eti.gt.0) then
           my_eti%instr_status=instr_scheduled  !insufficient resources (not an error)
          endif
          my_eti=>NULL()
         else
          xpcu_execute_eti=-999
         endif
         return
         end function xpcu_execute_eti
#endif

         integer function set_cleanup_flags(free_flags,tens_op_num,flag_to_add) !sets flags for <eti_task_cleanup>: MT only
         implicit none
         integer, intent(inout):: free_flags           !cumulative flags (for all tensor operands): 3 bits per tensor operand
         integer, intent(in):: tens_op_num,flag_to_add !tensor operand number and a specific flag to set
         integer:: j0,j1
         set_cleanup_flags=0
         if(tens_op_num.ge.0.and.tens_op_num.lt.max_tensor_operands) then
          if(flag_to_add.ge.0.and.flag_to_add.lt.8) then
           j1=tens_op_num*3; j0=0
           if(btest(free_flags,j1+2)) j0=j0+4
           if(btest(free_flags,j1+1)) j0=j0+2
           if(btest(free_flags,j1)) j0=j0+1
           free_flags=free_flags+(flag_to_add-j0)*(2**j1)
          else
           set_cleanup_flags=1
          endif
         else
          set_cleanup_flags=2
         endif
         return
         end function set_cleanup_flags

         integer function eti_task_cleanup(eti_num,free_flags) !cleans up after an ETI: MT only
!This function completely or partially destroys tensor arguments for a specific (dead) ETI.
!One can choose to keep some or all tensor operands in HAB/GAB/etc. For each tensor operand,
!one can keep the data in HAB (if any), or the data in GAB (if any), or both, or none.
!If for some tensor operand a flag is set to KEEP_NO_DATA, its AAR entry will be destroyed.
         implicit none
         integer, intent(in):: eti_num              !location of ETI in ETIQ (1..)
         integer, intent(in), optional:: free_flags !cumulative flags: 3 bits per tensor operand
         type(tens_instr_t), pointer:: my_eti=>NULL()
         type(tens_arg_t), pointer:: curr_arg=>NULL()
         integer jt,je,jf
         logical k_tbb,k_hab,k_gpu
!        integer(C_INT):: j0,j1,j2,j3,j4 !debug
         eti_task_cleanup=0
         if(eti_num.gt.0.and.eti_num.le.etiq%depth) then
          my_eti=>etiq%eti(eti_num)
          do jt=0,max_tensor_operands-1
           if(associated(my_eti%tens_op(jt)%op_aar_entry)) then
            curr_arg=>my_eti%tens_op(jt)%op_aar_entry
            if(present(free_flags)) then; jf=free_flags; else; jf=0; endif
            if(jf.ge.KEEP_GPU_DATA) then; jf=jf-KEEP_GPU_DATA; k_gpu=.true.; else; k_gpu=.false.; endif
            if(jf.ge.KEEP_HAB_DATA) then; jf=jf-KEEP_HAB_DATA; k_hab=.true.; else; k_hab=.false.; endif
            if(jf.ge.KEEP_TBB_DATA) then; jf=jf-KEEP_TBB_DATA; k_tbb=.true.; else; k_tbb=.false.; endif
            if(.not.k_hab) then !free the HAB entry of the tensor argument (if any)
             if(curr_arg%buf_entry_host.ge.0) then
              je=free_buf_entry_host(curr_arg%buf_entry_host); if(je.ne.0) eti_task_cleanup=eti_task_cleanup+1 !free HAB entry
              if(c_associated(curr_arg%tens_blck_c)) then
               if(associated(curr_arg%tens_blck_f)) then
                if(c_associated(curr_arg%tens_blck_c,c_loc(curr_arg%tens_blck_f))) nullify(curr_arg%tens_blck_f) !free F only if F->HAB
               endif
!               je=tensBlck_acc_id(curr_arg%tens_blck_c,j0,j1,j2,j3,j4) !debug
!               if(verbose) write(jo_cp,'("DEBUG(c_process::c_proc_life:eti_task_cleanup):",6(1x,i4))') je,j0,j1,j2,j3,j4 !debug
               je=tensBlck_hab_null(curr_arg%tens_blck_c); if(je.ne.0) eti_task_cleanup=eti_task_cleanup+10 !nullify HAB pointer
!               if(je.ne.0.and.verbose) write(jo_cp,'("ERROR(c_process::c_proc_life:eti_task_cleanup): tensBlck_hab_null error: ",i6)') je !debug
              else
               nullify(curr_arg%tens_blck_f)
              endif
             endif
            endif
            if(.not.k_gpu) then !destroy GPU data for the tensor argument (if any)
             if(c_associated(curr_arg%tens_blck_c)) then
              call tens_blck_dissoc(curr_arg%tens_blck_c,je); if(je.ne.0) eti_task_cleanup=eti_task_cleanup+100 !destroy tensBlck_t
             endif
            endif
            if(.not.(k_tbb.or.k_hab.or.k_gpu)) then
             nullify(curr_arg%tens_blck_f)
             je=aar_delete(my_eti%tens_op(jt)%tens_blck_id); if(je.ne.0) eti_task_cleanup=eti_task_cleanup+1000 !delete AA entry
             my_eti%tens_op(jt)%op_aar_entry=>NULL()
            endif
           endif
          enddo
          curr_arg=>NULL()
         else
          eti_task_cleanup=-999
         endif
         return
         end function eti_task_cleanup

         integer function aar_register(tkey,aar_p) !MT only
!Registers an argument in AAR. If the tensor block argument
!is already registered, the %times_needed field is increamented.
!INPUT:
! # tkey: tensor block identifier (key);
!OUTPUT:
! # aar_p: pointer to AAR entry.
         implicit none
         type(tens_blck_id_t), intent(in):: tkey
         type(tens_arg_t), pointer, intent(out):: aar_p
         type(tens_arg_t):: jtarg
         class(*), pointer:: jptr
         jtarg%mpi_tag=-1; jtarg%mpi_process=-1; jtarg%times_needed=0; jtarg%times_used=0
         jtarg%tens_blck_f=>NULL(); jtarg%tens_blck_c=C_NULL_PTR; jtarg%buf_entry_host=-1
         aar_register=aar%search(dict_add_if_not_found,tens_key_cmp,tkey,jtarg,value_out=jptr)
         if(aar_register.eq.dict_key_found.or.aar_register.eq.dict_key_not_found) then
          select type (jptr)
          type is (tens_arg_t)
           aar_p=>jptr; aar_p%times_needed=aar_p%times_needed+1
           aar_register=0
          class default
           aar_register=1
          end select
         endif
         return
         end function aar_register

         integer function aar_delete(tkey) !free AAR entry: MT only
         implicit none
         type(tens_blck_id_t), intent(in):: tkey
         integer:: ier
         aar_delete=aar%search(dict_delete_if_found,tens_key_cmp,tkey,destruct_key_func=cp_destructor)
         return
         end function aar_delete

         integer function eti_mark_aar_used(eti_num) !mark all tensor arguments of an ETI as used: MT only
         implicit none
         integer, intent(in):: eti_num !location of ETI in ETIQ
         type(tens_instr_t), pointer:: my_eti=>NULL()
         integer jt
         eti_mark_aar_used=0
         if(eti_num.gt.0.and.eti_num.le.etiq%depth) then
          my_eti=>etiq%eti(eti_num)
          do jt=0,max_tensor_operands-1
           if(associated(my_eti%tens_op(jt)%op_aar_entry)) then
!$OMP ATOMIC UPDATE
            my_eti%tens_op(jt)%op_aar_entry%times_used=my_eti%tens_op(jt)%op_aar_entry%times_used+1
           endif
          enddo
          my_eti=>NULL()
         else
          eti_mark_aar_used=-999
         endif
         return
         end function eti_mark_aar_used

         subroutine c_proc_quit(errc) !quit c_process
         implicit none
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
         j0=cp_destructor(etiq)
         if(j0.ne.0) then;  write(jo_cp,'("#ERROR(c_process::c_proc_life:c_proc_quit): ETIQ destruction failed!")'); ierr=ierr+20000; endif
 !ETIQ_STCU:
         j0=cp_destructor(etiq_stcu)
         if(j0.ne.0) then;  write(jo_cp,'("#ERROR(c_process::c_proc_life:c_proc_quit): ETIQ_STCU destruction failed!")'); ierr=ierr+50000; endif
         stcu_num_units=-1
 !ETIQ_NVCU:
         j0=cp_destructor(etiq_nvcu); j1=0; if(allocated(nvcu_tasks)) deallocate(nvcu_tasks,STAT=j1)
         if(j0.ne.0.or.j1.ne.0) then;  write(jo_cp,'("#ERROR(c_process::c_proc_life:c_proc_quit): ETIQ_NVCU destruction failed!")'); ierr=ierr+100000; endif
 !ETIQ_XPCU:
         j0=cp_destructor(etiq_xpcu); j1=0; if(allocated(xpcu_tasks)) deallocate(xpcu_tasks,STAT=j1)
         if(j0.ne.0.or.j1.ne.0) then;  write(jo_cp,'("#ERROR(c_process::c_proc_life:c_proc_quit): ETIQ_XPCU destruction failed!")'); ierr=ierr+300000; endif
!AAR:
         j0=aar%destroy(destruct_key_func=cp_destructor)
         if(j0.ne.0) then; write(jo_cp,'("#ERROR(c_process::c_proc_life:c_proc_quit): AAR destruction failed!")'); ierr=ierr+500000; endif
!TBB:
         j0=tbb%destroy(destruct_key_func=cp_destructor,destruct_val_func=cp_destructor)
         if(j0.ne.0) then; write(jo_cp,'("#ERROR(c_process::c_proc_life:c_proc_quit): TBB destruction failed!")'); ierr=ierr+1000000; endif
         return
         end subroutine c_proc_quit

        end subroutine c_proc_life
!----------------------------------------------------------
        recursive function cp_destructor(item) result(ierr) !universal destructor: destroys all allocated components
        implicit none
        class(*):: item !<item> itself is not destroyed, only its allocated components!
        integer:: ierr
        integer:: i,j
        ierr=0
        select type (item)
        type is (tens_blck_id_t)
!         if(verbose) write(jo_cp,'("#DEBUG(c_process::cp_destructor): tens_blck_id_t")') !debug
         if(allocated(item%tens_mlndx)) then; deallocate(item%tens_mlndx,STAT=i); if(i.ne.0) ierr=1; endif
         item%tens_name=' '
        type is (tbb_entry_t)
!         if(verbose) write(jo_cp,'("#DEBUG(c_process::cp_destructor): tbb_entry_t")') !debug
         call tensor_block_destroy(item%tens_blck,i); if(i.ne.0) ierr=2
         item%file_handle=0; item%file_offset=0_8; item%stored_size=0_8
        type is (tens_arg_t)
!         if(verbose) write(jo_cp,'("#DEBUG(c_process::cp_destructor): tens_arg_t")') !debug
         nullify(item%tens_blck_f)
         item%tens_blck_c=C_NULL_PTR; item%buf_entry_host=-1
         item%mpi_tag=-1; item%mpi_process=-1; item%times_needed=0; item%times_used=0
        type is (tens_operand_t)
!         if(verbose) write(jo_cp,'("#DEBUG(c_process::cp_destructor): tens_operand_t")') !debug
         i=cp_destructor(item%tens_blck_id); if(i.ne.0) ierr=3
         item%op_host=-1; item%op_pack_size=0_8; item%op_tag=0; item%op_price=0
         nullify(item%op_aar_entry)
        type is (tens_instr_t)
!         if(verbose) write(jo_cp,'("#DEBUG(c_process::cp_destructor): tens_instr_t")') !debug
         if(allocated(item%instr_aux)) then; deallocate(item%instr_aux,STAT=i); if(i.ne.0) ierr=4; endif
         do j=0,max_tensor_operands-1; i=cp_destructor(item%tens_op(j)); if(i.ne.0) ierr=4; enddo
         item%instr_code=instr_null; item%data_kind='  '; item%instr_priority=0
         item%instr_cost=0.0; item%instr_size=0.0; item%instr_status=instr_null
         item%instr_cu=cu_t(-1,-1); item%instr_handle=-1; item%args_ready=0
         item%time_touched=0.; item%time_data_ready=0.; item%time_issued=0.; item%time_completed=0.; item%time_uploaded=0.
        type is (etiq_cu_t)
!         if(verbose) write(jo_cp,'("#DEBUG(c_process::cp_destructor): etiq_cu_t")') !debug
         if(allocated(item%etiq_entry)) then; deallocate(item%etiq_entry,STAT=i); if(i.ne.0) ierr=5; endif
         if(allocated(item%te_conf)) then; deallocate(item%te_conf,STAT=i); if(i.ne.0) ierr=5; endif
         item%depth=0; item%scheduled=0; item%bp=0; item%ip=0
        type is (etiq_t)
!         if(verbose) write(jo_cp,'("#DEBUG(c_process::cp_destructor): etiq_t")') !debug
         if(allocated(item%ffe_stack)) then; deallocate(item%ffe_stack,STAT=i); if(i.ne.0) ierr=6; endif
         if(allocated(item%next)) then; deallocate(item%next,STAT=i); if(i.ne.0) ierr=6; endif
         if(allocated(item%eti)) then
          do j=lbound(item%eti,1),ubound(item%eti,1)
           i=cp_destructor(item%eti(j)); if(i.ne.0) ierr=6
          enddo
          deallocate(item%eti,STAT=i); if(i.ne.0) ierr=6
         endif
         item%depth=0; item%scheduled=0; item%ffe_sp=0
         item%last(:)=0; item%ip(:)=0; item%ic(:)=0
        class default
         write(jo_cp,'("#ERROR(c_process::cp_destructor): unknown type/class!")')
         ierr=-1
        end select
        return
        end function cp_destructor
!----------------------------------------------------------------
        integer function tens_blck_id_create(this,t_name,t_mlndx)
        implicit none
        class(tens_blck_id_t):: this
        character(*), intent(in), optional:: t_name
        integer, intent(in), optional:: t_mlndx(0:*) !0th element is the length of the multiindex [1:length]
        integer:: i,l
        tens_blck_id_create=0
        if(present(t_name)) then
         l=len_trim(t_name)
         if(l.gt.0) then
          i=1; do while(iachar(t_name(i:i)).le.32); i=i+1; if(i.gt.l) exit; enddo
          if(i.le.l) then; l=l-i+1; this%tens_name=t_name(i:i+min(tensor_name_len,l)-1); endif
         else
          this%tens_name=' '
         endif
        else
         this%tens_name=' '
        endif
        if(present(t_mlndx)) then
         l=t_mlndx(0) !length of the multiindex
         if(l.ge.0) then
          if(allocated(this%tens_mlndx)) then
           if(this%tens_mlndx(0).ne.l) then
            deallocate(this%tens_mlndx); allocate(this%tens_mlndx(0:l),STAT=i)
           else
            i=0
           endif
          else
           allocate(this%tens_mlndx(0:l),STAT=i)
          endif
          if(i.eq.0) then
           this%tens_mlndx(0)=l; do i=1,l; this%tens_mlndx(i)=t_mlndx(i); enddo
          else
           tens_blck_id_create=1
          endif
         else
          if(allocated(this%tens_mlndx)) deallocate(this%tens_mlndx)
         endif
        else
         if(allocated(this%tens_mlndx)) deallocate(this%tens_mlndx)
        endif
        return
        end function tens_blck_id_create
!--------------------------------------------------
        integer function tens_blck_id_destroy(this)
        implicit none
        class(tens_blck_id_t):: this
        integer:: i
        tens_blck_id_destroy=0
        if(allocated(this%tens_mlndx)) then
         deallocate(this%tens_mlndx,STAT=i); if(i.ne.0) tens_blck_id_destroy=1
        endif
        return
        end function tens_blck_id_destroy
!-----------------------------------------------
        integer function tens_key_cmp(key1,key2)
        implicit none
        class(tens_blck_id_t):: key1,key2
        integer:: i,l1,l2
        tens_key_cmp=dict_key_eq
        do i=1,tensor_name_len
         l1=iachar(key1%tens_name(i:i)); l2=iachar(key2%tens_name(i:i))
         if(l1.le.32.and.l2.le.32) exit !special symbols (including spaces) terminate the string
         if(l1.lt.l2) then
          tens_key_cmp=dict_key_lt; exit
         elseif(l1.gt.l2) then
          tens_key_cmp=dict_key_gt; exit
         endif
        enddo
        if(tens_key_cmp.eq.dict_key_eq) then
         if(allocated(key1%tens_mlndx)) then
          if(allocated(key2%tens_mlndx)) then
           l1=key1%tens_mlndx(0); l2=key2%tens_mlndx(0)
           if(l1.lt.l2) then
            tens_key_cmp=dict_key_lt
           elseif(l1.gt.l2) then
            tens_key_cmp=dict_key_gt
           else
            do i=1,l1
             if(key1%tens_mlndx(i).lt.key2%tens_mlndx(i)) then
              tens_key_cmp=dict_key_lt; exit
             elseif(key1%tens_mlndx(i).gt.key2%tens_mlndx(i)) then
              tens_key_cmp=dict_key_gt; exit
             endif
            enddo
           endif
          else
           tens_key_cmp=dict_key_gt
          endif
         else
          if(allocated(key2%tens_mlndx)) tens_key_cmp=dict_key_lt
         endif
        endif
        return
        end function tens_key_cmp
!----------------------------------------------------
        integer function etiq_init(this,queue_length) !SERIAL: MT only
!This function initializes an <etiq_t> object.
        implicit none
        class(etiq_t):: this
        integer(C_INT), intent(in):: queue_length
        integer(C_INT):: i
        etiq_init=0
        if(queue_length.gt.0) then
         allocate(this%eti(1:queue_length),this%next(1:queue_length),this%ffe_stack(1:queue_length),STAT=i)
         if(i.eq.0) then
          this%depth=queue_length; this%scheduled=0; this%last(:)=0; this%ip(:)=0; this%ic(:)=0
          do i=1,queue_length; this%eti(i)%instr_status=instr_null; enddo
          do i=1,queue_length; this%ffe_stack(i)=i; enddo; this%ffe_sp=1
         else
          etiq_init=1
         endif
        else
         etiq_init=2
        endif
        return
        end function etiq_init
!-------------------------------------------------------
        integer function etiq_cu_init(this,queue_length) !SERIAL: MT only
!This function initializes an <etiq_cu_t> object.
        implicit none
        class(etiq_cu_t):: this
        integer(C_INT), intent(in):: queue_length
        integer(C_INT):: i
        etiq_cu_init=0
        if(queue_length.gt.0) then
         allocate(this%etiq_entry(0:queue_length-1),this%te_conf(0:queue_length-1),STAT=i)
         if(i.eq.0) then
          do i=0,queue_length-1; this%etiq_entry(i)=0; enddo
          do i=0,queue_length-1; this%te_conf(i)%cu_id=cu_t(-1,-1); enddo
          this%depth=queue_length; this%scheduled=0; this%bp=0; this%ip=0
         else
          etiq_cu_init=1
         endif
        else
         etiq_cu_init=2
        endif
        return
        end function etiq_cu_init
!---------------------------------------------------------
        integer function tens_operand_pack(this,tens_data)
!This function packs the public part of <this> (tens_operand_t) into a plain byte packet <tens_data>.
!Format of the output packet:
! INTEGER(C_SIZE_T): packet size (in bytes);
! INTEGER(C_INT): length of the tensor name;
! INTEGER(1), DIMENSION(:): tensor name;
! INTEGER(C_INT): length of the tensor block key multiindex;
! INTEGER(C_INT), DIMENSION(:): tensor block key multiindex;
! INTEGER(C_INT): host MPI process of the tensor block (operand);
! INTEGER(C_SIZE_T): packed size of the tensor block (see <tens_blck_packet_size>);
! INTEGER(C_INT): MPI delivery tag;
! INTEGER(C_INT): price of the tensor block;
        implicit none
        class(tens_operand_t):: this
        type(C_PTR):: tens_data,c_addr
        integer(C_INT):: i
        integer(C_SIZE_T):: pack_size,i_s,i1_s
        integer(C_SIZE_T), pointer:: cz_p
        integer(C_INT), pointer:: i_p
        integer(C_INT), pointer:: id1_p(:)
        integer(1):: i1
        integer(1), pointer:: i1d1_p(:)

        tens_operand_pack=0
        if(c_associated(tens_data)) then
         i=0; i1=0; i_s=sizeof(i); i1_s=sizeof(i1)
         pack_size=0; pack_size=sizeof(pack_size)
!Tensor name:
         c_addr=ptr_offset(tens_data,pack_size); call c_f_pointer(c_addr,i_p)
         i_p=len_trim(this%tens_blck_id%tens_name); pack_size=pack_size+i_s !length
         if(i_p.gt.0) then
          c_addr=ptr_offset(tens_data,pack_size); call c_f_pointer(c_addr,i1d1_p,shape=[i_p])
          do i=1,i_p; i1d1_p(i)=iachar(this%tens_blck_id%tens_name(i:i)); enddo !tensor name
          pack_size=pack_size+i1_s*i_p
         else
          tens_operand_pack=1; return
         endif
!Tensor block key multiindex:
         c_addr=ptr_offset(tens_data,pack_size); call c_f_pointer(c_addr,i_p)
         if(allocated(this%tens_blck_id%tens_mlndx)) then
          if(lbound(this%tens_blck_id%tens_mlndx,1).eq.0) then
           if(ubound(this%tens_blck_id%tens_mlndx,1).eq.this%tens_blck_id%tens_mlndx(0)) then
            i_p=this%tens_blck_id%tens_mlndx(0); pack_size=pack_size+i_s
            c_addr=ptr_offset(tens_data,pack_size); call c_f_pointer(c_addr,id1_p,shape=[i_p])
            id1_p(1:i_p)=this%tens_blck_id%tens_mlndx(1:i_p); pack_size=pack_size+i_s*i_p
           else
            tens_operand_pack=2; return
           endif
          else
           tens_operand_pack=3; return
          endif
         else
          i_p=0; pack_size=pack_size+i_s
         endif
!Host MPI process:
         c_addr=ptr_offset(tens_data,pack_size); call c_f_pointer(c_addr,i_p)
         i_p=this%op_host; pack_size=pack_size+i_s
!Packed size of the tensor block:
         c_addr=ptr_offset(tens_data,pack_size); call c_f_pointer(c_addr,cz_p)
         cz_p=this%op_pack_size; pack_size=pack_size+sizeof(pack_size)
!MPI delivery tag:
         c_addr=ptr_offset(tens_data,pack_size); call c_f_pointer(c_addr,i_p)
         i_p=this%op_tag; pack_size=pack_size+i_s
!Price of the tensor block:
         c_addr=ptr_offset(tens_data,pack_size); call c_f_pointer(c_addr,i_p)
         i_p=this%op_price; pack_size=pack_size+i_s
!Packet size:
         call c_f_pointer(tens_data,cz_p); cz_p=pack_size
        else
         tens_operand_pack=999
        endif
        return
        end function tens_operand_pack
!-----------------------------------------------------------
        integer function tens_operand_unpack(this,tens_data)
!This function unpacks the public part of <this> (tens_operand_t) from a plain byte packet <tens_data>.
!The format of the packet is given in <tens_operand_pack>.
        implicit none
        class(tens_operand_t):: this
        type(C_PTR):: tens_data,c_addr
        integer(C_INT):: i
        integer(C_SIZE_T):: pack_size,s,i_s,i1_s
        integer(C_SIZE_T), pointer:: cz_p
        integer(C_INT), pointer:: i_p
        integer(C_INT), pointer:: id1_p(:)
        integer(1):: i1
        integer(1), pointer:: i1d1_p(:)

        tens_operand_unpack=0
        if(c_associated(tens_data)) then
         i=0; i1=0; i_s=sizeof(i); i1_s=sizeof(i1)
!Packet size:
         call c_f_pointer(tens_data,cz_p); pack_size=cz_p; s=sizeof(pack_size)
         if(pack_size.le.0) then; tens_operand_unpack=1; return; endif
!Tensor name:
         c_addr=ptr_offset(tens_data,s); call c_f_pointer(c_addr,i_p); s=s+i_s
         if(i_p.le.0.or.i_p.gt.tensor_name_len) then; tens_operand_unpack=2; return; endif
         c_addr=ptr_offset(tens_data,s); call c_f_pointer(c_addr,i1d1_p,shape=[i_p])
         do i=1,i_p; this%tens_blck_id%tens_name(i:i)=achar(i1d1_p(i)); enddo; s=s+i1_s*i_p
!Tensor block key multiindex:
         c_addr=ptr_offset(tens_data,s); call c_f_pointer(c_addr,i_p); s=s+i_s
         if(i_p.lt.0.or.i_p.gt.max_tensor_rank) then; tens_operand_unpack=3; return; endif
         c_addr=ptr_offset(tens_data,s); call c_f_pointer(c_addr,id1_p,shape=[i_p])
         if(.not.allocated(this%tens_blck_id%tens_mlndx)) then
          allocate(this%tens_blck_id%tens_mlndx(0:i_p),STAT=i); if(i.ne.0) then; tens_operand_unpack=4; return; endif
         endif
         if(size(this%tens_blck_id%tens_mlndx).ne.1+i_p) then
          deallocate(this%tens_blck_id%tens_mlndx,STAT=i); if(i.ne.0) then; tens_operand_unpack=5; return; endif
          allocate(this%tens_blck_id%tens_mlndx(0:i_p),STAT=i); if(i.ne.0) then; tens_operand_unpack=6; return; endif
         endif
         this%tens_blck_id%tens_mlndx(0)=i_p; this%tens_blck_id%tens_mlndx(1:i_p)=id1_p(1:i_p); s=s+i_s*i_p
!Host MPI process:
         c_addr=ptr_offset(tens_data,s); call c_f_pointer(c_addr,i_p); this%op_host=i_p; s=s+i_s
!Packed size of the tensor block:
         c_addr=ptr_offset(tens_data,s); call c_f_pointer(c_addr,cz_p); this%op_pack_size=cz_p; s=s+sizeof(pack_size)
!MPI delivery tag:
         c_addr=ptr_offset(tens_data,s); call c_f_pointer(c_addr,i_p); this%op_tag=i_p; s=s+i_s
!Price of the tensor block:
         c_addr=ptr_offset(tens_data,s); call c_f_pointer(c_addr,i_p); this%op_price=i_p; s=s+i_s
         if(s.ne.pack_size) tens_operand_unpack=7
        else
         tens_operand_unpack=999
        endif
        return
        end function tens_operand_unpack
!-----------------------------------------------
        integer function eti_pack(this,eti_data)
!This function packs the public part of an ETI <this> (tens_instr_t) into a plain byte packet <eti_data>.
!Format of the output packet:
! INTEGER(C_SIZE_T): packet size (in bytes);
! INTEGER(C_INT): tensor instruction code;
! INTEGER(C_INT): data kind;
! INTEGER(C_INT): length of the auxiliary instruction info;
! INTEGER(C_INT), DIMENSION(:): auxiliary instruction info;
! INTEGER(C_INT): tensor instruction priority;
! REAL(4): computational cost of the tensor instruction (Flops);
! REAL(4): total size of all tensor operands (Words of data kind);
! INTEGER: total number of tensor operands present in this packet;
! PACKET(tens_operand_t): operand number (C_INT) + packed tensor operand;
! PACKET(tens_operand_t): operand number (C_INT) + packed tensor operand;
! ...
        implicit none
        class(tens_instr_t):: this
        type(C_PTR):: eti_data,c_addr
        integer(C_INT):: i,ier
        integer(C_INT), pointer:: i_p,id1_p(:)
        integer(C_SIZE_T):: pack_size,i_s,r4_s
        integer(C_SIZE_T), pointer:: cz_p
        real(4):: rval4
        real(4), pointer:: r4_p
        integer(C_INT):: op_present(0:max_tensor_operands-1)

        eti_pack=0
        if(c_associated(eti_data)) then
         i=0; rval4=0.0; i_s=sizeof(i); r4_s=sizeof(rval4)
         pack_size=0; pack_size=sizeof(pack_size)
!Tensor instruction code:
         c_addr=ptr_offset(eti_data,pack_size); call c_f_pointer(c_addr,i_p)
         i_p=int(this%instr_code,C_INT); pack_size=pack_size+i_s
!Data kind:
         c_addr=ptr_offset(eti_data,pack_size); call c_f_pointer(c_addr,i_p)
         select case(this%data_kind)
         case('r4','R4')
          i_p=R4
         case('r8','R8')
          i_p=R8
         case('c8','C8')
          i_p=C8
         case default
          eti_pack=1; return !unknown data kind
         end select
         pack_size=pack_size+i_s
!Auxiliary instruction info:
         c_addr=ptr_offset(eti_data,pack_size); call c_f_pointer(c_addr,i_p)
         if(allocated(this%instr_aux)) then
          i_p=size(this%instr_aux); pack_size=pack_size+i_s
          c_addr=ptr_offset(eti_data,pack_size); call c_f_pointer(c_addr,id1_p,shape=[i_p])
          id1_p(1:i_p)=this%instr_aux(lbound(this%instr_aux,1):ubound(this%instr_aux,1))
          pack_size=pack_size+i_s*i_p
         else
          i_p=0; pack_size=pack_size+i_s
         endif
!Tensor instruction priority:
         c_addr=ptr_offset(eti_data,pack_size); call c_f_pointer(c_addr,i_p)
         i_p=this%instr_priority; pack_size=pack_size+i_s
!Computational cost:
         c_addr=ptr_offset(eti_data,pack_size); call c_f_pointer(c_addr,r4_p)
         r4_p=this%instr_cost; pack_size=pack_size+r4_s
!Memory demands:
         c_addr=ptr_offset(eti_data,pack_size); call c_f_pointer(c_addr,r4_p)
         r4_p=this%instr_size; pack_size=pack_size+r4_s
!Number of tensor operands:
         c_addr=ptr_offset(eti_data,pack_size); call c_f_pointer(c_addr,i_p); i_p=0
         do i=0,max_tensor_operands-1
          if(len_trim(this%tens_op(i)%tens_blck_id%tens_name).gt.0) then !tensor operand #i present
           i_p=i_p+1; op_present(i)=1
          else !tensor operand #i absent
           op_present(i)=0
          endif
         enddo
         pack_size=pack_size+i_s
!Packets for all present tensor operands:
         do i=0,max_tensor_operands-1
          if(op_present(i).ne.0) then
 !Operand number:
           c_addr=ptr_offset(eti_data,pack_size); call c_f_pointer(c_addr,i_p)
           i_p=i; pack_size=pack_size+i_s
 !Operand packet:
           c_addr=ptr_offset(eti_data,pack_size); ier=this%tens_op(i)%pack(c_addr)
           if(ier.eq.0) then
            call c_f_pointer(c_addr,cz_p); pack_size=pack_size+cz_p !cz_p now contains the size of the operand packet
           else
            eti_pack=2+i; return
           endif
          endif
         enddo
!Total packet size:
         call c_f_pointer(eti_data,cz_p); cz_p=pack_size
        else
         eti_pack=999
        endif
        return
        end function eti_pack
!-------------------------------------------------
        integer function eti_unpack(this,eti_data)
!This function unpacks a plain byte packet <eti_data> into <this> (tens_instr_t).
!The format of the packet is given in <eti_pack>.
        implicit none
        class(tens_instr_t):: this
        type(C_PTR):: eti_data,c_addr
        integer(C_INT):: i,n,ier
        integer(C_INT), pointer:: i_p,id1_p(:)
        integer(C_SIZE_T):: pack_size,s,i_s,r4_s
        integer(C_SIZE_T), pointer:: cz_p
        real(4):: rval4
        real(4), pointer:: r4_p
        integer(C_INT):: op_present(0:max_tensor_operands-1)

        eti_unpack=0
        if(c_associated(eti_data)) then
         i=0; rval4=0.0; i_s=sizeof(i); r4_s=sizeof(rval4)
!Total size of the packet:
         call c_f_pointer(eti_data,cz_p); pack_size=cz_p; s=sizeof(pack_size)
!Tensor instruction code:
         c_addr=ptr_offset(eti_data,s); call c_f_pointer(c_addr,i_p)
         this%instr_code=int(i_p,2); s=s+i_s
!Data kind:
         c_addr=ptr_offset(eti_data,s); call c_f_pointer(c_addr,i_p)
         select case(i_p)
         case(R4)
          this%data_kind='r4'
         case(R8)
          this%data_kind='r8'
         case(C8)
          this%data_kind='c8'
         case default
          eti_unpack=1; return !unknown data kind
         end select
         s=s+i_s
!Auxiliary instruction info:
         c_addr=ptr_offset(eti_data,s); call c_f_pointer(c_addr,i_p); s=s+i_s
         if(i_p.gt.0) then
          c_addr=ptr_offset(eti_data,s); call c_f_pointer(c_addr,id1_p,shape=[i_p])
          if(.not.allocated(this%instr_aux)) then
           allocate(this%instr_aux(0:i_p-1),STAT=i); if(i.ne.0) then; eti_unpack=2; return; endif
          endif
          if(size(this%instr_aux).ne.i_p) then
           deallocate(this%instr_aux,STAT=i); if(i.ne.0) then; eti_unpack=3; return; endif
           allocate(this%instr_aux(0:i_p-1),STAT=i); if(i.ne.0) then; eti_unpack=4; return; endif
          endif
          this%instr_aux(lbound(this%instr_aux,1):ubound(this%instr_aux,1))=id1_p(1:i_p)
          s=s+i_s*i_p
         elseif(i_p.lt.0) then
          eti_unpack=5; return
         endif
!Tensor instruction priority:
         c_addr=ptr_offset(eti_data,s); call c_f_pointer(c_addr,i_p)
         this%instr_priority=i_p; s=s+i_s
!Instruction cost:
         c_addr=ptr_offset(eti_data,s); call c_f_pointer(c_addr,r4_p)
         this%instr_cost=r4_p; s=s+r4_s
!Memory demands:
         c_addr=ptr_offset(eti_data,s); call c_f_pointer(c_addr,r4_p)
         this%instr_size=r4_p; s=s+r4_s
!Number of tensor operands:
         c_addr=ptr_offset(eti_data,s); call c_f_pointer(c_addr,i_p); n=i_p; s=s+i_s
         if(n.lt.0.or.n.gt.max_tensor_operands) then; eti_unpack=6; return; endif
         if(n.gt.0) then
          op_present(0:max_tensor_operands-1)=0
!Tensor operands:
          do i=1,n
           c_addr=ptr_offset(eti_data,s); call c_f_pointer(c_addr,i_p); s=s+i_s
           if(i_p.ge.0.and.i_p.lt.max_tensor_operands) then
            if(op_present(i_p).eq.0) then
             op_present(i_p)=i
             c_addr=ptr_offset(eti_data,s); call c_f_pointer(c_addr,cz_p); s=s+cz_p
             ier=this%tens_op(i_p)%unpack(c_addr); if(ier.ne.0) then; eti_unpack=10+i; return; endif
            else
             eti_unpack=100+i; return
            endif
           else
            eti_unpack=200+i; return
           endif
          enddo
         endif
!Check the total packet size:
         if(s.ne.pack_size) eti_unpack=7
        else
         eti_unpack=999
        endif
        return
        end function eti_unpack
!---------------------------------------------------------------
        integer(8) function tens_blck_packet_size(tens,dtk,ierr)
!Given an instance of tensor_block_t <tens> and required data kind <dtk>,
!this function returns the size (in bytes) of the corresponding packet.
!`Note that currently only the dense (DLF) tensor data layout is supported.
        implicit none
        type(tensor_block_t), intent(in):: tens
        character(2), intent(in):: dtk
        integer:: ierr
        integer(C_INT), parameter:: i=0
        integer(C_SIZE_T), parameter:: s=0
        real(4), parameter:: spn=0.0
        real(8), parameter:: dpn=0d0
        complex(8), parameter:: cdpn=cmplx(0d0,0d0,8)
        integer(8):: trank,i_size,s_size,cdpn_size
        ierr=0; tens_blck_packet_size=0_8
        if(tens%tensor_shape%num_dim.ge.0) then
         trank=int(tens%tensor_shape%num_dim,8)
         if(trank.eq.0_8.and.tens%tensor_block_size.ne.1_8) then; ierr=1; return; endif !trap
         if(trank.gt.0_8.and.tens%tensor_block_size.le.0_8) then; ierr=2; return; endif !trap
         i_size=sizeof(i); s_size=sizeof(s); cdpn_size=sizeof(cdpn)
         select case(dtk)
         case('r4')
          tens_blck_packet_size=cdpn_size+s_size*2_8+i_size*(2_8+trank*5_8)+sizeof(spn)*tens%tensor_block_size
         case('r8')
          tens_blck_packet_size=cdpn_size+s_size*2_8+i_size*(2_8+trank*5_8)+sizeof(dpn)*tens%tensor_block_size
         case('c8')
          tens_blck_packet_size=cdpn_size+s_size*2_8+i_size*(2_8+trank*5_8)+cdpn_size*tens%tensor_block_size
         case default
          ierr=3
         end select
        else
         ierr=4
        endif
        return
        end function tens_blck_packet_size
!--------------------------------------------------------------------------
        subroutine tens_blck_pack(tens,dtk,packet_size,pptr,entry_num,ierr) !SERIAL
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
!--------------------------------------------------
        subroutine tens_blck_unpack(tens,pptr,ierr) !SERIAL
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
! - This subroutine does NOT free the corresponding Host argument buffer (HAB) entry!
! - Packet structure is specified in <tens_blck_pack>.
! - The tensor block is always destroyed before re-creation and unpacking the data in it.
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
!--------------------------------------------------
        subroutine tens_blck_update(tens,pptr,ierr)
!This subroutine updates the <tensor_block_t> data
!according to the current data in the tensor block packet.
!INPUT:
! - tens - tensor block;
! - pptr - C_PTR to the HAB packet;
!OUTPUT:
! - ierr - error code (0:success).
        implicit none
        type(tensor_block_t), intent(inout):: tens
        type(C_PTR), intent(in):: pptr
        integer:: ierr
        integer(C_INT) i; integer(C_SIZE_T) offs
        integer(C_SIZE_T), pointer:: elems_count
        integer(C_INT), pointer:: trank,dtk
        real(4), pointer, contiguous:: elems_r4(:)
        real(8), pointer, contiguous:: elems_r8(:)
        complex(8), pointer, contiguous:: elems_c8(:)
        complex(8):: scalar_c8=dcmplx(0d0,0d0)
        type(C_PTR):: c_addr

        ierr=0
        offs=sizeof(offs); c_addr=ptr_offset(pptr,offs); call c_f_pointer(c_addr,dtk)
        offs=offs+sizeof(i); c_addr=ptr_offset(pptr,offs); call c_f_pointer(c_addr,elems_count)
        offs=offs+sizeof(offs); c_addr=ptr_offset(pptr,offs); call c_f_pointer(c_addr,trank)
        if(elems_count.gt.0.and.trank.ge.0) then
         offs=offs+5*trank*sizeof(i)+sizeof(scalar_c8); c_addr=ptr_offset(pptr,offs)
         if(elems_count.eq.tens%tensor_block_size) then
          select case(dtk)
          case(R4)
           call c_f_pointer(c_addr,elems_r4,shape=[elems_count])
           do offs=0,elems_count-1; tens%data_real4(offs)=elems_r4(offs+1); enddo
          case(R8)
           call c_f_pointer(c_addr,elems_r8,shape=[elems_count])
           do offs=0,elems_count-1; tens%data_real8(offs)=elems_r8(offs+1); enddo
          case(C8)
           call c_f_pointer(c_addr,elems_c8,shape=[elems_count])
           do offs=0,elems_count-1; tens%data_cmplx8(offs)=elems_c8(offs+1); enddo
          case default
           ierr=3
          end select
         else
          ierr=2
         endif
        else
         ierr=1
        endif
        return
        end subroutine tens_blck_update
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
!----------------------------------------------
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
!-------------------------------------------------------------------------------------------------------
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
