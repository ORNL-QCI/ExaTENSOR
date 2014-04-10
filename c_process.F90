!This module provides functionality for a QFORCE Computing Process (C-PROCESS).
!In essence, this is a (single-node) elementary tensor instruction scheduler (ETIS).
!AUTHOR: Dmitry I. Lyakh (Dmytro I. Liakh): quant4me@gmail.com
!REVISION: 2014/04/10
!NOTES:
! - Data synchronization in an instance of <tensor_block_t> (Fortran)
!   associated with a Host Argument Buffer entry can allocate only regular CPU memory
!   (the one outside the pinned Host Argument buffer). Hence that newly allocated
!   memory cannot be used with GPU.
!Acronyms:
! - ETI - Elementary Tensor Instruction, or simply TI.
! - HAB - Host Argument Buffer;
! - TIS - Tensor Instruction Scheduler;
! - TIQ - Tensor Instruction Queue;
! - LTBB - Local Tensor Block Bank;
! - AAL - Active Argument List;
! - TAL - Tensor Algebra Library;
       module c_process
        use, intrinsic:: ISO_C_BINDING
        use tensor_algebra
        use service
        use extern_names !`dependency to be removed
!PARAMETERS:
 !General:
        integer, parameter:: CZ=C_SIZE_T
 !Host Argument Buffer (HAB):
        integer(C_SIZE_T), parameter, private:: max_arg_buf_size=1024_CZ*1024_CZ*1024_CZ !max size in bytes of the HAB
        integer(C_SIZE_T), parameter, private:: min_arg_buf_size=32_CZ*1024_CZ*1024_CZ   !min size in bytes of the HAB
 !Tensor Instruction Scheduler (TIS):
  !Tensor naming:
        integer(C_INT), parameter:: tensor_name_len=32 !max number of characters used for tensor names
  !Tensor instruction status (any negative value will correspond to a failure, designating the error code):
        integer, parameter:: instr_null=0          !uninitialized instruction
        integer, parameter:: instr_data_wait=1     !instruction is waiting for data to arrive
        integer, parameter:: instr_ready_to_exec=2 !instruction is ready to be executed (data has arrived)
        integer, parameter:: instr_issued=3        !instruction is being executed
        integer, parameter:: instr_completed=4     !instruction has completed
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
 !Tensor Instruction Scheduler (TIS):
  !Locally present tensor argument that can reside either in HAB {tens_blck_c+buf_entry_host} or in LTBB {tens_blck_f}:
        type tens_arg_t
         type(tensor_block_t), pointer:: tens_blck_f=>NULL() !pointer to a tensor_block in LTBB
         type(C_PTR):: tens_blck_c       !C pointer to tensBlck_t in HAB (see "tensor_algebra_gpu_nvidia.h")
         integer(C_INT):: buf_entry_host !HAB entry number where the tensor block resides as a packet (-1: not there)
        end type tens_arg_t
  !Tensor operand:
        type tens_operand_t
         character(LEN=tensor_name_len):: tens_name    !tensor name
         integer(C_INT), pointer:: tens_mlndx(:)       !tensor block multiindex
         integer(C_SIZE_T):: tens_size                 !size of the tensor operand in bytes
         integer(C_INT):: tens_host                    !MPI process # where the tensor operand resides
         integer(C_INT):: tens_price                   !current price of the tensor operand (tensor block)
         type(tens_arg_t), pointer:: arg_entry=>NULL() !pointer to AAL entry (active tensor argument list)
        end type tens_operand_t
  !Dispatched elementary tensor instruction:
        type tens_instr_t
         integer:: instr_status                        !tensor instruction status (see above)
         integer:: instr_code                          !tensor instruction code (see above)
         integer:: instr_priority                      !tensor instruction priority
         real(8):: instr_cost                          !approx. instruction computational cost (FLOPs)
         real(8):: instr_size                          !approx. instruction memory demands (Bytes)
         real(4):: instr_time_beg                      !time the instruction scheduled
         real(4):: instr_time_end                      !time the instruction completed
         integer:: args_ready                          !each bit is set to 1 when the corresponding operand is in AAL
         type(tens_operand_t), pointer:: tens_op0=>NULL() !pointer to the tensor block operand #0
         type(tens_operand_t), pointer:: tens_op1=>NULL() !pointer to the tensor block operand #1
         type(tens_operand_t), pointer:: tens_op2=>NULL() !pointer to the tensor block operand #2
         type(tens_operand_t), pointer:: tens_op3=>NULL() !pointer to the tensor block opearnd #3
         type(C_PTR):: cuda_task                       !CUDA task handle associated with this tensor instruction (if any)
        end type tens_instr_t
!DATA:
 !Host Argument Buffer (HAB):
        integer(C_SIZE_T), private:: arg_buf_size_host=0 !actual size in bytes of the Host argument buffer
        integer(C_INT), private:: max_args_host=0 !max number of arguments (of lowest-size level) that can fit in HAB
 !Tensor Instruction Scheduler (TIS):
  !Tensor Instruction Queue (TIQ):
        integer(C_INT), private:: tens_instr_que_size=0 !actual size of the tensor instruction queue (TIQ), set by LR
        integer(C_INT), private:: tens_instr_que_lim=0  !current limit in the tensor instruction queue (TIQ)
        type(tens_instr_t), allocatable, private:: tens_instr_que(:) !tensor instruction queue
  !Active Tensor Argument List (AAL):
        integer(C_INT), private:: act_args_size=0 !actual size of the active tensor argument list (AAL), set by LR
        integer(C_INT), private:: act_args_lim=0  !current number of active (locally present) tensor arguments
        type(tens_arg_t), allocatable, private:: act_args(:) !list of active (locally present) tensor arguments
!--------------------------------------------------------------------------------------------------------------
       contains
!CODE:
        subroutine c_proc_life(ierr)
!This subroutine implements C-PROCESS life cycle.
!INPUT (external):
! - jo - output device (service.mod);
! - impir - MPI rank of the process (service.mod);
! - impis - global MPI comminicator size (service.mod);
! - max_threads - max number of threads available to the current MPI process (service.mod);
! - {gpu_start:gpu_start+gpu_count-1} - range of Nvidia GPUs assigned to the current process (service.mod);
!OUTPUT:
! - ierr - error code (0:success);
! - output of elementary tensor instructions executed.
        implicit none
        integer, intent(inout):: ierr
!------------------------------------
        integer(C_INT), parameter:: max_arg_buf_levels=256 !max number of argument buffer levels (do not exceed C values)
!------------------------------------------
        integer(C_INT) i,j,k,l,m,n,err_code
        integer(C_SIZE_T) blck_sizes(0:max_arg_buf_levels-1)
        real(8) tm

        ierr=0; write(jo,'("#MSG(c_process:c_proc_life): I am a C-process (Computing MPI Process), rank = ",i7)') impir
!Initialization:
        write(jo,'("#MSG(c_process:c_proc_life): Initialization:")')
 !Init TAL infrastructure (TAL buffers, cuBLAS, etc.):
        write(jo,'("#MSG(c_process:c_proc_life): Allocating TAL argument buffers ... ")',advance='no')
        tm=thread_wtime()
        arg_buf_size_host=max_arg_buf_size !desired (max) HAB size
        i=arg_buf_allocate(arg_buf_size_host,max_args_host,gpu_start,gpu_start+gpu_count-1); if(i.ne.0) then; write(jo,'("Failed!")'); call c_proc_quit(1); return; endif
        tm=thread_wtime()-tm; write(jo,'("Ok(",F4.1," sec): Host argument buffer size (B) = ",i11)') tm,arg_buf_size_host
        if(arg_buf_size_host.lt.min_arg_buf_size) then
         write(jo,'("#FATAL(c_process:c_proc_life): Host argument buffer size is lower than minimally allowed: ",i11)') min_arg_buf_size
         call c_proc_quit(2); return
        endif
  !Check Host argument buffer levels:
        i=get_blck_buf_sizes_host(blck_sizes)
        if(i.le.0.or.i.gt.max_arg_buf_levels) then
         write(jo,'("#ERROR(c_process:c_proc_life): Invalid number of Host argument buffer levels: ",i9,1x,i9)') max_arg_buf_levels,i
         call c_proc_quit(3); return
        else
         write(jo,'("#MSG(c_process:c_proc_life): Number of Host argument buffer levels = ",i4,":")') i
         do j=0,i-1
          write(jo,'("#MSG(c_process:c_proc_life): Level ",i4,": Size (B) = ",i11)') j,blck_sizes(j)
         enddo
        endif
        write(jo,'("#MSG(c_process:c_proc_life): Max number of arguments in the Host argument buffer = ",i7)') max_args_host
  !Check GPU argument buffer levels:
#ifndef NO_GPU
        do j=gpu_start,gpu_start+gpu_count-1
         if(gpu_is_mine(j).ne.NOT_REALLY) then
          i=get_blck_buf_sizes_gpu(j,blck_sizes)
          if(i.le.0.or.i.gt.max_arg_buf_levels) then
           write(jo,'("#ERROR(c_process:c_proc_life): Invalid number of GPU argument buffer levels: ",i9,1x,i9)') max_arg_buf_levels,i
           call c_proc_quit(4); return
          else
           write(jo,'("#MSG(c_process:c_proc_life): Number of GPU#",i2," argument buffer levels = ",i4,":")') j,i
           do k=0,i-1
            write(jo,'("#MSG(c_process:c_proc_life): Level ",i4,": Size (B) = ",i11)') k,blck_sizes(k)
           enddo
          endif
         endif
        enddo
#endif
 !Init AAL:
        write(jo,'("#MSG(c_process:c_proc_life): Allocating Active Argument List (AAL) ... ")',advance='no')
        tm=thread_wtime()
        act_args_size=max_args_host
        allocate(act_args(1:act_args_size),STAT=ierr); if(ierr.ne.0) then; write(jo,'("#ERROR(c_process:c_proc_life): Active argument list allocation failed!")'); call c_proc_quit(5); return; endif
        act_args_lim=0
        tm=thread_wtime()-tm; write(jo,'("Ok(",F4.1," sec):  AAL size = ",i7)') tm,act_args_size
 !Init TIQ:

!Test C-process functionality:
        call c_proc_test(ierr); if(ierr.ne.0) then; write(jo,'("#ERROR(c_process:c_proc_life): C-process functionality test failed: ",i7)') ierr; call c_proc_quit(6); return; endif
        call run_benchmarks(ierr); if(ierr.ne.0) then; write(jo,'("#ERROR(c_process:c_proc_life): C-process benchmarking failed: ",i7)') ierr; call c_proc_quit(7); return; endif
!Report to work to the local host:

!Receive the next batch of tensor instructions:

!Prioritize/reprioritize tensor instructions in the queue:

!If no new remote arguments have arrived, schedule some local instructions for execution (local hedge):

!Put MPI requests for remote arguments:

!Schedule those tensor instructions which have all their remote arguments delivered:

!Check the status of previously scheduled non-blocking instructions:

!Post MPI sends if remote tensor blocks have been computed:

        write(jo,'("#MSG(c_process:c_proc_life): Cleaning ... ")',advance='no')
        call c_proc_quit(0); if(ierr.eq.0) write(jo,'("Ok")')
        return

        contains

         subroutine c_proc_quit(errc)
         integer, intent(in):: errc
         integer(C_INT) j0,j1
         ierr=errc
         j0=arg_buf_clean_host(); if(j0.ne.0) then; write(jo,'("#WARNING(c_process:c_proc_life:c_proc_quit): Host Argument buffer is not clean!")'); ierr=ierr+100; endif
#ifndef NO_GPU
         do j1=gpu_start,gpu_start+gpu_count-1
          if(gpu_is_mine(j1).ne.NOT_REALLY) then
           j0=arg_buf_clean_gpu(j1); if(j0.ne.0) then; write(jo,'("#WARNING(c_process:c_proc_life:c_proc_quit): GPU#",i2," Argument buffer is not clean!")') j1; ierr=ierr+100*(2+j1); endif
          endif
         enddo
#endif
         j0=arg_buf_deallocate(gpu_start,gpu_start+gpu_count-1); arg_buf_size_host=0; max_args_host=0
         if(j0.ne.0) then; write(j0,'("#ERROR(c_process:c_proc_life:c_proc_quit): Deallocation of argument buffers failed!")'); ierr=ierr+10000; endif
         if(allocated(act_args)) deallocate(act_args); act_args_lim=0; act_args_size=0
         if(allocated(tens_instr_que)) deallocate(tens_instr_que); tens_instr_que_lim=0; tens_instr_que_size=0
         return
         end subroutine c_proc_quit

        end subroutine c_proc_life
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
!	write(jo,'("#DEBUG(c_process:tens_blck_pack): Creating a tensor block packet:")') !debug
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
!	write(jo,'(" Data kind    : ",i9,1x,i9)') off_data_kind,data_kind
!	write(jo,'(" Element count: ",i9,1x,i9)') off_elems_count,elems_count
!	write(jo,'(" Tensor rank  : ",i9,1x,i9)') off_trank,trank
!	write(jo,'(" Dims         :",i9,3x,32(1x,i4))') off_dims,dims(1:trank)
!	write(jo,'(" Divs         :",i9,3x,32(1x,i4))') off_divs,divs(1:trank)
!	write(jo,'(" Grps         :",i9,3x,32(1x,i4))') off_grps,grps(1:trank)
!	write(jo,'(" Base         :",i9,3x,32(1x,i4))') off_base,base(1:trank)
!	write(jo,'(" Scalar       : ",i9,3x,D25.14,1x,D25.14)') off_sclr,sclr
!	select case(dtk)
!	case('r4'); write(jo,'(" Elements(r4) : ",i9,1x,i9)') off_elems_r4
!	case('r8'); write(jo,'(" Elements(r8) : ",i9,1x,i9)') off_elems_r8
!	case('c8'); write(jo,'(" Elements(c8) : ",i9,1x,i9)') off_elems_c8
!	end select
!	write(jo,'(" Packet size  : ",i9,1x,i9)') off_packet_size,packet_size
!DEBUG end.
!Allocate an argument buffer space on Host:
        err_code=get_buf_entry_host(packet_size,pptr,entry_num); if(err_code.ne.0) then; ierr=5; return; endif
!	write(jo,'("#DEBUG(c_process:tens_blck_pack): Host argument buffer entry obtained: ",i7)') entry_num !debug
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
!	write(jo,'("#DEBUG(c_process:tens_blck_pack): packet created (size/entry): ",i12,1x,i7)') packet_size,entry_num !debug
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

        ierr=0; err_code=0
!	write(jo,'("#DEBUG(c_process:tens_blck_upack): Unpacking a tensor block packet:")') !debug
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
        if(associated(tens%tensor_shape%dim_extent)) then; deallocate(tens%tensor_shape%dim_extent,STAT=j); if(j.ne.0) nullify(tens%tensor_shape%dim_extent); endif
        if(trank.gt.0) then; allocate(tens%tensor_shape%dim_extent(1:trank),STAT=ierr); if(ierr.ne.0) then; ierr=3; return; endif; endif
        if(associated(tens%tensor_shape%dim_divider)) then; deallocate(tens%tensor_shape%dim_divider,STAT=j); if(j.ne.0) nullify(tens%tensor_shape%dim_divider); endif
        if(trank.gt.0) then; allocate(tens%tensor_shape%dim_divider(1:trank),STAT=ierr); if(ierr.ne.0) then; ierr=4; return; endif; endif
        if(associated(tens%tensor_shape%dim_group)) then; deallocate(tens%tensor_shape%dim_group,STAT=j); if(j.ne.0) nullify(tens%tensor_shape%dim_group); endif
        if(trank.gt.0) then; allocate(tens%tensor_shape%dim_group(1:trank),STAT=ierr); if(ierr.ne.0) then; ierr=5; return; endif; endif
        select case(data_kind)
        case(R4)
         if(associated(tens%data_real4)) then; deallocate(tens%data_real4,STAT=j); if(j.ne.0) nullify(tens%data_real4); endif
         if(trank.gt.0) then; allocate(tens%data_real4(0:elems_count-1),STAT=ierr); if(ierr.ne.0) then; ierr=6; return; endif; endif
        case(R8)
         if(associated(tens%data_real8)) then; deallocate(tens%data_real8,STAT=j); if(j.ne.0) nullify(tens%data_real8); endif
         if(trank.gt.0) then; allocate(tens%data_real8(0:elems_count-1),STAT=ierr); if(ierr.ne.0) then; ierr=7; return; endif; endif
        case(C8)
         if(associated(tens%data_cmplx8)) then; deallocate(tens%data_cmplx8,STAT=j); if(j.ne.0) nullify(tens%data_cmplx8); endif
         if(trank.gt.0) then; allocate(tens%data_cmplx8(0:elems_count-1),STAT=ierr); if(ierr.ne.0) then; ierr=8; return; endif; endif
        case default
         ierr=9; return
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
          if(off_elems_r4+elems_count*sizeof(elems_r4).ne.packet_size) then; ierr=10; return; endif
          call c_f_pointer(c_addr,ptr_elems_r4,shape=[elems_count])
          do s0=0,elems_count-1; tens%data_real4(s0)=ptr_elems_r4(1+s0); enddo
          nullify(ptr_elems_r4)
         endif
        case(R8)
         off_elems_r8=off_sclr+sizeof(elems_c8)
         if(trank.gt.0) then
          c_addr=ptr_offset(pptr,off_elems_r8)
          if(off_elems_r8+elems_count*sizeof(elems_r8).ne.packet_size) then; ierr=11; return; endif
          call c_f_pointer(c_addr,ptr_elems_r8,shape=[elems_count])
          do s0=0,elems_count-1; tens%data_real8(s0)=ptr_elems_r8(1+s0); enddo
          nullify(ptr_elems_r8)
         endif
        case(C8)
         off_elems_c8=off_sclr+sizeof(elems_c8)
         if(trank.gt.0) then
          c_addr=ptr_offset(pptr,off_elems_c8)
          if(off_elems_c8+elems_count*sizeof(elems_c8).ne.packet_size) then; ierr=12; return; endif
          call c_f_pointer(c_addr,ptr_elems_c8,shape=[elems_count])
          do s0=0,elems_count-1; tens%data_cmplx8(s0)=ptr_elems_c8(1+s0); enddo
          nullify(ptr_elems_c8)
         endif
        end select
!DEBUG begin:
!	write(jo,'(" Packet size  : ",i9,1x,i9)') off_packet_size,packet_size
!	write(jo,'(" Data kind    : ",i9,1x,i9)') off_data_kind,data_kind
!	write(jo,'(" Element count: ",i9,1x,i9)') off_elems_count,elems_count
!	write(jo,'(" Tensor rank  : ",i9,1x,i9)') off_trank,trank
!	write(jo,'(" Dims         :",i9,3x,32(1x,i4))') off_dims,tens%tensor_shape%dim_extent(1:trank)
!	write(jo,'(" Divs         :",i9,3x,32(1x,i4))') off_divs,tens%tensor_shape%dim_divider(1:trank)
!	write(jo,'(" Grps         :",i9,3x,32(1x,i4))') off_grps,tens%tensor_shape%dim_group(1:trank)
!	write(jo,'(" Base         :",i9,3x,32(1x,i4))') off_base,base(1:trank)
!	write(jo,'(" Scalar       : ",i9,3x,D25.14,1x,D25.14)') off_sclr,tens%scalar_value
!	select case(data_kind)
!	case(R4); write(jo,'(" Elements(r4) : ",i9,1x,i9)') off_elems_r4
!	case(R8); write(jo,'(" Elements(r8) : ",i9,1x,i9)') off_elems_r8
!	case(C8); write(jo,'(" Elements(c8) : ",i9,1x,i9)') off_elems_c8
!	end select
!	write(jo,'("#DEBUG(c_process:tens_blck_unpack): packet unpacked (size): ",i12)') packet_size !debug
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
!of the Host Argument Buffer cannot be freed during the lifetime of {tensor_block_t|tensBlck_t}.
!INPUT:
! - pptr - C pointer to a tensor block packet located in the Host argument buffer;
! - tens - an allocated (uninitialized) instance of tensor_block_t (F);
! - ctens - an allocated (uninitialized) instance of tensBlck_t (C);
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
!   For tensor_block_t (Fortran) one should still invoke <tensor_block_destroy>.
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
         write(jo,'("#FATAL(c_process:tens_blck_assoc): attempt to initialize a GPU-resident tensor block in GPU-free code compilation!")') !trap
         ierr=-1; return !attempt to initialize tensBlck_t in a GPU-free code compilation
#endif
        elseif(present(tens).and.(.not.present(ctens)).and.(.not.present(gpu_num))) then
 !Create an instance of tensor_block_t (Fortran):
         tens%tensor_block_size=elems_count
         tens%tensor_shape%num_dim=trank
         nullify(tens%tensor_shape%dim_extent); call c_f_pointer(addr_dims,tens%tensor_shape%dim_extent,shape=[trank])
         nullify(tens%tensor_shape%dim_divider); call c_f_pointer(addr_divs,tens%tensor_shape%dim_divider,shape=[trank])
         nullify(tens%tensor_shape%dim_group); call c_f_pointer(addr_grps,tens%tensor_shape%dim_group,shape=[trank])
         tens%scalar_value=elems_c8
         select case(data_kind)
         case(R4)
          call c_f_pointer(addr_host,ptr_elems_r4,shape=[elems_count])
          tens%data_real4(0:)=>ptr_elems_r4; nullify(ptr_elems_r4) !Is this portable?`
         case(R8)
          call c_f_pointer(addr_host,ptr_elems_r8,shape=[elems_count])
          tens%data_real8(0:)=>ptr_elems_r8; nullify(ptr_elems_r8) !Is this portable?`
         case(C8)
          call c_f_pointer(addr_host,ptr_elems_c8,shape=[elems_count])
          tens%data_cmplx8(0:)=>ptr_elems_c8; nullify(ptr_elems_c8) !Is this portable?`
         end select
        else
         ierr=99 !both <ctens> and <tens> cannot be absent/present simultaneously
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
        write(jo,'("#FATAL(c_process:tens_blck_dissoc): attempt to dissociate a GPU-resident tensor block in GPU-free code compilation!")') !trap
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

        ierr=0; write(jo,'("#MSG(c_process:c_proc_test): Testing C-process functionality ... ")',advance='no')
!TEST 0: Random contraction patterns:
!        write(jo,*)''
!        do i=1,16
!         call contr_pattern_rnd(8,100000000,shape0,shape1,shape2,cptrn,ierr)
!         if(ierr.ne.0) then; write(jo,*)'ERROR(c_process:c_proc_test): contr_pattern_rnd error ',ierr; return; endif
!         l=index(shape0,')'); k=index(shape1,')'); j=index(shape2,')')
!         call printl(jo,shape0(1:l)//'+='//shape1(1:k)//'*'//shape2(1:j)) !debug
!         m=tensor_shape_rank(shape1(1:k),ierr)+tensor_shape_rank(shape2(1:j),ierr)
!         write(jo,'("CPTRN:",64(1x,i2))') cptrn(1:m) !debug
!        enddo
!        return
!TEST 1 (CPU: tensor_block_contract):
        write(jo,'("1 ")',advance='no')
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
        gpu_id=gpu_busy_least(); call cudaSetDevice(gpu_id,err_code); if(err_code.ne.0) then; ierr=666; goto 999; endif
        call gpu_set_event_policy(EVENTS_ON)
!TEST 2 (GPU: matrix multiplication):
        write(jo,'("2 ")',advance='no')
        n=2; ext_beg(1:n)=1
        call tensor_block_create('(345,697)','r8',ftens(0),ierr,val_r8=0d0); if(ierr.ne.0) then; ierr=10; goto 999; endif
        call tensor_block_create('(37,345)','r8',ftens(1),ierr); if(ierr.ne.0) then; ierr=11; goto 999; endif
        call tensor_block_create('(37,697)','r8',ftens(2),ierr); if(ierr.ne.0) then; ierr=12; goto 999; endif
        call tensor_block_create('(345,697)','r4',ftens(3),ierr,val_r4=0.0); if(ierr.ne.0) then; ierr=13; goto 999; endif
        call tensor_block_contract((/-1,1,-1,2/),ftens(1),ftens(2),ftens(0),ierr,'r8'); if(ierr.ne.0) then; ierr=14; goto 999; endif
        call tensor_block_sync(ftens(0),'r8',ierr,'r4'); if(ierr.ne.0) then; ierr=15; goto 999; endif
        call tensor_block_sync(ftens(1),'r8',ierr,'r4'); if(ierr.ne.0) then; ierr=16; goto 999; endif
        call tensor_block_sync(ftens(2),'r8',ierr,'r4'); if(ierr.ne.0) then; ierr=17; goto 999; endif
!        call tensor_block_print(jo,'Matrix0:',ext_beg,ftens(0),ierr,'r4') !debug
        err_code=gpu_matrix_multiply_tn_r4(345_CZ,697_CZ,37_CZ,ftens(1)%data_real4,ftens(2)%data_real4,ftens(3)%data_real4); if(err_code.ne.0) then; write(jo,*)'GPU ERROR ',err_code; ierr=18; goto 999; endif
        i1=gpu_get_error_count(); if(i1.ne.0) then; write(jo,*)'GPU error count = ',i1; ierr=19; goto 999; endif
!        call tensor_block_print(jo,'Matrix3:',ext_beg,ftens(3),ierr,'r4') !debug
        cmp=tensor_block_cmp(ftens(0),ftens(3),ierr,'r4',cmp_thresh=1d-5); if(ierr.ne.0) then; ierr=20; goto 999; endif
        if(.not.cmp) then; ierr=21; goto 999; endif
        call tensor_block_destroy(ftens(3),ierr); if(ierr.ne.0) then; ierr=22; goto 999; endif
        call tensor_block_destroy(ftens(2),ierr); if(ierr.ne.0) then; ierr=23; goto 999; endif
        call tensor_block_destroy(ftens(1),ierr); if(ierr.ne.0) then; ierr=24; goto 999; endif
        call tensor_block_destroy(ftens(0),ierr); if(ierr.ne.0) then; ierr=25; goto 999; endif
!TEST 3 (GPU: tensor packing/unpacking/association, gpu_tensor_transpose):
        write(jo,'("3 ")',advance='no')
        n=4; ext_beg(1:n)=1; o2n(0:n)=(/+1,(j,j=1,n)/); ctrl=1
        do i=1,ifcl(n)
         call trng(ctrl,n,o2n,ngt); !write(jo,'(32(1x,i2))') o2n(0:n)
         call permutation_converter(.false.,n,n2o,o2n); !write(jo,'(32(1x,i2))') n2o(0:n)
         call tensor_block_create('(17,27,33,44)','r8',ftens(0),ierr); if(ierr.ne.0) then; ierr=26; goto 999; endif
         call tensor_block_create('(44,33,27,17)','r8',ftens(1),ierr,val_r8=0d0); if(ierr.ne.0) then; ierr=27; goto 999; endif
         tm0=thread_wtime()
         call tensor_block_copy(ftens(0),ftens(1),ierr,o2n); if(ierr.ne.0) then; ierr=28; goto 999; endif
!         write(jo,'("#DEBUG(c_process:c_proc_test): CPU tensor transpose time = ",F10.4)') thread_wtime()-tm0 !debug
         call tensor_block_sync(ftens(0),'r8',ierr,'r4'); if(ierr.ne.0) then; ierr=29; goto 999; endif
         call tensor_block_sync(ftens(1),'r8',ierr,'r4'); if(ierr.ne.0) then; ierr=30; goto 999; endif
         call tens_blck_pack(ftens(1),'r4',pack_size(1),entry_ptr(1),entry_num(1),ierr); if(ierr.ne.0) then; ierr=31; goto 999; endif
         call tens_blck_unpack(ftens(2),entry_ptr(1),ierr); if(ierr.ne.0) then; ierr=32; goto 999; endif
         cmp=tensor_block_cmp(ftens(1),ftens(2),ierr,'r4',cmp_thresh=1d-6); if(ierr.ne.0) then; ierr=33; goto 999; endif
         if(.not.cmp) then; ierr=34; goto 999; endif
         call tensor_block_destroy(ftens(2),ierr); if(ierr.ne.0) then; ierr=35; goto 999; endif
         call tensor_block_create('(17,27,33,44)','r4',ftens(2),ierr,val_r4=0.0); if(ierr.ne.0) then; ierr=36; goto 999; endif
         call tens_blck_pack(ftens(2),'r4',pack_size(2),entry_ptr(2),entry_num(2),ierr); if(ierr.ne.0) then; ierr=37; goto 999; endif
         call tens_blck_assoc(entry_ptr(1),ierr,ctens=ctens(1),gpu_num=gpu_id); if(ierr.ne.0) then; ierr=38; goto 999; endif
         call tens_blck_assoc(entry_ptr(2),ierr,ctens=ctens(2),gpu_num=gpu_id); if(ierr.ne.0) then; ierr=39; goto 999; endif
         tm0=thread_wtime()
         err_code=gpu_tensor_block_copy_dlf(n2o,ctens(1),ctens(2)); if(err_code.ne.0) then; write(jo,*)'GPU ERROR ',err_code; ierr=40; goto 999; endif
!         write(jo,'("#DEBUG(c_process:c_proc_test): GPU tensor transpose time = ",F10.4)') thread_wtime()-tm0 !debug
         i1=gpu_get_error_count(); if(i1.ne.0) then; write(jo,*)'GPU error count = ',i1; ierr=41; goto 999; endif
!         call print_gpu_debug_dump(jo) !debug
         if(tensor_block_norm2(ftens(2),ierr,'r4').ne.0d0) then; ierr=42; goto 999; endif
         call tens_blck_unpack(ftens(2),entry_ptr(2),ierr); if(ierr.ne.0) then; ierr=43; goto 999; endif
!         write(jo,*) tensor_block_norm2(ftens(0),ierr,'r4'),tensor_block_norm2(ftens(2),ierr,'r4') !debug
!         call tensor_block_print(jo,'Tensor0:',ext_beg,ftens(0),ierr,'r4') !debug
!         call tensor_block_print(jo,'Tensor1:',ext_beg,ftens(1),ierr,'r4') !debug
!         call tensor_block_print(jo,'Tensor2:',ext_beg,ftens(2),ierr,'r4') !debug
         cmp=tensor_block_cmp(ftens(0),ftens(2),ierr,'r4',cmp_thresh=1d-5,diff_count=diffc); if(ierr.ne.0) then; ierr=44; goto 999; endif
         if(.not.cmp) then; write(jo,*)'DIFF COUNT = ',diffc; ierr=45; goto 999; endif
         call tens_blck_dissoc(ctens(2),ierr); if(ierr.ne.0) then; ierr=46; goto 999; endif
         call tens_blck_dissoc(ctens(1),ierr); if(ierr.ne.0) then; ierr=47; goto 999; endif
         err_code=free_buf_entry_host(entry_num(2)); if(err_code.ne.0) then; ierr=48; goto 999; endif
         err_code=free_buf_entry_host(entry_num(1)); if(err_code.ne.0) then; ierr=49; goto 999; endif
         call tensor_block_destroy(ftens(2),ierr); if(ierr.ne.0) then; ierr=50; goto 999; endif
         call tensor_block_destroy(ftens(1),ierr); if(ierr.ne.0) then; ierr=51; goto 999; endif
         call tensor_block_destroy(ftens(0),ierr); if(ierr.ne.0) then; ierr=52; goto 999; endif
        enddo
!TEST 4 (GPU: gpu_tensor_contraction):
        write(jo,'("4 ")',advance='no')
        call tensor_block_create('(3,4,31,2,37,8)','r4',ftens(0),ierr,val_r4=0.0); if(ierr.ne.0) then; ierr=53; goto 999; endif
        call tensor_block_create('(2,31,7,8,11,29)','r4',ftens(1),ierr,val_r4=1.0); if(ierr.ne.0) then; ierr=54; goto 999; endif
        call tensor_block_create('(11,37,4,7,29,3)','r4',ftens(2),ierr,val_r4=1.0); if(ierr.ne.0) then; ierr=55; goto 999; endif
        call tensor_block_create('(3,4,31,2,37,8)','r4',ftens(3),ierr,val_r4=0.0); if(ierr.ne.0) then; ierr=56; goto 999; endif
        call tensor_block_create('(2,31,7,8,11,29)','r4',ftens(4),ierr,val_r4=0.0); if(ierr.ne.0) then; ierr=57; goto 999; endif
        call tensor_block_create('(11,37,4,7,29,3)','r4',ftens(5),ierr,val_r4=0.0); if(ierr.ne.0) then; ierr=58; goto 999; endif
        call tensor_block_create('(3,4,31,2,37,8)','r4',ftens(6),ierr,val_r4=0.0); if(ierr.ne.0) then; ierr=59; goto 999; endif
        cptrn(1:12)=(/4,3,-4,6,-1,-5,-5,5,2,-3,-6,1/)
        tm0=thread_wtime()
        call tensor_block_contract(cptrn,ftens(1),ftens(2),ftens(3),ierr,'r4'); if(ierr.ne.0) then; ierr=62; goto 999; endif
!       write(jo,'("#DEBUG(c_process:c_proc_test): CPU tensor contraction time = ",F10.4)') thread_wtime()-tm0 !debug
        call tens_blck_pack(ftens(0),'r4',pack_size(0),entry_ptr(0),entry_num(0),ierr); if(ierr.ne.0) then; ierr=64; goto 999; endif
        call tens_blck_pack(ftens(1),'r4',pack_size(1),entry_ptr(1),entry_num(1),ierr); if(ierr.ne.0) then; ierr=65; goto 999; endif
        call tens_blck_pack(ftens(2),'r4',pack_size(2),entry_ptr(2),entry_num(2),ierr); if(ierr.ne.0) then; ierr=66; goto 999; endif
        call tens_blck_pack(ftens(4),'r4',pack_size(4),entry_ptr(4),entry_num(4),ierr); if(ierr.ne.0) then; ierr=67; goto 999; endif
        call tens_blck_pack(ftens(5),'r4',pack_size(5),entry_ptr(5),entry_num(5),ierr); if(ierr.ne.0) then; ierr=68; goto 999; endif
        call tens_blck_pack(ftens(6),'r4',pack_size(6),entry_ptr(6),entry_num(6),ierr); if(ierr.ne.0) then; ierr=69; goto 999; endif
        call tens_blck_assoc(entry_ptr(0),ierr,ctens=ctens(0),gpu_num=gpu_id); if(ierr.ne.0) then; ierr=70; goto 999; endif
        call tens_blck_assoc(entry_ptr(1),ierr,ctens=ctens(1),gpu_num=gpu_id); if(ierr.ne.0) then; ierr=71; goto 999; endif
        call tens_blck_assoc(entry_ptr(2),ierr,ctens=ctens(2),gpu_num=gpu_id); if(ierr.ne.0) then; ierr=72; goto 999; endif
        call tens_blck_assoc(entry_ptr(4),ierr,ctens=ctens(4),gpu_num=gpu_id); if(ierr.ne.0) then; ierr=73; goto 999; endif
        call tens_blck_assoc(entry_ptr(5),ierr,ctens=ctens(5),gpu_num=gpu_id); if(ierr.ne.0) then; ierr=74; goto 999; endif
        call tens_blck_assoc(entry_ptr(6),ierr,ctens=ctens(6),gpu_num=gpu_id); if(ierr.ne.0) then; ierr=75; goto 999; endif
        o2n(0:6)=(/+1,1,2,3,4,5,6/);
        err_code=gpu_tensor_block_copy_dlf(o2n,ctens(1),ctens(4)); if(err_code.ne.0) then; ierr=76; goto 999; endif
        err_code=gpu_tensor_block_copy_dlf(o2n,ctens(2),ctens(5)); if(err_code.ne.0) then; ierr=77; goto 999; endif
        err_code=cuda_task_create(cuda_task(1)); if(err_code.ne.0) then; ierr=78; goto 999; endif
        err_code=cuda_task_create(cuda_task(2)); if(err_code.ne.0) then; ierr=79; goto 999; endif
        tm0=thread_wtime()
        call gpu_set_matmult_algorithm(BLAS_OFF)
        err_code=gpu_tensor_block_contract_dlf_(cptrn,ctens(1),ctens(2),ctens(0),COPY_BACK,cuda_task(1)); if(err_code.ne.0) then; write(jo,*)'GPU ERROR = ',err_code; ierr=80; goto 999; endif
        call gpu_set_matmult_algorithm(BLAS_ON)
        err_code=gpu_tensor_block_contract_dlf_(cptrn,ctens(4),ctens(5),ctens(6),COPY_BACK,cuda_task(2)); if(err_code.ne.0) then; write(jo,*)'GPU ERROR = ',err_code; ierr=81; goto 999; endif
        i1=0
        do while(i1.ne.3)
         if(cuda_task_complete(cuda_task(1)).ne.cuda_task_scheduled.and.(i1.eq.0.or.i1.eq.2)) then; tm1=thread_wtime()-tm0; i1=i1+1; endif
         if(cuda_task_complete(cuda_task(2)).ne.cuda_task_scheduled.and.(i1.eq.0.or.i1.eq.1)) then; tm2=thread_wtime()-tm0; i1=i1+2; endif
        enddo
!       write(jo,*)'GPU Times/Status = ',tm1,tm2,cuda_task_complete(cuda_task(1)),cuda_task_complete(cuda_task(2)) !debug
        err_code=cuda_task_destroy(cuda_task(2)); if(err_code.ne.0) then; ierr=82; goto 999; endif
        err_code=cuda_task_destroy(cuda_task(1)); if(err_code.ne.0) then; ierr=83; goto 999; endif
        call tens_blck_unpack(ftens(0),entry_ptr(0),ierr); if(ierr.ne.0) then; ierr=84; goto 999; endif
        call tens_blck_unpack(ftens(6),entry_ptr(6),ierr); if(ierr.ne.0) then; ierr=85; goto 999; endif
!        write(jo,*)'1-NORMS:',tensor_block_norm1(ftens(0),ierr,'r4'),tensor_block_norm1(ftens(3),ierr,'r4'),tensor_block_norm1(ftens(6),ierr,'r4') !debug
!        write(jo,*) ftens(0)%data_real4(0:5); write(jo,*) ftens(3)%data_real4(0:5); write(jo,*) ftens(6)%data_real4(0:5) !debug
        cmp=tensor_block_cmp(ftens(0),ftens(3),ierr,'r4',.true.,1d-3,diffc); if(ierr.ne.0) then; ierr=86; goto 999; endif
        if(.not.cmp) then; write(jo,*)'DIFF(1) COUNT = ',diffc; ierr=87; goto 999; endif
        cmp=tensor_block_cmp(ftens(6),ftens(3),ierr,'r4',.true.,1d-3,diffc); if(ierr.ne.0) then; ierr=88; goto 999; endif
        if(.not.cmp) then; write(jo,*)'DIFF(2) COUNT = ',diffc; ierr=89; goto 999; endif
        call tens_blck_dissoc(ctens(6),ierr); if(ierr.ne.0) then; ierr=90; goto 999; endif
        call tens_blck_dissoc(ctens(5),ierr); if(ierr.ne.0) then; ierr=91; goto 999; endif
        call tens_blck_dissoc(ctens(4),ierr); if(ierr.ne.0) then; ierr=92; goto 999; endif
        call tens_blck_dissoc(ctens(2),ierr); if(ierr.ne.0) then; ierr=93; goto 999; endif
        call tens_blck_dissoc(ctens(1),ierr); if(ierr.ne.0) then; ierr=94; goto 999; endif
        call tens_blck_dissoc(ctens(0),ierr); if(ierr.ne.0) then; ierr=95; goto 999; endif
        err_code=free_buf_entry_host(entry_num(6)); if(err_code.ne.0) then; ierr=96; goto 999; endif
        err_code=free_buf_entry_host(entry_num(5)); if(err_code.ne.0) then; ierr=97; goto 999; endif
        err_code=free_buf_entry_host(entry_num(4)); if(err_code.ne.0) then; ierr=98; goto 999; endif
        err_code=free_buf_entry_host(entry_num(2)); if(err_code.ne.0) then; ierr=99; goto 999; endif
        err_code=free_buf_entry_host(entry_num(1)); if(err_code.ne.0) then; ierr=100; goto 999; endif
        err_code=free_buf_entry_host(entry_num(0)); if(err_code.ne.0) then; ierr=101; goto 999; endif
        call tensor_block_destroy(ftens(6),ierr); if(ierr.ne.0) then; ierr=102; goto 999; endif
        call tensor_block_destroy(ftens(5),ierr); if(ierr.ne.0) then; ierr=103; goto 999; endif
        call tensor_block_destroy(ftens(4),ierr); if(ierr.ne.0) then; ierr=104; goto 999; endif
        call tensor_block_destroy(ftens(3),ierr); if(ierr.ne.0) then; ierr=105; goto 999; endif
        call tensor_block_destroy(ftens(2),ierr); if(ierr.ne.0) then; ierr=106; goto 999; endif
        call tensor_block_destroy(ftens(1),ierr); if(ierr.ne.0) then; ierr=107; goto 999; endif
        call tensor_block_destroy(ftens(0),ierr); if(ierr.ne.0) then; ierr=108; goto 999; endif
#endif
!Clean:
999     do i=0,test_args_lim; err_code=free_buf_entry_host(entry_num(i)); enddo
        do i=0,test_args_lim; call tensor_block_destroy(ftens(i),j); enddo
        if(ierr.eq.0) then; write(jo,'("Ok")'); else; write(jo,'("Failed: ERROR #",i4)') ierr; endif
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

        ierr=0; write(jo,'("#MSG(c_process:c_proc_test): Tensor algebra benchmarks:")')
#ifndef NO_GPU
        gpu_id=gpu_busy_least(); call cudaSetDevice(gpu_id,err_code); if(err_code.ne.0) then; ierr=666; goto 999; endif
        call gpu_set_event_policy(EVENTS_ON)
#endif
!TENSOR TRANSPOSE:
        write(jo,'(" TENSOR TRANSPOSE:")')
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
           call printl(jo,'  New Tensor Shape: '//tshape(1:tsl)//': ',.false.)
           write(jo,'(i10,1x,F16.4)') ftens(1)%tensor_block_size,tensor_block_norm1(ftens(1),ierr,'r4')
           write(jo,'(3x)',advance='no')
           tm=thread_wtime()
           call tensor_block_copy(ftens(1),ftens(0),ierr); if(ierr.ne.0) then; ierr=3; goto 999; endif
           write(jo,'("#DEBUG(tensor_algebra:tensor_block_copy_dlf): Direct time ",F10.6)') thread_wtime()-tm !debug
           call random_permutation(tens_rank,o2n)
!           o2n(1:tens_rank)=(/3,2,1/) !debug
           call permutation_converter(.false.,tens_rank,n2o,o2n)
           write(jo,'(3x,"Permutation:",32(1x,i2))') o2n(1:tens_rank)
           write(jo,'(3x,"Permutation:",32(1x,i2))') n2o(1:tens_rank)
           call set_transpose_algorithm(EFF_TRN_OFF) !scatter
           write(jo,'(3x)',advance='no')
           tm=thread_wtime()
           call tensor_block_copy(ftens(1),ftens(0),ierr,o2n); if(ierr.ne.0) then; ierr=4; goto 999; endif
           write(jo,'("#DEBUG(tensor_algebra:tensor_block_copy_scatter_dlf): Time ",F10.6)') thread_wtime()-tm !debug
           write(jo,'(3x)',advance='no')
           tm=thread_wtime()
           call tensor_block_copy(ftens(0),ftens(2),ierr,n2o); if(ierr.ne.0) then; ierr=5; goto 999; endif
           write(jo,'("#DEBUG(tensor_algebra:tensor_block_copy_scatter_dlf): Time ",F10.6)') thread_wtime()-tm !debug
           cmp=tensor_block_cmp(ftens(1),ftens(2),ierr,'r4',.true.,1d-4,diffc); if(ierr.ne.0) then; ierr=6; goto 999; endif
           write(jo,'(3x,l1,1x,i9,1x,F16.4)') cmp,diffc,tensor_block_norm1(ftens(2),ierr,'r4')
           if(.not.cmp) then; nfail=nfail+1; write(jo,'(3x,"Comparison Failed!")'); endif
           call set_transpose_algorithm(EFF_TRN_ON) !cache-efficient
           write(jo,'(3x)',advance='no')
           tm=thread_wtime()
           call tensor_block_copy(ftens(1),ftens(0),ierr,o2n); if(ierr.ne.0) then; ierr=7; goto 999; endif
           write(jo,'("#DEBUG(tensor_algebra:tensor_block_copy_dlf): Time ",F10.6)') thread_wtime()-tm !debug
           write(jo,'(3x)',advance='no')
           tm=thread_wtime()
           call tensor_block_copy(ftens(0),ftens(2),ierr,n2o); if(ierr.ne.0) then; ierr=8; goto 999; endif
           write(jo,'("#DEBUG(tensor_algebra:tensor_block_copy_dlf): Time ",F10.6)') thread_wtime()-tm !debug
           cmp=tensor_block_cmp(ftens(1),ftens(2),ierr,'r4',.true.,1d-4,diffc); if(ierr.ne.0) then; ierr=9; goto 999; endif
           write(jo,'(3x,l1,1x,i9,1x,F16.4)') cmp,diffc,tensor_block_norm1(ftens(2),ierr,'r4')
           if(.not.cmp) then; nfail=nfail+1; write(jo,'(3x,"Comparison Failed!")'); endif
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
           write(jo,'(3x)',advance='no')
           err_code=gpu_tensor_block_copy_dlf(o2n,ctens(1),ctens(0)); if(err_code.ne.0) then; write(jo,*)'GPU error ',err_code; call print_gpu_debug_dump(jo); ierr=18; goto 999; endif
           write(jo,'(3x)',advance='no')
           err_code=gpu_tensor_block_copy_dlf(n2o,ctens(0),ctens(2)); if(err_code.ne.0) then; write(jo,*)'GPU error ',err_code; call print_gpu_debug_dump(jo); ierr=19; goto 999; endif
           call tens_blck_unpack(ftens(2),entry_ptr(2),ierr); if(ierr.ne.0) then; ierr=20; goto 999; endif
           cmp=tensor_block_cmp(ftens(1),ftens(2),ierr,'r4',.true.,1d-4,diffc); if(ierr.ne.0) then; ierr=21; goto 999; endif
           write(jo,'(3x,l1,1x,i9,1x,F16.4)') cmp,diffc,tensor_block_norm1(ftens(2),ierr,'r4')
           if(.not.cmp) then; nfail=nfail+1; write(jo,'(3x,"Comparison Failed!")'); endif
           call gpu_set_transpose_algorithm(EFF_TRN_ON) !shared-memory on GPU
           err_code=cuda_task_create(cuda_task(1)); if(err_code.ne.0) then; ierr=22; goto 999; endif
           write(jo,'(3x)',advance='no')
           tm=thread_wtime()
           err_code=gpu_tensor_block_copy_dlf_(o2n,ctens(1),ctens(0),COPY_BACK,cuda_task(1)); if(err_code.ne.0) then; write(jo,*)'GPU error ',err_code; call print_gpu_debug_dump(jo); ierr=23; goto 999; endif
           err_code=cuda_task_wait(cuda_task(1)); if(err_code.ne.cuda_task_completed) then; ierr=24; goto 999; endif
           write(jo,'("#DEBUG(tensor_algebra_gpu_nvidia:gpu_tensor_block_copy_dlf_): Time ",F10.6)') thread_wtime()-tm
           write(jo,'(3x)',advance='no')
           err_code=gpu_tensor_block_copy_dlf(n2o,ctens(0),ctens(2)); if(err_code.ne.0) then; write(jo,*)'GPU error ',err_code; call print_gpu_debug_dump(jo); ierr=25; goto 999; endif
           err_code=cuda_task_destroy(cuda_task(1)); if(err_code.ne.0) then; ierr=26; goto 999; endif
           call tens_blck_unpack(ftens(2),entry_ptr(2),ierr); if(ierr.ne.0) then; ierr=27; goto 999; endif
           cmp=tensor_block_cmp(ftens(1),ftens(2),ierr,'r4',.true.,1d-4,diffc); if(ierr.ne.0) then; ierr=28; goto 999; endif
           write(jo,'(3x,l1,1x,i9,1x,F16.4)') cmp,diffc,tensor_block_norm1(ftens(2),ierr,'r4')
           if(.not.cmp) then; nfail=nfail+1; write(jo,'(3x,"Comparison Failed!")'); endif
           write(jo,'(3x)',advance='no')
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
         write(jo,'("Done: ",i6," failed transposes.")') nfail
        else
         write(jo,'("Failed: Error #",i7)') ierr
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
        write(ifh,*)' '; write(ifh,'("#DEBUG(c_process:print_gpu_debug_dump): GPU DEBUG DUMP:")')
        dump_size=gpu_get_debug_dump(dump)
        do i=1,dump_size,10; write(ifh,'(10(1x,i9))') dump(i:min(dump_size,i+9)); enddo
        return
        end subroutine print_gpu_debug_dump
#endif

       end module c_process
