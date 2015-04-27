!This module provides general services for MPI parallel programs.
       module service_mpi
        use, intrinsic:: ISO_C_BINDING, only: C_INT,C_SIZE_T
!If a function contains MPI calls it is classified as PARALLEL (PARALLEL: YES).
!FUNCTIONS:
! - subroutine file_handle(ch*:command,i:ifh,i:ierr): SERIAL: file handle manager.
! - subroutine quit(i:error_code,ch*:error_msg): PARALLEL: safe global exit.
! - get_memory_status(total_ram,free_ram,used_swap,ierr): SERIAL: Host RAM status.
!--------------------------------------------------------------------------------
!Parallel environment:
#ifdef USE_MPI_MOD
        use mpi          !MPI Fortran interface
        implicit none
        public
#else
        implicit none
        public
	include 'mpif.h' !MPI Fortran interface
#endif
!Parameters:
 !MPI kinds:
        integer(C_INT), parameter, public:: INT_MPI=4                 !default MPI integer kind
        integer(C_INT), parameter, public:: INT_ADDR=MPI_ADDRESS_KIND !default MPI address/size kind
 !File management:
	integer, parameter, private:: max_open_files=1024-16 !maximal amount of open files per process (first 16 file handles [0..15] are reserved)
!Types:
 !GPU info:
	type gpu_info_t
	 integer(C_SIZE_T):: totalGlobalMem
	 integer(C_SIZE_T):: sharedMemPerBlock
	 integer(C_INT):: regsPerBlock
	 integer(C_INT):: warpSize
	 integer(C_INT):: maxThreadsPerBlock
	 integer(C_INT):: maxThreadsDim(3)
	 integer(C_INT):: maxGridSize(3)
	 integer(C_INT):: clockRate
	 integer(C_SIZE_T):: totalConstMem
	 integer(C_INT):: major
	 integer(C_INT):: minor
	 integer(C_INT):: deviceOverlap
	 integer(C_INT):: multiProcessorCount
	 integer(C_INT):: concurrentKernels
	 integer(C_INT):: ECCEnabled
	 integer(C_INT):: asyncEngineCount
	 integer(C_INT):: memoryClockRate
	 integer(C_INT):: memoryBusWidth
	 integer(C_INT):: maxThreadsPerMultiProcessor
	end type gpu_info_t
!Internal Variables:
 !Private:
	integer, private:: j_
 !General:
	integer:: jo=6                             !default output device (6 is the screen)
	integer:: log_file                         !log-file handle (set by runtime)
	integer:: impis                            !MPI communicator size (set by runtime)
	integer:: impir                            !MPI rank of the process [0..impis-1] (set by runtime)
	integer:: max_threads=1                    !max number of threads per process (set by runtime)
	integer:: mpi_procs_per_node=0             !number of MPI processes per node (set by ENVIRONMENT)
	integer:: mpi_proc_id_on_node=0            !internal process ID within a node: [0..mpi_procs_per_node-1] (set by runtime)
	integer(C_INT):: gpus_found=0              !total number of GPUs found on the node (set by runtime)
	integer(C_INT):: gpu_count=0               !number of GPUs assigned to the current process (set by runtime)
	integer(C_INT):: gpu_start=0               !the number of the 1st GPU assigned to the current process: [gpu_start...gpu_start+gpu_count-1] (set by runtime)
	type(gpu_info_t),allocatable:: gpu_info(:) !information about available GPUs (set by runtime)
	integer:: my_role=-1                       !the role of the process (set by runtime)
	integer:: my_group=-1                      !group the process belongs to (set by runtime)
	real(8):: time_begin,time_end              !wall time for the MPI process
	integer:: exec_status=0                    !current execution status (0: success)
	character(MPI_MAX_PROCESSOR_NAME):: proc_name=' ' !processor name (set by runtime)
	integer:: proc_name_len                    !the length of the processor name
	integer:: mpi_thread_provided              !the level of multithreaded MPI service provided (set by runtime)
 !File Management:
	integer, private:: nof=0                          !current number of open files (local to each process)
	integer, private:: fhot(16:16+max_open_files-1)=0 !file handle occupancy table (first 16 file handles [0..15] are reserved)
	integer, private:: ffhs(0:max_open_files-1)=(/(j_,j_=16,16+max_open_files-1)/) !a stack of free file handles

       contains
!-----------------------------------------------
	subroutine file_handle(command,ifh,ierr)
!This subroutine provides the management of file handles.
!The subroutine and corrresponding tables are local to each process.
!INPUT:
! - command - 'get', 'free', 'stop';
! - ifh - file handle to be freed when command='free';
!OUTPUT:
! - ifh - file handle when command='get';
! - ierr - error code (0 - success);
!PARALLEL: NO.
!NOTES:
! - only the first 4 characters of the COMMAND really matter.
        use STSUBS
	implicit none
	character(*), intent(in):: command
	integer, intent(inout):: ifh
	integer, intent(out):: ierr
	integer i,j,k,l,m,n,k0,k1,k2,k3,k4,k5,k6,k7,ks,kf

	ierr=0
	l=len_trim(command)
	if(l.gt.0) then
!GET FILE HANDLE:
	 if(l.eq.3.and.command(1:3).eq.'get') then
	  if(nof.ge.max_open_files) then
	   call printl(jo,'#ERROR(service_mpi:file_handle): max amount of open files exceeded!')
	   ierr=1; goto 999
	  endif
	  ifh=ffhs(nof); fhot(ifh)=1; nof=nof+1
!FREE FILE HANDLE:
	 elseif(l.eq.4.and.command(1:4).eq.'free') then
	  if(nof.le.0) then
	   call printl(jo,'#ERROR(service_mpi:file_handle): attempt to free a file handle when no file handle is used!')
	   ierr=2; goto 999
	  endif
	  if(ifh.lt.16.or.ifh.ge.16+max_open_files) then
	   call printl(jo,'#ERROR(service_mpi:file_handle): invalid file handle!')
	   ierr=3; goto 999
	  endif
	  if(fhot(ifh).eq.0) then
	   call printl(jo,'#ERROR(service_mpi:file_handle): attempt to free an idle file handle!')
	   ierr=4; goto 999
	  endif
	  fhot(ifh)=0; nof=nof-1; ffhs(nof)=ifh
!CLOSE ALL OPENED FILES:
	 elseif(l.eq.4.and.command(1:4).eq.'stop') then
	  do n=16,16+max_open_files-1
	   if(nof.le.0) exit
	   if(fhot(n).ne.0) then; close(n); fhot(n)=0; nof=nof-1; ffhs(nof)=n; endif
	  enddo
	 else
	  call printl(jo,'#ERROR(service_mpi:file_handle): invalid command'//command(1:l)); ierr=5
	 endif
	else
	 call printl(jo,'#ERROR(service_mpi:file_handle): empty command passed!'); ierr=6
	endif
999	return
	end subroutine file_handle
!--------------------------------------------
	subroutine quit(error_code,error_msg)
!This subroutine prints the error message and safely terminates the parallel code execution.
!INPUT:
! - error_code - error code (-1 is the default non-specific error code);
! - error_msg - error message;
!OUTPUT:
! - error information printed by EACH process.
!PARALLEL: YES.
        use STSUBS
	implicit none
	integer i,j,k,l,m,n,k0,k1,k2,k3,k4,k5,k6,k7,ks,kf,ierr
	integer, intent(in):: error_code       !error code
	character(*), intent(in):: error_msg  !error message
!Time:
	time_end=MPI_WTIME()
!Error message:
	l=len_trim(error_msg); if(l.gt.0.and.(error_code.ne.0.or.impir.eq.0)) call printl(jo,error_msg(1:l))
	if(error_code.eq.0) then
	 write(jo,'(''###Process '',i8,'': Success: Wall Time (sec): '',f12.2)') impir,time_end-time_begin
	else
	 write(jo,'(''###Process '',i8,'': Error '',i12,'': Wall Time (sec): '',f12.2)') impir,error_code,time_end-time_begin
	endif
!Close all local files opened by the process:
	l=0; call file_handle('stop',l,ierr)
	if(ierr.ne.0) write(*,'(''Process # '',i8,'': Not all local files have been closed!'')') impir
	if(error_code.ne.0) then
	 l=0; call MPI_ABORT(MPI_COMM_WORLD,l,ierr)
	 if(ierr.ne.0) write(*,'(''Process # '',i8,'': MPI ABORT CODE: '',i12)') impir,ierr
	endif
	call MPI_FINALIZE(ierr)
	stop
	end subroutine quit
!----------------------------------------------------------------------
        subroutine get_memory_status(total_ram,free_ram,used_swap,ierr)
!This subroutine returns:
! - total_ram - total usable RAM available on the node in bytes;
! - free_ram - free usable RAM available on the node in bytes;
! - used_swap - current swap size in bytes;        
        implicit none
        interface
         integer(C_INT) function get_memory_stat(total_ram,free_ram,used_swap) bind(C)
          use, intrinsic:: ISO_C_BINDING
          integer(C_SIZE_T), intent(out):: total_ram
          integer(C_SIZE_T), intent(out):: free_ram
          integer(C_SIZE_T), intent(out):: used_swap
         end function get_memory_stat
        end interface
        integer(C_SIZE_T), intent(out):: total_ram,free_ram,used_swap
        integer, intent(inout):: ierr
        ierr=get_memory_stat(total_ram,free_ram,used_swap)
        return
        end subroutine get_memory_status

       end module service_mpi
