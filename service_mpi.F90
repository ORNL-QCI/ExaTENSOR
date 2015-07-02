!This module provides general services for MPI parallel programs.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2015/07/03
       module service_mpi
        use, intrinsic:: ISO_C_BINDING
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
        integer(C_INT), parameter, public:: INT_MPI=MPI_INTEGER_KIND   !default MPI integer kind
        integer(C_INT), parameter, public:: INT_ADDR=MPI_ADDRESS_KIND  !default MPI address/size kind
        integer(C_INT), parameter, public:: INT_OFFSET=MPI_OFFSET_KIND !default MPI offset kind
        integer(C_INT), parameter, public:: INT_COUNT=MPI_COUNT_KIND   !default MPI element count kind
 !File management:
        integer, parameter, private:: MAX_OPEN_FILES=1024-16 !max amount of open files per process (first 16 file handles [0..15] are reserved)
!Types:
 !NVidia GPU info:
        type gpu_info_t
         integer(C_SIZE_T):: totalGlobalMem !in bytes
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
 !Intel MIC info:
        type mic_info_t
         integer(C_SIZE_T):: totalGlobalMem !in bytes
        end type mic_info_t
!Configuration variables:
 !Private:
        integer, private:: j_
 !General:
        integer:: jo=6     !default output device (6 is the screen)
        integer:: log_file !log-file handle (set during runtime)
 !MPI/OMP:
        integer(INT_MPI):: GLOBAL_MPI_COMM=MPI_COMM_WORLD !defaults to MPI_COMM_WORLD
        integer(INT_MPI):: impis               !global MPI communicator size (set by runtime)
        integer(INT_MPI):: impir               !global MPI rank of the process [0..impis-1] (set by runtime)
        integer(INT_MPI):: mpi_thread_provided !level of multithreaded MPI support provided (set by runtime)
        integer:: max_threads=1         !max number of threads per MPI process (set during runtime)
        integer:: mpi_procs_per_node=0  !number of MPI processes per node (set by ENVIRONMENT)
        integer:: mpi_proc_id_on_node=0 !internal process ID within a node: [0..mpi_procs_per_node-1] (set during runtime)
 !Accelerators:
  !NVidia GPU:
        integer(C_INT):: gpus_found=0              !total number of NVidia GPUs found on the node (set during runtime)
        integer(C_INT):: gpu_count=0               !number of NVidia GPUs assigned to the current process (set during runtime)
        integer(C_INT):: gpu_start=0               !the number of the 1st NVidia GPU assigned to the current process: [gpu_start...gpu_start+gpu_count-1] (set during runtime)
        type(gpu_info_t),allocatable:: gpu_info(:) !information about available NVidia GPUs (set during runtime)
  !Intel MIC:
        integer(C_INT):: mics_found=0              !total number of Intel MICs found on the node (set during runtime)
        integer(C_INT):: mic_count=0               !number of Intel MICs assigned to the current process (set during runtime)
        integer(C_INT):: mic_start=0               !the number of the 1st Intel MIC assigned to the current process: [mic_start...mic_start+mic_count-1] (set during runtime)
        type(mic_info_t),allocatable:: mic_info(:) !information about available Intel MICs (set during runtime)
 !Process characteristics:
        real(8):: time_begin,time_end              !wall time for the MPI process
        integer:: exec_status=0                    !current execution status (0:success)
        character(MPI_MAX_PROCESSOR_NAME):: proc_name=' ' !processor name (set by runtime)
        integer:: proc_name_len                    !the length of the processor name
 !File Management:
        integer, private:: nof=0                          !current number of open files (local to each process)
        integer, private:: fhot(16:16+MAX_OPEN_FILES-1)=0 !file handle occupancy table (first 16 file handles [0..15] are reserved)
        integer, private:: ffhs(0:MAX_OPEN_FILES-1)=(/(j_,j_=16,16+MAX_OPEN_FILES-1)/) !a stack of free file handles
!Interfaces to Fortran wrappers to some MPI functions:
        interface
 !MPI_Get_address: Get the absolute MPI displacement of a local object (for remote accesses):
         subroutine MPI_Get_Displacement(location,disp,ierr) bind(c,name='MPI_Get_Displacement')
          import
          type(C_PTR), value, intent(in):: location !in: pointer to the local object
          integer(MPI_ADDRESS_KIND):: disp          !out: absolute MPI displacement
          integer(C_INT):: ierr                     !out: error code (0:success)
         end subroutine MPI_Get_Displacement

        end interface

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
        use STSUBS
        implicit none
        character(*), intent(in):: command
        integer, intent(inout):: ifh
        integer, intent(out):: ierr
        integer l,n

        ierr=0
        l=len_trim(command)
        if(l.gt.0) then
!GET FILE HANDLE:
         if(l.eq.3.and.(command(1:3).eq.'get'.or.command(1:3).eq.'GET')) then
          if(nof.lt.MAX_OPEN_FILES) then
           ifh=ffhs(nof); fhot(ifh)=1; nof=nof+1
          else
           call printl(jo,'#ERROR(service_mpi:file_handle): reached max amount of open files!')
           ierr=1
          endif
!FREE FILE HANDLE:
         elseif(l.eq.4.and.(command(1:4).eq.'free'.or.command(1:4).eq.'FREE')) then
          if(nof.gt.0) then
           if(ifh.ge.16.and.ifh.lt.16+MAX_OPEN_FILES) then
            if(fhot(ifh).ne.0) then
             fhot(ifh)=0; nof=nof-1; ffhs(nof)=ifh
            else
             call printl(jo,'#ERROR(service_mpi:file_handle): attempt to free an idle file handle!')
             ierr=4
            endif
           else
            call printl(jo,'#ERROR(service_mpi:file_handle): invalid file handle!')
            ierr=3
           endif
          else
           call printl(jo,'#ERROR(service_mpi:file_handle): attempt to free a file handle when no file handle is used!')
           ierr=2
          endif
!CLOSE ALL OPENED FILES:
         elseif(l.eq.4.and.(command(1:4).eq.'stop'.or.command(1:4).eq.'STOP')) then
          do n=16,16+MAX_OPEN_FILES-1
           if(nof.le.0) exit
           if(fhot(n).ne.0) then; close(n); fhot(n)=0; nof=nof-1; ffhs(nof)=n; endif
          enddo
         else
          call printl(jo,'#ERROR(service_mpi:file_handle): invalid command: '//command(1:l))
          ierr=5
         endif
        else
         call printl(jo,'#ERROR(service_mpi:file_handle): empty command passed!')
         ierr=6
        endif
        return
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
        integer, intent(in):: error_code     !in: error code
        character(*), intent(in):: error_msg !in: error message
        integer(INT_MPI):: i,errc
        integer:: l,ierr
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
        if(ierr.ne.0) write(*,'(''#WARNING: Process '',i8,'': Could not close all local files!'')') impir
        if(error_code.ne.0) then
         i=0; call MPI_ABORT(GLOBAL_MPI_COMM,i,errc)
         if(errc.ne.0) write(*,'(''#WARNING: Process '',i8,'': MPI ABORT CODE: '',i12)') impir,errc
        endif
        call MPI_FINALIZE(errc)
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
