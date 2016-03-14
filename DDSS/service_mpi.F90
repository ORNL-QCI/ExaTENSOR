!This module provides general services for MPI parallel programs.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2016/03/14

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

       module service_mpi
        use, intrinsic:: ISO_C_BINDING
        use dil_basic
        !Depends on <mpi_fort.c>
#ifdef USE_MPI_MOD
#ifdef FORTRAN2008
        use mpi_f08      !MPI Fortran 2008 interface `This will not work
#else
        use mpi          !MPI Fortran interface
#endif
        implicit none
        public           !must be PUBLIC because of sharing the MPI header (safe)
#else
        implicit none
        public           !must be PUBLIC because of sharing the MPI header (safe)
        include 'mpif.h' !MPI Fortran interface
#endif
!Parameters:
 !Internal:
        logical, private:: DEBUG=.false.
 !MPI kinds:
        integer(C_INT), parameter, public:: INT_MPI=MPI_INTEGER_KIND   !default MPI integer kind
        integer(C_INT), parameter, public:: INT_ADDR=MPI_ADDRESS_KIND  !default MPI address/size kind
        integer(C_INT), parameter, public:: INT_OFFSET=MPI_OFFSET_KIND !default MPI offset kind
        integer(C_INT), parameter, public:: INT_COUNT=MPI_COUNT_KIND   !default MPI element count kind
 !File management:
        integer, parameter, private:: MAX_OPEN_FILES=1024-16 !max amount of open files per process (first 16 file handles [0..15] are reserved)
!Types:
 !NVidia GPU info:
        type, private:: gpu_info_t
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
        type, private:: mic_info_t
         integer(C_SIZE_T):: totalGlobalMem !in bytes
        end type mic_info_t
!Configuration variables:
 !Private:
        integer, private:: j_
 !General:
        integer, protected:: jo=6             !default output device (6 is the screen), modified during runtime
        integer, protected:: log_file         !log-file handle (set during runtime)
        logical, private:: process_up=.false. !set to .true. when the process is successfully initialized
 !MPI/OMP:
        logical, private:: comm_imported=.false. !tells whether the global communicator was imported or created
        integer(INT_MPI), protected:: GLOBAL_MPI_COMM=MPI_COMM_WORLD !global MPI communicator (defaults to MPI_COMM_WORLD)
        integer(INT_MPI), protected:: impis !global MPI communicator size (set by runtime)
        integer(INT_MPI), protected:: impir !global MPI rank of the process [0..impis-1] (set by runtime)
        integer(INT_MPI), protected:: mpi_thread_provided !level of multithreaded MPI support provided (set by runtime)
        integer(C_INT), protected:: max_threads=1         !max number of threads per MPI process (set during runtime)
        integer(C_INT), protected:: mpi_procs_per_node=0  !number of MPI processes per node (set by ENVIRONMENT), 0 means undefined
        integer(C_INT), protected:: mpi_proc_id_on_node=0 !internal process ID within a node: [0..mpi_procs_per_node-1] (set during runtime)
 !Accelerators:
  !NVidia GPU:
        integer(C_INT), protected:: gpus_found=0               !total number of NVidia GPUs found on the node (set during runtime)
        integer(C_INT), protected:: gpu_count=0                !number of NVidia GPUs assigned to the current process (set during runtime)
        integer(C_INT), protected:: gpu_start=0                !the number of the 1st NVidia GPU assigned to the current process: [gpu_start...gpu_start+gpu_count-1] (set during runtime)
        type(gpu_info_t), allocatable, protected:: gpu_info(:) !information about available NVidia GPUs (set during runtime)
  !Intel MIC:
        integer(C_INT), protected:: mics_found=0               !total number of Intel MICs found on the node (set during runtime)
        integer(C_INT), protected:: mic_count=0                !number of Intel MICs assigned to the current process (set during runtime)
        integer(C_INT), protected:: mic_start=0                !the number of the 1st Intel MIC assigned to the current process: [mic_start...mic_start+mic_count-1] (set during runtime)
        type(mic_info_t), allocatable, protected:: mic_info(:) !information about available Intel MICs (set during runtime)
 !Process characteristics:
        real(8), protected:: time_begin,time_end                     !wall time for the MPI process
        integer, protected:: exec_status=0                           !current execution status (0:success)
        character(MPI_MAX_PROCESSOR_NAME), protected:: proc_name=' ' !processor name (set by runtime)
        integer(INT_MPI), protected:: proc_name_len                  !the length of the processor name
 !File Management:
        integer, private:: nof=0                          !current number of open files (local to each process)
        integer, private:: fhot(16:16+MAX_OPEN_FILES-1)=0 !file handle occupancy table (first 16 file handles [0..15] are reserved)
        integer, private:: ffhs(0:MAX_OPEN_FILES-1)=(/(j_,j_=16,16+MAX_OPEN_FILES-1)/) !a stack of free file handles
!Interfaces to Fortran wrappers to some MPI functions:
        interface
 !MPI_Get_address: Get the absolute MPI displacement of a local object (for remote accesses):
         subroutine MPI_Get_Displacement(location,disp,ierr) bind(c,name='MPI_Get_Displacement')
          import
          implicit none
          type(C_PTR), value, intent(in):: location !in: pointer to the local object
          integer(MPI_ADDRESS_KIND):: disp          !out: absolute MPI displacement
          integer(C_INT):: ierr                     !out: error code (0:success)
         end subroutine MPI_Get_Displacement
        end interface
!Visibility:
        public dil_process_start
        public dil_process_finish
        public dil_process_global_id
        public dil_global_comm_size
        public dil_global_comm_barrier
        public MPI_Get_Displacement
        private gpu_nvidia_probe
!       private intel_mic_probe
        public file_handle
        public quit
        public get_memory_status

       contains
!--------------------------------------------------
        subroutine dil_process_start(ierr,ext_comm)
!Starts the (MPI) process. By default the usual MPI_COMM_WORLD communicator
!is initialized, unless an existing MPI communicator <ext_comm> is provided.
!In the latter case, the MPI threading level must be at least MPI_THREAD_FUNNELED.
        use stsubs, only: numchar,icharnum,printl
#ifndef NO_OMP
#ifdef USE_OMP_MOD
        use omp_lib, only: omp_get_max_threads
        implicit none
#else
        implicit none
        integer, external:: omp_get_max_threads
#endif
#else
        implicit none
#endif
        integer(INT_MPI), intent(out):: ierr              !out: error code (0:success)
        integer(INT_MPI), intent(in), optional:: ext_comm !in: external communicator
        integer(INT_MPI):: errc
        integer:: k0
        character(1024):: str0

        ierr=0; if(process_up) then; ierr=1; return; endif !process already initialized
        process_up=.false.; exec_status=0
!Start MPI:
#ifdef NO_OMP
        if(present(ext_comm)) then
         GLOBAL_MPI_COMM=ext_comm; comm_imported=.true.
         if(DEBUG) write(jo,'("#DEBUG(service_mpi::dil_process_start): Imported a single-thread MPI communicator: ",i13)')&
                   &GLOBAL_MPI_COMM
        else
         GLOBAL_MPI_COMM=MPI_COMM_WORLD; comm_imported=.false.
         call MPI_INIT(errc)
         if(errc.ne.0) then
          ierr=2
          call quit(ierr,'#ERROR(ExaTensor::service_mpi::dil_process_start): MPI_INIT error!',comm_imported)
          return
         else
          if(DEBUG) write(jo,'("#DEBUG(service_mpi::dil_process_start): Created a single-thread MPI communicator: ",i13)')&
                    &GLOBAL_MPI_COMM
         endif
        endif
        mpi_thread_provided=0; max_threads=1
#else
        if(present(ext_comm)) then
         GLOBAL_MPI_COMM=ext_comm; comm_imported=.true.
         mpi_thread_provided=MPI_THREAD_MULTIPLE !assumes at least this level of MPI threading
         max_threads=omp_get_max_threads()
         if(DEBUG) write(jo,'("#DEBUG(service_mpi::dil_process_start): Imported a multi-threaded MPI communicator: ",i13,1x,i5)')&
                   &GLOBAL_MPI_COMM,max_threads
        else
         GLOBAL_MPI_COMM=MPI_COMM_WORLD; comm_imported=.false.
         call MPI_INIT_THREAD(MPI_THREAD_MULTIPLE,mpi_thread_provided,errc)
         if(errc.ne.0) then
          ierr=3
          call quit(ierr,'#ERROR(ExaTensor::service_mpi::dil_process_start): MPI_INIT_THREAD error!',comm_imported)
          return
         else
          if(mpi_thread_provided.gt.MPI_THREAD_SINGLE) then
           max_threads=omp_get_max_threads()
          else
           max_threads=1
          endif
          if(DEBUG)&
          &write(jo,'("#DEBUG(service_mpi::dil_process_start): Created a multi-threaded MPI communicator: ",i13,1x,i5,1x,i9)')&
                    &GLOBAL_MPI_COMM,max_threads,mpi_thread_provided
         endif
        endif
#endif
        call MPI_BARRIER(GLOBAL_MPI_COMM,errc)
        if(errc.ne.0) then
         ierr=4
         call quit(ierr,'#ERROR(ExaTensor::service_mpi::dil_process_start): Initial MPI_BARRIER failure!',comm_imported)
         return
        endif
!Start time:
        time_begin=MPI_WTIME() !walltime begin
!Check the size of the basic data types:
        if(C_INT.ne.4.or.C_LONG_LONG.ne.8.or.C_FLOAT.ne.4.or.C_DOUBLE.ne.8) then
         ierr=5
         call quit(ierr,'#ERROR(ExaTensor::service_mpi::dil_process_start): C/Fortran basic data types are not interoperable!',&
              &comm_imported)
         return
        endif
!Get the global communicator characteristics:
        call MPI_COMM_SIZE(GLOBAL_MPI_COMM,impis,errc)
        if(errc.ne.0) then
         ierr=6
         call quit(ierr,'#ERROR(ExaTensor::service_mpi::dil_process_start): MPI_COMM_SIZE failed!',comm_imported)
         return
        endif
        call MPI_COMM_RANK(GLOBAL_MPI_COMM,impir,errc)
        if(errc.ne.0) then
         ierr=7
         call quit(ierr,'#ERROR(ExaTensor::service_mpi::dil_process_start): MPI_COMM_RANK failed!',comm_imported)
         return
        endif
        if(DEBUG) write(jo,'("#DEBUG(service_mpi::dil_process_start): My MPI rank is ",i10," out of ",i10)') impir,impis
!Open the log file:
        call file_handle('get',log_file,ierr)
        if(ierr.ne.0) then
         ierr=8
         call quit(ierr,'#ERROR(ExaTensor::service_mpi::dil_process_start): unable to get a log-file handle!',comm_imported)
         return
        endif
        call numchar(impir,k0,str0)
        open(log_file,file='qforce.'//str0(1:k0)//'.log',form='FORMATTED',status='UNKNOWN',err=2000) !open a log file for each process
        if(impir.ne.0) then
         jo=log_file !redirect the standard output for slave processes to their log files
        else
         if(DEBUG)&
         &write(jo,'("#DEBUG(service_mpi::dil_process_start): MPI slave output has been redirected to individual log files!")')
        endif
!Greetings:
        write(jo,'("   *** ExaTensor v.15.11.16 by Dmitry I. Lyakh ***")')
!Info:
        write(jo,'("MPI number of processes            : ",i10)') impis
        write(jo,'("Current process rank               : ",i10)') impir
#ifndef NO_OMP
        select case(mpi_thread_provided)
        case(MPI_THREAD_SINGLE)
         write(jo,'("MPI threading level                :          ",i1," (Single)")') mpi_thread_provided
        case(MPI_THREAD_FUNNELED)
         write(jo,'("MPI threading level                :          ",i1," (Funneled)")') mpi_thread_provided
        case(MPI_THREAD_SERIALIZED)
         write(jo,'("MPI threading level                :          ",i1," (Serialized)")') mpi_thread_provided
        case(MPI_THREAD_MULTIPLE)
         write(jo,'("MPI threading level                :          ",i1," (Multiple)")') mpi_thread_provided
        case default
         write(jo,'("MPI threading level                :          ",i1," (Unknown)")') mpi_thread_provided
         ierr=9
         call quit(ierr,'#ERROR(ExaTensor::service_mpi::dil_process_start): Unknown/invalid MPI threading support!',comm_imported)
         return
        end select
#endif
        write(jo,'("Max number of threads/process      :       ",i4)') max_threads
#ifndef NO_BLAS
        write(jo,'("BLAS/LAPACK is not disabled.")')
#else
        write(jo,'("BLAS/LAPACK is disabled.")')
#endif
        call MPI_GET_PROCESSOR_NAME(proc_name,proc_name_len,errc)
        if(errc.ne.0) then
         ierr=10
         call quit(ierr,'#ERROR(ExaTensor::service_mpi::dil_process_start): MPI_GET_PROCESSOR_NAME failed!',comm_imported)
         return
        endif
        call printl(jo,'My processor name                  : '//proc_name(1:proc_name_len))
!Get environment variables:
        call get_environment(ierr)
        if(ierr.ne.0) then
         ierr=0
         write(jo,'("#WARNING(ExaTensor::service_mpi::dil_process_start): Unable to read environment variables! Ignored.")')
        endif
        if(mpi_procs_per_node.gt.0) write(jo,'("Number of MPI processes per node   :       ",i4)') mpi_procs_per_node
!Probe Nvidia GPU(s):
        call gpu_nvidia_probe(ierr)
        if(ierr.ne.0) then
         ierr=0
         write(jo,'("#WARNING(ExaTensor::service_mpi::dil_process_start): GPU(NVidia) probe failed!",'//&
                   &'" No GPU in use for this process! Ignored.")')
        endif
        call MPI_BARRIER(GLOBAL_MPI_COMM,errc)
        if(errc.ne.0) then
         ierr=11
         call quit(ierr,'#ERROR(ExaTensor::service_mpi::dil_process_start): Intermediate MPI_BARRIER failure!',comm_imported)
         return
        endif
!Clear the process status if no error:
        exec_status=ierr
        if(ierr.eq.0) process_up=.true.
        return
!------------------------------------------------------------------------------------------------------------------
2000    call quit(-1,'#ERROR(ExaTensor::service_mpi::dil_process_start): unable to open a log file!',comm_imported)
        ierr=12
        return

        contains

         subroutine get_environment(ier)
         integer, intent(inout):: ier
         integer:: j0
         character(64):: qppn
         ier=0; mpi_procs_per_node=0; mpi_proc_id_on_node=0
!QF_PROCS_PER_NODE:
         qppn=' '; call get_environment_variable("QF_PROCS_PER_NODE",qppn)
         j0=len_trim(qppn)
         if(j0.gt.0) then
          mpi_procs_per_node=icharnum(j0,qppn(1:j0))
          if(j0.gt.0) then
           mpi_proc_id_on_node=mod(impir,mpi_procs_per_node)
          else
           mpi_procs_per_node=0
           write(jo,'("#ERROR(service_mpi::dil_process_start:get_environment): ",'//&
                   &'"Environment variable QF_PROCS_PER_NODE is not a number!")')
           ier=1; return
          endif
         else
          write(jo,'("#ERROR(service_mpi::dil_process_start:get_environment): ",'//&
                  &'"Environment variable QF_PROCS_PER_NODE is not set!")')
          ier=2; return
         endif
         return
         end subroutine get_environment

        end subroutine dil_process_start
!------------------------------------------
        subroutine dil_process_finish(ierr)
!Terminates the (MPI) process after a successful execution.
!Module variable <exec_status> determines whether the execution was successful.
        implicit none
        integer(INT_MPI), intent(out):: ierr !out: error code (0:success)
        integer(INT_MPI):: errc
        integer:: erc

        ierr=0
        if(process_up) then
         call free_info_mem(ierr)
         if(ierr.ne.0) write(jo,'("#WARNING(ExaTensor::service_mpi::dil_process_finish): Unable to free info data!")')
         process_up=.false.
         if(exec_status.eq.0) then
          call MPI_BARRIER(GLOBAL_MPI_COMM,errc)
          if(errc.ne.0) then
           ierr=ierr+3
           call quit(ierr,'#ERROR(ExaTensor::service_mpi::dil_process_finish): Final MPI_BARRIER failure!',comm_imported)
          endif
         endif
         call quit(exec_status,'###ExaTensor exited:',comm_imported)
        else
         ierr=ierr+7
        endif
        return

        contains

         subroutine free_info_mem(ier) !This subroutine deallocates data structures used by QFORCE via its module.
          integer(INT_MPI), intent(inout):: ier
          ier=0
          if(allocated(gpu_info)) deallocate(gpu_info,STAT=ier); if(ier.ne.0) ier=ier+1
          if(allocated(mic_info)) deallocate(mic_info,STAT=ier); if(ier.ne.0) ier=ier+2
          return
         end subroutine free_info_mem

        end subroutine dil_process_finish
!-------------------------------------------------------
        function dil_process_global_id(ierr) result(gid)
!Returns the global (MPI) process rank.
        implicit none
        integer(INT_MPI):: gid                         !out: global rank of the process
        integer(INT_MPI), intent(out), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0; gid=-1
        if(process_up) then
         call MPI_COMM_RANK(GLOBAL_MPI_COMM,gid,errc)
        else
         errc=-1
        endif
        if(present(ierr)) ierr=errc
        return
        end function dil_process_global_id
!------------------------------------------------------
        function dil_global_comm_size(ierr) result(gcs)
!Returns the size of the global (MPI) communicator.
        implicit none
        integer(INT_MPI):: gcs                         !out: global communicator size
        integer(INT_MPI), intent(out), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        errc=0; gcs=-1
        if(process_up) then
         call MPI_COMM_SIZE(GLOBAL_MPI_COMM,gcs,errc)
        else
         errc=-1
        endif
        if(present(ierr)) ierr=errc
        return
        end function dil_global_comm_size
!-----------------------------------------------
        subroutine dil_global_comm_barrier(ierr)
!A barrier over the global (MPI) communicator.
        implicit none
        integer(INT_MPI), intent(out), optional:: ierr !out: error code (0:success)
        integer(INT_MPI):: errc

        if(process_up) then
         call MPI_BARRIER(GLOBAL_MPI_COMM,errc)
        else
         errc=-1
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine dil_global_comm_barrier
!----------------------------------------
        subroutine gpu_nvidia_probe(ierr) !SERIAL
!This subroutine distributes GPUs available on a node among the MPI processes running
!on that node. MPI processes running on the same node MUST have consecutive numbers.
!Some (or all) of the GPUs may fail to respond, leading to their ignorance.
!OUTPUT:
! - gpus_found - total number of GPUs found on a node;
! - gpu_start - first GPU ID assigned to the current MPI process;
! - gpu_count - number of consecutive GPU IDs assigned to the current MPI process;
! - ierr - error code (0:Success).
        use extern_names
        implicit none
        integer, intent(inout):: ierr
        integer i,j,k,l,m,n
        integer(C_INT):: active_gpu,err_code

        ierr=0; gpus_found=0; gpu_count=0; gpu_start=0
#ifndef NO_GPU
 !Init GPU/GPUs:
        call cudaGetDeviceCount(gpus_found,err_code)
        if(err_code.ne.0) then
         ierr=1; write(jo,'("#ERROR(gpu_nvidia_probe): Unable to get the GPU device count!")')
         gpus_found=0; return
        endif
        write(jo,'("Total amount of GPUs found on node : ",i10)') gpus_found
        if(gpus_found.gt.0) then
         if(allocated(gpu_info)) deallocate(gpu_info); allocate(gpu_info(0:gpus_found-1),STAT=ierr)
         if(ierr.ne.0) then
          ierr=2; write(jo,'("#ERROR(gpu_nvidia_probe): Unable to allocate GPU_INFO!")')
          gpu_count=0; return
         endif
         do active_gpu=0,gpus_found-1
          call cudaSetDevice(active_gpu,err_code)
          if(err_code.eq.0) then
           write(jo,'("Info on GPU device                 : ",i10,":")') active_gpu
           call cudaGetDeviceProperties(active_gpu,gpu_info(active_gpu)%totalGlobalMem,gpu_info(active_gpu)%sharedMemPerBlock,&
                                       &gpu_info(active_gpu)%regsPerBlock,gpu_info(active_gpu)%warpSize,&
                                       &gpu_info(active_gpu)%maxThreadsPerBlock,gpu_info(active_gpu)%maxThreadsDim,&
                                       &gpu_info(active_gpu)%maxGridSize,gpu_info(active_gpu)%clockRate,&
                                       &gpu_info(active_gpu)%totalConstMem,gpu_info(active_gpu)%major,gpu_info(active_gpu)%minor,&
                                       &gpu_info(active_gpu)%deviceOverlap,gpu_info(active_gpu)%multiProcessorCount,&
                                       &gpu_info(active_gpu)%concurrentKernels,gpu_info(active_gpu)%ECCEnabled,&
                                       &gpu_info(active_gpu)%asyncEngineCount,gpu_info(active_gpu)%memoryClockRate,&
                                       &gpu_info(active_gpu)%memoryBusWidth,gpu_info(active_gpu)%maxThreadsPerMultiProcessor,&
                                       &err_code)
           if(err_code.eq.0) then
            write(jo,'(1x,"GPU major revision number         : ",i10)') gpu_info(active_gpu)%major
            write(jo,'(1x,"GPU minor revision number         : ",i10)') gpu_info(active_gpu)%minor
            write(jo,'(1x,"GPU total global memory (B)       : ",i10)') gpu_info(active_gpu)%totalGlobalMem
            write(jo,'(1x,"GPU shared memory per SM (B)      : ",i10)') gpu_info(active_gpu)%sharedMemPerBlock
            write(jo,'(1x,"GPU constant memory (B)           : ",i10)') gpu_info(active_gpu)%totalConstMem
            write(jo,'(1x,"GPU registers per block           : ",i10)') gpu_info(active_gpu)%regsPerBlock
            write(jo,'(1x,"GPU warp size                     : ",i10)') gpu_info(active_gpu)%warpSize
            write(jo,'(1x,"GPU max threads per block         : ",i10)') gpu_info(active_gpu)%maxThreadsPerBlock
            write(jo,'(1x,"GPU max threads along each dim    : ",i10,1x,i10,1x,i10)') gpu_info(active_gpu)%maxThreadsDim(1:3)
            write(jo,'(1x,"GPU max blocks along each dim     : ",i10,1x,i10,1x,i10)') gpu_info(active_gpu)%maxGridSize(1:3)
            write(jo,'(1x,"GPU amount of multiprocessors     : ",i10)') gpu_info(active_gpu)%multiProcessorCount
            write(jo,'(1x,"GPU max threads per multiprocessor: ",i10)') gpu_info(active_gpu)%maxThreadsPerMultiProcessor
            write(jo,'(1x,"GPU copy/computation overlap      : ",i10)') gpu_info(active_gpu)%deviceOverlap
            write(jo,'(1x,"GPU number of transfer engines    : ",i10)') gpu_info(active_gpu)%asyncEngineCount
            write(jo,'(1x,"GPU kernel concurrency            : ",i10)') gpu_info(active_gpu)%concurrentKernels
            write(jo,'(1x,"GPU ECC status                    : ",i10)') gpu_info(active_gpu)%ECCEnabled
            write(jo,'(1x,"GPU clock rate (KHz)              : ",i10)') gpu_info(active_gpu)%clockRate
            write(jo,'(1x,"GPU memory clock rate (KHz)       : ",i10)') gpu_info(active_gpu)%memoryClockRate
            write(jo,'(1x,"GPU memory bus width (b)          : ",i10)') gpu_info(active_gpu)%memoryBusWidth
           else
            write(jo,'("#WARNING(gpu_nvidia_probe): Unable to get properties for GPU #",i3)') active_gpu
           endif
          else
           write(jo,'("#WARNING(gpu_nvidia_probe): Unable to set GPU device #",i3)') active_gpu
          endif
         enddo
         call restrict_gpu_amount(ierr) !restrict usable GPUs range when multiple MPI processes run on a node
         if(ierr.ne.0) then
          ierr=3; write(jo,'("#ERROR(gpu_nvidia_probe): Unable to restrict the amount of GPUs per MPI process!")')
          gpu_count=0; gpu_start=0
         endif
         write(jo,'("Range of GPUs assigned to process  : ",i4," -",i4)') gpu_start,gpu_start+gpu_count-1
         write(jo,'("Ok")')
        elseif(gpus_found.lt.0) then !invalid (negative) number of GPUs found
         gpus_found=0; gpu_count=0; ierr=4
         write(jo,'("#ERROR(gpu_nvidia_probe): Negative number of GPUs found! Ignored.")')
        endif
#endif
        return

        contains

         subroutine restrict_gpu_amount(ier) !SERIAL (Affects service::gpu_start, service::gpu_count)
!This subroutine distributes GPUs available on a node among the MPI processes residing on that node.
!It assumes that MPI processes are launched on each node consecutively!
         integer, intent(inout):: ier
         integer j0
         ier=0
         if(mpi_procs_per_node.gt.0) then
          j0=mod(gpus_found,mpi_procs_per_node)
          gpu_count=gpus_found/mpi_procs_per_node; gpu_start=mpi_proc_id_on_node*gpu_count
          if(mpi_proc_id_on_node.lt.j0) then
           gpu_start=gpu_start+mpi_proc_id_on_node; gpu_count=gpu_count+1
          else
           gpu_start=gpu_start+j0
          endif
         else
          ier=1 !mpi_procs_per_node has not been specified
         endif
         return
         end subroutine restrict_gpu_amount

        end subroutine gpu_nvidia_probe
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
        use stsubs, only: printl
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
!-----------------------------------------------------
        subroutine quit(error_code,error_msg,no_final)
!This subroutine prints the error message and safely terminates the parallel code execution.
!INPUT:
! - error_code - error code (-1 is the default non-specific error code);
! - error_msg - error message;
! - no_final - if .true., no MPI_Finalize and stop, but return;
!OUTPUT:
! - error information printed by EACH process.
!PARALLEL: YES.
        use stsubs, only: printl
        implicit none
        integer, intent(in):: error_code         !in: error code
        character(*), intent(in):: error_msg     !in: error message
        logical, intent(in), optional:: no_final !selects between Finalizing and just Returning
        integer(INT_MPI):: i,errc
        integer:: l,ierr
        logical:: nf
!Time:
        time_end=MPI_WTIME()
        if(present(no_final)) then; nf=no_final; else; nf=.false.; endif
!Error message:
        l=len_trim(error_msg); if(l.gt.0.and.(error_code.ne.0.or.impir.eq.0)) call printl(jo,error_msg(1:l))
        if(error_code.eq.0) then
         write(jo,'(''###ExaTensor Process '',i9,'': Success: Wall Time (sec): '',f12.2)') impir,time_end-time_begin
        else
         write(jo,'(''###ExaTensor Process '',i9,'': Error '',i12,'': Wall Time (sec): '',f12.2)')&
          &impir,error_code,time_end-time_begin
        endif
        flush(jo)
!Close all local files opened by the process:
        l=0; call file_handle('stop',l,ierr)
        if(ierr.ne.0) write(*,'(''#WARNING(ExaTensor::quit): Process '',i9,'': Could not close all local files!'')') impir
        if(error_code.ne.0.and.(.not.nf)) then
         i=0; call MPI_ABORT(GLOBAL_MPI_COMM,i,errc)
         if(errc.ne.0) write(*,'(''#WARNING(ExaTensor::quit): Process '',i9,'': MPI ABORT CODE: '',i12)') impir,errc
        endif
        if(.not.nf) call MPI_FINALIZE(errc)
        return
        end subroutine quit
!----------------------------------------------------------------------
        subroutine get_memory_status(total_ram,free_ram,used_swap,ierr)
!This subroutine returns:
! - total_ram - total usable RAM available on the node in bytes;
! - free_ram - free usable RAM available on the node in bytes;
! - used_swap - current swap size in bytes;
        use extern_names, only: get_memory_stat
        implicit none
        integer(C_SIZE_T), intent(out):: total_ram,free_ram,used_swap
        integer(C_INT), intent(inout):: ierr
        ierr=get_memory_stat(total_ram,free_ram,used_swap)
        return
        end subroutine get_memory_status

       end module service_mpi
