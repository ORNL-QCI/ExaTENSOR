!PROGRAM: Q-FORCE: Massively-Parallel Quantum Many-Body Methodology on Heterogeneous HPC systems.
!AUTHOR: Dmitry I. Lyakh (Dmytro I. Liakh): quant4me@gmail.com
!REVISION: 2014/05/20
!COMPILATION:
! - Fortran 2003 at least.
! - MPI 2.0 at least.
! - OpenMP 4.0 at least.
! - CUDA 5.0 at least.
! - GNU compiling flags: -c -O3 --free-line-length-none -x f95-cpp-input -fopenmp
!   GNU linking flags: -lgomp
!   GNU BLAS/LAPACK: -lblas -llapack
! - Intel compiling flags: -c -O3 -fpp -vec-threshold4 -vec-report2 -openmp -openmp-report2 -DUSE_MKL
!   Intel linking flags: -Bdynamic
!   Intel MKL: -lmkl_core -lmkl_intel_thread -lmkl_intel_lp64 -lmkl_blas95_lp64 -lmkl_lapack95_lp64 -liomp5
!MPI launch notes:
! - MPI processes launched on the same node MUST have consecutive numbers!
! - Environment variable <QF_PROCS_PER_NODE> MUST be defined when using accelerators (number of processes per node)!
!FPP directives:
! - NO_PHI - do not use Intel Xeon Phi (MIC);
! - NO_GPU - do not use Nvidia GPU (CUDA);
! - NO_BLAS - BLAS/LAPACK calls will be replaced by my own routines (D.I.L.);
! - USE_MKL - use Intel MKL library for BLAS/LAPACK;
! - NO_OMP - do not use OpenMP (single-threaded processes);`currently will not work
! - NO_GNU - Fortran compiler is not GNU (affects Fortran timers);
!OUTPUT DEVICE:
! - jo (@service.mod) - generic output device handle;
!ENUMERATION OF DEVICES ON A NODE:
! - device_id = 0: Host (SMP CPU node, may be NUMA);
! - device_id = [1:MAX_GPUS_PER_NODE]: Nvidia GPUs (GPU#=device_id-1);
! - device_id = [MAX_GPUS_PER_NODE+1:MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE]: Intel MICs (MIC#=device_id-1-MAX_GPUS_PER_NODE);
! - device_id = [MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE+1:MAX_GPUS_PER_NODE+MAX_MICS_PER_NODE+MAX_AMDS_PER_NODE]: AMD GPUs
        program main !PARALLEL
        use qforce
        implicit none
        integer i,j,k,l,m,n,k0,ierr
        character(1024) str0
        integer, external:: omp_get_max_threads

        ierr=0
!Check basic data types:
        if(C_INT.ne.4.or.C_LONG_LONG.ne.8.or.C_FLOAT.ne.4.or.C_DOUBLE.ne.8) then
         write(*,'("#FATAL(main): C/Fortran basic data types are not interoperable!")')
         stop
        endif
!Initialization of MPI/OpenMP/CUDA:
#ifdef NO_OMP
        call MPI_INIT(ierr); if(ierr.ne.0) call quit(ierr,'#ERROR(main): MPI_INIT error!')
        mpi_thread_provided=0; max_threads=1
#else
        call MPI_INIT_THREAD(MPI_THREAD_FUNNELED,mpi_thread_provided,ierr); if(ierr.ne.0) call quit(ierr,'#ERROR(main): MPI_INIT_THREAD error!')
        if(mpi_thread_provided.gt.MPI_THREAD_SINGLE) then
         max_threads=omp_get_max_threads()
        else
         max_threads=1
        endif
#endif
        time_begin=MPI_WTIME() !walltime begin
        call MPI_COMM_SIZE(MPI_COMM_WORLD,impis,ierr); if(ierr.ne.0) call quit(ierr,'#ERROR(main): MPI_COMM_SIZE error!')
        call MPI_COMM_RANK(MPI_COMM_WORLD,impir,ierr); if(ierr.ne.0) call quit(ierr,'#ERROR(main): MPI_COMM_RANK error!')
        call file_handle('get',log_file,ierr); if(ierr.ne.0) call quit(ierr,'#ERROR(main): unable to get a log-file handle!')
!       write(*,*) impis,impir,max_threads,jo,log_file; call quit(-1,'Test') !debug
        call numchar(impir,k0,str0); open(log_file,file='qforce.'//str0(1:k0)//'.log',form='FORMATTED',status='UNKNOWN',err=2000) !open the log file for each process
        if(impir.ne.0) jo=log_file !redirect the standard output for slave processes to their log-files
        write(jo,'("   *** Q-FORCE v.14.05.13 by Dmitry I. Lyakh ***")')
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
        end select
#endif
        write(jo,'("Max number of threads/process      :       ",i4)') max_threads
#ifndef NO_BLAS
        write(jo,'("BLAS/LAPACK is enabled.")')
#else
        write(jo,'("BLAS/LAPACK is disabled.")')
#endif
        call MPI_GET_PROCESSOR_NAME(proc_name,proc_name_len,ierr); if(ierr.ne.0) call quit(ierr,'#ERROR(main): MPI_GET_PROCESSOR_NAME error!')
        call printl(jo,'My processor name                  : '//proc_name(1:proc_name_len))
        call get_environment(ierr); if(ierr.ne.0) then; ierr=0; write(jo,'("#ERROR(main): Unable to read environment variables! Ignored.")'); endif
        if(mpi_procs_per_node.gt.0) write(jo,'("Number of MPI processes per node   :       ",i4)') mpi_procs_per_node
        call gpu_nvidia_init(ierr); if(ierr.ne.0) then; ierr=0; write(jo,'("#ERROR(main): GPU(Nvidia) initialization failed! No GPU in use for this process! Ignored.")'); endif
        call MPI_BARRIER(MPI_COMM_WORLD,ierr); if(ierr.ne.0) call quit(ierr,'#ERROR(main): MPI_BARRIER error (trap 1)')

!Run the process life:
        call proceed(ierr); if(ierr.ne.0) then; exec_status=ierr; write(jo,'("#ERROR(main): Process life dirty!")'); endif
        call qforce_free_memory(ierr); if(exec_status.eq.0.and.ierr.ne.0) then; exec_status=ierr; write(jo,'("#ERROR(main): Unable to free QFORCE data structures!")'); endif
        if(exec_status.eq.0) then; call MPI_BARRIER(MPI_COMM_WORLD,ierr); if(ierr.ne.0) call quit(ierr,'#ERROR(main): MPI_BARRIER error (trap 2)'); endif
!Terminate the process:
999     call quit(exec_status,'Program terminated. Enjoy the results! Bye!')
!---------------------------------------------------------------
2000    call quit(-1,'#ERROR(main): unable to open the log-file!')

        contains

         subroutine get_environment(ier)
         integer, intent(inout):: ier
         integer j0
         character(64) qppn
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
           write(jo,'("#ERROR(main:get_environment): Environment variable QF_PROCS_PER_NODE is not a number!")')
           ier=1; return
          endif
         else
          write(jo,'("#ERROR(main:get_environment): Environment variable QF_PROCS_PER_NODE is not set!")')
          ier=2; return
         endif
         return
         end subroutine get_environment

         subroutine qforce_free_memory(ier) !This subroutine deallocates data structures used by QFORCE via its module.
         integer, intent(inout):: ier
         integer(C_INT) err_code
         ier=0
         if(allocated(gpu_info)) deallocate(gpu_info,STAT=ier); if(ier.ne.0) ier=1
         return
         end subroutine qforce_free_memory

        end program main
!-----------------------------------------------------------------------
        subroutine gpu_nvidia_init(ierr) !SERIAL
!This subroutine distributes GPUs available on a node among MPI processes running on that node.
!MPI processes running on the same node MUST have consecutive numbers.
!Some (or all) of the GPUs may fail to respond, leading to their ignorance.
!OUTPUT:
! - gpus_found - total number of GPUs found on a node;
! - gpu_start - first GPU ID assigned to the current MPI process;
! - gpu_count - number of consecutive GPU IDs assigned to the current MPI process;
! - ierr - error code (0:Success).
        use qforce
        implicit none
        integer, intent(inout):: ierr
        integer i,j,k,l,m,n
        integer(C_INT):: active_gpu,err_code

        ierr=0; gpus_found=0; gpu_count=0; gpu_start=0
#ifndef NO_GPU
 !Init GPU/GPUs:
        call cudaGetDeviceCount(gpus_found,err_code)
        if(err_code.ne.0) then; ierr=1; write(jo,'("#ERROR(gpu_nvidia_init): unable to get the GPU device count!")'); gpus_found=0; return; endif
        write(jo,'("Total amount of GPUs found on node : ",i10)') gpus_found
        if(gpus_found.gt.0) then
         if(allocated(gpu_info)) deallocate(gpu_info); allocate(gpu_info(0:gpus_found-1),STAT=ierr)
         if(ierr.ne.0) then; ierr=2; write(jo,'("#ERROR(gpu_nvidia_init): unable to allocate GPU_INFO!")'); gpu_count=0; return; endif
         do active_gpu=0,gpus_found-1
          call cudaSetDevice(active_gpu,err_code)
          if(err_code.eq.0) then
           write(jo,'("Info on GPU device                 : ",i10,":")') active_gpu
           call cudaGetDeviceProperties(active_gpu,gpu_info(active_gpu)%totalGlobalMem,gpu_info(active_gpu)%sharedMemPerBlock, &
                                        gpu_info(active_gpu)%regsPerBlock,gpu_info(active_gpu)%warpSize, &
                                        gpu_info(active_gpu)%maxThreadsPerBlock,gpu_info(active_gpu)%maxThreadsDim,gpu_info(active_gpu)%maxGridSize, &
                                        gpu_info(active_gpu)%clockRate,gpu_info(active_gpu)%totalConstMem, &
                                        gpu_info(active_gpu)%major,gpu_info(active_gpu)%minor, &
                                        gpu_info(active_gpu)%deviceOverlap,gpu_info(active_gpu)%multiProcessorCount, &
                                        gpu_info(active_gpu)%concurrentKernels,gpu_info(active_gpu)%ECCEnabled, &
                                        gpu_info(active_gpu)%asyncEngineCount,gpu_info(active_gpu)%memoryClockRate, &
                                        gpu_info(active_gpu)%memoryBusWidth,gpu_info(active_gpu)%maxThreadsPerMultiProcessor,err_code)
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
            write(jo,'("#WARNING(gpu_nvidia_init): unable to get properties for GPU #",i3)') active_gpu
           endif
          else
           write(jo,'("#WARNING(gpu_nvidia_init): unable to set GPU device #",i3)') active_gpu
          endif
         enddo
         call restrict_gpu_amount(ierr) !restrict usable GPUs range when multiple MPI processes run on a node
         if(ierr.ne.0) then; ierr=5; write(jo,'("#ERROR(gpu_nvidia_init): unable to restrict the amount of GPUs per MPI process!")'); gpu_count=0; gpu_start=0; endif
         write(jo,'("Range of GPUs assigned to process  : ",i4,"-",i4)') gpu_start,gpu_start+gpu_count-1
         write(jo,'("Ok")')
        elseif(gpus_found.lt.0) then !invalid (negative) number of GPUs found
         gpus_found=0; gpu_count=0; ierr=6
         write(jo,'("#ERROR(gpu_nvidia_init): negative number of GPUs found! Ignored.")')
        endif
#endif
        return

        contains

         subroutine restrict_gpu_amount(ierr) !SERIAL (Affects service::gpu_start, service::gpu_count)
!This subroutine distributes GPUs available on a node among the MPI processes residing on that node.
!It assumes that MPI processes are launched on each node consecutively!
         integer, intent(inout):: ierr
         integer j0
         ierr=0
         if(mpi_procs_per_node.gt.0) then
          j0=mod(gpus_found,mpi_procs_per_node)
          gpu_count=gpus_found/mpi_procs_per_node; gpu_start=mpi_proc_id_on_node*gpu_count
          if(mpi_proc_id_on_node.lt.j0) then
           gpu_start=gpu_start+mpi_proc_id_on_node; gpu_count=gpu_count+1
          else
           gpu_start=gpu_start+j0
          endif
         else
          ierr=1 !mpi_procs_per_node has not been specified
         endif
         return
         end subroutine restrict_gpu_amount

        end subroutine gpu_nvidia_init
