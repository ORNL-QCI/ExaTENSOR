        subroutine proceed(ierr) !PARALLEL
!This subroutine assigns a role to each MPI process:
! - MPI process 0: Global Root.
! - The subsequent range of MPI process numbers is divided into
!   Computing Groups (contiguous segments of process numbers),
!   each group comprising a Local Root (1st process in the segment)
!   and Computing Processes (all other processes in the segment).
!   The number of computing processes per local root is determined
!   by qforce::c_procs_per_local_root. Thus, each computing group initially
!   contains (1+c_procs_per_local_root) processes. This number, as well
!   as process roles may change during the run time for load balancing purposes.
        use qforce
        implicit none
        integer, intent(inout):: ierr
        integer i,j,k,l,m,n,k0,k1,k2,k3

        ierr=0; write(jo,'("###LAUNCHED PROCESS ",i9)') impir
        call mpi_barrier(MPI_COMM_WORLD,ierr); if(ierr.ne.0) call quit(ierr,'#ERROR(proceed): MPI_BARRIER error (1)!')

!DEBUG begin:
	my_role=c_process_private
	call c_proc_life(ierr); if(ierr.ne.0)  then; write(jo,'("#ERROR(proceed): C_PROCESS failed!")'); ierr=3; endif
	call mpi_barrier(MPI_COMM_WORLD,ierr); if(ierr.ne.0) call quit(ierr,'#ERROR(proceed): MPI_BARRIER error (2)!')
	return
!DEBUG end.

        if(impir.eq.0) then !MPI root --> Global Root
         my_role=global_root; my_group=-1
!        call global_root_life(ierr); if(ierr.ne.0) then; write(jo,'("#ERROR(proceed): GLOBAL ROOT failed!")'); ierr=1; endif
        else !MPI slave --> Local Root OR Computing Process (C-process)
         i=impir-1; my_group=i/(1+c_procs_per_local_root)
         if(mod(i,(1+c_procs_per_local_root)).eq.0) then
          my_role=local_root
!         call local_root_life(ierr); if(ierr.ne.0) then; write(jo,'("#ERROR(proceed): LOCAL ROOT failed!")'); ierr=2; endif
         else
          my_role=c_process_private
          call c_proc_life(ierr); if(ierr.ne.0) then; write(jo,'("#ERROR(proceed): C_PROCESS failed!")'); ierr=3; endif
         endif
        endif
        call mpi_barrier(MPI_COMM_WORLD,ierr); if(ierr.ne.0) call quit(ierr,'#ERROR(proceed): MPI_BARRIER error (2)!')
        return
        end subroutine proceed
