        subroutine proceed(my_mpi_comm,ierr) !PARALLEL
!This subroutine assigns a role to each MPI process in <my_mpi_comm> communicator:
! * Role "Managing Process": Managing Processes (MP) form a tree with a single Level-0 Manager
!   as the Global Root (GR). Other Managers are indexed separately on each tree level.
!   Managers at the last tree level are called Local Roots (LR). Each Local Root
!   is assigned a group of Computing Processes (CP) and one or more Data Processes (DP).
! * Role "Computing Process": Computing Processes (CP) do the actual computations and
!   report to their respective Local Roots. Computing Processes are grouped into Computing
!   Groups (CG), each Computing Group being assigned a unique Local Root.
! * Role "Data Process": Data Processes (DP) facilitate data traffic and disk operations.
!This subroutine also creates local MPI communicators for computing groups and managers.
!       use g_process
!       use l_process
        use c_process
        use exatensor
        implicit none
        integer(INT_MPI), intent(in):: my_mpi_comm !in: global MPI communicator
        integer, intent(inout):: ierr              !out: error code (0:success)
        integer(INT_MPI):: errc
        integer(INTD):: i

        ierr=0
        GLOBAL_MPI_COMM=my_mpi_comm !this will be the global communication space for ExaTensor
!Preparing the global communication space:
        call MPI_COMM_SIZE(GLOBAL_MPI_COMM,impis,errc)
        if(errc.ne.0) then; ierr=1; call quit(ierr,'#ERROR(proceed): MPI_COMM_SIZE failed!'); endif
        call MPI_COMM_RANK(GLOBAL_MPI_COMM,impir,errc)
        if(errc.ne.0) then; ierr=2; call quit(ierr,'#ERROR(proceed): MPI_COMM_RANK failed!'); endif
        write(jo,'("###LAUNCHED PROCESS ",i9,": Syncing ... ")',ADVANCE='NO') impir
        call MPI_BARRIER(GLOBAL_MPI_COMM,errc)
        if(errc.eq.0) then
         write(jo,'("Ok")')
        else
         ierr=3; call quit(ierr,'#ERROR(proceed): MPI_BARRIER failed!')
        endif

!DEBUG begin:
	my_role=C_PROCESS_PRIVATE
	call c_proc_life(ierr); if(ierr.ne.0)  then; write(jo,'("#ERROR(proceed): C_PROCESS failed!")'); ierr=4; endif
	call MPI_BARRIER(GLOBAL_MPI_COMM,errc)
        if(errc.ne.0) then; ierr=5; call quit(ierr,'#ERROR(proceed): MPI_BARRIER error (1)!'); endif
	return
!DEBUG end.

        if(impir.eq.0) then !MPI root --> Global Root
         my_role=GLOBAL_ROOT; my_group=-1
!        call global_root_life(ierr); if(ierr.ne.0) then; write(jo,'("#ERROR(proceed): GLOBAL ROOT failed!")'); ierr=6; endif
        else !MPI slave --> Local Root OR Computing Process (C-process)
         i=impir-1; my_group=i/(1+C_PROCS_PER_LOCAL_ROOT); my_group_index=i-my_group*(1+C_PROCS_PER_LOCAL_ROOT)
         if(my_group_index.eq.0) then
          my_role=LOCAL_ROOT
!         call local_root_life(ierr); if(ierr.ne.0) then; write(jo,'("#ERROR(proceed): LOCAL ROOT failed!")'); ierr=7; endif
         else
          my_role=C_PROCESS_PRIVATE
          call c_proc_life(ierr); if(ierr.ne.0) then; write(jo,'("#ERROR(proceed): C_PROCESS failed!")'); ierr=8; endif
         endif
        endif
        call MPI_BARRIER(GLOBAL_MPI_COMM,errc)
        if(errc.ne.0) then; ierr=9; call quit(ierr,'#ERROR(proceed): MPI_BARRIER error (2)!'); endif
        return
        end subroutine proceed
