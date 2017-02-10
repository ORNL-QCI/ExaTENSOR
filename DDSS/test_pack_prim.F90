       program main
        use pack_prim_test
        use dil_basic, only: INTD,INTL
#ifdef USE_MPI_MOD
#ifdef FORTRAN2008
        use mpi_f08      !MPI Fortran 2008 interface `This will not work
#else
        use mpi          !MPI Fortran interface
#endif
        implicit none
#else
        implicit none
        include 'mpif.h' !MPI Fortran interface
#endif
        integer:: ierr,comm_size,my_rank
        integer(INTD):: errc

        call MPI_Init(ierr)
        if(ierr.eq.MPI_SUCCESS) then
         call MPI_Comm_Rank(MPI_COMM_WORLD,my_rank,ierr)
         call MPI_Comm_Size(MPI_COMM_WORLD,comm_size,ierr)
         if(my_rank.eq.0) write(*,'("MPI has been initialized: Communicator size = ",i11)') comm_size
!Packing primitives:
         if(my_rank.eq.0) write(*,'("Testing packing primitives:")')
         call MPI_Barrier(MPI_COMM_WORLD,ierr)
         ierr=test_pack_prim(errc)
         if(ierr.eq.0) then
          write(*,'(" Process ",i3," passed")') my_rank
         else
          write(*,'(" Process ",i3," failed: Error codes ",i11," -> ",i11)') my_rank,ierr,errc
          call MPI_Abort(MPI_COMM_WORLD,-1,ierr)
         endif
         call MPI_Finalize(ierr)
        else
         write(*,'("MPI initialization failed: Error ",i11)') ierr
        endif
        stop
       end program main
