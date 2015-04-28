!Distributed data storage primitives.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2015/04/27 (started 2015/03/18)
!All rights reserved!
        module distributed
!       use, intrinsic:: ISO_C_BINDING
        use service_mpi !includes ISO_C_BINDING
        implicit none
!PARAMETERS:
 !Output:
        integer, private:: CONS_OUT=6     !default output for this module
        logical, private:: VERBOSE=.true. !verbosity for errors
        logical, private:: DEBUG=.true.   !debugging mode
!TYPES:
 !Global data location descriptor:
        type, public:: DataLoc_t
         integer(INT_MPI), private:: PrRank=-1 !MPI rank on which the data resides
         integer(INT_MPI), private:: CommMPI   !MPI communicator the MPI rank refers to
         integer(INT_MPI), private:: Window    !MPI window the data is exposed with
         integer(INT_ADDR), private:: Offset   !Offset in the MPI window (in displacement units)
         integer(INT_ADDR), private:: Volume   !Number of elements (each element byte size = displacement unit)
         integer(INT_MPI), private:: DispUnit  !Displacement unit size in bytes
         contains
          procedure, public:: init=>DataLocInit
          procedure, public:: clean=>DataLocClean
        end type DataLoc_t
!DATA:

!FUNCTION VISIBILITY:
        private DataLocInit
        private DataLocClean

        contains
!METHODS:
!------------------------------------------------------------------------------
        subroutine DataLocInit(this,process_rank,mpi_window,offset,volume,ierr)
!DataLoc constructor.
        implicit none
        class(DataLoc_t), intent(out):: this           !out: DataLoc object
        integer(INT_MPI), intent(in):: process_rank    !in: process rank where the data resides
        integer(INT_MPI), intent(in):: mpi_window      !in: MPI window the data is exposed with
        integer(INT_ADDR), intent(in):: offset         !in: Offset in the MPI window (in displacement units)
        integer(INT_ADDR), intent(in):: volume         !in: Number of elements (each element byte size = displacement unit)
        integer(C_INT), intent(inout), optional:: ierr !out: error code (0:success)
        integer(C_INT):: errc

        errc=0
        if(process_rank.ge.0) then
         this%PrRank=process_rank
         this%CommMPI=MPI_COMM_WORLD
         this%Window=mpi_window
         if(offset.ge.0.and.volume.ge.0) then
          this%Offset=offset
          this%Volume=volume
         else
          errc=1
         endif
        else
         errc=2
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine DataLocInit
!-----------------------------------------
        subroutine DataLocClean(this,ierr)
!DataLoc cleaner.
        implicit none
        class(DataLoc_t), intent(out):: this           !out: DataLoc object
        integer(C_INT), intent(inout), optional:: ierr !out: error code (0:success)
        integer(C_INT):: errc

        errc=0
        this%PrRank=-1
        if(present(ierr)) ierr=errc
        return
        end subroutine DataLocClean

        end module distributed
