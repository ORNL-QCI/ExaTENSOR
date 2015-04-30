!Distributed data storage primitives.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2015/04/30 (started 2015/03/18)
!All rights reserved!
        module distributed
!       use, intrinsic:: ISO_C_BINDING
        use service_mpi !includes ISO_C_BINDING & MPI (must stay public)
        implicit none
        public !because of mpi.mod (mpif.h)
!PARAMETERS:
 !Output:
        integer, private:: CONS_OUT=6     !default output for this module
        logical, private:: VERBOSE=.true. !verbosity for errors
        logical, private:: DEBUG=.true.   !debugging mode
!TYPES:
 !Global data location descriptor:
        type, public:: DataDescr_t
         integer(INT_MPI), private:: PrRank=-1  !MPI rank on which the data resides
         integer(INT_MPI), private:: Window     !MPI window the data is exposed with
         integer(INT_MPI), private:: CommMPI    !MPI communicator the MPI rank refers to
         integer(INT_MPI), private:: DispUnit   !Displacement unit size in bytes
         integer(INT_ADDR), private:: Offset    !Offset in the MPI window (in displacement units)
         integer(INT_ADDR), private:: Volume    !Number of elements (each element size in bytes = displacement unit)
         contains
          procedure, public:: clean=>DataDescrClean      !clean a data descriptor
          procedure, public:: init=>DataDescrInit        !set up a data descriptor
!          procedure, public:: get_data=>DataDescrGetData !load data referred to by a data descriptor into a local buffer
!          procedure, public:: put_data=>DataDescrPutData !store data from a local buffer in the location specified by a data descriptor
        end type DataDescr_t
!DATA:

!FUNCTION VISIBILITY:
        private DataDescrClean
        private DataDescrInit

        contains
!METHODS:
!-------------------------------------------
        subroutine DataDescrClean(this,ierr)
!Data descriptor cleaner.
        implicit none
        class(DataDescr_t), intent(out):: this        !out: empty data descriptor
        integer(INTD), intent(inout), optional:: ierr !out: error code (0:success)
        integer(INTD):: errc

        errc=0
        this%PrRank=-1
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrClean
!--------------------------------------------------------------------------------------------------------
        subroutine DataDescrInit(this,process_rank,mpi_window,offset,volume,ierr,comm_mpi,check_thorough)
!Data descriptor constructor.
        implicit none
        class(DataDescr_t), intent(inout):: this             !inout: data descriptor
        integer(INT_MPI), intent(in):: process_rank          !in: process rank where the data resides
        integer(INT_MPI), intent(in):: mpi_window            !in: MPI window the data is exposed with
        integer(INT_ADDR), intent(in):: offset               !in: Offset in the MPI window (in displacement units)
        integer(INT_ADDR), intent(in):: volume               !in: Number of elements (each element byte size = displacement unit)
        integer(INTD), intent(inout), optional:: ierr        !out: error code (0:success)
        integer(INT_MPI), intent(in), optional:: comm_mpi    !in: MPI communicator (default: MPI_COMM_WORLD)
        logical(INTD), intent(in), optional:: check_thorough !in: if .true., enables the process order check (expensive): default=.false.
        integer(INT_MPI):: comm,my_group,win_group,res,errc
        integer(INT_ADDR):: attr
        logical(INTD):: flag

        errc=0
        if(process_rank.ge.0) then
         if(mpi_window.ne.MPI_WIN_NULL) then
          this%PrRank=process_rank
          this%Window=mpi_window
          if(present(comm_mpi)) then; comm=comm_mpi; else; comm=MPI_COMM_WORLD; endif
          call MPI_COMM_GROUP(comm,my_group,errc)
          if(errc.eq.0) then
           call MPI_WIN_GET_ATTR(mpi_window,MPI_WIN_DISP_UNIT,attr,flag,errc)
           if(errc.eq.0.and.flag) then
            this%DispUnit=int(attr,INT_MPI)
            call MPI_WIN_GET_GROUP(mpi_window,win_group,errc)
            if(errc.eq.0) then
             if(present(check_thorough)) then; flag=check_thorough; else; flag=.false.; endif
             res=MPI_UNEQUAL
             if(flag) then
              call MPI_WIN_GROUP_COMPARE(my_group,win_group,res,errc)
             else
              if(my_group.eq.win_group) res=MPI_IDENT
             endif
             if(errc.eq.0.and.res.eq.MPI_IDENT) then
              if(offset.ge.0.and.volume.ge.0) then
               this%Offset=offset
               this%Volume=volume
              else
               errc=1
              endif
             else
              errc=2
             endif
            else
             errc=3
            endif
           else
            errc=4
           endif
          else
           errc=5
          endif
         else
          errc=6
         endif
        else
         errc=7
        endif
        if(errc.ne.0) this%clean()
        if(present(ierr)) ierr=errc
        return
        end subroutine DataDescrInit
!-------------------------------------------------------------------

        end module distributed
