!TALSH::Fortran API testing.
        program main
        use, intrinsic:: ISO_C_BINDING
        use talsh
        implicit none

        interface
         subroutine test_talsh_c(ierr) bind(c)
          import
          integer(C_INT), intent(out):: ierr
         end subroutine test_talsh_c
        end interface

        integer(C_INT):: ierr

!Test C API interface:
        write(*,'("Testing TALSH C/C++ API ...")')
        call test_talsh_c(ierr)
        write(*,'("Done: Status ",i5)') ierr
!Test Fortran API interface:
        write(*,'("Testing TALSH Fortran API ...")')
!       call test_talsh_c(ierr)
        write(*,'("Done: Status ",i5)') ierr
        stop
        end program main
