!TALSH::Fortran API testing.
        program test_talsh_f
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
        call test_talsh_c(ierr)
!Test Fortran API interface:
!       ...
        stop
        end program test_talsh_f
