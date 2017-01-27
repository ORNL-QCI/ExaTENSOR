        program test
         use tensor_algebra
         use tensor_recursive_test
         implicit none
         integer(INTD):: ierr

         call test_tensor_recursive(ierr)
         print *,ierr
         stop
        end program test
