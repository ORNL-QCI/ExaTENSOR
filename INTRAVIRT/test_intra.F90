!INTRAVIRT tester
        program test_intra
         use tensor_algebra
         use tensor_recursive_test
         use parse_prim_test
         implicit none
         integer(INTD):: ierr

!Module tensor_recursive:
         call test_tensor_recursive(ierr)
!Module parse_prim:
         ierr=test_parse_prim()
         stop
        end program test_intra
