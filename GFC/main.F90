program main
 use dil_basic
 use stack_test
 use tree_test
 use dictionary_test
 implicit none
 real(8):: perf
 integer(INTD):: dev_out,ierr

 dev_out=6

 ierr=dil_test_tree(perf,dev_out)
 print *,'Tree testing status: ',ierr,': Performance: ',perf

 ierr=dil_test_stack(perf,dev_out)
 print *,'Stack testing status: ',ierr,': Performance: ',perf

 ierr=dil_test_dictionary(perf,dev_out)
 print *,'Dictionary testing status: ',ierr,': Performance: ',perf

 stop
end program main
