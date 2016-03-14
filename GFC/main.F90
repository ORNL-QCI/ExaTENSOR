!Copyright (C) 2014-2016 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2016 Oak Ridge National Laboratory (UT-Battelle)

!This file is part of ExaTensor.

!ExaTensor is free software: you can redistribute it and/or modify
!it under the terms of the GNU Lesser General Public License as published
!by the Free Software Foundation, either version 3 of the License, or
!(at your option) any later version.

!ExaTensor is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!GNU Lesser General Public License for more details.

!You should have received a copy of the GNU Lesser General Public License
!along with ExaTensor. If not, see <http://www.gnu.org/licenses/>.

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
