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
 use dictionary_test
 !use gfc_list_test
 use gfc_tree_test
 use gfc_dictionary_test
 implicit none
 real(8):: perf
 integer(INTD):: dev_out,ierr

 dev_out=6 !output device (defaults to screen)

!GFC containers:
 !List:
 !ierr=test_gfc_list(perf,dev_out)
 !if(ierr.eq.0) then
  !write(*,*) 'GFC::list testing status: ',ierr,'(PASSED): Performance: ',perf
 !else
  !write(*,*) 'GFC::list testing status: ',ierr,'(FAILED): Performance: ',perf
 !endif
 !Tree:
 ierr=test_gfc_tree(perf,dev_out)
 if(ierr.eq.0) then
  write(*,*) 'GFC::tree testing status: ',ierr,'(PASSED): Performance: ',perf
 else
  write(*,*) 'GFC::tree testing status: ',ierr,'(FAILED): Performance: ',perf
 endif
 !Dictionary:
 ierr=test_gfc_dictionary(perf,dev_out)
 if(ierr.eq.0) then
  write(*,*) 'GFC::dictionary testing status: ',ierr,'(PASSED): Performance: ',perf
 else
  write(*,*) 'GFC::dictionary testing status: ',ierr,'(FAILED): Performance: ',perf
 endif

!Legacy containers:
 !Stack:
 ierr=dil_test_stack(perf,dev_out)
 if(ierr.eq.0) then
  write(*,*) 'Legacy stack testing status: ',ierr,'(PASSED): Performance: ',perf
 else
  write(*,*) 'Legacy stack testing status: ',ierr,'(FAILED): Performance: ',perf
 endif
 !Dictionary:
 ierr=dil_test_dictionary(perf,dev_out)
 if(ierr.eq.0) then
  write(*,*) 'Legacy dictionary testing status: ',ierr,'(PASSED): Performance: ',perf
 else
  write(*,*) 'Legacy dictionary testing status: ',ierr,'(FAILED): Performance: ',perf
 endif

 stop
end program main
