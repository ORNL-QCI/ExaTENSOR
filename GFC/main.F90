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
 use gfc_vector_test
 !use gfc_list_test
 use gfc_tree_test
 use gfc_dictionary_test
 implicit none
 real(8):: perf
 integer(INTD):: dev_out,ierr
 real(8), external:: dil_test_infer_overhead

 dev_out=6 !output device (defaults to screen)

!GFC containers:
 !Vector:
 ierr=test_gfc_vector(perf,dev_out)
 if(ierr.eq.0) then
  write(*,*) 'gfc::vector testing status: ',ierr,'(PASSED): Performance: ',perf
 else
  write(*,*) 'gfc::vector testing status: ',ierr,'(FAILED): Performance: ',perf
 endif
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
  write(*,*) 'gfc::tree testing status: ',ierr,'(PASSED): Performance: ',perf
 else
  write(*,*) 'gfc::tree testing status: ',ierr,'(FAILED): Performance: ',perf
 endif
 !Dictionary:
 ierr=test_gfc_dictionary(perf,dev_out)
 if(ierr.eq.0) then
  write(*,*) 'gfc::dictionary testing status: ',ierr,'(PASSED): Performance: ',perf
 else
  write(*,*) 'gfc::dictionary testing status: ',ierr,'(FAILED): Performance: ',perf
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

!Dynamic type inferrence overhead:
 write(*,'(1x,"Dynamic type inferrence slowdown: ")',ADVANCE='NO')
 perf=dil_test_infer_overhead(2**23)
 write(*,*) perf

 stop
 contains

end program main

subroutine my_add(a,b,c)
 class(*), intent(in):: a
 class(*), intent(in):: b
 class(*), intent(inout):: c
 real(8), pointer:: ap,bp,cp

 select type(a); type is(real(8)); ap=>a; end select
 select type(b); type is(real(8)); bp=>b; end select
 select type(c); type is(real(8)); cp=>c; end select
 cp=cp+ap+bp
 return
end subroutine my_add

function dil_test_infer_overhead(n) result(slowdown)
 real(8):: slowdown
 integer, intent(in):: n
 integer:: i,j
 real(8), allocatable:: a(:),b(:)
 real(8):: c,tm,tm1,tm2

 interface
  subroutine my_add(a,b,c)
   class(*), intent(in):: a
   class(*), intent(in):: b
   class(*), intent(inout):: c
  end subroutine my_add
 end interface

 allocate(a(n),b(n))
 call random_number(a)
 call random_number(b)
!Direct:
 call cpu_time(tm)
 c=0d0
 do j=1,8
  do i=1,n
   c=c+a(i)+b(i)
  enddo
 enddo
 call cpu_time(tm1); tm1=tm1-tm
 write(*,'(F8.4,1x,F20.7,1x)',ADVANCE='NO') tm1,c
!Indirect:
 call cpu_time(tm)
 c=0d0
 do j=1,8
  do i=1,n
   call my_add(a(i),b(i),c)
  enddo
 enddo
 call cpu_time(tm2); tm2=tm2-tm
 write(*,'(F8.4,1x,F20.7,1x)',ADVANCE='NO') tm2,c
 slowdown=tm2/tm1
 deallocate(a,b)
 return
end function dil_test_infer_overhead
