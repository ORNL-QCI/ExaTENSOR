!TALSH::Fortran API testing.

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
        use, intrinsic:: ISO_C_BINDING
        implicit none
#ifndef NO_GPU
        interface
         subroutine test_nvtal_c(ierr) bind(c)
          import
          integer(C_INT), intent(out):: ierr
         end subroutine test_nvtal_c
        end interface
#endif
        integer(C_INT):: ierr

!Test NV-TAL C/C++ API interface:
#ifndef NO_GPU
!        write(*,'("Testing NV-TAL C/C++ API ...")')
!        call test_nvtal_c(ierr)
!        write(*,'("Done: Status ",i5)') ierr
!        if(ierr.ne.0) stop
!        write(*,*)''
#endif
!Test TAL-SH Fortran API interface:
        write(*,'("Testing TAL-SH Fortran API ...")')
        call test_talsh_f(ierr)
        write(*,'("Done: Status ",i5)') ierr
        if(ierr.ne.0) stop
        stop
        end program main
!------------------------------------
        subroutine test_talsh_f(ierr)
!Testing device-unified TAL-SH Fortran API.
        use, intrinsic:: ISO_C_BINDING
        use tensor_algebra
        use talsh
        implicit none
        integer(C_SIZE_T), parameter:: BUF_SIZE=8_8*1024_8*1024_8*1024_8 !desired Host argument buffer size in bytes
        integer(C_INT), parameter:: DIM_EXT=81 !tensor dimension extent
        integer(C_INT), parameter:: NUM_PASSES=2
        integer(C_INT), parameter:: MAX_TENSORS=24
        integer(C_INT), parameter:: MAX_TASKS=16
        integer(C_SIZE_T):: host_buf_size
        integer(C_INT):: i,j,m,n,ierr,num_gpus,host_arg_max
        integer(C_INT):: sts(MAX_TASKS)        !task statuses
        type(talsh_task_t):: tsks(MAX_TASKS)   !task handles
        type(talsh_tens_t):: tens(MAX_TENSORS) !tensors
        complex(8):: cval

        interface
         real(C_DOUBLE) function talshTensorImageNorm1_cpu(talsh_tens) bind(c,name='talshTensorImageNorm1_cpu')
          import
          implicit none
          type(talsh_tens_t), intent(in):: talsh_tens
         end function talshTensorImageNorm1_cpu
        end interface

        ierr=0
!Check GPU availability:
#ifndef NO_GPU
        write(*,'(1x,"Checking Nvidia GPU availability ... ")',ADVANCE='NO')
        ierr=cuda_get_device_count(num_gpus)
        write(*,'("Status ",i11,": Number of GPUs = ",i3)') ierr,num_gpus
        if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif
        !num_gpus=1
#else
        num_gpus=0
#endif
        if(MAX_TENSORS.lt.3*num_gpus) then; ierr=1; return; endif
!Initialize TALSH runtime:
        write(*,'(1x,"Initializing TALSH ... ")',ADVANCE='NO')
        host_buf_size=BUF_SIZE
#ifndef NO_GPU
        ierr=talsh_init(host_buf_size,host_arg_max,gpu_list=(/(i,i=0,num_gpus-1)/))
#else
        ierr=talsh_init(host_buf_size,host_arg_max)
#endif
        write(*,'("Status ",i11,": Size (Bytes) = ",i13,": Max args in HAB = ",i7)') ierr,host_buf_size,host_arg_max
        if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif

!Create tensors on Host and initialize them to a value:
        do i=1,MAX_TENSORS
         write(*,'(1x,"Constructing tensor block ",i2," ... ")',ADVANCE='NO') i
         select case(mod(i,3))
         case(1)
          cval=(0d0,0d0)
         case(2)
          cval=(1d-2,0d0)
         case(0)
          cval=(1d-3,0d0)
         end select
         ierr=talsh_tensor_construct(tens(i),R8,(/DIM_EXT,DIM_EXT,DIM_EXT,DIM_EXT/),init_val=cval)
         write(*,'("Status ",i11)') ierr; if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif
        enddo

!Tensor contractions on all GPUs:
        n=0 !number of tasks scheduled
        j=0
        do m=1,NUM_PASSES
         do i=0,num_gpus-1
          n=n+1
          write(*,'(1x,"Scheduling tensor contraction ",i2," on GPU ",i2,"... ")',ADVANCE='NO') n,i
          ierr=talsh_tensor_contract('D(a,b,i,j)+=L(j,c,k,a)*R(c,b,k,i)',tens(j+1),tens(j+2),tens(j+3),&
                                   &copy_ctrl=COPY_TTT,dev_id=talsh_flat_dev_id(DEV_NVIDIA_GPU,i),talsh_task=tsks(n))
          write(*,'("Status ",i11)') ierr; if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif
!         call talsh_task_print_info(tsks(n)) !debug
          j=j+3 !next set of tensor arguments
         enddo
        enddo
!Execute a tensor contraction on CPU (while the previous are running on GPUs):
!        n=n+1
!        write(*,'(1x,"Executing tensor contraction ",i2," on CPU ... ")',ADVANCE='NO') n
!        ierr=talsh_tensor_contract('D(a,b,i,j)+=L(j,c,k,a)*R(c,b,k,i)',tens(1),tens(2),tens(3),&
!                                  &dev_id=talsh_flat_dev_id(DEV_HOST,0),talsh_task=tsks(n))
!        write(*,'("Status ",i11)') ierr; if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif
!        call talsh_task_print_info(tsks(n)) !debug
!Synchronize and compare the results:
        write(*,'(1x,"Waiting upon completion of tensor contractions on all GPUs ... ")',ADVANCE='NO')
        ierr=talsh_tasks_wait(n,tsks,sts)
        write(*,'("Status ",i11," Completion =",8(1x,i8))') ierr,sts(1:n)
        if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif
!Printing results:
        do i=1,MAX_TENSORS,3
         call talsh_tensor_print_info(tens(i))
         print *,'TENSOR ',i,' NORM1 = ',talshTensorImageNorm1_cpu(tens(i))
        enddo

!Destruct TAL-SH task handles:
        do i=n,1,-1
         write(*,'(1x,"Destructing TAL-SH task handle ",i2," ... ")',ADVANCE='NO') i
         ierr=talsh_task_destruct(tsks(i))
         write(*,'("Status ",i11)') ierr; if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif
        enddo

!Destroy tensors:
        do i=MAX_TENSORS,1,-1
         write(*,'(1x,"Destructing tensor block ",i2," ... ")',ADVANCE='NO') i
         ierr=talsh_tensor_destruct(tens(i))
         write(*,'("Status ",i11)') ierr; if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif
        enddo
!Print run-time statistics:
        ierr=talsh_stats()
        if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif
!Shutdown TALSH:
        write(*,'(1x,"Shutting down TALSH ... ")',ADVANCE='NO')
        ierr=talsh_shutdown()
        write(*,'("Status ",i11)') ierr
        if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif
        return
        end subroutine test_talsh_f
