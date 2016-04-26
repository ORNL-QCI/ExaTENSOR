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
        if(ierr.ne.0) stop
        write(*,*)''
!Test Fortran API interface:
        write(*,'("Testing TALSH Fortran API ...")')
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
        integer(C_SIZE_T), parameter:: BUF_SIZE=1024*1024*1024 !desired Host argument buffer size in bytes
        integer(C_INT), parameter:: DIM_EXT=41 !tensor dimension extent
        integer(C_SIZE_T):: host_buf_size
        integer(C_INT):: i,ierr,host_arg_max,dev_gpu,dev_cpu
        type(talsh_tens_t):: tens(9) !three tensors for CPU, six for GPU
        type(talsh_task_t):: tsks(3) !three tasks (tensor contractions, three tensors per tensor contraction)

        ierr=0
!Init TALSH:
        write(*,'(1x,"Initializing TALSH ... ")',ADVANCE='NO')
        host_buf_size=BUF_SIZE
        ierr=talsh_init(host_buf_size,host_arg_max,gpu_list=(/0/))
        write(*,'("Status ",i11,": Size (Bytes) = ",i13,": Max Args in HAB = ",i7)') ierr,host_buf_size,host_arg_max
        if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif

!Create nine rank-4 tensors on Host and initialize them to value:
 !Tensor block 1:
        write(*,'(1x,"Constructing tensor block 1 ... ")',ADVANCE='NO')
        ierr=talsh_tensor_construct(tens(1),R8,(/DIM_EXT,DIM_EXT,DIM_EXT,DIM_EXT/),init_val=(0d0,0d0))
        write(*,'("Status ",i11)') ierr; if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif
 !Tensor block 2:
        write(*,'(1x,"Constructing tensor block 2 ... ")',ADVANCE='NO')
        ierr=talsh_tensor_construct(tens(2),R8,(/DIM_EXT,DIM_EXT,DIM_EXT,DIM_EXT/),init_val=(1d-2,0d0))
        write(*,'("Status ",i11)') ierr; if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif
 !Tensor block 3:
        write(*,'(1x,"Constructing tensor block 3 ... ")',ADVANCE='NO')
        ierr=talsh_tensor_construct(tens(3),R8,(/DIM_EXT,DIM_EXT,DIM_EXT,DIM_EXT/),init_val=(1d-3,0d0))
        write(*,'("Status ",i11)') ierr; if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif
 !Tensor block 4:
        write(*,'(1x,"Constructing tensor block 4 ... ")',ADVANCE='NO')
        ierr=talsh_tensor_construct(tens(4),R8,(/DIM_EXT,DIM_EXT,DIM_EXT,DIM_EXT/),init_val=(0d0,0d0))
        write(*,'("Status ",i11)') ierr; if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif
 !Tensor block 5:
        write(*,'(1x,"Constructing tensor block 5 ... ")',ADVANCE='NO')
        ierr=talsh_tensor_construct(tens(5),R8,(/DIM_EXT,DIM_EXT,DIM_EXT,DIM_EXT/),init_val=(1d-2,0d0))
        write(*,'("Status ",i11)') ierr; if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif
 !Tensor block 6:
        write(*,'(1x,"Constructing tensor block 6 ... ")',ADVANCE='NO')
        ierr=talsh_tensor_construct(tens(6),R8,(/DIM_EXT,DIM_EXT,DIM_EXT,DIM_EXT/),init_val=(1d-3,0d0))
        write(*,'("Status ",i11)') ierr; if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif
 !Tensor block 7:
        write(*,'(1x,"Constructing tensor block 7 ... ")',ADVANCE='NO')
        ierr=talsh_tensor_construct(tens(7),R8,(/DIM_EXT,DIM_EXT,DIM_EXT,DIM_EXT/),init_val=(0d0,0d0))
        write(*,'("Status ",i11)') ierr; if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif
 !Tensor block 8:
        write(*,'(1x,"Constructing tensor block 8 ... ")',ADVANCE='NO')
        ierr=talsh_tensor_construct(tens(8),R8,(/DIM_EXT,DIM_EXT,DIM_EXT,DIM_EXT/),init_val=(1d-2,0d0))
        write(*,'("Status ",i11)') ierr; if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif
 !Tensor block 9:
        write(*,'(1x,"Constructing tensor block 9 ... ")',ADVANCE='NO')
        ierr=talsh_tensor_construct(tens(9),R8,(/DIM_EXT,DIM_EXT,DIM_EXT,DIM_EXT/),init_val=(1d-3,0d0))
        write(*,'("Status ",i11)') ierr; if(ierr.ne.TALSH_SUCCESS) then; ierr=1; return; endif

!Choose execution devices:
        dev_cpu=talsh_flat_dev_id(DEV_HOST,0) !Host CPU
#ifndef NO_GPU
        dev_gpu=talsh_flat_dev_id(DEV_NVIDIA_GPU,0) !Nvidia GPU #0
#else
        dev_gpu=dev_cpu !fall back to CPU when no GPU in use
#endif
!Schedule two tensor contractions on GPU:

!Execute a tensor contraction on CPU (while the previous two are running on GPU):

!Synchronize and compare the results:


!Destroy tensors:
        do i=9,1,-1
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
