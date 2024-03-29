!This module contains interfaces to auxiliary C/CUDA functions callable from Fortran.

!Copyright (C) 2014-2022 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2022 Oak Ridge National Laboratory (UT-Battelle)

!LICENSE: BSD 3-Clause

       module extern_names
        use, intrinsic:: ISO_C_BINDING

        interface
!Auxiliary C functions:
 !Obtain an offset byte pointer:
         type(C_PTR) function ptr_offset(byte_ptr,byte_offset) bind(c,name='ptr_offset')
          import
          implicit none
          type(C_PTR), value, intent(in):: byte_ptr
          integer(C_SIZE_T), value, intent(in):: byte_offset
         end function ptr_offset
 !Converts a C pointer into an integer (address):
         function c_ptr_value(cptr) result(cpval) bind(c,name='c_ptr_value')
          import
          implicit none
          type(C_PTR), value:: cptr
          integer(C_SIZE_T):: cpval
         end function c_ptr_value
 !Converts an integer address into a C pointer:
         subroutine c_ptr_set(cpval,cptr) bind(c,name='c_ptr_set')
          import
          implicit none
          integer(C_SIZE_T), value:: cpval
          type(C_PTR):: cptr
         end subroutine c_ptr_set
 !Print a C pointer from Fortran:
         subroutine print_c_ptr(c_pointer) bind(c,name='print_c_ptr')
          import
          implicit none
          type(C_PTR), value:: c_pointer
         end subroutine print_c_ptr
 !Get the Linux memory status:
         integer(C_INT) function get_memory_stat(total_ram,free_ram,used_swap) bind(C)
          import
          integer(C_SIZE_T), intent(out):: total_ram
          integer(C_SIZE_T), intent(out):: free_ram
          integer(C_SIZE_T), intent(out):: used_swap
         end function get_memory_stat
 !Get an accurate C time:
         function accu_time() result(tm) bind(C,name='accu_time')
          import
          real(C_DOUBLE):: tm
         end function accu_time
!C wrappers for CUDA runtime functions:
#ifndef NO_GPU
 !Get the total number of available GPU(Nvidia) devices on a node:
         subroutine cudaGetDeviceCount(count,err_code) bind(c,name='cudagetdevicecount')
          import
          implicit none
          integer(C_INT), intent(out):: count
          integer(C_INT):: err_code
         end subroutine cudaGetDeviceCount
 !Set the active GPU(Nvidia) device:
         subroutine cudaSetDevice(device,err_code) bind(c,name='cudasetdevice')
          import
          implicit none
          integer(C_INT), value, intent(in):: device
          integer(C_INT):: err_code
         end subroutine cudaSetDevice
 !Get GPU(Nvidia) device properties:
         subroutine cudaGetDeviceProperties(device,totalGlobalMem_,sharedMemPerBlock_,regsPerBlock_,warpSize_,&
                                            &maxThreadsPerBlock_,maxThreadsDim_,maxGridSize_,clockRate_,totalConstMem_,&
                                            &major_,minor_,deviceOverlap_,multiProcessorCount_,concurrentKernels_,&
                                            &ECCEnabled_,asyncEngineCount_,memoryClockRate_,memoryBusWidth_,&
                                            &maxThreadsPerMultiProcessor_,err_code) bind(c,name='cudagetdeviceproperties')
          import
          implicit none
          integer(C_INT), value, intent(in):: device
          integer(C_SIZE_T), intent(out):: totalGlobalMem_
          integer(C_SIZE_T), intent(out):: sharedMemPerBlock_
          integer(C_INT), intent(out):: regsPerBlock_
          integer(C_INT), intent(out):: warpSize_
          integer(C_INT), intent(out):: maxThreadsPerBlock_
          integer(C_INT), intent(out):: maxThreadsDim_(3)
          integer(C_INT), intent(out):: maxGridSize_(3)
          integer(C_INT), intent(out):: clockRate_
          integer(C_SIZE_T), intent(out):: totalConstMem_
          integer(C_INT), intent(out):: major_
          integer(C_INT), intent(out):: minor_
          integer(C_INT), intent(out):: deviceOverlap_
          integer(C_INT), intent(out):: multiProcessorCount_
          integer(C_INT), intent(out):: concurrentKernels_
          integer(C_INT), intent(out):: ECCEnabled_
          integer(C_INT), intent(out):: asyncEngineCount_
          integer(C_INT), intent(out):: memoryClockRate_
          integer(C_INT), intent(out):: memoryBusWidth_
          integer(C_INT), intent(out):: maxThreadsPerMultiProcessor_
          integer(C_INT):: err_code
         end subroutine cudaGetDeviceProperties
 !GPU(Nvidia) device synchronization:
         subroutine cudaDeviceSynchronize(err_code) bind(c,name='cudadevicesynchronize')
          import
          implicit none
          integer(C_INT):: err_code
         end subroutine cudaDeviceSynchronize
#endif
        end interface

       end module extern_names
