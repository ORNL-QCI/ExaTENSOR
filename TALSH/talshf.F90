!ExaTensor::TAL-SH: Device-unified user-level API:
!REVISION: 2016/04/19

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
!--------------------------------------------------------------------------------
       module talsh
        use tensor_algebra_cpu_phi !device-specific tensor algebra API + basic
        implicit none
        private
!EXTERNAL PUBLIC:
        public tensor_shape_t    !tensor shape (Fortran)
        public tensor_block_t    !tensor block (Fortran)
        public MAX_SHAPE_STR_LEN !max length of a shape-defining string
!PARAMETERS:
 !Generic:
        integer(INTD), private:: CONS_OUT=6 !default output device for this module
        integer(INTD), private:: DEBUG=0    !debugging mode for this module
        logical, private:: VERBOSE=.true.   !verbosity for errors
 !Errors (keep consistent with "talsh.h"):
        integer(C_INT), parameter, public:: TALSH_SUCCESS=0                   !success
        integer(C_INT), parameter, public:: TALSH_FAILURE=-666                !generic failure
        integer(C_INT), parameter, public:: TALSH_NOT_AVAILABLE=-888          !information or feature not avaiable (in principle)
        integer(C_INT), parameter, public:: TALSH_NOT_IMPLEMENTED=-999        !feature not implemented yet
        integer(C_INT), parameter, public:: TALSH_NOT_INITIALIZED=1000000     !TALSH library has not been initialized yet
        integer(C_INT), parameter, public:: TALSH_ALREADY_INITIALIZED=1000001 !TALSH library has already been initialized
        integer(C_INT), parameter, public:: TALSH_INVALID_ARGS=1000002        !invalid arguments passed to a procedure
        integer(C_INT), parameter, public:: TALSH_INTEGER_OVERFLOW=1000003    !integer overflow occurred
        integer(C_INT), parameter, public:: TALSH_OBJECT_NOT_EMPTY=1000004    !object is not empty while expected so
        integer(C_INT), parameter, public:: TALSH_OBJECT_IS_EMPTY=1000005     !object is empty while not expected so
 !Host argument buffer:
        integer(C_SIZE_T), parameter, private:: HAB_SIZE_DEFAULT=1048576 !default size of the Host argument buffer in bytes
!DERIVED TYPES:
 !TAL-SH tensor block:
        type, bind(C):: talsh_tens_t
         type(C_PTR):: shape_p=C_NULL_PTR      !shape of the tensor block
         type(C_PTR):: dev_rsc=C_NULL_PTR      !list of device resources occupied by the tensor block body on each device
         type(C_PTR):: data_type=C_NULL_PTR    !list of data types for each device location occupied by the tensor body {R4,R8,C4,C8}
         type(C_PTR):: updated=C_NULL_PTR      !last update event number for each existing tensor block copy
         integer(C_INT):: dev_rsc_len=0        !capacity of .dev_rsc[], .data_type[], .updated[]
         integer(C_INT):: ndev=0               !number of devices the tensor block resides on: ndev <= dev_rsc_len
         integer(C_LONG_LONG):: last_update=-1 !last data update event number
         type(C_PTR):: tensF=C_NULL_PTR        !pointer to Fortran <tensor_block_t> (CPU, Intel Xeon Phi): Just a convenient alias to existing data
         type(C_PTR):: tensC=C_NULL_PTR        !pointer to C/C++ <tensBlck_t> (Nvidia GPU): Just a convenient alias to existing data
        end type talsh_tens_t
 !TAL-SH task handle:
        type, bind(C):: talsh_task_t
         type(C_PTR):: task_p=C_NULL_PTR    !pointer to the corresponding task object
         integer(C_INT):: dev_kind=DEV_NULL !device kind (DEV_NULL: uninitialized)
         real(C_DOUBLE):: flops=0d0         !number of floating point operations
         real(C_DOUBLE):: exec_time=0d0     !execution time in seconds
        end type talsh_task_t
!GLOBALS:

!INTERFACES:
        interface
 !TAL-SH control C/C++ API:
  !Initialize TAL-SH:
         integer(C_INT) function talshInit(host_buf_size,host_arg_max,ngpus,gpu_list,nmics,mic_list,namds,amd_list)&
                                          &bind(c,name='talshInit')
          import
          implicit none
          integer(C_SIZE_T), intent(inout):: host_buf_size
          integer(C_INT), intent(out):: host_arg_max
          integer(C_INT), value, intent(in):: ngpus
          integer(C_INT), intent(in):: gpu_list(*)
          integer(C_INT), value, intent(in):: nmics
          integer(C_INT), intent(in):: mic_list(*)
          integer(C_INT), value, intent(in):: namds
          integer(C_INT), intent(in):: amd_list(*)
         end function talshInit
  !Shutdown TAL-SH:
         integer(C_INT) function talshShutdown() bind(c,name='talshShutdown')
          import
          implicit none
         end function talshShutdown
  !Get the flat device Id:
         integer(C_INT) function talshFlatDevId(dev_kind,dev_num) bind(c,name='talshFlatDevId')
          import
          implicit none
          integer(C_INT), value, intent(in):: dev_kind
          integer(C_INT), value, intent(in):: dev_num
         end function talshFlatDevId
  !Get the kind-specific device Id:
         integer(C_INT) function talshKindDevId(dev_id,dev_kind) bind(c,name='talshKindDevId')
          import
          implicit none
          integer(C_INT), value, intent(in):: dev_id
          integer(C_INT), intent(out):: dev_kind
         end function talshKindDevId
  !Query the state of a device:
         integer(C_INT) function talshDeviceState_(dev_num,dev_kind) bind(c,name='talshDeviceState_')
          import
          implicit none
          integer(C_INT), value, intent(in):: dev_num
          integer(C_INT), value, intent(in):: dev_kind
         end function talshDeviceState_
  !Find the least busy device:
         integer(C_INT) function talshDeviceBusyLeast_(dev_kind) bind(c,name='talshDeviceBusyLeast_')
          import
          implicit none
          integer(C_INT), value, intent(in):: dev_kind
         end function talshDeviceBusyLeast_
  !Print run-time TAL-SH statistics for chosen devices:
         integer(C_INT) function talshStats_(dev_id,dev_kind) bind(c,name='talshStats_')
          import
          implicit none
          integer(C_INT), value, intent(in):: dev_id
          integer(C_INT), value, intent(in):: dev_kind
         end function talshStats_
 !TAL-SH tensor block C/C++ API:
  !Check whether a tensor block is empty (only be called on defined tensor blocks!):
         integer(C_INT) function talshTensorIsEmpty(tens_block) bind(c,name='talshTensorIsEmpty')
          import
          implicit none
          type(C_PTR), value, intent(in):: tens_block
         end function talshTensorIsEmpty
  !Construct a tensor block:
         integer(C_INT) function talshTensorConstruct_(tens_block,data_type,tens_rank,tens_dims,dev_id,&
                        ext_mem,in_hab,init_method,init_val_real,init_val_imag) bind(c,name='talshTensorConstruct_')
          import
          implicit none
          type(talsh_tens_t):: tens_block
          integer(C_INT), value, intent(in):: data_type
          integer(C_INT), value, intent(in):: tens_rank
          integer(C_INT), intent(in):: tens_dims(*)
          integer(C_INT), value, intent(in):: dev_id
          type(C_PTR), value:: ext_mem
          integer(C_INT), value, intent(in):: in_hab
          type(C_FUNPTR), value, intent(in):: init_method
          real(C_DOUBLE), value, intent(in):: init_val_real
          real(C_DOUBLE), value, intent(in):: init_val_imag
         end function talshTensorConstruct_
  !Destruct a tensor block:
         integer(C_INT) function talshTensorDestruct(tens_block) bind(c,name='talshTensorDestruct')
          import
          implicit none
          type(talsh_tens_t):: tens_block
         end function talshTensorDestruct
  !Get the volume of the tensor block:
         integer(C_SIZE_T) function talshTensorVolume(tens_block) bind(c,name='talshTensorVolume')
          import
          implicit none
          type(talsh_tens_t), intent(in):: tens_block
         end function talshTensorVolume
  !Get the shape of the tensor block:
         integer(C_INT) function talshTensorShape(tens_block,tens_shape) bind(c,name='talshTensorShape')
          import
          implicit none
          type(talsh_tens_t), intent(in):: tens_block
          type(talsh_tens_shape_t), intent(inout):: tens_shape
         end function talshTensorShape
  !Query the presence of the tensor block on device(s):
         integer(C_INT) function talshTensorPresence_(tens_block,ncopies,copies,data_types,dev_kind,dev_id)&
                                 bind(c,name='talshTensorPresence_')
          import
          implicit none
          type(talsh_tens_t), intent(in):: tens_block
          integer(C_INT), intent(out):: ncopies
          integer(C_INT), intent(inout):: copies(*)
          integer(C_INT), intent(inout):: data_types(*)
          integer(C_INT), value, intent(in):: dev_kind
          integer(C_INT), value, intent(in):: dev_id
         end function talshTensorPresence_

        end interface
!VISIBILITY:
 !TAL-SH device control API:
        public talsh_init
        public talsh_shutdown
        public talsh_flat_dev_id
        public talsh_kind_dev_id
        public talsh_device_state
        public talsh_device_busy_least
        public talsh_stats
 !TAL-SH tensor block API:
        public talsh_tensor_is_empty
        public talsh_tensor_construct
        public talsh_tensor_destruct
        public talsh_tensor_volume
        public talsh_tensor_shape
        public talsh_tensor_presence
 !TAL-SH task API:
!        public talsh_task_clean
!        public talsh_task_dev_id
!        public talsh_task_status
!        public talsh_task_completed
!        public talsh_task_wait
!        public talsh_tasks_wait
!        public talsh_task_time
 !TAL-SH tensor operations:
!        public talsh_tensor_place
!        public talsh_tensor_discard
!        public talsh_tensor_init
!        public talsh_tensor_scale
!        public talsh_tensor_norm1
!        public talsh_tensor_norm2
!        public talsh_tensor_copy
!        public talsh_tensor_add
!        public talsh_tensor_contract

       contains
!Fortran API definitions:
 !TAL-SH control API:
!----------------------------------------------------------------------------------------------
        function talsh_init(host_buf_size,host_arg_max,gpu_list,mic_list,amd_list) result(ierr)
         implicit none
         integer(C_INT):: ierr                                      !out: error code (0:success)
         integer(C_SIZE_T), intent(inout), optional:: host_buf_size !inout: desired size in bytes of the Host Argument Buffer (HAB).
                                                                    !       It will be replaced by the actual size.
         integer(C_INT), intent(out), optional:: host_arg_max       !out: max number of arguments the HAB can contain
         integer(C_INT), intent(in), optional:: gpu_list(1:)        !in: list of NVidia GPU's to use
         integer(C_INT), intent(in), optional:: mic_list(1:)        !in: list of Intel Xeon Phi's to use
         integer(C_INT), intent(in), optional:: amd_list(1:)        !in: list of AMD GPU's to use
         integer(C_INT):: ngpus,gpus(MAX_GPUS_PER_NODE)
         integer(C_INT):: nmics,mics(MAX_MICS_PER_NODE)
         integer(C_INT):: namds,amds(MAX_AMDS_PER_NODE)
         integer(C_SIZE_T):: hbuf_size
         integer(C_INT):: harg_max

         if(present(host_buf_size)) then; hbuf_size=host_buf_size; else; hbuf_size=HAB_SIZE_DEFAULT; endif
         if(present(gpu_list)) then; ngpus=size(gpu_list); gpus(1:ngpus)=gpu_list(1:ngpus); else; ngpus=0; endif
         if(present(mic_list)) then; nmics=size(mic_list); mics(1:nmics)=mic_list(1:nmics); else; nmics=0; endif
         if(present(amd_list)) then; namds=size(amd_list); amds(1:namds)=amd_list(1:namds); else; namds=0; endif
         ierr=talshInit(hbuf_size,harg_max,ngpus,gpus,nmics,mics,namds,amds)
         if(present(host_arg_max)) host_arg_max=harg_max
         if(present(host_buf_size)) host_buf_size=hbuf_size
         return
        end function talsh_init
!---------------------------------------------
        function talsh_shutdown() result(ierr)
         implicit none
         integer(C_INT):: ierr !out: error code (0:success)
         ierr=talshShutdown()
         return
        end function talsh_shutdown
!----------------------------------------------------------------
        function talsh_flat_dev_id(dev_kind,dev_num) result(res)
         implicit none
         integer(C_INT):: res                  !out: flat device Id [0..DEV_MAX-1]; Failure: DEV_MAX
         integer(C_INT), intent(in):: dev_kind !in: device kind
         integer(C_INT), intent(in):: dev_num  !in: device Id within its kind (0..MAX)
         res=talshFlatDevId(dev_kind,dev_num)
         return
        end function talsh_flat_dev_id
!--------------------------------------------------------------
        function talsh_kind_dev_id(dev_id,dev_kind) result(res)
         implicit none
         integer(C_INT):: res                   !out: kind-specific device Id [0..MAX]; Failure: DEV_NULL (negative)
         integer(C_INT), intent(in):: dev_id    !in: flat device Id
         integer(C_INT), intent(out):: dev_kind !out: device kind
         res=talshKindDevId(dev_id,dev_kind)
         return
        end function talsh_kind_dev_id
!----------------------------------------------------------------------
        function talsh_device_state(dev_num,dev_kind) result(dev_state)
         implicit none
         integer(C_INT):: dev_state                      !out: device state (Success:[DEV_OFF,DEV_ON,DEV_ON_BLAS])
         integer(C_INT), intent(in):: dev_num            !in: either a flat or kind specific (when <dev_kind> is present) device id
         integer(C_INT), intent(in), optional:: dev_kind !in: device kind (note that it changes the meaning of the <dev_num> argument)
         integer(C_INT):: devk

         if(present(dev_kind)) then; devk=dev_kind; else; devk=DEV_NULL; endif
         dev_state=talshDeviceState_(dev_num,devk)
         return
        end function talsh_device_state
!----------------------------------------------------------------
        function talsh_device_busy_least(dev_kind) result(dev_id)
         implicit none
         integer(C_INT):: dev_id                         !out: either a flat or kind specific device id
         integer(C_INT), intent(in), optional:: dev_kind !in: device kind (if absent, <dev_id> will return the flat device id)
         integer(C_INT):: devk

         if(present(dev_kind)) then; devk=dev_kind; else; devk=DEV_NULL; endif
         dev_id=talshDeviceBusyLeast_(devk)
         return
        end function talsh_device_busy_least
!---------------------------------------------------------
        function talsh_stats(dev_id,dev_kind) result(ierr)
         implicit none
         integer(C_INT):: ierr                           !out: error code (0:success)
         integer(C_INT), intent(in), optional:: dev_id   !in: device id (either flat or kind specific device id, see below)
         integer(C_INT), intent(in), optional:: dev_kind !in: device kind (if present, <dev_id> will be interpreted as kind specific)
         integer(C_INT):: devn,devk

         if(present(dev_id)) then; devn=dev_id; else; devn=-1; endif
         if(present(dev_kind)) then; devk=dev_kind; else; devk=DEV_NULL; endif
         ierr=talshStats_(devn,devk)
         return
        end function talsh_stats
!-------------------------------------------------------------
        function talsh_tensor_is_empty(tens_block) result(res)
         implicit none
         logical:: res                                       !out: .TRUE. if the tensor block is empty, .FALSE. otherwise
         type(talsh_tens_t), intent(in), target:: tens_block !in: tensor block

         res=.FALSE.
         if(talshTensorIsEmpty(c_loc(tens_block)).eq.YEP) res=.TRUE.
         return
        end function talsh_tensor_is_empty
!-------------------------------------------------------------------------------------------------------------------------------
        function talsh_tensor_construct(tens_block,data_type,tens_shape,dev_id,ext_mem,in_hab,init_method,init_val) result(ierr)
         implicit none
         integer(C_INT):: ierr
         type(talsh_tens_t), intent(inout):: tens_block       !inout: constructed tensor block (must be empty on entrance)
         integer(C_INT), intent(in):: data_type               !in: data type: {R4,R8,C4,C8,NO_TYPE}
         integer(C_INT), intent(in):: tens_shape(1:)          !in: tensor shape (volume = tensor rank)
         integer(C_INT), intent(in):: dev_id                  !in: flat device ID on which the tensor block will reside
         type(C_PTR), intent(in), optional:: ext_mem          !in: pointer to externally provided memory for tensor elements
         integer(C_INT), intent(in), optional:: in_hab        !in: if >=0, <ext_mem> points to the HAB entry #<in_hab>
         procedure(talsh_tens_init_i), optional:: init_method !in: user-defined initialization method (<init_val> must be absent)
         complex(8), intent(in), optional:: init_val          !in: initialization value (will be typecast to <data_type>, defaults to 0)
         type(C_PTR):: tens_body_p
         integer(C_INT):: hab_entry,tens_rank
         integer(C_INT), target:: tens_dims(1:MAX_TENSOR_RANK)
         type(C_FUNPTR):: init_method_p
         real(8):: val_real,val_imag

         ierr=TALSH_SUCCESS
         tens_rank=size(tens_shape) !tens_shape(1:) must have the exact volume = tensor rank
         if(tens_rank.ge.0.and.tens_rank.le.MAX_TENSOR_RANK) then
          if(tens_rank.gt.0) tens_dims(1:tens_rank)=tens_shape(1:tens_rank)
          if(present(ext_mem)) then; tens_body_p=ext_mem; else; tens_body_p=C_NULL_PTR; endif
          if(present(in_hab)) then; if(in_hab.ge.0) then; hab_entry=in_hab; else; hab_entry=-1; endif; else; hab_entry=-1; endif
          if(present(init_method)) then; init_method_p=c_funloc(init_method); else; init_method_p=C_NULL_FUNPTR; endif
          val_real=0d0; val_imag=0d0; if(present(init_val)) then; val_real=real(init_val); val_imag=aimag(init_val); endif
          ierr=talshTensorConstruct_(tens_block,data_type,tens_rank,tens_dims,dev_id,&
                                     tens_body_p,hab_entry,init_method_p,val_real,val_imag)
         else
          ierr=TALSH_INVALID_ARGS
         endif
         return
        end function talsh_tensor_construct
!--------------------------------------------------------------
        function talsh_tensor_destruct(tens_block) result(ierr)
         implicit none
         integer(C_INT):: ierr
         type(talsh_tens_t), intent(inout):: tens_block !inout: tensor block (will become empty on exit)
         ierr=talshTensorDestruct(tens_block)
         return
        end function talsh_tensor_destruct
!-----------------------------------------------------------
        function talsh_tensor_volume(tens_block) result(vol)
         implicit none
         integer(C_SIZE_T):: vol                     !out: number of elements in the tensor block (negative on error)
         type(talsh_tens_t), intent(in):: tens_block !in: tensor block
         vol=talshTensorVolume(tens_block)
         return
        end function talsh_tensor_volume
!----------------------------------------------------------------------
        function talsh_tensor_shape(tens_block,tens_shape) result(ierr)
         implicit none
         integer(C_INT):: ierr                                !out: error code (0:success)
         type(talsh_tens_t), intent(in):: tens_block          !in: tensor block
         type(talsh_tens_shape_t), intent(inout):: tens_shape !inout: tensor block shape (copy)
         ierr=talshTensorShape(tens_block,tens_shape)
         return
        end function talsh_tensor_shape
!--------------------------------------------------------------------------------------------------------
        function talsh_tensor_presence(tens_block,ncopies,copies,data_types,dev_kind,dev_id) result(ierr)
         integer(C_INT):: ierr                           !out: error code (0:success)
         type(talsh_tens_t), intent(in):: tens_block     !in: tensor block
         integer(C_INT), intent(out):: ncopies           !out: number of found copies of the tensor block
         integer(C_INT), intent(inout):: copies(1:*)     !out: copies found (list of flat device id's)
         integer(C_INT), intent(inout):: data_types(1:*) !out: data type of each copy
         integer(C_INT), intent(in), optional:: dev_kind !in: specific device kind of interest (defaults to All)
         integer(C_INT), intent(in), optional:: dev_id   !in: specific device of interest
         integer(C_INT):: devk,devnum

         if(present(dev_kind)) then; devk=dev_kind; else; devk=DEV_NULL; endif
         if(present(dev_id)) then; devnum=dev_id; else; devnum=-1; endif
         ierr=talshTensorPresence_(tens_block,ncopies,copies,data_types,devk,devnum)
         return
        end function talsh_tensor_presence

       end module talsh
