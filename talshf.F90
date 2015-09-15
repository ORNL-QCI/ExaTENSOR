!ExaTensor::TAL-SH: User-level API definition:
!REVISION: 2015/09/15
!Copyright (C) 2015 Dmitry I. Lyakh (email: quant4me@gmail.com)
!Copyright (C) 2015 Oak Ridge National Laboratory (UT-Battelle)

!This source file is free software; you can redistribute it and/or
!modify it under the terms of the GNU General Public License
!as published by the Free Software Foundation; either version 2
!of the License, or (at your option) any later version.

!This program is distributed in the hope that it will be useful,
!but WITHOUT ANY WARRANTY; without even the implied warranty of
!MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!GNU General Public License for more details.

!You should have received a copy of the GNU General Public License
!along with this program; if not, write to the Free Software
!Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
--------------------------------------------------------------------------------
       module talsh
        use tensor_algebra_cpu_phi
        implicit none
!       private
!PARAMETERS:
 !Generic:
        integer(INTD), private:: CONS_OUT=6   !default output device for this module
        logical(INTD), private:: DEBUG=.true. !debugging mode for this module
!DERIVED TYPES:

!GLOBALS:

!INTERFACES:

!PROCEDURE VISIBILITY:
 !Device control API:
        public talsh_init
        public talsh_shutdown
        public talsh_flat_dev_id
        public talsh_kind_dev_id
        public talsh_device_state
        public talsh_device_busy_least
        public talsh_stats
 !Tensor block API:
        public talsh_tensor_construct
        public talsh_tensor_destroy
        public talsh_tensor_volume
        public talsh_tensor_datatype
        public talsh_tensor_shape
        public talsh_tensor_presence
 !TAL-SH task API:
        public talsh_task_clean
        public talsh_task_dev_id
        public talsh_task_status
        public talsh_task_completed
        public talsh_task_wait
        public talsh_tasks_wait
        public talsh_task_time
 !Tensor operations:
        public talsh_tensor_place
        public talsh_tensor_discard
        public talsh_tensor_init
        public talsh_tensor_scale
        public talsh_tensor_norm1
        public talsh_tensor_norm2
        public talsh_tensor_copy
        public talsh_tensor_add
        public talsh_tensor_contract

       contains
!API DEFINTION:

       end module talsh
