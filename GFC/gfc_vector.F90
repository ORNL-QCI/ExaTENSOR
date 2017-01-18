!Generic Fortran Containers (GFC): Vector
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2017/01/18

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

       module gfc_vector
        use gfc_base
        use timers
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !output device
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
        integer(INTD), private:: DEBUG=0    !debugging level (0:none)
 !Vector length:
        integer(INTL), parameter, public:: GFC_VECTOR_MAX_LENGTH=(2_INTL)**36 !default vector length limit
!TYPES:
#if 0
 !Vector element:
        type, extends(gfc_cont_elem_t), public:: vector_elem_t
        contains
         procedure, public:: construct=>VectorElemConstruct !constructs the content (value) of the vector element
        end type vector_elem_t
 !Vector segment:
        type, private:: vector_seg_t
         integer(INTD), private:: max_elems=0                    !max length of the vector segment
         integer(INTD), private:: num_elems=0                    !current occupied length of the vector segment
         type(vector_elem_t), allocatable, private:: seg_elem(:) !elements of the vector segment: [0..max_elems-1]
         contains
          procedure, private:: construct=>VectorSegConstruct     !constructs a vector segment
          final:: VectorSegDestruct                              !destructs a vector segment
        end type vector_seg_t
 !Vector segment batch:
        type, private:: vector_batch_t
         integer(INTD), private:: max_segs=0                     !max number of segments in the batch
         integer(INTD), private:: num_segs=0                     !current number of active segments in the batch
         type(vector_seg_t), allocatable, private:: batch_seg(:) !segments of the batch: [0..max_segs-1]
         contains
          procedure, private:: construct=>VectorBatchConstruct   !constructs a vector batch
          final:: VectorBatchDestruct                            !destructs a vector batch
        end type vector_batch_t
 !Vector:
        type, extends(gfc_container_t), public:: vector_t
         integer(INTL), private:: lbnd=0_INTL                      !lower bound of the vector
         integer(INTL), private:: ubnd=-1_INTL                     !upper bound of the vector
         integer(INTD), private:: max_batches=0                    !max number of batches in the vector
         integer(INTD), private:: num_batches=0                    !current number of active batches in the vector
         type(vector_batch_t), allocatable, private:: vec_batch(:) !batches
         contains
          procedure, public:: capacity=>VectorCapacity      !maximal length of the vector
          procedure, public:: length=>VectorLength          !current vector length = (upper_bound - lower_bound + 1)
          procedure, public:: lower_bound=>VectorLowerBound !vector lower bound
          procedure, public:: upper_bound=>VectorUpperBound !vector upper bound
        end type vector_t
 !Vector iterator:
        type, extends(gfc_iter_t), public:: vector_iter_t
         class(vector_elem_t), pointer, private:: current=>NULL() !current element of the vector
         class(vector_t), pointer, private:: container=>NULL()    !vector associated with the iterator
         contains
          procedure, public:: init=>VectorIterInit                !initializes the iterator by associating it with a vector
          procedure, public:: reset=>VectorIterReset              !resets the iterator to the beginning of the vector
          procedure, public:: reset_back=>VectorIterResetBack     !resets the iterator to the end of the vector
          procedure, public:: release=>VectorIterRelease          !releases the iterator (dissocaites it from its container)
          procedure, public:: pointee=>VectorIterPointee          !returns the container element currently pointed to by the iterator
          procedure, public:: next=>VectorIterNext                !moves the iterator to the next vector element
          procedure, public:: previous=>VectorIterPrevious        !moves the iterator to the previous vector element
          procedure, public:: move_to=>VectorIterMoveTo           !moves the iterator to the specific vector element
          procedure, public:: append=>VectorIterAppend            !appends a new element at the end of the vector
          procedure, public:: insert=>VectorIterInsert            !inserts a new element at the current iterator position
          procedure, public:: delete=>VectorIterDelete            !deletes an element at the current iterator position
          procedure, public:: delete_all=>VectorIterDeleteAll     !deletes all elements of the vector
        end type vector_iter_t
!VISIBILITY:
 !vector_elem_t:
        private VectorElemConstruct
 !vector_seg_t:
        private VectorSegConstruct
 !vector_batch_t:
        private VectorBatchConstruct
 !vector_t:
        private VectorCapacity
        private VectorLength
        private VectorLowerBound
        private VectorUpperBound
 !vector_iter_t:
        private VectorIterInit
        private VectorIterReset
        private VectorIterResetBack
        private VectorIterRelease
        private VectorIterPointee
        private VectorIterNext
        private VectorIterPrevious
        private VectorIterMoveTo
        private VectorIterAppend
        private VectorIterInsert
        private VectorIterDelete
        private VectorIterDeleteAll

       contains
!IMPLEMENTATION:
![vector_elem_t]============================================================
#ifdef NO_GNU
        subroutine VectorElemConstruct(this,obj,ierr,assoc_only,copy_ctor_f) !`GCC has a bug with this line
#else
        subroutine VectorElemConstruct(this,obj,ierr,assoc_only)
#endif
!Constructs the content of the vector element.
         implicit none
         class(vector_elem_t), intent(inout):: this     !inout: vector element
         class(*), target, intent(in):: obj             !in: assigned value
         integer(INTD), intent(out), optional:: ierr    !out: error code
         logical, intent(in), optional:: assoc_only     !in: if TRUE, the value will be assigned by reference, otherwise by value (allocated): Defaults to FALSE
#ifdef NO_GNU
         procedure(gfc_copy_i), optional:: copy_ctor_f  !in: generic copy constructor
#endif
         integer(INTD):: errc

#ifdef NO_GNU
         if(present(copy_ctor_f)) then
          if(present(assoc_only)) then
           call this%construct_base(obj,errc,assoc_only,copy_ctor_f)
          else
           call this%construct_base(obj,errc,copy_ctor_f=copy_ctor_f)
          endif
         else
#endif
          if(present(assoc_only)) then
           call this%construct_base(obj,errc,assoc_only=assoc_only)
          else
           call this%construct_base(obj,errc)
          endif
#ifdef NO_GNU
         endif
#endif
         if(present(ierr)) ierr=errc
         return
        end subroutine VectorElemConstruct
!--------------------------------------------------------
        subroutine VectorSegConstruct(this,capacity,ierr)
         implicit none
         class(vector_seg_t), intent(inout):: this   !out: vector segment
         integer(INTD), intent(in):: capacity        !in: requested capacity
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
         if(this%max_elems.eq.0) then
          if(capacity.gt.0) then
           allocate(this%seg_elem(1:capacity),STAT=errc)
           if(errc.eq.0) then
            this%max_elems=capacity
            this%num_elems=0
           else
            errc=GFC_MEM_ALLOC_FAILED
           endif
          else
           errc=GFC_INVALID_ARGS
          endif
         else
          errc=GFC_ERROR
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine VectorSegConstruct


#endif
       end module gfc_vector
