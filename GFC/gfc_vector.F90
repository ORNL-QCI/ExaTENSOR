!Generic Fortran Containers (GFC): Vector
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2017/01/08

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
!TYPES:
#if 0
 !Vector element:
        type, extends(gfc_cont_elem_t), public:: vector_elem_t
        contains
         procedure, public:: construct=>VectorElemConstruct !constructs the content of the vector element
         procedure, public:: is_first=>VectorElemIsFirst    !returns GFC_TRUE if the element is the first
         procedure, public:: is_last=>VectorElemIsLast      !returns GFC_TRUE if the element is the last
        end type vector_elem_t
 !Vector segment:
        type, private:: vector_seg_t
         integer(INTL), private:: max_seg_len=0_INTL             !max length of the vector segment
         integer(INTL), private:: seg_len=0_INTL                 !current occupied length of the vector segment
         type(vector_elem_t), allocatable, private:: seg_elem(:) !elements of the vector segment
        end type vector_seg_t
 !Vector:
        type, extends(gfc_container_t), public:: vector_t
         integer(INTL), private:: num_segments=0_INTL      !number of segments in the vector
         integer(INTL), private:: segment_length=0_INTL    !length of each vector segment
         integer(INTL), private:: lower_bound=0_INTL       !lower bound of the vector
         integer(INTL), private:: upper_bound=-1_INTL      !upper bound of the vector
         type(vector_seg_t), allocatable, private:: seg(:) !vector segments
         contains
          procedure, public:: max_length=>VectorMaxLength   !maximal length of the vector
          procedure, public:: length=>VectorLength          !current vector length
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
          procedure, public:: insert=>VectorIterInsert            !inserts a new element anywhere in the vector (current iterator position)
          procedure, public:: delete=>VectorIterDelete            !deletes an element anywhere in the vector (current iterator position)
          procedure, public:: delete_all=>VectorIterDeleteAll     !deletes all elements of the vector
        end type vector_iter_t
#endif
       end module gfc_vector
