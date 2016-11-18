!Generic Fortran Containers (GFC): Dictionary (ordered map)
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2016/11/18

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

       module gfc_dictionary
        use gfc_base
        use timers
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !output device
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
        integer(INTD), private:: DEBUG=0    !debugging level (0:none)
 !Directions:
        logical, parameter, public:: GFC_DICT_LEFT=.FALSE.
        logical, parameter, public:: GFC_DICT_RIGHT=.TRUE.
        logical, parameter, public:: GFC_DICT_SUCCESSOR=.TRUE.
        logical, parameter, public:: GFC_DICT_PREDECESSOR=.FALSE.
 !Search:
        integer(INTD), parameter, public:: GFC_DICT_JUST_FIND=0
        integer(INTD), parameter, public:: GFC_DICT_DELETE_IF_FOUND=1
        integer(INTD), parameter, public:: GFC_DICT_REPLACE_IF_FOUND=2
        integer(INTD), parameter, public:: GFC_DICT_ADD_IF_NOT_FOUND=3
        integer(INTD), parameter, public:: GFC_DICT_ADD_OR_MODIFY=4
        integer(INTD), parameter, public:: GFC_DICT_FETCH_IF_FOUND=5
!TYPES:
 !Dictionary element:
        type, extends(gfc_cont_elem_t), public:: dict_elem_t
         class(*), pointer, private:: key=>NULL()                !dictionary element key
         class(dict_elem_t), pointer, private:: child_lt=>NULL() !lesser child element
         class(dict_elem_t), pointer, private:: child_gt=>NULL() !greater child element
         class(dict_elem_t), pointer, private:: parent=>NULL()   !parent element
         integer(INTD), private:: balance_fac                    !balance factor
         contains
          procedure, non_overridable, public:: is_root=>DictElemIsRoot !returns GFC_TRUE if the element is the root of the dictionary binary tree
          procedure, non_overridable, public:: is_leaf=>DictElemIsLeaf !returns GFC_TRUE if the element is a leaf of the dictionary binary tree
          procedure, non_overridable, public:: get_key=>DictElemGetKey !returns an unlimited polymorphic pointer to the element key
          procedure, public:: construct=>DictElemConstruct             !constructs a dictionary element
          procedure, public:: destruct=>DictElemDestruct               !destructs a dictionary element
          procedure, public:: predicate_key=>DictElemPredicateKey      !returns the result of predication on the element key
          procedure, public:: action_key=>DictElemActionKey            !performs a user-defined action on the element key
          procedure, public:: compare_key=>DictElemCompareKey          !compares the element key with another key
          procedure, public:: print_key=>DictElemPrintKey              !prints the element key
        end type dict_elem_t
 !Dictionary (all operations on the dictionary are performed via an iterator):
        type, extends(gfc_container_t), public:: dictionary_t
         class(dict_elem_t), pointer, private:: root=>NULL()       !root (boundary) element (beginning)
        end type dictionary_t
 !Dictionary iterator:
        type, extends(gfc_iter_t), public:: dictionary_iter_t
         class(dict_elem_t), pointer, private:: current=>NULL()    !currently pointed element of the container
         class(dictionary_t), pointer, private:: container=>NULL() !container
         contains
          procedure, public:: init=>DictionaryIterInit                 !associates the iterator with a container and positions it to the root element
          procedure, public:: reset=>DictionaryIterReset               !resets the iterator to the beginning of the container (root element)
          procedure, public:: release=>DictionaryIterRelease           !dissociates the iterator from its container
          procedure, public:: pointee=>DictionaryIterPointee           !returns a pointer to the container element currently in focus
          procedure, public:: next=>DictionaryIterNext                 !moves the iterator to the "next" element, if any (not necessarily in order)
          procedure, public:: previous=>DictionaryIterPrevious         !moves the iterator to the "previous" element, if any (not necessarily in order)
          procedure, public:: next_in_order=>DictionaryIterNextInOrder !moves the iterator to the next in order element, if any
          procedure, public:: prev_in_order=>DictionaryIterPrevInOrder !moves the iterator to the previous in order element, if any
          procedure, public:: search=>DictionaryIterSearch             !performs a key-based search in the dictionary
        end type dictionary_iter_t
!INTERFACES:

!VISIBILITY:
 !dict_elem_t:
        private DictElemIsRoot
        private DictElemIsLeaf
        private DictElemGetKey
        private DictElemConstruct
        private DictElemDestruct
        private DictElemPredicateKey
        private DictElemActionKey
        private DictElemCompareKey
        private DictElemPrintKey
 !dictionary_iter_t:
        private DictionaryIterInit
        private DictionaryIterReset
        private DictionaryIterRelease
        private DictionaryIterPointee
        private DictionaryIterNext
        private DictionaryIterPrevious
        private DictionaryIterNextInOrder
        private DictionaryIterPrevInOrder
        private DictionaryIterSearch
!DEFINITIONS:
       contains
![dict_elem_t]===================================
        function DictElemIsRoot(this) result(res)
         implicit none
         integer(INTD):: res                   !out: answer {GFC_TRUE,GFC_FALSE,GFC_ERROR}
         class(dict_elem_t), intent(in):: this !in: dictionary element

         if(associated(this%parent)) then
          res=GFC_FALSE
         else
          res=GFC_TRUE
         endif
         return
        end function DictElemIsRoot
!------------------------------------------------
        function DictElemIsLeaf(this) result(res)
         implicit none
         integer(INTD):: res                   !out: answer {GFC_TRUE,GFC_FALSE,GFC_ERROR}
         class(dict_elem_t), intent(in):: this !in: dictionary element

         if(associated(this%child_lt).or.associated(this%child_gt)) then
          res=GFC_FALSE
         else
          res=GFC_TRUE
         endif
         return
        end function DictElemIsLeaf
!-------------------------------------------------------
        function DictElemGetKey(this,ierr) result(key_p)
         implicit none
         class(*), pointer:: key_p                   !out: pointer to the element key
         class(dict_elem_t), intent(in):: this       !in: dictionary element
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS; key_p=>NULL()
         if(.not.this%is_empty()) then
          if(associated(this%key)) then
           key_p=>this%key
          else
           errc=GFC_CORRUPTED_CONT
          endif
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end function DictElemGetKey
![dictionary_t]=======================

![dictionary_iter_t]==================

       end module gfc_dictionary
