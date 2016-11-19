!Generic Fortran Containers (GFC): Dictionary (ordered map), AVL BST
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2016/11/19 (recycle of my old dictionary implementation)

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

!FOR DEVELOPERS ONLY:

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
        integer(INTD), parameter, public:: GFC_DICT_JUST_FIND=0        !action: just find the key in the dictionary, if it is there
        integer(INTD), parameter, public:: GFC_DICT_DELETE_IF_FOUND=1  !action: delete the dictionary element, if key found
        integer(INTD), parameter, public:: GFC_DICT_REPLACE_IF_FOUND=2 !action: replace the value of the dictionary element, if key found
        integer(INTD), parameter, public:: GFC_DICT_ADD_IF_NOT_FOUND=3 !action: create a new dictionary element, if key not found
        integer(INTD), parameter, public:: GFC_DICT_ADD_OR_MODIFY=4    !action: combines ADD_IF_NOT_FOUND and REPLACE_IF_FOUND
        integer(INTD), parameter, public:: GFC_DICT_FETCH_IF_FOUND=5   !action: simply fetch the dictionary element value, if key found
!TYPES:
 !Dictionary element:
        type, extends(gfc_cont_elem_t), public:: dict_elem_t
         class(*), pointer, private:: key=>NULL()                !dictionary element key
         class(dict_elem_t), pointer, private:: child_lt=>NULL() !lesser child element
         class(dict_elem_t), pointer, private:: child_gt=>NULL() !greater child element
         class(dict_elem_t), pointer, private:: parent=>NULL()   !parent element
         integer(INTD), private:: balance_fac                    !balance factor (for AVL BST)
         contains
          procedure, public:: construct=>DictElemConstruct             !constructs a dictionary element from a key-value pair
          procedure, non_overridable, public:: is_root=>DictElemIsRoot !returns GFC_TRUE if the element is the root of the dictionary binary tree
          procedure, non_overridable, public:: is_leaf=>DictElemIsLeaf !returns GFC_TRUE if the element is a leaf of the dictionary binary tree
          procedure, non_overridable, public:: get_key=>DictElemGetKey !returns an unlimited polymorphic pointer to the element key
#if 0
          procedure, public:: predicate_key=>DictElemPredicateKey      !returns the result of predication on the element key
          procedure, public:: DictElemActionKey                        !performs a user-defined action on the element key via a function
          procedure, public:: DictElemFunctorKey                       !performs a user-defined action on the element key via a functor
          procedure, public:: action_key=>DictElemActionKey,DictElemFunctorKey !performs a user-defined action: generic overload
          procedure, public:: compare_key=>DictElemCompareKey          !compares the element key with another key
          procedure, public:: print_key=>DictElemPrintKey              !prints the element key
#endif
        end type dict_elem_t
 !Dictionary (all operations on the dictionary are performed via an iterator):
        type, extends(gfc_container_t), public:: dictionary_t
         class(dict_elem_t), pointer, private:: root=>NULL()           !root (boundary) element (beginning)
        end type dictionary_t
 !Dictionary iterator:
#if 0
        type, extends(gfc_iter_t), public:: dictionary_iter_t
         class(dict_elem_t), pointer, private:: current=>NULL()        !currently pointed element of the container
         class(dictionary_t), pointer, private:: container=>NULL()     !container
         contains
          procedure, public:: init=>DictionaryIterInit                 !associates the iterator with a container and positions it to the root element
          procedure, public:: reset=>DictionaryIterReset               !resets the iterator to the beginning of the container (root element)
          procedure, public:: release=>DictionaryIterRelease           !dissociates the iterator from its container
          procedure, public:: pointee=>DictionaryIterPointee           !returns a pointer to the container element currently in focus
          procedure, public:: next=>DictionaryIterNext                 !moves the iterator to the "next" element, if any (not necessarily in order)
          procedure, public:: previous=>DictionaryIterPrevious         !moves the iterator to the "previous" element, if any (not necessarily in order)
          procedure, public:: next_in_order=>DictionaryIterNextInOrder !moves the iterator to the next-in-order element, if any
          procedure, public:: prev_in_order=>DictionaryIterPrevInOrder !moves the iterator to the previous-in-order element, if any
          procedure, public:: move_to_min=>DictionaryIterMoveToMin     !moves the iterator to the leftmost (minimal element)
          procedure, public:: move_to_max=>DictionaryIterMoveToMax     !moves the iterator to the rightmost (maximal element)
          procedure, public:: search=>DictionaryIterSearch             !performs a key-based search in the dictionary
        end type dictionary_iter_t
#endif
!INTERFACES:

!VISIBILITY:
 !dict_elem_t:
        private DictElemConstruct
        private DictElemIsRoot
        private DictElemIsLeaf
        private DictElemGetKey
#if 0
        private DictElemPredicateKey
        private DictElemActionKey
        private DictElemFunctorKey
        private DictElemCompareKey
        private DictElemPrintKey
#endif
 !dictionary_iter_t:
#if 0
        private DictionaryIterInit
        private DictionaryIterReset
        private DictionaryIterRelease
        private DictionaryIterPointee
        private DictionaryIterNext
        private DictionaryIterPrevious
        private DictionaryIterNextInOrder
        private DictionaryIterPrevInOrder
        private DictionaryIterMoveToMin
        private DictionaryIterMoveToMax
        private DictionaryIterSearch
#endif
!DEFINITIONS:
       contains
![dict_elem_t]===================================================================================
#ifdef NO_GNU
        subroutine DictElemConstruct(this,key,val,ierr,assoc_val,key_copy_ctor_f,val_copy_ctor_f) !`GCC has a bug with this line
#else
        subroutine DictElemConstruct(this,key,val,ierr,assoc_val)
#endif
!Given a key-value pair, constructs an element of a dictionary. Note
!that if construction fails, the element may be left underconstructed,
!requiring a separate call to the destructor after return.
         implicit none
         class(dict_elem_t), intent(inout):: this           !inout: element of the dictionary
         class(*), target, intent(in):: key                 !in: key to be stored (by value only)
         class(*), target, intent(in):: val                 !in: value to be stored (either by value or by reference)
         integer(INTD), intent(out), optional:: ierr        !out: error code
         logical, intent(in), optional:: assoc_val          !in: if TRUE, <val> will be stored by reference, otherwise by value (default)
#ifdef NO_GNU
         procedure(gfc_copy_i), optional:: key_copy_ctor_f  !in: user-defined generic copy constructor for the key
         procedure(gfc_copy_i), optional:: val_copy_ctor_f  !in: user-defined generic copy constructor for the value
#endif
         integer(INTD):: errc
         integer:: ier
         logical:: assoc

         errc=GFC_SUCCESS
         if(present(assoc_val)) then; assoc=assoc_val; else; assoc=.FALSE.; endif
         if(this%in_use(errc,set_lock=.TRUE.,report_refs=.FALSE.).eq.GFC_FALSE) then
          if(this%is_empty()) then
#ifdef NO_GNU
           if(present(val_copy_ctor_f)) then
            call this%construct_base(val,errc,assoc_only=assoc,copy_ctor_f=val_copy_ctor_f,locked=.TRUE.)
           else
#endif
            call this%construct_base(val,errc,assoc_only=assoc,locked=.TRUE.)
#ifdef NO_GNU
           endif
#endif
           if(errc.eq.GFC_SUCCESS) then !base constructor succeeded
#ifdef NO_GNU
            if(present(key_copy_ctor_f)) then
             this%key=>key_copy_ctor_f(key,errc)
            else
#endif
             allocate(this%key,SOURCE=key,STAT=ier)
             if(ier.ne.0) errc=GFC_MEM_ALLOC_FAILED
#ifdef NO_GNU
            endif
#endif
           endif
          else
           errc=GFC_ELEM_NOT_EMPTY
          endif
          if(errc.eq.GFC_SUCCESS) then
           call this%release_lock_(errc)
          else
           call this%release_lock_()
          endif
         else
          errc=GFC_IN_USE
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DictElemConstruct
!------------------------------------------------
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
![dictionary_t]=========================

![dictionary_iter_t]====================

       end module gfc_dictionary
