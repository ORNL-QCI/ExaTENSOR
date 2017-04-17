!Generic Fortran Containers (GFC): Dictionary (ordered map), AVL BST
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2017/04/16 (recycling my old dictionary implementation)

!Copyright (C) 2014-2017 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2017 Oak Ridge National Laboratory (UT-Battelle)

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
        use gfc_list
        use gfc_vector
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
          procedure, private:: DictElemConstruct
          generic, public:: dict_elem_ctor=>DictElemConstruct            !constructs a dictionary element from a key-value pair
          procedure, public:: destruct_keyval=>DictElemDestruct          !destructs a dictionary element (both key and value)
          procedure, non_overridable, public:: is_root=>DictElemIsRoot   !returns GFC_TRUE if the element is the root of the dictionary binary search tree
          procedure, non_overridable, public:: is_leaf=>DictElemIsLeaf   !returns GFC_TRUE if the element is a leaf of the dictionary binary search tree
          procedure, public:: get_key=>DictElemGetKey                    !returns an unlimited polymorphic pointer to the element key
          procedure, public:: predicate_key=>DictElemPredicateKey        !returns the result of predication on the element key
          procedure, public:: compare_key=>DictElemCompareKey            !compares the element key with another key
          procedure, public:: print_key=>DictElemPrintKey                !prints the dictionary element key
          procedure, public:: print_it=>DictElemPrintIt                  !prints the dictionary elemet (key,value)
        end type dict_elem_t
 !Dictionary (all operations on the dictionary are performed via an iterator):
        type, extends(gfc_container_t), public:: dictionary_t
         class(dict_elem_t), pointer, private:: root=>NULL()           !root of the AVL binary search tree (boundary element)
!        class(dict_elem_t), pointer, private:: first=>NULL()          !first element (boundary element) `Do I need this?
!        class(dict_elem_t), pointer, private:: last=>NULL()           !last element (boundary element) `Do I need this?
         logical, private:: key_storage=GFC_BY_VAL                     !dictionary key storage policy
         contains
          procedure, public:: is_empty=>DictionaryIsEmpty                 !returns GFC_TRUE if the dictionary is empty, GFC_FALSE otherwise (or error code)
          procedure, public:: is_subdictionary=>DictionaryIsSubdictionary !returns TRUE if the dictionary is subdictionary, FALSE otherwise
          procedure, public:: set_key_storage=>DictionarySetKeyStorage    !sets the key storage policy (by value or by reference), the dictionary must be empty
          procedure, private:: reroot_=>DictionaryReroot                  !PRIVATE: changes the root of the dictionary
        end type dictionary_t
 !Dictionary iterator:
        type, extends(gfc_iter_t), public:: dictionary_iter_t
         class(dict_elem_t), pointer, private:: current=>NULL()        !currently pointed element of the container
         class(dictionary_t), pointer, private:: container=>NULL()     !container
         contains
          procedure, private:: jump_=>DictionaryIterJump                 !PRIVATE: moves the iterator to an arbitrary position
          procedure, public:: init=>DictionaryIterInit                   !associates the iterator with a container and positions it to the root element
          procedure, public:: reset=>DictionaryIterReset                 !resets the iterator to the beginning of the container (root element)
          procedure, public:: release=>DictionaryIterRelease             !dissociates the iterator from its container
          procedure, public:: pointee=>DictionaryIterPointee             !returns a pointer to the container element currently in focus
          procedure, public:: next=>DictionaryIterNext                   !moves the iterator to the "next" element, if any (not necessarily in order)
          procedure, public:: previous=>DictionaryIterPrevious           !moves the iterator to the "previous" element, if any (not necessarily in order)
          procedure, public:: next_in_order=>DictionaryIterNextInOrder   !moves the iterator to the next-in-order element, if any
          procedure, public:: prev_in_order=>DictionaryIterPrevInOrder   !moves the iterator to the previous-in-order element, if any
          procedure, public:: move_in_order=>DictionaryIterMoveInOrder   !moves the iterator to either the next-in-order or previous-in-order element
          procedure, public:: move_to_min=>DictionaryIterMoveToMin       !moves the iterator to the minimal element
          procedure, public:: move_to_max=>DictionaryIterMoveToMax       !moves the iterator to the maximal element
          procedure, public:: move_up=>DictionaryIterMoveUp              !moves the iterator up the binary search tree (to the parent element)
          procedure, public:: move_down=>DictionaryIterMoveDown          !moves the iterator down the binary search tree, either left or right
          procedure, public:: get_key=>DictionaryIterGetKey              !returns a pointer to the key in the current iterator position
          procedure, public:: delete_all=>DictionaryIterDeleteAll        !deletes all elements of the dictionary
          procedure, public:: search=>DictionaryIterSearch               !performs a key-based search in the dictionary
          procedure, public:: sort_to_list=>DictionaryIterSortToList     !returns a list of references to dictionary elements in a sorted (by key) order
          procedure, public:: sort_to_vector=>DictionaryIterSortToVector !returns a vector of references to dictionary elements in a sorted (by key) order
        end type dictionary_iter_t
!INTERFACES:

!VISIBILITY:
 !dict_elem_t:
        private DictElemConstruct
        private DictElemDestruct
        private DictElemIsRoot
        private DictElemIsLeaf
        private DictElemGetKey
        private DictElemPredicateKey
        private DictElemCompareKey
        private DictElemPrintKey
        private DictElemPrintIt
 !dictionary_t:
        private DictionaryIsEmpty
        private DictionaryIsSubdictionary
        private DictionarySetKeyStorage
        private DictionaryReroot
 !dictionary_iter_t:
        private DictionaryIterJump
        private DictionaryIterInit
        private DictionaryIterReset
        private DictionaryIterRelease
        private DictionaryIterPointee
        private DictionaryIterNext
        private DictionaryIterPrevious
        private DictionaryIterNextInOrder
        private DictionaryIterPrevInOrder
        private DictionaryIterMoveInOrder
        private DictionaryIterMoveToMin
        private DictionaryIterMoveToMax
        private DictionaryIterMoveUp
        private DictionaryIterMoveDown
        private DictionaryIterGetKey
        private DictionaryIterDeleteAll
        private DictionaryIterSearch
        private DictionaryIterSortToList
        private DictionaryIterSortToVector
!DEFINITIONS:
       contains
![dict_elem_t]=============================================================================================
#ifdef NO_GNU
        subroutine DictElemConstruct(this,key,val,ierr,assoc_key,assoc_val,key_copy_ctor_f,val_copy_ctor_f) !`GCC has a bug with this line
#else
        subroutine DictElemConstruct(this,key,val,ierr,assoc_key,assoc_val)
#endif
!Given a key-value pair, constructs an element of a dictionary. Note
!that if construction fails, the element may be left underconstructed,
!requiring a separate call to the destructor after return.
         implicit none
         class(dict_elem_t), intent(inout):: this           !inout: element of the dictionary
         class(*), target, intent(in):: key                 !in: key to be stored (by value only)
         class(*), target, intent(in):: val                 !in: value to be stored (either by value or by reference)
         integer(INTD), intent(out), optional:: ierr        !out: error code
         logical, intent(in), optional:: assoc_key          !in: if TRUE, <key> will be stored by reference, otherwise by value (default)
         logical, intent(in), optional:: assoc_val          !in: if TRUE, <val> will be stored by reference, otherwise by value (default)
#ifdef NO_GNU
         procedure(gfc_copy_i), optional:: key_copy_ctor_f  !in: user-defined generic copy constructor for the key (by value)
         procedure(gfc_copy_i), optional:: val_copy_ctor_f  !in: user-defined generic copy constructor for the value (by value)
#endif
         integer(INTD):: errc
         integer:: ier
         logical:: assk,assv

         errc=GFC_SUCCESS
         if(present(assoc_key)) then; assk=assoc_key; else; assk=.FALSE.; endif
         if(present(assoc_val)) then; assv=assoc_val; else; assv=.FALSE.; endif
         if(this%in_use(errc,set_lock=.TRUE.,report_refs=.FALSE.).eq.GFC_FALSE) then
          if(this%is_empty()) then
#ifdef NO_GNU
           if(present(val_copy_ctor_f)) then
            call this%construct_base(val,errc,assoc_only=assv,copy_ctor_f=val_copy_ctor_f,locked=.TRUE.)
           else
#endif
            call this%construct_base(val,errc,assoc_only=assv,locked=.TRUE.)
#ifdef NO_GNU
           endif
#endif
           if(errc.eq.GFC_SUCCESS) then !base constructor succeeded
            if(assk) then
             this%key=>key
            else
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
           endif
          else
           errc=GFC_ELEM_NOT_EMPTY
          endif
          if(errc.eq.GFC_SUCCESS) then
           call this%release_lock(errc)
          else
           call this%release_lock()
          endif
         else
          errc=GFC_IN_USE
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DictElemConstruct
!---------------------------------------------------------------------
        subroutine DictElemDestruct(this,ierr,key_assoc,dtor_f,locked)
!Destructs the key-value pair inside the dictionary element.
!<dtor_f> provides an optional explicit destructor for the dictionary
!element value, if needed. Alternatively, the value may have the final
!subroutine defined. In contrast, the dictionary element key must
!have the final subroutine defined if it requires a non-trivial destruction.
         implicit none
         class(dict_elem_t), intent(inout):: this     !inout: element of a dictionary
         integer(INTD), intent(out), optional:: ierr  !out: error code
         logical, intent(in), optional:: key_assoc    !in: if TRUE, the key will be assumed stored by reference
         procedure(gfc_destruct_i), optional:: dtor_f !in: explicit destructor for the value
         logical, intent(in), optional:: locked       !in: if TRUE, the dictionary element will be assumed already locked (defaults to FALSE)
         integer(INTD):: errc
         logical:: assk,lckd,lck
         integer:: ier

         errc=GFC_SUCCESS
         if(present(key_assoc)) then; assk=key_assoc; else; assk=.FALSE.; endif
         if(present(locked)) then; lckd=locked; else; lckd=.FALSE.; endif; lck=lckd
         if(.not.lck) lck=(this%in_use(errc,set_lock=.TRUE.,report_refs=.FALSE.).eq.GFC_FALSE)
         if(lck) then
          if(errc.eq.GFC_SUCCESS) then
           if(associated(this%key)) then
            if(.not.assk) then
             if(present(dtor_f)) then
              call this%gfc_cont_elem_t%destruct(errc,dtor_f=dtor_f,locked=.TRUE.)
             else
              call this%gfc_cont_elem_t%destruct(errc,locked=.TRUE.)
             endif
             deallocate(this%key,STAT=ier)
             if(ier.ne.0.and.errc.eq.GFC_SUCCESS) errc=GFC_MEM_FREE_FAILED
            endif
            this%key=>NULL()
           else
            errc=GFC_CORRUPTED_CONT
           endif
          endif
          if(.not.lckd) then
           if(errc.eq.GFC_SUCCESS) then
            call this%release_lock(errc)
           else
            call this%release_lock()
           endif
          endif
         else
          if(errc.eq.GFC_SUCCESS) errc=GFC_IN_USE
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DictElemDestruct
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
!-----------------------------------------------------------------------
        function DictElemPredicateKey(this,predicat_f,ierr) result(pred)
!Evaluates a user-defined predicate on the key of a given dictionary element.
         implicit none
         integer(INTD):: pred                        !out: evaluated predicate value {GFC_TRUE,GFC_FALSE,GFC_ERROR}
         class(dict_elem_t), intent(in):: this       !in: element of a dictionary
         procedure(gfc_predicate_i):: predicat_f     !in: user-defined predicate
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS; pred=GFC_ERROR
         if(.not.this%is_empty()) then
          pred=predicat_f(this%key); if(pred.eq.GFC_ERROR) errc=GFC_ERROR
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end function DictElemPredicateKey
!-------------------------------------------------------------------------
        function DictElemCompareKey(this,other_key,cmp_f,ierr) result(cmp)
!Compares the dictionary element key (on the left) with another key (on the right).
         implicit none
         integer(INTD):: cmp                         !out: result of the comparison (see GFC_CMP_XXX in gfc_base.F90)
         class(dict_elem_t), intent(in):: this       !in: element of a dictionary whose key is being compared
         class(*), intent(in):: other_key            !in: the other value to be compared with
         procedure(gfc_cmp_i):: cmp_f                !in: user-defined comparison function
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS; cmp=GFC_CMP_NA
         if(.not.this%is_empty()) then
          cmp=cmp_f(this%key,other_key); if(cmp.eq.GFC_ERROR) errc=GFC_ERROR
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end function DictElemCompareKey
!------------------------------------------------------------
        subroutine DictElemPrintKey(this,print_f,ierr,dev_id)
!Prints the key of a dictionary element using a user-defined print function.
         implicit none
         class(dict_elem_t), intent(in):: this        !in: element of a dictionary
         procedure(gfc_print_i):: print_f             !in: user-defined printing function
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD), intent(in), optional:: dev_id !in: output device (default to screen, 6)
         integer(INTD):: dev,errc

         errc=GFC_SUCCESS
         if(.not.this%is_empty()) then
          if(present(dev_id)) then; dev=dev_id; else; dev=6; endif !defaults to screen
          write(dev,'("#GFC container element key:")')
          errc=print_f(this%key,dev); if(errc.ne.0) errc=GFC_ACTION_FAILED
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DictElemPrintKey
!---------------------------------------------------------------------------
        subroutine DictElemPrintIt(this,print_key_f,print_val_f,ierr,dev_id)
!Prints the dictionary element (key,value).
         implicit none
         class(dict_elem_t), intent(in):: this        !in: dictionary element
         procedure(gfc_print_i):: print_key_f         !in: key printing function
         procedure(gfc_print_i):: print_val_f         !in: value printing function
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD), intent(in), optional:: dev_id !in: output device (defaults to screen, 6)
         integer(INTD):: dev,errc

         errc=GFC_SUCCESS
         if(.not.this%is_empty()) then
          if(present(dev_id)) then; dev=dev_id; else; dev=6; endif !defaults to screen
          write(dev,'("#GFC dictionary element (key,value):")')
          call this%print_key(print_key_f,errc,dev)
          if(errc.eq.GFC_SUCCESS) call this%print_value(print_val_f,errc,dev)
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DictElemPrintIt
![dictionary_t]=====================================
        function DictionaryIsEmpty(this) result(res)
!Returns GFC_TRUE if the dictionary is empty, GFC_FALSE otherwise (or error code).
         implicit none
         integer(INTD):: res                    !out: result of query
         class(dictionary_t), intent(in):: this !in: dictionary

         if(associated(this%root)) then
          res=GFC_FALSE
         else
          res=GFC_TRUE
         endif
         return
        end function DictionaryIsEmpty
!----------------------------------------------------------------
        function DictionaryIsSubdictionary(this,ierr) result(res)
!Returns TRUE if the dictionary is a subdictionary of a larger dictionary.
         implicit none
         logical:: res                               !out: result
         class(dictionary_t), intent(in):: this      !in: dictionary
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS; res=.FALSE.
         if(associated(this%root)) then
          res=associated(this%root%parent)
         else
          errc=GFC_EMPTY_CONT
         endif
         if(present(ierr)) ierr=errc
         return
        end function DictionaryIsSubdictionary
!-----------------------------------------------------------
        subroutine DictionarySetKeyStorage(this,policy,ierr)
!Sets the key storage policy.
         implicit none
         class(dictionary_t), intent(inout):: this   !inout: dictionary (must be empty)
         logical, intent(in):: policy                !in: storage policy: {GFC_BY_VAL,GFC_BY_REF}
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
         if(.not.associated(this%root)) then
          this%key_storage=policy
         else
          errc=GFC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DictionarySetKeyStorage
!-------------------------------------------------
        subroutine DictionaryReroot(this,new_root)
!Changes the root of the dictionary.
         implicit none
         class(dictionary_t), intent(inout):: this             !inout: dictionary
         class(dict_elem_t), pointer, intent(inout):: new_root !in: pointer to the new root or NULL()

         if(associated(this%root)) call this%root%decr_ref_()
         this%root=>new_root
         if(associated(this%root)) call this%root%incr_ref_()
         return
        end subroutine DictionaryReroot
![dictionary_iter_t]================================
        subroutine DictionaryIterJump(this,new_elem)
!Moves the iterator to an arbitrary specified position.
         implicit none
         class(dictionary_iter_t), intent(inout):: this        !inout: dictionary iterator
         class(dict_elem_t), pointer, intent(inout):: new_elem !in: pointer to the new element or NULL()
         integer(INTD):: errc,sts

         if(associated(this%current)) call this%current%decr_ref_()
         this%current=>new_elem
         if(associated(this%current)) then
          call this%current%incr_ref_()
          errc=this%set_status_(GFC_IT_ACTIVE)
         else
          if(associated(this%container%root)) then
           errc=this%set_status_(GFC_IT_DONE)
          else
           errc=this%set_status_(GFC_IT_EMPTY)
          endif
         endif
         return
        end subroutine DictionaryIterJump
!----------------------------------------------------------
        function DictionaryIterInit(this,cont) result(ierr)
!Initializes an iterator and resets it to the beginning of the container.
         implicit none
         integer(INTD):: ierr                              !out: error code
         class(dictionary_iter_t), intent(inout):: this    !inout: iterator
         class(gfc_container_t), target, intent(in):: cont !in: container

         ierr=GFC_SUCCESS
         select type(cont)
         class is (dictionary_t)
          this%container=>cont
          ierr=this%reset()
         class default
          ierr=GFC_INVALID_ARGS
         end select
         return
        end function DictionaryIterInit
!------------------------------------------------------
        function DictionaryIterReset(this) result(ierr)
!Resets the iterator to the beginning (root element).
         implicit none
         integer(INTD):: ierr                           !out: error code
         class(dictionary_iter_t), intent(inout):: this !inout: iterator

         ierr=GFC_SUCCESS
         if(associated(this%container)) then
          if(associated(this%current)) call this%current%decr_ref_()
          this%current=>this%container%root
          if(associated(this%current)) then
           call this%current%incr_ref_()
           ierr=this%set_status_(GFC_IT_ACTIVE) !non-empty iterator/container
          else
           ierr=this%set_status_(GFC_IT_EMPTY) !empty iterator/container
          endif
          call this%reset_count() !reset all iteration counters
         else
          ierr=this%set_status_(GFC_IT_NULL)
          ierr=GFC_IT_NULL
         endif
         return
        end function DictionaryIterReset
!--------------------------------------------------------
        function DictionaryIterRelease(this) result(ierr)
!Dissociates the iterator from its container.
         implicit none
         integer(INTD):: ierr                           !out: error code
         class(dictionary_iter_t), intent(inout):: this !inout: iterator

         if(associated(this%current)) call this%current%decr_ref_()
         this%current=>NULL(); this%container=>NULL()
         call this%reset_count(); ierr=this%set_status_(GFC_IT_NULL)
         return
        end function DictionaryIterRelease
!--------------------------------------------------------------
        function DictionaryIterPointee(this,ierr) result(pntee)
!Returns the container element the iterator is currently pointing to.
         implicit none
         class(gfc_cont_elem_t), pointer:: pntee     !out: container element currently pointed to by the iterator
         class(dictionary_iter_t), intent(in):: this !in: iterator
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=this%get_status()
         if(errc.eq.GFC_IT_ACTIVE) then
          pntee=>this%current; errc=GFC_SUCCESS
         else
          pntee=>NULL()
         endif
         if(present(ierr)) ierr=errc
         return
        end function DictionaryIterPointee
!------------------------------------------------------------
        function DictionaryIterNext(this,elem_p) result(ierr)
!Moves the iterator position to the next element (in some sense).
!If <elem_p> is absent, the iterator moves to the next element, if any.
!If <elem_p> is present, the iterator simply returns the next element in <elem_p> without moving.
!Complexity: O(1)...O(logN). No additional memory is used.
         implicit none
         integer(INTD):: ierr                                            !out: error code
         class(dictionary_iter_t), intent(inout):: this                  !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element
         class(dict_elem_t), pointer:: dp

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           ierr=GFC_SUCCESS
           dp=>this%current%child_lt
           if(.not.associated(dp)) then
            dp=>this%current%child_gt
            if(.not.associated(dp)) then
             dp=>this%current
             do while(associated(dp))
              if(.not.associated(dp,this%container%root)) then !root of a subdictionary may have a parent
               if(associated(dp,dp%parent%child_lt)) then !moving up from the left
                if(associated(dp%parent%child_gt)) then
                 dp=>dp%parent%child_gt; exit
                endif
               endif
               dp=>dp%parent
              else
               dp=>NULL()
              endif
             enddo
            endif
           endif
           if(present(elem_p)) then
            elem_p=>dp
           else
            call this%current%decr_ref_()
            this%current=>dp
            if(associated(this%current)) then
             call this%current%incr_ref_()
            else
             ierr=this%set_status_(GFC_IT_DONE)
            endif
           endif
           if(.not.associated(dp)) ierr=GFC_IT_DONE
           dp=>NULL()
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function DictionaryIterNext
!----------------------------------------------------------------
        function DictionaryIterPrevious(this,elem_p) result(ierr)
!Moves the iterator position to the previous element (in some sense).
!If <elem_p> is absent, the iterator moves to the previous element, if any.
!If <elem_p> is present, the iterator simply returns the previous element in <elem_p> without moving.
!Complexity: O(1)...O(logN). No additional memory is used.
         implicit none
         integer(INTD):: ierr                                            !out: error code
         class(dictionary_iter_t), intent(inout):: this                  !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element
         class(dict_elem_t), pointer:: dp

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           ierr=GFC_SUCCESS
           if(associated(this%current,this%container%root)) then
            dp=>NULL()
           else
            dp=>this%current%parent
            if(associated(this%current,dp%child_gt).and.associated(dp%child_lt)) then
             dp=>dp%child_lt
             do
              if(associated(dp%child_gt)) then
               dp=>dp%child_gt
              else
               if(associated(dp%child_lt)) then
                dp=>dp%child_lt
               else
                exit
               endif
              endif
             enddo
            endif
           endif
           if(present(elem_p)) then
            elem_p=>dp
           else
            call this%current%decr_ref_()
            this%current=>dp
            if(associated(this%current)) then
             call this%current%incr_ref_()
            else
             ierr=this%set_status_(GFC_IT_DONE)
            endif
           endif
           if(.not.associated(dp)) ierr=GFC_IT_DONE
           dp=>NULL()
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function DictionaryIterPrevious
!-------------------------------------------------------------------
        function DictionaryIterNextInOrder(this,elem_p) result(ierr)
!Moves the iterator position to the next-in-order element.
!If <elem_p> is absent, the iterator moves to the next-in-order element, if any.
!If <elem_p> is present, the iterator simply returns the next-in-order element in <elem_p> without moving.
!Complexity: O(1)...O(logN). No additional memory is used.
         implicit none
         integer(INTD):: ierr                                            !out: error code
         class(dictionary_iter_t), intent(inout):: this                  !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element
         class(dict_elem_t), pointer:: dp

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           ierr=GFC_SUCCESS
           if(associated(this%current%child_gt)) then
            dp=>this%current%child_gt
            do while(associated(dp%child_lt)); dp=>dp%child_lt; enddo
           else
            ierr=GFC_IT_DONE; dp=>this%current
            do while(associated(dp%parent).and.(.not.associated(dp,this%container%root)))
             if(associated(dp,dp%parent%child_lt)) then
              dp=>dp%parent; ierr=GFC_SUCCESS; exit
             else
              dp=>dp%parent
             endif
            enddo
            if(ierr.eq.GFC_IT_DONE) dp=>NULL()
           endif
           if(present(elem_p)) then
            elem_p=>dp
           else
            call this%current%decr_ref_()
            this%current=>dp
            if(associated(this%current)) then
             call this%current%incr_ref_()
            else
             ierr=this%set_status_(GFC_IT_DONE)
            endif
           endif
           if(.not.associated(dp)) ierr=GFC_IT_DONE
           dp=>NULL()
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function DictionaryIterNextInOrder
!-------------------------------------------------------------------
        function DictionaryIterPrevInOrder(this,elem_p) result(ierr)
!Moves the iterator position to the previous-in-order element.
!If <elem_p> is absent, the iterator moves to the previous-in-order element, if any.
!If <elem_p> is present, the iterator simply returns the previous-in-order element in <elem_p> without moving.
!Complexity: O(1)...O(logN). No additional memory is used.
         implicit none
         integer(INTD):: ierr                                            !out: error code
         class(dictionary_iter_t), intent(inout):: this                  !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element
         class(dict_elem_t), pointer:: dp

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           ierr=GFC_SUCCESS
           if(associated(this%current%child_lt)) then
            dp=>this%current%child_lt
            do while(associated(dp%child_gt)); dp=>dp%child_gt; enddo
           else
            ierr=GFC_IT_DONE; dp=>this%current
            do while(associated(dp%parent).and.(.not.associated(dp,this%container%root)))
             if(associated(dp,dp%parent%child_gt)) then
              dp=>dp%parent; ierr=GFC_SUCCESS; exit
             else
              dp=>dp%parent
             endif
            enddo
            if(ierr.eq.GFC_IT_DONE) dp=>NULL()
           endif
           if(present(elem_p)) then
            elem_p=>dp
           else
            call this%current%decr_ref_()
            this%current=>dp
            if(associated(this%current)) then
             call this%current%incr_ref_()
            else
             ierr=this%set_status_(GFC_IT_DONE)
            endif
           endif
           if(.not.associated(dp)) ierr=GFC_IT_DONE
           dp=>NULL()
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function DictionaryIterPrevInOrder
!-----------------------------------------------------------------------------
        function DictionaryIterMoveInOrder(this,direction,elem_p) result(ierr)
!Moves the iterator position either to the next-in-order or previous-in-order element.
!If <elem_p> is absent, the iterator moves to the corresponding in-order element, if any.
!If <elem_p> is present, the iterator simply returns the corresponding in-order element
!in <elem_p> without moving. Complexity: O(1)...O(logN). No additional memory is used.
         implicit none
         integer(INTD):: ierr                                            !out: error code
         class(dictionary_iter_t), intent(inout):: this                  !inout: iterator
         logical, intent(in):: direction                                 !direction: {GFC_DICT_SUCCESSOR,GFC_DICT_PREDECESSOR}
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element

         if(present(elem_p)) then
          if(direction.eqv.GFC_DICT_SUCCESSOR) then
           ierr=this%next_in_order(elem_p)
          else
           ierr=this%prev_in_order(elem_p)
          endif
         else
          if(direction.eqv.GFC_DICT_SUCCESSOR) then
           ierr=this%next_in_order()
          else
           ierr=this%prev_in_order()
          endif
         endif
         return
        end function DictionaryIterMoveInOrder
!-----------------------------------------------------------------
        function DictionaryIterMoveToMin(this,elem_p) result(ierr)
!Moves the iterator position to the minimal element.
!If <elem_p> is absent, the iterator moves to the minimal element.
!If <elem_p> is present, the iterator simply returns the minimal element in <elem_p> without moving.
!Complexity: O(logN). No additional memory is used.
         implicit none
         integer(INTD):: ierr                                            !out: error code
         class(dictionary_iter_t), intent(inout):: this                  !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element
         class(dict_elem_t), pointer:: dp

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%container%root)) then
           ierr=GFC_SUCCESS; dp=>this%container%root
           do while(associated(dp%child_lt)); dp=>dp%child_lt; enddo
           if(present(elem_p)) then
            elem_p=>dp
           else
            call this%current%decr_ref_()
            this%current=>dp
            call this%current%incr_ref_()
           endif
           dp=>NULL()
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function DictionaryIterMoveToMin
!-----------------------------------------------------------------
        function DictionaryIterMoveToMax(this,elem_p) result(ierr)
!Moves the iterator position to the maximal element.
!If <elem_p> is absent, the iterator moves to the maximal element.
!If <elem_p> is present, the iterator simply returns the maximal element in <elem_p> without moving.
!Complexity: O(logN). No additional memory is used.
         implicit none
         integer(INTD):: ierr                                            !out: error code
         class(dictionary_iter_t), intent(inout):: this                  !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element
         class(dict_elem_t), pointer:: dp

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%container%root)) then
           ierr=GFC_SUCCESS; dp=>this%container%root
           do while(associated(dp%child_gt)); dp=>dp%child_gt; enddo
           if(present(elem_p)) then
            elem_p=>dp
           else
            call this%current%decr_ref_()
            this%current=>dp
            call this%current%incr_ref_()
           endif
           dp=>NULL()
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function DictionaryIterMoveToMax
!--------------------------------------------------------------
        function DictionaryIterMoveUp(this,elem_p) result(ierr)
!Moves the iterator position up the binary search tree (to the parental element).
!If <elem_p> is absent, the iterator moves up.
!If <elem_p> is present, the iterator simply returns the corresponding element in <elem_p> without moving.
!Complexity: O(1). No additional memory is used.
         implicit none
         integer(INTD):: ierr                                            !out: error code
         class(dictionary_iter_t), intent(inout):: this                  !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element
         class(dict_elem_t), pointer:: dp

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           ierr=GFC_SUCCESS; dp=>NULL()
           if(associated(this%current%parent)) then
            dp=>this%current%parent
           else
            ierr=GFC_NO_MOVE
           endif
           if(present(elem_p)) then
            elem_p=>dp
           else
            if(ierr.eq.GFC_SUCCESS) then
             call this%current%decr_ref_()
             this%current=>dp
             call this%current%incr_ref_()
            endif
           endif
           dp=>NULL()
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function DictionaryIterMoveUp
!--------------------------------------------------------------------------
        function DictionaryIterMoveDown(this,elem_p,direction) result(ierr)
!Moves the iterator position down the binary search tree, either left or right.
!If <elem_p> is absent, the iterator moves down.
!If <elem_p> is present, the iterator simply returns the corresponding element in <elem_p> without moving.
!Complexity: O(1). No additional memory is used.
         implicit none
         integer(INTD):: ierr                                            !out: error code
         class(dictionary_iter_t), intent(inout):: this                  !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element
         logical, intent(in), optional:: direction                       !in: direction: {GFC_DICT_LEFT,GFC_DICT_RIGHT}
         class(dict_elem_t), pointer:: dp

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           ierr=GFC_SUCCESS; dp=>NULL()
           if(present(direction)) then
            if(direction.eqv.GFC_DICT_LEFT) then !move down left
             if(associated(this%current%child_lt)) then
              dp=>this%current%child_lt
             else
              ierr=GFC_NO_MOVE
             endif
            else !move down right
             if(associated(this%current%child_gt)) then
              dp=>this%current%child_gt
             else
              ierr=GFC_NO_MOVE
             endif
            endif
           else !move left, if not left, move right
            if(associated(this%current%child_lt)) then
             dp=>this%current%child_lt
            else
             if(associated(this%current%child_gt)) then
              dp=>this%current%child_gt
             else
              ierr=GFC_NO_MOVE
             endif
            endif
           endif
           if(present(elem_p)) then
            elem_p=>dp
           else
            if(ierr.eq.GFC_SUCCESS) then
             call this%current%decr_ref_()
             this%current=>dp
             call this%current%incr_ref_()
            endif
           endif
           dp=>NULL()
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function DictionaryIterMoveDown
!-------------------------------------------------------------
        function DictionaryIterGetKey(this,ierr) result(key_p)
!Returns a pointer to the key in the current iterator position.
         implicit none
         class(*), pointer:: key_p                   !out: pointer to the current position key
         class(dictionary_iter_t), intent(in):: this !in: dictionary iterator
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         key_p=>NULL(); errc=this%get_status()
         if(errc.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           key_p=>this%current%get_key(errc)
          else
           errc=GFC_ERROR
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end function DictionaryIterGetKey
!-----------------------------------------------------------
        subroutine DictionaryIterDeleteAll(this,ierr,dtor_f)
!Deletes all dictionary elements, leaving dictionary empty at the end.
!Complexity: O(NLogN)
         implicit none
         class(dictionary_iter_t), intent(inout):: this !inout: dictionary iterator
         integer(INTD), intent(out), optional:: ierr    !out: error code
         procedure(gfc_destruct_i), optional:: dtor_f   !in: explicit destructor for the dictionary element value
         integer(INTD):: errc,errcode
         class(dict_elem_t), pointer:: dp
         integer:: ier

         errc=this%get_status()
         if(errc.eq.GFC_IT_ACTIVE) then
          errc=this%reset()
          if(errc.eq.GFC_SUCCESS) then
           dloop: do
            do while(errc.eq.GFC_SUCCESS); errc=this%move_down(); enddo
            if(errc.eq.GFC_NO_MOVE.and.associated(this%current)) then
             if(present(dtor_f)) then
              call this%current%destruct_keyval(errc,key_assoc=this%container%key_storage,dtor_f=dtor_f)
             else
              call this%current%destruct_keyval(errc,key_assoc=this%container%key_storage)
             endif
             if(errc.eq.GFC_SUCCESS) then
              if(associated(this%current,this%container%root)) then !last element (root)
               if(associated(this%current%parent)) then
                if(associated(this%current,this%current%parent%child_lt)) then
                 this%current%parent%child_lt=>NULL()
                elseif(associated(this%current,this%current%parent%child_gt)) then
                 this%current%parent%child_gt=>NULL()
                else
                 errc=GFC_CORRUPTED_CONT; exit dloop
                endif
               endif
               call this%current%decr_ref_()
               deallocate(this%current,STAT=ier)
               if(ier.ne.0) errc=GFC_MEM_FREE_FAILED
               this%current=>NULL(); errcode=this%set_status_(GFC_IT_EMPTY)
               if(errc.eq.GFC_SUCCESS) errc=errcode
               exit dloop
              else
               dp=>this%current
               call this%current%decr_ref_()
               this%current=>this%current%parent
               call this%current%incr_ref_()
               if(associated(this%current%child_lt,dp)) then
                this%current%child_lt=>NULL()
               elseif(associated(this%current%child_gt,dp)) then
                this%current%child_gt=>NULL()
               else
                errc=GFC_CORRUPTED_CONT; exit dloop
               endif
               dp%parent=>NULL()
               deallocate(dp,STAT=ier)
               if(ier.ne.0) then; errc=GFC_MEM_FREE_FAILED; exit dloop; endif
              endif
             else
              exit dloop
             endif
            else
             errc=GFC_ERROR; exit dloop
            endif
           enddo dloop
          else
           errc=GFC_ERROR
          endif
         elseif(errc.eq.GFC_IT_EMPTY) then
          errc=GFC_EMPTY_CONT
         else
          errc=GFC_NULL_CONT
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DictionaryIterDeleteAll
!-----------------------------------------------------------------------------------------------------------------------
#ifdef NO_GNU
        function DictionaryIterSearch(this,action,cmp_key_f,key,value_in,store_by,value_out,copy_ctor_val_f,dtor_val_f)& !`GCC bug
        &result(dict_search)
#else
        function DictionaryIterSearch(this,action,cmp_key_f,key,value_in,store_by,value_out,dtor_val_f)&
        &result(dict_search)
#endif
!Looks up a given key in the dictionary with optional actions. If the key is found (and not deleted)
!or newly added, then the current iterator position will point to that item, otherwise it will not change.
!The incoming iterator must have been initialized (either ACTIVE or EMPTY).
         implicit none
         integer(INTD):: dict_search                            !out: result: {GFC_FOUND,GFC_NOT_FOUND,specific errors}
         class(dictionary_iter_t), intent(inout):: this         !inout: dictionary iterator
         integer, intent(in):: action                           !in: requested action (see action parameters at the top of this module)
         procedure(gfc_cmp_i):: cmp_key_f                       !in: key comparison function, returns: {GFC_CMP_LT,GFC_CMP_GT,GFC_CMP_EQ,GFC_CMP_ERR}
         class(*), intent(in), target:: key                     !in: key being searched for
         class(*), intent(in), target, optional:: value_in      !in: an optional value to be stored with the key
         logical, intent(in), optional:: store_by               !in: storage type for newly added values: {GFC_BY_VAL,GFC_BY_REF}, defaults to GFC_BY_VAL
         class(*), pointer, intent(out), optional:: value_out   !out: when fetching, this will point to the value found by the key (NULL otherwise)
#ifdef NO_GNU
         procedure(gfc_copy_i), optional:: copy_ctor_val_f      !in: explicit copy constructor for the element value, if needed
#endif
         procedure(gfc_destruct_i), optional:: dtor_val_f       !in: explicit destructor for the dictionary element value, if needed
         class(dict_elem_t), pointer:: curr,old_cdp,leave,term,nullptr
         class(dictionary_t), pointer:: dict
         integer(INTD):: i,j,act,lev_p,grow,ierr

         dict_search=this%get_status()
         if(dict_search.eq.GFC_IT_DONE) then
          dict_search=this%reset(); if(dict_search.ne.GFC_SUCCESS) return
          dict_search=this%get_status()
         endif
         if(dict_search.ne.GFC_IT_ACTIVE.and.dict_search.ne.GFC_IT_EMPTY) return
         dict_search=GFC_NOT_FOUND; if(present(value_out)) value_out=>NULL(); nullptr=>NULL()
         if(associated(this%container)) then; dict=>this%container; else; dict_search=GFC_CORRUPTED_CONT; return; endif
!Look up the key:
         if(associated(dict%root)) then
          curr=>dict%root; lev_p=0
          sloop: do
           i=cmp_key_f(key,curr%key)
           if(i.eq.GFC_CMP_LT) then
            if(associated(curr%child_lt)) then; curr=>curr%child_lt; lev_p=lev_p+1; else; exit sloop; endif
           elseif(i.eq.GFC_CMP_GT) then
            if(associated(curr%child_gt)) then; curr=>curr%child_gt; lev_p=lev_p+1; else; exit sloop; endif
           elseif(i.eq.GFC_CMP_EQ) then
            dict_search=GFC_FOUND; exit sloop
           else
            if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): key comparison failed: error ",i11)') i
            dict_search=GFC_CMP_ERR; curr=>NULL(); return
           endif
          enddo sloop
         else
          i=0; curr=>NULL(); lev_p=-1 !empty dictionary
         endif
!Action:
         if(action.eq.GFC_DICT_ADD_OR_MODIFY) then
          if(dict_search.eq.GFC_NOT_FOUND) then
           act=GFC_DICT_ADD_IF_NOT_FOUND
          elseif(dict_search.eq.GFC_FOUND) then
           act=GFC_DICT_REPLACE_IF_FOUND
          endif
         else
          act=action
         endif
         select case(act) !process the action
         case(GFC_DICT_JUST_FIND) !no action
          if(dict_search.eq.GFC_FOUND) call this%jump_(curr)
         case(GFC_DICT_FETCH_IF_FOUND) !return the pointer to the stored <value> if found
          if(dict_search.eq.GFC_FOUND) then
           if(present(value_out)) then
            value_out=>curr%get_value(ierr)
            if(ierr.eq.GFC_SUCCESS) then; call this%jump_(curr); else; dict_search=ierr; endif
           else
            if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): fetch: absent value return pointer!")')
            dict_search=GFC_INVALID_ARGS; curr=>NULL(); return
           endif
          endif
         case(GFC_DICT_REPLACE_IF_FOUND) !replace the stored <value> if found
          if(dict_search.eq.GFC_FOUND) then
           if(present(dtor_val_f)) then
            call curr%gfc_cont_elem_t%destruct(ierr,dtor_f=dtor_val_f)
           else
            call curr%gfc_cont_elem_t%destruct(ierr)
           endif
           if(ierr.ne.GFC_SUCCESS) then
            if(VERBOSE)&
            &write(CONS_OUT,'("#ERROR(gfc::dictionary::search): replace: dictionary element value destruction failed!")')
            dict_search=GFC_MEM_FREE_FAILED; curr=>NULL(); return
           endif
           if(present(value_in)) then
            if(present(store_by)) then
             call curr%gfc_cont_elem_t%construct_base(value_in,ierr,assoc_only=store_by)
            else
             call curr%gfc_cont_elem_t%construct_base(value_in,ierr)
            endif
            if(ierr.ne.GFC_SUCCESS) then
             if(VERBOSE)&
             &write(CONS_OUT,'("#ERROR(gfc::dictionary::search): replace: dictionary element value construction failed!")')
             dict_search=GFC_MEM_ALLOC_FAILED; curr=>NULL(); return
            endif
           else
            if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): replace: absent input value!")')
            dict_search=GFC_INVALID_ARGS; curr=>NULL(); return
           endif
           if(present(value_out)) then; value_out=>curr%get_value(ierr); else; ierr=GFC_SUCCESS; endif
           if(ierr.eq.GFC_SUCCESS) then; call this%jump_(curr); else; dict_search=ierr; endif
          endif
         case(GFC_DICT_DELETE_IF_FOUND) !delete the item if found
          if(dict_search.eq.GFC_FOUND) then
           if(associated(this%current)) then
            if(associated(this%current,curr)) then; old_cdp=>NULL(); else; old_cdp=>this%current; endif
           else
            old_cdp=>NULL()
           endif
           call this%jump_(curr)
           if(associated(curr%child_lt).and.associated(curr%child_gt)) then !both subtrees are present
            if(curr%balance_fac.le.0) then !right subtree is taller or equal
             grow=-1; j=this%move_in_order(GFC_DICT_SUCCESSOR) !find in-order successor
             if(j.eq.GFC_SUCCESS) then
              if(associated(this%current%child_gt)) then
               leave=>this%current%child_gt
              else
               leave=>NULL()
              endif
             else
              if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): delete: next in-order not found!")')
              dict_search=GFC_CORRUPTED_CONT
             endif
            else !left subtree is taller
             grow=+1; j=this%move_in_order(GFC_DICT_PREDECESSOR) !find in-order predecessor
             if(j.eq.GFC_SUCCESS) then
              if(associated(this%current%child_lt)) then
               leave=>this%current%child_lt
              else
               leave=>NULL()
              endif
             else
              if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): delete: previous in-order not found!")')
              dict_search=GFC_CORRUPTED_CONT
             endif
            endif
            if(dict_search.eq.GFC_FOUND) then
             if(associated(this%current%parent,curr)) then
              term=>NULL()
             else
              term=>this%current%parent
             endif
             if(associated(curr%parent)) then
              if(associated(curr%parent%child_lt,curr)) then
               curr%parent%child_lt=>this%current
              elseif(associated(curr%parent%child_gt,curr)) then
               curr%parent%child_gt=>this%current
              else
               if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): delete: lost parent!")')
               dict_search=GFC_CORRUPTED_CONT
              endif
             else
              if(associated(dict%root,curr)) then
               call dict%reroot_(this%current)
              else
               if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): delete: lost root!")')
               dict_search=GFC_CORRUPTED_CONT
              endif
             endif
             if(dict_search.eq.GFC_FOUND) then
              if(grow.gt.0) then !reducing the left subtree
               if(associated(term)) then
                curr%child_lt%parent=>this%current
                this%current%child_lt=>curr%child_lt
                if(associated(leave)) then
                 term%child_gt=>leave; leave%parent=>term
                else
                 term%child_gt=>NULL()
                endif
               endif
               curr%child_gt%parent=>this%current
               this%current%child_gt=>curr%child_gt
              elseif(grow.lt.0) then !reducing the right subtree
               if(associated(term)) then
                curr%child_gt%parent=>this%current
                this%current%child_gt=>curr%child_gt
                if(associated(leave)) then
                 term%child_lt=>leave; leave%parent=>term
                else
                 term%child_lt=>NULL()
                endif
               endif
               curr%child_lt%parent=>this%current
               this%current%child_lt=>curr%child_lt
              endif
              if(associated(curr%parent)) then
               this%current%parent=>curr%parent
              else
               this%current%parent=>NULL()
              endif
              this%current%balance_fac=curr%balance_fac
              if(associated(term)) then
               call this%jump_(term); grow=-grow; term=>NULL()
              endif
              if(associated(leave)) leave=>NULL()
             endif
            endif
           else !at least one subtree is absent
            if(associated(curr%child_lt)) then !left subtree is present (a leave)
             if(associated(curr%parent)) then
              if(associated(curr%parent%child_lt,curr)) then
               curr%parent%child_lt=>curr%child_lt; call this%jump_(curr%parent); grow=+1
              elseif(associated(curr%parent%child_gt,curr)) then
               curr%parent%child_gt=>curr%child_lt; call this%jump_(curr%parent); grow=-1
              else
               if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): delete: lost parent!")')
               dict_search=GFC_CORRUPTED_CONT
              endif
              curr%child_lt%parent=>curr%parent
             else
              if(associated(dict%root,curr)) then
               call dict%reroot_(curr%child_lt)
               dict%root%parent=>NULL(); grow=0
              else
               if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): delete: lost root!")')
               dict_search=GFc_CORRUPTED_CONT
              endif
             endif
            elseif(associated(curr%child_gt)) then !right subtree is present (a leave)
             if(associated(curr%parent)) then
              if(associated(curr%parent%child_lt,curr)) then
               curr%parent%child_lt=>curr%child_gt; call this%jump_(curr%parent); grow=+1
              elseif(associated(curr%parent%child_gt,curr)) then
               curr%parent%child_gt=>curr%child_gt; call this%jump_(curr%parent); grow=-1
              else
               if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): delete: lost parent!")')
               dict_search=GFC_CORRUPTED_CONT
              endif
              curr%child_gt%parent=>curr%parent
             else
              if(associated(dict%root,curr)) then
               call dict%reroot_(curr%child_gt)
               dict%root%parent=>NULL(); grow=0
              else
               if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): delete: lost root!")')
               dict_search=GFC_CORRUPTED_CONT
              endif
             endif
            else !both subtrees are absent
             if(associated(curr%parent)) then
              if(associated(curr%parent%child_lt,curr)) then
               curr%parent%child_lt=>NULL(); call this%jump_(curr%parent); grow=+1
              elseif(associated(curr%parent%child_gt,curr)) then
               curr%parent%child_gt=>NULL(); call this%jump_(curr%parent); grow=-1
              else
               if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): delete: lost parent!")')
               dict_search=GFC_CORRUPTED_CONT
              endif
             else
              if(associated(dict%root,curr)) then
               call dict%reroot_(nullptr)
               call this%jump_(nullptr); grow=0
              else
               if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): delete: lost root!")')
               dict_search=GFC_CORRUPTED_CONT
              endif
             endif
            endif
           endif
           if(present(dtor_val_f)) then
            call curr%destruct_keyval(j,key_assoc=dict%key_storage,dtor_f=dtor_val_f)
           else
            call curr%destruct_keyval(j,key_assoc=dict%key_storage)
           endif
           if(j.ne.GFC_SUCCESS) then
            if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): delete: dictionary element destruction failed!")')
            dict_search=GFC_MEM_FREE_FAILED
           endif
           deallocate(curr,STAT=j)
           if(j.ne.0) then
            if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): delete: dictionary element deallocation failed!")')
            dict_search=GFC_MEM_FREE_FAILED
           endif
           if(dict_search.eq.GFC_FOUND.and.grow.ne.0) then
            ierr=-1; call rebalance(this%current,grow,ierr)
            if(ierr.ne.0) then
             if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): delete: rebalancing failed: ",i11)') ierr
             dict_search=GFC_CORRUPTED_CONT
            endif
           endif
           if(associated(old_cdp)) then
            call this%jump_(old_cdp); old_cdp=>NULL()
           else
            call this%jump_(nullptr)
           endif
          endif
         case(GFC_DICT_ADD_IF_NOT_FOUND) !add a new item if the key is not found
          if(dict_search.eq.GFC_NOT_FOUND) then
           if(lev_p.ge.0) then
            if(i.eq.GFC_CMP_LT) then
             allocate(curr%child_lt,STAT=j)
             if(j.ne.0) then
              if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): add: dictionary element allocation failed!")')
              dict_search=GFC_MEM_ALLOC_FAILED; curr=>NULL(); return
             endif
             curr%child_lt%parent=>curr
             grow=+1; curr=>curr%child_lt
            elseif(i.eq.GFC_CMP_GT) then
             allocate(curr%child_gt,STAT=j)
             if(j.ne.0) then
              if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): add: dictionary element allocation failed!")')
              dict_search=GFC_MEM_ALLOC_FAILED; curr=>NULL(); return
             endif
             curr%child_gt%parent=>curr
             grow=-1; curr=>curr%child_gt
            endif
            curr%child_lt=>NULL(); curr%child_gt=>NULL(); curr%balance_fac=0
            ierr=+1; call rebalance(curr%parent,grow,ierr)
            if(ierr.ne.0) then
             if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): add: rebalancing failed: ",i11)') ierr
             dict_search=GFC_CORRUPTED_CONT; curr=>NULL(); return
            endif
            if(present(value_in)) then
             if(present(store_by)) then
#ifdef NO_GNU
              if(present(copy_ctor_val_f)) then
               call curr%dict_elem_ctor(key,value_in,j,assoc_key=dict%key_storage,assoc_val=store_by,&
                                       &val_copy_ctor_f=copy_ctor_val_f)
              else
#endif
               call curr%dict_elem_ctor(key,value_in,j,assoc_key=dict%key_storage,assoc_val=store_by)
#ifdef NO_GNU
              endif
#endif
             else
#ifdef NO_GNU
              if(present(copy_ctor_val_f)) then
               call curr%dict_elem_ctor(key,value_in,j,assoc_key=dict%key_storage,val_copy_ctor_f=copy_ctor_val_f)
              else
#endif
               call curr%dict_elem_ctor(key,value_in,j,assoc_key=dict%key_storage)
#ifdef NO_GNU
              endif
#endif
             endif
             if(j.ne.GFC_SUCCESS) then
              if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): add: dictionary element construction failed!")')
              dict_search=GFC_MEM_ALLOC_FAILED; curr=>NULL(); return
             endif
            else
             if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): add: absent input value!")')
             dict_search=GFC_INVALID_ARGS; curr=>NULL(); return
            endif
            if(present(value_out)) then; value_out=>curr%get_value(ierr); else; ierr=GFC_SUCCESS; endif
            if(ierr.eq.GFC_SUCCESS) then; call this%jump_(curr); else; dict_search=ierr; endif
           else !empty dictionary
            allocate(curr,STAT=j)
            if(j.ne.0) then
             if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): add: dictionary root allocation failed!")')
             dict_search=GFC_MEM_ALLOC_FAILED; curr=>NULL(); return
            endif
            curr%parent=>NULL(); curr%child_lt=>NULL(); curr%child_gt=>NULL(); curr%balance_fac=0
            if(present(value_in)) then
             if(present(store_by)) then
#ifdef NO_GNU
              if(present(copy_ctor_val_f)) then
               call curr%dict_elem_ctor(key,value_in,j,assoc_key=dict%key_storage,assoc_val=store_by,&
                                       &val_copy_ctor_f=copy_ctor_val_f)
              else
#endif
               call curr%dict_elem_ctor(key,value_in,j,assoc_key=dict%key_storage,assoc_val=store_by)
#ifdef NO_GNU
              endif
#endif
             else
#ifdef NO_GNU
              if(present(copy_ctor_val_f)) then
               call curr%dict_elem_ctor(key,value_in,j,assoc_key=dict%key_storage,val_copy_ctor_f=copy_ctor_val_f)
              else
#endif
               call curr%dict_elem_ctor(key,value_in,j,assoc_key=dict%key_storage)
#ifdef NO_GNU
              endif
#endif
             endif
             if(j.ne.GFC_SUCCESS) then
              if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): add: dictionary root construction failed!")')
              dict_search=GFC_MEM_ALLOC_FAILED; curr=>NULL(); return
             endif
            else
             if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): add: absent input value for root!")')
             dict_search=GFC_INVALID_ARGS; curr=>NULL(); return
            endif
            call dict%reroot_(curr); ierr=this%reset()
            if(ierr.eq.GFC_SUCCESS) then
             if(present(value_out)) then; value_out=>curr%get_value(ierr); if(ierr.ne.GFC_SUCCESS) dict_search=ierr; endif
            else
             dict_search=ierr
            endif
           endif
          elseif(dict_search.eq.GFC_FOUND) then !return the found entry value (just in case)
           if(present(value_out)) then; value_out=>curr%get_value(ierr); else; ierr=GFC_SUCCESS; endif
           if(ierr.eq.GFC_SUCCESS) then; call this%jump_(curr); else; dict_search=ierr; endif
          endif
         case default
          if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search): unknown search action request: ",i11)') action
          dict_search=GFC_UNKNOWN_REQUEST; curr=>NULL(); return
         end select
         curr=>NULL()
         return

         contains

          subroutine rebalance(cr,gr,ier)
           class(dict_elem_t), pointer:: cr   !intermediate element on the way back to root
           integer(INTD), intent(in):: gr     !in(right/left subtree change): {-1;+1}
           integer(INTD), intent(inout):: ier !in(height decrease/increase): {-1;+1}; out: {0;{error_codes}}
           class(dict_elem_t), pointer:: cr_ptr
           integer(INTD):: jb,jc,jg,jd
!          if(DEBUG) write(CONS_OUT,'("Entered rebalance:")') !debug
           jd=ier; jg=gr; cr_ptr=>cr
           do while(jg.ne.0)
!           if(DEBUG) write(CONS_OUT,'(" Rebalance at level ",i3)') lev_p !debug
            cr_ptr%balance_fac=cr_ptr%balance_fac+jg*jd
            if(iabs(cr_ptr%balance_fac).ge.2) then !rotations needed
             if(cr_ptr%balance_fac.eq.-2) then
              jb=cr_ptr%child_gt%balance_fac
              if(jb.gt.0) then !{+1}
               jc=cr_ptr%child_gt%child_lt%balance_fac
               call rotate_double_left(cr_ptr)
               if(dict_search.ne.GFC_FOUND.and.dict_search.ne.GFC_NOT_FOUND) then; cr_ptr=>NULL(); ier=1; return; endif
               cr_ptr%parent%balance_fac=0; if(jd.gt.0) jg=0
               if(jc.gt.0) then
                cr_ptr%balance_fac=0; cr_ptr%parent%child_gt%balance_fac=-1
               elseif(jc.lt.0) then
                cr_ptr%balance_fac=1; cr_ptr%parent%child_gt%balance_fac=0
               else !jc=0
                cr_ptr%balance_fac=0; cr_ptr%parent%child_gt%balance_fac=0
               endif
              else !{-1;0}
               call rotate_simple_left(cr_ptr)
               if(dict_search.ne.GFC_FOUND.and.dict_search.ne.GFC_NOT_FOUND) then; cr_ptr=>NULL(); ier=2; return; endif
               if(jb.eq.0) then
                cr_ptr%balance_fac=-1; cr_ptr%parent%balance_fac=1; if(jd.lt.0) jg=0
               else
                cr_ptr%balance_fac=0; cr_ptr%parent%balance_fac=0; if(jd.gt.0) jg=0
               endif
              endif
              cr_ptr=>cr_ptr%parent
             elseif(cr_ptr%balance_fac.eq.2) then
              jb=cr_ptr%child_lt%balance_fac
              if(jb.lt.0) then !{-1}
               jc=cr_ptr%child_lt%child_gt%balance_fac
               call rotate_double_right(cr_ptr)
               if(dict_search.ne.GFC_FOUND.and.dict_search.ne.GFC_NOT_FOUND) then; cr_ptr=>NULL(); ier=3; return; endif
               cr_ptr%parent%balance_fac=0; if(jd.gt.0) jg=0
               if(jc.lt.0) then
                cr_ptr%balance_fac=0; cr_ptr%parent%child_lt%balance_fac=1
               elseif(jc.gt.0) then
                cr_ptr%balance_fac=-1; cr_ptr%parent%child_lt%balance_fac=0
               else !jc=0
                cr_ptr%balance_fac=0; cr_ptr%parent%child_lt%balance_fac=0
               endif
              else !{0;+1}
               call rotate_simple_right(cr_ptr)
               if(dict_search.ne.GFC_FOUND.and.dict_search.ne.GFC_NOT_FOUND) then; cr_ptr=>NULL(); ier=4; return; endif
               if(jb.eq.0) then
                cr_ptr%balance_fac=+1; cr_ptr%parent%balance_fac=-1; if(jd.lt.0) jg=0
               else
                cr_ptr%balance_fac=0; cr_ptr%parent%balance_fac=0; if(jd.gt.0) jg=0
               endif
              endif
              cr_ptr=>cr_ptr%parent
             else
              if(VERBOSE)&
              &write(CONS_OUT,'("#ERROR(gfc::dictionary::search): rebalance: invalid balance factor: ",i11)') cr_ptr%balance_fac
              cr_ptr=>NULL(); ier=5; return
             endif
            else !node balance factor changed to {-1;0;+1}
             if(jd.gt.0) then
              if(cr_ptr%balance_fac.eq.0) jg=0
             elseif(jd.lt.0) then
              if(cr_ptr%balance_fac.ne.0) jg=0
             endif
            endif
            if(associated(cr_ptr%parent)) then
             jg=iabs(jg); if(associated(cr_ptr%parent%child_gt,cr_ptr)) jg=-jg
             cr_ptr=>cr_ptr%parent; lev_p=lev_p-1
            else
             call dict%reroot_(cr_ptr)
             exit
            endif
           enddo
           ier=0; cr_ptr=>NULL()
!          if(DEBUG) write(CONS_OUT,'("Exited rebalance.")') !debug
           return
          end subroutine rebalance

          subroutine rotate_double_left(cr)
           class(dict_elem_t), pointer:: cr
           call rotate_simple_right(cr%child_gt)
           if(dict_search.ne.GFC_FOUND.and.dict_search.ne.GFC_NOT_FOUND) return
           call rotate_simple_left(cr)
           return
          end subroutine rotate_double_left

          subroutine rotate_double_right(cr)
           class(dict_elem_t), pointer:: cr
           call rotate_simple_left(cr%child_lt)
           if(dict_search.ne.GFC_FOUND.and.dict_search.ne.GFC_NOT_FOUND) return
           call rotate_simple_right(cr)
           return
          end subroutine rotate_double_right

          subroutine rotate_simple_left(cr)
           class(dict_elem_t), pointer:: cr
           class(dict_elem_t), pointer:: jp,jq,js
!          if(DEBUG) write(CONS_OUT,'("  Rotating left")') !debug
           jp=>cr; jq=>cr%child_gt; js=>jq%child_lt !js may be null
           jq%child_lt=>jp
           if(associated(jp%parent)) then
            if(associated(jp%parent%child_lt,jp)) then
             jp%parent%child_lt=>jq
            elseif(associated(jp%parent%child_gt,jp)) then
             jp%parent%child_gt=>jq
            else
             if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search:rotate_simple_left): broken parental link!")')
             dict_search=GFC_CORRUPTED_CONT; return
            endif
            jq%parent=>jp%parent
           else
            jq%parent=>NULL()
           endif
           jp%parent=>jq
           jp%child_gt=>js; if(associated(js)) js%parent=>jp
           return
          end subroutine rotate_simple_left

          subroutine rotate_simple_right(cr)
           class(dict_elem_t), pointer:: cr
           class(dict_elem_t), pointer:: jp,jq,js
!          if(DEBUG) write(CONS_OUT,'("  Rotating right")') !debug
           jp=>cr; jq=>cr%child_lt; js=>jq%child_gt !js may be null
           jq%child_gt=>jp
           if(associated(jp%parent)) then
            if(associated(jp%parent%child_lt,jp)) then
             jp%parent%child_lt=>jq
            elseif(associated(jp%parent%child_gt,jp)) then
             jp%parent%child_gt=>jq
            else
             if(VERBOSE) write(CONS_OUT,'("#ERROR(gfc::dictionary::search:rotate_simple_right): broken parental link!")')
             dict_search=GFC_CORRUPTED_CONT; return
            endif
            jq%parent=>jp%parent
           else
            jq%parent=>NULL()
           endif
           jp%parent=>jq
           jp%child_lt=>js; if(associated(js)) js%parent=>jp
           return
          end subroutine rotate_simple_right

        end function DictionaryIterSearch
!---------------------------------------------------------------------
        subroutine DictionaryIterSortToList(this,list_refs,ierr,order)
!Returns a bidirectional linked list of references to dictionary elements in a sorted order.
!The input dictionary iterator is allowed to be a subdictionary iterator.
!The current dictionary iterator position is kept intact, unless the
!iterator in the GFC_IT_DONE state, in which case it will be reset.
!The bidirectional linked list must be empty on entrance.
         implicit none
         class(dictionary_iter_t), intent(inout):: this !in: dictionary (or subdictionary) iterator
         class(list_bi_t), intent(inout):: list_refs    !out: bidirectional list of references to dictionary elements (in:empty,out:filled)
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD), intent(in), optional:: order    !in: sorting order: {GFC_ASCEND_ORDER,GFC_DESCEND_ORDER}, defaults to GFC_ASCEND_ORDER
         integer(INTD):: errc,ier,ord
         logical:: dir
         type(list_iter_t):: list_it
         class(dict_elem_t), pointer:: orig

         errc=list_it%init(list_refs)
         if(errc.eq.GFC_SUCCESS) then
          if(list_it%get_status().eq.GFC_IT_EMPTY) then
           errc=this%get_status()
           if(errc.eq.GFC_IT_DONE) then; errc=this%reset(); errc=this%get_status(); endif
           if(errc.eq.GFC_IT_ACTIVE) then
            if(present(order)) then; ord=order; else; ord=GFC_ASCEND_ORDER; endif
            orig=>this%current !save the original iterator position
            if(ord.eq.GFC_ASCEND_ORDER) then
             errc=this%move_to_min(); dir=GFC_DICT_SUCCESSOR
            else
             errc=this%move_to_max(); dir=GFC_DICT_PREDECESSOR
            endif
            if(errc.eq.GFC_SUCCESS) then
             sloop: do while(errc.eq.GFC_SUCCESS)
              errc=list_it%append(this%current,assoc_only=.TRUE.)
              if(errc.ne.GFC_SUCCESS) exit sloop
              errc=this%move_in_order(dir)
             enddo sloop
             if(errc.eq.GFC_IT_DONE) errc=this%reset()
            else
             errc=GFC_CORRUPTED_CONT
            endif
            call this%jump_(orig) !return to the original iterator position
           else
            if(errc.ne.GFC_IT_EMPTY) errc=GFC_NULL_CONT
           endif
          else
           errc=GFC_INVALID_ARGS
          endif
          ier=list_it%release()
          if(errc.eq.GFC_SUCCESS.and.ier.ne.GFC_SUCCESS) errc=ier
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DictionaryIterSortToList
!----------------------------------------------------------------------
        subroutine DictionaryIterSortToVector(this,vec_refs,ierr,order)
!Returns a vector of references to dictionary elements in a sorted order.
!The input dictionary iterator is allowed to be a subdictionary iterator.
!The current dictionary iterator position is kept intact, unless the
!iterator in the GFC_IT_DONE state, in which case it will be reset.
!The vector must be empty on entrance.
         implicit none
         class(dictionary_iter_t), intent(inout):: this !in: dictionary (or subdictionary) iterator
         class(vector_t), intent(inout):: vec_refs      !out: vector of references to dictionary elements (in:empty,out:filled)
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD), intent(in), optional:: order    !in: sorting order: {GFC_ASCEND_ORDER,GFC_DESCEND_ORDER}, defaults to GFC_ASCEND_ORDER
         integer(INTD):: errc,ier,ord
         logical:: dir
         type(vector_iter_t):: vec_it
         class(dict_elem_t), pointer:: orig

         errc=vec_it%init(vec_refs)
         if(errc.eq.GFC_SUCCESS) then
          if(vec_it%get_status().eq.GFC_IT_EMPTY) then
           errc=this%get_status()
           if(errc.eq.GFC_IT_DONE) then; errc=this%reset(); errc=this%get_status(); endif
           if(errc.eq.GFC_IT_ACTIVE) then
            if(present(order)) then; ord=order; else; ord=GFC_ASCEND_ORDER; endif
            orig=>this%current !save the original iterator position
            if(ord.eq.GFC_ASCEND_ORDER) then
             errc=this%move_to_min(); dir=GFC_DICT_SUCCESSOR
            else
             errc=this%move_to_max(); dir=GFC_DICT_PREDECESSOR
            endif
            if(errc.eq.GFC_SUCCESS) then
             sloop: do while(errc.eq.GFC_SUCCESS)
              errc=vec_it%append(this%current,assoc_only=.TRUE.)
              if(errc.ne.GFC_SUCCESS) exit sloop
              errc=this%move_in_order(dir)
             enddo sloop
             if(errc.eq.GFC_IT_DONE) errc=this%reset()
            else
             errc=GFC_CORRUPTED_CONT
            endif
            call this%jump_(orig) !return to the original iterator position
           else
            if(errc.ne.GFC_IT_EMPTY) errc=GFC_NULL_CONT
           endif
          else
           errc=GFC_INVALID_ARGS
          endif
          ier=vec_it%release()
          if(errc.eq.GFC_SUCCESS.and.ier.ne.GFC_SUCCESS) errc=ier
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine DictionaryIterSortToVector

       end module gfc_dictionary
!===============================
!TESTING:
!--------------------------------
       module gfc_dictionary_test
        use gfc_base
        use gfc_list
        use gfc_dictionary
        use timers, only: thread_wtime
        implicit none
        private
!PARAMETERS:
        integer(INTD), parameter, private:: KEY_LEN=6
        integer(INTD), parameter, private:: MAX_IND_VAL=7
!TYPES:
 !Key:
        type, private:: key_t
         integer(INTD):: rank=KEY_LEN
         integer(INTD):: dims(1:KEY_LEN)
        end type key_t
 !Value:
        type, private:: val_t
         real(8):: my_array(1:KEY_LEN)
         type(key_t):: key_stored
        end type val_t
!VISIBILITY:
        private cmp_key_test
        public test_gfc_dictionary

       contains
!-------------------------------------------------
        function cmp_key_test(up1,up2) result(cmp)
         implicit none
         integer(INTD):: cmp
         class(*), intent(in), target:: up1,up2
         integer(INTD):: i

         cmp=GFC_CMP_ERR
         select type(up1)
         class is(key_t)
          select type(up2)
          class is(key_t)
           if(up1%rank.lt.up2%rank) then
            cmp=GFC_CMP_LT
           elseif(up1%rank.gt.up2%rank) then
            cmp=GFC_CMP_GT
           else
            cmp=GFC_CMP_EQ
            do i=1,up1%rank
             if(up1%dims(i).lt.up2%dims(i)) then
              cmp=GFC_CMP_LT; exit
             elseif(up1%dims(i).gt.up2%dims(i)) then
              cmp=GFC_CMP_GT; exit
             endif
            enddo
           endif
          end select
         end select
         return
        end function cmp_key_test
!----------------------------------------------------
        function print_key(obj,dev_id) result(ierr)
         implicit none
         integer(INTD):: ierr
         class(*), intent(in):: obj
         integer(INTD), intent(in), optional:: dev_id
         integer(INTD):: dev

         ierr=GFC_SUCCESS
         if(present(dev_id)) then; dev=dev_id; else; dev=6; endif
         select type(obj)
         class is(key_t)
          write(dev,'(32(1x,i9))') obj%dims(1:obj%rank)
         class default
          ierr=GFC_ACTION_FAILED
         end select
         return
        end function print_key
!----------------------------------------------------
        function print_value(obj,dev_id) result(ierr)
         implicit none
         integer(INTD):: ierr
         class(*), intent(in):: obj
         integer(INTD), intent(in), optional:: dev_id
         integer(INTD):: dev

         ierr=GFC_SUCCESS
         if(present(dev_id)) then; dev=dev_id; else; dev=6; endif
         select type(obj)
         class is(val_t)
          write(dev,'(32(1x,D15.7))') obj%my_array(1:obj%key_stored%rank)
         class default
          ierr=GFC_ACTION_FAILED
         end select
         return
        end function print_value
!--------------------------------------------------------------
        function test_gfc_dictionary(perf,dev_out) result(ierr)
         implicit none
         integer(INTD):: ierr                          !out: error code (0:success)
         real(8), intent(out):: perf                   !out: performance index
         integer(INTD), intent(in), optional:: dev_out !in: default output device
!------------------------------------------------------
         integer(INTD), parameter:: MAX_ACTIONS=1000000
!------------------------------------------------------
         integer(INTD):: jo,i,j,nadd,ndel,nidl
         type(key_t):: key
         type(val_t):: val
         type(dictionary_t), target:: some_dict
         type(dictionary_iter_t):: dict_it
         type(list_bi_t):: list_refs
         type(list_iter_t):: list_it
         class(gfc_cont_elem_t), pointer:: pntee
         class(*), pointer:: uptr
         real(8):: tms,tm
         logical:: del

         ierr=GFC_SUCCESS; perf=0d0; uptr=>NULL()
         if(present(dev_out)) then; jo=dev_out; else; jo=6; endif
!Lookups/insertions:
         tms=thread_wtime()
         nadd=0; ndel=0; nidl=0
         j=dict_it%init(some_dict); if(j.ne.GFC_SUCCESS) then; call test_quit(1); return; endif
         do i=1,MAX_ACTIONS
          del=(mod(i,4).eq.3) !every fourth action will be a deletion
          key%rank=KEY_LEN; call get_rnd_key(key) !random key
          if(del) then !deletion (if found)
           j=dict_it%search(GFC_DICT_DELETE_IF_FOUND,cmp_key_test,key)
          else !insertion (if not found)
           val%my_array(1:KEY_LEN)=13.777d0; val%key_stored=key !value
           j=dict_it%search(GFC_DICT_ADD_IF_NOT_FOUND,cmp_key_test,key,val,GFC_BY_VAL,uptr)
          endif
          if(j.eq.GFC_FOUND) then
           if(del) then
            ndel=ndel+1
           else
            nidl=nidl+1
            if(.not.(associated(uptr))) then; call test_quit(2); return; endif
           endif
          elseif(j.eq.GFC_NOT_FOUND) then
           if(del) then
            nidl=nidl+1
           else
            nadd=nadd+1
            if(.not.(associated(uptr))) then; call test_quit(3); return; endif
           endif
          else
           write(jo,'("#DEBUG(gfc::dictionary:test): Dictionary search failed with error ",i11)') j !debug
           call test_quit(4); return
          endif
         enddo
         tm=thread_wtime(tms)
         perf=dble(MAX_ACTIONS)/tm
         !write(jo,'("Added ",i11,"; Deleted ",i11,"; Idle ",i11)') nadd,ndel,nidl !debug
!Traversing the final dictionary as a sorted list:
         tms=thread_wtime()
         call dict_it%sort_to_list(list_refs,j); if(j.ne.GFC_SUCCESS) then; call test_quit(5); return; endif
         j=list_it%init(list_refs); pntee=>NULL(); uptr=>NULL(); i=0
         do while(j.eq.GFC_SUCCESS)
          pntee=>list_it%pointee(); uptr=>pntee%get_value()
          select type(uptr)
          class is(dict_elem_t)
           i=i+1
           !call uptr%print_key(print_key,j,jo); if(j.ne.0) then; call test_quit(6); return; endif
          class default
           call test_quit(7); return
          end select
          j=list_it%next()
         enddo
         if(j.ne.GFC_IT_DONE) then; call test_quit(8); return; endif
         j=list_it%release(); if(j.ne.GFC_SUCCESS) then; call test_quit(9); return; endif
         tm=thread_wtime(tms)
         !write(jo,'("Traversing time for ",i11," elements is ",F10.4," sec")') i,tm !debug
!Success:
         call test_quit(0)
         return

         contains

          subroutine get_rnd_key(akey)
           type(key_t), intent(inout):: akey
           integer(INTD):: jj
           real(8):: jv(1:KEY_LEN)
           call random_number(jv(1:KEY_LEN))
           do jj=1,akey%rank
            akey%dims(jj)=nint(jv(jj)*dble(MAX_IND_VAL))
           enddo
           return
          end subroutine get_rnd_key

          subroutine test_quit(jerr)
           integer(INTD), intent(in):: jerr
           integer(INTD):: jj
           ierr=jerr
           if(ierr.ne.GFC_SUCCESS) then
            write(jo,'("#ERROR(gfc::dictionary::test): Test failed: Error code ",i13)') ierr
            write(jo,'("Please contact the developer at QUANT4ME@GMAIL.COM")')
           endif
           call dict_it%delete_all(jj); if(ierr.eq.GFC_SUCCESS.and.jj.ne.GFC_SUCCESS) ierr=10
           if(jj.ne.GFC_SUCCESS) write(jo,'("#ERROR(gfc::dictionary::test): Dictionary destruction failed: Error code ",i13)') jj
           jj=dict_it%release(); if(ierr.eq.GFC_SUCCESS.and.jj.ne.GFC_SUCCESS) ierr=11
           if(jj.ne.GFC_SUCCESS) write(jo,'("#ERROR(gfc::dictionary::test): Dictionary iterator release failed: Error code ",i13)')&
                                      &jj
           return
          end subroutine test_quit

        end function test_gfc_dictionary

       end module gfc_dictionary_test
