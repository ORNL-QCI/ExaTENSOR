!Generic Fortran Containers (GFC): Bi-directional linked list
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2017-11-09 (started 2016-02-28)

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

!NOTES:
! # A list is a linked derivative of an abstract (unlinked) GFC container.
!   A sublist is a list incorporated into another (composite) list.
! # All accesses, updates, scans, and actions are performed on a list
!   via the list iterator associated with the list. Elements of a sublist
!   can be accessed by both the original and the combined list iterators.
!   Multiple iterators can be associated with a list at the same time.

!FOR DEVELOPERS ONLY:
! # Boundary/current element reference counting needs to be implemented.

       module gfc_list
        use gfc_base
        use gfc_vector
        use timers
        implicit none
        private
!PARAMETERS:
!Basic:
        integer(INTD), private:: CONS_OUT=6 !output device
        logical, private:: VERBOSE=.TRUE.   !verbositiy for errors
        integer(INTD), private:: DEBUG=0    !debugging level (0:none)
!TYPES:
 !Linked list element:
        type, extends(gfc_cont_elem_t), public:: list_elem_t
         class(list_elem_t), pointer, private:: next_elem=>NULL()
         class(list_elem_t), pointer, private:: prev_elem=>NULL()
         contains
          procedure, private:: ListElemConstruct
          generic, public:: list_elem_ctor=>ListElemConstruct !constructs the content of the list element
          procedure, public:: is_first=>ListElemIsFirst       !returns GFC_TRUE if the element is the first in the list
          procedure, public:: is_last=>ListElemIsLast         !returns GFC_TRUE if the element is the last in the list
        end type list_elem_t
 !Iterator position in the linked list:
        type, public:: list_pos_t
         class(list_elem_t), pointer, private:: elem_p=>NULL() !pointer to a list element
         contains
          procedure, public:: is_set=>ListPosIsSet !returns TRUE if the list position is set
          procedure, public:: clean=>ListPosClean  !resets the list position to NULL
        end type list_pos_t
 !Linked list:
        type, extends(gfc_container_t), public:: list_bi_t
         class(list_elem_t), pointer, private:: first_elem=>NULL() !first element of the linked list
         class(list_elem_t), pointer, private:: last_elem=>NULL()  !last element of the linked list
         contains
          procedure, public:: is_empty=>ListIsEmpty                !returns GFC_TRUE if the list is empty, GFC_FALSE otherwise (or error code)
          procedure, public:: is_sublist=>ListIsSublist            !returns TRUE if the list is a sublist of a larger list, FALSE otherwise
          procedure, private:: ListDuplicateToList                 !duplicates a list into another list either by value or by reference
          procedure, private:: ListDuplicateToVector               !duplicates a list into a vector either by value or by reference
          generic, public:: duplicate=>ListDuplicateToList,ListDuplicateToVector
          procedure, private:: ListBiAssign                        !copy assignment
          generic, public:: assignment(=)=>ListBiAssign
          final:: list_bi_dtor                                     !dtor
        end type list_bi_t
 !List iterator:
        type, extends(gfc_iter_t), public:: list_iter_t
         class(list_elem_t), pointer, private:: current=>NULL() !current element of the linked list
         class(list_bi_t), pointer, private:: container=>NULL() !linked list associated with the iterator
         contains
          procedure, public:: init=>ListIterInit                   !initializes the iterator by associating it with a list and resetting to the beginning
          procedure, public:: reset=>ListIterReset                 !resets the iterator to the beginning of the list
          procedure, public:: reset_back=>ListIterResetBack        !resets the iterator to the end of the list
          procedure, public:: release=>ListIterRelease             !releases the iterator (dissociates it from the container)
          procedure, public:: pointee=>ListIterPointee             !returns the container element currently pointed to by the iterator
          procedure, public:: next=>ListIterNext                   !moves the iterator to the next list element
          procedure, public:: previous=>ListIterPrevious           !moves the iterator to the previous list element
          procedure, public:: on_first=>ListIterOnFirst            !returns TRUE if the iterator is positioned on the first element of the list
          procedure, public:: on_last=>ListIterOnLast              !returns TRUE if the iterator is positioned on the last element of the list
          procedure, public:: append=>ListIterAppend               !inserts a new element either at the beginning or at the end of the container
          procedure, public:: insert_elem=>ListIterInsertElem      !inserts a new element at the current position of the container
          procedure, public:: insert_list=>ListIterInsertList      !inserts another linked list at the current position of the container
          procedure, public:: move_elem=>ListIterMoveElem          !moves a list element from one iterator to another
          procedure, public:: move_list=>ListIterMoveList          !moves an entire list from one iterator to another
          procedure, public:: split=>ListIterSplit                 !splits the list into two parts at the current position
          procedure, public:: delete=>ListIterDelete               !deletes the list element in the current position
          procedure, public:: delete_sublist=>ListIterDeleteSublist!deletes all elements either prior or after the current iterator position
          procedure, public:: delete_all=>ListIterDeleteAll        !deletes all list elements
          procedure, public:: bookmark=>ListIterBookmark           !bookmarks the current iterator position
          procedure, public:: jump=>ListIterJump                   !jumps to a previously bookmarked iterator position
          procedure, public:: jump_=>ListIterJump_                 !PRIVATE: moves the iterator to a specified list element
          procedure, private:: pop_=>ListIterPop_                  !PRIVATE: Pops up the element at the current iterator position without destroying it
        end type list_iter_t
!INTERFACES:
!VISIBILITY:
 !list_elem_t:
        private ListElemConstruct
        private ListElemIsFirst
        private ListElemIsLast
 !list_pos_t:
        private ListPosIsSet
        private ListPosClean
 !list_bi_t:
        private ListIsEmpty
        private ListIsSublist
        private ListDuplicateToList
        private ListDuplicateToVector
        private ListBiAssign
        public list_bi_dtor
 !list_iter_t:
        private ListIterInit
        private ListIterReset
        private ListIterResetBack
        private ListIterRelease
        private ListIterPointee
        private ListIterNext
        private ListIterPrevious
        private ListIterOnFirst
        private ListIterOnLast
        private ListIterAppend
        private ListIterInsertElem
        private ListIterInsertList
        private ListIterMoveElem
        private ListIterMoveList
        private ListIterSplit
        private ListIterDelete
        private ListIterDeleteSublist
        private ListIterDeleteAll
        private ListIterBookmark
        private ListIterJump
        private ListIterJump_
        private ListIterPop_

       contains
!IMPLEMENTATION:
!-------------------------------------------------------------------------
#ifdef NO_GNU
        subroutine ListElemConstruct(this,obj,ierr,assoc_only,copy_ctor_f) !`GCC has a bug with this line
#else
        subroutine ListElemConstruct(this,obj,ierr,assoc_only)
#endif
!Constructs the content of the list element.
         implicit none
         class(list_elem_t), intent(inout):: this      !inout: element of a list
         class(*), target, intent(in):: obj            !in: value to be stored
         integer(INTD), intent(out), optional:: ierr   !out: error code
         logical, intent(in), optional:: assoc_only    !in: if TRUE, the value will be assigned by reference, otherwise by value (allocated): Defaults to FALSE
#ifdef NO_GNU
         procedure(gfc_copy_i), optional:: copy_ctor_f !in: generic copy constructor
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
        end subroutine ListElemConstruct
!------------------------------------------------------
        function ListElemIsFirst(this,ierr) result(res)
!Returns GFC_TRUE if this list element is the first in the list.
         implicit none
         integer(INTD):: res                         !out: result
         class(list_elem_t), intent(in):: this       !in: list element
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc

         errc=GFC_SUCCESS; res=GFC_TRUE
         if(associated(this%prev_elem)) res=GFC_FALSE
         if(present(ierr)) ierr=errc
         return
        end function ListElemIsFirst
!-----------------------------------------------------
        function ListElemIsLast(this,ierr) result(res)
!Returns GFC_TRUE if this list element is the last in the list.
         implicit none
         integer(INTD):: res                         !out: result
         class(list_elem_t), intent(in):: this       !in: list element
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc

         errc=GFC_SUCCESS; res=GFC_TRUE
         if(associated(this%next_elem)) res=GFC_FALSE
         if(present(ierr)) ierr=errc
         return
        end function ListElemIsLast
![list_pos_t]==================================
        function ListPosIsSet(this) result(ans)
!Returns TRUE if the list position is set, FALSE otherwise
         implicit none
         logical:: ans                        !out: answer
         class(list_pos_t), intent(in):: this !in: list position

         ans=associated(this%elem_p)
         return
        end function ListPosIsSet
!------------------------------------
        subroutine ListPosClean(this)
!Resets the list position to an empty state.
         implicit none
         class(list_pos_t), intent(inout):: this !inout: list position

         this%elem_p=>NULL()
         return
        end subroutine ListPosClean
!---------------------------------------------
        function ListIsEmpty(this) result(res)
!Returns GFC_TRUE if the list is empty, GFC_FALSE otherwise (or error code).
         implicit none
         integer(INTD):: res                 !out: result of query (or error code)
         class(list_bi_t), intent(in):: this !in: list

         if(associated(this%first_elem)) then
          if(associated(this%last_elem)) then
           res=GFC_FALSE
          else
           res=GFC_CORRUPTED_CONT
          endif
         else
          if(associated(this%last_elem)) then
           res=GFC_CORRUPTED_CONT
          else
           res=GFC_TRUE
          endif
         endif
         return
        end function ListIsEmpty
!----------------------------------------------------
        function ListIsSublist(this,ierr) result(res)
!Returns TRUE if the list is a sublist of another list, FALSE otherwise.
         implicit none
         logical:: res                               !out: result
         class(list_bi_t), intent(in):: this         !in: list
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc

         errc=GFC_SUCCESS; res=.FALSE.
         if(associated(this%first_elem).and.associated(this%last_elem)) then
          res=(associated(this%first_elem%prev_elem).or.associated(this%last_elem%next_elem))
         else
          if(associated(this%first_elem).or.associated(this%last_elem)) then
           errc=GFC_CORRUPTED_CONT
          else
           errc=GFC_EMPTY_CONT
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end function ListIsSublist
!----------------------------------------------------------------------------------
#ifdef NO_GNU
        function ListDuplicateToList(this,list,assoc_only,copy_ctor_f) result(ierr) !`GCC bug
#else
        function ListDuplicateToList(this,list,assoc_only) result(ierr)
#endif
!Duplicates a list into another list either by value or by reference.
         implicit none
         integer(INTD):: ierr                          !out: error code
         class(list_bi_t), intent(in):: this           !in: input list (must be non-empty on input)
         class(list_bi_t), intent(inout):: list        !out: output duplicate list (must be empty on input)
         logical, intent(in), optional:: assoc_only    !in: if TRUE, the list will be duplicated by reference, otherwise by value (default)
#ifdef NO_GNU
         procedure(gfc_copy_i), optional:: copy_ctor_f !in: optional copy constructor (may be needed when duplicating by value)
#endif
         integer(INTD):: errc
         type(list_iter_t):: ilit,olit
         class(*), pointer:: up
         logical:: assoc

         ierr=GFC_SUCCESS
         if(this%is_empty().eq.GFC_FALSE) then
          if(list%is_empty().eq.GFC_TRUE) then
           ierr=ilit%init(this)
           if(ierr.eq.GFC_SUCCESS) then
            ierr=olit%init(list)
            if(ierr.eq.GFC_SUCCESS) then
             assoc=.FALSE.; if(present(assoc_only)) assoc=assoc_only
             errc=GFC_SUCCESS
#ifdef NO_GNU
             if(present(copy_ctor_f)) then
              do while(errc.eq.GFC_SUCCESS)
               up=>ilit%get_value(ierr); if(ierr.ne.GFC_SUCCESS) exit
               ierr=olit%append(up,.FALSE.,assoc,copy_ctor_f)
               errc=ilit%next()
              enddo
             else
#endif
              do while(errc.eq.GFC_SUCCESS)
               up=>ilit%get_value(ierr); if(ierr.ne.GFC_SUCCESS) exit
               ierr=olit%append(up,.FALSE.,assoc)
               errc=ilit%next()
              enddo
#ifdef NO_GNU
             endif
#endif
             if(errc.eq.GFC_NO_MOVE) errc=GFC_SUCCESS
             if(errc.ne.GFC_SUCCESS) ierr=errc
             errc=olit%release(); if(errc.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=errc
            endif
            errc=ilit%release(); if(errc.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=errc
           endif
          else
           ierr=GFC_INVALID_ARGS
          endif
         else
          ierr=GFC_EMPTY_CONT
         endif
         return
        end function ListDuplicateToList
!--------------------------------------------------------------------------------------
#ifdef NO_GNU
        function ListDuplicateToVector(this,vector,assoc_only,copy_ctor_f) result(ierr) !`GCC bug
#else
        function ListDuplicateToVector(this,vector,assoc_only) result(ierr)
#endif
!Duplicates a list into a vector either by value or by reference.
         implicit none
         integer(INTD):: ierr                          !out: error code
         class(list_bi_t), intent(in):: this           !in: input list (must be non-empty on input)
         class(vector_t), intent(inout):: vector       !out: output vector (must be empty on input)
         logical, intent(in), optional:: assoc_only    !in: if TRUE, the list will be duplicated by reference, otherwise by value (default)
#ifdef NO_GNU
         procedure(gfc_copy_i), optional:: copy_ctor_f !in: optional copy constructor (may be needed when duplicating by value)
#endif
         integer(INTD):: errc
         type(list_iter_t):: ilit
         type(vector_iter_t):: ovit
         class(*), pointer:: up
         logical:: assoc

         ierr=GFC_SUCCESS
         if(this%is_empty().eq.GFC_FALSE) then
          if(vector%is_empty().eq.GFC_TRUE) then
           ierr=ilit%init(this)
           if(ierr.eq.GFC_SUCCESS) then
            ierr=ovit%init(vector)
            if(ierr.eq.GFC_SUCCESS) then
             assoc=.FALSE.; if(present(assoc_only)) assoc=assoc_only
             errc=GFC_SUCCESS
#ifdef NO_GNU
             if(present(copy_ctor_f)) then
              do while(errc.eq.GFC_SUCCESS)
               up=>ilit%get_value(ierr); if(ierr.ne.GFC_SUCCESS) exit
               ierr=ovit%append(up,assoc,copy_ctor_f)
               errc=ilit%next()
              enddo
             else
#endif
              do while(errc.eq.GFC_SUCCESS)
               up=>ilit%get_value(ierr); if(ierr.ne.GFC_SUCCESS) exit
               ierr=ovit%append(up,assoc)
               errc=ilit%next()
              enddo
#ifdef NO_GNU
             endif
#endif
             if(errc.eq.GFC_NO_MOVE) errc=GFC_SUCCESS
             if(errc.ne.GFC_SUCCESS) ierr=errc
             errc=ovit%release(); if(errc.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=errc
            endif
            errc=ilit%release(); if(errc.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=errc
           endif
          else
           ierr=GFC_INVALID_ARGS
          endif
         else
          ierr=GFC_EMPTY_CONT
         endif
         return
        end function ListDuplicateToVector
!----------------------------------------
        subroutine ListBiAssign(this,src)
!Copy assignment.
         implicit none
         class(list_bi_t), intent(out):: this !out: list
         class(list_bi_t), intent(in):: src   !in: source list
         integer(INTD):: errc

         errc=src%duplicate(this)
         return
        end subroutine ListBiAssign
!------------------------------------
        subroutine list_bi_dtor(this)
         implicit none
         type(list_bi_t):: this
         type(list_iter_t):: lit
         integer(INTD):: errc

         errc=lit%init(this)
         errc=lit%delete_all()
         errc=lit%release()
         return
        end subroutine list_bi_dtor
!----------------------------------------------------
        function ListIterInit(this,cont) result(ierr)
!Initializes the iterator and resets it to the beginning of the container.
         implicit none
         integer(INTD):: ierr                              !out: error code (0:success)
         class(list_iter_t), intent(inout):: this          !inout: iterator
         class(gfc_container_t), target, intent(in):: cont !in: container

         ierr=GFC_SUCCESS
         select type(cont)
         class is(list_bi_t)
          this%container=>cont
          ierr=this%reset()
         class default
          ierr=GFC_INVALID_ARGS
         end select
         return
        end function ListIterInit
!------------------------------------------------
        function ListIterReset(this) result(ierr)
!Resets the iterator to the beginning (first element).
         implicit none
         integer(INTD):: ierr                     !out: error code (0:success)
         class(list_iter_t), intent(inout):: this !inout: iterator

         ierr=GFC_SUCCESS
         if(associated(this%container)) then
          this%current=>this%container%first_elem
          if(associated(this%current)) then
           ierr=this%set_status_(GFC_IT_ACTIVE)
          else
           ierr=this%set_status_(GFC_IT_EMPTY)
          endif
          call this%reset_count() !reset all iteration counters
         else
          this%current=>NULL()
          ierr=this%set_status_(GFC_IT_NULL)
          ierr=GFC_IT_NULL
         endif
         return
        end function ListIterReset
!----------------------------------------------------
        function ListIterResetBack(this) result(ierr)
!Resets the iterator to the end (last element).
         implicit none
         integer(INTD):: ierr                     !out: error code (0:success)
         class(list_iter_t), intent(inout):: this !inout: iterator

         ierr=GFC_SUCCESS
         if(associated(this%container)) then
          this%current=>this%container%last_elem
          if(associated(this%current)) then
           ierr=this%set_status_(GFC_IT_ACTIVE)
          else
           ierr=this%set_status_(GFC_IT_EMPTY)
          endif
          call this%reset_count() !reset all iteration counters
         else
          this%current=>NULL()
          ierr=this%set_status_(GFC_IT_NULL)
          ierr=GFC_IT_NULL
         endif
         return
        end function ListIterResetBack
!--------------------------------------------------
        function ListIterRelease(this) result(ierr)
!Dissociates the iterator from its container.
         implicit none
         integer(INTD):: ierr                     !out: error code (0:success)
         class(list_iter_t), intent(inout):: this !inout: iterator

         this%current=>NULL(); this%container=>NULL()
         call this%reset_count(); ierr=this%set_status_(GFC_IT_NULL)
         return
        end function ListIterRelease
!--------------------------------------------------------
        function ListIterPointee(this,ierr) result(pntee)
!Returns the container element the iterator is currently pointing to.
         implicit none
         class(gfc_cont_elem_t), pointer:: pntee     !out: container element currently pointed to by the iterator
         class(list_iter_t), intent(in):: this       !in: iterator
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc

         errc=this%get_status()
         if(errc.eq.GFC_IT_ACTIVE) then
          pntee=>this%current; errc=GFC_SUCCESS
         else
          pntee=>NULL()
         endif
         if(present(ierr)) ierr=errc
         return
        end function ListIterPointee
!------------------------------------------------------
        function ListIterNext(this,elem_p) result(ierr)
!If <elem_p> is absent, the iterator moves to the next element, if any.
!If <elem_p> is present, the iterator simply returns the next element in <elem_p> without moving.
         implicit none
         integer(INTD):: ierr                                            !out: error code (0:success)
         class(list_iter_t), intent(inout):: this                        !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element
         class(list_elem_t), pointer:: lep

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           ierr=GFC_SUCCESS
           if(.not.associated(this%current,this%container%last_elem)) then
            lep=>this%current%next_elem
            if(present(elem_p)) then
             elem_p=>lep
            else
             this%current=>lep
            endif
           else
            this%current=>NULL()
            ierr=this%set_status_(GFC_IT_DONE)
            lep=>NULL()
           endif
           if(.not.associated(lep)) then; ierr=GFC_NO_MOVE; else; lep=>NULL(); endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function ListIterNext
!----------------------------------------------------------
        function ListIterPrevious(this,elem_p) result(ierr)
!If <elem_p> is absent, the iterator moves to the previous element, if any.
!If <elem_p> is present, the iterator simply returns the previous element in <elem_p> without moving.
         implicit none
         integer(INTD):: ierr                                            !out: error code (0:success)
         class(list_iter_t), intent(inout):: this                        !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element
         class(list_elem_t), pointer:: lep

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           ierr=GFC_SUCCESS
           if(.not.associated(this%current,this%container%first_elem)) then
            lep=>this%current%prev_elem
            if(present(elem_p)) then
             elem_p=>lep
            else
             this%current=>lep
            endif
           else
            this%current=>NULL()
            ierr=this%set_status_(GFC_IT_DONE)
            lep=>NULL()
           endif
           if(.not.associated(lep)) then; ierr=GFC_NO_MOVE; else; lep=>NULL(); endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function ListIterPrevious
!------------------------------------------------------
        function ListIterOnFirst(this,ierr) result(res)
!Returns TRUE if the iterator is positioned on the first element of the list
         logical:: res                               !out: result
         class(list_iter_t), intent(in):: this       !in: list iterator
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         res=.FALSE.; errc=this%get_status()
         if(errc.eq.GFC_IT_ACTIVE) then
          errc=GFC_SUCCESS
          res=associated(this%current,this%container%first_elem)
         endif
         if(present(ierr)) ierr=errc
         return
        end function ListIterOnFirst
!-----------------------------------------------------
        function ListIterOnLast(this,ierr) result(res)
!returns TRUE if the iterator is positioned on the last element of the list
         logical:: res                               !out: result
         class(list_iter_t), intent(in):: this       !in: list iterator
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         res=.FALSE.; errc=this%get_status()
         if(errc.eq.GFC_IT_ACTIVE) then
          errc=GFC_SUCCESS
          res=associated(this%current,this%container%last_elem)
         endif
         if(present(ierr)) ierr=errc
         return
        end function ListIterOnLast
!----------------------------------------------------------------------------------------
#ifdef NO_GNU
        function ListIterAppend(this,elem_val,at_top,assoc_only,copy_ctor_f) result(ierr) !`GCC/5.3.0 has a bug with this
#else
        function ListIterAppend(this,elem_val,at_top,assoc_only) result(ierr)
#endif
!Appends an element at the beginning or at the end of the list/sublist, either by value or by reference.
!For non-empty iterators, the iterator position is kept unchanged (even if it is GFC_IT_DONE).
!If the iterator is empty on input, its status will be changed to GFC_IT_ACTIVE via reset.
         implicit none
         integer(INTD):: ierr                       !out: error code (0:success)
         class(list_iter_t), intent(inout):: this   !inout: iterator
         class(*), target, intent(in):: elem_val    !in: value to be stored
         logical, intent(in), optional:: at_top     !in: TRUE:append at the top, FALSE:append at the end (default)
         logical, intent(in), optional:: assoc_only !in: storage type: TRUE:by reference, FALSE:by value (default)
#ifdef NO_GNU
         procedure(gfc_copy_i), optional:: copy_ctor_f !user-defined generic copy constructor (when storing by value only)
#endif
         logical:: assoc,top
         integer:: errc
         integer(INTL):: nelems
         class(list_elem_t), pointer:: lep,oep

         if(present(at_top)) then; top=at_top; else; top=.FALSE.; endif
         if(present(assoc_only)) then; assoc=assoc_only; else; assoc=.FALSE.; endif
         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE.or.ierr.eq.GFC_IT_DONE) then !non-empty container
          if(associated(this%container)) then
           if(associated(this%container%first_elem).and.associated(this%container%last_elem)) then
            ierr=GFC_SUCCESS
            if(top) then
             oep=>this%container%first_elem%prev_elem !not necessarily NULL for sublists
             allocate(this%container%first_elem%prev_elem,STAT=errc)
             if(errc.eq.0) lep=>this%container%first_elem%prev_elem
            else
             oep=>this%container%last_elem%next_elem !not necessarily NULL for sublists
             allocate(this%container%last_elem%next_elem,STAT=errc)
             if(errc.eq.0) lep=>this%container%last_elem%next_elem
            endif
            if(errc.eq.0) then
#ifdef NO_GNU
             if(present(copy_ctor_f)) then
              call lep%list_elem_ctor(elem_val,ierr,assoc_only=assoc,copy_ctor_f=copy_ctor_f)
             else
#endif
              call lep%list_elem_ctor(elem_val,ierr,assoc_only=assoc)
#ifdef NO_GNU
             endif
#endif
             if(ierr.eq.GFC_SUCCESS) then
              if(top) then
               lep%next_elem=>this%container%first_elem
               lep%prev_elem=>oep
               if(associated(oep)) oep%next_elem=>lep
               this%container%first_elem=>lep
              else
               lep%prev_elem=>this%container%last_elem
               lep%next_elem=>oep
               if(associated(oep)) oep%prev_elem=>lep
               this%container%last_elem=>lep
              endif
              if(this%container%num_elems_().ge.0) then !quick counting is on
               nelems=this%container%update_num_elems_(1_INTL,ierr); if(ierr.ne.GFC_SUCCESS) ierr=GFC_CORRUPTED_CONT
              endif
             else
              if(top) then
               deallocate(this%container%first_elem%prev_elem); this%container%first_elem%prev_elem=>oep
              else
               deallocate(this%container%last_elem%next_elem); this%container%last_elem%next_elem=>oep
              endif
              ierr=GFC_ERROR
             endif
             lep=>NULL()
            else
             if(top) then
              this%container%first_elem%prev_elem=>oep
             else
              this%container%last_elem%next_elem=>oep
             endif
             ierr=GFC_MEM_ALLOC_FAILED
            endif
            oep=>NULL()
           else
            ierr=GFC_CORRUPTED_CONT
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         elseif(ierr.eq.GFC_IT_EMPTY) then !empty container
          if(associated(this%container)) then
           ierr=GFC_SUCCESS
           if(.not.(associated(this%current).or.&
                    &associated(this%container%first_elem).or.associated(this%container%last_elem))) then
            allocate(this%container%first_elem,STAT=errc)
            if(errc.eq.0) then
#ifdef NO_GNU
             if(present(copy_ctor_f)) then
              call this%container%first_elem%list_elem_ctor(elem_val,ierr,assoc_only=assoc,copy_ctor_f=copy_ctor_f)
             else
#endif
              call this%container%first_elem%list_elem_ctor(elem_val,ierr,assoc_only=assoc)
#ifdef NO_GNU
             endif
#endif
             if(ierr.eq.GFC_SUCCESS) then
              this%container%first_elem%prev_elem=>NULL()
              this%container%first_elem%next_elem=>NULL()
              this%container%last_elem=>this%container%first_elem
              ierr=this%reset() !reset the iterator to GFC_IT_ACTIVE
              if(this%container%num_elems_().ge.0) then !quick counting is on
               nelems=this%container%update_num_elems_(1_INTL,ierr); if(ierr.ne.GFC_SUCCESS) ierr=GFC_CORRUPTED_CONT
              endif
             else
              deallocate(this%container%first_elem); this%container%first_elem=>NULL()
              ierr=GFC_ERROR
             endif
            else
             this%container%first_elem=>NULL()
             ierr=GFC_MEM_ALLOC_FAILED
            endif
           else
            ierr=GFC_CORRUPTED_CONT
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function ListIterAppend
!-----------------------------------------------------------------------------------------------------
#ifdef NO_GNU
        function ListIterInsertElem(this,elem_val,precede,assoc_only,no_move,copy_ctor_f) result(ierr) !`GCC/5.3.0 has a bug with this
#else
        function ListIterInsertElem(this,elem_val,precede,assoc_only,no_move) result(ierr)
#endif
!Inserts a new element in the current position (either right before or right after),
!either by value or by reference. If the iterator is empty it will be set to active
!after inserting the value as the first element of the container.
         implicit none
         integer(INTD):: ierr                          !out: error code (0:success)
         class(list_iter_t), intent(inout):: this      !inout: iterator
         class(*), target, intent(in):: elem_val       !in: value to be stored
         logical, intent(in), optional:: precede       !in: TRUE:insert before, FALSE:insert after (default)
         logical, intent(in), optional:: assoc_only    !in: storage type: TRUE:by reference, FALSE:by value (default)
         logical, intent(in), optional:: no_move       !in: if TRUE, the iterator position does not change, FALSE it does move to the newly added element (default)
#ifdef NO_GNU
         procedure(gfc_copy_i), optional:: copy_ctor_f !user-defined generic copy constructor
#endif
         logical:: assoc,before,move
         integer:: errc
         integer(INTL):: nelems
         class(list_elem_t), pointer:: lep

         if(present(precede)) then; before=precede; else; before=.FALSE.; endif
         if(present(assoc_only)) then; assoc=assoc_only; else; assoc=.FALSE.; endif
         if(present(no_move)) then; move=.not.no_move; else; move=.TRUE.; endif
         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%container)) then
           if(associated(this%current)) then
            ierr=GFC_SUCCESS
            allocate(lep,STAT=errc)
            if(errc.eq.0) then
#ifdef NO_GNU
             if(present(copy_ctor_f)) then
              call lep%list_elem_ctor(elem_val,ierr,assoc_only=assoc,copy_ctor_f=copy_ctor_f)
             else
#endif
              call lep%list_elem_ctor(elem_val,ierr,assoc_only=assoc)
#ifdef NO_GNU
             endif
#endif
             if(ierr.eq.GFC_SUCCESS) then
              if(before) then
               lep%prev_elem=>this%current%prev_elem
               lep%next_elem=>this%current
               if(associated(this%current%prev_elem)) this%current%prev_elem%next_elem=>lep
               this%current%prev_elem=>lep
               if(associated(this%current,this%container%first_elem)) this%container%first_elem=>lep
              else
               lep%next_elem=>this%current%next_elem
               lep%prev_elem=>this%current
               if(associated(this%current%next_elem)) this%current%next_elem%prev_elem=>lep
               this%current%next_elem=>lep
               if(associated(this%current,this%container%last_elem)) this%container%last_elem=>lep
              endif
              if(move) this%current=>lep
              if(this%container%num_elems_().ge.0) then !quick counting is on
               nelems=this%container%update_num_elems_(1_INTL,ierr); if(ierr.ne.GFC_SUCCESS) ierr=GFC_CORRUPTED_CONT
              endif
             else
              deallocate(lep); ierr=GFC_ERROR
             endif
            else
             ierr=GFC_MEM_ALLOC_FAILED
            endif
            lep=>NULL()
           else
            ierr=GFC_CORRUPTED_CONT
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         elseif(ierr.eq.GFC_IT_EMPTY) then
          if(associated(this%container)) then
           if(.not.(associated(this%current).or.&
                   &associated(this%container%first_elem).or.associated(this%container%last_elem))) then
            ierr=GFC_SUCCESS
            allocate(this%container%first_elem,STAT=errc)
            if(errc.eq.0) then
#ifdef NO_GNU
             if(present(copy_ctor_f)) then
              call this%container%first_elem%list_elem_ctor(elem_val,ierr,assoc_only=assoc,copy_ctor_f=copy_ctor_f)
             else
#endif
              call this%container%first_elem%list_elem_ctor(elem_val,ierr,assoc_only=assoc)
#ifdef NO_GNU
             endif
#endif
             if(ierr.eq.GFC_SUCCESS) then
              this%container%first_elem%prev_elem=>NULL()
              this%container%first_elem%next_elem=>NULL()
              this%container%last_elem=>this%container%first_elem
              ierr=this%reset() !reset the iterator to GFC_IT_ACTIVE (regardless of <no_move>)
              if(this%container%num_elems_().ge.0) then !quick counting is on
               nelems=this%container%update_num_elems_(1_INTL,ierr); if(ierr.ne.GFC_SUCCESS) ierr=GFC_CORRUPTED_CONT
              endif
             else
              deallocate(this%container%first_elem); this%container%first_elem=>NULL()
              ierr=GFC_ERROR
             endif
            else
             this%container%first_elem=>NULL()
             ierr=GFC_MEM_ALLOC_FAILED
            endif
           else
            ierr=GFC_CORRUPTED_CONT
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function ListIterInsertElem
!---------------------------------------------------------------------
        function ListIterInsertList(this,sublist,precede) result(ierr)
!Inserts another list at the current list iterator position (either before or after).
!The inserted list becomes a sublist after insertion. The iterator position does not change.
         implicit none
         integer(INTD):: ierr                      !out: error code (0:success)
         class(list_iter_t), intent(inout):: this  !inout: list iterator
         class(list_bi_t), intent(inout):: sublist !in: inserted list (must not be a sublist)
         logical, intent(in), optional:: precede   !in: if TRUE the sublist will be inserted prior to the current iterator position (defaults to FALSE)
         logical:: before
         integer(INTL):: nelems,new_elems

         if(present(precede)) then; before=precede; else; before=.FALSE.; endif
         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           if(associated(sublist%first_elem).and.associated(sublist%last_elem).and.(.not.sublist%is_sublist())) then
            ierr=GFC_SUCCESS
            if(before) then !insert prior to the current position
             if(associated(this%current%prev_elem)) then
              this%current%prev_elem%next_elem=>sublist%first_elem
              sublist%first_elem%prev_elem=>this%current%prev_elem
             else
              this%container%first_elem=>sublist%first_elem
             endif
             sublist%last_elem%next_elem=>this%current
             this%current%prev_elem=>sublist%last_elem
            else !insert after the current position
             if(associated(this%current%next_elem)) then
              this%current%next_elem%prev_elem=>sublist%last_elem
              sublist%last_elem%next_elem=>this%current%next_elem
             else
              this%container%last_elem=>sublist%last_elem
             endif
             sublist%first_elem%prev_elem=>this%current
             this%current%next_elem=>sublist%first_elem
            endif
            call this%container%quick_counting_off_() !turn off quick element counting
            call sublist%quick_counting_off_() !turn off quick element counting
           else
            ierr=GFC_INVALID_ARGS
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function ListIterInsertList
!-------------------------------------------------------------------
        function ListIterMoveElem(this,another,precede) result(ierr)
!Moves a list element at the current iterator position into
!another list iterator at its current position, either prior or after it.
!Another iterator will shift to the just added element (prior or after).
         implicit none
         integer(INTD):: ierr                        !out: error code
         class(list_iter_t), intent(inout):: this    !inout: source list iterator (from)
         class(list_iter_t), intent(inout):: another !inout: destination list iterator (to)
         logical, intent(in), optional:: precede     !in: if TRUE, the list element will be moved prior, otherwise after (default)
         class(list_elem_t), pointer:: list_elem
         logical:: before

         before=.FALSE.; if(present(precede)) before=precede
         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          ierr=another%get_status()
          if(ierr.eq.GFC_IT_ACTIVE) then
           list_elem=>this%pop_(ierr)
           if(ierr.eq.GFC_SUCCESS) then
            if(before) then
             list_elem%prev_elem=>another%current%prev_elem
             if(associated(list_elem%prev_elem)) list_elem%prev_elem%next_elem=>list_elem
             list_elem%next_elem=>another%current
             another%current%prev_elem=>list_elem
             if(associated(another%current,another%container%first_elem)) another%container%first_elem=>list_elem
            else
             list_elem%next_elem=>another%current%next_elem
             if(associated(list_elem%next_elem)) list_elem%next_elem%prev_elem=>list_elem
             list_elem%prev_elem=>another%current
             another%current%next_elem=>list_elem
             if(associated(another%current,another%container%last_elem)) another%container%last_elem=>list_elem
            endif
            call another%jump_(list_elem); list_elem=>NULL()
           endif
          elseif(ierr.eq.GFC_IT_EMPTY) then !another iterator was empty
           list_elem=>this%pop_(ierr)
           if(ierr.eq.GFC_SUCCESS) then
            another%container%first_elem=>list_elem
            another%container%last_elem=>list_elem
            call another%jump_(list_elem); list_elem=>NULL()
           endif
          else
           ierr=GFC_INVALID_ARGS
          endif
         else
          ierr=GFC_INVALID_REQUEST
         endif
         return
        end function ListIterMoveElem
!-------------------------------------------------------------------------------------
        function ListIterMoveList(this,another,max_elems,num_elems_moved) result(ierr)
!Moves the entire list or its part from one iterator to another by
!appending the moved elements right after another list's current iterator
!position. Another iterator will then reposition to the last appended element.
!If <max_elems> is absent, all list elements will be moved from <this> to
!<another> and <this> list will become empty at the end. If <max_elems>
!is present, only up to that number of list elements will be moved,
!starting from the current iterator position of <this>. If the end
!of <this> list is reached before <max_elems> elements have been moved,
!no further elements will be moved. An optional <num_elems_moved> will
!reflect the exact number of list elements moved.
         implicit none
         integer(INTD):: ierr                            !out: error code
         class(list_iter_t), intent(inout):: this        !inout: source list iterator (from)
         class(list_iter_t), intent(inout):: another     !inout: destination list iterator (to)
         integer(INTD), intent(in), optional:: max_elems !in: upper limit on the number of moved elements
         integer(INTD), intent(out), optional:: num_elems_moved !out: number of list elements moved
         integer(INTD):: errc,n,nelems,num_moved

         num_moved=0; ierr=another%get_status()
         if(ierr.eq.GFC_IT_ACTIVE.or.ierr.eq.GFC_IT_EMPTY) then
          nelems=-1; if(present(max_elems)) nelems=max_elems
          if(nelems.ne.0) then
           if(nelems.lt.0) then !move all list elements
            ierr=this%reset(); ierr=this%get_status()
            do while(ierr.eq.GFC_IT_ACTIVE)
             ierr=this%move_elem(another); if(ierr.ne.GFC_SUCCESS) exit
             num_moved=num_moved+1
             ierr=this%get_status()
            enddo
            if(ierr.eq.GFC_IT_EMPTY) ierr=GFC_SUCCESS
            errc=this%reset(); if(errc.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=errc
           else !move up to <nelems> list elements
            n=nelems; ierr=this%get_status()
            do while(ierr.eq.GFC_IT_ACTIVE)
             if(this%on_last()) n=1
             ierr=this%move_elem(another); if(ierr.ne.GFC_SUCCESS) exit
             num_moved=num_moved+1
             ierr=this%get_status(); n=n-1; if(n.eq.0) exit
            enddo
            if(ierr.eq.GFC_IT_EMPTY.or.ierr.eq.GFC_IT_ACTIVE) ierr=GFC_SUCCESS
           endif
          endif
         else
          ierr=GFC_INVALID_ARGS
         endif
         if(present(num_elems_moved)) num_elems_moved=num_moved
         return
        end function ListIterMoveList
!-------------------------------------------------------------------
        function ListIterSplit(this,new_list,keep_tail) result(ierr)
!Splits the list into two parts at the current iterator position.
!Depending on <keep_tail>, either the top or the bottom part will
!stay with the original iterator, including the current element.
!The other part will be returned as a separate list.
         implicit none
         integer(INTD):: ierr                         !out: error code (0:success)
         class(list_iter_t), intent(inout):: this     !inout: list iterator (cannot be a sublist iterator)
         class(list_bi_t), intent(inout):: new_list   !out: new list (cut), must be empty on entrance
         logical, intent(in), optional:: keep_tail    !in: if TRUE, the tail part will be associated with the iterator (defaults to FALSE)
         logical:: keep_top
         class(list_elem_t), pointer:: homo,lumo

         if(present(keep_tail)) then; keep_top=.not.keep_tail; else; keep_top=.TRUE.; endif
         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           if(.not.this%container%is_sublist(ierr)) then
            if(ierr.eq.GFC_SUCCESS) then
             if(new_list%is_empty().eq.GFC_TRUE) then
              if(keep_top) then
               homo=>this%current
               lumo=>this%current%next_elem
               if(associated(lumo)) then
                new_list%first_elem=>lumo; new_list%first_elem%prev_elem=>NULL()
                new_list%last_elem=>this%container%last_elem
                this%container%last_elem=>homo; this%container%last_elem%next_elem=>NULL()
               else
                ierr=GFC_INVALID_ARGS
               endif
              else
               homo=>this%current%prev_elem
               lumo=>this%current
               if(associated(homo)) then
                new_list%first_elem=>this%container%first_elem
                new_list%last_elem=>homo; new_list%last_elem%next_elem=>NULL()
                this%container%first_elem=>lumo; this%container%first_elem%prev_elem=>NULL()
               else
                ierr=GFC_INVALID_ARGS
               endif
              endif
              if(ierr.eq.GFC_SUCCESS) then
               call this%container%quick_counting_off_() !turn off quick element counting
               call new_list%quick_counting_off_() !turn off quick element counting
              endif
              homo=>NULL(); lumo=>NULL()
             else
              ierr=GFC_INVALID_ARGS
             endif
            else
             ierr=GFC_ERROR
            endif
           else
            ierr=GFC_INVALID_ARGS
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function ListIterSplit
!------------------------------------------------------------------
        function ListIterDelete(this,dtor_f,value_out) result(ierr)
!Deletes the element in the current iterator position. The current
!iterator position moves to the preious element, unless there is none.
!In the latter case, it moves to the next element, unless there is none.
!In the latter case, the iterator/container becomes empty.
!If <value_out> is present, the list element value will not be deleted
!but moved into <value_out>. In this case, <dtor_f> will be ignored.
         implicit none
         integer(INTD):: ierr                                 !out: error code (0:success)
         class(list_iter_t), intent(inout):: this             !inout: iterator
         procedure(gfc_destruct_i), optional:: dtor_f         !in: element value destructor
         class(*), pointer, intent(out), optional:: value_out !out: pointer to the saved value of the deleted list element
         class(list_elem_t), pointer:: lep
         integer(INTD):: errc
         integer(INTL):: nelems
         logical:: first,last
         character(256):: errmesg !debug

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           ierr=GFC_SUCCESS; lep=>this%current
           first=.FALSE.; if(associated(this%current,this%container%first_elem)) first=.TRUE.
           last=.FALSE.; if(associated(this%current,this%container%last_elem)) last=.TRUE.
           if(associated(this%current%prev_elem)) then
            this%current%prev_elem%next_elem=>this%current%next_elem
            if(associated(this%current%next_elem)) then
             this%current%next_elem%prev_elem=>this%current%prev_elem
            endif
           else
            if(associated(this%current%next_elem)) this%current%next_elem%prev_elem=>NULL()
           endif
           if(this%container%num_elems_().ge.0) then !quick counting is on
            nelems=this%container%update_num_elems_(-1_INTL,errc); if(errc.ne.GFC_SUCCESS) ierr=NOT_CLEAN
           endif
           if(first.and.last) then !the only element
            this%current=>NULL(); this%container%first_elem=>NULL(); this%container%last_elem=>NULL()
            errc=this%reset(); if(errc.ne.GFC_SUCCESS) ierr=NOT_CLEAN
           else
            if(first) then
             this%current=>this%current%next_elem
             this%container%first_elem=>this%current
            else
             this%current=>this%current%prev_elem
             if(last) this%container%last_elem=>this%current
            endif
           endif
           if(present(value_out)) then
            call lep%destruct(errc,value_out=value_out)
           else
            if(present(dtor_f)) then
             call lep%destruct(errc,dtor_f)
            else
             call lep%destruct(errc)
            endif
           endif
           if(errc.ne.GFC_SUCCESS) then
            if(VERBOSE) write(CONS_OUT,*)'GFC::list:ListIterDelete: element value destruction failed with error ',errc !debug
            if(ierr.eq.GFC_SUCCESS) ierr=NOT_CLEAN
           endif
           deallocate(lep,STAT=errc,ERRMSG=errmesg)
           if(errc.ne.0) then
            if(VERBOSE) write(CONS_OUT,*)'GFC::list:ListIterDelete: deallocate() failed: ',errc,': '//&
                        &errmesg(1:len_trim(errmesg)) !debug
            if(ierr.eq.GFC_SUCCESS) ierr=NOT_CLEAN
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function ListIterDelete
!-------------------------------------------------------------------------
        function ListIterDeleteSublist(this,dtor_f,backwards) result(ierr)
!Deletes all elements in the list starting from the current iterator position (either prior or after).
         implicit none
         integer(INTD):: ierr                         !out: error code (0:success)
         class(list_iter_t), intent(inout):: this     !inout: list iterator
         procedure(gfc_destruct_i), optional:: dtor_f !in: element value destructor
         logical, intent(in), optional:: backwards    !in: if TRUE, all preceding elements will be deleted instead of subsequent ones (defaults to FALSE)
         integer(INTD):: errc
         logical:: back
         class(list_elem_t), pointer:: lep

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          ierr=GFC_SUCCESS
          if(associated(this%current)) then
           if(present(backwards)) then; back=backwards; else; back=.FALSE.; endif
           lep=>this%current
           if(back) then
            this%current=>this%container%first_elem
           else
            this%current=>this%container%last_elem
           endif
           if(associated(this%current)) then
            if(present(dtor_f)) then
             do while(.not.associated(this%current,lep))
              errc=this%delete(dtor_f); if(errc.ne.GFC_SUCCESS) ierr=NOT_CLEAN
             enddo
             errc=this%delete(dtor_f); if(errc.ne.GFC_SUCCESS) ierr=NOT_CLEAN
            else
             do while(.not.associated(this%current,lep))
              errc=this%delete(); if(errc.ne.GFC_SUCCESS) ierr=NOT_CLEAN
             enddo
             errc=this%delete(); if(errc.ne.GFC_SUCCESS) ierr=NOT_CLEAN
            endif
            lep=>NULL()
           else
            ierr=GFC_CORRUPTED_CONT
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function ListIterDeleteSublist
!-----------------------------------------------------------
        function ListIterDeleteAll(this,dtor_f) result(ierr)
!Deletes all list elements leaving it in an empty state.
         implicit none
         integer(INTD):: ierr                         !out: error code
         class(list_iter_t), intent(inout):: this     !inout: list iterator
         procedure(gfc_destruct_i), optional:: dtor_f !in: element value destructor

         ierr=this%get_status()
         if(ierr.ne.GFC_IT_NULL) then
          ierr=this%reset()
          if(ierr.eq.GFC_SUCCESS) then
           if(present(dtor_f)) then
            ierr=this%delete_sublist(dtor_f)
           else
            ierr=this%delete_sublist()
           endif
           if(ierr.eq.GFC_IT_EMPTY) ierr=GFC_SUCCESS
          endif
         endif
         return
        end function ListIterDeleteAll
!------------------------------------------------------------
        function ListIterBookmark(this,bookmark) result(ierr)
!Bookmarks the current iterator position.
         implicit none
         integer(INTD):: ierr                        !out: error code
         class(list_iter_t), intent(in):: this       !in: list iterator
         class(list_pos_t), intent(out):: bookmark   !out: bookmarked list position

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          ierr=GFC_SUCCESS
          bookmark%elem_p=>this%current
         endif
         return
        end function ListIterBookmark
!--------------------------------------------------------
        function ListIterJump(this,bookmark) result(ierr)
!Moves the iterator to a previously bookmarked iterator position.
         implicit none
         integer(INTD):: ierr                     !out: error code
         class(list_iter_t), intent(inout):: this !inout: list iterator
         class(list_pos_t), intent(in):: bookmark !in: bookmark

         ierr=GFC_SUCCESS
         if(associated(bookmark%elem_p)) then
          call this%jump_(bookmark%elem_p)
         else
          ierr=GFC_NO_MOVE
         endif
         return
        end function ListIterJump
!----------------------------------------------
        subroutine ListIterJump_(this,new_elem)
!Moves the iterator to an arbitrary specified list element.
         implicit none
         class(list_iter_t), intent(inout):: this           !inout: list iterator
         class(list_elem_t), pointer, intent(in):: new_elem !in: pointer to the new element or NULL()
         integer(INTD):: errc,sts

         if(associated(this%current)) call this%current%decr_ref_()
         this%current=>new_elem
         if(associated(this%current)) then
          call this%current%incr_ref_()
          errc=this%set_status_(GFC_IT_ACTIVE)
         else
          if(this%container%is_empty().eq.GFC_TRUE) then
           errc=this%set_status_(GFC_IT_EMPTY)
          else
           errc=this%set_status_(GFC_IT_DONE)
          endif
         endif
         return
        end subroutine ListIterJump_
!------------------------------------------------------------
        function ListIterPop_(this,ierr,to_prev) result(elem)
!Pops up the element at the current iterator position without destroying it.
!By default, the current iterator position will try to move to the next element,
!unless <to_prev> says otherwise.
         implicit none
         class(list_elem_t), pointer:: elem          !out: owning pointer to the current element (detached from the list)
         class(list_iter_t), intent(inout):: this    !inout: list iterator
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: to_prev     !in: if TRUE, <this> iterator will move to the previous element (defaults to FALSE)
         class(list_elem_t), pointer:: list_elem_null=>NULL()
         integer(INTD):: errc
         logical:: prev

         elem=>NULL(); errc=this%get_status()
         if(errc.eq.GFC_IT_ACTIVE) then
          prev=.FALSE.; if(present(to_prev)) prev=to_prev
          errc=GFC_SUCCESS; elem=>this%current
          if(associated(elem,this%container%first_elem)) then
           if(.not.associated(elem,this%container%last_elem)) then
            this%container%first_elem=>elem%next_elem
            call this%jump_(elem%next_elem)
           else
            this%container%first_elem=>NULL(); this%container%last_elem=>NULL()
            call this%jump_(list_elem_null); errc=this%reset()
           endif
          else
           if(associated(elem,this%container%last_elem)) then
            this%container%last_elem=>elem%prev_elem
            call this%jump_(elem%prev_elem)
           else
            if(prev) then
             call this%jump_(elem%prev_elem)
            else
             call this%jump_(elem%next_elem)
            endif
           endif
          endif
          if(associated(elem%prev_elem)) elem%prev_elem%next_elem=>elem%next_elem
          if(associated(elem%next_elem)) elem%next_elem%prev_elem=>elem%prev_elem
          elem%prev_elem=>NULL(); elem%next_elem=>NULL()
         else
          errc=GFC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function ListIterPop_

       end module gfc_list
