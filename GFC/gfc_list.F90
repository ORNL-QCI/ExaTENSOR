!Generic Fortran Containers (GFC): Bi-directional linked list
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2017-05-11 (started 2016-02-28)

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
 !Linked list:
        type, extends(gfc_container_t), public:: list_bi_t
         class(list_elem_t), pointer, private:: first_elem=>NULL() !first element of the linked list
         class(list_elem_t), pointer, private:: last_elem=>NULL()  !last element of the linked list
        contains
         procedure, public:: is_empty=>ListIsEmpty                 !returns GFC_TRUE if the list is empty, GFC_FALSE otherwise (or error code)
         procedure, public:: is_sublist=>ListIsSublist             !returns TRUE if the list is a sublist of a larger list, FALSE otherwise
         procedure, private:: ListDuplicateToList                  !duplicates a list into another list either by value or by reference
         procedure, private:: ListDuplicateToVector                !duplicates a list into a vector either by value or by reference
         generic, public:: duplicate=>ListDuplicateToList,ListDuplicateToVector
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
         procedure, public:: append=>ListIterAppend               !inserts a new element either at the beginning or at the end of the container
         procedure, public:: insert_elem=>ListIterInsertElem      !inserts a new element at the current position of the container
         procedure, public:: insert_list=>ListIterInsertList      !inserts another linked list at the current position of the container
!        generic, public:: insert=>insert_elem,insert_list        !generic (`Ambiguous)
         procedure, public:: split=>ListIterSplit                 !splits the list into two parts at the current position
         procedure, public:: delete=>ListIterDelete               !deletes the list element in the current position
         procedure, public:: delete_sublist=>ListIterDeleteSublist!deletes all elements either prior or after the current iterator position
         procedure, public:: delete_all=>ListIterDeleteAll        !deletes all list elements
         procedure, public:: jump_=>ListIterJump                  !PRIVATE: moves the iterator to a specified list element
        end type list_iter_t
!INTERFACES:
!VISIBILITY:
 !list_elem_t:
        private ListElemConstruct
        private ListElemIsFirst
        private ListElemIsLast
 !list_bi_t:
        private ListIsEmpty
        private ListIsSublist
        private ListDuplicateToList
        private ListDuplicateToVector
 !list_iter_t:
        private ListIterInit
        private ListIterReset
        private ListIterResetBack
        private ListIterRelease
        private ListIterPointee
        private ListIterNext
        private ListIterPrevious
        private ListIterAppend
        private ListIterInsertElem
        private ListIterInsertList
        private ListIterSplit
        private ListIterDelete
        private ListIterDeleteSublist
        private ListIterDeleteAll
        private ListIterJump

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
        !class(*), target, intent(in):: obj            !in: value to be stored
         class(*), pointer, intent(in):: obj           !in: value to be stored
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
         type(list_iter_t):: ilit,ovit
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
!--------------------------------------------------------
        function ListIterDelete(this,dtor_f) result(ierr)
!Deletes the element in the current iterator position. The current
!iterator position moves to the preious element, unless there is none.
!In the latter case, it moves to the next element, unless there is none.
!In the latter case, the iterator/container becomes empty.
         implicit none
         integer(INTD):: ierr                         !out: error code (0:success)
         class(list_iter_t), intent(inout):: this     !inout: iterator
         procedure(gfc_destruct_i), optional:: dtor_f !in: element value destructor
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
            else
             this%current=>this%current%prev_elem
            endif
           endif
           if(present(dtor_f)) then
            call lep%destruct(errc,dtor_f)
           else
            call lep%destruct(errc)
           endif
           if(errc.ne.0) then
            write(CONS_OUT,*)'GFC::list:ListIterDelete: element value destruction failed with error ',errc !debug
            ierr=NOT_CLEAN
           endif
           deallocate(lep,STAT=errc,ERRMSG=errmesg)
           if(errc.ne.0) then
            write(CONS_OUT,*)'GFC::list:ListIterDelete: deallocate(lep) failed: ',errc,': '//errmesg(1:len_trim(errmesg)) !debug
            ierr=NOT_CLEAN
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
            this%current=>this%container%last_elem
           else
            this%current=>this%container%first_elem
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

         ierr=this%reset()
         if(ierr.eq.GFC_SUCCESS) then
          if(present(dtor_f)) then
           ierr=this%delete_sublist(dtor_f)
          else
           ierr=this%delete_sublist()
          endif
          if(ierr.eq.GFC_IT_EMPTY) ierr=GFC_SUCCESS
         endif
         return
        end function ListIterDeleteAll
!---------------------------------------------
        subroutine ListIterJump(this,new_elem)
!Moves the iterator to an arbitrary specified list element.
         implicit none
         class(list_iter_t), intent(inout):: this              !inout: list iterator
         class(list_elem_t), pointer, intent(inout):: new_elem !in: pointer to the new element or NULL()
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
        end subroutine ListIterJump

       end module gfc_list
