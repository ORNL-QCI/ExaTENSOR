!Generic Fortran Containers:: Linked list.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2016-03-07 (started 2016-02-28)
!Copyright (C) 2016 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2016 Oak Ridge National Laboratory (UT-Battelle)
!LICENSE: GNU GPL v.2 (or higher).
!NOTES:
! # A list is a linked derivative of an abstract (unlinked) GFC container.
!   A sublist is a list incorporated into another list.
! # All accesses, updates, and scans on a list are performed via
!   the list iterator associated with the list. Elements of a sublist
!   can be accessed by both the original and the combined iterators.
!   Multiple iterators can be associated with a list at the same time.
!FOR DEVELOPERS:
! * Currently, if an element is added/deleted via a sublist iterator
!   the total number of elements is not updated in the containing list,
!   and vice versa!
       module list
        use gfc_base
        use timers
        implicit none
        private
!PARAMETERS:
!Basic:
        integer(INTD), private:: CONS_OUT=6 !output device
        logical, private:: VERBOSE=.true.   !verbositiy for errors
        integer(INTD), private:: DEBUG=0    !debugging level (0:none)
!TYPES:
 !Linked list element:
        type, extends(gfc_cont_elem_t), public:: list_elem_t
         class(list_elem_t), pointer, private:: next_elem=>NULL()
         class(list_elem_t), pointer, private:: prev_elem=>NULL()
        contains
         procedure, public:: is_first=>ListElemIsFirst !returns GFC_TRUE if the element is the first in the list
         procedure, public:: is_last=>ListElemIsLast   !returns GFC_TRUE if the element is the last in the list
        end type list_elem_t
 !Linked list:
        type, extends(gfc_container_t), public:: list_bi_t
         class(list_elem_t), pointer, private:: first_elem=>NULL() !first element of the linked list
         class(list_elem_t), pointer, private:: last_elem=>NULL()  !last element of the linked list
        contains
         procedure, public:: is_sublist=>ListIsSublist             !returns TRUE if the list is a sublist of a larger list, FALSE otherwise
        end type list_bi_t
 !List iterator:
        type, extends(gfc_iter_t), public:: list_iter_t
         class(list_elem_t), pointer, private:: current=>NULL() !current element of the linked list
         class(list_bi_t), pointer, private:: container=>NULL() !linked list associated with the iterator
        contains
         procedure, public:: init=>ListIterInit                 !initializes the iterator by associating it with a list and resetting to the beginning
         procedure, public:: reset=>ListIterReset               !resets the iterator to the beginning of the list
         procedure, public:: reset_back=>ListIterResetBack      !resets the iterator to the end of the list
         procedure, public:: release=>ListIterRelease           !releases the iterator (dissociates it from the container)
         procedure, public:: pointee=>ListIterPointee           !returns the container element currently pointed to by the iterator
         procedure, public:: next=>ListIterNext                 !moves the iterator to the next list element
         procedure, public:: previous=>ListIterPrevious         !moves the iterator to the previous list element
         procedure, public:: append=>ListIterAppend             !inserts a new element either at the beginning or at the end of the container
         procedure, public:: insert_elem=>ListIterInsertElem    !inserts a new element at the current position of the container
         procedure, public:: insert_list=>ListIterInsertList    !inserts another linked list at the current position of the container
         generic, public:: insert=>ListIterInsertElem,ListIterInsertList !generic
         procedure, public:: split=>ListIterSplit               !splits the list into two parts at the current position of the container
         procedure, public:: delete=>ListIterDelete             !deletes an element or multiple elements starting from the current position
        end type list_iter_t
!INTERFACES:
!VISIBILITY:
        private ListElemIsFirst
        private ListElemIsLast
        private ListIsSublist
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

       contains
!IMPLEMENTATION:
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
!----------------------------------------------------
        function ListIsSublist(this,ierr) result(res)
!Returns TRUE if the list is a sublist of another list, FALSE otherwise.
         implicit none
         logical:: res                               !out: result
         class(list_bi_t), intent(in):: this         !in: list
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc

         errc=GFC_SUCCESS; res=.false.
         if(associated(this%first_elem).and.associated(this%last_elem)) then
          if(associated(this%first_elem%prev_elem).or.associated(this%last_elem%next_elem)) res=.true.
         else
          errc=GFC_EMPTY_CONT
         endif
         if(present(ierr)) ierr=errc
         return
        end function ListIsSublist
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
           ierr=this%set_status(GFC_IT_ACTIVE)
          else
           ierr=this%set_status(GFC_IT_EMPTY)
          endif
          call this%reset_count() !reset all iteration counters
         else
          ierr=this%set_status(GFC_IT_NULL)
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
           ierr=this%set_status(GFC_IT_ACTIVE)
          else
           ierr=this%set_status(GFC_IT_EMPTY)
          endif
          call this%reset_count() !reset all iteration counters
         else
          ierr=this%set_status(GFC_IT_NULL)
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
         call this%reset_count(); ierr=this%set_status(GFC_IT_NULL)
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
            ierr=this%set_status(GFC_IT_DONE)
            lep=>NULL()
           endif
           if(.not.associated(lep)) then; ierr=GFC_IT_DONE; else; lep=>NULL(); endif
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
            ierr=this%set_status(GFC_IT_DONE)
            lep=>NULL()
           endif
           if(.not.associated(lep)) then; ierr=GFC_IT_DONE; else; lep=>NULL(); endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function ListIterPrevious
!------------------------------------------------------------------------------------------
#ifdef NO_GNU
        function ListIterAppend(this,elem_val,at_top,assoc_only,copy_constr_f) result(ierr) !`GCC/5.3.0 has a bug with this
#else
        function ListIterAppend(this,elem_val,at_top,assoc_only) result(ierr)
#endif
!Appends an element at the beginning or at the end of the list, either by value or by reference.
!For non-empty iterators, the iterator position is kept unchanged (even if it is GFC_IT_DONE).
!If the iterator is empty on input, its status will be changed to GFC_IT_ACTIVE (via reset).
         implicit none
         integer(INTD):: ierr                            !out: error code (0:success)
         class(list_iter_t), intent(inout):: this        !inout: iterator
         class(*), target, intent(in):: elem_val         !in: value to be stored
         logical, intent(in), optional:: at_top          !in: TRUE:append at the top, FALSE:append at the end (default)
         logical, intent(in), optional:: assoc_only      !in: storage type: TRUE:by reference, FALSE:by value (default)
#ifdef NO_GNU
         procedure(gfc_copy_i), optional:: copy_constr_f !user-defined generic copy constructor
#endif
         logical:: assoc,top
         integer:: errc
         integer(INTL):: nelems
         class(list_elem_t), pointer:: lep

         if(present(at_top)) then; top=at_top; else; top=.false.; endif
         if(present(assoc_only)) then; assoc=assoc_only; else; assoc=.false.; endif
         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE.or.ierr.eq.GFC_IT_DONE) then !non-empty container
          if(associated(this%container)) then
           if(associated(this%container%first_elem).and.associated(this%container%last_elem)) then
            ierr=GFC_SUCCESS
            if(top) then
             allocate(this%container%first_elem%prev_elem,STAT=errc)
             if(errc.eq.0) lep=>this%container%first_elem%prev_elem
            else
             allocate(this%container%last_elem%next_elem,STAT=errc)
             if(errc.eq.0) lep=>this%container%last_elem%next_elem
            endif
            if(errc.eq.0) then
#ifdef NO_GNU
             if(present(copy_constr_f)) then
              call lep%construct(elem_val,ierr,assoc_only=assoc,copy_constr_func=copy_constr_f)
             else
              call lep%construct(elem_val,ierr,assoc_only=assoc)
             endif
#else
             call lep%construct(elem_val,ierr,assoc_only=assoc)
#endif
             if(ierr.eq.GFC_SUCCESS) then
              if(top) then
               lep%next_elem=>this%container%first_elem
               lep%prev_elem=>NULL()
               this%container%first_elem=>lep
              else
               lep%prev_elem=>this%container%last_elem
               lep%next_elem=>NULL()
               this%container%last_elem=>lep
              endif
              nelems=this%container%update_num_elems_(1_INTL,ierr); if(ierr.ne.GFC_SUCCESS) ierr=GFC_CORRUPTED_CONT
             else
              if(top) then
               deallocate(this%container%first_elem%prev_elem); this%container%first_elem%prev_elem=>NULL()
              else
               deallocate(this%container%last_elem%next_elem); this%container%last_elem%next_elem=>NULL()
              endif
              ierr=GFC_ERROR
             endif
             lep=>NULL()
            else
             if(top) then
              this%container%first_elem%prev_elem=>NULL()
             else
              this%container%last_elem%next_elem=>NULL()
             endif
             ierr=GFC_MEM_ALLOC_FAILED
            endif
           else
            ierr=GFC_CORRUPTED_CONT
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         elseif(ierr.eq.GFC_IT_EMPTY) then !empty container
          if(associated(this%container)) then
           ierr=GFC_SUCCESS
           if(.not.(associated(this%container%first_elem).or.associated(this%container%last_elem))) then
            allocate(this%container%first_elem,STAT=errc)
            if(errc.eq.0) then
#ifdef NO_GNU
             if(present(copy_constr_f)) then
              call this%container%first_elem%%construct(elem_val,ierr,assoc_only=assoc,copy_constr_func=copy_constr_f)
             else
              call this%container%first_elem%construct(elem_val,ierr,assoc_only=assoc)
             endif
#else
             call this%container%first_elem%construct(elem_val,ierr,assoc_only=assoc)
#endif
             if(ierr.eq.GFC_SUCCESS) then
              this%container%first_elem%prev_elem=>NULL()
              this%container%first_elem%next_elem=>NULL()
              this%container%last_elem=>this%container%first_elem
              ierr=this%reset() !reset the iterator to GFC_IT_ACTIVE
              nelems=this%container%update_num_elems_(1_INTL,ierr); if(ierr.ne.GFC_SUCCESS) ierr=GFC_CORRUPTED_CONT
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
!-------------------------------------------------------------------------------------------------------
#ifdef NO_GNU
        function ListIterInsertElem(this,elem_val,precede,assoc_only,no_move,copy_constr_f) result(ierr) !`GCC/5.3.0 has a bug with this
#else
        function ListIterInsertElem(this,elem_val,precede,assoc_only,no_move) result(ierr)
#endif
!Inserts a new element in the current position (either right before or right after),
!either by value or by reference. If the iterator is empty it will be set to active
!after inserting the value as the first element of the container.
         implicit none
         integer(INTD):: ierr                            !out: error code (0:success)
         class(list_iter_t), intent(inout):: this        !inout: iterator
         class(*), target, intent(in):: elem_val         !in: value to be stored
         logical, intent(in), optional:: precede         !in: TRUE:insert before, FALSE:insert after (default)
         logical, intent(in), optional:: assoc_only      !in: storage type: TRUE:by reference, FALSE:by value (default)
         logical, intent(in), optional:: no_move         !in: if TRUE, the iterator position does not change, FALSE it does move to the newly added element
#ifdef NO_GNU
         procedure(gfc_copy_i), optional:: copy_constr_f !user-defined generic copy constructor
#endif
         logical:: assoc,before,move
         integer:: errc
         integer(INTL):: nelems
         class(list_elem_t), pointer:: lep

         if(present(precede)) then; before=precede; else; before=.false.; endif
         if(present(assoc_only)) then; assoc=assoc_only; else; assoc=.false.; endif
         if(present(no_move)) then; move=.not.no_move; else; move=.true.; endif
         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%container)) then
           if(associated(this%current)) then
            ierr=GFC_SUCCESS
            allocate(lep,STAT=errc)
            if(errc.eq.0) then
#ifdef NO_GNU
             if(present(copy_constr_f)) then
              call lep%construct(elem_val,ierr,assoc_only=assoc,copy_constr_func=copy_constr_f)
             else
              call lep%construct(elem_val,ierr,assoc_only=assoc)
             endif
#else
             call lep%construct(elem_val,ierr,assoc_only=assoc)
#endif
             if(ierr.eq.GFC_SUCCESS) then
              if(before) then
               lep%prev_elem=>this%current%prev_elem
               lep%next_elem=>this%current
               if(associated(this%current%prev_elem)) this%current%prev_elem%next_elem=>lep
               this%current%prev_elem=>lep
              else
               lep%next_elem=>this%current%next_elem
               lep%prev_elem=>this%current
               if(associated(this%current%next_elem)) this%current%next_elem%prev_elem=>lep
               this%current%next_elem=>lep
              endif
              if(move) this%current=>lep
              nelems=this%container%update_num_elems_(1_INTL,ierr); if(ierr.ne.GFC_SUCCESS) ierr=GFC_CORRUPTED_CONT
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
             if(present(copy_constr_f)) then
              call this%container%first_elem%%construct(elem_val,ierr,assoc_only=assoc,copy_constr_func=copy_constr_f)
             else
              call this%container%first_elem%construct(elem_val,ierr,assoc_only=assoc)
             endif
#else
             call this%container%first_elem%construct(elem_val,ierr,assoc_only=assoc)
#endif
             if(ierr.eq.GFC_SUCCESS) then
              this%container%first_elem%prev_elem=>NULL()
              this%container%first_elem%next_elem=>NULL()
              this%container%last_elem=>this%container%first_elem
              ierr=this%reset() !reset the iterator to GFC_IT_ACTIVE (regardless of <no_move>)
              nelems=this%container%update_num_elems_(1_INTL,ierr); if(ierr.ne.GFC_SUCCESS) ierr=GFC_CORRUPTED_CONT
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
!Inserts another list at the current iterator position (either before or after).
!The inserted list becomes a sublist after insertion.
         implicit none
         integer(INTD):: ierr                     !out: error code (0:success)
         class(list_iter_t), intent(inout):: this !inout: list iterator
         class(list_bi_t), intent(in):: sublist   !in: inserted list
         logical, intent(in), optional:: precede  !in: if TRUE the sublist will be inserted prior to the current iterator position (defaults to FALSE)
         logical:: before
         integer(INTL):: nelems,new_elems

         if(present(precede)) then; before=precede; else; before=.false.; endif
         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           if(associated(sublist%first_elem).and.associated(sublist%last_elem)) then
            new_elems=sublist%num_elems(ierr)
            if(ierr.eq.GFC_SUCCESS.and.new_elems.gt.0) then
             if(before) then !insert prior to the current position
              if(associated(this%current%prev_elem)) then
               this%current%prev_elem%next_elem=>sublist%first_elem
               sublist%first_elem%prev_elem=>this%current%prev_elem
               sublist%last_elem%next_elem=>this%current
               this%current%prev_elem=>sublist%last_elem
              else
               this%container%first_elem=>sublist%first_elem
               sublist%last_elem%next_elem=>this%current
               this%current%prev_elem=>sublist%last_elem
              endif
             else !insert after the current position
              if(associated(this%current%next_elem)) then
               this%current%next_elem%prev_elem=>sublist%last_elem
               sublist%last_elem%next_elem=>this%current%next_elem
               sublist%first_elem%prev_elem=>this%current
               this%current%next_elem=>sublist%first_elem
              else
               this%container%last_elem=>sublist%last_elem
               sublist%first_elem%prev_elem=>this%current
               this%current%next_elem=>sublist%first_elem
              endif
             endif
             nelems=this%container%update_num_elems_(new_elems,ierr); if(ierr.ne.GFC_SUCCESS) ierr=GFC_CORRUPTED_CONT
            else
             ierr=GFC_ERROR
            endif
           else
            ierr=GFC_EMPTY_CONT
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function ListIterInsertList
!--------------------------------------------------------------------------------
        function ListIterSplit(this,new_list,keep_tail,exclude_iter) result(ierr)
!Splits the list into two parts at the current iterator position. Depending on
!the value of <exclude_iter>, the iterator will either stay in the first part
!of the list or in the second. One of the parts will be associated with the present
!iterator, another part will be returned as a separate list. After splitting,
!the iterator position is reset to the top of the list.
         implicit none
         integer(INTD):: ierr                         !out: error code (0:success)
         class(list_iter_t), intent(inout):: this     !inout: list iterator
         class(list_bi_t), intent(out):: new_list     !out: new list (cut)
         logical, intent(in), optional:: keep_tail    !in: if TRUE, the tail part will be associated with the iterator (defaults to FALSE)
         logical, intent(in), optional:: exclude_iter !in: if TRUE, the current element of the iterator will be excluded from the iterator (defaults to FALSE)
         logical:: keep_top,include_iter
         class(list_elem_t), pointer:: last1,first2

         if(present(keep_tail)) then; keep_top=.not.keep_tail; else; keep_top=.true.; endif
         if(present(exclude_iter)) then; include_iter=.not.exclude_iter; else; include_iter=.true.; endif
         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           if(exclude_iter) then
            last1=>this%current%prev_elem
            first2=>this%current
           else
            last1=>this%current
            first2=>this%current%next_elem
           endif
           if(keep_top) then
            if(
           else
            
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function ListIterSplit

       end module list
