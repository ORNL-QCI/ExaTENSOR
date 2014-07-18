       module lists
!Realizations of linked lists.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2014/07/18
!DESCRIPTION:
! CLASS(list_two_way_t):
!  NOTES:
!  # To avoid frequent memory allocates/deallocates, this class assumes that there exists
!    a known-size preallocated 1D array of objects that will contain the actualy data.
!    The class itself only supplies the infrastructure needed to organize the entries of that
!    1D array into a linked bidirectional list. Namely, it manipulates solely the entry numbers,
!    not the data stored in those entries. That is, an entry can be added to the linked list,
!    the next/previous linked entry can be retrieved, the current entry can be
!    set/reset/retrieved, an entry can be deleted from the linked list, etc.
!    Thus, this class does not deal with the actual data stored at all. It only helps to
!    organize an existing 1D array into a linked bidirectional list. Moreover, for each
!    existing 1D array multiple linked lists can be created, reflecting different
!    orderings of the array entries needed for specific purposes.
!  # Maximal length of a linked list can be zero (always empty list).
!    A negative max length means an uninitialized list.
        use, intrinsic:: ISO_C_BINDING, only: C_INT
!PARAMETERS:
        integer, parameter, private:: int_kind=C_INT                       !default integer in this module
        integer(int_kind), parameter, public:: list_err_invalid_arg=1      !invalid argument passed
        integer(int_kind), parameter, public:: list_err_mem_alloc_failed=2 !memory allocation failed
        integer(int_kind), parameter, public:: list_err_mem_free_failed=3  !memory deallocation failed
        integer(int_kind), parameter, public:: list_err_list_null=4        !list has not been initialized yet
        integer(int_kind), parameter, public:: list_err_list_exists=5      !list has already been initialized
        integer(int_kind), parameter, public:: list_err_list_corrupted=6   !list linking corrupted
        integer(int_kind), parameter, public:: list_err_entry_dead=7       !attempt to use a free (dead) entry (the one not in the list)
        integer(int_kind), parameter, public:: list_err_entry_busy=8       !attempt to add an entry which is already in the list
        integer(int_kind), parameter, public:: list_empty=9                !initialized list is empty
        integer(int_kind), parameter, public:: list_full=10                !initialized list is full
        integer(int_kind), parameter, public:: list_end=11                 !end of an initialized list (in any direction)
        logical, parameter, public:: go_to_the_top=.false.
        logical, parameter, public:: go_to_the_end=.true.
        logical, parameter, public:: append_before=.true.
        logical, parameter, public:: append_after=.false.
!TYPES:
 !Bidirectional linked list (linked in both ways):
        type, public:: list_two_way_t
         integer(int_kind), private:: list_max_length=-1   !max length of the list: range=[base_offset:base_offset+list_max_length-1]
         integer(int_kind), private:: base_offset          !lowest bound in the numeration of entries of the served 1D array
         integer(int_kind), private:: list_length          !current length of the linked list
         integer(int_kind), private:: first                !first entry number in the list (can be any from the above range)
         integer(int_kind), private:: last                 !last entry number in the list (can be any from the above range)
         integer(int_kind), private:: current              !current entry number in the list
         integer(int_kind), private:: free_ffe             !first free entry
         integer(int_kind), allocatable, private:: next(:) !next entry (negative: dead entry; positive: active entry)
         integer(int_kind), allocatable, private:: prev(:) !previous entry (negative: dead entry; positive: active entry)
         contains
          procedure, public:: create=>list_create_two_way      !create a bidirectional linked list
          procedure, public:: destroy=>list_destroy_two_way    !destroy a bidirectional linked list
          procedure, public:: max_length=>list_get_max_length  !get the max length of the list (length of the array served)
          procedure, public:: length=>list_get_length          !get the current length of the list
          procedure, public:: reset=>list_reset                !reset the current position in the list
          procedure, public:: test_first=>list_test_first      !test whether current == first
          procedure, public:: test_last=>list_test_last        !test whether current == last
          procedure, public:: set_position=>list_set_position  !set the current list position to a given entry
          procedure, public:: get_position=>list_get_position  !get the current position in the list
          procedure, public:: go_next=>list_move_to_next       !move to the next entry in the linked list
          procedure, public:: go_previous=>list_move_to_prev   !move to the previous entry in the linked list
          procedure, public:: add=>list_add_item_two_way       !add a new entry in the linked list
          procedure, public:: delete=>list_delete_item_two_way !delete an entry from the linked list
          procedure, public:: test_list=>list_test_two_way     !test the correctness of the linked list
        end type list_two_way_t
!MODULE PROCEDURES:
        private list_create_two_way
        private list_destroy_two_way
        private list_get_max_length
        private list_get_length
        private list_reset
        private list_test_first
        private list_test_last
        private list_set_position
        private list_get_position
        private list_move_to_next
        private list_move_to_prev
        private list_add_item_two_way
        private list_delete_item_two_way
        private list_test_two_way

       contains
!---------------------------------------------------------------------------------
        integer(int_kind) function list_create_two_way(this,list_max_len,list_beg)
!This function creates a bidirectional linked list infrastructure for some external 1D array.
!INPUT:
! # list_max_len - max length of the list (the length of the external 1D array served);
! # list_beg - lowest bound of the external 1D array: range of entry numbers = [list_beg:list_beg+list_max_len-1];
        implicit none
        class(list_two_way_t):: this
        integer(int_kind), intent(in):: list_max_len,list_beg
        integer(int_kind):: i
        list_create_two_way=0
        if(list_max_len.ge.0.and.list_max_len.lt.huge(i).and.huge(i)-list_max_len.ge.list_beg) then !check arguments
         if(this%list_max_length.lt.0) then !list is uninitialized
          this%base_offset=list_beg; this%list_length=0; this%first=-1; this%last=-1; this%current=-1
          if(list_max_len.gt.0) then
           allocate(this%next(1:list_max_len),STAT=i)
           if(i.ne.0) then; list_create_two_way=list_err_mem_alloc_failed; return; endif
           do i=1,list_max_len; this%next(i)=-(i+1); enddo; this%free_ffe=1
           allocate(this%prev(1:list_max_len),STAT=i)
           if(i.ne.0) then; deallocate(this%next); list_create_two_way=list_err_mem_alloc_failed; return; endif
           do i=2,list_max_len; this%prev(i)=-(i-1); enddo; this%prev(1)=-(list_max_len+1)
           this%list_max_length=list_max_len
          else !initialize a trivial zero-max-length list (always empty)
           this%list_max_length=0; this%free_ffe=-1
          endif
         else !list has already been initialized
          list_create_two_way=list_err_list_exists
         endif
        else
         list_create_two_way=list_err_invalid_arg
        endif
        return
        end function list_create_two_way
!------------------------------------------------------------
        integer(int_kind) function list_destroy_two_way(this)
!This function destroys a bidirectional linked list infrastructure.
        implicit none
        class(list_two_way_t):: this
        integer(int_kind):: i
        list_destroy_two_way=0
        if(this%list_max_length.ge.0) then !initialized list
         this%list_max_length=-1; this%base_offset=-1; this%list_length=-1
         this%first=-1; this%last=-1; this%current=-1; this%free_ffe=-1
         if(allocated(this%next)) then
          deallocate(this%next,STAT=i); if(i.ne.0) list_destroy_two_way=list_err_mem_free_failed
         endif
         if(allocated(this%prev)) then
          deallocate(this%prev,STAT=i); if(i.ne.0) list_destroy_two_way=list_err_mem_free_failed
         endif
        else !uninitialized list
         list_destroy_two_way=list_err_list_null
        endif
        return
        end function list_destroy_two_way
!-----------------------------------------------------------
        integer(int_kind) function list_get_max_length(this)
!This function returns the maximal length of a list (negative will mean an uninitialized list).
        implicit none
        class(list_two_way_t):: this
        list_get_max_length=this%list_max_length
        return
        end function list_get_max_length
!-------------------------------------------------------
        integer(int_kind) function list_get_length(this)
!This function returns the current length of a list (negative will mean an uninitialized list).
        implicit none
        class(list_two_way_t):: this
        if(this%list_max_length.ge.0) then !initialized list
         list_get_length=this%list_length
        else !uninitialized list
         list_get_length=-1
        endif
        return
        end function list_get_length
!------------------------------------------------------------
        integer(int_kind) function list_reset(this,go_to_end)
!This function resets the current entry of a list to its beginning or its end.
        implicit none
        class(list_two_way_t):: this
        logical, intent(in), optional:: go_to_end
        logical res
        list_reset=0
        if(this%list_max_length.ge.0) then !initialized list
         if(this%list_length.gt.0) then !non-empty list
          if(present(go_to_end)) then; res=go_to_end; else; res=go_to_the_top; endif
          if(res) then !go to the end of the list
           this%current=this%last
          else !go to the beginning of the list
           this%current=this%first
          endif
         else !empty list
          list_reset=list_empty
         endif
        else
         list_reset=list_err_list_null
        endif
        return
        end function list_reset
!-------------------------------------------------------------
        integer(int_kind) function list_test_first(this,first)
!This function returns .true. in <first> iff the current list entry is positioned at the top of the list.
        implicit none
        class(list_two_way_t):: this
        logical, intent(out):: first
        list_test_first=0; first=.false.
        if(this%list_max_length.ge.0) then
         if(this%list_length.gt.0) then
          if(this%current.eq.this%first) first=.true.
         else
          list_test_first=list_empty
         endif
        else
         list_test_first=list_err_list_null
        endif
        return
        end function list_test_first
!-----------------------------------------------------------
        integer(int_kind) function list_test_last(this,last)
!This function returns .true. in <last> iff the current list entry is positioned at the end of the list.
        implicit none
        class(list_two_way_t):: this
        logical, intent(out):: last
        list_test_last=0; last=.false.
        if(this%list_max_length.ge.0) then
         if(this%list_length.gt.0) then
          if(this%current.eq.this%last) last=.true.
         else
          list_test_last=list_empty
         endif
        else
         list_test_last=list_err_list_null
        endif
        return
        end function list_test_last
!--------------------------------------------------------------------
        integer(int_kind) function list_set_position(this,curr_entry)
!This function sets the current position in the list to a given value <curr_entry>.
!The current list position cannot point to dead entries (those not in the current list).
        implicit none
        class(list_two_way_t):: this
        integer(int_kind), intent(in):: curr_entry !array entry number (must belong to this linked list)
        integer(int_kind):: i
        list_set_position=0
        if(this%list_max_length.ge.0) then !initialized list
         if(this%list_length.gt.0) then !non-empty initialized list
          if(curr_entry.ge.this%base_offset.and.curr_entry.lt.this%base_offset+this%list_max_length) then !bounds are ok
           i=(curr_entry-this%base_offset)+1
           if(this%next(i).gt.0) then !active list entry (belongs to this linked list)
            this%current=i
           else !dead entry (does not belong to this linked list)
            list_set_position=list_err_entry_dead
           endif
          else
           list_set_position=list_err_invalid_arg
          endif
         else !empty initialized list
          list_set_position=list_empty
         endif
        else
         list_set_position=list_err_list_null
        endif
        return
        end function list_set_position
!--------------------------------------------------------------------
        integer(int_kind) function list_get_position(this,curr_entry)
!This function returns the current list entry in <curr_entry>.
        implicit none
        class(list_two_way_t):: this
        integer(int_kind), intent(out):: curr_entry
        list_get_position=0
        if(this%list_max_length.ge.0) then !initialized list
         if(this%list_length.gt.0) then !non-empty initialized list
          curr_entry=this%base_offset+this%current-1
         else !empty initialized list
          curr_entry=this%base_offset-1; list_get_position=list_empty
         endif
        else
         list_get_position=list_err_list_null
        endif
        return
        end function list_get_position
!-------------------------------------------------------------------
        integer(int_kind) function list_move_to_next(this,entry_num)
!This function moves the current list entry pointer to the next linked entry.
!The number of the next entry is (optionally) returned in <entry_num>.
        implicit none
        class(list_two_way_t):: this
        integer(int_kind), intent(out), optional:: entry_num
        integer(int_kind):: i
        list_move_to_next=0
        if(this%list_max_length.ge.0) then !initialized list
         if(this%list_length.gt.0) then !non-empty initialized list
          i=this%next(this%current)
          if(i.gt.0.and.i.le.this%list_max_length) then !link points to another active list entry
           this%current=i; if(present(entry_num)) entry_num=this%base_offset+this%current-1
          elseif(i.gt.this%list_max_length) then !end of the list reached
           list_move_to_next=list_end
          else !error: this%current pointed to a dead list entry
           list_move_to_next=list_err_entry_dead
          endif
         else !empty initialized list
          list_move_to_next=list_empty
         endif
        else
         list_move_to_next=list_err_list_null
        endif
        return
        end function list_move_to_next
!-------------------------------------------------------------------
        integer(int_kind) function list_move_to_prev(this,entry_num)
!This function moves the current list entry pointer to the previous linked entry.
!The number of the previous entry is (optionally) returned in <entry_num>.
        implicit none
        class(list_two_way_t):: this
        integer(int_kind), intent(out), optional:: entry_num
        integer(int_kind):: i
        list_move_to_prev=0
        if(this%list_max_length.ge.0) then !initialized list
         if(this%list_length.gt.0) then !non-empty initialized list
          i=this%prev(this%current)
          if(i.gt.0.and.i.le.this%list_max_length) then !link points to another active list entry
           this%current=i; if(present(entry_num)) entry_num=this%base_offset+this%current-1
          elseif(i.gt.this%list_max_length) then !end of the list reached
           list_move_to_prev=list_end
          else !error: this%current points to a dead list entry
           list_move_to_prev=list_err_entry_dead
          endif
         else !empty initialized list
          list_move_to_prev=list_empty
         endif
        else
         list_move_to_prev=list_err_list_null
        endif
        return
        end function list_move_to_prev
!--------------------------------------------------------------------------------------
        integer(int_kind) function list_add_item_two_way(this,new_entry_num,add_before)
!This function adds a new entry in the linked list. The new entry is linked either
!right after the current list entry or right before, depending on the <add_before>
!logical value. The current list pointer is moved to the new entry.
        implicit none
        class(list_two_way_t):: this
        integer(int_kind), intent(in):: new_entry_num !new entry number
        logical, intent(in), optional:: add_before !where to add (.true.: before; .false.: after)
        integer(int_kind):: i,j,k
        logical res
        list_add_item_two_way=0
        if(this%list_max_length.ge.0) then !initialized list
         if(this%list_length.lt.this%list_max_length) then !list is not full
          if(new_entry_num.ge.this%base_offset.and.new_entry_num.lt.this%base_offset+this%list_max_length) then !bounds ok
           i=new_entry_num-this%base_offset+1
           if(this%next(i).lt.0.and.this%prev(i).lt.0) then !the entry was not in the linked list
            j=-this%prev(i); k=-this%next(i)
            if(k.le.this%list_max_length) then; this%prev(k)=-j; if(i.eq.this%free_ffe) this%free_ffe=k; endif
            if(j.le.this%list_max_length) this%next(j)=-k
            if(this%list_length.gt.0) then !list already had at least one element
             if(this%current.gt.0.and.this%current.le.this%list_max_length) then
              j=this%prev(this%current); k=this%next(this%current)
              if(j.gt.0.and.k.gt.0) then !active entry
               if(present(add_before)) then; res=add_before; else; res=append_after; endif
               if(res) then !add before
                if(j.le.this%list_max_length) then; this%next(j)=i; else; this%first=i; endif
                this%prev(i)=j; this%next(i)=this%current; this%prev(this%current)=i
               else !add after
                if(k.le.this%list_max_length) then; this%prev(k)=i; else; this%last=i; endif
                this%next(i)=k; this%prev(i)=this%current; this%next(this%current)=i
               endif
               this%current=i
              else
               list_add_item_two_way=list_err_list_corrupted; return
              endif
             else
              list_add_item_two_way=list_err_list_corrupted; return
             endif
            else !list was empty
             this%next(i)=this%list_max_length+1; this%prev(i)=this%list_max_length+1
             this%first=i; this%last=i; this%current=i
            endif
            this%list_length=this%list_length+1
           else !the entry is already in the list: cannot be added again
            list_add_item_two_way=list_err_entry_busy
           endif
          else
           list_add_item_two_way=list_err_invalid_arg
          endif
         else
          list_add_item_two_way=list_full
         endif
        else
         list_add_item_two_way=list_err_list_null
        endif
        return
        end function list_add_item_two_way
!--------------------------------------------------------------------------
        integer(int_kind) function list_delete_item_two_way(this,entry_num)
!This function deletes the current list entry or the entry specified by user in <entry_num>.
!In the first case, if the previous linked entry exists, it becomes the current list entry;
!otherwise, if the next linked entry exists, it becomes the current list entry;
!otherwise the list will become empty after this deletion. In the second case (<entry_num>),
!the current list entry is redefined only if it coincides with the <entry_num>.
        implicit none
        class(list_two_way_t):: this
        integer(int_kind), intent(in), optional:: entry_num
        integer(int_kind):: j,k,l
        list_delete_item_two_way=0
        if(this%list_max_length.ge.0) then !initialized list
         if(this%list_length.gt.0) then !non-empty list
          l=this%current
          if(present(entry_num)) then
           if(entry_num.ge.this%base_offset.and.entry_num.lt.this%base_offset+this%list_max_length) then
            l=entry_num-this%base_offset+1
           else
            list_delete_item_two_way=list_err_invalid_arg; return
           endif
          endif
          if(l.gt.0.and.l.le.this%list_max_length) then !bounds
           if(this%next(l).gt.0.and.this%prev(l).gt.0) then !active entry
            j=this%prev(l); k=this%next(l)
            if(j.le.this%list_max_length) then
             this%next(j)=k; if(l.eq.this%current) this%current=j
            else
             if(k.le.this%list_max_length) then; this%first=k; else; this%first=-1; this%last=-1; this%current=-1; endif
            endif
            if(k.le.this%list_max_length) then
             this%prev(k)=j; if(l.eq.this%current) this%current=k
            else
             if(j.le.this%list_max_length) then; this%last=j; else; this%first=-1; this%last=-1; this%current=-1; endif
            endif
            this%prev(l)=-(this%list_max_length+1); this%next(l)=-(this%free_ffe); this%prev(this%free_ffe)=-l
            this%free_ffe=l; this%list_length=this%list_length-1
           else !entry already dead
            list_delete_item_two_way=list_err_entry_dead
           endif
          else
           list_delete_item_two_way=list_err_list_corrupted
          endif
         else
          list_delete_item_two_way=list_empty
         endif
        else
         list_delete_item_two_way=list_err_list_null
        endif
        return
        end function list_delete_item_two_way
!--------------------------------------------------------------------
        integer(int_kind) function list_test_two_way(this)
!The function checks the correctness of a linked list.
        implicit none
        class(list_two_way_t):: this
        integer(int_kind):: i,m,n
        list_test_two_way=0
        if(this%list_max_length.ge.0) then
         if(this%list_length.gt.0) then
          if(this%first.gt.0.and.this%first.le.this%list_max_length.and. &
             this%last.gt.0.and.this%last.le.this%list_max_length.and. &
             this%current.gt.0.and.this%current.le.this%list_max_length.and. &
             this%free_ffe.gt.0.and.this%free_ffe.le.this%list_max_length+1) then
 !Test forward:
           i=this%first; m=this%list_max_length+1; n=0
           do while(i.gt.0.and.i.le.this%list_max_length)
            if(this%prev(i).ne.m) then; list_test_two_way=list_err_list_corrupted; return; endif
            m=i; i=this%next(i); n=n+1
            if(n.gt.this%list_length) then; list_test_two_way=list_err_list_corrupted; return; endif
           enddo
           if(m.ne.this%last.or.i.ne.this%list_max_length+1.or.n.ne.this%list_length) then
            list_test_two_way=list_err_list_corrupted; return
           endif
 !Test backward:
           i=this%last; m=this%list_max_length+1; n=0
           do while(i.gt.0.and.i.le.this%list_max_length)
            if(this%next(i).ne.m) then; list_test_two_way=list_err_list_corrupted; return; endif
            m=i; i=this%prev(i); n=n+1
            if(n.gt.this%list_length) then; list_test_two_way=list_err_list_corrupted; return; endif
           enddo
           if(m.ne.this%first.or.i.ne.this%list_max_length+1.or.n.ne.this%list_length) then
            list_test_two_way=list_err_list_corrupted; return
           endif
 !Test free:
           i=this%free_ffe; m=this%list_max_length+1; n=0
           do while(i.gt.0.and.i.le.this%list_max_length)
            if(-this%prev(i).ne.m) then; list_test_two_way=list_err_list_corrupted; return; endif
            m=i; i=-this%next(i); n=n+1
            if(n.gt.this%list_max_length-this%list_length) then; list_test_two_way=list_err_list_corrupted; return; endif
           enddo
           if(i.ne.this%list_max_length+1.or.n.ne.this%list_max_length-this%list_length) then
            list_test_two_way=list_err_list_corrupted; return
           endif
          else
           list_test_two_way=list_err_list_corrupted
          endif
         else
          list_test_two_way=list_empty
         endif
        else
         list_test_two_way=list_err_list_null
        endif
        return
        end function list_test_two_way

       end module lists
