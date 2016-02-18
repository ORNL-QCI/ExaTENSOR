       module lists
!Special realization of linked lists.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2016/02/18
!DESCRIPTION:
! CLASS(list_two_way_t):
!  NOTES:
!  # To avoid frequent memory allocates/deallocates, this class assumes that there exists
!    a known-size preallocated 1D array of objects that will contain the actualy data.
!    The class itself only supplies the infrastructure needed to organize the entries of that
!    1D array into a linked bidirectional list. Namely, it manipulates solely the entry numbers,
!    not the data stored in those entries. That is, an entry can be added to the linked list,
!    the next/previous linked entry number can be retrieved, the current entry can be
!    set/reset/retrieved, an entry number can be excluded from the linked list, etc.
!    Thus, this class does not deal with the actual data stored at all. It only helps to
!    organize an existing 1D array into a linked bidirectional list. Moreover, for each
!    existing 1D array multiple linked lists can be imposed on it, reflecting different
!    orderings of the array entries needed for specific purposes.
!  # Maximal length of a linked list can be zero (always empty list).
!    A negative max length means an uninitialized list.
        use, intrinsic:: ISO_C_BINDING, only: C_INT
        implicit none
        private
!PARAMETERS:
        integer, parameter, private:: INT_KIND=C_INT                       !default integer in this module
        integer(INT_KIND), parameter, public:: LIST_ERR_INVALID_ARG=1      !invalid argument passed
        integer(INT_KIND), parameter, public:: LIST_ERR_MEM_ALLOC_FAILED=2 !memory allocation failed
        integer(INT_KIND), parameter, public:: LIST_ERR_MEM_FREE_FAILED=3  !memory deallocation failed
        integer(INT_KIND), parameter, public:: LIST_ERR_LIST_NULL=4        !list has not been initialized yet
        integer(INT_KIND), parameter, public:: LIST_ERR_LIST_EXISTS=5      !list has already been initialized
        integer(INT_KIND), parameter, public:: LIST_ERR_LIST_CORRUPTED=6   !list linking corrupted
        integer(INT_KIND), parameter, public:: LIST_ERR_ENTRY_DEAD=7       !attempt to use a free (dead) entry (the one not in the list)
        integer(INT_KIND), parameter, public:: LIST_ERR_ENTRY_BUSY=8       !attempt to add an entry which is already in the list
        integer(INT_KIND), parameter, public:: LIST_EMPTY=9                !initialized list is empty
        integer(INT_KIND), parameter, public:: LIST_FULL=10                !initialized list is full
        integer(INT_KIND), parameter, public:: LIST_END=11                 !end of an initialized list (in any direction)
        logical, parameter, public:: GO_TO_THE_TOP=.false.
        logical, parameter, public:: GO_TO_THE_END=.true.
        logical, parameter, public:: APPEND_BEFORE=.true.
        logical, parameter, public:: APPEND_AFTER=.false.
!TYPES:
 !Bidirectional linked list (linked in both ways):
        type, public:: list_two_way_t
         integer(INT_KIND), private:: list_max_length=-1   !max length of the list: range=[base_offset:base_offset+list_max_length-1]
         integer(INT_KIND), private:: base_offset          !lowest bound in the numeration of entries of the served 1D array
         integer(INT_KIND), private:: list_length          !current length of the linked list
         integer(INT_KIND), private:: first                !first entry number in the list (can be any from the above range)
         integer(INT_KIND), private:: last                 !last entry number in the list (can be any from the above range)
         integer(INT_KIND), private:: current              !current entry number in the list
         integer(INT_KIND), private:: free_ffe             !first free entry
         integer(INT_KIND), allocatable, private:: next(:) !next entry (negative: dead entry; positive: active entry)
         integer(INT_KIND), allocatable, private:: prev(:) !previous entry (negative: dead entry; positive: active entry)
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
        integer(INT_KIND) function list_create_two_way(this,list_max_len,list_beg)
!This function creates a bidirectional linked list infrastructure for some external 1D array.
!INPUT:
! # list_max_len - max length of the list (the length of the external 1D array served);
! # list_beg - lowest bound of the external 1D array: range of entry numbers = [list_beg:list_beg+list_max_len-1];
        implicit none
        class(list_two_way_t):: this
        integer(INT_KIND), intent(in):: list_max_len,list_beg
        integer(INT_KIND):: i
        list_create_two_way=0
        if(list_max_len.ge.0.and.list_max_len.lt.huge(i).and.huge(i)-list_max_len.ge.list_beg) then !check arguments
         if(this%list_max_length.lt.0) then !list is uninitialized
          this%base_offset=list_beg; this%list_length=0; this%first=-1; this%last=-1; this%current=-1
          if(list_max_len.gt.0) then
           allocate(this%next(1:list_max_len),STAT=i)
           if(i.ne.0) then; list_create_two_way=LIST_ERR_MEM_ALLOC_FAILED; return; endif
           do i=1,list_max_len; this%next(i)=-(i+1); enddo; this%free_ffe=1
           allocate(this%prev(1:list_max_len),STAT=i)
           if(i.ne.0) then; deallocate(this%next); list_create_two_way=LIST_ERR_MEM_ALLOC_FAILED; return; endif
           do i=2,list_max_len; this%prev(i)=-(i-1); enddo; this%prev(1)=-(list_max_len+1)
           this%list_max_length=list_max_len
          else !initialize a trivial zero-max-length list (always empty)
           this%list_max_length=0; this%free_ffe=-1
          endif
         else !list has already been initialized
          list_create_two_way=LIST_ERR_LIST_EXISTS
         endif
        else
         list_create_two_way=LIST_ERR_INVALID_ARG
        endif
        return
        end function list_create_two_way
!------------------------------------------------------------
        integer(INT_KIND) function list_destroy_two_way(this)
!This function destroys a bidirectional linked list infrastructure.
        implicit none
        class(list_two_way_t):: this
        integer(INT_KIND):: i
        list_destroy_two_way=0
        if(this%list_max_length.ge.0) then !initialized list
         this%list_max_length=-1; this%base_offset=-1; this%list_length=-1
         this%first=-1; this%last=-1; this%current=-1; this%free_ffe=-1
         if(allocated(this%next)) then
          deallocate(this%next,STAT=i); if(i.ne.0) list_destroy_two_way=LIST_ERR_MEM_FREE_FAILED
         endif
         if(allocated(this%prev)) then
          deallocate(this%prev,STAT=i); if(i.ne.0) list_destroy_two_way=LIST_ERR_MEM_FREE_FAILED
         endif
        else !uninitialized list
         list_destroy_two_way=LIST_ERR_LIST_NULL
        endif
        return
        end function list_destroy_two_way
!-----------------------------------------------------------
        integer(INT_KIND) function list_get_max_length(this)
!This function returns the maximal length of a list (negative will mean an uninitialized list).
        implicit none
        class(list_two_way_t):: this
        list_get_max_length=this%list_max_length
        return
        end function list_get_max_length
!-------------------------------------------------------
        integer(INT_KIND) function list_get_length(this)
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
        integer(INT_KIND) function list_reset(this,go_to_end)
!This function resets the current entry of a list to its beginning or its end.
        implicit none
        class(list_two_way_t):: this
        logical, intent(in), optional:: go_to_end
        logical res
        list_reset=0
        if(this%list_max_length.ge.0) then !initialized list
         if(this%list_length.gt.0) then !non-empty list
          if(present(go_to_end)) then; res=go_to_end; else; res=GO_TO_THE_TOP; endif
          if(res) then !go to the end of the list
           this%current=this%last
          else !go to the beginning of the list
           this%current=this%first
          endif
         else !empty list
          list_reset=LIST_EMPTY
         endif
        else
         list_reset=LIST_ERR_LIST_NULL
        endif
        return
        end function list_reset
!-------------------------------------------------------------
        integer(INT_KIND) function list_test_first(this,first)
!This function returns .true. in <first> iff the current list entry is positioned at the top of the list.
        implicit none
        class(list_two_way_t):: this
        logical, intent(out):: first
        list_test_first=0; first=.false.
        if(this%list_max_length.ge.0) then
         if(this%list_length.gt.0) then
          if(this%current.eq.this%first) first=.true.
         else
          list_test_first=LIST_EMPTY
         endif
        else
         list_test_first=LIST_ERR_LIST_NULL
        endif
        return
        end function list_test_first
!-----------------------------------------------------------
        integer(INT_KIND) function list_test_last(this,last)
!This function returns .true. in <last> iff the current list entry is positioned at the end of the list.
        implicit none
        class(list_two_way_t):: this
        logical, intent(out):: last
        list_test_last=0; last=.false.
        if(this%list_max_length.ge.0) then
         if(this%list_length.gt.0) then
          if(this%current.eq.this%last) last=.true.
         else
          list_test_last=LIST_EMPTY
         endif
        else
         list_test_last=LIST_ERR_LIST_NULL
        endif
        return
        end function list_test_last
!--------------------------------------------------------------------
        integer(INT_KIND) function list_set_position(this,curr_entry)
!This function sets the current position in the list to a given value <curr_entry>.
!The current list position cannot point to dead entries (those not in the current list).
        implicit none
        class(list_two_way_t):: this
        integer(INT_KIND), intent(in):: curr_entry !array entry number (must belong to this linked list)
        integer(INT_KIND):: i
        list_set_position=0
        if(this%list_max_length.ge.0) then !initialized list
         if(this%list_length.gt.0) then !non-empty initialized list
          if(curr_entry.ge.this%base_offset.and.curr_entry.lt.this%base_offset+this%list_max_length) then !bounds are ok
           i=(curr_entry-this%base_offset)+1
           if(this%next(i).gt.0) then !active list entry (belongs to this linked list)
            this%current=i
           else !dead entry (does not belong to this linked list)
            list_set_position=LIST_ERR_ENTRY_DEAD
           endif
          else
           list_set_position=LIST_ERR_INVALID_ARG
          endif
         else !empty initialized list
          list_set_position=LIST_EMPTY
         endif
        else
         list_set_position=LIST_ERR_LIST_NULL
        endif
        return
        end function list_set_position
!--------------------------------------------------------------------
        integer(INT_KIND) function list_get_position(this,curr_entry)
!This function returns the current list entry in <curr_entry>.
        implicit none
        class(list_two_way_t):: this
        integer(INT_KIND), intent(out):: curr_entry
        list_get_position=0
        if(this%list_max_length.ge.0) then !initialized list
         if(this%list_length.gt.0) then !non-empty initialized list
          curr_entry=this%base_offset+this%current-1
         else !empty initialized list
          curr_entry=this%base_offset-1; list_get_position=LIST_EMPTY
         endif
        else
         list_get_position=LIST_ERR_LIST_NULL
        endif
        return
        end function list_get_position
!-------------------------------------------------------------------
        integer(INT_KIND) function list_move_to_next(this,entry_num)
!This function moves the current list entry pointer to the next linked entry.
!The number of the next entry is (optionally) returned in <entry_num>.
        implicit none
        class(list_two_way_t):: this
        integer(INT_KIND), intent(out), optional:: entry_num
        integer(INT_KIND):: i
        list_move_to_next=0
        if(this%list_max_length.ge.0) then !initialized list
         if(this%list_length.gt.0) then !non-empty initialized list
          i=this%next(this%current)
          if(i.gt.0.and.i.le.this%list_max_length) then !link points to another active list entry
           this%current=i; if(present(entry_num)) entry_num=this%base_offset+this%current-1
          elseif(i.gt.this%list_max_length) then !end of the list reached
           list_move_to_next=LIST_END
          else !error: this%current pointed to a dead list entry
           list_move_to_next=LIST_ERR_ENTRY_DEAD
          endif
         else !empty initialized list
          list_move_to_next=LIST_EMPTY
         endif
        else
         list_move_to_next=LIST_ERR_LIST_NULL
        endif
        return
        end function list_move_to_next
!-------------------------------------------------------------------
        integer(INT_KIND) function list_move_to_prev(this,entry_num)
!This function moves the current list entry pointer to the previous linked entry.
!The number of the previous entry is (optionally) returned in <entry_num>.
        implicit none
        class(list_two_way_t):: this
        integer(INT_KIND), intent(out), optional:: entry_num
        integer(INT_KIND):: i
        list_move_to_prev=0
        if(this%list_max_length.ge.0) then !initialized list
         if(this%list_length.gt.0) then !non-empty initialized list
          i=this%prev(this%current)
          if(i.gt.0.and.i.le.this%list_max_length) then !link points to another active list entry
           this%current=i; if(present(entry_num)) entry_num=this%base_offset+this%current-1
          elseif(i.gt.this%list_max_length) then !end of the list reached
           list_move_to_prev=LIST_END
          else !error: this%current points to a dead list entry
           list_move_to_prev=LIST_ERR_ENTRY_DEAD
          endif
         else !empty initialized list
          list_move_to_prev=LIST_EMPTY
         endif
        else
         list_move_to_prev=LIST_ERR_LIST_NULL
        endif
        return
        end function list_move_to_prev
!--------------------------------------------------------------------------------------
        integer(INT_KIND) function list_add_item_two_way(this,new_entry_num,add_before)
!This function adds a new entry in the linked list. The new entry is linked either
!right after the current list entry or right before, depending on the <add_before>
!logical value. The current list pointer is moved to the new entry. Rules:
!If <new_entry_num> is within the valid range, it is considered as an explicit input,
!that is, an entry number to be added to the current linked list. If <new_entry_num>
!is outside the valid range, than the first free entry will be the entry to add.
!Its number will be returned in <new_entry_num>.
        implicit none
        class(list_two_way_t):: this
        integer(INT_KIND), intent(inout):: new_entry_num !new entry number (input or output)
        logical, intent(in), optional:: add_before !where to add (.true.: before; .false.: after)
        integer(INT_KIND):: i,j,k
        logical res
        list_add_item_two_way=0
        if(this%list_max_length.ge.0) then !initialized list
         if(this%list_length.lt.this%list_max_length) then !list is not full
          if(new_entry_num.ge.this%base_offset.and.new_entry_num.lt.this%base_offset+this%list_max_length) then
           i=new_entry_num-this%base_offset+1 !an explicit free entry number was specified by user: use it
          else
           i=this%free_ffe !user did not specify an explicit entry number: use the first free entry
          endif
          if(i.gt.0.and.i.le.this%list_max_length) then !bounds ok
           if(this%next(i).lt.0.and.this%prev(i).lt.0) then !the entry was not in the linked list
            j=-this%prev(i); k=-this%next(i)
            if(k.le.this%list_max_length) then; this%prev(k)=-j; if(i.eq.this%free_ffe) this%free_ffe=k; endif
            if(j.le.this%list_max_length) this%next(j)=-k
            if(this%list_length.gt.0) then !list already had at least one element
             if(this%current.gt.0.and.this%current.le.this%list_max_length) then
              j=this%prev(this%current); k=this%next(this%current)
              if(j.gt.0.and.k.gt.0) then !active entry
               if(present(add_before)) then; res=add_before; else; res=APPEND_AFTER; endif
               if(res) then !add before
                if(j.le.this%list_max_length) then; this%next(j)=i; else; this%first=i; endif
                this%prev(i)=j; this%next(i)=this%current; this%prev(this%current)=i
               else !add after
                if(k.le.this%list_max_length) then; this%prev(k)=i; else; this%last=i; endif
                this%next(i)=k; this%prev(i)=this%current; this%next(this%current)=i
               endif
               this%current=i
              else
               list_add_item_two_way=LIST_ERR_LIST_CORRUPTED; return
              endif
             else
              list_add_item_two_way=LIST_ERR_LIST_CORRUPTED; return
             endif
            else !list was empty
             this%next(i)=this%list_max_length+1; this%prev(i)=this%list_max_length+1
             this%first=i; this%last=i; this%current=i
            endif
            this%list_length=this%list_length+1
           else !the entry is already in the list: cannot be added again
            list_add_item_two_way=LIST_ERR_ENTRY_BUSY
           endif
          else
           list_add_item_two_way=LIST_ERR_LIST_CORRUPTED
          endif
         else
          list_add_item_two_way=LIST_FULL
         endif
        else
         list_add_item_two_way=LIST_ERR_LIST_NULL
        endif
        return
        end function list_add_item_two_way
!--------------------------------------------------------------------------
        integer(INT_KIND) function list_delete_item_two_way(this,entry_num)
!This function deletes the current list entry or the entry specified by user in <entry_num>.
!In the first case, if the previous linked entry exists, it becomes the current list entry;
!otherwise, if the next linked entry exists, it becomes the current list entry;
!otherwise the list will become empty after this deletion. In the second case (<entry_num>),
!the current list entry is redefined only if it coincides with the <entry_num>.
        implicit none
        class(list_two_way_t):: this
        integer(INT_KIND), intent(in), optional:: entry_num
        integer(INT_KIND):: j,k,l
        list_delete_item_two_way=0
        if(this%list_max_length.ge.0) then !initialized list
         if(this%list_length.gt.0) then !non-empty list
          l=this%current
          if(present(entry_num)) then
           if(entry_num.ge.this%base_offset.and.entry_num.lt.this%base_offset+this%list_max_length) then
            l=entry_num-this%base_offset+1
           else
            list_delete_item_two_way=LIST_ERR_INVALID_ARG; return
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
            list_delete_item_two_way=LIST_ERR_ENTRY_DEAD
           endif
          else
           list_delete_item_two_way=LIST_ERR_LIST_CORRUPTED
          endif
         else
          list_delete_item_two_way=LIST_EMPTY
         endif
        else
         list_delete_item_two_way=LIST_ERR_LIST_NULL
        endif
        return
        end function list_delete_item_two_way
!--------------------------------------------------------------------
        integer(INT_KIND) function list_test_two_way(this)
!The function checks the correctness of a linked list.
        implicit none
        class(list_two_way_t):: this
        integer(INT_KIND):: i,m,n
        list_test_two_way=0
        if(this%list_max_length.ge.0) then
         if(this%list_length.gt.0) then
          if(this%first.gt.0.and.this%first.le.this%list_max_length.and.&
             &this%last.gt.0.and.this%last.le.this%list_max_length.and.&
             &this%current.gt.0.and.this%current.le.this%list_max_length.and.&
             &this%free_ffe.gt.0.and.this%free_ffe.le.this%list_max_length+1) then
 !Test forward:
           i=this%first; m=this%list_max_length+1; n=0
           do while(i.gt.0.and.i.le.this%list_max_length)
            if(this%prev(i).ne.m) then; list_test_two_way=LIST_ERR_LIST_CORRUPTED; return; endif
            m=i; i=this%next(i); n=n+1
            if(n.gt.this%list_length) then; list_test_two_way=LIST_ERR_LIST_CORRUPTED; return; endif
           enddo
           if(m.ne.this%last.or.i.ne.this%list_max_length+1.or.n.ne.this%list_length) then
            list_test_two_way=LIST_ERR_LIST_CORRUPTED; return
           endif
 !Test backward:
           i=this%last; m=this%list_max_length+1; n=0
           do while(i.gt.0.and.i.le.this%list_max_length)
            if(this%next(i).ne.m) then; list_test_two_way=LIST_ERR_LIST_CORRUPTED; return; endif
            m=i; i=this%prev(i); n=n+1
            if(n.gt.this%list_length) then; list_test_two_way=LIST_ERR_LIST_CORRUPTED; return; endif
           enddo
           if(m.ne.this%first.or.i.ne.this%list_max_length+1.or.n.ne.this%list_length) then
            list_test_two_way=LIST_ERR_LIST_CORRUPTED; return
           endif
 !Test free:
           i=this%free_ffe; m=this%list_max_length+1; n=0
           do while(i.gt.0.and.i.le.this%list_max_length)
            if(-this%prev(i).ne.m) then; list_test_two_way=LIST_ERR_LIST_CORRUPTED; return; endif
            m=i; i=-this%next(i); n=n+1
            if(n.gt.this%list_max_length-this%list_length) then; list_test_two_way=LIST_ERR_LIST_CORRUPTED; return; endif
           enddo
           if(i.ne.this%list_max_length+1.or.n.ne.this%list_max_length-this%list_length) then
            list_test_two_way=LIST_ERR_LIST_CORRUPTED; return
           endif
          else
           list_test_two_way=LIST_ERR_LIST_CORRUPTED
          endif
         else
          list_test_two_way=LIST_EMPTY
         endif
        else
         list_test_two_way=LIST_ERR_LIST_NULL
        endif
        return
        end function list_test_two_way

       end module lists
