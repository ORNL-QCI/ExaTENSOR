       module lists
!Realization of linked lists.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2014/07/15
!NOTES:
! # CLASS(list_two_way_t):
!   To avoid frequent memory allocates/deallocates, this class assumes that
!   there exists a known-size array of objects which will contain the actualy data.
!   The class itself only supplies the infrastructure needed to organize that array of data
!   into a linked bidirectional list. Namely, it manipulates solely the entry numbers,
!   not the data stored in those entries. That is, an empty entry number can be provided,
!   the next/previous linked entry number can be returned, the current entry number can be
!   set/reset/retrieved, an entry number can be returned back to the stack of free entries.
!   Thus, this class does not deal with the actual data stored at all. It only helps to
!   organize an existing array into a linked bidirectional list.
! # Maximal length of a list can be zero (always empty list).
!   A negative max length means an uninitialized list.
        use, intrinsic:: ISO_C_BINDING
!PARAMETERS:
        integer, parameter, private:: int_kind=C_INT                       !default integer in this module
        integer(int_kind), parameter, public:: list_err_invalid_arg=1      !invalid argument passed
        integer(int_kind), parameter, public:: list_err_mem_alloc_failed=2 !memory allocation failed
        integer(int_kind), parameter, public:: list_err_mem_free_failed=3  !memory deallocation failed
        integer(int_kind), parameter, public:: list_err_list_null=4        !list has not been initialized yet
        integer(int_kind), parameter, public:: list_err_list_exists=5      !list has already been initialized
        integer(int_kind), parameter, public:: list_err_list_corrupted=6   !list linking corrupted
        integer(int_kind), parameter, public:: list_err_list_entry_null=7  !attempt to address a null (dead) entry
        integer(int_kind), parameter, public:: list_empty=7                !initialized list is empty
        integer(int_kind), parameter, public:: list_full=8                 !initialized list is empty
        integer(int_kind), parameter, public:: list_end=9                  !end of an initialized list (in any direction)
        integer(int_kind), parameter, private:: list_entry_null=-1         !marks null (dead) entries
        integer(int_kind), parameter, private:: list_link_null=-2          !marks null links on live entries
        logical, parameter, public:: go_to_the_top=.false.
        logical, parameter, public:: go_to_the_end=.true.
        logical, parameter, public:: append_before=.true.
        logical, parameter, public:: append_after=.false.
!TYPES:
 !Bidirectional linked list (linked in both ways):
        type, public:: list_two_way_t
         integer(int_kind), private:: list_max_length=-1          !max length of the list: range=[base_offset:base_offset+max_length-1]
         integer(int_kind), private:: base_offset                 !lowest bound in the numeration of entries: range=[base_offset:base_offset+max_length-1]
         integer(int_kind), private:: first                       !first entry number (can be any from the range)
         integer(int_kind), private:: last                        !last entry number (can be any from the range)
         integer(int_kind), private:: current                     !current entry number
         integer(int_kind), private:: free_p                      !free stack ponter = current length of the list
         integer(int_kind), allocatable, private:: free_stack(:)  !free entries stack
         integer(int_kind), allocatable, private:: next(:)        !next element
         integer(int_kind), allocatable, private:: prev(:)        !previous element
         contains
          procedure, public:: create=>list_create_two_way      !create a bidirectional linked list
          procedure, public:: destroy=>list_destroy_two_way    !destroy a bidirectional linked list
          procedure, public:: max_length=>list_get_max_length  !get max length of the list
          procedure, public:: length=>list_get_length          !get current length of the list
          procedure, public:: reset=>list_reset                !reset the current position of the list
          procedure, public:: test_first=>list_test_first      !test whether current == first
          procedure, public:: test_last=>list_test_last        !test whether current == last
          procedure, public:: set_position=>list_set_position  !set the current list position to a given entry
          procedure, public:: get_position=>list_get_position  !get the current position in the list
          procedure, public:: go_next=>list_move_to_next       !move to the next entry in the linked list
          procedure, public:: go_previous=>list_move_to_prev   !move to the previous entry in the linked list
          procedure, public:: add=>list_add_item_two_way       !register a new entry in the linked list
          procedure, public:: delete=>list_delete_item_two_way !delete an entry in the linked list
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

       contains
!---------------------------------------------------------------------------------
        integer(int_kind) function list_create_two_way(this,list_max_len,list_beg)
!This function creates a bidirectional linked list infrastructure for some external array.
!INPUT:
! # list_max_len - max length of the list (the length of the external array served);
! # list_beg - lowest bound of the external array: range of entry numbers = [list_beg:list_beg+list_max_len-1];
        implicit none
        class(list_two_way_t):: this
        integer(int_kind), intent(in):: list_max_len,list_beg
        integer(int_kind):: i
        list_create_two_way=0
        if(list_max_len.ge.0.and.(huge(i)-list_beg.ge.list_max_len-1)) then
         if(this%list_max_length.lt.0) then
          this%first=list_beg-1; this%last=list_beg-1; this%current=list_beg-1; this%free_p=-1; this%base_offset=list_beg
          if(list_max_len.gt.0) then
           allocate(this%free_stack(0:list_max_len-1),STAT=i)
           if(i.ne.0) then; list_create_two_way=list_err_mem_alloc_failed; return; endif
           do i=0,list_max_len-1; this%free_stack(i)=list_beg+i; enddo; this%free_p=0 !stack of free entry numbers
           allocate(this%next(0:list_max_len-1),STAT=i)
           if(i.ne.0) then; deallocate(this%free_stack); list_create_two_way=list_err_mem_alloc_failed; return; endif
           do i=0,list_max_len-1; this%next(i)=list_entry_null; enddo
           allocate(this%prev(0:list_max_len-1),STAT=i)
           if(i.ne.0) then; deallocate(this%free_stack,this%next); list_create_two_way=list_err_mem_alloc_failed; return; endif
           do i=0,list_max_len-1; this%prev(i)=list_entry_null; enddo
           this%list_max_length=list_max_len
          else !initialize a zero-max-length list (always empty)
           this%list_max_length=0
          endif
         else
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
         this%list_max_length=-1; this%base_offset=-1; this%first=-1; this%last=-1; this%current=-1; this%free_p=-1
         if(allocated(this%free_stack)) then 
          deallocate(this%free_stack,STAT=i); if(i.ne.0) list_destroy_two_way=list_err_mem_free_failed
         endif
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
        if(this%list_max_length.gt.0) then
         list_get_length=this%free_p
        elseif(this%list_max_length.eq.0) then
         list_get_length=0
        else !uninitialized list
         list_get_length=-1
        endif
        return
        end function list_get_length
!------------------------------------------------------------
        integer(int_kind) function list_reset(this,go_to_end)
!This function resets the current entry of a list to its beginning or end.
        implicit none
        class(list_two_way_t):: this
        logical, intent(in), optional:: go_to_end
        logical res
        list_reset=0
        if(this%list_max_length.ge.0) then !initialized list
         if(this%free_p.gt.0) then !non-empty list
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
         if(this%free_p.gt.0) then
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
         if(this%free_p.gt.0) then
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
!The current list position can point only to active list entries.
        implicit none
        class(list_two_way_t):: this
        integer(int_kind), intent(in):: curr_entry
        list_set_position=0
        if(this%list_max_length.ge.0) then !initialized list
         if(this%free_p.gt.0) then !non-empty initialized list
          if(curr_entry.ge.this%base_offset.and.curr_entry.lt.this%base_offset+this%list_max_length) then !bounds are ok
           if(this%next(curr_entry-this%base_offset).ne.list_entry_null) then !active list entry
            this%current=curr_entry
           else !dead entry
            list_set_position=list_err_list_entry_null
           endif
          else
           list_set_position=list_err_invalid_arg
          endif
         else
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
         if(this%free_p.gt.0) then !non-empty initialized list
          curr_entry=this%current
         else
          curr_entry=this%base_offset-1; list_get_position=list_empty
         endif
        else
         curr_entry=-1; list_get_position=list_err_list_null
        endif
        return
        end function list_get_position
!-------------------------------------------------------------------
        integer(int_kind) function list_move_to_next(this,entry_num)
!This function moves the current list pointer to the next linked entry.
!The number of the next entry is (optionally) returned in <entry_num>.
        implicit none
        class(list_two_way_t):: this
        integer(int_kind), intent(out), optional:: entry_num
        integer(int_kind):: i
        list_move_to_next=0
        if(this%list_max_length.ge.0) then !initialized list
         if(this%free_p.gt.0) then !non-empty initialized list
          i=this%next(this%current-this%base_offset)
          if(i.ge.0) then !link points to an active list entry
           this%current=this%base_offset+i
           if(present(entry_num)) entry_num=this%current
          elseif(i.eq.list_link_null) then !link does not point anywhere
           list_move_to_next=list_end
          else !the current entry was dead
           list_move_to_next=list_err_list_entry_null
          endif
         else
          list_move_to_next=list_empty
         endif
        else
         list_move_to_next=list_err_list_null
        endif
        return
        end function list_move_to_next
!-------------------------------------------------------------------
        integer(int_kind) function list_move_to_prev(this,entry_num)
!This function moves the current list pointer to the previous linked entry.
!The number of the previous entry is (optionally) returned in <entry_num>.
        implicit none
        class(list_two_way_t):: this
        integer(int_kind), intent(out), optional:: entry_num
        integer(int_kind):: i
        list_move_to_prev=0
        if(this%list_max_length.ge.0) then !initialized list
         if(this%free_p.gt.0) then !non-empty initialized list
          i=this%prev(this%current-this%base_offset)
          if(i.ge.0) then !link points to an active list entry
           this%current=this%base_offset+i
           if(present(entry_num)) entry_num=this%current
          elseif(i.eq.list_link_null) then !link does not point anywhere
           list_move_to_prev=list_end
          else !the current entry was dead
           list_move_to_prev=list_err_list_entry_null
          endif
         else
          list_move_to_prev=list_empty
         endif
        else
         list_move_to_prev=list_err_list_null
        endif
        return
        end function list_move_to_prev
!--------------------------------------------------------------------------------------
        integer(int_kind) function list_add_item_two_way(this,new_entry_num,add_before)
!This function registers a new entry in the linked list and returns its number in <new_entry_num>.
!The new entry is linked either right after the current list entry or right before.
!The current list pointer is set to the new entry.
        implicit none
        class(list_two_way_t):: this
        integer(int_kind), intent(out):: new_entry_num
        logical, intent(in), optional:: add_before
        integer(int_kind):: i,j
        logical res
        list_add_item_two_way=0
        if(this%list_max_length.ge.0) then !initialized list
         if(this%free_p.lt.this%list_max_length) then !list is not full
          new_entry_num=this%free_stack(this%free_p); this%free_p=this%free_p+1
          if(this%free_p.gt.1) then !add to a non-empty list: %current existed
           i=this%base_offset
           if(this%next(this%current-i).ne.list_entry_null) then !%current points to an active list entry
            if(present(add_before)) then; res=add_before; else; res=append_after; endif
            if(res) then !append before
             j=this%prev(this%current-i)
             this%prev(new_entry_num-i)=j; this%prev(this%current-i)=new_entry_num-i
             this%next(new_entry_num-i)=this%current-i; if(j.ge.0) this%next(j)=new_entry_num-i
             if(this%current.eq.this%first) this%first=new_entry_num
            else !append after
             j=this%next(this%current-i)
             this%next(new_entry_num-i)=j; this%next(this%current-i)=new_entry_num-i
             this%prev(new_entry_num-i)=this%current-i; if(j.ge.0) this%prev(j)=new_entry_num-i
             if(this%current.eq.this%last) this%last=new_entry_num
            endif
            this%current=new_entry_num
           else !%current pointed to a dead entry
            list_add_item_two_way=list_err_list_entry_null
           endif
          else !very first element is added
           this%first=new_entry_num; this%last=new_entry_num; this%current=new_entry_num
           this%next(new_entry_num-this%base_offset)=list_link_null; this%prev(new_entry_num-this%base_offset)=list_link_null
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
!This function deletes the current list entry (or the entry specified by user in <entry_num>).
!In the first case, if the previous linked entry exists, it becomes the current list entry;
!otherwise, if the next linked entry exists, it becomes the current list entry;
!otherwise the list will become empty after this deletion. In the second case (<entry_num>),
!the current list entry is redefined only if it coincides with the <entry_num>.
        implicit none
        class(list_two_way_t):: this
        integer(int_kind), intent(in), optional:: entry_num
        integer(int_kind):: i,j,b,pos
        list_delete_item_two_way=0
        if(this%list_max_length.ge.0) then !initialized list
         pos=this%current
         if(present(entry_num)) then
          if(entry_num.ge.this%base_offset.and.entry_num.lt.this%base_offset+this%list_max_length) then
           pos=entry_num
          else
           list_delete_item_two_way=list_err_invalid_arg; return
          endif
         endif
         if(this%free_p.gt.0) then !non-empty list
          b=this%base_offset; i=this%next(pos-b)
          if(i.ne.list_entry_null) then
           j=this%prev(pos-b); if(j.ge.0) this%next(j)=i; if(i.ge.0) this%prev(i)=j
           this%next(pos-b)=list_entry_null; this%prev(pos-b)=list_entry_null
           this%free_stack(this%free_p)=pos; this%free_p=this%free_p-1
           if(pos.eq.this%first) then; if(i.ge.0) then; this%first=b+i; else; this%first=b-1; endif; endif
           if(pos.eq.this%last) then; if(j.ge.0) then; this%last=b+j; else; this%last=b-1; endif; endif
           if(pos.eq.this%current) then
            if(j.ge.0) then
             this%current=b+j
            else
             if(i.ge.0) then
              this%current=b+i
             else
              if(this%free_p.eq.0) then !empty list left
               this%current=b-1
              else
               list_delete_item_two_way=list_err_list_corrupted
              endif
             endif
            endif
           endif
          else
           list_delete_item_two_way=list_err_list_entry_null
          endif
         else
          list_delete_item_two_way=list_empty
         endif
        else
         list_delete_item_two_way=list_err_list_null
        endif
        return
        end function list_delete_item_two_way

       end module lists
