       module lists
!Realization of linked lists.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2014/07/13
        use, intrinsic:: ISO_C_BINDING
!PARAMETERS:
        integer, parameter, private:: int_kind=C_INT
        integer(int_kind), parameter, public:: list_err_invalid_arg=1
        integer(int_kind), parameter, public:: list_err_mem_alloc_failed=2
        integer(int_kind), parameter, public:: list_err_mem_free_failed=3
        integer(int_kind), parameter, public:: list_err_list_null=4
        integer(int_kind), parameter, public:: list_err_list_exists=5
        integer(int_kind), parameter, public:: list_end=6
!TYPES:
        type, private:: list_t
         integer(int_kind), private:: list_max_length=-1
         integer(int_kind), private:: first
         integer(int_kind), private:: last
         integer(int_kind), private:: current
         integer(int_kind), private:: free_p
         integer(int_kind), allocatable, private:: free_stack(:)
         contains
          procedure, public, overridable:: create=>list_create
          procedure, public, overridable:: destroy=>list_destroy
          procedure, public, nonoverridable:: max_length=>list_get_max_length
          procedure, public, nonoverridable:: length=>list_get_length
          procedure, public, nonoverridable:: reset=>list_reset
        end type list_t
 !Unidirectional list (forward):
        type, extends(list_t), public:: list_one_way_t
         integer(int_kind), allocatable, private:: next(:)
         contains
          procedure, public, overridable:: create=>list_create_one_way
          procedure, public, overridable:: destroy=>list_destroy_one_way
          procedure, public, overridable:: add=>list_add_item_one_way
          procedure, public, overridable:: delete=>list_delete_item_one_way
          procedure, public, nonoverridable:: next=>list_move_to_next
        end type list_one_way_t
 !Bidirectional list (both ways):
        type, extends(list_one_way_t), public:: list_two_way_t
         integer(int_kind), allocatable, private:: prev(:)
         contains
          procedure, public, overridable:: create=>list_create_two_way
          procedure, public, overridable:: destroy=>list_destroy_two_way
          procedure, public, overridable:: add=>list_add_item_two_way
          procedure, public, overridable:: delete=>list_delete_item_two_way
          procedure, public, nonoverridable:: previous=>list_move_to_prev
        end type list_two_way_t
!FUNCTIONS:
        private list_create
        private list_destroy
        private list_create_one_way
        private list_destroy_one_way
        private list_create_two_way
        private list_destroy_two_way
        private list_get_max_length
        private list_get_length
        private list_reset
        private list_add_item_one_way
        private list_delete_item_one_way
        private list_add_item_two_way
        private list_delete_item_two_way
        private list_move_to_next
        private list_move_to_prev

       contains
!-------------------------------------------------------------------------
        integer(int_kind) function list_create(this,list_max_len,list_beg)
!This function creates a list.
!INPUT:
! # list_max_len - max length of the list;
! # list _beg - first entry number: range of entry numbers = [list_beg..(list_beg+list_max_len-1)];
        implicit none
        class(list_t):: this
        integer(int_kind), intent(in):: list_max_len,list_beg
        integer(int_kind):: i
        list_create=0
        if(list_max_len.ge.0) then
         if(this%list_max_length.lt.0) then
          this%list_max_length=0; this%first=-1; this%last=-1; this%current=-1; this%free_p=-1
          if(list_max_len.gt.0) then
           allocate(this%free_stack(0:list_max_len-1),STAT=i)
           if(i.ne.0) then; list_create=list_err_mem_alloc_failed; return; endif
           do i=0,list_max_len-1; free_stack(i)=list_beg+i; enddo !stack of free entry numbers
           this%free_p=0; this%list_max_length=list_max_len
          endif
         else
          list_create=list_err_list_exists
         endif
        else
         list_create=list_err_invalid_arg
        endif
        return
        end function list_create
!---------------------------------------------------------------------------------
        integer(int_kind) function list_create_one_way(this,list_max_len,list_beg)
!This function creates a unidirectional linked list.
!INPUT:
! # list_max_len - max length of the list;
! # list _beg - first entry number: range of entry numbers = [list_beg..(list_beg+list_max_len-1)];
        implicit none
        class(list_one_way_t):: this
        integer(int_kind), intent(in):: list_max_len,list_beg
        integer(int_kind):: i
        list_create_one_way=0
        if(list_max_len.ge.0) then
         if(this%list_max_length.lt.0) then
          this%list_max_length=0; this%first=-1; this%last=-1; this%current=-1; this%free_p=-1
          if(list_max_len.gt.0) then
           allocate(this%free_stack(0:list_max_len-1),STAT=i)
           if(i.ne.0) then; list_create_one_way=list_err_mem_alloc_failed; return; endif
           do i=0,list_max_len-1; free_stack(i)=list_beg+i; enddo; this%free_p=0 !stack of free entry numbers
           allocate(this%next(list_beg:list_beg+list_max_len-1),STAT=i)
           if(i.ne.0) then; deallocate(this%free_stack); list_create_one_way=list_err_mem_alloc_failed; return; endif
           this%list_max_length=list_max_len
          endif
         else
          list_create_one_way=list_err_list_exists
         endif
        else
         list_create_one_way=list_err_invalid_arg
        endif
        return
        end function list_create_one_way
!---------------------------------------------------------------------------------
        integer(int_kind) function list_create_two_way(this,list_max_len,list_beg)
!This function creates a bidirectional linked list.
!INPUT:
! # list_max_len - max length of the list;
! # list _beg - first entry number: range of entry numbers = [list_beg..(list_beg+list_max_len-1)];
        implicit none
        class(list_two_way_t):: this
        integer(int_kind), intent(in):: list_max_len,list_beg
        integer(int_kind):: i
        list_create_two_way=0
        if(list_max_len.ge.0) then
         if(this%list_max_length.lt.0) then
          this%list_max_length=0; this%first=-1; this%last=-1; this%current=-1; this%free_p=-1
          if(list_max_len.gt.0) then
           allocate(this%free_stack(0:list_max_len-1),STAT=i)
           if(i.ne.0) then; list_create_two_way=list_err_mem_alloc_failed; return; endif
           do i=0,list_max_len-1; free_stack(i)=list_beg+i; enddo; this%free_p=0 !stack of free entry numbers
           allocate(this%next(list_beg:list_beg+list_max_len-1),STAT=i)
           if(i.ne.0) then; deallocate(this%free_stack); list_create_two_way=list_err_mem_alloc_failed; return; endif
           allocate(this%prev(list_beg:list_beg+list_max_len-1),STAT=i)
           if(i.ne.0) then; deallocate(this%free_stack,this%next); list_create_two_way=list_err_mem_alloc_failed; return; endif
           this%list_max_length=list_max_len
          endif
         else
          list_create_two_way=list_err_list_exists
         endif
        else
         list_create_two_way=list_err_invalid_arg
        endif
        return
        end function list_create_two_way
!----------------------------------------------------
        integer(int_kind) function list_destroy(this)
!This function destroys a list.
        implicit none
        class(list_t):: this
        integer(int_kind):: i
        list_destroy=0
        if(this%list_max_length.ge.0) then
         this%list_max_length=-1; this%first=-1; this%last=-1; this%current=-1; this%free_p=-1
         if(allocated(this%free_stack)) then 
          deallocate(this%free_stack,STAT=i); if(i.ne.0) list_destroy=list_err_mem_free_failed
         endif
        else
         list_destroy=list_err_list_null
        endif
        return
        end function list_destroy
!------------------------------------------------------------
        integer(int_kind) function list_destroy_one_way(this)
!This function destroys a unidirectional linked list.
        implicit none
        class(list_one_way_t):: this
        integer(int_kind):: i
        list_destroy_one_way=0
        if(this%list_max_length.ge.0) then
         this%list_max_length=-1; this%first=-1; this%last=-1; this%current=-1; this%free_p=-1
         if(allocated(this%free_stack)) then 
          deallocate(this%free_stack,STAT=i); if(i.ne.0) list_destroy_one_way=list_err_mem_free_failed
         endif
         if(allocated(this%next)) then 
          deallocate(this%next,STAT=i); if(i.ne.0) list_destroy_one_way=list_err_mem_free_failed
         endif
        else
         list_destroy_one_way=list_err_list_null
        endif
        return
        end function list_destroy_one_way
!------------------------------------------------------------
        integer(int_kind) function list_destroy_two_way(this)
!This function destroys a bidirectional linked list.
        implicit none
        class(list_two_way_t):: this
        integer(int_kind):: i
        list_destroy_two_way=0
        if(this%list_max_length.ge.0) then
         this%list_max_length=-1; this%first=-1; this%last=-1; this%current=-1; this%free_p=-1
         if(allocated(this%free_stack)) then 
          deallocate(this%free_stack,STAT=i); if(i.ne.0) list_destroy_two_way=list_err_mem_free_failed
         endif
         if(allocated(this%next)) then 
          deallocate(this%next,STAT=i); if(i.ne.0) list_destroy_two_way=list_err_mem_free_failed
         endif
         if(allocated(this%prev)) then 
          deallocate(this%prev,STAT=i); if(i.ne.0) list_destroy_two_way=list_err_mem_free_failed
         endif
        else
         list_destroy_two_way=list_err_list_null
        endif
        return
        end function list_destroy_two_way
!-----------------------------------------------------------
        integer(int_kind) function list_get_max_length(this)
!This function returns the maximal length of a list (negative will mean an uninitialized list).
        implicit none
        class(list_t):: this
        list_get_max_length=this%list_max_length
        return
        end function list_get_max_length
!-------------------------------------------------------
        integer(int_kind) function list_get_length(this)
!This function returns the current length of a list (negative will mean an uninitialized list).
        implicit none
        class(list_t):: this
        if(this%list_max_length.gt.0) then
         list_get_length=this%free_p
        elseif(this%list_max_length.eq.0) then
         list_get_length=0
        else
         list_get_length=-1
        endif
        return
        end function list_get_length
!------------------------------------------------------------
        integer(int_kind) function list_reset(this,go_to_end)
!This function resets the current entry of a list to its beginning (or end).
        implicit none
        class(list_t):: this
        logical, intent(in), optional:: go_to_end
        logical res
        list_reset=0
        if(this%list_max_length.ge.0) then
         if(present(go_to_end)) then; res=go_to_end; else; res=.false.; endif
         if(res) then !go to the end of the list
          this%current=this%last
         else !go to the beginning of the list
          this%current=this%first
         endif
        else
         list_reset=list_err_list_null
        endif
        return
        end function list_reset
!-------------------------------------------------------------------
        integer(int_kind) function list_move_to_next(this,entry_num)
!This function moves the current entry pointer to the next entry.
!The number of the next entry is returned in <entry_num>.
        implicit none
        class(list_one_way_t):: this
        integer, intent(out):: entry_num
        list_move_to_next=0; entry_num=-1
        if(this%list_max_length.ge.0) then
         
        else
         list_move_to_next=list_err_list_null
        endif
        return
        end function list_move_to_next

       end module lists
