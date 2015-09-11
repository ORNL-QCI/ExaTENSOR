!Generic implementation of a stack based on OOP Fortran 2008.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2015/09/11
!Copyright (C) 2015 Dmitry I. Lyakh (email: quant4me@gmail.com)
!Copyright (C) 2015 Oak Ridge National Laboratory (UT-Battelle)
       module stack
        implicit none
        private
!PARAMETERS:
 !Generic:
        integer, private:: CONS_OUT=6        !default output
        logical, private:: DEBUG=.true.      !debugging mode
        integer, parameter, private:: INTD=4 !default integer kind
        integer, parameter, private:: INTL=8 !long integer kind
 !Error codes:
        integer(INTD), parameter, public:: STACK_SUCCESS=0         !success
        integer(INTD), parameter, public:: STACK_ERROR=-1          !unclassified error
        integer(INTD), parameter, public:: STACK_NOT_INITIALIZED=1 !error: stack is not initialized
        integer(INTD), parameter, public:: STACK_EMPTY=2           !error: stack is empty
        integer(INTD), parameter, public:: STACK_NOT_EMPTY=3       !error: stack is not empty
        integer(INTD), parameter, public:: STACK_MEM_ALLOC_FAIL=4  !error: memory allocation failed
        integer(INTD), parameter, public:: STACK_MEM_FREE_FAIL=5   !error: memory deallocation failed
        integer(INTD), parameter, public:: STACK_CORRUPTED=6       !error: stack got corrupted
!DERIVED TYPES:
 !Polymorphic stack entry:
        type, private:: stack_entry_t
         class(*), pointer, private:: item                     !stored item
         class(stack_entry_t), pointer, private:: prev=>NULL() !pointer to the previous stored item
         logical, private:: pointed_only                       !flag tht distinguish allocated items from associated
        end type stack_entry_t
 !Stack:
        type, public:: stack_t
         integer(INTL), private:: num_items=0_INTL             !total number of items in the stack
         class(stack_entry_t), pointer, private:: last=>NULL() !pointer to the last item
         logical, private:: corrupted=.false.                  !set to .TRUE. when the stack is corrupted
         contains
          procedure, public:: init=>stack_init   !initializes a new stack
          procedure, public:: state=>stack_state !returns the current state of the stack
          procedure, public:: push=>stack_push   !adds an item to the stack
          procedure, public:: pull=>stack_pull   !extracts the last item from the stack
        end type stack_t
!VISIBILITY:
        private stack_init
        private stack_state
        private stack_push
        private stack_pull

       contains
!---------------------------------------
        subroutine stack_init(this,ierr)
!Initializes a new stack.
        implicit none
        class(stack_t), intent(inout):: this        !inout: stack
        integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
        integer(INTD):: errc

        if(this%num_items.eq.0.and.(.not.(this%corrupted.or.associated(this%last)))) then
         errc=STACK_SUCCESS
        else
         if(this%corrupted) then
          errc=STACK_CORRUPTED
         else
          if(this%num_items.gt.0.and.associated(this%last)) then
           errc=STACK_NOT_EMPTY
          else
           this%corrupted=.true.; errc=STACK_CORRUPTED
          endif
         endif
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine stack_init
!----------------------------------------------------------
        subroutine stack_state(this,num_items,ierr,test_it)
!Returns the current state of a stack.
        implicit none
        class(stack_t), intent(in):: this           !in: stack
        integer(INTL), intent(out):: num_items      !out: total number of items in the stack
        integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
        logical, intent(in), optional:: test_it     !in: if .true., the stack linking will be tested
        integer(INTD):: errc
        integer(INTL):: i
        class(stack_entry_t), pointer:: sp

        num_items=-1
        if(.not.this%corrupted) then
         errc=STACK_SUCCESS
         if(this%num_items.gt.0.and.associated(this%last)) then
          if(present(test_it)) then
           if(test_it) then
            sp=>this%last
            do i=1,this%num_items-1
             if(associated(sp%prev)) then
              sp=>sp%prev
             else
              this%corrupted=.true.; errc=STACK_CORRUPTED
              exit
             endif
            enddo
            if(errc.eq.STACK_SUCCESS.and.associated(sp%prev)) then
             this%corrupted=.true.; errc=STACK_CORRUPTED
            endif
            nullify(sp)
           endif
          endif
          if(errc.eq.STACK_SUCCESS) num_items=this%num_items
         elseif(this%num_items.lt.0.or.associated(this%last)) then
          this%corrupted=.true.; errc=STACK_CORRUPTED
         endif
        else
         errc=STACK_CORRUPTED
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine stack_state
!-------------------------------------------------------
        subroutine stack_push(this,item,ierr,point_only)
!Adds an item to a stack.
        implicit none
        class(stack_t), intent(inout):: this        !in: stack
        class(*), intent(in):: item                 !in: item
        integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
        logical, intent(in), optional:: point_only  !in: if .true., the %item field will be just associated with <item> (default: false)
        integer(INTD):: errc
        logical:: only_point
        class(stack_entry_t), pointer:: sp

        if(.not.this%corrupted) then
         errc=STACK_SUCCESS
         if(present(point_only)) then; only_point=point_only; else; only_point=.false.; endif
         if(this%num_items.gt.0.and.associated(this%last)) then
          sp=>this%last; call add_new_entry(errc)
          if(errc.eq.STACK_SUCCESS.or.errc.eq.STACK_CORRUPTED) then
           this%last%prev=>sp
          else
           this%last=>sp
          endif
          sp=>NULL()
         elseif(this%num_items.eq.0.and.(.not.associated(this%last))) then
          call add_new_entry(errc)
          if(errc.eq.STACK_SUCCESS.or.errc.eq.STACK_CORRUPTED) then
           this%last%prev=>NULL()
          else
           this%last=>NULL()
          endif
         else
          this%corrupted=.true.; errc=STACK_CORRUPTED
         endif
        else
         errc=STACK_CORRUPTED
        endif
        if(present(ierr)) ierr=errc
        return

        contains

         subroutine add_new_entry(erc)
          integer(INTD), intent(out):: erc
          erc=0; allocate(this%last,STAT=erc) !allocate a new stack entry
          if(erc.eq.0) then
           if(only_point) then
            this%last%item=>item
           else
            allocate(this%last%item,source=item,STAT=erc)
            if(erc.ne.0) then
             deallocate(this%last,STAT=erc)
             if(erc.eq.0) then
              erc=STACK_MEM_ALLOC_FAIL
             else
              this%corrupted=.true.; erc=STACK_CORRUPTED
             endif
            endif
           endif
          else
           this%last=>NULL(); erc=STACK_MEM_ALLOC_FAIL
          endif
          if(erc.eq.0) erc=STACK_SUCCESS
          return
         end subroutine add_new_entry

        end subroutine stack_push

       end module stack
