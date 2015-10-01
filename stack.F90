!Generic implementation of a stack (OO Fortran 2003).
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2015/10/01
!Copyright (C) 2015 Dmitry I. Lyakh (email: quant4me@gmail.com)
!Copyright (C) 2015 Oak Ridge National Laboratory (UT-Battelle)
!LICENSE: GPL v.2
!DESCRIPTION:
! 1) Allocate a stack object via Fortran allocate;
! 2) Clean the stack at any point via .clean;
! 3) Push/pop any objects to the stack via .push and .pop;
! 4) Clean the stack at the end via .clean (if it is not empty);
! 5) Deallocate the stack object via Fortran deallocate;
! 6) While the stack is in scope, .state can return info on it.
!NOTES:
! a) An object pushed to the stack can be of any scalar type/class.
! b) It is illegal to push a NULL object.
! c) So far, it is illegal to push non-scalar objects (arrays).
       module stack
        use dil_kinds
        implicit none
        private
!PARAMETERS:
 !Generic:
        integer, private:: CONS_OUT=6        !default output
        logical, private:: VERBOSE=.true.    !verbosity for errors
        logical, private:: DEBUG=.true.      !debugging mode
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
         logical, private:: pointed_only                       !flag that distinguishes allocated items from associated
        end type stack_entry_t
 !Stack:
        type, public:: stack_t
         integer(INTL), private:: num_items=0_INTL             !total number of items in the stack
         class(stack_entry_t), pointer, private:: last=>NULL() !pointer to the last item
         logical, private:: corrupted=.false.                  !set to .TRUE. when the stack is corrupted
         contains
          procedure, public:: clean=>stack_clean !cleans the stack (pops out all items)
          procedure, public:: state=>stack_state !returns the current state of the stack
          procedure, public:: push=>stack_push   !adds an item to the stack
          procedure, public:: pop=>stack_pop     !extracts the last item from the stack
        end type stack_t
!VISIBILITY:
        private stack_clean
        private stack_state
        private stack_push
        private stack_pop

       contains
!----------------------------------------
        subroutine stack_clean(this,ierr)
!Cleans a stack (pops out and deletes all items).
        implicit none
        class(stack_t), intent(inout):: this        !inout: stack
        integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
        integer(INTD):: i,errc
        class(stack_entry_t), pointer:: sp

        if(this%num_items.eq.0.and.(.not.(this%corrupted.or.associated(this%last)))) then
         errc=STACK_SUCCESS
        else
         if(this%corrupted) then
          errc=STACK_CORRUPTED
         else
          if(this%num_items.gt.0.and.associated(this%last)) then
           errc=STACK_SUCCESS; sp=>this%last
           do while(associated(sp))
            this%last=>sp%prev; nullify(sp%prev)
            if(.not.sp%pointed_only) then
             deallocate(sp%item,STAT=i); if(i.ne.0) errc=STACK_MEM_FREE_FAIL
            else
             nullify(sp%item)
            endif
            deallocate(sp,STAT=i); if(i.ne.0) errc=STACK_MEM_FREE_FAIL
            this%num_items=this%num_items-1
            sp=>this%last
           enddo
           if(this%num_items.ne.0) then
            this%corrupted=.true.; if(errc.eq.STACK_SUCCESS) errc=STACK_CORRUPTED
           endif
          else
           this%corrupted=.true.; errc=STACK_CORRUPTED
          endif
         endif
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine stack_clean
!----------------------------------------------------------
        subroutine stack_state(this,num_items,ierr,test_it)
!Returns the current state of the stack.
        implicit none
        class(stack_t), intent(inout):: this        !in: stack
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
!Adds an item to a stack. An item can be an arbitrary object.
!It can either be cloned (default) or pointed to (point_only=TRUE).
        implicit none
        class(stack_t), intent(inout):: this        !in: stack
        class(*), intent(in), target:: item         !in: item
        integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
        logical, intent(in), optional:: point_only  !in: if TRUE, the %item field will be just associated with <item> (default: false)
        integer(INTD):: errc
        logical:: only_point
        class(stack_entry_t), pointer:: sp

        if(.not.this%corrupted) then
         errc=STACK_SUCCESS
         if(present(point_only)) then; only_point=point_only; else; only_point=.false.; endif
         if(this%num_items.gt.0.and.associated(this%last)) then !not the first item
          sp=>this%last; call add_new_entry(errc)
          if(errc.eq.STACK_SUCCESS) then
           this%last%prev=>sp
          elseif(errc.eq.STACK_CORRUPTED) then
           if(associated(this%last)) this%last%prev=>sp
          else
           this%last=>sp
          endif
          sp=>NULL()
         elseif(this%num_items.eq.0.and.(.not.associated(this%last))) then !first item
          call add_new_entry(errc)
          if(errc.eq.STACK_SUCCESS) then
           this%last%prev=>NULL()
          elseif(errc.eq.STACK_CORRUPTED) then
           if(associated(this%last)) this%last%prev=>NULL()
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
          erc=0; this%last=>NULL(); allocate(this%last,STAT=erc) !allocate a new stack entry as the last one
          if(erc.eq.0) then
           if(only_point) then
            this%last%item=>item; this%last%pointed_only=.true.
            this%num_items=this%num_items+1
           else
            allocate(this%last%item,source=item,STAT=erc)
            if(erc.eq.0) then
             this%last%pointed_only=.false.
             this%num_items=this%num_items+1
            else
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
!-------------------------------------------
        subroutine stack_pop(this,item,ierr)
!Pops up the last item out of the stack.
        implicit none
        class(stack_t), intent(inout):: this        !in: stack
        class(*), intent(out), pointer:: item       !out: item
        integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
        integer(INTD):: i,errc
        class(stack_entry_t), pointer:: sp

        errc=0; item=>NULL()
        if(.not.this%corrupted) then
         if(this%num_items.gt.0.and.associated(this%last)) then
          if(associated(this%last%item)) then
           item=>this%last%item; this%last%item=>NULL()
           sp=>this%last%prev; this%last%prev=>NULL()
           deallocate(this%last,STAT=i)
           if(i.eq.0) then
            this%num_items=this%num_items-1; this%last=>sp
            errc=STACK_SUCCESS
           else
            this%num_items=this%num_items-1; this%last=>sp
            errc=STACK_MEM_FREE_FAIL
           endif
           sp=>NULL()
          else !NULL pointers cannot be stored in the stack as items
           this%corrupted=.true.; errc=STACK_CORRUPTED
          endif
         elseif(this%num_items.lt.0.or.associated(this%last)) then
          this%corrupted=.true.; errc=STACK_CORRUPTED
         else
          errc=STACK_EMPTY
         endif
        else
         errc=STACK_CORRUPTED
        endif
        if(present(ierr)) ierr=errc
        return
        end subroutine stack_pop

       end module stack
!----------------------------------------------------------------------
!TESTING PART:
!----------------------------------------------------------------------
       module stack_test
        use dil_kinds
        use stack
        use timers, only: thread_wtime
        implicit none
        private
!VISIBILITY:
        public dil_test_stack

       contains

        function dil_test_stack(perf,dev_out) result(ierr)
         implicit none
         integer(INTD):: ierr                          !out: error code (0:success)
         real(8), intent(out):: perf                   !out: performance index
         integer(INTD), intent(in), optional:: dev_out !in: default output device
!-----------------------------------------------------
         integer(INTD), parameter:: MAX_ACTIONS=1000000

         type base_t
          integer(INTD):: int_field
          real(8):: real8_field
         end type base_t

         type, extends(base_t):: derv_t
          logical:: lgc_field
          type(base_t), pointer:: base_ptr
         end type derv_t

         integer(INTD):: jo,i,pushed,popped
         integer(INTL):: n
         type(base_t), target:: base
         type(derv_t):: derv
         type(stack_t), pointer:: my_stack
         class(*), pointer:: item
         class(base_t), pointer:: bp
         real(8):: rn,tm,sm

         ierr=0_INTD; perf=0d0
         if(present(dev_out)) then; jo=dev_out; else; jo=6_INTD; endif
         base=base_t(3_INTD,-1.134d0)
         derv%base_t=base_t(-11_INTD,-1.134d0); derv%lgc_field=.true.; derv%base_ptr=>base
         allocate(my_stack,STAT=ierr); if(ierr.ne.0) then; call test_quit(1_INTD); return; endif
         call my_stack%clean(ierr); if(ierr.ne.STACK_SUCCESS) then; call test_quit(2_INTD); return; endif
         pushed=0; popped=0; sm=0d0; tm=thread_wtime()
         do i=1,MAX_ACTIONS
          nullify(item)
          call random_number(rn)
          if(rn.gt.2d-1) then !push (more often)
           if(rn.gt.6d-1) then
            call my_stack%push(base,ierr,point_only=.false.)
           else
            call my_stack%push(derv,ierr,point_only=.false.)
           endif
           if(ierr.ne.STACK_SUCCESS) then; call test_quit(3_INTD); return; endif
           pushed=pushed+1
          else !pop (less often)
           call my_stack%state(n,ierr,test_it=.false.)
           if(ierr.ne.STACK_SUCCESS) then; call test_quit(4_INTD); return; endif
           call my_stack%pop(item,ierr)
           if(.not.((n.gt.0.and.ierr.eq.STACK_SUCCESS.and.associated(item)).or.(n.le.0.and.ierr.eq.STACK_EMPTY))) then
            call test_quit(5_INTD); return
           endif
           if(ierr.eq.STACK_SUCCESS) then
            select type(item)
            class is(base_t)
             bp=>item; sm=sm+bp%real8_field
            class default
             call test_quit(6_INTD); return
            end select
            popped=popped+1
           endif
          endif
         enddo
         ierr=STACK_SUCCESS
         do while(ierr.eq.STACK_SUCCESS)
          call my_stack%pop(item,ierr)
          if(ierr.eq.STACK_SUCCESS) then
           popped=popped+1; if(.not.associated(item)) then; call test_quit(7_INTD); return; endif
           select type(item)
           class is(base_t)
            bp=>item; sm=sm+bp%real8_field
           class default
            call test_quit(8_INTD); return
           end select
          endif
         enddo
         if(ierr.ne.STACK_EMPTY) then; call test_quit(9_INTD); return; endif
         if(pushed.ne.popped) then; call test_quit(10_INTD); return; endif
         tm=thread_wtime(tm); perf=dble(MAX_ACTIONS)/tm
!        write(jo,*) popped,sm !debug
         call my_stack%clean()
         deallocate(my_stack,STAT=ierr); if(ierr.ne.0) then; call test_quit(11_INTD); return; endif
         return

         contains

          subroutine test_quit(jerr)
           integer(INTD), intent(in):: jerr
           ierr=jerr
           if(ierr.ne.0) then
            write(jo,'("#ERROR(stack::dil_test_stack): Test failed: Error code ",i13)') ierr
            write(jo,'("Please contact the developer at QUANT4ME@GMAIL.COM")')
           endif
           if(associated(my_stack)) then
            call my_stack%clean(); deallocate(my_stack)
           endif
           return
          end subroutine test_quit

        end function dil_test_stack

       end module stack_test
