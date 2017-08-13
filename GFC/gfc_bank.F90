!Generic Fortran Containers (GFC): Bank of preallocated reusable objects of a specific type
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2017/08/13

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

       module gfc_bank
        use gfc_base
        use timers
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !output device
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
        integer(INTD), private:: DEBUG=0    !debugging level (0:none)
!TYPES:
 !Reusable object:
        type, public:: borrowed_object_t
         class(*), pointer, private:: object_p=>NULL() !pointer to a reusable object
         integer(INTD), private:: offset=-1            !object entry position in the object bank: [0..max]
         contains
          procedure, private:: reset=>BorrowedObjectReset          !resets the pointer either to another target or NULL
          procedure, public:: is_set=>BorrowedObjectIsSet          !returns TRUE if the object is set, FALSE otherwise
          procedure, public:: get_object=>BorrowedObjectGetObject  !returns a pointer to the stored (reusable) object for client's use
          procedure, private:: get_offset=>BorrowedObjectGetOffset !returns the offset of the borrowed object in the bank
          final:: borrowed_object_dtor
        end type borrowed_object_t
 !Object bank entry:
        type, private:: obj_bank_entry_t
         class(*), allocatable, public:: object !stored object
        end type obj_bank_entry_t
 !Bank of reusable objects:
        type, public:: object_bank_t
         type(obj_bank_entry_t), allocatable, private:: entries(:)                      !object bank entries with objects to loan
         integer(INTD), allocatable, private:: free_entries(:)                          !offsets of currently available entries
         integer(INTD), private:: free_left=0                                           !number of currently available entries left in the bank
         integer(INTD), private:: capacity=0                                            !total bank capacity
         procedure(gfc_allocate_scalar_i), nopass, pointer, private:: allocator=>NULL() !non-member generic allocator for each object
         procedure(gfc_destruct_i), nopass, pointer, private:: constructor=>NULL()      !non-member generic constructor for each object
         procedure(gfc_destruct_i), nopass, pointer, private:: destructor=>NULL()       !non-member generic destructor for each object
         contains
          procedure, public:: init=>ObjectBankInit          !initializes the object bank
          procedure, public:: borrow=>ObjectBankBorrow      !returns a resuable (borrowed) object from the bank
          procedure, public:: give_back=>ObjectBankGiveBack !puts the borrowed object back in the bank
          final:: object_bank_dtor                          !dtor
        end type object_bank_t
!VISIBILITY:
 !borrowed_object_t:
        private BorrowedObjectReset
        private BorrowedObjectIsSet
        private BorrowedObjectGetObject
        private BorrowedObjectGetOffset
        public borrowed_object_dtor
 !object_bank_t:
        private ObjectBankInit
        private ObjectBankBorrow
        private ObjectBankGiveBack
        public object_bank_dtor

       contains
!IMPLEMENTATION:
![borrowed_object_t]==========================================
        subroutine BorrowedObjectReset(this,ierr,new_p,offset)
         implicit none
         class(borrowed_object_t), intent(inout):: this   !inout: borrowed object
         integer(INTD), intent(out), optional:: ierr      !out: error code
         class(*), intent(in), pointer, optional:: new_p  !in: pointer to the newly borrowed bank object
         integer(INTD), intent(in), optional:: offset     !in: offset of the newly borrowed object in the bank
         integer(INTD):: errc

         errc=GFC_SUCCESS
         this%object_p=>NULL(); this%offset=-1
         if(present(new_p)) then
          if(present(offset)) then
           this%object_p=>new_p
           this%offset=offset
          else
           errc=GFC_INVALID_ARGS
          endif
         else
          if(present(offset)) errc=GFC_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine BorrowedObjectReset
!----------------------------------------------------------
        function BorrowedObjectIsSet(this,ierr) result(res)
         implicit none
         logical:: res                               !out: result
         class(borrowed_object_t), intent(in):: this !in: borrowed object
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS; res=associated(this%object_p)
         if(present(ierr)) ierr=errc
         return
        end function BorrowedObjectIsSet
!----------------------------------------------------------------
        function BorrowedObjectGetObject(this,ierr) result(obj_p)
         implicit none
         class(*), pointer:: obj_p                   !out: pointer to the borrowed object value
         class(borrowed_object_t), intent(in):: this !in: borrowed object
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS; obj_p=>this%object_p
         if(.not.associated(obj_p)) errc=GFC_ELEM_EMPTY
         if(present(ierr)) ierr=errc
         return
        end function BorrowedObjectGetObject
!-----------------------------------------------------------------
        function BorrowedObjectGetOffset(this,ierr) result(offset)
         implicit none
         integer(INTD):: offset                      !out: offset: >=0
         class(borrowed_object_t), intent(in):: this !in: borrowed object
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         offset=-1
         if(this%is_set(errc)) then
          offset=this%offset
         else
          if(errc.eq.GFC_SUCCESS) errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end function BorrowedObjectGetOffset
!--------------------------------------------
        subroutine borrowed_object_dtor(this)
         implicit none
         type(borrowed_object_t):: this

         this%object_p=>NULL()
         this%offset=-1
         return
        end subroutine borrowed_object_dtor
![object_bank_t]===============================================================
        subroutine ObjectBankInit(this,capacity,allocator_f,ierr,ctor_f,dtor_f)
!Initializes a bank of reusable entries to some capacity.
         implicit none
         class(object_bank_t), intent(inout):: this     !inout: bank of reusable objects
         integer(INTD), intent(in):: capacity           !in: desired capacity of the bank
         procedure(gfc_allocate_scalar_i):: allocator_f !in: non-member generic allocator procedure
         integer(INTD), intent(out), optional:: ierr    !out: error code
         procedure(gfc_destruct_i), optional:: ctor_f   !in: non-member generic constructor procedure
         procedure(gfc_destruct_i), optional:: dtor_f   !in: non-member generic destructor procedure
         integer(INTD):: errc,i
         logical:: fc

         errc=GFC_SUCCESS
         if(capacity.gt.0) then
          fc=.FALSE.; if(present(ctor_f)) fc=.TRUE.
          allocate(this%free_entries(0:capacity-1),STAT=errc)
          if(errc.eq.0) then
           allocate(this%entries(0:capacity-1),STAT=errc)
           if(errc.eq.0) then
            do i=0,capacity-1
             errc=allocator_f(this%entries(i)%object); if(errc.ne.0) exit
             if(fc) errc=ctor_f(this%entries(i)%object); if(errc.ne.0) exit
            enddo
            if(errc.eq.0) then
             this%capacity=0
             do i=capacity-1,0,-1
              this%free_entries(i)=this%capacity
              this%capacity=this%capacity+1
             enddo
             this%free_left=capacity
             this%allocator=>allocator_f
             if(fc) this%constructor=>ctor_f
             if(present(dtor_f)) this%destructor=>dtor_f
            else
             errc=GFC_MEM_ALLOC_FAILED
            endif
           else
            errc=GFC_MEM_ALLOC_FAILED
           endif
          else
           errc=GFC_MEM_ALLOC_FAILED
          endif
          if(errc.ne.GFC_SUCCESS) call object_bank_dtor(this)
         else
          errc=GFC_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ObjectBankInit
!----------------------------------------------------------------
        function ObjectBankBorrow(this,ierr) result(borrowed_obj)
!Borrows an object from the bank. The object is stored inside <borrowed_object_t>.
         implicit none
         type(borrowed_object_t):: borrowed_obj             !out: borrowed object (contains the stored object)
         class(object_bank_t), intent(inout), target:: this !inout: bank of reusable objects
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc,i
         class(*), pointer:: uptr

         errc=GFC_SUCCESS
         if(this%free_left.gt.0) then
          this%free_left=this%free_left-1
          i=this%free_entries(this%free_left)
          if(allocated(this%entries(i)%object)) then
           if(associated(this%constructor)) errc=this%constructor(this%entries(i)%object)
           if(errc.eq.0) then
            uptr=>this%entries(i)%object
            call borrowed_obj%reset(errc,uptr,i)
            uptr=>NULL()
           else
            errc=GFC_ACTION_FAILED
           endif
          else
           errc=GFC_CORRUPTED_CONT
          endif
         else
          errc=GFC_OVERFLOW
         endif
         if(present(ierr)) ierr=errc
         return
        end function ObjectBankBorrow
!------------------------------------------------------------
        subroutine ObjectBankGiveBack(this,borrowed_obj,ierr)
!Returns a borrowed object back to the bank.
         implicit none
         class(object_bank_t), intent(inout):: this            !inout: bank of reusable objects
         type(borrowed_object_t), intent(inout):: borrowed_obj !inout: defined borrowed object
         integer(INTD), intent(out), optional:: ierr
         integer(INTD):: errc,i

         if(borrowed_obj%is_set(errc)) then
          if(errc.eq.GFC_SUCCESS) then
           if(this%free_left.lt.this%capacity) then
            i=borrowed_obj%get_offset(errc)
            if(errc.eq.GFC_SUCCESS) then
             if(i.ge.0.and.i.lt.this%capacity) then
              this%free_entries(this%free_left)=i
              this%free_left=this%free_left+1
              call borrowed_obj%reset()
              if(associated(this%destructor)) errc=this%destructor(this%entries(i)%object)
             else
              errc=GFC_INVALID_ARGS
             endif
            endif
           else
            errc=GFC_INVALID_REQUEST
           endif
          endif
         else
          if(errc.eq.GFC_SUCCESS) errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ObjectBankGiveBack
!----------------------------------------
        subroutine object_bank_dtor(this)
         implicit none
         type(object_bank_t):: this
         integer(INTD):: i,errc

         if(allocated(this%entries)) then
          if(associated(this%destructor)) then
           do i=this%capacity-1,0,-1
            errc=this%destructor(this%entries(i)%object)
           enddo
          endif
          do i=this%capacity-1,0,-1
           deallocate(this%entries(i)%object)
          enddo
          deallocate(this%entries)
         endif
         if(allocated(this%free_entries)) deallocate(this%free_entries)
         this%capacity=0; this%free_left=0
         this%allocator=>NULL(); this%destructor=>NULL()
         return
        end subroutine object_bank_dtor

       end module gfc_bank
![TESTING]=============================================================
       module gfc_bank_test
        use gfc_base
        use gfc_bank
        use timers, only: thread_wtime
        implicit none
        private

        type, private:: my_type
         integer:: my_int
         real:: my_real
         character(9):: my_string='My string'
        end type my_type

        public test_gfc_bank

       contains

        function allocator(obj) result(ierr)
         implicit none
         integer(INTD):: ierr
         class(*), allocatable, intent(out):: obj

         allocate(my_type::obj,STAT=ierr)
         return
        end function allocator

        function test_gfc_bank(perf,dev_out) result(ierr)
         implicit none
         integer(INTD):: ierr
         real(8), intent(out):: perf
         integer(INTD), intent(in), optional:: dev_out
         integer(INTD), parameter:: BANK_SIZE=1024
         type(object_bank_t):: obank
         type(borrowed_object_t):: bob(BANK_SIZE)
         class(my_type), pointer:: my_ptr
         class(*), pointer:: uptr
         integer(INTD):: jo,i
         real(8):: tms,tm

         ierr=0; perf=0d0
         jo=6; if(present(dev_out)) jo=dev_out
         tms=thread_wtime()
         call obank%init(BANK_SIZE,allocator,ierr); if(ierr.ne.0) then; ierr=1; return; endif
         do i=1,BANK_SIZE
          bob(i)=obank%borrow(ierr); if(ierr.ne.0) then; ierr=2; return; endif
          uptr=>bob(i)%get_object(ierr); if(ierr.ne.0) then; ierr=3; return; endif
          my_ptr=>NULL(); select type(uptr); class is(my_type); my_ptr=>uptr; end select
          if(.not.associated(my_ptr)) then; ierr=4; return; endif
          my_ptr%my_int=i; my_ptr%my_real=dble(i)
         enddo
         do i=1,BANK_SIZE
          call obank%give_back(bob(i),ierr); if(ierr.ne.0) then; ierr=5; return; endif
         enddo
         call object_bank_dtor(obank)
         tm=thread_wtime(tms); perf=dble(BANK_SIZE)/tm
         return
        end function test_gfc_bank

       end module gfc_bank_test
