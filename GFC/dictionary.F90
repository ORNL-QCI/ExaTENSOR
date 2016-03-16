!Generic dictionary implementation (OO Fortran) based on AVL BST.
!So far it is not derived from GFC::Base.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2016/03/16

!Copyright (C) 2014-2016 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2016 Oak Ridge National Laboratory (UT-Battelle)

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

       module dictionary
!DESCRIPTION:
!#Dictionary items:
!  In order to store an item ({key;value}) in the dictionary,
!  one must call <search> with the corresponding action parameter
!  and provide the <key> and <value>. If the <value> is of derived data type,
!  its components may stay uninitialized because <search> can return
!  a pointer to the stored <value> that can be used for accessing/allocating its components.
!  In other words, <search> simply clones the <key> and <value> in the state
!  they are passed to <search>. After that, the stored <key> cannot be modified,
!  whereas the stored <value> can. The stored <value> can be retrieved by its <key>.
!  All <keys> in the same dictionary must be of the same class whereas
!  <values> can be of different classes/types.
!  Every dictionary has an internal iterator which is set every time a key is
!  found (and not deleted) or newly added. There are explicit moving methods
!  which move the position of the iterator up/down the binary tree. The subtree
!  rooted at the current iterator position can also be traversed or printed.
!#Key comparison function:
!  The key comparison function must be supplied to <search> (see cmp_key_func_i interface below).
!  The key comparison function must operate on unlimited polymorphic entities (<keys>)!
!#Item destructor function:
!  If <key> or <value> is of derived type with pointer/allocatable components,
!  a destructor function may need to be supplied (see destruct_func_i interface below).
!  The destructor function must operate on an unlimited polymorphic entity!
!  The destructor function must free all dynamic components of <key> or <value>,
!  but not the <key> or <value> themselves (but if it does, that should not cause a problem).
!#Printing dictionary items:
!  If one wants to print dictionary items, the item printing function
!  must be supplied (see print_func_i interface below). The item printing
!  function must operate on an unlimited polymorphic entity (item)!
        implicit none
        private
!PARAMETERS:
 !General:
        integer, private:: jo_dict=6      !default output device
        logical, private:: debug=.true.   !debug mode
        logical, private:: verbose=.true. !turns on/off verbose mode (errors, warnings, etc)
        integer, parameter, private:: LONGINT=8 !long integer length
 !Input named parameters:
  !Positioning:
        logical, parameter, public:: DICT_LEFT=.true.
        logical, parameter, public:: DICT_RIGHT=.false.
        logical, parameter, public:: DICT_SUCCESSOR=.true.
        logical, parameter, public:: DICT_PREDECESSOR=.false.
  !Search:
        integer, parameter, public:: DICT_JUST_SEARCH=0
        integer, parameter, public:: DICT_DELETE_IF_FOUND=1
        integer, parameter, public:: DICT_REPLACE_IF_FOUND=2
        integer, parameter, public:: DICT_ADD_IF_NOT_FOUND=3
        integer, parameter, public:: DICT_ADD_OR_MODIFY=4
        integer, parameter, public:: DICT_FETCH_IF_FOUND=5
 !Output named parameters:
  !General:
        integer, parameter, public:: DICT_SUCCESS=0
        integer, parameter, public:: DICT_UNKNOWN_ERR=-666
  !Search:
        integer, parameter, public:: DICT_KEY_FOUND=DICT_SUCCESS
        integer, parameter, public:: DICT_KEY_NOT_FOUND=1
  !Specific errors:
        integer, parameter, public:: DICT_CURRENT_ENTRY_NULL=-1
        integer, parameter, public:: DICT_ENTRY_NOT_EXIST=-2
        integer, parameter, public:: DICT_KEY_TYPE_MISMATCH=-3
        integer, parameter, public:: DICT_UNKNOWN_REQUEST=-4
        integer, parameter, public:: DICT_ABSENT_ARGUMENT=-5
        integer, parameter, public:: DICT_ALLOC_FAILED=-6
        integer, parameter, public:: DICT_FREE_FAILED=-7
        integer, parameter, public:: DICT_CORRUPTED=-8
        integer, parameter, public:: DICT_ERR_MAX_HEIGHT=-9
  !Key comparison results:
        integer, parameter, public:: DICT_KEY_EQ=0
        integer, parameter, public:: DICT_KEY_LT=-1
        integer, parameter, public:: DICT_KEY_GT=+1
        integer, parameter, public:: DICT_KEY_ERR=DICT_UNKNOWN_ERR
!TYPES:
 !Dictionary entry:
        type, private:: dict_entry_t
         class(*), allocatable, private:: entry_key
         class(*), allocatable, private:: entry_val
         class(dict_entry_t), pointer, private:: child_lt !left child (points only to an allocated target!)
         class(dict_entry_t), pointer, private:: child_gt !right child (points only to an allocated target!)
         class(dict_entry_t), pointer, private:: parent   !parent (points only to an allocated target!)
         integer, private:: balance_fac
        end type dict_entry_t
 !Dictionary:
        type, public:: dict_t
         class(dict_entry_t), pointer, private:: root=>NULL() !dictionary root (points only to an allocated target!)
         class(dict_entry_t), pointer, private:: curr_entry=>NULL() !current dictionary position (CDP): iterator
         integer(LONGINT), private:: num_entries=0 !total number of entries in the dictionary
         contains
          procedure, public:: volume=>dict_volume !returns the total number of entries in the dictionary
          procedure, public:: reset=>dict_reset !resets the current dictionary position (CDP) to the root entry
          procedure, public:: move_down=>dict_move_down !the current dictionary position (CDP) is moved downward the tree
          procedure, public:: move_up=>dict_move_up !the current dictionary position (CDP) is moved upward the tree
          procedure, public:: next_in_order=>dict_next_in_order !moves CDP to the in-order successor/predecessor
          procedure, public:: traverse_subtree=>dict_traverse !traverses the dictionary subtreee starting from the current position
          procedure, public:: destroy=>dict_destroy !destroys the dictionary
          procedure, public:: print_subtree=>dict_print !prints a subtree of the dictionary starting at CDP
          procedure, public:: print_entry=>dict_print_entry !print the current (CDP) entry
          procedure, public:: search=>dict_search !searches for an item by a given key with further optional actions
        end type dict_t
!INTERFACES:
        abstract interface
 !Generic key comparison function:
         integer function cmp_key_func_i(key1,key2)
          class(*), intent(in):: key1
          class(*), intent(in):: key2
         end function cmp_key_func_i
 !Generic printing function:
         integer function print_func_i(dev_id,item)
          integer, intent(in):: dev_id !device ID (6: screen)
          class(*), intent(in):: item !item to print
         end function print_func_i
 !Generic destructor:
         recursive function destruct_func_i(item) result(ierr)
          class(*), intent(inout):: item
          integer:: ierr
         end function destruct_func_i
        end interface
!OBJECT VISIBILITY:
 !Procedures:
        public dict_verbose
        private dict_volume
        private dict_reset
        private dict_move_down
        private dict_move_up
        private dict_next_in_order
        private dict_traverse
        private dict_destroy
        private dict_print
        private dict_print_entry
        private dict_search
        private print_key
 !Interfaces:
        public cmp_key_func_i
        public print_func_i
        public destruct_func_i
!-------------------------------------------
       contains
!SERVICE PROCEDURES:
        subroutine dict_verbose(verb,dev_id)
!Sets the verbosity mode for this module.
!INPUT:
! # verb: verbosity mode (.true./.false);
! # dev_id: optional output device id.
         implicit none
         logical, intent(in):: verb
         integer, intent(in), optional:: dev_id
         verbose=verb; if(present(dev_id)) jo_dict=dev_id
         return
        end subroutine dict_verbose
!METHODS:
!--------------------------------------------------
        integer(LONGINT) function dict_volume(this)
         implicit none
         class(dict_t), intent(in):: this
         dict_volume=this%num_entries
         return
        end function dict_volume
!----------------------------------------
        integer function dict_reset(this)
!Resets the current dictionary position to the root entry: CDP=>ROOT.
         implicit none
         class(dict_t):: this
         dict_reset=DICT_SUCCESS
         if(associated(this%root)) then
          this%curr_entry=>this%root
         else
          this%curr_entry=>NULL()
         endif
         return
        end function dict_reset
!-------------------------------------------------
        integer function dict_move_down(this,left)
!The current dictionary position is moved downward the tree (left or right): CDP=>CHILD(left,right)
         implicit none
         class(dict_t):: this
         logical, intent(in):: left !chooses between the left (.true.) and right (.false.) children
         dict_move_down=DICT_SUCCESS
         if(associated(this%curr_entry)) then
          if(left) then !move to the left
           if(associated(this%curr_entry%child_lt)) then
            this%curr_entry=>this%curr_entry%child_lt
           else
            dict_move_down=DICT_ENTRY_NOT_EXIST
           endif
          else !move to the right
           if(associated(this%curr_entry%child_gt)) then
            this%curr_entry=>this%curr_entry%child_gt
           else
            dict_move_down=DICT_ENTRY_NOT_EXIST
           endif
          endif
         else
          dict_move_down=DICT_CURRENT_ENTRY_NULL
         endif
         return
        end function dict_move_down
!------------------------------------------
        integer function dict_move_up(this)
!The current dictionary position is moved upward the tree (to the parental entry): CDP=>PARENT
         implicit none
         class(dict_t):: this
         dict_move_up=DICT_SUCCESS
         if(associated(this%curr_entry)) then
          if(associated(this%curr_entry%parent)) then
           this%curr_entry=>this%curr_entry%parent
          else
           dict_move_up=DICT_ENTRY_NOT_EXIST
          endif
         else
          dict_move_up=DICT_CURRENT_ENTRY_NULL
         endif
         return
        end function dict_move_up
!-----------------------------------------------------
        integer function dict_next_in_order(this,succ)
!Moves CDP to the in-order successor/predecessor: CDP=>successor/predecessor
!If there is none, CDP is left unchanged.
         implicit none
         class(dict_t):: this
         logical, intent(in):: succ !switches between successor (.true.) and predecessor (.false.)
         class(dict_entry_t), pointer:: curr

         dict_next_in_order=DICT_SUCCESS
         if(associated(this%curr_entry)) then
          curr=>this%curr_entry
          if(succ) then !find in-order successor
           if(associated(curr%child_gt)) then
            curr=>curr%child_gt
            do while(associated(curr%child_lt)); curr=>curr%child_lt; enddo
            this%curr_entry=>curr; curr=>NULL(); return
           else
            dict_next_in_order=DICT_ENTRY_NOT_EXIST
            do while(associated(curr%parent))
             if(associated(curr%parent%child_lt,curr)) then
              this%curr_entry=>curr%parent; curr=>NULL(); dict_next_in_order=DICT_SUCCESS; return
             else
              curr=>curr%parent
             endif
            enddo
           endif
          else !find in-order predecessor
           if(associated(curr%child_lt)) then
            curr=>curr%child_lt
            do while(associated(curr%child_gt)); curr=>curr%child_gt; enddo
            this%curr_entry=>curr; curr=>NULL(); return
           else
            dict_next_in_order=DICT_ENTRY_NOT_EXIST
            do while(associated(curr%parent))
             if(associated(curr%parent%child_gt,curr)) then
              this%curr_entry=>curr%parent; curr=>NULL(); dict_next_in_order=DICT_SUCCESS; return
             else
              curr=>curr%parent
             endif
            enddo
           endif
          endif
         else
          dict_next_in_order=DICT_CURRENT_ENTRY_NULL
         endif
         curr=>NULL()
         return
        end function dict_next_in_order
!------------------------------------------------------------------------
        integer function dict_traverse(this,subtree_size,max_height,iter)
!Traverses the dictionary subtree starting from the current dictionary position (CDP).
!The current dictionary position (CDP) is left unchanged at the end.
!If <iter> is present, then this function returns after each new dictionary entry,
!keeping the current dictionary position in <iter> (for internal use only).
!INPUT:
! # this: dictionary (passed-object dummy);
! # subtree_size: If <iter> is present, <subtree_size> MUST be initialized to ZERO in the first call;
!OUTPUT:
! # subtree_size: total number of dictionary entries traversed (at the end);
! # iter: iterator (opaque object for internal use only).
         implicit none
         class(dict_t):: this
         integer(LONGINT), intent(inout):: subtree_size
         integer, intent(inout):: max_height
         class(dict_entry_t), pointer, optional:: iter
         integer, save:: drct,lev_p,left
         class(dict_entry_t), pointer, save:: curr=>NULL()
         logical:: no_iter,no_go

!$OMP THREADPRIVATE(drct,lev_p,left,curr)

!        if(debug) write(jo_dict,'("Entered traverse:")') !debug
         dict_traverse=DICT_SUCCESS; if(present(iter)) then; no_iter=.false.; else; no_iter=.true.; endif
         if(no_iter) subtree_size=0_LONGINT
         if(subtree_size.eq.0_LONGINT) then !no iterator or first call with iterator
          if(associated(this%curr_entry)) then
           curr=>this%curr_entry; drct=1; lev_p=0; left=-1; max_height=1; no_go=.false.
          else
           curr=>NULL(); max_height=0; no_go=.true.
          endif
         else
          if(.not.associated(curr)) then; dict_traverse=DICT_UNKNOWN_REQUEST; return; endif
          no_go=.false.
         endif
         if(.not.no_go) then
          do while(lev_p.ge.0)
           if(.not.allocated(curr%entry_key)) then
            if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_traverse): item with no key detected!")')
            dict_traverse=DICT_CORRUPTED; return
           endif
!          if(debug) then; write(jo_dict,'("Lev/Drct ",i3,"/",i2": ")',advance='no') lev_p,drct; call print_key(curr); endif !debug
           if(drct.gt.0) then !going downwards
            if(left.lt.0) then !going downwards left
             if(associated(curr%child_lt)) then !proceed one level down to the left
!             if(debug) write(jo_dict,'("Jump left")') !debug
              curr=>curr%child_lt; lev_p=lev_p+1; max_height=max(lev_p+1,max_height)
             else
              left=-left !switch to the right branch on the same level
             endif
            else !going downwards right
             if(associated(curr%child_gt)) then !proceed one level down to the right
!             if(debug) write(jo_dict,'("Jump right")') !debug
              curr=>curr%child_gt; left=-1; lev_p=lev_p+1; max_height=max(lev_p+1,max_height)
             else !switch the direction to upwards
              drct=-drct
             endif
            endif
           else !going upwards
            subtree_size=subtree_size+1_LONGINT
            if(iabs(curr%balance_fac).gt.1) then !trap
             if(verbose) write(jo_dict,'("#WARNING(dictionary::dict_traverse): unbalanced entry: ",i13)') curr%balance_fac
            endif
!           if(debug) write(jo_dict,'("Jump up")') !debug
            lev_p=lev_p-1; if(.not.no_iter) iter=>curr
            if(lev_p.ge.0) then
             if(associated(curr%parent)) then 
              if(associated(curr%parent%child_lt,curr)) then
               left=+1; drct=-drct !switch to the right branch and change direction to downwards
              endif
              curr=>curr%parent
             else
              if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_traverse): absent parent!")')
              dict_traverse=DICT_CORRUPTED; return
             endif
            endif
            if(.not.no_iter) return
           endif
          enddo
          curr=>NULL()
          if(no_iter.and.associated(this%curr_entry,this%root)) then !trap
           if(subtree_size.ne.this%num_entries) then
            if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_traverse): wrong total number of entries: ",i13,1x,i13)')&
             &subtree_size,this%num_entries
            dict_traverse=DICT_CORRUPTED; return
           endif
          endif
         else
          dict_traverse=DICT_CURRENT_ENTRY_NULL
         endif
         if((.not.no_iter).and.dict_traverse.eq.DICT_SUCCESS) dict_traverse=DICT_ENTRY_NOT_EXIST !end of iterations
!        if(debug) write(jo_dict,'("Exited traverse: max height = ",i11)') max_height !debug
         return
        end function dict_traverse
!------------------------------------------------------------------------------
        integer function dict_destroy(this,destruct_key_func,destruct_val_func)
!Destroys a dictionary.
!INPUT:
! # this: dictionary;
! # destruct_key_func: key destructor;
! # destruct_val_func: value destructor;
!OUTPUT:
! # this: empty dictionary (and freed memory);
!NOTES:
! # This function tries to destroy all dictionary entries, regardless of errors.
         implicit none
         class(dict_t):: this
         procedure(destruct_func_i), optional:: destruct_key_func,destruct_val_func
         class(dict_entry_t), pointer:: iter
         integer(LONGINT):: subtree_size
         integer i,max_height

         dict_destroy=DICT_SUCCESS
         if(associated(this%root)) then
          i=this%reset(); subtree_size=0_LONGINT; iter=>NULL()
          do while(this%traverse_subtree(subtree_size,max_height,iter).eq.DICT_SUCCESS)
           if(associated(iter)) then
            if(present(destruct_key_func)) then
             i=destruct_key_func(iter%entry_key)
             if(i.ne.DICT_SUCCESS) then
              if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_destroy): key destructor failed: error ",i11)') i
              dict_destroy=DICT_FREE_FAILED
             endif
            endif
            if(allocated(iter%entry_key)) then
             deallocate(iter%entry_key,STAT=i)
             if(i.ne.DICT_SUCCESS) then
              if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_destroy): key deallocation failed: error ",i11)') i
              dict_destroy=DICT_FREE_FAILED
             endif
            endif
            if(present(destruct_val_func)) then
             i=destruct_val_func(iter%entry_val)
             if(i.ne.DICT_SUCCESS) then
              if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_destroy): value destructor failed: error ",i11)') i
              dict_destroy=DICT_FREE_FAILED
             endif
            endif
            if(allocated(iter%entry_val)) then
             deallocate(iter%entry_val,STAT=i)
             if(i.ne.DICT_SUCCESS) then
              if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_destroy): value deallocation failed: error ",i11)') i
              dict_destroy=DICT_FREE_FAILED
             endif
            endif
            deallocate(iter,STAT=i)
            if(i.ne.DICT_SUCCESS) then
             if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_destroy): item deallocation failed: error ",i11)') i
             dict_destroy=DICT_FREE_FAILED
            endif
           else
            if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_destroy): null iterator!")')
            dict_destroy=DICT_CORRUPTED
           endif
          enddo
          iter=>NULL()
         endif
         this%root=>NULL(); this%curr_entry=>NULL(); this%num_entries=0_LONGINT
         return
        end function dict_destroy
!----------------------------------------------------------
        integer function dict_print(this,dev_id,print_func)
!Prints a dictionary subtree starting from the current dictionary position (CDP).
         implicit none
         class(dict_t):: this
         integer, intent(in):: dev_id
         procedure(print_func_i):: print_func
         class(dict_entry_t), pointer:: iter
         integer(LONGINT):: subtree_size
         integer i,max_height

         dict_print=DICT_SUCCESS
         write(dev_id,'("#PRINTING DICTIONARY SUBTREE:")')
         if(associated(this%curr_entry)) then
          subtree_size=0_LONGINT; iter=>NULL()
          do while(this%traverse_subtree(subtree_size,max_height,iter).eq.DICT_SUCCESS)
           write(dev_id,'("<item>")')
           if(debug) then !debug begin
            write(dev_id,'(1x,"Balance: ",i6)') iter%balance_fac !debug
            if(associated(iter%parent)) then
             write(dev_id,'(1x,"Parent: ")',advance='no'); i=print_func(dev_id,iter%parent%entry_key) !debug
            endif
            if(associated(iter%child_lt)) then
             write(dev_id,'(1x,"Left: ")',advance='no'); i=print_func(dev_id,iter%child_lt%entry_key) !debug
            endif
            if(associated(iter%child_gt)) then
             write(dev_id,'(1x,"Right: ")',advance='no'); i=print_func(dev_id,iter%child_gt%entry_key) !debug
            endif
           endif !debug end
           write(dev_id,'(1x,"<key>")')
           i=print_func(dev_id,iter%entry_key) !print key
           if(i.ne.0) then
            if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_print): key printing failed: error ",i11)') i
            dict_print=DICT_UNKNOWN_ERR
           endif
           write(dev_id,'(1x,"</key>")')
           write(dev_id,'(1x,"<value>")')
           i=print_func(dev_id,iter%entry_val) !print value
           if(i.ne.0) then
            if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_print): value printing failed: error ",i11)') i
            dict_print=DICT_UNKNOWN_ERR
           endif
           write(dev_id,'(1x,"</value>")')
           write(dev_id,'("</item>")')
          enddo
          iter=>NULL()
          if(dict_print.eq.DICT_SUCCESS) then
           write(dev_id,'("#DONE SUCCESSFULLY: Total number of entries = ",i13)') subtree_size
          else
           write(dev_id,'("#DONE WITH ERRORS: Total number of entries = ",i13)') subtree_size
          endif
         else
          write(dev_id,'("#DONE SUCCESSFULLY: Empty.")')
          dict_print=DICT_CURRENT_ENTRY_NULL
         endif
         return
        end function dict_print
!----------------------------------------------------------------
        integer function dict_print_entry(this,dev_id,print_func)
!Prints a dictionary entry at the current dictionary position (CDP).
         implicit none
         class(dict_t):: this
         integer, intent(in):: dev_id
         procedure(print_func_i):: print_func
         integer i

         dict_print_entry=DICT_SUCCESS
         write(dev_id,'("<item_current>")')
         if(debug) then !debug begin
          write(dev_id,'(1x,"Balance: ",i6)') this%curr_entry%balance_fac !debug
          if(associated(this%curr_entry%parent)) then
           write(dev_id,'(1x,"Parent: ")',advance='no'); i=print_func(dev_id,this%curr_entry%parent%entry_key) !debug
          endif
          if(associated(this%curr_entry%child_lt)) then
           write(dev_id,'(1x,"Left: ")',advance='no'); i=print_func(dev_id,this%curr_entry%child_lt%entry_key) !debug
          endif
          if(associated(this%curr_entry%child_gt)) then
           write(dev_id,'(1x,"Right: ")',advance='no'); i=print_func(dev_id,this%curr_entry%child_gt%entry_key) !debug
          endif
         endif !debug end
         write(dev_id,'(1x,"<key>")')
         i=print_func(dev_id,this%curr_entry%entry_key) !print key
         if(i.ne.0) then
          if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_print_entry): key printing failed: error ",i11)') i
          dict_print_entry=DICT_UNKNOWN_ERR
         endif
         write(dev_id,'(1x,"</key>")')
         write(dev_id,'(1x,"<value>")')
         i=print_func(dev_id,this%curr_entry%entry_val) !print value
         if(i.ne.0) then
          if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_print_entry): value printing failed: error ",i11)') i
          dict_print_entry=DICT_UNKNOWN_ERR
         endif
         write(dev_id,'(1x,"</value>")')
         write(dev_id,'("</item_current>")')
         return
        end function dict_print_entry
!-----------------------------------------------------------------------------------------------------------------------------
        integer function dict_search(this,action,cmp_key_func,item_key,value_in,value_out,destruct_key_func,destruct_val_func)
!Looks up a given key in the dictionary with optional actions.
!If the key is found (and not deleted) or newly added, then the current dictionary position (CDP) will point to that item,
!otherwise CDP will not change.
!INPUT:
! # this: dictionary (passed-object dummy);
! # action: requested action (see action parameters at the top of this module);
! # cmp_key_func: key comparison function (must return: {DICT_KEY_LT, DICT_KEY_GT, DICT_KEY_EQ, DICT_KEY_ERR});
! # item_key: scalar key;
! # value_in: optional (scalar) item value (to be stored with the key by cloning);
! # destruct_key_func: key destructor;
! # destruct_val_func: value destructor;
!OUTPUT:
! # dict_search: {DICT_KEY_FOUND, DICT_KEY_NOT_FOUND, specific errors} (see output parameters at the top of this module);
! # value_out: when fetching, the <value_out> poly-pointer will point to the value found by the key (NULL otherwise);
! # this: possibly modified dictionary.
         implicit none
         class(dict_t), target:: this
         integer, intent(in):: action
         procedure(cmp_key_func_i):: cmp_key_func
         class(*):: item_key
         class(*), optional:: value_in
         class(*), pointer, optional:: value_out
         procedure(destruct_func_i), optional:: destruct_key_func,destruct_val_func
         class(dict_entry_t), pointer:: curr,old_cdp,leave,term
         integer:: i,j,act,lev_p,grow,ierr

         dict_search=DICT_KEY_NOT_FOUND; if(present(value_out)) value_out=>NULL()
!Look up the key:
         if(associated(this%root)) then
          curr=>this%root; lev_p=0
          sloop: do
           i=cmp_key_func(item_key,curr%entry_key)
           if(i.eq.DICT_KEY_LT) then
            if(associated(curr%child_lt)) then; curr=>curr%child_lt; lev_p=lev_p+1; else; exit sloop; endif
           elseif(i.eq.DICT_KEY_GT) then
            if(associated(curr%child_gt)) then; curr=>curr%child_gt; lev_p=lev_p+1; else; exit sloop; endif
           elseif(i.eq.DICT_KEY_EQ) then
            dict_search=DICT_KEY_FOUND; exit sloop
           else
            if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): key comparison failed: error ",i11)') i
            dict_search=DICT_KEY_TYPE_MISMATCH; curr=>NULL(); return
           endif
          enddo sloop
         else
          i=0; curr=>NULL(); lev_p=-1
         endif
!Action:
         if(action.eq.DICT_ADD_OR_MODIFY) then
          if(dict_search.eq.DICT_KEY_NOT_FOUND) then
           act=DICT_ADD_IF_NOT_FOUND
          elseif(dict_search.eq.DICT_KEY_FOUND) then
           act=DICT_REPLACE_IF_FOUND
          endif
         else
          act=action
         endif
         select case(act) !Process the action
         case(DICT_JUST_SEARCH) !no action
          if(dict_search.eq.DICT_KEY_FOUND) this%curr_entry=>curr
         case(DICT_FETCH_IF_FOUND) !return the pointer to the stored <value> if found
          if(dict_search.eq.DICT_KEY_FOUND) then
           this%curr_entry=>curr
           if(present(value_out)) then
            value_out=>curr%entry_val
           else
            if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): fetch: absent value return pointer!")')
            dict_search=DICT_ABSENT_ARGUMENT; curr=>NULL(); return
           endif
          endif
         case(DICT_REPLACE_IF_FOUND) !replace the stored <value> if found
          if(dict_search.eq.DICT_KEY_FOUND) then
           if(present(destruct_val_func)) then
            j=destruct_val_func(curr%entry_val)
            if(j.ne.0) then
             if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): replace: dictionary entry value destruction failed!")')
             dict_search=DICT_FREE_FAILED; curr=>NULL(); return
            endif
           endif
           deallocate(curr%entry_val,STAT=j)
           if(j.ne.0) then
            if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): replace: dictionary entry value deallocation failed!")')
            dict_search=DICT_FREE_FAILED; curr=>NULL(); return
           endif
           if(present(value_in)) then
            allocate(curr%entry_val,source=value_in,STAT=j)
            if(j.ne.0) then
             if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): replace: dictionary entry value allocation failed!")')
             dict_search=DICT_ALLOC_FAILED; curr=>NULL(); return
            endif
           else
            if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): replace: absent value!")')
            dict_search=DICT_ABSENT_ARGUMENT; curr=>NULL(); return
           endif
           this%curr_entry=>curr
           if(present(value_out)) value_out=>curr%entry_val
          endif
         case(DICT_DELETE_IF_FOUND) !delete the item if found
          if(dict_search.eq.DICT_KEY_FOUND) then
           if(associated(this%curr_entry)) then
            if(associated(this%curr_entry,curr)) then; old_cdp=>NULL(); else; old_cdp=>this%curr_entry; endif
           else
            old_cdp=>NULL()
           endif
           this%curr_entry=>curr
           if(associated(curr%child_lt).and.associated(curr%child_gt)) then !both subtrees are present
            if(curr%balance_fac.le.0) then !right subtree is taller or equal
             grow=-1; j=this%next_in_order(.true.) !find in-order successor
             if(j.eq.DICT_SUCCESS) then
              if(associated(this%curr_entry%child_gt)) then
               leave=>this%curr_entry%child_gt
              else
               leave=>NULL()
              endif
             else
              if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): delete: next in-order not found!")')
              dict_search=DICT_CORRUPTED
             endif
            else !left subtree is taller
             grow=+1; j=this%next_in_order(.false.) !find in-order predecessor
             if(j.eq.DICT_SUCCESS) then
              if(associated(this%curr_entry%child_lt)) then
               leave=>this%curr_entry%child_lt
              else
               leave=>NULL()
              endif
             else
              if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): delete: next in-order not found!")')
              dict_search=DICT_CORRUPTED
             endif
            endif
            if(dict_search.eq.DICT_KEY_FOUND) then
             if(associated(this%curr_entry%parent,curr)) then
              term=>NULL()
             else
              term=>this%curr_entry%parent
             endif
             if(associated(curr%parent)) then
              if(associated(curr%parent%child_lt,curr)) then
               curr%parent%child_lt=>this%curr_entry
              elseif(associated(curr%parent%child_gt,curr)) then
               curr%parent%child_gt=>this%curr_entry
              else
               if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): delete: lost parent!")')
               dict_search=DICT_CORRUPTED
              endif
             else
              if(associated(this%root,curr)) then
               this%root=>this%curr_entry
              else
               if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): delete: lost root!")')
               dict_search=DICT_CORRUPTED
              endif
             endif
             if(dict_search.eq.DICT_KEY_FOUND) then
              if(grow.gt.0) then !reducing the left subtree
               if(associated(term)) then
                curr%child_lt%parent=>this%curr_entry
                this%curr_entry%child_lt=>curr%child_lt
                if(associated(leave)) then
                 term%child_gt=>leave; leave%parent=>term
                else
                 term%child_gt=>NULL()
                endif
               endif
               curr%child_gt%parent=>this%curr_entry
               this%curr_entry%child_gt=>curr%child_gt
              elseif(grow.lt.0) then !reducing the right subtree
               if(associated(term)) then
                curr%child_gt%parent=>this%curr_entry
                this%curr_entry%child_gt=>curr%child_gt
                if(associated(leave)) then
                 term%child_lt=>leave; leave%parent=>term
                else
                 term%child_lt=>NULL()
                endif
               endif
               curr%child_lt%parent=>this%curr_entry
               this%curr_entry%child_lt=>curr%child_lt
              endif
              if(associated(curr%parent)) then
               this%curr_entry%parent=>curr%parent
              else
               this%curr_entry%parent=>NULL()
              endif
              this%curr_entry%balance_fac=curr%balance_fac
              if(associated(term)) then; this%curr_entry=>term; grow=-grow; term=>NULL(); endif
              if(associated(leave)) leave=>NULL()
             endif
            endif
           else !at least one subtree is absent
            if(associated(curr%child_lt)) then !left subtree is present (a leave)
             if(associated(curr%parent)) then
              if(associated(curr%parent%child_lt,curr)) then
               curr%parent%child_lt=>curr%child_lt; this%curr_entry=>curr%parent; grow=+1
              elseif(associated(curr%parent%child_gt,curr)) then
               curr%parent%child_gt=>curr%child_lt; this%curr_entry=>curr%parent; grow=-1
              else
               if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): delete: lost parent!")')
               dict_search=DICT_CORRUPTED
              endif
              curr%child_lt%parent=>curr%parent
             else
              if(associated(this%root,curr)) then
               this%root=>curr%child_lt; this%root%parent=>NULL(); grow=0
              else
               if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): delete: lost root!")')
               dict_search=DICT_CORRUPTED
              endif
             endif
            elseif(associated(curr%child_gt)) then !right subtree is present (a leave)
             if(associated(curr%parent)) then
              if(associated(curr%parent%child_lt,curr)) then
               curr%parent%child_lt=>curr%child_gt; this%curr_entry=>curr%parent; grow=+1
              elseif(associated(curr%parent%child_gt,curr)) then
               curr%parent%child_gt=>curr%child_gt; this%curr_entry=>curr%parent; grow=-1
              else
               if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): delete: lost parent!")')
               dict_search=DICT_CORRUPTED
              endif
              curr%child_gt%parent=>curr%parent
             else
              if(associated(this%root,curr)) then
               this%root=>curr%child_gt; this%root%parent=>NULL(); grow=0
              else
               if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): delete: lost root!")')
               dict_search=DICT_CORRUPTED
              endif
             endif
            else !both subtrees are absent
             if(associated(curr%parent)) then
              if(associated(curr%parent%child_lt,curr)) then
               curr%parent%child_lt=>NULL(); this%curr_entry=>curr%parent; grow=+1
              elseif(associated(curr%parent%child_gt,curr)) then
               curr%parent%child_gt=>NULL(); this%curr_entry=>curr%parent; grow=-1
              else
               if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): delete: lost parent!")')
               dict_search=DICT_CORRUPTED
              endif
             else
              if(associated(this%root,curr)) then
               this%root=>NULL(); this%curr_entry=>NULL(); grow=0
              else
               if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): delete: lost root!")')
               dict_search=DICT_CORRUPTED
              endif
             endif
            endif
           endif
           if(present(destruct_key_func)) then
            j=destruct_key_func(curr%entry_key)
            if(j.ne.0) then
             if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): delete: key destructor failed!")')
             dict_search=DICT_FREE_FAILED
            endif
           endif
           if(allocated(curr%entry_key)) then
            deallocate(curr%entry_key,STAT=j)
            if(j.ne.0) then
             if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): delete: key deallocation failed!")')
             dict_search=DICT_FREE_FAILED
            endif
           endif
           if(present(destruct_val_func)) then
            j=destruct_val_func(curr%entry_val)
            if(j.ne.0) then
             if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): delete: value destructor failed!")')
             dict_search=DICT_FREE_FAILED
            endif
           endif
           if(allocated(curr%entry_val)) then
            deallocate(curr%entry_val,STAT=j)
            if(j.ne.0) then
             if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): delete: value deallocation failed!")')
             dict_search=DICT_FREE_FAILED
            endif
           endif
           deallocate(curr,STAT=j)
           if(j.ne.0) then
            if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): delete: entry deallocation failed!")')
            dict_search=DICT_FREE_FAILED
           endif
           this%num_entries=this%num_entries-1_LONGINT
           if(dict_search.eq.DICT_KEY_FOUND.and.grow.ne.0) then
            ierr=-1; call rebalance(this%curr_entry,grow,ierr)
            if(ierr.ne.0) then
             if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): delete: rebalancing failed: ",i11)') ierr
             dict_search=DICT_CORRUPTED
            endif
           endif
           if(associated(old_cdp)) then; this%curr_entry=>old_cdp; old_cdp=>NULL(); else; this%curr_entry=>NULL(); endif
          endif
         case(DICT_ADD_IF_NOT_FOUND) !add a new item if the key is not found
          if(dict_search.eq.DICT_KEY_NOT_FOUND) then
           if(lev_p.ge.0) then
            if(i.eq.DICT_KEY_LT) then
             allocate(curr%child_lt,STAT=j)
             if(j.ne.0) then
              if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): add: dictionary entry allocation failed!")')
              dict_search=DICT_ALLOC_FAILED; curr=>NULL(); return
             endif
             curr%child_lt%parent=>curr
             grow=+1; curr=>curr%child_lt
            elseif(i.eq.DICT_KEY_GT) then
             allocate(curr%child_gt,STAT=j)
             if(j.ne.0) then
              if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): add: dictionary entry allocation failed!")')
              dict_search=DICT_ALLOC_FAILED; curr=>NULL(); return
             endif
             curr%child_gt%parent=>curr
             grow=-1; curr=>curr%child_gt
            endif
            curr%child_lt=>NULL(); curr%child_gt=>NULL(); curr%balance_fac=0
            ierr=+1; call rebalance(curr%parent,grow,ierr)
            if(ierr.ne.0) then
             if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): add: rebalancing failed: ",i11)') ierr
             dict_search=DICT_CORRUPTED; curr=>NULL(); return
            endif
            allocate(curr%entry_key,source=item_key,STAT=j)
            if(j.ne.0) then
             if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): add: dictionary item key allocation failed!")')
             dict_search=DICT_ALLOC_FAILED; curr=>NULL(); return
            endif
            if(present(value_in)) then
             allocate(curr%entry_val,source=value_in,STAT=j)
             if(j.ne.0) then
              if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): add: dictionary item value allocation failed!")')
              dict_search=DICT_ALLOC_FAILED; curr=>NULL(); return
             endif
            else
             if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): add: absent value!")')
             dict_search=DICT_ABSENT_ARGUMENT; curr=>NULL(); return
            endif
            this%num_entries=this%num_entries+1_LONGINT
            this%curr_entry=>curr
            if(present(value_out)) value_out=>curr%entry_val
           else !empty dictionary
            allocate(this%root,STAT=j)
            if(j.ne.0) then
             if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): add: dictionary root allocation failed!")')
             dict_search=DICT_ALLOC_FAILED; curr=>NULL(); return
            endif
            this%root%parent=>NULL(); this%root%child_lt=>NULL(); this%root%child_gt=>NULL(); this%root%balance_fac=0
            allocate(this%root%entry_key,source=item_key,STAT=j)
            if(j.ne.0) then
             if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): add: dictionary root key allocation failed!")')
             dict_search=DICT_ALLOC_FAILED; curr=>NULL(); return
            endif
            if(present(value_in)) then
             allocate(this%root%entry_val,source=value_in,STAT=j)
             if(j.ne.0) then
              if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): add: dictionary root value allocation failed!")')
              dict_search=DICT_ALLOC_FAILED; curr=>NULL(); return
             endif
            else
             if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): add: absent value!")')
             dict_search=DICT_ABSENT_ARGUMENT; curr=>NULL(); return
            endif
            this%num_entries=1_LONGINT
            this%curr_entry=>this%root
            if(present(value_out)) value_out=>this%root%entry_val
           endif
          elseif(dict_search.eq.DICT_KEY_FOUND) then !return the found entry value (just in case)
           this%curr_entry=>curr
           if(present(value_out)) value_out=>curr%entry_val
          endif
         case default
          if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search): unknown search action request: ",i11)') action
          dict_search=DICT_UNKNOWN_REQUEST; curr=>NULL(); return
         end select
         curr=>NULL()
         return

         contains

          subroutine rebalance(cr,gr,ier)
           class(dict_entry_t), pointer:: cr !intermediate entry on the way back to root
           integer, intent(in):: gr !in(right/left subtree change): {-1;+1}
           integer, intent(inout):: ier !in(height decrease/increase): {-1;+1}; out: {0;{error_codes}}
           class(dict_entry_t), pointer:: cr_ptr
           integer:: jb,jc,jg,jd
!          if(debug) write(jo_dict,'("Entered rebalance:")') !debug
           jd=ier; jg=gr; cr_ptr=>cr
           do while(jg.ne.0)
!           if(debug) write(jo_dict,'(" Rebalance at level ",i3)') lev_p !debug
            cr_ptr%balance_fac=cr_ptr%balance_fac+jg*jd
            if(iabs(cr_ptr%balance_fac).ge.2) then !rotations needed
             if(cr_ptr%balance_fac.eq.-2) then
              jb=cr_ptr%child_gt%balance_fac
              if(jb.gt.0) then !{+1}
               jc=cr_ptr%child_gt%child_lt%balance_fac
               call rotate_double_left(cr_ptr)
               if(dict_search.ne.DICT_KEY_FOUND.and.dict_search.ne.DICT_KEY_NOT_FOUND) then; cr_ptr=>NULL(); ier=1; return; endif
               cr_ptr%parent%balance_fac=0; if(jd.gt.0) jg=0
               if(jc.gt.0) then
                cr_ptr%balance_fac=0; cr_ptr%parent%child_gt%balance_fac=-1
               elseif(jc.lt.0) then
                cr_ptr%balance_fac=1; cr_ptr%parent%child_gt%balance_fac=0
               else !jc=0
                cr_ptr%balance_fac=0; cr_ptr%parent%child_gt%balance_fac=0
               endif
              else !{-1;0}
               call rotate_simple_left(cr_ptr)
               if(dict_search.ne.DICT_KEY_FOUND.and.dict_search.ne.DICT_KEY_NOT_FOUND) then; cr_ptr=>NULL(); ier=2; return; endif
               if(jb.eq.0) then
                cr_ptr%balance_fac=-1; cr_ptr%parent%balance_fac=1; if(jd.lt.0) jg=0
               else
                cr_ptr%balance_fac=0; cr_ptr%parent%balance_fac=0; if(jd.gt.0) jg=0
               endif
              endif
              cr_ptr=>cr_ptr%parent
             elseif(cr_ptr%balance_fac.eq.2) then
              jb=cr_ptr%child_lt%balance_fac
              if(jb.lt.0) then !{-1}
               jc=cr_ptr%child_lt%child_gt%balance_fac
               call rotate_double_right(cr_ptr)
               if(dict_search.ne.DICT_KEY_FOUND.and.dict_search.ne.DICT_KEY_NOT_FOUND) then; cr_ptr=>NULL(); ier=3; return; endif
               cr_ptr%parent%balance_fac=0; if(jd.gt.0) jg=0
               if(jc.lt.0) then
                cr_ptr%balance_fac=0; cr_ptr%parent%child_lt%balance_fac=1
               elseif(jc.gt.0) then
                cr_ptr%balance_fac=-1; cr_ptr%parent%child_lt%balance_fac=0
               else !jc=0
                cr_ptr%balance_fac=0; cr_ptr%parent%child_lt%balance_fac=0
               endif
              else !{0;+1}
               call rotate_simple_right(cr_ptr)
               if(dict_search.ne.DICT_KEY_FOUND.and.dict_search.ne.DICT_KEY_NOT_FOUND) then; cr_ptr=>NULL(); ier=4; return; endif
               if(jb.eq.0) then
                cr_ptr%balance_fac=+1; cr_ptr%parent%balance_fac=-1; if(jd.lt.0) jg=0
               else
                cr_ptr%balance_fac=0; cr_ptr%parent%balance_fac=0; if(jd.gt.0) jg=0
               endif
              endif
              cr_ptr=>cr_ptr%parent
             else
              if(verbose)&
               &write(jo_dict,'("#ERROR(dictionary::dict_search): rebalance: invalid balance factor: ",i11)') cr_ptr%balance_fac
              cr_ptr=>NULL(); ier=5; return
             endif
            else !node balance factor changed to {-1;0;+1}
             if(jd.gt.0) then
              if(cr_ptr%balance_fac.eq.0) jg=0
             elseif(jd.lt.0) then
              if(cr_ptr%balance_fac.ne.0) jg=0
             endif
            endif
            if(associated(cr_ptr%parent)) then
             jg=iabs(jg); if(associated(cr_ptr%parent%child_gt,cr_ptr)) jg=-jg
             cr_ptr=>cr_ptr%parent; lev_p=lev_p-1
            else
             this%root=>cr_ptr
             exit
            endif
           enddo
           ier=0; cr_ptr=>NULL()
!          if(debug) write(jo_dict,'("Exited rebalance.")') !debug
           return
          end subroutine rebalance

          subroutine rotate_double_left(cr)
           class(dict_entry_t), pointer:: cr
           call rotate_simple_right(cr%child_gt)
           if(dict_search.ne.DICT_KEY_FOUND.and.dict_search.ne.DICT_KEY_NOT_FOUND) return
           call rotate_simple_left(cr)
           return
          end subroutine rotate_double_left

          subroutine rotate_double_right(cr)
           class(dict_entry_t), pointer:: cr
           call rotate_simple_left(cr%child_lt)
           if(dict_search.ne.DICT_KEY_FOUND.and.dict_search.ne.DICT_KEY_NOT_FOUND) return
           call rotate_simple_right(cr)
           return
          end subroutine rotate_double_right

          subroutine rotate_simple_left(cr)
           class(dict_entry_t), pointer:: cr
           class(dict_entry_t), pointer:: jp,jq,js
!          if(debug) write(jo_dict,'("  Rotating left")') !debug
           jp=>cr; jq=>cr%child_gt; js=>jq%child_lt !js may be null
           jq%child_lt=>jp
           if(associated(jp%parent)) then
            if(associated(jp%parent%child_lt,jp)) then
             jp%parent%child_lt=>jq
            elseif(associated(jp%parent%child_gt,jp)) then
             jp%parent%child_gt=>jq
            else
             if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search:rotate_simple_left): broken parental link!")')
             dict_search=DICT_CORRUPTED; return
            endif
            jq%parent=>jp%parent
           else
            jq%parent=>NULL()
           endif
           jp%parent=>jq
           jp%child_gt=>js; if(associated(js)) js%parent=>jp
           return
          end subroutine rotate_simple_left

          subroutine rotate_simple_right(cr)
           class(dict_entry_t), pointer:: cr
           class(dict_entry_t), pointer:: jp,jq,js
!          if(debug) write(jo_dict,'("  Rotating right")') !debug
           jp=>cr; jq=>cr%child_lt; js=>jq%child_gt !js may be null
           jq%child_gt=>jp
           if(associated(jp%parent)) then
            if(associated(jp%parent%child_lt,jp)) then
             jp%parent%child_lt=>jq
            elseif(associated(jp%parent%child_gt,jp)) then
             jp%parent%child_gt=>jq
            else
             if(verbose) write(jo_dict,'("#ERROR(dictionary::dict_search:rotate_simple_right): broken parental link!")')
             dict_search=DICT_CORRUPTED; return
            endif
            jq%parent=>jp%parent
           else
            jq%parent=>NULL()
           endif
           jp%parent=>jq
           jp%child_lt=>js; if(associated(js)) js%parent=>jp
           return
          end subroutine rotate_simple_right

        end function dict_search
!---------------------------------
        subroutine print_key(item) !debug procedure
         implicit none
         class(dict_entry_t), pointer:: item
         integer pr,cl,cr
         logical apr,acl,acr
         pr=0; cl=0; cr=0; apr=.false.; acl=.false.; acr=.false.
         select type (ke=>item%entry_key)
         type is (integer)
          if(associated(item%parent)) then
           apr=.true.
           select type (p=>item%parent%entry_key)
           type is (integer)
            pr=p
           end select
          endif
          if(associated(item%child_lt)) then
           acl=.true.
           select type (l=>item%child_lt%entry_key)
           type is (integer)
            cl=l
           end select
          endif
          if(associated(item%child_gt)) then
           acr=.true.
           select type (r=>item%child_gt%entry_key)
           type is (integer)
            cr=r
           end select
          endif
          write(jo_dict,'("Key ",i11,": Parent(",l1,")",i11,"; Left(",l1,")",i11,"; Right(",l1,")",i11)') ke,apr,pr,acl,cl,acr,cr
         class default
          if(verbose) write(jo_dict,'("#ERROR(dictionary::print_key): non-integer key:")')
          if(associated(item%parent)) then
           apr=.true.
           select type (p=>item%parent%entry_key)
           type is (integer)
            pr=p
           end select
          endif
          if(associated(item%child_lt)) then
           acl=.true.
           select type (l=>item%child_lt%entry_key)
           type is (integer)
            cl=l
           end select
          endif
          if(associated(item%child_gt)) then
           acr=.true.
           select type (r=>item%child_gt%entry_key)
           type is (integer)
            cr=r
           end select
          endif
          write(jo_dict,'("Key ?: Parent(",l1,")",i11,"; Left(",l1,")",i11,"; Right(",l1,")",i11)') apr,pr,acl,cl,acr,cr
         end select
         return
        end subroutine print_key

       end module dictionary
!-----------------------------------------------
!TESTING:
!-----------------------------------------------
       module dictionary_test
        use dictionary
        use timers, only: thread_wtime
        implicit none
        private
!PARAMETERS:
        integer, parameter, private:: KEY_LEN=16
        integer, parameter, private:: MAX_IND_VAL=7
!TYPES:
 !Key:
        type, private:: key_t
         integer:: rank=KEY_LEN
         integer:: dims(1:KEY_LEN)
        end type key_t
 !Values:
        type, private:: val_t
         real(8):: my_array(1:KEY_LEN)
         type(key_t):: key_stored
        end type val_t
!VISIBILITY:
        public dil_test_dictionary
        private cmp_key_test

       contains
!-------------------------------------------------
        function cmp_key_test(up1,up2) result(cmp)
         implicit none
         integer:: cmp
         class(*), intent(in):: up1,up2
         integer:: i

         cmp=DICT_KEY_ERR
         select type(up1)
         class is(key_t)
          select type(up2)
          class is(key_t)
           if(up1%rank.lt.up2%rank) then
            cmp=DICT_KEY_LT
           elseif(up1%rank.gt.up2%rank) then
            cmp=DICT_KEY_GT
           else
            cmp=DICT_KEY_EQ
            do i=1,up1%rank
             if(up1%dims(i).lt.up2%dims(i)) then
              cmp=DICT_KEY_LT; exit
             elseif(up1%dims(i).gt.up2%dims(i)) then
              cmp=DICT_KEY_GT; exit
             endif
            enddo
           endif
          end select
         end select
         return
        end function cmp_key_test
!--------------------------------------------------------------
        function dil_test_dictionary(perf,dev_out) result(ierr)
         implicit none
         integer:: ierr                          !out: error code (0:success)
         real(8), intent(out):: perf             !out: performance index
         integer, intent(in), optional:: dev_out !in: default output device
!-----------------------------------------------
         integer, parameter:: MAX_ACTIONS=1000000
         integer:: jo,i,j,fnd,nfnd
         type(key_t), target:: key
         type(val_t):: val
         type(dict_t):: my_dict
         class(*), pointer:: uptr
         real(8):: tm

         ierr=0; perf=0d0; uptr=>NULL()
         if(present(dev_out)) then; jo=dev_out; else; jo=6; endif
         fnd=0; nfnd=0; key%rank=KEY_LEN
         tm=thread_wtime()
         do i=1,MAX_ACTIONS
          call get_rnd_key(key) !random key
          val%my_array(1:KEY_LEN)=13.777d0; val%key_stored=key !value
          j=my_dict%search(DICT_ADD_IF_NOT_FOUND,cmp_key_test,key,val,uptr)
          if(j.eq.DICT_KEY_FOUND) then
           fnd=fnd+1
           if(.not.associated(uptr)) then; call test_quit(1); return; endif
          elseif(j.eq.DICT_KEY_NOT_FOUND) then
           nfnd=nfnd+1
          else
           call test_quit(2)
          endif
         enddo
         tm=thread_wtime(tm)
         perf=dble(MAX_ACTIONS)/tm
         return

         contains

          subroutine get_rnd_key(akey)
           type(key_t), intent(inout):: akey
           integer:: jj
           real(8):: jv(1:KEY_LEN)
           call random_number(jv(1:KEY_LEN))
           do jj=1,akey%rank
            akey%dims(jj)=nint(jv(jj)*dble(MAX_IND_VAL))
           enddo
           return
          end subroutine get_rnd_key

          subroutine test_quit(jerr)
           integer, intent(in):: jerr
           integer:: jj
           ierr=jerr
           if(ierr.ne.0) then
            write(jo,'("#ERROR(dictionary::dil_test_dictionary): Test failed: Error code ",i13)') ierr
            write(jo,'("Please contact the developer at QUANT4ME@GMAIL.COM")')
           endif
           jj=my_dict%destroy()
           return
          end subroutine test_quit

        end function dil_test_dictionary

       end module dictionary_test
