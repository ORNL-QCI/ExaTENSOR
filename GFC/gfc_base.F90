!Generic Fortran Containers (GFC): Base
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2016-11-19 (started 2016-02-17)

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

!FOR DEVELOPERS ONLY:
! # The base SCAN method as well as the concrete container constructors need to update
!   the reference count in the corresponding container element (in use). An element
!   of a container is considered IN USE if either it is a boundary element of some
!   container/iterator or it is associated with the current position of some
!   iterator or it is explicitly locked via the .in_use(errc,.TRUE.) method.
!   If the element of a container is IN USE, it cannot be deleted or relinked.
! # Inconsistency: When multiple iterators are associated with the same container,
!   a deletion of a container element via one of the iterators may result in an
!   undefined value for the boundary and/or current pointer in another iterator
!   in case the deleted element is one of the boundary or the current element
!   in that other iterator (also the other iterator's status may change to EMPTY).
!   Also, some derived containers (e.g., list_bi_t) allow splitting the container
!   into parts that can break other iterators associated with the same list.
!   So far the solution is to make user liable for such illegal use of GFC.
! # Quick counting does not work with composite containers and subcontainers
!   and probably it should not be used at all. Currently gfc_container_t::num_elems_()
!   will not return the total number of elements without quick counting. However, one
!   can always traverse the container in order to count the total number of elements.
!   Since .num_elems_() is currently marked with a trailing underscore, it shall not
!   be used by a user.

       module gfc_base
        use dil_basic
        use timers
#ifndef NO_OMP
        use omp_lib
#endif
        implicit none
        public
!PARAMETERS:
 !Error codes:
        integer(INTD), parameter:: GFC_SUCCESS=0            !success
        integer(INTD), parameter:: GFC_ERROR=-666           !unspecified error
        integer(INTD), parameter:: GFC_INVALID_ARGS=-1      !invalid arguments
        integer(INTD), parameter:: GFC_NULL_CONT=-2         !uninitialized container
        integer(INTD), parameter:: GFC_EMPTY_CONT=-3        !empty container
        integer(INTD), parameter:: GFC_MEM_ALLOC_FAILED=-4  !memory allocation failed
        integer(INTD), parameter:: GFC_MEM_FREE_FAILED=-5   !memory deallocation failed
        integer(INTD), parameter:: GFC_CORRUPTED_CONT=-6    !corrupted container
        integer(INTD), parameter:: GFC_ELEM_EMPTY=-7        !element of the container is empty
        integer(INTD), parameter:: GFC_ELEM_NOT_EMPTY=-8    !element of the container is not empty
        integer(INTD), parameter:: GFC_ACTION_FAILED=-9     !user-defined action failed on an element
        integer(INTD), parameter:: GFC_METHOD_UNDEFINED=-10 !undefined method called on an object
        integer(INTD), parameter:: GFC_IN_USE=-11           !object is in use by others
 !Predicates (GFC_ERROR applies here as well):
        integer(INTD), parameter:: GFC_TRUE=1  !TRUE value
        integer(INTD), parameter:: GFC_FALSE=0 !FALSE value
 !Search results (GFC_ERROR applies here as well):
        integer(INTD), parameter:: GFC_FOUND=GFC_SUCCESS !element found
        integer(INTD), parameter:: GFC_NOT_FOUND=1       !element not found
 !Comparison/relation (GFC_ERROR applies here as well):
        integer(INTD), parameter:: GFC_CMP_EQ=0      !equivalent objects
        integer(INTD), parameter:: GFC_CMP_LT=-1     !object1 < object2
        integer(INTD), parameter:: GFC_CMP_GT=+1     !object1 > object2
        integer(INTD), parameter:: GFC_CMP_CHILD=-2  !object1 is a child of object2
        integer(INTD), parameter:: GFC_CMP_PARENT=+2 !object1 is a parent of object2
        integer(INTD), parameter:: GFC_CMP_NA=-6     !objects are not comparable
 !GFC iterator status:
        integer(INTD), parameter:: GFC_IT_NULL=1000   !uninitialized iterator
        integer(INTD), parameter:: GFC_IT_EMPTY=1001  !empty initialized iterator
        integer(INTD), parameter:: GFC_IT_ACTIVE=1002 !active (non-empty) iterator
        integer(INTD), parameter:: GFC_IT_DONE=1003   !pass the end of the container (done), needs to be reset to continue
        integer(INTD), parameter:: GFC_NO_MOVE=1004   !no move possible (for custom moves)
!TYPES:
 !Element of a container:
        type, public:: gfc_cont_elem_t
         class(*), pointer, private:: value_p=>NULL() !element value (data): either associated (by reference) or allocated (by value)
         integer(INTD), private:: alloc=GFC_FALSE     !GFC_FALSE: value is stored by reference or null; GFC_TRUE: value is stored by value
#ifndef NO_OMP
         integer(INTD), private:: ref_count=0 !reference count (incremented when the container boundary or an iterator point to the element)
         integer(INTD), private:: lock=0      !update lock (for concurrent updates): Can be set by user to ensure an exclusive access
#endif
         contains
          procedure, public:: construct_base=>ContElemConstruct !constructs a new container element, either by reference or by value
          procedure, public:: destruct=>ContElemDestruct        !destructs an existing container element (releases memory occupied by its value)
          procedure, public:: get_value=>ContElemGetValue       !returns a pointer to the element value (unlimited polymorphic)
          procedure, public:: is_empty=>ContElemIsEmpty         !returns TRUE if the element of the container is empty, FALSE otherwise
          procedure, public:: stored_by_value=>ContElemStoredByValue !returns TRUE if the element value is stored by value, FALSE otherwise (stored by reference)
          procedure, public:: predicate=>ContElemPredicate      !returns the value of a user-given predicate applied to the element
          procedure, public:: ContElemAction                    !acts on the element with a user-defined action function
          procedure, public:: ContElemFunctor                   !acts on the element with a user-defined functor
          generic, public:: action=>ContElemAction,ContElemFunctor !generic overload (acts on the container element value)
          procedure, public:: compare=>ContElemCompare          !compares the value of the element with the value of another element
          procedure, public:: print_it=>ContElemPrintIt         !prints the value of the element with a user-defined print function
          procedure, public:: in_use=>ContElemInUse             !returns TRUE if the element of the container is currently in use, hence cannot be deleted
          procedure, public:: release_lock_=>ContElemReleaseLock !PRIVATE: releases the lock on the container element
          procedure, public:: incr_ref_=>ContElemIncrRef         !PRIVATE: increments the reference count for the container element
          procedure, public:: decr_ref_=>ContElemDecrRef         !PRIVATE: decrements the reference count for the container element
        end type gfc_cont_elem_t
 !Base container:
        type, abstract, public:: gfc_container_t
         integer(INTL), private:: volume=0_INTL !volume of the container (total number of elements when quick counting is on), -1 means quick counting is off
         contains
          procedure, non_overridable, public:: num_elems_=>ContNumElems !PRIVATE: returns the total number of elements in the container
          procedure, non_overridable, public:: update_num_elems_=>ContUpdateNumElems !PRIVATE: updates the number of elements
          procedure, non_overridable, public:: quick_counting_off_=>ContQuickCountingOff !PRIVATE: turns off quick element counting
        end type gfc_container_t
 !Base iterator:
        type, abstract, public:: gfc_iter_t
         integer(INTD), private:: state=GFC_IT_NULL !current state of the iterator (see possible iterator states above)
         integer(INTL), private:: tot_count=0_INTL  !total number of iterations after the last reset
         integer(INTL), private:: pred_count=0_INTL !number of iterations with a TRUE predicate after the last reset
         contains
          procedure, non_overridable, public:: get_status=>IterGetStatus  !returns the status of the iterator
          procedure, non_overridable, public:: set_status_=>IterSetStatus !PRIVATE: sets the status of the iterator
          procedure, public:: reset_count=>IterResetCount                 !resets all iteration counters to zero
          procedure, public:: total_count=>IterTotalCount                 !returns the total iteration count since the last reset
          procedure, public:: predicated_count=>IterPredicatedCount       !returns the TRUE predicated iteration count since the last reset
          procedure, public:: scan=>IterScan                              !traverses the container with an optional action
          procedure(gfc_it_init_i), deferred, public:: init       !initializes the iterator (associates it with a container and positions it on the root)
          procedure(gfc_it_reset_i), deferred, public:: reset     !resets the iterator to the beginning (root element)
          procedure(gfc_it_reset_i), deferred, public:: release   !dissociates the iterator from its container
          procedure(gfc_it_pointee_i), deferred, public:: pointee !returns the element currently pointed to
          procedure(gfc_it_next_i), deferred, public:: next       !proceeds to the next element of the container
          procedure(gfc_it_next_i), deferred, public:: previous   !proceeds to the previous element of the container
        end type gfc_iter_t
 !Base functor:
        type, abstract, public:: gfc_func_t
         contains
          procedure(gfc_func_act_i), deferred, public:: act !performs an action on an unlimited polymorphic object
        end type gfc_func_t
!ABSTRACT INTERFACES:
        abstract interface
 !Generics:
  !GFC generic predicate:
         function gfc_predicate_i(obj) result(pred)
          import:: INTD
          integer(INTD):: pred       !out: predicate value: {GFC_TRUE, GFC_FALSE, GFC_ERROR, other integers}
          class(*), intent(in):: obj !in: arbitrary object
         end function gfc_predicate_i
  !GFC generic relation:
         function gfc_cmp_i(obj1,obj2) result(cmp)
          import:: INTD
          integer(INTD):: cmp         !out: result of the comparison (relation): See CMP parameters above
          class(*), intent(in):: obj1 !in: object 1
          class(*), intent(in):: obj2 !in: object 2
         end function gfc_cmp_i
  !GFC generic copy constructor:
         function gfc_copy_i(obj,ierr) result(clone)
          import:: INTD
          class(*), pointer:: clone                   !out: clone
          class(*), intent(in):: obj                  !in: original object
          integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         end function gfc_copy_i
  !GFC generic destructor (destructs the inner state of an object, but does not DEALLOCATE the object):
         function gfc_destruct_i(obj) result(ierr)
          import:: INTD
          integer(INTD):: ierr          !out: error code (0:success)
          class(*), intent(inout):: obj !inout: object for destruction
         end function gfc_destruct_i
  !GFC generic action:
         function gfc_action_i(obj) result(ierr)
          import:: INTD
          integer(INTD):: ierr          !out: error code (0:success)
          class(*), intent(inout):: obj !inout: object to be acted on
         end function gfc_action_i
  !GFC generic printing:
         function gfc_print_i(obj,dev_id) result(ierr)
          import:: INTD
          integer(INTD):: ierr                         !out: error code (0:success)
          class(*), intent(in):: obj                   !in: arbitrary object
          integer(INTD), intent(in), optional:: dev_id !in: output device id (defaults to screen: 6)
         end function gfc_print_i
 !Deferred:
  !Deferred: GFC iterator: .init:
         function gfc_it_init_i(this,cont) result(ierr)
          import:: gfc_iter_t,gfc_container_t,INTD
          integer(INTD):: ierr                              !out: error code
          class(gfc_iter_t), intent(inout):: this           !inout: GFC iterator
          class(gfc_container_t), target, intent(in):: cont !in: GFC container
         end function gfc_it_init_i
  !Deferred: GFC iterator: .reset .release:
         function gfc_it_reset_i(this) result(ierr)
          import:: gfc_iter_t,INTD
          integer(INTD):: ierr                    !out: error code
          class(gfc_iter_t), intent(inout):: this !inout: GFC iterator
         end function gfc_it_reset_i
  !Deferred: GFC iterator: .pointee:
         function gfc_it_pointee_i(this,ierr) result(pntee)
          import:: gfc_iter_t,gfc_cont_elem_t,INTD
          class(gfc_cont_elem_t), pointer:: pntee     !out: container element currently pointed to by the iterator
          class(gfc_iter_t), intent(in):: this        !in: GFC iterator
          integer(INTD), intent(out), optional:: ierr !out: error code
         end function gfc_it_pointee_i
  !Deferred: GFC iterator: .next .previous:
         function gfc_it_next_i(this,elem_p) result(ierr)
          import:: gfc_iter_t,gfc_cont_elem_t,INTD
          integer(INTD):: ierr                                            !out: error code
          class(gfc_iter_t), intent(inout):: this                         !inout: GFC iterator
          class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element
         end function gfc_it_next_i
  !Deferred: GFC functor action: .act:
         function gfc_func_act_i(this,obj) result(ierr)
          import:: gfc_func_t,INTD
          integer(INTD):: ierr                    !out: error code
          class(gfc_func_t), intent(inout):: this !inout: GFC functor (may change its state)
          class(*), intent(inout):: obj           !inout: arbitrary object the functor is acting upon
         end function gfc_func_act_i
        end interface
!VISIBILITY:
        private ContElemConstruct
        private ContElemDestruct
        private ContElemGetValue
        private ContElemIsEmpty
        private ContElemStoredByValue
        private ContElemPredicate
        private ContElemAction
        private ContElemFunctor
        private ContElemCompare
        private ContElemPrintIt
        private ContNumElems
        private ContUpdateNumElems
        private ContQuickCountingOff
        private IterGetStatus
        private IterSetStatus
        private IterResetCount
        private IterTotalCount
        private IterPredicatedCount
        private IterScan

       contains
!IMPLEMENTATION:
!--------------------------------------------------------------------------------
#ifdef NO_GNU
        subroutine ContElemConstruct(this,obj,ierr,assoc_only,copy_ctor_f,locked) !`GCC/5.3.0 has a bug with this line
#else
        subroutine ContElemConstruct(this,obj,ierr,assoc_only,locked)
#endif
!Constructs the base part of the element of a container (sets its value).
!The constructor is allowed to construct the value of a boundary element
!as well as an element pointed to by a container iterator.
         implicit none
         class(gfc_cont_elem_t), intent(inout):: this  !inout: element of a container
         class(*), target, intent(in):: obj            !in: value to be stored in this element
         integer(INTD), intent(out), optional:: ierr   !out: error code (0:success)
         logical, intent(in), optional:: assoc_only    !in: if TRUE, <obj> will be stored by reference, otherwise by value (default)
#ifdef NO_GNU
         procedure(gfc_copy_i), optional:: copy_ctor_f !in: user-defined generic copy constructor
#endif
         logical, intent(in), optional:: locked        !in: if TRUE, the container element will be assumed already locked (defaults to FALSE)
         integer(INTD):: errc
         integer:: errcode
         logical:: assoc,lckd,lck

         errc=GFC_SUCCESS
         if(present(assoc_only)) then; assoc=assoc_only; else; assoc=.FALSE.; endif
         if(present(locked)) then; lckd=locked; else; lckd=.FALSE.; endif; lck=lckd
         if(this%is_empty()) then
          if(.not.lck) lck=(this%in_use(errc,set_lock=.TRUE.,report_refs=.FALSE.).eq.GFC_FALSE)
          if(lck) then
           if(errc.eq.GFC_SUCCESS) then
            if(assoc) then
             this%value_p=>obj
             this%alloc=GFC_FALSE
            else
#ifdef NO_GNU
             if(present(copy_ctor_f)) then
              this%value_p=>copy_ctor_f(obj,errc)
             else
#endif
              allocate(this%value_p,SOURCE=obj,STAT=errcode)
              if(errcode.ne.0) errc=GFC_MEM_ALLOC_FAILED
#ifdef NO_GNU
             endif
#endif
             if(errc.eq.GFC_SUCCESS) then
              this%alloc=GFC_TRUE
             else
              this%value_p=>NULL()
             endif
            endif
           endif
           if(.not.lckd) then
            if(errc.eq.GFC_SUCCESS) then
             call this%release_lock_(errc)
            else
             call this%release_lock_()
            endif
           endif
          else
           if(errc.eq.GFC_SUCCESS) errc=GFC_IN_USE
          endif
         else
          errc=GFC_ELEM_NOT_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContElemConstruct
!-----------------------------------------------------------
        subroutine ContElemDestruct(this,ierr,dtor_f,locked)
!Destructs the value of an element of a container (not the element itself).
!On return, the element will be empty, unless an error occurs. An optional
!explicit destructor, if provided, will only be invoked for allocated values.
!Alternatively, the value may also have the final procedure defined.
!The container element being destructed is allowed to be a boundary
!element or an element pointed to by a container iterator.
         implicit none
         class(gfc_cont_elem_t), intent(inout):: this !inout: element of a container
         integer(INTD), intent(out), optional:: ierr  !out: error code (0:success)
         procedure(gfc_destruct_i), optional:: dtor_f !in: explicit destructor for the value of this element (for allocated only)
         logical, intent(in), optional:: locked       !in: if TRUE, the container element will be assumed already locked (defaults to FALSE)
         integer(INTD):: errc
         integer:: errcode
         logical:: lckd,lck

         errc=GFC_SUCCESS
         if(present(locked)) then; lckd=locked; else; lckd=.FALSE.; endif; lck=lckd
         if(.not.this%is_empty()) then
          if(.not.lck) lck=(this%in_use(errc,set_lock=.TRUE.,report_refs=.FALSE.).eq.GFC_FALSE)
          if(lck) then
           if(errc.eq.GFC_SUCCESS) then
            if(this%alloc.eq.GFC_TRUE) then
             if(present(dtor_f)) then
              errc=dtor_f(this%value_p)
              if(errc.ne.0) errc=GFC_MEM_FREE_FAILED
             endif
             deallocate(this%value_p,STAT=errcode)
             if(errcode.ne.0.and.errc.eq.GFC_SUCCESS) errc=GFC_MEM_FREE_FAILED
            endif
            this%value_p=>NULL()
            this%alloc=GFC_FALSE
           endif
           if(.not.lckd) then
            if(errc.eq.GFC_SUCCESS) then
             call this%release_lock_(errc)
            else
             call this%release_lock_()
            endif
           endif
          else
           if(errc.eq.GFC_SUCCESS) errc=GFC_IN_USE
          endif
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContElemDestruct
!---------------------------------------------------------
        function ContElemGetValue(this,ierr) result(val_p)
!Returns an unlimited polymorphic pointer to the element value.
         implicit none
         class(*), pointer:: val_p                   !out: pointer to the element value
         class(gfc_cont_elem_t), intent(in):: this   !in: element of a container
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc

         errc=GFC_SUCCESS; val_p=>NULL()
         if(.not.this%is_empty()) then
          val_p=>this%value_p
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end function ContElemGetValue
!--------------------------------------------------
        function ContElemIsEmpty(this) result(empt)
!Returns TRUE if the element of a container is empty, FALSE otherwise.
         implicit none
         logical:: empt                            !out: empty (true) or not (false)
         class(gfc_cont_elem_t), intent(in):: this !in: element of a container

         empt=.not.associated(this%value_p)
         return
        end function ContElemIsEmpty
!------------------------------------------------------------
        function ContElemStoredByValue(this,ierr) result(res)
!Returns TRUE if the element value is stored by value, FALSE otherwise (stored by reference).
         implicit none
         logical:: res                               !out: answer
         class(gfc_cont_elem_t), intent(in):: this   !in: element of a container
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS; res=.FALSE.
         if(.not.this%is_empty()) then
          res=(this%alloc.eq.GFC_TRUE)
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end function ContElemStoredByValue
!--------------------------------------------------------------------
        function ContElemPredicate(this,predicat_f,ierr) result(pred)
!Evaluates a user-defined predicate on the value of a given container element.
         implicit none
         integer(INTD):: pred                        !out: evaluated predicate value {GFC_TRUE,GFC_FALSE,GFC_ERROR}
         class(gfc_cont_elem_t), intent(in):: this   !in: element of a container
         procedure(gfc_predicate_i):: predicat_f     !in: user-defined predicate
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc

         errc=GFC_SUCCESS; pred=GFC_ERROR
         if(.not.this%is_empty()) then
          pred=predicat_f(this%value_p); if(pred.eq.GFC_ERROR) errc=GFC_ERROR
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end function ContElemPredicate
!----------------------------------------------------
        subroutine ContElemAction(this,action_f,ierr)
!Performs a user-defined action on the value of a given element of a container
!via a function. It is the user responsibility to avoid race conditions
!when updating container element values!
         implicit none
         class(gfc_cont_elem_t), intent(inout):: this !inout: element of a container
         procedure(gfc_action_i):: action_f           !in: user-defined action
         integer(INTD), intent(out), optional:: ierr  !out: error code (0:success)
         integer(INTD):: errc

         errc=GFC_SUCCESS
         if(.not.this%is_empty()) then
          errc=action_f(this%value_p); if(errc.ne.0) errc=GFC_ACTION_FAILED
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContElemAction
!----------------------------------------------------
        subroutine ContElemFunctor(this,functor,ierr)
!Performs a user-defined action on the value of a given container element
!via a functor. It is the user responsibility to avoid race conditions
!when updating container element values!
         implicit none
         class(gfc_cont_elem_t), intent(inout):: this !inout: element of a container
         class(gfc_func_t), intent(inout):: functor   !inout: functor (may change its state)
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
         if(.not.this%is_empty()) then
          errc=functor%act(this); if(errc.ne.0) errc=GFC_ACTION_FAILED
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContElemFunctor
!------------------------------------------------------------------------
        function ContElemCompare(this,other_value,cmp_f,ierr) result(cmp)
!Compares the value of a given container element (on the left) with another
!value (on the right) using an appropriate comparison function.
         implicit none
         integer(INTD):: cmp                         !out: result of the comparison (see GFC_CMP_XXX above)
         class(gfc_cont_elem_t), intent(in):: this   !in: element of a container whose value is being compared
         class(*), intent(in):: other_value          !in: the other value to be compared with
         procedure(gfc_cmp_i):: cmp_f                !in: user-defined comparison function
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc

         errc=GFC_SUCCESS; cmp=GFC_CMP_NA
         if(.not.this%is_empty()) then
          cmp=cmp_f(this%value_p,other_value); if(cmp.eq.GFC_ERROR) errc=GFC_ERROR
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end function ContElemCompare
!-----------------------------------------------------------
        subroutine ContElemPrintIt(this,print_f,ierr,dev_id)
!Prints the given element of a container using a user-defined print function.
         implicit none
         class(gfc_cont_elem_t), intent(in):: this    !in: element of a container
         procedure(gfc_print_i):: print_f             !in: user-defined printing function
         integer(INTD), intent(out), optional:: ierr  !out: error code (0:success)
         integer(INTD), intent(in), optional:: dev_id !in: output device (default to screen, 6)
         integer(INTD):: dev
         integer(INTD):: errc

         errc=GFC_SUCCESS
         if(.not.this%is_empty()) then
          if(present(dev_id)) then; dev=dev_id; else; dev=6; endif !defaults to screen
          errc=print_f(this%value_p,dev); if(errc.ne.0) errc=GFC_ACTION_FAILED
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContElemPrintIt
!----------------------------------------------------------------------------
        function ContElemInUse(this,ierr,set_lock,report_refs) result(in_use)
!Returns GFC_TRUE if the element of the container is currently in use,
!hence cannot be deleted. The element of the container is considered
!IN USE if it is currently locked. Additionally, if <report_refs> = TRUE,
!an element of a container is considered IN USE if it is a boundary element
!of some container or it is associated with the current iterator position in
!some iterator. If the element of the container is not in use,
!a status GFC_FALSE is returned. In the latter case, setting <set_lock> to TRUE
!will set the lock, thus providing an exclusive access to the container element.
         implicit none
         integer(INTD):: in_use                       !out: GFC_TRUE, GFC_FALSE, or some other (error) status
         class(gfc_cont_elem_t), intent(inout):: this !inout: element of the container
         integer(INTD), intent(out), optional:: ierr  !out: error code
         logical, intent(in), optional:: set_lock     !in: if TRUE and the element is not in use, the lock will be set (default FALSE)
         logical, intent(in), optional:: report_refs  !in: if TRUE, being referred to by an iterator is considered IN USE as well (default FALSE)
         integer(INTD):: errc
         logical:: refs

         errc=GFC_SUCCESS
#ifndef NO_OMP
         in_use=GFC_TRUE
         if(present(report_refs)) then; refs=report_refs; else; refs=.FALSE.; endif
!$OMP CRITICAL (GFC_LOCK)
         if(this%lock.eq.0) then
          if((.not.refs).or.this%ref_count.eq.0) then
           in_use=GFC_FALSE
           if(present(set_lock)) then
            if(set_lock) this%lock=1
           endif
          endif
         endif
!$OMP END CRITICAL (GFC_LOCK)
#else
         in_use=GFC_FALSE
#endif
         if(present(ierr)) ierr=errc
         return
        end function ContElemInUse
!------------------------------------------------
        subroutine ContElemReleaseLock(this,ierr) !INTERNAL USE ONLY
!Releases the lock in a container element.
         implicit none
         class(gfc_cont_elem_t), intent(inout):: this !inout: element of the container
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
#ifndef NO_OMP
!$OMP CRITICAL (GFC_LOCK)
         this%lock=0
!$OMP END CRITICAL (GFC_LOCK)
#endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContElemReleaseLock
!--------------------------------------------
        subroutine ContElemIncrRef(this,ierr) !INTERNAL USE ONLY!
!Increments the reference count.
         implicit none
         class(gfc_cont_elem_t), intent(inout):: this !inout: container element
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
#ifndef NO_OMP
!$OMP CRITICAL (GFC_LOCK)
         if(this%ref_count.ge.0) then
          this%ref_count=this%ref_count+1
         else
          errc=GFC_CORRUPTED_CONT
         endif
!$OMP END CRITICAL (GFC_LOCK)
#endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContElemIncrRef
!--------------------------------------------
        subroutine ContElemDecrRef(this,ierr) !INTERNAL USE ONLY!
!Decrements the reference count.
         implicit none
         class(gfc_cont_elem_t), intent(inout):: this !inout: container element
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
#ifndef NO_OMP
!$OMP CRITICAL (GFC_LOCK)
         if(this%ref_count.gt.0) then
          this%ref_count=this%ref_count-1
         else
          errc=GFC_CORRUPTED_CONT
         endif
!$OMP END CRITICAL (GFC_LOCK)
#endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContElemDecrRef
!------------------------------------------------------
        function ContNumElems(this,ierr) result(nelems) !INTERNAL USE ONLY!
!Returns the total number of elements stored in the container,
!unless quick counting has been deactivated (returns -1 then).
         implicit none
         integer(INTL):: nelems                      !out: total number of elements in the container or -1
         class(gfc_container_t), intent(in):: this   !in: container
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc

         errc=GFC_SUCCESS
         if(this%volume.ge.0) then
          nelems=this%volume
         else
          nelems=-1; errc=GFC_METHOD_UNDEFINED
          !`Count somehow else (quick counting deactivated)
         endif
         if(present(ierr)) ierr=errc
         return
        end function ContNumElems
!----------------------------------------------------------------------
        function ContUpdateNumElems(this,new_elems,ierr) result(nelems) !INTERNAL USE ONLY!
!Updates the total number of elements in the container,
!unless quick counting has been deactivated (returns -1 then).
         implicit none
         integer(INTL):: nelems                       !out: total number of elements after the update
         class(gfc_container_t), intent(inout):: this !inout: container
         integer(INTL), intent(in):: new_elems        !in: number of elements added (+) or deleted (-)
         integer(INTD), intent(out), optional:: ierr  !out: error code (0:success)
         integer(INTD):: errc

         errc=GFC_SUCCESS
         if(this%volume.ge.0) then
          nelems=this%volume+new_elems
          if(nelems.ge.0) then
           this%volume=nelems
          else
           errc=GFC_INVALID_ARGS
          endif
         else
          nelems=-1; errc=GFC_METHOD_UNDEFINED !quick counting deactivated
         endif
         if(present(ierr)) ierr=errc
         return
        end function ContUpdateNumElems
!-------------------------------------------------
        subroutine ContQuickCountingOff(this,ierr) !INTERNAL USE ONLY!
!Turns off quick element counting.
         implicit none
         class(gfc_container_t), intent(inout):: this !inout: container
         integer(INTD), intent(out), optional:: ierr
         integer(INTD):: errc

         errc=GFC_SUCCESS
         this%volume=-1
         if(present(ierr)) ierr=errc
         return
        end subroutine ContQuickCountingOff
!----------------------------------------------------
        function IterGetStatus(this,ierr) result(sts)
!Returns the status of the iterator.
         implicit none
         integer(INTD):: sts                         !out: current status of the iterator
         class(gfc_iter_t), intent(in):: this        !in: iterator
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc

         errc=GFC_SUCCESS
         sts=this%state
         if(present(ierr)) ierr=errc
         return
        end function IterGetStatus
!----------------------------------------------------
        function IterSetStatus(this,sts) result(ierr) !INTERNAL USE ONLY!
!Sets the status of the iterator.
         implicit none
         integer(INTD):: ierr                    !out: error code (0:success)
         class(gfc_iter_t), intent(inout):: this !in: iterator
         integer(INTD), intent(in):: sts         !in: status

         ierr=GFC_SUCCESS
         select case(sts)
         case(GFC_IT_NULL,GFC_IT_EMPTY,GFC_IT_ACTIVE,GFC_IT_DONE)
          this%state=sts
         case default
          ierr=GFC_INVALID_ARGS
         end select
         return
        end function IterSetStatus
!-------------------------------------------
        subroutine IterResetCount(this,ierr)
!Resets all iteration counters.
         implicit none
         class(gfc_iter_t), intent(inout):: this     !inout: iterator
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc

         errc=GFC_SUCCESS
         this%tot_count=0; this%pred_count=0
         if(present(ierr)) ierr=errc
         return
        end subroutine IterResetCount
!-----------------------------------------------------
        function IterTotalCount(this,ierr) result(cnt)
!Returns the total iteration count since the last reset.
         implicit none
         integer(INTL):: cnt                         !out: total iteration count
         class(gfc_iter_t), intent(in):: this        !in: iterator
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc

         cnt=-1_INTL; errc=this%get_status()
         if(errc.ne.GFC_IT_NULL) then
          cnt=this%tot_count; errc=GFC_SUCCESS
         endif
         if(present(ierr)) ierr=errc
         return
        end function IterTotalCount
!----------------------------------------------------------
        function IterPredicatedCount(this,ierr) result(cnt)
!Returns the predicated iteration count since the last reset.
         implicit none
         integer(INTL):: cnt                         !out: predicated iteration count
         class(gfc_iter_t), intent(in):: this        !in: iterator
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc

         cnt=-1_INTL; errc=this%get_status()
         if(errc.ne.GFC_IT_NULL) then
          cnt=this%pred_count; errc=GFC_SUCCESS
         endif
         if(present(ierr)) ierr=errc
         return
        end function IterPredicatedCount
!-------------------------------------------------------------------------------------------------------------
        function IterScan(this,return_each,predicate_f,action_f,backward,skip_current,time_limit) result(ierr)
!Traverses the container via an associated iterator beginning from the current position of the iterator.
!Returns GFC_IT_DONE upon reaching the end of the container; returns GFC_IT_ACTIVE in intermediate returns.
!If the active scan is time limited, at least one move of the iterator will be done before returning.
         implicit none
         integer(INTD):: ierr !out: error code (GFC_IT_ACTIVE:intermediate return, GFC_IT_DONE:done, Other:empty or error)
         class(gfc_iter_t), intent(inout):: this !inout: iterator
         logical, intent(in), optional:: return_each !in: if TRUE, each successful match will be returned (defaults to FALSE)
         procedure(gfc_predicate_i), optional:: predicate_f !in: predicating function
         procedure(gfc_action_i), optional:: action_f !in: action function
         logical, intent(in), optional:: backward !in: if TRUE, the container will be traversed in the backward direction (defaults to FALSE)
         logical, intent(in), optional:: skip_current !in: if TRUE, the current element of the container will be skipped (forces a move) (defaults to FALSE)
         real(8), intent(in), optional:: time_limit !in: if specified, the active scan will be interrupted after this time limit (sec)
         logical:: ret,pred,act,bkw,moved,skip
         integer(INTD):: pred_val
         class(*), pointer:: elem_val
         class(gfc_cont_elem_t), pointer:: curr
         real(8):: tml,tms

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(present(return_each)) then; ret=return_each; else; ret=.false.; endif
          if(present(predicate_f)) then; pred=.true.; else; pred=.false.; endif
          if(present(action_f)) then; act=.true.; else; act=.false.; endif
          if(present(backward)) then; bkw=backward; else; bkw=.false.; endif
          if(present(skip_current)) then; skip=skip_current; else; skip=.false.; endif
          if(present(time_limit).and.(.not.ret)) then; tml=time_limit; else; tml=-1d0; endif
          if(tml.gt.0d0) tms=thread_wtime()
          ierr=GFC_SUCCESS; moved=.false.
          do while(ierr.eq.GFC_SUCCESS)
           curr=>this%pointee(ierr)
           if(ierr.eq.GFC_SUCCESS.and.associated(curr)) then
            elem_val=>curr%get_value(ierr)
            if(ierr.eq.GFC_SUCCESS.and.associated(elem_val)) then
             if(skip) then
              pred_val=GFC_FALSE
             else
              pred_val=GFC_TRUE; if(pred) pred_val=curr%predicate(predicate_f)
             endif
             if(pred_val.eq.GFC_TRUE) then
              this%pred_count=this%pred_count+1
              if(act) then
               call curr%action(action_f,ierr); if(ierr.ne.0) then; ierr=GFC_ACTION_FAILED; exit; endif
              endif
              if(ret) then; ierr=GFC_IT_ACTIVE; exit; endif !intermediate return on passive scans
             endif
             if(tml.gt.0d0) then
              if(thread_wtime(tms).gt.tml.and.moved) then; ierr=GFC_IT_ACTIVE; exit; endif !time limit exceeded (but at least one move done)
             endif
             if(bkw) then; ierr=this%previous(); else; ierr=this%next(); endif !move to the next/previous element
             moved=.true.; skip=.false.; this%tot_count=this%tot_count+1 !register a move of the iterator
            else
             ierr=GFC_CORRUPTED_CONT
            endif
           else
            ierr=GFC_CORRUPTED_CONT
           endif
          enddo
         endif
         return
        end function IterScan

       end module gfc_base
