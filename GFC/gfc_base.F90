!Generic Fortran Containers (GFC): Base
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2016-02-20 (started 2016-02-17)
!Copyright (C) 2016 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2016 Oak Ridge National Laboratory (UT-Battelle)
!LICENSE: GNU GPL v.2
!DESIGN:
! # A GFC container is a structured collection of objects of any class;
!   Objects of different classes can be stored in the same container.
! # Objects can be added to the container either by value or by reference;
!   The same container may have elements added both ways.
! # Each container has an associated iterator for scanning over its elements.
!   The structure of the container determines the scanning sequence, that is,
!   the way the elements of the container are traversed over.
!   There are four major ways of scanning over a container:
!    1) Unpredicated passive scan: The iterator returns each new encountered
!       element of the container.
!    2) Predicated passive scan: The iterator returns each new encountered
!       element of the container that satisfies a certain condition.
!    3) Unpredicated active scan: The iterator traverses the entire container
!       or its part and applies a user-defined action to each element.
!    4) Predicated active scan: The iterator traverses the entire container
!       or its part and applies a user-defined action to each element that
!       satisfies a certain condition.
!   Additionally, active scans allow for a time-limited execution, in which
!   the scan is interrupted after a given time interval.
!   Each specific class of containers has its own iterator class because
!   of different linkage between the elements of different containers.
!   All insertion, deletion, and search operations are done via iterators.
! # The container element deletion operation may require a user-supplied
!   destructor which will release all resources occupied by the object
!   stored in that element, unless the object has FINAL methods defined.
!NOTES:
! # This implementation of Generic Fortran containers heavily relies on
!   dynamic polymorphism, thus introducing certain memory overhead which
!   can be significant when storing small objects! Also dynamic type
!   inferrence may decrease the efficiency when storing and operating
!   on small objects.
       module gfc_base
        use dil_basic
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
        integer(INTD), parameter:: GFC_END_CONT=1           !end of container
 !Predicates (GFC_ERROR applies here as well):
        integer(INTD), parameter:: GFC_TRUE=1  !TRUE value
        integer(INTD), parameter:: GFC_FALSE=0 !FALSE value
 !Comparison/relation (GFC_ERROR applies here as well):
        integer(INTD), parameter:: GFC_CMP_EQ=0      !equivalent objects
        integer(INTD), parameter:: GFC_CMP_LT=-1     !object1 < object2
        integer(INTD), parameter:: GFC_CMP_GT=+1     !object1 > object2
        integer(INTD), parameter:: GFC_CMP_CHILD=-2  !object1 is a child of object2
        integer(INTD), parameter:: GFC_CMP_PARENT=+2 !object1 is a parent of object2
        integer(INTD), parameter:: GFC_CMP_NA=-6     !objects are not comparable
 !GFC iterator status:
        integer(INTD), parameter:: GFC_IT_NULL=1000    !uninitialized iterator
        integer(INTD), parameter:: GFC_IT_EMPTY=1001   !empty initialized iterator
        integer(INTD), parameter:: GFC_IT_ACTIVE=1002  !active (non-empty) iterator
        integer(INTD), parameter:: GFC_IT_DONE=1003    !pass the end of the container (done)
!TYPES:
 !Element of a container:
        type, public:: gfc_cont_elem_t
         class(*), pointer, private:: value_p=>NULL() !element value (data): either associated (by reference) or allocated (by value)
         integer(INTD), private:: alloc=GFC_FALSE     !FALSE: value is stored by reference or null; TRUE: value is stored by value
         contains
          procedure, public:: construct=>ContElemConstruct !constructs a new container element, either by reference or by value
          procedure, public:: destruct=>ContElemDestruct   !destructs an existing container element (releases memory occupied by value)
          procedure, public:: get_value=>ContElemGetValue  !returns a pointer to the element value
          procedure, public:: is_empty=>ContElemIsEmpty    !returns TRUE if the element of the container is empty, FALSE otherwise
          procedure, public:: predicate=>ContElemPredicate !returns the value of a user-given predicate applied to the element
          procedure, public:: action=>ContElemAction       !acts on the element with a user-defined action
          procedure, public:: compare=>ContElemCompare     !compares the value of the element with the value of another element
          procedure, public:: print_it=>ContElemPrintIt    !prints the value of the element with a user-defined print function
        end type gfc_cont_elem_t
 !Base container:
        type, abstract, public:: gfc_container_t
         integer(INTL), private:: volume=0_INTL !volume of the container (total number of elements)
#ifndef NO_OMP
         integer(omp_lock_kind), private:: lock
#endif
         contains
          procedure, non_overridable, public:: num_elems=>ContNumElems !returns the total number of elements in the container
        end type gfc_container_t
 !Base iterator:
        type, abstract, public:: gfc_iter_t
         integer(INTD), private:: state=GFC_IT_NULL !current state of the iterator
         contains
          procedure, non_overridable, public:: get_status=>IterGetStatus !returns the status of the iterator
          procedure, non_overridable, public:: set_status=>IterSetStatus !sets the status of the iterator
          procedure(gfc_it_init_i), deferred, public:: init       !initializes the iterator (associates it with a container and sets it to the root)
          procedure(gfc_it_pointee_i), deferred, public:: pointee !returns the element currently pointed to
          procedure(gfc_it_next_i), deferred, public:: next       !proceeds to the next element of the container
          procedure(gfc_it_next_i), deferred, public:: previous   !proceeds to the previous element of the container
          procedure(gfc_it_query_i), deferred, public:: on_first  !returns GFC_TRUE if the iterator is positioned at the first element, GFC_FALSE otherwise
          procedure(gfc_it_query_i), deferred, public:: on_last   !returns GFC_TRUE if the iterator is positioned at the last element, GFC_FALSE otherwise
        end type gfc_iter_t
!ABSTRACT INTERFACES:
        abstract interface
 !GFC predicate:
         function gfc_predicate_i(obj) result(pred)
          import:: INTD
          integer(INTD):: pred       !out: predicate value: {GFC_TRUE, GFC_FALSE, GFC_ERROR, other integers}
          class(*), intent(in):: obj !in: arbitrary object
         end function gfc_predicate_i
 !GFC relational:
         function gfc_cmp_i(obj1,obj2) result(cmp)
          import:: INTD
          integer(INTD):: cmp         !out: result of the comparison (relation): See CMP parameters above
          class(*), intent(in):: obj1 !in: object 1
          class(*), intent(in):: obj2 !in: object 2
         end function gfc_cmp_i
 !GFC destructor:
         function gfc_destruct_i(obj) result(ierr)
          import:: INTD
          integer(INTD):: ierr          !out: error code (0:success)
          class(*), intent(inout):: obj !inout: object for destruction
         end function gfc_destruct_i
 !GFC action:
         function gfc_action_i(obj) result(ierr)
          import:: INTD
          integer(INTD):: ierr          !out: error code (0:success)
          class(*), intent(inout):: obj !inout: object to be acted on
         end function gfc_action_i
 !GFC printing:
         function gfc_print_i(obj,dev_id) result(ierr)
          import:: INTD
          integer(INTD):: ierr                         !out: error code (0:success)
          class(*), intent(in):: obj                   !in: arbitrary object
          integer(INTD), intent(in), optional:: dev_id !in: output device id (defaults to screen: 6)
         end function gfc_print_i
 !Deferred: GFC iterator: .init:
         function gfc_it_init_i(this,cont) result(ierr)
          import:: gfc_iter_t,gfc_container_t,INTD
          integer(INTD):: ierr                              !error code
          class(gfc_iter_t), intent(inout):: this           !GFC iterator
          class(gfc_container_t), target, intent(in):: cont !GFC container
         end function gfc_it_init_i
 !Deferred: GFC iterator: .pointee:
         function gfc_it_pointee_i(this,ierr) result(pntee)
          import:: gfc_iter_t,gfc_cont_elem_t,INTD
          class(gfc_cont_elem_t), pointer:: pntee     !container element currently pointed to by the iterator
          class(gfc_iter_t), intent(in):: this        !iterator
          integer(INTD), intent(out), optional:: ierr !error code
         end function gfc_it_pointee_i
 !Deferred: GFC iterator: .next .previous:
         function gfc_it_next_i(this,elem_p) result(ierr)
          import:: gfc_iter_t,gfc_cont_elem_t,INTD
          integer(INTD):: ierr                                            !error code
          class(gfc_iter_t), intent(inout):: this                         !GFC iterator
          class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !pointer to the container element
         end function gfc_it_next_i
 !Deferred: GFC iterator: query:
         function gfc_it_query_i(this) result(res)
          import:: gfc_iter_t,INTD
          integer(INTD):: res                  !result of the specific query
          class(gfc_iter_t), intent(in):: this !GFC iterator
         end function gfc_it_query_i
        end interface
!VISIBILITY:
        private ContElemConstruct
        private ContElemDestruct
        private ContElemGetValue
        private ContElemIsEmpty
        private ContElemPredicate
        private ContElemAction
        private ContElemCompare
        private ContElemPrintIt
        private ContNumElems
        private IterGetStatus
        private IterSetStatus

       contains
!IMPLEMENTATION:
!-------------------------------------------------------------
        subroutine ContElemConstruct(this,obj,ierr,assoc_only)
!Constructs an element of a container (fills in its contents).
         implicit none
         class(gfc_cont_elem_t), intent(inout):: this !inout: element of a container
         class(*), target, intent(in):: obj           !in: value to be stored in this element
         integer(INTD), intent(out), optional:: ierr  !out: error code (0:success)
         logical, intent(in), optional:: assoc_only   !if TRUE, <obj> will be stored by reference, otherwise by value (default)
         integer:: errc
         logical:: assoc

         errc=GFC_SUCCESS
         if(present(assoc_only)) then; assoc=assoc_only; else; assoc=.false.; endif
         if(this%is_empty()) then
          if(assoc) then
           this%value_p=>obj
           this%alloc=GFC_FALSE
          else
           allocate(this%value_p,SOURCE=obj,STAT=errc); if(errc.ne.0) errc=GFC_MEM_ALLOC_FAILED
           if(errc.eq.GFC_SUCCESS) then; this%alloc=GFC_TRUE; else; this%value_p=>NULL(); endif
          endif
         else
          errc=GFC_ELEM_NOT_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContElemConstruct
!--------------------------------------------------------
        subroutine ContElemDestruct(this,ierr,destructor)
!Destructs the contents of an element of a container (not the element itself).
!On return, the element will be empty, unless an error occurs.
!An optional explicit destructor will only be invoked for allocated values.
         implicit none
         class(gfc_cont_elem_t), intent(inout):: this     !inout: element of a container
         integer(INTD), intent(out), optional:: ierr      !out: error code (0:success)
         procedure(gfc_destruct_i), optional:: destructor !in: explicit destructor for the value of this element (for allocated only)
         integer:: errc

         errc=GFC_SUCCESS
         if(.not.this%is_empty()) then
          if(this%alloc.eq.GFC_TRUE) then
           if(present(destructor)) then
            errc=destructor(this%value_p)
            if(errc.ne.0) errc=GFC_MEM_FREE_FAILED
           endif
           deallocate(this%value_p,STAT=errc); if(errc.ne.0) errc=GFC_MEM_FREE_FAILED
          endif
          this%value_p=>NULL()
          this%alloc=GFC_FALSE
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContElemDestruct
!---------------------------------------------------------
        function ContElemGetValue(this,ierr) result(val_p)
!Returns a pointer to the element value.
         implicit none
         class(gfc_cont_elem_t), intent(in):: this   !in: element of a container
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         class(*), pointer:: val_p                   !out: pointer to the element value
         integer(INTD):: errc

         val_p=>NULL()
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
!------------------------------------------------------------------
        function ContElemPredicate(this,predicat,ierr) result(pred)
!Returns the value of a user-defined predicate evaluated on
!the given element of a container.
         implicit none
         integer(INTD):: pred                        !out: evaluated predicate value
         class(gfc_cont_elem_t), intent(in):: this   !in: element of a container
         procedure(gfc_predicate_i):: predicat       !in: user-defined predicate
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer:: errc

         errc=GFC_SUCCESS; pred=GFC_ERROR
         if(.not.this%is_empty()) then
          pred=predicat(this%value_p); if(pred.eq.GFC_ERROR) errc=GFC_ERROR
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end function ContElemPredicate
!----------------------------------------------------
        subroutine ContElemAction(this,act_func,ierr)
!Performs a user-defined action on a given element of a container.
         implicit none
         class(gfc_cont_elem_t), intent(inout):: this !inout: element of a container
         procedure(gfc_action_i):: act_func           !in: user-defined action
         integer(INTD), intent(out), optional:: ierr  !out: error code (0:success)
         integer:: errc

         errc=GFC_SUCCESS
         if(.not.this%is_empty()) then
          errc=act_func(this%value_p); if(errc.ne.0) errc=GFC_ACTION_FAILED
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContElemAction
!---------------------------------------------------------------------------
        function ContElemCompare(this,other_value,cmp_func,ierr) result(cmp)
!Compares the value of a given container element with another value
!using an appropriate comparison function.
         implicit none
         integer(INTD):: cmp                         !out: result of the comparison (see GFC_CMP_XXX above)
         class(gfc_cont_elem_t), intent(in):: this   !in: element of a container whose value is being compared
         class(*), intent(in):: other_value          !in: the other value to be compared with
         procedure(gfc_cmp_i):: cmp_func             !in: appropriate comparison function
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer:: errc

         errc=GFC_SUCCESS; cmp=GFC_CMP_NA
         if(.not.this%is_empty()) then
          cmp=cmp_func(this%value_p,other_value); if(cmp.eq.GFC_ERROR) errc=GFC_ERROR
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end function ContElemCompare
!--------------------------------------------------------------
        subroutine ContElemPrintIt(this,print_func,ierr,dev_id)
!Prints the given element of a container.
         implicit none
         class(gfc_cont_elem_t), intent(in):: this    !in: element of a container
         procedure(gfc_print_i):: print_func          !in: appropriate printing function
         integer(INTD), intent(out), optional:: ierr  !out: error code (0:success)
         integer(INTD), intent(in), optional:: dev_id !in: output device (default to screen, 6)
         integer(INTD):: dev
         integer:: errc

         errc=GFC_SUCCESS
         if(.not.this%is_empty()) then
          if(present(dev_id)) then; dev=dev_id; else; dev=6; endif
          errc=print_func(this%value_p,dev); if(errc.ne.0) errc=GFC_ACTION_FAILED
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContElemPrintIt
!------------------------------------------------------
        function ContNumElems(this,ierr) result(nelems)
!Returns the total number of elements stored in the container.
         implicit none
         integer(INTL):: nelems                      !out: total number of elements in the container
         class(gfc_container_t), intent(in):: this   !in: container
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc

         errc=GFC_SUCCESS
         nelems=this%volume; if(nelems.lt.0) errc=GFC_CORRUPTED_CONT
         if(present(ierr)) ierr=errc
         return
        end function ContNumElems
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
        function IterSetStatus(this,sts) result(ierr)
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

       end module gfc_base
