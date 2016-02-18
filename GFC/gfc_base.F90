!Generic Fortran Containers (GFC): Base
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2016-02-18 (started 2016-02-17)
!Copyright (C) 2016 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2016 Oak Ridge National Laboratory (UT-Battelle)
!LICENSE: GNU GPL v.2
!DESIGN:
! # A GFC container is a structured collection of objects of any class;
!   Objects of different classes can be stored in the same container.
! # Objects can be added to the container either by value or by reference;
!   The same container may have elements added both ways.
! # Each container has an associated iterator for scanning over its elements;
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
!   Each specific class of containers has its own iterator class because
!   of different linkage between the elements of container.
! # The element deletion operation may require a user-supplied destructor which
!   will release all resources occupied by the object except the object itself.
!   In order to completely destroy the object stored by value in the container,
!   a DEALLOCATE(object) operation will be required, which is done automatically
!   on all objects stored by value.
!NOTES:
! # This implementation of Generic Fortran containers heavily relies on
!   dynamic polymorphism, thus introducing some memory overhead which
!   can be significant when storing small objects!
       module gfc_base
        use dil_basic
        implicit none
        public
!PARAMETERS:
 !Error codes:
        integer(INTD), parameter:: GFC_SUCCESS=0           !success
        integer(INTD), parameter:: GFC_ERROR=-666          !unspecified error
        integer(INTD), parameter:: GFC_INVALID_ARGS=-1     !invalid arguments
        integer(INTD), parameter:: GFC_NULL_CONT=-2        !uninitialized container
        integer(INTD), parameter:: GFC_EMPTY_CONT=-3       !empty container
        integer(INTD), parameter:: GFC_MEM_ALLOC_FAILED=-4 !memory allocation failed
        integer(INTD), parameter:: GFC_MEM_FREE_FAILED=-5  !memory deallocation failed
        integer(INTD), parameter:: GFC_CORRUPTED_CONT=-6   !corrupted container
        integer(INTD), parameter:: GFC_ELEM_EMPTY=-7       !element of the container is empty
        integer(INTD), parameter:: GFC_ELEM_NOT_EMPTY=-8   !element of the container is not empty
        integer(INTD), parameter:: GFC_ACTION_FAILED=-9    !user-defined action failed on an element
        integer(INTD), parameter:: GFC_END_CONT=1          !end of container
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
        integer(INTD), parameter:: GFC_IT_EMPTY=-1   !empty iterator
        integer(INTD), parameter:: GFC_IT_ACTIVE=0   !active iterator
        integer(INTD), parameter:: GFC_IT_DONE=1     !pass the end of the container (done)
!TYPES:
 !Element of a container:
        type, public:: gfc_cont_elem_t
         class(*), pointer, private:: value_p=>NULL() !element value (data)
         integer(INTD), private:: alloc=GFC_FALSE     !FALSE: value is stored by reference or null; TRUE: value is stored by value
         contains
          procedure, public:: construct=>ContElemConstruct !constructs a new container element, either by reference or by value
          procedure, public:: destruct=>ContElemDestruct   !destructs an existing container element (releases memory occupied by value)
          procedure, public:: is_empty=>ContElemEmpty      !returns TRUE if the element of the container is empty, FALSE otherwise
          procedure, public:: predicate=>ContElemPredicate !returns the value of a user-given predicate applied to the element
          procedure, public:: action=>ContElemAction       !acts on the element with a user-defined action
          procedure, public:: compare=>ContElemCompare     !compares the value of the element with the value of another element
          procedure, public:: print_it=>ContElemPrintIt    !prints the value of the element with a user-defined print function
        end type gfc_cont_elem_t
 !Base container:
        type, abstract, public:: gfc_container_t
         class(gfc_cont_elem_t), pointer, private:: root=>NULL() !root element of the container (defines the real type of the container)
         integer(INTL), private:: volume=0_INTL                  !volume of the container (total number of elements)
         contains
          procedure(gfc_cont_num_elems_i), deferred, public:: num_elems !returns the total number of elements in the container
        end type gfc_container_t
 !Base iterator:
        type, abstract, public:: gfc_iter_t
         class(gfc_cont_elem_t), pointer, private:: curr_elem=>NULL() !container element currently pointed to
         class(gfc_container_t), pointer, private:: container=>NULL() !associated container
         integer(INTD), private:: state=GFC_IT_EMPTY                  !state of the iterator
         contains
          procedure(gfc_it_init_i), deferred, public:: init      !initializes the iterator to the root of the container
          procedure(gfc_it_next_i), deferred, public:: next      !proceeds to the next element of the container
          procedure(gfc_it_next_i), deferred, public:: previous  !proceeds to the previous element of the container
          procedure(gfc_it_query_i), deferred, public:: is_empty !returns GFC_TRUE if the iterator is empty, GFC_FALSE otherwise
          procedure(gfc_it_query_i), deferred, public:: is_done  !returns GFC_TRUE if the iterator is beyond the end of the container, GFC_FALSE otherwise
          procedure(gfc_it_query_i), deferred, public:: on_first !returns GFC_TRUE if the iterator is positioned at the beginning, GFC_FALSE otherwise
          procedure(gfc_it_query_i), deferred, public:: on_last  !returns GFC_TRUE if the iterator is positioned at the last element, GFC_FALSE otherwise
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
 !Deferred: GFC container: .num_elems:
         function gfc_cont_num_elems_i(this) result(nl)
          import:: gfc_container_t,INTL
          class(gfc_container_t), intent(in):: this !GFC container
          integer(INTL):: nl                        !total number of elements or error (negative)
         end function gfc_cont_num_elems_i
 !Deferred: GFC iterator: .init:
         function gfc_it_init_i(this,cont) result(ierr)
          import:: gfc_iter_t,gfc_container_t,INTD
          class(gfc_iter_t), intent(inout):: this           !GFC iterator
          class(gfc_container_t), target, intent(in):: cont !GFC container
          integer(INTD):: ierr                              !error code
         end function gfc_it_init_i
 !Deferred: GFC iterator: .next .previous:
         function gfc_it_next_i(this) result(ierr)
          import:: gfc_iter_t,INTD
          class(gfc_iter_t), intent(inout):: this !GFC iterator
          integer(INTD):: ierr                    !error code
         end function gfc_it_next_i
 !Deferred: GFC iterator: query:
         function gfc_it_query_i(this) result(res)
          import:: gfc_iter_t,INTD
          class(gfc_iter_t), intent(in):: this !GFC iterator
          integer(INTD):: res                  !result
         end function gfc_it_query_i
        end interface
!VISIBILITY:
        private ContElemConstruct
        private ContElemDestruct
        private ContElemEmpty
        private ContElemPredicate
        private ContElemAction
        private ContElemCompare
        private ContElemPrintIt

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
         if(this%empty()) then
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
         if(.not.this%empty()) then
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
!------------------------------------------------
        function ContElemEmpty(this) result(empt)
!Returns TRUE if the element of a container is empty, FALSE otherwise.
         implicit none
         logical:: empt                            !out: empty (true) or not (false)
         class(gfc_cont_elem_t), intent(in):: this !in: element of a container

         empt=.not.associated(this%value_p)
         return
        end function ContElemEmpty
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
         if(.not.this%empty()) then
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
         if(.not.this%empty()) then
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
         if(.not.this%empty()) then
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
         if(.not.this%empty()) then
          if(present(dev_id)) then; dev=dev_id; else; dev=6; endif
          errc=print_func(this%value_p,dev); if(errc.ne.0) errc=GFC_ACTION_FAILED
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContElemPrintIt

       end module gfc_base
