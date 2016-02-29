!Generic Fortran Containers:: Linked list.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2016-02-29 (started 2016-02-28)
!Copyright (C) 2016 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2016 Oak Ridge National Laboratory (UT-Battelle)
!LICENSE: GNU GPL v.2 (or higher).
!NOTES:
! # A list is a linked derivative of an abstract (unlinked) GFC container.
! # All accesses, updates, and scans on a list are performed via
!   a list iterator associated with the list. Multiple iterators
!   can be associated with a list at the same time.
       module list
        use gfc_base
        use timers
        implicit none
        private
!PARAMETERS:
!Basic:
        integer(INTD), private:: CONS_OUT=6 !output device
        logical, private:: VERBOSE=.true.   !verbositiy for errors
        integer(INTD), private:: DEBUG=0    !debugging level (0:none)
!TYPES:
 !Linked list element:
        type, extends(gfc_cont_elem_t), public:: list_elem_t
         class(list_elem_t), pointer, private:: next_elem=>NULL()
         class(list_elem_t), pointer, private:: prev_elem=>NULL()
        contains
         procedure, public:: is_first=>ListElemIsFirst !returns GFC_TRUE if the element is the first in the list
         procedure, public:: is_last=>ListElemIsLast   !returns GFC_TRUE if the element is the last in the list
        end type list_elem_t
 !Linked list:
        type, extends(gfc_container_t), public:: list_bi_t
         class(list_elem_t), pointer, private:: first_elem=>NULL() !first element of the linked list
         class(list_elem_t), pointer, private:: last_elem=>NULL()  !last element of the linked list
        end type list_bi_t
 !List iterator:
        type, extends(gfc_iter_t), public:: list_iter_t
         class(list_elem_t), pointer, private:: current=>NULL() !current element of the linked list
         class(list_bi_t), pointer, private:: container=>NULL() !linked list associated with the iterator
        contains
         procedure, public:: init=>ListIterInit                 !initializes the iterator by associating it with a list and resetting to the beginning
         procedure, public:: reset=>ListIterReset               !resets the iterator to the beginning of the list
         procedure, public:: release=>ListIterRelease           !releases the iterator (dissociates it from the container)
         procedure, public:: pointee=>ListIterPointee           !returns the container element currently pointed to by the iterator
         procedure, public:: next=>ListIterNext                 !moves the iterator to the next list element
         procedure, public:: previous=>ListIterPrevious         !moves the iterator to the previous list element
         procedure, public:: push_top=>ListIterPushTop          !inserts a new element at the beginning of the container
         procedure, public:: push_end=>ListIterPushEnd          !inserts a new element at the end of the container
         procedure, public:: insert_elem=>ListIterInsertElem    !inserts a new element at the current position of the container
         procedure, public:: insert_list=>ListIterInsertList    !inserts another linked list at the current position of the container
         generic, public:: insert=>ListIterInsertElem,ListIterInsertList !generic
         procedure, public:: extract=>ListIterExtract           !extracts an element from the list (deletion is optional)
         procedure, public:: split=>ListIterSplit               !splits the list into two parts at the current position of the container
         procedure, public:: delete=>ListIterDelete             !deletes an element or multiple elements starting from the current position
        end type list_iter_t
!INTERFACES:
!VISIBILITY:
        private ListElemIsFirst
        private ListElemIsLast
        private ListIterInit
        private ListIterReset
        private ListIterRelease
        private ListIterPointee
        private ListIterNext
        private ListIterPrevious
        private ListIterPushTop
        private ListIterPushEnd
        private ListIterInsertElem
        private ListIterInsertList
        private ListIterExtract
        private ListIterSplit
        private ListIterDelete

       contains
!IMPLEMENTATION:
!----------------------------

       end module list
