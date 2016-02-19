!Generic Fortran Containers:: Tree.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2016-02-19 (started 2016-02-17)
!Copyright (C) 2016 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2016 Oak Ridge National Laboratory (UT-Battelle)
!LICENSE: GNU GPL v.2
       module tree
        use gfc_base
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !output device
        logical, private:: VERBOSE=.true.   !verbositiy for errors
        integer(INTD), private:: DEBUG=0    !debugging level (0:none)
!TYPES:
 !Tree vertex with a bi-directional sibling linked list:
        type, extends(gfc_cont_elem_t), public:: tree_vertex_t
         class(tree_vertex_t), pointer, private:: parent=>NULL()
         class(tree_vertex_t), pointer, private:: next_sibling=>NULL()
         class(tree_vertex_t), pointer, private:: prev_sibling=>NULL()
         class(tree_vertex_t), pointer, private:: first_child=>NULL()
         integer(INTD), private:: num_child=0
         contains
          procedure, public:: num_children=>TreeVertexNumChildren !returns the total number of children
          procedure, public:: num_siblings=>TreeVertexNumSiblings !returns the total number of siblings
        end type tree_vertex_t
 !Tree:
        type, extends(gfc_container_t), public:: tree_t
         class(tree_vertex_t), pointer, private:: root=>NULL() !root element
         contains
          procedure, public:: is_null=>TreeIsNull     !returns TRUE if container is uninitialized, FALSE otherwise
          procedure, public:: num_elems=>TreeNumElems !returns the total number of elements in the container
        end type tree_t
 !Tree iterator:
        type, extends(gfc_iter_t), public:: tree_iter_t
         class(tree_vertex_t), pointer, private:: current=>NULL()
         class(tree_t), pointer, private:: container=>NULL()
         contains
          procedure, public:: init=>TreeIterInit
          procedure, public:: pointee=>TreeIterPointee
          procedure, public:: next=>TreeIterNext
          procedure, public:: previous=>TreeIterPrevious
          procedure, public:: on_first=>TreeIterOnFirst
          procedure, public:: on_last=>TreeIterOnLast
          procedure, public:: is_null=>TreeIterIsNull
          procedure, public:: is_empty=>TreeIterIsEmpty
          procedure, public:: is_done=>TreeIterIsDone
          procedure, public:: scan=>TreeIterScan
          procedure, public:: add_subtree=>TreeIterAddSubtree
          procedure, public:: delete_subtree=>TreeIterDeleteSubtree
          procedure, public:: move_subtree=>TreeIterMoveSubtree
        end type tree_iter_t
!GLOBAL DATA:
!VISIBILITY:
       contains
!IMPLEMENTATION:
!------------------------------

       end module tree
