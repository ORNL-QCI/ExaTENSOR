!Generic Fortran Containers:: Tree.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2016-02-18 (started 2016-02-17)
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
 !Tree vertex with a one-directional sibling linked list:
        type, extends(gfc_cont_elem_t), public:: tree_vertex_uni_t
         class(tree_vertex_uni_t), pointer, private:: parent=>NULL()
         class(tree_vertex_uni_t), pointer, private:: next_sibling=>NULL()
         class(tree_vertex_uni_t), pointer, private:: first_child=>NULL()
         integer(INTD), private:: num_child=0
         contains
          procedure, public:: add_subtree=>TreeUniAddSubtree       !adds a subtree to the given vertex in the last child position
          procedure, public:: delete_subtree=>TreeUniDeleteSubtree !deletes the subtree originating at the given vertex
          procedure, public:: move_subtree=>TreeUniMoveSubtree     !moves the subtree originating at the given vertex to another vertex as the last child of the latter
        end type tree_vertex_uni_t
 !Tree vertex with a bi-directional sibling linked list:
        type, extends(tree_vertex_uni_t), public:: tree_vertex_bi_t
         class(tree_vertex_uni_t), pointer, private:: prev_sibling=>NULL()
        end type tree_vertex_bi_t
 !Tree:
        type, extends(gfc_container_t), public:: tree_t
         contains
          procedure, public:: num_elems=>TreeNumElems
        end type tree_t
 !Tree iterator:
        type, extends(gfc_iter_t), public:: tree_iter_t
         contains
          
        end type tree_iter_t

!GLOBAL DATA:

!VISIBILITY:

       contains
!IMPLEMENTATION:
!------------------------------

       end module tree
