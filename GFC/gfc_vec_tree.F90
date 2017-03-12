!Generic Fortran Containers (GFC): Vector tree, combines vector and tree:
!The elements are initially inserted in a vector with an option to be
!later added in a tree, thus imposing a tree relationship on them.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2017/03/12

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

       module gfc_vec_tree
        use gfc_base
        use gfc_vector
        use gfc_tree
        use timers
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !output device
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
        integer(INTD), private:: DEBUG=0    !debugging level (0:none)
!TYPES:
 !Tree vertex with an explicit numeric id:
        type, extends(tree_vertex_t), private:: vec_tree_elem_t
         integer(INTL), private:: vertex_id=-1_INTL                     !tree vertex id (position in the vector): [0..max]
        end type vec_tree_elem_t
 !Vector element with a pointer to vec_tree_elem_t:
        type, extends(vector_elem_t), private:: tree_vec_elem_t
         class(vec_tree_elem_t), pointer, private:: tree_elem_p=>NULL() !pointer to an element of the vector tree
        end type tree_vec_elem_t
 !Vector tree:
        type, public:: vec_tree_t
         type(vector_t), private:: vertices !vector of vertices
         type(tree_t), private:: tree       !tree imposed on the vector of vertices (tree vertex values are stored by reference)
        end type vec_tree_t
 !Vector tree iterator for the vector part:
        type, extends(vector_iter_t), public:: tree_vec_iter_t
        end type tree_vec_iter_t
 !Vector tree iterator for the tree part:
        type, extends(tree_iter_t), public:: vec_tree_iter_t
        end type vec_tree_iter_t

       contains
!IMPLEMENTATION:

       end module gfc_vec_tree
