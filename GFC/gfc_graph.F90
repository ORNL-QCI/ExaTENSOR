!Generic Fortran Containers (GFC): Graph
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2017/02/06

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

       module gfc_graph
        use gfc_base
        use gfc_vector
        use timers
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !output device
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
        integer(INTD), private:: DEBUG=0    !debugging level (0:none)
! Edge direction:
        integer(INTD), parameter, public:: GFC_GRAPH_DIR_OUTWARD=-1
        integer(INTD), parameter, public:: GFC_GRAPH_DIR_INWARD=+1
!TYPES:
#if 0
 !Graph edge:
        type, public:: graph_edge_t
         integer(INTL), private:: vertex=-1_INTL
         contains
          procedure, private:: GraphEdgeConstruct
          generic, public:: graph_edge_ctor=>GraphEdgeConstruct
          procedure, public:: get_conn_vertex=>GraphEdgeGetConnVertex
        end type graph_edge_t
 !Directed colored graph edge (general graph edge):
        type, extends(graph_edge_t), public:: graph_edge_gen_t
         integer(INTL), private:: color=0_INTL
         contains
          procedure, private:: GraphEdgeGenConstruct
          generic, public:: graph_edge_gen_ctor=>GraphEdgeGenConstruct
          procedure, public:: get_color=>GraphEdgeGetColor
          procedure, public:: get_direction=>GraphEdgeGetDirection
        end type graph_edge_gen_t
 !Graph vertex:
        type, public:: graph_vertex_t
         integer(INTL), private:: color=0_INTL
         integer(INTL), private:: num_edges=0_INTL
         type(list_bi_t), private:: edges
         contains
          procedure, private:: GraphVertexConstruct
          generic, public:: graph_vertex_ctor=>GraphVertexConstruct
          final:: graph_vertex_dtor
        end type graph_vertex_t
 !Graph:
        type, extends(gfc_container_t), public:: graph_t
         type(vector_t), private:: vertices
         integer(INTL), private:: num_edges=0_INTL
         contains
          procedure, public:: is_empty=>GraphIsEmpty
          procedure, public:: num_vertices=>GraphNumVertices
          final:: graph_dtor
        end type graph_t
 !Graph iterator:
        type, extends(gfc_iter_t), public:: graph_iter_t
         class(graph_vertex_t), pointer, private:: current=>NULL()
         class(graph_t), pointer, private:: container=>NULL()
         contains
          procedure, public:: init=>GraphIterInit
          procedure, public:: reset=>GraphIterReset
          procedure, public:: release=>GraphIterRelease
          procedure, public:: pointee=>GraphIterPointee
          procedure, public:: next=>GraphIterNext
          procedure, public:: previous=>GraphIterPrevious
          procedure, public:: get_num_vertices=>GraphIterGetNumVertices
          procedure, public:: get_num_edges=>GraphIterGetNumEdges
          procedure, public:: move_to=>GraphIterMoveTo
          procedure, public:: get_edges=>GraphIterGetEdges
          procedure, public:: append_vertex=>GraphIterAppendVertex
          procedure, public:: append_edge=>GraphIterAppendEdge
          procedure, public:: merge_vertex=>GraphIterMergeVertex
          procedure, public:: split_vertex=>GraphIterSplitVertex
          procedure, public:: delete_vertex=>GraphIterDeleteVertex
          procedure, public:: delete_edge=>GraphIterDeleteEdge
          procedure, public:: delete_all=>GraphIterDeleteAll
          final:: graph_iter_dtor
        end type graph_iter_t
#endif
       end module gfc_graph
