!Generic Fortran Containers (GFC): Graph
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2017/03/15

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
 !Link (edge) direction:
        integer(INTD), parameter, public:: GFC_GRAPH_DIR_NONE=0
        integer(INTD), parameter, public:: GFC_GRAPH_DIR_OUTWARD=-1
        integer(INTD), parameter, public:: GFC_GRAPH_DIR_INWARD=+1
!TYPES:
#if 0
 !Link between vertices:
        type, public:: graph_link_t
         logical, private:: directed=.FALSE.            !directed or not (if directed, the order of vertices will matter)
         integer(INTD), private:: num_verts=0           !number of vertices in the link (2 for ordinary graphs, >2 for hypergraphs)
         integer(INTD), allocatable, private:: verts(:) !vertices participating in the link
         contains
          procedure, private:: GraphLinkCtor
          generic, public:: graph_link_ctor=>GraphLinkCtor !ctor
          final:: graph_link_dtor                          !dtor
        end type graph_link_t
 !Graph vertex:
        type, public:: graph_vertex_t
         integer(INTL), private:: color=0_INTL                 !vertex color
         integer(INTL), private:: num_links=0_INTL             !number of links
         type(list_bi_t), private:: links                      !links
         contains
          procedure, private:: GraphVertexCtor
          generic, public:: graph_vertex_ctor=>GraphVertexCtor !ctor
          final:: graph_vertex_dtor                            !dtor
        end type graph_vertex_t
 !Graph container:
        type, extends(gfc_container_t), public:: graph_t
         type(vector_t), private:: vertices                    !vertices: [0..max]
         integer(INTL), private:: num_links=0_INTL             !total number of links in the graph
         contains
          procedure, public:: is_empty=>GraphIsEmpty           !returns TRUE if the graph is empty
          procedure, public:: num_vertices=>GraphNumVertices   !returns the number of vertices in the graph (graph cardinality)
          final:: graph_dtor                                   !dtor
        end type graph_t
 !Graph iterator:
        type, extends(gfc_iter_t), public:: graph_iter_t
         integer(INTL), private:: curr_vertex=-1_INTL              !current vertex number
         class(graph_vertex_t), pointer, private:: current=>NULL() !pointer to the current vertex
         class(graph_t), pointer, private:: container=>NULL()      !pointer to the graph container
         contains
          procedure, public:: init=>GraphIterInit                       !initializes the iterator by associating it with a graph container
          procedure, public:: reset=>GraphIterReset                     !resets the iterator to the first graph vertex (vertex 0)
          procedure, public:: release=>GraphIterRelease                 !releases the iterator by dissociating it from its container
          procedure, public:: pointee=>GraphIterPointee                 !returns a pointer to the current graph container element
          procedure, public:: next=>GraphIterNext                       !moves the iterator to the next graph vertex
          procedure, public:: previous=>GraphIterPrevious               !moves the iterator to the previous graph vertex
          procedure, public:: get_num_vertices=>GraphIterGetNumVertices !returns the total number of vertices in the graph
          procedure, public:: get_num_links=>GraphIterGetNumLinks       !returns the total number of links in the graph
          procedure, public:: get_links=>GraphIterGetLinks              !returns the list of links for specific graph vertices
          procedure, public:: move_to=>GraphIterMoveTo                  !moves the iterator to the specific graph vertex
          procedure, public:: append_vertex=>GraphIterAppendVertex      !appends a new vertex to the graph
          procedure, public:: append_link=>GraphIterAppendLink          !appends a new link between two or more graph vertices
          procedure, public:: merge_vertices=>GraphIterMergeVertices    !merges two or more graph vertices into a single vertex
          procedure, public:: split_vertex=>GraphIterSplitVertex        !splits a graph vertex into two or more vertices
          procedure, public:: delete_vertex=>GraphIterDeleteVertex      !deletes a specific graph vertex
          procedure, public:: delete_link=>GraphIterDeleteLink          !deletes a specific graph link
          procedure, public:: delete_all=>GraphIterDeleteAll            !deletes all graph vertices
          final:: graph_iter_dtor
        end type graph_iter_t
#endif

       end module gfc_graph
