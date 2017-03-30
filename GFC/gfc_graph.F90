!Generic Fortran Containers (GFC): Graph
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2017/03/30

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
        use gfc_list
        use gfc_dictionary
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
 !Graph vertex:
        type, public:: graph_vertex_t
         integer(INTL), private:: color=0_INTL !vertex color
         contains
          procedure, private:: GraphVertexCtor
          generic, public:: graph_vertex_ctor=>GraphVertexCtor !ctor
          procedure, public:: get_color=>GraphVertexGetColor   !returns the color of the graph vertex
          procedure, public:: compare=>GraphVertexCompare      !comparator
          final:: graph_vertex_dtor                            !dtor
        end type graph_vertex_t
#if 0
 !Link between vertices (directed/undirected edge or hyperedge):
        type, public:: graph_link_t
         integer(INTL), private:: color=0_INTL          !link color
         logical, private:: directed=.FALSE.            !directed or not (if directed, the order of vertices will matter)
         integer(INTD), private:: rank=0                !number of vertices in the link (2 for ordinary graphs, >2 for hypergraphs)
         integer(INTD), allocatable, private:: verts(:) !vertices participating in the link
         contains
          procedure, private:: GraphLinkCtor
          generic, public:: graph_link_ctor=>GraphLinkCtor !ctor
          procedure, public:: get_color=>GraphLinkGetColor !returns the link color
          procedure, public:: get_rank=>GraphLinkGetRank   !returns the number of vertices in the link (2 for ordinary edges, >2 for hyperedges)
          procedure, public:: compare=>GraphLinkCompare    !comparator
          final:: graph_link_dtor                          !dtor
        end type graph_link_t
 !Graph vertex links:
        type, private:: vert_link_ref_t
         integer(INTL), private:: num_links=0_INTL  !number of links the vertex participates in
         type(dictionary_t), private:: vert_links   !vertex link references organized as an ordered associative container
         contains
          procedure, private:: VertLinkRefCtor                      !ctor
          generic, public:: vert_link_ref_ctor=>VertLinkRefCtor     !ctors
          procedure, public:: get_num_links=>VertLinkRefGetNumLinks !returns the number of links
          procedure, public:: find_link=>VertLinkRefFindLink        !finds a specific vertex link, adds new links, deletes links (all by reference)
          procedure, public:: add_link=>VertLinkRefAddLink          !adds a link reference to the vertex
          procedure, public:: delete_link=>VertLinkRefDeleteLink    !deletes a link reference
          procedure, public:: delete_all=>VertLinkRefDeleteAll      !deletes all link references
          final: vert_link_ref_dtor                                 !dtor
        end type vert_link_ref_t
 !Graph container:
        type, extends(gfc_container_t), public:: graph_t
         type(vector_t), private:: vertices                !graph vertices stored by value: [0..N-1], N is graph cardinality
         type(list_bi_t), private:: links                  !links between vertices stored by value: Each unique (generally ordered) combination of vertices may only have one link
         type(vector_t), private:: link_ref                !vector of <vert_link_ref_t> objects for all vertices
         integer(INTL), private:: num_links=0_INTL         !total number of links in the graph
         contains
          procedure, public:: is_empty=>GraphIsEmpty                !returns TRUE if the graph is empty
          procedure, public:: get_num_vertices=>GraphGetNumVertices !returns the total number of vertices in the graph (graph cardinality)
          procedure, public:: get_num_links=>GraphGetNumLinks       !returns the total number of links in the graph
          final:: graph_dtor                                        !dtor
        end type graph_t
 !Graph iterator:
        type, extends(gfc_iter_t), public:: graph_iter_t
         integer(INTL), private:: curr_vertex=-1_INTL              !current vertex number: [0..max]
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
          procedure, public:: get_links=>GraphIterGetLinks              !returns the list of links (by reference) for specific graph vertices
          procedure, public:: move_to=>GraphIterMoveTo                  !moves the iterator to the specific graph vertex
          procedure, public:: find_link=>GraphIterFindLink              !finds a specific link in the graph
          procedure, public:: append_vertex=>GraphIterAppendVertex      !appends a new vertex to the graph
          procedure, public:: append_link=>GraphIterAppendLink          !appends a new link between two or more graph vertices
          procedure, public:: delete_vertex=>GraphIterDeleteVertex      !deletes a specific graph vertex
          procedure, public:: delete_link=>GraphIterDeleteLink          !deletes a specific graph link
          procedure, public:: delete_all=>GraphIterDeleteAll            !deletes all graph vertices and links
          procedure, public:: merge_vertices=>GraphIterMergeVertices    !merges two or more graph vertices into a single vertex
         !procedure, public:: split_vertex=>GraphIterSplitVertex        !splits a graph vertex into two or more vertices
        end type graph_iter_t
#endif
!VISIBILITY:
 !non-member:
 !graph_vertex_t:
        private GraphVertexCtor
        private GraphVertexGetColor
        private GraphVertexCompare
        public graph_vertex_dtor
#if 0
 !graph_link_t:
        private GraphLinkCtor
        private GraphLinkGetColor
        private GraphLinkGetRank
        private GraphLinkCompare
        public graph_link_dtor
 !vert_link_ref_t
        private VertLinkRefCtor
        private VertLinkRefGetNumLinks
        private VertLinkRefFindLink
        private VertLinkRefAddLink
        private VertLinkRefDeleteLink
        private VertLinkRefDeleteAll
        public vert_link_ref_dtor
 !graph_t:
        private GraphIsEmpty
        private GraphGetNumVertices
        private GraphGetNumLinks
        public graph_dtor
 !graph_iter_t:
        private GraphIterInit
        private GraphIterReset
        private GraphIterRelease
        private GraphIterPointee
        private GraphIterNext
        private GraphIterPrevious
        private GraphIterGetNumVertices
        private GraphIterGetNumLinks
        private GraphIterGetLinks
        private GraphIterMoveTo
        private GraphIterFindLink
        private GraphIterAppendVertex
        private GraphIterAppendLink
        private GraphIterDeleteVertex
        private GraphIterDeleteLink
        private GraphIterDeleteAll
        private GraphIterMergeVertices
#endif
!IMPLEMENTATION:
       contains
![non-member]=========================

![graph_vertex_t]=====================
        subroutine GraphVertexCtor(this,color,ierr)
!Ctor.
         implicit none
         class(graph_vertex_t), intent(out):: this   !out: graph vertex
         integer(INTL), intent(in):: color           !in: vertex color
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
         this%color=color
         if(present(ierr)) ierr=errc
         return
        end subroutine GraphVertexCtor
!------------------------------------------------------------
        function GraphVertexGetColor(this,ierr) result(color)
!Returns the color of the graph vertex.
         implicit none
         integer(INTL):: color                       !out: vertex color
         class(graph_vertex_t), intent(in):: this    !in: graph vertex
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
         color=this%color
         if(present(ierr)) ierr=errc
         return
        end function GraphVertexGetColor
!------------------------------------------------------------
        function GraphVertexCompare(this,another) result(cmp)
!Compares two graph vertices (comparator).
         implicit none
         integer(INTD):: cmp                         !out: result: {CMP_EQ,CMP_LT,CMP_GT,CMP_ERR}
         class(graph_vertex_t), intent(in):: this    !in: graph vertex 1
         class(graph_vertex_t), intent(in):: another !in: graph vertex 2

         if(this%color.lt.another%color) then
          cmp=CMP_LT
         elseif(this%color.gt.another%color) then
          cmp=CMP_GT
         else !equ
          cmp=CMP_EQ
         endif
         return
        end function GraphVertexCompare
!-----------------------------------------
        subroutine graph_vertex_dtor(this)
!Dtor.
         implicit none
         type(graph_vertex_t):: this

         this%color=0_INTL
         return
        end subroutine graph_vertex_dtor
![graph_link_t]=======================

![vert_link_ref_t]====================

![graph_t]============================

![graph_iter_t]=======================

       end module gfc_graph
