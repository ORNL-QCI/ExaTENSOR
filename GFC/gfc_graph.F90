!Generic Fortran Containers (GFC): Graph
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2017/03/31

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
        use combinatoric
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
 !Link between vertices (directed/undirected edge or hyperedge):
        type, public:: graph_link_t
         integer(INTL), private:: color=0_INTL          !link color
         logical, private:: directed=.FALSE.            !directed or not (if directed, the order of vertices will matter)
         integer(INTD), private:: rank=0                !number of vertices in the link (2 for ordinary graphs, >2 for hypergraphs)
         integer(INTL), allocatable, private:: verts(:) !vertices participating in the link
         contains
          procedure, private:: GraphLinkCtor
          generic, public:: graph_link_ctor=>GraphLinkCtor !ctor
          procedure, public:: is_set=>GraphLinkIsSet       !returns TRUE if the graph link is set
          procedure, public:: get_color=>GraphLinkGetColor !returns the link color
          procedure, public:: get_rank=>GraphLinkGetRank   !returns the number of vertices in the link (2 for ordinary edges, >2 for hyperedges)
          procedure, public:: compare=>GraphLinkCompare    !comparator
          final:: graph_link_dtor                          !dtor
        end type graph_link_t
 !Graph vertex links:
        type, private:: vert_link_refs_t
         integer(INTL), private:: num_links=0_INTL  !number of links the vertex participates in
         type(dictionary_t), private:: vert_links   !vertex link references organized as an ordered associative container
         contains
          procedure, private:: VertLinkRefsCtor                      !ctor
          generic, public:: vert_link_refs_ctor=>VertLinkRefsCtor    !ctors
          procedure, public:: get_num_links=>VertLinkRefsGetNumLinks !returns the number of links
          procedure, public:: find_link=>VertLinkRefsFindLink        !finds a specific vertex link
          procedure, public:: add_link=>VertLinkRefsAddLink          !adds a link reference to the vertex
          procedure, public:: delete_link=>VertLinkRefsDeleteLink    !deletes a link reference from the vertex
          procedure, public:: delete_all=>VertLinkRefsDeleteAll      !deletes all link references on the vertex
          final:: vert_link_refs_dtor                                !dtor
        end type vert_link_refs_t
#if 0
 !Graph container:
        type, extends(gfc_container_t), public:: graph_t
         type(vector_t), private:: vertices                !graph vertices stored by value: [0..N-1], N is graph cardinality
         type(list_bi_t), private:: links                  !links between vertices stored by value: Each unique (generally ordered) combination of vertices may only have one link
         type(vector_t), private:: link_refs               !vector of <vert_link_refs_t> objects for all vertices
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
        public cmp_graph_links
 !graph_vertex_t:
        private GraphVertexCtor
        private GraphVertexGetColor
        private GraphVertexCompare
        public graph_vertex_dtor
 !graph_link_t:
        private GraphLinkCtor
        private GraphLinkIsSet
        private GraphLinkGetColor
        private GraphLinkGetRank
        private GraphLinkCompare
        public graph_link_dtor
 !vert_link_refs_t
        private VertLinkRefsCtor
        private VertLinkRefsGetNumLinks
        private VertLinkRefsFindLink
        private VertLinkRefsAddLink
        private VertLinkRefsDeleteLink
        private VertLinkRefsDeleteAll
        public vert_link_refs_dtor
#if 0
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
![non-member]==========================================
        function cmp_graph_links(obj1,obj2) result(cmp)
!Non-member comparator for graph links.
         implicit none
         integer(INTD):: cmp                 !out: comparison result: {CMP_EQ,CMP_LT,CMP_GT,CMP_ER}
         class(*), intent(in), target:: obj1 !in: graph link 1
         class(*), intent(in), target:: obj2 !in: graph link 2
         class(graph_link_t), pointer:: glp1,glp2

         glp1=>NULL(); select type(obj1); class is(graph_link_t); glp1=>obj1; end select
         glp2=>NULL(); select type(obj2); class is(graph_link_t); glp2=>obj2; end select
         if(associated(glp1).and.associated(glp2)) then
          cmp=glp1%compare(glp2)
         else
          cmp=GFC_CMP_ERR
         endif
         return
        end function cmp_graph_links
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
![graph_link_t]====================================================
        subroutine GraphLinkCtor(this,vertices,ierr,directed,color)
!Ctor.
         implicit none
         class(graph_link_t), intent(out):: this     !out: link between graph vertices
         integer(INTL), intent(in):: vertices(1:)    !in: vertices participating in the link (2 or more)
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: directed    !in: TRUE for directed links (the order of vertices will matter)
         integer(INTL), intent(in), optional:: color !in: color of the link (for multigraphs)
         integer(INTD):: sgn,errc
         integer(INTL):: ax

         errc=GFC_SUCCESS
         this%rank=size(vertices)
         if(this%rank.ge.2) then !2 for graphs, >2 for hypergraphs
          allocate(this%verts(1:this%rank),STAT=errc)
          if(errc.eq.0) then
           this%verts(:)=vertices(:)
           if(present(color)) then; this%color=color; else; this%color=0_INTL; endif
           if(present(directed)) then; this%directed=directed; else; this%directed=.FALSE.; endif
           if(.not.this%directed) then !undirected links need to be ordered
            if(this%rank.eq.2) then
             if(this%verts(1).gt.this%verts(2)) then
              ax=this%verts(1); this%verts(1)=this%verts(2); this%verts(2)=ax
             endif
            else
             call merge_sort(int(this%rank,INTL),this%verts(1:this%rank),sgn)
            endif
           endif
          else
           errc=GFC_MEM_ALLOC_FAILED
          endif
         else
          errc=GFC_INVALID_ARGS
         endif
         if(errc.ne.GFC_SUCCESS) call graph_link_dtor(this)
         if(present(ierr)) ierr=errc
         return
        end subroutine GraphLinkCtor
!-----------------------------------------------------
        function GraphLinkIsSet(this,ierr) result(res)
!Returns TRUE if the graph link is set.
         implicit none
         logical:: res                               !out: result
         class(graph_link_t), intent(in):: this      !in: graph link
         integer(INTD), intent(out), optional:: ierr !out: error code

         res=(this%rank.ge.2)
         if(present(ierr)) ierr=GFC_SUCCESS
         return
        end function GraphLinkIsSet
!----------------------------------------------------------
        function GraphLinkGetColor(this,ierr) result(color)
!Returns the color of the graph link.
         implicit none
         integer(INTL):: color                       !out: link color
         class(graph_link_t), intent(in):: this      !in: graph link
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS; color=this%color
         if(this%rank.lt.2) errc=GFC_ELEM_EMPTY
         if(present(ierr)) ierr=errc
         return
        end function GraphLinkGetColor
!-----------------------------------------------------------------
        function GraphLinkGetRank(this,ierr,directed) result(rank)
!Returns the rank of the graph link (2 for graphs, >2 for hypergraphs).
         implicit none
         integer(INTD):: rank                        !out: graph link rank
         class(graph_link_t), intent(in):: this      !in: graph link
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(out), optional:: directed   !out: TRUE if the graph link is directed
         integer(INTD):: errc

         errc=GFC_SUCCESS; rank=this%rank
         if(this%rank.ge.2) then
          if(present(directed)) directed=this%directed
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end function GraphLinkGetRank
!----------------------------------------------------------
        function GraphLinkCompare(this,another) result(cmp)
!Compares two graph links (comparator).
         implicit none
         integer(INTD):: cmp                       !out: result: {CMP_EQ,CMP_LT,CMP_GT,CMP_ERR}
         class(graph_link_t), intent(in):: this    !in: graph vertex 1
         class(graph_link_t), intent(in):: another !in: graph vertex 2
         integer(INTD):: i

         cmp=CMP_EQ
         if(this%rank.ge.2.and.another%rank.ge.2) then
          if(this%rank.lt.another%rank) then
           cmp=CMP_LT
          elseif(this%rank.gt.another%rank) then
           cmp=CMP_GT
          else
           if((.not.this%directed).and.another%directed) then
            cmp=CMP_LT
           elseif(this%directed.and.(.not.another%directed)) then
            cmp=CMP_GT
           else
            do i=1,this%rank
             if(this%verts(i).lt.another%verts(i)) then
              cmp=CMP_LT; exit
             elseif(this%verts(i).gt.another%verts(i)) then
              cmp=CMP_GT; exit
             endif
            enddo
            if(cmp.eq.CMP_EQ) then
             if(this%color.lt.another%color) then
              cmp=CMP_LT
             elseif(this%color.gt.another%color) then
              cmp=CMP_GT
             endif
            endif
           endif
          endif
         else
          cmp=CMP_ER
         endif
         return
        end function GraphLinkCompare
!---------------------------------------
        subroutine graph_link_dtor(this)
!Dtor.
         implicit none
         type(graph_link_t):: this

         if(allocated(this%verts)) deallocate(this%verts)
         this%rank=0; this%color=0_INTL; this%directed=.FALSE.
         return
        end subroutine graph_link_dtor
![vert_link_refs_t]===========================
        subroutine VertLinkRefsCtor(this,ierr)
!Ctor.
         implicit none
         class(vert_link_refs_t), intent(out):: this
         integer(INTD), intent(out), optional:: ierr

         if(present(ierr)) ierr=GFC_SUCCESS
         return
        end subroutine VertLinkRefsCtor
!--------------------------------------------------------------------
        function VertLinkRefsGetNumLinks(this,ierr) result(num_links)
!Returns the number of graph links attached to the vertex.
         implicit none
         integer(INTD):: num_links                   !out: number of links
         class(vert_link_refs_t), intent(in):: this  !in: vertex info
         integer(INTD), intent(out), optional:: ierr !out: error code

         num_links=this%num_links
         if(present(ierr)) ierr=GFC_SUCCESS
         return
        end function VertLinkRefsGetNumLinks
!------------------------------------------------------------------
        function VertLinkRefsFindLink(this,link,ierr) result(found)
!Searches for a specific link in the vertex info object.
         implicit none
         logical:: found                             !out: TRUE if found
         class(vert_link_refs_t), intent(in):: this  !in: vertex info
         class(graph_link_t), intent(in):: link      !in: graph link
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: i,errc
         type(dictionary_iter_t):: dit

         errc=GFC_SUCCESS; found=.FALSE.
         if(link%is_set()) then
          if(this%num_links.gt.0) then
           errc=dit%init(this%vert_links)
           if(errc.eq.GFC_SUCCESS) then
            errc=dit%search(GFC_DICT_JUST_FIND,cmp_graph_links,link)
            found=(errc.eq.GFC_FOUND); if(errc.eq.GFC_NOT_FOUND.or.found) errc=GFC_SUCCESS
            i=dit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=i
           endif
          endif
         else
          errc=GFC_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end function VertLinkRefsFindLink
!-----------------------------------------------------------
        function VertLinkRefsAddLink(this,link) result(ierr)
!Adds a graph link to the vertex info object.
         implicit none
         integer(INTD):: ierr                          !out: error code
         class(vert_link_refs_t), intent(inout):: this !inout: vertex info
         class(graph_link_t), intent(in):: link        !in: graph link
         integer(INTD):: i
         type(dictionary_iter_t):: dit

         if(link%is_set(ierr)) then
          ierr=dit%init(this%vert_links)
          if(ierr.eq.GFC_SUCCESS) then
           i=0 !`dummy value, because a gfc::dictionary is used instead of gfc::set
           ierr=dit%search(GFC_DICT_ADD_IF_NOT_FOUND,cmp_graph_links,link,i) !`the key <link> should be added by reference in future (when gfc::set is avaliable)
           if(ierr.ne.GFC_NOT_FOUND) then
            if(ierr.eq.GFC_FOUND) ierr=GFC_INVALID_REQUEST
           else
            ierr=GFC_SUCCESS
           endif
           i=dit%release(); if(i.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=i
          endif
         else
          ierr=GFC_INVALID_ARGS
         endif
         return
        end function VertLinkRefsAddLink
!--------------------------------------------------------------
        function VertLinkRefsDeleteLink(this,link) result(ierr)
!Deletes a graph link from the vertex info object.
         implicit none
         integer(INTD):: ierr                          !out: error code
         class(vert_link_refs_t), intent(inout):: this !inout: vertex info
         class(graph_link_t), intent(in):: link        !in: graph link
         integer(INTD):: i
         type(dictionary_iter_t):: dit

         if(link%is_set(ierr)) then
          if(this%num_links.gt.0) then
           ierr=dit%init(this%vert_links)
           if(ierr.eq.GFC_SUCCESS) then
            ierr=dit%search(GFC_DICT_DELETE_IF_FOUND,cmp_graph_links,link)
            if(ierr.ne.GFC_FOUND) then
             if(ierr.eq.GFC_NOT_FOUND) ierr=GFC_INVALID_REQUEST
            else
             ierr=GFC_SUCCESS
            endif
            i=dit%release(); if(i.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=i
           endif
          else
           ierr=GFC_INVALID_REQUEST
          endif
         else
          ierr=GFC_INVALID_ARGS
         endif
         return
        end function VertLinkRefsDeleteLink
!--------------------------------------------------------
        function VertLinkRefsDeleteAll(this) result(ierr)
!Deletes all link references on the vertex.
         implicit none
         integer(INTD):: ierr                          !out: error code
         class(vert_link_refs_t), intent(inout):: this !inout: vertex info
         type(dictionary_iter_t):: dit
         integer(INTD):: i

         ierr=dit%init(this%vert_links)
         if(ierr.eq.GFC_SUCCESS) then
          call dit%delete_all(ierr)
          i=dit%release(); if(i.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=i
         endif
         return
        end function VertLinkRefsDeleteAll
!-------------------------------------------
        subroutine vert_link_refs_dtor(this)
         implicit none
         type(vert_link_refs_t):: this
         integer(INTD):: errc

         errc=this%delete_all(); this%num_links=0_INTL
         return
        end subroutine vert_link_refs_dtor
![graph_t]============================

![graph_iter_t]=======================

       end module gfc_graph
