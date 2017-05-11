!Generic Fortran Containers (GFC): Graph
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2017/05/11

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

!NOTES:
! # GFC::graph contains graph vertices (<graph_vertex_t>) and graph links (<graph_link_t>).
!   Graph vertices are internally stored in a gfc::vector, making the dynamic type of the
!   pointer returned by the graph iterator (method .pointee) to be <vector_elem_t>.

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
          procedure, private:: GraphVertexCtor                 !ctor
          generic, public:: graph_vertex_ctor=>GraphVertexCtor !ctors
          procedure, public:: get_color=>GraphVertexGetColor   !returns the color of the graph vertex
          procedure, public:: compare=>GraphVertexCompare      !comparator
          procedure, public:: print_it=>GraphVertexPrintIt     !prints the graph vertex info
          final:: graph_vertex_dtor                            !dtor
        end type graph_vertex_t
 !Link between vertices (directed/undirected edge or hyperedge):
        type, public:: graph_link_t
         integer(INTL), private:: color=0_INTL             !link color
         logical, private:: directed=.FALSE.               !directed or not (if directed, the order of vertices will matter)
         integer(INTD), private:: rank=0                   !number of vertices in the link (2 for ordinary graphs, >2 for hypergraphs)
         integer(INTL), allocatable, private:: vertices(:) !vertices participating in the link
         contains
          procedure, private:: GraphLinkCtor               !ctor
          generic, public:: graph_link_ctor=>GraphLinkCtor !ctors
          procedure, public:: is_set=>GraphLinkIsSet       !returns TRUE if the graph link is set
          procedure, public:: get_color=>GraphLinkGetColor !returns the link color
          procedure, public:: get_rank=>GraphLinkGetRank   !returns the number of vertices in the link (2 for ordinary edges, >2 for hyperedges)
          procedure, public:: compare=>GraphLinkCompare    !comparator
          procedure, public:: print_it=>GraphLinkPrintIt   !prints the graph link info
          final:: graph_link_dtor                          !dtor
        end type graph_link_t
 !Graph vertex links:
        type, private:: vert_link_refs_t
         integer(INTL), private:: num_links=0_INTL  !number of links the vertex participates in
         type(dictionary_t), private:: vert_links   !vertex link references organized as an ordered associative container
         contains
          procedure, private:: VertLinkRefsCtor                        !ctor
          generic, public:: vert_link_refs_ctor=>VertLinkRefsCtor      !ctors
          procedure, public:: get_num_links=>VertLinkRefsGetNumLinks   !returns the number of links attached to the corresponding vertex
          procedure, public:: get_links_iter=>VertLinkRefsGetLinksIter !returns a dictionary iterator to the vertex links
          procedure, public:: find_link=>VertLinkRefsFindLink          !finds a specific vertex link
          procedure, public:: add_link=>VertLinkRefsAddLink            !adds a link reference to the vertex
          procedure, public:: delete_link=>VertLinkRefsDeleteLink      !deletes a link reference from the vertex
          procedure, public:: delete_all=>VertLinkRefsDeleteAll        !deletes all link references on the vertex
          final:: vert_link_refs_dtor                                  !dtor
        end type vert_link_refs_t
 !Graph container:
        type, extends(gfc_container_t), public:: graph_t
         type(vector_t), private:: vertices                !graph vertices stored by value: [0..N-1], N is graph cardinality
         type(vector_t), private:: link_refs               !vector of <vert_link_refs_t> objects containing link references for all vertices
         type(list_bi_t), private:: links                  !links (stored by value): Each unique (generally ordered) combination of vertices may only have one link
         integer(INTL), private:: num_links=0_INTL         !total number of links in the graph
         contains
          procedure, private:: incr_num_links_=>GraphIncrNumLinks   !PRIVATE: increments the number of graph links by one
          procedure, private:: decr_num_links_=>GraphDecrNumLinks   !PRIVATE: decrements the number of graph links by one
          procedure, public:: is_empty=>GraphIsEmpty                !returns TRUE if the graph is empty
          procedure, public:: is_set=>GraphIsSet                    !returns GFC_TRUE if the graph is set, plus additional info
          procedure, public:: get_num_vertices=>GraphGetNumVertices !returns the total number of vertices in the graph (graph cardinality)
          procedure, public:: get_num_links=>GraphGetNumLinks       !returns the total number of links in the graph
          procedure, public:: print_it=>GraphPrintIt                !prints the graph
        end type graph_t
 !Graph iterator:
        type, extends(gfc_iter_t), public:: graph_iter_t
         class(graph_t), pointer, private:: container=>NULL()       !pointer to the graph container
         type(vector_iter_t), private:: vert_it                     !vertex iterator
         type(vector_iter_t), private:: vert_ln_it                  !vertex links iterator
         type(list_iter_t), private:: link_it                       !graph links iterator
         contains
          procedure, private:: update_status_=>GraphIterUpdateStatus    !PRIVATE: updates graph iterator status ater component updates
          procedure, private:: renumber_=>GraphIterRenumber             !PRIVATE: renumbers graph vertices in graph links in case of vertex deletion
          procedure, public:: init=>GraphIterInit                       !initializes the iterator by associating it with a graph container
          procedure, public:: reset=>GraphIterReset                     !resets the iterator to the first graph vertex (vertex 0)
          procedure, public:: release=>GraphIterRelease                 !releases the iterator by dissociating it from its container
          procedure, public:: pointee=>GraphIterPointee                 !returns a pointer to the current graph container element (<vector_elem_t>)
          procedure, public:: next=>GraphIterNext                       !moves the iterator to the next graph vertex
          procedure, public:: previous=>GraphIterPrevious               !moves the iterator to the previous graph vertex
          procedure, public:: get_vertex=>GraphIterGetVertex            !returns a pointer to a specified graph vertex (<graph_vertex_t>)
          procedure, public:: get_num_vertices=>GraphIterGetNumVertices !returns the total number of vertices in the graph
          procedure, public:: get_num_links=>GraphIterGetNumLinks       !returns the total number of links in the graph
          procedure, public:: get_links=>GraphIterGetLinks              !returns the list of links (by reference) for specific graph vertices
          procedure, public:: move_to=>GraphIterMoveTo                  !moves the iterator to the specific graph vertex
          procedure, public:: find_link=>GraphIterFindLink              !finds a specific link in the graph
          procedure, public:: append_vertex=>GraphIterAppendVertex      !appends a new vertex to the graph
          procedure, public:: append_link=>GraphIterAppendLink          !appends a new link between two or more graph vertices
          procedure, public:: delete_link=>GraphIterDeleteLink          !deletes a specific graph link
          procedure, public:: delete_vertex=>GraphIterDeleteVertex      !deletes a specific graph vertex
          procedure, public:: delete_all=>GraphIterDeleteAll            !deletes all graph vertices and links
          procedure, public:: merge_vertices=>GraphIterMergeVertices    !merges two or more graph vertices into a single vertex
        end type graph_iter_t
!VISIBILITY:
 !non-member:
        public cmp_graph_links
        public print_graph_link
        public print_graph_link_from_list
 !graph_vertex_t:
        private GraphVertexCtor
        private GraphVertexGetColor
        private GraphVertexCompare
        private GraphVertexPrintIt
        public graph_vertex_dtor
 !graph_link_t:
        private GraphLinkCtor
        private GraphLinkIsSet
        private GraphLinkGetColor
        private GraphLinkGetRank
        private GraphLinkCompare
        private GraphLinkPrintIt
        public graph_link_dtor
 !vert_link_refs_t
        private VertLinkRefsCtor
        private VertLinkRefsGetNumLinks
        private VertLinkRefsGetLinksIter
        private VertLinkRefsFindLink
        private VertLinkRefsAddLink
        private VertLinkRefsDeleteLink
        private VertLinkRefsDeleteAll
        public vert_link_refs_dtor
 !graph_t:
        private GraphIncrNumLinks
        private GraphDecrNumLinks
        private GraphIsEmpty
        private GraphIsSet
        private GraphGetNumVertices
        private GraphGetNumLinks
        private GraphPrintIt
 !graph_iter_t:
        private GraphIterUpdateStatus
        private GraphIterRenumber
        private GraphIterInit
        private GraphIterReset
        private GraphIterRelease
        private GraphIterPointee
        private GraphIterNext
        private GraphIterPrevious
        private GraphIterGetVertex
        private GraphIterGetNumVertices
        private GraphIterGetNumLinks
        private GraphIterGetLinks
        private GraphIterMoveTo
        private GraphIterFindLink
        private GraphIterAppendVertex
        private GraphIterAppendLink
        private GraphIterDeleteLink
        private GraphIterDeleteVertex
        private GraphIterDeleteAll
        private GraphIterMergeVertices
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
!---------------------------------------------------------
        function print_graph_link(obj,dev_id) result(ierr)
!Prints graph_link_t.
         implicit none
         integer(INTD):: ierr                         !out: error code
         class(*), intent(in), target:: obj           !in: graph_link_t object
         integer(INTD), intent(in), optional:: dev_id !in: output device (default to screen)
         integer(INTD):: devo

         ierr=GFC_SUCCESS
         if(present(dev_id)) then; devo=dev_id; else; devo=6; endif
         select type(obj)
         class is(graph_link_t)
          call obj%print_it(ierr,devo)
         class default
          ierr=GFC_ACTION_FAILED
         end select
         return
        end function print_graph_link
!-------------------------------------------------------------------
        function print_graph_link_from_list(obj,dev_id) result(ierr)
!Prints graph_link_t from the list of links.
         implicit none
         integer(INTD):: ierr                         !out: error code
         class(*), intent(in), target:: obj           !in: list_elem_t<graph_link_t> object
         integer(INTD), intent(in), optional:: dev_id !in: output device (default to screen)
         class(*), pointer:: up
         integer(INTD):: devo

         ierr=GFC_SUCCESS
         if(present(dev_id)) then; devo=dev_id; else; devo=6; endif
         select type(obj)
         class is(list_elem_t)
          up=>obj%get_value(ierr)
          if(ierr.eq.GFC_SUCCESS) ierr=print_graph_link(up,devo)
         class default
          ierr=GFC_ACTION_FAILED
         end select
         return
        end function print_graph_link_from_list
![graph_vertex_t]==================================
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
!-----------------------------------------------
        subroutine GraphVertexPrintIt(this,ierr)
!Prints the graph vertex info.
         implicit none
         class(graph_vertex_t), intent(in):: this    !in: graph vertex
         integer(INTD), intent(out), optional:: ierr !out: error code

         write(*,'("Vertex color = ",i13)') this%color
         if(present(ierr)) ierr=GFC_SUCCESS
         return
        end subroutine GraphVertexPrintIt
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
          allocate(this%vertices(1:this%rank),STAT=errc)
          if(errc.eq.0) then
           this%vertices(:)=vertices(:)
           if(present(color)) then; this%color=color; else; this%color=0_INTL; endif
           if(present(directed)) then; this%directed=directed; else; this%directed=.FALSE.; endif
           if(.not.this%directed) then !undirected links need to be ordered
            if(this%rank.eq.2) then
             if(this%vertices(1).gt.this%vertices(2)) then
              ax=this%vertices(1); this%vertices(1)=this%vertices(2); this%vertices(2)=ax
             endif
            else
             call merge_sort(int(this%rank,INTL),this%vertices(1:this%rank),sgn)
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
         integer(INTD):: errc

         errc=GFC_SUCCESS; res=(this%rank.ge.2)
         if(this%rank.eq.1) errc=GFC_CORRUPTED_CONT
         if(present(ierr)) ierr=errc
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
             if(this%vertices(i).lt.another%vertices(i)) then
              cmp=CMP_LT; exit
             elseif(this%vertices(i).gt.another%vertices(i)) then
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
!----------------------------------------------------
        subroutine GraphLinkPrintIt(this,ierr,dev_id)
!Prints the graph link.
         implicit none
         class(graph_link_t), intent(in):: this       !in: graph link
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD), intent(in), optional:: dev_id !in: output device (defaults to screen)
         integer(INTD):: errc,devo

         if(present(dev_id)) then; devo=dev_id; else; devo=6; endif
         if(this%is_set(errc)) then
          if(this%directed) then
           write(devo,'("Link (directed, color = ",i13,"):",16(1x,i13))') this%color,this%vertices(1:this%rank)
          else
           write(devo,'("Link (color = ",i13,"):",16(1x,i13))') this%color,this%vertices(1:this%rank)
          endif
         else
          errc=GFC_ELEM_EMPTY
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine GraphLinkPrintIt
!---------------------------------------
        subroutine graph_link_dtor(this)
!Dtor.
         implicit none
         type(graph_link_t):: this

         if(allocated(this%vertices)) deallocate(this%vertices)
         this%rank=0; this%color=0_INTL; this%directed=.FALSE.
         return
        end subroutine graph_link_dtor
![vert_link_refs_t]===========================
        subroutine VertLinkRefsCtor(this,ierr)
!Ctor.
         implicit none
         class(vert_link_refs_t), intent(out):: this !out: vertex info
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         this%num_links=0_INTL
         call this%vert_links%set_key_storage(GFC_BY_REF,errc)
         if(present(ierr)) ierr=errc
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
!---------------------------------------------------------------------
        function VertLinkRefsGetLinksIter(this,dict_iter) result(ierr)
!Returns a dictionary iterator to the vertex links.
         implicit none
         integer(INTD):: ierr                                !out: error code
         class(vert_link_refs_t), intent(in):: this          !in: vertex info
         class(dictionary_iter_t), intent(inout):: dict_iter !inout: dictionary iterator

         ierr=dict_iter%get_status()
         if(ierr.eq.GFC_IT_NULL) ierr=dict_iter%init(this%vert_links)
         return
        end function VertLinkRefsGetLinksIter
!------------------------------------------------------------------------------
        function VertLinkRefsFindLink(this,link,ierr,list_link_p) result(found)
!Searches for a specific graph link in the vertex info object.
         implicit none
         logical:: found                             !out: TRUE if found
         class(vert_link_refs_t), intent(in):: this  !in: vertex info
         class(graph_link_t), intent(in):: link      !in: graph link
         integer(INTD), intent(out), optional:: ierr !out: error code
         class(list_elem_t), intent(out), pointer, optional:: list_link_p !out: pointer to the list element containing the graph link
         integer(INTD):: i,errc
         type(dictionary_iter_t):: dit
         class(*), pointer:: up

         errc=GFC_SUCCESS; found=.FALSE.
         if(link%is_set()) then
          if(this%num_links.gt.0) then
           errc=dit%init(this%vert_links)
           if(errc.eq.GFC_SUCCESS) then
            up=>NULL()
            errc=dit%search(GFC_DICT_FETCH_IF_FOUND,cmp_graph_links,link,value_out=up)
            found=(errc.eq.GFC_FOUND)
            if(errc.eq.GFC_NOT_FOUND.or.found) then
             errc=GFC_SUCCESS
             if(found.and.present(list_link_p)) then
              if(associated(up)) then
               select type(up)
               class is(list_elem_t)
                list_link_p=>up
               class default
                errc=GFC_CORRUPTED_CONT
               end select
              else
               errc=GFC_ERROR
              endif
             endif
            endif
            up=>NULL()
            i=dit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=i
           endif
          endif
         else
          errc=GFC_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end function VertLinkRefsFindLink
!--------------------------------------------------------------------
        function VertLinkRefsAddLink(this,link,link_ref) result(ierr)
!Adds a graph link to the vertex info object.
         implicit none
         integer(INTD):: ierr                                !out: error code
         class(vert_link_refs_t), intent(inout):: this       !inout: vertex info
         class(graph_link_t), intent(in), target:: link      !in: graph link
         class(list_elem_t), intent(in), pointer:: link_ref  !in: pointer to the list element where value <link> is stored by value
         integer(INTD):: i
         type(dictionary_iter_t):: dit

         if(link%is_set(ierr)) then
          if(ierr.eq.GFC_SUCCESS) then
           ierr=dit%init(this%vert_links)
           if(ierr.eq.GFC_SUCCESS) then
            ierr=dit%search(GFC_DICT_ADD_IF_NOT_FOUND,cmp_graph_links,link,link_ref,GFC_BY_REF)
            if(ierr.eq.GFC_NOT_FOUND) then
             this%num_links=this%num_links+1_INTL
             ierr=GFC_SUCCESS
            else
             if(ierr.eq.GFC_FOUND) then
              ierr=GFC_INVALID_REQUEST !link already exists
             endif
            endif
            i=dit%release(); if(i.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=i
           endif
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
          if(ierr.eq.GFC_SUCCESS) then
           if(this%num_links.gt.0) then
            ierr=dit%init(this%vert_links)
            if(ierr.eq.GFC_SUCCESS) then
             ierr=dit%search(GFC_DICT_DELETE_IF_FOUND,cmp_graph_links,link)
             if(ierr.eq.GFC_FOUND) then
              this%num_links=this%num_links-1_INTL
              ierr=GFC_SUCCESS
             else
              if(ierr.eq.GFC_NOT_FOUND) ierr=GFC_INVALID_REQUEST !link not found at the vertex
             endif
             i=dit%release(); if(i.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=i
            endif
           else
            ierr=GFC_INVALID_REQUEST
           endif
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
          call dit%delete_all(ierr); if(ierr.eq.GFC_SUCCESS) this%num_links=0_INTL
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
![graph_t]================================
        subroutine GraphIncrNumLinks(this)
         implicit none
         class(graph_t), intent(inout):: this !inout: graph

         this%num_links=this%num_links+1_INTL
         return
        end subroutine GraphIncrNumLinks
!-----------------------------------------
        subroutine GraphDecrNumLinks(this)
         implicit none
         class(graph_t), intent(inout):: this !inout: graph

         this%num_links=this%num_links-1_INTL
         return
        end subroutine GraphDecrNumLinks
!----------------------------------------------
        function GraphIsEmpty(this) result(res)
!Returns TRUE if the graph is empty.
         implicit none
         integer(INTD):: res               !out: result: {GFC_TRUE,GFC_FALSE,error}
         class(graph_t), intent(in):: this !in: graph
         integer(INTL):: nv

         nv=this%get_num_vertices(res)
         if(res.eq.GFC_SUCCESS) then
          if(nv.gt.0_INTL) then; res=GFC_FALSE; else; res=GFC_TRUE; endif
         endif
         return
        end function GraphIsEmpty
!------------------------------------------------------------------------
        function GraphIsSet(this,ierr,num_vertices,num_links) result(res)
!Returns TRUE if the graph is set, plus additional info.
         implicit none
         logical:: res                                       !out: result
         class(graph_t), intent(in):: this                   !in: graph
         integer(INTD), intent(out), optional:: ierr         !out: error code
         integer(INTL), intent(out), optional:: num_vertices !out: number of vertices in the graph
         integer(INTL), intent(out), optional:: num_links    !out: number of (unique) links in the graph
         integer(INTD):: errc
         integer(INTL):: nv

         res=.FALSE.; nv=this%get_num_vertices(errc)
         if(errc.eq.GFC_SUCCESS) then
          res=(nv.gt.0_INTL)
          if(present(num_vertices)) num_vertices=nv
          if(present(num_links)) num_links=this%num_links
         endif
         if(present(ierr)) ierr=errc
         return
        end function GraphIsSet
!-------------------------------------------------------------------
        function GraphGetNumVertices(this,ierr) result(num_vertices)
!Returns the total number of vertices in the graph.
         implicit none
         integer(INTL):: num_vertices                !out: number of vertices
         class(graph_t), intent(in):: this           !in: graph
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         num_vertices=this%vertices%length(errc)
         if(present(ierr)) ierr=errc
         return
        end function GraphGetNumVertices
!-------------------------------------------------------------
        function GraphGetNumLinks(this,ierr) result(num_links)
!Returns the total number of (unique) links in the graph.
!Each (hyper-)link between two or more vertices is counted only once.
         implicit none
         integer(INTL):: num_links                   !out: number of links
         class(graph_t), intent(in):: this           !in: graph
         integer(INTD), intent(out), optional:: ierr !out: error code

         num_links=this%num_links
         if(present(ierr)) ierr=GFC_SUCCESS
         return
        end function GraphGetNumLinks
!-----------------------------------------
        subroutine GraphPrintIt(this,ierr)
!Prints the graph.
         implicit none
         class(graph_t), intent(in), target:: this   !in: graph
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: i,errc
         integer(INTL):: nv,nl,l
         type(graph_iter_t):: git
         class(graph_vertex_t), pointer:: gvp
         class(graph_link_t), pointer:: glp
         type(list_bi_t):: links
         type(list_iter_t):: lit
         class(*), pointer:: up

         write(*,'("#Printing graph:")')
         errc=git%init(this)
         if(errc.eq.GFC_SUCCESS) then
          nv=git%get_num_vertices(errc)
          if(errc.eq.GFC_SUCCESS.and.nv.gt.0) then
           do l=0_INTL,nv-1_INTL !loop over vertices
            gvp=>git%get_vertex(l,errc); if(errc.ne.GFC_SUCCESS) exit
            write(*,'("Vertex ",i13,1x)',ADVANCE='NO') l; call gvp%print_it()
            call git%get_links((/l/),links,nl,errc); if(errc.ne.GFC_SUCCESS) exit
            if(nl.gt.0) then
             errc=lit%init(links)
             do while(errc.eq.GFC_SUCCESS)
              up=>lit%get_value(errc); if(errc.ne.GFC_SUCCESS) exit
              glp=>NULL(); select type(up); class is(graph_link_t); glp=>up; end select
              if(.not.associated(glp)) then; errc=GFC_ERROR; exit; endif
              write(*,'(1x)',ADVANCE='NO'); call glp%print_it(errc); if(errc.ne.GFC_SUCCESS) exit
              errc=lit%next()
             enddo
             if(errc.eq.GFC_NO_MOVE) errc=GFC_SUCCESS
             i=lit%delete_all(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=i
             if(errc.eq.GFC_SUCCESS) then
              errc=lit%release(); if(errc.ne.GFC_SUCCESS) exit
             else
              exit
             endif
            endif
           enddo
          endif
          i=git%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=i
         endif
         write(*,'("#End of printing")')
         if(present(ierr)) ierr=errc
         return
        end subroutine GraphPrintIt
![graph_iter_t]===============================
        subroutine GraphIterUpdateStatus(this)
!Updates the iterator status after component updates.
         implicit none
         class(graph_iter_t), intent(inout):: this !inout: iterator
         integer(INTD):: errc

         errc=this%set_status_(this%vert_it%get_status())
         return
        end subroutine GraphIterUpdateStatus
!--------------------------------------------------------
        subroutine GraphIterRenumber(this,vertex_id,ierr)
!Renumbers graph vertices in graph links in case of vertex deletion.
!All vertex id's larger than <vertex_id> will be decremented by one.
         implicit none
         class(graph_iter_t), intent(inout):: this   !inout: graph iterator
         integer(INTL), intent(in):: vertex_id       !in: deleted vertex id
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: i,errc
         integer(INTL):: vlv,vll,vid
         class(*), pointer:: up
         class(graph_link_t), pointer:: glp

         if(vertex_id.ge.0_INTL) then
          vlv=this%vert_it%get_length(errc)
          if(errc.eq.GFC_SUCCESS) then
           vll=this%vert_ln_it%get_length(errc)
           if(errc.eq.GFC_SUCCESS) then
            if(vlv.eq.vll) then
             if(vlv.gt.0_INTL.and.vertex_id.lt.vlv.and.this%get_num_links().gt.0_INTL) then
              errc=this%link_it%reset()
              lloop: do while(errc.eq.GFC_SUCCESS)
               up=>this%link_it%get_value(errc); if(errc.ne.GFC_SUCCESS) exit lloop
               if(.not.associated(up)) then; errc=GFC_ERROR; exit lloop; endif
               select type(up); class is(graph_link_t); glp=>up; end select
               if(.not.associated(glp)) then; errc=GFC_CORRUPTED_CONT; exit lloop; endif
               do i=1,glp%rank
                vid=glp%vertices(i); if(vid.eq.vertex_id) then; errc=GFC_ERROR; exit lloop; endif
                if(vid.gt.vertex_id) glp%vertices(i)=vid-1_INTL
               enddo
               errc=this%link_it%next()
              enddo lloop
              if(errc.eq.GFC_NO_MOVE) errc=this%link_it%reset()
              glp=>NULL(); up=>NULL()
             endif
            else
             errc=GFC_CORRUPTED_CONT
            endif
           endif
          endif
         else
          errc=GFC_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine GraphIterRenumber
!-----------------------------------------------------
        function GraphIterInit(this,cont) result(ierr)
!Initializes the iterator with its container and resets it to the beginning.
         implicit none
         integer(INTD):: ierr                              !out: error code
         class(graph_iter_t), intent(inout):: this         !inout: iterator
         class(gfc_container_t), target, intent(in):: cont !in: container

         select type(cont)
         class is(graph_t)
          this%container=>cont
          ierr=this%vert_it%init(this%container%vertices)
          if(ierr.eq.GFC_SUCCESS) then
           ierr=this%vert_ln_it%init(this%container%link_refs)
           if(ierr.eq.GFC_SUCCESS) then
            ierr=this%link_it%init(this%container%links)
            if(ierr.eq.GFC_SUCCESS) ierr=this%reset()
           endif
          endif
         class default
          ierr=GFC_INVALID_ARGS
         end select
         return
        end function GraphIterInit
!-------------------------------------------------
        function GraphIterReset(this) result(ierr)
!Resets the iterator to the beginning of the container.
         implicit none
         integer(INTD):: ierr                      !out: error code
         class(graph_iter_t), intent(inout):: this !inout: iterator

         ierr=this%vert_it%reset()
         if(ierr.eq.GFC_SUCCESS) then
          ierr=this%vert_ln_it%reset()
          if(ierr.eq.GFC_SUCCESS) ierr=this%link_it%reset()
          call this%reset_count()
         endif
         call this%update_status_()
         return
        end function GraphIterReset
!---------------------------------------------------
        function GraphIterRelease(this) result(ierr)
!Dissociates the iterator from its container.
         implicit none
         integer(INTD):: ierr                      !out: error code
         class(graph_iter_t), intent(inout):: this !inout: iterator
         integer(INTD):: errc

         call this%reset_count()
         ierr=this%link_it%release()
         errc=this%vert_ln_it%release(); if(ierr.eq.GFC_SUCCESS.and.errc.ne.GFC_SUCCESS) ierr=errc
         errc=this%link_it%release(); if(ierr.eq.GFC_SUCCESS.and.errc.ne.GFC_SUCCESS) ierr=errc
         this%container=>NULL()
         errc=this%set_status_(GFC_IT_NULL); if(ierr.eq.GFC_SUCCESS.and.errc.ne.GFC_SUCCESS) ierr=errc
         return
        end function GraphIterRelease
!---------------------------------------------------------
        function GraphIterPointee(this,ierr) result(pntee)
!Returns a pointer to the current container element.
!Note that the dynamic type of the returned pointer is <vector_elem_t>
!since the graph vertices (<graph_vertex_t>) are stored in a vector.
         implicit none
         class(gfc_cont_elem_t), pointer:: pntee     !out: container element currently pointed to by the iterator: <vector_elem_t>
         class(graph_iter_t), intent(in):: this      !in: iterator
         integer(INTD), intent(out), optional:: ierr !out: error code

         pntee=>NULL()
         if(present(ierr)) then
          pntee=>this%vert_it%pointee(ierr)
         else
          pntee=>this%vert_it%pointee()
         endif
         return
        end function GraphIterPointee
!-------------------------------------------------------
        function GraphIterNext(this,elem_p) result(ierr)
!Moves the iterator to the next element of the graph. If <elem_p> is present,
!the next element will be returned in <elem_p> without moving the iterator.
         implicit none
         integer(INTD):: ierr                                            !out: error code
         class(graph_iter_t), intent(inout):: this                       !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element

         if(present(elem_p)) then
          ierr=this%vert_it%next(elem_p)
         else
          ierr=this%vert_it%next()
          if(ierr.eq.GFC_SUCCESS) ierr=this%vert_ln_it%next()
         endif
         call this%update_status_()
         return
        end function GraphIterNext
!-----------------------------------------------------------
        function GraphIterPrevious(this,elem_p) result(ierr)
!Moves the iterator to the previous element of the graph. If <elem_p> is present,
!the previous element will be returned in <elem_p> without moving the iterator.
         implicit none
         integer(INTD):: ierr                                            !out: error code
         class(graph_iter_t), intent(inout):: this                       !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element

         if(present(elem_p)) then
          ierr=this%vert_it%previous(elem_p)
         else
          ierr=this%vert_it%previous()
          if(ierr.eq.GFC_SUCCESS) ierr=this%vert_ln_it%previous()
         endif
         call this%update_status_()
         return
        end function GraphIterPrevious
!----------------------------------------------------------------------
        function GraphIterGetVertex(this,vert_id,ierr) result(vertex_p)
!Returns a pointer to a specific graph vertex.
         implicit none
         class(graph_vertex_t), pointer:: vertex_p   !out: pointer to the specified graph vertex
         class(graph_iter_t), intent(in):: this      !in: iterator
         integer(INTL), intent(in):: vert_id         !in: vertex id: [0..max]
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         class(*), pointer:: up

         vertex_p=>NULL()
         up=>this%vert_it%element_value(vert_id,errc)
         if(errc.eq.GFC_SUCCESS) then
          if(associated(up)) then
           select type(up); class is(graph_vertex_t); vertex_p=>up; end select
           if(.not.associated(vertex_p)) errc=GFC_CORRUPTED_CONT
           up=>NULL()
          else
           errc=GFC_ERROR
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end function GraphIterGetVertex
!-----------------------------------------------------------------------
        function GraphIterGetNumVertices(this,ierr) result(num_vertices)
!Returns the total number of vertices in the graph.
         implicit none
         integer(INTL):: num_vertices                !out: number of vertices
         class(graph_iter_t), intent(in):: this      !in: iterator
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         num_vertices=0_INTL; errc=this%get_status()
         if(errc.eq.GFC_IT_ACTIVE.or.errc.eq.GFC_IT_DONE) then
          num_vertices=this%container%get_num_vertices(errc)
         elseif(errc.eq.GFC_IT_EMPTY) then
          errc=GFC_SUCCESS
         endif
         if(present(ierr)) ierr=errc
         return
        end function GraphIterGetNumVertices
!-----------------------------------------------------------------
        function GraphIterGetNumLinks(this,ierr) result(num_links)
!Returns the total number of (unique) links in the graph.
         implicit none
         integer(INTL):: num_links                   !out: number of links
         class(graph_iter_t), intent(in):: this      !in: iterator
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         num_links=0_INTL; errc=this%get_status()
         if(errc.eq.GFC_IT_ACTIVE.or.errc.eq.GFC_IT_DONE) then
          num_links=this%container%get_num_links(errc)
         elseif(errc.eq.GFC_IT_EMPTY) then
          errc=GFC_SUCCESS
         endif
         if(present(ierr)) ierr=errc
         return
        end function GraphIterGetNumLinks
!---------------------------------------------------------------------------
        subroutine GraphIterGetLinks(this,vertices,link_list,num_links,ierr)
!Returns a list of graph links (by reference) for specific graph vertices.
!<link_list> is allowed to be non-empty on input (will be appended).
!The iterator position is kept intact.
         implicit none
         class(graph_iter_t), intent(in):: this      !in: graph iterator
         integer(INTL), intent(in):: vertices(1:)    !in: vertices of interest (vertex id's): vertex id = [0..max]
         type(list_bi_t), intent(inout):: link_list  !inout: list of graph links (by reference) attached to the specified graph vertices
         integer(INTL), intent(out):: num_links      !out: number of links extracted from the graph
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: i,n,errc
         class(*), pointer:: up
         class(vert_link_refs_t), pointer:: vlrp
         type(list_iter_t):: lit
         type(dictionary_iter_t):: dit

         num_links=0_INTL; errc=this%get_status()
         if(errc.eq.GFC_IT_ACTIVE.or.errc.eq.GFC_IT_DONE) then
          errc=lit%init(link_list)
          if(errc.eq.GFC_SUCCESS) then
           n=size(vertices)
           aloop: do i=1,n
            up=>this%vert_ln_it%element_value(vertices(i),errc); if(errc.ne.GFC_SUCCESS) exit aloop !up => <vert_link_refs_t>
            if(.not.associated(up)) then; errc=GFC_ERROR; exit aloop; endif
            select type(up); class is(vert_link_refs_t); vlrp=>up; end select
            if(.not.associated(vlrp)) then; errc=GFC_CORRUPTED_CONT; exit aloop; endif
            if(vlrp%num_links.gt.0_INTL) then
             errc=dit%init(vlrp%vert_links); if(errc.ne.GFC_SUCCESS) exit aloop
             do while(errc.eq.GFC_SUCCESS)
              up=>dit%get_key(errc); if(.not.associated(up)) errc=GFC_ERROR
              if(errc.eq.GFC_SUCCESS) then
               num_links=num_links+1_INTL
               errc=lit%append(up,assoc_only=.TRUE.) !up => <graph_link_t>
               if(errc.eq.GFC_SUCCESS) errc=dit%move_in_order(GFC_DICT_SUCCESSOR)
              endif
             enddo
             if(errc.ne.GFC_NO_MOVE) exit aloop
             errc=dit%release(); if(errc.ne.GFC_SUCCESS) exit aloop
            endif
           enddo aloop
           i=lit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=i
           vlrp=>NULL(); up=>NULL()
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine GraphIterGetLinks
!------------------------------------------------------------
        function GraphIterMoveTo(this,vertex_id) result(ierr)
!Moves the graph iterator to a specific graph vertex.
         implicit none
         integer(INTD):: ierr                      !out: error code: {GFC_SUCCESS,errors}
         class(graph_iter_t), intent(inout):: this !inout: graph iterator
         integer(INTL), intent(in):: vertex_id     !in: vertex id: [0..max]

         ierr=this%vert_it%move_to(vertex_id)
         if(ierr.eq.GFC_SUCCESS) ierr=this%vert_ln_it%move_to(vertex_id)
         call this%update_status_()
         return
        end function GraphIterMoveTo
!---------------------------------------------------------------
        function GraphIterFindLink(this,link,ierr) result(found)
!Finds a specific link in the graph. Passing an empty graph here
!will not raise an error.
         implicit none
         logical:: found                             !out: found or not: {TRUE|FALSE}
         class(graph_iter_t), intent(in):: this      !in: graph iterator
         class(graph_link_t), intent(in):: link      !in: graph link being searched for
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer(INTL):: vid
         class(*), pointer:: up
         class(vert_link_refs_t), pointer:: vlrp

         found=.FALSE.; errc=this%get_status()
         if(errc.eq.GFC_IT_ACTIVE.or.errc.eq.GFC_IT_DONE) then
          if(link%is_set(errc)) then
           if(errc.eq.GFC_SUCCESS) then
            vid=link%vertices(1) !first vertex of the link
            up=>this%vert_ln_it%element_value(vid,errc)
            if(errc.eq.GFC_SUCCESS) then
             if(associated(up)) then
              select type(up); class is(vert_link_refs_t); vlrp=>up; end select
              if(associated(vlrp)) then
               found=vlrp%find_link(link,errc)
              else
               errc=GFC_CORRUPTED_CONT
              endif
             else
              errc=GFC_ERROR
             endif
            endif
            vlrp=>NULL(); up=>NULL()
           endif
          else
           errc=GFC_INVALID_ARGS
          endif
         elseif(errc.eq.GFC_IT_EMPTY) then
          if(.not.link%is_set(errc)) errc=GFC_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end function GraphIterFindLink
!---------------------------------------------------------------
        function GraphIterAppendVertex(this,vertex) result(ierr)
!Appends a new vertex to the graph.
         implicit none
         integer(INTD):: ierr                               !out: error code
         class(graph_iter_t), intent(inout):: this          !inout: graph iterator
         class(graph_vertex_t), intent(in), target:: vertex !in: new vertex
         type(vert_link_refs_t):: vlr

         ierr=this%vert_it%append(vertex)
         if(ierr.eq.GFC_SUCCESS) then
          call vlr%vert_link_refs_ctor(ierr)
          if(ierr.eq.GFC_SUCCESS) ierr=this%vert_ln_it%append(vlr)
         endif
         call this%update_status_()
         return
        end function GraphIterAppendVertex
!-----------------------------------------------------------
        function GraphIterAppendLink(this,link) result(ierr)
!Appends a new link to the graph.
         implicit none
         integer(INTD):: ierr                      !out: error code
         class(graph_iter_t), intent(inout):: this !inout: graph iterator
         class(graph_link_t), intent(in):: link    !in: new graph link
         integer(INTD):: i
         integer(INTL):: vid
         class(gfc_cont_elem_t), pointer:: gcp
         class(list_elem_t), pointer:: lep
         class(*), pointer:: up
         class(graph_link_t), pointer:: glp
         class(vert_link_refs_t), pointer:: vlrp

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE.or.ierr.eq.GFC_IT_DONE) then
          if(link%is_set(ierr)) then
           if(ierr.eq.GFC_SUCCESS) then
            ierr=this%link_it%append(link) !append the new graph link to the list of graph links (by value)
            !write(*,'("GraphIterAppendLink:link_list.append:ierr = ",i9,1x)',ADVANCE='NO') ierr; call link%print_it() !debug[-666]
            if(ierr.eq.GFC_SUCCESS) then
             call this%container%incr_num_links_()
             ierr=this%link_it%reset_back()
             if(ierr.eq.GFC_SUCCESS) then
              gcp=>this%link_it%pointee(ierr)
              if(.not.associated(gcp)) ierr=GFC_ERROR
              if(ierr.eq.GFC_SUCCESS) then
               up=>gcp%get_value(ierr) !graph_link_t
               if(.not.associated(up).and.ierr.eq.GFC_SUCCESS) ierr=GFC_ERROR
               if(ierr.eq.GFC_SUCCESS) then
                glp=>NULL(); select type(up); class is(graph_link_t); glp=>up; end select
                if(associated(glp)) then
                 select type(gcp)
                 class is(list_elem_t)
                  lep=>gcp
                  do i=1,glp%rank !loop over participating vertices
                   vid=glp%vertices(i)
                   up=>this%vert_ln_it%element_value(vid,ierr); if(ierr.ne.GFC_SUCCESS) exit
                   if(associated(up)) then
                    select type(up); class is(vert_link_refs_t); vlrp=>up; end select
                    if(associated(vlrp)) then
                     ierr=vlrp%add_link(glp,lep) !add a reference to the graph link to each participating vertex
                     if(ierr.ne.GFC_SUCCESS) exit
                    else
                     ierr=GFC_CORRUPTED_CONT; exit
                    endif
                   else
                    ierr=GFC_ERROR; exit
                   endif
                  enddo
                  vlrp=>NULL(); up=>NULL(); lep=>NULL()
                 class default
                  ierr=GFC_CORRUPTED_CONT
                 end select
                 glp=>NULL()
                else
                 ierr=GFC_CORRUPTED_CONT
                endif
               endif
              endif
              gcp=>NULL()
             endif
            endif
           endif
          else
           ierr=GFC_INVALID_ARGS
          endif
         endif
         call this%update_status_()
         return
        end function GraphIterAppendLink
!-----------------------------------------------------------
        function GraphIterDeleteLink(this,link) result(ierr)
!Deletes a specific graph link.
         implicit none
         integer(INTD):: ierr                      !out: error code
         class(graph_iter_t), intent(inout):: this !inout: graph iterator
         class(graph_link_t), intent(in):: link    !in: graph link to delete
         integer(INTD):: i
         integer(INTL):: vid
         class(*), pointer:: up
         class(vert_link_refs_t), pointer:: vlrp
         class(list_elem_t), pointer:: lep
         logical:: found

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE.or.ierr.eq.GFC_IT_DONE) then
          if(link%is_set(ierr)) then
           if(ierr.eq.GFC_SUCCESS) then
            lep=>NULL()
            vloop: do i=link%rank,1,-1 !loop over the vertices participating in the link
             vid=link%vertices(i) !vertex id
             up=>this%vert_ln_it%element_value(vid,ierr)
             if(.not.associated(up)) ierr=GFC_ERROR
             if(ierr.eq.GFC_SUCCESS) then
              vlrp=>NULL(); select type(up); class is(vert_link_refs_t); vlrp=>up; end select
              if(associated(vlrp)) then
               found=vlrp%find_link(link,ierr,lep)
               if(found.and.ierr.eq.GFC_SUCCESS) then
                ierr=vlrp%delete_link(link); if(ierr.ne.GFC_SUCCESS) exit vloop !delete the link reference from each vertex
                if(i.eq.1) then !last vertex
                 call this%link_it%jump_(lep) !jump to the list element to be deleted
                 ierr=this%link_it%delete()   !delete that list element from the list
                 if(ierr.ne.GFC_SUCCESS) exit vloop
                 call this%container%decr_num_links_()
                endif
               else
                ierr=GFC_CORRUPTED_CONT; exit vloop
               endif
              else
               ierr=GFC_CORRUPTED_CONT; exit vloop
              endif
             else
              exit vloop
             endif
            enddo vloop
            lep=>NULL(); vlrp=>NULL(); up=>NULL()
           endif
          else
           ierr=GFC_INVALID_ARGS
          endif
         endif
         call this%update_status_()
         return
        end function GraphIterDeleteLink
!------------------------------------------------------------------
        function GraphIterDeleteVertex(this,vertex_id) result(ierr)
!Deletes a specific graph vertex (and all its links). If the vertex
!being deleted is not the last vertex, all subsequent graph vertices
!will be renumbered by one less. The graph links will be renumbered as well.
!The iterator position will not change, unless the iterator is positioned
!on the vertex being deleted, in which case the iterator will acquire
!a status of GFC_IT_DONE.
         implicit none
         integer(INTD):: ierr                      !out: error code
         class(graph_iter_t), intent(inout):: this !inout: graph iterator
         integer(INTL), intent(in):: vertex_id     !in: id of the vertex to be deleted: [0..max]
         integer(INTD):: i
         integer(INTL):: vid
         class(*), pointer:: up
         class(vert_link_refs_t), pointer:: vlrp
         type(dictionary_iter_t):: dit
         class(graph_link_t), pointer:: glp

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE.or.ierr.eq.GFC_IT_DONE) then
          vid=this%vert_it%get_offset() !save current iterator position, -1 for GFC_IT_DONE
          ierr=this%vert_it%move_to(vertex_id)
          if(ierr.eq.GFC_SUCCESS) then
           ierr=this%vert_ln_it%move_to(vertex_id)
           if(ierr.eq.GFC_SUCCESS) then
            up=>this%vert_ln_it%get_value(ierr)
            if(ierr.eq.GFC_SUCCESS) then
             if(associated(up)) then
              select type(up); class is(vert_link_refs_t); vlrp=>up; end select
              if(associated(vlrp)) then
               ierr=vlrp%get_links_iter(dit)
               dloop: do while(ierr.eq.GFC_SUCCESS)
                up=>dit%get_key(ierr); if(ierr.ne.GFC_SUCCESS) exit dloop
                if(.not.associated(up)) then; ierr=GFC_ERROR; exit dloop; endif
                select type(up); class is(graph_link_t); glp=>up; end select
                ierr=this%delete_link(glp); if(ierr.ne.GFC_SUCCESS) exit dloop
                ierr=dit%next()
               enddo dloop
               if(ierr.eq.GFC_NO_MOVE) then
                ierr=dit%release()
                if(ierr.eq.GFC_SUCCESS) ierr=this%vert_ln_it%delete()
                if(ierr.eq.GFC_SUCCESS) ierr=this%vert_it%delete()
                if(ierr.eq.GFC_SUCCESS) call this%renumber_(vertex_id,ierr)
               else
                i=dit%release()
               endif
              else
               ierr=GFC_CORRUPTED_CONT
              endif
             else
              ierr=GFC_ERROR
             endif
            endif
            glp=>NULL(); vlrp=>NULL(); up=>NULL()
           endif
           if(vid.ge.0_INTL.and.vid.ne.vertex_id) then
            if(vid.gt.vertex_id) vid=vid-1_INTL !vertex numeration shift by one
            i=this%vert_it%move_to(vid); if(i.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=i
            i=this%vert_ln_it%move_to(vid); if(i.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=i
           else
            i=this%vert_it%set_status_(GFC_IT_DONE); if(i.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=i
            i=this%vert_ln_it%set_status_(GFC_IT_DONE); if(i.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=i
           endif
          endif
         endif
         call this%update_status_()
         return
        end function GraphIterDeleteVertex
!-----------------------------------------------------
        function GraphIterDeleteAll(this) result(ierr)
!Deletes all graph vertices (and links).
         implicit none
         integer(INTD):: ierr                      !out: error code
         class(graph_iter_t), intent(inout):: this !inout: graph iterator
         class(*), pointer:: up
         class(graph_link_t), pointer:: glp
         integer(INTD):: errc

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_DONE) then
          ierr=this%reset(); if(ierr.eq.GFC_SUCCESS) ierr=this%get_status()
         endif
         if(ierr.eq.GFC_IT_ACTIVE) then
          ierr=GFC_SUCCESS
          if(this%get_num_links().gt.0_INTL) then
           ierr=this%link_it%reset()
           do while(ierr.eq.GFC_SUCCESS)
            up=>this%link_it%get_value(ierr); if(ierr.ne.GFC_SUCCESS) exit
            if(.not.associated(up)) then; ierr=GFC_ERROR; exit; endif
            glp=>NULL(); select type(up); class is(graph_link_t); glp=>up; end select
            if(.not.associated(glp)) then; ierr=GFC_CORRUPTED_CONT; exit; endif
            ierr=this%delete_link(glp); if(ierr.ne.GFC_SUCCESS) exit
            ierr=this%link_it%next()
           enddo
           if(ierr.eq.GFC_NO_MOVE) ierr=GFC_SUCCESS
           glp=>NULL(); up=>NULL()
          endif
          errc=this%vert_ln_it%delete_all(); if(errc.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=errc
          errc=this%vert_it%delete_all(); if(errc.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=errc
         elseif(ierr.eq.GFC_IT_EMPTY) then
          ierr=GFC_SUCCESS
         endif
         call this%update_status_()
         return
        end function GraphIterDeleteAll
!-------------------------------------------------------------------
        subroutine GraphIterMergeVertices(this,vertex1,vertex2,ierr)
!Merges two graph vertices into a single one. The vertex numeration
!may change, in which case all relevant graph links will be updated.
         class(graph_iter_t), intent(inout):: this   !inout: graph iterator
         integer(INTL), intent(in):: vertex1         !in: vertex id 1
         integer(INTL), intent(in):: vertex2         !in: vertex id 2
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
         !`Write
         if(present(ierr)) ierr=errc
         return
        end subroutine GraphIterMergeVertices

       end module gfc_graph
!TESTING====================
       module gfc_graph_test
        use gfc_base
        use gfc_graph
        use timers, only: thread_wtime
        implicit none
        private

        public test_gfc_graph

       contains

        function test_gfc_graph(perf,dev_out) result(ierr)
         implicit none
         integer(INTD):: ierr
         real(8), intent(out):: perf
         integer(INTD), intent(in), optional:: dev_out
         integer(INTD), parameter:: MAX_VERTICES=100000
         integer(INTD):: jo,i
         integer(INTL):: v0,v1,v2
         type(graph_vertex_t):: vrt
         type(graph_link_t):: lnk
         type(graph_t):: graph
         type(graph_iter_t):: git
         real(8):: tms,tm

         if(present(dev_out)) then; jo=dev_out; else; jo=6; endif
         perf=0d0; tms=thread_wtime()
         ierr=git%init(graph); if(ierr.ne.GFC_SUCCESS) then; ierr=1; return; endif
         call vrt%graph_vertex_ctor(5_INTL,ierr); if(ierr.ne.GFC_SUCCESS) then; ierr=2; return; endif
         do i=0,MAX_VERTICES-1
          ierr=git%append_vertex(vrt); if(ierr.ne.GFC_SUCCESS) then; ierr=3; return; endif
         enddo
         do i=0,MAX_VERTICES-1
          v0=int(i,INTL)
          v1=int(mod(i+1,MAX_VERTICES),INTL)
          v2=int(mod(i+2,MAX_VERTICES),INTL)
          call lnk%graph_link_ctor((/v0,v1/),ierr); if(ierr.ne.GFC_SUCCESS) then; ierr=4; return; endif
          !write(*,*); write(*,'("New Vertex",i7,": link: ")',ADVANCE='NO') i; call lnk%print_it() !debug
          ierr=git%append_link(lnk); if(ierr.ne.GFC_SUCCESS) then; print *,ierr; ierr=5; return; endif
          call lnk%graph_link_ctor((/v0,v2/),ierr); if(ierr.ne.GFC_SUCCESS) then; ierr=6; return; endif
          !write(*,*); write(*,'("New Vertex",i7,": link: ")',ADVANCE='NO') i; call lnk%print_it() !debug
          ierr=git%append_link(lnk); if(ierr.ne.GFC_SUCCESS) then; print *,ierr; ierr=7; return; endif
         enddo
         !call graph%print_it(ierr); if(ierr.ne.GFC_SUCCESS) then; ierr=8; return; endif !debug
         ierr=git%delete_all(); if(ierr.ne.GFC_SUCCESS) then; ierr=9; return; endif
         ierr=git%release(); if(ierr.ne.GFC_SUCCESS) then; ierr=10; return; endif
         tm=thread_wtime(tms); perf=dble(MAX_VERTICES)/tm
         return
        end function test_gfc_graph

       end module gfc_graph_test
