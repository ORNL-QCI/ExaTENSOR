!Generic Fortran Containers:: Tree.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2016-02-20 (started 2016-02-17)
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
         integer(INTL), private:: num_child=0
         contains
          procedure, non_overridable, public:: num_children=>TreeVertexNumChildren !returns the total number of children
          procedure, public:: num_siblings=>TreeVertexNumSiblings !returns the total number of siblings
        end type tree_vertex_t
 !Tree:
        type, extends(gfc_container_t), public:: tree_t
         class(tree_vertex_t), pointer, private:: root=>NULL() !root element
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
          procedure, public:: scan=>TreeIterScan
          procedure, public:: add_subtree=>TreeIterAddSubtree
          procedure, public:: delete_subtree=>TreeIterDeleteSubtree
          procedure, public:: move_subtree=>TreeIterMoveSubtree
        end type tree_iter_t
!GLOBAL DATA:
!VISIBILITY:
 !Interfaces:
        public gfc_predicate_i
        public gfc_cmp_i
        public gfc_destruct_i
        public gfc_action_i
        public gfc_print_i
 !Procedures:
        private TreeVertexNumChildren
        private TreeVertexNumSiblings
        private TreeIterInit
        private TreeIterPointee
        private TreeIterNext
        private TreeIterPrevious
        private TreeIterOnFirst
        private TreeIterOnLast
        private TreeIterScan
        private TreeIterAddSubtree
        private TreeIterDeleteSubtree
        private TreeIterMoveSubtree

       contains
!IMPLEMENTATION:
!---------------------------------------------------------------
        function TreeVertexNumChildren(this,ierr) result(nchild)
!Returns the total number of children attached to the tree vertex: O(1) worst.
         implicit none
         integer(INTL):: nchild                      !out: number of children
         class(tree_vertex_t), intent(in):: this     !in: tree vertex
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc

         errc=GFC_SUCCESS
         nchild=this%num_child; if(nchild.lt.0) errc=GFC_CORRUPTED_CONT
         if(present(ierr)) ierr=errc
         return
        end function TreeVertexNumChildren
!--------------------------------------------------------------
        function TreeVertexNumSiblings(this,ierr) result(nsibl)
!Returns the total number of siblings for a tree vertex: O(N) worst.
         implicit none
         integer(INTL):: nsibl                       !out: number of children
         class(tree_vertex_t), intent(in):: this     !in: tree vertex
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc
         class(tree_vertex_t), pointer:: tvp

         errc=GFC_SUCCESS; nsibl=0
         tvp=>this%prev_sibling
         do while(associated(tvp))
          nsibl=nsibl+1; tvp=>tvp%prev_sibling
         enddo
         tvp=>this%next_sibling
         do while(associated(tvp))
          nsibl=nsibl+1; tvp=>tvp%next_sibling
         enddo
         if(present(ierr)) ierr=errc
         return
        end function TreeVertexNumSiblings
!----------------------------------------------------
        function TreeIterInit(this,cont) result(ierr)
!Initializes an iterator.
         implicit none
         integer(INTD):: ierr                              !out: error code (0:success)
         class(tree_iter_t), intent(inout):: this          !inout: iterator
         class(gfc_container_t), target, intent(in):: cont !in: container

         ierr=GFC_SUCCESS
         select type(cont)
         class is (tree_t)
          this%container=>cont
          this%current=>cont%root
          if(associated(cont%root)) then
           ierr=this%set_status(GFC_IT_ACTIVE) !non-empty iterator
          else
           ierr=this%set_status(GFC_IT_EMPTY) !empty iterator
          endif
         class default
          ierr=GFC_INVALID_ARGS
         end select
         return
        end function TreeIterInit
!--------------------------------------------------------
        function TreeIterPointee(this,ierr) result(pntee)
!Returns the container element the iterator is currently pointing to.
         implicit none
         class(gfc_cont_elem_t), pointer:: pntee     !out: container element currently pointed to by the iterator
         class(tree_iter_t), intent(in):: this       !in: iterator
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc

         errc=this%get_status()
         if(errc.eq.GFC_IT_ACTIVE) then
          pntee=>this%current; errc=GFC_SUCCESS
         else
          pntee=>NULL()
         endif
         if(present(ierr)) ierr=errc
         return
        end function TreeIterPointee
!------------------------------------------------------
        function TreeIterNext(this,elem_p) result(ierr)
!If <elem_p> is absent, the iterator is moved to the next element, if any.
!If <elem_p> is present, the iterator simply returns the next element in <elem_p>.
         implicit none
         integer(INTD):: ierr                                            !out:error code
         class(tree_iter_t), intent(inout):: this                        !inout: GFC iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element

         ierr=GFC_SUCCESS
         return
        end function TreeIterNext
!----------------------------------------------------------
        function TreeIterPrevious(this,elem_p) result(ierr)
!If <elem_p> is absent, the iterator is moved to the previous element, if any.
!If <elem_p> is present, the iterator simply returns the previous element in <elem_p>.
         implicit none
         integer(INTD):: ierr                                            !out:error code
         class(tree_iter_t), intent(inout):: this                        !inout: GFC iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element

         ierr=GFC_SUCCESS
         return
        end function TreeIterPrevious
!-------------------------------------------------
        function TreeIterOnFirst(this) result(res)
!Returns GFC_TRUE if the iterator is positioned on the first element of the container.
         implicit none
         integer(INTD):: res                      !out: result
         class(tree_iter_t), intent(inout):: this !in: iterator
         class(gfc_cont_elem_t), pointer:: tvp

         res=this%get_status()
         if(res.eq.GFC_IT_ACTIVE) then
          res=this%previous(tvp); tvp=>NULL()
          if(res.eq.GFC_SUCCESS) then
           res=GFC_FALSE
          elseif(res.eq.GFC_END_CONT) then
           res=GFC_TRUE
          endif
         endif
         return
        end function TreeIterOnFirst
!------------------------------------------------
        function TreeIterOnLast(this) result(res)
!Returns GFC_TRUE if the iterator is positioned on the last element of the container.
         implicit none
         integer(INTD):: res                      !out: result
         class(tree_iter_t), intent(inout):: this !in: iterator
         class(gfc_cont_elem_t), pointer:: tvp

         res=this%get_status()
         if(res.eq.GFC_IT_ACTIVE) then
          res=this%next(tvp); tvp=>NULL()
          if(res.eq.GFC_SUCCESS) then
           res=GFC_FALSE
          elseif(res.eq.GFC_END_CONT) then
           res=GFC_TRUE
          endif
         endif
         return
        end function TreeIterOnLast

       end module tree
