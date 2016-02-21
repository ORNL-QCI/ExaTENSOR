!Generic Fortran Containers:: Tree.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2016-02-21 (started 2016-02-17)
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
 !Tree iterator:
        integer(INTD), parameter, private:: TREE_IT_DOWN=1  !downward direction
        integer(INTD), parameter, private:: TREE_IT_RIGHT=2 !right direction
        integer(INTD), parameter, private:: TREE_IT_UP=3    !upward direction
        integer(INTD), parameter, private:: TREE_IT_LEFT=4  !left direction
!TYPES:
 !Tree vertex with a bi-directional sibling ring (linked) list:
        type, extends(gfc_cont_elem_t), public:: tree_vertex_t
         class(tree_vertex_t), pointer, private:: parent=>NULL()       !link to the parent vertex
         class(tree_vertex_t), pointer, private:: next_sibling=>NULL() !link to the next sibling vertex (in a ring)
         class(tree_vertex_t), pointer, private:: prev_sibling=>NULL() !link to the previous sibling vertex (in a ring)
         class(tree_vertex_t), pointer, private:: first_child=>NULL()  !link to the first child vertex
         integer(INTL), private:: num_child=0                          !total number of children vertices
         contains
          procedure, non_overridable, public:: num_children=>TreeVertexNumChildren   !returns the total number of children
          procedure, non_overridable, public:: num_siblings=>TreeVertexNumSiblings   !returns the total number of siblings
          procedure, non_overridable, public:: first_sibling=>TreeVertexFirstSibling !returns GFC_TRUE if the vertex is the first in the sibling (ring) list
          procedure, non_overridable, public:: last_sibling=>TreeVertexLastSibling   !returns GFC_TRUE if the vertex is the last in the sibling (ring) list
        end type tree_vertex_t
 !Tree (all operations on the tree are performend via an iterator):
        type, extends(gfc_container_t), public:: tree_t
         class(tree_vertex_t), pointer, private:: root=>NULL() !root element (beginning)
        end type tree_t
 !Tree iterator:
        type, extends(gfc_iter_t), public:: tree_iter_t
         class(tree_vertex_t), pointer, private:: current=>NULL()   !currently pointed element of the container
         class(tree_t), pointer, private:: container=>NULL()        !container
         contains
          procedure, public:: init=>TreeIterInit                    !associates the iterator with a container and sets its position to the root element
          procedure, public:: reset=>TreeIterReset                  !resets the iterator to the beginning of the container (root element)
          procedure, public:: pointee=>TreeIterPointee              !returns a pointer to the container element currently in focus
          procedure, public:: next=>TreeIterNext                    !moves the iterator to the next element
          procedure, public:: previous=>TreeIterPrevious            !moves the iterator to the previous element
!          procedure, public:: scan=>TreeIterScan                    !traverses the container with an optional action
!          procedure, public:: add_element=>TreeIterAddElement       !adds a new (child) element to the element of the container currently pointed to
!          procedure, public:: add_subtree=>TreeIterAddSubtree       !adds a subtree to the element of the container currently pointed to
!          procedure, public:: delete_subtree=>TreeIterDeleteSubtree !deletes a subtree beginning from the currently pointed element of the container
!          procedure, public:: move_subtree=>TreeIterMoveSubtree     !moves the subtree beginning from the currently pointed element of the container to another location
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
        private TreeVertexFirstSibling
        private TreeVertexLastSibling
        private TreeIterInit
        private TreeIterReset
        private TreeIterPointee
        private TreeIterNext
        private TreeIterPrevious
!        private TreeIterScan
!        private TreeIterAddElement
!        private TreeIterAddSubtree
!        private TreeIterDeleteSubtree
!        private TreeIterMoveSubtree

       contains
!IMPLEMENTATION:
!---------------------------------------------------------------
        function TreeVertexNumChildren(this,ierr) result(nchild)
!Returns the total number of children attached to the tree vertex.
!Complexity: O(1) worst.
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
!Returns the total number of siblings for a tree vertex.
!Complexity: O(N) worst.
         implicit none
         integer(INTL):: nsibl                       !out: number of children
         class(tree_vertex_t), intent(in):: this     !in: tree vertex
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc
         class(tree_vertex_t), pointer:: tvp,ftvp

         errc=GFC_SUCCESS; nsibl=0
         tvp=>this%next_sibling; ftvp=>tvp%prev_sibling
         do while(.not.associated(tvp,ftvp))
          nsibl=nsibl+1; tvp=>tvp%next_sibling
         enddo
         if(present(ierr)) ierr=errc
         return
        end function TreeVertexNumSiblings
!--------------------------------------------------------
        function TreeVertexFirstSibling(this) result(res)
!Returns GFC_TRUE if the tree vertex is the first vertex in the sibling list.
         implicit none
         integer(INTD):: res                             !out: result
         class(tree_vertex_t), target, intent(in):: this !in: tree vertex

         res=GFC_FALSE
         if(associated(this%parent)) then
          if(associated(this%parent%first_child,this)) res=GFC_TRUE
         else !root vertex (only one)
          res=GFC_TRUE
         endif
         return
        end function TreeVertexFirstSibling
!-------------------------------------------------------
        function TreeVertexLastSibling(this) result(res)
!Returns GFC_TRUE if the tree vertex is the last vertex in the sibling list.
         implicit none
         integer(INTD):: res                     !out: result
         class(tree_vertex_t), intent(in):: this !in: tree vertex

         res=GFC_FALSE
         if(associated(this%parent)) then
          if(associated(this%parent%first_child,this%next_sibling)) res=GFC_TRUE
         else !root vertex (only one)
          res=GFC_TRUE
         endif
         return
        end function TreeVertexLastSibling
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
          ierr=this%reset()
         class default
          ierr=GFC_INVALID_ARGS
         end select
         return
        end function TreeIterInit
!------------------------------------------------
        function TreeIterReset(this) result(ierr)
!Resets the iterator to the beginning (root element).
         implicit none
         integer(INTD):: ierr                     !out: error code (0:success)
         class(tree_iter_t), intent(inout):: this !inout: iterator

         ierr=GFC_SUCCESS
         if(associated(this%container)) then
          this%current=>this%container%root
          if(associated(this%current)) then
           ierr=this%set_status(GFC_IT_ACTIVE) !non-empty iterator
          else
           ierr=this%set_status(GFC_IT_EMPTY) !empty iterator
          endif
         else
          ierr=this%set_status(GFC_IT_NULL)
          ierr=GFC_IT_NULL
         endif
         return
        end function TreeIterReset
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
!If <elem_p> is absent, the iterator moves to the next element, if any.
!If <elem_p> is present, the iterator simply returns the next element in <elem_p> without moving.
!Complexity: O(1)...O(N).
         implicit none
         integer(INTD):: ierr                                            !out: error code (0:success)
         class(tree_iter_t), intent(inout):: this                        !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element
         class(tree_vertex_t), pointer:: tvp

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           tvp=>this%current%first_child
           if(.not.associated(tvp)) then
            tvp=>this%current
            do while(associated(tvp))
             if(tvp%last_sibling().eq.GFC_FALSE) then !not the last sibling
              tvp=>tvp%next_sibling; exit
             else
              tvp=>tvp%parent
             endif
            enddo
           endif
           if(present(elem_p)) then
            elem_p=>tvp
           else
            this%current=>tvp
            if(.not.associated(tvp)) ierr=this%set_status(GFC_IT_DONE)
           endif
           if(.not.associated(tvp)) ierr=GFC_IT_DONE
           tvp=>NULL()
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function TreeIterNext
!----------------------------------------------------------
        function TreeIterPrevious(this,elem_p) result(ierr)
!If <elem_p> is absent, the iterator moves to the previous element, if any.
!If <elem_p> is present, the iterator simply returns the previous element in <elem_p> without moving.
!Complexity: O(1)..O(N).
         implicit none
         integer(INTD):: ierr                                            !out: error code (0:success)
         class(tree_iter_t), intent(inout):: this                        !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element
         class(tree_vertex_t), pointer:: tvp

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           tvp=>this%current
           if(tvp%first_sibling().eq.GFC_TRUE) then
            tvp=>tvp%parent
           else
            tvp=>tvp%prev_sibling
            do while(associated(tvp%first_child))
             tvp=>tvp%first_child%prev_sibling !last sibling (because of ring)
            enddo
           endif
           if(present(elem_p)) then
            elem_p=>tvp
           else
            this%current=>tvp
            if(.not.associated(tvp)) ierr=this%set_status(GFC_IT_DONE)
           endif
           if(.not.associated(tvp)) ierr=GFC_IT_DONE
           tvp=>NULL()
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function TreeIterPrevious

       end module tree
