!Generic Fortran Containers (GFC): Tree
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2016-03-16 (started 2016-02-17)

!Copyright (C) 2014-2016 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2016 Oak Ridge National Laboratory (UT-Battelle)

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
! # A tree is a linked derivative of an abstract (unlinked) GFC container.
!   A subtree is a tree object that is linked as a part of a larger tree,
!   thus having its root element linked to other elements of the larger tree.
! # All accesses, updates, scans, and actions on a tree are performed via
!   a tree iterator associated with the tree. When attaching a tree
!   to another tree, the attached tree elements can be accessed either
!   via its own iterator or via the combined tree iterator. Multiple
!   iterators can be associated with a tree at a time.

       module tree
        use gfc_base
        use timers
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !output device
        logical, private:: VERBOSE=.true.   !verbositiy for errors
        integer(INTD), private:: DEBUG=0    !debugging level (0:none)
 !Tree iterator directions:
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
          procedure, non_overridable, public:: is_root=>TreeVertexIsRoot             !returns GFC_TRUE if the vertex is the root of the tree, GFC_FALSE otherwise
          procedure, non_overridable, public:: is_leaf=>TreeVertexIsLeaf             !returns GFC_TRUE if the vertex is a leaf, GFC_FALSE otherwise
        end type tree_vertex_t
 !Tree (all operations on the tree are performend via an iterator):
        type, extends(gfc_container_t), public:: tree_t
         class(tree_vertex_t), pointer, private:: root=>NULL() !root element (beginning)
         contains
          procedure, public:: is_subtree=>TreeIsSubtree !returns TRUE if the tree is a subtree of a larger tree, FALSE otherwise
        end type tree_t
 !Tree iterator:
        type, extends(gfc_iter_t), public:: tree_iter_t
         class(tree_vertex_t), pointer, private:: current=>NULL()   !currently pointed element of the container
         class(tree_t), pointer, private:: container=>NULL()        !container
         contains
          procedure, public:: init=>TreeIterInit                    !associates the iterator with a container and sets its position to the root element
          procedure, public:: reset=>TreeIterReset                  !resets the iterator to the beginning of the container (root element)
          procedure, public:: release=>TreeIterRelease              !dissociates the iterator from its container
          procedure, public:: pointee=>TreeIterPointee              !returns a pointer to the container element currently in focus
          procedure, public:: next=>TreeIterNext                    !moves the iterator to the next element
          procedure, public:: previous=>TreeIterPrevious            !moves the iterator to the previous element
          procedure, public:: to_sibling=>TreeIterToSibling         !moves the iterator to one of the siblings
          procedure, public:: to_child=>TreeIterToChild             !moves the iterator to one of the children
          procedure, public:: to_parent=>TreeIterToParent           !moves the iterator to the parent
          procedure, public:: my_parent=>TreeIterMyParent           !returns the parent of the current vertex
          procedure, public:: add_leaf=>TreeIterAddLeaf             !adds a new leaf element to the element of the container currently pointed to
          procedure, public:: delete_leaf=>TreeIterDeleteLeaf       !deletes the leaf pointed to by the iterator (if it is actually a leaf)
          procedure, public:: attach_subtree=>TreeIterAttachSubtree !attaches a subtree to the element of the container currently pointed to as the last child
          procedure, public:: detach_subtree=>TreeIterDetachSubtree !detaches a subtree beginning from the currently pointed element of the container
          procedure, public:: delete_subtree=>TreeIterDeleteSubtree !deletes a subtree beginning from the currently pointed element of the container
        end type tree_iter_t
!GLOBAL DATA:
!VISIBILITY:
 !Procedures:
        private TreeVertexNumChildren
        private TreeVertexNumSiblings
        private TreeVertexFirstSibling
        private TreeVertexLastSibling
        private TreeVertexIsRoot
        private TreeVertexIsLeaf
        private TreeIsSubtree
        private TreeIterInit
        private TreeIterReset
        private TreeIterRelease
        private TreeIterPointee
        private TreeIterNext
        private TreeIterPrevious
        private TreeIterToSibling
        private TreeIterToChild
        private TreeIterToParent
        private TreeIterMyParent
        private TreeIterAddLeaf
        private TreeIterDeleteLeaf
        private TreeIterAttachSubtree
        private TreeIterDetachSubtree
        private TreeIterDeleteSubtree

       contains
!IMPLEMENTATION:
!---------------------------------------------------------------
        function TreeVertexNumChildren(this,ierr) result(nchild)
!Returns the total number of children attached to the tree vertex.
!Complexity: O(1).
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
!Complexity: O(1).
         implicit none
         integer(INTL):: nsibl                       !out: number of siblings
         class(tree_vertex_t), intent(in):: this     !in: tree vertex
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc

         errc=GFC_SUCCESS
         if(associated(this%parent)) then
          nsibl=this%parent%num_child-1
         else !root vertex cannot have siblings
          nsibl=0
         endif
         if(nsibl.lt.0) errc=GFC_CORRUPTED_CONT
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
!--------------------------------------------------
        function TreeVertexIsRoot(this) result(res)
!Returns GFC_TRUE if the vertex is the root of the tree.
         implicit none
         integer(INTD):: res                     !out: result {GFC_TRUE,GFC_FALSE}
         class(tree_vertex_t), intent(in):: this !in: tree vertex

         if(associated(this%parent)) then
          res=GFC_FALSE
         else
          res=GFC_TRUE
         endif
         return
        end function TreeVertexIsRoot
!--------------------------------------------------
        function TreeVertexIsLeaf(this) result(res)
!Returns GFC_TRUE if the vertex is a leaf.
         implicit none
         integer(INTD):: res                     !out: result
         class(tree_vertex_t), intent(in):: this !in: tree vertex

         res=GFC_FALSE
         if(.not.associated(this%first_child)) res=GFC_TRUE
         return
        end function TreeVertexIsLeaf
!-----------------------------------------------
        function TreeIsSubtree(this) result(res)
!Returns TRUE if the tree is a subtree of a larger tree.
         implicit none
         logical:: res                    !out: result
         class(tree_t), intent(in):: this !in: tree

         res=.false.
         if(associated(this%root)) then
          if(associated(this%root%parent)) res=.true.
         endif
         return
        end function TreeIsSubtree
!----------------------------------------------------
        function TreeIterInit(this,cont) result(ierr)
!Initializes an iterator and resets it to the beginning of the container.
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
           ierr=this%set_status_(GFC_IT_ACTIVE) !non-empty iterator/container
          else
           ierr=this%set_status_(GFC_IT_EMPTY) !empty iterator/container
          endif
          call this%reset_count() !reset all iteration counters
         else
          ierr=this%set_status_(GFC_IT_NULL)
          ierr=GFC_IT_NULL
         endif
         return
        end function TreeIterReset
!--------------------------------------------------
        function TreeIterRelease(this) result(ierr)
!Dissociates the iterator from its container.
         implicit none
         integer(INTD):: ierr                     !out: error code (0:success)
         class(tree_iter_t), intent(inout):: this !inout: iterator

         this%current=>NULL(); this%container=>NULL()
         call this%reset_count(); ierr=this%set_status_(GFC_IT_NULL)
         return
        end function TreeIterRelease
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
!Complexity: O(1)...O(N). No additional memory is used.
         implicit none
         integer(INTD):: ierr                                            !out: error code (0:success)
         class(tree_iter_t), intent(inout):: this                        !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element
         class(tree_vertex_t), pointer:: tvp

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           ierr=GFC_SUCCESS
           tvp=>this%current%first_child
           if(.not.associated(tvp)) then
            tvp=>this%current
            do while(associated(tvp))
             if(.not.associated(tvp,this%container%root)) then !root of a subtree may have a parent
              if(tvp%last_sibling().eq.GFC_FALSE) then !not the last sibling
               tvp=>tvp%next_sibling; exit
              else
               tvp=>tvp%parent
              endif
             else
              tvp=>NULL()
             endif
            enddo
           endif
           if(present(elem_p)) then
            elem_p=>tvp
           else
            this%current=>tvp
            if(.not.associated(tvp)) ierr=this%set_status_(GFC_IT_DONE)
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
!Complexity: O(1)..O(N). No additional memory is used.
         implicit none
         integer(INTD):: ierr                                            !out: error code (0:success)
         class(tree_iter_t), intent(inout):: this                        !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element
         class(tree_vertex_t), pointer:: tvp
         logical:: on_root

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           ierr=GFC_SUCCESS
           tvp=>this%current; on_root=associated(tvp,this%container%root)
           if((.not.on_root).and.tvp%first_sibling().eq.GFC_TRUE) then
            tvp=>tvp%parent
           else
            if(.not.on_root) tvp=>tvp%prev_sibling
            if((.not.on_root).or.associated(tvp%first_child)) then
             do while(associated(tvp%first_child))
              tvp=>tvp%first_child%prev_sibling !last sibling among the children (because of ring linking)
             enddo
            else
             tvp=>NULL()
            endif
           endif
           if(present(elem_p)) then
            elem_p=>tvp
           else
            this%current=>tvp
            if(.not.associated(tvp)) ierr=this%set_status_(GFC_IT_DONE)
           endif
           if(.not.associated(tvp)) ierr=GFC_IT_DONE
           tvp=>NULL()
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function TreeIterPrevious
!----------------------------------------------------------------
        function TreeIterToSibling(this,to_previous) result(ierr)
!Moves the iterator either to the next or to the previous sibling.
         implicit none
         integer(INTD):: ierr                        !out: error code (0:success)
         class(tree_iter_t), intent(inout):: this    !inout: iterator
         logical, intent(in), optional:: to_previous !in: if TRUE, the iterator will move to the previous sibling (defaults to FALSE)
         class(tree_vertex_t), pointer:: tvp
         logical:: to_next

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           if(.not.associated(this%current,this%container%root)) then
            if(present(to_previous)) then; to_next=.not.to_previous; else; to_next=.true.; endif
            ierr=GFC_SUCCESS
            if(to_next) then
             tvp=>this%current%next_sibling
             if(associated(tvp,this%current%parent%first_child)) then
              ierr=GFC_NO_MOVE
             else
              this%current=>tvp
             endif
             tvp=>NULL()
            else
             if(associated(this%current,this%current%parent%first_child)) then
              ierr=GFC_NO_MOVE
             else
              this%current=>this%current%prev_sibling
             endif
            endif
           else
            ierr=GFC_NO_MOVE !tree/subtree root does not have siblings within its iterator
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function TreeIterToSibling
!--------------------------------------------------
        function TreeIterToChild(this) result(ierr)
!Moves the iterator to the first child, if any.
         implicit none
         integer(INTD):: ierr                        !out: error code (0:success)
         class(tree_iter_t), intent(inout):: this    !inout: iterator

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           if(associated(this%current%first_child)) then
            this%current=>this%current%first_child
            ierr=GFC_SUCCESS
           else
            ierr=GFC_NO_MOVE
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function TreeIterToChild
!---------------------------------------------------
        function TreeIterToParent(this) result(ierr)
!Moves the iterator to the parent, if any.
         implicit none
         integer(INTD):: ierr                        !out: error code (0:success)
         class(tree_iter_t), intent(inout):: this    !inout: iterator

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           if(associated(this%current,this%container%root)) then
            ierr=GFC_NO_MOVE
           else
            this%current=>this%current%parent
            ierr=GFC_SUCCESS
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function TreeIterToParent
!----------------------------------------------------------
        function TreeIterMyParent(this,ierr) result(parent)
!Returns the parent of the current vertex.
         implicit none
         class(gfc_cont_elem_t), pointer:: parent    !out: parent element of the current vertex
         class(tree_iter_t), intent(in):: this       !in: iterator
         integer(INTD), intent(out), optional:: ierr !out: error code (0:success)
         integer(INTD):: errc

         errc=this%get_status()
         if(errc.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           parent=>this%current%parent; errc=GFC_SUCCESS
          else
           parent=>NULL(); errc=GFC_CORRUPTED_CONT
          endif
         else
          parent=>NULL()
         endif
         if(present(ierr)) ierr=errc
         return
        end function TreeIterMyParent
!--------------------------------------------------------------------------------------------
#ifdef NO_GNU
        function TreeIterAddLeaf(this,elem_val,assoc_only,no_move,copy_constr_f) result(ierr) !`GCC/5.3.0 has a bug with this
#else
        function TreeIterAddLeaf(this,elem_val,assoc_only,no_move) result(ierr)
#endif
!Creates a new container element (leaf) as the last child of the currently pointed
!element and stores the value <elem_val> in it, either by value or by reference.
         implicit none
         integer(INTD):: ierr                       !out: error code (0:success)
         class(tree_iter_t), intent(inout):: this   !inout: iterator
         class(*), target, intent(in):: elem_val    !in: value to store in the container
         logical, intent(in), optional:: assoc_only !in: TRUE: store by reference, FALSE: store by value (defaults to FALSE)
         logical, intent(in), optional:: no_move    !in: if TRUE, the iterator will not move to the newly added element (defaults to FALSE)
#ifdef NO_GNU
         procedure(gfc_copy_i), optional:: copy_constr_f !in: user-defined generic copy constructor
#endif
         class(tree_vertex_t), pointer:: tvp
         integer:: errc
         logical:: assoc,nomo
         integer(INTL):: nelems

         if(present(assoc_only)) then; assoc=assoc_only; else; assoc=.false.; endif
         if(present(no_move)) then; nomo=no_move; else; nomo=.false.; endif
         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          ierr=GFC_SUCCESS
          if(associated(this%container).and.associated(this%current)) then
           if(associated(this%current%first_child)) then
            tvp=>this%current%first_child%prev_sibling !last sibling among children
            allocate(tvp%next_sibling,STAT=errc)
            if(errc.eq.0) then
#ifdef NO_GNU
             if(present(copy_constr_f)) then
              call tvp%next_sibling%construct(elem_val,ierr,assoc_only=assoc,copy_constr_func=copy_constr_f)
             else
              call tvp%next_sibling%construct(elem_val,ierr,assoc_only=assoc)
             endif
#else
             call tvp%next_sibling%construct(elem_val,ierr,assoc_only=assoc)
#endif
             if(ierr.eq.GFC_SUCCESS) then
              tvp%next_sibling%prev_sibling=>tvp
              tvp=>tvp%next_sibling
             else
              deallocate(tvp%next_sibling)
             endif
            else
             tvp%next_sibling=>this%current%first_child
             ierr=GFC_MEM_ALLOC_FAILED
            endif
           else
            allocate(this%current%first_child,STAT=errc)
            if(errc.eq.0) then
#ifdef NO_GNU
             if(present(copy_constr_f)) then
              call this%current%first_child%construct(elem_val,ierr,assoc_only=assoc,copy_constr_func=copy_constr_f)
             else
              call this%current%first_child%construct(elem_val,ierr,assoc_only=assoc)
             endif
#else
             call this%current%first_child%construct(elem_val,ierr,assoc_only=assoc)
#endif
             if(ierr.eq.GFC_SUCCESS) then
              tvp=>this%current%first_child
             else
              deallocate(this%current%first_child)
             endif
            else
             this%current%first_child=>NULL()
             ierr=GFC_MEM_ALLOC_FAILED
            endif
           endif
           if(ierr.eq.GFC_SUCCESS) then
            tvp%parent=>this%current
            tvp%next_sibling=>this%current%first_child
            this%current%first_child%prev_sibling=>tvp
            this%current%num_child=this%current%num_child+1
            nelems=this%container%update_num_elems_(1_INTL,ierr); if(ierr.ne.GFC_SUCCESS) ierr=GFC_CORRUPTED_CONT
            if(.not.nomo) this%current=>tvp
           endif
           tvp=>NULL()
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         elseif(ierr.eq.GFC_IT_EMPTY) then !very first element of the container
          ierr=GFC_SUCCESS
          if(associated(this%container)) then
           if(.not.(associated(this%container%root).or.associated(this%current))) then
            allocate(this%container%root,STAT=errc)
            if(errc.eq.0) then
#ifdef NO_GNU
             if(present(copy_constr_f)) then
              call this%container%root%construct(elem_val,ierr,assoc_only=assoc,copy_constr_func=copy_constr_f)
             else
              call this%container%root%construct(elem_val,ierr,assoc_only=assoc)
             endif
#else
             call this%container%root%construct(elem_val,ierr,assoc_only=assoc)
#endif
             if(ierr.eq.GFC_SUCCESS) then
              ierr=this%reset() !move to the just added first element regardless of <no_move> and change the EMPTY status to ACTIVE
              nelems=this%container%update_num_elems_(1_INTL,ierr); if(ierr.ne.GFC_SUCCESS) ierr=GFC_CORRUPTED_CONT
             else
              deallocate(this%container%root)
              this%container%root=>NULL()
             endif
            else
             this%container%root=>NULL()
             ierr=GFC_MEM_ALLOC_FAILED
            endif
           else
            ierr=GFC_CORRUPTED_CONT
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function TreeIterAddLeaf
!-------------------------------------------------------------------
        function TreeIterDeleteLeaf(this,destruct_func) result(ierr)
!Deletes a leaf from a tree and moves the iterator to its parent.
         implicit none
         integer(INTD):: ierr                                !out: error code (0:success)
         class(tree_iter_t), intent(inout):: this            !inout: iterator
         procedure(gfc_destruct_i), optional:: destruct_func !in: value destructor
         integer(INTL):: totelems
         integer(INTD):: errc
         class(tree_vertex_t), pointer:: tvp

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           ierr=GFC_SUCCESS; tvp=>this%current
           if(tvp%is_leaf().eq.GFC_TRUE) then
            if(present(destruct_func)) then !destructs the value of the current element
             call tvp%destruct(errc,destruct_func)
            else
             call tvp%destruct(errc)
            endif
            if(errc.ne.GFC_SUCCESS) ierr=NOT_CLEAN
            if(tvp%num_siblings(errc).gt.0) then
             if(errc.eq.GFC_SUCCESS) then
              if(associated(tvp%parent)) then
               if(tvp%first_sibling().eq.GFC_TRUE) tvp%parent%first_child=>tvp%next_sibling
              endif
              tvp%prev_sibling%next_sibling=>tvp%next_sibling
              tvp%next_sibling%prev_sibling=>tvp%prev_sibling
             else
              ierr=GFC_CORRUPTED_CONT
             endif
            else
             if(errc.eq.GFC_SUCCESS) then
              if(associated(tvp%parent)) tvp%parent%first_child=>NULL()
             else
              ierr=GFC_CORRUPTED_CONT
             endif
            endif
            if(ierr.eq.GFC_SUCCESS.or.ierr.eq.NOT_CLEAN) then
             totelems=this%container%update_num_elems_(-1_INTL,errc)
             if(errc.eq.GFC_SUCCESS) then
              if(associated(tvp,this%container%root)) then
               this%current=>NULL(); this%container%root=>NULL()
               errc=this%set_status_(GFC_IT_EMPTY); if(errc.ne.GFC_SUCCESS) ierr=GFC_CORRUPTED_CONT
              else
               tvp%parent%num_child=tvp%parent%num_child-1
               this%current=>tvp%parent
              endif
              deallocate(tvp,STAT=errc); if(errc.ne.0) ierr=NOT_CLEAN
             else
              ierr=GFC_CORRUPTED_CONT
             endif
            endif
           else
            ierr=GFC_INVALID_ARGS !the current element is not a leaf
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function TreeIterDeleteLeaf
!----------------------------------------------------------------
        function TreeIterAttachSubtree(this,subtree) result(ierr)
!Attaches a tree as the last child to the current iterator position,
!thus making the attached tree a subtree (root will have a parent/siblings).
!The attached tree root shall not have a parent/siblings at the beginning.
!The iterator position does not change.
         implicit none
         integer(INTD):: ierr                     !out: error code (0:success)
         class(tree_iter_t), intent(inout):: this !inout: iterator
         class(tree_t), intent(inout):: subtree   !inout: tree at input, subtree at output
         class(tree_vertex_t), pointer:: tvp
         integer(INTL):: nelems,totelems

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           if(associated(subtree%root)) then
            ierr=GFC_SUCCESS
!           nelems=subtree%num_elems_(ierr) !this is irrelevant
            if(.not.associated(subtree%root%parent)) then !the attached tree cannot be a subtree prior to that
             if(associated(this%current%first_child)) then
              tvp=>this%current%first_child%prev_sibling !tvp => last sibling
              tvp%next_sibling=>subtree%root
              subtree%root%prev_sibling=>tvp
              subtree%root%next_sibling=>this%current%first_child
              this%current%first_child%prev_sibling=>subtree%root
             else
              this%current%first_child=>subtree%root
              subtree%root%next_sibling=>subtree%root
              subtree%root%prev_sibling=>subtree%root
             endif
             subtree%root%parent=>this%current
             this%current%num_child=this%current%num_child+1
!            totelems=this%container%update_num_elems_(nelems,ierr)
!            if(ierr.ne.GFC_SUCCESS) ierr=GFC_CORRUPTED_CONT
             call this%container%quick_counting_off_() !turn off quick counting in the combined container
             call subtree%quick_counting_off_() !turn off quick counting in the subcontainer
            else
             ierr=GFC_INVALID_ARGS
            endif
           else
            ierr=GFC_INVALID_ARGS
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function TreeIterAttachSubtree
!----------------------------------------------------------------
        function TreeIterDetachSubtree(this,subtree) result(ierr)
!Detaches a subtree beginning at the current iterator position and returns it as a tree.
!The iterator is moved to the parental vertex.
         implicit none
         integer(INTD):: ierr                     !out: error code (0:success)
         class(tree_iter_t), intent(inout):: this !inout: iterator
         class(tree_t), intent(inout):: subtree   !inout: subtree (must be empty at entrance)
         class(tree_vertex_t), pointer:: psib,nsib
         type(tree_iter_t):: subtree_it
         integer(INTL):: nelems,totelems

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           if(.not.associated(subtree%root)) then !subtree must be empty on entrance
            ierr=GFC_SUCCESS
            psib=>this%current%prev_sibling; nsib=>this%current%next_sibling
            subtree%root=>this%current
            if(associated(this%current%parent)) then
             this%current%parent%num_child=this%current%parent%num_child-1
             if(this%current%num_siblings(ierr).gt.0) then
              if(ierr.eq.GFC_SUCCESS) this%current%parent%first_child=>this%current%next_sibling
             else
              if(ierr.eq.GFC_SUCCESS) this%current%parent%first_child=>NULL()
             endif
            endif
            if(ierr.eq.GFC_SUCCESS) then
             if(associated(this%current,this%container%root)) then
              this%current=>NULL(); this%container%root=>NULL()
              ierr=this%set_status_(GFC_IT_EMPTY); if(ierr.ne.GFC_SUCCESS) ierr=GFC_CORRUPTED_CONT
             else
              this%current=>this%current%parent
             endif
             if(ierr.eq.GFC_SUCCESS) then
              subtree%root%parent=>NULL()
              psib%next_sibling=>nsib; nsib%prev_sibling=>psib
              subtree%root%prev_sibling=>subtree%root
              subtree%root%next_sibling=>subtree%root
              if(this%container%num_elems_().gt.0) then !quick counting is still on
               ierr=subtree_it%init(subtree)
               if(ierr.eq.GFC_SUCCESS) then
                ierr=subtree_it%scan()
                if(ierr.eq.GFC_SUCCESS) then
                 nelems=subtree_it%total_count()
                 totelems=subtree%update_num_elems_(nelems)
                 totelems=this%container%update_num_elems_(-nelems,ierr)
                else
                 call this%container%quick_counting_off_()
                 call subtree%quick_counting_off_()
                endif
               else
                call this%container%quick_counting_off_()
                call subtree%quick_counting_off_()
               endif
               ierr=subtree_it%release()
              else
               call subtree%quick_counting_off_()
              endif
             endif
            else
             ierr=GFC_CORRUPTED_CONT
            endif
            psib=>NULL(); nsib=>NULL()
           else
            ierr=GFC_INVALID_ARGS
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function TreeIterDetachSubtree
!----------------------------------------------------------------------
        function TreeIterDeleteSubtree(this,destruct_func) result(ierr)
!Completely deletes a subtree starting from the current iterator position.
!The iterator is moved to the parent at the end. A return status NOT_CLEAN
!indicates that some memory deallocation and/or object destruction failed.
         implicit none
         integer(INTD):: ierr                                !out: error code (0:success)
         class(tree_iter_t), intent(inout):: this            !inout: iterator
         procedure(gfc_destruct_i), optional:: destruct_func !in: value destruction function
         class(tree_vertex_t), pointer:: tvp
         logical:: subtree,dsf,ntcl

         ierr=this%get_status(); ntcl=.false.
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current)) then
           dsf=.false.; if(present(destruct_func)) dsf=.true.
           if(associated(this%current%parent)) then
            tvp=>this%current%parent; subtree=.true.
           else
            tvp=>NULL(); subtree=.false.
           endif
           do while(associated(this%current))
            do while(associated(this%current%first_child)) !find a leaf
             this%current=>this%current%first_child
            enddo
            if(dsf) then
             ierr=this%delete_leaf(destruct_func)
            else
             ierr=this%delete_leaf()
            endif
            if(ierr.eq.NOT_CLEAN) then; ntcl=.true.; ierr=GFC_SUCCESS; endif
            if(ierr.ne.GFC_SUCCESS) exit
            if(subtree) then; if(associated(this%current,tvp)) exit; endif
           enddo
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         if(ntcl.and.ierr.eq.GFC_SUCCESS) ierr=NOT_CLEAN
         return
        end function TreeIterDeleteSubtree

       end module tree
!=========================================
!TESTING:
       module tree_test
        use gfc_base
        use tree
        use timers, only: thread_wtime
        implicit none
        private

        public dil_test_tree

        type, private:: some_t
         real(8):: some_real=0d0
         integer(INTD):: some_int=0
         integer(INTL):: some_long=0
         real(8), pointer:: some_arr(:)=>NULL()
        end type some_t

        type, private:: vertex_ptr_t
         type(gfc_cont_elem_t), pointer:: ptr=>NULL()
        end type vertex_ptr_t

       contains

        function some_destructor(obj) result(ierr)
         implicit none
         integer(INTD):: ierr
         class(*), intent(inout):: obj
         ierr=1
         select type(obj)
         class is(some_t)
          if(associated(obj%some_arr)) deallocate(obj%some_arr)
          obj%some_real=0d0
          obj%some_long=0
          ierr=0
         end select
         return
        end function some_destructor

        function some_action(obj) result(ierr)
         implicit none
         integer(INTD):: ierr
         class(*), intent(inout):: obj
         ierr=1
         select type(obj)
         class is(some_t)
          if(.not.associated(obj%some_arr)) allocate(obj%some_arr(16))
          obj%some_long=size(obj%some_arr)
          ierr=0
         end select
         return
        end function some_action

        function print_action(obj) result(ierr)
         implicit none
         integer(INTD):: ierr
         class(*), intent(inout):: obj
         ierr=1
         select type(obj)
         class is(some_t)
          write(*,'("vertex ",i9," (",F6.4,")")') obj%some_int,obj%some_real
          ierr=0
         end select
         return
        end function print_action

        function some_predicate(obj) result(pred)
         implicit none
         integer(INTD):: pred
         class(*), intent(in):: obj
         pred=GFC_FALSE
         select type(obj)
         class is(some_t)
          if(obj%some_real.gt.5d-1) pred=GFC_TRUE
         end select
         return
        end function some_predicate

        function dil_test_tree(perf,dev_out) result(ierr)
         implicit none
         integer(INTD):: ierr
         real(8), intent(out):: perf
         integer(INTD), intent(in), optional:: dev_out
         integer(INTD), parameter:: MAX_TREE_ELEMS=1000000
         class(gfc_cont_elem_t), pointer:: tvp
         class(*), pointer:: val_p
         integer(INTD):: jo,i,j,m
         type(some_t):: some_val
         type(tree_t):: some_tree
         type(tree_iter_t):: some_iter
         real(8):: tms,tm

         if(present(dev_out)) then; jo=dev_out; else; jo=6; endif
         perf=0d0; tms=thread_wtime()
         ierr=some_iter%init(some_tree); if(ierr.ne.GFC_SUCCESS) then; ierr=1; return; endif
!Add elements to the tree:
         some_val%some_real=1d0; some_val%some_int=1 !root
         ierr=some_iter%add_leaf(some_val); if(ierr.ne.GFC_SUCCESS) then; ierr=2; return; endif
         i=MAX_TREE_ELEMS-1
         do while(i.gt.0)
          ierr=some_iter%reset(); if(ierr.ne.GFC_SUCCESS) then; ierr=3; return; endif
          do while(i.gt.0)
           ierr=some_iter%scan(.TRUE.,some_predicate,some_action); if(ierr.ne.GFC_IT_ACTIVE) exit
           call random_number(some_val%some_real); some_val%some_int=MAX_TREE_ELEMS-i+1
!          tvp=>some_iter%pointee(); val_p=>tvp%get_value(); ierr=print_action(val_p) !debug: parent
           ierr=some_iter%add_leaf(some_val); if(ierr.ne.GFC_SUCCESS) then; ierr=4; return; endif
!          tvp=>some_iter%pointee(); val_p=>tvp%get_value(); ierr=print_action(val_p) !debug: newly added
           i=i-1
          enddo
          if(ierr.ne.GFC_SUCCESS.and.ierr.ne.GFC_IT_DONE) then; ierr=5; return; endif
         enddo
!        write(jo,'("Total number of elements in the tree = ",i9)') some_tree%num_elems_(ierr) !debug
         if(ierr.ne.GFC_SUCCESS) then; ierr=6; return; endif
         ierr=some_iter%reset(); if(ierr.ne.GFC_SUCCESS) then; ierr=7; return; endif
         ierr=some_iter%scan(); if(ierr.ne.GFC_IT_DONE) then; ierr=8; return; endif
!        write(jo,'("Total number of traversed elements   = ",i9)') some_iter%total_count() !debug
!        if(some_tree%num_elems_().ne.some_iter%total_count()) then; ierr=9; return; endif
!Delete the tree:
         ierr=some_iter%reset(); if(ierr.ne.GFC_SUCCESS) then; ierr=10; return; endif
         ierr=some_iter%delete_subtree(some_destructor); if(ierr.ne.GFC_SUCCESS) then; ierr=11; return; endif
         tm=thread_wtime(tms); perf=dble(MAX_TREE_ELEMS)/tm
         return
        end function dil_test_tree

       end module tree_test
