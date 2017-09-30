!Generic Fortran Containers (GFC): Vector tree, combines vector and tree:
!The elements are initially inserted in a vector with an option to be
!later added in a tree, thus imposing a tree relationship on them.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2017/09/30

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

!DESCRIPTION:
! # Vector tree is a generic vector with an additional tree structure imposed on it.
!   Initially the container elements are appended into the vector. During this process
!   or after it, each existing vector element can also be added to the associated tree.
!   Each existing vector element can only be added once to the tree.
!   The vector-tree iterator can move to the next/previous element in the vector, as well
!   as perform tree-like moves (to child, to parent, to cousin, etc.). If the next/previous
!   element of the vector has not been yet added to the tree, the corresponding next/previous
!   move of the iterator will fail. Deletion of an internal vector element, that is, not the
!   last one, will require renumeration of the vector elements, with the O(N) worst cost.

       module gfc_vec_tree
        use gfc_base
        use gfc_vector
        use gfc_tree
        use timers
        implicit none
        private
        public LEFT_SIBLING,RIGHT_SIBLING
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !output device
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
        integer(INTD), private:: DEBUG=0    !debugging level (0:none)
!TYPES:
 !Tree position:
        type, private:: tree_pos_t
         class(tree_vertex_t), pointer, private:: pos=>NULL() !pointer to a tree vertex
        end type tree_pos_t
 !Vector tree:
        type, extends(gfc_container_t), public:: vec_tree_t
         type(vector_t), private:: vertices !vector of values of vertices
         type(vector_t), private:: tree_pos !vector of positions of vertices in the imposed tree
         type(tree_t), private:: tree       !tree imposed on the vector of vertices (tree vertex values are IDs of the vector vertices)
         contains
          procedure, public:: is_empty=>VecTreeIsEmpty       !returns GFC_TRUE if the container is empty, GFC_FALSE otherwise (or error code)
          procedure, public:: is_full=>VecTreeIsFull         !returns GFC_TRUE if the container is full, GFC_FALSE otherwise (or error code)
          procedure, public:: capacity=>VecTreeCapacity      !returns maximal length of the container
          procedure, public:: length=>VecTreeLength          !returns current length of the container
          procedure, public:: lower_bound=>VecTreeLowerBound !returns the lower bound of the container
          procedure, public:: upper_bound=>VecTreeUpperBound !returns the upper bound of the container
        end type vec_tree_t
 !Vector tree iterator:
        type, extends(gfc_iter_t), public:: vec_tree_iter_t
         class(vec_tree_t), pointer, private:: container=>NULL() !vector tree container
         type(vector_iter_t), private:: val_it                   !iterator for vertex values
         type(vector_iter_t), private:: pos_it                   !iterator for vertex tree positions
         type(tree_iter_t), private:: tree_it                    !iterator for the tree
         contains
          procedure, private:: update_status_=>VecTreeIterUpdateStatus  !PRIVATE: updates vector tree iterator status after component updates
          procedure, public:: init=>VecTreeIterInit                     !associates the iterator with a container and sets its position to the beginning
          procedure, public:: reset=>VecTreeIterReset                   !resets the iterator to the beginning of the container (only the vector part if tree is not fully built)
          procedure, public:: reset_back=>VecTreeIterResetBack          !resets the iterator to the end of the container (only the vector part if tree is not fully built)
          procedure, public:: release=>VecTreeIterRelease               !releases the iterator (dissociates it from its container)
          procedure, public:: pointee=>VecTreeIterPointee               !returns a pointer to the container element currently pointed to
          procedure, public:: next=>VecTreeIterNext                     !moves the iterator to the next container element (only the vector part if tree is not fully built)
          procedure, public:: previous=>VecTreeIterPrevious             !moves the iterator to the previous container element (only the vector part if tree is not fully built)
          procedure, public:: get_length=>VecTreeIterGetLength          !returns the current length of the container
          procedure, public:: get_offset=>VecTreeIterGetOffset          !returns the current offset in the container
          procedure, public:: get_level=>VecTreeIterGetLevel            !returns the distance from the tree root for the current iterator position
          procedure, public:: get_root_id=>VecTreeIterGetRootId         !returns the offset of the tree root node
          procedure, public:: get_num_children=>VecTreeIterGetNumChildren !returns the total number of children in the current (tree) iterator position
          procedure, public:: get_num_siblings=>VecTreeIterGetNumSiblings !returns the total number of siblings in the current (tree) iterator position
          procedure, public:: on_first_sibling=>VecTreeIterOnFirstSibling !returns GFC_TRUE if tree is positioned on the first sibling
          procedure, public:: on_last_sibling=>VecTreeIterOnLastSibling   !returns GFC_TRUE if tree is positioned on the last sibling
          procedure, public:: element_value=>VecTreeIterElementValue    !returns a pointer to the value of a specific element
          procedure, public:: move_to=>VecTreeIterMoveTo                !moves the iterator to the specific container element by its sequential number
          procedure, public:: move_to_sibling=>VecTreeIterMoveToSibling !moves the iterator to the next/previous sibling of the current element in the tree
          procedure, public:: move_to_child=>VecTreeIterMoveToChild     !moves the iterator to the first child of the current element in the tree
          procedure, public:: move_to_parent=>VecTreeIterMoveToParent   !moves the iterator to the parent of the current element in the tree
          procedure, public:: move_up=>VecTreeIterMoveUp                !moves the iterator towrards the tree root a specific number of hops
          procedure, public:: move_to_cousin=>VecTreeIterMoveToCousin   !moves the iterator to the next/previous cousin of the current element in the tree (same tree level)
          procedure, public:: append=>VecTreeIterAppend                 !appends a new element at the end of the container
          procedure, public:: add_leaf=>VecTreeIterAddLeaf              !appends a new child element to the tree vertex currently pointed to
          procedure, public:: delete_all=>VecTreeIterDeleteAll          !deletes all container elements
        end type vec_tree_iter_t
!VISIBILITY:
 !vec_tree_t:
        private VecTreeIsEmpty
        private VecTreeIsFull
        private VecTreeCapacity
        private VecTreeLength
        private VecTreeLowerBound
        private VecTreeUpperBound
 !vec_tree_iter_t:
        private VecTreeIterUpdateStatus
        private VecTreeIterInit
        private VecTreeIterReset
        private VecTreeIterResetBack
        private VecTreeIterRelease
        private VecTreeIterPointee
        private VecTreeIterNext
        private VecTreeIterPrevious
        private VecTreeIterGetLength
        private VecTreeIterGetOffset
        private VecTreeIterGetLevel
        private VecTreeIterGetRootId
        private VecTreeIterGetNumChildren
        private VecTreeIterGetNumSiblings
        private VecTreeIterOnFirstSibling
        private VecTreeIterOnLastSibling
        private VecTreeIterElementValue
        private VecTreeIterMoveTo
        private VecTreeIterMoveToSibling
        private VecTreeIterMoveToChild
        private VecTreeIterMoveToParent
        private VecTreeIterMoveUp
        private VecTreeIterMoveToCousin
        private VecTreeIterAppend
        private VecTreeIterAddLeaf
        private VecTreeIterDeleteAll
!IMPLEMENTATION:
       contains
![vec_tree_t]====================================
        function VecTreeIsEmpty(this) result(res)
!Returns GFC_TRUE if the container is empty, GFC_FALSE otherwise (or error code).
         implicit none
         integer(INTD):: res                  !out: result of the query (or error code)
         class(vec_tree_t), intent(in):: this !in: vector tree

         res=this%vertices%is_empty()
         return
        end function VecTreeIsEmpty
!-----------------------------------------------
        function VecTreeIsFull(this) result(res)
!Returns GFC_TRUE if the container is full, GFC_FALSE otherwise (or error code).
         implicit none
         integer(INTD):: res                  !out: result of the query (or error code)
         class(vec_tree_t), intent(in):: this !in: vector tree

         res=this%vertices%is_full()
         return
        end function VecTreeIsFull
!-----------------------------------------------------------
        function VecTreeCapacity(this,ierr) result(capacity)
!Returns the capacity (max length) of the container.
         implicit none
         integer(INTL):: capacity                    !out: container capacity (max length)
         class(vec_tree_t), intent(in):: this        !in: vector tree
         integer(INTD), intent(out), optional:: ierr !out: error code

         if(present(ierr)) then
          capacity=this%vertices%capacity(ierr)
         else
          capacity=this%vertices%capacity()
         endif
         return
        end function VecTreeCapacity
!-------------------------------------------------------
        function VecTreeLength(this,ierr) result(length)
!Returns the current length of the container.
         implicit none
         integer(INTL):: length                      !out: container length
         class(vec_tree_t), intent(in):: this        !in: vector tree
         integer(INTD), intent(out), optional:: ierr !out: error code

         if(present(ierr)) then
          length=this%vertices%length(ierr)
         else
          length=this%vertices%length()
         endif
         return
        end function VecTreeLength
!----------------------------------------------------------
        function VecTreeLowerBound(this,ierr) result(lower)
!Returns the current lower bound of the container.
         implicit none
         integer(INTL):: lower                       !out: container lower bound
         class(vec_tree_t), intent(in):: this        !in: vector tree
         integer(INTD), intent(out), optional:: ierr !out: error code

         if(present(ierr)) then
          lower=this%vertices%lower_bound(ierr)
         else
          lower=this%vertices%lower_bound()
         endif
         return
        end function VecTreeLowerBound
!----------------------------------------------------------
        function VecTreeUpperBound(this,ierr) result(upper)
!Returns the current upper bound of the container.
         implicit none
         integer(INTL):: upper                       !out: container upper bound
         class(vec_tree_t), intent(in):: this        !in: vector tree
         integer(INTD), intent(out), optional:: ierr !out: error code

         if(present(ierr)) then
          upper=this%vertices%upper_bound(ierr)
         else
          upper=this%vertices%upper_bound()
         endif
         return
        end function VecTreeUpperBound
![vec_tree_iter_t]==============================
        subroutine VecTreeIterUpdateStatus(this)
!Updates the iterator status after component updates.
         implicit none
         class(vec_tree_iter_t), intent(inout):: this !inout: iterator
         integer(INTD):: errc

         errc=this%set_status_(this%val_it%get_status())
         return
        end subroutine VecTreeIterUpdateStatus
!-------------------------------------------------------
        function VecTreeIterInit(this,cont) result(ierr)
!Initializes the iterator with its container and resets it to the beginning.
         implicit none
         integer(INTD):: ierr                              !out: error code
         class(vec_tree_iter_t), intent(inout):: this      !inout: iterator
         class(gfc_container_t), target, intent(in):: cont !in: container

         select type(cont)
         class is(vec_tree_t)
          this%container=>cont
          ierr=this%val_it%init(this%container%vertices)
          if(ierr.eq.GFC_SUCCESS) then
           ierr=this%pos_it%init(this%container%tree_pos)
           if(ierr.eq.GFC_SUCCESS) then
            ierr=this%tree_it%init(this%container%tree)
            if(ierr.eq.GFC_SUCCESS) ierr=this%reset()
           endif
          endif
         class default
          ierr=GFC_INVALID_ARGS
         end select
         return
        end function VecTreeIterInit
!---------------------------------------------------
        function VecTreeIterReset(this) result(ierr)
!Resets the iterator to the beginning of the container (vector part).
!If the first element of the vector part is not in the tree yet,
!the tree subiterator will be set to GFC_IT_DONE.
         implicit none
         integer(INTD):: ierr                         !out: error code
         class(vec_tree_iter_t), intent(inout):: this !inout: iterator
         class(*), pointer:: up
         class(tree_pos_t), pointer:: tpp

         ierr=this%val_it%reset()
         if(ierr.eq.GFC_SUCCESS) then
          ierr=this%pos_it%reset()
          if(ierr.eq.GFC_SUCCESS) then
           if(this%pos_it%get_status().eq.GFC_IT_ACTIVE) then
            up=>this%pos_it%get_value(ierr)
            if(ierr.eq.GFC_SUCCESS.and.associated(up)) then
             select type(up); class is(tree_pos_t); tpp=>up; class default; ierr=GFC_ERROR; end select
             if(ierr.eq.GFC_SUCCESS) call this%tree_it%jump_(tpp%pos)
            else
             ierr=GFC_CORRUPTED_CONT
            endif
           else
            ierr=this%tree_it%set_status_(GFC_IT_EMPTY)
           endif
           call this%tree_it%reset_count()
           call this%reset_count() !reset all iteration counters
          endif
         endif
         call this%update_status_()
         return
        end function VecTreeIterReset
!-------------------------------------------------------
        function VecTreeIterResetBack(this) result(ierr)
!Resets the iterator to the end of the container (vector part).
!If the last element of the vector part is not in the tree yet,
!the tree subiterator will be set to GFC_IT_DONE.
         implicit none
         integer(INTD):: ierr                         !out: error code
         class(vec_tree_iter_t), intent(inout):: this !inout: iterator
         class(*), pointer:: up
         class(tree_pos_t), pointer:: tpp

         ierr=this%val_it%reset_back()
         if(ierr.eq.GFC_SUCCESS) then
          ierr=this%pos_it%reset_back()
          if(ierr.eq.GFC_SUCCESS) then
           if(this%pos_it%get_status().eq.GFC_IT_ACTIVE) then
            up=>this%pos_it%get_value(ierr)
            if(ierr.eq.GFC_SUCCESS.and.associated(up)) then
             select type(up); class is(tree_pos_t); tpp=>up; class default; ierr=GFC_ERROR; end select
             if(ierr.eq.GFC_SUCCESS) call this%tree_it%jump_(tpp%pos)
            else
             ierr=GFC_CORRUPTED_CONT
            endif
           else
            ierr=this%tree_it%set_status_(GFC_IT_EMPTY)
           endif
           call this%tree_it%reset_count()
           call this%reset_count() !reset all iteration counters
          endif
         endif
         call this%update_status_()
         return
        end function VecTreeIterResetBack
!-----------------------------------------------------
        function VecTreeIterRelease(this) result(ierr)
!Dissociates the iterator from its container.
         implicit none
         integer(INTD):: ierr                         !out: error code
         class(vec_tree_iter_t), intent(inout):: this !inout: iterator
         integer(INTD):: errc

         call this%reset_count()
         ierr=this%tree_it%release()
         errc=this%pos_it%release(); if(ierr.eq.GFC_SUCCESS.and.errc.ne.GFC_SUCCESS) ierr=errc
         errc=this%val_it%release(); if(ierr.eq.GFC_SUCCESS.and.errc.ne.GFC_SUCCESS) ierr=errc
         this%container=>NULL()
         errc=this%set_status_(GFC_IT_NULL); if(ierr.eq.GFC_SUCCESS.and.errc.ne.GFC_SUCCESS) ierr=errc
         return
        end function VecTreeIterRelease
!-----------------------------------------------------------
        function VecTreeIterPointee(this,ierr) result(pntee)
!Returns a pointer to the current container element.
         implicit none
         class(gfc_cont_elem_t), pointer:: pntee     !out: container element currently pointed to by the iterator
         class(vec_tree_iter_t), intent(in):: this   !in: iterator
         integer(INTD), intent(out), optional:: ierr !out: error code

         pntee=>NULL()
         if(present(ierr)) then
          pntee=>this%val_it%pointee(ierr)
         else
          pntee=>this%val_it%pointee()
         endif
         return
        end function VecTreeIterPointee
!---------------------------------------------------------
        function VecTreeIterNext(this,elem_p) result(ierr)
!Moves the iterator to the next element of the vector part. If the next element is not
!in the tree yet, the tree subiterator will be set to GFC_IT_DONE. If <elem_p> is present,
!the next element will be returned in <elem_p> without moving the iterator.
         implicit none
         integer(INTD):: ierr                                            !out: error code
         class(vec_tree_iter_t), intent(inout):: this                    !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element
         class(*), pointer:: up
         class(tree_pos_t), pointer:: tpp

         if(present(elem_p)) then
          ierr=this%val_it%next(elem_p)
         else
          ierr=this%val_it%next()
          if(ierr.eq.GFC_SUCCESS) then
           ierr=this%pos_it%next()
           if(ierr.eq.GFC_SUCCESS) then
            up=>this%pos_it%get_value(ierr)
            if(ierr.eq.GFC_SUCCESS.and.associated(up)) then
             select type(up); class is(tree_pos_t); tpp=>up; class default; ierr=GFC_ERROR; end select
             if(ierr.eq.GFC_SUCCESS) call this%tree_it%jump_(tpp%pos)
            else
             ierr=GFC_CORRUPTED_CONT
            endif
           endif
          endif
         endif
         call this%update_status_()
         return
        end function VecTreeIterNext
!-------------------------------------------------------------
        function VecTreeIterPrevious(this,elem_p) result(ierr)
!Moves the iterator to the previous element of the vector part. If the previous element is not
!in the tree yet, the tree subiterator will be set to GFC_IT_DONE. If <elem_p> is present,
!the previous element will be returned in <elem_p> without moving the iterator.
         implicit none
         integer(INTD):: ierr                                            !out: error code
         class(vec_tree_iter_t), intent(inout):: this                    !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element
         class(*), pointer:: up
         class(tree_pos_t), pointer:: tpp

         if(present(elem_p)) then
          ierr=this%val_it%previous(elem_p)
         else
          ierr=this%val_it%previous()
          if(ierr.eq.GFC_SUCCESS) then
           ierr=this%pos_it%previous()
           if(ierr.eq.GFC_SUCCESS) then
            up=>this%pos_it%get_value(ierr)
            if(ierr.eq.GFC_SUCCESS.and.associated(up)) then
             select type(up); class is(tree_pos_t); tpp=>up; class default; ierr=GFC_ERROR; end select
             if(ierr.eq.GFC_SUCCESS) call this%tree_it%jump_(tpp%pos)
            else
             ierr=GFC_CORRUPTED_CONT
            endif
           endif
          endif
         endif
         call this%update_status_()
         return
        end function VecTreeIterPrevious
!--------------------------------------------------------------
        function VecTreeIterGetLength(this,ierr) result(length)
!Returns the current length of the associated vector tree container.
         implicit none
         integer(INTL):: length                      !out: length of the vector tree
         class(vec_tree_iter_t), intent(in):: this   !in: vector tree iterator
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(associated(this%container)) then
          length=this%container%length(errc)
         else
          length=-1_INTL; errc=GFC_NULL_CONT
         endif
         if(present(ierr)) ierr=errc
         return
        end function VecTreeIterGetLength
!--------------------------------------------------------------
        function VecTreeIterGetOffset(this,ierr) result(offset)
!Returns the current offset (position) in the vector part of a vector tree.
         implicit none
         integer(INTL):: offset                      !out: offset
         class(vec_tree_iter_t), intent(in):: this   !in: vector tree iterator
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         offset=this%val_it%get_offset(errc)
         if(present(ierr)) ierr=errc
         return
        end function VecTreeIterGetOffset
!------------------------------------------------------------
        function VecTreeIterGetLevel(this,ierr) result(level)
!Returns the distance from the tree root for the current iterator position.
         implicit none
         integer(INTD):: level                        !out: distance from the tree root
         class(vec_tree_iter_t), intent(inout):: this !in: vector tree iterator
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         level=this%tree_it%get_level(errc)
         if(present(ierr)) ierr=errc
         return
        end function VecTreeIterGetLevel
!--------------------------------------------------------------
        function VecTreeIterGetRootId(this,ierr) result(offset)
!Returns the offset of the tree root node.
         implicit none
         integer(INTL):: offset                      !out: offset of the tree root node
         class(vec_tree_iter_t), intent(in):: this   !in: vector tree iterator
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         class(gfc_cont_elem_t), pointer:: gep
         class(*), pointer:: up

         offset=-1_INTL; gep=>this%tree_it%get_root(errc)
         if(errc.eq.GFC_SUCCESS.and.associated(gep)) then
          up=>gep%get_value(errc)
          if(errc.eq.GFC_SUCCESS) then
           select type(up); type is(integer(INTL)); offset=up; end select
           if(offset.lt.0_INTL) errc=GFC_CORRUPTED_CONT
          endif
         else
          if(errc.eq.GFC_SUCCESS) errc=GFC_EMPTY_CONT
         endif
         if(present(ierr)) ierr=errc
         return
        end function VecTreeIterGetRootId
!-------------------------------------------------------------------
        function VecTreeIterGetNumChildren(this,ierr) result(nchild)
!Returns the total number of children at the current (tree) iterator position.
         implicit none
         integer(INTL):: nchild                      !out: number of children
         class(vec_tree_iter_t), intent(in):: this   !in: vector tree iterator
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         nchild=this%tree_it%get_num_children(errc)
         if(present(ierr)) ierr=errc
         return
        end function VecTreeIterGetNumChildren
!------------------------------------------------------------------
        function VecTreeIterGetNumSiblings(this,ierr) result(nsibl)
!Returns the total number of siblings at the current (tree) iterator position.
         implicit none
         integer(INTL):: nsibl                       !out: number of siblings
         class(vec_tree_iter_t), intent(in):: this   !in: vector tree iterator
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         nsibl=this%tree_it%get_num_siblings(errc)
         if(present(ierr)) ierr=errc
         return
        end function VecTreeIterGetNumSiblings
!----------------------------------------------------------------
        function VecTreeIterOnFirstSibling(this,ierr) result(res)
!Returns GFC_TRUE if the iterator is positioned on the first sibling.
         implicit none
         integer(INTD):: res                         !out: {GFC_TRUE,GFC_FALSE,GFC_ERROR}
         class(vec_tree_iter_t), intent(in):: this   !in: vector tree iterator
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         res=this%tree_it%on_first_sibling(errc)
         if(present(ierr)) ierr=errc
         return
        end function VecTreeIterOnFirstSibling
!---------------------------------------------------------------
        function VecTreeIterOnLastSibling(this,ierr) result(res)
!Returns GFC_TRUE if the iterator is positioned on the last sibling.
         implicit none
         integer(INTD):: res                         !out: {GFC_TRUE,GFC_FALSE,GFC_ERROR}
         class(vec_tree_iter_t), intent(in):: this   !in: vector tree iterator
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         res=this%tree_it%on_last_sibling(errc)
         if(present(ierr)) ierr=errc
         return
        end function VecTreeIterOnLastSibling
!-----------------------------------------------------------------------
        function VecTreeIterElementValue(this,offset,ierr) result(val_p)
!Returns a pointer to the value of a specific element.
         implicit none
         class(*), pointer:: val_p                   !out: pointer to the value of a specific element
         class(vec_tree_iter_t), intent(in):: this   !in: vector tree iterator
         integer(INTL), intent(in):: offset          !in: offset of the element in the vector
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         val_p=>this%val_it%element_value(offset,errc)
         if(present(ierr)) ierr=errc
         return
        end function VecTreeIterElementValue
!-----------------------------------------------------------
        function VecTreeIterMoveTo(this,offset) result(ierr)
!Moves the iterator to a specific vector position.
         implicit none
         integer(INTD):: ierr                         !out: error code
         class(vec_tree_iter_t), intent(inout):: this !inout: vector tree iterator
         integer(INTL), intent(in):: offset           !in: offset
         class(*), pointer:: up
         class(tree_pos_t), pointer:: tpp

         ierr=this%val_it%move_to(offset)
         if(ierr.eq.GFC_SUCCESS) then
          ierr=this%pos_it%move_to(offset)
          if(ierr.eq.GFC_SUCCESS) then
           up=>this%pos_it%get_value(ierr)
           if(ierr.eq.GFC_SUCCESS.and.associated(up)) then
            select type(up); class is(tree_pos_t); tpp=>up; class default; ierr=GFC_ERROR; end select
            if(ierr.eq.GFC_SUCCESS) call this%tree_it%jump_(tpp%pos)
           else
            ierr=GFC_ERROR
           endif
          endif
         endif
         call this%update_status_()
         return
        end function VecTreeIterMoveTo
!-----------------------------------------------------------------------
        function VecTreeIterMoveToSibling(this,to_previous) result(ierr)
!Moves the iterator either to the next or to the previous sibling in the tree.
         implicit none
         integer(INTD):: ierr                         !out: error code
         class(vec_tree_iter_t), intent(inout):: this !inout: vector tree iterator
         logical, intent(in), optional:: to_previous  !in: if TRUE, the iterator will move to the previous sibling (defaults to FALSE)
         class(*), pointer:: up
         integer(INTL):: offset

         if(present(to_previous)) then
          ierr=this%tree_it%move_to_sibling(to_previous)
         else
          ierr=this%tree_it%move_to_sibling()
         endif
         if(ierr.eq.GFC_SUCCESS) then
          up=>this%tree_it%get_value(ierr)
          if(ierr.eq.GFC_SUCCESS.and.associated(up)) then
           select type(up); type is(integer(INTL)); offset=up; class default; ierr=GFC_ERROR; end select
           if(ierr.eq.GFC_SUCCESS) then
            ierr=this%val_it%move_to(offset)
            if(ierr.eq.GFC_SUCCESS) ierr=this%pos_it%move_to(offset)
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         call this%update_status_()
         return
        end function VecTreeIterMoveToSibling
!---------------------------------------------------------
        function VecTreeIterMoveToChild(this) result(ierr)
!Moves the iterator to the first child in the tree.
         implicit none
         integer(INTD):: ierr                         !out: error code
         class(vec_tree_iter_t), intent(inout):: this !inout: vector tree iterator
         class(*), pointer:: up
         integer(INTL):: offset

         ierr=this%tree_it%move_to_child()
         if(ierr.eq.GFC_SUCCESS) then
          up=>this%tree_it%get_value(ierr)
          if(ierr.eq.GFC_SUCCESS.and.associated(up)) then
           select type(up); type is(integer(INTL)); offset=up; class default; ierr=GFC_ERROR; end select
           if(ierr.eq.GFC_SUCCESS) then
            ierr=this%val_it%move_to(offset)
            if(ierr.eq.GFC_SUCCESS) ierr=this%pos_it%move_to(offset)
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         call this%update_status_()
         return
        end function VecTreeIterMoveToChild
!----------------------------------------------------------
        function VecTreeIterMoveToParent(this) result(ierr)
!Moves the iterator to the parent in the tree.
         implicit none
         integer(INTD):: ierr                         !out: error code
         class(vec_tree_iter_t), intent(inout):: this !inout: vector tree iterator
         class(*), pointer:: up
         integer(INTL):: offset

         ierr=this%tree_it%move_to_parent()
         if(ierr.eq.GFC_SUCCESS) then
          up=>this%tree_it%get_value(ierr)
          if(ierr.eq.GFC_SUCCESS.and.associated(up)) then
           select type(up); type is(integer(INTL)); offset=up; class default; ierr=GFC_ERROR; end select
           if(ierr.eq.GFC_SUCCESS) then
            ierr=this%val_it%move_to(offset)
            if(ierr.eq.GFC_SUCCESS) ierr=this%pos_it%move_to(offset)
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         call this%update_status_()
         return
        end function VecTreeIterMoveToParent
!-------------------------------------------------------------
        function VecTreeIterMoveUp(this,num_hops) result(ierr)
!Moves the iterator towards the tree root a specific number of hops.
         implicit none
         integer(INTD):: ierr                         !out: error code
         class(vec_tree_iter_t), intent(inout):: this !inout: vector tree iterator
         integer(INTD), intent(in):: num_hops         !in: number of hops to move up
         integer(INTD):: n

         ierr=GFC_SUCCESS
         if(num_hops.ge.0) then
          n=num_hops
          do while(n.gt.0)
           ierr=this%move_to_parent(); if(ierr.ne.GFC_SUCCESS) exit
           n=n-1
          enddo
         else
          ierr=GFC_INVALID_ARGS
         endif
         return
        end function VecTreeIterMoveUp
!----------------------------------------------------------------------
        function VecTreeIterMoveToCousin(this,to_previous) result(ierr)
!Moves the iterator either to the next or to the previous cousin in the tree.
!A cousin is a tree vertex at the same tree level.
         implicit none
         integer(INTD):: ierr                         !out: error code
         class(vec_tree_iter_t), intent(inout):: this !inout: vector tree iterator
         logical, intent(in), optional:: to_previous  !in: if TRUE, the iterator will move to the previous sibling (defaults to FALSE)
         class(*), pointer:: up
         integer(INTL):: offset

         if(present(to_previous)) then
          ierr=this%tree_it%move_to_cousin(to_previous)
         else
          ierr=this%tree_it%move_to_cousin()
         endif
         if(ierr.eq.GFC_SUCCESS) then
          up=>this%tree_it%get_value(ierr)
          if(ierr.eq.GFC_SUCCESS.and.associated(up)) then
           select type(up); type is(integer(INTL)); offset=up; class default; ierr=GFC_ERROR; end select
           if(ierr.eq.GFC_SUCCESS) then
            ierr=this%val_it%move_to(offset)
            if(ierr.eq.GFC_SUCCESS) ierr=this%pos_it%move_to(offset)
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         call this%update_status_()
         return
        end function VecTreeIterMoveToCousin
!------------------------------------------------------------------------------------
#ifdef NO_GNU
        function VecTreeIterAppend(this,elem_val,assoc_only,copy_ctor_f) result(ierr) !`GCC bug
#else
        function VecTreeIterAppend(this,elem_val,assoc_only) result(ierr) !`GCC bug
#endif
!Appends a new element to the end of the vector. The iterator position is kept unchanged,
!unless the container is empty in which case it will be reset to the first element.
         implicit none
         integer(INTD):: ierr                         !out: error code
         class(vec_tree_iter_t), intent(inout):: this !inout: vector tree iterator
         class(*), target, intent(in):: elem_val      !in: value to be stored
         logical, intent(in), optional:: assoc_only   !in: storage type: TRUE:by reference, FALSE:by value (default)
#ifdef NO_GNU
         procedure(gfc_copy_i), optional:: copy_ctor_f !user-defined generic copy constructor (when storing by value only)
#endif
         type(tree_pos_t):: trpos
         logical:: assoc

         if(present(assoc_only)) then; assoc=assoc_only; else; assoc=.FALSE.; endif
#ifdef NO_GNU
         if(present(copy_ctor_f)) then
          ierr=this%val_it%append(elem_val,assoc,copy_ctor_f)
         else
#endif
          ierr=this%val_it%append(elem_val,assoc)
#ifdef NO_GNU
         endif
#endif
         if(ierr.eq.GFC_SUCCESS) ierr=this%pos_it%append(trpos)
         call this%update_status_()
         return
        end function VecTreeIterAppend
!----------------------------------------------------------------------
        function VecTreeIterAddLeaf(this,elem_num,no_move) result(ierr)
!Adds a specific vector element as a leaf to the tree. The iterator position
!will move to the newly added element, unless <no_move>=TRUE.
         implicit none
         integer(INTD):: ierr                         !out: error code
         class(vec_tree_iter_t), intent(inout):: this !inout: vector tree iterator
         integer(INTL), intent(in):: elem_num         !in: vector element number: [0..max]
         logical, intent(in), optional:: no_move      !in: if TRUE, the iterator position will not change (defaults to FALSE)
         logical:: nmove
         class(*), pointer:: up
         class(tree_pos_t), pointer:: tpp
         class(gfc_cont_elem_t), pointer:: gep

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_DONE) then; ierr=this%reset(); ierr=this%get_status(); endif
         if(ierr.eq.GFC_IT_ACTIVE) then
          ierr=GFC_SUCCESS
          if(present(no_move)) then; nmove=no_move; else; nmove=.FALSE.; endif
          if(elem_num.ge.0_INTL.and.elem_num.lt.this%get_length()) then
           up=>this%pos_it%element_value(elem_num,ierr)
           if(ierr.eq.GFC_SUCCESS) then
            if(associated(up)) then
             select type(up); class is(tree_pos_t); tpp=>up; class default; ierr=GFC_ERROR; end select
             if(ierr.eq.GFC_SUCCESS) then
              if(.not.associated(tpp%pos)) then !vector element is not in the tree
               ierr=this%tree_it%get_status(); if(ierr.eq.GFC_IT_EMPTY) nmove=.FALSE.
               ierr=this%tree_it%add_leaf(elem_num,no_move=nmove)
               if(ierr.eq.GFC_SUCCESS) then
                if(nmove) then
                 gep=>this%tree_it%get_child(-1,ierr) !-1:last child
                else
                 gep=>this%tree_it%pointee(ierr)
                endif
                if(ierr.eq.GFC_SUCCESS) then
                 select type(gep); class is(tree_vertex_t); tpp%pos=>gep; class default; ierr=GFC_ERROR; end select
                 if(ierr.eq.GFC_SUCCESS.and.(.not.nmove)) then
                  ierr=this%pos_it%move_to(elem_num)
                  if(ierr.eq.GFC_SUCCESS) ierr=this%val_it%move_to(elem_num)
                 endif
                endif
                gep=>NULL()
               endif
              else
               ierr=GFC_INVALID_ARGS
              endif
             endif
             tpp=>NULL()
            else
             ierr=GFC_ERROR
            endif
           else
            ierr=GFC_ERROR
           endif
          else
           ierr=GFC_INVALID_ARGS
          endif
         endif
         call this%update_status_()
         return
        end function VecTreeIterAddLeaf
!--------------------------------------------------------------
        function VecTreeIterDeleteAll(this,dtor_f) result(ierr)
!Deletes all elements of the vector tree, resetting it to an empty state.
         implicit none
         integer(INTD):: ierr                         !out: error code
         class(vec_tree_iter_t), intent(inout):: this !inout: vector tree iterator
         procedure(gfc_destruct_i), optional:: dtor_f !in: explicit container value destructor
         integer(INTD):: errc

         ierr=this%reset()
         if(ierr.eq.GFC_SUCCESS) then
          errc=this%tree_it%delete_all(); if(errc.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=errc
          errc=this%pos_it%delete_all(); if(errc.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=errc
          if(present(dtor_f)) then
           errc=this%val_it%delete_all(dtor_f); if(errc.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=errc
          else
           errc=this%val_it%delete_all(); if(errc.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=errc
          endif
         endif
         call this%update_status_()
         return
        end function VecTreeIterDeleteAll

       end module gfc_vec_tree
