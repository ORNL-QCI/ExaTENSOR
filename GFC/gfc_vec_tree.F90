!Generic Fortran Containers (GFC): Vector tree, combines vector and tree:
!The elements are initially inserted in a vector with an option to be
!later added in a tree, thus imposing a tree relationship on them.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2017/03/13

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
          procedure, public:: init=>VecTreeIterInit                     !associates the iterator with a container and sets its position to the beginning
          procedure, public:: reset=>VecTreeIterReset                   !resets the iterator to the beginning of the container (only the vector part if tree is not fully built)
          procedure, public:: reset_back=>VecTreeIterResetBack          !resets the iterator to the end of the container (only the vector part if tree is not fully built)
          procedure, public:: release=>VecTreeIterRelease               !releases the iterator (dissociates it from its container)
          procedure, public:: pointee=>VecTreeIterPointee               !returns a pointer to the container element currently pointed to
          procedure, public:: next=>VecTreeIterNext                     !moves the iterator to the next container element (only the vector part if tree is not fully built)
          procedure, public:: previous=>VecTreeIterPrevious             !moves the iterator to the previous container element (only the vector part if tree is not fully built)
#if 0
          procedure, public:: get_length=>VecTreeIterGetLength          !returns the current length of the container
          procedure, public:: get_offset=>VecTreeIterGetOffset          !returns the current offset in the container
          procedure, public:: move_to=>VecTreeIterMoveTo                !moves the iterator to the specific container element by its sequential number
          procedure, public:: move_to_sibling=>VecTreeIterMoveToSibling !moves the iterator to the next/previous sibling of the current element in the tree
          procedure, public:: move_to_child=>VecTreeIterMoveToChild     !moves the iterator to the first child of the current element in the tree
          procedure, public:: move_to_parent=>VecTreeIterMoveToParent   !moves the iterator to the parent of the current element in the tree
          procedure, public:: move_to_cousin=>VecTreeIterMoveToCousin   !moves the iterator to the next/previous cousin of the current element in the tree (same tree level)
          procedure, public:: append=>VecTreeIterAppend                 !appends a new element at the end of the container
          procedure, public:: delete=>VecTreeIterDelete                 !deletes an element at the current iterator position
          procedure, public:: add_leaf=>VecTreeIterAddLeaf              !appends a new child element to the tree vertex currently pointed to
          procedure, public:: delete_leaf=>VecTreeIterDeleteLeaf        !deletes the tree leaf element currently pointed to by the iterator
          procedure, public:: delete_all=>VecTreeIterDeleteAll          !deletes all container elements
#endif
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
        private VecTreeIterInit
        private VecTreeIterReset
        private VecTreeIterResetBack
        private VecTreeIterRelease
        private VecTreeIterPointee
        private VecTreeIterNext
        private VecTreeIterPrevious
#if 0
        private VecTreeIterGetLength
        private VecTreeIterGetOffset
        private VecTreeIterMoveTo
        private VecTreeIterMoveToSibling
        private VecTreeIterMoveToChild
        private VecTreeIterMoveToParent
        private VecTreeIterMoveToCousin
        private VecTreeIterAppend
        private VecTreeIterDelete
        private VecTreeIterAddLeaf
        private VecTreeIterDeleteLeaf
        private VecTreeIterDeleteAll
#endif
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
![vec_tree_iter_t]======================================
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
           up=>this%pos_it%get_value(ierr)
           if(ierr.eq.GFC_SUCCESS.and.associated(up)) then
            select type(up); class is(tree_pos_t); tpp=>up; class default; ierr=GFC_ERROR; end select
            if(ierr.eq.GFC_SUCCESS) then
             if(associated(tpp%pos)) then
              call this%tree_it%jump_(tpp%pos)
             else !the first vector element is not in the tree yet
              ierr=this%tree_it%set_status_(GFC_IT_DONE)
             endif
            endif
           else
            ierr=GFC_CORRUPTED_CONT
           endif
          endif
         endif
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
           up=>this%pos_it%get_value(ierr)
           if(ierr.eq.GFC_SUCCESS.and.associated(up)) then
            select type(up); class is(tree_pos_t); tpp=>up; class default; ierr=GFC_ERROR; end select
            if(ierr.eq.GFC_SUCCESS) then
             if(associated(tpp%pos)) then
              call this%tree_it%jump_(tpp%pos)
             else !the last vector element is not in the tree yet
              ierr=this%tree_it%set_status_(GFC_IT_DONE)
             endif
            endif
           else
            ierr=GFC_CORRUPTED_CONT
           endif
          endif
         endif
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
         class(gfc_cont_elem_t), pointer:: pntee     !out: container element currently pointed to be the iterator
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
             if(ierr.eq.GFC_SUCCESS) then
              if(associated(tpp%pos)) then
               call this%tree_it%jump_(tpp%pos)
              else !the vector element is not in the tree yet
               ierr=this%tree_it%set_status_(GFC_IT_DONE)
              endif
             endif
            else
             ierr=GFC_CORRUPTED_CONT
            endif
           endif
          endif
         endif
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
             if(ierr.eq.GFC_SUCCESS) then
              if(associated(tpp%pos)) then
               call this%tree_it%jump_(tpp%pos)
              else !the vector element is not in the tree yet
               ierr=this%tree_it%set_status_(GFC_IT_DONE)
              endif
             endif
            else
             ierr=GFC_CORRUPTED_CONT
            endif
           endif
          endif
         endif
         return
        end function VecTreeIterPrevious

       end module gfc_vec_tree
