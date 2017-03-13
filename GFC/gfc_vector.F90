!Generic Fortran Containers (GFC): Vector (non-contiguous)
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

       module gfc_vector
        use gfc_base
        use timers
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !output device
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
        integer(INTD), private:: DEBUG=0    !debugging level (0:none)
 !Vector length:
        integer(INTL), parameter, public:: GFC_VECTOR_SEG_LENGTH=(2_INTL)**9              !default vector segment length
        integer(INTL), parameter, public:: GFC_VECTOR_MAX_LENGTH=GFC_VECTOR_SEG_LENGTH**4 !default max vector capacity
        integer(INTL), parameter, public:: GFC_VECTOR_OVERHEAD=GFC_VECTOR_SEG_LENGTH*4    !default memory overhead
!TYPES:
 !Vector element:
        type, extends(gfc_cont_elem_t), public:: vector_elem_t
        contains
         procedure, private:: VectorElemConstruct
         generic, public:: vector_elem_ctor=>VectorElemConstruct !constructs the content (value) of the vector element
        end type vector_elem_t
 !Vector segment:
        type, private:: vector_seg_t
         type(vector_elem_t), allocatable, private:: seg_elem(:) !elements of the vector segment: [0..max_elems-1]
         contains
          procedure, private:: VectorSegConstruct
          generic:: vector_seg_ctor=>VectorSegConstruct          !constructs a vector segment
          final:: vector_seg_dtor                                !destructs a vector segment
        end type vector_seg_t
 !Vector segment batch:
        type, private:: vector_batch_t
         type(vector_seg_t), allocatable, private:: batch_seg(:) !segments of the batch: [0..max_segs-1]
         contains
          procedure, private:: VectorBatchConstruct
          generic:: vector_batch_ctor=>VectorBatchConstruct      !constructs a vector batch
          final:: vector_batch_dtor                              !destructs a vector batch
        end type vector_batch_t
 !Vector tile:
        type, private:: vector_tile_t
         type(vector_batch_t), allocatable, private:: tile_batch(:) !batches of the tile: [0..max_batches-1]
         contains
          procedure, private:: VectorTileConstruct
          generic:: vector_tile_ctor=>VectorTileConstruct           !constructs a vector tile
          final:: vector_tile_dtor                                  !destructs a vector tile
        end type vector_tile_t
 !Vector:
        type, extends(gfc_container_t), public:: vector_t
         integer(INTL), private:: max_length=GFC_VECTOR_MAX_LENGTH !vector capacity (max possible number of elements)
         integer(INTL), private:: lbnd=0_INTL                      !lower bound of the vector
         integer(INTL), private:: ubnd=-1_INTL                     !upper bound of the vector
         type(vector_tile_t), allocatable, private:: vec_tile(:)   !vector tiles: [0..max_tiles-1]
         contains
          procedure, public:: is_empty=>VectorIsEmpty       !returns GFC_TRUE if the vector is empty, GFC_FALSE otherwise (or error code)
          procedure, public:: is_full=>VectorIsFull         !returns GFC_TRUE if the vector is full, GFC_FALSE otherwise (or error code)
          procedure, public:: capacity=>VectorCapacity      !returns maximal length of the vector
          procedure, public:: length=>VectorLength          !returns current vector length = (upper_bound - lower_bound + 1)
          procedure, public:: lower_bound=>VectorLowerBound !returns vector lower bound
          procedure, public:: upper_bound=>VectorUpperBound !returns vector upper bound
          procedure, private:: incr_len_=>VectorIncrLen
          procedure, private:: decr_len_=>VectorDecrLen
          procedure, private:: adjust_structure_=>VectorAdjustStructure
        end type vector_t
 !Vector iterator:
        type, extends(gfc_iter_t), public:: vector_iter_t
         class(vector_elem_t), pointer, private:: current=>NULL() !current element of the vector
         class(vector_t), pointer, private:: container=>NULL()    !vector associated with the iterator
         integer(INTL), private:: curr_offset=-1_INTL             !current element offset in the vector
         contains
          procedure, public:: init=>VectorIterInit                  !initializes the iterator by associating it with a vector
          procedure, public:: reset=>VectorIterReset                !resets the iterator to the beginning of the vector
          procedure, public:: reset_back=>VectorIterResetBack       !resets the iterator to the end of the vector
          procedure, public:: release=>VectorIterRelease            !releases the iterator (dissocaites it from its container)
          procedure, public:: pointee=>VectorIterPointee            !returns the container element currently pointed to by the iterator
          procedure, public:: next=>VectorIterNext                  !moves the iterator to the next vector element
          procedure, public:: previous=>VectorIterPrevious          !moves the iterator to the previous vector element
          procedure, public:: get_length=>VectorIterGetLength       !returns the current length of the vector
          procedure, public:: get_offset=>VectorIterGetOffset       !returns the offset of the current iterator position: [0..MAX]
          procedure, public:: element=>VectorIterElement            !returns a pointer to the specific vector element
          procedure, public:: element_value=>VectorIterElementValue !returns an unlimited polymorphic pointer to the value of a specific vector element
          procedure, public:: move_to=>VectorIterMoveTo             !moves the iterator to the specific vector element
          procedure, public:: append=>VectorIterAppend              !appends a new element at the end of the vector
          procedure, public:: insert=>VectorIterInsert              !inserts a new element at the current iterator position
          procedure, public:: swap_last=>VectorIterSwapLast         !swaps the currently pointed to element with the last element of the vector
          procedure, public:: delete=>VectorIterDelete              !deletes an element at the current iterator position
          procedure, public:: delete_all=>VectorIterDeleteAll       !deletes all elements of the vector
        end type vector_iter_t
!VISIBILITY:
        private flat2quadruplet
        private quadruplet2flat
 !vector_elem_t:
        private VectorElemConstruct
 !vector_seg_t:
        private VectorSegConstruct
 !vector_batch_t:
        private VectorBatchConstruct
 !vector_tile_t:
        private VectorTileConstruct
 !vector_t:
        private VectorIsEmpty
        private VectorIsFull
        private VectorCapacity
        private VectorLength
        private VectorLowerBound
        private VectorUpperBound
        private VectorIncrLen
        private VectorDecrLen
        private VectorAdjustStructure
 !vector_iter_t:
        private VectorIterInit
        private VectorIterReset
        private VectorIterResetBack
        private VectorIterRelease
        private VectorIterPointee
        private VectorIterNext
        private VectorIterPrevious
        private VectorIterGetLength
        private VectorIterGetOffset
        private VectorIterElement
        private VectorIterElementValue
        private VectorIterMoveTo
        private VectorIterAppend
        private VectorIterInsert
        private VectorIterSwapLast
        private VectorIterDelete
        private VectorIterDeleteAll

       contains
!IMPLEMENTATION:
!==================================================
        subroutine flat2quadruplet(flat,quadruplet)
         implicit none
         integer(INTL), intent(in):: flat               !in: flat offset (0..MAX)
         integer(INTD), intent(inout):: quadruplet(1:4) !out: 4-component representation with seniority 1->4
         integer(INTL):: i,k
         integer(INTD):: j

         i=flat
         do j=1,4
          k=i/GFC_VECTOR_SEG_LENGTH
          quadruplet(j)=i-k*GFC_VECTOR_SEG_LENGTH
          i=k
         enddo
         if(i.ne.0_INTL) then
          write(CONS_OUT,'("#FATAL(GFC::vector): index out of range: ",i18)') flat
          stop
         endif
         return
        end subroutine flat2quadruplet
!--------------------------------------------------------
        function quadruplet2flat(quadruplet) result(flat)
         implicit none
         integer(INTL):: flat                        !out: flat offset
         integer(INTD), intent(in):: quadruplet(1:4) !in: 4-component representation with seniority 1->4
         integer(INTD):: j

         flat=quadruplet(4)
         do j=3,1,-1
          flat=flat*GFC_VECTOR_SEG_LENGTH+quadruplet(j)
         enddo
         return
        end function quadruplet2flat
![vector_elem_t]============================================================
#ifdef NO_GNU
        subroutine VectorElemConstruct(this,obj,ierr,assoc_only,copy_ctor_f) !`GCC has a bug with this line
#else
        subroutine VectorElemConstruct(this,obj,ierr,assoc_only)
#endif
!Constructs the content of the vector element.
         implicit none
         class(vector_elem_t), intent(inout):: this     !inout: vector element (empty on input)
         class(*), target, intent(in):: obj             !in: assigned value
         integer(INTD), intent(out), optional:: ierr    !out: error code
         logical, intent(in), optional:: assoc_only     !in: if TRUE, the value will be assigned by reference, otherwise by value (allocated): Defaults to FALSE
#ifdef NO_GNU
         procedure(gfc_copy_i), optional:: copy_ctor_f  !in: generic copy constructor
#endif
         integer(INTD):: errc

#ifdef NO_GNU
         if(present(copy_ctor_f)) then
          if(present(assoc_only)) then
           call this%construct_base(obj,errc,assoc_only,copy_ctor_f)
          else
           call this%construct_base(obj,errc,copy_ctor_f=copy_ctor_f)
          endif
         else
#endif
          if(present(assoc_only)) then
           call this%construct_base(obj,errc,assoc_only=assoc_only)
          else
           call this%construct_base(obj,errc)
          endif
#ifdef NO_GNU
         endif
#endif
         if(present(ierr)) ierr=errc
         return
        end subroutine VectorElemConstruct
![vector_seg_t]==========================================
        subroutine VectorSegConstruct(this,capacity,ierr)
         implicit none
         class(vector_seg_t), intent(out):: this     !out: vector segment
         integer(INTD), intent(in):: capacity        !in: requested capacity
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
         if(capacity.gt.0) then
          allocate(this%seg_elem(0:capacity-1),STAT=errc)
          if(errc.ne.0) errc=GFC_MEM_ALLOC_FAILED
         else
          errc=GFC_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine VectorSegConstruct
!---------------------------------------
        subroutine vector_seg_dtor(this)
         implicit none
         type(vector_seg_t):: this

         if(allocated(this%seg_elem)) deallocate(this%seg_elem)
         return
        end subroutine vector_seg_dtor
![vector_batch_t]==========================================
        subroutine VectorBatchConstruct(this,capacity,ierr)
         implicit none
         class(vector_batch_t), intent(out):: this   !out: vector batch
         integer(INTD), intent(in):: capacity        !in: requested capacity
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
         if(capacity.gt.0) then
          allocate(this%batch_seg(0:capacity-1),STAT=errc)
          if(errc.ne.0) errc=GFC_MEM_ALLOC_FAILED
         else
          errc=GFC_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine VectorBatchConstruct
!-----------------------------------------
        subroutine vector_batch_dtor(this)
         implicit none
         type(vector_batch_t):: this

         if(allocated(this%batch_seg)) deallocate(this%batch_seg)
         return
        end subroutine vector_batch_dtor
![vector_tile_t]==========================================
        subroutine VectorTileConstruct(this,capacity,ierr)
         implicit none
         class(vector_tile_t), intent(out):: this    !out: vector tile
         integer(INTD), intent(in):: capacity        !in: requested capacity
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
         if(capacity.gt.0) then
          allocate(this%tile_batch(0:capacity-1),STAT=errc)
          if(errc.ne.0) errc=GFC_MEM_ALLOC_FAILED
         else
          errc=GFC_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine VectorTileConstruct
!----------------------------------------
        subroutine vector_tile_dtor(this)
         implicit none
         type(vector_tile_t):: this

         if(allocated(this%tile_batch)) deallocate(this%tile_batch)
         return
        end subroutine vector_tile_dtor
![vector_t]=====================================
        function VectorIsEmpty(this) result(res)
!Returns GFC_TRUE if the vector is empty, GFC_FALSE otherwise (or error code).
         implicit none
         integer(INTD):: res                !out: result of the query (or error code)
         class(vector_t), intent(in):: this !in: vector

         res=GFC_FALSE
         if(this%ubnd.lt.this%lbnd) res=GFC_TRUE
         return
        end function VectorIsEmpty
!----------------------------------------------
        function VectorIsFull(this) result(res)
!Returns GFC_TRUE if the vector is full, GFC_FALSE otherwise (or error code).
         implicit none
         integer(INTD):: res                !out: result of the query (or error code)
         class(vector_t), intent(in):: this !in: vector

         res=GFC_FALSE
         if(this%ubnd-this%lbnd+1_INTL.ge.this%max_length) res=GFC_TRUE
         return
        end function VectorIsFull
!----------------------------------------------------------
        function VectorCapacity(this,ierr) result(capacity)
!Returns the max capacity of the vector.
         implicit none
         integer(INTL):: capacity                    !out: vector capacity (max length)
         class(vector_t), intent(in):: this          !in: vector
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
         capacity=this%max_length
         if(present(ierr)) ierr=errc
         return
        end function VectorCapacity
!------------------------------------------------------
        function VectorLength(this,ierr) result(length)
!Returns the current length of the vector.
         implicit none
         integer(INTL):: length                      !out: current vector length
         class(vector_t), intent(in):: this          !in: vector
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
         length=this%ubnd-this%lbnd+1_INTL
         if(present(ierr)) ierr=errc
         return
        end function VectorLength
!---------------------------------------------------------
        function VectorLowerBound(this,ierr) result(lower)
!Returns the lower bound of the vector.
         implicit none
         integer(INTL):: lower                       !out: vector lower bound
         class(vector_t), intent(in):: this          !in: vector
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
         lower=this%lbnd
         if(present(ierr)) ierr=errc
         return
        end function VectorLowerBound
!---------------------------------------------------------
        function VectorUpperBound(this,ierr) result(upper)
!Returns the upper bound of the vector.
         implicit none
         integer(INTL):: upper                       !out: vector upper bound
         class(vector_t), intent(in):: this          !in: vector
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
         upper=this%ubnd
         if(present(ierr)) ierr=errc
         return
        end function VectorUpperBound
!-------------------------------------
        subroutine VectorIncrLen(this)
         implicit none
         class(vector_t), intent(inout):: this !inout: vector

         this%ubnd=this%ubnd+1_INTL
         return
        end subroutine VectorIncrLen
!-------------------------------------
        subroutine VectorDecrLen(this)
         implicit none
         class(vector_t), intent(inout):: this !inout: vector

         this%ubnd=this%ubnd-1_INTL
         return
        end subroutine VectorDecrLen
!----------------------------------------------------------------------------
        function VectorAdjustStructure(this,offset,for_deletion) result(ierr)
!Adjusts the vector structure either prior to insertion (for_deletion=.FALSE.)
!or after a deletion (for_deletion=.TRUE.).
         implicit none
         integer(INTD):: ierr                         !out: error code
         class(vector_t), intent(inout):: this        !inout: vector
         integer(INTL), intent(in):: offset           !in: offset in the vector (for insertion or deletion)
         logical, intent(in), optional:: for_deletion !in: if TRUE, activates adjustment after a deletion, defaults to FALSE (adjustment prior to insertion)
         logical:: fordel
         integer(INTL):: i,l
         integer(INTD):: q(1:4),w(1:4),j
         integer:: errc

         ierr=GFC_SUCCESS; errc=0
         if(present(for_deletion)) then; fordel=for_deletion; else; fordel=.FALSE.; endif
         if(offset.ge.0_INTL) then
          l=this%length() !original vector length (prior to insertion/deletion)
          if(fordel) then !adjustment after a deletion of a vector element at position <offset>
           if(offset.lt.l) then
            if(offset.lt.l-1_INTL) then !inside deletion: shift all subsequent elements, O(N) in general
             do i=offset,l-2_INTL
              call flat2quadruplet(i,q); call flat2quadruplet(i+1_INTL,w)
              this%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))=&
              &this%vec_tile(w(4))%tile_batch(w(3))%batch_seg(w(2))%seg_elem(w(1))
             enddo
            endif
            call flat2quadruplet(l-1_INTL,q) !last vector element is now empty
            call this%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))%clean_(ierr)
            if(q(2).eq.0.and.q(1).eq.0) then !batch has become empty, deallocate it
             do j=GFC_VECTOR_SEG_LENGTH-1,0,-1
              if(allocated(this%vec_tile(q(4))%tile_batch(q(3))%batch_seg(j)%seg_elem))&
               &deallocate(this%vec_tile(q(4))%tile_batch(q(3))%batch_seg(j)%seg_elem,STAT=errc)
              if(errc.ne.0) ierr=GFC_MEM_FREE_FAILED
              !print *,'Deallocating ',q(4),q(3),j,' --> ',errc !debug
             enddo
             if(q(3).eq.0) then
              do j=GFC_VECTOR_SEG_LENGTH-1,0,-1
               if(allocated(this%vec_tile(q(4))%tile_batch(j)%batch_seg))&
                &deallocate(this%vec_tile(q(4))%tile_batch(j)%batch_seg,STAT=errc)
               if(errc.ne.0) ierr=GFC_MEM_FREE_FAILED
               !print *,'Deallocating ',q(4),j,' --> ',errc !debug
              enddo
              if(q(4).eq.0) then !last element has been deleted: vector is empty
               do j=GFC_VECTOR_SEG_LENGTH-1,0,-1
                if(allocated(this%vec_tile(j)%tile_batch))&
                 &deallocate(this%vec_tile(j)%tile_batch,STAT=errc)
                if(errc.ne.0) ierr=GFC_MEM_FREE_FAILED
                !print *,'Deallocating ',j,' --> ',errc !debug
               enddo
               deallocate(this%vec_tile,STAT=errc); if(errc.ne.0) ierr=GFC_MEM_FREE_FAILED
               !print *,'Deallocating vec_tile --> ',errc !debug
              endif
             endif
            endif
           else
            ierr=GFC_INVALID_ARGS
           endif
          else !adjustment prior to insertion of an element at position <offset>
           if(l.lt.this%max_length) then
            call flat2quadruplet(l,q)
            if(.not.allocated(this%vec_tile)) allocate(this%vec_tile(0:GFC_VECTOR_SEG_LENGTH-1),STAT=errc)
            if(errc.eq.0) then
             if(.not.allocated(this%vec_tile(q(4))%tile_batch))&
               &allocate(this%vec_tile(q(4))%tile_batch(0:GFC_VECTOR_SEG_LENGTH-1),STAT=errc)
             if(errc.eq.0) then
              if(.not.allocated(this%vec_tile(q(4))%tile_batch(q(3))%batch_seg))&
                &allocate(this%vec_tile(q(4))%tile_batch(q(3))%batch_seg(0:GFC_VECTOR_SEG_LENGTH-1),STAT=errc)
              if(errc.eq.0) then
               if(.not.allocated(this%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem))&
                &allocate(this%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(0:GFC_VECTOR_SEG_LENGTH-1),STAT=errc)
               if(errc.eq.0) then
                if(offset.lt.l) then !inside insertion: shift all subsequent elements, O(N) in general
                 do i=l,offset+1_INTL,-1_INTL
                  call flat2quadruplet(i,q); call flat2quadruplet(i-1_INTL,w)
                  this%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))=&
                  &this%vec_tile(w(4))%tile_batch(w(3))%batch_seg(w(2))%seg_elem(w(1))
                 enddo
                endif
                call flat2quadruplet(offset,q) !element #<offset> is empty now
                call this%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))%clean_(ierr)
               else
                ierr=GFC_MEM_ALLOC_FAILED
               endif
              else
               ierr=GFC_MEM_ALLOC_FAILED
              endif
             else
              ierr=GFC_MEM_ALLOC_FAILED
             endif
            else
             ierr=GFC_MEM_ALLOC_FAILED
            endif
           else
            ierr=GFC_OVERFLOW
           endif
          endif
         else
          ierr=GFC_INVALID_ARGS
         endif
         return
        end function VectorAdjustStructure
![vector_iter_t]=======================================
        function VectorIterInit(this,cont) result(ierr)
!Initializes the iterator and resets it to the beginning of the container.
         implicit none
         integer(INTD):: ierr                              !out: error code
         class(vector_iter_t), intent(inout):: this        !inout: iterator
         class(gfc_container_t), target, intent(in):: cont !in: container

         ierr=GFC_SUCCESS
         select type(cont)
         class is(vector_t)
          this%container=>cont
          ierr=this%reset()
         class default
          ierr=GFC_INVALID_ARGS
         end select
         return
        end function VectorIterInit
!--------------------------------------------------
        function VectorIterReset(this) result(ierr)
!Resets the iterator to the beginning (first element).
         implicit none
         integer(INTD):: ierr                       !out: error code
         class(vector_iter_t), intent(inout):: this !inout: iterator

         ierr=GFC_SUCCESS
         if(associated(this%container)) then
          if(this%container%length().gt.0_INTL) then
           this%current=>this%container%vec_tile(0)%tile_batch(0)%batch_seg(0)%seg_elem(0)
           this%curr_offset=0_INTL
           ierr=this%set_status_(GFC_IT_ACTIVE)
          else
           this%current=>NULL(); this%curr_offset=-1_INTL
           ierr=this%set_status_(GFC_IT_EMPTY)
          endif
          call this%reset_count() !reset all iteration counters
         else
          this%current=>NULL(); this%curr_offset=-1_INTL
          ierr=this%set_status_(GFC_IT_NULL)
          ierr=GFC_IT_NULL
         endif
         return
        end function VectorIterReset
!------------------------------------------------------
        function VectorIterResetBack(this) result(ierr)
!Resets the iterator to the end (last element).
         implicit none
         integer(INTD):: ierr                       !out: error code
         class(vector_iter_t), intent(inout):: this !inout: iterator
         integer(INTD):: q(1:4)

         ierr=GFC_SUCCESS
         if(associated(this%container)) then
          this%curr_offset=this%container%length()-1_INTL
          if(this%curr_offset.ge.0_INTL) then
           call flat2quadruplet(this%curr_offset,q)
           this%current=>this%container%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))
           ierr=this%set_status_(GFC_IT_ACTIVE)
          else
           this%current=>NULL()
           ierr=this%set_status_(GFC_IT_EMPTY)
          endif
          call this%reset_count() !reset all iteration counters
         else
          this%current=>NULL(); this%curr_offset=-1_INTL
          ierr=this%set_status_(GFC_IT_NULL)
          ierr=GFC_IT_NULL
         endif
         return
        end function VectorIterResetBack
!----------------------------------------------------
        function VectorIterRelease(this) result(ierr)
!Dissociates the iterator from its container.
         implicit none
         integer(INTD):: ierr                       !out: error code
         class(vector_iter_t), intent(inout):: this !inout: iterator

         this%current=>NULL(); this%container=>NULL(); this%curr_offset=-1_INTL
         call this%reset_count(); ierr=this%set_status_(GFC_IT_NULL)
         return
        end function VectorIterRelease
!----------------------------------------------------------
        function VectorIterPointee(this,ierr) result(pntee)
!Returns the container element the iterator is currently pointing to.
         implicit none
         class(gfc_cont_elem_t), pointer:: pntee     !out: container element currently pointed to by the iterator
         class(vector_iter_t), intent(in):: this     !in: iterator
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
        end function VectorIterPointee
!--------------------------------------------------------
        function VectorIterNext(this,elem_p) result(ierr)
!If <elem_p> is absent, the iterator moves to the next element, if any.
!If <elem_p> is present, the iterator simply returns the next element in <elem_p> without moving.
         implicit none
         integer(INTD):: ierr                                            !out: error code (0:success)
         class(vector_iter_t), intent(inout):: this                      !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element
         class(vector_elem_t), pointer:: vep
         integer(INTL):: l,k
         integer(INTD):: q(1:4)

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current).and.this%curr_offset.ge.0_INTL) then
           ierr=GFC_SUCCESS
           l=this%container%length()
           k=this%curr_offset+1_INTL
           if(k.lt.l) then
            call flat2quadruplet(k,q)
            vep=>this%container%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))
            if(present(elem_p)) then
             elem_p=>vep
            else
             this%curr_offset=k
             this%current=>vep
            endif
           else
            this%current=>NULL(); this%curr_offset=-1_INTL
            ierr=this%set_status_(GFC_IT_DONE)
            vep=>NULL()
           endif
           if(.not.associated(vep)) then; ierr=GFC_IT_DONE; else; vep=>NULL(); endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function VectorIterNext
!------------------------------------------------------------
        function VectorIterPrevious(this,elem_p) result(ierr)
!If <elem_p> is absent, the iterator moves to the previous element, if any.
!If <elem_p> is present, the iterator simply returns the previous element in <elem_p> without moving.
         implicit none
         integer(INTD):: ierr                                            !out: error code (0:success)
         class(vector_iter_t), intent(inout):: this                      !inout: iterator
         class(gfc_cont_elem_t), pointer, intent(out), optional:: elem_p !out: pointer to the container element
         class(vector_elem_t), pointer:: vep
         integer(INTL):: l,k
         integer(INTD):: q(1:4)

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current).and.this%curr_offset.ge.0_INTL) then
           ierr=GFC_SUCCESS
           k=this%curr_offset-1_INTL
           if(k.ge.0_INTL) then
            call flat2quadruplet(k,q)
            vep=>this%container%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))
            if(present(elem_p)) then
             elem_p=>vep
            else
             this%curr_offset=k
             this%current=>vep
            endif
           else
            this%current=>NULL(); this%curr_offset=-1_INTL
            ierr=this%set_status_(GFC_IT_DONE)
            vep=>NULL()
           endif
           if(.not.associated(vep)) then; ierr=GFC_IT_DONE; else; vep=>NULL(); endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function VectorIterPrevious
!-------------------------------------------------------------
        function VectorIterGetLength(this,ierr) result(length)
!Returns the current length of the vector.
         implicit none
         integer(INTL):: length                      !out: length
         class(vector_iter_t), intent(in):: this     !in: vector iterator
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(associated(this%container)) then
          length=this%container%length(errc)
         else
          length=-1_INTL; errc=GFC_NULL_CONT
         endif
         if(present(ierr)) ierr=errc
         return
        end function VectorIterGetLength
!-------------------------------------------------------------
        function VectorIterGetOffset(this,ierr) result(offset)
!Returns the offset of the current vector element.
         implicit none
         integer(INTL):: offset                      !out: offset
         class(vector_iter_t), intent(in):: this     !in: vector iterator
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         offset=-1_INTL; errc=this%get_status()
         if(errc.eq.GFC_IT_ACTIVE) then
          errc=GFC_SUCCESS; offset=this%curr_offset
          if(offset.lt.0_INTL) errc=GFC_IT_DONE
         endif
         if(present(ierr)) ierr=errc
         return
        end function VectorIterGetOffset
!------------------------------------------------------------------
        function VectorIterElement(this,offset,ierr) result(elem_p)
!Returns a pointer to the specific vector element.
         implicit none
         class(gfc_cont_elem_t), pointer:: elem_p    !out: pointer to the specific vector element
         class(vector_iter_t), intent(in):: this     !in: vector iterator
         integer(INTL), intent(in):: offset          !in: vector element offset
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: q(4),errc

         elem_p=>NULL(); errc=this%get_status()
         if(errc.eq.GFC_IT_ACTIVE.or.errc.eq.GFC_IT_DONE) then
          if(offset.ge.0_INTL.and.offset.lt.this%container%length()) then
           errc=GFC_SUCCESS
           call flat2quadruplet(offset,q)
           elem_p=>this%container%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))
          else
           errc=GFC_INVALID_ARGS
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end function VectorIterElement
!----------------------------------------------------------------------
        function VectorIterElementValue(this,offset,ierr) result(val_p)
!Returns a pointer to the value of a specific vector element.
         implicit none
         class(*), pointer:: val_p                   !out: pointer to the value of a specific vector element
         class(vector_iter_t), intent(in):: this     !in: vector iterator
         integer(INTL), intent(in):: offset          !in: vector element offset
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         class(gfc_cont_elem_t), pointer:: cep

         cep=>this%element(offset,errc)
         if(errc.eq.GFC_SUCCESS) val_p=>cep%get_value(errc)
         if(present(ierr)) ierr=errc
         return
        end function VectorIterElementValue
!----------------------------------------------------------
        function VectorIterMoveTo(this,offset) result(ierr)
!Moves the iterator to the given vector position.
         implicit none
         integer(INTD):: ierr                       !out: error code
         class(vector_iter_t), intent(inout):: this !in: vector iterator
         integer(INTL), intent(in):: offset         !in: offset
         integer(INTD):: q(4)

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE.or.ierr.eq.GFC_IT_DONE) then
          if(offset.ge.0_INTL.and.offset.lt.this%container%length()) then
           this%curr_offset=offset
           call flat2quadruplet(offset,q)
           this%current=>this%container%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))
           if(ierr.eq.GFC_IT_DONE) ierr=this%set_status_(GFC_IT_ACTIVE)
          else
           ierr=GFC_NO_MOVE
          endif
         endif
         return
        end function VectorIterMoveTo
!-----------------------------------------------------------------------------------
#ifdef NO_GNU
        function VectorIterAppend(this,elem_val,assoc_only,copy_ctor_f) result(ierr) !`GCC has a bug with this
#else
        function VectorIterAppend(this,elem_val,assoc_only) result(ierr)
#endif
!Appends a vector element to the end of the vector. The iterator position is kept unchanged,
!unless the container is empty in which case it will be reset to the first element.
         implicit none
         integer(INTD):: ierr                       !out: error code
         class(vector_iter_t), intent(inout):: this !inout: iterator
         class(*), target, intent(in):: elem_val    !in: value to be stored
         logical, intent(in), optional:: assoc_only !in: storage type: TRUE:by reference, FALSE:by value (default)
#ifdef NO_GNU
         procedure(gfc_copy_i), optional:: copy_ctor_f !user-defined generic copy constructor (when storing by value only)
#endif
         logical:: assoc
         integer(INTL):: offset,nelems
         integer(INTD):: q(1:4)

         if(present(assoc_only)) then; assoc=assoc_only; else; assoc=.FALSE.; endif
         ierr=this%get_status()
         if(ierr.eq.GFC_IT_DONE) then; ierr=this%reset(); ierr=this%get_status(); endif
         if(ierr.eq.GFC_IT_ACTIVE) then !non-empty container
          if(associated(this%container)) then
           offset=this%container%length()
           ierr=this%container%adjust_structure_(offset)
           if(ierr.eq.GFC_SUCCESS) then
            call flat2quadruplet(offset,q)
#ifdef NO_GNU
            if(present(copy_ctor_f)) then
             call this%container%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))%&
                 &vector_elem_ctor(elem_val,ierr,assoc_only=assoc,copy_ctor_f=copy_ctor_f)
            else
#endif
             call this%container%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))%&
                 &vector_elem_ctor(elem_val,ierr,assoc_only=assoc)
#ifdef NO_GNU
            endif
#endif
            if(ierr.eq.GFC_SUCCESS) then
             call this%container%incr_len_()
             if(this%container%num_elems_().ge.0) then !quick counting is on
              nelems=this%container%update_num_elems_(1_INTL,ierr); if(ierr.ne.GFC_SUCCESS) ierr=GFC_CORRUPTED_CONT
             endif
            endif
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         elseif(ierr.eq.GFC_IT_EMPTY) then !empty container
          if(associated(this%container)) then
           if(this%container%is_empty().eq.GFC_TRUE) then
            offset=0_INTL
            ierr=this%container%adjust_structure_(offset)
            if(ierr.eq.GFC_SUCCESS) then
             call flat2quadruplet(offset,q)
#ifdef NO_GNU
             if(present(copy_ctor_f)) then
              call this%container%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))%&
                  &vector_elem_ctor(elem_val,ierr,assoc_only=assoc,copy_ctor_f=copy_ctor_f)
             else
#endif
              call this%container%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))%&
                  &vector_elem_ctor(elem_val,ierr,assoc_only=assoc)
#ifdef NO_GNU
             endif
#endif
             if(ierr.eq.GFC_SUCCESS) then
              call this%container%incr_len_()
              ierr=this%reset() !reset the iterator to GFC_IT_ACTIVE
              if(this%container%num_elems_().ge.0) then !quick counting is on
               nelems=this%container%update_num_elems_(1_INTL,ierr); if(ierr.ne.GFC_SUCCESS) ierr=GFC_CORRUPTED_CONT
              endif
             endif
            endif
           else
            ierr=GFC_CORRUPTED_CONT
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function VectorIterAppend
!-----------------------------------------------------------------------------------
#ifdef NO_GNU
        function VectorIterInsert(this,elem_val,assoc_only,copy_ctor_f) result(ierr) !`GCC has a bug with this
#else
        function VectorIterInsert(this,elem_val,assoc_only) result(ierr)
#endif
!Inserts a vector element at the current iterator position. The iterator position is kept unchanged,
!unless the container is empty in which case it will be reset to the first element.
         implicit none
         integer(INTD):: ierr                       !out: error code
         class(vector_iter_t), intent(inout):: this !inout: iterator
         class(*), target, intent(in):: elem_val    !in: value to be stored
         logical, intent(in), optional:: assoc_only !in: storage type: TRUE:by reference, FALSE:by value (default)
#ifdef NO_GNU
         procedure(gfc_copy_i), optional:: copy_ctor_f !user-defined generic copy constructor (when storing by value only)
#endif
         logical:: assoc
         integer(INTL):: offset,nelems
         integer(INTD):: q(1:4)

         if(present(assoc_only)) then; assoc=assoc_only; else; assoc=.FALSE.; endif
         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then !non-empty container
          if(associated(this%container).and.this%curr_offset.ge.0_INTL) then
           offset=this%curr_offset
           ierr=this%container%adjust_structure_(offset)
           if(ierr.eq.GFC_SUCCESS) then
            call flat2quadruplet(offset,q)
#ifdef NO_GNU
            if(present(copy_ctor_f)) then
             call this%container%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))%&
                 &vector_elem_ctor(elem_val,ierr,assoc_only=assoc,copy_ctor_f=copy_ctor_f)
            else
#endif
             call this%container%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))%&
                 &vector_elem_ctor(elem_val,ierr,assoc_only=assoc)
#ifdef NO_GNU
            endif
#endif
            if(ierr.eq.GFC_SUCCESS) then
             this%current=>this%container%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))
             call this%container%incr_len_()
             if(this%container%num_elems_().ge.0) then !quick counting is on
              nelems=this%container%update_num_elems_(1_INTL,ierr); if(ierr.ne.GFC_SUCCESS) ierr=GFC_CORRUPTED_CONT
             endif
            endif
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         elseif(ierr.eq.GFC_IT_EMPTY) then !empty container
          if(associated(this%container)) then
           if(this%container%is_empty().eq.GFC_TRUE) then
            offset=0_INTL
            ierr=this%container%adjust_structure_(offset)
            if(ierr.eq.GFC_SUCCESS) then
             call flat2quadruplet(offset,q)
#ifdef NO_GNU
             if(present(copy_ctor_f)) then
              call this%container%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))%&
                  &vector_elem_ctor(elem_val,ierr,assoc_only=assoc,copy_ctor_f=copy_ctor_f)
             else
#endif
              call this%container%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))%&
                  &vector_elem_ctor(elem_val,ierr,assoc_only=assoc)
#ifdef NO_GNU
             endif
#endif
             if(ierr.eq.GFC_SUCCESS) then
              call this%container%incr_len_()
              ierr=this%reset() !reset the iterator to GFC_IT_ACTIVE
              if(this%container%num_elems_().ge.0) then !quick counting is on
               nelems=this%container%update_num_elems_(1_INTL,ierr); if(ierr.ne.GFC_SUCCESS) ierr=GFC_CORRUPTED_CONT
              endif
             endif
            endif
           else
            ierr=GFC_CORRUPTED_CONT
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function VectorIterInsert
!-----------------------------------------------------
        function VectorIterSwapLast(this) result(ierr)
!Swaps the currently pointed to element with the last element of the vector.
!The iterator position is unchanged.
         implicit none
         integer(INTD):: ierr                       !out: error code
         class(vector_iter_t), intent(inout):: this !inout: (active) iterator
         integer(INTL):: last
         integer(INTD):: q(1:4),w(1:4)
         type(vector_elem_t):: vec_elem

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          ierr=GFC_SUCCESS; last=this%get_length()-1_INTL
          if(last.ge.0_INTL.and.this%curr_offset.ge.0_INTL) then
           if(this%curr_offset.ne.last) then
            call flat2quadruplet(this%curr_offset,q); call flat2quadruplet(last,w)
            vec_elem=this%container%vec_tile(w(4))%tile_batch(w(3))%batch_seg(w(2))%seg_elem(w(1))
            this%container%vec_tile(w(4))%tile_batch(w(3))%batch_seg(w(2))%seg_elem(w(1))=&
            &this%container%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))
            this%container%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))=vec_elem
            this%current=>this%container%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         return
        end function VectorIterSwapLast
!--------------------------------------------------------------
        function VectorIterDelete(this,destruct_f) result(ierr)
!Deletes the element at the current iterator position.
!The current iterator position is kept unchanged, pointing to the next
!element. If there is no next element, it will point to the previous
!element. If none, then the iterator will become empty.
         integer(INTD):: ierr                             !out: error code
         class(vector_iter_t), intent(inout):: this       !inout: iterator
         procedure(gfc_destruct_i), optional:: destruct_f !in: element value destructor
         integer(INTD):: errc,q(1:4)
         integer(INTL):: l,nelems

         ierr=this%get_status()
         if(ierr.eq.GFC_IT_ACTIVE) then
          if(associated(this%current).and.this%curr_offset.ge.0_INTL.and.this%curr_offset.lt.this%get_length()) then
           !call flat2quadruplet(this%curr_offset,q); write(*,'("Deleting ",i7,4(1x,i4))') this%curr_offset,q(4:1:-1) !debug
           if(present(destruct_f)) then
            call this%current%destruct(ierr,destruct_f)
           else
            call this%current%destruct(ierr)
           endif
           if(ierr.ne.GFC_SUCCESS) ierr=NOT_CLEAN
           errc=this%container%adjust_structure_(this%curr_offset,for_deletion=.TRUE.)
           if(errc.ne.GFC_SUCCESS) ierr=GFC_ERROR
           call this%container%decr_len_(); l=this%get_length()
           if(l.gt.0_INTL) then
            if(this%curr_offset.ge.l) this%curr_offset=l-1_INTL !reposition to the last element
            call flat2quadruplet(this%curr_offset,q)
            this%current=>this%container%vec_tile(q(4))%tile_batch(q(3))%batch_seg(q(2))%seg_elem(q(1))
           else
            this%current=>NULL(); this%curr_offset=-1_INTL
            errc=this%reset()
           endif
           if(this%container%num_elems_().ge.0) then !quick counting is on
            nelems=this%container%update_num_elems_(-1_INTL,errc); if(errc.ne.GFC_SUCCESS.and.ierr.eq.GFC_SUCCESS) ierr=NOT_CLEAN
           endif
          else
           ierr=GFC_CORRUPTED_CONT
          endif
         endif
         !if(ierr.ne.GFC_SUCCESS) print *,'D: ',ierr !debug
         return
        end function VectorIterDelete
!-----------------------------------------------------------------
        function VectorIterDeleteAll(this,destruct_f) result(ierr)
!Deletes all elements of the vector, leaving it empty.
         implicit none
         integer(INTD):: ierr                             !out: error code
         class(vector_iter_t), intent(inout):: this       !inout: iterator
         procedure(gfc_destruct_i), optional:: destruct_f !in: element value destructor
         integer(INTD):: errc

         ierr=this%reset_back()
         if(ierr.eq.GFC_SUCCESS) then
          ierr=this%get_status()
          if(ierr.eq.GFC_IT_ACTIVE) then
           ierr=GFC_SUCCESS
           if(present(destruct_f)) then
            do while(this%get_length().gt.0_INTL)
             errc=this%delete(destruct_f); if(errc.ne.GFC_SUCCESS) ierr=NOT_CLEAN
            enddo
           else
            do while(this%get_length().gt.0_INTL)
             errc=this%delete(); if(errc.ne.GFC_SUCCESS) ierr=NOT_CLEAN
            enddo
           endif
          else
           if(ierr.eq.GFC_IT_EMPTY) ierr=GFC_SUCCESS
          endif
         endif
         !if(ierr.ne.GFC_SUCCESS) print *,'DA: ',ierr
         return
        end function VectorIterDeleteAll

       end module gfc_vector
!TESTING=====================
       module gfc_vector_test
        use gfc_base
        use gfc_vector
        use timers, only: thread_wtime
        implicit none
        private

        public test_gfc_vector

        type, private:: some_t
         real(8):: some_real=0d0
         integer(INTD):: some_int=0
         integer(INTL):: some_long=0
         real(8), pointer:: some_arr(:)=>NULL()
        end type some_t

      contains

       function test_gfc_vector(perf,dev_out) result(ierr)
        implicit none
        integer(INTD):: ierr
        real(8), intent(out):: perf
        integer(INTD), intent(in), optional:: dev_out
        integer(INTD), parameter:: MAX_VEC_ELEMS=1000000
        integer(INTD):: jo,i
        type(some_t):: some_val
        type(vector_t):: some_vector
        type(vector_iter_t):: some_iter
        real(8):: tms,tm

        if(present(dev_out)) then; jo=dev_out; else; jo=6; endif
        perf=0d0; tms=thread_wtime()
        ierr=some_iter%init(some_vector); if(ierr.ne.GFC_SUCCESS) then; ierr=1; return; endif
!Append elements to the vector:
        do i=1,MAX_VEC_ELEMS
         some_val%some_int=i
         ierr=some_iter%append(some_val); if(ierr.ne.GFC_SUCCESS) then; ierr=2; return; endif
        enddo
        if(some_iter%get_length(ierr).ne.MAX_VEC_ELEMS) then; ierr=3; return; endif
        if(ierr.ne.GFC_SUCCESS) then; ierr=4; return; endif
        ierr=some_iter%delete_all(); if(ierr.ne.GFC_SUCCESS) then; ierr=5; return; endif
        if(some_iter%get_length(ierr).ne.0) then; ierr=6; return; endif
        if(ierr.ne.GFC_SUCCESS) then; ierr=7; return; endif
        ierr=some_iter%release(); if(ierr.ne.GFC_SUCCESS) then; ierr=8; return; endif
        tm=thread_wtime(tms); perf=dble(MAX_VEC_ELEMS)/tm
        return
       end function test_gfc_vector

      end module gfc_vector_test
