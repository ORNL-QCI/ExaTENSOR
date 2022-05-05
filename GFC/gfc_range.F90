!Generic Fortran Containers (GFC): Range
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2018/09/24

!Copyright (C) 2014-2022 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2022 Oak Ridge National Laboratory (UT-Battelle)

!LICENSE: BSD 3-Clause

       module gfc_range
        use gfc_base
        use timers
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !output device
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
        integer(INTD), private:: DEBUG=0    !debugging level (0:none)
!TYPES:
 !GFC multi-integer range:
        type, public:: gfc_range_t
         integer(INTD), private:: num_dims=0             !number of dimensions in the range (>0): [1:num_dims]
         integer(INTL), allocatable, private:: start(:)  !beginning of each dimension: [1:num_dims]
         integer(INTL), allocatable, private:: finish(:) !end of each dimension: [1:num_dims]
         integer(INTL), allocatable, private:: stride(:) !stride for each dimension (defaults to dimension extent): [0:num_dims]
         contains
          procedure, private:: GFCRangeCtor4                             !ctor for integer(4) ranges
          procedure, private:: GFCRangeCtor8                             !ctor for integer(8) ranges
          generic, public:: gfc_range_ctor=>GFCRangeCtor4,GFCRangeCtor8  !ctor
          procedure, public:: get_num_dims=>GFCRangeGetNumDims           !returns the number of dimensions in the range
          procedure, public:: get_dim_extent=>GFCRangeGetDimExtent       !returns the extent of a specific dimension
          procedure, public:: get_dim_extents=>GFCRangeGetDimExtents     !returns the extents of all dimensions
          procedure, public:: get_dim_beginning=>GFCRangeGetDimStart     !returns the beginning of a specific dimension
          procedure, public:: get_dim_beginnings=>GFCRangeGetDimStarts   !returns the beginnings of all dimensions
          procedure, public:: get_dim_end=>GFCRangeGetDimFinish          !returns the end of a specific dimension
          procedure, public:: get_dim_ends=>GFCRangeGetDimFinishes       !returns the ends of all dimensions
          procedure, public:: get_dim_stride=>GFCRangeGetDimStride       !returns the stride of a specific dimension
          procedure, public:: get_dim_strides=>GFCRangeGetDimStrides     !returns the strides of all dimensions
          procedure, public:: get_volume=>GFCRangeGetVolume              !returns the range volume (product of all extents)
          procedure, public:: get_global_volume=>GFCRangeGetGlobalVolume !returns the global volume (product of all strides)
          final:: gfc_range_dtor                                         !dtor
        end type gfc_range_t
 !GFC multi-integer range iterator:
        type, public:: gfc_range_iter_t
         class(gfc_range_t), pointer, private:: range=>NULL() !defining range
         integer(INTL), allocatable, private:: offset(:)      !current offset [0..max] for each dimension within the defining range
         integer(INTL), private:: position=-1_INTL            !linear position of the iterator within its range [0..volume-1]
         integer(INTL), private:: global_position=-1_INTL     !absolute linear position of the iterator which accounts for strides [0..max]
         contains
          procedure, public:: init=>GFCRangeIterInit                !initializes the range iterator by associating it with a range
          procedure, public:: reset=>GFCRangeIterReset              !resets the range iterator
          procedure, public:: next=>GFCRangeIterNext                !moves the range iterator to the next position
          procedure, public:: previous=>GFCRangeIterPrevious        !moves the range iterator to the previous position
          procedure, public:: get_offsets=>GFCRangeIterGetOffsets   !returns the current absolute offsets for all dimensions [start:finish]
          procedure, public:: get_position=>GFCRangeIterGetPosition !returns the current linear position of the range iterator
          procedure, public:: get_global_position=>GFCRangeIterGetGlobalPosition !returns the global linear position of the range iterator that includes strides
          procedure, public:: get_volume=>GFCRangeIterGetVolume     !returns the range volume (product of all extents)
          procedure, public:: get_global_volume=>GFCRangeIterGetGlobalVolume !returns the global volume (product of all strides)
          procedure, public:: release=>GFCRangeIterRelease          !releases the range iterator by dissociating it from its range
          final:: gfc_range_iter_dtor
        end type gfc_range_iter_t

       contains
![gfc_range_t]=================================================
        subroutine GFCRangeCtor4(this,extent,ierr,start,stride)
!Constructs a range.
         implicit none
         class(gfc_range_t), intent(out):: this           !out: range
         integer(INTD), intent(in):: extent(1:)           !in: extent of each dimension (>0)
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD), intent(in), optional:: start(1:)  !in: beginning of each dimension (defaults to zero)
         integer(INTD), intent(in), optional:: stride(1:) !in: stride for each dimension (defaults to the dimension extent)
         integer(INTD):: errc,n,i

         errc=GFC_SUCCESS; n=size(extent)
         if(n.gt.0) then
          allocate(this%start(1:n),this%finish(1:n),this%stride(0:n),STAT=errc)
          if(errc.eq.0) then
           if(present(start)) then
            this%start(1:n)=start(1:n)
           else
            this%start(1:n)=0_INTL
           endif
           this%finish(1:n)=this%start(1:n)+extent(1:n)-1
           this%stride(0)=1
           if(present(stride)) then
            this%stride(1:n)=stride(1:n)
            do i=1,n
             if(this%stride(i).lt.extent(i)) then; errc=GFC_INVALID_ARGS; exit; endif
            enddo
           else
            this%stride(1:n)=extent(1:n)
           endif
           if(errc.eq.GFC_SUCCESS) this%num_dims=n
          else
           errc=GFC_MEM_ALLOC_FAILED
          endif
         else
          errc=GFC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine GFCRangeCtor4
!--------------------------------------------------------------
        subroutine GFCRangeCtor8(this,extent,ierr,start,stride)
!Constructs a range.
         implicit none
         class(gfc_range_t), intent(out):: this           !out: range
         integer(INTL), intent(in):: extent(1:)           !in: extent of each dimension (>0)
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTL), intent(in), optional:: start(1:)  !in: beginning of each dimension (defaults to zero)
         integer(INTL), intent(in), optional:: stride(1:) !in: stride for each dimension (defaults to the dimension extent)
         integer(INTD):: errc,n,i

         errc=GFC_SUCCESS; n=size(extent)
         if(n.gt.0) then
          allocate(this%start(1:n),this%finish(1:n),this%stride(0:n),STAT=errc)
          if(errc.eq.0) then
           if(present(start)) then
            this%start(1:n)=start(1:n)
           else
            this%start(1:n)=0_INTL
           endif
           this%finish(1:n)=this%start(1:n)+extent(1:n)-1_INTL
           this%stride(0)=1_INTL
           if(present(stride)) then
            this%stride(1:n)=stride(1:n)
            do i=1,n
             if(this%stride(i).lt.extent(i)) then; errc=GFC_INVALID_ARGS; exit; endif
            enddo
           else
            this%stride(1:n)=extent(1:n)
           endif
           if(errc.eq.GFC_SUCCESS) this%num_dims=n
          else
           errc=GFC_MEM_ALLOC_FAILED
          endif
         else
          errc=GFC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine GFCRangeCtor8
!--------------------------------------------------------------
        function GFCRangeGetNumDims(this,ierr) result(num_dims)
!Returns the number of dimensions in the range, or zero if not set.
         implicit none
         integer(INTD):: num_dims                    !out: number of dimensions, or zero
         class(gfc_range_t), intent(in):: this       !in: range
         integer(INTD), intent(out), optional:: ierr !out: error code

         num_dims=this%num_dims
         if(present(ierr)) ierr=GFC_SUCCESS
         return
        end function GFCRangeGetNumDims
!--------------------------------------------------------------------
        function GFCRangeGetDimExtent(this,dimsn,ierr) result(extent)
!Returns the extent of a specific dimension of the range.
         implicit none
         integer(INTL):: extent                      !out: dimension extent
         class(gfc_range_t), intent(in):: this       !in: range
         integer(INTD), intent(in):: dimsn           !in: dimension number [1..rank]
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
         if(this%num_dims.gt.0) then
          if(dimsn.ge.1.and.dimsn.le.this%num_dims) then
           extent=this%finish(dimsn)-this%start(dimsn)+1_INTL
          else
           errc=GFC_INVALID_ARGS
          endif
         else
          errc=GFC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function GFCRangeGetDimExtent
!------------------------------------------------------------------
        subroutine GFCRangeGetDimExtents(this,extent,num_dims,ierr)
!Returns the extents of all dimensions of the range.
         implicit none
         class(gfc_range_t), intent(in):: this       !in: range
         integer(INTL), intent(inout):: extent(1:)   !out: extents
         integer(INTD), intent(out):: num_dims       !out: number of dimensions in the range
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS; num_dims=this%num_dims
         if(num_dims.gt.0) then
          extent(1:num_dims)=this%finish(1:num_dims)-this%start(1:num_dims)+1_INTL
         else
          errc=GFC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine GFCRangeGetDimExtents
!------------------------------------------------------------------
        function GFCRangeGetDimStart(this,dimsn,ierr) result(start)
!Returns the beginning of a specific dimension of the range.
         implicit none
         integer(INTL):: start                       !out: dimension beginning
         class(gfc_range_t), intent(in):: this       !in: range
         integer(INTD), intent(in):: dimsn           !in: dimension number [1..rank]
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
         if(this%num_dims.gt.0) then
          if(dimsn.ge.1.and.dimsn.le.this%num_dims) then
           start=this%start(dimsn)
          else
           errc=GFC_INVALID_ARGS
          endif
         else
          errc=GFC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function GFCRangeGetDimStart
!----------------------------------------------------------------
        subroutine GFCRangeGetDimStarts(this,start,num_dims,ierr)
!Returns the beginnings of all dimensions of the range.
         implicit none
         class(gfc_range_t), intent(in):: this       !in: range
         integer(INTL), intent(inout):: start(1:)    !out: beginnings
         integer(INTD), intent(out):: num_dims       !out: number of dimensions in the range
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS; num_dims=this%num_dims
         if(num_dims.gt.0) then
          start(1:num_dims)=this%start(1:num_dims)
         else
          errc=GFC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine GFCRangeGetDimStarts
!--------------------------------------------------------------------
        function GFCRangeGetDimFinish(this,dimsn,ierr) result(finish)
!Returns the end of a specific dimension of the range.
         implicit none
         integer(INTL):: finish                      !out: dimension end
         class(gfc_range_t), intent(in):: this       !in: range
         integer(INTD), intent(in):: dimsn           !in: dimension number [1..rank]
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
         if(this%num_dims.gt.0) then
          if(dimsn.ge.1.and.dimsn.le.this%num_dims) then
           finish=this%finish(dimsn)
          else
           errc=GFC_INVALID_ARGS
          endif
         else
          errc=GFC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function GFCRangeGetDimFinish
!-------------------------------------------------------------------
        subroutine GFCRangeGetDimFinishes(this,finish,num_dims,ierr)
!Returns the ends of all dimensions of the range.
         implicit none
         class(gfc_range_t), intent(in):: this       !in: range
         integer(INTL), intent(inout):: finish(1:)   !out: ends
         integer(INTD), intent(out):: num_dims       !out: number of dimensions in the range
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS; num_dims=this%num_dims
         if(num_dims.gt.0) then
          finish(1:num_dims)=this%finish(1:num_dims)
         else
          errc=GFC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine GFCRangeGetDimFinishes
!--------------------------------------------------------------------
        function GFCRangeGetDimStride(this,dimsn,ierr) result(stride)
!Returns the stride for a specific dimension of the range.
         implicit none
         integer(INTL):: stride                      !out: dimension stride
         class(gfc_range_t), intent(in):: this       !in: range
         integer(INTD), intent(in):: dimsn           !in: dimension number [1..rank]
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
         if(this%num_dims.gt.0) then
          if(dimsn.ge.1.and.dimsn.le.this%num_dims) then
           stride=this%stride(dimsn)
          else
           errc=GFC_INVALID_ARGS
          endif
         else
          errc=GFC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function GFCRangeGetDimStride
!------------------------------------------------------------------
        subroutine GFCRangeGetDimStrides(this,stride,num_dims,ierr)
!Returns the strides of all dimensions of the range.
         implicit none
         class(gfc_range_t), intent(in):: this       !in: range
         integer(INTL), intent(inout):: stride(1:)   !out: strides
         integer(INTD), intent(out):: num_dims       !out: number of dimensions in the range
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS; num_dims=this%num_dims
         if(num_dims.gt.0) then
          stride(1:num_dims)=this%stride(1:num_dims)
         else
          errc=GFC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine GFCRangeGetDimStrides
!--------------------------------------------------------
        function GFCRangeGetVolume(this,ierr) result(vol)
!Returns the volume of the range (total number of combinations).
         implicit none
         integer(INTL):: vol                         !out: volume
         class(gfc_range_t), intent(in):: this       !in: range
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,i

         errc=GFC_SUCCESS; vol=0_INTL
         if(this%num_dims.gt.0) then
          vol=this%finish(1)-this%start(1)+1_INTL
          do i=2,this%num_dims
           vol=vol*(this%finish(i)-this%start(i)+1_INTL)
          enddo
         else
          errc=GFC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function GFCRangeGetVolume
!--------------------------------------------------------------
        function GFCRangeGetGlobalVolume(this,ierr) result(vol)
!Returns the global volume (product of strides).
         implicit none
         integer(INTL):: vol                         !out: volume
         class(gfc_range_t), intent(in):: this       !in: range
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,i

         errc=GFC_SUCCESS; vol=0_INTL
         if(this%num_dims.gt.0) then
          vol=this%stride(1)
          do i=2,this%num_dims
           vol=vol*this%stride(i)
          enddo
         else
          errc=GFC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function GFCRangeGetGlobalVolume
!--------------------------------------
        subroutine gfc_range_dtor(this)
         implicit none
         type(gfc_range_t):: this

         if(allocated(this%stride)) deallocate(this%stride)
         if(allocated(this%finish)) deallocate(this%finish)
         if(allocated(this%start)) deallocate(this%start)
         this%num_dims=0
        end subroutine gfc_range_dtor
![gfc_range_iter_t]=================================
        subroutine GFCRangeIterInit(this,range,ierr)
!Initializes a range iterator by associating it with a range.
         implicit none
         class(gfc_range_iter_t), intent(out):: this    !out: range iterator
         class(gfc_range_t), intent(in), target:: range !in: range
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD):: errc,n

         n=range%get_num_dims(errc)
         if(errc.eq.GFC_SUCCESS.and.n.gt.0) then
          allocate(this%offset(1:n),STAT=errc)
          if(errc.eq.0) then
           this%range=>range
           this%offset(:)=0_INTL
           this%position=0_INTL
           this%global_position=0_INTL
          else
           errc=GFC_MEM_ALLOC_FAILED
          endif
         else
          if(errc.eq.GFC_SUCCESS) errc=GFC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine GFCRangeIterInit
!----------------------------------------------
        subroutine GFCRangeIterReset(this,ierr)
!Resets the range iterator.
         implicit none
         class(gfc_range_iter_t), intent(inout):: this  !inout: range iterator
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD):: errc

         errc=GFC_SUCCESS
         if(allocated(this%offset)) then
          this%offset(:)=0_INTL
          this%position=0_INTL
          this%global_position=0_INTL
         else
          errc=GFC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine GFCRangeIterReset
!---------------------------------------------------------------------
        function GFCRangeIterNext(this,pos,gl_pos,offset) result(ierr)
!Moves the range iterator to the next position. If the range iterator
!is positioned at the end of the range, GFC_NO_MOVE will be returned
!and the range iterator will change its state to finished.
!If <offset(1:)> is present, it will contain absolute offsets
!of the current range iterator position in its associated range.
         implicit none
         integer(INTD):: ierr                                !out: error code
         class(gfc_range_iter_t), intent(inout):: this       !inout: range iterator
         integer(INTL), intent(out), optional:: pos          !out: iterator position after the move
         integer(INTL), intent(out), optional:: gl_pos       !out: global iterator position after the move
         integer(INTL), intent(inout), optional:: offset(1:) !out: absolute position in each dimension of the range after the move
         integer(INTD):: n,i

         ierr=GFC_NO_MOVE
         if(associated(this%range)) then
          if(this%position.ge.0_INTL) then !valid iterator
           n=this%range%get_num_dims()
           do i=1,n
            if(this%range%start(i)+this%offset(i).lt.this%range%finish(i)) then
             this%offset(i)=this%offset(i)+1_INTL
             this%position=this%position+1_INTL
             this%global_position=this%global_position+this%range%stride(i-1)
             ierr=GFC_SUCCESS; exit
            endif
            this%global_position=this%global_position-this%offset(i)*this%range%stride(i-1)
            this%offset(i)=0_INTL
           enddo
           if(ierr.eq.GFC_SUCCESS) then
            if(present(pos)) pos=this%position
            if(present(gl_pos)) gl_pos=this%global_position
            if(present(offset)) offset(1:n)=this%range%start(1:n)+this%offset(1:n)
           else !iterator is over
            this%position=-1_INTL
            this%global_position=-1_INTL
           endif
          endif
         else
          ierr=GFC_INVALID_REQUEST
         endif
         return
        end function GFCRangeIterNext
!-------------------------------------------------------------------------
        function GFCRangeIterPrevious(this,pos,gl_pos,offset) result(ierr)
!Moves the range iterator to the previous position. If the range iterator
!is positioned at the beginning of the range, GFC_NO_MOVE will be returned
!and the range iterator will change its state to finished.
!If <offset(1:)> is present, it will contain absolute offsets
!of the current range iterator position in its associated range.
         implicit none
         integer(INTD):: ierr                                !out: error code
         class(gfc_range_iter_t), intent(inout):: this       !inout: range iterator
         integer(INTL), intent(out), optional:: pos          !out: iterator position after the move
         integer(INTL), intent(out), optional:: gl_pos       !out: global iterator position after the move
         integer(INTL), intent(inout), optional:: offset(1:) !out: absolute position in each dimension of the range after the move
         integer(INTD):: n,i

         ierr=GFC_NO_MOVE
         if(associated(this%range)) then
          if(this%position.ge.0_INTL) then !valid iterator
           n=this%range%get_num_dims()
           do i=1,n
            if(this%offset(i).gt.0_INTL) then
             this%offset(i)=this%offset(i)-1_INTL
             this%position=this%position-1_INTL
             this%global_position=this%global_position-this%range%stride(i-1)
             ierr=GFC_SUCCESS; exit
            endif
            this%offset(i)=this%range%finish(i)-this%range%start(i)
            this%global_position=this%global_position+this%offset(i)*this%range%stride(i-1)
           enddo
           if(ierr.eq.GFC_SUCCESS) then
            if(present(pos)) pos=this%position
            if(present(gl_pos)) gl_pos=this%global_position
            if(present(offset)) offset(1:n)=this%range%start(1:n)+this%offset(1:n)
           else !iterator is over
            this%position=-1_INTL
            this%global_position=-1_INTL
           endif
          endif
         else
          ierr=GFC_INVALID_REQUEST
         endif
         return
        end function GFCRangeIterPrevious
!-------------------------------------------------------------------
        subroutine GFCRangeIterGetOffsets(this,offset,num_dims,ierr)
!Returns the current absolute offsets for all range dimensions.
         implicit none
         class(gfc_range_iter_t), intent(in):: this  !in: range iterator
         integer(INTL), intent(inout):: offset(1:)   !out: absolute offsets
         integer(INTD), intent(out):: num_dims       !out: number of dimensions
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,n

         errc=GFC_SUCCESS
         if(this%position.ge.0_INTL) then !valid iterator
          num_dims=this%range%get_num_dims()
          offset(1:num_dims)=this%range%start(1:num_dims)+this%offset(1:num_dims)
         else !finished or empty iterator
          errc=GFC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine GFCRangeIterGetOffsets
!--------------------------------------------------------------
        function GFCRangeIterGetPosition(this,ierr) result(pos)
!Returns the current linear iterator position, negative in case
!the range iterator is finished or empty.
         implicit none
         integer(INTL):: pos                         !out: position
         class(gfc_range_iter_t), intent(in):: this  !in: range iterator
         integer(INTD), intent(out), optional:: ierr !out: error code

         pos=this%position
         if(present(ierr)) ierr=GFC_SUCCESS
         return
        end function GFCRangeIterGetPosition
!--------------------------------------------------------------------
        function GFCRangeIterGetGlobalPosition(this,ierr) result(pos)
!Returns the current global linear iterator position, negative in case
!the range iterator is finished or empty.
         implicit none
         integer(INTL):: pos                         !out: global position
         class(gfc_range_iter_t), intent(in):: this  !in: range iterator
         integer(INTD), intent(out), optional:: ierr !out: error code

         pos=this%global_position
         if(present(ierr)) ierr=GFC_SUCCESS
         return
        end function GFCRangeIterGetGlobalPosition
!------------------------------------------------------------
        function GFCRangeIterGetVolume(this,ierr) result(vol)
!Returns the volume of the associated range (product of dimension extents).
         implicit none
         integer(INTL):: vol                         !out: volume
         class(gfc_range_iter_t), intent(in):: this  !in: range iterator
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(associated(this%range)) then
          vol=this%range%get_volume(errc)
         else
          errc=GFC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function GFCRangeIterGetVolume
!------------------------------------------------------------------
        function GFCRangeIterGetGlobalVolume(this,ierr) result(vol)
!Returns the global volume of the associated range (product of dimension strides).
         implicit none
         integer(INTL):: vol                         !out: volume
         class(gfc_range_iter_t), intent(in):: this  !in: range iterator
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(associated(this%range)) then
          vol=this%range%get_global_volume(errc)
         else
          errc=GFC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function GFCRangeIterGetGlobalVolume
!------------------------------------------------
        subroutine GFCRangeIterRelease(this,ierr)
!Releases the range iterator by dissociating it from its range.
         implicit none
         class(gfc_range_iter_t), intent(inout):: this !inout: range iterator
         integer(INTD), intent(out), optional:: ierr   !out: error code

         this%position=-1_INTL
         this%global_position=-1_INTL
         if(allocated(this%offset)) deallocate(this%offset)
         this%range=>NULL()
         if(present(ierr)) ierr=GFC_SUCCESS
         return
        end subroutine GFCRangeIterRelease
!-------------------------------------------
        subroutine gfc_range_iter_dtor(this)
         implicit none
         type(gfc_range_iter_t):: this

         call this%release()
        end subroutine gfc_range_iter_dtor

       end module gfc_range
