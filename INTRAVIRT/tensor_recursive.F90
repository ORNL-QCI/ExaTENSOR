!ExaTENSOR: Recursive tensors
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/01/27

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

       module tensor_recursive
!Acronyms:
! # SAT: Subspace Aggregation Tree.
! # MUD: Maximally Uniform Distribution.
        use tensor_algebra !includes dil_basic
        use stsubs
        use gfc_base
        use gfc_list
        use gfc_tree
        use subspaces
        use distributed
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !output device
        integer(INTD), private:: DEBUG=0    !debugging level
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
 !Tensor naming:
        integer(INTD), parameter, public:: TEREC_MAX_TENSOR_NAME_LEN=32 !max length of an alphanumeric_ tensor name
 !Error codes:
        integer(INTD), parameter, public:: TEREC_SUCCESS=SUCCESS
        integer(INTD), parameter, public:: TEREC_ERROR=GENERIC_ERROR
        integer(INTD), parameter, public:: TEREC_INVALID_ARGS=-1
        integer(INTD), parameter, public:: TEREC_INVALID_REQUEST=-2
        integer(INTD), parameter, public:: TEREC_MEM_ALLOC_FAILED=-3
        integer(INTD), parameter, public:: TEREC_MEM_FREE_FAILED=-4
        integer(INTD), parameter, public:: TEREC_UNABLE_COMPLETE=-5
        integer(INTD), parameter, public:: TEREC_OBJ_CORRUPTED=-6
 !Storage layout for tensor blocks (locally stored tensors):
        integer(INTD), parameter, public:: TEREC_LAY_NONE=0  !none
        integer(INTD), parameter, public:: TEREC_LAY_RECUR=1 !storage layout is inferred from that of individual constituent tensors (recursive)
        integer(INTD), parameter, public:: TEREC_LAY_FDIMS=2 !straightforward dimension lead storage layout (Fortran style)
        integer(INTD), parameter, public:: TEREC_LAY_CDIMS=3 !straightforward dimension lead storage layout (C style)
        integer(INTD), parameter, public:: TEREC_LAY_DSYMM=4 !dimension led storage layout with permutational symmetry
        integer(INTD), parameter, public:: TEREC_LAY_BRICK=5 !bricked storage layout (all bricks of the same size, padding if needed)
        integer(INTD), parameter, public:: TEREC_LAY_BSYMM=6 !bricked storage layour with permutational symmetry on the brick level (all bricks of the same size, padding if needed)
        integer(INTD), parameter, public:: TEREC_LAY_SPARS=7 !sparse tensor storage layout
 !Index restriction kinds:
        integer(INTD), parameter, public:: TEREC_IND_RESTR_NONE=0 !no restrictions
        integer(INTD), parameter, public:: TEREC_IND_RESTR_LT=1   !indices within the group are < ordered: i1 < i2 < i3
        integer(INTD), parameter, public:: TEREC_IND_RESTR_GT=2   !indices within the group are > ordered: i1 > i2 > i3
        integer(INTD), parameter, public:: TEREC_IND_RESTR_LE=3   !indices within the group are <= ordered: i1 <= i2 <= i3
        integer(INTD), parameter, public:: TEREC_IND_RESTR_GE=4   !indices within the group are >= ordered i1 >= i2 >= i3
!TYPES:
 !Tensor signature (unique tensor identifier):
        type, public:: tens_signature_t
         character(:), allocatable, private:: char_name     !character tensor name (alphanumeric_)
         integer(INTD), private:: num_dim=-1                !number of tensor dimensions (tensor order, tensor rank)
         integer(INTL), allocatable, private:: space_idx(:) !subspace id for each tensor dimension
         contains
          !initial:: tens_signature_ctor                     !ctor
          procedure, public:: is_set=>TensSignatureIsSet     !returns .TRUE. if the tensor signature is set
          procedure, public:: get_name=>TensSignatureGetName !returns the alphanumeric_ tensor name
          procedure, public:: get_rank=>TensSignatureGetRank !returns the rank of the tensor (number of dimensions)
          procedure, public:: get_spec=>TensSignatureGetSpec !returns the tensor subspace multi-index (specification)
          final:: tens_signature_dtor                        !dtor
        end type tens_signature_t
 !Tensor shape:
        type, public:: tens_shape_t
         integer(INTL), allocatable, private:: dim_extent(:) !tensor dimension extents (resolution)
         integer(INTD), allocatable, private:: dim_group(:)  !tensor dimension groups (group 0 is default with no restrictions)
         integer(INTD), allocatable, private:: group_spec(:) !specification of the restriction kind for each index restriction group
         integer(INTD), private:: num_dim=-1                 !number of tensor dimensions (tensor order, tensor rank)
         integer(INTD), private:: num_grps=0                 !number of defined (non-trivial) index restriction groups
#if 0
         contains
          procedure, public:: set_dims=>TensShapeSetDims     !sets tensor dimension extents and tensor rank (ctor)
          procedure, public:: set_group=>TensShapeSetGroup   !creates a new index restriction group
          procedure, public:: is_set=>TensShapeIsSet         !returns .TRUE. if the tensor shape is set
          procedure, public:: get_dims=>TensShapeGetDims     !returns tensor dimension extents
          procedure, public:: get_rank=>TensShapeGetRank     !return the rank of the tensor (number of dimensions)
          procedure, public:: get_group=>TensShapeGetGroup   !returns the index restriction group number for given dimensions
          procedure, public:: same_group=>TensShapeSameGroup !checks whether specific tensor dimensions belong to the same group
          procedure, public:: num_groups=>TensShapeNumGroups !returns the total number of non-trivial index groups set on the shape
          procedure, public:: clean=>TensShapeClean          !cleans the object (dtor)
          final:: TensShapeFinal                             !dtor
#endif
        end type tens_shape_t
 !Tensor header (signature+shape):
        type, public:: tens_header_t
         type(tens_signature_t), private:: signature         !tensor signature
         type(tens_shape_t), private:: shape                 !tensor shape
#if 0
         contains
          procedure, public:: TensHeaderSetAll                       !sets the tensor header by specifying all underlying components (ctor)
          procedure, public:: TensHeaderSetParts                     !sets the tensor heeader by providing defined signature and shape (ctor)
          generic, public:: set=>TensHeaderSetAll,TensHeaderSetParts !sets the tensor header (ctor)
          procedure, public:: is_set=>TensHeaderIsSet                !returns .TRUE. if the tensor header is set
          procedure, public:: get_name=>TensHeaderGetName            !returns the alphanumeric_ tensor name
          procedure, public:: get_rank=>TensHeaderGetRank            !returns the rank of the tensor (number of dimensions)
          procedure, public:: get_spec=>TensHeaderGetSpec            !returns the tensor subspace multi-index (specification)
          procedure, public:: get_dims=>TensHeaderGetDims            !returns tensor dimension extents
          procedure, public:: get_signature=>TensHeaderGetSignature  !returns the pointer to the tensor signature
          procedure, public:: get_shape=>TensHeaderGetShape          !returns the pointer the the tensor shape
          final:: TensHeaderFinal                                    !dtor
#endif
        end type tens_header_t
 !Simple (dense) tensor block (part):
        type, public:: tens_simple_part_t
         type(tens_header_t), private:: header               !header of the constituent simple tensor block
         integer(INTL), private:: offset=-1_INTL             !offset of the constituent simple tensor block in the parental composite tensor block (locally stored)
         integer(INTD), private:: layout=TEREC_LAY_NONE      !simple storage layout: {TEREC_LAY_FDIMS,TEREC_LAY_CDIMS} only
#if 0
         contains
          procedure, public:: set=>TensSimplePartSet              !sets the tensor simple part (ctor)
          procedure, public:: is_set=>TensSimplePartIsSet         !return .TRUE. if the simple part is set
          procedure, public:: get_offset=>TensSimplePartGetOffset !returns the offset of the simple part in the parental tensor block
          procedure, public:: get_layout=>TensSimplePartGetLayout !returns the simple layout of the simple part: {TEREC_LAY_FDIMS,TEREC_LAY_CDIMS} only
          procedure, public:: get_header=>TensSimplePartGetHeader !returns a pointer to the header of the simple tensor part
          final:: TensSimplePartFinal                             !dtor
#endif
        end type tens_simple_part_t
 !Storage layout for locally stored blocks (abstract base):
        type, abstract, public:: tens_layout_t
         integer(INTD), private:: layout=TEREC_LAY_NONE      !tensor block storage layout (see above), either simple or composite
         contains
          procedure(tens_layout_map_i), deferred, public:: map                       !maps a specific element of the tensor block (layout specific)
          procedure(tens_layout_extract_i), deferred, public:: extract_simple_blocks !creates a list of constituent simple (dense) parts of the tensor block
#if 0
          procedure, public:: get_layout=>TensLayoutGetLayout                        !returns the tensor storage layout kind
#endif
        end type tens_layout_t
 !Tensor body:
        type, public:: tens_body_t
         type(DataDescr_t), allocatable, private:: data_descr !DDSS data descriptor (if physically stored as a whole)
         class(tens_layout_t), allocatable, private:: layout  !tensor block storage layout (if physically stored as a whole)
         type(list_bi_t), allocatable, private:: subtensors   !list of constituent tensors in terms of tensor headers
#if 0
         contains
          procedure, public:: set_location=>TensBodySetLocation     !sets the data location and storage layout if physically stored as a whole (builder ctor)
          procedure, public:: add_subtensor=>TensBodyAddSubtensor   !registers a constituent subtensor (builder ctor)
          procedure, public:: get_data_descr=>TensBodyGetDataDescr  !returns a pointer to the DDSS data descriptor
          procedure, public:: get_layout=>TensBodyGetLayout         !returns a pointer to the storage layout
          procedure, public:: get_subtensors=>TensBodyGetSubtensors !returns a pointer to the list of constituent subtensors (each subtensor is represented by a tensor header)
          procedure, public:: is_mapped=>TensBodyIsMapped           !returns .TRUE. if the tensor body is mapped (stored as a whole)
          procedure, public:: clean=>TensBodyClean                  !resets the object to an empty state (dtor)
          final:: TensBodyFinal                                     !dtor
#endif
        end type tens_body_t
 !Recursive tensor:
        type, public:: tens_rcrsv_t
         type(tens_header_t), private:: header                !tensor header (signature + shape)
         type(tens_body_t), private:: body                    !tensor body (data location and storage layout)
#if 0
         contains
          procedure, public:: set_header=>TensRcrsvSetHeader  !sets the tensor header (builder ctor)
          procedure, public:: set_body=>TensRcrsvSetBody      !sets the tensor body (builder ctor)
          procedure, public:: set=>TensRcrsvSet               !sets the tensor (ctor)
          procedure, public:: is_set=>TensRcrsvIsSet          !returns .TRUE. is the tensor is set
          procedure, public:: get_header=>TensRcrsvGetHeader  !returns a pointer to the tensor header
          procedure, public:: get_body=>TensRcrsvGetBody      !returns a pointer to the tensor body
          procedure, public:: clean=>TensRcrsvClean           !resets the object to an empty state (dtor)
          final:: TensRcrsvFinal                              !dtor
#endif
        end type tens_rcrsv_t
!INTERFACES:
 !Abstract:
        abstract interface
  !tens_layout_t: .map():
         function tens_layout_map_i(this,mlndx,ierr) result(offset)
          import:: INTD,INTL,tens_layout_t
          implicit none
          integer(INTL):: offset                      !out: offset of the requested tensor element
          class(tens_layout_t), intent(in):: this     !in: tensor block storage layout
          integer(INTL), intent(in):: mlndx(1:)       !in: input multi-index specifying the individual tensor element
          integer(INTD), intent(out), optional:: ierr !out: error code
         end function tens_layout_map_i
  !tens_layout_t: .extract_simple_blocks():
         subroutine tens_layout_extract_i(this,num_parts,part,ierr)
          import:: INTD,INTL,list_bi_t,tens_layout_t
          implicit none
          class(tens_layout_t), intent(in):: this     !in: tensor block storage layout
          integer(INTL), intent(out):: num_parts      !out: number of constituent simple (dense) blocks
          type(list_bi_t), intent(inout):: part       !out: list of the constituent simple (dense) blocks with their signatures and pointers
          integer(INTD), intent(out), optional:: ierr !out: error code
         end subroutine tens_layout_extract_i
        end interface
!VISIBILITY:
 !tens_signature_t:
        public tens_signature_ctor
        private TensSignatureIsSet
        private TensSignatureGetName
        private TensSignatureGetRank
        private TensSignatureGetSpec
        public tens_signature_dtor
#if 0
 !tens_shape_t:
        private TensShapeSetDims
        private TensShapeSetGroup
        private TensShapeIsSet
        private TensShapeGetDims
        private TensShapeGetRank
        private TensShapeGetGroup
        private TensShapeSameGroup
        private TensShapeNumGroups
        private TensShapeClean
 !tens_header_t:
        private TensHeaderSetAll
        private TensHeaderSetParts
        private TensHeaderIsSet
        private TensHeaderGetName
        private TensHeaderGetRank
        private TensHeaderGetSpec
        private TensHeaderGetDims
        private TensHeaderGetSignature
        private TensHeaderGetShape
 !tens_simple_part_t:
        private TensSimplePartSet
        private TensSimplePartIsSet
        private TensSimplePartGetOffset
        private TensSimplePartGetLayout
        private TensSimplePartGetHeader
 !tens_layout_t:
        private TensLayoutGetLayout
 !tens_body_t:
        private TensBodySetLocation
        private TensBodyAddSubtensor
        private TensBodyGetDataDescr
        private TensBodyGetLayout
        private TensBosyGetSubtensors
        private TensBodyIsMapped
        private TensBodyClean
 !tens_rcrsv_t:
        private TensRcrsvSetHeader
        private TensRcrsvSetBody
        private TensRcrsvSet
        private TensRcrsvGetHeader
        private TensRcrsvGetBody
        private TensRcrsvClean
#endif
!DATA:

       contains
!IMPLEMENTATION:
![tens_signature_t]--------------------------------------------------
        subroutine tens_signature_ctor(this,subspaces,ierr,tens_name)
!CTOR for tens_signature_t.
         implicit none
         type(tens_signature_t), intent(out):: this          !out: tensor signature
         integer(INTL), intent(in), optional:: subspaces(1:) !in: multi-index of subspaces
         integer(INTD), intent(out), optional:: ierr         !out: error code
         character(*), intent(in), optional:: tens_name      !in: alphanumeric_ tensor name (no spaces allowed!)
         integer(INTD):: errc,n

         errc=TEREC_SUCCESS
         n=0; if(present(subspaces)) n=size(subspaces)
         if(n.gt.0) then
          allocate(this%space_idx(1:n),STAT=errc)
          if(errc.eq.0) then
           this%num_dim=n !true tensor
          else
           errc=TEREC_MEM_ALLOC_FAILED
          endif
         else
          this%num_dim=0 !scalar
         endif
         if(errc.eq.TEREC_SUCCESS.and.present(tens_name)) then
          if(alphanumeric_string(tens_name)) then
           n=len(tens_name)
           allocate(character(len=n)::this%char_name,STAT=errc)
           if(errc.eq.0) then
            this%char_name=tens_name
           else
            errc=TEREC_MEM_ALLOC_FAILED
           endif
          else
           errc=TEREC_INVALID_ARGS
          endif
         endif
         if(errc.ne.TEREC_SUCCESS) call tens_signature_dtor(this)
         if(present(ierr)) ierr=errc
         return
        end subroutine tens_signature_ctor
!---------------------------------------------------------
        function TensSignatureIsSet(this,ierr) result(res)
!Returns TRUE of the tensor_signature_t object is set, FALSE otherwise.
         implicit none
         logical:: res                               !out: result
         class(tens_signature_t), intent(in):: this  !in: tensor signature
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         res=(this%num_dim.ge.0)
         if(present(ierr)) ierr=errc
         return
        end function TensSignatureIsSet
!--------------------------------------------------------------------
        subroutine TensSignatureGetName(this,tens_name,name_len,ierr)
!Returns the alphanumeric_ name of the tensor.
         implicit none
         class(tens_signature_t), intent(in):: this  !in: tensor signature
         character(*), intent(inout):: tens_name     !out: tensor name
         integer(INTD), intent(out):: name_len       !out: length of the tensor name
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         if(this%is_set()) then
          if(allocated(this%char_name)) then
           name_len=len(this%char_name)
           if(len(tens_name).ge.name_len) then
            tens_name(1:name_len)=this%char_name(1:name_len)
           else
            name_len=0; errc=TEREC_UNABLE_COMPLETE
           endif
          else
           errc=TEREC_INVALID_REQUEST
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensSignatureGetName
!-----------------------------------------------------------
        function TensSignatureGetRank(this,ierr) result(res)
!Returns the rank of the tensor (number of dimensions).
         implicit none
         integer(INTD):: res                         !out: result
         class(tens_signature_t), intent(in):: this  !in: tensor signature
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         if(this%is_set()) then
          res=this%num_dim
         else
          res=-1; errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensSignatureGetRank
!-------------------------------------------------------------------
        subroutine TensSignatureGetSpec(this,subspaces,num_dim,ierr)
!Returns the defining subspaces of the tensor (subspace multi-index).
         implicit none
         class(tens_signature_t), intent(in):: this   !in: tensor signature
         integer(INTL), intent(inout):: subspaces(1:) !out: defining subspaces (their IDs)
         integer(INTD), intent(out):: num_dim         !out: number of tensor dimensions
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         if(this%is_set()) then
          num_dim=this%num_dim
          if(size(subspaces).ge.num_dim) then
           subspaces(1:num_dim)=this%space_idx(1:num_dim)
          else
           num_dim=-1; errc=TEREC_UNABLE_COMPLETE
          endif
         else
          num_dim=-1; errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensSignatureGetSpec
!----------------------------------------------
        subroutine tens_signature_dtor(this)
!DTOR for tens_signature_t.
         implicit none
         type(tens_signature_t):: this

         if(allocated(this%char_name)) deallocate(this%char_name)
         if(allocated(this%space_idx)) deallocate(this%space_idx)
         this%num_dim=-1
         return
        end subroutine tens_signature_dtor

       end module tensor_recursive
!==================================
       module tensor_recursive_test
        use tensor_algebra
        use tensor_recursive
        implicit none
        private
        public test_tensor_recursive
       contains
!---------------------------------------------
        subroutine test_tensor_recursive(ierr)
         implicit none
         integer(INTD), intent(out):: ierr

         call test_tens_signature(ierr); if(ierr.ne.0) return
        end subroutine test_tensor_recursive
!-------------------------------------------
        subroutine test_tens_signature(ierr)
         implicit none
         integer(INTD), intent(out):: ierr
         integer(INTD):: errc,n
         integer(INTL):: mlndx(1:MAX_TENSOR_RANK)
         character(32):: tname
         type(tens_signature_t):: tsigna

         ierr=0
         call tens_signature_ctor(tsigna,(/3_INTL,4_INTL,2_INTL/),ierr,'Tensor')
         if(ierr.eq.TEREC_SUCCESS.and.tsigna%is_set()) then
          n=tsigna%get_rank(ierr)
          if(ierr.eq.TEREC_SUCCESS.and.n.eq.3) then
           call tsigna%get_spec(mlndx,n,ierr)
           if(ierr.eq.TEREC_SUCCESS.and.n.eq.3.and.&
             &mlndx(1).eq.3.and.mlndx(2).eq.4.and.mlndx(3).eq.2) then
            call tsigna%get_name(tname,n,ierr)
            if(tname(1:n).ne.'Tensor') ierr=4
           else
            ierr=3
           endif
          else
           ierr=2
          endif
         else
          ierr=1
         endif
         call tens_signature_dtor(tsigna)
         return
        end subroutine test_tens_signature

       end module tensor_recursive_test
