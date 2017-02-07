!ExaTENSOR: Recursive tensors
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/02/07

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
 !Storage layout for tensor blocks for locally stored tensors (must span a contiguous integer range!):
        integer(INTD), parameter, public:: TEREC_LAY_NONE=0    !none
        integer(INTD), parameter, public:: TEREC_LAY_RECUR=1   !storage layout is inferred from that of individual constituent tensors (recursive)
        integer(INTD), parameter, public:: TEREC_LAY_FDIMS=2   !straightforward dimension lead storage layout (Fortran style)
        integer(INTD), parameter, public:: TEREC_LAY_CDIMS=3   !straightforward dimension lead storage layout (C style)
        integer(INTD), parameter, public:: TEREC_LAY_DSYMM=4   !dimension led storage layout with permutational symmetry
        integer(INTD), parameter, public:: TEREC_LAY_BRICK=5   !bricked storage layout (all bricks of the same size, padding if needed)
        integer(INTD), parameter, public:: TEREC_LAY_BSYMM=6   !bricked storage layour with permutational symmetry on the brick level (all bricks of the same size, padding if needed)
        integer(INTD), parameter, public:: TEREC_LAY_SPARS=7   !sparse tensor storage layout
        integer(INTD), parameter, public:: TEREC_NUM_LAYOUTS=8 !total number of tensor layouts
 !Index restriction kinds (must span a contiguous integer range!):
        integer(INTD), parameter, public:: TEREC_IND_RESTR_NONE=0 !no restrictions
        integer(INTD), parameter, public:: TEREC_IND_RESTR_LT=1   !indices within the group are < ordered: i1 < i2 < i3
        integer(INTD), parameter, public:: TEREC_IND_RESTR_GT=2   !indices within the group are > ordered: i1 > i2 > i3
        integer(INTD), parameter, public:: TEREC_IND_RESTR_LE=3   !indices within the group are <= ordered: i1 <= i2 <= i3
        integer(INTD), parameter, public:: TEREC_IND_RESTR_GE=4   !indices within the group are >= ordered i1 >= i2 >= i3
        integer(INTD), parameter, public:: TEREC_NUM_IND_RESTR=5  !total number of index restrictions (0..max)
!TYPES:
 !Tensor signature (unique tensor identifier):
        type, public:: tens_signature_t
         character(:), allocatable, private:: char_name         !character tensor name (alphanumeric_)
         integer(INTD), private:: num_dims=-1                   !number of tensor dimensions (tensor order, tensor rank)
         integer(INTL), allocatable, private:: space_idx(:)     !subspace id for each tensor dimension
         class(h_space_t), pointer, private:: h_space_p=>NULL() !pointer to the underlying hierarchical vector space
         contains
          !initial:: tens_signature_ctor                     !ctor
          procedure, public:: is_set=>TensSignatureIsSet     !returns .TRUE. if the tensor signature is set
          procedure, public:: get_name=>TensSignatureGetName !returns the alphanumeric_ tensor name
          procedure, public:: get_rank=>TensSignatureGetRank !returns the rank of the tensor (number of dimensions)
          procedure, public:: get_spec=>TensSignatureGetSpec !returns the tensor subspace multi-index (specification)
          procedure, public:: compare=>TensSignatureCompare  !compares the tensor signature with another tensor signature
          procedure, public:: print_it=>TensSignaturePrintIt !prints the tensor signature
          final:: tens_signature_dtor                        !dtor
        end type tens_signature_t
 !Tensor shape:
        type, public:: tens_shape_t
         integer(INTL), allocatable, private:: dim_extent(:) !tensor dimension extents (resolution)
         integer(INTD), allocatable, private:: dim_group(:)  !tensor dimension groups (group 0 is default with no restrictions)
         integer(INTD), allocatable, private:: group_spec(:) !specification of the restriction kind for each index restriction group (not every defined group needs to be present in dim_group(:))
         integer(INTD), private:: num_dims=-1                !number of tensor dimensions (tensor order, tensor rank)
         integer(INTD), private:: num_grps=0                 !number of defined (non-trivial) index restriction groups
         contains
          !initial:: tens_shape_ctor                         !ctor
          procedure, public:: set_groups=>TensShapeSetGroups !creates a new index restriction group
          procedure, public:: is_set=>TensShapeIsSet         !returns .TRUE. if the tensor shape is set
          procedure, public:: get_dims=>TensShapeGetDims     !returns tensor dimension extents
          procedure, public:: get_rank=>TensShapeGetRank     !return the rank of the tensor (number of dimensions)
          procedure, public:: get_group=>TensShapeGetGroup   !returns a restricted index group (specific dimensions belonging to the specified group)
          procedure, public:: same_group=>TensShapeSameGroup !checks whether specific tensor dimensions belong to the same group
          procedure, public:: num_groups=>TensShapeNumGroups !returns the total number of non-trivial index groups defined in the tensor shape
          procedure, public:: compare=>TensShapeCompare      !compares the tensor shape with another tensor shape
          procedure, public:: print_it=>TensShapePrintIt     !prints the tensor shape
          final:: tens_shape_dtor                            !dtor
        end type tens_shape_t
 !Tensor header (signature+shape):
        type, public:: tens_header_t
         type(tens_signature_t), private:: signature         !tensor signature
         type(tens_shape_t), private:: shape                 !tensor shape
         contains
          !initial:: tens_header_ctor                                !ctor
          procedure, public:: add_shape=>TensHeaderAddShape          !ctor for a deferred tensor shape specification
          procedure, public:: is_set=>TensHeaderIsSet                !returns .TRUE. if the tensor header is set (with or without shape)
          procedure, public:: get_name=>TensHeaderGetName            !returns the alphanumeric_ tensor name
          procedure, public:: get_rank=>TensHeaderGetRank            !returns the rank of the tensor (number of dimensions)
          procedure, public:: get_spec=>TensHeaderGetSpec            !returns the tensor subspace multi-index (specification)
          procedure, public:: get_dims=>TensHeaderGetDims            !returns tensor dimension extents
          procedure, public:: num_groups=>TensHeaderNumGroups        !returns the total number of non-trivial index groups defined in the tensor shape
          procedure, public:: get_group=>TensHeaderGetGroup          !returns a restricted index group (specific dimensions belonging to the specified group)
          procedure, public:: same_group=>TensHeaderSameGroup        !checks whether specific tensor dimensions belong to the same group
          procedure, public:: get_signature=>TensHeaderGetSignature  !returns the pointer to the tensor signature
          procedure, public:: get_shape=>TensHeaderGetShape          !returns the pointer the the tensor shape
          procedure, public:: compare=>TensHeaderCompare             !compares the tensor header with another tensor header
          procedure, public:: print_it=>TensHeaderPrintIt            !prints the tensor header
#ifdef NO_GNU
          final:: tens_header_dtor                                   !dtor `GCC/5.4.0 bug
#endif
        end type tens_header_t
 !Simple (dense) tensor block (part):
        type, public:: tens_simple_part_t
         type(tens_header_t), private:: header               !header of the constituent simple tensor block
         integer(INTL), private:: offset=-1_INTL             !offset of the constituent simple tensor block in the parental composite tensor block (locally stored)
         integer(INTD), private:: layout=TEREC_LAY_NONE      !simple storage layout: {TEREC_LAY_FDIMS,TEREC_LAY_CDIMS} only
         contains
          !initial:: tens_simple_part_ctor                        !ctor
          procedure, public:: is_set=>TensSimplePartIsSet         !return .TRUE. if the simple part is set (signature, shape, layout, offset)
          procedure, public:: get_offset=>TensSimplePartGetOffset !returns the offset of the simple part in the parental tensor block
          procedure, public:: get_layout=>TensSimplePartGetLayout !returns the simple layout of the simple tensor part: {TEREC_LAY_FDIMS,TEREC_LAY_CDIMS} only
          procedure, public:: get_header=>TensSimplePartGetHeader !returns a pointer to the header of the simple tensor part
          final:: tens_simple_part_dtor                           !dtor
        end type tens_simple_part_t
 !Storage layout for locally stored blocks (abstract base):
        type, abstract, public:: tens_layout_t
         integer(INTD), private:: layout=TEREC_LAY_NONE      !tensor block storage layout (see above), either simple or composite
         contains
          procedure(tens_layout_map_i), deferred, public:: map                       !maps a specific element of the tensor block (layout specific)
          procedure(tens_layout_extract_i), deferred, public:: extract_simple_blocks !creates a list of constituent simple (dense) parts of the tensor block
          procedure, public:: is_set=>TensLayoutIsSet                                !returns TRUE if the tensor layout is set
          procedure, public:: get_layout=>TensLayoutGetLayout                        !returns the tensor storage layout kind
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
        public valid_tensor_layout
        public cmp_tens_headers
 !tens_signature_t:
        public tens_signature_ctor
        private TensSignatureIsSet
        private TensSignatureGetName
        private TensSignatureGetRank
        private TensSignatureGetSpec
        private TensSignatureCompare
        private TensSignaturePrintIt
        public tens_signature_dtor
 !tens_shape_t:
        public tens_shape_ctor
        private TensShapeSetGroups
        private TensShapeIsSet
        private TensShapeGetDims
        private TensShapeGetRank
        private TensShapeGetGroup
        private TensShapeSameGroup
        private TensShapeNumGroups
        private TensShapeCompare
        private TensShapePrintIt
        public tens_shape_dtor
 !tens_header_t:
        public tens_header_ctor
        private TensHeaderAddShape
        private TensHeaderIsSet
        private TensHeaderGetName
        private TensHeaderGetRank
        private TensHeaderGetSpec
        private TensHeaderGetDims
        private TensHeaderNumGroups
        private TensHeaderGetGroup
        private TensHeaderSameGroup
        private TensHeaderGetSignature
        private TensHeaderGetShape
        private TensHeaderCompare
        private TensHeaderPrintIt
        public tens_header_dtor
 !tens_simple_part_t:
        public tens_simple_part_ctor
        private TensSimplePartIsSet
        private TensSimplePartGetOffset
        private TensSimplePartGetLayout
        private TensSimplePartGetHeader
        public tens_simple_part_dtor
 !tens_layout_t:
        private TensLayoutIsSet
        private TensLayoutGetLayout
#if 0
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
!=======================================================
        function valid_tensor_layout(layout) result(res)
!Returns TRUE if the tensor layout is valid.
         logical:: res
         integer(INTD), intent(in):: layout

         res=(layout.ge.0.and.layout.lt.TEREC_NUM_LAYOUTS)
         return
        end function valid_tensor_layout
!-----------------------------------------------------
        function cmp_tens_headers(th1,th2) result(cmp)
!Compares two tensor headers.
         implicit none
         integer(INTD):: cmp                !out: result of comparison (see GFC/gfc_base.F90)
         class(*), intent(in), target:: th1 !in: tensor header 1
         class(*), intent(in), target:: th2 !in: tensor header 2
         class(tens_header_t), pointer:: thp1,thp2

         thp1=>NULL(); thp2=>NULL()
         select type(th1); class is(tens_header_t); thp1=>th1; end select
         select type(th2); class is(tens_header_t); thp2=>th2; end select
         if(associated(thp1).and.associated(thp2)) then
          cmp=thp1%compare(thp2)
         else
          cmp=GFC_CMP_ERR
         endif
         return
        end function cmp_tens_headers
![tens_signature_t]==========================================================
        subroutine tens_signature_ctor(this,ierr,subspaces,tens_name,h_space)
!CTOR for tens_signature_t.
         implicit none
         type(tens_signature_t), intent(out):: this               !out: tensor signature
         integer(INTD), intent(out), optional:: ierr              !out: error code
         integer(INTL), intent(in), optional:: subspaces(1:)      !in: multi-index of subspaces
         character(*), intent(in), optional:: tens_name           !in: alphanumeric_ tensor name (no spaces allowed!)
         class(h_space_t), intent(in), target, optional:: h_space !in: underlying hierarchical vector space
         integer(INTD):: errc,n

         errc=TEREC_SUCCESS
         n=0; if(present(subspaces)) n=size(subspaces)
         if(n.gt.0) then
          if(present(h_space)) then
           allocate(this%space_idx(1:n),STAT=errc)
           if(errc.eq.0) then
            this%space_idx(1:n)=subspaces(1:n)
            this%num_dims=n !true tensor
            this%h_space_p=>h_space
           else
            errc=TEREC_MEM_ALLOC_FAILED
           endif
          else
           errc=TEREC_INVALID_ARGS
          endif
         else
          this%num_dims=0 !scalar
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
!------------------------------------------------------------------
        function TensSignatureIsSet(this,ierr,num_dims) result(res)
!Returns TRUE if the tensor_signature_t object is set, FALSE otherwise.
         implicit none
         logical:: res                                   !out: result
         class(tens_signature_t), intent(in):: this      !in: tensor signature
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD), intent(out), optional:: num_dims !out: tensor rank (if set)
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         res=(this%num_dims.ge.0)
         if(present(num_dims)) num_dims=this%num_dims
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
          res=this%num_dims
         else
          res=-1; errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensSignatureGetRank
!--------------------------------------------------------------------
        subroutine TensSignatureGetSpec(this,subspaces,num_dims,ierr)
!Returns the defining subspaces of the tensor (subspace multi-index).
         implicit none
         class(tens_signature_t), intent(in):: this   !in: tensor signature
         integer(INTL), intent(inout):: subspaces(1:) !out: defining subspaces (their IDs)
         integer(INTD), intent(out):: num_dims        !out: number of tensor dimensions
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         if(this%is_set()) then
          num_dims=this%num_dims
          if(size(subspaces).ge.num_dims) then
           subspaces(1:num_dims)=this%space_idx(1:num_dims)
          else
           errc=TEREC_UNABLE_COMPLETE
          endif
         else
          num_dims=-1; errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensSignatureGetSpec
!--------------------------------------------------------------
        function TensSignatureCompare(this,another) result(cmp)
!Compares the tensor signature with another tensor signature.
         implicit none
         integer(INTD):: cmp                           !out: comparison result: {GFC_CMP_EQ,GFC_CMP_LT,GFC_CMP_GT,GFC_CMP_PARENT,GFC_CMP_CHILD,GFC_CMP_NA,GFC_CMP_ERR}
         class(tens_signature_t), intent(in):: this    !in: tensor signature 1
         class(tens_signature_t), intent(in):: another !in: tensor signature 2
         integer(INTD):: nl1,nl2,ch1,ch2,i,rel,errc
         integer(INTL):: s1,s2
         class(subspace_t), pointer:: sb1,sb2

         errc=0
         if(this%is_set().and.another%is_set()) then
          cmp=GFC_CMP_EQ
!Compare names:
          nl1=len(this%char_name); nl2=len(another%char_name)
          if(nl1.lt.nl2) then; cmp=GFC_CMP_LT; elseif(nl1.gt.nl2) then; cmp=GFC_CMP_GT; endif
          if(cmp.eq.GFC_CMP_EQ) then
           do i=1,nl1
            ch1=iachar(this%char_name(i:i))
            ch2=iachar(another%char_name(i:i))
            if(ch1.lt.ch2) then; cmp=GFC_CMP_LT; exit; elseif(ch1.gt.ch2) then; cmp=GFC_CMP_GT; exit; endif
           enddo
           if(cmp.eq.GFC_CMP_EQ) then
            if(associated(this%h_space_p,another%h_space_p)) then
!Compare specs:
             if(this%num_dims.lt.another%num_dims) then
              cmp=GFC_CMP_LT
             elseif(this%num_dims.gt.another%num_dims) then
              cmp=GFC_CMP_GT
             else
              do i=1,this%num_dims
               s1=this%space_idx(i); s2=another%space_idx(i)
               sb1=>this%h_space_p%get_subspace(s1,errc); if(errc.ne.0) exit
               sb2=>this%h_space_p%get_subspace(s2,errc); if(errc.ne.0) exit
               rel=sb1%relate(sb2,this%h_space_p)
               if(rel.eq.GFC_CMP_LT.or.rel.eq.GFC_CMP_GT) then
                cmp=rel; exit
               else
                if(rel.eq.GFC_CMP_PARENT.or.rel.eq.GFC_CMP_CHILD) then
                 if(cmp.eq.GFC_CMP_EQ) then
                  cmp=rel
                 else
                  if(rel.ne.cmp) then; cmp=GFC_CMP_NA; exit; endif
                 endif
                else
                 if(rel.ne.GFC_CMP_EQ) then; cmp=rel; exit; endif
                endif
               endif
              enddo
              if(errc.ne.0) cmp=GFC_CMP_ERR
              sb1=>NULL(); sb2=>NULL()
             endif
            else
             cmp=GFC_CMP_ERR !tensors with the same name cannot reside in differe hierarchical spaces
            endif
           endif
          endif
         else
          cmp=GFC_CMP_ERR
         endif
         return
        end function TensSignatureCompare
!----------------------------------------------------------------
        subroutine TensSignaturePrintIt(this,ierr,dev_id,nspaces)
!Prints the tensor signature.
         implicit none
         class(tens_signature_t), intent(in):: this    !in: tensor signature
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD), intent(in), optional:: dev_id  !in: output device handle (6:screen)
         integer(INTD), intent(in), optional:: nspaces !in: number of preceding spaces (left alignment)
         integer(INTD):: errc,i,dev

         errc=TEREC_SUCCESS
         if(present(dev_id)) then; dev=dev_id; else; dev=6; endif
         if(present(nspaces)) then; do i=1,nspaces; write(dev,'(" ")',ADVANCE='NO'); enddo; endif
         if(this%is_set()) then
          call printl(dev,this%char_name,adv=.FALSE.)
          write(dev,'("(")',ADVANCE='NO')
          do i=1,this%num_dims
           write(dev,'(i9)',ADVANCE='NO') this%space_idx(i)
           if(i.lt.this%num_dims) write(dev,'(",")',ADVANCE='NO')
          enddo
          write(dev,'(")")')
         else
          write(dev,'("EMPTY TENSOR SIGNATURE")')
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensSignaturePrintIt
!-------------------------------------------
        subroutine tens_signature_dtor(this)
!DTOR for tens_signature_t.
         implicit none
         type(tens_signature_t):: this

         !write(*,'("#DEBUG: Entered tens_signature_dtor")') !debug
         if(allocated(this%char_name)) deallocate(this%char_name)
         if(allocated(this%space_idx)) deallocate(this%space_idx)
         this%num_dims=-1; this%h_space_p=>NULL()
         !write(*,'("#DEBUG: Exited tens_signature_dtor")') !debug
         return
        end subroutine tens_signature_dtor
![tens_shape_t]==============================================================
        subroutine tens_shape_ctor(this,ierr,dim_extent,dim_group,group_spec)
!CTOR for tens_shape_t.
         implicit none
         type(tens_shape_t), intent(out):: this               !out: tensor shape
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTL), intent(in), optional:: dim_extent(1:) !in: tensor dimension extents
         integer(INTD), intent(in), optional:: dim_group(1:)  !in: dimension grouping: dim_group(x)=y means dimension x belongs to group y>0 (group 0 is default)
         integer(INTD), intent(in), optional:: group_spec(1:) !in: group specification: group_spec(x)=y means group x has restriction y (see on top)
         integer(INTD):: errc,i,j,k,m,n
         logical:: pr1,pr2

         errc=TEREC_SUCCESS
         n=0; if(present(dim_extent)) n=size(dim_extent)
         if(n.gt.0) then !true tensor
          do i=1,n; if(dim_extent(i).le.0_INTL) then; errc=TEREC_INVALID_ARGS; exit; endif; enddo
          if(errc.eq.TEREC_SUCCESS) then
           allocate(this%dim_extent(1:n),STAT=errc)
           if(errc.eq.0) then
            this%dim_extent(1:n)=dim_extent(1:n)
            this%num_dims=n
           else
            errc=TEREC_MEM_ALLOC_FAILED
           endif
          endif
         else !scalar
          this%num_dims=0
         endif
         if(errc.eq.TEREC_SUCCESS) then
          this%num_grps=0
          pr1=present(dim_group); pr2=present(group_spec)
          if(pr1.and.pr2) then
           if(size(dim_group).eq.n) then
            m=size(group_spec) !total number of non-trivial index groups
            allocate(this%dim_group(1:n),STAT=errc)
            if(errc.eq.0) then
             allocate(this%group_spec(1:m),STAT=errc)
             if(errc.eq.0) then
              this%group_spec(1:m)=group_spec(1:m)
              this%num_grps=m
              do i=1,n
               j=dim_group(i)
               if(j.gt.0.and.j.le.m) then
                k=group_spec(j) !restriction kind (see top)
                if(k.lt.0.or.k.ge.TEREC_NUM_IND_RESTR) then; errc=TEREC_INVALID_ARGS; exit; endif
               else
                if(j.ne.0) then; errc=TEREC_INVALID_ARGS; exit; endif
               endif
               this%dim_group(i)=j
              enddo
             else
              errc=TEREC_MEM_ALLOC_FAILED
             endif
            else
             errc=TEREC_MEM_ALLOC_FAILED
            endif
           else
            errc=TEREC_INVALID_ARGS
           endif
          else
           if(pr1.or.pr2) errc=TEREC_INVALID_ARGS
          endif
         endif
         if(errc.ne.TEREC_SUCCESS) call tens_shape_dtor(this)
         if(present(ierr)) ierr=errc
         return
        end subroutine tens_shape_ctor
!--------------------------------------------------------------------
        subroutine TensShapeSetGroups(this,dim_group,group_spec,ierr)
!Sets group restrictions on tensor dimensions. If group restrictions
!have already been previously set, it will return an error.
         implicit none
         class(tens_shape_t), intent(inout):: this   !inout: tensor shape
         integer(INTD), intent(in):: dim_group(1:)   !in: dimension groups: dim_group(x)=y means dimension x belongs to group y>0 (0 is the default group)
         integer(INTD), intent(in):: group_spec(1:)  !in: group specification: group_spec(x)=y means group x has restriction y (see on top)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,i,j,k,n,m

         errc=TEREC_SUCCESS
         if(this%is_set()) then
          n=this%num_dims
          if(n.gt.0) then
           if(.not.allocated(this%group_spec)) then !groups have not been previously set up
            if(size(dim_group).eq.n) then
             m=size(group_spec)
             do i=1,n
              j=dim_group(i)
              if(j.gt.0.and.j.le.m) then
               k=group_spec(j)
               if(k.lt.0.or.k.ge.TEREC_NUM_IND_RESTR) then; errc=TEREC_INVALID_ARGS; exit; endif
              else
               if(j.ne.0) then; errc=TEREC_INVALID_ARGS; exit; endif
              endif
             enddo
             if(.not.allocated(this%dim_group)) then
              allocate(this%dim_group(1:n),STAT=errc)
              if(errc.eq.0) then; this%dim_group(1:n)=0; else; errc=TEREC_MEM_ALLOC_FAILED; endif
             endif
             if(errc.eq.TEREC_SUCCESS) then
              allocate(this%group_spec(1:m),STAT=errc)
              if(errc.eq.0) then
               this%dim_group(1:n)=dim_group(1:n)
               this%group_spec(1:m)=group_spec(1:m)
               this%num_grps=m
              else
               errc=TEREC_MEM_ALLOC_FAILED
              endif
             endif
            else
             errc=TEREC_INVALID_ARGS
            endif
           else
            errc=TEREC_INVALID_REQUEST
           endif
          else
           errc=TEREC_INVALID_REQUEST
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensShapeSetGroups
!------------------------------------------------------------------------
        function TensShapeIsSet(this,ierr,num_dims,num_groups) result(res)
!Returns TRUE if the tensor shape is set.
         implicit none
         logical:: res                                     !out: result
         class(tens_shape_t), intent(in):: this            !in: tensor shape
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD), intent(out), optional:: num_dims   !out: number of dimensions
         integer(INTD), intent(out), optional:: num_groups !out: number of dimension groups
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         res=(this%num_dims.ge.0)
         if(present(num_dims)) num_dims=this%num_dims
         if(present(num_groups)) num_groups=this%num_grps
         if(present(ierr)) ierr=errc
         return
        end function TensShapeIsSet
!-----------------------------------------------------------
        subroutine TensShapeGetDims(this,dims,num_dims,ierr)
!Returns tensor dimension extents together with the tensor rank.
         implicit none
         class(tens_shape_t), intent(in):: this      !in: tensor shape
         integer(INTL), intent(inout):: dims(1:)     !out: tensor dimension extents
         integer(INTD), intent(out):: num_dims       !out: number of tensor dimensions
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         if(this%is_set(num_dims=num_dims)) then
          if(num_dims.gt.0) then
           if(size(dims).ge.num_dims) then
            dims(1:num_dims)=this%dim_extent(1:num_dims)
           else
            errc=TEREC_UNABLE_COMPLETE
           endif
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensShapeGetDims
!-------------------------------------------------------
        function TensShapeGetRank(this,ierr) result(res)
!Returns the tensor rank.
         implicit none
         integer(INTD):: res                         !out: tensor rank
         class(tens_shape_t), intent(in):: this      !in: tensor shape
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         if(.not.this%is_set(num_dims=res)) errc=TEREC_INVALID_REQUEST
         if(present(ierr)) ierr=errc
         return
        end function TensShapeGetRank
!------------------------------------------------------------------------------------
        subroutine TensShapeGetGroup(this,group,group_dims,num_dims,ierr,group_restr)
!Returns the index restriction group (specific dimensions belonging to the specified group).
         implicit none
         class(tens_shape_t), intent(in):: this             !in: tensor shape
         integer(INTD), intent(in):: group                  !in: requested dimension group
         integer(INTD), intent(inout):: group_dims(1:)      !out: dimensions which belong to the requested dimension group
         integer(INTD), intent(out):: num_dims              !out: number of dimensions in the group
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD), intent(out), optional:: group_restr !out: group restriction kind
         integer(INTD):: errc,i,n,m

         errc=TEREC_SUCCESS; num_dims=0
         if(this%is_set(num_dims=n)) then
          if(group.ge.0.and.group.le.this%num_grps) then
           if(present(group_restr)) then
            if(group.gt.0) then
             group_restr=this%group_spec(group)
            else
             group_restr=TEREC_IND_RESTR_NONE !group 0 is default with no restrictions
            endif
           endif
           m=size(group_dims)
           do i=1,n
            if(this%dim_group(i).eq.group) then
             num_dims=num_dims+1
             if(num_dims.le.m) then
              group_dims(num_dims)=i
             else
              errc=TEREC_UNABLE_COMPLETE; exit
             endif
            endif
           enddo
          else
           errc=TEREC_INVALID_ARGS
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensShapeGetGroup
!----------------------------------------------------------------------------
        function TensShapeSameGroup(this,dims,ierr,group_restr) result(group)
!Returns the group number if the specified tensor dimensions belong to
!the same group, -1 otherwise.
         implicit none
         integer(INTD):: group                              !out: group number (>=0)
         class(tens_shape_t), intent(in):: this             !in: tensor shape
         integer(INTD), intent(in):: dims(1:)               !in: tensor dimensions to check
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD), intent(out), optional:: group_restr !out: group restriction kind
         integer(INTD):: errc,i,m,n

         errc=TEREC_SUCCESS; group=-1
         if(this%is_set(num_dims=n)) then
          m=size(dims)
          if(m.le.n) then
           if(m.gt.0) then
            do i=1,m; if(dims(i).lt.1.or.dims(i).gt.n) then; errc=TEREC_INVALID_ARGS; exit; endif; enddo
            if(errc.eq.TEREC_SUCCESS) then
             if(this%num_grps.gt.0) then !non-trivial grouping present
              group=this%dim_group(dims(1))
              do i=2,m
               if(this%dim_group(dims(i)).ne.group) then; group=-1; exit; endif
              enddo
              if(present(group_restr)) then
               if(group.eq.0) then
                group_restr=TEREC_IND_RESTR_NONE
               elseif(group.gt.0) then
                if(group.le.this%num_grps) then
                 group_restr=this%group_spec(group)
                else
                 errc=TEREC_OBJ_CORRUPTED
                endif
               endif
              endif
             else !only default group 0 (no restrictions)
              group=0; if(present(group_restr)) group_restr=TEREC_IND_RESTR_NONE
             endif
            endif
           else !dims(:) is empty
            if(n.eq.0) then !scalar (0-dimensional tensor)
             group=0; if(present(group_restr)) group_restr=TEREC_IND_RESTR_NONE
            else
             errc=TEREC_INVALID_ARGS
            endif
           endif
          else
           errc=TEREC_INVALID_ARGS
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensShapeSameGroup
!----------------------------------------------------------------
        function TensShapeNumGroups(this,ierr) result(num_groups)
!Returns the total number of non-trivial index groups defined in the tensor shape.
         implicit none
         integer(INTD):: num_groups                  !out: number of non-trivial dimension groups
         class(tens_shape_t), intent(in):: this      !in: tensor shape
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS; num_groups=0
         if(.not.this%is_set(num_groups=num_groups)) errc=TEREC_INVALID_REQUEST
         if(present(ierr)) ierr=errc
         return
        end function TensShapeNumGroups
!-------------------------------------------------------------------------
        function TensShapeCompare(this,another,compare_groups) result(cmp)
!Compares the given tensor shape with another tensor shape.
         implicit none
         integer(INTD):: cmp                            !out: comparison result: {GFC_CMP_EQ,GFC_CMP_LT,GFC_CMP_GT,GFC_CMP_ERR}
         class(tens_shape_t), intent(in):: this         !in: tensor shape 1
         class(tens_shape_t), intent(in):: another      !in: tensor shape 2
         logical, intent(in), optional:: compare_groups !in: if FALSE, dimension groups will not be taken into account (defaults to TRUE)
         integer(INTD):: i,g1,g2,gmap(MAX_TENSOR_RANK)
         logical:: comp_grps

         if(present(compare_groups)) then; comp_grps=compare_groups; else; comp_grps=.TRUE.; endif
         if(this%is_set().and.another%is_set()) then
          if(this%num_dims.lt.another%num_dims) then
           cmp=GFC_CMP_LT
          elseif(this%num_dims.gt.another%num_dims) then
           cmp=GFC_CMP_GT
          else
           cmp=GFC_CMP_EQ
           do i=1,this%num_dims
            if(this%dim_extent(i).lt.another%dim_extent(i)) then
             cmp=GFC_CMP_LT; exit
            elseif(this%dim_extent(i).gt.another%dim_extent(i)) then
             cmp=GFC_CMP_GT; exit
            endif
           enddo
           if(cmp.eq.GFC_CMP_EQ.and.comp_grps) then
            if(this%num_grps.lt.another%num_grps) then
             cmp=GFC_CMP_LT
            elseif(this%num_grps.gt.another%num_grps) then
             cmp=GFC_CMP_GT
            else
             if(this%num_grps.gt.0) then
              if(this%num_grps.le.MAX_TENSOR_RANK) then
               gmap(1:this%num_grps)=0
               do i=1,this%num_dims
                g1=this%dim_group(i); g2=another%dim_group(i)
                if(g1.gt.0.and.g2.gt.0) then !both groups are non-trivial
                 if(this%group_spec(g1).lt.another%group_spec(g2)) then
                  cmp=GFC_CMP_LT; exit
                 elseif(this%group_spec(g1).gt.another%group_spec(g2)) then
                  cmp=GFC_CMP_GT; exit
                 else
                  if(gmap(g1).gt.0) then
                   if(gmap(g1).lt.g2) then
                    cmp=GFC_CMP_LT; exit
                   elseif(gmap(g1).gt.g2) then
                    cmp=GFC_CMP_GT; exit
                   endif
                  else
                   gmap(g1)=g2
                  endif
                 endif
                else
                 if(g1.lt.g2) then
                  cmp=GFC_CMP_LT; exit
                 elseif(g1.gt.g2) then
                  cmp=GFC_CMP_GT; exit
                 endif
                endif
               enddo
              else
               cmp=GFC_CMP_ERR
              endif
             endif
            endif
           endif
          endif
         else
          cmp=GFC_CMP_ERR
         endif
         return
        end function TensShapeCompare
!------------------------------------------------------------
        subroutine TensShapePrintIt(this,ierr,dev_id,nspaces)
!Prints the tensor shape.
         implicit none
         class(tens_shape_t), intent(in):: this        !in: tensor shape
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD), intent(in), optional:: dev_id  !in: output device handle (6:screen)
         integer(INTD), intent(in), optional:: nspaces !in: number of preceding spaces (left alignment)
         character(2), parameter:: rel_sign(0:TEREC_NUM_IND_RESTR-1)=(/'  ','< ','> ','<=','>='/)
         integer(INTD):: errc,i,j,dev,ng,gr,gdims(1:MAX_TENSOR_RANK)

         errc=TEREC_SUCCESS
         if(present(dev_id)) then; dev=dev_id; else; dev=6; endif
         if(present(nspaces)) then; do i=1,nspaces; write(dev,'(" ")',ADVANCE='NO'); enddo; endif
         if(this%is_set()) then
          write(dev,'("SHAPE")',ADVANCE='NO')
          write(dev,'("(")',ADVANCE='NO')
          if(this%num_groups().gt.0) then
           do i=1,this%num_dims
            write(dev,'(i9,":",i2)',ADVANCE='NO') this%dim_extent(i),this%dim_group(i)
            if(i.lt.this%num_dims) write(dev,'(",")',ADVANCE='NO')
           enddo
           write(dev,'("):")',ADVANCE='NO')
           do i=1,this%num_grps
            write(dev,'(1x,i2,"[")',ADVANCE='NO') i
            call this%get_group(i,gdims,ng,errc,gr); if(errc.ne.TEREC_SUCCESS) exit
            do j=1,ng
             write(dev,'(i2)',ADVANCE='NO') gdims(j)
             if(j.lt.ng) write(dev,'(A2)',ADVANCE='NO') rel_sign(gr)
            enddo
            write(dev,'("]")',ADVANCE='NO')
           enddo
           write(dev,'()')
          else
           do i=1,this%num_dims
            write(dev,'(i9)',ADVANCE='NO') this%dim_extent(i)
            if(i.lt.this%num_dims) write(dev,'(",")',ADVANCE='NO')
           enddo
           write(dev,'(")")')
          endif
         else
          write(dev,'("EMPTY TENSOR SHAPE")')
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensShapePrintIt
!---------------------------------------
        subroutine tens_shape_dtor(this)
!DTOR for tens_shape_t.
         implicit none
         type(tens_shape_t):: this

         !write(*,'("#DEBUG: Entered tens_shape_dtor")') !debug
         if(allocated(this%group_spec)) deallocate(this%group_spec)
         if(allocated(this%dim_group)) deallocate(this%dim_group)
         if(allocated(this%dim_extent)) deallocate(this%dim_extent)
         this%num_dims=-1; this%num_grps=0
         !write(*,'("#DEBUG: Exited tens_shape_dtor")') !debug
         return
        end subroutine tens_shape_dtor
![tens_header_t]==================================================================================
        subroutine tens_header_ctor(this,ierr,tens_name,subspaces,dim_extent,dim_group,group_spec)
!CTOR for tens_header_t. Each subsequent optional argument implies the existence of all preceding
!optional arguments, except <ierr> and <tens_name>. If no optional arguments are present, except
!maybe <tens_name> and/or <ierr>, a scalar header will be constructed. <dim_group> and <group_spec>
!must either be both present or both absent. More specifically:
! # Constructing a scalar tensor header: Do not pass any optional arguments except <tens_name> and/or <ierr>;
! # Constructing a true tensor header without shape: Pass only <subspaces> and, optionally, <tens_name> and/or <ierr>;
! # Constructing a true tensor header with a shape: Pass <subspaces> and <dim_extent> with all other arguments optional.
         implicit none
         type(tens_header_t), intent(out):: this              !out: tensor header
         integer(INTD), intent(out), optional:: ierr          !out: error code
         character(*), intent(in), optional:: tens_name       !in: alphanumeric_ tensor name
         integer(INTL), intent(in), optional:: subspaces(1:)  !in: subspace multi-index (specification): Length = tensor rank
         integer(INTL), intent(in), optional:: dim_extent(1:) !in: dimension extents: Length = tensor rank
         integer(INTD), intent(in), optional:: dim_group(1:)  !in: dimension restriction groups: Length = tensor rank
         integer(INTD), intent(in), optional:: group_spec(1:) !in: dimension restriction group specification
         integer(INTD):: errc
         logical:: pr_nam,pr_sub,pr_dim,pr_grp,pr_grs

         errc=TEREC_SUCCESS
         pr_nam=present(tens_name)
         pr_sub=present(subspaces)
         pr_dim=present(dim_extent)
         pr_grp=present(dim_group)
         pr_grs=present(group_spec)
         if(pr_dim.and.(.not.pr_sub)) errc=TEREC_INVALID_ARGS
         if((pr_grp.or.pr_grs).and.(.not.pr_dim)) errc=TEREC_INVALID_ARGS
         if((pr_grp.and.(.not.pr_grs)).or.(pr_grs.and.(.not.pr_grp))) errc=TEREC_INVALID_ARGS
         if(pr_sub.and.pr_dim) then
          if(size(subspaces).ne.size(dim_extent)) errc=TEREC_INVALID_ARGS
         endif
         if(errc.eq.TEREC_SUCCESS) then
 !tensor signature:
          if(pr_sub) then !explicit tensor
           if(pr_nam) then
            call tens_signature_ctor(this%signature,errc,subspaces,tens_name)
           else
            call tens_signature_ctor(this%signature,errc,subspaces)
           endif
          else !scalar
           if(pr_nam) then
            call tens_signature_ctor(this%signature,errc,tens_name=tens_name)
           else
            call tens_signature_ctor(this%signature,errc)
           endif
           if(errc.eq.TEREC_SUCCESS) call tens_shape_ctor(this%shape,errc)
          endif
 !tensor shape (optional):
          if(errc.eq.TEREC_SUCCESS) then
           if(pr_dim) then
            if(pr_grp) then
             call tens_shape_ctor(this%shape,errc,dim_extent,dim_group,group_spec)
            else
             call tens_shape_ctor(this%shape,errc,dim_extent)
            endif
           endif
          endif
         endif
         if(errc.ne.TEREC_SUCCESS) call tens_header_dtor(this)
         if(present(ierr)) ierr=errc
         return
        end subroutine tens_header_ctor
!-----------------------------------------------------------------------------------------
        subroutine TensHeaderAddShape(this,ierr,dim_extent,dim_group,group_spec,overwrite)
!Sets the tensor header shape in case it has not been set initially.
         implicit none
         class(tens_header_t), intent(inout):: this           !inout: tensor header
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTL), intent(in), optional:: dim_extent(1:) !in: dimension extents: Length = tensor rank
         integer(INTD), intent(in), optional:: dim_group(1:)  !in: dimension restriction groups: Length = tensor rank
         integer(INTD), intent(in), optional:: group_spec(1:) !in: dimension restriction group specification
         logical, intent(in), optional:: overwrite            !in: if TRUE, the exsiting shape will be overwritten, defaults to FALSE
         integer(INTD):: errc,n
         logical:: overwr

         errc=TEREC_SUCCESS
         if(this%is_set(num_dims=n)) then
          if(present(overwrite)) then; overwr=overwrite; else; overwr=.FALSE.; endif
          if((.not.this%shape%is_set()).or.overwr) then
           if(present(dim_extent)) then
            if(size(dim_extent).eq.n) then
             if(n.gt.0) then
              if(present(dim_group)) then
               if(present(group_spec)) then
                call tens_shape_ctor(this%shape,errc,dim_extent,dim_group,group_spec)
               else
                errc=TEREC_INVALID_ARGS
               endif
              else
               if(present(group_spec)) then
                errc=TEREC_INVALID_ARGS
               else
                call tens_shape_ctor(this%shape,errc,dim_extent)
               endif
              endif
             else !scalar shape
              call tens_shape_ctor(this%shape,errc)
             endif
            else
             errc=TEREC_INVALID_ARGS
            endif
           else !scalar shape
            if(n.ne.0.or.present(dim_group).or.present(group_spec)) then
             errc=TEREC_INVALID_ARGS
            else
             call tens_shape_ctor(this%shape,errc)
            endif
           endif
          else
           errc=TEREC_INVALID_REQUEST
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensHeaderAddShape
!---------------------------------------------------------------------------------
        function TensHeaderIsSet(this,ierr,num_dims,num_groups,shaped) result(res)
!Returns TRUE if the tensor header is set (with or without shape).
         implicit none
         logical:: res                                     !out: result
         class(tens_header_t), intent(in):: this           !in: tensor header
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD), intent(out), optional:: num_dims   !out: number of dimensions
         integer(INTD), intent(out), optional:: num_groups !out: number of restricted dimension groups
         logical, intent(out), optional:: shaped           !out: TRUE if the tensor shape is set
         integer(INTD):: errc,nd,ng
         logical:: shpd

         res=this%signature%is_set(errc,num_dims=nd)
         if(present(num_dims)) num_dims=nd
         if(res.and.errc.eq.TEREC_SUCCESS) then
          shpd=this%shape%is_set(errc,num_dims=nd,num_groups=ng)
         else
          shpd=.FALSE.; ng=0
         endif
         if(present(num_groups)) num_groups=ng
         if(present(shaped)) shaped=shpd
         if(present(ierr)) ierr=errc
         return
        end function TensHeaderIsSet
!-----------------------------------------------------------------
        subroutine TensHeaderGetName(this,tens_name,name_len,ierr)
!Returns the alphanumeric_ name of the tensor.
         implicit none
         class(tens_header_t), intent(in):: this     !in: tensor header
         character(*), intent(inout):: tens_name     !out: tensor name
         integer(INTD), intent(out):: name_len       !out: length of the tensor name
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         call this%signature%get_name(tens_name,name_len,errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensHeaderGetName
!--------------------------------------------------------
        function TensHeaderGetRank(this,ierr) result(res)
!Returns the rank of the tensor (number of dimensions).
         implicit none
         integer(INTD):: res                         !out: result
         class(tens_header_t), intent(in):: this     !in: tensor header
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         res=this%signature%get_rank(errc)
         if(present(ierr)) ierr=errc
         return
        end function TensHeaderGetRank
!-----------------------------------------------------------------
        subroutine TensHeaderGetSpec(this,subspaces,num_dims,ierr)
!Returns the defining subspaces of the tensor (subspace multi-index).
         implicit none
         class(tens_header_t), intent(in):: this      !in: tensor header
         integer(INTL), intent(inout):: subspaces(1:) !out: defining subspaces (their IDs)
         integer(INTD), intent(out):: num_dims        !out: number of tensor dimensions
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         call this%signature%get_spec(subspaces,num_dims,errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensHeaderGetSpec
!------------------------------------------------------------
        subroutine TensHeaderGetDims(this,dims,num_dims,ierr)
!Returns tensor dimension extents together with the tensor rank.
         implicit none
         class(tens_header_t), intent(in):: this     !in: tensor header
         integer(INTL), intent(inout):: dims(1:)     !out: tensor dimension extents
         integer(INTD), intent(out):: num_dims       !out: number of tensor dimensions
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         call this%shape%get_dims(dims,num_dims,errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensHeaderGetDims
!-----------------------------------------------------------------
        function TensHeaderNumGroups(this,ierr) result(num_groups)
!Returns the total number of non-trivial index groups defined in the tensor header shape.
         implicit none
         integer(INTD):: num_groups                  !out: number of non-trivial dimension groups
         class(tens_header_t), intent(in):: this     !in: tensor header
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         num_groups=this%shape%num_groups(errc)
         if(present(ierr)) ierr=errc
         return
        end function TensHeaderNumGroups
!-------------------------------------------------------------------------------------
        subroutine TensHeaderGetGroup(this,group,group_dims,num_dims,ierr,group_restr)
!Returns the index restriction group (specific dimensions belonging to the specified group).
         implicit none
         class(tens_header_t), intent(in):: this            !in: tensor header
         integer(INTD), intent(in):: group                  !in: requested dimension group
         integer(INTD), intent(inout):: group_dims(1:)      !out: dimensions which belong to the requested dimension group
         integer(INTD), intent(out):: num_dims              !out: number of dimensions in the group
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD), intent(out), optional:: group_restr !out: group restriction kind
         integer(INTD):: errc,gr

         call this%shape%get_group(group,group_dims,num_dims,errc,gr)
         if(present(group_restr)) group_restr=gr
         if(present(ierr)) ierr=errc
         return
        end subroutine TensHeaderGetGroup
!-----------------------------------------------------------------------------
        function TensHeaderSameGroup(this,dims,ierr,group_restr) result(group)
!Returns the group number if the specified tensor dimensions belong to the same group, -1 otherwise.
         implicit none
         integer(INTD):: group                              !out: group number (>=0)
         class(tens_header_t), intent(in):: this            !in: tensor header
         integer(INTD), intent(in):: dims(1:)               !in: tensor dimensions to check
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD), intent(out), optional:: group_restr !out: group restriction kind
         integer(INTD):: errc,gr

         group=this%shape%same_group(dims,errc,gr)
         if(present(group_restr)) group_restr=gr
         if(present(ierr)) ierr=errc
         return
        end function TensHeaderSameGroup
!---------------------------------------------------------------------
        function TensHeaderGetSignature(this,ierr) result(signature_p)
!Returns a pointer to the tensor header signature.
         implicit none
         type(tens_signature_t), pointer:: signature_p   !out: pointer to the tensor signature
         class(tens_header_t), intent(in), target:: this !in: tensor header
         integer(INTD), intent(out), optional:: ierr     !out: error code

         signature_p=>this%signature
         return
        end function TensHeaderGetSignature
!-------------------------------------------------------------
        function TensHeaderGetShape(this,ierr) result(shape_p)
!Returns a pointer to the tensor header shape.
         implicit none
         type(tens_shape_t), pointer:: shape_p           !out: pointer to the tensor shape
         class(tens_header_t), intent(in), target:: this !in: tensor header
         integer(INTD), intent(out), optional:: ierr     !out: error code

         shape_p=>this%shape
         return
        end function TensHeaderGetShape
!--------------------------------------------------------------------------
        function TensHeaderCompare(this,another,compare_groups) result(cmp)
!Compares the given tensor header with another tensor header.
         implicit none
         integer(INTD):: cmp                            !out: comparison result: {GFC_CMP_EQ,GFC_CMP_LT,GFC_CMP_GT,GFC_CMP_PARENT,GFC_CMP_CHILD,GFC_CMP_NA,GFC_CMP_ERR}
         class(tens_header_t), intent(in):: this        !in: tensor header 1
         class(tens_header_t), intent(in):: another     !in: tensor header 2
         logical, intent(in), optional:: compare_groups !in: if FALSE, the shape dimension groups will not be taken into account (defaults to TRUE)

         cmp=this%signature%compare(another%signature)
         if(cmp.eq.GFC_CMP_EQ) then
          if(present(compare_groups)) then
           cmp=this%shape%compare(another%shape,compare_groups=compare_groups)
          else
           cmp=this%shape%compare(another%shape)
          endif
         endif
         return
        end function TensHeaderCompare
!-------------------------------------------------------------
        subroutine TensHeaderPrintIt(this,ierr,dev_id,nspaces)
!Prints the tensor header.
         implicit none
         class(tens_header_t), intent(in):: this       !in: tensor header
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD), intent(in), optional:: dev_id  !in: output device id (6:screen)
         integer(INTD), intent(in), optional:: nspaces !out: left alignment
         integer(INTD):: errc,dev,i

         errc=TEREC_SUCCESS
         if(present(dev_id)) then; dev=dev_id; else; dev=6; endif
         if(present(nspaces)) then
          do i=1,nspaces; write(dev,'(" ")',ADVANCE='NO'); enddo
          write(dev,'("TENSOR HEADER{")')
          call this%signature%print_it(errc,dev,nspaces+1)
          call this%shape%print_it(errc,dev,nspaces+1)
          do i=1,nspaces; write(dev,'(" ")',ADVANCE='NO'); enddo
          write(dev,'("}")')
         else
          write(dev,'("TENSOR HEADER{")')
          call this%signature%print_it(errc,dev,1)
          call this%shape%print_it(errc,dev,1)
          write(dev,'("}")')
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensHeaderPrintIt
!----------------------------------------
        subroutine tens_header_dtor(this)
!DTOR for tens_header_t.
         implicit none
         type(tens_header_t):: this

         return
        end subroutine tens_header_dtor
![tens_simple_part_t]===================================================
        subroutine tens_simple_part_ctor(this,header,layout,offset,ierr)
!CTOR for tens_simple_part_t.
         implicit none
         type(tens_simple_part_t), intent(out):: this !out: tensor simple part
         class(tens_header_t), intent(in):: header    !in: tensor header
         integer(INTD), intent(in):: layout           !in: simple storage layout: {TEREC_LAY_FDIMS,TEREC_LAY_CDIMS}
         integer(INTL), intent(in):: offset           !in: offset in the parental tensor block
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc
         logical:: shpd

         errc=TEREC_SUCCESS
         if(header%is_set(errc,shaped=shpd)) then
          if(errc.eq.TEREC_SUCCESS.and.shpd) then
           if(layout.eq.TEREC_LAY_FDIMS.or.layout.eq.TEREC_LAY_CDIMS) then
            this%header=header
            this%offset=offset
            this%layout=layout
           else
            errc=TEREC_INVALID_ARGS
           endif
          else
           errc=TEREC_INVALID_ARGS
          endif
         else
          errc=TEREC_INVALID_ARGS
         endif
         if(errc.ne.TEREC_SUCCESS) call tens_simple_part_dtor(this)
         if(present(ierr)) ierr=errc
         return
        end subroutine tens_simple_part_ctor
!----------------------------------------------------------
        function TensSimplePartIsSet(this,ierr) result(res)
!Returns TRUE of the tensor simple part is set.
         implicit none
         logical:: res                                !out: result
         class(tens_simple_part_t), intent(in):: this !in: tensor simple part
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         res=(this%layout.ne.TEREC_LAY_NONE)
         if(present(ierr)) ierr=errc
         return
        end function TensSimplePartIsSet
!-----------------------------------------------------------------
        function TensSimplePartGetOffset(this,ierr) result(offset)
!Returns the offset of the tensor simple part.
         implicit none
         integer(INTL):: offset                       !out: offset
         class(tens_simple_part_t), intent(in):: this !in: tensor simple part
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         if(this%is_set()) then
          offset=this%offset
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensSimplePartGetOffset
!-----------------------------------------------------------------
        function TensSimplePartGetLayout(this,ierr) result(layout)
!Returns the layout of the tensor simple part.
         implicit none
         integer(INTL):: layout                       !out: layout
         class(tens_simple_part_t), intent(in):: this !in: tensor simple part
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         if(this%is_set()) then
          layout=this%layout
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensSimplePartGetLayout
!-----------------------------------------------------------------
        function TensSimplePartGetHeader(this,ierr) result(header)
!Returns a pointer to the header of the tensor simple part.
         implicit none
         type(tens_header_t), pointer:: header                !out: pointer to the header
         class(tens_simple_part_t), intent(in), target:: this !in: tensor simple part
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         if(this%is_set()) then
          header=>this%header
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensSimplePartGetHeader
!---------------------------------------------
        subroutine tens_simple_part_dtor(this)
!DTOR for tens_simple_part_t.
         implicit none
         type(tens_simple_part_t):: this

         this%offset=-1_INTL
         this%layout=TEREC_LAY_NONE
         return
        end subroutine tens_simple_part_dtor
![tens_layout_t]===============================
        function TensLayoutIsSet(this,ierr) result(res)
!Returns TRUE if the tensor layout is set.
         implicit none
         logical:: res                               !out: result
         class(tens_layout_t), intent(in):: this     !in: tensor layout
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         res=(this%layout.ne.TEREC_LAY_NONE)
         if(present(ierr)) ierr=errc
         return
        end function TensLayoutIsSet
!-------------------------------------------------------------
        function TensLayoutGetLayout(this,ierr) result(layout)
!Returns the tensor layout.
         implicit none
         integer(INTD):: layout                      !out: tensor layout
         class(tens_layout_t), intent(in):: this     !in: tensor layout
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         layout=this%layout
         if(present(ierr)) ierr=errc
         return
        end function TensLayoutGetLayout

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

         write(*,'("Testing class tens_signature_t ... ")',ADVANCE='NO')
         call test_tens_signature(ierr)
         if(ierr.eq.0) then; write(*,'("PASSED")'); else; write(*,'("FAILED: Error ",i11)') ierr; return; endif
         write(*,'("Testing class tens_shape_t ... ")',ADVANCE='NO')
         call test_tens_shape(ierr)
         if(ierr.eq.0) then; write(*,'("PASSED")'); else; write(*,'("FAILED: Error ",i11)') ierr; return; endif
         write(*,'("Testing class tens_header_t ... ")',ADVANCE='NO')
         call test_tens_header(ierr)
         if(ierr.eq.0) then; write(*,'("PASSED")'); else; write(*,'("FAILED: Error ",i11)') ierr; return; endif
         write(*,'("Testing class tens_simple_part_t ... ")',ADVANCE='NO')
         call test_tens_simple_part(ierr)
         if(ierr.eq.0) then; write(*,'("PASSED")'); else; write(*,'("FAILED: Error ",i11)') ierr; return; endif
         return
        end subroutine test_tensor_recursive
!-------------------------------------------
        subroutine test_tens_signature(ierr)
         implicit none
         integer(INTD), intent(out):: ierr
         integer(INTD):: n
         integer(INTL):: mlndx(1:MAX_TENSOR_RANK)
         character(32):: tname
         type(tens_signature_t):: tsigna

         ierr=0
         call tens_signature_ctor(tsigna,ierr,(/3_INTL,4_INTL,2_INTL/),'Tensor')
         if(ierr.eq.TEREC_SUCCESS.and.tsigna%is_set()) then
          n=tsigna%get_rank(ierr)
          if(ierr.eq.TEREC_SUCCESS.and.n.eq.3) then
           call tsigna%get_spec(mlndx,n,ierr)
           if(ierr.eq.TEREC_SUCCESS.and.n.eq.3.and.&
             &mlndx(1).eq.3.and.mlndx(2).eq.4.and.mlndx(3).eq.2) then
            call tsigna%get_name(tname,n,ierr)
            if(.not.(n.eq.len('Tensor').and.tname(1:n).eq.'Tensor')) ierr=4
           else
            ierr=3
           endif
          else
           ierr=2
          endif
         else
          ierr=1
         endif
         !call tens_signature_dtor(tsigna)
         return
        end subroutine test_tens_signature
!---------------------------------------
        subroutine test_tens_shape(ierr)
         implicit none
         integer(INTD), intent(out):: ierr
         integer(INTD):: i,m,n
         integer(INTL):: dims(1:MAX_TENSOR_RANK)
         integer(INTD):: grps(1:MAX_TENSOR_RANK)
         integer(INTD):: grp_spec(1:MAX_TENSOR_RANK)
         type(tens_shape_t):: tshape

         ierr=0
         n=6; dims(1:n)=(/128_INTL,64_INTL,256_INTL,64_INTL,128_INTL,64_INTL/)
         m=2; grps(1:n)=(/1,2,0,2,1,2/); grp_spec(1:m)=(/TEREC_IND_RESTR_LT,TEREC_IND_RESTR_GE/)
         call tens_shape_ctor(tshape,ierr,dims(1:n),grps(1:n),grp_spec(1:m))
         if(ierr.eq.TEREC_SUCCESS) then
          !write(*,'()'); call tshape%print_it() !debug
          grps(1:3)=(/2,4,6/)
          if(tshape%same_group(grps(1:3),ierr,i).eq.2) then
           if(ierr.eq.TEREC_SUCCESS.and.i.eq.TEREC_IND_RESTR_GE) then
            if(tshape%num_groups(ierr).eq.m) then
             if(ierr.ne.TEREC_SUCCESS) ierr=5
            else
             ierr=4
            endif
           else
            ierr=3
           endif
          else
           ierr=2
          endif
         else
          ierr=1
         endif
         !call tens_shape_dtor(tshape)
         return
        end subroutine test_tens_shape
!----------------------------------------
        subroutine test_tens_header(ierr)
         implicit none
         integer(INTD), intent(out):: ierr
         integer(INTD):: i,l,m,n
         integer(INTL):: dims(1:MAX_TENSOR_RANK)
         integer(INTD):: grps(1:MAX_TENSOR_RANK)
         integer(INTD):: grp_spec(1:MAX_TENSOR_RANK)
         type(tens_signature_t):: tsigna
         type(tens_shape_t):: tshape
         type(tens_header_t):: thead
         character(32):: tens_name

         ierr=0
         n=6; dims(1:n)=(/128_INTL,64_INTL,256_INTL,64_INTL,128_INTL,64_INTL/)
         m=2; grps(1:n)=(/1,2,0,2,1,2/); grp_spec(1:m)=(/TEREC_IND_RESTR_LT,TEREC_IND_RESTR_GE/)
         call tens_header_ctor(thead,ierr,'Tensor',(/1_INTL,2_INTL,3_INTL,2_INTL,1_INTL,2_INTL/))
         if(ierr.eq.TEREC_SUCCESS) then
          call thead%add_shape(ierr,dims(1:n),grps(1:n),grp_spec(1:m))
          if(ierr.eq.TEREC_SUCCESS) then
           !call thead%print_it(ierr) !debug
           if(ierr.eq.TEREC_SUCCESS) then
            call thead%get_name(tens_name,l,ierr)
            if(ierr.eq.TEREC_SUCCESS.and.tens_name(1:l).eq.'Tensor') then
             if(thead%get_rank(ierr).eq.6) then
              call thead%get_dims(dims,n,ierr)
              if(ierr.eq.TEREC_SUCCESS.and.n.eq.6.and.dims(1).eq.128.and.dims(2).eq.64.and.&
                                                     &dims(3).eq.256.and.dims(4).eq.64.and.&
                                                     &dims(5).eq.128.and.dims(6).eq.64) then
              else
               ierr=6
              endif
             else
              ierr=5
             endif
            else
             ierr=4
            endif
           else
            ierr=3
           endif
          else
           ierr=2
          endif
         else
          ierr=1
         endif
         !call tens_header_dtor(thead)
         return
        end subroutine test_tens_header
!---------------------------------------------
        subroutine test_tens_simple_part(ierr)
         implicit none
         integer(INTD), intent(out):: ierr
         integer(INTD):: i,l,m,n
         integer(INTL):: dims(1:MAX_TENSOR_RANK)
         integer(INTD):: grps(1:MAX_TENSOR_RANK)
         integer(INTD):: grp_spec(1:MAX_TENSOR_RANK)
         type(tens_signature_t):: tsigna
         type(tens_shape_t):: tshape
         type(tens_header_t):: thead
         type(tens_header_t), pointer:: thp
         character(32):: tens_name
         type(tens_simple_part_t):: tpart

         ierr=0
         n=6; dims(1:n)=(/128_INTL,64_INTL,256_INTL,64_INTL,128_INTL,64_INTL/)
         m=2; grps(1:n)=(/1,2,0,2,1,2/); grp_spec(1:m)=(/TEREC_IND_RESTR_LT,TEREC_IND_RESTR_GE/)
         call tens_header_ctor(thead,ierr,'Tensor',(/1_INTL,2_INTL,3_INTL,2_INTL,1_INTL,2_INTL/))
         if(ierr.eq.TEREC_SUCCESS) then
          call thead%add_shape(ierr,dims(1:n),grps(1:n),grp_spec(1:m))
          if(ierr.eq.TEREC_SUCCESS) then
           !call thead%print_it(ierr) !debug
           if(ierr.eq.TEREC_SUCCESS) then
            call thead%get_name(tens_name,l,ierr)
            if(ierr.eq.TEREC_SUCCESS.and.tens_name(1:l).eq.'Tensor') then
             if(thead%get_rank(ierr).eq.6) then
              dims(:)=0_INTL
              call thead%get_dims(dims,n,ierr)
              if(ierr.eq.TEREC_SUCCESS.and.n.eq.6.and.dims(1).eq.128.and.dims(2).eq.64.and.&
                                                     &dims(3).eq.256.and.dims(4).eq.64.and.&
                                                     &dims(5).eq.128.and.dims(6).eq.64) then
               call tens_simple_part_ctor(tpart,thead,TEREC_LAY_FDIMS,4096_INTL,ierr)
               if(ierr.eq.TEREC_SUCCESS) then
                if(tpart%get_offset(ierr).eq.4096.and.tpart%get_layout(i).eq.TEREC_LAY_FDIMS) then
                 if(ierr.eq.TEREC_SUCCESS.and.i.eq.TEREC_SUCCESS) then
                  thp=>tpart%get_header(ierr)
                  if(ierr.eq.TEREC_SUCCESS.and.associated(thp)) then
                   tens_name=' '
                   call thp%get_name(tens_name,l,ierr)
                   if(ierr.eq.TEREC_SUCCESS.and.tens_name(1:l).eq.'Tensor') then
                    if(thp%get_rank(ierr).eq.6) then
                     dims(:)=0_INTL
                     call thp%get_dims(dims,n,ierr)
                     if(ierr.eq.TEREC_SUCCESS.and.n.eq.6.and.dims(1).eq.128.and.dims(2).eq.64.and.&
                                                     &dims(3).eq.256.and.dims(4).eq.64.and.&
                                                     &dims(5).eq.128.and.dims(6).eq.64) then
                     else
                      ierr=13
                     endif
                    else
                     ierr=12
                    endif
                   else
                    ierr=11
                   endif
                  else
                   ierr=10
                  endif
                 else
                  ierr=9
                 endif
                else
                 ierr=8
                endif
               else
                ierr=7
               endif
              else
               ierr=6
              endif
             else
              ierr=5
             endif
            else
             ierr=4
            endif
           else
            ierr=3
           endif
          else
           ierr=2
          endif
         else
          ierr=1
         endif
         !call tens_simple_part_dtor(tpart)
         !call tens_header_dtor(thead)
         return
        end subroutine test_tens_simple_part

       end module tensor_recursive_test
