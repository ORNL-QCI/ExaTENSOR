!ExaTENSOR: Recursive tensors
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/04/10

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
        use gfc_vec_tree
        use subspaces
        use pack_prim
        use distributed
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !output device
        integer(INTD), private:: DEBUG=0    !debugging level
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
 !Error codes:
        integer(INTD), parameter, public:: TEREC_SUCCESS=SUCCESS
        integer(INTD), parameter, public:: TEREC_ERROR=GENERIC_ERROR
        integer(INTD), parameter, public:: TEREC_INVALID_ARGS=-1
        integer(INTD), parameter, public:: TEREC_INVALID_REQUEST=-2
        integer(INTD), parameter, public:: TEREC_MEM_ALLOC_FAILED=-3
        integer(INTD), parameter, public:: TEREC_MEM_FREE_FAILED=-4
        integer(INTD), parameter, public:: TEREC_UNABLE_COMPLETE=-5
        integer(INTD), parameter, public:: TEREC_OBJ_CORRUPTED=-6
 !Symbolic tensor name:
        integer(INTD), parameter, public:: TEREC_MAX_TENS_NAME_LEN=64 !max length of the tensor name (alphanumeric_ string)
 !Storage layout for tensor blocks for locally stored tensors (must span a contiguous integer range!):
        integer(INTD), parameter, public:: TEREC_LAY_NONE=0    !none
        integer(INTD), parameter, public:: TEREC_LAY_RECUR=1   !storage layout is inferred from that of individual constituent tensors (recursive)
        integer(INTD), parameter, public:: TEREC_LAY_FDIMS=2   !straightforward dimension led storage layout (Fortran style)
        integer(INTD), parameter, public:: TEREC_LAY_CDIMS=3   !straightforward dimension led storage layout (C style)
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
         integer(INTD), private:: num_dims=-1                   !number of tensor dimensions (aka tensor order in math or tensor rank in physics)
         integer(INTL), allocatable, private:: space_idx(:)     !subspace id for each tensor dimension
         !`Not packed:
         class(h_space_t), pointer, private:: h_space_p=>NULL() !pointer to the underlying hierarchical vector space specification (external target!)
         contains
          procedure, private:: TensSignatureCtor                   !ctor
          procedure, private:: TensSignatureCtorUnpack             !ctor by unpacking
          generic, public:: tens_signature_ctor=>TensSignatureCtor,TensSignatureCtorUnpack
          procedure, public:: pack=>TensSignaturePack              !packs the object into a packet
          procedure, public:: is_set=>TensSignatureIsSet           !returns .TRUE. if the tensor signature is set
          procedure, public:: get_name=>TensSignatureGetName       !returns the alphanumeric_ tensor name
          procedure, public:: get_rank=>TensSignatureGetRank       !returns the rank of the tensor (number of dimensions)
          procedure, public:: get_spec=>TensSignatureGetSpec       !returns the tensor subspace multi-index (specification)
          procedure, public:: relate=>TensSignatureRelate          !relates the tensor signature to another tensor signature: {CMP_EQ,CMP_CN,CMP_IN,CMP_OV,CMP_NC}
          procedure, public:: compare=>TensSignatureCompare        !compares the tensor signature with another tensor signature: {CMP_EQ,CMP_LT,CMP_GT,CMP_ER}
          procedure, public:: print_it=>TensSignaturePrintIt       !prints the tensor signature
          final:: tens_signature_dtor                              !dtor
        end type tens_signature_t
 !Tensor shape:
        type, public:: tens_shape_t
         integer(INTL), allocatable, private:: dim_extent(:) !tensor dimension extents (resolution)
         integer(INTD), allocatable, private:: dim_group(:)  !tensor dimension groups (group 0 is default with no restrictions)
         integer(INTD), allocatable, private:: group_spec(:) !specification of the restriction kind for each index restriction group (not every defined group needs to be present in dim_group(:))
         integer(INTD), private:: num_dims=-1                !number of tensor dimensions (aka tensor order or tensor rank)
         integer(INTD), private:: num_grps=0                 !number of defined (non-trivial) index restriction groups
         contains
          procedure, private:: TensShapeCtor                 !ctor
          procedure, private:: TensShapeCtorUnpack           !ctor by unpacking
          generic, public:: tens_shape_ctor=>TensShapeCtor,TensShapeCtorUnpack
          procedure, public:: pack=>TensShapePack            !packs the object into a packet
          procedure, public:: set_dims=>TensShapeSetDims     !sets dimension extents (if they have not been set previously)
          procedure, public:: set_groups=>TensShapeSetGroups !creates new index restriction groups
          procedure, public:: is_set=>TensShapeIsSet         !returns .TRUE. if the tensor shape is set
          procedure, public:: get_dims=>TensShapeGetDims     !returns tensor dimension extents
          procedure, public:: get_rank=>TensShapeGetRank     !returns the rank of the tensor (number of dimensions)
          procedure, public:: get_dim_group=>TensShapeGetDimGroup !returns the restriction group for a specific tensor dimension (0: no restrictions)
          procedure, public:: get_group=>TensShapeGetGroup   !returns a restricted index group (specific dimensions belonging to the specified group)
          procedure, public:: same_group=>TensShapeSameGroup !checks whether specific tensor dimensions belong to the same group
          procedure, public:: num_groups=>TensShapeNumGroups !returns the total number of non-trivial index groups defined in the tensor shape
          procedure, public:: compare=>TensShapeCompare      !compares the tensor shape with another tensor shape: {CMP_EQ,CMP_LT,CMP_GT,CMP_ER}
          procedure, public:: print_it=>TensShapePrintIt     !prints the tensor shape
          final:: tens_shape_dtor                            !dtor
        end type tens_shape_t
 !Tensor header (signature+shape):
        type, public:: tens_header_t
         type(tens_signature_t), private:: signature         !tensor signature
         type(tens_shape_t), private:: shape                 !tensor shape
         contains
          procedure, private:: TensHeaderCtor                       !ctor
          procedure, private:: TensHeaderCtorUnpack                 !ctor by unpacking
          generic, public:: tens_header_ctor=>TensHeaderCtor,TensHeaderCtorUnpack
          procedure, public:: pack=>TensHeaderPack                  !packs the object into a packet
          procedure, public:: add_shape=>TensHeaderAddShape         !ctor for a deferred tensor shape specification
          procedure, public:: set_dims=>TensHeaderSetDims           !sets dimension extents (if they have not been set previously)
          procedure, public:: set_groups=>TensHeaderSetGroups       !sets index restriction groups if they have not been previously set
          procedure, public:: is_set=>TensHeaderIsSet               !returns .TRUE. if the tensor header is set (with or without shape)
          procedure, public:: get_name=>TensHeaderGetName           !returns the alphanumeric_ tensor name
          procedure, public:: get_rank=>TensHeaderGetRank           !returns the rank of the tensor (number of dimensions)
          procedure, public:: get_spec=>TensHeaderGetSpec           !returns the tensor subspace multi-index (specification)
          procedure, public:: get_dims=>TensHeaderGetDims           !returns tensor dimension extents
          procedure, public:: num_groups=>TensHeaderNumGroups       !returns the total number of non-trivial index groups defined in the tensor shape
          procedure, public:: get_dim_group=>TensHeaderGetDimGroup  !returns the restriction group for a specific tensor dimension (0: no restrictions)
          procedure, public:: get_group=>TensHeaderGetGroup         !returns a restricted index group (specific dimensions belonging to the specified group)
          procedure, public:: same_group=>TensHeaderSameGroup       !checks whether specific tensor dimensions belong to the same group
          procedure, public:: get_signature=>TensHeaderGetSignature !returns the pointer to the tensor signature
          procedure, public:: get_shape=>TensHeaderGetShape         !returns the pointer the the tensor shape
          procedure, public:: compare=>TensHeaderCompare            !compares the tensor header with another tensor header: {CMP_EQ,CMP_LT,CMP_GT,CMP_ER}
          procedure, public:: print_it=>TensHeaderPrintIt           !prints the tensor header
#ifdef NO_GNU
          final:: tens_header_dtor                                  !dtor `GCC/5.4.0 bug
#endif
        end type tens_header_t
 !Simple (dense) tensor block (part):
        type, public:: tens_simple_part_t
         type(tens_header_t), private:: header               !header of the constituent simple tensor block
         integer(INTL), private:: offset=-1_INTL             !offset of the constituent simple tensor block in the parental composite tensor block (locally stored)
         integer(INTD), private:: layout=TEREC_LAY_NONE      !simple storage layout: {TEREC_LAY_FDIMS,TEREC_LAY_CDIMS} only
         contains
          procedure, private:: TensSimplePartCtor                     !ctor
          procedure, private:: TensSimplePartCtorUnpack               !ctor by unpacking
          generic, public:: tens_simple_part_ctor=>TensSimplePartCtor,TensSimplePartCtorUnpack
          procedure, public:: pack=>TensSimplePartPack                !packs the object into a packet
          procedure, public:: is_set=>TensSimplePartIsSet             !return TRUE if the simple part is set (signature, shape, layout, offset)
          procedure, public:: get_offset=>TensSimplePartGetOffset     !returns the offset of the simple part in the parental tensor block
          procedure, public:: get_layout=>TensSimplePartGetLayout     !returns the simple layout of the simple tensor part: {TEREC_LAY_FDIMS,TEREC_LAY_CDIMS} only
          procedure, public:: get_header=>TensSimplePartGetHeader     !returns a pointer to the header of the simple tensor part
          final:: tens_simple_part_dtor                               !dtor
        end type tens_simple_part_t
 !Storage layout for locally stored blocks (abstract base):
        type, abstract, public:: tens_layout_t
         integer(INTD), private:: layout=TEREC_LAY_NONE          !tensor block storage layout (see above), either simple or composite
         integer(INTD), private:: data_type=NO_TYPE              !data type of tensor elements: {R4,R8,C4,C8}
         class(DataDescr_t), allocatable, private:: data_descr   !DDSS data descriptor for physically stored tensor body (tensor elements)
         contains
          procedure, private:: set_location=>TensLayoutSetLocation                  !sets the phyiscal location of the data via a DDSS data descriptor
          procedure, public:: is_set=>TensLayoutIsSet                               !returns TRUE if the tensor layout is set
          procedure, public:: get_data_type=>TensLayoutGetDataType                  !returns the data type of the stored tensor elements
          procedure, public:: get_layout_kind=>TensLayoutGetLayoutKind              !returns the tensor storage layout kind
          procedure, public:: get_body_ptr=>TensLayoutGetBodyPtr                    !returns a C pointer to the tensor body
          procedure, public:: get_body_size=>TensLayoutGetBodySize                  !returns the size of the stored tensor body in bytes
          procedure(tens_layout_volume_i), deferred, public:: get_volume            !returns the physical volume of the tensor block (number of physically stored elements)
          procedure(tens_layout_map_i), deferred, public:: map                      !maps a specific element of the tensor block (layout specific)
          procedure(tens_layout_extract_i), deferred, public:: extract_simple_parts !creates a list of constituent simple (dense) parts of the tensor block
        end type tens_layout_t
 !Concrete storage layout "Fortran-dimension-led":
        type, extends(tens_layout_t), public:: tens_layout_fdims_t
         class(tens_header_t), pointer, private:: header=>NULL() !pointer to the defining tensor header
         contains
          procedure, private:: TensLayoutFdimsCtor
          generic, public:: tens_layout_fdims_ctor=>TensLayoutFdimsCtor
          procedure, public:: get_volume=>TensLayoutFdimsGetVolume
          procedure, public:: map=>TensLayoutFdimsMap
          procedure, public:: extract_simple_parts=>TensLayoutFdimsExtract
          final:: tens_layout_fdims_dtor
        end type tens_layout_fdims_t
 !Tensor body:
        type, public:: tens_body_t
         integer(INTD), private:: num_subtensors=0            !number of subtensors in the list
         type(list_bi_t), private:: subtensors                !list of constituent tensors in terms of tensor headers
         class(tens_layout_t), allocatable, private:: layout  !tensor block storage layout (if physically stored as a whole)
         contains
          procedure, private:: TensBodyCtorBase                     !basic ctor (layout + data type)
          generic, public:: tens_body_ctor=>TensBodyCtorBase        !ctors
          procedure, public:: is_set=>TensBodyIsSet                 !returns TRUE if the tensor body is set (plus additional info)
          procedure, public:: add_subtensor=>TensBodyAddSubtensor   !registers a constituent subtensor by providing its tensor header
          procedure, public:: set_layout=>TensBodySetLayout         !sets the tensor body storage layout if physically stored as a whole
          procedure, public:: set_location=>TensBodySetLocation     !sets the tensor body data location if physically stored as a whole (via a DDSS data descriptor)
          procedure, public:: get_layout=>TensBodyGetLayout         !returns a pointer to the tensor body storage layout
          procedure, public:: get_num_subtensors=>TensBodyGetNumSubtensors !returns the total number of constituent subtensors
          procedure, public:: get_subtensors=>TensBodyGetSubtensors !returns a pointer to the list of constituent subtensors (each subtensor is represented by a tensor header)
          final:: tens_body_dtor                                    !dtor
        end type tens_body_t
 !Recursive tensor:
        type, public:: tens_rcrsv_t
         type(tens_header_t), private:: header                !tensor header (signature + shape)
         type(tens_body_t), private:: body                    !tensor body (recursive composition, data location and storage layout)
         contains
          procedure, private:: TensRcrsvCtorSigna                    !ctor by tensor signature and optionally tensor shape
          procedure, private:: TensRcrsvCtorHead                     !ctor by tensor header
          generic, public:: tens_rcrsv_ctor=>TensRcrsvCtorSigna,TensRcrsvCtorHead
          procedure, public:: is_set=>TensRcrsvIsSet                 !returns TRUE if the tensor is set (signature defined) plus other info
          procedure, public:: add_subtensor=>TensRcrsvAddSubtensor   !registers a constituent subtensor by providing its tensor header
          procedure, public:: add_subtensors=>TensRcrsvAddSubtensors !registers constituent subtensors by providing a list of their tensor headers
          procedure, public:: set_shape=>TensRcrsvSetShape           !sets the tensor shape (if it has not been set yet)
          procedure, public:: set_layout=>TensRcrsvSetLayout         !sets the tensor body storage layout
          procedure, public:: set_location=>TensRcrsvSetLocation     !sets the physical location of the tensor body data
          procedure, public:: get_header=>TensRcrsvGetHeader         !returns a pointer to the tensor header
          procedure, public:: get_body=>TensRcrsvGetBody             !returns a pointer to the tensor body
          procedure, public:: split=>TensRcrsvSplit                  !splits the tensor into subtensors (a list of subtensors by their headers)
          final:: tens_rcrsv_dtor                                    !dtor
        end type tens_rcrsv_t
!INTERFACES:
 !Abstract:
        abstract interface
  !tens_layout_t: .volume():
         function tens_layout_volume_i(this) result(vol)
          import:: INTL,tens_layout_t
          implicit none
          integer(INTL):: vol                         !out: physical volume (number of tensor elements physically stored)
          class(tens_layout_t), intent(in):: this     !in: tensor block storage layout
         end function tens_layout_volume_i
  !tens_layout_t: .map():
         function tens_layout_map_i(this,mlndx) result(offset)
          import:: INTL,tens_layout_t
          implicit none
          integer(INTL):: offset                      !out: offset of the requested tensor element
          class(tens_layout_t), intent(in):: this     !in: tensor block storage layout
          integer(INTL), intent(in):: mlndx(1:)       !in: input multi-index specifying the individual tensor element
         end function tens_layout_map_i
  !tens_layout_t: .extract_simple_parts():
         subroutine tens_layout_extract_i(this,num_parts,parts,ierr)
          import:: INTD,INTL,list_bi_t,tens_layout_t
          implicit none
          class(tens_layout_t), intent(in):: this     !in: tensor block storage layout
          integer(INTL), intent(out):: num_parts      !out: number of constituent simple (dense) blocks
          type(list_bi_t), intent(inout):: parts      !out: list of the constituent simple (dense) blocks with their headers and offsets
          integer(INTD), intent(out), optional:: ierr !out: error code
         end subroutine tens_layout_extract_i
        end interface
!VISIBILITY:
        public valid_tensor_layout
        public cmp_tens_signatures
        public cmp_tens_headers
        public print_tens_header_f
 !tens_signature_t:
        private TensSignatureCtor
        private TensSignatureCtorUnpack
        private TensSignaturePack
        private TensSignatureIsSet
        private TensSignatureGetName
        private TensSignatureGetRank
        private TensSignatureGetSpec
        private TensSignatureRelate
        private TensSignatureCompare
        private TensSignaturePrintIt
        public tens_signature_dtor
 !tens_shape_t:
        private TensShapeCtor
        private TensShapeCtorUnpack
        private TensShapePack
        private TensShapeSetDims
        private TensShapeSetGroups
        private TensShapeIsSet
        private TensShapeGetDims
        private TensShapeGetRank
        private TensShapeGetDimGroup
        private TensShapeGetGroup
        private TensShapeSameGroup
        private TensShapeNumGroups
        private TensShapeCompare
        private TensShapePrintIt
        public tens_shape_dtor
 !tens_header_t:
        private TensHeaderCtor
        private TensHeaderCtorUnpack
        private TensHeaderPack
        private TensHeaderAddShape
        private TensHeaderSetDims
        private TensHeaderSetGroups
        private TensHeaderIsSet
        private TensHeaderGetName
        private TensHeaderGetRank
        private TensHeaderGetSpec
        private TensHeaderGetDims
        private TensHeaderNumGroups
        private TensHeaderGetDimGroup
        private TensHeaderGetGroup
        private TensHeaderSameGroup
        private TensHeaderGetSignature
        private TensHeaderGetShape
        private TensHeaderCompare
        private TensHeaderPrintIt
        public tens_header_dtor
 !tens_simple_part_t:
        private TensSimplePartCtor
        private TensSimplePartCtorUnpack
        private TensSimplePartPack
        private TensSimplePartIsSet
        private TensSimplePartGetOffset
        private TensSimplePartGetLayout
        private TensSimplePartGetHeader
        public tens_simple_part_dtor
 !tens_layout_t:
        private TensLayoutSetLocation
        private TensLayoutIsSet
        private TensLayoutGetDataType
        private TensLayoutGetLayoutKind
        private TensLayoutGetBodyPtr
        private TensLayoutGetBodySize
 !tens_layout_fdims_t:
        private TensLayoutFdimsCtor
        private TensLayoutFdimsGetVolume
        private TensLayoutFdimsMap
        private TensLayoutFdimsExtract
        public tens_layout_fdims_dtor
 !tens_body_t:
        private TensBodyCtorBase
        private TensBodyIsSet
        private TensBodyAddSubtensor
        private TensBodySetLayout
        private TensBodySetLocation
        private TensBodyGetLayout
        private TensBodyGetNumSubtensors
        private TensBodyGetSubtensors
        public tens_body_dtor
 !tens_rcrsv_t:
        private TensRcrsvCtorSigna
        private TensRcrsvCtorHead
        private TensRcrsvIsSet
        private TensRcrsvAddSubtensor
        private TensRcrsvAddSubtensors
        private TensRcrsvSetShape
        private TensRcrsvSetLayout
        private TensRcrsvSetLocation
        private TensRcrsvGetHeader
        private TensRcrsvGetBody
        private TensRcrsvSplit
        public tens_rcrsv_dtor
!DATA:

       contains
!IMPLEMENTATION:
![Non-member]===========================================
        function valid_tensor_layout(layout) result(res)
!Returns TRUE if the tensor layout is valid.
         logical:: res
         integer(INTD), intent(in):: layout

         res=(layout.ge.0.and.layout.lt.TEREC_NUM_LAYOUTS)
         return
        end function valid_tensor_layout
!--------------------------------------------------------
        function cmp_tens_signatures(ts1,ts2) result(cmp)
!Comparator for tensor signatures.
         implicit none
         integer(INTD):: cmp                !out: result of comparison: {CMP_EQ,CMP_LT,CMP_GT,CMP_ER}
         class(*), intent(in), target:: ts1 !in: tensor signature 1
         class(*), intent(in), target:: ts2 !in: tensor signature 2
         class(tens_signature_t), pointer:: tsp1,tsp2

         tsp1=>NULL(); tsp2=>NULL()
         select type(ts1); class is(tens_signature_t); tsp1=>ts1; end select
         select type(ts2); class is(tens_signature_t); tsp2=>ts2; end select
         if(associated(tsp1).and.associated(tsp2)) then
          cmp=tsp1%compare(tsp2)
         else
          cmp=CMP_ER
         endif
         return
        end function cmp_tens_signatures
!-----------------------------------------------------
        function cmp_tens_headers(th1,th2) result(cmp)
!Comparator for tensor headers.
         implicit none
         integer(INTD):: cmp                !out: result of comparison: {CMP_EQ,CMP_LT,CMP_GT,CMP_ER}
         class(*), intent(in), target:: th1 !in: tensor header 1
         class(*), intent(in), target:: th2 !in: tensor header 2
         class(tens_header_t), pointer:: thp1,thp2

         thp1=>NULL(); thp2=>NULL()
         select type(th1); class is(tens_header_t); thp1=>th1; end select
         select type(th2); class is(tens_header_t); thp2=>th2; end select
         if(associated(thp1).and.associated(thp2)) then
          cmp=thp1%compare(thp2)
         else
          cmp=CMP_ER
         endif
         return
        end function cmp_tens_headers
!-----------------------------------------------------
        function print_tens_header_f(obj) result(ierr)
!Prints a tensor header (for GFC use).
         integer(INTD):: ierr
         class(*), intent(inout):: obj

         ierr=GFC_SUCCESS
         select type(obj)
         class is(tens_header_t)
          call obj%print_it(ierr)
         class default
          ierr=GFC_ACTION_FAILED
         end select
         return
        end function print_tens_header_f
![tens_signature_t]========================================================
        subroutine TensSignatureCtor(this,ierr,subspaces,tens_name,h_space)
!CTOR for tens_signature_t.
         implicit none
         class(tens_signature_t), intent(out):: this              !out: tensor signature
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
          if(len(tens_name).gt.0) then
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
         !else
          !errc=TEREC_INVALID_ARGS
          endif
         endif
         if(errc.ne.TEREC_SUCCESS) call tens_signature_dtor(this)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensSignatureCtor
!-----------------------------------------------------------
        subroutine TensSignatureCtorUnpack(this,packet,ierr)
!Ctor by unpacking.
         implicit none
         class(tens_signature_t), intent(out):: this     !out: tensor signature
         class(obj_pack_t), intent(inout):: packet       !in: packet
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: i,nd,errc
         integer(INTL):: sidx(1:MAX_TENSOR_RANK)
         character(TEREC_MAX_TENS_NAME_LEN):: tname
         logical:: pcn

         tname=' '
         call unpack_builtin(packet,nd,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,pcn,errc)
         if(errc.eq.PACK_SUCCESS) then
          if(nd.gt.0) then
           do i=1,nd
            call unpack_builtin(packet,sidx(i),errc); if(errc.ne.PACK_SUCCESS) exit
           enddo
          endif
          if(errc.eq.PACK_SUCCESS.and.pcn) call unpack_builtin(packet,tname,errc)
          if(errc.eq.PACK_SUCCESS) then
           if(nd.gt.0) then
            if(pcn) then
             call this%tens_signature_ctor(errc,sidx(1:nd),tname(1:len_trim(tname)))
            else
             call this%tens_signature_ctor(errc,sidx(1:nd))
            endif
           elseif(nd.eq.0) then
            if(pcn) then
             call this%tens_signature_ctor(errc,tens_name=tname(1:len_trim(tname)))
            else
             call this%tens_signature_ctor(errc)
            endif
           else
            errc=TEREC_ERROR
           endif
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensSignatureCtorUnpack
!-----------------------------------------------------
        subroutine TensSignaturePack(this,packet,ierr)
!Packs the object into a packet:
! + this%num_dims;
! + logical {TRUE/FALSE}: whether or not the tensor is named;
! + [OPTIONAL]: this%space_idx(1:this%num_dims);
! + [OPTIONAL]: this%char_name;
         implicit none
         class(tens_signature_t), intent(in):: this  !in: tensor signature
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: i,errc
         logical:: pcn

         if(this%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) then
           call pack_builtin(packet,this%num_dims,errc)
           if(errc.eq.PACK_SUCCESS) then
            pcn=allocated(this%char_name)
            call pack_builtin(packet,pcn,errc)
            if(errc.eq.PACK_SUCCESS) then
             if(this%num_dims.gt.0) then
              if(allocated(this%space_idx)) then
               do i=1,this%num_dims
                call pack_builtin(packet,this%space_idx(i),errc); if(errc.ne.PACK_SUCCESS) exit
               enddo
              else
               errc=TEREC_ERROR
              endif
             endif
             if(errc.eq.PACK_SUCCESS.and.pcn) call pack_builtin(packet,this%char_name,errc)
            endif
           endif
          endif
         else
          if(errc.eq.TEREC_SUCCESS) errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensSignaturePack
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
!------------------------------------------------------------------------------
        subroutine TensSignatureGetSpec(this,subspaces,num_dims,ierr,h_space_p)
!Returns the defining subspaces of the tensor (subspace multi-index).
         implicit none
         class(tens_signature_t), intent(in):: this                   !in: tensor signature
         integer(INTL), intent(inout):: subspaces(1:)                 !out: defining subspaces (their IDs)
         integer(INTD), intent(out):: num_dims                        !out: number of tensor dimensions
         integer(INTD), intent(out), optional:: ierr                  !out: error code
         class(h_space_t), intent(out), pointer, optional:: h_space_p !out: pointer to the underlying hierarchical vector space
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         if(this%is_set()) then
          num_dims=this%num_dims
          if(size(subspaces).ge.num_dims) then
           subspaces(1:num_dims)=this%space_idx(1:num_dims)
           if(present(h_space_p)) h_space_p=>this%h_space_p
          else
           errc=TEREC_UNABLE_COMPLETE
          endif
         else
          num_dims=-1; errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensSignatureGetSpec
!------------------------------------------------------------------
        function TensSignatureRelate(this,another) result(relation)
!Relates the tensor signature to another tensor signature.
         implicit none
         integer(INTD):: relation                      !out: relation: {CMP_EQ,CMP_CN,CMP_IN,CMP_OV,CMP_NC}
         class(tens_signature_t), intent(in):: this    !in: tensor signature 1
         class(tens_signature_t), intent(in):: another !in: tensor signature 2
         integer(INTD):: nl1,nl2,ch1,ch2,i,cmp,errc
         integer(INTL):: s1,s2

         errc=0
         if(this%is_set().and.another%is_set()) then
          relation=CMP_EQ
!Compare names:
          nl1=len(this%char_name); nl2=len(another%char_name)
          if(nl1.ne.nl2) relation=CMP_NC
          if(relation.eq.CMP_EQ) then
           do i=1,nl1
            ch1=iachar(this%char_name(i:i))
            ch2=iachar(another%char_name(i:i))
            if(ch1.ne.ch2) then; relation=CMP_NC; exit; endif
           enddo
           if(relation.eq.CMP_EQ) then
            if(associated(this%h_space_p)) then
             if(associated(this%h_space_p,another%h_space_p)) then
!Compare specs:
              if(this%num_dims.ne.another%num_dims) then
               relation=CMP_NC
              else
               do i=1,this%num_dims
                s1=this%space_idx(i); s2=another%space_idx(i)
                cmp=this%h_space_p%relate_subspaces(s1,s2,errc); if(errc.ne.0) exit
                if(cmp.eq.CMP_ER.or.cmp.eq.CMP_NC) then; relation=cmp; exit; endif
                if(cmp.ne.CMP_EQ) then
                 if(relation.eq.CMP_EQ) then
                  relation=cmp
                 else
                  if(cmp.ne.relation) then; relation=CMP_NC; exit; endif
                 endif
                endif
               enddo
               if(errc.ne.0) relation=CMP_ER
              endif
             else
              relation=CMP_ER !tensors with the same name cannot reside in differe hierarchical spaces
             endif
            else
             relation=CMP_ER
            endif
           endif
          endif
         else
          cmp=CMP_ER
         endif
         return
        end function TensSignatureRelate
!--------------------------------------------------------------
        function TensSignatureCompare(this,another) result(cmp)
!Compares the tensor signature with another tensor signature.
         implicit none
         integer(INTD):: cmp                           !out: comparison result: {CMP_EQ,CMP_LT,CMP_GT,CMP_ER}
         class(tens_signature_t), intent(in):: this    !in: tensor signature 1
         class(tens_signature_t), intent(in):: another !in: tensor signature 2
         integer(INTD):: nl1,nl2,ch1,ch2,i,errc
         integer(INTL):: s1,s2

         errc=0
         if(this%is_set().and.another%is_set()) then
          cmp=CMP_EQ
!Compare names:
          nl1=len(this%char_name); nl2=len(another%char_name)
          if(nl1.lt.nl2) then; cmp=CMP_LT; elseif(nl1.gt.nl2) then; cmp=CMP_GT; endif
          if(cmp.eq.CMP_EQ) then
           do i=1,nl1
            ch1=iachar(this%char_name(i:i))
            ch2=iachar(another%char_name(i:i))
            if(ch1.lt.ch2) then; cmp=CMP_LT; exit; elseif(ch1.gt.ch2) then; cmp=CMP_GT; exit; endif
           enddo
           if(cmp.eq.CMP_EQ) then
            if(associated(this%h_space_p)) then
             if(associated(this%h_space_p,another%h_space_p)) then
!Compare specs:
              if(this%num_dims.lt.another%num_dims) then
               cmp=CMP_LT
              elseif(this%num_dims.gt.another%num_dims) then
               cmp=CMP_GT
              else
               do i=1,this%num_dims
                s1=this%space_idx(i); s2=another%space_idx(i)
                cmp=this%h_space_p%compare_subspaces(s1,s2,errc)
                if(cmp.ne.CMP_EQ.or.errc.ne.0) exit
               enddo
               if(errc.ne.0) cmp=CMP_ER
              endif
             else
              cmp=CMP_ER !tensors with the same name cannot reside in differe hierarchical spaces
             endif
            else
             cmp=CMP_ER
            endif
           endif
          endif
         else
          cmp=CMP_ER
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
![tens_shape_t]============================================================
        subroutine TensShapeCtor(this,ierr,dim_extent,dim_group,group_spec)
!CTOR for tens_shape_t.
         implicit none
         class(tens_shape_t), intent(out):: this              !out: tensor shape
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTL), intent(in), optional:: dim_extent(1:) !in: tensor dimension extents (dimension extent of 0 means an unresolved dimension)
         integer(INTD), intent(in), optional:: dim_group(1:)  !in: dimension grouping: dim_group(x)=y means dimension x belongs to group y>0 (group 0 is default)
         integer(INTD), intent(in), optional:: group_spec(1:) !in: group specification: group_spec(x)=y means group x has restriction y (see on top)
         integer(INTD):: errc,i,j,k,m,n
         logical:: pr1,pr2

         errc=TEREC_SUCCESS
         n=0; if(present(dim_extent)) n=size(dim_extent)
         if(n.gt.0) then !true tensor
          do i=1,n; if(dim_extent(i).lt.0_INTL) then; errc=TEREC_INVALID_ARGS; exit; endif; enddo
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
        end subroutine TensShapeCtor
!-------------------------------------------------------
        subroutine TensShapeCtorUnpack(this,packet,ierr)
!Ctor by unpacking.
         implicit none
         class(tens_shape_t), intent(out):: this     !out: tensor shape
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: nd,ng,i,errc
         integer(INTL):: dim_ext(1:MAX_TENSOR_RANK)
         integer(INTD):: dim_grp(1:MAX_TENSOR_RANK),grp_spc(1:MAX_TENSOR_RANK)
         logical:: dpr

         call unpack_builtin(packet,nd,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,ng,errc)
         if(errc.eq.PACK_SUCCESS) then
          if(nd.gt.0) then
           call unpack_builtin(packet,dpr,errc)
           if(errc.eq.PACK_SUCCESS) then
            if(dpr) then
             do i=1,nd
              call unpack_builtin(packet,dim_ext(i),errc); if(errc.ne.PACK_SUCCESS) exit
             enddo
            else
             dim_ext(1:nd)=0_INTL !deferred dimension extents
            endif
            if(errc.eq.PACK_SUCCESS.and.ng.gt.0) then
             do i=1,nd
              call unpack_builtin(packet,dim_grp(i),errc); if(errc.ne.PACK_SUCCESS) exit
             enddo
             if(errc.eq.PACK_SUCCESS) then
              do i=1,ng
               call unpack_builtin(packet,grp_spc(i),errc); if(errc.ne.PACK_SUCCESS) exit
              enddo
             endif
            endif
            if(errc.eq.PACK_SUCCESS) then
             if(ng.gt.0) then
              call this%tens_shape_ctor(errc,dim_ext(1:nd),dim_grp(1:nd),grp_spc(1:ng))
             elseif(ng.eq.0) then
              call this%tens_shape_ctor(errc,dim_ext(1:nd))
             else
              errc=TEREC_ERROR
             endif
            endif
           endif
          elseif(nd.eq.0) then
           call this%tens_shape_ctor(errc)
          else
           errc=TEREC_ERROR
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensShapeCtorUnpack
!-------------------------------------------------
        subroutine TensShapePack(this,packet,ierr)
!Packs the object into a packet:
! + this%num_dims;
! + this%num_grps;
! + [OPTIONAL]: logical {TRUE/FALSE}: TRUE if this%dim_extent(:) is allocated;
! + [OPTIONAL]: this%dim_extent(1:this%num_dims);
! + [OPTIONAL]: this%dim_group(1:this%num_dims);
! + [OPTIONAL]: this%group_spec(1:this%num_grps).
         implicit none
         class(tens_shape_t), intent(in):: this      !in: tensor shape
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: i,errc,nd,ng
         logical:: dpr

         if(this%is_set(errc,nd,ng)) then
          if(errc.eq.TEREC_SUCCESS) then
           call pack_builtin(packet,nd,errc)
           if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,ng,errc)
           if(errc.eq.PACK_SUCCESS) then
            if(nd.gt.0) then
             dpr=allocated(this%dim_extent)
             call pack_builtin(packet,dpr,errc) !flag: presence of dimension extents
             if(errc.eq.PACK_SUCCESS) then
              if(dpr) then
               do i=1,nd
                call pack_builtin(packet,this%dim_extent(i),errc); if(errc.ne.PACK_SUCCESS) exit
               enddo
              endif
              if(errc.eq.PACK_SUCCESS.and.ng.gt.0) then
               do i=1,nd
                call pack_builtin(packet,this%dim_group(i),errc); if(errc.ne.PACK_SUCCESS) exit
               enddo
               if(errc.eq.PACK_SUCCESS) then
                do i=1,ng
                 call pack_builtin(packet,this%group_spec(i),errc); if(errc.ne.PACK_SUCCESS) exit
                enddo
               endif
              endif
             endif
            endif
           endif
          endif
         else
          if(errc.eq.TEREC_SUCCESS) errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensShapePack
!--------------------------------------------------------
        subroutine TensShapeSetDims(this,dim_extent,ierr)
!Sets dimension extents (if they have not been set previously).
!An attempt to set an already resolved tensor dimension will raise an error.
         implicit none
         class(tens_shape_t), intent(inout):: this   !inout: tensor shape
         integer(INTL), intent(in):: dim_extent(1:)  !in: dimension extents (those equal to 0 will be ignored)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: i,n,errc

         if(this%is_set(errc,n)) then
          if(size(dim_extent).eq.n) then
           do i=1,n
            if(dim_extent(i).gt.0) then
             if(this%dim_extent(i).eq.0) then !unresolved dimension
              this%dim_extent(i)=dim_extent(i)
             else
              errc=TEREC_INVALID_REQUEST; exit
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
        end subroutine TensShapeSetDims
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

         if(this%is_set(errc,n)) then
          if(errc.eq.TEREC_SUCCESS) then
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
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensShapeSetGroups
!------------------------------------------------------------------------------------
        function TensShapeIsSet(this,ierr,num_dims,num_groups,unresolved) result(res)
!Returns TRUE if the tensor shape is set, plus additional info.
         implicit none
         logical:: res                                     !out: result
         class(tens_shape_t), intent(in):: this            !in: tensor shape
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD), intent(out), optional:: num_dims   !out: number of dimensions
         integer(INTD), intent(out), optional:: num_groups !out: number of dimension groups
         integer(INTD), intent(out), optional:: unresolved !number of unresolved tensor dimensions
         integer(INTD):: i,errc

         errc=TEREC_SUCCESS
         res=(this%num_dims.ge.0)
         if(present(num_dims)) num_dims=this%num_dims
         if(present(num_groups)) num_groups=this%num_grps
         if(present(unresolved)) then
          unresolved=0
          if(res) then
           do i=1,this%num_dims
            if(this%dim_extent(i).eq.0_INTL) unresolved=unresolved+1 !unresolved tensor dimension
           enddo
          endif
         endif
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
!-------------------------------------------------------------------------------
        function TensShapeGetDimGroup(this,dimsn,group_restr,ierr) result(group)
!Returns the restriction group number and type of restriction for a specific tensor dimension.
         implicit none
         integer(INTD):: group                       !out: restriction group number (0:default group with no restrictions)
         class(tens_shape_t), intent(in):: this      !in: tensor shape
         integer(INTD), intent(in):: dimsn           !in: tensor dimension
         integer(INTD), intent(out):: group_restr    !out: type of index restriction (see top for TEREC_IND_RESTR_XXX)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: n,errc

         group=0; group_restr=TEREC_IND_RESTR_NONE
         if(this%is_set(errc,num_dims=n)) then
          if(errc.eq.TEREC_SUCCESS) then
           if(dimsn.gt.0.and.dimsn.le.n) then
            if(this%num_grps.gt.0) then
             group=this%dim_group(dimsn)
             if(group.gt.0) group_restr=this%group_spec(group)
            endif
           else
            errc=TEREC_INVALID_ARGS
           endif
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensShapeGetDimGroup
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
         integer(INTD):: cmp                            !out: comparison result: {CMP_EQ,CMP_LT,CMP_GT,CMP_ER}
         class(tens_shape_t), intent(in):: this         !in: tensor shape 1
         class(tens_shape_t), intent(in):: another      !in: tensor shape 2
         logical, intent(in), optional:: compare_groups !in: if FALSE, dimension groups will not be taken into account (defaults to TRUE)
         integer(INTD):: i,g1,g2,gmap(MAX_TENSOR_RANK)
         logical:: comp_grps

         if(present(compare_groups)) then; comp_grps=compare_groups; else; comp_grps=.TRUE.; endif
         if(this%is_set().and.another%is_set()) then
          if(this%num_dims.lt.another%num_dims) then
           cmp=CMP_LT
          elseif(this%num_dims.gt.another%num_dims) then
           cmp=CMP_GT
          else
           cmp=CMP_EQ
           do i=1,this%num_dims
            if(this%dim_extent(i).lt.another%dim_extent(i)) then
             cmp=CMP_LT; exit
            elseif(this%dim_extent(i).gt.another%dim_extent(i)) then
             cmp=CMP_GT; exit
            endif
           enddo
           if(cmp.eq.CMP_EQ.and.comp_grps) then
            if(this%num_grps.lt.another%num_grps) then
             cmp=CMP_LT
            elseif(this%num_grps.gt.another%num_grps) then
             cmp=CMP_GT
            else
             if(this%num_grps.gt.0) then
              if(this%num_grps.le.MAX_TENSOR_RANK) then
               gmap(1:this%num_grps)=0
               do i=1,this%num_dims
                g1=this%dim_group(i); g2=another%dim_group(i)
                if(g1.gt.0.and.g2.gt.0) then !both groups are non-trivial
                 if(this%group_spec(g1).lt.another%group_spec(g2)) then
                  cmp=CMP_LT; exit
                 elseif(this%group_spec(g1).gt.another%group_spec(g2)) then
                  cmp=CMP_GT; exit
                 else
                  if(gmap(g1).gt.0) then
                   if(gmap(g1).lt.g2) then
                    cmp=CMP_LT; exit
                   elseif(gmap(g1).gt.g2) then
                    cmp=CMP_GT; exit
                   endif
                  else
                   gmap(g1)=g2
                  endif
                 endif
                else
                 if(g1.lt.g2) then
                  cmp=CMP_LT; exit
                 elseif(g1.gt.g2) then
                  cmp=CMP_GT; exit
                 endif
                endif
               enddo
              else
               cmp=CMP_ER
              endif
             endif
            endif
           endif
          endif
         else
          cmp=CMP_ER
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
![tens_header_t]========================================================================================
        subroutine TensHeaderCtor(this,ierr,tens_name,subspaces,h_space,dim_extent,dim_group,group_spec)
!CTOR for tens_header_t. Each subsequent optional argument implies the existence of all preceding
!optional arguments, except <ierr> and <tens_name>. If no optional arguments are present, except
!maybe <tens_name> and/or <ierr>, a scalar header will be constructed. <dim_group> and <group_spec>
!must either be both present or both absent. More specifically:
! # Constructing a scalar tensor header: Do not pass any optional arguments except <tens_name> and/or <ierr>;
! # Constructing a true tensor header without shape: Pass only <subspaces>, <h_space>, and optionally <tens_name> and/or <ierr>;
! # Constructing a true tensor header with a shape: Pass <subspaces>, <h_space>, and <dim_extent> with all other arguments optional.
!   Note that it is ok to pass dimension extents equal to 0 for unresolved tensor dimensions (to be set later).
         implicit none
         class(tens_header_t), intent(out):: this             !out: tensor header
         integer(INTD), intent(out), optional:: ierr          !out: error code
         character(*), intent(in), optional:: tens_name       !in: alphanumeric_ tensor name
         integer(INTL), intent(in), optional:: subspaces(1:)  !in: subspace multi-index (specification): Length = tensor rank
         class(h_space_t), intent(in), target, optional:: h_space !in: underlying hierarchical vector space
         integer(INTL), intent(in), optional:: dim_extent(1:) !in: dimension extents: Length = tensor rank
         integer(INTD), intent(in), optional:: dim_group(1:)  !in: dimension restriction groups: Length = tensor rank
         integer(INTD), intent(in), optional:: group_spec(1:) !in: dimension restriction group specification
         integer(INTD):: errc,m
         logical:: pr_nam,pr_sub,pr_hsp,pr_dim,pr_grp,pr_grs

         errc=TEREC_SUCCESS
         pr_nam=present(tens_name)
         pr_sub=present(subspaces)
         pr_hsp=present(h_space)
         pr_dim=present(dim_extent)
         pr_grp=present(dim_group)
         pr_grs=present(group_spec)
         if(pr_sub.and.(.not.pr_hsp)) errc=TEREC_INVALID_ARGS
         if(pr_dim.and.(.not.pr_sub)) errc=TEREC_INVALID_ARGS
         if((pr_grp.or.pr_grs).and.(.not.pr_dim)) errc=TEREC_INVALID_ARGS
         if((pr_grp.and.(.not.pr_grs)).or.(pr_grs.and.(.not.pr_grp))) errc=TEREC_INVALID_ARGS
         if(pr_sub.and.pr_dim) then
          m=size(subspaces); if(m.ne.size(dim_extent)) errc=TEREC_INVALID_ARGS
          if(m.le.0) then; pr_sub=.FALSE.; pr_dim=.FALSE.; endif
         endif
         if(errc.eq.TEREC_SUCCESS) then
 !tensor signature:
          if(pr_sub) then !explicit tensor
           if(pr_nam) then
            if(pr_hsp) then
             call this%signature%tens_signature_ctor(errc,subspaces,tens_name,h_space)
            else
             call this%signature%tens_signature_ctor(errc,subspaces,tens_name) !`this will never happen
            endif
           else
            if(pr_hsp) then
             call this%signature%tens_signature_ctor(errc,subspaces,h_space=h_space)
            else
             call this%signature%tens_signature_ctor(errc,subspaces) !`this will never happen
            endif
           endif
          else !scalar
           if(pr_nam) then
            call this%signature%tens_signature_ctor(errc,tens_name=tens_name)
           else
            call this%signature%tens_signature_ctor(errc)
           endif
           if(errc.eq.TEREC_SUCCESS) call this%shape%tens_shape_ctor(errc)
          endif
 !tensor shape (optional):
          if(errc.eq.TEREC_SUCCESS) then
           if(pr_dim) then !only explicit tensors
            if(pr_grp) then
             call this%shape%tens_shape_ctor(errc,dim_extent,dim_group,group_spec)
            else
             call this%shape%tens_shape_ctor(errc,dim_extent)
            endif
           endif
          endif
         endif
         if(errc.ne.TEREC_SUCCESS) call tens_header_dtor(this)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensHeaderCtor
!--------------------------------------------------------
        subroutine TensHeaderCtorUnpack(this,packet,ierr)
!Ctor by unpacking.
         implicit none
         class(tens_header_t), intent(out):: this    !out: tensor header
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         logical:: shaped

         call unpack_builtin(packet,shaped,errc)
         if(errc.eq.PACK_SUCCESS) call this%signature%tens_signature_ctor(packet,errc)
         if(errc.eq.PACK_SUCCESS.and.shaped) call this%shape%tens_shape_ctor(packet,errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensHeaderCtorUnpack
!--------------------------------------------------
        subroutine TensHeaderPack(this,packet,ierr)
!Packs the object into a packet:
! + logical {TRUE|FALSE}: whether the tensor header shape is set or not
! + this%signature;
! + [OPTIONAL]: this%shape;
         implicit none
         class(tens_header_t), intent(in):: this     !in: tensor header
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         logical:: shpd

         if(this%is_set(errc,shaped=shpd)) then
          if(errc.eq.TEREC_SUCCESS) then
           call pack_builtin(packet,shpd,errc)
           if(errc.eq.PACK_SUCCESS) call this%signature%pack(packet,errc)
           if(errc.eq.PACK_SUCCESS.and.shpd) call this%shape%pack(packet,errc)
          endif
         else
          if(errc.eq.TEREC_SUCCESS) errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensHeaderPack
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
                call this%shape%tens_shape_ctor(errc,dim_extent,dim_group,group_spec)
               else
                errc=TEREC_INVALID_ARGS
               endif
              else
               if(present(group_spec)) then
                errc=TEREC_INVALID_ARGS
               else
                call this%shape%tens_shape_ctor(errc,dim_extent)
               endif
              endif
             else !scalar shape
              call this%shape%tens_shape_ctor(errc)
             endif
            else
             errc=TEREC_INVALID_ARGS
            endif
           else !scalar shape
            if(n.ne.0.or.present(dim_group).or.present(group_spec)) then
             errc=TEREC_INVALID_ARGS
            else
             call this%shape%tens_shape_ctor(errc)
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
!---------------------------------------------------------
        subroutine TensHeaderSetDims(this,dim_extent,ierr)
!Sets dimension extents (if they have not been set previously).
!An attempt to reset an already resolved tensor dimension will raise an error.
         implicit none
         class(tens_header_t), intent(inout):: this  !inout: tensor header
         integer(INTL), intent(in):: dim_extent(1:)  !in: dimension extents (those equal to 0 will be ignored)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer:: errc

         if(this%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) call this%shape%set_dims(dim_extent,errc)
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensHeaderSetDims
!---------------------------------------------------------------------
        subroutine TensHeaderSetGroups(this,dim_group,group_spec,ierr)
!Sets index restriction groups if they have not been previously set.
!An attempt to redefine existing groups will raise an error.
         implicit none
         class(tens_header_t), intent(inout):: this  !inout: tensor header
         integer(INTD), intent(in):: dim_group(1:)   !in: index restriction groups
         integer(INTD), intent(in):: group_spec(1:)  !in: restriction groups specification
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer:: errc

         if(this%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) call this%shape%set_groups(dim_group,group_spec,errc)
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensHeaderSetGroups
!--------------------------------------------------------------------------------------------
        function TensHeaderIsSet(this,ierr,num_dims,num_groups,shaped,unresolved) result(res)
!Returns TRUE if the tensor header is set (with or without shape), plus additional info.
         implicit none
         logical:: res                                     !out: result
         class(tens_header_t), intent(in):: this           !in: tensor header
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD), intent(out), optional:: num_dims   !out: number of dimensions
         integer(INTD), intent(out), optional:: num_groups !out: number of restricted dimension groups
         logical, intent(out), optional:: shaped           !out: TRUE if the tensor shape is set
         integer(INTD), intent(out), optional:: unresolved !out: number of unresolved tensor dimensions
         integer(INTD):: errc,nd,ng,unres
         logical:: shpd

         res=this%signature%is_set(errc,num_dims=nd)
         if(present(num_dims)) num_dims=nd
         if(res.and.errc.eq.TEREC_SUCCESS) then
          shpd=this%shape%is_set(errc,num_dims=nd,num_groups=ng,unresolved=unres)
         else
          shpd=.FALSE.; ng=0; unres=-1
         endif
         if(present(num_groups)) num_groups=ng
         if(present(shaped)) shaped=shpd
         if(present(unresolved)) unresolved=unres
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
!---------------------------------------------------------------------------
        subroutine TensHeaderGetSpec(this,subspaces,num_dims,ierr,h_space_p)
!Returns the defining subspaces of the tensor (subspace multi-index).
         implicit none
         class(tens_header_t), intent(in):: this                      !in: tensor header
         integer(INTL), intent(inout):: subspaces(1:)                 !out: defining subspaces (their IDs)
         integer(INTD), intent(out):: num_dims                        !out: number of tensor dimensions
         integer(INTD), intent(out), optional:: ierr                  !out: error code
         class(h_space_t), intent(out), pointer, optional:: h_space_p !out: pointer to the underlying hierarchical vector space
         integer(INTD):: errc

         if(present(h_space_p)) then
          call this%signature%get_spec(subspaces,num_dims,errc,h_space_p)
         else
          call this%signature%get_spec(subspaces,num_dims,errc)
         endif
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
!--------------------------------------------------------------------------------
        function TensHeaderGetDimGroup(this,dimsn,group_restr,ierr) result(group)
!Returns the restriction group for a specific tensor dimension.
         implicit none
         integer(INTD):: group                       !out: restriction group
         class(tens_header_t), intent(in):: this     !in: tensor header
         integer(INTD), intent(in):: dimsn           !in: specific tensor dimension
         integer(INTD), intent(out):: group_restr    !out: type of restriction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         group=this%shape%get_dim_group(dimsn,group_restr,errc)
         if(present(ierr)) ierr=errc
         return
        end function TensHeaderGetDimGroup
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
         integer(INTD):: cmp                            !out: comparison result: {CMP_EQ,CMP_LT,CMP_GT,CMP_ER}
         class(tens_header_t), intent(in):: this        !in: tensor header 1
         class(tens_header_t), intent(in):: another     !in: tensor header 2
         logical, intent(in), optional:: compare_groups !in: if FALSE, the shape dimension groups will not be taken into account (defaults to TRUE)

         cmp=this%signature%compare(another%signature)
         if(cmp.eq.CMP_EQ) then
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
![tens_simple_part_t]================================================
        subroutine TensSimplePartCtor(this,header,layout,offset,ierr)
!CTOR for tens_simple_part_t.
         implicit none
         class(tens_simple_part_t), intent(out):: this !out: tensor simple part
         class(tens_header_t), intent(in):: header     !in: tensor header
         integer(INTD), intent(in):: layout            !in: simple storage layout: {TEREC_LAY_FDIMS,TEREC_LAY_CDIMS}
         integer(INTL), intent(in):: offset            !in: offset in the parental tensor block
         integer(INTD), intent(out), optional:: ierr   !out: error code
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
        end subroutine TensSimplePartCtor
!------------------------------------------------------------
        subroutine TensSimplePartCtorUnpack(this,packet,ierr)
!Ctor by unpacking.
         implicit none
         class(tens_simple_part_t), intent(out):: this !out: tensor simple part
         class(obj_pack_t), intent(inout):: packet     !inout: packet
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: errc

         call this%header%tens_header_ctor(packet,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%offset,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%layout,errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensSimplePartCtorUnpack
!------------------------------------------------------
        subroutine TensSimplePartPack(this,packet,ierr)
!Packs the object into a packet:
! + this%header;
! + this%offset;
! + this%layout.
         implicit none
         class(tens_simple_part_t), intent(in):: this !in: tensor simple part
         class(obj_pack_t), intent(inout):: packet    !inout: packet
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         if(this%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) then
           call this%header%pack(packet,errc)
           if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%offset,errc)
           if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%layout,errc)
          endif
         else
          if(errc.eq.TEREC_SUCCESS) errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensSimplePartPack
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
![tens_layout_t]==============================================
        subroutine TensLayoutSetLocation(this,data_descr,ierr)
!Sets the phyiscal location of the tensor body.
         implicit none
         class(tens_layout_t), intent(inout):: this  !inout: tensor body layout
         class(DataDescr_t), intent(in):: data_descr !in: DDSS data descriptor for the tensor body
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         if(this%layout.ne.TEREC_LAY_NONE) then
          if(data_descr%is_set().and.(.not.allocated(this%data_descr))) then
           allocate(this%data_descr,SOURCE=data_descr,STAT=errc); if(errc.ne.0) errc=TEREC_MEM_ALLOC_FAILED
          else
           errc=TEREC_INVALID_ARGS
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensLayoutSetLocation
!--------------------------------------------------------------
        function TensLayoutIsSet(this,ierr,located) result(res)
!Returns TRUE if the tensor layout is set.
         implicit none
         logical:: res                               !out: result
         class(tens_layout_t), intent(in):: this     !in: tensor layout
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(out), optional:: located    !out: returns TRUE if the tensor body has physical location
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         res=(this%layout.ne.TEREC_LAY_NONE)
         if(present(located)) located=allocated(this%data_descr)
         if(present(ierr)) ierr=errc
         return
        end function TensLayoutIsSet
!------------------------------------------------------------------
        function TensLayoutGetDataType(this,ierr) result(data_type)
!Returns the data type of stored tensor elements.
         implicit none
         integer(INTD):: data_type                   !out: data type
         class(tens_layout_t), intent(in):: this     !in: tensor layout
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%is_set(errc)) then
          data_type=this%data_type
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensLayoutGetDataType
!-----------------------------------------------------------------
        function TensLayoutGetLayoutKind(this,ierr) result(layout)
!Returns the tensor layout kind.
         implicit none
         integer(INTD):: layout                      !out: tensor layout kind
         class(tens_layout_t), intent(in):: this     !in: tensor layout
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         layout=this%layout
         if(present(ierr)) ierr=errc
         return
        end function TensLayoutGetLayoutKind
!--------------------------------------------------------------
        function TensLayoutGetBodyPtr(this,ierr) result(body_p)
!Returns a C pointer to the stored tensor body.
         implicit none
         type(C_PTR):: body_p                        !out: C pointer to the tensor body
         class(tens_layout_t), intent(in):: this     !in: tensor layout
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         logical:: locd

         body_p=C_NULL_PTR
         if(this%is_set(errc,locd)) then
          if(locd) then
           if(this%data_descr%is_set()) then
            body_p=this%data_descr%get_data_ptr(errc)
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
        end function TensLayoutGetBodyPtr
!------------------------------------------------------------------
        function TensLayoutGetBodySize(this,ierr) result(body_size)
!Returns the size of the stored tensor body in bytes.
         implicit none
         integer(INTL):: body_size                   !out: tensor body size in bytes
         class(tens_layout_t), intent(in):: this     !in: tensor layout
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         logical:: locd

         body_size=0_INTL
         if(this%is_set(errc,locd)) then
          if(locd) then
           if(this%data_descr%is_set()) then
            body_size=this%data_descr%data_size(errc)
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
        end function TensLayoutGetBodySize
![tens_layout_fdims_t]=================================================
        subroutine TensLayoutFdimsCtor(this,tens_header,data_type,ierr)
!Constructs the "Fortran dimension led" tensor body layout.
         implicit none
         class(tens_layout_fdims_t), intent(out):: this         !out: tensor body layout
         class(tens_header_t), intent(in), target:: tens_header !in: tensor header (logical tensor spec for which the physical layout is constructed)
         integer(INTD), intent(in):: data_type                  !in: data type for tensor elements: {R4,R8,C4,C8}
         integer(INTD), intent(out), optional:: ierr            !out: error code
         integer(INTD):: errc,ds,unres
         logical:: shpd

         if(tens_header%is_set(errc,shaped=shpd,unresolved=unres)) then
          if(errc.eq.TEREC_SUCCESS) then
           if(shpd.and.unres.eq.0) then
            if(tens_valid_data_kind(data_type,ds).eq.YEP) then
             if(ds.gt.0) then
              this%layout=TEREC_LAY_FDIMS
              this%data_type=data_type
              this%header=>tens_header
             else
              errc=TEREC_INVALID_ARGS
             endif
            else
             errc=TEREC_INVALID_ARGS
            endif
           else
            errc=TEREC_INVALID_REQUEST
           endif
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensLayoutFdimsCtor
!----------------------------------------------------------
        function TensLayoutFdimsGetVolume(this) result(vol)
         implicit none
         integer(INTL):: vol                           !out: physical volume of the tensor body (number of stored tensor elements)
         class(tens_layout_fdims_t), intent(in):: this !in: tensor layout
         integer(INTL):: dims(1:MAX_TENSOR_RANK)
         integer(INTD):: i,n,errc

         vol=0_INTL
         if(associated(this%header)) then
          call this%header%get_dims(dims,n,errc)
          if(errc.eq.TEREC_SUCCESS) then
           vol=1_INTL; do i=1,n; vol=vol*dims(i); enddo
          else
           vol=-1_INTL
          endif
         endif
         return
        end function TensLayoutFdimsGetVolume
!-------------------------------------------------------------
        function TensLayoutFdimsMap(this,mlndx) result(offset)
!Given a multi-index position of the tensor element inside tensor body,
!returns its linear offset in the tensor body. The multi-index position
!is specified relative to the tensor body and index numeration starts from 1.
         implicit none
         integer(INTL):: offset                        !out: linear tensor element offset
         class(tens_layout_fdims_t), intent(in):: this !in: tensor layout
         integer(INTL), intent(in):: mlndx(1:)         !in: multi-index position of the tensor element
         integer(INTL):: dims(1:MAX_TENSOR_RANK)
         integer(INTD):: i,n,errc

         offset=-1_INTL
         if(associated(this%header)) then
          call this%header%get_dims(dims,n,errc)
          if(errc.eq.TEREC_SUCCESS) then
           if(n.gt.0) then
            offset=(mlndx(n)-1_INTL)
            do i=n-1,1,-1
             offset=offset*dims(i)+(mlndx(i)-1_INTL)
            enddo
           else
            offset=0_INTL
           endif
          endif
         endif
         return
        end function TensLayoutFdimsMap
!-------------------------------------------------------------------
        subroutine TensLayoutFdimsExtract(this,num_parts,parts,ierr)
         implicit none
         class(tens_layout_fdims_t), intent(in):: this !in: tensor layout
         integer(INTL), intent(out):: num_parts        !out: number of simple parts extracted from the tensor layout
         type(list_bi_t), intent(inout):: parts        !list of the simple parts extracted from the tensor layout
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: i,errc
         logical:: locd
         type(list_iter_t):: lit
         type(tens_simple_part_t):: tsp

         num_parts=0_INTL
         if(this%is_set(errc,locd)) then
          if(locd) then
           errc=lit%init(parts)
           if(errc.eq.GFC_SUCCESS) then
            errc=lit%get_status()
            if(errc.eq.GFC_IT_EMPTY) then
             call tsp%tens_simple_part_ctor(this%header,TEREC_LAY_FDIMS,0_INTL,errc)
             if(errc.eq.TEREC_SUCCESS) then
              errc=lit%append(tsp)
              if(errc.eq.GFC_SUCCESS) then
               num_parts=num_parts+1_INTL
              else
               errc=TEREC_UNABLE_COMPLETE
              endif
             endif
             call tens_simple_part_dtor(tsp)
            else
             errc=TEREC_INVALID_ARGS
            endif
            i=lit%release()
           else
            errc=TEREC_UNABLE_COMPLETE
           endif
          else
           errc=TEREC_INVALID_REQUEST
          endif
         else
          errc=TEREC_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensLayoutFdimsExtract
!---------------------------------------------
        subroutine tens_layout_fdims_dtor(this)
         implicit none
         type(tens_layout_fdims_t):: this

         if(allocated(this%data_descr)) deallocate(this%data_descr)
         this%header=>NULL()
         this%data_type=NO_TYPE
         this%layout=TEREC_LAY_NONE
         return
        end subroutine tens_layout_fdims_dtor
![tens_body_t]================================
        subroutine TensBodyCtorBase(this,ierr)
!Default ctor.
         implicit none
         class(tens_body_t), intent(out):: this      !out: tensor body
         integer(INTD), intent(out), optional:: ierr !out: error code

         if(present(ierr)) ierr=TEREC_SUCCESS
         return
        end subroutine TensBodyCtorBase
!------------------------------------------------------------------
        function TensBodyIsSet(this,ierr,layed,located) result(res)
!Returns TRUE if the tensor body is set (plus additional info), that is,
!if it consists of at least one subtensor.
         implicit none
         logical:: res                               !out: result
         class(tens_body_t), intent(in):: this       !in: tensor body
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(out), optional:: layed      !out: TRUE if the tensor physical layout has been set
         logical, intent(out), optional:: located    !out: TRUE if the tensor body has been physically mapped
         integer(INTD):: errc
         logical:: layd,locd

         errc=TEREC_SUCCESS; layd=.FALSE.; locd=.FALSE.
         res=(this%num_subtensors.gt.0)
         if(res) then
          if(allocated(this%layout)) layd=this%layout%is_set(errc,locd)
         endif
         if(present(layed)) layed=layd
         if(present(located)) located=locd
         if(present(ierr)) ierr=errc
         return
        end function TensBodyIsSet
!-----------------------------------------------------------
        subroutine TensBodyAddSubtensor(this,subtensor,ierr)
!Registers a constituent subtensor by providing its tensor header.
         implicit none
         class(tens_body_t), intent(inout):: this     !inout: tensor body
         class(tens_header_t), intent(in):: subtensor !in: constituent subtensor
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: i,errc
         type(list_iter_t):: lit

         errc=lit%init(this%subtensors)
         if(errc.eq.GFC_SUCCESS) then
          errc=lit%append(subtensor)
          if(errc.eq.GFC_SUCCESS) then
           this%num_subtensors=this%num_subtensors+1
          else
           errc=TEREC_UNABLE_COMPLETE
          endif
          i=lit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=TEREC_ERROR
         else
          errc=TEREC_UNABLE_COMPLETE
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensBodyAddSubtensor
!--------------------------------------------------------------------
        subroutine TensBodySetLayout(this,layout_kind,data_type,ierr)
!Sets tensor body storage layout.
         implicit none
         class(tens_body_t), intent(inout):: this    !inout: tensor body
         integer(INTD), intent(in):: layout_kind     !in: layout kind
         integer(INTD), intent(in):: data_type       !in: numeric data type: {R4,R8,C4,C8}
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         logical:: layd,locd
         type(list_iter_t):: lit
         class(tens_header_t), pointer:: thp
         class(*), pointer:: up

         if(this%is_set(errc,layd,locd)) then
          if(errc.eq.TEREC_SUCCESS.and.(.not.layd)) then
           lselect: select case(layout_kind)
           case(TEREC_LAY_RECUR) !all constituent subtensors are mapped sequentially to a contiguous chunk of local memory
            stop !`Implement
           case(TEREC_LAY_FDIMS) !a single subtensor is mapped as "Fortran-dimension-led"
            if(this%num_subtensors.eq.1) then
             errc=lit%init(this%subtensors); if(errc.ne.GFC_SUCCESS) then; errc=TEREC_UNABLE_COMPLETE; exit lselect; endif
             up=>lit%get_value(errc); if(errc.ne.GFC_SUCCESS) then; errc=TEREC_UNABLE_COMPLETE; exit lselect; endif
             select type(up); class is(tens_header_t); thp=>up; end select
             if(.not.associated(thp)) then; errc=TEREC_UNABLE_COMPLETE; exit lselect; endif
             errc=lit%release(); if(errc.ne.GFC_SUCCESS) then; errc=TEREC_UNABLE_COMPLETE; exit lselect; endif
             allocate(tens_layout_fdims_t::this%layout,STAT=errc)
             if(errc.eq.0) then
              select type(layout=>this%layout)
              class is(tens_layout_fdims_t)
               call layout%tens_layout_fdims_ctor(thp,data_type,errc)
              class default
               errc=TEREC_ERROR
              end select
             else
              errc=TEREC_MEM_ALLOC_FAILED
             endif
            else
             errc=TEREC_INVALID_REQUEST
            endif
           case(TEREC_LAY_CDIMS)
            if(this%num_subtensors.eq.1) then
             stop !`Implement
            else
             errc=TEREC_INVALID_REQUEST
            endif
           case(TEREC_LAY_DSYMM)
            if(this%num_subtensors.eq.1) then
             stop !`Implement
            else
             errc=TEREC_INVALID_REQUEST
            endif
           case(TEREC_LAY_BRICK)
            if(this%num_subtensors.eq.1) then
             stop !`Implement
            else
             errc=TEREC_INVALID_REQUEST
            endif
           case(TEREC_LAY_BSYMM)
            if(this%num_subtensors.eq.1) then
             stop !`Implement
            else
             errc=TEREC_INVALID_REQUEST
            endif
           case(TEREC_LAY_SPARS)
            if(this%num_subtensors.eq.1) then
             stop !`Implement
            else
             errc=TEREC_INVALID_REQUEST
            endif
           case(TEREC_LAY_NONE) !destroy existing layout
            if(allocated(this%layout)) then
             deallocate(this%layout,STAT=errc); if(errc.ne.0) errc=TEREC_MEM_FREE_FAILED
            endif
           case default
            errc=TEREC_INVALID_ARGS
           end select lselect
          else
           errc=TEREC_INVALID_REQUEST
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensBodySetLayout
!-----------------------------------------------------------
        subroutine TensBodySetLocation(this,data_descr,ierr)
!Sets the physical location of the tensor body via a DDSS data descriptor.
         implicit none
         class(tens_body_t), intent(inout):: this    !inout: tensor body
         class(DataDescr_t), intent(in):: data_descr !in: DDSS data descriptor for tensor body data
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         logical:: layd,locd

         if(this%is_set(errc,layd,locd)) then
          if(errc.eq.TEREC_SUCCESS.and.layd.and.(.not.locd)) then
           call this%layout%set_location(data_descr,errc)
          else
           errc=TEREC_INVALID_REQUEST
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensBodySetLocation
!-------------------------------------------------------------
        function TensBodyGetLayout(this,ierr) result(layout_p)
!Returns a pointer to the tensor body storage layout.
         implicit none
         class(tens_layout_t), pointer:: layout_p      !out: pointer to the tensor body layout
         class(tens_body_t), intent(in), target:: this !in: tensor body
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: errc
         logical:: layd,locd

         layout_p=>NULL()
         if(this%is_set(errc,layd,locd)) then
          if(errc.eq.TEREC_SUCCESS.and.layd) then
           layout_p=>this%layout
          else
           errc=TEREC_INVALID_REQUEST
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensBodyGetLayout
!--------------------------------------------------------------------------
        function TensBodyGetNumSubtensors(this,ierr) result(num_subtensors)
!Returns the total number of constituent subtensors.
         implicit none
         integer(INTD):: num_subtensors              !out: total number of constituent subtensors
         class(tens_body_t), intent(in):: this       !in: tensor body
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         num_subtensors=this%num_subtensors
         if(.not.this%is_set(errc)) errc=TEREC_INVALID_REQUEST
         if(present(ierr)) ierr=errc
         return
        end function TensBodyGetNumSubtensors
!-------------------------------------------------------------------------
        function TensBodyGetSubtensors(this,ierr) result(subtensor_list_p)
!Returns a pointer to the list of constituent subtensors (each subtensor is represented by a tensor header)
         implicit none
         type(list_bi_t), pointer:: subtensor_list_p   !out: pointer to the subtensor list
         class(tens_body_t), intent(in), target:: this !in: tensor body
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: errc

         subtensor_list_p=>NULL()
         if(this%is_set(errc)) then
          subtensor_list_p=>this%subtensors
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensBodyGetSubtensors
!--------------------------------------
        subroutine tens_body_dtor(this)
         implicit none
         type(tens_body_t):: this
         type(list_iter_t):: lit
         integer(INTD):: errc

         if(allocated(this%layout)) deallocate(this%layout)
         errc=lit%init(this%subtensors); if(errc.eq.GFC_SUCCESS) errc=lit%delete_all()
         errc=lit%release()
         this%num_subtensors=0
         return
        end subroutine tens_body_dtor
![tens_rcrsv_t]=============================================================================================
        subroutine TensRcrsvCtorSigna(this,tens_name,subspaces,h_space,ierr,dim_extent,dim_group,group_spec)
!Constructs a tensor by specifying a tensor signature and optionally a shape.
!See TensHeaderCtor for restrictions.
         implicit none
         class(tens_rcrsv_t), intent(out):: this              !out: tensor
         character(*), intent(in):: tens_name                 !in: alphanumeric_ tensor name
         integer(INTL), intent(in):: subspaces(1:)            !in: subspace multi-index (signature)
         class(h_space_t), intent(in), target:: h_space       !in: hierarchical vector space (externally persistent)
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTL), intent(in), optional:: dim_extent(1:) !in: dimension extents (those equal to 0 are unresolved)
         integer(INTD), intent(in), optional:: dim_group(1:)  !in: dimension restriction groups
         integer(INTD), intent(in), optional:: group_spec(1:) !in: restriction groups specification
         integer(INTD):: errc

         if(present(dim_extent)) then !signature + shape (possibly with unresolved dimensions, those equal to 0)
          if(present(dim_group)) then
           if(present(group_spec)) then
            call this%header%tens_header_ctor(errc,tens_name,subspaces,h_space,dim_extent,dim_group,group_spec)
           else
            errc=TEREC_INVALID_ARGS
           endif
          else
           if(.not.present(group_spec)) then
            call this%header%tens_header_ctor(errc,tens_name,subspaces,h_space,dim_extent)
           else
            errc=TEREC_INVALID_ARGS
           endif
          endif
         else !signature only
          if(.not.(present(dim_group).or.present(group_spec))) then
           call this%header%tens_header_ctor(errc,tens_name,subspaces,h_space)
          else
           errc=TEREC_INVALID_ARGS
          endif
         endif
         if(errc.ne.TEREC_SUCCESS) call tens_rcrsv_dtor(this)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensRcrsvCtorSigna
!-----------------------------------------------------
        subroutine TensRcrsvCtorHead(this,header,ierr)
!Constructs a tensor by providing a preset tensor header.
         implicit none
         class(tens_rcrsv_t), intent(out):: this     !out: tensor
         type(tens_header_t), intent(in):: header    !in: existing preset tensor header
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(header%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) this%header=header
         else
          errc=TEREC_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensRcrsvCtorHead
!-------------------------------------------------------------------------------------
        function TensRcrsvIsSet(this,ierr,shaped,unresolved,layed,located) result(res)
!Returns TRUE if the tensor is set, plus additional info.
         implicit none
         logical:: res                                     !out: result
         class(tens_rcrsv_t), intent(in):: this            !in: tensor
         integer(INTD), intent(out), optional:: ierr       !out: error code
         logical, intent(out), optional:: shaped           !out: TRUE if tensor shape is set (even with unresolved dimensions)
         integer(INTD), intent(out), optional:: unresolved !out: number of unresolved tensor dimensions
         logical, intent(out), optional:: layed            !out: TRUE if the tensor body storage layout is set
         logical, intent(out), optional:: located          !out: TRUE if the physical location for tensor body data is set
         integer(INTD):: errc,unres
         logical:: shpd,layd,locd

         res=this%header%is_set(errc,shaped=shpd,unresolved=unres)
         if(present(shaped)) shaped=shpd
         if(present(unresolved)) unresolved=unres
         if(res) then
          shpd=this%body%is_set(errc,layd,locd)
          if(present(layed)) layed=layd
          if(present(located)) located=locd
         else
          if(present(layed)) layed=.FALSE.
          if(present(located)) located=.FALSE.
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensRcrsvIsSet
!------------------------------------------------------------
        subroutine TensRcrsvAddSubtensor(this,subtensor,ierr)
!Registers a constituent subtensor by providing its tensor header.
         implicit none
         class(tens_rcrsv_t), intent(inout):: this    !inout: tensor
         class(tens_header_t), intent(in):: subtensor !in: subtensor
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         if(this%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) call this%body%add_subtensor(subtensor,errc)
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensRcrsvAddSubtensor
!--------------------------------------------------------------
        subroutine TensRcrsvAddSubtensors(this,subtensors,ierr)
!Registers constituent subtensors by providing a list of their tensor headers.
         implicit none
         class(tens_rcrsv_t), intent(inout):: this   !inout: tensor
         type(list_bi_t), intent(in):: subtensors    !in: list of subtensors (their tensor headers)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: i,errc
         type(list_iter_t):: lit
         class(tens_header_t), pointer:: thp
         class(*), pointer:: up

         if(this%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) then
           errc=lit%init(subtensors)
           if(errc.eq.GFC_SUCCESS) then
            errc=GFC_IT_ACTIVE; thp=>NULL()
            do while(errc.eq.GFC_IT_ACTIVE)
             up=>lit%get_value(errc); if(errc.ne.GFC_SUCCESS) exit
             select type(up); class is(tens_header_t); thp=>up; end select
             if(.not.associated(thp)) then; errc=TEREC_ERROR; exit; endif
             call this%body%add_subtensor(thp,errc); if(errc.ne.TEREC_SUCCESS) exit
             thp=>NULL(); up=>NULL()
             errc=lit%scanp(return_each=.TRUE.,skip_current=.TRUE.)
            enddo
            if(errc.eq.GFC_IT_DONE) then
             errc=lit%release()
            else
             errc=lit%release()
             errc=TEREC_UNABLE_COMPLETE
            endif
           else
            errc=TEREC_UNABLE_COMPLETE
           endif
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensRcrsvAddSubtensors
!------------------------------------------------------------------------------
        subroutine TensRcrsvSetShape(this,dim_extent,ierr,dim_group,group_spec)
!Sets the tensor shape. If the shape is already set, it will try to set unresolved
!tensor dimensions and optionally index restriction groups (if not yet set).
!If the shape is unset, it will be constructed anew.
         implicit none
         class(tens_rcrsv_t), intent(inout):: this            !inout: tensor
         integer(INTL), intent(in):: dim_extent(1:)           !in: dimension extents (those equal to 0 are unresolved)
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTD), intent(in), optional:: dim_group(1:)  !in: index restriction groups
         integer(INTD), intent(in), optional:: group_spec(1:) !in: restriction groups specification
         integer(INTD):: errc
         logical:: shpd

         if(this%is_set(errc,shaped=shpd)) then
          if(errc.eq.TEREC_SUCCESS) then
           if(shpd) then
            call this%header%set_dims(dim_extent,errc)
            if(errc.eq.TEREC_SUCCESS) then
             if(present(dim_group)) then
              if(present(group_spec)) then
               call this%header%set_groups(dim_group,group_spec,errc)
              else
               errc=TEREC_INVALID_ARGS
              endif
             else
              if(present(group_spec)) errc=TEREC_INVALID_ARGS
             endif
            endif
           else
            if(present(dim_group)) then
             if(present(group_spec)) then
              call this%header%add_shape(errc,dim_extent,dim_group,group_spec)
             else
              errc=TEREC_INVALID_ARGS
             endif
            else
             if(.not.present(group_spec)) then
              call this%header%add_shape(errc,dim_extent)
             else
              errc=TEREC_INVALID_ARGS
             endif
            endif
           endif
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensRcrsvSetShape
!---------------------------------------------------------------------
        subroutine TensRcrsvSetLayout(this,layout_kind,data_type,ierr)
!Sets the tensor body storage layout.
         implicit none
         class(tens_rcrsv_t), intent(inout):: this   !inout: tensor
         integer(INTD), intent(in):: layout_kind     !in: tensor body storage layout kind
         integer(INTD), intent(in):: data_type       !in: tensor body data type: {R4,R8,C4,C8}
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,unres
         logical:: shpd,layd

         if(this%is_set(errc,shaped=shpd,unresolved=unres,layed=layd)) then
          if(errc.eq.TEREC_SUCCESS) then
           if(shpd.and.unres.eq.0.and.(.not.layd)) then
            call this%body%set_layout(layout_kind,data_type,errc)
           else
            errc=TEREC_INVALID_REQUEST
           endif
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensRcrsvSetLayout
!------------------------------------------------------------
        subroutine TensRcrsvSetLocation(this,data_descr,ierr)
!Sets the physical location of the tensor body data via a DDSS data descriptor.
         implicit none
         class(tens_rcrsv_t), intent(inout):: this   !inout: tensor
         class(DataDescr_t), intent(in):: data_descr !in: DDSS data descriptor for tensor body data
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,unres
         logical:: shpd,layd,locd

         if(this%is_set(errc,shaped=shpd,unresolved=unres,layed=layd,located=locd)) then
          if(errc.eq.TEREC_SUCCESS) then
           if(shpd.and.unres.eq.0.and.layd.and.(.not.locd)) then
            call this%body%set_location(data_descr,errc)
           else
            errc=TEREC_INVALID_REQUEST
           endif
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensRcrsvSetLocation
!--------------------------------------------------------------
        function TensRcrsvGetHeader(this,ierr) result(header_p)
!Returns a pointer to the tensor header.
         implicit none
         class(tens_header_t), pointer:: header_p       !out: pointer to the tensor header
         class(tens_rcrsv_t), intent(in), target:: this !in: tensor
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD):: errc

         header_p=>NULL()
         if(this%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) header_p=>this%header
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensRcrsvGetHeader
!----------------------------------------------------------
        function TensRcrsvGetBody(this,ierr) result(body_p)
!Returns a pointer to the tensor body.
         implicit none
         class(tens_body_t), pointer:: body_p           !out: pointer to the tensor body
         class(tens_rcrsv_t), intent(in), target:: this !in: tensor
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD):: errc

         body_p=>NULL()
         if(this%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) body_p=>this%body
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensRcrsvGetBody
!--------------------------------------------------------------------------------
        subroutine TensRcrsvSplit(this,split_dims,subtensors,ierr,num_subtensors)
!Splits the given tensor into subtensors and appends those to a list of subtensors (by their headers).
!If the input parental tensor is shaped with concrete dimension extents, the children subtensors
!will not carry concrete dimension extents, but deferred dimension extents instead. However,
!they will obey the parental dimension restriction grouping with the following assumption:
!ASSUMPTION: Restricted tensor dimensions belonging to the same group are supposed to depend on
!each other from right to left, e.g., in {T(a,b,c,d,e):[a<c<e],[b<d]}, a dependent index
!on the right always depends on some of the previous indices on the left. The first
!restricted index on the left in each group is actually independent (e.g., "a" and "b" above).
         implicit none
         class(tens_rcrsv_t), intent(in):: this                !in: parental tensor (either shaped or unshaped)
         integer(INTD), intent(in):: split_dims(1:)            !in: tensor dimensions to be split (at least one)
         type(list_bi_t), intent(inout):: subtensors           !out: list of subtensors specified by their tensor headers
         integer(INTD), intent(out), optional:: ierr           !out: error code
         integer(INTD), intent(out), optional:: num_subtensors !out: number of subtensors generated
         integer(INTD):: i,j,n,tnl,nsb,nd,sd,nsubt,ngroups,errc
         character(128):: tens_name
         integer(INTL):: sidx(1:MAX_TENSOR_RANK)       !parental subspace multi-index
         integer(INTL), allocatable:: sbuf(:)          !temporary buffer for holding children subspace id's
         integer(INTD):: midx(1:MAX_TENSOR_RANK)       !subspace iterator register
         integer(INTD):: firo(0:MAX_TENSOR_RANK)       !base offset in sbuf() for each tensor dimension
         integer(INTD):: swid(0:MAX_TENSOR_RANK)       !number of children subspaces for each dimension in sbuf()
         integer(INTD):: deps(1:MAX_TENSOR_RANK)       !dimension dependencies (grouping)
         integer(INTD):: depk(1:MAX_TENSOR_RANK)       !dimension dependency kinds
         integer(INTD):: tmpd(1:MAX_TENSOR_RANK)       !temporary (reduced subtensor dimension dependencies)
         integer(INTD):: dim_group(1:MAX_TENSOR_RANK)  !dimension groups
         integer(INTD):: group_spec(1:MAX_TENSOR_RANK) !dimension group restriction specs
         class(h_space_t), pointer:: hsp
         type(list_iter_t):: lit
         type(tens_header_t):: thead
         logical:: shpd

         nsubt=0 !number of generated subtensors
         if(this%is_set(errc,shaped=shpd)) then !shaped tensors with deferred dimension extents are expected
          if(errc.eq.TEREC_SUCCESS) then
           call this%header%get_name(tens_name,tnl,errc)
           if(errc.eq.TEREC_SUCCESS) then
            hsp=>NULL()
            call this%header%get_spec(sidx,nd,errc,hsp) !nd: total number of tensor dimensions; sidx(1:nd): subspace id's
            if(errc.eq.TEREC_SUCCESS) then
             errc=lit%init(subtensors)
             if(errc.eq.GFC_SUCCESS) then
              if(nd.gt.0.and.associated(hsp)) then !true tensor on hierarchical vector space
               sd=size(split_dims) !sd: number of tensor dimensions to split
               if(sd.gt.0.and.sd.le.nd) then !true splitting
                call extract_subspaces_to_sbuf(errc) !sbuf(1:nsb),firo(1:nd),swid(1:nd)
                if(errc.eq.TEREC_SUCCESS) then
                 call setup_index_dependencies(errc)
                 if(errc.eq.TEREC_SUCCESS) then
                  midx(1:nd)=-1; n=2 !n=2 is a special setting to start iterating midx(:)
                  do
                   call get_next_midx(n,errc); if(errc.ne.TEREC_SUCCESS.or.n.le.0) exit
                   call construct_subtensor_header(thead,errc); if(errc.ne.TEREC_SUCCESS) exit
                   errc=lit%append(thead); if(errc.eq.GFC_SUCCESS) then; nsubt=nsubt+1; else; exit; endif
                  enddo
                  if(errc.ne.TEREC_SUCCESS) errc=TEREC_UNABLE_COMPLETE
                 endif
                endif
               elseif(sd.eq.0) then !no splitting, return the original header
                errc=lit%append(this%header)
                if(errc.eq.GFC_SUCCESS) then; nsubt=nsubt+1; else; errc=TEREC_UNABLE_COMPLETE; endif
               else
                errc=TEREC_INVALID_ARGS
               endif
              else
               if(nd.eq.0) then
                errc=TEREC_INVALID_REQUEST !scalars cannot be split further
               else
                errc=TEREC_ERROR !unable to retrieve the hierarhical vector space info
               endif
              endif
              i=lit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.TEREC_SUCCESS) errc=TEREC_ERROR
             else
              errc=TEREC_ERROR
             endif
            endif
           endif
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(allocated(sbuf)) deallocate(sbuf)
         if(present(num_subtensors)) num_subtensors=nsubt
         if(present(ierr)) ierr=errc
         return

         contains

          subroutine extract_subspaces_to_sbuf(jerr) !sets sbuf(1:nsb),firo(1:nd),swid(1:nd)
           implicit none
           integer(INTD), intent(out):: jerr
           integer(INTD):: jj,js,jd
           type(vec_tree_iter_t):: vt_it
           class(*), pointer:: jup
           class(subspace_t), pointer:: jssp

           jerr=vt_it%init(hsp%get_aggr_tree(nsb))
           if(jerr.eq.GFC_SUCCESS.and.nsb.eq.0) then
 !Count:
            nsb=0; firo(0:nd)=0; swid(0:nd)=1
            do jj=1,sd !loop over the dimensions to split
             jd=split_dims(jj)
             jerr=vt_it%move_to(sidx(jd)); if(jerr.ne.GFC_SUCCESS) exit
             js=vt_it%get_num_children(jerr); if(jerr.ne.GFC_SUCCESS) exit
             firo(jd)=nsb+1
             if(js.gt.0) then !splitting will occur
              swid(jd)=-js; nsb=nsb+js !place for children subspaces
             else !no children -> splitting is impossible
              swid(jd)=-1; nsb=nsb+1 !place for the parental subspace
             endif
            enddo
            if(jerr.eq.GFC_SUCCESS) then
 !Insert unsplit dimensions:
             jj=0
             do jd=1,nd
              if(firo(jd).gt.0) then
               firo(jd)=firo(jd)+jj
              else
               firo(jd)=firo(jd-1)+abs(swid(jd-1)); jj=jj+1
              endif
             enddo
             nsb=nsb+jj
 !Set up:
             allocate(sbuf(1:nsb),STAT=jerr)
             if(jerr.eq.0) then
              sloop: do jd=1,nd
               jj=firo(jd)
               if(swid(jd).lt.0) then !split dimension
                swid(jd)=-swid(jd); js=swid(jd)
                jerr=vt_it%move_to(sidx(jd)); if(jerr.ne.GFC_SUCCESS) exit sloop
                if(js.gt.0) then
                 jerr=vt_it%move_to_child()
                 do while(js.gt.0.and.jerr.eq.GFC_SUCCESS)
                  jup=>vt_it%get_value(jerr); if(jerr.ne.GFC_SUCCESS) exit
                  select type(jup); class is(subspace_t); jssp=>jup; end select
                  sbuf(jj)=jssp%get_id(jerr); if(jerr.ne.0) exit
                  jj=jj+1; js=js-1; if(js.gt.0) jerr=vt_it%move_to_sibling()
                 enddo
                 if(jerr.ne.GFC_SUCCESS) exit sloop
                else !no children subspaces
                 sbuf(jj)=sidx(jd)
                endif
               else !unsplit dimension
                sbuf(jj)=sidx(jd)
               endif
              enddo sloop
              if(jerr.ne.GFC_SUCCESS) jerr=TEREC_UNABLE_COMPLETE
             else
              jerr=TEREC_MEM_ALLOC_FAILED
             endif
            else
             jerr=TEREC_ERROR
            endif
            jj=vt_it%release(); if(jj.ne.GFC_SUCCESS.and.jerr.eq.TEREC_SUCCESS) jerr=TEREC_ERROR
           else
            jerr=TEREC_ERROR
           endif
           return
          end subroutine extract_subspaces_to_sbuf

          subroutine setup_index_dependencies(jerr) !sets ngroups,deps(1:nd),depk(1:nd)
           implicit none
           integer(INTD), intent(out):: jerr
           integer(INTD):: jj,ji,jn,jr

           jerr=TEREC_SUCCESS
           deps(1:nd)=0; depk(1:nd)=TEREC_IND_RESTR_NONE !no dependencies by default
           if(shpd) then
            ngroups=this%header%num_groups(jerr) !number of non-trivial index restriction groups
            if(jerr.eq.TEREC_SUCCESS) then
             if(ngroups.gt.0) then
              do jj=1,ngroups
               call this%header%get_group(jj,dim_group,jn,jerr,jr); if(jerr.ne.TEREC_SUCCESS) exit
               if(jn.gt.0) then
                depk(dim_group(1))=jr
                do ji=2,jn
                 deps(dim_group(ji))=dim_group(ji-1); depk(dim_group(ji))=jr
                enddo
               endif
              enddo
             endif
            endif
           else
            ngroups=0
           endif
           return
          end subroutine setup_index_dependencies

          subroutine get_next_midx(np,jerr) !iterates midx(1:nd)
           implicit none
           integer(INTD), intent(inout):: np
           integer(INTD), intent(out):: jerr

           jerr=TEREC_SUCCESS; np=np-1
           mloop: do while(np.gt.0.and.np.le.nd)
            if(midx(np).ge.0) then
             if(midx(np)+1.lt.swid(np)) then !not the last value
              midx(np)=midx(np)+1
              if(deps(np).gt.0) then !dependent dimension: acceptance check
               if(restr_move_rejected(np,jerr)) cycle mloop
               if(jerr.ne.TEREC_SUCCESS) exit mloop
              endif
              np=np+1
             else !the last value
              midx(np)=-1; np=np-1
             endif
            else
             midx(np)=0 !first value
             if(deps(np).gt.0) then !dependent dimension: acceptance check
              if(restr_move_rejected(np,jerr)) cycle mloop
              if(jerr.ne.TEREC_SUCCESS) exit mloop
             endif
             np=np+1
            endif
           enddo mloop
           return
          end subroutine get_next_midx

          function restr_move_rejected(np,jerr) result(rejected) !rejects non-conforming midx(:) values
           implicit none
           logical:: rejected                !out: rejection decision
           integer(INTD), intent(in):: np    !in: specific dimension that has just changed its value
           integer(INTD), intent(out):: jerr !out: error code
           integer(INTD):: jd,jcmp
           integer(INTL):: js1,js2

           rejected=.FALSE.; jerr=TEREC_SUCCESS
           jd=deps(np) !dimension which dimension np depends on (must be on the left of np)
           if(jd.lt.np) then
            if(jd.gt.0) then
             js1=sbuf(firo(jd)+midx(jd)) !left subspace id
             js2=sbuf(firo(np)+midx(np)) !right subspace id (dependent)
             jcmp=hsp%compare_subranges(js1,js2,jerr)
             if(jerr.eq.0) then
              if(depk(np).eq.TEREC_IND_RESTR_LT.or.depk(np).eq.TEREC_IND_RESTR_LE) then
               if(jcmp.eq.CMP_GT) rejected=.TRUE.
              elseif(depk(np).eq.TEREC_IND_RESTR_GT.or.depk(np).eq.TEREC_IND_RESTR_GE) then
               if(jcmp.eq.CMP_LT) rejected=.TRUE.
              else
               jerr=TEREC_UNABLE_COMPLETE
              endif
             else
              jerr=TEREC_ERROR
             endif
            endif
           else
            jerr=TEREC_ERROR
           endif
           return
          end function restr_move_rejected

          subroutine construct_subtensor_header(thd,jerr) !constructs the next subtensor header, builds tmpd(1:nd), dim_group(1:nd), group_spec(1:) prior to that
           implicit none
           type(tens_header_t), intent(inout):: thd !out: subtensor header
           integer(INTD), intent(out):: jerr        !out: error code
           integer(INTD):: jd,jp,jgr,jcmp
           integer(INTL):: js1,js2

           jerr=TEREC_SUCCESS
 !Set children subspace multi-index (reuse sidx(:)):
           do jd=1,nd; sidx(jd)=sbuf(firo(jd)+midx(jd)); enddo
 !Set index restriction groups, if needed:
           if(ngroups.gt.0) then
  !Remove unneeded dimension dependencies:
            do jd=1,nd
             jp=deps(jd); tmpd(jd)=jp
             if(jp.gt.0) then !dependent dimension in the parental tensor: jd->jp
              if(jp.ge.jd) then; jerr=TEREC_ERROR; exit; endif !strict right-on-the-left dependency order
              js1=sbuf(firo(jp)+midx(jp)) !left subspace id
              js2=sbuf(firo(jd)+midx(jd)) !right subspace id
              jcmp=hsp%compare_subranges(js1,js2,jerr); if(jerr.ne.0) then; jerr=TEREC_ERROR; exit; endif
              if(depk(jd).eq.TEREC_IND_RESTR_LT.or.depk(jd).eq.TEREC_IND_RESTR_LE) then
               if(jcmp.eq.CMP_LT) tmpd(jd)=0 !dependency is automatically satisfied on (disjoint) children subspaces
              elseif(depk(jd).eq.TEREC_IND_RESTR_GT.or.depk(jd).eq.TEREC_IND_RESTR_GE) then
               if(jcmp.eq.CMP_GT) tmpd(jd)=0 !dependency is automatically satisfied on (disjoint) children subspaces
              else
               jerr=TEREC_UNABLE_COMPLETE; exit
              endif
             endif
            enddo
  !Mark first (independent) dimension in each restriction group:
            if(jerr.eq.TEREC_SUCCESS) then
             do jd=1,nd
              jp=tmpd(jd)
              if(jp.gt.0) then
               if(tmpd(jp).eq.0) tmpd(jp)=-1
              endif
             enddo
  !Create dim_group(:) and group_spec(:):
             jgr=0 !number of index restriction groups in the current subtensor
             do jd=1,nd
              jp=tmpd(jd)
              if(jp.gt.0) then
               dim_group(jd)=dim_group(jp)
              elseif(jp.lt.0) then !the 1st dimension of a new restriction group
               jgr=jgr+1; dim_group(jd)=jgr; group_spec(jgr)=depk(jd)
              else
               dim_group(jd)=0
              endif
             enddo
  !Construct a new subtensor with index restrictions:
             call thd%tens_header_ctor(jerr,tens_name(1:tnl),sidx(1:nd),hsp,(/(0_INTL,jd=1,nd)/),dim_group(1:nd),group_spec(1:jgr))
            endif
           else
  !Construct a new subtensor without index restrictions:
            call thd%tens_header_ctor(jerr,tens_name(1:tnl),sidx(1:nd),hsp,(/(0_INTL,jd=1,nd)/))
           endif
           return
          end subroutine construct_subtensor_header

        end subroutine TensRcrsvSplit
!---------------------------------------
        subroutine tens_rcrsv_dtor(this)
         implicit none
         type(tens_rcrsv_t):: this

         return
        end subroutine tens_rcrsv_dtor

       end module tensor_recursive
!==================================
       module tensor_recursive_test
        use tensor_algebra
        use gfc_base
        use gfc_list
        use subspaces
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
         write(*,'("Testing class tens_rcrsv_t ... ")',ADVANCE='NO')
         call test_tens_rcrsv(ierr)
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
         type(h_space_t):: h_space

         ierr=0
         call tsigna%tens_signature_ctor(ierr,(/3_INTL,4_INTL,2_INTL/),'Tensor',h_space)
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
         call tshape%tens_shape_ctor(ierr,dims(1:n),grps(1:n),grp_spec(1:m))
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
         type(h_space_t):: h_space
         character(32):: tens_name

         ierr=0
         n=6; dims(1:n)=(/128_INTL,64_INTL,256_INTL,64_INTL,128_INTL,64_INTL/)
         m=2; grps(1:n)=(/1,2,0,2,1,2/); grp_spec(1:m)=(/TEREC_IND_RESTR_LT,TEREC_IND_RESTR_GE/)
         call thead%tens_header_ctor(ierr,'Tensor',(/1_INTL,2_INTL,3_INTL,2_INTL,1_INTL,2_INTL/),h_space)
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
         type(h_space_t):: h_space
         character(32):: tens_name
         type(tens_simple_part_t):: tpart

         ierr=0
         n=6; dims(1:n)=(/128_INTL,64_INTL,256_INTL,64_INTL,128_INTL,64_INTL/)
         m=2; grps(1:n)=(/1,2,0,2,1,2/); grp_spec(1:m)=(/TEREC_IND_RESTR_LT,TEREC_IND_RESTR_GE/)
         call thead%tens_header_ctor(ierr,'Tensor',(/1_INTL,2_INTL,3_INTL,2_INTL,1_INTL,2_INTL/),h_space)
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
               call tpart%tens_simple_part_ctor(thead,TEREC_LAY_FDIMS,4096_INTL,ierr)
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
!---------------------------------------
        subroutine test_tens_rcrsv(ierr)
         implicit none
         integer(INTD), intent(out):: ierr
         integer(INTD), parameter:: tens_rank=4          !tensor rank
         integer(INTL), parameter:: TEST_SPACE_DIM=33    !vector space dimension
         type(spher_symmetry_t):: symm(1:TEST_SPACE_DIM) !symmetry of each basis vector
         type(subspace_basis_t):: full_basis             !full vector space basis
         type(h_space_t):: hspace                        !hierarchical representation of the vector space
         integer(INTL):: spcx(1:MAX_TENSOR_RANK),dims(1:MAX_TENSOR_RANK),space_id,max_res
         integer(INTD):: dimg(1:MAX_TENSOR_RANK),grps(1:MAX_TENSOR_RANK),num_subtensors
         type(tens_header_t), pointer:: thp
         type(tens_rcrsv_t):: tensor
         type(list_bi_t):: subtensors
         type(list_iter_t):: lit
         class(subspace_t), pointer:: ssp
         class(*), pointer:: up

 !Build a hierarchical representation for a test vector space:
         call register_test_space(ierr); if(ierr.ne.TEREC_SUCCESS) then; ierr=1; return; endif
 !Create the full tensor (over the full space):
  !Get full space id and its max resolution:
         space_id=hspace%get_common_subspace(0_INTL,TEST_SPACE_DIM-1_INTL,ierr); if(ierr.ne.0) then; ierr=2; return; endif
         ssp=>hspace%get_subspace(space_id,ierr); if(ierr.ne.0) then; ierr=3; return; endif
         if(.not.associated(ssp)) then; ierr=4; return; endif
         max_res=ssp%get_max_resolution(ierr); if(ierr.ne.0) then; ierr=5; return; endif
         !write(*,*) 'Space ID = ',space_id,': Max resolution = ',max_res !debug
  !Create a tensor over the full space:
         spcx(1:tens_rank)=space_id
         dims(1:tens_rank)=max_res
         dimg(1:tens_rank)=(/1,1,2,2/); grps(1:2)=(/TEREC_IND_RESTR_LT,TEREC_IND_RESTR_GT/)
         call tensor%tens_rcrsv_ctor('T2',spcx(1:tens_rank),hspace,ierr,dims(1:tens_rank),dimg(1:tens_rank),grps(1:2))
         if(ierr.ne.0) then; ierr=6; return; endif
 !Split the tensor into subtensors:
         call tensor%split((/1,2,3,4/),subtensors,ierr,num_subtensors); if(ierr.ne.0) then; ierr=7; return; endif
         !write(*,*) 'Number of subtensors generated = ',num_subtensors !debug
         ierr=lit%init(subtensors); if(ierr.ne.0) then; ierr=8; return; endif
         !ierr=lit%scanp(action_f=print_tens_header_f); if(ierr.eq.GFC_IT_DONE) ierr=0 !debug
         ierr=lit%delete_all(); if(ierr.ne.0) then; ierr=9; return; endif
         ierr=lit%release(); if(ierr.ne.0) then; ierr=10; return; endif

         return

         contains

          subroutine register_test_space(jerr)
           implicit none
           integer(INTD), intent(out):: jerr
           integer(INTL):: jj

           jerr=0
           call full_basis%subspace_basis_ctor(TEST_SPACE_DIM,jerr); if(jerr.ne.0) return
           do jj=1,TEST_SPACE_DIM
            call symm(jj)%spher_symmetry_ctor(int((jj-1)/5,INTD),0,jerr); if(jerr.ne.0) return
            call full_basis%set_basis_func(jj,BASIS_ABSTRACT,jerr,symm=symm(jj)); if(jerr.ne.0) return
           enddo
           call full_basis%finalize(jerr); if(jerr.ne.0) return
           call hspace%h_space_ctor(full_basis,jerr); if(jerr.ne.0) return
           !call hspace%print_it() !debug
           return
          end subroutine register_test_space

        end subroutine test_tens_rcrsv

       end module tensor_recursive_test
