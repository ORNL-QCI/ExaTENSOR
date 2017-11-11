!ExaTENSOR: Recursive (hierarchical) tensors
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/11/11

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
!NOTES:
! # Tensor definition stages:
!   a) Null: No tensor;
!   b) Defined: Tensor signature is defined;
!   c) Resolved: Tensor shape is fully defined (all dimensions resolved);
!   d) Laid-Out: Tensor layout is defined;
!   e) Mapped: Tensor body is physically mapped.
! # Tensor body value definition stages:
!   a) TEREC_BODY_UNDEF: Tensor body value is undefined;
!   b) TEREC_BODY_DEF: Tensor body value is defined and can be used;
!   c) TEREC_BODY_UPDATE: Tensor body value is currently being updated and cannot be used.
        use tensor_algebra !includes dil_basic
        use stsubs
        use timers
        use combinatoric
        use gfc_base
        use gfc_list
        use gfc_vector
        use gfc_vec_tree
        use gfc_dictionary
        use multords, only: multord_i8e
        use subspaces
        use pack_prim
        use distributed, only: DataDescr_t
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
 !Tensor contraction generator:
        integer(INTD), parameter, public:: TEREC_TCG_BUF_SIZE=1048576 !max number of subtensors per parental tensor split and potential symmetry unfolding
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
 !Tensor body value state:
        integer(INTD), parameter, public:: TEREC_BODY_UNDEF=0     !tensor body value is undefined
        integer(INTD), parameter, public:: TEREC_BODY_DEF=1       !tensor body value is defined and can be used
        integer(INTD), parameter, public:: TEREC_BODY_UPDATE=2    !tensor body value is currently being updated
!TYPES:
 !Register of hierarchical spaces:
        type, private:: hspace_register_t
         logical, private:: initialized=.FALSE.        !initialization status
         type(dictionary_t), private:: name2id         !symbolic name -> id map
         type(vector_t), private:: hspaces             !vector of h_space_t objects
         type(dictionary_iter_t), private:: name2id_it !name2id iterator
         type(vector_iter_t), private:: hspaces_it     !hspaces iterator
         contains
          procedure, private:: init=>HspaceRegisterInit                   !initializes the register
          procedure, public:: register_space=>HspaceRegisterRegisterSpace !registers a new hierarchical vector space
          procedure, public:: get_space_id=>HspaceRegisterGetSpaceId      !returns the registered id of the space by its name
          procedure, private:: HspaceRegisterGetSpaceByName               !returns a pointer to the hierarchical space by its name
          procedure, private:: HspaceRegisterGetSpaceById                 !returns a pointer to the hierarchical space by its id
          generic, public:: get_space=>HspaceRegisterGetSpaceByName,HspaceRegisterGetSpaceById !returns a pointer to the hierarchical space
          final:: hspace_register_dtor                                    !dtor
        end type hspace_register_t
 !Registered hierarchical space:
        type, public:: hspace_reg_t
         integer(INTD), private:: space_id=-1                  !registered space id: [0..max]
         class(h_space_t), pointer, private:: hspace_p=>NULL() !non-owning pointer to a defined (persistent) hierarchical vector space object
         contains
          procedure, private:: HspaceRegCtor                    !ctor
          procedure, private:: HspaceRegCtorUnpack              !ctor by unpacking
          generic, public:: hspace_reg_ctor=>HspaceRegCtor,HspaceRegCtorUnpack
          procedure, public:: pack=>HspaceRegPack               !packs the object into a packet
          procedure, public:: is_set=>HspaceRegIsSet            !returns TRUE if the object is set
          procedure, public:: get_space_id=>HspaceRegGetSpaceId !returns the registered space id: [0..max]
          procedure, public:: get_space=>HspaceRegGetSpace      !returns a pointer to the hierarchical space definition
          final:: hspace_reg_dtor                               !dtor
        end type hspace_reg_t
 !Tensor signature (unique tensor identifier):
        type, public:: tens_signature_t
         character(:), allocatable, private:: char_name       !character tensor name (alphanumeric_)
         integer(INTD), private:: num_dims=-1                 !number of tensor dimensions (aka tensor order in math or tensor rank in physics)
         integer(INTL), allocatable, private:: space_idx(:)   !subspace id for each tensor dimension
         type(hspace_reg_t), allocatable, private:: hspace(:) !hierarchical vector space id for each tensor dimension (optional)
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
          procedure, public:: get_groups=>TensShapeGetGroups !returns dimension groups and group restrictions
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
          procedure, public:: is_valid=>TensHeaderIsValid           !returns .TRUE. if the set tensor header is valid
          procedure, public:: get_name=>TensHeaderGetName           !returns the alphanumeric_ tensor name
          procedure, public:: get_rank=>TensHeaderGetRank           !returns the rank of the tensor (number of dimensions)
          procedure, public:: get_spec=>TensHeaderGetSpec           !returns the tensor subspace multi-index (specification)
          procedure, public:: get_dims=>TensHeaderGetDims           !returns tensor dimension extents
          procedure, public:: num_groups=>TensHeaderNumGroups       !returns the total number of non-trivial index groups defined in the tensor shape
          procedure, public:: get_dim_group=>TensHeaderGetDimGroup  !returns the restriction group for a specific tensor dimension (0: no restrictions)
          procedure, public:: get_groups=>TensHeaderGetGroups       !returns dimension groups and group restrictions
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
          procedure, public:: unpack_base=>TensLayoutUnpackBase                     !unpacks the object from a packet
          procedure, public:: pack_base=>TensLayoutPackBase                         !packs the object into a packet
          procedure, public:: get_data_type=>TensLayoutGetDataType                  !returns the data type of the stored tensor elements
          procedure, public:: get_layout_kind=>TensLayoutGetLayoutKind              !returns the tensor storage layout kind
          procedure, public:: get_body_ptr=>TensLayoutGetBodyPtr                    !returns a C pointer to the tensor body
          procedure, public:: get_body_size=>TensLayoutGetBodySize                  !returns the size of the stored tensor body in bytes
          procedure, public:: get_data_descr=>TensLayoutGetDataDescr                !returns a pointer to the DDSS data descriptor
          procedure(tens_layout_volume_i), deferred, public:: get_volume            !returns the physical volume of the tensor block (number of physically stored elements)
          procedure(tens_layout_map_i), deferred, public:: map                      !maps a specific element of the tensor block (layout specific)
          procedure(tens_layout_extract_i), deferred, public:: extract_simple_parts !creates a list of constituent simple (dense) parts of the tensor block
          procedure(tens_layout_update_i), deferred, private:: update_header        !updates the tensor header the layout is associated with
        end type tens_layout_t
 !Concrete storage layout "Fortran-dimension-led":
        type, extends(tens_layout_t), public:: tens_layout_fdims_t
         class(tens_header_t), pointer, private:: header=>NULL() !non-owning pointer to the defining tensor header
         contains
          procedure, private:: TensLayoutFdimsCtor                          !ctor
          procedure, private:: TensLayoutFdimsCtorUnpack                    !ctor by unpacking
          generic, public:: tens_layout_fdims_ctor=>TensLayoutFdimsCtor,TensLayoutFdimsCtorUnpack
          procedure, public:: pack=>TensLayoutFdimsPack                     !packs the object into a packet
          procedure, public:: get_volume=>TensLayoutFdimsGetVolume          !returns the physical tensor volume (number of elements stored)
          procedure, public:: map=>TensLayoutFdimsMap                       !addresses a specific tensor element
          procedure, public:: extract_simple_parts=>TensLayoutFdimsExtract  !extracts simpe dense tensor parts (bricks) from the tensor block
          procedure, private:: update_header=>TensLayoutFdimsUpdateHeader   !updates the tensor header the layout is associated with
          final:: tens_layout_fdims_dtor
        end type tens_layout_fdims_t
 !Tensor body:
        type, public:: tens_body_t
         integer(INTD), private:: state=TEREC_BODY_UNDEF      !tensor body value state: {TEREC_BODY_UNDEF,TEREC_BODY_DEF,TEREC_BODY_UPDATE}
         integer(INTD), private:: num_subtensors=0            !number of subtensors in the list
         type(list_bi_t), private:: subtensors                !list of constituent tensors in terms of tensor headers
         class(tens_layout_t), allocatable, private:: layout  !tensor block storage layout (if physically stored as a whole)
         contains
          procedure, private:: TensBodyCtorBase                     !basic ctor
          procedure, private:: TensBodyCtorUnpack                   !ctor by unpacking
          generic, public:: tens_body_ctor=>TensBodyCtorBase,TensBodyCtorUnpack
          procedure, public:: pack=>TensBodyPack                    !packs the object into a packet
          procedure, public:: is_set=>TensBodyIsSet                 !returns TRUE if the tensor body is set (plus additional info)
          procedure, public:: add_subtensor=>TensBodyAddSubtensor   !registers a constituent subtensor by providing its tensor header
          procedure, public:: set_layout=>TensBodySetLayout         !sets the tensor body storage layout if physically stored as a whole
          procedure, public:: set_location=>TensBodySetLocation     !sets the tensor body data location if physically stored as a whole (via a DDSS data descriptor)
          procedure, public:: reset_state=>TensBodyResetState       !resets the tensor body value state
          procedure, public:: get_state=>TensBodyGetState           !returns the tensor body value state
          procedure, public:: get_layout=>TensBodyGetLayout         !returns a pointer to the tensor body storage layout
          procedure, public:: get_num_subtensors=>TensBodyGetNumSubtensors !returns the total number of constituent subtensors
          procedure, public:: get_subtensors=>TensBodyGetSubtensors !returns a pointer to the list of constituent subtensors (each subtensor is represented by a tensor header)
          procedure, public:: print_it=>TensBodyPrintIt             !prints the tensor body info
          procedure, private:: update_header=>TensBodyUpdateHeader  !updates the tensor header the body is associated with
          final:: tens_body_dtor                                    !dtor
        end type tens_body_t
 !Recursive tensor:
        type, public:: tens_rcrsv_t
         type(tens_header_t), private:: header !tensor header (signature + shape)
         type(tens_body_t), private:: body     !tensor body (recursive composition, data location and storage layout)
         contains
          procedure, private:: TensRcrsvCtorSigna                    !ctor by tensor signature and optionally tensor shape
          procedure, private:: TensRcrsvCtorHead                     !ctor by tensor header
          procedure, private:: TensRcrsvCtorClone                    !copy ctor
          procedure, private:: TensRcrsvCtorUnpack                   !ctor by unpacking
          generic, public:: tens_rcrsv_ctor=>TensRcrsvCtorSigna,TensRcrsvCtorHead,TensRcrsvCtorClone,TensRcrsvCtorUnpack
          procedure, public:: pack=>TensRcrsvPack                    !packs the object into a packet
          procedure, public:: is_set=>TensRcrsvIsSet                 !returns TRUE if the tensor is set (signature defined) plus other info
          procedure, public:: get_name=>TensRcrsvGetName             !returns the alphanumeric_ tensor name
          procedure, public:: get_rank=>TensRcrsvGetRank             !returns the rank of the tensor (number of dimensions)
          procedure, public:: get_spec=>TensRcrsvGetSpec             !returns the tensor subspace multi-index (specification)
          procedure, public:: get_dims=>TensRcrsvGetDims             !returns tensor dimension extents
          procedure, public:: add_subtensor=>TensRcrsvAddSubtensor   !registers a constituent subtensor by providing its tensor header
          procedure, public:: add_subtensors=>TensRcrsvAddSubtensors !registers constituent subtensors by providing a list of their tensor headers
          procedure, public:: set_shape=>TensRcrsvSetShape           !sets the tensor shape (if it has not been set yet)
          procedure, public:: set_layout=>TensRcrsvSetLayout         !sets the tensor body storage layout
          procedure, public:: set_location=>TensRcrsvSetLocation     !sets the physical location of the tensor body data
          procedure, public:: update=>TensRcrsvUpdate                !updates the tensor information (new resolution -> new layout -> new body)
          procedure, public:: reset_state=>TensRcrsvResetState       !resets the tensor body value state
          procedure, public:: get_state=>TensRcrsvGetState           !returns the tensor body value state
          procedure, public:: get_header=>TensRcrsvGetHeader         !returns a pointer to the tensor header
          procedure, public:: get_body=>TensRcrsvGetBody             !returns a pointer to the tensor body
          procedure, public:: get_descriptor=>TensRcrsvGetDescriptor !returns a tensor descriptor uniquely characterizing tensor signature, shape, layout, and location
          procedure, private:: TensRcrsvSplitList                    !splits the tensor into subtensors (a list of subtensors by their headers)
          procedure, private:: TensRcrsvSplitVector                  !splits the tensor into subtensors (a vector of subtensors by their headers)
          generic, public:: split=>TensRcrsvSplitList,TensRcrsvSplitVector
          procedure, public:: print_it=>TensRcrsvPrintIt             !prints the tensor info
#ifdef NO_GNU
          final:: tens_rcrsv_dtor                                    !dtor `GCC/5.4.0 bug
#endif
        end type tens_rcrsv_t
 !Tensor descriptor:
        type, public:: tens_descr_t
         character(:), allocatable, private:: char_name !symbolic tensor name
         integer(INTL), allocatable, private:: info(:)  !combined: space_idx,subspace_idx,dim_extent,dim_group,dim_group_restr,layout_kind,volume,location
         integer(INTD), private:: rank=-1               !tensor rank (number of dimensions)
         contains
          procedure, public:: print_it=>TensDescrPrintIt !prints the object
          procedure, public:: compare=>TensDescrCompare  !compares with another instance
          final:: tens_descr_dtor                        !dtor
        end type tens_descr_t
 !Tensor argument (reference to a recursive tensor):`Must not be cloned if .alloc=TRUE
        type, private:: tens_argument_t
         class(tens_rcrsv_t), pointer, private:: tens_p=>NULL() !pointer to a persistent tensor
         logical, private:: alloc=.FALSE.                       !TRUE if the tensor argument was allocated, FALSE if associated
         contains
          procedure, private:: set_tensor=>TensArgumentSetTensor           !sets up the tensor argument (ctor) by pointer association
          procedure, private:: allocate_tensor=>TensArgumentAllocateTensor !allocates an empty tensor for a subsequent definition
          procedure, private:: is_set=>TensArgumentIsSet                   !returns TRUE if the tensor argument is set, plus additional info
          procedure, private:: free_tensor=>TensArgumentFreeTensor         !frees the tensor (either by deallocation or by dissociation only)
          final:: tens_argument_dtor                                       !dtor
        end type tens_argument_t
 !Tensor operation (abstract):
        type, abstract, public:: tens_operation_t
         integer(INTD), private:: num_args=0                                !number of tensor arguments
         type(tens_argument_t), private:: tens_arg(0:MAX_TENSOR_OPERANDS-1) !tensor arguments: [0..num_args-1], argument 0 is always the destination tensor
         contains
          procedure(tens_operation_query_i), deferred, public:: is_set    !returns TRUE if the tensor operation is fully set
          procedure(tens_operation_query_i), deferred, public:: args_full !returns TRUE if all required tensor arguments are set
          procedure, public:: clean=>TensOperationClean                   !cleans the tensor operation to an empty state
          procedure, public:: set_argument=>TensOperationSetArgument      !sets up the next tensor argument by pointer association
          procedure, public:: reset_argument=>TensOperationResetArgument  !resets an already set argument
          procedure, public:: get_num_args=>TensOperationGetNumArgs       !returns the number of set arguments
          procedure, public:: get_argument=>TensOperationGetArgument      !returns a pointer to the specific tensor argument (tens_rcrsv_t)
          procedure, private:: allocate_argument=>TensOperationAllocateArgument !allocates the next tensor argument for a subsequent setup
          procedure, private:: free_arguments=>TensOperationFreeArguments       !deallocates/dissociates all arguments
        end type tens_operation_t
 !Tensor dimension permutation:
        type, public:: permutation_t
         integer(INTD), private:: length=0
         integer(INTD), allocatable, private:: prm(:) !prm(1:length) is the permutation itself, prm(0) is the current sign of the permutation
         contains
          procedure, public:: reset=>PermutationReset          !resets (constructs) the permutation
          procedure, public:: get_access=>PermutationGetAccess !returns a pointer to the permutation array, with or without the sign
          procedure, public:: get_sign=>PermutationGetSign     !returns the current sign of the permutation
          procedure, public:: set_sign=>PermutationSetSign     !sets a new sign to the permutation
          procedure, public:: invert=>PermutationInvert        !inverts the permutation (in-place)
          final:: permutation_dtor                             !dtor
        end type permutation_t
 !Extended digital tensor contraction pattern:
        type, public:: contr_ptrn_ext_t
         integer(INTD), private:: ddim=-1                     !destination tensor rank
         integer(INTD), private:: ldim=-1                     !left tensor rank
         integer(INTD), private:: rdim=-1                     !right tensor rank
         logical, private:: ind_restr_set=.FALSE.             !TRUE if index permutational restrictions were set
         integer(INTD), private:: dind_pos(1:MAX_TENSOR_RANK) !corresponding positions for indices of the destination tensor
         integer(INTD), private:: lind_pos(1:MAX_TENSOR_RANK) !corresponding positions for indices of the left tensor
         integer(INTD), private:: rind_pos(1:MAX_TENSOR_RANK) !corresponding positions for indices of the right tensor
         integer(INTD), private:: dind_res(1:MAX_TENSOR_RANK) !permutational dependencies for the destination tensor indices (y=dind_res(x): position x depends on position y on the left, -1 first position)
         integer(INTD), private:: lind_res(1:MAX_TENSOR_RANK) !permutational dependencies for the left tensor indices (y=lind_res(x): position x depends on position y on the left, -1 first position)
         integer(INTD), private:: rind_res(1:MAX_TENSOR_RANK) !permutational dependencies for the right tensor indices (y=rind_res(x): position x depends on position y on the left, -1 first position)
         contains
          procedure, public:: set_index_corr=>ContrPtrnExtSetIndexCorr !sets index correspondence pattern (contraction pattern): basic ctor
          procedure, public:: set_store_symm=>ContrPtrnExtSetStoreSymm !sets index permutational symmetry restrictions due to tensor storage: post-ctor
          procedure, public:: set_operl_symm=>ContrPtrnExtSetOperlSymm !sets index permutational symmetry restrictions due to tensor operation: post-ctor
          procedure, public:: break_dim_symm=>ContrPtrnExtBreakDimSymm !breaks the dimension symmetry for a specific tensor dimension
          procedure, public:: unpack=>ContrPtrnExtUnpack               !unpacks the object from a packet
          procedure, public:: is_set=>ContrPtrnExtIsSet                !returns TRUE of the tensor contraction pattern is set
          procedure, public:: pack=>ContrPtrnExtPack                   !packs the object into a packet
          procedure, public:: get_contr_ptrn=>ContrPtrnExtGetContrPtrn !returns the classical (basic) digital contraction pattern used by TAL-SH for example
          procedure, public:: get_dim_symmetry=>ContrPtrnExtGetDimSymmetry !returns the symmetric restrictions for a specific tensor argument
          procedure, public:: print_it=>ContrPtrnExtPrintIt            !prints the extended tensor contraction
        end type contr_ptrn_ext_t
 !Tensor contraction:
        type, extends(tens_operation_t), public:: tens_contraction_t
         type(contr_ptrn_ext_t), private:: contr_ptrn                     !extended tensor contraction pattern
         complex(8), private:: alpha=(1d0,0d0)                            !alpha prefactor
         contains
          procedure, private:: TensContractionAssign                       !copy assignment
          generic, public:: assignment(=)=>TensContractionAssign
          procedure, public:: is_set=>TensContractionIsSet                 !returns TRUE if the tensor contraction is fully set
          procedure, public:: args_full=>TensContractionArgsFull           !returns TRUE if all tensor contraction arguments have been set
          procedure, public:: set_contr_ptrn=>TensContractionSetContrPtrn  !sets the tensor contraction pattern (all tensor arguments must have been set already)
          procedure, public:: set_operl_symm=>TensContractionSetOperlSymm  !sets index permutational symmetry restrictions due to tensor operation (both contraction pattern and arguments must have been set already)
          procedure, public:: unpack=>TensContractionUnpack                !unpacks the object from a packet
          procedure, public:: pack=>TensContractionPack                    !packs the object into a packet
          procedure, public:: get_prefactor=>TensContractionGetPrefactor   !returns the scalar prefactor
          procedure, public:: get_ext_contr_ptrn=>TensContractionGetExtContrPtrn !returns a pointer to the extended tensor contraction pattern
          procedure, public:: get_contr_ptrn=>TensContractionGetContrPtrn  !returns the classical (basic) digital contraction pattern used by TAL-SH for example
          procedure, private:: import_replace=>TensContractionImportReplace!creates a new tensor contraction by replacing tensor arguments in an existing tensor contraction (plus symmetry adjustment)
          procedure, public:: split=>TensContractionSplit                  !splits the tensor contraction into a list of subtensor contractions based on the pre-existing lists of argument subtensors
          procedure, public:: print_it=>TensContractionPrintIt             !prints the tensor contraction info
        end type tens_contraction_t
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
  !tens_layout_t: .update_header():
         subroutine tens_layout_update_i(this,header,ierr)
          import:: INTD,tens_header_t,tens_layout_t
          implicit none
          class(tens_layout_t), intent(inout):: this        !inout: tensor block storage layout
          class(tens_header_t), intent(in), target:: header !in: tensor header the layout is associated with
          integer(INTD), intent(out), optional:: ierr       !out: error code
         end subroutine tens_layout_update_i
  !tens_operation_t: .is_set():
         function tens_operation_query_i(this,ierr) result(ans)
          import:: INTD,tens_operation_t
          implicit none
          logical:: ans                               !out: answer
          class(tens_operation_t), intent(in):: this  !in: tensor operation
          integer(INTD), intent(out), optional:: ierr !out: error code
         end function tens_operation_query_i
  !tens_rcrsv_t: split():
         function tens_rcrsv_split_i(this,subtensors,num_subtensors) result(ierr)
          import:: INTD,tens_rcrsv_t,vector_t
          implicit none
          integer(INTD):: ierr                        !out: error code
          class(tens_rcrsv_t), intent(in):: this      !in: parental tensor
          type(vector_t), intent(inout):: subtensors  !inout: vector of subtensors (in general, can be non-empty on entrance)
          integer(INTD), intent(out):: num_subtensors !out: number of subtensors generated from the parental tensor
         end function tens_rcrsv_split_i
        end interface
!VISIBILITY:
 !non-member:
        public valid_tensor_layout
        public cmp_integers
        public cmp_real
        public cmp_strings
        public cmp_tens_signatures
        public cmp_tens_headers
        public cmp_tens_descriptors
        public build_test_hspace
        public print_tens_header_f
        public print_tcg_buffer
 !hspace_register_t:
        private HspaceRegisterInit
        private HspaceRegisterRegisterSpace
        private HspaceRegisterGetSpaceId
        private HspaceRegisterGetSpaceByName
        private HspaceRegisterGetSpaceById
        public hspace_register_dtor
 !hspace_reg_t:
        private HspaceRegCtor
        private HspaceRegCtorUnpack
        private HspaceRegPack
        private HspaceRegIsSet
        private HspaceRegGetSpaceId
        private HspaceRegGetSpace
        public hspace_reg_dtor
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
        private TensShapeGetGroups
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
        private TensHeaderIsValid
        private TensHeaderGetName
        private TensHeaderGetRank
        private TensHeaderGetSpec
        private TensHeaderGetDims
        private TensHeaderNumGroups
        private TensHeaderGetDimGroup
        private TensHeaderGetGroups
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
        private TensLayoutUnpackBase
        private TensLayoutPackBase
        private TensLayoutGetDataType
        private TensLayoutGetLayoutKind
        private TensLayoutGetBodyPtr
        private TensLayoutGetBodySize
        private TensLayoutGetDataDescr
 !tens_layout_fdims_t:
        private TensLayoutFdimsCtor
        private TensLayoutFdimsCtorUnpack
        private TensLayoutFdimsPack
        private TensLayoutFdimsGetVolume
        private TensLayoutFdimsMap
        private TensLayoutFdimsExtract
        private TensLayoutFdimsUpdateHeader
        public tens_layout_fdims_dtor
 !tens_body_t:
        private TensBodyCtorBase
        private TensBodyCtorUnpack
        private TensBodyPack
        private TensBodyIsSet
        private TensBodyAddSubtensor
        private TensBodySetLayout
        private TensBodySetLocation
        private TensBodyResetState
        private TensBodyGetState
        private TensBodyGetLayout
        private TensBodyGetNumSubtensors
        private TensBodyGetSubtensors
        private TensBodyPrintIt
        private TensBodyUpdateHeader
        public tens_body_dtor
 !tens_rcrsv_t:
        private TensRcrsvCtorSigna
        private TensRcrsvCtorHead
        private TensRcrsvCtorClone
        private TensRcrsvCtorUnpack
        private TensRcrsvPack
        private TensRcrsvIsSet
        private TensRcrsvGetName
        private TensRcrsvGetRank
        private TensRcrsvGetSpec
        private TensRcrsvGetDims
        private TensRcrsvAddSubtensor
        private TensRcrsvAddSubtensors
        private TensRcrsvSetShape
        private TensRcrsvSetLayout
        private TensRcrsvSetLocation
        private TensRcrsvUpdate
        private TensRcrsvResetState
        private TensRcrsvGetState
        private TensRcrsvGetHeader
        private TensRcrsvGetBody
        private TensRcrsvGetDescriptor
        private TensRcrsvSplitList
        private TensRcrsvSplitVector
        private TensRcrsvPrintIt
        public tens_rcrsv_dtor
 !tens_descr_t:
        private TensDescrPrintIt
        private TensDescrCompare
        public tens_descr_dtor
 !tens_argument_t:
        private TensArgumentSetTensor
        private TensArgumentAllocateTensor
        private TensArgumentIsSet
        private TensArgumentFreeTensor
        private tens_argument_dtor
 !tens_operation_t:
        private TensOperationClean
        private TensOperationSetArgument
        private TensOperationResetArgument
        private TensOperationGetNumArgs
        private TensOperationGetArgument
        private TensOperationAllocateArgument
        private TensOperationFreeArguments
 !permutation_t:
        private PermutationReset
        private PermutationGetAccess
        private PermutationGetSign
        private PermutationSetSign
        private PermutationInvert
        public permutation_dtor
 !contr_ptrn_ext_t:
        private ContrPtrnExtSetIndexCorr
        private ContrPtrnExtSetStoreSymm
        private ContrPtrnExtSetOperlSymm
        private ContrPtrnExtBreakDimSymm
        private ContrPtrnExtUnpack
        private ContrPtrnExtIsSet
        private ContrPtrnExtPack
        private ContrPtrnExtGetContrPtrn
        private ContrPtrnExtGetDimSymmetry
        private ContrPtrnExtPrintIt
 !tens_contraction_t:
        private TensContractionAssign
        private TensContractionIsSet
        private TensContractionArgsFull
        private TensContractionSetContrPtrn
        private TensContractionSetOperlSymm
        private TensContractionUnpack
        private TensContractionPack
        private TensContractionGetPrefactor
        private TensContractionGetExtContrPtrn
        private TensContractionGetContrPtrn
        private TensContractionImportReplace
        private TensContractionSplit
        private TensContractionPrintIt
!DATA:
 !Register of hierarchical vector spaces (only these spaces can be used in tensors):
        type(hspace_register_t), public:: hspace_register
 ![TESTING]: prerequisites for building a hierarchical space:
        integer(INTL), parameter, public:: HSPACE_DIM_=20                    !vector space dimension
        type(spher_symmetry_t), target, public:: hspace_symm_(1:HSPACE_DIM_) !symmetry of basis vectors
        type(subspace_basis_t), target, public:: hspace_basis_               !vector space basis
 !Tensor contraction generator: subtensor buffer:
        integer(INTL), allocatable, target, private:: tcg_ind_buf(:,:) !subtensor index buffer (private to each OpenMP thread)
        integer(INTL), allocatable, target, private:: tcg_num_buf(:)   !subtensor number buffer (private to each OpenMP thread)
!$OMP THREADPRIVATE(tcg_ind_buf,tcg_num_buf)

       contains
!IMPLEMENTATION:
![Non-member]===========================================
        function valid_tensor_layout(layout) result(res)
!Returns TRUE if the tensor layout is valid.
         implicit none
         logical:: res
         integer(INTD), intent(in):: layout

         res=(layout.ge.0.and.layout.lt.TEREC_NUM_LAYOUTS)
         return
        end function valid_tensor_layout
!---------------------------------------------------
        function cmp_integers(int1,int2) result(cmp)
!Generic comparator for integers or arbitrary kinds.
         implicit none
         integer(INTD):: cmp                 !out: result of comparison: {CMP_EQ,CMP_LT,CMP_GT,CMP_ER}
         class(*), intent(in), target:: int1 !in: integer 1
         class(*), intent(in), target:: int2 !in: integer 2
         integer(INTL):: i1,i2

         cmp=CMP_ER
         select type(int1)
         type is(integer(1))
          i1=int1; cmp=CMP_EQ
         type is(integer(2))
          i1=int1; cmp=CMP_EQ
         type is(integer(4))
          i1=int1; cmp=CMP_EQ
         type is(integer(8))
          i1=int1; cmp=CMP_EQ
         end select
         if(cmp.eq.CMP_EQ) then
          cmp=CMP_ER
          select type(int2)
          type is(integer(1))
           i2=int2; cmp=CMP_EQ
          type is(integer(2))
           i2=int2; cmp=CMP_EQ
          type is(integer(4))
           i2=int2; cmp=CMP_EQ
          type is(integer(8))
           i2=int2; cmp=CMP_EQ
          end select
          if(cmp.eq.CMP_EQ) then
           if(i1.lt.i2) then
            cmp=CMP_LT
           elseif(i1.gt.i2) then
            cmp=CMP_GT
           endif
          endif
         endif
         return
        end function cmp_integers
!-------------------------------------------------
        function cmp_real(real1,real2) result(cmp)
!Generic comparator for reals or arbitrary kinds with finite zero threshold.
         implicit none
         integer(INTD):: cmp                  !out: result of comparison: {CMP_EQ,CMP_LT,CMP_GT,CMP_ER}
         class(*), intent(in), target:: real1 !in: real 1
         class(*), intent(in), target:: real2 !in: real 2
         real(8):: r1,r2,diff

         cmp=CMP_ER
         select type(real1)
         type is(real(4))
          r1=real1; cmp=CMP_EQ
         type is(real(8))
          r1=real1; cmp=CMP_EQ
         end select
         if(cmp.eq.CMP_EQ) then
          cmp=CMP_ER
          select type(real2)
          type is(real(4))
           r2=real2; cmp=CMP_EQ
          type is(real(8))
           r2=real2; cmp=CMP_EQ
          end select
          if(cmp.eq.CMP_EQ) then
           diff=r1-r2
           if(diff.lt.-DP_ZERO_THRESH) then
            cmp=CMP_LT
           elseif(diff.gt.DP_ZERO_THRESH) then
            cmp=CMP_GT
           endif
          endif
         endif
         return
        end function cmp_real
!--------------------------------------------------
        function cmp_strings(str1,str2) result(cmp)
!Generic comparator for strings.
         implicit none
         integer(INTD):: cmp                 !out: result of comparison: {CMP_EQ,CMP_LT,CMP_GT,CMP_ER}
         class(*), intent(in), target:: str1 !in: string 1
         class(*), intent(in), target:: str2 !in: string 2

         cmp=CMP_ER
         select type(str1)
         type is(character(*))
          select type(str2)
          type is(character(*))
           cmp=str_cmp(str1,str2)
          end select
         end select
         return
        end function cmp_strings
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
!---------------------------------------------------------
        function cmp_tens_descriptors(td1,td2) result(cmp)
!Comparator for tensor descriptors.
         implicit none
         integer(INTD):: cmp                !out: result of comparison: {CMP_EQ,CMP_LT,CMP_GT,CMP_ER}
         class(*), intent(in), target:: td1 !in: tensor descriptor 1
         class(*), intent(in), target:: td2 !in: tensor descriptor 2
         class(tens_descr_t), pointer:: tdp1,tdp2

         tdp1=>NULL(); tdp2=>NULL()
         select type(td1); class is(tens_descr_t); tdp1=>td1; end select
         select type(td2); class is(tens_descr_t); tdp2=>td2; end select
         if(associated(tdp1).and.associated(tdp2)) then
          cmp=tdp1%compare(tdp2)
         else
          cmp=CMP_ER
         endif
         return
        end function cmp_tens_descriptors
!---------------------------------------------------------------------------
        function build_test_hspace(space_name,ierr,space_p) result(space_id)
!Builds a hierarchical vector space for testing/debugging.
         implicit none
         integer(INTD):: space_id                                   !out: registered space id
         character(*), target, intent(in):: space_name              !in: vector space name
         integer(INTD), intent(out), optional:: ierr                !out: error code
         class(h_space_t), pointer, intent(out), optional:: space_p !out: pointer to the hierarchical vector space
         class(h_space_t), pointer:: hsp
         integer(INTD):: errc
         integer(INTL):: l

         hsp=>NULL()
         space_id=hspace_register%register_space(space_name,errc,hsp)
         if(errc.eq.TEREC_SUCCESS.and.associated(hsp)) then
          call hspace_basis_%subspace_basis_ctor(HSPACE_DIM_,errc)
          if(errc.eq.0) then
           do l=1_INTL,HSPACE_DIM_
            call hspace_symm_(l)%spher_symmetry_ctor(int((l-1)/5,INTD),0,errc); if(errc.ne.0) exit
            call hspace_basis_%set_basis_func(l,BASIS_ABSTRACT,errc,symm=hspace_symm_(l)); if(errc.ne.0) exit
           enddo
           if(errc.eq.0) then
            call hspace_basis_%finalize(errc)
            if(errc.eq.0) then
             !write(CONS_OUT,'("#DEBUG(build_test_hspace): Building test H-space from a basis ... ")') !debug
             call hsp%h_space_ctor(hspace_basis_,errc)
             !write(CONS_OUT,'(i9,"-dimensional full space -> ",i9," subspaces:")') hsp%get_space_dim(),hsp%get_num_subspaces() !debug
             !call hsp%print_it() !debug
            endif
           endif
          endif
         else
          if(errc.eq.TEREC_SUCCESS) errc=TEREC_ERROR
         endif
         if(present(space_p)) space_p=>hsp
         if(present(ierr)) ierr=errc
         return
        end function build_test_hspace
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
!-------------------------------------------------------
        subroutine print_tcg_buffer(start,finish,length)
!Dumps the current content of the TCG buffer.
         implicit none
         integer(INTD), intent(in):: start
         integer(INTD), intent(in):: finish
         integer(INTD), intent(in):: length
         integer(INTD):: i

         write(CONS_OUT,'("PRINTING TCG BUFFER: ",i8," - ",i8)') start,finish
         if(allocated(tcg_num_buf)) then
          if(start.ge.lbound(tcg_num_buf,1).and.finish.le.ubound(tcg_num_buf,1)) then
           do i=start,finish
            write(CONS_OUT,'(i8,3x,32(1x,i5))') tcg_num_buf(i),tcg_ind_buf(1:length,i)
           enddo
          endif
         endif
         write(CONS_OUT,'("END OF PRINTING")')
         return
        end subroutine print_tcg_buffer
![hspace_register_t]=========================================================================
        subroutine HspaceRegisterInit(this,ierr)
!Initializes the register of hierarchical vector spaces.
         implicit none
         class(hspace_register_t), intent(inout):: this !inout: register of hierarchical vector spaces
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD):: errc

         if(this%initialized) then
          errc=TEREC_INVALID_REQUEST
         else
          errc=this%name2id_it%init(this%name2id)
          if(errc.eq.GFC_SUCCESS) then
           errc=this%hspaces_it%init(this%hspaces)
           if(errc.eq.GFC_SUCCESS) this%initialized=.TRUE.
          endif
         endif
         if(errc.ne.TEREC_SUCCESS) call hspace_register_dtor(this)
         if(present(ierr)) ierr=errc
         return
        end subroutine HspaceRegisterInit
!--------------------------------------------------------------------------------------------
        function HspaceRegisterRegisterSpace(this,space_name,ierr,hspace_p) result(hspace_id)
!Registers an empty hierarchical vector space under the name <space_name> and, optionally,
!returns a pointer to the just registered space for a subsequent definition.
         implicit none
         integer(INTD):: hspace_id                                   !out: registered id of the hierarchical vector space: [0..max]
         class(hspace_register_t), intent(inout):: this              !inout: register of hierarchical vector spaces
         character(*), target, intent(in):: space_name               !in: space name
         integer(INTD), intent(out), optional:: ierr                 !out: error code
         class(h_space_t), pointer, intent(out), optional:: hspace_p !out: pointer to the just registered empty hierarchical vector space (for further construction)
         integer(INTD):: errc
         type(h_space_t), target:: hspace_empty
         class(*), pointer:: up
         !type(h_space_t), pointer:: hptr !debug

         errc=TEREC_SUCCESS
         if(.not.this%initialized) call this%init(errc)
         if(errc.eq.TEREC_SUCCESS) then
          hspace_id=this%hspaces_it%get_length(errc)
          if(errc.eq.GFC_SUCCESS) then
           errc=this%name2id_it%search(GFC_DICT_ADD_IF_NOT_FOUND,cmp_strings,space_name,hspace_id)
           if(errc.eq.GFC_NOT_FOUND) then
            !write(*,'("#DEBUG(hspace_register_t.register_space): Appending a local empty h_space_t ...")') !debug
            !hptr=>hspace_empty; call dump_bytes(c_loc(hptr),size_of(hspace_empty),'dump0') !debug
            errc=this%hspaces_it%append(hspace_empty)
            if(errc.eq.GFC_SUCCESS) then
             if(present(hspace_p)) then
              hspace_p=>NULL()
              up=>this%hspaces_it%element_value(int(hspace_id,INTL),errc)
              if(errc.eq.GFC_SUCCESS.and.associated(up)) then
               !select type(up); type is(h_space_t); hptr=>up; end select !debug
               !call dump_bytes(c_loc(hptr),size_of(hspace_empty),'dump1') !debug
               !hptr=>hspace_empty; call dump_bytes(c_loc(hptr),size_of(hspace_empty),'dump2') !debug
               select type(up); type is(h_space_t); hspace_p=>up; end select
               if(.not.associated(hspace_p)) errc=TEREC_ERROR
              else
               if(errc.eq.GFC_SUCCESS) errc=TEREC_ERROR
              endif
             endif
            endif
           else
            if(errc.eq.GFC_FOUND) errc=TEREC_INVALID_REQUEST
           endif
          endif
         endif
         if(errc.ne.TEREC_SUCCESS) then
          if(present(hspace_p)) hspace_p=>NULL()
          hspace_id=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function HspaceRegisterRegisterSpace
!--------------------------------------------------------------------------------
        function HspaceRegisterGetSpaceId(this,space_name,ierr) result(hspace_id)
!Given the name of a registered hierarchical space, returns its id.
         implicit none
         integer(INTD):: hspace_id                      !out: space id
         class(hspace_register_t), intent(inout):: this !in: register of hierarchical spaces
         character(*), intent(in), target:: space_name  !in: name of the hierarchical space
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD):: errc
         class(*), pointer:: up

         errc=TEREC_SUCCESS; hspace_id=-1
         if(.not.this%initialized) call this%init(errc)
         if(errc.eq.TEREC_SUCCESS) then
          errc=this%name2id_it%search(GFC_DICT_JUST_FIND,cmp_strings,space_name,value_out=up)
          if(errc.eq.GFC_FOUND) then
           !if(associated(up)) then
            select type(up); type is(integer); hspace_id=int(up,INTD); class default; errc=TEREC_ERROR; end select
           !else
            !errc=TEREC_ERROR
           !endif
          else
           errc=TEREC_INVALID_ARGS
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end function HspaceRegisterGetSpaceId
!---------------------------------------------------------------------------------------------
        function HspaceRegisterGetSpaceByName(this,space_name,ierr,hspace_id) result(hspace_p)
!Given the name of a registered hierarchical space, returns a pointer to the stored space itself.
         implicit none
         class(h_space_t), pointer:: hspace_p             !out: pointer to the stored hierarchical vector space
         class(hspace_register_t), intent(inout):: this   !in: register of hierarchical vector spaces
         character(*), intent(in), target:: space_name    !in: space name
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD), intent(out), optional:: hspace_id !out: registered space id: [0..max], -1:unregistered
         integer(INTD):: errc,hid
         class(*), pointer:: up

         hspace_p=>NULL(); hid=this%get_space_id(space_name,errc)
         if(errc.eq.TEREC_SUCCESS) then
          up=>this%hspaces_it%element_value(int(hid,INTL),errc)
          if(errc.eq.GFC_SUCCESS) then
           select type(up); class is(h_space_t); hspace_p=>up; end select
           if(.not.associated(hspace_p)) errc=TEREC_ERROR
          endif
         endif
         if(present(hspace_id)) hspace_id=hid
         if(present(ierr)) ierr=errc
         return
        end function HspaceRegisterGetSpaceByName
!-------------------------------------------------------------------------------
        function HspaceRegisterGetSpaceById(this,space_id,ierr) result(hspace_p)
!Returns a pointer to a hierarchical vector space stored in the register with id <space_id>.
         implicit none
         class(h_space_t), pointer:: hspace_p             !out: pointer to the stored hierarchical vector space
         class(hspace_register_t), intent(inout):: this   !in: register of hierarchical vector spaces
         integer(INTD), intent(in):: space_id             !in: registered space id
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD):: errc
         class(*), pointer:: up

         errc=TEREC_SUCCESS; hspace_p=>NULL()
         if(.not.this%initialized) call this%init(errc)
         if(errc.eq.TEREC_SUCCESS) then
          up=>this%hspaces_it%element_value(int(space_id,INTL),errc)
          if(errc.eq.GFC_SUCCESS) then
           select type(up); class is(h_space_t); hspace_p=>up; end select
           if(.not.associated(hspace_p)) errc=TEREC_ERROR
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end function HspaceRegisterGetSpaceById
!--------------------------------------------
        subroutine hspace_register_dtor(this)
         implicit none
         type(hspace_register_t):: this
         integer(INTD):: errc

         errc=this%name2id_it%get_status()
         if(errc.ne.GFC_IT_NULL) then
          call this%name2id_it%delete_all(errc)
          errc=this%name2id_it%release()
         endif
         errc=this%hspaces_it%get_status()
         if(errc.ne.GFC_IT_NULL) then
          errc=this%hspaces_it%delete_all()
          errc=this%hspaces_it%release()
         endif
         this%initialized=.FALSE.
         return
        end subroutine hspace_register_dtor
![hspace_reg_t]=====================================
        subroutine HspaceRegCtor(this,space_id,ierr)
         implicit none
         class(hspace_reg_t), intent(out):: this     !out: registered hierarchical space
         integer(INTD), intent(in):: space_id        !in: registered space id
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         this%space_id=space_id
         this%hspace_p=>hspace_register%get_space(space_id,errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine HspaceRegCtor
!-------------------------------------------------------
        subroutine HspaceRegCtorUnpack(this,packet,ierr)
!Unpacks an object from the packet.
         implicit none
         class(hspace_reg_t), intent(out):: this     !out: registered hierarchical space
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,hid

         call unpack_builtin(packet,hid,errc)
         if(errc.eq.PACK_SUCCESS) call this%hspace_reg_ctor(hid,errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine HspaceRegCtorUnpack
!-------------------------------------------------
        subroutine HspaceRegPack(this,packet,ierr)
!Packs the object into a packet.
         implicit none
         class(hspace_reg_t), intent(in):: this      !in: registered hierarchical space
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         call pack_builtin(packet,this%space_id,errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine HspaceRegPack
!-----------------------------------------------------
        function HspaceRegIsSet(this,ierr) result(ans)
!Returns TRUE if the object is set.
         implicit none
         logical:: ans                               !out: answer
         class(hspace_reg_t), intent(in):: this      !in: registered hierarchical space
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS; ans=(this%space_id.ge.0)
         if(present(ierr)) ierr=errc
         return
        end function HspaceRegIsSet
!----------------------------------------------------------------
        function HspaceRegGetSpaceId(this,ierr) result(hspace_id)
!Returns the registered id of the hierarchical space.
         implicit none
         integer(INTD):: hspace_id                   !out: space id
         class(hspace_reg_t), intent(in):: this      !in: registered hierarchical space
         integer(INTD), intent(out), optional:: ierr !out: error code

         hspace_id=this%space_id
         if(present(ierr).and.hspace_id.lt.0) ierr=TEREC_INVALID_REQUEST
         return
        end function HspaceRegGetSpaceId
!-------------------------------------------------------------
        function HspaceRegGetSpace(this,ierr) result(hspace_p)
!Returns a pointer to the stored hierarchical vector space.
         implicit none
         class(h_space_t), pointer:: hspace_p        !out: pointer to the hierarchical space
         class(hspace_reg_t), intent(in):: this      !in: registered hierarchical space
         integer(INTD), intent(out), optional:: ierr !out: error code

         hspace_p=>this%hspace_p
         if(present(ierr).and.(.not.associated(hspace_p))) ierr=TEREC_INVALID_REQUEST
         return
        end function HspaceRegGetSpace
!---------------------------------------
        subroutine hspace_reg_dtor(this)
         implicit none
         type(hspace_reg_t):: this

         this%space_id=-1; this%hspace_p=>NULL()
         return
        end subroutine hspace_reg_dtor
![tens_signature_t]========================================================
        subroutine TensSignatureCtor(this,ierr,subspaces,tens_name,hspaces)
!CTOR for tens_signature_t.
         implicit none
         class(tens_signature_t), intent(out):: this              !out: tensor signature
         integer(INTD), intent(out), optional:: ierr              !out: error code
         integer(INTL), intent(in), optional:: subspaces(1:)      !in: multi-index of subspaces
         character(*), intent(in), optional:: tens_name           !in: alphanumeric_ tensor name (no spaces allowed!)
         integer(INTD), intent(in), optional:: hspaces(1:)        !in: hierarchical vector spaces id's
         integer(INTD):: errc,n,i

         errc=TEREC_SUCCESS
         n=0; if(present(subspaces)) n=size(subspaces)
         if(n.gt.0) then !true tensor
          allocate(this%space_idx(1:n),STAT=errc)
          if(errc.eq.0) then
           this%space_idx(1:n)=subspaces(1:n)
           this%num_dims=n
          else
           errc=TEREC_MEM_ALLOC_FAILED
          endif
          if(errc.eq.TEREC_SUCCESS.and.present(hspaces)) then
           allocate(this%hspace(1:n),STAT=errc)
           if(errc.eq.0) then
            do i=1,n
             call this%hspace(i)%hspace_reg_ctor(hspaces(i),errc); if(errc.ne.TEREC_SUCCESS) exit
            enddo
           else
            errc=TEREC_MEM_ALLOC_FAILED
           endif
          endif
         else !scalar tensor
          this%num_dims=0
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
         integer(INTD):: i,nd,errc,hidx(1:MAX_TENSOR_RANK)
         integer(INTL):: sidx(1:MAX_TENSOR_RANK)
         character(TEREC_MAX_TENS_NAME_LEN):: tname
         logical:: pcn,hsn

         tname=' '
         call unpack_builtin(packet,nd,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,pcn,errc)
         if(errc.eq.PACK_SUCCESS) then
          if(nd.gt.0) then
           call unpack_builtin(packet,hsn,errc)
           if(errc.eq.PACK_SUCCESS) then
            do i=1,nd
             call unpack_builtin(packet,sidx(i),errc); if(errc.ne.PACK_SUCCESS) exit
            enddo
            if(errc.eq.PACK_SUCCESS.and.hsn) then
             do i=1,nd
              call unpack_builtin(packet,hidx(i),errc); if(errc.ne.PACK_SUCCESS) exit
             enddo
            endif
           endif
          endif
          if(errc.eq.PACK_SUCCESS.and.pcn) call unpack_builtin(packet,tname,errc)
          if(errc.eq.PACK_SUCCESS) then
           if(nd.gt.0) then
            if(pcn) then
             if(hsn) then
              call this%tens_signature_ctor(errc,sidx(1:nd),tname(1:len_trim(tname)),hidx(1:nd))
             else
              call this%tens_signature_ctor(errc,sidx(1:nd),tname(1:len_trim(tname)))
             endif
            else
             if(hsn) then
              call this%tens_signature_ctor(errc,sidx(1:nd),hspaces=hidx(1:nd))
             else
              call this%tens_signature_ctor(errc,sidx(1:nd))
             endif
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
! + logical {TRUE|FALSE}: whether or not the tensor is named;
! + [OPTIONAL]: logical {TRUE|FALSE}: whether or not .hspace(:) is present;
! + [OPTIONAL]: this%space_idx(1:this%num_dims);
! + [OPTIONAL]: this%hspace(1:this%num_dims)%space_id;
! + [OPTIONAL]: this%char_name;
         implicit none
         class(tens_signature_t), intent(in):: this  !in: tensor signature
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: i,errc
         logical:: pcn,hsn

         if(this%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) then
           call pack_builtin(packet,this%num_dims,errc)
           if(errc.eq.PACK_SUCCESS) then
            pcn=allocated(this%char_name)
            call pack_builtin(packet,pcn,errc)
            if(errc.eq.PACK_SUCCESS) then
             if(this%num_dims.gt.0) then
              hsn=allocated(this%hspace)
              call pack_builtin(packet,hsn,errc)
              if(errc.eq.PACK_SUCCESS) then
               if(allocated(this%space_idx)) then
                do i=1,this%num_dims
                 call pack_builtin(packet,this%space_idx(i),errc); if(errc.ne.PACK_SUCCESS) exit
                enddo
               else
                errc=TEREC_ERROR
               endif
               if(errc.eq.PACK_SUCCESS.and.hsn) then
                do i=1,this%num_dims
                 call pack_builtin(packet,this%hspace(i)%space_id,errc); if(errc.ne.PACK_SUCCESS) exit
                enddo
               endif
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
!--------------------------------------------------------------------------
        function TensSignatureIsSet(this,ierr,num_dims,hspaced) result(res)
!Returns TRUE if the tensor_signature_t object is set, FALSE otherwise.
         implicit none
         logical:: res                                   !out: result
         class(tens_signature_t), intent(in):: this      !in: tensor signature
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD), intent(out), optional:: num_dims !out: tensor rank (if set)
         logical, intent(out), optional:: hspaced        !out: TRUE if the tensor signature is defined over hierarchical vector spaces
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         res=(this%num_dims.ge.0)
         if(present(num_dims)) num_dims=this%num_dims
         if(present(hspaced)) hspaced=allocated(this%hspace)
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
!----------------------------------------------------------------------------
        subroutine TensSignatureGetSpec(this,subspaces,num_dims,ierr,hspaces)
!Returns the defining subspaces of the tensor (subspace multi-index).
         implicit none
         class(tens_signature_t), intent(in):: this                !in: tensor signature
         integer(INTL), intent(inout):: subspaces(1:)              !out: defining subspaces (their IDs)
         integer(INTD), intent(out):: num_dims                     !out: number of tensor dimensions
         integer(INTD), intent(out), optional:: ierr               !out: error code
         type(hspace_reg_t), intent(inout), optional:: hspaces(1:) !out: hierarchical vector spaces for each dimension
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         if(this%is_set()) then
          num_dims=this%num_dims
          if(size(subspaces).ge.num_dims) then
           subspaces(1:num_dims)=this%space_idx(1:num_dims)
           if(present(hspaces)) then
            if(size(hspaces).ge.num_dims) then
             hspaces(1:num_dims)=this%hspace(1:num_dims)
            else
             errc=TEREC_UNABLE_COMPLETE
            endif
           endif
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
!Both signatures must be associated with hierarchical vector spaces.
         implicit none
         integer(INTD):: relation                      !out: relation: {CMP_EQ,CMP_CN,CMP_IN,CMP_OV,CMP_NC}
         class(tens_signature_t), intent(in):: this    !in: tensor signature 1
         class(tens_signature_t), intent(in):: another !in: tensor signature 2
         integer(INTD):: n1,n2,nl1,nl2,ch1,ch2,i,cmp,errc
         integer(INTL):: s1,s2
         logical:: hs1,hs2

         errc=0
         if(this%is_set(num_dims=n1,hspaced=hs1).and.another%is_set(num_dims=n2,hspaced=hs2)) then
          if(hs1.and.hs2.and.n1.eq.n2) then
           relation=CMP_EQ
!Compare names:
           nl1=len(this%char_name); nl2=len(another%char_name)
           if(nl1.eq.nl2) then
            do i=1,nl1
             ch1=iachar(this%char_name(i:i)); ch2=iachar(another%char_name(i:i))
             if(ch1.ne.ch2) then; relation=CMP_NC; exit; endif
            enddo
            if(relation.eq.CMP_EQ) then
!Compare specs:
 !Compare spaces:
             do i=1,n1
              if(this%hspace(i)%space_id.lt.0) then; relation=CMP_ER; exit; endif
              if(.not.associated(this%hspace(i)%hspace_p)) then; relation=CMP_ER; exit; endif !trap
              if(this%hspace(i)%space_id.ne.another%hspace(i)%space_id) then
               relation=CMP_NC; exit
              endif
             enddo
 !Compare subspaces:
             if(relation.eq.CMP_EQ) then
              do i=1,n1
               s1=this%space_idx(i); s2=another%space_idx(i)
               cmp=this%hspace(i)%hspace_p%relate_subspaces(s1,s2,errc)
               if(errc.ne.0) then; cmp=CMP_ER; else; if(cmp.eq.CMP_OV) cmp=CMP_NC; endif
               if(cmp.eq.CMP_ER.or.cmp.eq.CMP_NC) then; relation=cmp; exit; endif
               if(cmp.ne.CMP_EQ) then !{CMP_CN,CMP_IN}
                if(relation.eq.CMP_EQ) then
                 relation=cmp
                else
                 if(cmp.ne.relation) relation=CMP_OV !overlap
                endif
               endif
              enddo
             endif
            endif
           else
            relation=CMP_NC
           endif
          else
           relation=CMP_NC
          endif
         else
          relation=CMP_ER
         endif
         return
        end function TensSignatureRelate
!--------------------------------------------------------------
        function TensSignatureCompare(this,another) result(cmp)
!Compares the tensor signature with another tensor signature (formal comparator).
         implicit none
         integer(INTD):: cmp                           !out: comparison result: {CMP_EQ,CMP_LT,CMP_GT,CMP_ER}
         class(tens_signature_t), intent(in):: this    !in: tensor signature 1
         class(tens_signature_t), intent(in):: another !in: tensor signature 2
         integer(INTD):: n1,n2,nl1,nl2,ch1,ch2,i,errc
         integer(INTL):: s1,s2
         logical:: hs1,hs2

         errc=0
         if(this%is_set(num_dims=n1,hspaced=hs1).and.another%is_set(num_dims=n2,hspaced=hs2)) then
          cmp=CMP_EQ
!Compare names:
          nl1=len(this%char_name); nl2=len(another%char_name)
          if(nl1.lt.nl2) then; cmp=CMP_LT; elseif(nl1.gt.nl2) then; cmp=CMP_GT; endif
          if(cmp.eq.CMP_EQ) then
           do i=1,nl1
            ch1=iachar(this%char_name(i:i)); ch2=iachar(another%char_name(i:i))
            if(ch1.lt.ch2) then; cmp=CMP_LT; exit; elseif(ch1.gt.ch2) then; cmp=CMP_GT; exit; endif
           enddo
           if(cmp.eq.CMP_EQ) then
!Compare specs:
            if(n1.lt.n2) then
             cmp=CMP_LT
            elseif(n1.gt.n2) then
             cmp=CMP_GT
            else
             if(hs1.and.hs2) then !both signatures are over hierarchical spaces
              do i=1,this%num_dims
               if(this%hspace(i)%space_id.lt.another%hspace(i)%space_id) then
                cmp=CMP_LT; exit
               elseif(this%hspace(i)%space_id.gt.another%hspace(i)%space_id) then
                cmp=CMP_GT; exit
               endif
              enddo
              if(cmp.eq.CMP_EQ) then
               do i=1,this%num_dims
                s1=this%space_idx(i); s2=another%space_idx(i)
                cmp=this%hspace(i)%hspace_p%compare_subspaces(s1,s2,errc); if(errc.ne.0) cmp=CMP_ER
                if(cmp.ne.CMP_EQ) exit
               enddo
              endif
             elseif(.not.(hs1.or.hs2)) then !both signatures are over the default space
              do i=1,this%num_dims
               if(this%space_idx(i).lt.another%space_idx(i)) then
                cmp=CMP_LT; exit
               elseif(this%space_idx(i).gt.another%space_idx(i)) then
                cmp=CMP_GT; exit
               endif
              enddo
             else
              if(hs1) then; cmp=CMP_LT; else; cmp=CMP_GT; endif
             endif
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
         if(allocated(this%hspace)) deallocate(this%hspace)
         if(allocated(this%char_name)) deallocate(this%char_name)
         if(allocated(this%space_idx)) deallocate(this%space_idx)
         this%num_dims=-1
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
!------------------------------------------------------------------------------------------
        subroutine TensShapeGetGroups(this,num_dims,num_groups,dim_groups,ierr,group_restr)
!Returns symmetric dimension groups and group restrictions.
         implicit none
         class(tens_shape_t), intent(in):: this                   !in: tensor shape
         integer(INTD), intent(out):: num_dims                    !out: tensor rank
         integer(INTD), intent(out):: num_groups                  !out: number of index groups
         integer(INTD), intent(inout):: dim_groups(1:)            !out: dimension groups
         integer(INTD), intent(out), optional:: ierr              !out: error code
         integer(INTD), intent(inout), optional:: group_restr(1:) !out: group restrictions
         integer(INTD):: errc

         if(this%is_set(errc,num_dims=num_dims)) then
          if(errc.eq.TEREC_SUCCESS) then
           num_groups=this%num_grps
           dim_groups(1:num_dims)=this%dim_group(1:num_dims)
           if(present(group_restr).and.num_groups.gt.0) group_restr=this%group_spec(1:num_groups)
          endif
         else
          if(errc.eq.TEREC_SUCCESS) errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensShapeGetGroups
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
        subroutine TensHeaderCtor(this,ierr,tens_name,subspaces,hspaces,dim_extent,dim_group,group_spec)
!CTOR for tens_header_t. Each subsequent optional argument implies the existence of all preceding
!optional arguments, except <ierr> and <tens_name>. If no optional arguments are present, except
!maybe <tens_name> and/or <ierr>, a scalar header will be constructed. <dim_group> and <group_spec>
!must either be both present or both absent. More specifically:
! # Constructing a scalar tensor header: Do not pass any optional arguments except <tens_name> and/or <ierr>;
! # Constructing a true tensor header without shape: Pass only <subspaces>, <hspaces>, and optionally <tens_name> and/or <ierr>;
! # Constructing a true tensor header with a shape: Pass <subspaces>, <hspaces>, and <dim_extent> with all other arguments optional.
!   Note that it is ok to pass dimension extents equal to 0 for unresolved tensor dimensions (to be set later).
         implicit none
         class(tens_header_t), intent(out):: this             !out: tensor header
         integer(INTD), intent(out), optional:: ierr          !out: error code
         character(*), intent(in), optional:: tens_name       !in: alphanumeric_ tensor name
         integer(INTL), intent(in), optional:: subspaces(1:)  !in: subspace multi-index (specification): Length = tensor rank
         integer(INTD), intent(in), optional:: hspaces(1:)    !in: hierarchical vector space id for each tensor dimension
         integer(INTL), intent(in), optional:: dim_extent(1:) !in: dimension extents: Length = tensor rank
         integer(INTD), intent(in), optional:: dim_group(1:)  !in: dimension restriction groups: Length = tensor rank
         integer(INTD), intent(in), optional:: group_spec(1:) !in: dimension restriction group specification
         integer(INTD):: errc,m
         logical:: pr_nam,pr_sub,pr_hsp,pr_dim,pr_grp,pr_grs

         errc=TEREC_SUCCESS
         pr_nam=present(tens_name)
         pr_sub=present(subspaces)
         pr_hsp=present(hspaces)
         pr_dim=present(dim_extent)
         pr_grp=present(dim_group)
         pr_grs=present(group_spec)
         if(pr_hsp.and.(.not.pr_sub)) errc=TEREC_INVALID_ARGS
         if(pr_dim.and.(.not.pr_sub)) errc=TEREC_INVALID_ARGS
         if((pr_grp.or.pr_grs).and.(.not.pr_dim)) errc=TEREC_INVALID_ARGS
         if((pr_grp.and.(.not.pr_grs)).or.(pr_grs.and.(.not.pr_grp))) errc=TEREC_INVALID_ARGS
         if(pr_sub) then
          m=size(subspaces)
          if(m.le.0) then; pr_sub=.FALSE.; pr_hsp=.FALSE.; pr_dim=.FALSE.; pr_grp=.FALSE.; pr_grs=.FALSE.; endif
          if(pr_dim) then; if(m.ne.size(dim_extent)) errc=TEREC_INVALID_ARGS; endif
         endif
         if(errc.eq.TEREC_SUCCESS) then
 !tensor signature:
          if(pr_sub) then !explicit tensor
           if(pr_nam) then
            if(pr_hsp) then
             call this%signature%tens_signature_ctor(errc,subspaces,tens_name,hspaces)
            else
             call this%signature%tens_signature_ctor(errc,subspaces,tens_name) !`this will never happen
            endif
           else
            if(pr_hsp) then
             call this%signature%tens_signature_ctor(errc,subspaces,hspaces=hspaces)
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
!----------------------------------------------------------------------------------------------------
        function TensHeaderIsSet(this,ierr,num_dims,num_groups,shaped,unresolved,hspaced) result(res)
!Returns TRUE if the tensor header is set (with or without shape), plus additional info.
         implicit none
         logical:: res                                     !out: result
         class(tens_header_t), intent(in):: this           !in: tensor header
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD), intent(out), optional:: num_dims   !out: number of dimensions
         integer(INTD), intent(out), optional:: num_groups !out: number of restricted dimension groups
         logical, intent(out), optional:: shaped           !out: TRUE if the tensor shape is set
         integer(INTD), intent(out), optional:: unresolved !out: number of unresolved tensor dimensions
         logical, intent(out), optional:: hspaced          !out: TRUE if the tensor dimensions are over hierarchical spaces
         integer(INTD):: errc,nd,ng,unres
         logical:: shpd,hspc

         res=this%signature%is_set(errc,num_dims=nd,hspaced=hspc)
         if(present(num_dims)) num_dims=nd
         if(res.and.errc.eq.TEREC_SUCCESS) then
          shpd=this%shape%is_set(errc,num_dims=nd,num_groups=ng,unresolved=unres)
         else
          shpd=.FALSE.; hspc=.FALSE.; ng=0; unres=-1
         endif
         if(present(num_groups)) num_groups=ng
         if(present(shaped)) shaped=shpd
         if(present(unresolved)) unresolved=unres
         if(present(hspaced)) hspaced=hspc
         if(present(ierr)) ierr=errc
         return
        end function TensHeaderIsSet
!--------------------------------------------------------
        function TensHeaderIsValid(this,ierr) result(ans)
!Returns TRUE if the set tensor header is valid.
         implicit none
         logical:: ans                               !out: answer
         class(tens_header_t), intent(in):: this     !in: tensor header
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,nd,ng,unres,i,j,fl(1:MAX_TENSOR_RANK)
         logical:: shpd,hspc

         ans=.FALSE.
         if(this%is_set(errc,nd,ng,shpd,unres,hspc)) then
          if(errc.eq.TEREC_SUCCESS) then
           if(nd.gt.0) then
 !Signature check:
            if(allocated(this%signature%space_idx)) then
             do i=1,nd; if(this%signature%space_idx(i).lt.0_INTL) return; enddo
             if(allocated(this%signature%hspace)) then
              do i=1,nd
               if(this%signature%hspace(i)%space_id.lt.0) return
               if(.not.associated(this%signature%hspace(i)%hspace_p)) return
              enddo
             endif
 !Shape:
             if(shpd) then
              do i=1,nd; if(this%shape%dim_extent(i).lt.0) return; enddo
 !Shape/signature consistency:
              if(ng.gt.0) then
               if(ng.gt.MAX_TENSOR_RANK) return
               do i=1,nd; if(this%shape%dim_group(i).lt.0) return; enddo
               if(hspc) then
                fl(1:ng)=-1
                do i=1,nd
                 j=this%shape%dim_group(i)
                 if(j.gt.0) then
                  if(fl(j).lt.0) then
                   fl(j)=this%signature%hspace(i)%space_id
                  else
                   if(fl(j).ne.this%signature%hspace(i)%space_id) return
                  endif
                 endif
                enddo
               endif
              endif
             endif
            endif
           else !scalar
            if(ng.le.0.and.unres.le.0) ans=.TRUE.
           endif
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensHeaderIsValid
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
!-------------------------------------------------------------------------
        subroutine TensHeaderGetSpec(this,subspaces,num_dims,ierr,hspaces)
!Returns the defining subspaces of the tensor (subspace multi-index).
         implicit none
         class(tens_header_t), intent(in):: this                   !in: tensor header
         integer(INTL), intent(inout):: subspaces(1:)              !out: defining subspaces (their IDs)
         integer(INTD), intent(out):: num_dims                     !out: number of tensor dimensions
         integer(INTD), intent(out), optional:: ierr               !out: error code
         type(hspace_reg_t), intent(inout), optional:: hspaces(1:) !out: hierarchical vector space for each tensor dimension
         integer(INTD):: errc

         if(present(hspaces)) then
          call this%signature%get_spec(subspaces,num_dims,errc,hspaces)
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
!-------------------------------------------------------------------------------------------
        subroutine TensHeaderGetGroups(this,num_dims,num_groups,dim_groups,ierr,group_restr)
!Returns symmetric dimension groups and group restrictions.
         implicit none
         class(tens_header_t), intent(in):: this                  !in: tensor header
         integer(INTD), intent(out):: num_dims                    !out: tensor rank
         integer(INTD), intent(out):: num_groups                  !out: number of groups
         integer(INTD), intent(inout):: dim_groups(1:)            !out: dimension groups
         integer(INTD), intent(out), optional:: ierr              !out: error code
         integer(INTD), intent(inout), optional:: group_restr(1:) !out: group restrictions
         integer(INTD):: errc

         if(present(group_restr)) then
          call this%shape%get_groups(num_dims,num_groups,dim_groups,errc,group_restr)
         else
          call this%shape%get_groups(num_dims,num_groups,dim_groups,errc)
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensHeaderGetGroups
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
         class(DataDescr_t), intent(in):: data_descr !in: DDSS data descriptor for the tensor body (will be cloned)
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
!--------------------------------------------------------
        subroutine TensLayoutUnpackBase(this,packet,ierr)
!Unpacks the object from a packet.
         implicit none
         class(tens_layout_t), intent(out):: this    !out: tensor body layout
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         logical:: dda

         call unpack_builtin(packet,this%layout,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%data_type,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,dda,errc)
         if(errc.eq.PACK_SUCCESS.and.dda) then
          if(.not.allocated(this%data_descr)) allocate(this%data_descr)
          call this%data_descr%unpack(packet,errc)
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensLayoutUnpackBase
!------------------------------------------------------
        subroutine TensLayoutPackBase(this,packet,ierr)
!Packs the object into a packet.
         implicit none
         class(tens_layout_t), intent(in):: this     !in: tensor body layout
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         logical:: dda

         call pack_builtin(packet,this%layout,errc)
         if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%data_type,errc)
         if(errc.eq.PACK_SUCCESS) then
          dda=allocated(this%data_descr)
          call pack_builtin(packet,dda,errc)
          if(errc.eq.PACK_SUCCESS.and.dda) call this%data_descr%pack(packet,errc)
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensLayoutPackBase
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
         integer(INTD):: errc,dtk,dts
         logical:: locd

         body_size=0_INTL
         if(this%is_set(errc,locd)) then
          if(locd) then !located tensor
           if(this%data_descr%is_set()) then
            body_size=this%data_descr%data_size(errc)
           else
            errc=TEREC_INVALID_REQUEST
           endif
          else !not located (only laid out)
           dtk=this%get_data_type(errc)
           if(errc.eq.TEREC_SUCCESS) then
            if(tens_valid_data_kind(dtk,dts).eq.YEP) then
             body_size=this%get_volume()*dts
            else
             errc=TEREC_OBJ_CORRUPTED
            endif
           endif
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensLayoutGetBodySize
!-----------------------------------------------------------------
        function TensLayoutGetDataDescr(this,ierr) result(descr_p)
!Returns a pointer to the DDSS data descriptor.
         implicit none
         class(DataDescr_t), pointer:: descr_p           !out: pointer to the DDSS data descriptor
         class(tens_layout_t), intent(in), target:: this !in: tensor layout
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc
         logical:: locd

         errc=TEREC_SUCCESS; descr_p=>NULL()
         if(this%is_set(errc,locd)) then
          if(locd) then
           if(errc.eq.TEREC_SUCCESS) descr_p=>this%data_descr
          else
           errc=TEREC_INVALID_REQUEST
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensLayoutGetDataDescr
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
!---------------------------------------------------------------------------
        subroutine TensLayoutFdimsCtorUnpack(this,packet,ierr,tens_header_p)
!Unpacks the object from a packet.
         implicit none
         class(tens_layout_fdims_t), intent(out):: this                      !out: tensor body layout
         class(obj_pack_t), intent(inout):: packet                           !inout: packet
         integer(INTD), intent(out), optional:: ierr                         !out: error code
         class(tens_header_t), pointer, intent(in), optional:: tens_header_p !in: pointer to the corresponding tensor header
         integer(INTD):: errc

         call this%unpack_base(packet,errc)
         this%header=>NULL(); if(present(tens_header_p)) this%header=>tens_header_p
         if(present(ierr)) ierr=errc
         return
        end subroutine TensLayoutFdimsCtorUnpack
!-------------------------------------------------------
        subroutine TensLayoutFdimsPack(this,packet,ierr)
!Packs the object into a packet.
         implicit none
         class(tens_layout_fdims_t), intent(in):: this !in: tensor body layout
         class(obj_pack_t), intent(inout):: packet     !inout: packet
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: errc

         call this%pack_base(packet,errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensLayoutFdimsPack
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
!---------------------------------------------------------------
        subroutine TensLayoutFdimsUpdateHeader(this,header,ierr)
!Updates the tensor header the layout is associated with.
         implicit none
         class(tens_layout_fdims_t), intent(inout):: this  !inout: tensor body layout
         class(tens_header_t), intent(in), target:: header !in: tensor header
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         this%header=>header
         if(present(ierr)) ierr=errc
         return
        end subroutine TensLayoutFdimsUpdateHeader
!----------------------------------------------
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
!--------------------------------------------------------------------
        subroutine TensBodyCtorUnpack(this,packet,ierr,tens_header_p)
!Unpacks the object from a packet.
         implicit none
         class(tens_body_t), intent(out):: this      !out: tensor body
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         class(tens_header_t), pointer, intent(in), optional:: tens_header_p !in: pointer to the corresponding tensor header
         integer(INTD):: i,n,lay,errc
         type(list_iter_t):: lit
         type(tens_header_t):: thp
         class(tens_layout_fdims_t), pointer:: fl
         logical:: laid

         call unpack_builtin(packet,laid,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%state,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%num_subtensors,errc)
         if(errc.eq.PACK_SUCCESS.and.laid) then
          call unpack_builtin(packet,lay,errc)
          if(errc.eq.PACK_SUCCESS) then
           if(.not.allocated(this%layout)) then
            select case(lay)
            case(TEREC_LAY_FDIMS)
             allocate(tens_layout_fdims_t::this%layout,STAT=i)
             if(i.eq.0) then
              select type(lat=>this%layout)
              type is(tens_layout_fdims_t)
               if(present(tens_header_p)) then
                call lat%tens_layout_fdims_ctor(packet,errc,tens_header_p)
               else
                call lat%tens_layout_fdims_ctor(packet,errc)
               endif
              end select
             else
              errc=TEREC_MEM_ALLOC_FAILED
             endif
            case default
             errc=TEREC_ERROR
            end select
           else
            errc=TEREC_ERROR
           endif
          endif
         endif
         if(errc.eq.PACK_SUCCESS.and.this%num_subtensors.gt.0) then
          errc=lit%init(this%subtensors); n=this%num_subtensors
          do while(n.gt.0)
           call thp%tens_header_ctor(packet,errc); if(errc.ne.PACK_SUCCESS) exit
           errc=lit%append(thp); if(errc.ne.GFC_SUCCESS) exit
           n=n-1
          enddo
          i=lit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=i
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensBodyCtorUnpack
!------------------------------------------------
        subroutine TensBodyPack(this,packet,ierr)
!Packs the object into a packet.
         implicit none
         class(tens_body_t), intent(in):: this       !in: tensor body
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: i,errc
         type(list_iter_t):: lit
         class(tens_header_t), pointer:: thp
         class(*), pointer:: up
         logical:: laid

         laid=allocated(this%layout); call pack_builtin(packet,laid,errc)
         if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%state,errc)
         if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%num_subtensors,errc)
         if(errc.eq.PACK_SUCCESS.and.laid) then
          call pack_builtin(packet,this%layout%layout,errc)
          if(errc.eq.PACK_SUCCESS) then
           select type(lat=>this%layout)
           type is(tens_layout_fdims_t)
            call lat%pack(packet,errc)
           class default
            errc=TEREC_ERROR
           end select
          endif
         endif
         if(errc.eq.PACK_SUCCESS.and.this%num_subtensors.gt.0) then
          errc=lit%init(this%subtensors)
          do while(errc.eq.GFC_SUCCESS)
           up=>lit%get_value(errc); if(errc.ne.GFC_SUCCESS) exit
           select type(up); class is(tens_header_t); thp=>up; end select
           if(.not.associated(thp)) then; errc=GFC_ERROR; exit; endif
           call thp%pack(packet,errc); if(errc.ne.PACK_SUCCESS) then; errc=GFC_ERROR; exit; endif
           errc=lit%scanp(return_each=.TRUE.,skip_current=.TRUE.)
          enddo
          if(errc.eq.GFC_NO_MOVE) errc=GFC_SUCCESS
          i=lit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=i
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensBodyPack
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
           select case(layout_kind)
           case(TEREC_LAY_RECUR) !all constituent subtensors are mapped sequentially to a contiguous chunk of local memory
            stop !`Implement
           case(TEREC_LAY_FDIMS) !a single subtensor is mapped as "Fortran-dimension-led"
            if(this%num_subtensors.eq.1) then
             errc=lit%init(this%subtensors)
             if(errc.eq.GFC_SUCCESS) then
              up=>lit%get_value(errc)
              if(errc.eq.GFC_SUCCESS) then
               thp=>NULL(); select type(up); class is(tens_header_t); thp=>up; end select
               if(associated(thp)) then
                errc=lit%release()
                if(errc.eq.GFC_SUCCESS) then
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
                 errc=TEREC_ERROR
                endif
               else
                errc=TEREC_UNABLE_COMPLETE
               endif
              else
               errc=TEREC_UNABLE_COMPLETE
              endif
             else
              errc=TEREC_UNABLE_COMPLETE
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
           end select
          else
           errc=TEREC_INVALID_REQUEST
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensBodySetLayout
!----------------------------------------------------------------------
        subroutine TensBodySetLocation(this,data_descr,ierr,body_state)
!Sets the physical location of the tensor body via a DDSS data descriptor.
         implicit none
         class(tens_body_t), intent(inout):: this         !inout: tensor body
         class(DataDescr_t), intent(in):: data_descr      !in: DDSS data descriptor for tensor body (will be cloned)
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD), intent(in), optional:: body_state !in: tensor body value state
         integer(INTD):: errc
         logical:: layd,locd

         if(this%is_set(errc,layd,locd)) then
          if(errc.eq.TEREC_SUCCESS.and.layd.and.(.not.locd)) then
           call this%layout%set_location(data_descr,errc)
           if(errc.eq.TEREC_SUCCESS.and.present(body_state)) this%state=body_state
          else
           errc=TEREC_INVALID_REQUEST
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensBodySetLocation
!----------------------------------------------------------
        subroutine TensBodyResetState(this,ierr,body_state)
!Resets the tensor body value state.
         implicit none
         class(tens_body_t), intent(inout):: this         !inout: tensor body
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD), intent(in), optional:: body_state !in: new body state (none will reset to TEREC_BODY_UNDEF)
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         if(present(body_state)) then
          this%state=body_state
         else
          this%state=TEREC_BODY_UNDEF
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensBodyResetState
!--------------------------------------------------------------
        function TensBodyGetState(this,ierr) result(body_state)
!Returns the tensor body value state.
         implicit none
         integer(INTD):: body_state                  !out: tensor body value state
         class(tens_body_t), intent(in):: this       !in: tensor body
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         body_state=this%state
         if(present(ierr)) ierr=errc
         return
        end function TensBodyGetState
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
!-----------------------------------------------------------
        subroutine TensBodyPrintIt(this,ierr,dev_id,nspaces)
!Prints the tensor body info.
         implicit none
         class(tens_body_t), intent(in):: this         !in: tensor body
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD), intent(in), optional:: dev_id  !in: output device id (6:screen)
         integer(INTD), intent(in), optional:: nspaces !out: left alignment
         integer(INTD):: errc,dev,i
         class(tens_layout_t), pointer:: tens_layout

         errc=TEREC_SUCCESS
         if(present(dev_id)) then; dev=dev_id; else; dev=6; endif
         if(present(nspaces)) then
          do i=1,nspaces; write(dev,'(" ")',ADVANCE='NO'); enddo
          write(dev,'("TENSOR BODY{")')
          do i=1,nspaces; write(dev,'(" ")',ADVANCE='NO'); enddo
          write(dev,'(" Number of subtensors    = ",i7)') this%num_subtensors
          do i=1,nspaces; write(dev,'(" ")',ADVANCE='NO'); enddo
          write(dev,'(" Tensor body value state = ",i7)') this%state
          tens_layout=>this%get_layout()
          if(associated(tens_layout)) then
           do i=1,nspaces; write(dev,'(" ")',ADVANCE='NO'); enddo
           write(dev,'(" Tensor body layout      = ",i7)') tens_layout%layout
           do i=1,nspaces; write(dev,'(" ")',ADVANCE='NO'); enddo
           write(dev,'(" Tensor body data kind   = ",i7)') tens_layout%data_type
          else
           do i=1,nspaces; write(dev,'(" ")',ADVANCE='NO'); enddo
           write(dev,'(" Tensot body has no layout yet")')
          endif
          do i=1,nspaces; write(dev,'(" ")',ADVANCE='NO'); enddo
          write(dev,'("}")')
         else
          write(dev,'("TENSOR BODY{")')
          write(dev,'(" Number of subtensors    = ",i7)') this%num_subtensors
          write(dev,'(" Tensor body value state = ",i7)') this%state
          tens_layout=>this%get_layout()
          if(associated(tens_layout)) then
           write(dev,'(" Tensor body layout      = ",i7)') tens_layout%layout
           write(dev,'(" Tensor body data kind   = ",i7)') tens_layout%data_type
          else
           write(dev,'(" Tensot body has no layout yet")')
          endif
          write(dev,'("}")')
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensBodyPrintIt
!--------------------------------------------------------
        subroutine TensBodyUpdateHeader(this,header,ierr)
!Updates the tensor header the body is associated with.
         implicit none
         class(tens_body_t), intent(inout):: this          !inout: tensor body
         class(tens_header_t), intent(in), target:: header !in: tensor header the body is associated with
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         if(allocated(this%layout)) then
          call this%layout%update_header(header,errc)
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensBodyUpdateHeader
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
         this%state=TEREC_BODY_UNDEF
         return
        end subroutine tens_body_dtor
![tens_rcrsv_t]=============================================================================================
        subroutine TensRcrsvCtorSigna(this,tens_name,subspaces,hspaces,ierr,dim_extent,dim_group,group_spec)
!Constructs a tensor by specifying a tensor signature and optionally a shape.
!See TensHeaderCtor for restrictions.
         implicit none
         class(tens_rcrsv_t), intent(out):: this              !out: tensor
         character(*), intent(in):: tens_name                 !in: alphanumeric_ tensor name
         integer(INTL), intent(in):: subspaces(1:)            !in: subspace multi-index (signature)
         integer(INTD), intent(in):: hspaces(1:)              !in: hierarchical vector space id for each tensor dimension
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTL), intent(in), optional:: dim_extent(1:) !in: dimension extents (those equal to 0 are unresolved)
         integer(INTD), intent(in), optional:: dim_group(1:)  !in: dimension restriction groups
         integer(INTD), intent(in), optional:: group_spec(1:) !in: restriction groups specification
         integer(INTD):: errc

         if(present(dim_extent)) then !signature + shape (possibly with unresolved dimensions, those equal to 0)
          if(present(dim_group)) then
           if(present(group_spec)) then
            call this%header%tens_header_ctor(errc,tens_name,subspaces,hspaces,dim_extent,dim_group,group_spec)
           else
            errc=TEREC_INVALID_ARGS
           endif
          else
           if(.not.present(group_spec)) then
            call this%header%tens_header_ctor(errc,tens_name,subspaces,hspaces,dim_extent)
           else
            errc=TEREC_INVALID_ARGS
           endif
          endif
         else !signature only
          if(.not.(present(dim_group).or.present(group_spec))) then
           call this%header%tens_header_ctor(errc,tens_name,subspaces,hspaces)
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
!-------------------------------------------------
        subroutine TensRcrsvCtorClone(this,source)
!Copy ctor.
         implicit none
         class(tens_rcrsv_t), intent(out):: this  !out: tensor
         class(tens_rcrsv_t), intent(in):: source !in: source tensor

         this%header=source%header
         this%body=source%body
         call this%body%update_header(this%header)
         return
        end subroutine TensRcrsvCtorClone
!-------------------------------------------------------
        subroutine TensRcrsvCtorUnpack(this,packet,ierr)
!Unpacks the object from a packet.
         implicit none
         class(tens_rcrsv_t), intent(out):: this     !out: tensor
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         class(tens_header_t), pointer:: tens_header_p

         call this%header%tens_header_ctor(packet,errc)
         if(errc.eq.PACK_SUCCESS) then
          tens_header_p=>this%get_header(errc)
          if(errc.eq.TEREC_SUCCESS) call this%body%tens_body_ctor(packet,errc,tens_header_p)
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensRcrsvCtorUnpack
!-------------------------------------------------
        subroutine TensRcrsvPack(this,packet,ierr)
!Packs the object into a packet.
         implicit none
         class(tens_rcrsv_t), intent(in):: this      !in: tensor
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         call this%header%pack(packet,errc)
         if(errc.eq.PACK_SUCCESS) call this%body%pack(packet,errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensRcrsvPack
!----------------------------------------------------------------------------------------------------------------
        function TensRcrsvIsSet(this,ierr,num_dims,shaped,unresolved,hspaced,layed,located,symmetric) result(res)
!Returns TRUE if the tensor is set, plus additional info.
         implicit none
         logical:: res                                     !out: result
         class(tens_rcrsv_t), intent(in):: this            !in: tensor
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD), intent(out), optional:: num_dims   !out: number of tensor dimensions
         logical, intent(out), optional:: shaped           !out: TRUE if tensor shape is set (even with unresolved dimensions)
         integer(INTD), intent(out), optional:: unresolved !out: number of unresolved tensor dimensions
         logical, intent(out), optional:: hspaced          !out: TRUE if the tensor dimensions are over hierarchical spaces
         logical, intent(out), optional:: layed            !out: TRUE if the tensor body storage layout is set
         logical, intent(out), optional:: located          !out: TRUE if the physical location for tensor body data is set
         logical, intent(out), optional:: symmetric        !out: TRUE if the tensor has symmetric dimensions, FALSE otherwise
         integer(INTD):: errc,nd,ng,unres
         logical:: shpd,hspc,layd,locd

         res=this%header%is_set(errc,num_dims=nd,num_groups=ng,shaped=shpd,unresolved=unres,hspaced=hspc)
         if(present(num_dims)) num_dims=nd
         if(present(shaped)) shaped=shpd
         if(present(unresolved)) unresolved=unres
         if(present(hspaced)) hspaced=hspc
         if(present(symmetric)) symmetric=(ng.gt.0)
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
!----------------------------------------------------------------
        subroutine TensRcrsvGetName(this,tens_name,name_len,ierr)
!Returns the alphanumeric_ name of the tensor.
         implicit none
         class(tens_rcrsv_t), intent(in):: this      !in: tensor
         character(*), intent(inout):: tens_name     !out: tensor name
         integer(INTD), intent(out):: name_len       !out: length of the tensor name
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         call this%header%get_name(tens_name,name_len,errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensRcrsvGetName
!--------------------------------------------------------
        function TensRcrsvGetRank(this,ierr) result(rank)
!Returns the rank of the tensor (number of dimensions).
         implicit none
         integer(INTD):: rank                        !out: tensor rank
         class(tens_rcrsv_t), intent(in):: this      !in: tensor
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         rank=this%header%get_rank(errc)
         if(present(ierr)) ierr=errc
         return
        end function TensRcrsvGetRank
!------------------------------------------------------------------------
        subroutine TensRcrsvGetSpec(this,subspaces,num_dims,ierr,hspaces)
!Returns the defining subspaces of the tensor (subspace multi-index).
         implicit none
         class(tens_rcrsv_t), intent(in):: this                    !in: tensor
         integer(INTL), intent(inout):: subspaces(1:)              !out: defining subspaces (their IDs)
         integer(INTD), intent(out):: num_dims                     !out: number of tensor dimensions
         integer(INTD), intent(out), optional:: ierr               !out: error code
         type(hspace_reg_t), intent(inout), optional:: hspaces(1:) !out: hierarchical vector space for each tensor dimension
         integer(INTD):: errc

         if(present(hspaces)) then
          call this%header%get_spec(subspaces,num_dims,errc,hspaces)
         else
          call this%header%get_spec(subspaces,num_dims,errc)
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensRcrsvGetSpec
!-----------------------------------------------------------
        subroutine TensRcrsvGetDims(this,dims,num_dims,ierr)
!Returns tensor dimension extents together with the tensor rank.
         implicit none
         class(tens_rcrsv_t), intent(in):: this      !in: tensor
         integer(INTL), intent(inout):: dims(1:)     !out: tensor dimension extents
         integer(INTD), intent(out):: num_dims       !out: number of tensor dimensions
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         call this%header%get_dims(dims,num_dims,errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensRcrsvGetDims
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
         class(DataDescr_t), intent(in):: data_descr !in: DDSS data descriptor for tensor body (will be cloned)
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
!----------------------------------------------------
        subroutine TensRcrsvUpdate(this,another,ierr)
!Updates the tensor information (new resolution -> new layout -> new body).
!<this> tensor must at least be defined (have its signature set).
!<another> tensor must have the same signature, and possibly, shape.
!Then, any additional information will be transferred from <another>
!tensor to <this> tensor, thus updating it.
         implicit none
         class(tens_rcrsv_t), intent(inout):: this   !inout: tensor being updated
         class(tens_rcrsv_t), intent(in):: another   !in: same tensor with more complete information (later stage of definition)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,this_nd,this_unres,an_nd,an_unres
         logical:: this_shp,this_lay,this_loc,an_shp,an_lay,an_loc

         if(this%is_set(errc,num_dims=this_nd,shaped=this_shp,unresolved=this_unres,layed=this_lay,located=this_loc)) then
          if(errc.eq.TEREC_SUCCESS) then
           if(another%is_set(errc,num_dims=an_nd,shaped=an_shp,unresolved=an_unres,layed=an_lay,located=an_loc)) then
            if(errc.eq.TEREC_SUCCESS) then
             if(this%header%signature%compare(another%header%signature).eq.CMP_EQ) then
 !Update shape:
              if(an_shp) then
               if(.not.this_shp) then !copy shape from <another>
                this%header%shape=another%header%shape
               else
                if(this%header%shape%compare(another%header%shape).ne.CMP_EQ) errc=TEREC_INVALID_REQUEST
               endif
              endif
 !Update layout (copy the full body info):
              if(an_lay.and.errc.eq.TEREC_SUCCESS) then
               if(.not.this_lay) then
                this%body=another%body
                call this%body%update_header(this%header,errc)
               else
 !Update location (actually copy the full body info):
                if(an_loc.and.(.not.this_loc)) then
                 this%body=another%body
                 call this%body%update_header(this%header,errc)
                endif
               endif
              endif
             else
              errc=TEREC_INVALID_REQUEST
             endif
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
        end subroutine TensRcrsvUpdate
!-----------------------------------------------------------
        subroutine TensRcrsvResetState(this,ierr,body_state)
!Resets the tensor body value state.
         implicit none
         class(tens_rcrsv_t), intent(inout):: this        !inout: tensor
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD), intent(in), optional:: body_state !in: tensor body value state
         integer(INTD):: errc

         if(present(body_state)) then
          call this%body%reset_state(errc,body_state)
         else
          call this%body%reset_state(errc)
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensRcrsvResetState
!---------------------------------------------------------------
        function TensRcrsvGetState(this,ierr) result(body_state)
!Returns the tensor body value state.
         implicit none
         integer(INTD):: body_state                  !out: tensor body value state
         class(tens_rcrsv_t), intent(in):: this      !in: tensor
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         body_state=this%body%get_state(errc)
         if(present(ierr)) ierr=errc
         return
        end function TensRcrsvGetState
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
!--------------------------------------------------------------------------------------------
        function TensRcrsvGetDescriptor(this,ierr,skip_body,skip_location) result(tens_descr)
!Returns a tensor descriptor <tens_descr_t> object uniquely characterizing
!the tensor signature, shape, layout kind, data type, and location.
!Essentially, this function is an indirect CTOR for <tens_descr_t>.
!tens_descr.info(:) format:
! 1. hspace_idx(1:rank): space id for each tensor dimension;
! 2. subspace_idx(1:rank): subspace id for each tensor dimension;
! 3. dim_extent(1:rank): extent for each tensor dimension;
! 4. dim_group(1:rank): previous dimension for each (symmetric) tensor dimension (0:none);
! 5. dim_group_restriction(1:rank): group restriction for each symmetric tensor dimension;
! 6. layout kind;
! 7. data type;
! 8. body size in bytes;
! 9. location (DDSS process id): Skipped when <skip_location>=TRUE.
!TOTAL size = 5*rank + 1 + 1 + 1 + 1 = 5*rank + 4 [elements]
!`Note: This function violates object encapsulation by
! directly accessing data members of the data members
! of the <tens_rcrsv_t> class. However, since it is a read-only
! access and all accessed data members are defined in the same
! module, it should not cause a problem.
         implicit none
         type(tens_descr_t):: tens_descr               !out: tensor descriptor
         class(tens_rcrsv_t), intent(in):: this        !in: tensor
         integer(INTD), intent(out), optional:: ierr   !out: error code
         logical, intent(in), optional:: skip_body     !in: if TRUE, the tensor body info (layout, location) will be omitted (defaults to FALSE)
         logical, intent(in), optional:: skip_location !in: if TRUE, the tensor location will be omitted (defaults to FALSE)
         integer(INTD):: errc,num_dims,unresolved,i,j,k,ng,gres
         logical:: shaped,hspaced,layed,located,symmetric,skiploc,skipbody
         class(DataDescr_t), pointer:: descr_p
         integer(INTD):: grp_pos(0:MAX_TENSOR_RANK*2) !the factor of 2 is just used to make it big enough
         !real(8):: tms

         !tms=thread_wtime()
         skiploc=.FALSE.; if(present(skip_location)) skiploc=skip_location
         skipbody=.FALSE.; if(present(skip_body)) then; skipbody=skip_body; skiploc=((.not.skipbody).and.skiploc); endif
         if(this%is_set(errc,num_dims,shaped,unresolved,hspaced,layed,located,symmetric)) then
          if(errc.eq.TEREC_SUCCESS) then
           if(unresolved.eq.0.and.shaped.and.(layed.or.skipbody).and.(located.or.skiploc)) then !sufficiently defined tensor
            tens_descr%rank=num_dims
            tens_descr%char_name=this%header%signature%char_name
            allocate(tens_descr%info(num_dims*5+4)) !see the format right above
            i=0
            if(num_dims.gt.0) then
             tens_descr%info(i+1:i+num_dims)=this%header%signature%hspace(1:num_dims)%space_id; i=i+num_dims
             tens_descr%info(i+1:i+num_dims)=this%header%signature%space_idx(1:num_dims); i=i+num_dims
             tens_descr%info(i+1:i+num_dims)=this%header%shape%dim_extent(1:num_dims); i=i+num_dims
             if(symmetric) then !dimension symmetry may be present
              ng=this%header%num_groups(errc) !number of non-trivial symmetric dimension groups
              if(errc.eq.TEREC_SUCCESS) then
               if(ng.gt.0) then !dimension symmetry indeed present
                if(ng.le.size(grp_pos)) then
                 grp_pos(0:ng)=0
                 do j=1,num_dims
                  k=this%header%get_dim_group(j,gres,errc); if(errc.ne.TEREC_SUCCESS) exit
                  tens_descr%info(i+j)=grp_pos(k); tens_descr%info(i+num_dims+j)=gres
                  if(k.gt.0) grp_pos(k)=j
                 enddo
                 i=i+num_dims*2
                else
                 errc=TEREC_ERROR
                endif
               else !no dimension symmetry
                tens_descr%info(i+1:i+num_dims)=0; i=i+num_dims
                tens_descr%info(i+1:i+num_dims)=TEREC_IND_RESTR_NONE; i=i+num_dims
               endif
              endif
             else !no dimension symmetry
              tens_descr%info(i+1:i+num_dims)=0; i=i+num_dims
              tens_descr%info(i+1:i+num_dims)=TEREC_IND_RESTR_NONE; i=i+num_dims
             endif
            endif
            if(errc.eq.TEREC_SUCCESS) then
             if(.not.skipbody) then
              i=i+1; tens_descr%info(i)=this%body%layout%get_layout_kind()
              i=i+1; tens_descr%info(i)=this%body%layout%get_data_type()
              i=i+1; tens_descr%info(i)=this%body%layout%get_body_size()
              if(.not.skiploc) then
               descr_p=>this%body%layout%get_data_descr(errc)
               if(errc.eq.TEREC_SUCCESS) then
                if(descr_p%is_set(errc,proc_rank=j)) then
                 if(errc.eq.0) then
                  i=i+1; tens_descr%info(i)=j
                 else
                  errc=TEREC_OBJ_CORRUPTED
                 endif
                else
                 errc=TEREC_OBJ_CORRUPTED
                endif
               endif
               descr_p=>NULL()
              else
               tens_descr%info(i+1:)=-1
              endif
             else
              tens_descr%info(i+1:)=-1
             endif
            endif
           else
            errc=TEREC_INVALID_ARGS
           endif
          endif
         else
          if(errc.eq.TEREC_SUCCESS) errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         !write(CONS_OUT,*)'#MSG: tens_rcrsv_t.get_descriptor() timing: ',thread_wtime(tms) !timing
         return
        end function TensRcrsvGetDescriptor
!-------------------------------------------------------------------------------------------------
        subroutine TensRcrsvSplitList(this,split_dims,subtensors,ierr,num_subtensors,headers_only)
!Splits the given tensor into subtensors and appends those to a list of subtensors, either as
!tensors or as headers only. If the input parental tensor is shaped with concrete dimension extents,
!the children subtensors in contrast will not carry concrete dimension extents, but deferred dimension
!extents instead. However, they will obey the parental dimension restriction grouping with the following
!assumption: Restricted tensor dimensions belonging to the same group are supposed to depend on
!each other from right to left, e.g., in {T(a,b,c,d,e):[a<c<e],[b<d]}, a dependent index
!on the right always depends on some of the previous indices on the left. The first
!restricted index on the left in each group is actually independent (e.g., "a" and "b" above).
         implicit none
         class(tens_rcrsv_t), intent(in):: this                !in: parental tensor (either shaped or unshaped)
         integer(INTD), intent(in):: split_dims(1:)            !in: tensor dimensions to be split (at least one)
         type(list_bi_t), intent(inout):: subtensors           !out: list of subtensors
         integer(INTD), intent(out), optional:: ierr           !out: error code
         integer(INTD), intent(out), optional:: num_subtensors !out: number of subtensors generated
         logical, intent(in), optional:: headers_only          !in: if TRUE, the output list will contain tensor headers only (not full tensors)
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
         type(hspace_reg_t):: hsreg(1:MAX_TENSOR_RANK) !hierarchical space for each tensor dimension
         type(list_iter_t):: lit
         type(tens_header_t):: thead
         type(tens_rcrsv_t):: subtens
         logical:: shpd,hspc,head_only

         nsubt=0 !number of generated subtensors
         head_only=.FALSE.; if(present(headers_only)) head_only=headers_only
         if(this%is_set(errc,shaped=shpd,hspaced=hspc)) then !shaped tensors over hierarchical vector spaces are expected
          if(errc.eq.TEREC_SUCCESS) then
           call this%header%get_name(tens_name,tnl,errc)
           if(errc.eq.TEREC_SUCCESS) then
            call this%header%get_spec(sidx,nd,errc,hsreg) !nd: total number of tensor dimensions; sidx(1:nd): subspace id's
            if(errc.eq.TEREC_SUCCESS) then
             errc=lit%init(subtensors)
             if(errc.eq.GFC_SUCCESS) then
              if(nd.gt.0) then !true tensor over hierarchical vector spaces
               sd=size(split_dims) !sd: number of tensor dimensions to split
               if(sd.gt.0.and.sd.le.nd) then !true splitting
                if(hspc) then
                 call extract_subspaces_to_sbuf(errc) !sbuf(1:nsb),firo(1:nd),swid(1:nd)
                 if(errc.eq.TEREC_SUCCESS) then
                  call setup_index_dependencies(errc)
                  if(errc.eq.TEREC_SUCCESS) then
                   midx(1:nd)=-1; n=2 !n=2 is a special setting to start iterating midx(:)
                   do
                    call get_next_midx(n,errc); if(errc.ne.TEREC_SUCCESS.or.n.le.0) exit
                    call construct_subtensor_header(thead,errc); if(errc.ne.TEREC_SUCCESS) exit
                    if(head_only) then
                     errc=lit%append(thead)
                    else
                     call subtens%tens_rcrsv_ctor(thead,errc)
                     if(errc.eq.TEREC_SUCCESS) errc=lit%append(subtens)
                    endif
                    if(errc.eq.GFC_SUCCESS) then; nsubt=nsubt+1; else; exit; endif
                   enddo
                   if(errc.ne.TEREC_SUCCESS) errc=TEREC_UNABLE_COMPLETE
                  endif
                 endif
                else
                 errc=TEREC_INVALID_REQUEST
                endif
               elseif(sd.eq.0) then !no splitting, return the original header
                if(head_only) then
                 errc=lit%append(this%header)
                else
                 call subtens%tens_rcrsv_ctor(this%header,errc)
                 if(errc.eq.TEREC_SUCCESS) errc=lit%append(subtens)
                endif
                if(errc.eq.GFC_SUCCESS) then; nsubt=nsubt+1; else; errc=TEREC_UNABLE_COMPLETE; endif
               else
                errc=TEREC_INVALID_ARGS
               endif
              elseif(nd.eq.0) then !scalar
               if(head_only) then
                errc=lit%append(this%header)
               else
                call subtens%tens_rcrsv_ctor(this%header,errc)
                if(errc.eq.TEREC_SUCCESS) errc=lit%append(subtens)
               endif
               if(errc.eq.GFC_SUCCESS) then; nsubt=nsubt+1; else; errc=TEREC_UNABLE_COMPLETE; endif
              else
               errc=TEREC_OBJ_CORRUPTED
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
           integer(INTD):: jj,js,jd,je
           type(vec_tree_iter_t):: vt_it(1:MAX_TENSOR_RANK)
           class(*), pointer:: jup
           class(subspace_t), pointer:: jssp
 !Init:
           do jj=1,sd
            jd=split_dims(jj)
            jerr=vt_it(jd)%init(hsreg(jd)%hspace_p%get_aggr_tree(nsb))
            if(jerr.ne.GFC_SUCCESS.or.nsb.ne.0) exit
           enddo
 !Count:
           if(jerr.eq.GFC_SUCCESS.and.nsb.eq.0) then
            nsb=0; firo(0:nd)=0; swid(0:nd)=1
            do jj=1,sd !loop over the dimensions to split
             jd=split_dims(jj)
             jerr=vt_it(jd)%move_to(sidx(jd)); if(jerr.ne.GFC_SUCCESS) exit
             js=vt_it(jd)%get_num_children(jerr); if(jerr.ne.GFC_SUCCESS) exit
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
                jerr=vt_it(jd)%move_to(sidx(jd)); if(jerr.ne.GFC_SUCCESS) exit sloop
                if(js.gt.0) then
                 jerr=vt_it(jd)%move_to_child()
                 do while(js.gt.0.and.jerr.eq.GFC_SUCCESS)
                  jup=>vt_it(jd)%get_value(jerr); if(jerr.ne.GFC_SUCCESS) exit
                  select type(jup); class is(subspace_t); jssp=>jup; end select
                  sbuf(jj)=jssp%get_id(jerr); if(jerr.ne.0) exit
                  jj=jj+1; js=js-1; if(js.gt.0) jerr=vt_it(jd)%move_to_sibling()
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
            do jj=1,sd
             jd=split_dims(jj)
             je=vt_it(jd)%release(); if(je.ne.GFC_SUCCESS.and.jerr.eq.TEREC_SUCCESS) jerr=TEREC_ERROR
            enddo
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
             jcmp=hsreg(np)%hspace_p%compare_subranges(js1,js2,jerr)
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
              jcmp=hsreg(jd)%hspace_p%compare_subranges(js1,js2,jerr); if(jerr.ne.0) then; jerr=TEREC_ERROR; exit; endif
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
             call thd%tens_header_ctor(jerr,tens_name(1:tnl),sidx(1:nd),hsreg(1:nd)%space_id,&
                                     &(/(0_INTL,jd=1,nd)/),dim_group(1:nd),group_spec(1:jgr))
            endif
           else
  !Construct a new subtensor without index restrictions:
            call thd%tens_header_ctor(jerr,tens_name(1:tnl),sidx(1:nd),hsreg(1:nd)%space_id,(/(0_INTL,jd=1,nd)/))
           endif
           return
          end subroutine construct_subtensor_header

        end subroutine TensRcrsvSplitList
!---------------------------------------------------------------------------------------------------
        subroutine TensRcrsvSplitVector(this,split_dims,subtensors,ierr,num_subtensors,headers_only)
!Splits the given tensor into subtensors and appends those to a vector of subtensors, either as
!tensors or as headers only. If the input parental tensor is shaped with concrete dimension extents,
!the children subtensors in contrast will not carry concrete dimension extents, but deferred dimension
!extents instead. However, they will obey the parental dimension restriction grouping with the following
!assumption: Restricted tensor dimensions belonging to the same group are supposed to depend on
!each other from right to left, e.g., in {T(a,b,c,d,e):[a<c<e],[b<d]}, a dependent index
!on the right always depends on some of the previous indices on the left. The first
!restricted index on the left in each group is actually independent (e.g., "a" and "b" above).
         implicit none
         class(tens_rcrsv_t), intent(in):: this                !in: parental tensor (either shaped or unshaped)
         integer(INTD), intent(in):: split_dims(1:)            !in: tensor dimensions to be split (at least one)
         type(vector_t), intent(inout):: subtensors            !out: vector of subtensors
         integer(INTD), intent(out), optional:: ierr           !out: error code
         integer(INTD), intent(out), optional:: num_subtensors !out: number of subtensors generated
         logical, intent(in), optional:: headers_only          !in: if TRUE, the output list will contain tensor headers only (not full tensors)
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
         type(hspace_reg_t):: hsreg(1:MAX_TENSOR_RANK) !hierarchical space for each tensor dimension
         type(vector_iter_t):: vit
         type(tens_header_t):: thead
         type(tens_rcrsv_t):: subtens
         logical:: shpd,hspc,head_only

         nsubt=0 !number of generated subtensors
         head_only=.FALSE.; if(present(headers_only)) head_only=headers_only
         if(this%is_set(errc,shaped=shpd,hspaced=hspc)) then !shaped tensors over hierarchical vector spaces are expected
          if(errc.eq.TEREC_SUCCESS) then
           call this%header%get_name(tens_name,tnl,errc)
           if(errc.eq.TEREC_SUCCESS) then
            call this%header%get_spec(sidx,nd,errc,hsreg) !nd: total number of tensor dimensions; sidx(1:nd): subspace id's
            if(errc.eq.TEREC_SUCCESS) then
             errc=vit%init(subtensors)
             if(errc.eq.GFC_SUCCESS) then
              if(nd.gt.0) then !true tensor over hierarchical vector spaces
               sd=size(split_dims) !sd: number of tensor dimensions to split
               if(sd.gt.0.and.sd.le.nd) then !true splitting
                if(hspc) then
                 call extract_subspaces_to_sbuf(errc) !sbuf(1:nsb),firo(1:nd),swid(1:nd)
                 if(errc.eq.TEREC_SUCCESS) then
                  call setup_index_dependencies(errc)
                  if(errc.eq.TEREC_SUCCESS) then
                   midx(1:nd)=-1; n=2 !n=2 is a special setting to start iterating midx(:)
                   do
                    call get_next_midx(n,errc); if(errc.ne.TEREC_SUCCESS.or.n.le.0) exit
                    call construct_subtensor_header(thead,errc); if(errc.ne.TEREC_SUCCESS) exit
                    if(head_only) then
                     errc=vit%append(thead)
                    else
                     call subtens%tens_rcrsv_ctor(thead,errc)
                     if(errc.eq.TEREC_SUCCESS) errc=vit%append(subtens)
                    endif
                    if(errc.eq.GFC_SUCCESS) then; nsubt=nsubt+1; else; exit; endif
                   enddo
                   if(errc.ne.TEREC_SUCCESS) errc=TEREC_UNABLE_COMPLETE
                  endif
                 endif
                else
                 errc=TEREC_INVALID_REQUEST
                endif
               elseif(sd.eq.0) then !no splitting, return the original header
                if(head_only) then
                 errc=vit%append(this%header)
                else
                 call subtens%tens_rcrsv_ctor(this%header,errc)
                 if(errc.eq.TEREC_SUCCESS) errc=vit%append(subtens)
                endif
                if(errc.eq.GFC_SUCCESS) then; nsubt=nsubt+1; else; errc=TEREC_UNABLE_COMPLETE; endif
               else
                errc=TEREC_INVALID_ARGS
               endif
              elseif(nd.eq.0) then !scalar
               if(head_only) then
                errc=vit%append(this%header)
               else
                call subtens%tens_rcrsv_ctor(this%header,errc)
                if(errc.eq.TEREC_SUCCESS) errc=vit%append(subtens)
               endif
               if(errc.eq.GFC_SUCCESS) then; nsubt=nsubt+1; else; errc=TEREC_UNABLE_COMPLETE; endif
              else
               errc=TEREC_OBJ_CORRUPTED
              endif
              i=vit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.TEREC_SUCCESS) errc=TEREC_ERROR
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
           integer(INTD):: jj,js,jd,je
           type(vec_tree_iter_t):: vt_it(1:MAX_TENSOR_RANK)
           class(*), pointer:: jup
           class(subspace_t), pointer:: jssp
 !Init:
           do jj=1,sd
            jd=split_dims(jj)
            jerr=vt_it(jd)%init(hsreg(jd)%hspace_p%get_aggr_tree(nsb))
            if(jerr.ne.GFC_SUCCESS.or.nsb.ne.0) exit
           enddo
 !Count:
           if(jerr.eq.GFC_SUCCESS.and.nsb.eq.0) then
            nsb=0; firo(0:nd)=0; swid(0:nd)=1
            do jj=1,sd !loop over the dimensions to split
             jd=split_dims(jj)
             jerr=vt_it(jd)%move_to(sidx(jd)); if(jerr.ne.GFC_SUCCESS) exit
             js=vt_it(jd)%get_num_children(jerr); if(jerr.ne.GFC_SUCCESS) exit
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
                jerr=vt_it(jd)%move_to(sidx(jd)); if(jerr.ne.GFC_SUCCESS) exit sloop
                if(js.gt.0) then
                 jerr=vt_it(jd)%move_to_child()
                 do while(js.gt.0.and.jerr.eq.GFC_SUCCESS)
                  jup=>vt_it(jd)%get_value(jerr); if(jerr.ne.GFC_SUCCESS) exit
                  select type(jup); class is(subspace_t); jssp=>jup; end select
                  sbuf(jj)=jssp%get_id(jerr); if(jerr.ne.0) exit
                  jj=jj+1; js=js-1; if(js.gt.0) jerr=vt_it(jd)%move_to_sibling()
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
            do jj=1,sd
             jd=split_dims(jj)
             je=vt_it(jd)%release(); if(je.ne.GFC_SUCCESS.and.jerr.eq.TEREC_SUCCESS) jerr=TEREC_ERROR
            enddo
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
             jcmp=hsreg(np)%hspace_p%compare_subranges(js1,js2,jerr)
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
              jcmp=hsreg(jd)%hspace_p%compare_subranges(js1,js2,jerr); if(jerr.ne.0) then; jerr=TEREC_ERROR; exit; endif
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
             call thd%tens_header_ctor(jerr,tens_name(1:tnl),sidx(1:nd),hsreg(1:nd)%space_id,&
                                     &(/(0_INTL,jd=1,nd)/),dim_group(1:nd),group_spec(1:jgr))
            endif
           else
  !Construct a new subtensor without index restrictions:
            call thd%tens_header_ctor(jerr,tens_name(1:tnl),sidx(1:nd),hsreg(1:nd)%space_id,(/(0_INTL,jd=1,nd)/))
           endif
           return
          end subroutine construct_subtensor_header

        end subroutine TensRcrsvSplitVector
!------------------------------------------------------------
        subroutine TensRcrsvPrintIt(this,ierr,dev_id,nspaces)
!Prints the tensor info.
         implicit none
         class(tens_rcrsv_t), intent(in):: this        !in: tensor
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD), intent(in), optional:: dev_id  !in: output device (defaults to screen)
         integer(INTD), intent(in), optional:: nspaces !in: number of leading spaces
         integer(INTD):: errc,devo,nsp,i

         devo=6; if(present(dev_id)) devo=dev_id
         nsp=0; if(present(nspaces)) nsp=nspaces
         do i=1,nsp; write(devo,'(" ")',ADVANCE='NO'); enddo
         write(devo,'("TENSOR{")')
         call this%header%print_it(errc,devo,nsp+1)
         call this%body%print_it(errc,devo,nsp+1)
         do i=1,nsp; write(devo,'(" ")',ADVANCE='NO'); enddo
         write(devo,'("}")')
         if(present(ierr)) ierr=errc
         return
        end subroutine TensRcrsvPrintIt
!---------------------------------------
        subroutine tens_rcrsv_dtor(this)
         implicit none
         type(tens_rcrsv_t):: this

         return
        end subroutine tens_rcrsv_dtor
![tens_descr_t]======================================
        subroutine TensDescrPrintIt(this,ierr,dev_id)
         implicit none
         class(tens_descr_t), intent(in):: this        !in: tensor descriptor
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD), intent(in), optional:: dev_id  !in: output device (defaults to screen)
         integer(INTD):: errc,devo,i,l,n

         errc=TEREC_SUCCESS
         devo=6; if(present(dev_id)) devo=dev_id
         write(devo,'("TENSOR_DESCRIPTOR{")')
         write(devo,'(1x,i4)') this%rank
         write(devo,*) this%char_name
         l=0; n=this%rank
         do i=1,5
          write(devo,'(64(1x,i6))') this%info(l+1:l+n); l=l+n
         enddo
         do i=1,4
          l=l+1; write(devo,'(1x,i11)') this%info(l)
         enddo
         write(devo,'("}")')
         if(present(ierr)) ierr=errc
         return
        end subroutine TensDescrPrintIt
!----------------------------------------------------------
        function TensDescrCompare(this,another) result(cmp)
         implicit none
         integer(INTD):: cmp                               !out: comparison result: {CMP_EQ,CMP_LT,CMP_GT,CMP_ER}
         class(tens_descr_t), intent(in), target:: this    !in: tensor descriptor 1
         class(tens_descr_t), intent(in), target:: another !in: tensor descriptor 2
         integer(INTD):: l1,l2,i1,i2

         cmp=CMP_EQ
         if(this%rank.lt.another%rank) then
          cmp=CMP_LT
         elseif(this%rank.gt.another%rank) then
          cmp=CMP_GT
         else
          l1=len(this%char_name); l2=len(another%char_name)
          if(l1.lt.l2) then
           cmp=CMP_LT
          elseif(l1.gt.l2) then
           cmp=CMP_GT
          else
           do l2=1,l1
            i1=iachar(this%char_name(l2:l2)); i2=iachar(another%char_name(l2:l2))
            if(i1.lt.i2) then; cmp=CMP_LT; exit; elseif(i1.gt.i2) then; cmp=CMP_GT; exit; endif
           enddo
           if(cmp.eq.CMP_EQ) then
            l1=size(this%info)
            if(size(another%info).eq.l1) then !trap
             do l2=1,l1
              if(this%info(l2).lt.another%info(l2)) then
               cmp=CMP_LT; exit
              elseif(this%info(l2).gt.another%info(l2)) then
               cmp=CMP_GT; exit
              endif
             enddo
            else
             cmp=CMP_ER
            endif
           endif
          endif
         endif
         return
        end function TensDescrCompare
!---------------------------------------
        subroutine tens_descr_dtor(this)
         implicit none
         type(tens_descr_t):: this

         if(allocated(this%info)) deallocate(this%info)
         if(allocated(this%char_name)) deallocate(this%char_name)
         this%rank=-1
         return
        end subroutine tens_descr_dtor
![tens_argument_t]========================================
        subroutine TensArgumentSetTensor(this,tensor,ierr)
!Sets the tensor argument by pointer association.
         implicit none
         class(tens_argument_t), intent(inout):: this     !inout: tensor argument
         class(tens_rcrsv_t), intent(in), target:: tensor !in: tensor target (tens_rcrsv_t)
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD):: errc

         if(.not.this%alloc) then
          if(tensor%is_set(errc)) then
           if(errc.eq.TEREC_SUCCESS) this%tens_p=>tensor
          else
           if(errc.eq.TEREC_SUCCESS) errc=TEREC_INVALID_ARGS
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensArgumentSetTensor
!-------------------------------------------------------
        subroutine TensArgumentAllocateTensor(this,ierr)
!Allocates an empty tensor for a subsequent definition.
         implicit none
         class(tens_argument_t), intent(inout):: this !inout: tensor argument
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         if(.not.this%alloc) then
          allocate(this%tens_p,STAT=errc)
          if(errc.eq.0) then; this%alloc=.TRUE.; else; errc=TEREC_MEM_ALLOC_FAILED; endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensArgumentAllocateTensor
!------------------------------------------------------------------------------
        function TensArgumentIsSet(this,ierr,num_dims,tens_p,alloc) result(ans)
!Returns TRUE if the tensor argument is set, plus additional info.
         implicit none
         logical:: ans                                   !out: answer
         class(tens_argument_t), intent(in):: this       !in: tensor argument
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD), intent(out), optional:: num_dims !out: tensor rank
         class(tens_rcrsv_t), pointer, intent(inout), optional:: tens_p !out: pointer to the tensor
         logical, intent(out), optional:: alloc          !out: allocation status of the tensor pointer (TRUE:allocated; FALSE:associated)
         integer(INTD):: errc
         logical:: alcd

         errc=TEREC_SUCCESS; ans=associated(this%tens_p); alcd=this%alloc
         if(ans.and.present(num_dims)) then
          if(.not.this%tens_p%is_set(errc,num_dims=num_dims)) then
           if(errc.eq.TEREC_SUCCESS) errc=TEREC_OBJ_CORRUPTED
          endif
         endif
         if(present(tens_p)) tens_p=>this%tens_p
         if(present(alloc)) alloc=alcd
         if(present(ierr)) ierr=errc
         return
        end function TensArgumentIsSet
!---------------------------------------------------
        subroutine TensArgumentFreeTensor(this,ierr)
!Frees the tensor (deallocation or dissociation).
         implicit none
         class(tens_argument_t), intent(inout):: this !inout: tensor argument
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         if(associated(this%tens_p)) then
          if(this%alloc) then
           deallocate(this%tens_p,STAT=errc); if(errc.ne.0) errc=TEREC_MEM_FREE_FAILED
           this%alloc=.FALSE.
          endif
          this%tens_p=>NULL()
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensArgumentFreeTensor
!------------------------------------------
        subroutine tens_argument_dtor(this)
         implicit none
         type(tens_argument_t):: this
         integer(INTD):: errc

         call this%free_tensor(errc)
         return
        end subroutine tens_argument_dtor
![tens_operation_t]=============================
        subroutine TensOperationClean(this,ierr)
!Cleans the tensor operation.
         implicit none
         class(tens_operation_t), intent(inout):: this !out: empty (clean)tensor operation
         integer(INTD), intent(out), optional:: ierr   !out: error code

         this%num_args=0
         if(present(ierr)) ierr=TEREC_SUCCESS
         return
        end subroutine TensOperationClean
!------------------------------------------------------------
        subroutine TensOperationSetArgument(this,tensor,ierr)
!Sets the next tensor operation argument.
         implicit none
         class(tens_operation_t), intent(inout):: this    !inout: tensor operation
         class(tens_rcrsv_t), intent(in), target:: tensor !in: tensor to be set as an argument
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD):: errc

         if(.not.this%args_full(errc)) then
          if(errc.eq.TEREC_SUCCESS) then
           if(this%num_args.lt.MAX_TENSOR_OPERANDS) then
            call this%tens_arg(this%num_args)%set_tensor(tensor,errc)
            if(errc.eq.TEREC_SUCCESS) this%num_args=this%num_args+1
           else
            errc=TEREC_INVALID_REQUEST
           endif
          endif
         else
          if(errc.eq.TEREC_SUCCESS) errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOperationSetArgument
!----------------------------------------------------------------------
        subroutine TensOperationResetArgument(this,tensor,arg_num,ierr)
!Resets an already set tensor argument.
         implicit none
         class(tens_operation_t), intent(inout):: this    !inout: tensor operation
         class(tens_rcrsv_t), intent(in), target:: tensor !in: tensor to be set as the new argument
         integer(INTD), intent(in):: arg_num              !in: argument number: [0..max]
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         if(arg_num.ge.0.and.arg_num.lt.this%num_args) then
          call this%tens_arg(arg_num)%set_tensor(tensor,errc)
         else
          errc=TEREC_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOperationResetArgument
!-------------------------------------------------------------------
        function TensOperationGetNumArgs(this,ierr) result(num_args)
!Returns the current number of set tensor arguments.
         implicit none
         integer(INTD):: num_args                    !out: number of tensor arguments set
         class(tens_operation_t), intent(in):: this  !in: tensor operation
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS; num_args=this%num_args
         if(num_args.lt.0.or.num_args.gt.MAX_TENSOR_OPERANDS) errc=TEREC_OBJ_CORRUPTED
         if(present(ierr)) ierr=errc
         return
        end function TensOperationGetNumArgs
!--------------------------------------------------------------------------
        function TensOperationGetArgument(this,arg_num,ierr) result(tens_p)
!Returns a pointer to the specific tensor argument.
         implicit none
         class(tens_rcrsv_t), pointer:: tens_p       !out: pointer to the tensor argument (tens_rcrsv_t)
         class(tens_operation_t), intent(in):: this  !in: tensor operation
         integer(INTD), intent(in):: arg_num         !in: argument number: [0..max]
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS; tens_p=>NULL()
         if(arg_num.ge.0.and.arg_num.lt.this%num_args) then
          tens_p=>this%tens_arg(arg_num)%tens_p
         else
          errc=TEREC_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensOperationGetArgument
!----------------------------------------------------------
        subroutine TensOperationAllocateArgument(this,ierr)
!Allocates the next argument in a tensor operation.
         implicit none
         class(tens_operation_t), intent(inout):: this !inout: tensor operation
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         if(this%num_args.lt.MAX_TENSOR_OPERANDS) then
          call this%tens_arg(this%num_args)%allocate_tensor(errc)
         else
          errc=TEREC_UNABLE_COMPLETE
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOperationAllocateArgument
!-------------------------------------------------------
        subroutine TensOperationFreeArguments(this,ierr)
!Frees all arguments in the tensor operation.
         implicit none
         class(tens_operation_t), intent(inout):: this !inout: tensor operation
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: errc,ier

         errc=TEREC_SUCCESS
         do while(this%num_args.gt.0)
          this%num_args=this%num_args-1
          call this%tens_arg(this%num_args)%free_tensor(ier)
          if(ier.ne.TEREC_SUCCESS.and.errc.eq.TEREC_SUCCESS) errc=ier
         enddo
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOperationFreeArguments
![permutation_t]=====================================
        subroutine PermutationReset(this,length,ierr)
!Resets the permutation (ctor).
         implicit none
         class(permutation_t), intent(inout):: this  !inout: permutation
         integer(INTD), intent(in):: length          !in: new length of the permutation
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         if(length.gt.0) then
          if(allocated(this%prm)) then
           if(size(this%prm).lt.1+length) then
            deallocate(this%prm)
            allocate(this%prm(0:length),STAT=errc)
           endif
          else
           allocate(this%prm(0:length),STAT=errc)
          endif
          if(errc.eq.0) then
           this%length=length; this%prm(0)=0
          else
           this%length=0; errc=TEREC_MEM_ALLOC_FAILED
          endif
         elseif(length.eq.0) then
          if(allocated(this%prm)) deallocate(this%prm)
          this%length=0
         else
          errc=TEREC_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine PermutationReset
!-------------------------------------------------------------------------------
        function PermutationGetAccess(this,length,ierr,with_sign) result(perm_p)
!Returns a pointer to the permutation body, with or without the sign.
         implicit none
         integer(INTD), pointer:: perm_p(:)              !out: ponter to the permutation body
         class(permutation_t), intent(in), target:: this !in: permutation
         integer(INTD), intent(out):: length             !out: length of the permutation
         integer(INTD), intent(out), optional:: ierr     !out: error code
         logical, intent(in), optional:: with_sign       !in: with or without sign (default is without)
         integer(INTD):: errc
         logical:: ws

         errc=TEREC_SUCCESS; perm_p=>NULL(); length=0
         if(this%length.gt.0) then
          ws=.FALSE.; if(present(with_sign)) ws=with_sign
          if(ws) then
           perm_p(0:this%length)=>this%prm(0:this%length) !perm_p(0) is the sign, perm_p(1:length) is the permutation itself
          else
           perm_p(1:this%length)=>this%prm(1:this%length)
          endif
          length=this%length
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function PermutationGetAccess
!---------------------------------------------------------------
        function PermutationGetSign(this,ierr) result(perm_sign)
!Returns the currently set sign of the permutation. If the sign
!has not been set previously, it will be recomputed, which
!is O(NlogN) operation.
         implicit none
         integer(INTD):: perm_sign                   !out: permutation sign: {-1,+1}
         class(permutation_t), intent(inout):: this  !inout: permutation
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS; perm_sign=0
         if(this%length.gt.0) then
          if(abs(this%prm(0)).eq.1) then
           perm_sign=this%prm(0)
          else
           call this%set_sign(errc)
           if(errc.eq.TEREC_SUCCESS) perm_sign=this%prm(0)
          endif
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function PermutationGetSign
!--------------------------------------------------------
        subroutine PermutationSetSign(this,ierr,new_sign)
!Sets either a user-defined sign or computes one from the permutation.
         implicit none
         class(permutation_t), intent(inout):: this     !inout: permutation
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD), intent(in), optional:: new_sign !in: new sign (by user)
         integer(INTD):: errc,i
         integer(INTD), allocatable:: trn(:)

         errc=TEREC_SUCCESS
         if(this%length.gt.1) then
          allocate(trn(0:this%length),STAT=errc)
          if(errc.eq.0) then
           trn(0:this%length)=(/+1,(i,i=1,this%length)/)
           call merge_sort_key_int(this%length,this%prm(1:this%length),trn(0:this%length))
           if(abs(trn(0)).eq.1) then
            this%prm(0)=trn(0)
           else
            errc=TEREC_UNABLE_COMPLETE
           endif
           deallocate(trn)
          else
           errc=TEREC_MEM_ALLOC_FAILED
          endif
         elseif(this%length.eq.1) then
          this%prm(0)=+1
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine PermutationSetSign
!----------------------------------------------
        subroutine PermutationInvert(this,ierr)
!Inverts the permutation.
         implicit none
         class(permutation_t), intent(inout):: this  !inout: permutation (in:original, out:inverted)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: trn(1:this%length)
         integer(INTD):: errc,i,j

         errc=TEREC_SUCCESS
         if(this%length.gt.0) then
          trn(1:this%length)=this%prm(1:this%length)
          this%prm(1:this%length)=0
          do i=1,this%length
           j=trn(i)
           if(j.gt.0.and.j.le.this%length) then
            if(this%prm(j).eq.0) then
             this%prm(j)=i
            else
             errc=TEREC_OBJ_CORRUPTED; exit
            endif
           else
            errc=TEREC_OBJ_CORRUPTED; exit
           endif
          enddo
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine PermutationInvert
!----------------------------------------
        subroutine permutation_dtor(this)
         implicit none
         type(permutation_t):: this

         if(allocated(this%prm)) deallocate(this%prm)
         this%length=0
         return
        end subroutine permutation_dtor
![contr_ptrn_ext_t]=======================================================
        subroutine ContrPtrnExtSetIndexCorr(this,nd,nl,nr,contr_ptrn,ierr)
!Sets tensor dimension correspondence in a tensor contraction (contraction pattern).
!No strict validity check for <contr_ptrn>!
         implicit none
         class(contr_ptrn_ext_t), intent(inout):: this   !out: extended contraction pattern spec
         integer(INTD), intent(in):: nd                  !in: destination tensor rank
         integer(INTD), intent(in):: nl                  !in: left tensor rank
         integer(INTD), intent(in):: nr                  !in: right tensor rank
         integer(INTD), intent(in):: contr_ptrn(1:nl+nr) !in: valid basic digital tensor contraction pattern (see TAL-SH specs)
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,i,j

         errc=TEREC_SUCCESS
         if(nd.ge.0.and.nd.le.MAX_TENSOR_RANK.and.nl.ge.0.and.nl.le.MAX_TENSOR_RANK.and.nr.ge.0.and.nr.le.MAX_TENSOR_RANK) then
          do i=1,nl
           j=contr_ptrn(i)
           if(j.gt.0) then !uncontracted left index
            this%lind_pos(i)=-j; this%dind_pos(j)=-i
           elseif(j.lt.0) then !contracted index
            this%lind_pos(i)=-j; this%rind_pos(-j)=i
           else
            errc=TEREC_INVALID_ARGS; exit
           endif
          enddo
          if(errc.eq.TEREC_SUCCESS) then
           do i=1,nr
            j=contr_ptrn(nl+i)
            if(j.gt.0) then !uncontracted right index
             this%rind_pos(i)=-j; this%dind_pos(j)=i
            elseif(j.lt.0) then !contracted index
             this%rind_pos(i)=-j; this%lind_pos(-j)=i
            else
             errc=TEREC_INVALID_ARGS; exit
            endif
           enddo
           if(errc.eq.TEREC_SUCCESS) then
            this%ddim=nd; this%ldim=nl; this%rdim=nr
            this%dind_res(1:nd)=0; this%lind_res(1:nl)=0; this%rind_res(1:nr)=0 !clear restrictions
            this%ind_restr_set=.FALSE.
           endif
          endif
         else
          errc=TEREC_INVALID_ARGS
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContrPtrnExtSetIndexCorr
!-------------------------------------------------------------------------
        subroutine ContrPtrnExtSetStoreSymm(this,tens_num,restr_inds,ierr)
!Sets index permutational symmetries due to tensor storage. At least two
!tensor dimensions must be passed here to set up a dimension symmetry dependency.
!Previously set dimension symmetries are kept intact.
         implicit none
         class(contr_ptrn_ext_t), intent(inout), target:: this !inout: extended contraction pattern spec
         integer(INTD), intent(in):: tens_num                  !in: tensor number (0:D,1:L,2:R)
         integer(INTD), intent(in):: restr_inds(1:)            !in: restricted index positions in tensor <tens_num>
         integer(INTD), intent(out), optional:: ierr           !out: error code
         integer(INTD):: errc,n,m,i,j0,j1
         integer(INTD), pointer:: ind_res(:)

         if(this%is_set(errc)) then
          if(tens_num.eq.0) then !destination tensor
           n=this%ddim; ind_res=>this%dind_res
          elseif(tens_num.eq.1) then !left input tensor
           n=this%ldim; ind_res=>this%lind_res
          elseif(tens_num.eq.2) then !right input tensor
           n=this%rdim; ind_res=>this%rind_res
          else
           errc=TEREC_INVALID_ARGS
          endif
          if(errc.eq.TEREC_SUCCESS) then
           m=size(restr_inds)
           if(m.ge.2) then !at least two index positions are expected
            if(m.le.n) then
             do i=2,m
              j0=restr_inds(i-1); j1=restr_inds(i)
              if(j0.gt.0.and.j0.lt.n.and.j1.gt.1.and.j1.le.n.and.j0.lt.j1) then
               ind_res(j1)=j0 !index j1 on the right depends on index j0 on the left
              else
               errc=TEREC_INVALID_ARGS; exit
              endif
             enddo
             ind_res(restr_inds(1))=-1 !first restricted index mark (-1)
             if(errc.eq.TEREC_SUCCESS) this%ind_restr_set=.TRUE.
            else
             errc=TEREC_INVALID_ARGS
            endif
           elseif(m.eq.1) then
            errc=TEREC_INVALID_REQUEST
           endif
          endif
         else
          if(errc.eq.TEREC_SUCCESS) errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContrPtrnExtSetStoreSymm
!--------------------------------------------------------------------------------------
        subroutine ContrPtrnExtSetOperlSymm(this,tens_num,restr_inds,ierr,destroy_symm)
!Sets index permutational symmetries due to tensor operation. At least one tensor
!dimension must be passed here. Provided that <destroy_symm>=TRUE, if a tensor
!dimension is a part of a pre-exsiting symmetry group, the group will be destroyed.
!Otherwise, the symmetry may be imposed only on non-symmetric tensor dimensions.
!To set up a new symmetry group, at least two tensor dimensions must be passed here.
         implicit none
         class(contr_ptrn_ext_t), intent(inout), target:: this !inout: extended contraction pattern spec
         integer(INTD), intent(in):: tens_num                  !in: tensor argument number (0:D,1:L,2:R)
         integer(INTD), intent(in):: restr_inds(1:)            !in: restricted index positions in tensor <tens_num>
         integer(INTD), intent(out), optional:: ierr           !out: error code
         logical, intent(in), optional:: destroy_symm          !in: if TRUE, an attempt to reset up an existing symmetric group will destroy it first (defaults to FALSE)
         integer(INTD):: errc,n,m,i,j0,j1
         integer(INTD), pointer:: ind_res(:)
         logical:: destr

         if(this%is_set(errc)) then
          if(tens_num.eq.0) then !destination tensor
           n=this%ddim; ind_res=>this%dind_res
          elseif(tens_num.eq.1) then !left input tensor
           n=this%ldim; ind_res=>this%lind_res
          elseif(tens_num.eq.2) then !right input tensor
           n=this%rdim; ind_res=>this%rind_res
          else
           errc=TEREC_INVALID_ARGS
          endif
          if(errc.eq.TEREC_SUCCESS) then
           destr=.FALSE.; if(present(destroy_symm)) destr=destroy_symm
           m=size(restr_inds)
           if(m.ge.1) then
            if(m.le.n) then
 !Check:
             do i=2,m
              j0=restr_inds(i-1); j1=restr_inds(i)
              if(.not.(j0.ge.1.and.j0.lt.n.and.j1.ge.2.and.j1.le.n.and.j0.lt.j1)) then
               errc=TEREC_INVALID_ARGS; exit
              endif
             enddo
             if(errc.eq.TEREC_SUCCESS.and.(.not.destr)) then
              if(m.ge.2) then
               do i=1,m
                if(ind_res(restr_inds(i)).ne.0) then; errc=TEREC_INVALID_REQUEST; exit; endif
               enddo
              else
               errc=TEREC_INVALID_REQUEST
              endif
             else
 !Destroy previous symmetry group:
              if(errc.eq.TEREC_SUCCESS) then
               do i=1,m
                j0=restr_inds(i)
                do while(ind_res(j0).gt.0)
                 j1=ind_res(j0); ind_res(j0)=0; j0=j1
                enddo
                ind_res(j0)=0
                j0=restr_inds(i)
                do j1=restr_inds(i)+1,n
                 if(ind_res(j1).eq.j0) then; ind_res(j1)=0; j0=j1; endif
                enddo
               enddo
              endif
             endif
             if(errc.eq.TEREC_SUCCESS) then
 !Construct new symmetry group:
              do i=2,m
               ind_res(restr_inds(i))=restr_inds(i-1)
              enddo
              if(m.gt.1) then
               ind_res(restr_inds(1))=-1 !first restricted index mark (-1)
               this%ind_restr_set=.TRUE.
              endif
             endif
            else
             errc=TEREC_INVALID_ARGS
            endif
           endif
          endif
         else
          if(errc.eq.TEREC_SUCCESS) errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContrPtrnExtSetOperlSymm
!----------------------------------------------------------------------
        subroutine ContrPtrnExtBreakDimSymm(this,tens_num,dim_num,ierr)
!Breaks permutational symmetry on a specific tensor dimension.
         implicit none
         class(contr_ptrn_ext_t), intent(inout), target:: this !inout: extended tensor contraction pattern
         integer(INTD), intent(in):: tens_num                  !in: tensor argument number (0:D,1:L,2:R)
         integer(INTD), intent(in):: dim_num                   !in: tensor dimension number: [1:rank]
         integer(INTD), intent(out), optional:: ierr           !out: error code
         integer(INTD):: errc,n,i
         integer(INTD), pointer:: ind_res(:)
         logical:: last

         if(this%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) then
           if(tens_num.eq.0) then !destination tensor
            n=this%ddim; ind_res=>this%dind_res
           elseif(tens_num.eq.1) then !left input tensor
            n=this%ldim; ind_res=>this%lind_res
           elseif(tens_num.eq.2) then !right input tensor
            n=this%rdim; ind_res=>this%rind_res
           else
            errc=TEREC_INVALID_ARGS
           endif
           if(errc.eq.TEREC_SUCCESS) then
            if(dim_num.gt.0.and.dim_num.le.n) then
             if(this%ind_restr_set) then
              i=ind_res(dim_num)
              if(i.gt.0) then
               if(i.lt.dim_num) then
                if(ind_res(i).lt.0) ind_res(i)=0 !first dimension from a symmetric group
               else
                errc=TEREC_OBJ_CORRUPTED
               endif
              endif
              if(errc.eq.TEREC_SUCCESS) then
               last=.TRUE.
               do i=dim_num+1,n
                if(ind_res(i).eq.dim_num) then; last=.FALSE.; exit; endif
               enddo
               if(last) then; ind_res(dim_num)=0; else; ind_res(dim_num)=-1; endif
              endif
             endif
            else
             errc=TEREC_INVALID_ARGS
            endif
           endif
          endif
         else
          if(errc.eq.TEREC_SUCCESS) errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContrPtrnExtBreakDimSymm
!------------------------------------------------------
        subroutine ContrPtrnExtUnpack(this,packet,ierr)
!Unpacks an object from a packet.
         implicit none
         class(contr_ptrn_ext_t), intent(inout):: this !out: extended tensor contraction pattern
         class(obj_pack_t), intent(inout):: packet     !inout: packet
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: errc

         call unpack_builtin(packet,this%ddim,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%ldim,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%rdim,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%ind_restr_set,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%dind_pos,this%ddim,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%lind_pos,this%ldim,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%rind_pos,this%rdim,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%dind_res,this%ddim,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%lind_res,this%ldim,errc)
         if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%rind_res,this%rdim,errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine ContrPtrnExtUnpack
!--------------------------------------------------------------
        function ContrPtrnExtIsSet(this,ierr,restr) result(ans)
!Returns TRUE if the extended tensor contraction pattern is set.
         implicit none
         logical:: ans                               !out: answer
         class(contr_ptrn_ext_t), intent(in):: this  !in: extended tensor contraction pattern
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(out), optional:: restr      !out: TRUE of index permutational symmetries are present
         integer(INTD):: errc

         errc=TEREC_SUCCESS
         ans=(this%ddim.ge.0.and.this%ldim.ge.0.and.this%rdim.ge.0)
         if(present(restr)) restr=ans.and.this%ind_restr_set
         if(present(ierr)) ierr=errc
         return
        end function ContrPtrnExtIsSet
!----------------------------------------------------
        subroutine ContrPtrnExtPack(this,packet,ierr)
!Packs the object into a packet.
         implicit none
         class(contr_ptrn_ext_t), intent(in):: this  !in: extended tensor contraction pattern
         class(obj_pack_t), intent(inout):: packet   !inout: packet
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         call pack_builtin(packet,this%ddim,errc)
         if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%ldim,errc)
         if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%rdim,errc)
         if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%ind_restr_set,errc)
         if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%dind_pos,this%ddim,errc)
         if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%lind_pos,this%ldim,errc)
         if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%rind_pos,this%rdim,errc)
         if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%dind_res,this%ddim,errc)
         if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%lind_res,this%ldim,errc)
         if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%rind_res,this%rdim,errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine ContrPtrnExtPack
!----------------------------------------------------------------------
        subroutine ContrPtrnExtGetContrPtrn(this,nl,nr,contr_ptrn,ierr)
!Returns the basic digital tensor contraction pattern used by TAL-SH.
         implicit none
         class(contr_ptrn_ext_t), intent(in):: this     !in: extended tensor contraction pattern
         integer(INTD), intent(out):: nl                !out: left tensor rank
         integer(INTD), intent(out):: nr                !out: right tensor rank
         integer(INTD), intent(inout):: contr_ptrn(1:*) !out: basic digital tensor contraction pattern
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD):: errc,i

         if(this%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) then
           nl=this%ldim; nr=this%rdim
           do i=1,nl; contr_ptrn(i)=-this%lind_pos(i); enddo
           do i=1,nr; contr_ptrn(nl+i)=-this%rind_pos(i); enddo
          endif
         else
          if(errc.eq.TEREC_SUCCESS) errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContrPtrnExtGetContrPtrn
!---------------------------------------------------------------------------------
        subroutine ContrPtrnExtGetDimSymmetry(this,tens_num,num_dim,dim_symm,ierr)
!Returns the dimension symmetry restrictions for a specific tensor argument.
         implicit none
         class(contr_ptrn_ext_t), intent(in):: this  !in: extended tensor contraction pattern
         integer(INTD), intent(in):: tens_num        !in: tensor argument number (0:D,1:L,2:R)
         integer(INTD), intent(out):: num_dim        !out: number of tensor dimensions
         integer(INTD), intent(inout):: dim_symm(1:) !out: tensor dimension symmetry restrictions
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         num_dim=-1
         if(this%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) then
           select case(tens_num)
           case(0)
            num_dim=this%ddim
            dim_symm(1:num_dim)=this%dind_res(1:num_dim)
           case(1)
            num_dim=this%ldim
            dim_symm(1:num_dim)=this%lind_res(1:num_dim)
           case(2)
            num_dim=this%rdim
            dim_symm(1:num_dim)=this%rind_res(1:num_dim)
           case default
            errc=TEREC_INVALID_ARGS
           end select
          endif
         else
          if(errc.eq.TEREC_SUCCESS) errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContrPtrnExtGetDimSymmetry
!---------------------------------------------------------------
        subroutine ContrPtrnExtPrintIt(this,ierr,dev_id,nspaces)
!Prints the extended tensor contraction.
         implicit none
         class(contr_ptrn_ext_t), intent(in):: this    !in: extended tensor contraction pattern
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD), intent(in), optional:: dev_id  !in: output device id (defaults to screeen)
         integer(INTD), intent(in), optional:: nspaces !in: number of leading spaces
         integer(INTD):: errc,devo,nsp,i,l
         character(4):: num_pos

         devo=6; if(present(dev_id)) devo=dev_id
         nsp=0; if(present(nspaces)) nsp=nspaces
         if(this%is_set(errc)) then
          do i=1,nsp; write(devo,'(" ")',ADVANCE='NO'); enddo
          if(this%ddim.gt.0) then
           call numchar(this%ddim,l,num_pos)
           write(devo,'("{",'//num_pos(1:l)//'(1x,i3),"}=")',ADVANCE='NO') this%dind_pos(1:this%ddim)
          else
           write(devo,'("{}=")',ADVANCE='NO')
          endif
          if(this%ldim.gt.0) then
           call numchar(this%ldim,l,num_pos)
           write(devo,'("{",'//num_pos(1:l)//'(1x,i3),"}*")',ADVANCE='NO') this%lind_pos(1:this%ldim)
          else
           write(devo,'("{}*")',ADVANCE='NO')
          endif
          if(this%rdim.gt.0) then
           call numchar(this%rdim,l,num_pos)
           write(devo,'("{",'//num_pos(1:l)//'(1x,i3),"}")',ADVANCE='NO') this%rind_pos(1:this%rdim)
          else
           write(devo,'("{}")',ADVANCE='NO')
          endif
          if(this%ind_restr_set) then
           write(devo,'(":     ")',ADVANCE='NO')
           if(this%ddim.gt.0) then
            call numchar(this%ddim,l,num_pos)
            write(devo,'("{",'//num_pos(1:l)//'(1x,i2),"}=")',ADVANCE='NO') this%dind_res(1:this%ddim)
           else
            write(devo,'("{}=")',ADVANCE='NO')
           endif
           if(this%ldim.gt.0) then
            call numchar(this%ldim,l,num_pos)
            write(devo,'("{",'//num_pos(1:l)//'(1x,i2),"}*")',ADVANCE='NO') this%lind_res(1:this%ldim)
           else
            write(devo,'("{}*")',ADVANCE='NO')
           endif
           if(this%rdim.gt.0) then
            call numchar(this%rdim,l,num_pos)
            write(devo,'("{",'//num_pos(1:l)//'(1x,i2),"}")',ADVANCE='NO') this%rind_res(1:this%rdim)
           else
            write(devo,'("{}")',ADVANCE='NO')
           endif
          endif
          write(devo,'()')
         else
          do i=1,nsp; write(devo,'(" ")',ADVANCE='NO'); enddo
          write(devo,'("{Empty/invalid extended tensor contraction}")')
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine ContrPtrnExtPrintIt
![tens_contraction_t]=============================
        subroutine TensContractionAssign(this,src)
!Copy assignment.
         implicit none
         class(tens_contraction_t), intent(out):: this !out: copy
         class(tens_contraction_t), intent(in):: src   !in: source
         integer(INTD):: i

         this%num_args=src%num_args
         do i=0,src%num_args-1
          this%tens_arg(i)=src%tens_arg(i)
         enddo
         this%contr_ptrn=src%contr_ptrn
         this%alpha=src%alpha
         return
        end subroutine TensContractionAssign
!-----------------------------------------------------------
        function TensContractionIsSet(this,ierr) result(ans)
!Returns TRUE if the tensor contraction is fully set.
         implicit none
         logical:: ans                                !out: answer
         class(tens_contraction_t), intent(in):: this !in: tensor contraction
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc,n
         class(tens_rcrsv_t), pointer:: trp

         ans=.FALSE.
         if(this%get_num_args(errc).eq.3) then !binary tensor contraction has 3 tensor arguments: D=L*R
          if(errc.eq.TEREC_SUCCESS) then
           trp=>this%get_argument(0,errc) !destination tensor
           if(errc.eq.TEREC_SUCCESS) then
            if(trp%is_set(errc,num_dims=n)) then
             if(errc.eq.TEREC_SUCCESS.and.n.eq.this%contr_ptrn%ddim) then
              trp=>this%get_argument(1,errc) !left tensor
              if(errc.eq.TEREC_SUCCESS) then
               if(trp%is_set(errc,num_dims=n)) then
                if(errc.eq.TEREC_SUCCESS.and.n.eq.this%contr_ptrn%ldim) then
                 trp=>this%get_argument(2,errc) !right tensor
                 if(errc.eq.TEREC_SUCCESS) then
                  if(trp%is_set(errc,num_dims=n)) then
                   if(errc.eq.TEREC_SUCCESS.and.n.eq.this%contr_ptrn%rdim) then
                    ans=.TRUE.
                   endif
                  endif
                 endif
                endif
               endif
              endif
             endif
            endif
           endif
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensContractionIsSet
!--------------------------------------------------------------
        function TensContractionArgsFull(this,ierr) result(ans)
!Returns TRUE if all tensor contraction arguments have been set.
         implicit none
         logical:: ans                                !out: answer
         class(tens_contraction_t), intent(in):: this !in: tensor contraction
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         ans=(this%get_num_args(errc).eq.3) !binary tensor contraction has 3 arguments: D=L*R
         ans=ans.and.(errc.eq.TEREC_SUCCESS)
         if(present(ierr)) ierr=errc
         return
        end function TensContractionArgsFull
!-------------------------------------------------------------------------
        subroutine TensContractionSetContrPtrn(this,contr_ptrn,ierr,alpha)
!Sets the extended tensor contraction pattern (all tensor arguments must have been set already).
!Additionally, if tensor arguments have permutational symmetries, they will be incorporated into
!the extended tensor contraction pattern. Later on, additional permutational symmetries can be
!imposed onto non-symmetric indices. An optional tensor contraction prefactor <alpha> can be supplied.
         implicit none
         class(tens_contraction_t), intent(inout):: this !inout: tensor contraction
         integer(INTD), intent(in):: contr_ptrn(1:*)     !in: basic digital tensor contraction pattern (see TAL-SH specs)
         integer(INTD), intent(out), optional:: ierr     !out: error code
         complex(8), intent(in), optional:: alpha        !in: complex tensor contraction perfactor (defaults to 1.0)
         integer(INTD):: errc,nd,nl,nr,m,i,l,grs,restr_dims(1:MAX_TENSOR_RANK)
         class(tens_rcrsv_t), pointer:: dtrp,ltrp,rtrp
         class(tens_header_t), pointer:: thp

         if(this%args_full(errc)) then !all tensor arguments must have been set already
          if(errc.eq.TEREC_SUCCESS) then
           dtrp=>this%get_argument(0,errc)
           if(errc.eq.TEREC_SUCCESS) then
            if(dtrp%is_set(errc,num_dims=nd)) then
             ltrp=>this%get_argument(1,errc)
             if(errc.eq.TEREC_SUCCESS) then
              if(ltrp%is_set(errc,num_dims=nl)) then
               rtrp=>this%get_argument(2,errc)
               if(errc.eq.TEREC_SUCCESS) then
                if(rtrp%is_set(errc,num_dims=nr)) then
 !Basic contraction pattern:
                 call this%contr_ptrn%set_index_corr(nd,nl,nr,contr_ptrn,errc)
 !Destination tensor storage symmetries:
                 if(errc.eq.TEREC_SUCCESS) then
                  thp=>dtrp%get_header(errc)
                  if(errc.eq.TEREC_SUCCESS) then
                   m=thp%num_groups(errc)
                   if(errc.eq.TEREC_SUCCESS.and.m.gt.0) then
                    do i=1,m
                     call thp%get_group(i,restr_dims,l,errc,grs); if(errc.ne.TEREC_SUCCESS) exit
                     if(l.ge.2) then
                      call this%contr_ptrn%set_store_symm(0,restr_dims(1:l),errc); if(errc.ne.TEREC_SUCCESS) exit
                     endif
                    enddo
                   endif
                  endif
                 endif
 !Left tensor storage symmetries:
                 if(errc.eq.TEREC_SUCCESS) then
                  thp=>ltrp%get_header(errc)
                  if(errc.eq.TEREC_SUCCESS) then
                   m=thp%num_groups(errc)
                   if(errc.eq.TEREC_SUCCESS.and.m.gt.0) then
                    do i=1,m
                     call thp%get_group(i,restr_dims,l,errc,grs); if(errc.ne.TEREC_SUCCESS) exit
                     if(l.ge.2) then
                      call this%contr_ptrn%set_store_symm(1,restr_dims(1:l),errc); if(errc.ne.TEREC_SUCCESS) exit
                     endif
                    enddo
                   endif
                  endif
                 endif
 !Right tensor storage symmetries:
                 if(errc.eq.TEREC_SUCCESS) then
                  thp=>rtrp%get_header(errc)
                  if(errc.eq.TEREC_SUCCESS) then
                   m=thp%num_groups(errc)
                   if(errc.eq.TEREC_SUCCESS.and.m.gt.0) then
                    do i=1,m
                     call thp%get_group(i,restr_dims,l,errc,grs); if(errc.ne.TEREC_SUCCESS) exit
                     if(l.ge.2) then
                      call this%contr_ptrn%set_store_symm(2,restr_dims(1:l),errc); if(errc.ne.TEREC_SUCCESS) exit
                     endif
                    enddo
                   endif
                  endif
                 endif
 !Tensor contraction prefactor:
                 if(errc.eq.TEREC_SUCCESS) then
                  if(present(alpha)) then; this%alpha=alpha; else; this%alpha=(1d0,0d0); endif
                 endif
                else
                 if(errc.eq.TEREC_SUCCESS) errc=TEREC_INVALID_REQUEST
                endif
               endif
              else
               if(errc.eq.TEREC_SUCCESS) errc=TEREC_INVALID_REQUEST
              endif
             endif
            else
             if(errc.eq.TEREC_SUCCESS) errc=TEREC_INVALID_REQUEST
            endif
           endif
          endif
         else
          if(errc.eq.TEREC_SUCCESS) errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensContractionSetContrPtrn
!----------------------------------------------------------------------------
        subroutine TensContractionSetOperlSymm(this,tens_num,restr_inds,ierr)
!Sets index permutational symmetry restrictions due to tensor operation
!(both contraction pattern and arguments must have been set already).
!If any of the given tensor dimensions already belongs to a symmetric group,
!an error will be returned.
         implicit none
         class(tens_contraction_t), intent(inout):: this !inout: tensor contraction
         integer(INTD), intent(in):: tens_num            !in: specific tensor argument (0:D,1:L,2:R)
         integer(INTD), intent(in):: restr_inds(1:)      !in: restricted tensor dimensions
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc

         if(this%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) call this%contr_ptrn%set_operl_symm(tens_num,restr_inds,errc)
         else
          if(errc.eq.TEREC_SUCCESS) errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensContractionSetOperlSymm
!---------------------------------------------------------
        subroutine TensContractionUnpack(this,packet,ierr)
!Unpacks the object from a packet.
         implicit none
         class(tens_contraction_t), intent(inout):: this !in: tensor contraction
         class(obj_pack_t), intent(inout):: packet       !inout: packet
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,i
         class(tens_rcrsv_t), pointer:: tens_p

         call unpack_builtin(packet,this%num_args,errc)
         if(errc.eq.PACK_SUCCESS) then
          do i=0,this%num_args-1
           call this%allocate_argument(errc); if(errc.ne.TEREC_SUCCESS) exit
           tens_p=>this%get_argument(i,errc); if(errc.ne.TEREC_SUCCESS) exit
           call tens_p%tens_rcrsv_ctor(packet,errc); if(errc.ne.TEREC_SUCCESS) exit
          enddo
          if(errc.eq.PACK_SUCCESS) call this%contr_ptrn%unpack(packet,errc)
          if(errc.eq.PACK_SUCCESS) call unpack_builtin(packet,this%alpha,errc)
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensContractionUnpack
!-------------------------------------------------------
        subroutine TensContractionPack(this,packet,ierr)
!Packs the object into a packet.
         implicit none
         class(tens_contraction_t), intent(in):: this !in: tensor contraction
         class(obj_pack_t), intent(inout):: packet    !inout: packet
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc,i
         class(tens_rcrsv_t), pointer:: tens_p

         call pack_builtin(packet,this%num_args,errc)
         if(errc.eq.PACK_SUCCESS) then
          do i=0,this%num_args-1
           tens_p=>this%get_argument(i,errc); if(errc.ne.TEREC_SUCCESS) exit
           call tens_p%pack(packet,errc); if(errc.ne.TEREC_SUCCESS) exit
          enddo
          if(errc.eq.PACK_SUCCESS) call this%contr_ptrn%pack(packet,errc)
          if(errc.eq.PACK_SUCCESS) call pack_builtin(packet,this%alpha,errc)
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensContractionPack
!------------------------------------------------------------------------
        function TensContractionGetPrefactor(this,ierr) result(prefactor)
!Returns the scalar prefactor.
         implicit none
         complex(8):: prefactor                       !out: tensor contraction prefactor
         class(tens_contraction_t), intent(in):: this !in: tensor contraction
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         prefactor=this%alpha
         if(.not.this%is_set(errc)) errc=TEREC_INVALID_REQUEST
         if(present(ierr)) ierr=errc
         return
        end function TensContractionGetPrefactor
!------------------------------------------------------------------------------
        function TensContractionGetExtContrPtrn(this,ierr) result(contr_ptrn_p)
!Returns a pointer to the extended tensor contraction pattern.
         implicit none
         class(contr_ptrn_ext_t), pointer:: contr_ptrn_p      !out: pointer to the extended tensor contraction pattern
         class(tens_contraction_t), target, intent(in):: this !in: tensor contraction
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTD):: errc

         errc=TEREC_SUCCESS; contr_ptrn_p=>NULL()
         if(this%is_set()) then
          contr_ptrn_p=>this%contr_ptrn
         else
          errc=TEREC_INVALID_REQUEST
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensContractionGetExtContrPtrn
!-------------------------------------------------------------------------
        subroutine TensContractionGetContrPtrn(this,nl,nr,contr_ptrn,ierr)
!Returns the basic digital tensor contraction pattern used by TAL-SH.
         implicit none
         class(tens_contraction_t), intent(in):: this   !in: tensor contraction
         integer(INTD), intent(out):: nl                !out: left tensor rank
         integer(INTD), intent(out):: nr                !out: right tensor rank
         integer(INTD), intent(inout):: contr_ptrn(1:*) !out: basic digital tensor contraction pattern
         integer(INTD), intent(out), optional:: ierr    !out: error code
         integer(INTD):: errc

         call this%contr_ptrn%get_contr_ptrn(nl,nr,contr_ptrn,errc)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensContractionGetContrPtrn
!--------------------------------------------------------------------------------------
        subroutine TensContractionImportReplace(this,tens_contr,ierr,dtens,ltens,rtens)
!Creates a new tensor contraction by replacing tensor arguments in
!an existing tensor contraction (plus symmetry adjustment).
         implicit none
         class(tens_contraction_t), intent(out):: this            !out: derived tensor contraction
         class(tens_contraction_t), intent(in):: tens_contr       !in: parental tensor contraction
         integer(INTD), intent(out), optional:: ierr              !out: error code
         type(tens_rcrsv_t), intent(in), target, optional:: dtens !in: new destination tensor argument
         type(tens_rcrsv_t), intent(in), target, optional:: ltens !in: new left tensor argument
         type(tens_rcrsv_t), intent(in), target, optional:: rtens !in: new right tensor argument
         integer(INTD):: i,errc,nd,cmp
         integer(INTL):: sidx(1:MAX_TENSOR_RANK)
         type(tens_header_t), pointer:: thp
         type(hspace_reg_t):: ths(1:MAX_TENSOR_RANK)
         type(h_space_t), pointer:: hsp

         if(tens_contr%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) then
           this=tens_contr !clone the parental tensor contraction
 !Destination tensor replacement:
           if(present(dtens)) then
            call this%reset_argument(dtens,0,errc)
            if(errc.eq.TEREC_SUCCESS) then
             thp=>dtens%get_header(errc)
             if(errc.eq.TEREC_SUCCESS) then
              call thp%get_spec(sidx,nd,errc,ths)
              if(errc.eq.TEREC_SUCCESS) then
               do i=2,nd
                if(ths(i-1)%get_space_id().eq.ths(i)%get_space_id()) then
                 hsp=>ths(i)%get_space(errc); if(errc.ne.TEREC_SUCCESS) exit
                 cmp=hsp%compare_subranges(sidx(i-1),sidx(i))
                 if(cmp.eq.CMP_ER) then; errc=TEREC_ERROR; exit; endif
                 if(cmp.eq.CMP_LT) then !`I assume LT index ordering in a tensor contraction
                  call this%contr_ptrn%break_dim_symm(0,i,errc); if(errc.ne.TEREC_SUCCESS) exit
                 endif
                endif
               enddo
              endif
             endif
            endif
           endif
 !Left tensor replacement:
           if(present(ltens).and.errc.eq.TEREC_SUCCESS) then
            call this%reset_argument(ltens,1,errc)
            if(errc.eq.TEREC_SUCCESS) then
             thp=>ltens%get_header(errc)
             if(errc.eq.TEREC_SUCCESS) then
              call thp%get_spec(sidx,nd,errc,ths)
              if(errc.eq.TEREC_SUCCESS) then
               do i=2,nd
                if(ths(i-1)%get_space_id().eq.ths(i)%get_space_id()) then
                 hsp=>ths(i)%get_space(errc); if(errc.ne.TEREC_SUCCESS) exit
                 cmp=hsp%compare_subranges(sidx(i-1),sidx(i))
                 if(cmp.eq.CMP_ER) then; errc=TEREC_ERROR; exit; endif
                 if(cmp.eq.CMP_LT) then !`I assume LT index ordering in a tensor contraction
                  call this%contr_ptrn%break_dim_symm(1,i,errc); if(errc.ne.TEREC_SUCCESS) exit
                 endif
                endif
               enddo
              endif
             endif
            endif
           endif
 !Right tensor replacement:
           if(present(rtens).and.errc.eq.TEREC_SUCCESS) then
            call this%reset_argument(rtens,2,errc)
            if(errc.eq.TEREC_SUCCESS) then
             thp=>rtens%get_header(errc)
             if(errc.eq.TEREC_SUCCESS) then
              call thp%get_spec(sidx,nd,errc,ths)
              if(errc.eq.TEREC_SUCCESS) then
               do i=2,nd
                if(ths(i-1)%get_space_id().eq.ths(i)%get_space_id()) then
                 hsp=>ths(i)%get_space(errc); if(errc.ne.TEREC_SUCCESS) exit
                 cmp=hsp%compare_subranges(sidx(i-1),sidx(i))
                 if(cmp.eq.CMP_ER) then; errc=TEREC_ERROR; exit; endif
                 if(cmp.eq.CMP_LT) then !`I assume LT index ordering in a tensor contraction
                  call this%contr_ptrn%break_dim_symm(2,i,errc); if(errc.ne.TEREC_SUCCESS) exit
                 endif
                endif
               enddo
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
        end subroutine TensContractionImportReplace
!--------------------------------------------------------------------------------
        subroutine TensContractionSplit(this,tens_split_f,subops,ierr,num_subops)
!Splits a defined tensor contraction into a list of subcontractions. The splitting
!is driven by the tensor argument decomposition into unique subtensors. More precisely,
!each tensor argument is represented as a direct sum of its unique children subtensors.
!Then all possible non-zero combinations of those subtensors that match the tensor
!contraction pattern will form the list of subcontractions.
         implicit none
         class(tens_contraction_t), intent(in):: this      !in: parental tensor contraction
         procedure(tens_rcrsv_split_i):: tens_split_f      !in: tensor splitting function (splits a tensor into a vector of unique subtensors)
         type(list_bi_t), intent(inout):: subops           !inout: list of subtensor contractions
         integer(INTD), intent(out), optional:: ierr       !out: error code
         integer(INTD), intent(out), optional:: num_subops !out: number of subcontractions generated from the parental tensor contraction here
         integer(INTD):: i,errc,nsub,tcgl
         integer(INTD):: drank,lrank,rrank,dsl,lsl,rsl,dstart,dfinish,lstart,lfinish,rstart,rfinish,pstart,pfinish
         integer(INTD):: nci,cptrn(1:MAX_TENSOR_RANK*2)  !basic tensor contraction pattern (as specified by TAL-SH)
         integer(INTD):: ord(1:MAX_TENSOR_RANK,0:2)      !N2O dimension order (0:destination,1:left,2:right)
         integer(INTD):: adj(1:MAX_TENSOR_RANK,0:2)      !tree level adjustment for dimensions of tensor arguments (0:destination,1:left,2:right)
         integer(INTL):: dmsi,lmsi,rmsi                  !max subspace id in the sorting lists for D,L,R
         type(hspace_reg_t):: ths(1:MAX_TENSOR_RANK,0:2) !h_space of each dimension of each tensor argument (0:destination,1:left,2:right)
         type(vector_t):: dsubs,lsubs,rsubs              !vector of subtensors for each tensor argument
         type(vector_iter_t):: dvit,lvit,rvit            !vector iterator for each tensor argument
         type(list_iter_t):: slit                        !list iterator for the list of subcontractions
         real(8):: tm(0:6),tmf

         tm(0)=thread_wtime(); tm(1:)=tm(0)
         nsub=0 !number of generated subcontractions
         if(this%is_set(errc)) then
          if(errc.eq.TEREC_SUCCESS) then
           if(check_allocate_buffers()) then !checks/allocates sorting buffers tcg_ind_buf/tcg_num_buf
            tm(1)=thread_wtime(tm(0))
            call generate_subtensors(errc) !generates dsubs, lsubs, and rsubs, and associates them with iterators (dvit,lvit,rvit)
            if(errc.eq.TEREC_SUCCESS) then
             tm(2)=thread_wtime(tm(0))
             call align_levels(errc) !sets cptrn(:), ths(:,:), and adj(:,:) for all tensor arguments (SAT level adjustment for tensor dimensions)
             if(errc.eq.TEREC_SUCCESS) then
              tm(3)=thread_wtime(tm(0))
              call build_descriptors(errc) !generates lists of subtensor descriptors in tcg_ind_buf/tcg_num_buf for each tensor argument
              if(errc.eq.TEREC_SUCCESS) then
               tm(4)=thread_wtime(tm(0))
               errc=slit%init(subops) !iterator for the list of subcontractions
               if(errc.eq.GFC_SUCCESS) then
                call generate_subcontractions(errc) !generates subtensor contractions by matching descriptor lists
                i=slit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.TEREC_SUCCESS) errc=TEREC_ERROR
                tm(5)=thread_wtime(tm(0))
               endif
               i=rvit%delete_all(); if(i.ne.GFC_SUCCESS.and.errc.eq.TEREC_SUCCESS) errc=TEREC_ERROR
               i=rvit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.TEREC_SUCCESS) errc=TEREC_ERROR
               i=lvit%delete_all(); if(i.ne.GFC_SUCCESS.and.errc.eq.TEREC_SUCCESS) errc=TEREC_ERROR
               i=lvit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.TEREC_SUCCESS) errc=TEREC_ERROR
               i=dvit%delete_all(); if(i.ne.GFC_SUCCESS.and.errc.eq.TEREC_SUCCESS) errc=TEREC_ERROR
               i=dvit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.TEREC_SUCCESS) errc=TEREC_ERROR
               tm(6)=thread_wtime(tm(0))
              endif
             endif
            endif
           else
            errc=TEREC_MEM_ALLOC_FAILED
           endif
          endif
         else
          if(errc.eq.TEREC_SUCCESS) errc=TEREC_INVALID_ARGS
         endif
         if(present(num_subops)) then
          if(errc.eq.TEREC_SUCCESS) then; num_subops=nsub; else; num_subops=-1; endif
         endif
         if(present(ierr)) ierr=errc
         tmf=thread_wtime(tm(0))
         !write(CONS_OUT,'("#DEBUG(TensContractionSplit): Timings (msec): ",F8.3,":",6(1x,F8.3))')&
              !&tmf*1d3,tm(1)*1d3,(tm(2)-tm(1))*1d3,(tm(3)-tm(2))*1d3,(tm(4)-tm(3))*1d3,(tm(5)-tm(4))*1d3,(tm(6)-tm(5))*1d3 !debug
         return

         contains

          function check_allocate_buffers() result(jsts) !allocates thread private work buffers
           implicit none
           logical:: jsts
           logical:: jfi,jfn
           integer(INTD):: jerr

           jfi=allocated(tcg_ind_buf); jfn=allocated(tcg_num_buf)
 !Allocate tcg_ind_buf(:,:), if not allocated:
           if(.not.jfi) then
            allocate(tcg_ind_buf(1:MAX_TENSOR_RANK,1:TEREC_TCG_BUF_SIZE*2),STAT=jerr) !twice memory for sorting
            if(jerr.eq.0) then; jfi=.TRUE.; tcg_ind_buf(:,:)=0; endif
           endif
 !Allocate tcg_num_buf(:), if not allocated:
           if(.not.jfn) then
            if(jfi) then
             allocate(tcg_num_buf(1:TEREC_TCG_BUF_SIZE*2),STAT=jerr) !twice memory for sorting
             if(jerr.eq.0) then; jfn=.TRUE.; tcg_num_buf(:)=0; endif
             if(.not.jfn) then; deallocate(tcg_ind_buf); jfi=.FALSE.; endif
            endif
           else
            if(.not.jfi) then; deallocate(tcg_num_buf); jfn=.FALSE.; endif
           endif
           jsts=jfi.and.jfn
           return
          end function check_allocate_buffers

          subroutine generate_subtensors(jerr)
           implicit none
           integer(INTD), intent(out):: jerr
           class(tens_rcrsv_t), pointer:: jtrp

           jerr=TEREC_SUCCESS
 !Destination tensor argument:
           jtrp=>this%get_argument(0,jerr)
           if(jerr.eq.TEREC_SUCCESS) then
            jerr=tens_split_f(jtrp,dsubs,dsl); if(dsl.le.0.and.jerr.eq.TEREC_SUCCESS) jerr=TEREC_ERROR
            !write(CONS_OUT,'("#DEBUG(TensContractionSplit:generate_subtensors): Status ",i10,": D length ",i10)') jerr,dsl !debug
           endif
 !Left tensor argument:
           if(jerr.eq.TEREC_SUCCESS) then
            jtrp=>this%get_argument(1,jerr)
            if(jerr.eq.TEREC_SUCCESS) then
             jerr=tens_split_f(jtrp,lsubs,lsl); if(lsl.le.0.and.jerr.eq.TEREC_SUCCESS) jerr=TEREC_ERROR
             !write(CONS_OUT,'("#DEBUG(TensContractionSplit:generate_subtensors): Status ",i10,": L length ",i10)') jerr,lsl !debug
            endif
           endif
 !Right tensor argument:
           if(jerr.eq.TEREC_SUCCESS) then
            jtrp=>this%get_argument(2,jerr)
            if(jerr.eq.TEREC_SUCCESS) then
             jerr=tens_split_f(jtrp,rsubs,rsl); if(rsl.le.0.and.jerr.eq.TEREC_SUCCESS) jerr=TEREC_ERROR
             !write(CONS_OUT,'("#DEBUG(TensContractionSplit:generate_subtensors): Status ",i10,": R length ",i10)') jerr,rsl !debug
            endif
           endif
 !Init list iterators:
           if(jerr.eq.TEREC_SUCCESS) jerr=dvit%init(dsubs)
           if(jerr.eq.TEREC_SUCCESS) jerr=lvit%init(lsubs)
           if(jerr.eq.TEREC_SUCCESS) jerr=rvit%init(rsubs)
           return
          end subroutine generate_subtensors

          subroutine align_levels(jerr)
           implicit none
           integer(INTD), intent(out):: jerr
           class(*), pointer:: jup
           class(tens_rcrsv_t), pointer:: ltrp,rtrp
           class(tens_header_t), pointer:: lthp,rthp
           class(h_space_t), pointer:: lhsp,rhsp
           integer(INTD):: j1,j2,jl1,jl2
           integer(INTL):: jts(1:MAX_TENSOR_RANK,0:2)

           drank=0; lrank=0; rrank=0
 !Destination tensor argument (just get h_spaces):
           jerr=dvit%reset()
           if(jerr.eq.GFC_SUCCESS) then
            jup=>dvit%get_value(jerr)
            if(jerr.eq.GFC_SUCCESS) then
             ltrp=>NULL(); select type(jup); class is(tens_rcrsv_t); ltrp=>jup; end select
             if(associated(ltrp)) then
              lthp=>ltrp%get_header()
              call lthp%get_spec(jts(:,0),drank,jerr,ths(:,0))
             else
              jerr=TEREC_OBJ_CORRUPTED
             endif
            endif
           endif
 !Left and right tensor arguments (get h_spaces and set SAT level adjustment for all):
           if(jerr.eq.TEREC_SUCCESS) then
            jerr=lvit%reset()
            if(jerr.eq.GFC_SUCCESS) then
             jup=>lvit%get_value(jerr)
             if(jerr.eq.GFC_SUCCESS) then
              ltrp=>NULL(); select type(jup); class is(tens_rcrsv_t); ltrp=>jup; end select
              if(associated(ltrp)) then
               lthp=>ltrp%get_header()
               call lthp%get_spec(jts(:,1),lrank,jerr,ths(:,1))
               if(jerr.eq.TEREC_SUCCESS) then
                jerr=rvit%reset()
                if(jerr.eq.GFC_SUCCESS) then
                 jup=>rvit%get_value(jerr)
                 if(jerr.eq.GFC_SUCCESS) then
                  rtrp=>NULL(); select type(jup); class is(tens_rcrsv_t); rtrp=>jup; end select
                  if(associated(rtrp)) then
                   rthp=>rtrp%get_header()
                   call rthp%get_spec(jts(:,2),rrank,jerr,ths(:,2))
                   if(jerr.eq.TEREC_SUCCESS) then
                    adj(1:drank,0)=0; adj(1:lrank,1)=0; adj(1:rrank,2)=0
                    call this%get_contr_ptrn(j1,j2,cptrn,jerr) !get basic contraction pattern
                    if(jerr.eq.TEREC_SUCCESS) then
                     if(lrank+rrank.gt.0) then
                      do j1=1,lrank !dimensions of the left tensor
                       j2=cptrn(j1)
                       if(j2.lt.0) then !contracted dimension: abs(j2) = position in the right tensor
                        j2=-j2 !corresponding dimension in the right tensor
                        lhsp=>ths(j1,1)%get_space(jerr); if(jerr.ne.TEREC_SUCCESS) exit
                        jl1=lhsp%get_subspace_level(jts(j1,1),jerr); if(jerr.ne.TEREC_SUCCESS) exit
                        rhsp=>ths(j2,2)%get_space(jerr); if(jerr.ne.TEREC_SUCCESS) exit
                        jl2=rhsp%get_subspace_level(jts(j2,2),jerr); if(jerr.ne.TEREC_SUCCESS) exit
                        if(jl1.lt.jl2) then; adj(j2,2)=jl2-jl1; elseif(jl1.gt.jl2) then; adj(j1,1)=jl1-jl2; endif
                       else !uncontracted dimension: j2 = position in the destination tensor
                        lhsp=>ths(j1,1)%get_space(jerr); if(jerr.ne.TEREC_SUCCESS) exit
                        jl1=lhsp%get_subspace_level(jts(j1,1),jerr); if(jerr.ne.TEREC_SUCCESS) exit
                        rhsp=>ths(j2,0)%get_space(jerr); if(jerr.ne.TEREC_SUCCESS) exit
                        jl2=rhsp%get_subspace_level(jts(j2,0),jerr); if(jerr.ne.TEREC_SUCCESS) exit
                        if(jl1.lt.jl2) then; adj(j2,0)=jl2-jl1; elseif(jl1.gt.jl2) then; adj(j1,1)=jl1-jl2; endif
                       endif
                      enddo
                      do j1=1,rrank !dimensions of the right tensor
                       j2=cptrn(lrank+j1)
                       if(j2.gt.0) then !uncontracted dimension: j2 = position in the destination tensor
                        lhsp=>ths(j1,2)%get_space(jerr); if(jerr.ne.TEREC_SUCCESS) exit
                        jl1=lhsp%get_subspace_level(jts(j1,2),jerr); if(jerr.ne.TEREC_SUCCESS) exit
                        rhsp=>ths(j2,0)%get_space(jerr); if(jerr.ne.TEREC_SUCCESS) exit
                        jl2=rhsp%get_subspace_level(jts(j2,0),jerr); if(jerr.ne.TEREC_SUCCESS) exit
                        if(jl1.lt.jl2) then; adj(j2,0)=jl2-jl1; elseif(jl1.gt.jl2) then; adj(j1,2)=jl1-jl2; endif
                       endif
                      enddo
                     endif
                    endif
                   endif
                  else
                   jerr=TEREC_OBJ_CORRUPTED
                  endif
                 endif
                endif
               endif
              else
               jerr=TEREC_OBJ_CORRUPTED
              endif
             endif
            endif
           endif
           nci=(lrank+rrank-drank)/2 !number of contracted indices
           return
          end subroutine align_levels

          subroutine build_descriptors(jerr)
           implicit none
           integer(INTD), intent(out):: jerr
           integer(INTL):: sidx(1:MAX_TENSOR_RANK)
           class(tens_rcrsv_t), pointer:: jtrp
           class(tens_header_t), pointer:: jthp
           class(*), pointer:: jup
           integer(INTD):: ji,ja,jl,jc,jnd,dim_restr(1:MAX_TENSOR_RANK)
           class(h_space_t), pointer:: jhsp
           logical:: approved

           jerr=TEREC_SUCCESS; tcgl=0 !tcgl: current length of the tcg_ind_buf(:)/tcg_num_buf(:)
 !Left subtensors:
           lstart=tcgl+1; lmsi=0_INTL !start offset of the tensor descriptors and maximum subspace id
           call this%contr_ptrn%get_dim_symmetry(1,jnd,dim_restr,jerr)
           if(jerr.eq.TEREC_SUCCESS) then
            jl=0; jerr=lvit%reset()
  !Iterate over subtensors:
            lloop: do while(jerr.eq.GFC_SUCCESS)
   !Get subtensor header:
             jup=>lvit%get_value(jerr); if(jerr.ne.GFC_SUCCESS) exit lloop
             jtrp=>NULL(); select type(jup); class is(tens_rcrsv_t); jtrp=>jup; end select
             if(.not.associated(jtrp)) then; jerr=TEREC_OBJ_CORRUPTED; exit lloop; endif !trap
             jthp=>jtrp%get_header()
             if(.not.associated(jthp)) then; jerr=TEREC_OBJ_CORRUPTED; exit lloop; endif !trap
             call jthp%get_spec(sidx,lrank,jerr); if(jerr.ne.TEREC_SUCCESS) exit lloop
   !Append the subtensor multi-index (descriptor) into the sorting list (adjust dimension SAT level, if needed):
             approved=.TRUE.
             tcgl=tcgl+1; tcg_num_buf(tcgl)=int(jl,INTL) !subtensor number: [0..max]
             do ji=1,lrank
              ja=dim_restr(ji)
              if(ja.gt.0.and.ja.le.lrank) then
               jhsp=>ths(ji,1)%get_space(jerr); if(jerr.ne.TEREC_SUCCESS) exit lloop
               jc=jhsp%compare_subranges(sidx(ja),sidx(ji))
               if(jc.eq.CMP_ER) then; jerr=TEREC_ERROR; exit lloop; endif
               if(jc.eq.CMP_GT) then; approved=.FALSE.; exit; endif !`I assume LT index ordering for symmetric tensor dimensions in a tensor contraction
              endif
              if(adj(ji,1).gt.0) then !promotion to an ancestor SAT level is needed
               jhsp=>ths(ji,1)%get_space(jerr); if(jerr.ne.TEREC_SUCCESS) exit lloop
               tcg_ind_buf(ji,tcgl)=jhsp%get_ancestor_id(sidx(ji),adj(ji,1),jerr); if(jerr.ne.TEREC_SUCCESS) exit lloop
              else
               tcg_ind_buf(ji,tcgl)=sidx(ji)
              endif
             enddo
             if(approved) then
              do ji=1,lrank; lmsi=max(lmsi,tcg_ind_buf(ji,tcgl)); enddo
             else
              tcgl=tcgl-1
             endif
             jl=jl+1; jerr=lvit%next() !next subtensor
            enddo lloop
            if(jerr.eq.GFC_NO_MOVE) jerr=TEREC_SUCCESS
           endif
           lfinish=tcgl !lfinish: end offset of the left tensor descriptors
           !call print_tcg_buffer(lstart,lfinish,lrank) !debug
 !Right subtensors:
           rstart=tcgl+1; rmsi=0_INTL !start offset of the tensor descriptors and maximum subspace id
           call this%contr_ptrn%get_dim_symmetry(2,jnd,dim_restr,jerr)
           if(jerr.eq.TEREC_SUCCESS) then
            jl=0; jerr=rvit%reset()
  !Iterate over subtensors:
            rloop: do while(jerr.eq.GFC_SUCCESS)
   !Get subtensor header:
             jup=>rvit%get_value(jerr); if(jerr.ne.GFC_SUCCESS) exit rloop
             jtrp=>NULL(); select type(jup); class is(tens_rcrsv_t); jtrp=>jup; end select
             if(.not.associated(jtrp)) then; jerr=TEREC_OBJ_CORRUPTED; exit rloop; endif !trap
             jthp=>jtrp%get_header()
             if(.not.associated(jthp)) then; jerr=TEREC_OBJ_CORRUPTED; exit rloop; endif !trap
             call jthp%get_spec(sidx,rrank,jerr); if(jerr.ne.TEREC_SUCCESS) exit rloop
   !Append the subtensor multi-index (descriptor) into the sorting list (adjust dimension SAT level, if needed):
             approved=.TRUE.
             tcgl=tcgl+1; tcg_num_buf(tcgl)=int(jl,INTL) !subtensor number: [0..max]
             do ji=1,rrank
              ja=dim_restr(ji)
              if(ja.gt.0.and.ja.le.rrank) then
               jhsp=>ths(ji,2)%get_space(jerr); if(jerr.ne.TEREC_SUCCESS) exit rloop
               jc=jhsp%compare_subranges(sidx(ja),sidx(ji))
               if(jc.eq.CMP_ER) then; jerr=TEREC_ERROR; exit rloop; endif
               if(jc.eq.CMP_GT) then; approved=.FALSE.; exit; endif !`I assume LT index ordering for symmetric tensor dimensions in a tensor contraction
              endif
              if(adj(ji,2).gt.0) then !promotion to an ancestor SAT level is needed
               jhsp=>ths(ji,2)%get_space(jerr); if(jerr.ne.TEREC_SUCCESS) exit rloop
               tcg_ind_buf(ji,tcgl)=jhsp%get_ancestor_id(sidx(ji),adj(ji,2),jerr); if(jerr.ne.TEREC_SUCCESS) exit rloop
              else
               tcg_ind_buf(ji,tcgl)=sidx(ji)
              endif
             enddo
             if(approved) then
              do ji=1,rrank; rmsi=max(rmsi,tcg_ind_buf(ji,tcgl)); enddo
             else
              tcgl=tcgl-1
             endif
             jl=jl+1; jerr=rvit%next() !next subtensor
            enddo rloop
            if(jerr.eq.GFC_NO_MOVE) jerr=TEREC_SUCCESS
           endif
           rfinish=tcgl !rfinish: end offset of the right tensor descriptors
           !call print_tcg_buffer(rstart,rfinish,rrank) !debug
 !Destination subtensors:
           dstart=tcgl+1; dmsi=0_INTL !start offset of the tensor descriptors and maximum subspace id
           call this%contr_ptrn%get_dim_symmetry(0,jnd,dim_restr,jerr)
           if(jerr.eq.TEREC_SUCCESS) then
            jl=0; jerr=dvit%reset()
  !Iterate over subtensors:
            dloop: do while(jerr.eq.GFC_SUCCESS)
   !Get subtensor header:
             jup=>dvit%get_value(jerr); if(jerr.ne.GFC_SUCCESS) exit dloop
             jtrp=>NULL(); select type(jup); class is(tens_rcrsv_t); jtrp=>jup; end select
             if(.not.associated(jtrp)) then; jerr=TEREC_OBJ_CORRUPTED; exit dloop; endif !trap
             jthp=>jtrp%get_header()
             if(.not.associated(jthp)) then; jerr=TEREC_OBJ_CORRUPTED; exit dloop; endif !trap
             call jthp%get_spec(sidx,drank,jerr); if(jerr.ne.TEREC_SUCCESS) exit dloop
   !Append the subtensor multi-index (descriptor) into the sorting list (adjust dimension SAT level, if needed):
             approved=.TRUE.
             tcgl=tcgl+1; tcg_num_buf(tcgl)=int(jl,INTL) !subtensor number: [0..max]
             do ji=1,drank
              ja=dim_restr(ji)
              if(ja.gt.0.and.ja.le.drank) then
               jhsp=>ths(ji,0)%get_space(jerr); if(jerr.ne.TEREC_SUCCESS) exit dloop
               jc=jhsp%compare_subranges(sidx(ja),sidx(ji))
               if(jc.eq.CMP_ER) then; jerr=TEREC_ERROR; exit dloop; endif
               if(jc.eq.CMP_GT) then; approved=.FALSE.; exit; endif !`I assume LT index ordering for symmetric tensor dimensions in a tensor contraction
              endif
              if(adj(ji,0).gt.0) then !promotion to an ancestor SAT level is needed
               jhsp=>ths(ji,0)%get_space(jerr); if(jerr.ne.TEREC_SUCCESS) exit dloop
               tcg_ind_buf(ji,tcgl)=jhsp%get_ancestor_id(sidx(ji),adj(ji,0),jerr); if(jerr.ne.TEREC_SUCCESS) exit dloop
              else
               tcg_ind_buf(ji,tcgl)=sidx(ji)
              endif
             enddo
             if(approved) then
              do ji=1,drank; dmsi=max(dmsi,tcg_ind_buf(ji,tcgl)); enddo
             else
              tcgl=tcgl-1
             endif
             jl=jl+1; jerr=dvit%next() !next subtensor
            enddo dloop
            if(jerr.eq.GFC_NO_MOVE) jerr=TEREC_SUCCESS
           endif
           dfinish=tcgl !dfinish: end offset of the destination tensor descriptors
           !call print_tcg_buffer(dstart,dfinish,drank) !debug
           return
          end subroutine build_descriptors

          subroutine generate_subcontractions(jerr)
           implicit none
           integer(INTD), intent(out):: jerr
           integer(INTL), pointer, contiguous:: ext_buf(:),iv(:,:),v(:)
           integer(INTD):: j1,j2,jn,ji,jj,jif,jjf,jnc,jnu,jprm(1:MAX_TENSOR_RANK)
           integer(INTL):: jub

           jerr=TEREC_SUCCESS
 !Determine dimension order:
           jnc=0; jnu=0
           do ji=1,lrank
            jj=cptrn(ji)
            if(jj.lt.0) then !contracted index
             jj=-jj; jnc=jnc+1; ord(jnc,1)=ji; ord(jnc,2)=jj
            else !uncontracted index
             jnu=jnu+1; ord(nci+jnu,1)=ji; ord(jnu,0)=jj
            endif
           enddo
           jnc=0
           do ji=1,rrank
            jj=cptrn(lrank+ji)
            if(jj.gt.0) then !uncontracted index
             jnc=jnc+1; ord(nci+jnc,2)=ji; jnu=jnu+1; ord(jnu,0)=jj
            endif
           enddo
           !write(CONS_OUT,'("D dimension order:",32(1x,i2))') ord(1:drank,0) !debug: position in D
           !write(CONS_OUT,'("L dimension order:",32(1x,i2))') ord(1:lrank,1) !debug: N2O for L
           !write(CONS_OUT,'("R dimension order:",32(1x,i2))') ord(1:rrank,2) !debug: N2O for R
           !write(CONS_OUT,'("Current legnth of the TCG buffer = ",i9)') tcgl !debug
 !Sort tensor descriptors:
  !Use the rest of the TCG buffer as an external buffer:
           jub=(int(ubound(tcg_ind_buf,2),INTL)-tcgl)*int(size(tcg_ind_buf,1),INTL)
           ext_buf(1:jub)=>tcg_ind_buf(:,tcgl+1:) !`Does this introduce a temporary copy?
  !Sort left descriptors:
           if(lrank.gt.0) then
            iv(1:,1:)=>tcg_ind_buf(:,lstart:lfinish)
            v(1:)=>tcg_num_buf(lstart:lfinish)
            !write(CONS_OUT,'("Sorting segment ",i10,1x,i10)') lstart,lfinish !debug
            !do j1=1,lfinish-lstart+1; write(CONS_OUT,'(2x,i10,64(1x,i4))') v(j1),iv(1:lrank,j1); enddo !debug
            !write(CONS_OUT,'("Calling multord_i8e() ...")') !debug
            call multord_i8e(lrank,lmsi,int(lfinish-lstart+1,INTL),ord(1:lrank,1),iv,v,ext_buf)
            !call print_tcg_buffer(lstart,lfinish,lrank) !debug
           endif
  !Sort right descriptors:
           if(rrank.gt.0) then
            iv(1:,1:)=>tcg_ind_buf(:,rstart:rfinish)
            v(1:)=>tcg_num_buf(rstart:rfinish)
            !write(CONS_OUT,'("Sorting segment ",i10,1x,i10)') rstart,rfinish !debug
            !do j1=1,rfinish-rstart+1; write(CONS_OUT,'(2x,i10,64(1x,i4))') v(j1),iv(1:rrank,j1); enddo !debug
            !write(CONS_OUT,'("Calling multord_i8e() ...")') !debug
            call multord_i8e(rrank,rmsi,int(rfinish-rstart+1,INTL),ord(1:rrank,2),iv,v,ext_buf)
            !call print_tcg_buffer(rstart,rfinish,rrank) !debug
           endif
  !Sort destination descriptors:
           if(drank.gt.0) then
            iv(1:,1:)=>tcg_ind_buf(:,dstart:dfinish)
            v(1:)=>tcg_num_buf(dstart:dfinish)
            !write(CONS_OUT,'("Sorting segment ",i10,1x,i10)') dstart,dfinish !debug
            !do j1=1,dfinish-dstart+1; write(CONS_OUT,'(2x,i10,64(1x,i4))') v(j1),iv(1:drank,j1); enddo !debug
            !write(CONS_OUT,'("Calling multord_i8e() ...")') !debug
            call multord_i8e(drank,dmsi,int(dfinish-dstart+1,INTL),(/(ji,ji=1,drank)/),iv,v,ext_buf)
            !call print_tcg_buffer(dstart,dfinish,drank) !debug
           endif
  !Match contracted multi-indices:
           pstart=tcgl+1; ji=lstart; jj=rstart
           do while(ji.le.lfinish.and.jj.le.rfinish)
   !Check first match:
            if(nci.gt.0) then
             jnc=multindx_cmp(nci,tcg_ind_buf(ord(1:nci,1),ji),nci,tcg_ind_buf(ord(1:nci,2),jj))
            else
             jnc=0
            endif
            if(jnc.eq.0) then
             if(nci.gt.0) then
   !Determine the size of the left block:
              jif=ji+1
              do while(jif.le.lfinish)
               jnu=multindx_cmp(nci,tcg_ind_buf(ord(1:nci,1),ji),nci,tcg_ind_buf(ord(1:nci,1),jif))
               if(jnu.ne.0) exit
               jif=jif+1
              enddo
              jif=jif-1
   !Determine the size of the right block:
              jjf=jj+1
              do while(jjf.le.rfinish)
               jnu=multindx_cmp(nci,tcg_ind_buf(ord(1:nci,2),jj),nci,tcg_ind_buf(ord(1:nci,2),jjf))
               if(jnu.ne.0) exit
               jjf=jjf+1
              enddo
              jjf=jjf-1
             else
              jif=lfinish; jjf=rfinish
             endif
   !Take Cartesian product of the blocks:
             jn=0
             do j2=jj,jjf
              jn=tcg_num_buf(j2)*lsl
              if(drank.gt.0) then
               do j1=ji,jif
                tcgl=tcgl+1; tcg_num_buf(tcgl)=jn+tcg_num_buf(j1)
                tcg_ind_buf(1:drank,tcgl)=(/tcg_ind_buf(ord(nci+1:lrank,1),j1),tcg_ind_buf(ord(nci+1:rrank,2),j2)/)
               enddo
              else
               do j1=ji,jif
                tcgl=tcgl+1; tcg_num_buf(tcgl)=jn+tcg_num_buf(j1)
               enddo
              endif
             enddo
   !Proceed further:
             ji=jif+1; jj=jjf+1
            else
             if(jnc.lt.0) then; ji=ji+1; else; jj=jj+1; endif
            endif
           enddo
           pfinish=tcgl
           !call print_tcg_buffer(pstart,pfinish,drank) !debug
  !Filter with the destination multi-indices:
           if(drank.gt.0) then
   !Use the rest of the TCG buffer as an external buffer:
            jub=(int(ubound(tcg_ind_buf,2),INTL)-tcgl)*int(size(tcg_ind_buf,1),INTL)
            ext_buf(1:jub)=>tcg_ind_buf(:,tcgl+1:) !`Does this introduce a temporary copy?
   !Sort the Cartesian products:
            jprm(ord(1:drank,0))=(/(ji,ji=1,drank)/)
            iv=>tcg_ind_buf(:,pstart:pfinish)
            v=>tcg_num_buf(pstart:pfinish)
            call multord_i8e(drank,max(lmsi,rmsi),int(pfinish-pstart+1,INTL),jprm(1:drank),iv,v,ext_buf)
            !call print_tcg_buffer(pstart,pfinish,drank) !debug
   !Filter with the destination multi-index:
            ji=dstart; jj=pstart
            do while(ji.le.dfinish.and.jj.le.pfinish)
    !Check first match:
             jnc=multindx_cmp(drank,tcg_ind_buf(1:drank,ji),drank,tcg_ind_buf(jprm(1:drank),jj))
             if(jnc.eq.0) then
    !Determine the size of the block:
              jjf=jj+1
              do while(jjf.le.pfinish)
               jnu=multindx_cmp(drank,tcg_ind_buf(1:drank,jj),drank,tcg_ind_buf(1:drank,jjf))
               if(jnu.ne.0) exit
               jjf=jjf+1
              enddo
              jjf=jjf-1
    !Record subtensor contractions:
              call record_subcontractions(ji,jj,jjf,jerr); if(jerr.ne.TEREC_SUCCESS) return
              ji=ji+1; jj=jjf+1
             else
              if(jnc.lt.0) then; ji=ji+1; else; jj=jj+1; endif
             endif
            enddo
           else !scalar destination
            call record_subcontractions(dstart,pstart,pfinish,jerr)
           endif
           return
          end subroutine generate_subcontractions

          subroutine record_subcontractions(ds,ps,pf,jerr)
           implicit none
           integer(INTD), intent(in):: ds
           integer(INTD), intent(in):: ps
           integer(INTD), intent(in):: pf
           integer(INTD), intent(out):: jerr
           integer(INTD):: jdn,jln,jrn,ji
           type(tens_contraction_t):: tcontr
           type(tens_contraction_t), pointer:: jtcp
           class(tens_rcrsv_t), pointer:: dtrp,ltrp,rtrp
           class(*), pointer:: jup

           jerr=TEREC_SUCCESS
 !Associate destination tensor header:
           jdn=tcg_num_buf(ds) !destination subtensor number
           jup=>dvit%element_value(int(jdn,INTL),jerr); if(jerr.ne.GFC_SUCCESS) return
           dtrp=>NULL(); select type(jup); class is(tens_rcrsv_t); dtrp=>jup; end select
           if(.not.associated(dtrp)) then; jerr=TEREC_ERROR; return; endif !trap
           do ji=ps,pf
            jrn=tcg_num_buf(ji)/lsl     !right subtensor number
            jln=tcg_num_buf(ji)-jrn*lsl !left subtensor number
            !write(CONS_OUT,'("New subcontraction: ",i6," = ",i6," * ",i6)') jdn,jln,jrn !debug
 !Associate left tensor header:
            jup=>lvit%element_value(int(jln,INTL),jerr); if(jerr.ne.GFC_SUCCESS) exit
            ltrp=>NULL(); select type(jup); class is(tens_rcrsv_t); ltrp=>jup; end select
            if(.not.associated(ltrp)) then; jerr=TEREC_ERROR; exit; endif !trap
 !Associate right tensor header:
            jup=>rvit%element_value(int(jrn,INTL),jerr); if(jerr.ne.GFC_SUCCESS) exit
            rtrp=>NULL(); select type(jup); class is(tens_rcrsv_t); rtrp=>jup; end select
            if(.not.associated(rtrp)) then; jerr=TEREC_ERROR; exit; endif !trap
 !Construct and record a subcontraction:
            jerr=slit%append(tcontr); if(jerr.ne.GFC_SUCCESS) exit
            jerr=slit%reset_back(); if(jerr.ne.GFC_SUCCESS) exit
            jup=>slit%get_value(jerr); if(jerr.ne.GFC_SUCCESS) exit
            jtcp=>NULL(); select type(jup); class is(tens_contraction_t); jtcp=>jup; end select
            if(.not.associated(jtcp)) then; jerr=TEREC_ERROR; exit; endif !trap
            call jtcp%import_replace(this,jerr,dtrp,ltrp,rtrp); if(jerr.ne.TEREC_SUCCESS) exit
            nsub=nsub+1; jtcp=>NULL()
           enddo
           return
          end subroutine record_subcontractions

        end subroutine TensContractionSplit
!------------------------------------------------------------------
        subroutine TensContractionPrintIt(this,ierr,dev_id,nspaces)
!Prints the tensor contraction info.
         implicit none
         class(tens_contraction_t), intent(in):: this  !in: tensor contraction
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD), intent(in), optional:: dev_id  !in: output device (defaults to screen)
         integer(INTD), intent(in), optional:: nspaces !in: number of leading spaces
         integer(INTD):: errc,devo,nsp,n,i,j
         class(tens_rcrsv_t), pointer:: trp

         devo=6; if(present(dev_id)) devo=dev_id
         nsp=0; if(present(nspaces)) nsp=nspaces
         do i=1,nsp; write(devo,'(" ")',ADVANCE='NO'); enddo
         write(devo,'("TENSOR CONTRACTION{")')
         if(this%is_set(errc)) then
          call this%contr_ptrn%print_it(j,devo,nsp+1); if(errc.eq.TEREC_SUCCESS.and.j.ne.TEREC_SUCCESS) errc=j
          do i=1,nsp+1; write(devo,'(" ")',ADVANCE='NO'); enddo
          write(devo,'("Complex Prefactor = (",D21.14,",",D21.14,")")') this%alpha
          n=this%get_num_args(j); if(errc.eq.TEREC_SUCCESS.and.j.ne.TEREC_SUCCESS) errc=j
          do i=0,n-1
           trp=>NULL(); trp=>this%get_argument(i,j)
           if(errc.eq.TEREC_SUCCESS.and.j.ne.TEREC_SUCCESS) errc=j
           call trp%print_it(j,devo,nsp+1)
           if(errc.eq.TEREC_SUCCESS.and.j.ne.TEREC_SUCCESS) errc=j
          enddo
         else
          do i=1,nsp+1; write(devo,'(" ")',ADVANCE='NO'); enddo
          write(devo,'("Empty/invalid tensor contraction")')
         endif
         do i=1,nsp; write(devo,'(" ")',ADVANCE='NO'); enddo
         write(devo,'("}")')
         if(present(ierr)) ierr=errc
         return
        end subroutine TensContractionPrintIt

       end module tensor_recursive
!==================================
       module tensor_recursive_test
        use tensor_algebra
        use gfc_base
        use gfc_list
        use gfc_vector
        use subspaces
        use tensor_recursive
        use distributed, only: DataDescr_t,data_descr_rnd_
        implicit none
        private
        public test_tensor_recursive
!GLOBAL DATA:
        type(vector_t), private:: subtensor_storage !persistent storage for subtensors

       contains
!---------------------------------------------
        subroutine test_tensor_recursive(ierr)
         implicit none
         integer(INTD), intent(out):: ierr
         logical, parameter:: FTEST_TENS_SIGNATURE=.TRUE.
         logical, parameter:: FTEST_TENS_SHAPE=.TRUE.
         logical, parameter:: FTEST_TENS_HEADER=.TRUE.
         logical, parameter:: FTEST_TENS_SIMPLE_PART=.TRUE.
         logical, parameter:: FTEST_TENS_RCRSV=.TRUE.
         logical, parameter:: FTEST_TENS_CONTRACTION=.TRUE.
         logical, parameter:: FTEST_CMP_INTEGERS=.TRUE.

         if(FTEST_TENS_SIGNATURE) then
          write(*,'("Testing class tens_signature_t ... ")',ADVANCE='NO')
          call test_tens_signature(ierr)
          if(ierr.eq.0) then; write(*,'("PASSED")'); else; write(*,'("FAILED: Error ",i11)') ierr; return; endif
         endif
         if(FTEST_TENS_SHAPE) then
          write(*,'("Testing class tens_shape_t ... ")',ADVANCE='NO')
          call test_tens_shape(ierr)
          if(ierr.eq.0) then; write(*,'("PASSED")'); else; write(*,'("FAILED: Error ",i11)') ierr; return; endif
         endif
         if(FTEST_TENS_HEADER) then
          write(*,'("Testing class tens_header_t ... ")',ADVANCE='NO')
          call test_tens_header(ierr)
          if(ierr.eq.0) then; write(*,'("PASSED")'); else; write(*,'("FAILED: Error ",i11)') ierr; return; endif
         endif
         if(FTEST_TENS_SIMPLE_PART) then
          write(*,'("Testing class tens_simple_part_t ... ")',ADVANCE='NO')
          call test_tens_simple_part(ierr)
          if(ierr.eq.0) then; write(*,'("PASSED")'); else; write(*,'("FAILED: Error ",i11)') ierr; return; endif
         endif
         if(FTEST_TENS_RCRSV) then
          write(*,'("Testing class tens_rcrsv_t ... ")',ADVANCE='NO')
          call test_tens_rcrsv(ierr)
          if(ierr.eq.0) then; write(*,'("PASSED")'); else; write(*,'("FAILED: Error ",i11)') ierr; return; endif
         endif
         if(FTEST_TENS_CONTRACTION) then
          write(*,'("Testing class tens_contraction_t ... ")',ADVANCE='NO')
          call test_tens_contraction(ierr)
          if(ierr.eq.0) then; write(*,'("PASSED")'); else; write(*,'("FAILED: Error ",i11)') ierr; return; endif
         endif
         if(FTEST_CMP_INTEGERS) then
          write(*,'("Testing generic integer comparator ... ")',ADVANCE='NO')
          call test_cmp_integers(ierr)
          if(ierr.eq.0) then; write(*,'("PASSED")'); else; write(*,'("FAILED: Error ",i11)') ierr; return; endif
         endif
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
         call tsigna%tens_signature_ctor(ierr,(/3_INTL,4_INTL,2_INTL/),'Tensor')
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
         character(32):: tens_name

         ierr=0
         n=6; dims(1:n)=(/128_INTL,64_INTL,256_INTL,64_INTL,128_INTL,64_INTL/)
         m=2; grps(1:n)=(/1,2,0,2,1,2/); grp_spec(1:m)=(/TEREC_IND_RESTR_LT,TEREC_IND_RESTR_GE/)
         call thead%tens_header_ctor(ierr,'Tensor',(/1_INTL,2_INTL,3_INTL,2_INTL,1_INTL,2_INTL/))
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
         call thead%tens_header_ctor(ierr,'Tensor',(/1_INTL,2_INTL,3_INTL,2_INTL,1_INTL,2_INTL/))
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
         integer(INTD), parameter:: tens_rank=4 !tensor rank
         class(h_space_t), pointer:: hspace     !hierarchical representation of the vector space
         integer(INTL):: spcx(1:MAX_TENSOR_RANK),dims(1:MAX_TENSOR_RANK),space_id,max_res
         integer(INTD):: dimg(1:MAX_TENSOR_RANK),grps(1:MAX_TENSOR_RANK),num_subtensors,hsid,j
         type(tens_header_t), pointer:: thp
         type(tens_rcrsv_t):: tensor
         type(list_bi_t):: subtensors
         type(list_iter_t):: lit
         class(subspace_t), pointer:: ssp
         type(tens_descr_t):: tdescr
         class(*), pointer:: up

 !Build a hierarchical representation for a test vector space:
         hspace=>NULL()
         !write(*,'("#DEBUG(test_tens_rcrsv): Building test H-space ... ")') !debug
         hsid=build_test_hspace('TestSpace0',ierr,hspace); if(ierr.ne.TEREC_SUCCESS) then; ierr=1; return; endif
         !write(*,'("#DEBUG(test_tens_rcrsv): H-space has been built.")') !debug
 !Create the full tensor (over the full space):
  !Get full space id and its max resolution:
         space_id=hspace%get_root_id(ierr); if(ierr.ne.0) then; ierr=3; return; endif
         ssp=>hspace%get_subspace(space_id,ierr); if(ierr.ne.0) then; ierr=4; return; endif
         if(.not.associated(ssp)) then; ierr=5; return; endif
         max_res=ssp%get_max_resolution(ierr); if(ierr.ne.0) then; ierr=6; return; endif
         !write(*,*) 'Space ID = ',space_id,': Max resolution = ',max_res !debug
  !Create a tensor over the full space:
         spcx(1:tens_rank)=space_id
         dims(1:tens_rank)=max_res
         dimg(1:tens_rank)=(/1,1,2,2/); grps(1:2)=(/TEREC_IND_RESTR_LT,TEREC_IND_RESTR_GT/)
         call tensor%tens_rcrsv_ctor('T2',spcx(1:tens_rank),(/(hsid,j=1,tens_rank)/),ierr,&
                                    &dims(1:tens_rank),dimg(1:tens_rank),grps(1:2))
         if(ierr.ne.0) then; ierr=7; return; endif
 !Check the tensor descriptor:
         thp=>tensor%get_header(ierr); if(ierr.ne.0) then; ierr=8; return; endif
         call tensor%add_subtensor(thp,ierr); thp=>NULL(); if(ierr.ne.0) then; ierr=9; return; endif
         call tensor%set_layout(TEREC_LAY_FDIMS,R8,ierr); if(ierr.ne.0) then; ierr=10; return; endif
         call tensor%set_location(data_descr_rnd_,ierr); if(ierr.ne.0) then; ierr=11; return; endif
         tdescr=tensor%get_descriptor(ierr); if(ierr.ne.0) then; ierr=12; return; endif
         !call tdescr%print_it() !debug
 !Split the tensor into subtensors:
         call tensor%split((/1,2,3,4/),subtensors,ierr,num_subtensors,.TRUE.); if(ierr.ne.0) then; ierr=13; return; endif
         !write(*,*) 'Number of subtensors generated = ',num_subtensors !debug
         ierr=lit%init(subtensors); if(ierr.ne.0) then; ierr=14; return; endif
         !ierr=lit%scanp(action_f=print_tens_header_f); if(ierr.eq.GFC_IT_DONE) ierr=lit%reset() !debug
         !if(ierr.ne.0) then; ierr=15; return; endif !debug
         ierr=lit%delete_all(); if(ierr.ne.0) then; ierr=16; return; endif
         ierr=lit%release(); if(ierr.ne.0) then; ierr=17; return; endif
         return
        end subroutine test_tens_rcrsv
!------------------------------------------------------------------------------
        function tens_split_func(tensor,subtensors,num_subtensors) result(ierr)
!Testing only: Returns a list of subtensors for a given tensor.
         implicit none
         integer(INTD):: ierr                        !out: error code
         class(tens_rcrsv_t), intent(in):: tensor    !in: tensor
         type(vector_t), intent(inout):: subtensors  !out: vector of subtensors
         integer(INTD), intent(out):: num_subtensors !out: number of generated subtensors
         integer(INTD):: i,nd
         integer(INTL):: first,last
         type(vector_iter_t):: vit,rvit
         class(*), pointer:: up

         num_subtensors=0
         if(tensor%is_set(ierr,num_dims=nd)) then
          if(ierr.eq.TEREC_SUCCESS) then
           ierr=vit%init(subtensor_storage)
           if(ierr.eq.GFC_SUCCESS) then
            first=vit%get_length()
            call tensor%split((/(i,i=1,nd)/),subtensor_storage,ierr,num_subtensors) !subtensor_storage is a global vector
            if(ierr.eq.TEREC_SUCCESS) then
             ierr=vit%reset() !update iterator status
             if(ierr.eq.GFC_SUCCESS) then
              last=vit%get_length()-1_INTL
              if(first.le.last) then
               ierr=vit%move_to(first)
               if(ierr.eq.GFC_SUCCESS) then
                ierr=rvit%init(subtensors)
                if(ierr.eq.GFC_SUCCESS) then
                 do
                  up=>NULL(); up=>vit%get_value(); if(.not.associated(up)) then; ierr=TEREC_OBJ_CORRUPTED; exit; endif
                  ierr=rvit%append(up,assoc_only=.TRUE.); if(ierr.ne.GFC_SUCCESS) exit
                  ierr=vit%next(); if(ierr.ne.GFC_SUCCESS) exit
                 enddo
                 if(ierr.eq.GFC_NO_MOVE) ierr=GFC_SUCCESS
                 i=rvit%release(); if(i.ne.GFC_SUCCESS.and.ierr.eq.TEREC_SUCCESS) ierr=TEREC_ERROR
                endif
               endif
              else
               ierr=TEREC_ERROR
              endif
             endif
            endif
            i=vit%release(); if(i.ne.GFC_SUCCESS.and.ierr.eq.TEREC_SUCCESS) ierr=TEREC_ERROR
           endif
          endif
         else
          ierr=TEREC_INVALID_REQUEST
         endif
         return
        end function tens_split_func
!---------------------------------------------
        subroutine test_tens_contraction(ierr)
         implicit none
         integer(INTD), intent(out):: ierr
         !-------------------------------------
         integer(INTD), parameter:: tens_rank=4 !tensor rank
         !-------------------------------------
         class(h_space_t), pointer:: hspace     !hierarchical representation of the vector space
         integer(INTL):: spcx(1:MAX_TENSOR_RANK),dims(1:MAX_TENSOR_RANK),space_id,max_res
         integer(INTD):: dimg(1:MAX_TENSOR_RANK),grps(1:MAX_TENSOR_RANK)
         integer(INTD):: j,hsid,num_subcontractions
         class(subspace_t), pointer:: ssp
         type(tens_rcrsv_t):: dtens,ltens,rtens,stens
         type(tens_contraction_t):: tens_contr
         type(tens_contraction_t), pointer:: subcontr_p
         type(list_bi_t):: subcontractions
         type(list_iter_t):: lit
         type(vector_iter_t):: vit
         class(*), pointer:: up

!Build a hierarchical representation for a test vector space:
         hsid=build_test_hspace('TestSpace1',ierr,hspace); if(ierr.ne.TEREC_SUCCESS) then; ierr=1; return; endif

!Create tensor arguments (over the full space):
 !Get full space id and its max resolution:
         space_id=hspace%get_root_id(ierr); if(ierr.ne.0) then; ierr=3; return; endif
         ssp=>hspace%get_subspace(space_id,ierr); if(ierr.ne.0) then; ierr=4; return; endif
         if(.not.associated(ssp)) then; ierr=5; return; endif
         max_res=ssp%get_max_resolution(ierr); if(ierr.ne.0) then; ierr=6; return; endif
         !write(*,*) 'Space ID = ',space_id,': Max resolution = ',max_res !debug
 !Create a tensor over the full space:
         spcx(1:tens_rank)=space_id
         dims(1:tens_rank)=max_res
         call dtens%tens_rcrsv_ctor('Z2',spcx(1:tens_rank),(/(hsid,j=1,tens_rank)/),ierr,&
                                   &dims(1:tens_rank))
         if(ierr.ne.0) then; ierr=7; return; endif
 !Create a tensor over the full space:
         spcx(1:tens_rank)=space_id
         dims(1:tens_rank)=max_res
         dimg(1:tens_rank)=(/1,1,2,2/); grps(1:2)=(/TEREC_IND_RESTR_LT,TEREC_IND_RESTR_LT/)
         call ltens%tens_rcrsv_ctor('H2',spcx(1:tens_rank),(/(hsid,j=1,tens_rank)/),ierr,&
                                   &dims(1:tens_rank),dimg(1:tens_rank),grps(1:2))
         if(ierr.ne.0) then; ierr=8; return; endif
 !Create a tensor over the full space:
         spcx(1:tens_rank)=space_id
         dims(1:tens_rank)=max_res
         dimg(1:tens_rank)=(/1,1,2,2/); grps(1:2)=(/TEREC_IND_RESTR_LT,TEREC_IND_RESTR_LT/)
         call rtens%tens_rcrsv_ctor('T2',spcx(1:tens_rank),(/(hsid,j=1,tens_rank)/),ierr,&
                                   &dims(1:tens_rank),dimg(1:tens_rank),grps(1:2))
         if(ierr.ne.0) then; ierr=9; return; endif
  !Create a scalar tensor:
         call stens%tens_rcrsv_ctor('dE',spcx(1:0),(/(hsid,j=1,0)/),ierr)
         if(ierr.ne.0) then; ierr=10; return; endif
#if 1
!Tensor contraction 1:
 !Create the full tensor contraction specification for Z2(a,b,i,j)+=H2(i,k,a,c)*T2(b,c,j,k): [a<b] [i<j]:
         call tens_contr%clean(ierr); if(ierr.ne.0) then; ierr=11; return; endif
  !Set tensor contraction arguments:
         call tens_contr%set_argument(dtens,ierr); if(ierr.ne.0) then; ierr=12; return; endif
         call tens_contr%set_argument(ltens,ierr); if(ierr.ne.0) then; ierr=13; return; endif
         call tens_contr%set_argument(rtens,ierr); if(ierr.ne.0) then; ierr=14; return; endif
         !print *,'Tensor contraction (args_set,fully set): ',tens_contr%args_full(),tens_contr%is_set() !debug
  !Set the tensor contraction pattern:
         call tens_contr%set_contr_ptrn((/3,-4,1,-2, 2,-4,4,-2/),ierr); if(ierr.ne.0) then; ierr=15; return; endif
  !Impose additional symmetries, if needed:
         call tens_contr%set_operl_symm(0,(/1,2/),ierr); if(ierr.ne.0) then; ierr=16; return; endif ![a<b]
         call tens_contr%set_operl_symm(0,(/3,4/),ierr); if(ierr.ne.0) then; ierr=17; return; endif ![i<j]
         !call tens_contr%print_it() !debug
 !Split the tensor contraction into a list of subtensor contractions:
         call tens_contr%split(tens_split_func,subcontractions,ierr,num_subcontractions)
         if(ierr.ne.0) then; ierr=18; return; endif
         !write(*,*) 'Number of subtensor contractions = ',num_subcontractions !debug
         if(num_subcontractions.ne.16) then; ierr=19; return; endif
         ierr=lit%init(subcontractions); if(ierr.ne.GFC_SUCCESS) then; ierr=20; return; endif
  !Check the GFC list (debug):
         ierr=lit%reset(); if(ierr.ne.GFC_SUCCESS) then; ierr=21; return; endif
         ierr=lit%scanp(); if(ierr.ne.GFC_IT_DONE) then; ierr=22; return; endif
         if(lit%total_count(ierr).ne.num_subcontractions) then; ierr=23; return; endif
         ierr=lit%reset(); if(ierr.ne.GFC_SUCCESS) then; ierr=21; return; endif
 !Print subcontractions (debug):
#if 0
         do j=1,num_subcontractions
          write(*,'("Subtensor contraction ",i7,":")') j
          up=>lit%get_value(ierr); if(ierr.ne.GFC_SUCCESS) then; ierr=24; return; endif
          subcontr_p=>NULL(); select type(up); class is(tens_contraction_t); subcontr_p=>up; end select
          if(.not.associated(subcontr_p)) then; ierr=25; return; endif
          call subcontr_p%print_it(ierr,nspaces=1); if(ierr.ne.TEREC_SUCCESS) then; ierr=26; return; endif
          ierr=lit%next(); if(ierr.ne.GFC_SUCCESS.and.j.lt.num_subcontractions) then; ierr=27; return; endif
         enddo
         if(ierr.ne.GFC_NO_MOVE) then; ierr=28; return; endif
#endif
 !Destroy subcontractions:
         ierr=lit%delete_all(); if(ierr.ne.GFC_SUCCESS) then; ierr=29; return; endif
         ierr=lit%release(); if(ierr.ne.GFC_SUCCESS) then; ierr=30; return; endif
#endif
#if 1
!Tensor contraction 2:
 !Create the full tensor contraction specification for dE()+=H2(i,k,a,c)*T2(a,c,i,k):
         call tens_contr%clean(ierr); if(ierr.ne.0) then; ierr=31; return; endif
  !Set tensor contraction arguments:
         call tens_contr%set_argument(stens,ierr); if(ierr.ne.0) then; ierr=32; return; endif
         call tens_contr%set_argument(ltens,ierr); if(ierr.ne.0) then; ierr=33; return; endif
         call tens_contr%set_argument(rtens,ierr); if(ierr.ne.0) then; ierr=34; return; endif
         !print *,'Tensor contraction (args_set,fully set): ',tens_contr%args_full(),tens_contr%is_set() !debug
  !Set the tensor contraction pattern:
         call tens_contr%set_contr_ptrn((/-3,-4,-1,-2, -3,-4,-1,-2/),ierr); if(ierr.ne.0) then; ierr=35; return; endif
  !Impose additional symmetries, if needed:
         !call tens_contr%print_it() !debug
 !Split the tensor contraction into a list of subtensor contractions:
         call tens_contr%split(tens_split_func,subcontractions,ierr,num_subcontractions)
         if(ierr.ne.0) then; ierr=36; return; endif
         !write(*,*) 'Number of subtensor contractions = ',num_subcontractions !debug
         if(num_subcontractions.ne.9) then; ierr=37; return; endif
         ierr=lit%init(subcontractions); if(ierr.ne.GFC_SUCCESS) then; ierr=38; return; endif
  !Check the GFC list (debug):
         ierr=lit%reset(); if(ierr.ne.GFC_SUCCESS) then; ierr=39; return; endif
         ierr=lit%scanp(); if(ierr.ne.GFC_IT_DONE) then; ierr=40; return; endif
         if(lit%total_count(ierr).ne.num_subcontractions) then; ierr=41; return; endif
         ierr=lit%reset(); if(ierr.ne.GFC_SUCCESS) then; ierr=39; return; endif
 !Print subcontractions (debug):
#if 0
         do j=1,num_subcontractions
          write(*,'("Subtensor contraction ",i7,":")') j
          up=>lit%get_value(ierr); if(ierr.ne.GFC_SUCCESS) then; ierr=42; return; endif
          subcontr_p=>NULL(); select type(up); class is(tens_contraction_t); subcontr_p=>up; end select
          if(.not.associated(subcontr_p)) then; ierr=43; return; endif
          call subcontr_p%print_it(ierr,nspaces=1); if(ierr.ne.TEREC_SUCCESS) then; ierr=44; return; endif
          ierr=lit%next(); if(ierr.ne.GFC_SUCCESS.and.j.lt.num_subcontractions) then; ierr=45; return; endif
         enddo
         if(ierr.ne.GFC_NO_MOVE) then; ierr=46; return; endif
#endif
 !Destroy subcontractions:
         ierr=lit%delete_all(); if(ierr.ne.GFC_SUCCESS) then; ierr=47; return; endif
         ierr=lit%release(); if(ierr.ne.GFC_SUCCESS) then; ierr=48; return; endif
#endif

 !Release global resources:
         ierr=vit%init(subtensor_storage); if(ierr.ne.GFC_SUCCESS) then; ierr=100; return; endif
         ierr=vit%delete_all(); if(ierr.ne.GFC_SUCCESS) then; ierr=101; return; endif
         ierr=vit%release(); if(ierr.ne.GFC_SUCCESS) then; ierr=102; return; endif
         return
        end subroutine test_tens_contraction
!-----------------------------------------
        subroutine test_cmp_integers(ierr)
         implicit none
         integer(INTD), intent(out):: ierr
         integer(INTD), parameter:: num_times=1000000
         integer(INTD):: i,j,m,n
         integer(INTD), allocatable:: i1(:)
         integer(INTL), allocatable:: i2(:)
         real(8), allocatable:: rn(:)
         real(8):: tm1,tm2

         ierr=0
         allocate(rn(num_times),i1(num_times),i2(num_times))
         call random_number(rn); rn(:)=rn(:)*real(num_times,8)
         do i=1,num_times; i1(i)=int(rn(i),INTD); enddo
         call random_number(rn); rn(:)=rn(:)*real(num_times,8)
         do i=1,num_times; i2(i)=int(rn(i),INTD); enddo
 !Generic cross-kind comparison:
         m=0; tm2=thread_wtime()
         do i=1,num_times
          j=cmp_integers(i1(i),i2(i))
          if(j.eq.CMP_LT) then; m=m-1; elseif(j.eq.CMP_GT) then; m=m+1; endif
         enddo
         tm2=thread_wtime(tm2); write(*,'(1x,F7.4)',ADVANCE='NO') tm2
 !Regular cross-kind comparison:
         n=0; tm1=thread_wtime()
         do i=1,num_times
          if(i1(i).lt.i2(i)) then; n=n-1; elseif(i1(i).gt.i2(i)) then; n=n+1; endif
         enddo
         tm1=thread_wtime(tm1); write(*,'(1x,F7.4)',ADVANCE='NO') tm1
 !Print slowdown:
         write(*,'(1x,"Slowdown = ",F8.1,"X",1x)',ADVANCE='NO') tm2/tm1
         if(m.ne.n) ierr=1
         deallocate(rn,i1,i2)
         return
        end subroutine test_cmp_integers

       end module tensor_recursive_test
