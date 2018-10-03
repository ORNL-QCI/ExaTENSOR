!ExaTENSOR: Infrastructure for a recursive adaptive vector space decomposition
!and hierarchical vector space representation.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2018/10/03

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
! # Input: A vector space of dimension N with a specific basis.
! # Subspacing: Each original basis function is a subspace.
!   Existing subspaces are aggregated together to form larger subspaces,
!   up to the full space, which can be represented by the subspace
!   aggregation tree (SAT). In this process, the subspace dimension,
!   called "max_resolution", is a sum of the dimensions of the children
!   subspaces, and the subspace basis is a direct sum of the bases
!   from the children subspaces. Additionally, each subspace may
!   define one or more reduced basis sets of lower dimension/quality.
!   The reduced basis sets can be defined in multiple ways:
!    a) In an explicit form: Each basis function is specified analytically;
!    b) In a superposition form: Each basis function is represented as a
!       linear combination of the original basis vectors.
! # Subspace basis categories:
!    a) Abstract: Only the basis function kind is specified;
!    b) Real space supported: A support in the real space is specified as well;
!    c) Symmetric: A symmetry is specified as well.
!   Without any further specialization, a basis function is assumed
!   to be some abstract basis function of a specific basis set kind.
! # Output: The subspace aggregation tree is constructed in which each
!   subspace has a unique non-negative ID. A subspace overlap matrix
!   shows the overlap of the subspace supports.
! # Use:
!    1) Construct an empty subspace_basis_t object of specified dimension.
!    2) OPTIONAL: Define an array of specific basis functions (basis_func_supp_t).
!    3) OPTIONAL: Define an array of symmetries for the basis functions.
!    4) Set up each basis function in the subspace_basis_t object using
!       the two arrays created above.
!    5) Finalize the subspace_basis_t object via the .finalize() method.
!    6) Construct the hierarchical vector space object (h_space_t) using
!       the subspace_basis_t object.

       module subspaces
        use dil_basic
        use stsubs
        use gfc_base
        use gfc_list
        use gfc_vec_tree
        implicit none
        private
!PARAMETERS:
 !Output:
        integer, private:: CONS_OUT=6     !default output for this module
        integer, private:: DEBUG=0        !debugging mode level (0:none)
        logical, private:: VERBOSE=.TRUE. !verbosity for errors
 !Basis set kind:
        integer(INTD), parameter, public:: BASIS_NONE=0      !no basis set
        integer(INTD), parameter, public:: BASIS_ABSTRACT=1  !abstract basis set
        integer(INTD), parameter, public:: BASIS_GAUSS=2     !Gaussian basis set
        integer(INTD), parameter, public:: BASIS_SLATER=3    !Slater basis set
        integer(INTD), parameter, public:: BASIS_PLANEWAVE=4 !planewave basis set
        integer(INTD), parameter, public:: BASIS_HARMONIC=5  !harmonic basis set
        integer(INTD), parameter, public:: BASIS_POLYNOM=6   !polynomial basis set
        integer(INTD), parameter, public:: BASIS_WAVELET=7   !wavelet basis set
 !Symmetry:
        integer(INTD), parameter, public:: SYMMETRY_NONE=-1  !no symmetry
 !Subspace relationship:
        integer(INTD), parameter, public:: SUBSPACE_SAME=CMP_EQ      !subspace 1 and subspace 2 are the same
        integer(INTD), parameter, public:: SUBSPACE_LESS=CMP_LT      !subspace 1 is on the left of subspace 2
        integer(INTD), parameter, public:: SUBSPACE_GREATER=CMP_GT   !subspace 1 is on the right of subspace 2
        integer(INTD), parameter, public:: SUBSPACE_PARENT=CMP_CN    !subspace 1 is a parent of subspace 2
        integer(INTD), parameter, public:: SUBSPACE_CHILD=CMP_IN     !subspace 1 is a child of subspace 2
        integer(INTD), parameter, public:: SUBSPACE_SIBLING=CMP_OV   !subspace 1 and subspace 2 have the same parent
        integer(INTD), parameter, public:: SUBSPACE_UNRELATED=CMP_NC !subspace 1 and subspace 2 are unrelated
!TYPES:
 !Real space vector:
        type, public:: real_vec_t
         integer(INTL), private:: num_dim=0      !number of dimensions (0 means empty)
         real(8), allocatable, public:: coord(:) !components of the real space vector
         contains
          procedure, private:: RealVecCtor             !real space vector ctor
          generic, public:: real_vec_ctor=>RealVecCtor
          procedure, public:: dimsn=>RealVecDimsn      !returns dimension of the vector
          procedure, public:: norm2=>RealVecNorm2      !returns the 2-norm of the vector
          procedure, public:: scale=>RealVecScale      !vector scaling by a scalar
          procedure, public:: add=>RealVecAdd          !computes the sum of two tensors
          procedure, public:: average=>RealVecAverage  !computes an average of two real space vectors
          final:: real_vec_dtor                        !real space vector dtor
        end type real_vec_t
 !Real space 1d range  = segment[min:max]:
        type, public:: range1d_t
         real(8), private:: min_coord=0d0 !minimum coordinate (lower bound)
         real(8), private:: max_coord=0d0 !maximum coordinate (upper bound)
         contains
          procedure, public:: set=>Range1dSet                !sets the real range (ctor)
          procedure, public:: lower_bound=>Range1dLowerBound !returns the real range lower bound
          procedure, public:: upper_bound=>Range1dUpperBound !returns the real range upper bound
          procedure, public:: length=>Range1dLength          !returns the real range length = (upper - lower)
          procedure, public:: overlap=>Range1dOverlap        !returns the overlap of two real ranges
          procedure, public:: union=>Range1dUnion            !returns the minimal real range containing two given real ranges
          procedure, public:: split=>Range1dSplit            !splits the real range
        end type range1d_t
 !Integer semi-interval (min:max] = Integer range [min+1:max]:
        type, public:: seg_int_t
         integer(INTL), private:: min_coord=0 !minimum coordinate (lower bound): does not belong to the integer range
         integer(INTL), private:: max_coord=0 !maximum coordinate (upper bound): belongs to the integer range
         contains
          procedure, public:: set=>SegIntSet                !sets the integer range (ctor)
          procedure, public:: lower_bound=>SegIntLowerBound !returns the integer range lower bound
          procedure, public:: upper_bound=>SegIntUpperBound !returns the integer range upper bound
          procedure, public:: length=>SegIntLength          !returns the integer range length = (upper - lower)
          procedure, public:: overlap=>SegIntOverlap        !returns the overlap of two integer ranges
          procedure, public:: union=>SegIntUnion            !returns the minimal integer range containing two given integer ranges
          procedure, public:: split=>SegIntSplit            !splits the integer range
          procedure, public:: print_range=>SegIntPrintRange !prints the integer range
        end type seg_int_t
 !Real space rectangular hypercube (orthotope):
        type, public:: orthotope_t
         integer(INTL), private:: num_dim=0                    !number of dimensions
         type(range1d_t), allocatable, private:: extent(:)     !extent of each dimension (min,max)
         contains
          procedure, private:: OrthotopeCtor                   !orthotope ctor
          generic, public:: orthotope_ctor=>OrthotopeCtor
          procedure, public:: dimsn=>OrthotopeDimsn            !returns the real space dimension orthotope resides in
          procedure, public:: set_extent=>OrthotopeSetExtent   !sets the extent along a specific dimension
          procedure, public:: get_extent=>OrthotopeGetExtent   !returns the extent of a specific dimension
          procedure, public:: lower_bound=>OrthotopeLowerBound !returns the lower bound of a specific extent
          procedure, public:: upper_bound=>OrthotopeUpperBound !returns the upper bound of a specific extent
          procedure, public:: length=>OrthotopeLength          !returns the length along a specific extent
          procedure, public:: volume=>OrthotopeVolume          !returns the volume of the orthotope
          procedure, public:: overlap=>OrthotopeOverlap        !returns the overlap of two orthotopes
          procedure, public:: union=>OrthotopeUnion            !returns the minimal orthotope containing two given orthotopes
          final:: orthotope_dtor                               !orthotope dtor
        end type orthotope_t
 !Symmetry (abstract):
        type, abstract, public:: symmetry_t
         contains
          procedure(symm_compare_i), deferred, public:: compare !compares two symmetries
          procedure(symm_combine_i), deferred, public:: combine !combines two symmetries
        end type symmetry_t
 !Generic symmetry (color):
        type, extends(symmetry_t), public:: color_symmetry_t
         integer(INTD), private:: color=0                          !color (-inf;+inf)
         contains
          procedure, private:: ColorSymmetryCtor                   !ctor
          generic, public:: color_symmetry_ctor=>ColorSymmetryCtor
          procedure, public:: compare=>ColorSymmetryCompare        !compares two color symmetries
          procedure, public:: combine=>ColorSymmetryCombine        !combines two color symmetries
          final:: color_symmetry_dtor
        end type color_symmetry_t
 !Spherical symmetry (orbital momentum):
        type, extends(symmetry_t), public:: spher_symmetry_t
         integer(INTD), private:: orb_moment=SYMMETRY_NONE   !total orbital momentum: [0,1,2,...)
         integer(INTD), private:: orb_z_proj=0               !Z-axis projection of the total orbital momentum: [-L...+L]
         contains
          procedure, private:: SpherSymmetryCtor                            !constructor
          generic, public:: spher_symmetry_ctor=>SpherSymmetryCtor
          procedure, public:: get_orb_momentum=>SpherSymmetryGetOrbMomentum !returns the total orbital momentum and its Z-axis projection
          procedure, public:: compare=>SpherSymmetryCompare                 !compares two spherical symmetries
          procedure, public:: combine=>SpherSymmetryCombine                 !combines two spherical symmetries (common lower irrep, if any)
          procedure, public:: print_it=>SpherSymmetryPrintIt                !prints the symmetry information
          final:: spher_symmetry_dtor                                       !destructor
        end type spher_symmetry_t
 !Abstract basis function (basis function support only):
        type, public:: basis_func_supp_t
         integer(INTD), private:: supp_dim=-1     !dimensionality of the real space support on which the basis function resides
         type(real_vec_t), private:: center       !center of the effective function support in the real space
         type(orthotope_t), private:: supp_box    !supporting orthotope (multidimensional real space support)
         contains
          procedure, private:: BasisFuncSuppCtorEmpty               !constructs a trivial basis function support (ctor)
          procedure, private:: BasisFuncSuppCtorReal                !constructs a non-trivial basis function support (ctor)
          generic, public:: basis_func_supp_ctor=>BasisFuncSuppCtorEmpty,BasisFuncSuppCtorReal
          procedure, public:: is_set=>BasisFuncSuppIsSet            !returns .TRUE. if the basis function support is set
          procedure, public:: supp_dimsn=>BasisFuncSuppDimsn        !returns the support dimension (>=0), 0 is trivial (no real support)
#if !(defined(__GNUC__) && __GNUC__ < 8)
          final:: basis_func_supp_dtor
#endif
        end type basis_func_supp_t
 !Gaussian basis function:
        type, extends(basis_func_supp_t), public:: basis_func_gauss_t
         integer(INTD), private:: num_prims=0        !number of contracted primitives
         integer(INTD), private:: orb_moment=-1      !orbital momentum (0,1,2,3,...)
         real(8), allocatable, private:: exponent(:) !primitive exponents
         complex(8), allocatable, private:: coef(:)  !primitive contraction coefficients
         contains
          procedure, private:: BasisFuncGaussCtor                     !sets up the basis function (ctor)
          generic, public:: basis_func_gauss_ctor=>BasisFuncGaussCtor
          final:: basis_func_gauss_dtor                               !destructs the basis function (dtor)
        end type basis_func_gauss_t
 !Typed basis function:
        type, public:: basis_func_t
         integer(INTD), private:: basis_kind=BASIS_NONE                    !specific basis kind (mandatory)
         class(basis_func_supp_t), pointer, private:: basis_func_p=>NULL() !non-owning pointer to a persistent basis function of this basis kind (optional)
         class(symmetry_t), pointer, private:: symm_p=>NULL()              !non-owning pointer to a persistent symmetry object of the basis function (optional)
         contains
          procedure, private:: BasisFuncCtor                         !sets up the basis function (ctor)
          generic, public:: basis_func_ctor=>BasisFuncCtor
          procedure, private:: get_basis_func=>BasisFuncGetBasisFunc !returns basis function description (attributes)
          final:: basis_func_dtor                                    !destructs the basis function (dtor)
        end type basis_func_t
 !Subspace basis:
        type, public:: subspace_basis_t
         integer(INTL), private:: space_dim=0                     !number of basis functions
         integer(INTD), private:: supp_dim=0                      !dimensionality of the real space support on which the basis functions reside
         type(real_vec_t), private:: center                       !center of the effective subspace basis support in real space
         type(orthotope_t), private:: supp_box                    !effective subspace basis support in real space (multidimensional orthotope)
         class(symmetry_t), allocatable, private:: symm           !symmetry of the subspace basis (if any, for all basis functions)
         type(basis_func_t), allocatable, private:: basis_func(:) !basis functions specified by reference: [1..space_dim]
         contains
          procedure, private:: SubspaceBasisCtor                        !creates an abstract subspace without specific basis functions (ctor)
          generic, public:: subspace_basis_ctor=>SubspaceBasisCtor
          procedure, public:: dimsn=>SubspaceBasisDimsn                 !returns the dimension of the subspace
          procedure, public:: supp_dimsn=>SubspaceBasisSuppDimsn        !returns the support space dimension
          procedure, public:: set_basis_func=>SubspaceBasisSetBasisFunc !sets a specific basis function (builder)
          procedure, public:: get_basis_func=>SubspaceBasisGetBasisFunc !returns a pointer to a specific basis function from the subspace basis
          procedure, public:: finalize=>SubspaceBasisFinalize           !finalizes the subspace basis (sets up the support and overall symmetry)
          procedure, public:: get_symmetry=>SubspaceBasisGetSymmetry    !returns a pointer to the subspace basis symmetry object
          procedure, public:: get_center=>SubspaceBasisGetCenter        !returns a pointer to the center of the subspace basis in the real space
          procedure, public:: get_support=>SubspaceBasisGetSupport      !returns a pointer to the supporting orthotope of the subspace basis
          final:: subspace_basis_dtor                                   !destructs the subspace basis (dtor)
        end type subspace_basis_t
 !Subspace:
        type, public:: subspace_t
         integer(INTL), private:: subspace_id=-1   !unique subspace ID (registered ID): must be non-negative, -1 means undefined
         type(seg_int_t), private:: basis_subrange !subrange of the basis vectors defining this subspace, specified as a semi-interval
         integer(INTD), private:: supp_dim=0       !dimensionality of the real space support on which the basis functions reside
         integer(INTL), private:: max_resolution=0 !max resolution level (max subspace dimension): 0 means undefined
         type(list_bi_t), private:: bases          !basis sets (subspace_basis_t) for each registered resolution level
         contains
          procedure, private:: SubspaceCtorBase                            !initializes a subspace with id and real space support dimension only (ctor)
          generic, public:: subspace_ctor=>SubspaceCtorBase
          procedure, public:: is_set=>SubspaceIsSet                        !returns TRUE if the subspace is set (constructed)
          procedure, public:: get_id=>SubspaceGetId                        !returns the subspace id
          procedure, public:: get_basis_subrange=>SubspaceGetBasisSubrange !returns the subspace basis subrange
          procedure, public:: get_supp_dim=>SubspaceGetSuppDim             !returns the dimensionality of the subspace support
          procedure, public:: get_max_resolution=>SubspaceGetMaxResolution !returns the max resolution (dimension) of the subspace
          procedure, public:: compare_range=>SubspaceCompareRange          !compares ranges of two subspaces: {CMP_EQ,CMP_LT,CMP_GT,CMP_OV,CMP_ER}
          procedure, public:: register_basis=>SubspaceRegisterBasis        !registers a specific basis of the subspace
          procedure, public:: resolve=>SubspaceResolve                     !resolves the subspace with a specific basis (based on some condition)
          procedure, public:: print_it=>SubspacePrintIt                    !prints the subspace information
          final:: subspace_dtor                                            !destroys the subspace (dtor)
        end type subspace_t
 !Hierarchical composite index:
        type, public:: h_index_t
         integer(INTL), public:: subspace_id=-1 !subspace ID (registered ID): must be non-negative, -1 means undefined
         integer(INTL), public:: resolution=0   !subspace resolution level 1<=resolution<=max_resolution: 0 means undefined
         integer(INTL), public:: component=0    !subspace component number at the given level of resolution: [1..resolution], 0 means undefined
        end type h_index_t
 !Hierarchical vector space representation:
        type, public:: h_space_t
         integer(INTL), private:: space_dim=0                 !dimension of the vector space
         integer(INTL), private:: num_subspaces=0             !number of subspaces defined in the vector space: [0..space_dim-1] are original basis functions, [space_dim..num_subspaces-1] are their aggregates
         type(vec_tree_t), private:: subspaces                !subspaces defined in the vector space: [0..space_dim-1] are original basis functions, [space_dim..num_subspaces-1] are their aggregates
         complex(8), pointer, private:: metric_p(:,:)=>NULL() !non-owning pointer to the original (persistent) metric tensor: g12=<bf1|bf2>
         real(8), allocatable, private:: overlap(:,:)         !subspace support overlap matrix (extent of support overlap between all subspaces)
         contains
          procedure, private:: HSpaceCtorSimple                         !constructs a simple hierarchical representation of a vector space (ctor)
          generic, public:: h_space_ctor=>HSpaceCtorSimple              !ctors
          procedure, public:: is_set=>HSpaceIsSet                       !returns TRUE if the hierarchical vector space is set
          procedure, public:: get_space_dim=>HSpaceGetSpaceDim          !returns the dimension of the vector space
          procedure, public:: get_num_subspaces=>HSpaceGetNumSubspaces  !returns the total number of defined subspaces in the vector space
          procedure, public:: get_root_id=>HSpaceGetRootId              !returns the id of the full space (root of the subspace aggregation tree)
          procedure, public:: get_ancestor_id=>HSpaceGetAncestorId      !returns the id of a specific ancestor subspace
          procedure, public:: get_subspace_level=>HSpaceGetSubspaceLevel!returns the distance from the root for the specific subspace
          procedure, public:: get_subspace=>HSpaceGetSubspace           !returns a pointer to the requested subspace of the hierarchical vector space
          procedure, public:: get_aggr_tree=>HSpaceGetAggrTree          !returns a pointer to the subspace aggregation tree (->subspaces)
          procedure, public:: get_level_composition=>HSpaceGetLevelComposition !returns an ordered list of subspace ids forming a given level of the subspace aggregation tree (direct sum decomposition)
          procedure, public:: compare_subspaces=>HSpaceCompareSubspaces !compares two subspaces from the hierarchical vector space: {CMP_EQ,CMP_LT,CMP_GT,CMP_ER}
          procedure, public:: relate_subspaces=>HSpaceRelateSubspaces   !relates two subspaces from the hierarchical vector space: {CMP_EQ,CMP_CN,CMP_IN,CMP_OV,CMP_NC}
          procedure, public:: compare_subranges=>HSpaceCompareSubranges !compares basis subranges of two subspaces from the hierarchical vector space: {CMP_EQ,CMP_LT,CMP_GT,CMP_OV,CMP_ER}
          procedure, public:: get_common_subspace=>HSpaceGetCommonSubspace !returns the minimal (registered) subspace containing two given subspaces
          procedure, public:: print_it=>HSpacePrintIt                   !prints the hierarchical vector space
          final:: h_space_dtor                                          !destructs the hierarchical representation of a vector space
        end type h_space_t
!INTERFACES:
 !symmetry_t:
        abstract interface
  !Deferred: .compare:
         function symm_compare_i(this,symm) result(cmp)
          import:: symmetry_t,INTD
          implicit none
          integer(INTD):: cmp                  !out: comparison result (see module dil_basic)
          class(symmetry_t), intent(in):: this !in: symmetry 1
          class(symmetry_t), intent(in):: symm !in: symmetry 2
         end function symm_compare_i
  !Deferred: .combine:
         subroutine symm_combine_i(this,symm,ierr)
          import:: symmetry_t,INTD
          implicit none
          class(symmetry_t), intent(inout):: this     !inout: symmetry 1 (updated)
          class(symmetry_t), intent(in):: symm        !in: symmetry 2
          integer(INTD), intent(out), optional:: ierr !out: error code
         end subroutine symm_combine_i
        end interface
!VISIBILITY:
 !Non-member:
        public build_basis_hierarchy_abstract   !establishes a hierarchy for an abstract basis with possible symmetries
        public build_basis_hierarchy_real_space !establishes a hierarchy for a real space supported basis with possible symmetries
 !real_vec_t:
        private RealVecCtor
        private RealVecDimsn
        private RealVecNorm2
        private RealVecScale
        private RealVecAdd
        private RealVecAverage
        public real_vec_dtor
 !range1d_t:
        private Range1dSet
        private Range1dLowerBound
        private Range1dUpperBound
        private Range1dLength
        private Range1dOverlap
        private Range1dUnion
        private Range1dSplit
 !seg_int_t:
        private SegIntSet
        private SegIntLowerBound
        private SegIntUpperBound
        private SegIntLength
        private SegIntOverlap
        private SegIntUnion
        private SegIntSplit
        private SegIntPrintRange
        public seg_int_print_range !debug
 !orthotope_t:
        private OrthotopeCtor
        private OrthotopeDimsn
        private OrthotopeSetExtent
        private OrthotopeGetExtent
        private OrthotopeLowerBound
        private OrthotopeUpperBound
        private OrthotopeLength
        private OrthotopeVolume
        private OrthotopeOverlap
        private OrthotopeUnion
        public orthotope_dtor
 !color_symmetry_t:
        private ColorSymmetryCtor
        private ColorSymmetryCompare
        private ColorSymmetryCombine
        public color_symmetry_dtor
 !spher_symmetry_t:
        private SpherSymmetryCtor
        private SpherSymmetryGetOrbMomentum
        private SpherSymmetryCompare
        private SpherSymmetryCombine
        private SpherSymmetryPrintIt
        public spher_symmetry_dtor
 !basis_func_supp_t:
        private BasisFuncSuppCtorEmpty
        private BasisFuncSuppCtorReal
        private BasisFuncSuppIsSet
        private BasisFuncSuppDimsn
        public basis_func_supp_dtor
 !basis_func_gauss_t:
        private BasisFuncGaussCtor
        public basis_func_gauss_dtor
 !basis_func_t:
        private BasisFuncCtor
        private BasisFuncGetBasisFunc
        public basis_func_dtor
 !subspace_basis_t:
        private SubspaceBasisCtor
        private SubspaceBasisDimsn
        private SubspaceBasisSuppDimsn
        private SubspaceBasisSetBasisFunc
        private SubspaceBasisGetBasisFunc
        private SubspaceBasisFinalize
        private SubspaceBasisGetSymmetry
        private SubspaceBasisGetCenter
        private SubspaceBasisGetSupport
        public subspace_basis_dtor
 !subspace_t:
        private SubspaceCtorBase
        private SubspaceIsSet
        private SubspaceGetId
        private SubspaceGetBasisSubrange
        private SubspaceGetSuppDim
        private SubspaceGetMaxResolution
        private SubspaceCompareRange
        private SubspaceRegisterBasis
        private SubspaceResolve
        private SubspacePrintIt
        public subspace_dtor
 !h_space_t:
        private HSpaceCtorSimple
        private HSpaceIsSet
        private HSpaceGetSpaceDim
        private HSpaceGetNumSubspaces
        private HSpaceGetRootId
        private HSpaceGetAncestorId
        private HSpaceGetSubspaceLevel
        private HSpaceGetSubspace
        private HSpaceGetAggrTree
        private HSpaceGetLevelComposition
        private HSpaceCompareSubspaces
        private HSpaceRelateSubspaces
        private HSpaceCompareSubranges
        private HSpaceGetCommonSubspace
        private HSpacePrintIt
        public h_space_dtor

       contains
!IMPLEMENTATION:
![real_vec_t]==================================
        subroutine RealVecCtor(this,dimsn,ierr)
!Creates an empty real space vector. If the vector is defined on input,
!it will be automatically destructed prior to the re-initialization.
         implicit none
         class(real_vec_t), intent(out):: this       !out: empty real space vector
         integer(INTL), intent(in):: dimsn           !in: desired vector dimension
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(.not.allocated(this%coord)) then
          if(dimsn.gt.0) then
           allocate(this%coord(1:dimsn),STAT=errc)
           if(errc.eq.0) then; this%num_dim=dimsn; else; errc=3; endif
          else
           errc=2
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine RealVecCtor
!------------------------------------------------
        function RealVecDimsn(this) result(dimsn)
!Returns the dimension of the vector.
         implicit none
         integer(INTL):: dimsn                       !out: vector dimension
         class(real_vec_t), intent(in):: this        !in: real space vector

         dimsn=this%num_dim
         return
        end function RealVecDimsn
!-----------------------------------------------------
        function RealVecNorm2(this,ierr) result(norm2)
!Returns the 2-norm of the vector.
         implicit none
         real(8):: norm2                             !out: vector 2-norm
         class(real_vec_t), intent(in):: this        !in: real space vector
         integer(INTD), intent(out), optional:: ierr !out: error code
         !------------------------------------------
         integer(INTL), parameter:: LARGE_VECTOR=(2_INTL)**20
         !---------------------------------------------------
         integer(INTD):: errc
         integer(INTL):: i

         errc=0; norm2=0d0
         if(this%num_dim.ge.LARGE_VECTOR) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED) REDUCTION(+:norm2)
          do i=1,this%num_dim
           norm2=norm2+(this%coord(i))*(this%coord(i))
          enddo
!$OMP END PARALLEL DO
         else
          if(this%num_dim.gt.0) then
           do i=1,this%num_dim
            norm2=norm2+(this%coord(i))*(this%coord(i))
           enddo
          else
           errc=1
          endif
         endif
         if(errc.eq.0) norm2=dsqrt(norm2)
         if(present(ierr)) ierr=errc
         return
        end function RealVecNorm2
!------------------------------------------------
        subroutine RealVecScale(this,scalar,ierr)
!Multiplies a real space vector by a scalar.
         implicit none
         class(real_vec_t), intent(inout):: this     !inout: real vector
         real(8), intent(in):: scalar                !in: scalar factor
         integer(INTD), intent(out), optional:: ierr !out: error code
         !---------------------------------------------------
         integer(INTL), parameter:: LARGE_VECTOR=(2_INTL)**20
         !---------------------------------------------------
         integer(INTD):: errc
         integer(INTL):: i,n

         errc=0; n=this%num_dim
         if(n.ge.LARGE_VECTOR) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
          do i=1,n; this%coord(i)=this%coord(i)*scalar; enddo
!$OMP END PARALLEL DO
         elseif(n.le.0) then
          errc=1
         else
          do i=1,n; this%coord(i)=this%coord(i)*scalar; enddo
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine RealVecScale
!--------------------------------------------------------
        subroutine RealVecAdd(this,real_vec,ierr,res_vec)
!Adds two vectors. If <res_vec> is present, it will contain the result,
!otherwise <this> will contain the result.
         implicit none
         class(real_vec_t), intent(inout):: this             !inout: real space vector 1 (can be updated)
         class(real_vec_t), intent(in):: real_vec            !in: real space vector 2
         integer(INTD), intent(out), optional:: ierr         !out: error code
         type(real_vec_t), intent(inout), optional:: res_vec !out: resulting vector
         !---------------------------------------------------
         integer(INTL), parameter:: LARGE_VECTOR=(2_INTL)**20
         !---------------------------------------------------
         integer(INTD):: errc
         integer(INTL):: i,n

         errc=0; n=this%num_dim
         if(n.gt.0.and.n.eq.real_vec%num_dim) then
          if(present(res_vec)) then
           if(res_vec%num_dim.ne.n) call res_vec%real_vec_ctor(n,errc)
           if(errc.eq.0) then
            if(n.ge.LARGE_VECTOR) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
             do i=1,n; res_vec%coord(i)=this%coord(i)+real_vec%coord(i); enddo
!$OMP END PARALLEL DO
            else
             do i=1,n; res_vec%coord(i)=this%coord(i)+real_vec%coord(i); enddo
            endif
           else
            errc=2
           endif
          else
           if(n.ge.LARGE_VECTOR) then
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) SCHEDULE(GUIDED)
            do i=1,n; this%coord(i)=this%coord(i)+real_vec%coord(i); enddo
!$OMP END PARALLEL DO
           else
            do i=1,n; this%coord(i)=this%coord(i)+real_vec%coord(i); enddo
           endif
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine RealVecAdd
!------------------------------------------------------------
        subroutine RealVecAverage(this,real_vec,ierr,res_vec)
!Computes the average of two vectors. If <res_vec> is present,
!it will contain the result, otherwise <this> will contain the result.
         implicit none
         class(real_vec_t), intent(inout):: this             !inout: real space vector 1 (can be updated)
         class(real_vec_t), intent(in):: real_vec            !in: real space vector 2
         integer(INTD), intent(out), optional:: ierr         !out: error code
         type(real_vec_t), intent(inout), optional:: res_vec !out: resulting vector
         integer(INTD):: errc

         errc=0
         if(present(res_vec)) then
          call this%add(real_vec,errc,res_vec)
          if(errc.eq.0) call res_vec%scale(0.5d0,errc)
         else
          call this%add(real_vec,errc)
          if(errc.eq.0) call this%scale(0.5d0,errc)
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine RealVecAverage
!-------------------------------------
        subroutine real_vec_dtor(this)
         implicit none
         type(real_vec_t):: this

         if(allocated(this%coord)) deallocate(this%coord)
         this%num_dim=0
         return
        end subroutine real_vec_dtor
![range1d_t]========================================
        subroutine Range1dSet(this,lower,upper,ierr)
!Range ctor (defines or redefines the range).
         implicit none
         class(range1d_t), intent(inout):: this      !inout: range
         real(8), intent(in), optional:: lower       !in: new lower bound
         real(8), intent(in), optional:: upper       !in: new upper bound
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         real(8):: l,u

         errc=0
         if(present(lower)) then; l=lower; else; l=this%min_coord; endif
         if(present(upper)) then; u=upper; else; u=this%max_coord; endif
         if(l.le.u) then
          this%min_coord=l; this%max_coord=u
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine Range1dSet
!--------------------------------------------------
        function Range1dLowerBound(this) result(lb)
!Returns the range lower bound.
         implicit none
         real(8):: lb                        !out: lower bound
         class(range1d_t), intent(in):: this !in: range

         lb=this%min_coord
         return
        end function Range1dLowerBound
!--------------------------------------------------
        function Range1dUpperBound(this) result(ub)
!Returns the range upper bound.
         implicit none
         real(8):: ub                        !out: upper bound
         class(range1d_t), intent(in):: this !in: range

         ub=this%max_coord
         return
        end function Range1dUpperBound
!--------------------------------------------------
        function Range1dLength(this) result(length)
!Returns the length of the range.
         implicit none
         real(8):: length                    !out: length of the range
         class(range1d_t), intent(in):: this !in: range

         length=this%max_coord-this%min_coord
         return
        end function Range1dLength
!------------------------------------------------------
        subroutine Range1dOverlap(this,range,res_range)
!Returns the overlap of two ranges. If <res_range> is present,
!it will contain the overlap. Otherwise, <this> will contain the result.
         implicit none
         class(range1d_t), intent(inout):: this                !inout: range 1 (can be updated)
         class(range1d_t), intent(in):: range                  !in: range 2
         type(range1d_t), intent(inout), optional:: res_range  !out: resulting range

         if(this%max_coord.le.range%min_coord.or.this%min_coord.ge.range%max_coord) then !no overlap
          if(present(res_range)) then
           call res_range%set(0d0,0d0)
          else
           call this%set(0d0,0d0)
          endif
         else !overlap
          if(present(res_range)) then
           call res_range%set(max(this%min_coord,range%min_coord),min(this%max_coord,range%max_coord))
          else
           call this%set(max(this%min_coord,range%min_coord),min(this%max_coord,range%max_coord))
          endif
         endif
         return
        end subroutine Range1dOverlap
!----------------------------------------------------
        subroutine Range1dUnion(this,range,res_range)
!Returns the union of two ranges. If <res_range> is present,
!it will contain the union. Otherwise, <this> will contain the result.
         implicit none
         class(range1d_t), intent(inout):: this               !inout: range 1 (can be updated)
         class(range1d_t), intent(in):: range                 !in: range 2
         type(range1d_t), intent(inout), optional:: res_range !out: resulting range

         if(present(res_range)) then
          call res_range%set(min(this%min_coord,range%min_coord),max(this%max_coord,range%max_coord))
         else
          call this%set(min(this%min_coord,range%min_coord),max(this%max_coord,range%max_coord))
         endif
         return
        end subroutine Range1dUnion
!-------------------------------------------------------------
        subroutine Range1dSplit(this,num_segs,segs,ierr,align)
!Splits the range in a number of contiguous segments with an optional alignment.
!If present, align(:) must contain strictly positive real offsets within the range
!ordered in an ascending order.
         implicit none
         class(range1d_t), intent(in):: this          !in: input range
         integer(INTD), intent(in):: num_segs         !in: number of segments to split the range into
         class(range1d_t), intent(inout):: segs(1:)   !out: segments (subranges)
         integer(INTD), intent(out), optional:: ierr  !out: error code
         real(8), intent(inout), optional:: align(1:) !in: inner alignment boundaries (>0) relative to the beginning of the range, excluding beginning and end of the range
         integer(INTD):: errc,nbnd,nchnk,i,j,k,l,left,left2,right
         real(8):: rl,incr,lb,ub,rb
         logical:: next

         errc=0
         if(num_segs.gt.0.and.num_segs.le.size(segs)) then
          rl=this%length()
          if(rl.gt.0d0) then
           if(present(align)) then; nbnd=size(align); else; nbnd=0; endif
           if(nbnd.gt.0) then !aligned splitting
            if(num_segs.le.nbnd) then
             ub=this%upper_bound()
             nchnk=nbnd+1 !number of alignment chunks
             do while(nchnk.gt.num_segs)
 !Find the smallest alignment chunk:
              left=0; left2=0; right=0
              j=0; k=0; l=0; lb=0d0; next=.TRUE.
              do i=1,nbnd
               if(align(i).gt.0d0) then !boundary is still active
                if(.not.next) then; right=i; next=.TRUE.; endif
                if(j.eq.0) then
                 incr=align(i)-lb; j=i; left=k; left2=l; next=.FALSE.
                else
                 if(align(i)-lb.lt.incr) then; incr=align(i)-lb; j=i; left=k; left2=l; next=.FALSE.; endif
                endif
                lb=align(i); l=k; k=i
               endif
              enddo
              if(.not.next) then; right=nbnd+1; next=.TRUE.; endif
              if(ub-lb.lt.incr) then; incr=ub-lb; j=nbnd+1; left=k; left2=l; right=0; endif
 !Merge the smallest alignment chunk with its left or right neighbor:
              if(left2.gt.0.and.right.gt.0) then
               if(right.le.nbnd) then; rb=align(right); else; rb=ub; endif
               if(align(left)-align(left2).lt.rb-align(j)) then
                align(left)=-align(left) !deactivate the boundary due to merge
               else
                align(j)=-align(j) !deactivate the boundary due to merge
               endif
              else
               if(right.gt.0) then
                align(j)=-align(j) !deactivate the boundary due to merge
               else
                align(left)=-align(left) !deactivate the boundary due to merge
               endif
              endif
              nchnk=nchnk-1
             enddo
 !Set the final boundaries:
             lb=this%lower_bound(); ub=lb; j=0
             do i=1,nbnd
              if(align(i).gt.0d0) then
               j=j+1; call segs(j)%set(ub,lb+align(i),errc); if(errc.ne.0) exit
               ub=lb+align(i)
              else
               align(i)=-align(i)
              endif
             enddo
             if(errc.eq.0) call segs(num_segs)%set(ub,this%upper_bound(),errc)
            elseif(num_segs.eq.nbnd+1) then
             lb=this%lower_bound(); call segs(1)%set(lb,lb+align(1),errc)
             if(errc.eq.0) then
              do i=2,num_segs-1
               call segs(i)%set(lb+align(i-1),lb+align(i),errc); if(errc.ne.0) exit
              enddo
              if(errc.eq.0) call segs(num_segs)%set(lb+align(nbnd),this%upper_bound(),errc)
             endif
            else
             errc=3
            endif
           else !unaligned splitting
            lb=this%lower_bound(); incr=rl/real(num_segs,8)
            do i=1,num_segs-1
             ub=lb+incr
             call segs(i)%set(lb,ub,errc); if(errc.ne.0) exit
             lb=ub
            enddo
            if(errc.eq.0) call segs(num_segs)%set(lb,this%upper_bound(),errc)
           endif
          else
           errc=2
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine Range1dSplit
![seg_int_t]=======================================
        subroutine SegIntSet(this,lower,upper,ierr)
!seg_int_t ctor (defines or redefines the integer range).
         implicit none
         class(seg_int_t), intent(inout):: this      !inout: integer range
         integer(INTL), intent(in), optional:: lower !in: new lower bound
         integer(INTL), intent(in), optional:: upper !in: new upper bound
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer(INTL):: l,u

         errc=0
         if(present(lower)) then; l=lower; else; l=this%min_coord; endif
         if(present(upper)) then; u=upper; else; u=this%max_coord; endif
         if(l.le.u) then
          this%min_coord=l; this%max_coord=u
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine SegIntSet
!-------------------------------------------------
        function SegIntLowerBound(this) result(lb)
!Returns the lower bound of the integer range.
         implicit none
         integer(INTL):: lb                  !out: lower bound
         class(seg_int_t), intent(in):: this !in: integer range

         lb=this%min_coord
         return
        end function SegIntLowerBound
!-------------------------------------------------
        function SegIntUpperBound(this) result(ub)
!Returns the upper bound of the integer range.
         implicit none
         integer(INTL):: ub                  !out: upper bound
         class(seg_int_t), intent(in):: this !in: integer range

         ub=this%max_coord
         return
        end function SegIntUpperBound
!-------------------------------------------------
        function SegIntLength(this) result(length)
!Returns the length of the integer range.
         implicit none
         integer(INTL):: length              !out: length of the integer range
         class(seg_int_t), intent(in):: this !in: integer range

         length=this%max_coord-this%min_coord
         return
        end function SegIntLength
!-----------------------------------------------------
        subroutine SegIntOverlap(this,range,res_range)
!Returns the overlap of two integer ranges. If <res_range> is present,
!it will contain the overlap. Otherwise, <this> will contain the result.
         implicit none
         class(seg_int_t), intent(inout):: this                !inout: integer range 1 (can be updated)
         class(seg_int_t), intent(in):: range                  !in: integer range 2
         type(seg_int_t), intent(inout), optional:: res_range  !out: resulting integer range

         if(this%max_coord.le.range%min_coord.or.this%min_coord.ge.range%max_coord) then !no overlap
          if(present(res_range)) then
           call res_range%set(0_INTL,0_INTL)
          else
           call this%set(0_INTL,0_INTL)
          endif
         else !overlap
          if(present(res_range)) then
           call res_range%set(max(this%min_coord,range%min_coord),min(this%max_coord,range%max_coord))
          else
           call this%set(max(this%min_coord,range%min_coord),min(this%max_coord,range%max_coord))
          endif
         endif
         return
        end subroutine SegIntOverlap
!---------------------------------------------------
        subroutine SegIntUnion(this,range,res_range)
!Returns the union of two integer ranges. If <res_range> is present,
!it will contain the union. Otherwise, <this> will contain the result.
         implicit none
         class(seg_int_t), intent(inout):: this               !inout: integer range 1 (can be updated)
         class(seg_int_t), intent(in):: range                 !in: integer range 2
         type(seg_int_t), intent(inout), optional:: res_range !out: resulting integer range

         if(present(res_range)) then
          call res_range%set(min(this%min_coord,range%min_coord),max(this%max_coord,range%max_coord))
         else
          call this%set(min(this%min_coord,range%min_coord),max(this%max_coord,range%max_coord))
         endif
         return
        end subroutine SegIntUnion
!------------------------------------------------------------
        subroutine SegIntSplit(this,num_segs,segs,ierr,align)
!Splits the integer range in a number of integer segments with an optional alignment.
!If present, align(:) must contain strictly positive integer offsets within the range
!ordered in an ascending order.
         implicit none
         class(seg_int_t), intent(in):: this                !in: input integer range
         integer(INTD), intent(in):: num_segs               !in: number of segments to split the integer range into
         class(seg_int_t), intent(inout):: segs(1:)         !out: integer segments (subranges)
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTL), intent(inout), optional:: align(1:) !in: inner alignment boundaries (>0) relative to the beginning of the range, excluding beginning and end of the range
         integer(INTD):: errc,nbnd,nchnk,i,j,k,l,left,left2,right
         integer(INTL):: rl,incr,rem,lb,ub
         logical:: next

         errc=0; rl=this%length()
         if(num_segs.gt.0.and.num_segs.le.size(segs).and.num_segs.le.rl) then
          if(rl.gt.0) then
           if(present(align)) then; nbnd=size(align); else; nbnd=0; endif
           !write(CONS_OUT,'("#DEBUG(seg_int_t:split): align(:) array size = ",i11)') nbnd !debug
           if(nbnd.gt.0) then !aligned splitting
            if(num_segs.le.nbnd) then
             ub=this%upper_bound()
             nchnk=nbnd+1 !number of alignment chunks
             do while(nchnk.gt.num_segs)
 !Find the smallest alignment chunk:
              left=0; left2=0; right=0
              j=0; k=0; l=0; lb=0; next=.TRUE.
              do i=1,nbnd
               if(align(i).gt.0) then !boundary is still active
                if(.not.next) then; right=i; next=.TRUE.; endif
                if(j.eq.0) then
                 incr=align(i)-lb; j=i; left=k; left2=l; next=.FALSE.
                else
                 if(align(i)-lb.lt.incr) then; incr=align(i)-lb; j=i; left=k; left2=l; next=.FALSE.; endif
                endif
                lb=align(i); l=k; k=i
               endif
              enddo
              if(.not.next) then; right=nbnd+1; next=.TRUE.; endif
              if(ub-lb.lt.incr) then; incr=ub-lb; j=nbnd+1; left=k; left2=l; right=0; endif
 !Merge the smallest alignment chunk with its left or right neighbor:
              if(left2.gt.0.and.right.gt.0) then
               if(right.le.nbnd) then; i=align(right); else; i=ub; endif
               if(align(left)-align(left2).lt.i-align(j)) then
                align(left)=-align(left) !deactivate the boundary due to merge
               else
                align(j)=-align(j) !deactivate the boundary due to merge
               endif
              else
               if(right.gt.0) then
                align(j)=-align(j) !deactivate the boundary due to merge
               else
                align(left)=-align(left) !deactivate the boundary due to merge
               endif
              endif
              nchnk=nchnk-1
             enddo
 !Set the final boundaries:
             lb=this%lower_bound(); ub=lb; j=0
             do i=1,nbnd
              if(align(i).gt.0) then
               j=j+1; call segs(j)%set(ub,lb+align(i),errc); if(errc.ne.0) exit
               ub=lb+align(i)
              else
               align(i)=-align(i)
              endif
             enddo
             if(errc.eq.0) call segs(num_segs)%set(ub,this%upper_bound(),errc)
            elseif(num_segs.eq.nbnd+1) then
             lb=this%lower_bound(); call segs(1)%set(lb,lb+align(1),errc)
             if(errc.eq.0) then
              do i=2,num_segs-1
               call segs(i)%set(lb+align(i-1),lb+align(i),errc); if(errc.ne.0) exit
              enddo
              if(errc.eq.0) call segs(num_segs)%set(lb+align(nbnd),this%upper_bound(),errc)
             endif
            else
             errc=3
            endif
           else !unaligned splitting
            lb=this%lower_bound(); incr=rl/int(num_segs,INTL); rem=rl-incr*int(num_segs,INTL)
            do i=1,num_segs
             ub=lb+incr; if(i.le.rem) ub=ub+1_INTL
             call segs(i)%set(lb,ub,errc); if(errc.ne.0) exit
             lb=ub
            enddo
           endif
          else
           errc=2
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine SegIntSplit
!------------------------------------------------
        subroutine SegIntPrintRange(this,dev_out)
!Prints the integer range.
         implicit none
         class(seg_int_t), intent(in):: this           !in: integer range
         integer(INTD), intent(in), optional:: dev_out !in: output device (defaults to screen)
         integer(INTD):: devo

         if(present(dev_out)) then; devo=dev_out; else; devo=6; endif
         write(devo,*) 'Range [',this%lower_bound()+1_INTL,this%upper_bound(),']'
         return
        end subroutine SegIntPrintRange
!-----------------------------------------------------
        function seg_int_print_range(obj) result(ierr) !debug
!Non-member generic printing action for use in GFC.
         implicit none
         integer(INTD):: ierr                  !out: error code
         class(*), intent(inout), target:: obj !in: seg_int_t object
         class(seg_int_t), pointer:: obp

         ierr=0; obp=>NULL()
         select type(obj); class is(seg_int_t); obp=>obj; end select
         if(associated(obp)) then
          call obp%print_range()
         else
          ierr=-1
         endif
         return
        end function seg_int_print_range
![orthotope_t]===================================
        subroutine OrthotopeCtor(this,dimsn,ierr)
!Creates an empty orthotope. If the orthotope is defined on input,
!it will automatically be destructed prior to redefinition.
         implicit none
         class(orthotope_t), intent(out):: this      !out: empty orthotope
         integer(INTL), intent(in):: dimsn           !in: real space dimension
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(.not.allocated(this%extent)) then
          if(dimsn.gt.0) then
           allocate(this%extent(1:dimsn),STAT=errc)
           if(errc.eq.0) then; this%num_dim=dimsn; else; errc=3; endif
          else
           errc=2
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine OrthotopeCtor
!--------------------------------------------------
        function OrthotopeDimsn(this) result(dimsn)
!Returns the real space dimension the orthotope lives in.
         implicit none
         integer(INTL):: dimsn                       !out: real space dimension (0 means empty orthotope)
         class(orthotope_t), intent(in):: this       !in: orthotope

         dimsn=this%num_dim
         return
        end function OrthotopeDimsn
!------------------------------------------------------------
        subroutine OrthotopeSetExtent(this,dimsn,extent,ierr)
!Sets the orthotope extent along a specific dimension.
         implicit none
         class(orthotope_t), intent(inout):: this    !inout: orthotope
         integer(INTL), intent(in):: dimsn           !in: specific dimension to set extent over
         type(range1d_t), intent(in):: extent        !in: extent
         integer(INTD), intent(out), optional:: ierr !error code
         integer(INTD):: errc

         errc=0
         if(dimsn.gt.0.and.dimsn.le.this%num_dim) then
          this%extent(dimsn)=extent
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine OrthotopeSetExtent
!------------------------------------------------------------------
        function OrthotopeGetExtent(this,dimsn,ierr) result(extent)
!Returns the orthotope extent along a specific dimension.
         implicit none
         type(range1d_t):: extent                    !out: extent of the specific dimension
         class(orthotope_t), intent(in):: this       !in: orthotope
         integer(INTL), intent(in):: dimsn           !in: specific dimension
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(dimsn.gt.0.and.dimsn.le.this%num_dim) then
          extent=this%extent(dimsn)
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function OrthotopeGetExtent
!---------------------------------------------------------------
        function OrthotopeLowerBound(this,dimsn,ierr) result(lb)
!Returns the orthotope extent lower bound over a specific dimension.
         implicit none
         real(8):: lb                                !out: extent lower bound
         class(orthotope_t), intent(in):: this       !in: orthotope
         integer(INTL), intent(in):: dimsn           !in: specific dimension
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; lb=0d0
         if(dimsn.gt.0.and.dimsn.le.this%num_dim) then
          lb=this%extent(dimsn)%lower_bound()
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function OrthotopeLowerBound
!---------------------------------------------------------------
        function OrthotopeUpperBound(this,dimsn,ierr) result(ub)
!Returns the orthotope extent upper bound over a specific dimension.
         implicit none
         real(8):: ub                                !out: extent upper bound
         class(orthotope_t), intent(in):: this       !in: orthotope
         integer(INTL), intent(in):: dimsn           !in: specific dimension
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; ub=0d0
         if(dimsn.gt.0.and.dimsn.le.this%num_dim) then
          ub=this%extent(dimsn)%upper_bound()
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function OrthotopeUpperBound
!---------------------------------------------------------------
        function OrthotopeLength(this,dimsn,ierr) result(length)
!Returns the orthotope length over a specific dimension.
         implicit none
         real(8):: length                            !out: length over a specific dimension
         class(orthotope_t), intent(in):: this       !in: orthotope
         integer(INTL), intent(in):: dimsn           !in: specific dimension
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; length=0d0
         if(dimsn.gt.0.and.dimsn.le.this%num_dim) then
          length=this%extent(dimsn)%length()
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function OrthotopeLength
!---------------------------------------------------------
        function OrthotopeVolume(this,ierr) result(volume)
!Returns the volume of the orthotope.
         implicit none
         real(8):: volume                            !out: volume of the orthotope
         class(orthotope_t), intent(in):: this       !in: orthotope
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer(INTL):: i

         errc=0; volume=0d0
         if(this%num_dim.gt.0) then
          volume=1d0
          do i=1,this%num_dim; volume=volume*this%extent(i)%length(); enddo
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function OrthotopeVolume
!---------------------------------------------------------------
        subroutine OrthotopeOverlap(this,orthotope,ierr,overlap)
!Returns the overlap of two orthotopes. If <overlap> is present,
!it will contain the result, otherwise <this> will contain the result.
         implicit none
         class(orthotope_t), intent(inout):: this             !inout: orthotope 1 (can be updated)
         class(orthotope_t), intent(in):: orthotope           !in: orthotope 2
         integer(INTD), intent(out), optional:: ierr          !out: error code
         type(orthotope_t), intent(inout), optional:: overlap !out: overlap (orthotope)
         integer(INTD):: errc
         integer(INTL):: i,n

         errc=0; n=this%num_dim
         if(n.gt.0.and.n.eq.orthotope%num_dim) then
          if(present(overlap)) then
           if(overlap%num_dim.ne.n) call overlap%orthotope_ctor(n,errc)
           if(errc.eq.0) then
            do i=1,n; call this%extent(i)%overlap(orthotope%extent(i),overlap%extent(i)); enddo
           else
            errc=2
           endif
          else
           do i=1,n; call this%extent(i)%overlap(orthotope%extent(i)); enddo
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine OrthotopeOverlap
!-----------------------------------------------------------
        subroutine OrthotopeUnion(this,orthotope,ierr,union)
!Returns the union of two orthotopes. If <union> is present,
!it will contain the result, otherwise <this> will contain the result.
         implicit none
         class(orthotope_t), intent(inout):: this           !inout: orthotope 1 (can be updated)
         class(orthotope_t), intent(in):: orthotope         !in: orthotope 2
         integer(INTD), intent(out), optional:: ierr        !out: error code
         type(orthotope_t), intent(inout), optional:: union !out: union (orthotope)
         integer(INTD):: errc
         integer(INTL):: i,n

         errc=0; n=this%num_dim
         if(n.gt.0.and.n.eq.orthotope%num_dim) then
          if(present(union)) then
           if(union%num_dim.ne.n) call union%orthotope_ctor(n,errc)
           if(errc.eq.0) then
            do i=1,n; call this%extent(i)%union(orthotope%extent(i),union%extent(i)); enddo
           else
            errc=2
           endif
          else
           do i=1,n; call this%extent(i)%union(orthotope%extent(i)); enddo
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine OrthotopeUnion
!--------------------------------------
        subroutine orthotope_dtor(this)
         implicit none
         type(orthotope_t):: this

         if(allocated(this%extent)) deallocate(this%extent)
         this%num_dim=0
         return
        end subroutine orthotope_dtor
![color_symmetry_t]==================================
        subroutine ColorSymmetryCtor(this,color,ierr)
!CTOR for color_symmetry_t.
         implicit none
         class(color_symmetry_t), intent(out):: this !out: generic symmetry object
         integer(INTD), intent(in):: color           !in: color (generic symmetry number)
         integer(INTD), intent(out), optional:: ierr !out: error code

         this%color=color
         if(present(ierr)) ierr=0
         return
        end subroutine ColorSymmetryCtor
!-----------------------------------------------------------
        function ColorSymmetryCompare(this,symm) result(cmp)
!Compares two symmetries.
         implicit none
         integer(INTD):: cmp                        !out: comparison result (see module dil_basic)
         class(color_symmetry_t), intent(in):: this !inout: generic symmetry object
         class(symmetry_t), intent(in):: symm       !in: another generic symmetry object

         select type(symm)
         class is(color_symmetry_t)
          if(this%color.lt.symm%color) then
           cmp=CMP_LT
          elseif(this%color.gt.symm%color) then
           cmp=CMP_GT
          else
           cmp=CMP_EQ
          endif
         class default
          cmp=CMP_ER
         end select
         return
        end function ColorSymmetryCompare
!------------------------------------------------------
        subroutine ColorSymmetryCombine(this,symm,ierr)
!Combines two symmetries.
         implicit none
         class(color_symmetry_t), intent(inout):: this !inout: generic symmetry object
         class(symmetry_t), intent(in):: symm          !in: another generic symmetry object
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: errc

         errc=0
         select type(symm)
         class is(color_symmetry_t)
          this%color=this%color+symm%color
         class default
          errc=1
         end select
         if(present(ierr)) ierr=errc
         return
        end subroutine ColorSymmetryCombine
!-------------------------------------------
        subroutine color_symmetry_dtor(this)
!DTOR for color_symmetry_t.
         implicit none
         type(color_symmetry_t):: this

         this%color=0
         return
        end subroutine color_symmetry_dtor
![spher_symmetry_t]==================================================
        subroutine SpherSymmetryCtor(this,orb_mom_tot,orb_mom_z,ierr)
!CTOR for spher_symmetry_t.
         implicit none
         class(spher_symmetry_t), intent(out):: this !out: spherical symmetry object
         integer(INTD), intent(in):: orb_mom_tot     !in: total orbital momentum (>=0), SYMMETRY_NONE means no symmetry
         integer(INTD), intent(in):: orb_mom_z       !in: Z-axis projection of the total orbital momentum, if defined
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(orb_mom_tot.ge.0.and.abs(orb_mom_z).le.orb_mom_tot) then !symmetry
          this%orb_moment=orb_mom_tot
          this%orb_z_proj=orb_mom_z
         elseif(orb_mom_tot.eq.SYMMETRY_NONE) then !no symmetry
          this%orb_moment=orb_mom_tot
          this%orb_z_proj=0
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine SpherSymmetryCtor
!--------------------------------------------------------------------
        function SpherSymmetryGetOrbMomentum(this,z_proj) result(res)
!Returns the total orbital momentum and its Z-axis projection.
         implicit none
         integer(INTD):: res                           !out: total orbital momentum
         class(spher_symmetry_t), intent(in):: this    !in: spherical symmetry object
         integer(INTD), intent(out), optional:: z_proj !out: Z-axis projection of the total orbital momentum

         res=this%orb_moment
         if(res.ge.0.and.present(z_proj)) z_proj=this%orb_z_proj
         return
        end function SpherSymmetryGetOrbMomentum
!-----------------------------------------------------------
        function SpherSymmetryCompare(this,symm) result(cmp)
!Compares two symmetries.
         implicit none
         integer(INTD):: cmp                        !out: comparison result (see module dil_basic)
         class(spher_symmetry_t), intent(in):: this !inout: spherical symmetry object
         class(symmetry_t), intent(in):: symm       !in: another spherical symmetry object
         integer(INTD):: om1,om2,mz1,mz2

         cmp=CMP_NC
         select type(symm)
         class is(spher_symmetry_t)
          om1=this%get_orb_momentum(mz1); om2=symm%get_orb_momentum(mz2)
          if(om1.lt.om2) then
           cmp=CMP_LT
          elseif(om1.gt.om2) then
           cmp=CMP_GT
          else
           if(mz1.lt.mz2) then
            cmp=CMP_LT
           elseif(mz1.gt.mz2) then
            cmp=CMP_GT
           else
            cmp=CMP_EQ
           endif
          endif
         class default
          cmp=CMP_ER
         end select
         return
        end function SpherSymmetryCompare
!------------------------------------------------------
        subroutine SpherSymmetryCombine(this,symm,ierr) !`Reimplement properly
!Combines two symmetries by finding a common lower symmetry.
         implicit none
         class(spher_symmetry_t), intent(inout):: this !inout: spherical symmetry object
         class(symmetry_t), intent(in):: symm          !in: another spherical symmetry object
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: errc

         errc=0
         select type(symm)
         class is(spher_symmetry_t)
          if(symm%get_orb_momentum().ne.this%get_orb_momentum()) call this%spher_symmetry_ctor(SYMMETRY_NONE,0,errc)
         class default
          errc=1
         end select
         if(present(ierr)) ierr=errc
         return
        end subroutine SpherSymmetryCombine
!----------------------------------------------------
        subroutine SpherSymmetryPrintIt(this,dev_out)
!Prints the symmetry information.
         implicit none
         class(spher_symmetry_t), intent(in):: this    !in: spherical symmetry object
         integer(INTD), intent(in), optional:: dev_out !in: output device (defaults to screen)
         integer(INTD):: devo

         devo=6; if(present(dev_out)) devo=dev_out
         write(devo,'("Spherical Symmetry (total,Z): ",i5,1x,i5)') this%orb_moment,this%orb_z_proj
         return
        end subroutine SpherSymmetryPrintIt
!-------------------------------------------
        subroutine spher_symmetry_dtor(this)
!DTOR for spher_symmetry_t.
         implicit none
         type(spher_symmetry_t):: this

         this%orb_moment=SYMMETRY_NONE; this%orb_z_proj=0
         return
        end subroutine spher_symmetry_dtor
![basis_func_supp_t]================================
        subroutine BasisFuncSuppCtorEmpty(this,ierr)
!Constructs an empty basis function support.
         implicit none
         class(basis_func_supp_t), intent(out):: this !out: basis function support
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         errc=0; this%supp_dim=0
         if(present(ierr)) ierr=errc
         return
        end subroutine BasisFuncSuppCtorEmpty
!---------------------------------------------------------------------
        subroutine BasisFuncSuppCtorReal(this,center,support_box,ierr)
!Constructs a real basis function support.
         implicit none
         class(basis_func_supp_t), intent(out):: this !out: basis function support
         type(real_vec_t), intent(in):: center        !in: effective center of the basis function support
         type(orthotope_t), intent(in):: support_box  !in: containing box (support)
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc,n

         n=center%dimsn()
         if(n.ge.0.and.support_box%dimsn().eq.n) then
          this%supp_dim=n
          this%center=center
          this%supp_box=support_box
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine BasisFuncSuppCtorReal
!---------------------------------------------------------
        function BasisFuncSuppIsSet(this,ierr) result(ans)
!Returns TRUE if the basis function support is set.
         implicit none
         logical:: ans                               !out: answer
         class(basis_func_supp_t), intent(in):: this !in: basis function support
         integer(INTD), intent(out), optional:: ierr !out: error code

         ans=(this%supp_dim.ge.0)
         if(present(ierr)) ierr=0
         return
        end function BasisFuncSuppIsSet
!-----------------------------------------------------------
        function BasisFuncSuppDimsn(this,ierr) result(dimsn)
!Returns the support dimension (>=0).
         implicit none
         integer(INTD):: dimsn                       !out: support space dimension
         class(basis_func_supp_t), intent(in):: this !in: basis function support
         integer(INTD), intent(out), optional:: ierr !out: error code

         dimsn=this%supp_dim
         if(present(ierr)) then
          if(dimsn.lt.0) ierr=1
         endif
         return
        end function BasisFuncSuppDimsn
!--------------------------------------------
        subroutine basis_func_supp_dtor(this)
!DTOR for basis_func_supp_t.
         implicit none
         type(basis_func_supp_t):: this

         this%supp_dim=-1
         return
        end subroutine basis_func_supp_dtor
![basis_func_gauss_t]======================================================
        subroutine BasisFuncGaussCtor(this,orb_moment,exponents,coefs,ierr)
!Sets up a Gaussian basis function. If it is already defined on input,
!it will automatically be destructed prior to redefinition.
         implicit none
         class(basis_func_gauss_t), intent(out):: this !out: basis function
         integer(INTD), intent(in):: orb_moment        !in: orbital momentum
         real(8), intent(in):: exponents(1:)           !in: exponents
         complex(8), intent(in):: coefs(1:)            !in: contraction coefficients
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: np,errc

         errc=0
         if(orb_moment.ge.0) then
          np=size(exponents) !number of primitives
          if(np.eq.size(coefs)) then
           allocate(this%exponent(1:np),STAT=errc)
           if(errc.eq.0) then
            allocate(this%coef(1:np),STAT=errc)
            if(errc.eq.0) then
             this%num_prims=np
             this%orb_moment=orb_moment
            else
             deallocate(this%exponent)
             errc=4
            endif
           else
            errc=3
           endif
          else
           errc=2
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine BasisFuncGaussCtor
!---------------------------------------------
        subroutine basis_func_gauss_dtor(this)
         implicit none
         type(basis_func_gauss_t):: this

         if(allocated(this%coef)) deallocate(this%coef)
         if(allocated(this%exponent)) deallocate(this%exponent)
         this%num_prims=0; this%orb_moment=-1
         return
        end subroutine basis_func_gauss_dtor
![basis_func_t]=======================================================
        subroutine BasisFuncCtor(this,basis_kind,ierr,basis_func,symm)
!Sets up a basis function of a given kind. If the basis function is already set,
!it will be redefined (no non-trivial destruction is assumed). If no basis
!function is passed here, an abstract basis function of <basis_kind> is assumed.
         implicit none
         class(basis_func_t), intent(out):: this                             !out: basis function
         integer(INTD), intent(in):: basis_kind                              !in: basis kind
         integer(INTD), intent(out), optional:: ierr                         !out: error code
         class(basis_func_supp_t), intent(in), target, optional:: basis_func !in: specification of the basis function (external)
         class(symmetry_t), intent(in), target, optional:: symm              !in: optional basis function symmetry (external)
         integer(INTD):: errc

         errc=0
         this%basis_func_p=>NULL(); this%symm_p=>NULL()
         if(basis_kind.ne.BASIS_NONE) then
          this%basis_kind=basis_kind
          if(present(basis_func)) this%basis_func_p=>basis_func
          if(present(symm)) this%symm_p=>symm
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine BasisFuncCtor
!---------------------------------------------------------------------------------------
        function BasisFuncGetBasisFunc(this,ierr,basis_func_p,symm_p) result(basis_kind)
!Returns basis function description. Note that the only guaranteed
!information to be returned is the <basis_kind>.
         implicit none
         integer(INTD):: basis_kind                                 !out: basis function kind (see top)
         class(basis_func_t), intent(in):: this                     !in: basis function
         integer(INTD), intent(out), optional:: ierr                !out: error code
         class(basis_func_supp_t), pointer, optional:: basis_func_p !out: polymorphic pointer to the basis function description
         class(symmetry_t), pointer, optional:: symm_p              !out: polymorphic pointer to the basis function symmetry (if any)
         integer(INTD):: errc

         errc=0; basis_kind=this%basis_kind
         if(basis_kind.ne.BASIS_NONE) then
          if(present(basis_func_p)) basis_func_p=>this%basis_func_p !does not have to be defined
          if(present(symm_p)) symm_p=>this%symm_p                   !does not have to be defined
         else
          if(present(basis_func_p)) basis_func_p=>NULL()
          if(present(symm_p)) symm_p=>NULL()
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function BasisFuncGetBasisFunc
!---------------------------------------
        subroutine basis_func_dtor(this)
!DTOR for basis_func_t.
         implicit none
         type(basis_func_t):: this

         this%basis_kind=BASIS_NONE
         this%basis_func_p=>NULL()
         this%symm_p=>NULL()
         return
        end subroutine basis_func_dtor
![subspace_basis_t]==================================
        subroutine SubspaceBasisCtor(this,dimsn,ierr)
!Creates an empty subspace basis. If the subspace basis is defined on input,
!it will automatically be destructed prior to redefinition.
         implicit none
         class(subspace_basis_t), intent(out):: this !out: empty subspace basis
         integer(INTL), intent(in):: dimsn           !in: dimension of the subspace
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(dimsn.gt.0) then
          allocate(this%basis_func(1:dimsn),STAT=errc)
          if(errc.eq.0) then
           this%space_dim=dimsn
          else
           errc=2
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine SubspaceBasisCtor
!------------------------------------------------------
        function SubspaceBasisDimsn(this) result(dimsn)
!Returns the dimension of the subspace basis.
         implicit none
         integer(INTL):: dimsn                      !out: subspace basis dimension
         class(subspace_basis_t), intent(in):: this !in: subspace basis

         dimsn=this%space_dim
         return
        end function SubspaceBasisDimsn
!------------------------------------------------------------------
        function SubspaceBasisSuppDimsn(this,ierr) result(supp_dim)
!Returns the dimension of the supporting real space.
         implicit none
         integer(INTD):: supp_dim                    !out: support dimension
         class(subspace_basis_t), intent(in):: this  !in: subspace basis
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(this%space_dim.gt.0) then
          supp_dim=this%supp_dim
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function SubspaceBasisSuppDimsn
!------------------------------------------------------------------------------------------
        subroutine SubspaceBasisSetBasisFunc(this,func_num,basis_kind,ierr,basis_func,symm)
!Sets up a specific basis function in the subspace basis.
         implicit none
         class(subspace_basis_t), intent(inout):: this                       !inout: subspace basis
         integer(INTL), intent(in):: func_num                                !in: basis function number in the subspace basis
         integer(INTD), intent(in):: basis_kind                              !in: basis kind
         integer(INTD), intent(out), optional:: ierr                         !out: error code
         class(basis_func_supp_t), intent(in), target, optional:: basis_func !in: specific basis function (external)
         class(symmetry_t), intent(in), target, optional:: symm              !in: basis function symmetry (external)
         integer(INTD):: errc
         integer(INTL):: n
         logical:: specific

         errc=0; n=this%dimsn(); specific=present(basis_func)
         if(n.gt.0) then
          if(func_num.gt.0.and.func_num.le.n) then
           if(this%supp_dim.eq.0.and.specific) then
            if(basis_func%supp_dim.gt.0) this%supp_dim=basis_func%supp_dim
           endif
           if(specific) then
            if(this%supp_dim.eq.basis_func%supp_dim) then
             if(present(symm)) then
              call this%basis_func(func_num)%basis_func_ctor(basis_kind,errc,basis_func,symm)
             else
              call this%basis_func(func_num)%basis_func_ctor(basis_kind,errc,basis_func)
             endif
             if(errc.ne.0) errc=5
            else
             errc=4
            endif
           else
            if(present(symm)) then
             call this%basis_func(func_num)%basis_func_ctor(basis_kind,errc,symm=symm)
            else
             call this%basis_func(func_num)%basis_func_ctor(basis_kind,errc)
            endif
            if(errc.ne.0) errc=3
           endif
          else
           errc=2
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine SubspaceBasisSetBasisFunc
!--------------------------------------------------------------------------------
        function SubspaceBasisGetBasisFunc(this,func_num,ierr) result(basis_func)
!Returns a pointer to a specific basis function in the subspace basis.
         implicit none
         class(basis_func_t), pointer:: basis_func          !out: pointer to a specific basis function
         class(subspace_basis_t), intent(in), target:: this !in: subspace basis
         integer(INTL), intent(in):: func_num               !in: basis function number
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc
         integer(INTL):: n

         errc=0; n=this%dimsn()
         if(n.gt.0) then
          if(func_num.gt.0.and.func_num.le.n) then
           basis_func=>this%basis_func(func_num)
          else
           errc=2
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function SubspaceBasisGetBasisFunc
!------------------------------------------------------------
        subroutine SubspaceBasisFinalize(this,ierr,num_undef)
!Finalizes the subspace basis setup after all basis functions have been set.
!Computes the average center and support in the real space in case of
!real space supported basis sets. Also computes the overall symmetry group irrep.
!If not all basis functions were set, the <num_undef> argument will return
!the number of unset basis functions and an error code.
         implicit none
         class(subspace_basis_t), intent(inout):: this    !inout: subspace basis
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTL), intent(out), optional:: num_undef !out: number of undefined basis functions
         integer(INTD):: errc
         integer(INTL):: i,n,nun
         class(basis_func_supp_t), pointer:: bas_func
         class(symmetry_t), pointer:: symm
         logical:: initb

         errc=0; nun=-1; n=this%dimsn()
         if(n.gt.0) then
          nun=0; initb=.FALSE.
          bloop: do i=1,n
           if(this%basis_func(i)%basis_kind.ne.BASIS_NONE) then
            bas_func=>this%basis_func(i)%basis_func_p
            symm=>this%basis_func(i)%symm_p
            if(associated(bas_func)) then
             if(.not.initb) then
              this%supp_dim=bas_func%supp_dim
              if(this%supp_dim.gt.0) then
               call this%center%real_vec_ctor(int(this%supp_dim,INTL),errc)
               if(errc.ne.0) then; errc=9; exit bloop; endif
               this%center=bas_func%center
               call this%supp_box%orthotope_ctor(int(this%supp_dim,INTL),errc)
               if(errc.ne.0) then; errc=8; exit bloop; endif
               this%supp_box=bas_func%supp_box
              endif
              initb=.TRUE.
             else
              if(bas_func%supp_dim.eq.this%supp_dim) then
               call this%center%add(bas_func%center,errc)
               if(errc.ne.0) then; errc=7; exit bloop; endif
               call this%supp_box%union(bas_func%supp_box,errc)
               if(errc.ne.0) then; errc=6; exit bloop; endif
              else
               errc=5; exit bloop !support dimension mismatch
              endif
             endif
            endif
            if(associated(symm)) then
             if(.not.allocated(this%symm)) then
              allocate(this%symm,source=symm,STAT=errc)
              if(errc.ne.0) then; errc=4; exit bloop; endif
             else
              call this%symm%combine(symm,errc)
              if(errc.ne.0) then; errc=3; exit bloop; endif
             endif
            endif
           else !undefined basis function
            nun=nun+1
           endif
          enddo bloop
          if(errc.eq.0) then
           if(nun.lt.n.and.this%supp_dim.gt.0) call this%center%scale(1d0/real(n-nun,8),errc)
           if(nun.ne.0) errc=2 !some basis functions are still undefined
          endif
         else
          errc=1
         endif
         if(present(num_undef)) num_undef=nun
         if(present(ierr)) ierr=errc
         return
        end subroutine SubspaceBasisFinalize
!----------------------------------------------------------------
        function SubspaceBasisGetSymmetry(this,ierr) result(symm)
!Returns a pointer to the subspace basis symmetry object.
         implicit none
         class(symmetry_t), pointer:: symm                  !out: subspace basis symmetry
         class(subspace_basis_t), intent(in), target:: this !in: subspace basis
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc

         errc=0
         if(this%space_dim.gt.0) then
          symm=>NULL(); if(allocated(this%symm)) symm=>this%symm
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function SubspaceBasisGetSymmetry
!----------------------------------------------------------------
        function SubspaceBasisGetCenter(this,ierr) result(center)
!Returns a pointer to the subspace basis center in real space.
         implicit none
         class(real_vec_t), pointer:: center                !out: effective center of the subspace basis support
         class(subspace_basis_t), intent(in), target:: this !in: subspace basis
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc

         errc=0
         if(this%space_dim.gt.0.and.this%supp_dim.gt.0) then
          center=>this%center
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function SubspaceBasisGetCenter
!-------------------------------------------------------------------
        function SubspaceBasisGetSupport(this,ierr) result(supp_box)
!Returns a pointer to the subspace basis support box.
         implicit none
         class(orthotope_t), pointer:: supp_box             !out: subspace basis support (orthotope)
         class(subspace_basis_t), intent(in), target:: this !in: subspace basis
         integer(INTD), intent(out), optional:: ierr        !out: error code
         integer(INTD):: errc

         errc=0
         if(this%space_dim.gt.0.and.this%supp_dim.gt.0) then
          supp_box=>this%supp_box
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function SubspaceBasisGetSupport
!-------------------------------------------
        subroutine subspace_basis_dtor(this)
         implicit none
         type(subspace_basis_t):: this

         if(allocated(this%basis_func)) deallocate(this%basis_func)
         if(allocated(this%symm)) deallocate(this%symm)
         this%space_dim=0; this%supp_dim=0
         return
        end subroutine subspace_basis_dtor
![subspace_t]=============================================
        subroutine SubspaceCtorBase(this,id,ierr,subrange)
!Initializes a subspace. If the subspace is defined on input,
!it will automatically be destructed prior to re-initialization.
         implicit none
         class(subspace_t), intent(out):: this            !out: empty subspace
         integer(INTL), intent(in):: id                   !in: subspace id (must be non-negative)
         integer(INTD), intent(out), optional:: ierr      !out: error code
         type(seg_int_t), intent(in), optional:: subrange !in: defining subrange of basis vectors, specified as a semi-interval
         integer(INTD):: errc

         errc=0
         if(id.ge.0) then
          this%subspace_id=id   !unique subspace id
          this%supp_dim=0       !will be inferred from the first registered basis supported in real space
          this%max_resolution=0 !will be inferred from the first registered basis
          if(present(subrange)) this%basis_subrange=subrange !defining basis subrange
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine SubspaceCtorBase
!----------------------------------------------------
        function SubspaceIsSet(this,ierr) result(ans)
!Returns TRUE if the subspace is set.
         implicit none
         logical:: ans                               !out: answer
         class(subspace_t), intent(in):: this        !in: subspace
         integer(INTD), intent(out), optional:: ierr !out: error code

         ans=(this%subspace_id.ge.0)
         if(present(ierr)) ierr=0
         return
        end function SubspaceIsSet
!---------------------------------------------------
        function SubspaceGetId(this,ierr) result(id)
         implicit none
         integer(INTL):: id                          !out: subspace id
         class(subspace_t), intent(in):: this        !in: subspace
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         id=this%subspace_id
         if(present(ierr)) ierr=errc
         return
        end function SubspaceGetId
!--------------------------------------------------------------------
        function SubspaceGetBasisSubrange(this,ierr) result(subrange)
!Returns the subspace basis subrange.
         implicit none
         type(seg_int_t):: subrange                  !out: basis subrange
         class(subspace_t), intent(in):: this        !in: subspace
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(this%is_set()) then
          subrange=this%basis_subrange
         else
          call subrange%set(0_INTL,0_INTL)
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function SubspaceGetBasisSubrange
!--------------------------------------------------------------
        function SubspaceGetSuppDim(this,ierr) result(supp_dim)
         implicit none
         integer(INTD):: supp_dim                    !out: subspace support dimension
         class(subspace_t), intent(in):: this        !in: subspace
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(this%subspace_id.ge.0) then
          supp_dim=this%supp_dim
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function SubspaceGetSuppDim
!-------------------------------------------------------------------
        function SubspaceGetMaxResolution(this,ierr) result(max_res)
         implicit none
         integer(INTL):: max_res                     !out: max resolution
         class(subspace_t), intent(in):: this        !in: subspace
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(this%subspace_id.ge.0) then
          max_res=this%max_resolution
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function SubspaceGetMaxResolution
!--------------------------------------------------------------
        function SubspaceCompareRange(this,another) result(cmp)
!Compares the basis range of the current subspace with another subspace.
         implicit none
         integer(INTD):: cmp                     !out: comparison result: {CMP_EQ,CMP_LT,CMP_GT,CMP_OV,CMP_ER}
         class(subspace_t), intent(in):: this    !in: subspace 1
         class(subspace_t), intent(in):: another !in: subspace 2
         integer(INTL):: l1,u1,l2,u2

         if(this%is_set().and.another%is_set()) then
          l1=this%basis_subrange%lower_bound()+1_INTL; u1=this%basis_subrange%upper_bound()
          l2=another%basis_subrange%lower_bound()+1_INTL; u2=another%basis_subrange%upper_bound()
          if(l1.lt.l2) then
           if(u1.lt.l2) then; cmp=CMP_LT; else; cmp=CMP_OV; endif
          elseif(l1.gt.l2) then
           if(u2.lt.l1) then; cmp=CMP_GT; else; cmp=CMP_OV; endif
          else
           if(u1.eq.u2) then; cmp=CMP_EQ; else; cmp=CMP_OV; endif
          endif
         else
          cmp=CMP_ER
         endif
         return
        end function SubspaceCompareRange
!--------------------------------------------------------
        subroutine SubspaceRegisterBasis(this,basis,ierr)
!Registers a new subspace basis. The new subspace basis must have
!the same support dimension as any previously registered basis.
         implicit none
         class(subspace_t), intent(inout):: this             !inout: subspace
         class(subspace_basis_t), intent(in), target:: basis !in: new subspace basis
         integer(INTD), intent(out), optional:: ierr         !out: error code
         integer(INTD):: errc,ier,sd
         integer(INTL):: n
         type(list_iter_t):: basis_it
         !type(subspace_basis_t), pointer:: dptr !debug

         if(this%subspace_id.ge.0) then
          errc=basis_it%init(this%bases)
          if(errc.eq.GFC_SUCCESS) then
           errc=basis_it%reset()
           if(errc.eq.GFC_SUCCESS) then
            n=basis%dimsn()
            if(n.gt.0) then
             sd=basis%supp_dimsn()
             if(sd.gt.0.and.this%supp_dim.eq.0) this%supp_dim=sd
             if(sd.eq.this%supp_dim) then
              !select type(basis); type is(subspace_basis_t); size_=size_of(basis); dptr=>basis; ptr_=c_loc(dptr); end select !debug
              !call dump_bytes(ptr_,size_,'dump1') !debug
              errc=basis_it%append(basis)
              if(errc.eq.GFC_SUCCESS) then
               this%max_resolution=max(this%max_resolution,n)
              else
               errc=6
              endif
             else
              errc=5
             endif
            else
             errc=4
            endif
           else
            errc=3
           endif
           ier=basis_it%release()
          else
           errc=2
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine SubspaceRegisterBasis
!---------------------------------------------------------------------
        function SubspaceResolve(this,ierr,bas_pred_f) result(basis_p)
!Returns a pointer to the subspace basis satisfying a certain (optional) condition.
!The condition is specified via a GFC predicate object. If no condition is supplied,
!the very first basis will be returned (max resolution basis).
         implicit none
         class(subspace_basis_t), pointer:: basis_p                   !out: pointer to the subspace basis
         class(subspace_t), intent(in):: this                         !in: subspace
         integer(INTD), intent(out), optional:: ierr                  !out: error code
         class(gfc_predicate_t), intent(inout), optional:: bas_pred_f !in: optional predicate for search
         integer(INTD):: errc,ier
         type(list_iter_t):: basis_it
         class(gfc_cont_elem_t), pointer:: cont_elem
         class(*), pointer:: elem_value

         if(this%subspace_id.ge.0) then
          if(this%max_resolution.gt.0) then
           errc=basis_it%init(this%bases)
           if(errc.eq.GFC_SUCCESS) then
            errc=basis_it%reset()
            if(errc.eq.GFC_SUCCESS) then
             if(present(bas_pred_f)) then
              ier=basis_it%scanf(return_each=.TRUE.,predicate_f=bas_pred_f)
             else
              ier=basis_it%scanf(return_each=.TRUE.)
             endif
             if(ier.eq.GFC_IT_ACTIVE) then !found
              cont_elem=>basis_it%pointee(errc)
              if(errc.eq.GFC_SUCCESS) then
               elem_value=>cont_elem%get_value(errc)
               if(errc.eq.GFC_SUCCESS) then
                select type(elem_value)
                class is(subspace_basis_t)
                 basis_p=>elem_value
                class default
                 errc=7
                end select
               else
                errc=6
               endif
              else
               errc=5
              endif
             else
              errc=2 !not found
             endif
            else
             errc=4
            endif
            ier=basis_it%release()
           else
            errc=3
           endif
          else
           errc=2 !not found
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function SubspaceResolve
!-----------------------------------------------
        subroutine SubspacePrintIt(this,dev_out)
!Prints the subspace information.
         implicit none
         class(subspace_t), intent(in):: this          !in: subspace
         integer(INTD), intent(in), optional:: dev_out !in: output device (defaults to screen)
         integer(INTD):: devo

         devo=6; if(present(dev_out)) devo=dev_out
         write(devo,'("Subspace ",i10,": [",i9,":",i9,"] -> ",i9," (dim max)")') this%subspace_id,&
              &this%basis_subrange%lower_bound()+1_INTL,this%basis_subrange%upper_bound(),this%max_resolution
         return
        end subroutine SubspacePrintIt
!-------------------------------------
        subroutine subspace_dtor(this)
         implicit none
         type(subspace_t):: this
         type(list_iter_t):: basis_it
         integer(INTD):: errc

         errc=basis_it%init(this%bases)
         if(errc.eq.GFC_SUCCESS) then
          errc=basis_it%delete_all()
          errc=basis_it%release()
         endif
         this%max_resolution=0
         this%supp_dim=0
         this%subspace_id=-1
         return
        end subroutine subspace_dtor
![h_space_t]===============================================================
        subroutine HSpaceCtorSimple(this,full_basis,ierr,branch_fac,metric)
!Constructs a simple hierarchical vector space representation with a subspace aggregation tree.
!The original basis functions will be hierarchically aggregated into larger subspaces,
!up to the full space. If symmetry is present, the subspaces will be aligned on the
!symmetry boundaries. Each subspace in the subspace aggregation tree will have its
!unique id and an associated max-resolution basis, the latter consisting of a subset
!of the original full basis (later on, reduced basis sets can be defined for any subspace).
!Each subspace can also be directly accessed by its id in the vector tree <this%subspaces>.
!Storage complexity O(NlogN), where N is the original full space dimension.
         implicit none
         class(h_space_t), intent(out), target:: this             !out: hierarchical representation of the vector space
         class(subspace_basis_t), intent(in), target:: full_basis !in: full basis of the vector space, {Psi_i}
         integer(INTD), intent(out), optional:: ierr              !out: error code
         integer(INTD), intent(in), optional:: branch_fac         !in: aggregation tree branching factor (defaults to 2)
         complex(8), intent(in), optional, target:: metric(:,:)   !in: metric tensor: g_ij=<Psi_i|Psi_j>: Hermitian complex matrix
         logical:: split
         integer(INTD):: i,l,m,brf,tlevel,errc
         integer(INTL):: nbnd
         integer(INTL), allocatable:: bndr(:)
         type(vec_tree_iter_t):: vt_it
         type(seg_int_t):: seg
         type(seg_int_t), allocatable:: segs(:)
         type(subspace_t):: subspace
         class(subspace_t), pointer:: subsp
         type(subspace_basis_t):: basis
         class(subspace_basis_t), pointer:: basp
         class(*), pointer:: up

         errc=0
         this%num_subspaces=0_INTL; this%space_dim=full_basis%dimsn()
         if(this%space_dim.gt.0) then
!Construct the subspace aggregation tree (SAT) by recursive basis splitting:
          allocate(bndr(1:this%space_dim),STAT=errc)
          if(errc.eq.0) then
           if(present(branch_fac)) then; brf=branch_fac; else; brf=2; endif
           if(brf.ge.2) then
            allocate(segs(1:brf),STAT=errc)
            if(errc.eq.0) then
             errc=vt_it%init(this%subspaces)
             if(errc.eq.GFC_SUCCESS) then
 !Add individual 1-dimensional subspaces first (span of each basis vector):
              do while(this%num_subspaces.lt.this%space_dim) !subspaces [0:N-1] are spans of individual basis vectors [1..N]
               call seg%set(this%num_subspaces,this%num_subspaces+1_INTL,errc); if(errc.ne.0) exit
               call construct_subspace_basis(full_basis,seg,basis,errc); if(errc.ne.0) exit
               call subspace%subspace_ctor(this%num_subspaces,errc,seg); if(errc.ne.0) exit
               errc=vt_it%append(subspace); if(errc.ne.GFC_SUCCESS) exit
               errc=vt_it%reset_back(); if(errc.ne.GFC_SUCCESS) exit
               up=>vt_it%get_value(errc); if(errc.ne.GFC_SUCCESS) exit
               subsp=>NULL(); select type(up); class is(subspace_t); subsp=>up; end select
               if(.not.associated(subsp)) then; errc=1; exit; endif
               !write(*,'("#DEBUG(h_space_t.ctor): Registering 1d subspace basis ...")') !debug
               call subsp%register_basis(basis,errc); if(errc.ne.0) exit
               !write(*,'("#DEBUG(h_space_t.ctor): Registered.")'); stop !debug
               this%num_subspaces=this%num_subspaces+1_INTL
              enddo
 !Add the tree root (full vector space defined by the provided full basis):
              if(errc.eq.0) then
               if(this%space_dim.gt.1) then !more than one basis function in the full space
                call seg%set(0_INTL,this%space_dim,errc) !full space: (0:N] = [1:N]
                if(errc.eq.0) call subspace%subspace_ctor(this%num_subspaces,errc,seg)
                if(errc.eq.0) errc=vt_it%append(subspace)
                if(errc.eq.GFC_SUCCESS) errc=vt_it%reset_back();
                if(errc.eq.GFC_SUCCESS) up=>vt_it%get_value(errc)
                if(errc.eq.GFC_SUCCESS) then
                 subsp=>NULL(); select type(up); class is(subspace_t); subsp=>up; end select
                 if(associated(subsp)) then
                  call subsp%register_basis(full_basis,errc)
                  if(errc.eq.0) errc=vt_it%add_leaf(this%num_subspaces)
                  if(errc.eq.GFC_SUCCESS) this%num_subspaces=this%num_subspaces+1_INTL
                 else
                  errc=2
                 endif
                endif
 !Recursively split the full space into subspaces while respecting symmetry boundaries, if any:
                if(errc.eq.0) then
                 tlevel=-1; split=.TRUE.
                 tloop: do while(split)
                  split=.FALSE.; tlevel=tlevel+1
                  !write(*,'("Splitting Tree Level ",i4)') tlevel !debug
                  do while(errc.eq.GFC_SUCCESS)
  !Process the current tree vertex:
                   up=>vt_it%get_value(errc); if(errc.ne.GFC_SUCCESS) exit tloop
                   subsp=>NULL(); select type(up); class is(subspace_t); subsp=>up; end select
                   if(.not.associated(subsp)) then; errc=3; exit tloop; endif
                   m=int(min(subsp%get_max_resolution(errc),int(brf,INTL)),INTD); if(errc.ne.0) exit tloop
                   if(m.gt.1) then !aggregate subspace, thus can be split further
                    seg=subsp%get_basis_subrange(errc); if(errc.ne.0) exit tloop
                    basp=>subsp%resolve(errc); if(errc.ne.0) exit tloop
                    if(associated(basp)) then
                     call set_symmetry_boundaries(basp,nbnd,bndr,errc); if(errc.ne.0) exit tloop
                    else
                     errc=4; exit tloop
                    endif
                    !write(*,'("Range ")',ADVANCE='NO'); call seg%print_range() !debug
                    !write(*,'(" to be split into ",i6," segments")') m !debug
                    !print *,'with boundaries ',bndr(1:nbnd) !debug
                    if(nbnd.gt.0) then
                     call seg%split(m,segs,errc,bndr(1:nbnd))
                    else
                     call seg%split(m,segs,errc)
                    endif
                    do l=1,m
                     !write(*,'("Subrange ",i3,": ")',ADVANCE='NO') l; call segs(l)%print_range() !debug
                     if(segs(l)%length().gt.1_INTL) then !aggregate subspace
                      call construct_subspace_basis(full_basis,segs(l),basis,errc); if(errc.ne.0) exit tloop
                      call subspace%subspace_ctor(this%num_subspaces,errc,segs(l)); if(errc.ne.0) exit tloop
                      errc=vt_it%append(subspace); if(errc.ne.GFC_SUCCESS) exit tloop
                      up=>vt_it%element_value(this%num_subspaces,errc); if(errc.ne.GFC_SUCCESS) exit tloop
                      subsp=>NULL(); select type(up); class is(subspace_t); subsp=>up; end select
                      if(.not.associated(subsp)) then; errc=5; exit tloop; endif
                      call subsp%register_basis(basis,errc); if(errc.ne.0) exit tloop
                      errc=vt_it%add_leaf(this%num_subspaces,no_move=.TRUE.); if(errc.ne.GFC_SUCCESS) exit tloop
                      this%num_subspaces=this%num_subspaces+1_INTL
                     else !1-dimensional subspace
                      errc=vt_it%add_leaf(segs(l)%lower_bound(),no_move=.TRUE.); if(errc.ne.GFC_SUCCESS) exit tloop
                     endif
                    enddo
                    split=.TRUE.
                   endif
  !Move to the right cousin (within the current tree level):
                   errc=vt_it%move_to_cousin()
                  enddo
                  if(errc.eq.GFC_NO_MOVE) then; errc=GFC_SUCCESS; else; exit tloop; endif
  !Move to the children level:
                  do while(errc.eq.GFC_SUCCESS); errc=vt_it%move_to_cousin(to_previous=.TRUE.); enddo
                  if(errc.eq.GFC_NO_MOVE) then; errc=GFC_SUCCESS; else; exit tloop; endif
                  cloop: do
                   errc=vt_it%move_to_child(); if(errc.eq.GFC_SUCCESS) exit cloop
                   if(errc.eq.GFC_NO_MOVE) then
                    errc=vt_it%move_to_cousin(); if(errc.eq.GFC_SUCCESS) cycle cloop
                    if(errc.eq.GFC_NO_MOVE.and.(.not.split)) errc=GFC_SUCCESS
                   endif
                   exit tloop
                  enddo cloop
                 enddo tloop
                else
                 errc=6
                endif
               endif
              else
               errc=7
              endif
              i=vt_it%release(); if(errc.eq.0.and.i.ne.GFC_SUCCESS) errc=8
             else
              errc=9
             endif
             if(allocated(segs)) deallocate(segs)
            else
             errc=10
            endif
           else
            errc=11
           endif
           if(allocated(bndr)) deallocate(bndr)
          else
           errc=12
          endif
!Construct the support overlap matrix between all subspaces:
          if(errc.eq.0.and.present(metric)) then
           this%metric_p=>metric
           !`Construct overlaps
          endif
         else
          errc=13
         endif
         if(errc.ne.0) call h_space_dtor(this)
         if(present(ierr)) ierr=errc
         return

        contains

         subroutine set_symmetry_boundaries(bas,jnbnd,jbnd,jerr)
          class(subspace_basis_t), intent(in):: bas
          integer(INTL), intent(out):: jnbnd
          integer(INTL), intent(inout):: jbnd(1:)
          integer(INTD), intent(out):: jerr
          class(basis_func_t), pointer:: bfp
          class(basis_func_supp_t), pointer:: bfsp
          class(symmetry_t), pointer:: curr_symm,next_symm
          integer(INTL):: jdim,ji
          integer(INTD):: curr_bfk,next_bfk

          jdim=bas%dimsn(); jnbnd=0_INTL !will be the number of boundaries
          if(jdim.gt.0) then
           do ji=1,jdim-1
            bfp=>bas%get_basis_func(ji,jerr); if(jerr.ne.0) exit
            curr_bfk=bfp%get_basis_func(jerr,bfsp,curr_symm); if(jerr.ne.0) exit
            bfp=>bas%get_basis_func(ji+1_INTL,jerr); if(jerr.ne.0) exit
            next_bfk=bfp%get_basis_func(jerr,bfsp,next_symm); if(jerr.ne.0) exit
            !select type(curr_symm); class is(spher_symmetry_t); call curr_symm%print_it(); end select !debug
            !select type(next_symm); class is(spher_symmetry_t); call next_symm%print_it(); end select !debug
            if(curr_bfk.eq.next_bfk) then
             if(associated(curr_symm)) then
              if(associated(next_symm)) then
               if(curr_symm%compare(next_symm).ne.CMP_EQ) then; jnbnd=jnbnd+1_INTL; jbnd(jnbnd)=ji; endif !symmetry change boundary
              else
               jnbnd=jnbnd+1_INTL; jbnd(jnbnd)=ji !symmetry change boundary
              endif
             else
              if(associated(next_symm)) then; jnbnd=jnbnd+1_INTL; jbnd(jnbnd)=ji; endif !symmetry change boundary
             endif
            else
             jnbnd=jnbnd+1_INTL; jbnd(jnbnd)=ji !basis kind change boundary
            endif
           enddo
          else
           jerr=1
          endif
          return
         end subroutine set_symmetry_boundaries

         subroutine construct_subspace_basis(parbas,sgs,bas,jerr)
          implicit none
          class(subspace_basis_t), intent(in), target:: parbas !in: parental basis
          type(seg_int_t), intent(in):: sgs                    !in: defining integer semi-interval (subrange of basis functions)
          type(subspace_basis_t), intent(inout):: bas          !out: subspace basis being constructed based on the <sgs> subrange
          integer(INTD), intent(out):: jerr                    !out: error code
          integer(INTD):: jbfk
          integer(INTL):: jdim,jl,ju,jj
          class(basis_func_t), pointer:: jbfp
          class(basis_func_supp_t), pointer:: jbfsp
          class(symmetry_t), pointer:: jsmp

          jerr=0; jdim=sgs%length()
          if(jdim.gt.0) then
           jl=sgs%lower_bound()+1_INTL; ju=sgs%upper_bound()
           call bas%subspace_basis_ctor(jdim,jerr)
           if(jerr.eq.0) then
            do jj=jl,ju
             jbfp=>parbas%get_basis_func(jj,jerr); if(jerr.ne.0) exit
             jbfk=jbfp%get_basis_func(jerr,jbfsp,jsmp); if(jerr.ne.0) exit
             if(associated(jbfsp)) then
              if(associated(jsmp)) then
               call bas%set_basis_func(jj-jl+1_INTL,jbfk,jerr,jbfsp,jsmp)
              else
               call bas%set_basis_func(jj-jl+1_INTL,jbfk,jerr,jbfsp)
              endif
             else
              if(associated(jsmp)) then
               call bas%set_basis_func(jj-jl+1_INTL,jbfk,jerr,symm=jsmp)
              else
               call bas%set_basis_func(jj-jl+1_INTL,jbfk,jerr)
              endif
             endif
             if(jerr.ne.0) exit
            enddo
            if(jerr.eq.0) call bas%finalize(jerr)
           endif
          else
           jerr=1
          endif
          if(jerr.ne.0) call subspace_basis_dtor(bas)
          return
         end subroutine construct_subspace_basis

        end subroutine HSpaceCtorSimple
!--------------------------------------------------
        function HSpaceIsSet(this,ierr) result(res)
!Returns TRUE if the hierarchical vector space is set.
         implicit none
         logical:: res                               !out: result
         class(h_space_t), intent(in):: this         !in: hierarchical vector space representation
         integer(INTD), intent(out), optional:: ierr !out: eror code
         integer(INTD):: errc

         errc=0; res=(this%space_dim.gt.0)
         if(present(ierr)) ierr=errc
         return
        end function HSpaceIsSet
!--------------------------------------------------------
        function HSpaceGetSpaceDim(this,ierr) result(res)
!Returns the dimension of the vector space.
         implicit none
         integer(INTL):: res                         !out: vector space dimension
         class(h_space_t), intent(in):: this         !in: hierarchical vector space representation
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(this%space_dim.gt.0) then
          res=this%space_dim
         else
          res=0_INTL; errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function HSpaceGetSpaceDim
!------------------------------------------------------------
        function HSpaceGetNumSubspaces(this,ierr) result(res)
!Returns the total number of defined subspaces in the vector space.
         implicit none
         integer(INTL):: res                         !out: total number of subspaces, including the 1-dimensional basis subspaces
         class(h_space_t), intent(in):: this         !in: hierarchical vector space representation
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(this%num_subspaces.gt.0) then
          res=this%num_subspaces
         else
          res=0_INTL; errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function HSpaceGetNumSubspaces
!----------------------------------------------------------
        function HSpaceGetRootId(this,ierr) result(root_id)
!Returns the id of the full space (root of the subspace aggregation tree).
         implicit none
         integer(INTL):: root_id                     !out: root id (sequential offset)
         class(h_space_t), intent(in):: this         !in: hierarchical vector space
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,i
         type(vec_tree_iter_t):: vtit

         root_id=-1_INTL; errc=vtit%init(this%subspaces)
         if(errc.eq.GFC_SUCCESS) then
          root_id=vtit%get_root_id(errc)
          i=vtit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=i
         endif
         if(present(ierr)) ierr=errc
         return
        end function HSpaceGetRootId
!------------------------------------------------------------------------------------------------
        function HSpaceGetAncestorId(this,subspace_id,ancestor_distance,ierr) result(ancestor_id)
!Returns the id of a specific ancestor subspace.
         implicit none
         integer(INTL):: ancestor_id                   !out: ancestor id (sequential offset)
         class(h_space_t), intent(in):: this           !in: hierarchical vector space
         integer(INTL), intent(in):: subspace_id       !in: subspace id
         integer(INTD), intent(in):: ancestor_distance !in: distance to the requested ancestor (1:parent,2,3,...)
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: errc,i
         type(vec_tree_iter_t):: vtit

         ancestor_id=-1_INTL
         if(ancestor_distance.gt.0) then
          errc=vtit%init(this%subspaces)
          if(errc.eq.GFC_SUCCESS) then
           errc=vtit%move_to(subspace_id)
           if(errc.eq.GFC_SUCCESS) then
            errc=vtit%move_up(ancestor_distance)
            if(errc.eq.GFC_SUCCESS) ancestor_id=vtit%get_offset(errc)
           endif
           i=vtit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=i
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function HSpaceGetAncestorId
!---------------------------------------------------------------------------
        function HSpaceGetSubspaceLevel(this,subspace_id,ierr) result(level)
!Returns the distance from the root for a specific subspace.
         implicit none
         integer(INTD):: level                       !out: distance from the root (in hops)
         class(h_space_t), intent(inout):: this      !in: hierarchical vector space
         integer(INTL), intent(in):: subspace_id     !in: subspace id
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,i
         type(vec_tree_iter_t):: vtit

         level=-1; errc=vtit%init(this%subspaces)
         if(errc.eq.GFC_SUCCESS) then
          errc=vtit%move_to(subspace_id)
          if(errc.eq.GFC_SUCCESS) level=vtit%get_level(errc)
          i=vtit%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=i
         endif
         if(present(ierr)) ierr=errc
         return
        end function HSpaceGetSubspaceLevel
!---------------------------------------------------------------------------
        function HSpaceGetSubspace(this,subspace_id,ierr) result(subspace_p)
!Returns a pointer to the requested subspace from the hierarchical vector space.
         implicit none
         class(subspace_t), pointer:: subspace_p     !out: pointer to the requested subspace
         class(h_space_t), intent(in), target:: this !in: hierarchical vector space representation
         integer(INTL), intent(in):: subspace_id     !in: requested subspace id: [0..max]
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,i
         class(*), pointer:: up
         type(vec_tree_iter_t):: vt_it

         errc=0; subspace_p=>NULL()
         if(this%is_set()) then
          if(subspace_id.ge.0.and.subspace_id.lt.this%num_subspaces) then
           errc=vt_it%init(this%subspaces)
           if(errc.eq.GFC_SUCCESS) then
            errc=vt_it%move_to(subspace_id)
            if(errc.eq.GFC_SUCCESS) then
             up=>vt_it%get_value(errc)
             if(errc.eq.GFC_SUCCESS) then
              select type(up)
              class is(subspace_t)
               subspace_p=>up
              class default
               errc=7
              end select
             else
              errc=6
             endif
            else
             errc=5
            endif
            i=vt_it%release(); if(errc.eq.0.and.i.ne.GFC_SUCCESS) errc=4
           else
            errc=3
           endif
          else
           errc=2
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function HSpaceGetSubspace
!----------------------------------------------------------------
        function HSpaceGetAggrTree(this,ierr) result(aggr_tree_p)
!Returns a pointer to the subspace aggregation tree.
         implicit none
         class(vec_tree_t), pointer:: aggr_tree_p    !out: pointer to the subspace aggregation tree
         class(h_space_t), intent(in), target:: this !in: hierarchical vector space
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; aggr_tree_p=>NULL()
         if(this%is_set()) then
          aggr_tree_p=>this%subspaces
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function HSpaceGetAggrTree
!------------------------------------------------------------------------------------
        subroutine HSpaceGetLevelComposition(this,level,subspaces,num_subspaces,ierr)
!Returns an ordered list of subspaces (their ids) forming a given level of the subspace
!aggregation tree (direct sum decomposition of the complete space at the given level).
!If the <subspaces> array is already allocated, its length must be large enough.
         implicit none
         class(h_space_t), intent(in):: this                      !in: hierarchical vector space
         integer(INTD), intent(in):: level                        !in: subspace aggregation tree level (0 = root)
         integer(INTL), intent(inout), allocatable:: subspaces(:) !out: an ordered list of subspaces (their ids) forming the given tree level
         integer(INTL), intent(out):: num_subspaces               !out: number of subspaces in the list
         integer(INTD), intent(out), optional:: ierr              !out: error code
         integer(INTD):: errc,ier
         integer(INTL):: ln
         class(vec_tree_t), pointer:: aggr_tree
         type(vec_tree_iter_t):: tree
         class(*), pointer:: subspace

         errc=0; num_subspaces=0
         if(level.ge.0) then
 !Get the subspace aggregation tree:
          aggr_tree=>this%get_aggr_tree(errc)
          if(errc.eq.0) then
           errc=tree%init(aggr_tree)
           if(errc.eq.GFC_SUCCESS) then
 !Check the allocation status of the output array and allocate it if needed:
            ln=0; if(allocated(subspaces)) ln=size(subspaces)
            if(ln.le.0) then
 !Count the number of subspaces at the given tree level:
             errc=tree%find_first_of_level(level)
             if(errc.eq.GFC_SUCCESS) then
              ln=0
              do while(errc.eq.GFC_SUCCESS); ln=ln+1; errc=tree%move_to_cousin(); enddo
              if(errc.eq.GFC_NO_MOVE) then
               allocate(subspaces(1:ln),STAT=errc); if(errc.ne.0) errc=-12
              else
               errc=-11
              endif
             else
              errc=-10
             endif
            endif
 !Form the list of subspaces (their ids):
            if(errc.eq.0) then
             errc=tree%find_first_of_level(level)
             if(errc.eq.GFC_SUCCESS) then
              do while(errc.eq.GFC_SUCCESS)
               subspace=>tree%get_value(errc); if(errc.ne.GFC_SUCCESS) then; errc=-9; exit; endif
               select type(subspace)
               class is(subspace_t)
                num_subspaces=num_subspaces+1
                subspaces(num_subspaces)=subspace%get_id(errc)
                if(errc.ne.0) then; errc=-8; exit; endif
               class default
                errc=-7; exit
               end select
               errc=tree%move_to_cousin()
              enddo
              if(errc.eq.GFC_NO_MOVE) then; errc=0; else; errc=-6; endif
             else
              errc=-5
             endif
            endif
            ier=tree%release(); if(ier.ne.GFC_SUCCESS.and.errc.eq.0) errc=-4
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine HSpaceGetLevelComposition
!---------------------------------------------------------------------
        function HSpaceCompareSubspaces(this,id1,id2,ierr) result(cmp)
!Compares two subspaces from the hierarchical vector space.
         implicit none
         integer(INTD):: cmp                         !out: comparison result: {CMP_EQ,CMP_LT,CMP_GT,CMP_ER}
         class(h_space_t), intent(in):: this         !in: hierarchical vector space
         integer(INTL), intent(in):: id1             !in: subspace 1 id
         integer(INTL), intent(in):: id2             !in: subspace 2 id
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; cmp=CMP_ER
         if(id1.ge.0.and.id1.lt.this%num_subspaces.and.id2.ge.0.and.id2.lt.this%num_subspaces) then
          if(id1.eq.id2) then
           cmp=CMP_EQ
          else
           if(id1.lt.this%space_dim.and.id2.lt.this%space_dim) then
            if(id1.lt.id2) then; cmp=CMP_LT; else; cmp=CMP_GT; endif
           elseif(id1.ge.this%space_dim.and.id2.ge.this%space_dim) then
            if(id1.lt.id2) then; cmp=CMP_GT; else; cmp=CMP_LT; endif
           elseif(id1.lt.this%space_dim.and.id2.ge.this%space_dim) then
            cmp=CMP_LT
           else
            cmp=CMP_GT
           endif
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function HSpaceCompareSubspaces
!-------------------------------------------------------------------------
        function HSpaceRelateSubspaces(this,id1,id2,ierr) result(relation)
!Relates two subspaces from the hierarchical vector space.
         implicit none
         integer(INTD):: relation                    !out: relation of subspace 1 to subspace 2: {CMP_EQ,CMP_CN,CMP_IN,CMP_OV,CMP_NC}
         class(h_space_t), intent(in):: this         !in: hierarhical vector space
         integer(INTL), intent(in):: id1             !in: subspace 1 id
         integer(INTL), intent(in):: id2             !in: subspace 2 id
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,i
         integer(INTL):: sid,pid1,pid2
         type(vec_tree_iter_t):: vt_it

         errc=0; relation=CMP_NC
         if(id1.ge.0.and.id1.lt.this%num_subspaces.and.id2.ge.0.and.id2.lt.this%num_subspaces) then
          if(id1.eq.id2) then
           relation=CMP_EQ
          else
           errc=vt_it%init(this%subspaces)
           if(errc.eq.GFC_SUCCESS) then
            pid1=-1_INTL; pid2=-2_INTL
 !Try "subspace 1 is a parent of subspace 2":
            errc=vt_it%move_to(id2)
            if(errc.eq.GFC_SUCCESS) then
             do while(errc.eq.GFC_SUCCESS)
              errc=vt_it%move_to_parent(); if(errc.ne.GFC_SUCCESS) exit
              sid=vt_it%get_offset(errc); if(errc.ne.GFC_SUCCESS) exit
              if(sid.eq.id1) then; relation=CMP_CN; exit; endif
              if(pid2.lt.0_INTL) pid2=sid
             enddo
             if(errc.eq.GFC_NO_MOVE) errc=GFC_SUCCESS
 !Try "subspace 1 is a child of subspace 2":
             if(errc.eq.GFC_SUCCESS.and.relation.eq.CMP_NC) then
              errc=vt_it%move_to(id1)
              if(errc.eq.GFC_SUCCESS) then
               do while(errc.eq.GFC_SUCCESS)
                errc=vt_it%move_to_parent(); if(errc.ne.GFC_SUCCESS) exit
                sid=vt_it%get_offset(errc); if(errc.ne.GFC_SUCCESS) exit
                if(sid.eq.id2) then; relation=CMP_IN; exit; endif
                if(pid1.lt.0_INTL) pid1=sid
               enddo
               if(errc.eq.GFC_NO_MOVE) errc=GFC_SUCCESS
 !Try "subspace 1 and subspace 2 are siblings":
               if(errc.eq.GFC_SUCCESS.and.relation.eq.CMP_NC) then
                if(pid1.eq.pid2) relation=CMP_OV !siblings
               endif
              else
               errc=5
              endif
             endif
            else
             errc=4
            endif
            i=vt_it%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=3
           else
            errc=2
           endif
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function HSpaceRelateSubspaces
!---------------------------------------------------------------------
        function HSpaceCompareSubranges(this,id1,id2,ierr) result(cmp)
!Compares the basis subranges of two subspaces from the hierarhical vector space.
         implicit none
         integer(INTD):: cmp                         !out: comparison result: {CMP_EQ,CMP_LT,CMP_GT,CMP_OV,CMP_ER}
         class(h_space_t), intent(in):: this         !in: hierarchical vector space
         integer(INTL), intent(in):: id1             !in: subspace 1 id
         integer(INTL), intent(in):: id2             !in: subspace 2 id
         integer(INTD), intent(out), optional:: ierr !out: error code
         type(vec_tree_iter_t):: vt_it
         integer(INTD):: i,errc
         class(*), pointer:: up
         class(subspace_t), pointer:: sp1,sp2

         errc=0; cmp=CMP_ER
         if(id1.ge.0.and.id1.lt.this%num_subspaces.and.id2.ge.0.and.id2.lt.this%num_subspaces) then
          errc=vt_it%init(this%subspaces)
          if(errc.eq.GFC_SUCCESS) then
           up=>vt_it%element_value(id1,errc)
           if(errc.eq.GFC_SUCCESS) then
            sp1=>NULL(); select type(up); class is(subspace_t); sp1=>up; end select
            if(associated(sp1)) then
             up=>vt_it%element_value(id2,errc)
             if(errc.eq.GFC_SUCCESS) then
              sp2=>NULL(); select type(up); class is(subspace_t); sp2=>up; end select
              if(associated(sp2)) then
               cmp=sp1%compare_range(sp2); if(cmp.eq.CMP_ER) errc=4
              else
               errc=3
              endif
             endif
            else
             errc=2
            endif
           endif
           i=vt_it%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.GFC_SUCCESS) errc=i
          endif
         else
          errc=1
         endif
         up=>NULL(); sp1=>NULL(); sp2=>NULL()
         if(present(ierr)) ierr=errc
         return
        end function HSpaceCompareSubranges
!------------------------------------------------------------------------------
        function HSpaceGetCommonSubspace(this,id1,id2,ierr) result(subspace_id)
!Returns the minimal (registered) subspace containing two given subspaces.
         implicit none
         integer(INTL):: subspace_id                 !out: common (registered) subspace id
         class(h_space_t), intent(in):: this         !in: hierarchical vector space
         integer(INTL), intent(in):: id1             !in: input subspace id 1
         integer(INTL), intent(in):: id2             !in: input subspace id 2
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: i,cmp,errc
         type(vec_tree_iter_t):: vt_it

         if(this%is_set(errc)) then
          if(errc.eq.0) then
           if(id1.ge.0.and.id1.lt.this%num_subspaces.and.id2.ge.0.and.id2.lt.this%num_subspaces) then
            errc=vt_it%init(this%subspaces)
            if(errc.eq.GFC_SUCCESS) then
             subspace_id=id1; errc=vt_it%move_to(subspace_id)
             if(errc.eq.GFC_SUCCESS) then
              cmp=this%relate_subspaces(subspace_id,id2,errc)
              if(errc.eq.0) then
               do while(cmp.ne.CMP_CN)
                errc=vt_it%move_to_parent(); if(errc.ne.GFC_SUCCESS) exit
                subspace_id=vt_it%get_offset(errc); if(errc.ne.GFC_SUCCESS) exit
                cmp=this%relate_subspaces(subspace_id,id2,errc); if(errc.ne.GFC_SUCCESS) exit
               enddo
               if(errc.ne.GFC_SUCCESS) errc=8
              else
               errc=7
              endif
             else
              errc=6
             endif
             i=vt_it%release(); if(i.ne.GFC_SUCCESS.and.errc.eq.0) errc=5
            else
             errc=4
            endif
           else
            errc=3
           endif
          else
           errc=2
          endif
         else
          errc=1
         endif
         if(errc.ne.0) subspace_id=-1_INTL
         if(present(ierr)) ierr=errc
         return
        end function HSpaceGetCommonSubspace
!--------------------------------------------------
        subroutine HSpacePrintIt(this,ierr,dev_out)
!Prints the hierarchical vector space representation.
         implicit none
         class(h_space_t), intent(in):: this           !in: hierarchical vector space
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD), intent(in), optional:: dev_out !in: output device (defaults to screen)
         integer(INTD):: i,errc,devo
         class(*), pointer:: up
         class(subspace_t), pointer:: ssp
         type(vec_tree_iter_t):: vt_it

         devo=6; if(present(dev_out)) devo=dev_out
         errc=vt_it%init(this%subspaces)
         if(errc.eq.GFC_SUCCESS) then
          errc=GFC_IT_ACTIVE; ssp=>NULL()
          do while(errc.eq.GFC_IT_ACTIVE)
           up=>vt_it%get_value(errc); if(errc.ne.GFC_SUCCESS) exit
           select type(up); class is(subspace_t); ssp=>up; end select
           if(.not.associated(ssp)) then; errc=2; exit; endif
           call ssp%print_it(devo); ssp=>NULL()
           errc=vt_it%scanp(return_each=.TRUE.,skip_current=.TRUE.)
          enddo
          if(errc.eq.GFC_IT_DONE) then
           errc=vt_it%release()
          else
           i=vt_it%release()
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine HSpacePrintIt
!------------------------------------
        subroutine h_space_dtor(this)
         implicit none
         type(h_space_t):: this
         integer(INTD):: errc
         type(vec_tree_iter_t):: vt_it

         this%metric_p=>NULL()
         if(allocated(this%overlap)) deallocate(this%overlap)
         errc=vt_it%init(this%subspaces)
         if(errc.eq.GFC_SUCCESS) then
          errc=vt_it%delete_all()
          errc=vt_it%release()
         endif
         this%num_subspaces=0; this%space_dim=0
         return
        end subroutine h_space_dtor
!-----------------------------------------------------------------------------
        subroutine build_basis_hierarchy_abstract(basis,order,boundaries,ierr)
!Sorts and recursively aggregates bases into a hierarchical representatation.
         implicit none
         class(subspace_basis_t), intent(in):: basis             !in: vector space basis
         integer(INTL), intent(out), allocatable:: order(:)      !out: sorted order of basis vectors (old numbers)
         integer(INTL), intent(out), allocatable:: boundaries(:) !out: subspace boundaries
         integer(INTD), intent(out), optional:: ierr             !out: error code
         integer(INTD):: errc
         integer(INTL):: n,nsegs

         errc=0; n=basis%dimsn()
         if(n.gt.0) then
          !`Finish
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine build_basis_hierarchy_abstract
!-------------------------------------------------------------------------------
        subroutine build_basis_hierarchy_real_space(basis,order,boundaries,ierr)
!Sorts and recursively aggregates bases into a hierarchical representatation
!based on the support locality criteria.
         implicit none
         class(subspace_basis_t), intent(in):: basis             !in: vector space basis supported on real space
         integer(INTL), intent(out), allocatable:: order(:)      !out: sorted order of basis vectors (old numbers)
         integer(INTL), intent(out), allocatable:: boundaries(:) !out: subspace boundaries
         integer(INTD), intent(out), optional:: ierr             !out: error code
         integer(INTD):: errc
         integer(INTL):: n,nsegs

         errc=0; n=basis%dimsn()
         if(n.gt.0) then
          !`Finish
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine build_basis_hierarchy_real_space

       end module subspaces
