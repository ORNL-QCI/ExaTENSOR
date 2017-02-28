!Infrastructure for a recursive adaptive vector space decomposition.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/02/28

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
!   shows which subspaces overlap. For orthogonal basis sets, two
!   subspaces overlap if and only if they are releated in SAT.
!   For non-orthogonal subspaces, two subspaces overlap if and
!   only if their last level children subspaces (leaves) overlap.

       module subspaces
        use dil_basic
        use gfc_base
        use gfc_list
        use gfc_tree
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
!TYPES:
 !Real space vector:
        type, public:: real_vec_t
         integer(INTL), private:: num_dim=0      !number of dimensions (0 means empty)
         real(8), allocatable, public:: coord(:) !components of the real space vector
         contains
          procedure, private:: RealVecCtor
          generic, public:: real_vec_ctor=>RealVecCtor !real space vector ctor
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
 !Integer range = semi-interval(min:max]:
        type, public:: seg_int_t
         integer(INTL), private:: min_coord=0 !minimum coordinate (lower bound)
         integer(INTL), private:: max_coord=0 !maximum coordinate (upper bound)
         contains
          procedure, public:: set=>SegIntSet                !sets the integer range (ctor)
          procedure, public:: lower_bound=>SegIntLowerBound !returns the integer range lower bound
          procedure, public:: upper_bound=>SegIntUpperBound !returns the integer range upper bound
          procedure, public:: length=>SegIntLength          !returns the integer range length = (upper - lower)
          procedure, public:: overlap=>SegIntOverlap        !returns the overlap of two integer ranges
          procedure, public:: union=>SegIntUnion            !returns the minimal integer range containing two given integer ranges
          procedure, public:: split=>SegIntSplit            !splits the integer range
        end type seg_int_t
 !Real space rectangular hypercube (orthotope):
        type, public:: orthotope_t
         integer(INTL), private:: num_dim=0                    !number of dimensions
         type(range1d_t), allocatable, private:: extent(:)     !extent of each dimension (min,max)
         contains
          procedure, private:: OrthotopeCtor
          generic, public:: orthotope_ctor=>OrthotopeCtor      !orthotope ctor
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
 !Symmetry:
        type, abstract, public:: space_symmetry_t
         contains
          procedure(space_symm_combine_i), deferred, public:: combine !combines symmetry group irreps
        end type space_symmetry_t
 !Abstract basis function (basis function support only):
        type, public:: basis_func_supp_t
         integer(INTD), private:: supp_dim=0      !dimensionality of the real space support on which the basis function resides
         type(real_vec_t), private:: center       !center of the effective function support in the real space
         type(orthotope_t), private:: supp_box    !supporting orthotope (multidimensional real space support)
        end type basis_func_supp_t
 !Gaussian basis function:
        type, extends(basis_func_supp_t), public:: basis_func_gauss_t
         integer(INTD), private:: num_prims=0        !number of contracted primitives
         integer(INTD), private:: orb_moment=-1      !orbital momentum (0,1,2,3,...)
         real(8), allocatable, private:: exponent(:) !primitive exponents
         complex(8), allocatable, private:: coef(:)  !primitive contraction coefficients
         contains
          procedure, public:: set=>BasisFuncGaussSet !sets up the basis function
          final:: BasisFuncGaussDestroy              !destructs the basis function
        end type basis_func_gauss_t
 !Typed basis function:
        type, public:: basis_func_t
         integer(INTD), private:: basis_kind=BASIS_NONE                    !specific basis kind
         class(basis_func_supp_t), pointer, private:: basis_func_p=>NULL() !pointer to a basis function of this kind
         class(space_symmetry_t), pointer, private:: symm=>NULL()          !symmetry of the basis function (if any)
         contains
          procedure, public:: set=>BasisFuncSet !sets up the basis function
        end type basis_func_t
 !Subspace basis:
        type, public:: subspace_basis_t
         integer(INTL), private:: space_dim=0                     !number of basis functions
         integer(INTD), private:: supp_dim=0                      !dimensionality of the real space support on which the basis functions reside
         type(real_vec_t), private:: center                       !center of the effective subspace basis support in real space
         type(orthotope_t), private:: supp_box                    !effective subspace basis support in real space (multidimensional orthotope)
         class(space_symmetry_t), pointer, private:: symm=>NULL() !symmetry of the subspace basis (if any)
         type(basis_func_t), allocatable, private:: basis_func(:) !basis functions [1..space_dim]
         contains
          procedure, public:: create=>SubspaceBasisCreate               !creates an empty subspace
          procedure, public:: dimsn=>SubspaceBasisDimsn                 !returns the dimension of the subspace
          procedure, public:: supp_dimsn=>SubspaceBasisSuppDimsn        !returns the support space dimension
          procedure, public:: set_basis_func=>SubspaceBasisSetBasisFunc !sets a specific basis function
          procedure, public:: get_basis_func=>SubspaceBasisGetBasisFunc !returns a pointer to a specific basis function
          procedure, public:: finalize=>SubspaceBasisFinalize           !finalizes the subspace basis (sets up the support and symmetry)
          procedure, public:: get_symmetry=>SubspaceBasisGetSymmetry    !returns a pointer to the subspace basis symmetry object
          procedure, public:: get_center=>SubspaceBasisGetCenter        !returns a pointer to the center of the subspace basis in the real space
          procedure, public:: get_support=>SubspaceBasisGetSupport      !returns a pointer to the supporting orthotope of the subspace basis
          final:: SubspaceBasisDestroy                                  !destructs the subspace basis
        end type subspace_basis_t
 !Subspace:
        type, public:: subspace_t
         integer(INTL), private:: subspace_id=-1   !unique subspace ID (registered ID): must be non-negative, -1 means undefined
         integer(INTD), private:: supp_dim=0       !dimensionality of the real space support on which the basis functions reside
         integer(INTL), private:: max_resolution=0 !max resolution level (max subspace dimension): 0 means undefined
         type(list_bi_t), private:: bases          !basis sets (subspace_basis_t) for each registered resolution level
         contains
          procedure, public:: init=>SubspaceInit                           !initializes the subspace (id and support dimension only)
          procedure, public:: get_id=>SubspaceGetId                        !returns the subspace id
          procedure, public:: get_supp_dim=>SubspaceGetSuppDim             !returns the dimensionality of the subspace support
          procedure, public:: get_max_resolution=>SubspaceGetMaxResolution !returns the max resolution (dimension) of the subspace
          procedure, public:: register_basis=>SubspaceRegisterBasis        !registers a specific basis of the subspace
          procedure, public:: resolve=>SubspaceResolve                     !resolves the subspace with a specific basis (based on some condition)
          procedure, public:: relate=>SubspaceRelate                       !relates the subspaces to another subspace
          final:: SubspaceDestroy                                          !destroys the subspace
        end type subspace_t
 !Hierarchical composite index:
        type, public:: h_index_t
         integer(INTL), public:: subspace_id=-1 !subspace ID (registered ID): must be non-negative, -1 means undefined
         integer(INTL), public:: resolution=0   !subspace resolution level 1<=resolution<=max_resolution: 0 means undefined
         integer(INTL), public:: component=0    !subspace component number at the given level of resolution: [1..resolution], 0 means undefined
        end type h_index_t
 !Hierarchical space representation:
        type, public:: h_space_t
         integer(INTL), private:: space_dim=0                 !dimension of the vector space
         integer(INTL), private:: num_subspaces=0             !number of subspaces defined in the vector space
         type(subspace_t), allocatable, private:: subspace(:) !hierarchical subspaces defined in the vector space: [0..num_subspaces-1]
         type(tree_t), private:: aggr_tree                    !subspace aggregation tree (SAT): Hierarchical representation
         complex(8), pointer, private:: metric(:,:)=>NULL()   !pointer to the original metric tensor: g12=<bf1|bf2>
         real(8), allocatable, private:: overlap(:,:)         !subspace overlap matrix (extent of overlap between all subspaces)
         contains
          procedure, public:: construct=>HSpaceConstruct      !constructs a hierarchical representation of a vector space
          procedure, public:: is_set=>HSpaceIsSet             !return TRUE if the hierarchical vector space is set
          procedure, public:: get_subspace=>HSpaceGetSubspace !returns a pointer to the required subspace of the hierarchical vector space
          final:: HSpaceDestruct                              !destructs the hierarchical representation of a vector space
        end type h_space_t
!INTERFACES:
 !space_symmetry_t:
        abstract interface
  !Deferred: .combine:
         subroutine space_symm_combine_i(this,symm,ierr)
          import:: space_symmetry_t,INTD
          implicit none
          class(space_symmetry_t), intent(inout):: this       !inout: symmetry 1 (updated)
          class(space_symmetry_t), intent(in):: symm          !in: symmetry 2
          integer(INTD), intent(out), optional:: ierr         !out: error code
         end subroutine space_symm_combine_i
        end interface
!VISIBILITY:
 !Non-member:
        public build_basis_hierarchy_abstract   !establishes a hierarchy for an abstract basis with possible symmetries
        public build_basis_hierarchy_real_space !establishes a hierarchy for a real space supported basis with possible symmetries
 !real_vec_t:
        private RealVecCreate
        private RealVecDimsn
        private RealVecNorm2
        private RealVecScale
        private RealVecAdd
        private RealVecAverage
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
 !orthotope_t:
        private OrthotopeCreate
        private OrthotopeDimsn
        private OrthotopeSetExtent
        private OrthotopeGetExtent
        private OrthotopeLowerBound
        private OrthotopeUpperBound
        private OrthotopeLength
        private OrthotopeVolume
        private OrthotopeOverlap
        private OrthotopeUnion
 !basis_func_gauss_t:
        private BasisFuncGaussSet
 !basis_func_t:
        private BasisFuncSet
 !subspace_basis_t:
        private SubspaceBasisCreate
        private SubspaceBasisDimsn
        private SubspaceBasisSuppDimsn
        private SubspaceBasisSetBasisFunc
        private SubspaceBasisGetBasisFunc
        private SubspaceBasisFinalize
        private SubspaceBasisGetSymmetry
        private SubspaceBasisGetCenter
        private SubspaceBasisGetSupport
 !subspace_t:
        private SubspaceInit
        private SubspaceGetId
        private SubspaceGetSuppDim
        private SubspaceGetMaxResolution
        private SubspaceRegisterBasis
        private SubspaceResolve
        private SubspaceRelate
 !h_space_t:
        private HSpaceConstruct
        private HSpaceIsSet
        private HSpaceGetSubspace

       contains
!IMPLEMENTATION:
![real_vec_t]====================================
        subroutine RealVecCreate(this,dimsn,ierr)
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
        end subroutine RealVecCreate
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
           if(res_vec%num_dim.ne.n) call res_vec%create(n,errc)
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
!--------------------------------------
        subroutine RealVecDestroy(this)
         implicit none
         type(real_vec_t):: this

         if(allocated(this%coord)) deallocate(this%coord)
         this%num_dim=0
         return
        end subroutine RealVecDestroy
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
         real(8):: rl,incr,lb,ub
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
               if(align(left)-align(left2).lt.align(right)-align(j)) then
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
           if(nbnd.gt.0) then !aligned splitting
            if(num_segs.le.nbnd) then
             ub=this%upper_bound()
             nchnk=nbnd+1 !number of alignment chunks
             do while(nchnk.gt.num_segs)
 !Find the smallest alignment chunk:
              left=0; left2=0; right=0
              j=0; k=0; l=0; lb=0d0; next=.TRUE.
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
               if(align(left)-align(left2).lt.align(right)-align(j)) then
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
![orthotope_t]=====================================
        subroutine OrthotopeCreate(this,dimsn,ierr)
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
        end subroutine OrthotopeCreate
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
           if(overlap%num_dim.ne.n) call overlap%create(n,errc)
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
           if(union%num_dim.ne.n) call union%create(n,errc)
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
!----------------------------------------
        subroutine OrthotopeDestroy(this)
         implicit none
         type(orthotope_t):: this

         if(allocated(this%extent)) deallocate(this%extent)
         this%num_dim=0
         return
        end subroutine OrthotopeDestroy
![basis_func_gauss_t]=====================================================
        subroutine BasisFuncGaussSet(this,orb_moment,exponents,coefs,ierr)
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
        end subroutine BasisFuncGaussSet
!---------------------------------------------
        subroutine BasisFuncGaussDestroy(this)
         implicit none
         type(basis_func_gauss_t):: this

         if(allocated(this%coef)) deallocate(this%coef)
         if(allocated(this%exponent)) deallocate(this%exponent)
         this%num_prims=0; this%orb_moment=-1
         return
        end subroutine BasisFuncGaussDestroy
![basis_func_t]======================================================
        subroutine BasisFuncSet(this,basis_kind,ierr,basis_func,symm)
!Sets up a basis function of a given kind. If the basis function is already set,
!it will be redefined (no non-trivial destruction is assumed). If no basis
!function is passed here, an abstract basis function of <basis_kind> is assumed.
         implicit none
         class(basis_func_t), intent(out):: this                             !out: basis function
         integer(INTD), intent(in):: basis_kind                              !in: basis kind
         integer(INTD), intent(out), optional:: ierr                         !out: error code
         class(basis_func_supp_t), intent(in), target, optional:: basis_func !in: specific basis function
         class(space_symmetry_t), intent(in), target, optional:: symm        !in: basis function symmetry
         integer(INTD):: errc

         errc=0
         this%basis_func_p=>NULL(); this%symm=>NULL()
         if(basis_kind.ne.BASIS_NONE) then
          this%basis_kind=basis_kind
          if(present(basis_func)) this%basis_func_p=>basis_func
          if(present(symm)) this%symm=>symm
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine BasisFuncSet
![subspace_basis_t]====================================
        subroutine SubspaceBasisCreate(this,dimsn,ierr)
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
        end subroutine SubspaceBasisCreate
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
         class(basis_func_supp_t), intent(in), target, optional:: basis_func !in: specific basis function
         class(space_symmetry_t), intent(in), target, optional:: symm        !in: basis function symmetry
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
              call this%basis_func(func_num)%set(basis_kind,errc,basis_func,symm)
             else
              call this%basis_func(func_num)%set(basis_kind,errc,basis_func)
             endif
             if(errc.ne.0) errc=5
            else
             errc=4
            endif
           else
            if(present(symm)) then
             call this%basis_func(func_num)%set(basis_kind,errc,symm=symm)
            else
             call this%basis_func(func_num)%set(basis_kind,errc)
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
!real space supported basis sets.  Also computes the overall symmetry group irrep.
!If not all basis functions were set, the <num_undef> argument will return
!the number of unset basis functions and an error code.
         implicit none
         class(subspace_basis_t), intent(inout):: this    !inout: subspace basis
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTL), intent(out), optional:: num_undef !out: number of undefined basis functions
         integer(INTD):: errc
         integer(INTL):: i,n,nun
         class(basis_func_supp_t), pointer:: bas_func
         class(space_symmetry_t), pointer:: symm
         logical:: initb,inits

         errc=0; nun=-1; n=this%dimsn()
         if(n.gt.0) then
          nun=0; initb=.FALSE.; inits=.FALSE.
          bloop: do i=1,n
           if(this%basis_func(i)%basis_kind.ne.BASIS_NONE) then
            bas_func=>this%basis_func(i)%basis_func_p
            symm=>this%basis_func(i)%symm
            if(associated(bas_func)) then
             if(.not.initb) then
              this%supp_dim=bas_func%supp_dim
              if(this%supp_dim.gt.0) then
               call this%center%create(int(this%supp_dim,INTL),errc)
               if(errc.ne.0) then; errc=8; exit bloop; endif
               this%center=bas_func%center
               call this%supp_box%create(int(this%supp_dim,INTL),errc)
               if(errc.ne.0) then; errc=7; exit bloop; endif
               this%supp_box=bas_func%supp_box
              endif
              initb=.TRUE.
             else
              if(bas_func%supp_dim.eq.this%supp_dim) then
               call this%center%add(bas_func%center,errc)
               if(errc.ne.0) then; errc=6; exit; exit bloop; endif
               call this%supp_box%union(bas_func%supp_box,errc)
               if(errc.ne.0) then; errc=5; exit; exit bloop; endif
              else
               errc=4; exit bloop !support dimension mismatch
              endif
             endif
            endif
            if(associated(symm)) then
             if(.not.inits) then
              this%symm=>symm; inits=.TRUE.
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
           if(nun.lt.n) call this%center%scale(1d0/real(n-nun,8),errc)
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
         class(space_symmetry_t), pointer:: symm     !out: subspace basis symmetry
         class(subspace_basis_t), intent(in):: this  !in: subspace basis
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(this%space_dim.gt.0) then
          symm=>this%symm
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
!--------------------------------------------
        subroutine SubspaceBasisDestroy(this)
         implicit none
         type(subspace_basis_t):: this

         if(allocated(this%basis_func)) deallocate(this%basis_func)
         this%symm=>NULL()
         this%space_dim=0; this%supp_dim=0
         return
        end subroutine SubspaceBasisDestroy
![subspace_t]================================
        subroutine SubspaceInit(this,id,ierr)
!Initializes a subspace. If the subspace is defined on input,
!it will automatically be destructed prior to re-initialization.
         implicit none
         class(subspace_t), intent(out):: this       !out: empty subspace
         integer(INTL), intent(in):: id              !in: subspace id (must be non-negative)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(id.ge.0) then
          this%subspace_id=id
          this%supp_dim=0       !will be inferred from the first registered basis supported in real space
          this%max_resolution=0 !will be inferred from the first registered basis
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine SubspaceInit
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
!--------------------------------------------------------
        subroutine SubspaceRegisterBasis(this,basis,ierr)
!Registers a new subspace basis. The new subspace basis must have
!the same support dimension as any previously registered basis.
         implicit none
         class(subspace_t), intent(inout):: this     !inout: subspace
         class(subspace_basis_t), intent(in):: basis !in: new subspace basis
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,ier,sd
         integer(INTL):: n
         type(list_iter_t):: basis_it

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
!the very first basis will be returned.
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
!----------------------------------------------------------------
        function SubspaceRelate(this,another,h_space) result(cmp)
!Relates the given subspace with another subspace in the given subspace hierarchy.
         implicit none
         integer(INTD):: cmp                     !out: relation: {GFC_CMP_EQ,GFC_CMP_LT,GFC_CMP_GT,GFC_CMP_CHILD,GFC_CMP_PARENT,GFC_CMP_ERR}
         class(subspace_t), intent(in):: this    !in: subspace 1
         class(subspace_t), intent(in):: another !in: subspace 2
         class(h_space_t), intent(in):: h_space  !in: hierarchical subspace representation

         cmp=GFC_CMP_EQ
         !`Finish
         return
        end function SubspaceRelate
!---------------------------------------
        subroutine SubspaceDestroy(this)
         implicit none
         type(subspace_t):: this
         type(list_iter_t):: basis_it
         integer(INTD):: errc

         errc=basis_it%init(this%bases)
         if(errc.eq.GFC_SUCCESS) then
          errc=basis_it%reset(); if(errc.eq.GFC_SUCCESS) errc=basis_it%delete_all()
          errc=basis_it%release()
         endif
         this%max_resolution=0
         this%supp_dim=0
         this%subspace_id=-1
         return
        end subroutine SubspaceDestroy
![h_space_t]===================================================
        subroutine HSpaceConstruct(this,full_basis,ierr,metric)
!Constructs a hierarchical vector space representation with a subspace aggregation tree.
!The original basis functions will be hierarchically aggregated into larger subspaces,
!up to the full space. Each subspace will have a unique id.
         implicit none
         class(h_space_t), intent(out):: this                   !out: hierarchical representation of the vector space
         class(subspace_basis_t), intent(in):: full_basis       !in: full basis of the vector space
         integer(INTD), intent(out), optional:: ierr            !out: error code
         complex(8), intent(in), optional, target:: metric(:,:) !in: metric tensor: g12=<bf1|bf2>: Hermitian matrix
         integer(INTD):: errc
         integer(INTL):: n
         type(tree_iter_t):: sat_it

         errc=0
!Construct the subspace aggregation tree (SAT) by recursion:
         n=full_basis%dimsn()
         if(n.gt.0) then
          errc=sat_it%init(this%aggr_tree)
          if(errc.eq.GFC_SUCCESS) then
 !Add the root (full space):
           !`Finish
 !Recursively split the full space into subspaces:
           
           errc=sat_it%release()
!Construct the overlap matrix between all subspaces:
           if(present(metric)) then
            this%metric=>metric
            !`Write
           endif
          else
           errc=2
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine HSpaceConstruct
!--------------------------------------------------
        function HSpaceIsSet(this,ierr) result(res)
!Returns TRUE if the hierarchical vector space is set.
         implicit none
         logical:: res                               !out: result
         class(h_space_t), intent(in):: this         !in: hierarchical vector space
         integer(INTD), intent(out), optional:: ierr !out: eror code
         integer(INTD):: errc

         errc=0; res=(this%space_dim.gt.0)
         if(present(ierr)) ierr=errc
         return
        end function HSpaceIsSet
!---------------------------------------------------------------------------
        function HSpaceGetSubspace(this,subspace_id,ierr) result(subspace_p)
!Returns a pointer to the requested subspace from the hierarchical vector space.
         implicit none
         class(subspace_t), pointer:: subspace_p      !out: pointer to the requested subspace
         class(h_space_t), intent(in), target:: this  !in: hierarchical vector space
         integer(INTL), intent(in):: subspace_id      !in: requested subspace id
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         errc=0
         if(this%is_set()) then
          if(subspace_id.ge.0.and.subspace_id.lt.this%num_subspaces) then
           subspace_p=>this%subspace(subspace_id)
          else
           errc=2
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function HSpaceGetSubspace
!--------------------------------------
        subroutine HSpaceDestruct(this)
         implicit none
         type(h_space_t):: this
         integer(INTD):: errc
         type(tree_iter_t):: tree_it

         this%metric=>NULL()
         if(allocated(this%overlap)) deallocate(this%overlap)
         if(allocated(this%subspace)) deallocate(this%subspace)
         errc=tree_it%init(this%aggr_tree)
         if(errc.eq.GFC_SUCCESS) then
          errc=tree_it%reset(); if(errc.eq.GFC_SUCCESS) errc=tree_it%delete_subtree()
          errc=tree_it%release()
         endif
         this%space_dim=0; this%num_subspaces=0
         return
        end subroutine HSpaceDestruct
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
         class(subspace_basis_t), intent(in):: basis             !in: vector space basis supported in real space
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
