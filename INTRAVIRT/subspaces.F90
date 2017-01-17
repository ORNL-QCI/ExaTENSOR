!Infrastructure for a recursive adaptive vector space decomposition.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/01/17

!Copyright (C) 2014-2016 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2016 Oak Ridge National Laboratory (UT-Battelle)

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
        integer(INTD), parameter, public:: BASIS_NONE=0    !No basis set
        integer(INTD), parameter, public:: BASIS_GAUSS=1   !Gaussian basis set
        integer(INTD), parameter, public:: BASIS_WAVELET=2 !Wavelet basis set
!TYPES:
 !Real space vector:
        type, public:: real_vec_t
         integer(INTL), private:: num_dim=0      !number of dimensions (0 means empty)
         real(8), allocatable, public:: coord(:) !components of the real space vector
         contains
          procedure, public:: create=>RealVecCreate !creates an empty real space vector
          procedure, public:: dimsn=>RealVecDimsn   !returns dimension of the vector
          procedure, public:: norm2=>RealVecNorm2   !returns the 2-norm of the vector
          final:: RealVecDestroy                    !destroys the vector
        end type real_vec_t
 !Real space 1d extent (segment[min:max]):
        type, public:: extent1d_t
         real(8), private:: min_coord=0d0 !minimum coordinate (lower bound)
         real(8), private:: max_coord=0d0 !maximum coordinate (upper bound)
         contains
          procedure, public:: set=>Extent1dSet                !sets the extent (ctor)
          procedure, public:: lower_bound=>Extent1dLowerBound !returns the extent lower bound
          procedure, public:: upper_bound=>Extent1dUpperBound !returns the extent upper bound
          procedure, public:: length=>Extent1dLength          !returns the extent length
          procedure, public:: overlap=>Extent1dOverlap        !returns the overlap of two extents
          procedure, public:: union=>Extent1dUnion            !returns the union of two extents
        end type extent1d_t
 !Real space rectangular hypercube (orthotope):
        type, public:: orthotope_t
         integer(INTL), private:: num_dim=0                   !number of dimensions
         type(extent1d_t), allocatable, private:: extent(:)   !extent of each dimension (min,max)
         contains
          procedure, public:: create=>OrthotopeCreate          !creates an empty orthotope
          procedure, public:: dimsn=>OrthotopeDimsn            !returns the real space dimension orthotope resides in
          procedure, public:: set_extent=>OrthotopeSetExtent   !sets the extent along a specific dimension
          procedure, public:: get_extent=>OrthotopeGetExtent   !returns the extent of a specific dimension
          procedure, public:: lower_bound=>OrthotopeLowerBound !returns the lower bound of a specific extent
          procedure, public:: upper_bound=>OrthotopeUpperBound !returns the upper bound of a specific extent
          procedure, public:: length=>OrthotopeLength          !returns the length along a specific extent
          procedure, public:: volume=>OrthotopeVolume          !returns the volume of the orthotope
          procedure, public:: overlap=>OrthotopeOverlap        !returns the overlap of two orthotopes
          procedure, public:: union=>OrthotopeUnion            !returns the minimal orthotope containing two orthotopes
          final:: OrthotopeDestroy                             !destroys the orthotope
        end type orthotope_t
 !Symmetry:
        type, abstract, public:: space_symmetry_t
         contains
          procedure(space_symm_combine_i), deferred, public:: get_combined !combines group irreps (symmetry group multiplication)
        end type space_symmetry_t
 !Abstract basis function:
        type, abstract, public:: basis_func_abs_t
         integer(INTD), private:: supp_dim=0      !dimensionality of the real space support on which the basis function resides
         type(real_vec_t), private:: center       !center of the effective function support in the real space
         type(orthotope_t), private:: supp_box    !supporting orthotope (multidimensional real space support)
        end type basis_func_abs_t
 !Gaussian basis function:
        type, extends(basis_func_abs_t), public:: basis_func_gauss_t
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
         integer(INTD), private:: basis_kind=BASIS_NONE                   !specific basis kind
         class(space_symmetry_t), pointer, private:: symm=>NULL()         !symmetry of the basis function (if any)
         class(basis_func_abs_t), pointer, private:: basis_func_p=>NULL() !pointer to a basis function of this kind
         contains
          procedure, public:: set=>BasisFuncSet !sets up the basis function
        end type basis_func_t
 !Subspace basis:
        type, public:: subspace_basis_t
         integer(INTL), private:: space_dim=0                     !number of basis functions
         integer(INTD), private:: supp_dim=0                      !dimensionality of the real space support on which the basis functions reside
         class(space_symmetry_t), pointer, private:: symm=>NULL() !symmetry of the subspace basis (if any)
         type(real_vec_t), private:: center                       !center of the effective subspace basis support in real space
         type(orthotope_t), private:: supp_box                    !effective subspace basis support in real space (multidimensional orthotope)
         type(basis_func_t), allocatable, private:: basis_func(:) !basis functions [1..space_dim]
         contains
          procedure, public:: create=>SubspaceBasisCreate               !creates an empty subspace
          procedure, public:: dimsn=>SubspaceBasisDimsn                 !returns the dimension of the subspace
          procedure, public:: set_basis_func=>SubspaceBasisSetBasisFunc !sets a specific basis function
          procedure, public:: get_basis_func=>SubspaceBasisGetBasisFunc !returns a pointer to a specific basis function
          procedure, public:: finalize=>SubspaceBasisFinalize           !finalizes the subspace basis (sets up the support and symmetry)
          final:: SubspaceBasisDestroy                                  !destructs the subspace basis
        end type subspace_basis_t
 !Subspace:
        type, public:: subspace_t
         integer(INTL), private:: subspace_id=-1   !unique subspace ID (registered ID): must be non-negative, -1 means undefined
         integer(INTD), private:: supp_dim=0       !dimensionality of the real space support on which the basis functions reside
         integer(INTL), private:: max_resolution=0 !max resolution level (max subspace dimension): 0 means undefined
         class(subspace_basis_t), pointer, private:: curr_basis=>NULL() !currently set subspace basis
         type(list_bi_t), private:: bases          !basis sets (subspace_basis_t) for each registered resolution level
         contains
          procedure, public:: init=>SubspaceInit                           !initializes the subspace (id and support dimension only)
          procedure, public:: get_id=>SubspaceGetId                        !returns the subspace id
          procedure, public:: get_supp_dim=>SubspaceGetSuppDim             !returns the dimensionality of the subspace support
          procedure, public:: get_max_resolution=>SubspaceGetMaxResolution !returns the max resolution (dimension) of the subspace
          procedure, public:: get_symmetry=>SubspaceGetSymmetry            !returns the current symmetry of the subspace
          procedure, public:: get_center=>SubspaceGetCenter                !returns the current center of the subspace in the real space
          procedure, public:: get_support=>SubspaceGetSupport              !returns the current supporting orthotope
          procedure, public:: register_basis=>SubspaceRegisterBasis        !registers a specific basis of the subspace
          procedure, public:: set_basis=>SubspaceSetBasis                  !sets a specific subspace basis as current (based on some condition)
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
         type(subspace_t), allocatable, private:: subspace(:) !hierarchical subspaces defined in the vector space
         type(tree_t), private:: aggr_tree                    !subspace aggregation tree (SAT): Hierarchical representation
         contains
          procedure, public:: construct=>HSpaceConstruct      !constructs a hierarchical representation of a vector space
          final:: HSpaceDestruct                              !destructs the hierarchical representation of a vector space
        end type h_space_t
!INTERFACES:
 !space_symmetry_t:
        abstract interface
  !Deferred: .get_combined:
         subroutine space_symm_combine_i(this,other_symm,ierr)
          implicit none
          import:: space_symmetry_t,INTD
          class(space_symmetry_t), intent(inout):: this
          class(space_symmetry_t), intent(in):: other_symm
          integer(INTD), intent(out), optional:: ierr
         end subroutine space_symm_combine_i
        end interface
!VISIBILITY:
 !real_vec_t:
        private RealVecCreate
        private RealVecDimsn
        private RealVecNorm2
 !extent1d_t:
        private Extent1dSet
        private Extent1dLowerBound
        private Extent1dUpperBound
        private Extent1dLength
        private Extent1dOverlap
        private Extent1dUnion
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
        private SubspaceBasisSetBasisFunc
        private SubspaceBasisGetBasisFunc
        private SubspaceBasisFinalize
 !subspace_t:
        private SubspaceInit
        private SubspaceGetId
        private SubspaceGetSuppDim
        private SubspaceGetMaxResolution
        private SubspaceGetSymmetry
        private SubspaceGetCenter
        private SubspaceGetSupport
        private SubspaceRegisterBasis
        private SubspaceGetBasis
        private SubspaceUpdateSupport
 !h_space_t:
        private HSpaceConstruct

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
         implicit none
         integer(INTL):: dimsn                       !out: vector dimension
         class(real_vec_t), intent(in):: this        !in: real space vector

         dimsn=this%num_dim
         return
        end function RealVecDimsn
!-----------------------------------------------------
        function RealVecNorm2(this,ierr) result(norm2)
         implicit none
         real(8):: norm2                             !out: vector 2-norm
         class(real_vec_t), intent(in):: this        !in: real space vector
         integer(INTD), intent(out), optional:: ierr !out: error code
         !------------------------------------------
         integer(INTL), parameter:: LARGE_VECTOR=(2_INTL)**20
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
         if(present(ierr)) ierr=errc
         return
        end function RealVecNorm2
!--------------------------------------
        subroutine RealVecDestroy(this)
         implicit none
         type(real_vec_t):: this

         if(allocated(this%coord)) deallocate(this%coord)
         this%num_dim=0
         return
        end subroutine RealVecDestroy
![extent1d_t]========================================
        subroutine Extent1dSet(this,lower,upper,ierr)
         implicit none
         class(extent1d_t), intent(inout):: this     !inout: extent
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
        end subroutine Extent1dSet
!---------------------------------------------------
        function Extent1dLowerBound(this) result(lb)
         implicit none
         real(8):: lb
         class(extent1d_t), intent(in):: this

         lb=this%min_coord
         return
        end function Extent1dLowerBound
!---------------------------------------------------
        function Extent1dUpperBound(this) result(ub)
         implicit none
         real(8):: ub
         class(extent1d_t), intent(in):: this

         ub=this%max_coord
         return
        end function Extent1dUpperBound
!--------------------------------------------------
        function Extent1dLength(this) result(length)
         implicit none
         real(8):: length                     !out: extent length
         class(extent1d_t), intent(in):: this !in: extent

         length=this%max_coord-this%min_coord
         return
        end function Extent1dLength
!---------------------------------------------------------------
        function Extent1dOverlap(this,extent) result(res_extent)
         implicit none
         type(extent1d_t):: res_extent          !out: resulting extent (overlap)
         class(extent1d_t), intent(in):: this   !in: extent 1
         class(extent1d_t), intent(in):: extent !in: extent 2

         if(this%max_coord.le.extent%min_coord.or.this%min_coord.ge.extent%max_coord) then
          call res_extent%set(0d0,0d0) !no overlap
         else
          call res_extent%set(max(this%min_coord,extent%min_coord),min(this%max_coord,extent%max_coord))
         endif
         return
        end function Extent1dOverlap
!-------------------------------------------------------------
        function Extent1dUnion(this,extent) result(res_extent)
         implicit none
         type(extent1d_t):: res_extent          !out: resulting extent (union)
         class(extent1d_t), intent(in):: this   !in: extent 1
         class(extent1d_t), intent(in):: extent !in: extent 2

         call res_extent%set(min(this%min_coord,extent%min_coord),max(this%max_coord,extent%max_coord))
         return
        end function Extent1dUnion
![orthotope_t]=====================================
        subroutine OrthotopeCreate(this,dimsn,ierr)
!Creates an empty orthotope. If the orthotope is defined in input,
!it will be automatically destructed prior to the re-initialization.
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
         implicit none
         integer(INTL):: dimsn                       !out: real space dimension (0 means empty orthotope)
         class(orthotope_t), intent(in):: this       !in: orthotope

         dimsn=this%num_dim
         return
        end function OrthotopeDimsn
!------------------------------------------------------------
        subroutine OrthotopeSetExtent(this,dimsn,extent,ierr)
         implicit none
         class(orthotope_t), intent(inout):: this    !inout: orthotope
         integer(INTL), intent(in):: dimsn           !in: specific dimension to set extent over
         type(extent1d_t), intent(in):: extent       !in: extent
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
         implicit none
         type(extent1d_t):: extent                   !out: extent of the specific dimension
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
!---------------------------------------------------------------------
        function OrthotopeOverlap(this,orthotope,ierr) result(overlap)
         implicit none
         type(orthotope_t):: overlap                 !out: overlap (orthotope)
         class(orthotope_t), intent(in):: this       !in: orthotope 1
         class(orthotope_t), intent(in):: orthotope  !in: orthotope 2
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer(INTL):: i,n

         errc=0; n=this%num_dim
         if(n.gt.0.and.n.eq.orthotope%num_dim) then
          call overlap%create(n,errc)
          if(errc.eq.0) then
           do i=1,n
            overlap%extent(i)=this%extent(i)%overlap(orthotope%extent(i))
           enddo
          else
           errc=2
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function OrthotopeOverlap
!-----------------------------------------------------------------
        function OrthotopeUnion(this,orthotope,ierr) result(union)
         implicit none
         type(orthotope_t):: union                   !out: union (orthotope)
         class(orthotope_t), intent(in):: this       !in: orthotope 1
         class(orthotope_t), intent(in):: orthotope  !in: orthotope 2
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer(INTL):: i,n

         errc=0; n=this%num_dim
         if(n.gt.0.and.n.eq.orthotope%num_dim) then
          call union%create(n,errc)
          if(errc.eq.0) then
           do i=1,n
            union%extent(i)=this%extent(i)%union(orthotope%extent(i))
           enddo
          else
           errc=2
          endif
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function OrthotopeUnion
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
!Sets the Gaussian basis function. If it is already defined, it will be reset.
         implicit none
         class(basis_func_gauss_t), intent(out):: this !out: basis function
         integer(INTD), intent(in):: orb_moment        !in: orbital momentum
         real(8), intent(in):: exponents(1:)           !in: exponents
         complex(8), intent(in):: coefs(1:)            !in: contraction coefficients
         integer(INTD), intent(out), optional:: ierr   !out: error code
         integer(INTD):: np,errc

         errc=0
         if(orb_moment.ge.0) then
          np=size(exponents)
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
        subroutine BasisFuncSet(this,basis_kind,basis_func,ierr,symm)
         implicit none
         class(basis_func_t), intent(out):: this                  !out: basis function
         integer(INTD), intent(in):: basis_kind                   !in: basis kind
         class(basis_func_abs_t), intent(in), target:: basis_func !in: specific basis function
         integer(INTD), intent(out), optional:: ierr              !out: error code
         type(space_symmetry_t), intent(in), optional:: symm      !in: basis function symmetry
         integer(INTD):: errc

         errc=0
         if(basis_kind.ne.BASIS_NONE) then
          this%basis_kind=basis_kind
          this%basis_func_p=>basis_func
          if(present(symm)) this%symm=symm
         else
          this%basis_func=>NULL()
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine BasisFuncSet
![subspace_basis_t]====================================
        subroutine SubspaceBasisCreate(this,dimsn,ierr)
!Creates an empty subspace basis. If the subspace basis is defined on input,
!it will be automatically destructed prior to the re-initialization.
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
         implicit none
         integer(INTL):: dimsn                      !out: subspace basis dimension
         class(subspace_basis_t), intent(in):: this !in: subspace basis

         dimsn=this%space_dim
         return
        end function SubspaceBasisDimsn
!------------------------------------------------------------------------------------------
        subroutine SubspaceBasisSetBasisFunc(this,func_num,basis_kind,basis_func,ierr,symm)
         implicit none
         class(subspace_basis_t), intent(inout):: this            !inout: subspace basis
         integer(INTL), intent(in):: func_num                     !in: basis function number in the subspace basis
         integer(INTD), intent(in):: basis_kind                   !in: basis kind
         class(basis_func_abs_t), intent(in), target:: basis_func !in: specific basis function
         integer(INTD), intent(out), optional:: ierr              !out: error code
         type(space_symmetry_t), intent(in), optional:: symm      !in: basis function symmetry
         integer(INTD):: errc
         integer(INTL):: n

         errc=0; n=this%dimsn()
         if(n.gt.0) then
          if(func_num.gt.0.and.func_num.le.n) then
           if(present(symm)) then
            call this%basis_func(func_num)%set(basis_kind,basis_func,errc,symm)
           else
            call this%basis_func(func_num)%set(basis_kind,basis_func,errc)
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
!--------------------------------------------
        subroutine SubspaceBasisDestroy(this)
         implicit none
         type(subspace_basis_t):: this

         if(allocated(this%basis_func)) deallocate(this%basis_func)
         this%space_dim=0
         return
        end subroutine SubspaceBasisDestroy
![subspace_t]=========================================
        subroutine SubspaceInit(this,id,supp_dim,ierr)
!Sets up a subspace. If the subspace is defined on input,
!it will be automatically destructed prior to re-initialization.
         implicit none
         class(subspace_t), intent(out):: this       !out: empty subspace
         integer(INTL), intent(in):: id              !in: subspace id (must be non-negative)
         integer(INTD), intent(in):: supp_dim        !in: subspace support dimension
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(id.ge.0.and.supp_dim.gt.0) then
          this%subspace_id=id
          this%supp_dim=supp_dim
          this%max_resolution=0
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
!-----------------------------------------------------------
        function SubspaceGetSymmetry(this,ierr) result(symm)
         implicit none
         type(space_symmetry_t):: symm               !out: subspace symmetry
         class(subspace_t), intent(in):: this        !in: subspace
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(this%subspace_id.ge.0.and.this%max_resolution.gt.0) then
          symm=this%symm
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function SubspaceGetSymmetry
!-----------------------------------------------------------
        function SubspaceGetCenter(this,ierr) result(center)
         implicit none
         type(real_vec_t):: center                   !out: effective center of the subspace support
         class(subspace_t), intent(in):: this        !in: subspace
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(this%subspace_id.ge.0.and.this%max_resolution.gt.0) then
          center=this%center
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function SubspaceGetCenter
!--------------------------------------------------------------
        function SubspaceGetSupport(this,ierr) result(supp_box)
         implicit none
         type(orthotope_t):: supp_box                !out: subspace support (orthotope)
         class(subspace_t), intent(in):: this        !in: subspace
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(this%subspace_id.ge.0.and.this%max_resolution.gt.0) then
          supp_box=this%supp_box
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function SubspaceGetSupport
!--------------------------------------------------------
        subroutine SubspaceRegisterBasis(this,basis,ierr)
!Registers a new subspace basis.
         implicit none
         class(subspace_t), intent(inout):: this     !inout: subspace
         class(subspace_basis_t), intent(in):: basis !in: new subspace basis
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,ier
         type(list_iter_t):: basis_it

         if(this%subspace_id.ge.0) then
          errc=basis_it%init(this%bases)
          if(errc.eq.GFC_SUCCESS) then
           errc=basis_it%reset()
           if(errc.eq.GFC_SUCCESS) then
            errc=basis_it%append(basis)
            if(errc.eq.GFC_SUCCESS) then
             this%max_resolution=max(this%max_resolution,basis%space_dim)
             call this%update_support(errc)
            else
             errc=5
            endif
           else
            errc=4
           endif
           ier=basis_it%release()
          else
           errc=3
          endif
          ier=basis_it%release(); if(ier.ne.GFC_SUCCESS) errc=2
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine SubspaceRegisterBasis
!----------------------------------------------------------------------
        function SubspaceGetBasis(this,ierr,bas_pred_f) result(basis_p)
!Returns a pointer to the subspace basis satisfying a certain (optional) condition.
!The condition is specified via a GFC predicate object. If no condition is supplied,
!the first basis will be returned.
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
             else !not found (or error)
              errc=2
             endif
            else
             errc=4
            endif
            ier=basis_it%release()
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
        end function SubspaceGetBasis
!--------------------------------------------------
        subroutine SubspaceUpdateSupport(this,ierr)
         implicit none
         class(subspace_t), intent(inout):: this     !inout: subspace
         class(subspace_basis_t), intent(in):: basis !in: subspace basis
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,ier
         integer(INTL):: bf
         type(list_iter_t):: basis_it
         class(gfc_cont_elem_t), pointer:: cont_elem
         class(*), pointer:: elem_value
         class(subspace_basis_t), pointer:: basis_p

         if(this%subspace_id.ge.0) then
          if(this%max_resolution.gt.0) then
           errc=basis_it%init(this%bases)
           if(errc.eq.GFC_SUCCESS) then
            errc=basis_it%reset()
            if(errc.eq.GFC_SUCCESS) then
             ier=basis_it%scanp(return_each=.TRUE.)
             sloop: do while(ier.eq.GFC_IT_ACTIVE)
              cont_elem=>basis_it%pointee(errc)
              if(errc.eq.GFC_SUCCESS) then
               elem_value=>cont_elem%get_value(errc)
               if(errc.eq.GFC_SUCCESS) then
                select type(elem_value)
                class is(subspace_basis_t)
                 basis_p=>elem_value
                 !`
                class default
                 errc=7; exit sloop
                end select
               else
                errc=6; exit sloop
               endif
              else
               errc=5; exit sloop
              endif
              ier=basis_it%scanp(return_each=.TRUE.,skip_current=.TRUE.)
             enddo sloop
            else
             errc=4
            endif
            ier=basis_it%release()
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
        end subroutine SubspaceUpdateSupport
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
![h_space_t]===========================================
        subroutine HSpaceConstruct(this,space_dim,ierr)
!Constructs a hierarchical space representation with a subspace aggregation tree.
         implicit none
         class(h_space_t), intent(inout):: this      !out: hierarchical representation of a vector space (SAT)
         integer(INTL), intent(in):: space_dim       !in: space dimension
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         !`Write
         if(present(ierr)) ierr=errc
         return
        end subroutine HSpaceConstruct
!--------------------------------------
        subroutine HSpaceDestruct(this)
         implicit none
         type(h_space_t):: this
         integer(INTD):: errc
         type(tree_iter_t):: tree_it

         if(allocated(this%subspace)) deallocate(this%subspace)
         errc=tree_it%init(this%aggr_tree)
         if(errc.eq.GFC_SUCCESS) then
          errc=tree_it%reset(); if(errc.eq.GFC_SUCCESS) errc=tree_it%delete_subtree()
          errc=tree_it%release()
         endif
         this%space_dim=0; this%num_subspaces=0
         return
        end subroutine HSpaceDestruct

       end module subspaces
