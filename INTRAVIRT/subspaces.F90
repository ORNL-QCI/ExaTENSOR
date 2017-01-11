!Infrastructure for a recursive adaptive vector space decomposition.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/01/11

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
        end type extent1d_t
 !Real space rectangular hypercube (orthotope):
        type, public:: orthotope_t
         integer(INTL), private:: num_dim=0                   !number of dimensions
         type(extent1d_t), allocatable, private:: extent(:)   !extent of each dimension (min,max)
         !contains
          !procedure, public:: create=>OrthotopeCreate          !creates an empty orthotope
          !procedure, public:: set_extent=>OrthotopeSetExtent   !sets extent along each dimension
          !procedure, public:: lower_bound=>OrthotopeLowerBound !returns the lower bound of a specific extent
          !procedure, public:: upper_bound=>OrthotopeUpperBound !returns the upper bound of a specific extent
          !procedure, public:: length=>OrthotopeLength          !returns the length along a specific extent
          !procedure, public:: volume=>OrthotopeVolume          !returns the volume of the orthotope
          !procedure, public:: overlap=>OrthotopeOverlap        !returns the overlap of two orthotopes
          !procedure, public:: union=>OrthotopeUnion            !returns the minimal orthotope containing two orthotopes
          !final:: OrthotopeDestroy                             !destroys the orthotope
        end type orthotope_t
 !Symmetry:
        type, public:: space_symmetry_t
         integer(INTD), public:: irrep=0          !irreducible representation (non-negative, 0 means no symmetry)
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
        end type basis_func_gauss_t
 !Typed basis function:
        type, public:: basis_func_t
         integer(INTD), private:: basis_kind=BASIS_NONE                 !specific basis kind
         class(basis_func_abs_t), pointer, private:: basis_func=>NULL() !pointer to a basis function of this kind
        end type basis_func_t
 !Subspace basis:
        type, public:: subspace_basis_t
         integer(INTL), private:: space_dim=0                     !number of basis functions
         type(basis_func_t), allocatable, private:: basis_func(:) !basis functions [1..space_dim]
        end type subspace_basis_t
 !Subspace:
        type, public:: subspace_t
         integer(INTL), private:: subspace_id=-1                 !unique subspace ID (registered ID): must be non-negative, -1 means undefined
         integer(INTD), private:: supp_dim=0                     !dimensionality of the real space support on which the basis functions reside
         integer(INTL), private:: max_resolution=0               !max resolution level (max subspace dimension): 0 means undefined
         type(space_symmetry_t), private:: symm                  !symmetry of the subspace (if any)
         type(real_vec_t), private:: center                      !center of the effective subspace support in real space
         type(orthotope_t), private:: supp_box                   !effective subspace support in real space (multidimensional orthotope)
         type(subspace_basis_t), allocatable, private:: basis(:) !basis sets for each resolution level [1..max_resolution]: Defined only for terminal subspaces!
        end type subspace_t
 !Hierarchical composite index:
        type, public:: h_index_t
         integer(INTL), public:: subspace_id=-1 !subspace ID (registered ID): must be non-negative, -1 means undefined
         integer(INTL), public:: resolution=0   !subspace resolution level 1<=resolution<=max_resolution: 0 means undefined
         integer(INTL), public:: component=0    !subspace component number at the given level of resolution: [1..resolution], 0 means undefined
        end type h_index_t
 !Hierarchical space representation:
        type, public:: h_space_t
         integer(INTL), private:: space_dim=0                 !original dimension of the vector space
         integer(INTL), private:: num_subspaces=0             !number of subspaces in the direct-sum decomposition of the vector space
         type(subspace_t), allocatable, private:: subspace(:) !subspaces
         type(tree_t), private:: aggr_tree                    !subspace aggregation tree (SAT)
         !contains
          !procedure, public:: construct=>HSpaceConstruct
          !final:: HSpaceDestruct
        end type h_space_t
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
 !orthotope_t:
 !h_space_t:
        !private HSpaceConstruct

       contains
!IMPLEMENTATION:
![real_vec_t]====================================
        subroutine RealVecCreate(this,dimsn,ierr)
         implicit none
         class(real_vec_t), intent(inout):: this     !out: empty real space vector
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
         type(extent1d_t):: res_extent          !resulting extent (overlap)
         class(extent1d_t), intent(in):: this   !in: extent 1
         class(extent1d_t), intent(in):: extent !in: extent 2

         if(this%max_coord.le.extent%min_coord.or.this%min_coord.ge.extent%max_coord) then
          call res_extent%set(0d0,0d0) !no overlap
         else
          call res_extent%set(max(this%min_coord,extent%min_coord),min(this%max_coord,extent%max_coord))
         endif
         return
        end function Extent1dOverlap

       end module subspaces
