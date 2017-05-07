!Linear-scaling sorting subroutines operating with multi-index keys.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/05/06 (origin 2005 PhD work, used in Mol.Phys.2007)

!Copyright (C) 2005 Dmitry I. Lyakh (Liakh)

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
!The following subroutines sort a list of items according to their unsigned integer multi-index keys.
!Each integer index in a multi-index key must be non-negative!
!The formal scaling of the sorting algorithm is O(L*N), where N is the number of items
!and L is the length of the integer multi-index key (number of integer indices in the key).
!An additional RAM of the same size as the input is allocated: Total RAM ~ O(2*(L+1)*N),
!unless an external buffer is provided.
!There are two versions of the sorting subroutines:
! 1) Disk version (suffix "d"): The input/output is stored on disk:
!    The input file contains the entries to be sorted in the format:
!    <Entry VALUE (integer/real)>  <Entry KEY (integer multi-index)>
!    The (sorted) output will overwrite the input in the same file.
! 2) RAM version: The input/output is stored in RAM:
!    A 1d array of entry values (integer/real);
!    A 2d array of entry keys (each key is an integer multi-index).
!    The (sorted) output will overwrite the input in the memory.
!Additionally, suffix "e" means an external buffer can be provided.

!INTERFACES:
! - multord_r8d(i:ifh,i:n,i:nl,i:mov,i[1]:ip1);
! - multord_r8(i:n,i:nl,i:mov,i[1]:ip1,i[2]:iv,r8[1]:v);
! - multord_id(i:ifh,i:n,i:nl,i:mov,i[1]:ip1);
! - multord_i(i:n,i:nl,i:mov,i[1]:ip1,i[2]:iv,i[1]:v);
! - multord_i8(i8:n,i8:nl,i8:mov,i8[1]:ip1,i8[2]:iv,i8[1]:v);
! - multord_i8e(i8:n,i8:nl,i8:mov,i8[1]:ip1,i8[2]:iv,i8[1]:v,i8[1]:ext_buf);
        module multords
        implicit none
        private
!PARAMETERS:
 !Key length:
        integer, parameter, private:: MAX_MLNDX_LEN=1024    !max length of a multi-index key (keep consistent with MLNDX_FMT)
        character(4), parameter, private:: MLNDX_FMT='1024' !multi-index key read/write format (keep consistent with MAX_MLNDX_LEN)
 !Output:
        integer, private:: CONS_OUT=6     !default output device
        logical, private:: VERBOSE=.TRUE. !verbosity for errors
        logical, private:: DEBUG=.FALSE.  !debugging
!VISIBILITY:
        public multord_r8d !<real(8) value> + <integer multi-index key>: Disk version
        public multord_r8  !<real(8) value> + <integer multi-index key>: RAM version
        public multord_id  !<integer value> + <integer multi-index key>: Disk version
        public multord_i   !<integer value> + <integer multi-index key>: RAM version
        public multord_i8  !<integer(8) value> + <integer(8) multi-index key>: RAM version
        public multord_i8e !<integer(8) value> + <integer(8) multi-index key>: RAM version with an optional buffer

        contains
!-----------------------------------------------
	subroutine multord_r8d(ifh,n,nl,mov,ip1)
!INPUT:
!	ifh - file handle (file with items to be sorted);
!	n   - number of indices in the multiindex keys;
!	nl  - the largest index value encountered in the multi-index keys (the smallest is zero);
!	mov - number of items to be sorted;
!	ip1(1:n) - index place priorities (first encountered 0 terminates sorting, thus allowing partial key comparison);
!OUTPUT:
!       File <ifh> with the sorted item list.
!NOTES:
!	The items are read from the begining of file <ifh> with a subsequent overwriting.
	implicit none
	integer, intent(in):: ifh,n,nl,mov,ip1(1:n)
	integer, parameter:: index_base=2**6 !a base for index splitting: MUST BE EVEN!
	integer i,j,k,l,m,k1,k2,k3,ks,kf,ierr
	integer isp(0:index_base-1),ibp(0:index_base-1),m1,m2,np,kp,fml,index_split,index_part
	character(128) fm
	integer, allocatable:: ifv(:),iv(:,:)
	real(8), allocatable:: v(:)

	index_part(i,j)=mod(i/(index_base**(j-1)),index_base)

	if(n.gt.0.and.nl.gt.0.and.mov.gt.1) then !at least one index with a non-zero limit, and two items
	 if(n.gt.MAX_MLNDX_LEN) then
	  if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_r8d): unsupported multi-index key length: ',n
	  stop
	 endif
!Determine index splitting:
	 index_split=0; k=nl; do while(k.gt.0); k=k/index_base; index_split=index_split+1; enddo
!Create the output format specification:
	 if(nl.gt.0.and.nl.lt.10**9) then
	  k=int(dlog10(dble(nl)))+1
	  fml=16+len(MLNDX_FMT); fm(1:fml)='(d25.15,'//MLNDX_FMT//'(1x,i'//achar(iachar('0')+k)//'))'
	 elseif(nl.eq.0) then
	  fml=16+len(MLNDX_FMT); fm(1:fml)='(d25.15,'//MLNDX_FMT//'(1x,i1))'
	 elseif(nl.ge.10**9.and.nl.le.huge(1_4)) then
	  fml=17+len(MLNDX_FMT); fm(1:fml)='(d25.15,'//MLNDX_FMT//'(1x,i10))'
	 else
	  if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_r8d): invalid key-index upper bound: ',nl
	  stop
	 endif
!Array allocation:
	 allocate(v(mov+mov+1),iv(n,mov+mov+1),ifv(mov+1),STAT=ierr) !work structures
	 if(ierr.eq.0) then
!Read items:
	  rewind(ifh)
	  do k=1,mov
	   read(ifh,*,err=2000,end=2001) v(k),iv(1:n,k)     !get an item from the file: {real number, integer multi-index key}
	  enddo
!Init:
	  ifv(:)=0; ifv(1)=1; ifv(mov+1)=1 !special bound values
	  isp(:)=0
!Begin:
	  m1=0; m2=mov !flip/flop bases
	  iloop: do np=1,n
	   kp=ip1(np)
	   if(kp.gt.0.and.kp.le.n) then
	    do l=index_split,1,-1
	     k=1
	     do while(k.le.mov)
	      ibp(0)=k
 !count:
	      do
	       k2=index_part(iv(kp,k+m1),l); isp(k2)=isp(k2)+1
	       k=k+1; if(ifv(k).ne.0) exit
	      enddo
 !mark:
!	      do k1=1,index_base-1
!	       ibp(k1)=ibp(k1-1)+isp(k1-1); isp(k1-1)=0; ifv(ibp(k1))=1
!	      enddo
!	      ifv(ibp(index_base-1)+isp(index_base-1))=1; isp(index_base-1)=0
	      do k1=1,index_base-3,2
	       j=ibp(k1-1)+isp(k1-1)
	       ibp(k1)=j; ibp(k1+1)=j+isp(k1)
	       ifv(j)=1; ifv(ibp(k1+1))=1
	      enddo
	      j=ibp(index_base-2)+isp(index_base-2); ibp(index_base-1)=j
	      ifv(j)=1; ifv(j+isp(index_base-1))=1
	      isp(0:index_base-1)=0
 !resort:
	      do k1=ibp(0),k-1
	       k2=index_part(iv(kp,k1+m1),l); k3=ibp(k2)+m2; ibp(k2)=ibp(k2)+1
	       v(k3)=v(k1+m1); iv(1:n,k3)=iv(1:n,k1+m1)
	      enddo
	     enddo !next k
	     m=m1; m1=m2; m2=m !switch buffers
	    enddo !next l: index part being analyzed
	   elseif(kp.le.0) then
	    exit iloop
	   else
	    if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_r8d): invalid index priority: ',kp,np
	    stop
	   endif
	  enddo iloop !next np: prioritized index-place
!Save results:
	  rewind(ifh)
	  do k=m1+1,m1+mov
	   write(ifh,fm(1:fml)) v(k),iv(1:n,k)
	  enddo
	  deallocate(v,iv,ifv)
	 else
	  if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_r8d): memory allocation failed: ',mov
	  stop
	 endif
	endif
	return
!------------------
2000	if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_r8d): invalid format of the input file!'; stop
2001	if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_r8d): unexpected end of file while reading the items: ',mov; stop
	end subroutine multord_r8d
!-----------------------------------------------
	subroutine multord_r8(n,nl,mov,ip1,iv,v)
!INPUT:
!	n   - number of indices in the multi-index keys;
!	nl  - the largest index value encountered in the multi-index keys (the smallest is zero);
!	mov - number of items to be sorted;
!	ip1(1:n) - index place priorities (first encountered 0 terminates sorting, thus allowing partial key comparison);
!	iv(1:n,1:mov) - multi-index keys;
!	v(1:mov) - items.
!OUTPUT:
!       iv/v sorted.
	implicit none
	integer, intent(in):: n,nl,mov,ip1(1:n)
	integer, intent(inout), target:: iv(1:n,1:*)
	real(8), intent(inout), target:: v(1:*)
	integer, parameter:: index_base=2**6 !a base for index splitting: MUST BE EVEN!
	integer i,j,k,l,m,k1,k2,k3,ks,kf,ierr
	integer isp(0:index_base-1),ibp(0:index_base-1),m1,m2,np,kp,index_split,index_part
	logical final_copy
	integer, allocatable, target:: ivb(:,:),ifv(:)
	real(8), allocatable, target:: vb(:)
	integer, pointer, contiguous:: iv1(:,:),iv2(:,:)
	real(8), pointer, contiguous:: v1(:),v2(:)

	index_part(i,j)=mod(i/(index_base**(j-1)),index_base)

	if(n.gt.0.and.nl.gt.0.and.mov.gt.1) then !at least one index with a non-zero limit, and two items
	 if(n.gt.MAX_MLNDX_LEN) then
	  if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_r8): unsupported multi-index key length: ',n
	  stop
	 endif
!Determine index splitting:
	 index_split=0; k=nl; do while(k.gt.0); k=k/index_base; index_split=index_split+1; enddo
!Array allocation:
	 allocate(ifv(mov+1),ivb(n,mov),vb(mov),STAT=ierr) !work structures
	 if(ierr.eq.0) then
!Init:
	  ifv(:)=0; ifv(1)=1; ifv(mov+1)=1 !special bound values
	  isp(:)=0
!Begin:
	  iv1=>iv(:,1:mov); iv2=>ivb(:,1:mov); v1=>v(1:mov); v2=>vb(1:mov); final_copy=.false.
	  iloop: do np=1,n
	   kp=ip1(np)
	   if(kp.gt.0.and.kp.le.n) then
	    do l=index_split,1,-1
	     k=1
	     do while(k.le.mov)
	      ibp(0)=k
 !count:
	      do
	       k2=index_part(iv1(kp,k),l); isp(k2)=isp(k2)+1
	       k=k+1; if(ifv(k).ne.0) exit
	      enddo
 !mark:
	      do k1=1,index_base-3,2
	       j=ibp(k1-1)+isp(k1-1)
	       ibp(k1)=j; ibp(k1+1)=j+isp(k1)
	       ifv(j)=1; ifv(ibp(k1+1))=1
	      enddo
	      j=ibp(index_base-2)+isp(index_base-2); ibp(index_base-1)=j
	      ifv(j)=1; ifv(j+isp(index_base-1))=1
	      isp(0:index_base-1)=0
 !resort:
	      do k1=ibp(0),k-1
	       k2=index_part(iv1(kp,k1),l); k3=ibp(k2); ibp(k2)=ibp(k2)+1
	       v2(k3)=v1(k1); iv2(1:n,k3)=iv1(1:n,k1)
	      enddo
	     enddo !next k
	     if(.not.final_copy) then !switch buffers
	      iv1=>ivb(:,1:mov); iv2=>iv(:,1:mov); v1=>vb(1:mov); v2=>v(1:mov); final_copy=.true.
	     else
	      iv1=>iv(:,1:mov); iv2=>ivb(:,1:mov); v1=>v(1:mov); v2=>vb(1:mov); final_copy=.false.
	     endif
	    enddo !next l: index part being analyzed
	   elseif(kp.le.0) then
	    exit iloop
	   else
	    if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_r8): invalid index priority: ',kp,np
	    stop
	   endif
	  enddo iloop !next np: prioritized index-place
!Save results:
	  if(final_copy) then
	   iv(1:n,1:mov)=ivb(1:n,1:mov); v(1:mov)=vb(1:mov)
	  endif
	  nullify(iv1,iv2,v1,v2); deallocate(ifv,ivb,vb)
	 else
	  if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_r8): memory allocation failed: ',mov
	  stop
	 endif
	endif
	return
	end subroutine multord_r8
!----------------------------------------------
	subroutine multord_id(ifh,n,nl,mov,ip1)
!INPUT:
!	ifh - file handle (file with items to be sorted);
!	n   - number of indices in the multiindex keys;
!	nl  - the largest index value encountered in the multi-index keys (the smallest is zero);
!	mov - number of items to be sorted;
!	ip1(1:n) - index place priorities (first encountered 0 terminates sorting, thus allowing partial key comparison);
!OUTPUT:
!       File <ifh> with the sorted item list.
!NOTES:
!	The items are read from the begining of file <ifh> with a subsequent overwriting.
	implicit none
	integer, intent(in):: ifh,n,nl,mov,ip1(1:n)
	integer, parameter:: index_base=2**6 !a base for index splitting: MUST BE EVEN!
	integer i,j,k,l,m,k1,k2,k3,ks,kf,ierr
	integer isp(0:index_base-1),ibp(0:index_base-1),m1,m2,np,kp,fml,index_split,index_part
	character(128) fm
	integer, allocatable:: ifv(:),iv(:,:)
	integer, allocatable:: v(:)

	index_part(i,j)=mod(i/(index_base**(j-1)),index_base)

	if(n.gt.0.and.nl.gt.0.and.mov.gt.1) then !at least one index with a non-zero limit, and two items
	 if(n.gt.MAX_MLNDX_LEN) then
	  if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_id): unsupported multi-index key length: ',n
	  stop
	 endif
!Determine index splitting:
	 index_split=0; k=nl; do while(k.gt.0); k=k/index_base; index_split=index_split+1; enddo
!Create the output format specification:
	 if(nl.gt.0.and.nl.lt.10**9) then
	  k=int(dlog10(dble(nl)))+1
	  fml=13+len(MLNDX_FMT); fm(1:fml)='(i11,'//MLNDX_FMT//'(1x,i'//achar(iachar('0')+k)//'))'
	 elseif(nl.eq.0) then
	  fml=13+len(MLNDX_FMT); fm(1:fml)='(i11,'//MLNDX_FMT//'(1x,i1))'
	 elseif(nl.ge.10**9.and.nl.le.huge(1_4)) then
	  fml=14+len(MLNDX_FMT); fm(1:fml)='(i11,'//MLNDX_FMT//'(1x,i10))'
	 else
	  if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_id): invalid key-index upper bound: ',nl
	  stop
	 endif
!Array allocation:
	 allocate(v(mov+mov+1),iv(n,mov+mov+1),ifv(mov+1),STAT=ierr) !work structures
	 if(ierr.eq.0) then
!Read items:
	  rewind(ifh)
	  do k=1,mov
	   read(ifh,*,err=2000,end=2001) v(k),iv(1:n,k)     !get an item from the file: {real number, integer multi-index key}
	  enddo
!Init:
	  ifv(:)=0; ifv(1)=1; ifv(mov+1)=1 !special bound values
	  isp(:)=0
!Begin:
	  m1=0; m2=mov !flip/flop bases
	  iloop: do np=1,n
	   kp=ip1(np)
	   if(kp.gt.0.and.kp.le.n) then
	    do l=index_split,1,-1
	     k=1
	     do while(k.le.mov)
	      ibp(0)=k
 !count:
	      do
	       k2=index_part(iv(kp,k+m1),l); isp(k2)=isp(k2)+1
	       k=k+1; if(ifv(k).ne.0) exit
	      enddo
 !mark:
	      do k1=1,index_base-3,2
	       j=ibp(k1-1)+isp(k1-1)
	       ibp(k1)=j; ibp(k1+1)=j+isp(k1)
	       ifv(j)=1; ifv(ibp(k1+1))=1
	      enddo
	      j=ibp(index_base-2)+isp(index_base-2); ibp(index_base-1)=j
	      ifv(j)=1; ifv(j+isp(index_base-1))=1
	      isp(0:index_base-1)=0
 !resort:
	      do k1=ibp(0),k-1
	       k2=index_part(iv(kp,k1+m1),l); k3=ibp(k2)+m2; ibp(k2)=ibp(k2)+1
	       v(k3)=v(k1+m1); iv(1:n,k3)=iv(1:n,k1+m1)
	      enddo
	     enddo !next k
	     m=m1; m1=m2; m2=m !switch buffers
	    enddo !next l: index part being analyzed
	   elseif(kp.le.0) then
	    exit iloop
	   else
	    if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_id): invalid index priority: ',kp,np
	    stop
	   endif
	  enddo iloop !next np: prioritized index-place
!Save results:
	  rewind(ifh)
	  do k=m1+1,m1+mov
	   write(ifh,fm(1:fml)) v(k),iv(1:n,k)
	  enddo
	  deallocate(v,iv,ifv)
	 else
	  if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_id): memory allocation failed: ',mov
	  stop
	 endif
	endif
	return
!------------------
2000	if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_id): invalid format of the input file!'; stop
2001	if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_id): unexpected end of file while reading the items: ',mov; stop
	end subroutine multord_id
!----------------------------------------------
	subroutine multord_i(n,nl,mov,ip1,iv,v)
!INPUT:
!	n   - number of indices in the multi-index keys;
!	nl  - the largest index value encountered in the multi-index keys (the smallest is zero);
!	mov - number of items to be sorted;
!	ip1(1:n) - index place priorities (first encountered 0 terminates sorting, thus allowing partial key comparison);
!	iv(1:n,1:mov) - multi-index keys;
!	v(1:mov) - items.
!OUTPUT:
!       iv/v sorted.
	implicit none
	integer, intent(in):: n,nl,mov,ip1(1:n)
	integer, intent(inout), target:: iv(1:n,1:*)
	integer, intent(inout), target:: v(1:*)
	integer, parameter:: index_base=2**6 !a base for index splitting: MUST BE EVEN!
	integer i,j,k,l,m,k1,k2,k3,ks,kf,ierr
	integer isp(0:index_base-1),ibp(0:index_base-1),m1,m2,np,kp,index_split,index_part
	logical final_copy
	integer, allocatable, target:: ivb(:,:),ifv(:)
	integer, allocatable, target:: vb(:)
	integer, pointer, contiguous:: iv1(:,:),iv2(:,:)
	integer, pointer, contiguous:: v1(:),v2(:)

	index_part(i,j)=mod(i/(index_base**(j-1)),index_base)

	if(n.gt.0.and.nl.gt.0.and.mov.gt.1) then !at least one index with a non-zero limit, and two items
	 if(n.gt.MAX_MLNDX_LEN) then
	  if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_i): unsupported multi-index key length: ',n
	  stop
	 endif
!Determine index splitting:
	 index_split=0; k=nl; do while(k.gt.0); k=k/index_base; index_split=index_split+1; enddo
!Array allocation:
	 allocate(ifv(mov+1),ivb(n,mov),vb(mov),STAT=ierr) !work structures
	 if(ierr.eq.0) then
!Init:
	  ifv(:)=0; ifv(1)=1; ifv(mov+1)=1 !special bound values
	  isp(:)=0
!Begin:
	  iv1=>iv(:,1:mov); iv2=>ivb(:,1:mov); v1=>v(1:mov); v2=>vb(1:mov); final_copy=.false.
	  iloop: do np=1,n
	   kp=ip1(np)
	   if(kp.gt.0.and.kp.le.n) then
	    do l=index_split,1,-1
	     k=1
	     do while(k.le.mov)
	      ibp(0)=k
 !count:
	      do
	       k2=index_part(iv1(kp,k),l); isp(k2)=isp(k2)+1
	       k=k+1; if(ifv(k).ne.0) exit
	      enddo
 !mark:
	      do k1=1,index_base-3,2
	       j=ibp(k1-1)+isp(k1-1)
	       ibp(k1)=j; ibp(k1+1)=j+isp(k1)
	       ifv(j)=1; ifv(ibp(k1+1))=1
	      enddo
	      j=ibp(index_base-2)+isp(index_base-2); ibp(index_base-1)=j
	      ifv(j)=1; ifv(j+isp(index_base-1))=1
	      isp(0:index_base-1)=0
 !resort:
	      do k1=ibp(0),k-1
	       k2=index_part(iv1(kp,k1),l); k3=ibp(k2); ibp(k2)=ibp(k2)+1
	       v2(k3)=v1(k1); iv2(1:n,k3)=iv1(1:n,k1)
	      enddo
	     enddo !next k
	     if(.not.final_copy) then !switch buffers
	      iv1=>ivb(:,1:mov); iv2=>iv(:,1:mov); v1=>vb(1:mov); v2=>v(1:mov); final_copy=.true.
	     else
	      iv1=>iv(:,1:mov); iv2=>ivb(:,1:mov); v1=>v(1:mov); v2=>vb(1:mov); final_copy=.false.
	     endif
	    enddo !next l: index part being analyzed
	   elseif(kp.le.0) then
	    exit iloop
	   else
	    if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_i): invalid index priority: ',kp,np
	    stop
	   endif
	  enddo iloop !next np: prioritized index-place
!Save results:
	  if(final_copy) then
	   iv(1:n,1:mov)=ivb(1:n,1:mov); v(1:mov)=vb(1:mov)
	  endif
	  nullify(iv1,iv2,v1,v2); deallocate(ifv,ivb,vb)
	 else
	  if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_i): memory allocation failed: ',mov
	  stop
	 endif
	endif
	return
	end subroutine multord_i
!-----------------------------------------------
	subroutine multord_i8(n,nl,mov,ip1,iv,v)
!INPUT:
!	n   - number of indices in the multi-index keys;
!	nl  - the largest index value encountered in the multi-index keys (the smallest is zero);
!	mov - number of items to be sorted;
!	ip1(1:n) - index place priorities (first encountered 0 terminates sorting, thus allowing partial key comparison);
!	iv(1:n,1:mov) - multi-index keys;
!	v(1:mov) - items.
!OUTPUT:
!       iv/v sorted.
	implicit none
	integer(8), intent(in):: n,nl,mov,ip1(1:n)
	integer(8), intent(inout), target:: iv(1:n,1:*)
	integer(8), intent(inout), target:: v(1:*)
	integer(8), parameter:: index_base=2**6 !a base for index splitting: MUST BE EVEN!
	integer(8) i,j,k,l,m,k1,k2,k3,ks,kf
        integer(4):: ierr
	integer(8) isp(0:index_base-1),ibp(0:index_base-1),m1,m2,np,kp,index_split,index_part
	logical final_copy
	integer(8), allocatable, target:: ivb(:,:),ifv(:)
	integer(8), allocatable, target:: vb(:)
	integer(8), pointer, contiguous:: iv1(:,:),iv2(:,:)
	integer(8), pointer, contiguous:: v1(:),v2(:)

	index_part(i,j)=mod(i/(index_base**(j-1)),index_base)

	if(n.gt.0_8.and.nl.gt.0_8.and.mov.gt.1_8) then !at least one index with a non-zero limit, and two items
	 if(n.gt.MAX_MLNDX_LEN) then
	  if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_i8): unsupported multi-index key length: ',n
	  stop
	 endif
!Determine index splitting:
	 index_split=0_8; k=nl; do while(k.gt.0_8); k=k/index_base; index_split=index_split+1_8; enddo
!Array allocation:
	 allocate(ifv(mov+1_8),ivb(n,mov),vb(mov),STAT=ierr) !work structures
	 if(ierr.eq.0) then
!Init:
	  ifv(:)=0_8; ifv(1_8)=1_8; ifv(mov+1_8)=1_8 !special bound values
	  isp(:)=0_8
!Begin:
	  iv1=>iv(:,1_8:mov); iv2=>ivb(:,1_8:mov); v1=>v(1_8:mov); v2=>vb(1_8:mov); final_copy=.false.
	  iloop: do np=1_8,n
	   kp=ip1(np)
	   if(kp.gt.0_8.and.kp.le.n) then
	    do l=index_split,1_8,-1_8
	     k=1_8
	     do while(k.le.mov)
	      ibp(0_8)=k
 !count:
	      do
	       k2=index_part(iv1(kp,k),l); isp(k2)=isp(k2)+1_8
	       k=k+1_8; if(ifv(k).ne.0_8) exit
	      enddo
 !mark:
	      do k1=1_8,index_base-3_8,2_8
	       j=ibp(k1-1_8)+isp(k1-1_8)
	       ibp(k1)=j; ibp(k1+1_8)=j+isp(k1)
	       ifv(j)=1_8; ifv(ibp(k1+1_8))=1_8
	      enddo
	      j=ibp(index_base-2_8)+isp(index_base-2_8); ibp(index_base-1_8)=j
	      ifv(j)=1_8; ifv(j+isp(index_base-1_8))=1_8
	      isp(0_8:index_base-1_8)=0_8
 !resort:
	      do k1=ibp(0_8),k-1_8
	       k2=index_part(iv1(kp,k1),l); k3=ibp(k2); ibp(k2)=ibp(k2)+1_8
	       v2(k3)=v1(k1); iv2(1_8:n,k3)=iv1(1_8:n,k1)
	      enddo
	     enddo !next k
	     if(.not.final_copy) then !switch buffers
	      iv1=>ivb(:,1_8:mov); iv2=>iv(:,1_8:mov); v1=>vb(1_8:mov); v2=>v(1_8:mov); final_copy=.true.
	     else
	      iv1=>iv(:,1_8:mov); iv2=>ivb(:,1_8:mov); v1=>v(1_8:mov); v2=>vb(1_8:mov); final_copy=.false.
	     endif
	    enddo !next l: index part being analyzed
	   elseif(kp.le.0_8) then
	    exit iloop
	   else
	    if(VERBOSE) write(CONS_OUT,*)'ERRROR(multord_i8): invalid index priority: ',kp,np
	    stop
	   endif
	  enddo iloop !next np: prioritized index-place
!Save results:
	  if(final_copy) then
	   iv(1_8:n,1_8:mov)=ivb(1_8:n,1_8:mov); v(1_8:mov)=vb(1_8:mov)
	  endif
	  nullify(iv1,iv2,v1,v2); deallocate(ifv,ivb,vb)
	 else
	  if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_i8): memory allocation failed: ',mov
	  stop
	 endif
	endif
	return
	end subroutine multord_i8
!--------------------------------------------------------
        subroutine multord_i8e(n,nl,mov,ip1,iv,v,ext_buf)
        implicit none
        integer(4), intent(in):: n                      !in: number of indices in the multi-index keys
        integer(8), intent(in):: nl                     !in: the largest index value encountered in the multi-index keys (the smallest is zero)
        integer(8), intent(in):: mov                    !in: number of items to be sorted
        integer(4), intent(in):: ip1(1:n)               !in: index place priorities (first encountered 0 terminates sorting, thus allowing partial key comparison)
        integer(8), intent(inout), target:: iv(1:n,1:*) !inout: multi-index keys
        integer(8), intent(inout), target:: v(1:*)      !inout: items
        integer(8), intent(in), contiguous, target, optional:: ext_buf(1:) !in: external buffer (to avoid memory allocation)
        !-----------------------------------------------------------------
        integer(8), parameter:: INDEX_BASE=2**6 !base for index splitting: MUST BE EVEN!
        !--------------------------------------
        integer(4):: ierr,np,kp
        integer(8):: i,j,k,l,m,k1,k2,k3,ks,kf
        integer(8):: isp(0:INDEX_BASE-1),ibp(0:INDEX_BASE-1),m1,m2,index_split,index_part
        integer(8), pointer, contiguous:: ivb(:,:),vb(:),ifv(:),iv1(:,:),iv2(:,:),v1(:),v2(:)
        logical:: mem_allocated,final_copy

        index_part(i,j)=mod(i/(INDEX_BASE**(j-1)),INDEX_BASE)

        if(n.gt.0.and.nl.gt.0_8.and.mov.gt.1_8) then !at least one index with a non-zero limit, and two items
!Determine index splitting:
         index_split=0_8; k=nl; do while(k.gt.0_8); k=k/INDEX_BASE; index_split=index_split+1_8; enddo
!Work buffer allocation:
         mem_allocated=.FALSE.; ierr=0
         if(present(ext_buf)) then
          if(size(ext_buf).gt.mov*int(n+2,8)) then
           l=0_8
           ivb(1:n,1_8:mov)=>ext_buf(l+1_8:l+mov*int(n,8))
           l=l+mov*int(n,8)
           vb(1_8:mov)=>ext_buf(l+1_8:l+mov)
           l=l+mov
           ifv(1_8:mov+1_8)=>ext_buf(l+1_8:l+mov+1_8)
          else
           ierr=-1
          endif
         else
          allocate(ivb(n,mov),vb(mov),ifv(mov+1_8),STAT=ierr) !work structures
          if(ierr.eq.0) mem_allocated=.TRUE.
         endif
         if(ierr.eq.0) then
!Init:
          ifv(:)=0_8; ifv(1_8)=1_8; ifv(mov+1_8)=1_8 !special bound values
          isp(:)=0_8
!Begin:
          iv1=>iv(:,1_8:mov); iv2=>ivb(:,1_8:mov); v1=>v(1_8:mov); v2=>vb(1_8:mov); final_copy=.FALSE.
          iloop: do np=1,n
           kp=ip1(np)
           if(kp.gt.0.and.kp.le.n) then
            do l=index_split,1_8,-1_8
             k=1_8
             do while(k.le.mov)
              ibp(0_8)=k
 !count:
              do
               k2=index_part(iv1(kp,k),l); isp(k2)=isp(k2)+1_8
               k=k+1_8; if(ifv(k).ne.0_8) exit
              enddo
 !mark:
              do k1=1_8,INDEX_BASE-3_8,2_8
               j=ibp(k1-1_8)+isp(k1-1_8)
               ibp(k1)=j; ibp(k1+1_8)=j+isp(k1)
               ifv(j)=1_8; ifv(ibp(k1+1_8))=1_8
              enddo
              j=ibp(INDEX_BASE-2_8)+isp(INDEX_BASE-2_8); ibp(INDEX_BASE-1_8)=j
              ifv(j)=1_8; ifv(j+isp(INDEX_BASE-1_8))=1_8
              isp(0_8:INDEX_BASE-1_8)=0_8
 !resort:
              do k1=ibp(0_8),k-1_8
               k2=index_part(iv1(kp,k1),l); k3=ibp(k2); ibp(k2)=ibp(k2)+1_8
               v2(k3)=v1(k1); iv2(1:n,k3)=iv1(1:n,k1)
              enddo
             enddo !next k
             if(.not.final_copy) then !switch buffers
              iv1=>ivb(:,1_8:mov); iv2=>iv(:,1_8:mov); v1=>vb(1_8:mov); v2=>v(1_8:mov); final_copy=.TRUE.
             else
              iv1=>iv(:,1_8:mov); iv2=>ivb(:,1_8:mov); v1=>v(1_8:mov); v2=>vb(1_8:mov); final_copy=.FALSE.
             endif
            enddo !next l: index part being analyzed
           elseif(kp.le.0) then
            exit iloop
           else
            if(VERBOSE) write(CONS_OUT,*)'ERRROR(multord_i8e): invalid index priority: ',kp,np
            stop
           endif
          enddo iloop !next np: prioritized index-place
!Save results:
          if(final_copy) then
           iv(1:n,1_8:mov)=ivb(1:n,1_8:mov); v(1_8:mov)=vb(1_8:mov)
          endif
          nullify(iv1,iv2,v1,v2)
          if(mem_allocated) deallocate(ifv,vb,ivb)
          nullify(ifv,vb,ivb)
         else
          if(VERBOSE) write(CONS_OUT,*)'ERROR(multord_i8e): insufficient memory: ',mov
          stop
         endif
        endif
        return
        end subroutine multord_i8e

        end module multords
!==========================
!TESTING:
!===========================
        module multords_test
        use multords
        use timers, only: thread_wtime
        implicit none
        private
        public test_multord_sort
        contains
!------------------------------------------------------------
        function test_multord_sort(perf,dev_out) result(ierr)
         integer(4):: ierr
         real(8), intent(out):: perf
         integer(4), intent(in), optional:: dev_out
         integer(4), parameter:: KEYLEN=8
         integer(8), parameter:: NITEMS=1000000_8,UBND=9999_8
         integer(8), allocatable:: iv(:,:),v(:),buf(:)
         integer(4):: ipr(1:KEYLEN),devo,i
         integer(8):: l
         real(8):: tms,tm,rndv(1:KEYLEN)

         ierr=0; perf=0d0; devo=6; if(present(dev_out)) devo=dev_out
         allocate(iv(1:KEYLEN,NITEMS),v(NITEMS),buf(NITEMS*(KEYLEN+3)),STAT=ierr)
         if(ierr.ne.0) then; ierr=1; return; endif
         do l=1_8,NITEMS
          v(l)=l
          call random_number(rndv)
          do i=1,KEYLEN; iv(i,l)=int(rndv(i)*dble(UBND),8); enddo
         enddo
         ipr(1:KEYLEN)=(/(i,i=1,KEYLEN)/)
         tms=thread_wtime()
         call multord_i8e(KEYLEN,UBND,NITEMS,ipr,iv,v,buf)
         tm=thread_wtime(tms); perf=dble(NITEMS)/tm
         cloop: do l=2_8,NITEMS
          do i=1,KEYLEN
           if(iv(i,l).lt.iv(i,l-1_8)) then
            ierr=2; exit cloop
           elseif(iv(i,l).gt.iv(i,l-1_8)) then
            exit
           endif
          enddo
         enddo cloop
         if(allocated(buf)) deallocate(buf)
         if(allocated(v)) deallocate(v)
         if(allocated(iv)) deallocate(iv)
         return
        end function test_multord_sort

        end module multords_test
