	module combinatoric
!Combinatoric Procedures.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!Revision: 2014/09/24
!All rights reserved! No single part can be taken or reproduced!
!PROCEDURES:
! - TRNG(i:ctrl,i:ni,i[1]:trn,i[1]:ngt): permutation generator which returns each new permutation.
! - TRSIGN(i:n,i[1]:itr): reorders indices in an ascending order and determines the sign of the corresponding permutation (Bubble). Use MERGE_SORT for fast.
! - i8:FACTORIAL(i:n): factorial of a number.
! - i:NOID(i:m,i:n): this function returns a binomial coefficient (integer).
! - i8:NOID8(i8:m,i8:n): this function returns a binomial coefficient (integer*8).
! - GPGEN(i:ctrl,i:ni,i[1]:vh,i[1]:trn,i[2]:cil): permutation generator with partial-ordering restrictions, returns each new permutation.
! - TR_CYCLE(i:ni,i[1]:trn,i:nc,i[2]:cyc): decomposes a permutation into permutation cycles and determines the sign of the permutation.
! - l:PERM_TRIVIAL(i:ni,i[1]:trn): checks whether the given permutation is trivial or not.
! - l:PERM_TRIVIAL_INT8(i8:ni,i8[1]:trn): checks whether the given permutation is trivial or not.
! - l:PERM_OK(i:ni,i[1]:trn): checks whether the given permutation is legitimate or not.
! - PERM2TRANS(i:ni,i[1]:trn1,i[1]:trn2,i:ntrp,i[2]:trp): creates a list of elementary transpositions relating one permutation to another.
! - PERMUTATION_CONVERTER(l:seq2pos,i:ni,i[1]:n2o,i[1]:o2n): converts between two permutation representations (N2O and O2N).
! - i:HASH_ARR_INT(i:hash_range,i:ni,i[1]:arr): returns a hash-mask for an integer array.
! - i:HASH_ARR_INT8(i:hash_range,i:ni,i8[1]:arr): returns a hash-mask for an integer*8 array.
! - i:CMP_MULTINDS(i:ml1,i[1]:m1,i:ml2,i[1]:m2): compares two integer multiindices.
! - i:CMP_MULTINDS_INT8(i:ml1,i8[1]:m1,i:ml2,i8[1]:m2): compares two integer*8 multiindices.
! - MULTINDX_MERGE(i:ml1,i[1]:m1,i:ml2,i[1]:m2,i:mlr,i[1]:mr,i:sign_corr): merges two multiindices, providing a potential sign correction.
! - i:CMP_ARRAYS_INT(l:preorder,i:ml1,i[1]:m1,i:ml2,i[1]:m2,i[1]:trn): compares two integer arrays with an optional preodering.
! - i:CMP_ARRAYS_INT8(l:preorder,i:ml1,i8[1]:m1,i:ml2,i8[1]:m2,i[1]:trn): compares two integer(8) arrays with an optional preodering.
! - CLANAL(i:nt,r8:ctol,r8[1]:dta,i:ncl,i[1]:cdta,r8[1]:cmv): cluster analysis of a one-dimensional set of points.
! - RANDOM_PERMUTATION(i:ni,i[1]:trn,l:no_trivial): returns a random permutation (default integer).
! - RANDOM_PERMUTATION_INT8(i8:ni,i8[1]:trn,l:no_trivial): returns a random permutation (integer*8 variant).
! - RANDOM_COMPOSITION(l:ordered,i:irange,i:ni,i[1]:trn): returns a random sequence of ni natural numbers from the range [1..irange] without repeats.
! - MERGE_SORT_INT(i:ni,i[1]:trn): fast sorting algorithm for an integer array (default integer).
! - MERGE_SORT_KEY_INT(i:ni,i[1]:key,i[1]:trn): fast sorting algorithm for an integer array, based on integer keys (default integer).
! - MERGE_SORT_INT8(i8:ni,i8[1]:trn): fast sorting algorithm for an integer*8 array.
! - MERGE_SORT_KEY_INT8(i8:ni,i8[1]:key,i8[1]:trn): fast sorting algorithm for an integer*8 array, based on integer keys (integer*8).
! - MERGE_SORT_REAL8(i:ni,r8[1]:trn): fast sorting algorithm for a real*8 array.
! - MERGE_SORT_KEY_REAL8(i:ni,r8[1]:key,i[1]:trn): fast sorting algorithm for an integer array, based on real*8 keys.
! - MERGE_SORT_CMPLX8(i:ni,c8[1]:trn): fast sorting algorithm for a complex*8 array.
! - MERGE_SORT_KEY_CMPLX8(i:ni,c8[1]:key,i[1]:trn): fast sorting algorithm for an integer array, based on comlex*8 keys.
! - SORT_SLOTS(i:ni,slot[1]:keys,i[1]:trn): fast sorting algorithm for an integer permutation using keys of type slot.
! - NULLIFY_SPARSE_GRAPH(graph_sparse:gr): cleans an object of type graph_sparse.
! - COPY_SPARSE_GRAPH(l:clean,graph_sparse:gr1,graph_sparse:gr2): copies gr1 to gr2 (gr2 is pre-cleaned if clean is .true.).
! - i:COMPARE_DESCR_SLOTS(slot:sl1,slot:sl2): compares two descriptor slots: -1:sl1<sl2, 0:sl1=sl2, +1:sl1>sl2.
! - SORT_VERTEX_CONNECTIONS(graph_sparse:gr): rearranges vertex slots (connections) in an ordered form.
! - GRAPH_STANDARD(graph_sparse:gr,i[1]:trn,i:auto,i:ierr): returns the standard permutation for an arbitrary graph (CRA algorithm of D.I.Lyakh,2011).
! - l:GRAPH_ISOMORPHIC(graph_sparse:g1,graph_sparse:g2,i:ierr): this function returns true if the two graphs are isomorphic (CRA algorithm of D.I.Lyakh,2011).
! - r8:DINT2FRAC(r8:dpi,i:prec): converts a double-precision integer into a fractional number (<1d0) by reflecting its decimal digits agains the decimal point (e.g., 28560. --> .06582).
! - r8:FRAC2DINT(r8:dpf,i:prec): inverse of DINT2FRAC (e.g., .06582 --> 28560.).
! - GRAPH2MATRIX(graph_sparse:gr,i:nv,c8[2]:adj): returns the adjacency matrix adj(1:nv,1:nv) for an object gr of type(graph_sparse).
! - MATRIX2GRAPH(i:nv,c8[2]:adj,graph_sparse:gr): returns an object gr of type(graph_sparse) for the adjacency matrix adj(1:nv,1:nv).
! - r8:CHECK_HERMITICITY(i:nv,c8[2]:adj,l:correct): computes the maximal deviation from hermiticity for an arbitrary complex(8) matrix (+optional correction).
! - STAR_GRAPH(i:nconn,i:gdim,c8[2]:graph): returns a star-graph of cardinality gdim (as an adjacency matrix) with a half-width nconn.
! - RANDOM_GRAPH(i:color_range,i:valence_range,i:bond_range,i:gdim,c8[2]:graph): returns a random graph adjacency matrix graph(1:gdim,1:gdim) of cardinality gdim.
! - REMOVE_RANDOM_EDGES(i:edge_num,i:gd,c8[2]:graph): randomly removes edges from an adjacency matrix graph(1:gd,1:gd).
!NOTES:
! - As a rule, a permutation is passed as an integer array (0:n), where the 0th element is the sign of the permutation {1..n}.
! - TRSIGN uses the buble algorithm => will be very slow for long permutations.
! - For some subroutines, the permutation {1..n} must contain all integers from the segment [1..n].
! - Arbitrary (directed, colored) (multi)-graphs can be represented either by a complex(8) hermitian adjacency matrix, or by the type(graph_sparse).
!   Any given (directed) (multiple) edge/loop (connection/self-connection) is represented by a complex number encoding its particular kind.
!   All procedures working with graphs use only comparisons (subtractions). No other operations, like additions, multiplications, etc. are performed.
!   It is EXTREMELY important to ensure stability of the double-precision number representation, say, up to 12 decimal digits of precision
!   (because of the digit-wise comparisons)! Hermiticity of the adjacency matrix must be STRICTLY enforced (digit-wise)!
!   Wherever appropriate, the zero threshold for numerical comparisons is 1d-13 (which is an order of magnitude lower than the double-precision epsilon).
!   A double-precision representation restricts the sum of the number of colors and the number of possible undirected connection kinds to 12 decimal digits,
!   since the diagonal elements of the adjacency matrix are of the form: <color.number_of_loops> ~ 12 decimal digits of guaranteed precision.
!   Complex off-diagonal elements have the form: real=<number_of_undirected_edges.number_of_compensated_directed_edge_pairs> ~ 12 digits;
!                                                imag=<number_of_uncompensated_directed_edges.0> ~ 12 digits.
!   Above, the fractional part of each real number is generated as an inverted sequence of digits of the corresponding integer connection_kind,
!   whereas the integer part of each real number is translated directly. Since the total maximal number of digits in a real number is set to 12,
!   the integer and fractional parts are given 6 digits each, thus restricting the max number of colors and the max number of connection kinds to 1000000.
! - On some machines, DINT2FRAC and FRAC2DINT can be numerically unstable: 1572/10000 = 0.1572 --> 0.15719999999999.
!   This is, however, being checked, leading to an exception and termination of the program.

!GLOBAL PARAMETERS:
	integer, parameter:: dp_codon_len=6 !number of decimal digits in a double-precision codon (codon is a sequence of decimal digits representing either the integer or the fractional part of the real number)
	real(8), parameter:: dp_zero_thresh=1d-13 !certified guaranteed precision of a double-precision real number (one order lower than the epsilon)
!DATA TYPES:
        type slot !describes a particular (possbily multiple, possibly directed) connection of a certain vertex
         integer:: vert_attr=0           !attribute of the connected vertex (can be either VERTEX_ID or VERTEX_CLASS, depending on the purpose)
	 complex(8):: conn_typ=(0d0,0d0) !connection type: off-diag: real=x.y: x is the amount of undirected edges, y is the amount of directed compensated edge-pairs; imag=+/-v.w: v is the amount of uncompensated directed edges (sign shows which), w is reserved
	end type slot                    !                     diag: real=0.y: y is the amount of loops (color is taken out); imag=0.0

	type vertex                         !describes a vertex (color + connections to other vertices)
         real(8):: color=0d0                !vertex color (should be integer, except some special cases)
         integer:: conn_num=0               !number of connections of the current vertex (possibly including the only possible self-connection)
         type(slot), allocatable:: slots(:) !vertex connections (VERT_ATTR=VERTEX_ID): SLOTS(1:CONN_NUM)
        end type vertex

	type descriptor                     !vertex descriptor
	 integer:: vertex_class=0           !vertex class
         integer:: conn_num=0               !length of the descriptor (number of connections of the current vertex)
         type(slot), allocatable:: slots(:) !slots (VERT_ATTR=VERTEX_CLASS): SLOTS(1:CONN_NUM)
        end type descriptor

        type classification                 !describes a classification
         integer:: length=0                 !length of the classification
	 integer:: num_classes=0            !number of classes: [1..num_classes]
	 integer, allocatable:: class(:)    !classes of the vertices (by vertex ID): CLASS# = class(VERTEX_ID)
	 integer, allocatable:: new_id(:)   !new positions of the vertices (by vertex ID): NEW# = new_id(VERTEX_ID)
        end type classification

        type graph_sparse                         !desribes an arbitrary graph (colored directed multigraph). Vertex numeration should begin with 1.
         integer:: gcard=0                        !graph cardinality (number of vertices)
	 integer:: graph_typ=0                    !graph type (default is 0: non-specific)
	 logical:: descr_allocated=.false.        !.true. if descr(:) is allocated
	 type(vertex), allocatable:: vertices(:)  !description of the vertices: a vertex should not have more than one (complex) connection to another vertex
	 type(descriptor), allocatable:: descr(:) !vertex descriptors: auxiliary: allocated whenever needed
	end type graph_sparse
!----------------------------------------------------------------------------------------------------------
	contains
!METHODS:
	subroutine trng(ctrl,ni,trn,ngt)
!Permutation generator: Returns each subsequent permutation.
! CTRL - control argument. The first call must have CTRL<>0 and initialized TRN!
!        Then, on a regular output, CTRL will be 0 unless permutations are over (CTRL=-1).
! NI - number of indices to permute;
! TRN - current permutation (TRN(0) is the sign of the permutation);
! NGT - work structure. Do not change it outside this subroutine!
!NOTES:
! - Permutation index values are not restricted (any array of integer numbers).
	implicit none
	integer, intent(in):: ni
	integer, intent(inout):: ctrl,trn(0:*),ngt(0:*)
	integer j,k

	if(ctrl.ne.0) then !first call: NGT initialization. The original permutation is returned.
	 ngt(0)=0; do j=1,ni; ngt(j)=j-1; enddo
	 ctrl=0
	else !subsequent call: a new permutation is returned.
	 k=1 !maybe k=2 will accelerate
	 do while(k.le.ni)
	  if(ngt(k).ne.k-1) call transp(k,ngt(k)+1)
	  if(ngt(k).ne.0) then
	   call transp(k,ngt(k))
	   ngt(k)=ngt(k)-1
	   return
	  else
	   ngt(k)=k-1
	   k=k+1
	  endif
	 enddo
	 ctrl=-1
	endif
	return

	contains

	 subroutine transp(m,n)
	  implicit none
	  integer m,n,l
	  l=trn(m); trn(m)=trn(n); trn(n)=l; trn(0)=-trn(0)
	  return
	 end subroutine transp

	end subroutine trng
!------------------------------------------------------------------
	subroutine trsign(n,itr) !ITR - permutation of the length N
!This subroutine orders integers in ITR in an ascending order.
! ITR - integers. ITR(0) will be the sign of the ordered permutation.
! N - number of integers to be ordered.
!NOTES:
! - Permutation index values are not restricted (any array of integers).
! - Bubble algorithm (only for small arrays).
	implicit none
	integer, intent(in):: n
	integer, intent(inout):: itr(0:*)
	integer k,l,isgn
	k=1
	isgn=+1
	do while(k.lt.n)
	 if(itr(k).gt.itr(k+1)) then
	  l=itr(k); itr(k)=itr(k+1); itr(k+1)=l; isgn=-isgn
	  if(k.gt.1) then; k=k-1; else; k=2; endif
	 else
	  k=k+1
	 endif
	enddo
	itr(0)=isgn
	return
	end subroutine trsign
!-------------------------------------------------------------------
	integer(8) function factorial(n) !returns N! for N>=0, and -1 otherwise
	implicit none
	integer, intent(in):: n
	integer(8) k

	if(n.ge.0) then
	 factorial=1_8
	 do k=2_8,int(n,8)
	  factorial=factorial*k
	  if(factorial.lt.0_8) then; write(*,*)'ERROR(combinatoric:factorial): integer(8) overflow!'; stop; endif !trap
	 enddo
	else
	 factorial=-1_8
	endif
	return
	end function factorial
!------------------------------------------------------------------------------------------------------
	integer function noid(m,n) !returns the number of unique distributions of N objects on M places
	implicit none
	integer, intent(in):: m,n
	integer k,l

	if(n.gt.m.or.n.lt.0.or.m.lt.0) then
	 noid=0
	 return
	elseif(n.eq.m.or.n.eq.0) then
	 noid=1
	 return
	endif
	noid=1; l=m
	do k=1,n; noid=noid*l/k; l=l-1; enddo
	if(noid.le.0) then; write(*,*)'ERROR(combinatoric:noid): integer overflow:',m,n,noid; stop; endif !trap
	return
	end function noid
!----------------------------------------------------------------------------------------------------------
	integer(8) function noid8(m,n) !returns the number of unique distributions of N objects on M places
	implicit none
	integer(8), intent(in):: m,n
	integer(8) k,l

	if(n.gt.m.or.n.lt.0.or.m.lt.0) then
	 noid8=0_8
	 return
	elseif(n.eq.m.or.n.eq.0) then
	 noid8=1_8
	 return
	endif
	noid8=1_8; l=m
	do k=1,n; noid8=noid8*l/k; l=l-1; enddo
	if(noid8.le.0_8) then; write(*,*)'ERROR(combinatoric:noid8): integer*8 overflow: ',m,n,noid8; stop; endif !trap
	return
	end function noid8
!-------------------------------------------
	subroutine gpgen(ctrl,ni,vh,trn,cil)
!This subroutine generates all unique permutations of NI items,
!in which the items belonging to the same host are always ordered.
!INPUT:
! - ctrl - control argument: at the begining must be <>0; 0 - next permutation; -1 - permutations are over.
! - ni - number of items;
! - vh(1:ni) - index hosts;
!INPUT(dummy)/OUTPUT:
! - trn(0:ni) - current permutation (trn(0) is the sign);
! - cil(0:1,0:ni) - connected list (for internal use only, do not set or change it outside!).
!OUTPUT:
! - trn(0:ni) - current permutation (trn(0) is the sign).
!NOTES:
! - The first permutation is created here (fully ordered one).
!   Its sign is +1. Each next permutation is generated from the previous one.
	implicit none
!------------------------------------------
	integer, parameter:: max_item=16384 !max allowed number of items (because of the permutation sign determination, see below)
!------------------------------------------
	integer i,j,k,l,m,n,k1,k2,k3,k4,k5,k6,ks,kf,ierr
	integer, intent(in):: ni,vh(1:ni)
	integer, intent(inout):: ctrl,trn(0:ni),cil(0:1,0:ni)

	if(ni.gt.0) then
	 if(ni.gt.max_item) then
	  write(*,*)'ERROR(gpgen): legnth of the permutation exceeds the maximal value: ',max_item,ni
	  stop
	 endif
	 if(ctrl.ne.0) then
	  call first_call
	 else
!free the last box:
	  m=vh(ni); n=ni
	  do while(vh(n).eq.m)
	   call free_item(trn(n)); n=n-1
	   if(n.eq.0) exit
	  enddo
!get the next composition:
	  if(n.gt.0) then
	   ks=-1
	   do while(n.gt.0)
!	    write(*,'(''DEBUG1: '',128(i1,1x))') trn(1:n) !debug
!	    write(*,'(4x,128(i2,1x))') (k6,k6=0,ni); write(*,'(4x,128(i2,1x))') cil(0,0:ni) !debug
!           write(*,'(4x,128(i2,1x))') cil(1,0:ni) !debug
	    if(ks.lt.0) then
	     i=trn(n); call free_item(i); j=cil(1,i)
	     if(j.gt.0) then
	      m=cil(1,j); call engage_item(j); trn(n)=j; ks=+1
	     endif
	    else
	     if(vh(n).eq.vh(n-1)) then
	      j=m
	      if(j.gt.0) then
	       m=cil(1,j); call engage_item(j); trn(n)=j
	      else
	       ks=-1
	      endif
	     else
	      j=cil(0,0); m=cil(1,j); call engage_item(j); trn(n)=j
	     endif
	    endif
!	    write(*,'(''DEBUG2: '',128(i1,1x))') trn(1:n) !debug
!	    write(*,'(4x,128(i2,1x))') (k6,k6=0,ni); write(*,'(4x,128(i2,1x))') cil(0,0:ni) !debug
!           write(*,'(4x,128(i2,1x))') cil(1,0:ni) !debug
	    n=n+ks
	    if(n.gt.ni) then !success
	     call determine_trn_sign
	     return
	    endif
	   enddo !n
	   ctrl=-1 !end
	  else
	   ctrl=-1
	  endif
	 endif
	else
	 ctrl=-1
	endif
	return

	contains

	 subroutine first_call
	  implicit none
	  integer j1
	  trn(0)=+1; do j1=1,ni; trn(j1)=j1; enddo
	  cil(0:1,1:ni)=-1; cil(0:1,0)=(/ni+1,0/)
	  ctrl=0
	  return
	 end subroutine first_call

	 subroutine free_item(it)
	  implicit none
	  integer, intent(in):: it
	  integer j1,j2
	  if(it.lt.cil(0,0)) then
	   if(cil(0,0).le.ni) then
	    cil(0:1,it)=(/0,cil(0,0)/); cil(0,cil(0,0))=it; cil(0,0)=it
	   else
	    cil(0:1,it)=(/0,0/); cil(0:1,0)=it
	   endif
	  elseif(it.gt.cil(1,0)) then
	   if(cil(1,0).ge.1) then
	    cil(0:1,it)=(/cil(1,0),0/); cil(1,cil(1,0))=it; cil(1,0)=it
	   else
	    cil(0:1,it)=(/0,0/); cil(0:1,0)=it
	   endif
	  else
	   if(it-cil(0,0).le.cil(1,0)-it) then !moving up: insert the new item between the minimal and maximal elements
	    j1=cil(0,0); do while(cil(1,j1).lt.it); j1=cil(1,j1); enddo
	    j2=cil(1,j1); cil(0:1,it)=(/j1,j2/); cil(1,j1)=it; cil(0,j2)=it
	   else !moving down: insert the new item between the minimal and maximal elements
	    j1=cil(1,0); do while(cil(0,j1).gt.it); j1=cil(0,j1); enddo
	    j2=cil(0,j1); cil(0:1,it)=(/j2,j1/); cil(1,j2)=it; cil(0,j1)=it
	   endif
	  endif
	  return
	 end subroutine free_item

	 subroutine engage_item(it)
	  implicit none
	  integer, intent(in):: it
	  integer j1,jd,ju
	  jd=1; ju=1
	  if(it.eq.cil(0,0)) then
	   j1=cil(1,it)
	   if(j1.gt.0) then; cil(0,j1)=0; cil(0,0)=j1; else; cil(0,0)=ni+1; endif
	   jd=0
	  endif
	  if(it.eq.cil(1,0)) then
	   j1=cil(0,it)
	   if(j1.gt.0) then; cil(1,j1)=0; cil(1,0)=j1; else; cil(1,0)=0; endif
	   ju=0
	  endif
	  if(jd*ju.eq.1) then; cil(1,cil(0,it))=cil(1,it); cil(0,cil(1,it))=cil(0,it); endif
	  cil(0:1,it)=-1
	  return
	 end subroutine engage_item

	 subroutine determine_trn_sign
	  implicit none
	  integer occ(max_item),j1,j2,j3,js
	  js=+1; occ(1:ni)=0; j1=1; j2=0; j3=0
	  do
	   if(occ(j1).eq.0) then
	    j3=j3+1; j2=j2+1; occ(j1)=1; j1=trn(j1)
	   else
	    if(mod(j2,2).eq.0) js=-js
	    if(j3.eq.ni) exit
	    j2=0; j1=1; do while(occ(j1).ne.0); j1=j1+1; enddo
	   endif
	  enddo
	  trn(0)=js
	  return
	 end subroutine determine_trn_sign

	end subroutine gpgen
!-----------------------------------------
	subroutine tr_cycle(ni,trn,nc,cyc)
!This subroutine extracts permutation cycles and the sign from a given permutation.
!INPUT:
! - ni - number of indices;
! - trn(0:ni) - permutation, in which trn(0) is the sign returned;
!OUTPUT:
! - trn(0) - sign of the permutation;
! - nc - number of permutation cycles;
! - cyc(0:1,1:ni) - permutation cycles: cyc(1,:) is an index of a cycle; cyc(0,:) is the number of the cycle the index belongs to.
!NOTE:
! - nc=-666 - means ERROR.
! - permutation index values must lie within the range [1..ni].
	implicit none
	integer, intent(in):: ni
	integer, intent(inout):: trn(0:*)
	integer, intent(out):: nc,cyc(0:1,*)
	integer, parameter:: max_in_mem=1024
	integer i,j,k,l,m,n,k1,k2,k3,k4,k5,k6,ks,kf,ierr
	integer, target:: ibuss(1:max_in_mem)
	integer, allocatable, target:: ibusa(:)
	integer, pointer:: ibus(:)

	nc=0
	if(ni.gt.0) then
	 if(ni.gt.max_in_mem) then; allocate(ibusa(1:ni)); ibus=>ibusa; else; ibus=>ibuss; endif
	 if(trn_ok()) then
!	  ibus(1:ni)=0 !busy flags
	  trn(0)=+1; n=0; m=0
	  do while(n.lt.ni)
	   i=m+1; do while(ibus(i).ne.0); i=i+1; enddo; m=i
	   nc=nc+1; l=i; j=0
	   do
	    n=n+1; j=j+1; cyc(0:1,n)=(/nc,trn(i)/); ibus(i)=nc
	    if(trn(i).eq.l) then; exit; else; i=trn(i); endif
	   enddo
	   if(mod(j,2).eq.0) trn(0)=-trn(0)
	  enddo
	 else
	  trn(0)=-667; nc=-666
	 endif
	 nullify(ibus); if(ni.gt.max_in_mem) deallocate(ibusa)
	else
	 trn(0)=-666; nc=-666
	endif
	return

	contains

	 logical function trn_ok()
	  integer j1
	  ibus(1:ni)=0
	  do j1=1,ni
	   if(trn(j1).le.0.or.trn(j1).gt.ni) then; trn_ok=.false.; return; endif
	   if(ibus(trn(j1)).ne.0) then; trn_ok=.false.; return; endif
	   ibus(trn(j1))=j1
	  enddo
	  ibus(1:ni)=0
	  trn_ok=.true.
	  return
	 end function trn_ok

	end subroutine tr_cycle
!--------------------------------------------
	logical function perm_trivial(ni,trn)
!Checks whether the given permutation trn(0:ni) is trivial or not.
	implicit none
	integer, intent(in):: ni,trn(0:*)
	integer i
	perm_trivial=.true.
	do i=1,ni; if(trn(i).ne.i) then; perm_trivial=.false.; exit; endif; enddo
	return
	end function perm_trivial
!-------------------------------------------------
	logical function perm_trivial_int8(ni,trn)
!Checks whether the given permutation trn(0:ni) is trivial or not.
	implicit none
	integer(8), intent(in):: ni,trn(0:*)
	integer(8) i
	perm_trivial_int8=.true.
	do i=1,ni; if(trn(i).ne.i) then; perm_trivial_int8=.false.; exit; endif; enddo
	return
	end function perm_trivial_int8
!---------------------------------------
	logical function perm_ok(ni,trn)
!Checks whether the given permutation trn(0:ni) is legitimate or not.
!NOTE: keep the permutation fit into the stack!
	implicit none
	integer, intent(in):: ni,trn(0:*)
	integer i,j,ibus(1:ni)
	ibus(1:ni)=0
	do i=1,ni
	 j=trn(i)
	 if(j.le.0.or.j.gt.ni) then; perm_ok=.false.; return; endif
	 if(ibus(j).ne.0) then; perm_ok=.false.; return; else; ibus(j)=i; endif
	enddo
	perm_ok=.true.
	return
	end function perm_ok
!---------------------------------------------------
	subroutine perm2trans(ni,trn1,trn2,ntrp,trp)
!This subroutine creates a sequence of elementary pair transpositions
!which connect one permutation (old) with another one (new).
!The sign of the second (new) permutation is changed accordingly.
!Schematically: TRN1 ---list_of_transpositions---> TRN2 with a modified sign.
!INPUT:
! - ni - number of indices;
! - trn1(0:ni) - old permutation (&0 - sign);
! - trn2(0:ni) - new permutation (&0 - sign);
!OUTPUT:
! - ntrp - number of elementary transpositions;
! - trp(1:2,1:ntrp) - elementary transpositions in terms of old index numbers.
!                     Each transposition interchanges two old indices according to their original numbers;
! - trn2(0) - sign of the permutation TRN2 with respect to the permutation TRN1.
!NOTES:
! - Permutation index values must span the range [1..ni].
! - Because AUTOMATIC arrays are employed, the number of indices cannot be too large.
!   Switch to pointers if you want to process larger permutations.
	implicit none
	integer, intent(in):: ni,trn1(0:*),trn2(0:*)
	integer, intent(out):: ntrp,trp(2,*)
	integer i,j,k,l,m,n,k1,k2,k3,k4,k5,k6,ks,kf,ierr
	integer trn(1:ni),ipos(1:ni)
	integer isgn

	ntrp=0
	if(ni.gt.1) then
	 if(trn_ok()) then
	  isgn=+1; do j=1,ni; trn(j)=trn1(j); enddo; do j=1,ni; ipos(trn1(j))=j; enddo
	  do i=1,ni
	   if(trn(i).ne.trn2(i)) then
	    ntrp=ntrp+1; trp(1:2,ntrp)=(/trn(i),trn2(i)/)
	    j=trn(i); trn(i)=trn2(i); trn(ipos(trn2(i)))=j
	    ipos(j)=ipos(trn2(i)); ipos(trn2(i))=i
	    isgn=-isgn
	   endif
	  enddo
	 else
	  write(*,*)'ERROR(perm2trans): invalid input permutation: ',ni,trn1(1:ni),trn2(1:ni)
	  stop
	 endif
	endif
	return

	contains

	 logical function trn_ok()
	  integer j1
	  trn_ok=.true.
	  do j1=1,ni
	   if(trn1(j1).lt.1.or.trn1(j1).gt.ni.or.trn2(j1).lt.1.or.trn2(j1).gt.ni) then; trn_ok=.false.; exit; endif
	  enddo
	  return
	 end function trn_ok

	end subroutine perm2trans
!-----------------------------------------------------------
	subroutine permutation_converter(seq2pos,ni,n2o,o2n)
!This subroutine converts between two permutation representations:
!(n2o) new_to_old: a sequence of original item numbers (position --> ID);
!(o2n) old_to_new: new positions of items (ID --> new position).
!Briefly: n2o(POSITION)=ID; o2n(ID)=POSITION. They are inverse permutations.
!INPUT:
! - seq2pos - .TRUE. means a conversion from n2o to o2n, .FALSE. o2n to n2o;
! - ni - number of items;
! - n2o(0:ni)/o2n(0:ni);
!OUTPUT:
! - n2o(0:ni)/o2n(0:ni).
!NOTES:
! - The elements of n2o(1:ni)/o2n(1:ni) must span the range [1..ni] (no argument-validity check is done here!).
	implicit none
	logical, intent(in):: seq2pos
	integer, intent(in):: ni
	integer, intent(inout):: n2o(0:ni),o2n(0:ni)
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf,ierr

	if(seq2pos) then
	 o2n(0)=n2o(0); do i=1,ni; o2n(n2o(i))=i; enddo !loop over positions
	else
	 n2o(0)=o2n(0); do i=1,ni; n2o(o2n(i))=i; enddo !loop over item IDs
	endif
	return
	end subroutine permutation_converter
!-------------------------------------------------------
	integer function hash_arr_int(hash_range,ni,arr)
!This function returns a hash-mask for a given integer array.
!INPUT:
! - hash_range - defines the range of this function [0..hash_range-1];
! - ni - number of items in the array;
! - arr(0:ni-1) - items;
!OUTPUT:
! - hash_arr_int - hash-mask.
	implicit none
	integer, intent(in):: hash_range,ni
	integer, intent(in):: arr(0:*)
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf,ierr

	hash_arr_int=0
	if(hash_range.gt.0.and.ni.ge.0) then
	 do i=0,ni-1
	  hash_arr_int=hash_arr_int+mod(arr(i),hash_range)
	  hash_arr_int=mod(hash_arr_int,hash_range)
	 enddo
	else
	 write(*,*)'ERROR(combinatoric:hash_arr_int): invalid arguments: ',hash_range,ni
	 stop
	endif
	return
	end function hash_arr_int
!--------------------------------------------------------
	integer function hash_arr_int8(hash_range,ni,arr)
!This function returns a hash-mask for a given integer*8 array.
!INPUT:
! - hash_range - defines the range of this function [0..hash_range-1];
! - ni - number of items in the array;
! - arr(0:ni-1) - items (integer*8);
!OUTPUT:
! - hash_arr_int - hash-mask.
	implicit none
	integer, intent(in):: hash_range,ni
	integer(8), intent(in):: arr(0:*)
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf,ierr
	integer(8) hr

	hash_arr_int8=0
	if(hash_range.gt.0.and.ni.ge.0) then
	 hr=int(hash_range,8)
	 do i=0,ni-1
	  hash_arr_int8=hash_arr_int8+int(mod(arr(i),hr),4)
	  hash_arr_int8=mod(hash_arr_int8,hash_range)
	 enddo
	else
	 write(*,*)'ERROR(combinatoric:hash_arr_int8): invalid arguments: ',hash_range,ni
	 stop
	endif
	return
	end function hash_arr_int8
!---------------------------------------------------
	integer function cmp_multinds(ml1,m1,ml2,m2)
!This function compares two integer multiindices.
!INPUT:
! - m1(1:ml1) - 1st multiindex;
! - m2(1:ml2) - 2nd multiindex.
!OUTPUT:
! - cmp_multinds: -X (1<=X<=ml1,ml1=ml2): first difference in two multiindices occurs at position X (left-to-right) and m1(X)<m2(X);
!                 +X (1<=X<=ml1,ml1=ml2): first difference in two multiindices occurs at position X (left-to-right) and m1(X)>m2(X);
!                 -X (X>ml1,X>ml2): first multiindex is shorter than the second one (ml1<ml2);
!                 +X (X>ml1,X>ml2): first multiindex is longer than the second one (ml1>ml2);
!                  0 (ml1=ml2): the two multiindices are equal.
	implicit none
	integer, intent(in):: ml1,ml2,m1(1:*),m2(1:*)
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf,ierr

	if(ml1.ge.0.and.ml2.ge.0) then
	 if(ml1.eq.ml2) then
	  cmp_multinds=0
	  do i=1,ml1
	   if(m1(i).ne.m2(i)) then; cmp_multinds=sign(i,m1(i)-m2(i)); exit; endif
	  enddo
	 elseif(ml1.lt.ml2) then
	  cmp_multinds=-(max(ml1,ml2)+1)
	 else
	  cmp_multinds=+(max(ml1,ml2)+1)
	 endif
	else
	 write(*,'("ERROR(combinatoric:cmp_multinds): negative multiindex length passed: ",i10,1x,i10)') ml1,ml2
	 stop
	endif
	return
	end function cmp_multinds
!--------------------------------------------------------
	integer function cmp_multinds_int8(ml1,m1,ml2,m2)
!This function compares two integer*8 multiindices.
!INPUT:
! - m1(1:ml1) - 1st multiindex (integer*8);
! - m2(1:ml2) - 2nd multiindex (integer*8).
!OUTPUT:
! - cmp_multinds: -X (1<=X<=ml1,ml1=ml2): first difference in two multiindices occurs at position X (left-to-right) and m1(X)<m2(X);
!                 +X (1<=X<=ml1,ml1=ml2): first difference in two multiindices occurs at position X (left-to-right) and m1(X)>m2(X);
!                 -X (X>ml1,X>ml2): first multiindex is shorter than the second one (ml1<ml2);
!                 +X (X>ml1,X>ml2): first multiindex is longer than the second one (ml1>ml2);
!                  0 (ml1=ml2): the two multiindices are equal.
	implicit none
	integer, intent(in):: ml1,ml2
	integer(8), intent(in):: m1(1:*),m2(1:*)
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf,ierr

	if(ml1.ge.0.and.ml2.ge.0) then
	 if(ml1.eq.ml2) then
	  cmp_multinds_int8=0
	  do i=1,ml1
	   if(m1(i).ne.m2(i)) then; cmp_multinds_int8=int(sign(int(i,8),m1(i)-m2(i)),4); exit; endif
	  enddo
	 elseif(ml1.lt.ml2) then
	  cmp_multinds_int8=-(max(ml1,ml2)+1)
	 else
	  cmp_multinds_int8=+(max(ml1,ml2)+1)
	 endif
	else
	 write(*,'("ERROR(combinatoric:cmp_multinds_int8): negative multiindex length passed: ",i10,1x,i10)') ml1,ml2
	 stop
	endif
	return
	end function cmp_multinds_int8
!----------------------------------------------------------------
	subroutine multindx_merge(ml1,m1,ml2,m2,mlr,mr,sign_corr)
!This subroutine merges two multiindices (with index reodering).
!INPUT:
! - ml1 - length of the 1st multiindex;
! - m1(1:ml1) - 1st multiindex (left);
! - ml2 - length of the 2nd multiindex;
! - m2(1:ml2) - 2nd multiindex (right);
!OUTPUT:
! - mlr - length of the multiindex-result (ml1+ml2);
! - mr(1:mlr) - multiindex-result;
! - sign_corr - sign correction (+1/-1/0): 0 when a repeated index is present.
!NOTES:
! - No error checks.
	implicit none
	integer, intent(in):: ml1,m1(1:*),ml2,m2(1:*)
	integer, intent(out):: mlr,mr(1:*),sign_corr
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf,ierr

	mlr=0; sign_corr=+1
!merge:
	if(ml1.gt.0.and.ml2.gt.0) then
	 k1=1; k2=1
	 mloop: do
	  mlr=mlr+1
	  if(m1(k1).le.m2(k2)) then
	   mr(mlr)=m1(k1)
	   if(k1.lt.ml1) then
	    k1=k1+1
	   else
	    l=ml2-k2+1; mr(mlr+1:mlr+l)=m2(k2:ml2); mlr=mlr+l
	    exit mloop
	   endif
	  else
	   mr(mlr)=m2(k2); sign_corr=sign_corr*(-1+2*mod(ml1-k1,2))
	   if(k2.lt.ml2) then
	    k2=k2+1
	   else
	    l=ml1-k1+1; mr(mlr+1:mlr+l)=m1(k1:ml1); mlr=mlr+l
	    exit mloop
	   endif
	  endif
	 enddo mloop
	elseif(ml1.gt.0.and.ml2.le.0) then
	 mlr=ml1; mr(1:ml1)=m1(1:ml1)
	elseif(ml1.le.0.and.ml2.gt.0) then
	 mlr=ml2; mr(1:ml2)=m2(1:ml2)
	endif
!check index repeats:
	do i=1,mlr-1; if(mr(i).eq.mr(i+1)) then; sign_corr=0; exit; endif; enddo
	return
	end subroutine multindx_merge
!------------------------------------------------------------------
	integer function cmp_arrays_int(preorder,ml1,m1,ml2,m2,trn)
!This function compares two integer arrays with an optional preodering.
!INPUT:
! - preorder - if .true., both arrays will be ordered before the comparison;
! - ml1 - length of the 1st array;
! - m1(1:ml1) - 1st array;
! - ml2 - length of the 2nd array;
! - m2(1:ml2) - 2nd array;
!OUTPUT:
! - cmp_arrays_int: -X (1<=X<=ml1,ml1=ml2): first difference in two arrays occurs at position X (left-to-right) and m1(X)<m2(X);
!                   +X (1<=X<=ml1,ml1=ml2): first difference in two arrays occurs at position X (left-to-right) and m1(X)>m2(X);
!                   -X (X>ml1,X>ml2): first array is shorter than the second one (ml1<ml2);
!                   +X (X>ml1,X>ml2): first array is longer than the second one (ml1>ml2);
!                    0 (ml1=ml2): the two arrays are equal.
! - trn(0:) - if preorder=.true. and the arrays are equal (cmp_arrays_int=0),
!             this (optional) output array will contain the permutation matching the two arrays:
!             A new order of old elements of the 2nd array (N2O) that matches the 1st array (both with the original ordering).
	implicit none
	logical, intent(in):: preorder
	integer, intent(in):: ml1,ml2
	integer, intent(in):: m1(1:ml1),m2(1:ml2)
	integer, intent(out), optional:: trn(0:*)
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf,ierr
	integer, allocatable:: prm1(:),prm2(:)

	if(ml1.ge.0.and.ml2.ge.0) then
	 if(ml1.lt.ml2) then
	  cmp_arrays_int=-(max(ml1,ml2)+1)
	 elseif(ml1.gt.ml2) then
	  cmp_arrays_int=(max(ml1,ml2)+1)
	 else !ml1=ml2
	  cmp_arrays_int=0
	  if(preorder) then
	   allocate(prm1(0:ml1),prm2(0:ml2),STAT=ierr)
	   if(ierr.ne.0) then; write(*,*)'ERROR(combinatoric:cmp_arrays_int): memory allocation failed!'; stop; endif
	   prm1(0)=+1; do i=1,ml1; prm1(i)=i; enddo
	   prm2(0)=+1; do i=1,ml2; prm2(i)=i; enddo
	   call merge_sort_key_int(ml1,m1,prm1) !prm1 is N2O
	   call merge_sort_key_int(ml2,m2,prm2) !prm2 is N2O
	   do i=1,ml1
	    if(m1(prm1(i)).ne.m2(prm2(i))) then; cmp_arrays_int=sign(i,m1(prm1(i))-m2(prm2(i))); exit; endif
	   enddo
	   if(cmp_arrays_int.eq.0.and.present(trn)) then
	    trn(0)=prm1(0)*prm2(0); do i=1,ml1; trn(prm1(i))=prm2(i); enddo
	   endif
	   deallocate(prm1,prm2)
	  else
	   do i=1,ml1
	    if(m1(i).ne.m2(i)) then; cmp_arrays_int=sign(i,m1(i)-m2(i)); exit; endif
	   enddo
	  endif
	 endif
	else !invalid arguments
	 write(*,*)'ERROR(combinatoric:cmp_arrays_int): invalid arguments: ',ml1,ml2
	 stop
	endif
	return
	end function cmp_arrays_int
!-------------------------------------------------------------------
	integer function cmp_arrays_int8(preorder,ml1,m1,ml2,m2,trn)
!This function compares two integer(8) arrays with an optional preodering.
!INPUT:
! - preorder - if .true., both arrays will be ordered before the comparison;
! - ml1 - length of the 1st array;
! - m1(1:ml1) - 1st array;
! - ml2 - length of the 2nd array;
! - m2(1:ml2) - 2nd array;
!OUTPUT:
! - cmp_arrays_int8: -X (1<=X<=ml1,ml1=ml2): first difference in two arrays occurs at position X (left-to-right) and m1(X)<m2(X);
!                    +X (1<=X<=ml1,ml1=ml2): first difference in two arrays occurs at position X (left-to-right) and m1(X)>m2(X);
!                    -X (X>ml1,X>ml2): first array is shorter than the second one (ml1<ml2);
!                    +X (X>ml1,X>ml2): first array is longer than the second one (ml1>ml2);
!                     0 (ml1=ml2): the two arrays are equal.
! - trn(0:) - if preorder=.true. and the arrays are equal (cmp_arrays_int8=0),
!             this (optional) output array will contain the permutation matching the two arrays:
!             A new order of old elements of the 2nd array (N2O) that matches the 1st array (both with the original ordering).
	implicit none
	logical, intent(in):: preorder
	integer, intent(in):: ml1,ml2
	integer(8), intent(in):: m1(1:ml1),m2(1:ml2)
	integer, intent(out), optional:: trn(0:*)
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf,ierr
	integer(8), allocatable:: prm1(:),prm2(:)
	integer(8) ml

	if(ml1.ge.0.and.ml2.ge.0) then
	 if(ml1.lt.ml2) then
	  cmp_arrays_int8=-(max(ml1,ml2)+1)
	 elseif(ml1.gt.ml2) then
	  cmp_arrays_int8=(max(ml1,ml2)+1)
	 else !ml1=ml2
	  cmp_arrays_int8=0; ml=int(ml1,8)
	  if(preorder) then
	   allocate(prm1(0:ml),prm2(0:ml),STAT=ierr)
	   if(ierr.ne.0) then; write(*,*)'ERROR(combinatoric:cmp_arrays_int8): memory allocation failed!'; stop; endif
	   prm1(0)=+1_8; do i=1,ml1; prm1(i)=int(i,8); enddo
	   prm2(0)=+1_8; do i=1,ml2; prm2(i)=int(i,8); enddo
	   call merge_sort_key_int8(ml,m1,prm1) !prm1 is N2O
	   call merge_sort_key_int8(ml,m2,prm2) !prm2 is N2O
	   do i=1,ml1
	    if(m1(prm1(i)).ne.m2(prm2(i))) then; cmp_arrays_int8=int(sign(int(i,8),m1(prm1(i))-m2(prm2(i))),4); exit; endif
	   enddo
	   if(cmp_arrays_int8.eq.0.and.present(trn)) then
	    trn(0)=int(prm1(0)*prm2(0),4); do i=1,ml1; trn(prm1(i))=int(prm2(i),4); enddo
	   endif
	   deallocate(prm1,prm2)
	  else
	   do i=1,ml1
	    if(m1(i).ne.m2(i)) then; cmp_arrays_int8=int(sign(int(i,8),m1(i)-m2(i)),4); exit; endif
	   enddo
	  endif
	 endif
	else !invalid arguments
	 write(*,*)'ERROR(combinatoric:cmp_arrays_int8): invalid arguments: ',ml1,ml2
	 stop
	endif
	return
	end function cmp_arrays_int8
!--------------------------------------------------
	subroutine clanal(nt,ctol,dta,ncl,cdta,cmv)
!This subroutine clusterizes a one-dimensional set of points.
!INPUT:
! - nt - number of points;
! - ctol - clusterization tolerance;
! - dta(1:nt) - points;
!OUTPUT:
! - ncl - number of classes;
! - cdta(1:nt) - class of each point;
! - cmv(1:ncl) - class mean values.
!NOTES:
! - if no points has been passed, NCL is set to -1.
	implicit none
	integer i,j,k,l,m,n,k1,k2,k3,k4,k5,k6,ks,kf,ierr
	integer, intent(in):: nt        !number of points
	real(8), intent(in):: ctol      !clusterization tolerance (0..1)
	real(8), intent(in):: dta(*)    !data (points)
	integer, intent(out):: ncl      !number of classes
	integer, intent(out):: cdta(*)  !class# the point belongs to
	real(8), intent(out):: cmv(*)   !class mean values
	real(8) wh,sw,minv,maxv,val
	integer nfup,nfdp

	if(nt.gt.0) then
!calculate range of the values:
	 minv=dta(1); maxv=dta(1)
	 do k=2,nt
	  minv=min(dta(k),minv); maxv=max(dta(k),maxv)
	 enddo
	 wh=maxv-minv   !width of the distribution
	 sw=wh/dble(nt) !specific width
!init cdta:
	 cdta(1:nt)=0   !all points are not classified at the begining
!clusterize:
	 ncl=1         !current class #
	 cdta(1)=1     !element 1 is assigned to class 1
	 nfup=2        !number of the first unclassified point
	 n=1           !will be the number of points in the current class
	 cmv(1)=dta(1)
	 do            !until all the points are classified
	  nfdp=0                     !number of the 1st declined point
	  do k=nfup,nt
	   if(cdta(k).eq.0) then
	    val=dist(k,ncl)
	    if(val.lt.ctol*sw) then  !accepted to the current class
	     cdta(k)=ncl; cmv(ncl)=cmv(ncl)+dta(k)
	     n=n+1
	    else                     !declined
	     if(nfdp.eq.0) nfdp=k
	    endif
	   endif
	  enddo
	  cmv(ncl)=cmv(ncl)/dble(n)  !mean value of the current class
	  if(nfdp.gt.0) then
	   ncl=ncl+1; cdta(nfdp)=ncl
	   nfup=nfdp+1; cmv(ncl)=dta(nfdp)
	   n=1
	  else
	   exit  !all points have been classified
	  endif
	 enddo
	else
	 ncl=-1  !no input data found
	endif   !nt>0: points exist

	return

	contains

	 real(8) function dist(i_point,class_num)
	  integer, intent(in):: i_point,class_num
	  integer npc,lp
	  dist=0d0
	  npc=0      !number of points found for the class#class_num
	  do lp=1,nt
	   if(cdta(lp).eq.class_num) then
	    dist=dist+abs(dta(i_point)-dta(lp))
	    npc=npc+1
	   endif
	  enddo
	  dist=dist/dble(npc)
	  return
	 end function dist

	end subroutine clanal
!-------------------------------------------------------
	subroutine random_permutation(ni,trn,no_trivial)
!This subroutine returns a random permutation of NI items [1..NI].
!INPUT:
! - ni - number of items, range [1:ni], if ni<=0 nothing will be done;
!OUTPUT:
! - trn(0:ni) - generated permutation (trn(0) - sign of the permutation);
	implicit none
	integer, intent(in):: ni
	integer, intent(out):: trn(0:ni)
	logical, intent(in), optional:: no_trivial
!----------------------------------------------
	integer, parameter:: random_chunk=2**10 !size of the chunk of random numbers generated in one call
	integer, parameter:: num_repeats=5      !the bigger the number, the better the generator quality (more expensive)
!----------------------------------------------
	integer i,j,k,l,m,n,nr,ierr
	real(8):: ra(1:random_chunk)

	if(ni.gt.0) then
	 trn(0)=+1; do i=1,ni; trn(i)=i; enddo !initial permutation
	 ploop: do
	  do nr=1,num_repeats
	   do i=1,ni,random_chunk
	    l=min(i+random_chunk-1,ni)-i+1
	    call random_number(ra(1:l))
	    ra(1:l)=ra(1:l)*2d0
	    do k=1,l
	     if(ra(k).lt.1d0) then
	      j=i+k-1
	      n=int(ra(k)*dble(ni))+1; if(n.gt.ni) n=ni
	      if(n.ne.j) then; m=trn(j); trn(j)=trn(n); trn(n)=m; trn(0)=-trn(0); endif
	     endif
	    enddo
	   enddo
	  enddo
	  if(present(no_trivial)) then
	   if(no_trivial) then
	    if(.not.perm_trivial(ni,trn)) exit ploop
	   else
	    exit ploop
	   endif
	  else
	   exit ploop
	  endif
	 enddo ploop
!	else
!	 write(*,*)'ERROR(random_permutation): negative or zero number of items: ',ni
!	 stop
	endif
	return
	end subroutine random_permutation
!------------------------------------------------------------
	subroutine random_permutation_int8(ni,trn,no_trivial)
!This subroutine returns a random permutation of NI items [1..NI].
!INPUT:
! - ni - number of items, range [1:ni], if ni<=0 nothing will be done;
!OUTPUT:
! - trn(0:ni) - generated permutation (trn(0) - sign of the permutation);
	implicit none
	integer(8), intent(in):: ni
	integer(8), intent(out):: trn(0:ni)
	logical, intent(in), optional:: no_trivial
!-------------------------------------------------
	integer(8), parameter:: random_chunk=2**10 !size of the chunk of random numbers generated in one call
	integer(8), parameter:: num_repeats=5_8    !the bigger the number, the better the generator quality (more expensive)
!-------------------------------------------------
	integer(8) i,j,k,l,m,n,nr,ierr
	real(8):: ra(1:random_chunk)

	if(ni.gt.0_8) then
	 trn(0)=+1_8; do i=1_8,ni; trn(i)=i; enddo !initial permutation
	 ploop: do
	  do nr=1_8,num_repeats
	   do i=1_8,ni,random_chunk
	    l=min(i+random_chunk-1_8,ni)-i+1_8
	    call random_number(ra(1_8:l))
	    ra(1_8:l)=ra(1_8:l)*2d0
	    do k=1_8,l
	     if(ra(k).lt.1d0) then
	      j=i+k-1_8
	      n=int(ra(k)*dble(ni),8)+1_8; if(n.gt.ni) n=ni
	      if(n.ne.j) then; m=trn(j); trn(j)=trn(n); trn(n)=m; trn(0)=-trn(0); endif
	     endif
	    enddo
	   enddo
	  enddo
	  if(present(no_trivial)) then
	   if(no_trivial) then
	    if(.not.perm_trivial_int8(ni,trn)) exit ploop
	   else
	    exit ploop
	   endif
	  else
	   exit ploop
	  endif
	 enddo ploop
!	else
!	 write(*,*)'ERROR(random_permutation_int8): negative or zero number of items: ',ni
!	 stop
	endif
	return
	end subroutine random_permutation_int8
!-----------------------------------------------------------
	subroutine random_composition(ordered,irange,ni,trn)
!This subroutine returns a random sequence of ni natural numbers from the range [1..ni] without repeats.
!INPUT:
! - ordered - if .true. the sequence will be ordered;
! - irange - range of natural numbers to be used: [1..irange];
! - ni - number of elements in the sequence (length of the sequence);
!OUTPUT:
! - trn(0:ni) - the sequence generated, trn(0) is the sign of the permutation if unordered.
	implicit none
	logical, intent(in):: ordered
	integer, intent(in):: irange
	integer, intent(in):: ni
	integer, intent(out):: trn(0:ni)
	integer, parameter:: rnd_chunk=2**12
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf,ierr
	integer, allocatable:: prm(:)
	real(8) rnd_buf(1:rnd_chunk),rn,accept_thresh

	accept_thresh(k,l)=dble(k-l)/dble(k) !k>=l, k!=0: k - amount of objects left for selection; l - number of items to select.

	if(ni.ge.0.and.irange.ge.1.and.irange.ge.ni) then
	 if(ni.gt.0) then
	  if(ordered) then
	   k=irange; l=ni; k0=0; n=0
	   do i=1,irange
	    if(k0.eq.0) then; k0=min(rnd_chunk,k); call random_number(rnd_buf(1:k0)); endif
	    if(rnd_buf(k0).ge.accept_thresh(k,l)) then
	     n=n+1; trn(n)=i
	     l=l-1; if(l.eq.0) exit
	    endif
	    k=k-1; k0=k0-1
	   enddo
	   if(n.eq.ni) then
	    trn(0)=+1
	   else
	    write(*,*)'ERROR(combinatoric:random_composition): trap: invalid number of items: ',n,ni,irange,ordered
	    stop
	   endif
	  else
	   allocate(prm(0:ni),STAT=ierr)
	   if(ierr.ne.0) then; write(*,*)'ERROR(combinatoric:random_composition): allocation failed!'; stop; endif
	   prm(0)=+1; do i=1,ni; prm(i)=i; enddo
	   call random_permutation(ni,prm)
	   k=irange; l=ni; k0=0; n=0
	   do i=1,irange
	    if(k0.eq.0) then; k0=min(rnd_chunk,k); call random_number(rnd_buf(1:k0)); endif
	    if(rnd_buf(k0).ge.accept_thresh(k,l)) then
	     n=n+1; trn(prm(n))=i
	     l=l-1; if(l.eq.0) exit
	    endif
	    k=k-1; k0=k0-1
	   enddo
	   deallocate(prm)
	   if(n.eq.ni) then
	    trn(0)=prm(0)
	   else
	    write(*,*)'ERROR(combinatoric:random_composition): trap: invalid number of items: ',n,ni,irange,ordered
	    stop
	   endif
	  endif
	 endif
	else
	 write(*,*)'ERROR(combinatoric:random_composition): incompatible or invalid arguments: ',ni,irange
	 stop
	endif
	return
	end subroutine random_composition
!----------------------------------------
	subroutine merge_sort_int(ni,trn)
!This subroutine sorts an array of NI items in a non-descending order.
!The algorithm was suggested by Johann von Neumann.
!INPUT:
! - ni - number of items;
! - trn(1:ni) - items (an array of arbitrary integers);
!OUTPUT:
! - trn(0:ni) - sorted items and sign in trn(0).
!NOTES:
! - In order to accelerate the procedure use flip/flop for TRN(:)||PRM(:).
	implicit none
	integer, intent(in):: ni
	integer, intent(inout):: trn(0:ni)
	integer i,j,k,l,m,n,k0,k1,k2,k3,k4,k5,k6,k7,ks,kf,ierr
	integer, allocatable:: prm(:)

	trn(0)=+1
	if(ni.gt.1) then
	 allocate(prm(1:ni))
	 n=1
	 do while(n.lt.ni)
	  m=n*2
	  do i=1,ni,m
	   k1=i; k2=i+n
	   if(k2.gt.ni) then
	    k2=ni+1; k3=0; k4=0 !no right block, only left block
	   else
	    k3=i+n; k4=min(ni+1,i+m) !right block present
	   endif
	   kf=min(ni+1,i+m)-i; l=0
	   do while(l.lt.kf)
	    if(k3.ge.k4) then !right block is over
	     prm(i+l:i+kf-1)=trn(k1:k2-1); l=kf
	    elseif(k1.ge.k2) then !left block is over
	     prm(i+l:i+kf-1)=trn(k3:k4-1); l=kf
	    else
	     if(trn(k1)-trn(k3).gt.0) then
	      prm(i+l)=trn(k3); k3=k3+1; trn(0)=(1-2*mod(k2-k1,2))*trn(0)
	     else
	      prm(i+l)=trn(k1); k1=k1+1
	     endif
	     l=l+1
	    endif
	   enddo
	  enddo
	  trn(1:ni)=prm(1:ni)
	  n=m
	 enddo
	 deallocate(prm)
	endif
	return
	end subroutine merge_sort_int
!------------------------------------------------
	subroutine merge_sort_key_int(ni,key,trn)
!This subroutine sorts an array of NI items in a non-descending order according to their keys.
!The algorithm is due to Johann von Neumann.
!INPUT:
! - ni - number of items;
! - key(1:ni) - item keys (retrieved by old item numbers!): arbitrary integers;
! - trn(0:ni) - initial permutation of ni items (a sequence of old numbers), trn(0) is the initial sign;
!OUTPUT:
! - trn(0:ni) - sorted permutation (new sequence of old numbers) of ni items (according to their keys), the sign is in trn(0).
	implicit none
	integer, intent(in):: ni,key(1:ni)
	integer, intent(inout):: trn(0:ni)
	integer, parameter:: max_in_mem=1024
	integer i,j,k,l,m,n,k0,k1,k2,k3,k4,k5,k6,k7,ks,kf,ierr
	integer, target:: prms(1:max_in_mem)
	integer, allocatable,target:: prma(:)
	integer, pointer:: prm(:)

	if(ni.gt.1) then
	 if(ni.le.max_in_mem) then; prm=>prms; else; allocate(prma(1:ni)); prm=>prma; endif
	 n=1
	 do while(n.lt.ni)
	  m=n*2
	  do i=1,ni,m
	   k1=i; k2=i+n
	   if(k2.gt.ni) then
	    k2=ni+1; k3=0; k4=0 !no right block, only left block
	   else
	    k3=i+n; k4=min(ni+1,i+m) !right block present
	   endif
	   kf=min(ni+1,i+m)-i; l=0
	   do while(l.lt.kf)
	    if(k3.ge.k4) then !right block is over
	     prm(i+l:i+kf-1)=trn(k1:k2-1); l=kf
	    elseif(k1.ge.k2) then !left block is over
	     prm(i+l:i+kf-1)=trn(k3:k4-1); l=kf
	    else
	     if(key(trn(k1))-key(trn(k3)).gt.0) then
	      prm(i+l)=trn(k3); k3=k3+1; trn(0)=(1-2*mod(k2-k1,2))*trn(0)
	     else
	      prm(i+l)=trn(k1); k1=k1+1
	     endif
	     l=l+1
	    endif
	   enddo
	  enddo
	  trn(1:ni)=prm(1:ni)
	  n=m
	 enddo
	 nullify(prm); if(ni.gt.max_in_mem) deallocate(prma)
	endif
	return
	end subroutine merge_sort_key_int
!-----------------------------------------
	subroutine merge_sort_int8(ni,trn)
!This subroutine sorts an array of NI items in a non-descending order.
!The algorithm was suggested by Johann von Neumann.
!INPUT:
! - ni - number of items;
! - trn(1:ni) - items (arbitrary integer*8 numbers);
!OUTPUT:
! - trn(0:ni) - sorted items and sign in trn(0).
!NOTES:
! - In order to accelerate the procedure use flip/flop for TRN(:)||PRM(:).
	implicit none
	integer(8), intent(in):: ni
	integer(8), intent(inout):: trn(0:ni)
	integer(8) i,j,k,l,m,n,k0,k1,k2,k3,k4,k5,k6,k7,ks,kf,ierr
	integer(8), allocatable:: prm(:)

	trn(0)=+1_8
	if(ni.gt.1_8) then
	 allocate(prm(1_8:ni))
	 n=1_8
	 do while(n.lt.ni)
	  m=n*2_8
	  do i=1_8,ni,m
	   k1=i; k2=i+n
	   if(k2.gt.ni) then
	    k2=ni+1_8; k3=0_8; k4=0_8 !no right block, only left block
	   else
	    k3=i+n; k4=min(ni+1_8,i+m) !right block present
	   endif
	   kf=min(ni+1_8,i+m)-i; l=0_8
	   do while(l.lt.kf)
	    if(k3.ge.k4) then !right block is over
	     prm(i+l:i+kf-1_8)=trn(k1:k2-1_8); l=kf
	    elseif(k1.ge.k2) then !left block is over
	     prm(i+l:i+kf-1_8)=trn(k3:k4-1_8); l=kf
	    else
	     if(trn(k1)-trn(k3).gt.0_8) then
	      prm(i+l)=trn(k3); k3=k3+1_8; trn(0_8)=(1_8-2_8*mod(k2-k1,2_8))*trn(0_8)
	     else
	      prm(i+l)=trn(k1); k1=k1+1_8
	     endif
	     l=l+1_8
	    endif
	   enddo
	  enddo
	  trn(1_8:ni)=prm(1_8:ni)
	  n=m
	 enddo
	 deallocate(prm)
	endif
	return
	end subroutine merge_sort_int8
!-------------------------------------------------
	subroutine merge_sort_key_int8(ni,key,trn)
!This subroutine sorts an array of NI items in a non-descending order according to their keys.
!The algorithm is due to Johann von Neumann.
!INPUT:
! - ni - number of items;
! - key(1:ni) - item keys (retrieved by old item numbers!): arbitrary integer*8 numbers;
! - trn(0:ni) - initial permutation of ni items (a sequence of old numbers), trn(0) is the initial sign;
!OUTPUT:
! - trn(0:ni) - sorted permutation of ni items (according to their keys), the sign is in trn(0).
	implicit none
	integer(8), intent(in):: ni,key(1:ni)
	integer(8), intent(inout):: trn(0:ni)
	integer(8), parameter:: max_in_mem=1024
	integer(8) i,j,k,l,m,n,k0,k1,k2,k3,k4,k5,k6,k7,ks,kf,ierr
	integer(8), target:: prms(1:max_in_mem)
	integer(8), allocatable,target:: prma(:)
	integer(8), pointer:: prm(:)

	if(ni.gt.1_8) then
	 if(ni.le.max_in_mem) then; prm=>prms; else; allocate(prma(1_8:ni)); prm=>prma; endif
	 n=1_8
	 do while(n.lt.ni)
	  m=n*2_8
	  do i=1_8,ni,m
	   k1=i; k2=i+n
	   if(k2.gt.ni) then
	    k2=ni+1_8; k3=0_8; k4=0_8 !no right block, only left block
	   else
	    k3=i+n; k4=min(ni+1_8,i+m) !right block present
	   endif
	   kf=min(ni+1_8,i+m)-i; l=0_8
	   do while(l.lt.kf)
	    if(k3.ge.k4) then !right block is over
	     prm(i+l:i+kf-1_8)=trn(k1:k2-1_8); l=kf
	    elseif(k1.ge.k2) then !left block is over
	     prm(i+l:i+kf-1_8)=trn(k3:k4-1_8); l=kf
	    else
	     if(key(trn(k1))-key(trn(k3)).gt.0_8) then
	      prm(i+l)=trn(k3); k3=k3+1_8; trn(0_8)=(1_8-2_8*mod(k2-k1,2_8))*trn(0_8)
	     else
	      prm(i+l)=trn(k1); k1=k1+1_8
	     endif
	     l=l+1_8
	    endif
	   enddo
	  enddo
	  trn(1_8:ni)=prm(1_8:ni)
	  n=m
	 enddo
	 nullify(prm); if(ni.gt.max_in_mem) deallocate(prma)
	endif
	return
	end subroutine merge_sort_key_int8
!------------------------------------------
	subroutine merge_sort_real8(ni,trn)
!This subroutine sorts an array of NI items in a non-descending order.
!The algorithm was suggested by Johann von Neumann.
!INPUT:
! - ni - number of items;
! - trn(1:ni) - items (arbitrary real*8 numbers);
!OUTPUT:
! - trn(0:ni) - sorted items and sign in trn(0).
	implicit none
	integer, intent(in):: ni
	real(8), intent(inout):: trn(0:ni)
	integer i,j,k,l,m,n,k0,k1,k2,k3,k4,k5,k6,k7,ks,kf,ierr
	real(8), parameter:: ds(0:1)=(/+1d0,-1d0/)
	real(8), allocatable:: prm(:)

	trn(0)=+1d0
	if(ni.gt.1) then
	 allocate(prm(1:ni))
	 n=1
	 do while(n.lt.ni)
	  m=n*2
	  do i=1,ni,m
	   k1=i; k2=i+n
	   if(k2.gt.ni) then
	    k2=ni+1; k3=0; k4=0 !no right block, only left block
	   else
	    k3=i+n; k4=min(ni+1,i+m) !right block present
	   endif
	   kf=min(ni+1,i+m)-i; l=0
	   do while(l.lt.kf)
	    if(k3.ge.k4) then !right block is over
	     prm(i+l:i+kf-1)=trn(k1:k2-1); l=kf
	    elseif(k1.ge.k2) then !left block is over
	     prm(i+l:i+kf-1)=trn(k3:k4-1); l=kf
	    else
	     if(trn(k1)-trn(k3).gt.0d0) then
	      prm(i+l)=trn(k3); k3=k3+1; trn(0)=ds(mod(k2-k1,2))*trn(0)
	     else
	      prm(i+l)=trn(k1); k1=k1+1
	     endif
	     l=l+1
	    endif
	   enddo
	  enddo
	  trn(1:ni)=prm(1:ni)
	  n=m
	 enddo
	 deallocate(prm)
	endif
	return
	end subroutine merge_sort_real8
!--------------------------------------------------
	subroutine merge_sort_key_real8(ni,key,trn)
!This subroutine sorts an array of NI items in a non-descending order according to their keys.
!The algorithm is due to Johann von Neumann.
!INPUT:
! - ni - number of items;
! - key(1:ni) - item keys (retrieved by old item numbers!): arbitrary real*8 numbers;
! - trn(0:ni) - initial permutation of ni items (a sequence of old numbers), trn(0) is the initial sign;
!OUTPUT:
! - trn(0:ni) - sorted permutation of ni items (according to their keys), the sign is in trn(0).
	implicit none
	integer, intent(in):: ni
	real(8), intent(in):: key(1:ni)
	integer, intent(inout):: trn(0:ni)
	integer, parameter:: max_in_mem=1024
	integer i,j,k,l,m,n,k0,k1,k2,k3,k4,k5,k6,k7,ks,kf,ierr
	integer, target:: prms(1:max_in_mem)
	integer, allocatable,target:: prma(:)
	integer, pointer:: prm(:)

	if(ni.gt.1) then
	 if(ni.le.max_in_mem) then; prm=>prms; else; allocate(prma(1:ni)); prm=>prma; endif
	 n=1
	 do while(n.lt.ni)
	  m=n*2
	  do i=1,ni,m
	   k1=i; k2=i+n
	   if(k2.gt.ni) then
	    k2=ni+1; k3=0; k4=0 !no right block, only left block
	   else
	    k3=i+n; k4=min(ni+1,i+m) !right block present
	   endif
	   kf=min(ni+1,i+m)-i; l=0
	   do while(l.lt.kf)
	    if(k3.ge.k4) then !right block is over
	     prm(i+l:i+kf-1)=trn(k1:k2-1); l=kf
	    elseif(k1.ge.k2) then !left block is over
	     prm(i+l:i+kf-1)=trn(k3:k4-1); l=kf
	    else
	     if(key(trn(k1))-key(trn(k3)).gt.0d0) then
	      prm(i+l)=trn(k3); k3=k3+1; trn(0)=(1-2*mod(k2-k1,2))*trn(0)
	     else
	      prm(i+l)=trn(k1); k1=k1+1
	     endif
	     l=l+1
	    endif
	   enddo
	  enddo
	  trn(1:ni)=prm(1:ni)
	  n=m
	 enddo
	 nullify(prm); if(ni.gt.max_in_mem) deallocate(prma)
	endif
	return
	end subroutine merge_sort_key_real8
!-------------------------------------------
	subroutine merge_sort_cmplx8(ni,trn)
!This subroutine sorts an array of NI items in a non-descending order.
!The algorithm was suggested by Johann von Neumann.
!COMPLEX(8) comparison (non-standard): (x1,y1)>(x2,y2) iff [[x1>x2].OR.[x1=x2.AND.y1>y2]]
!INPUT:
! - ni - number of items;
! - trn(1:ni) - items (arbitrary complex*8 numbers);
!OUTPUT:
! - trn(0:ni) - sorted items and sign in trn(0).
	implicit none
	integer, intent(in):: ni
	complex(8), intent(inout):: trn(0:ni)
	integer i,j,k,l,m,n,k0,k1,k2,k3,k4,k5,k6,k7,ks,kf,ierr
	real(8), parameter:: ds(0:1)=(/+1d0,-1d0/)
	complex(8), allocatable:: prm(:)
	real(8) sgn

	sgn=+1d0
	if(ni.gt.1) then
	 allocate(prm(1:ni))
	 n=1
	 do while(n.lt.ni)
	  m=n*2
	  do i=1,ni,m
	   k1=i; k2=i+n
	   if(k2.gt.ni) then
	    k2=ni+1; k3=0; k4=0 !no right block, only left block
	   else
	    k3=i+n; k4=min(ni+1,i+m) !right block present
	   endif
	   kf=min(ni+1,i+m)-i; l=0
	   do while(l.lt.kf)
	    if(k3.ge.k4) then !right block is over
	     prm(i+l:i+kf-1)=trn(k1:k2-1); l=kf
	    elseif(k1.ge.k2) then !left block is over
	     prm(i+l:i+kf-1)=trn(k3:k4-1); l=kf
	    else
	     if(dble(trn(k1))-dble(trn(k3)).gt.0d0.or. &
	        (dble(trn(k1)).eq.dble(trn(k3)).and.dimag(trn(k1))-dimag(trn(k3)).gt.0d0)) then
	      prm(i+l)=trn(k3); k3=k3+1; sgn=ds(mod(k2-k1,2))*sgn
	     else
	      prm(i+l)=trn(k1); k1=k1+1
	     endif
	     l=l+1
	    endif
	   enddo
	  enddo
	  trn(1:ni)=prm(1:ni)
	  n=m
	 enddo
	 deallocate(prm)
	endif
	trn(0)=dcmplx(sgn,0d0)
	return
	end subroutine merge_sort_cmplx8
!---------------------------------------------------
	subroutine merge_sort_key_cmplx8(ni,key,trn)
!This subroutine sorts an array of NI items in a non-descending order according to their keys.
!The algorithm is due to Johann von Neumann.
!COMPLEX(8) comparison (non-standard): (x1,y1)>(x2,y2) iff [[x1>x2].OR.[x1=x2.AND.y1>y2]]
!INPUT:
! - ni - number of items;
! - key(1:ni) - item keys (retrieved by old item numbers!): arbitrary complex*8 numbers;
! - trn(0:ni) - initial permutation of ni items (a sequence of old numbers), trn(0) is the initial sign;
!OUTPUT:
! - trn(0:ni) - sorted permutation of ni items (according to their keys), the sign is in trn(0).
	implicit none
	integer, intent(in):: ni
	complex(8), intent(in):: key(1:ni)
	integer, intent(inout):: trn(0:ni)
	integer, parameter:: max_in_mem=1024
	integer i,j,k,l,m,n,k0,k1,k2,k3,k4,k5,k6,k7,ks,kf,ierr
	integer, target:: prms(1:max_in_mem)
	integer, allocatable,target:: prma(:)
	integer, pointer:: prm(:)

	if(ni.gt.1) then
	 if(ni.le.max_in_mem) then; prm=>prms; else; allocate(prma(1:ni)); prm=>prma; endif
	 n=1
	 do while(n.lt.ni)
	  m=n*2
	  do i=1,ni,m
	   k1=i; k2=i+n
	   if(k2.gt.ni) then
	    k2=ni+1; k3=0; k4=0 !no right block, only left block
	   else
	    k3=i+n; k4=min(ni+1,i+m) !right block present
	   endif
	   kf=min(ni+1,i+m)-i; l=0
	   do while(l.lt.kf)
	    if(k3.ge.k4) then !right block is over
	     prm(i+l:i+kf-1)=trn(k1:k2-1); l=kf
	    elseif(k1.ge.k2) then !left block is over
	     prm(i+l:i+kf-1)=trn(k3:k4-1); l=kf
	    else
	     if(dble(key(trn(k1)))-dble(key(trn(k3))).gt.0d0.or. &
	        (dble(key(trn(k1))).eq.dble(key(trn(k3))).and.dimag(key(trn(k1)))-dimag(key(trn(k3))).gt.0d0)) then
	      prm(i+l)=trn(k3); k3=k3+1; trn(0)=(1-2*mod(k2-k1,2))*trn(0)
	     else
	      prm(i+l)=trn(k1); k1=k1+1
	     endif
	     l=l+1
	    endif
	   enddo
	  enddo
	  trn(1:ni)=prm(1:ni)
	  n=m
	 enddo
	 nullify(prm); if(ni.gt.max_in_mem) deallocate(prma)
	endif
	return
	end subroutine merge_sort_key_cmplx8
!-------------------------------------------
	subroutine sort_slots(ni,keys,trn)
!This subroutine sorts an array of NI items of type SLOT putting them in a non-descending order
!by their key-values in accordance with a special relation of order.
!Stable MERGE_SORT algorithm of Johann von Neumann is used.
!INPUT:
! - ni - number of items;
! - keys(1:ni) - item keys: given an original item ID, KEYS(ID) returns the item key of type SLOT;
! - trn(0:ni) - an input permutation of items given as a sequence of their original numbers (IDs), where TRN(0) is the current sign;
!OUTPUT:
! - trn(0:ni) - new permutation of items corresponding to sorted keys; sign in TRN(0).
!NOTES:
! - KEYS(:) does not change since its argument is the original item number (item ID).
	implicit none
	integer, intent(in):: ni
	type(slot), intent(in):: keys(1:ni)
	integer, intent(inout):: trn(0:ni)
!--------------------------------------------
	integer, parameter:: max_dim_mem=1024 !maximum number of items that does not need an additional memory allocation
!--------------------------------------------
	integer i,j,k,l,m,n,k0,k1,k2,k3,k4,k5,k6,k7,ks,kf,ierr
	integer, target:: prms(1:max_dim_mem)
	integer, allocatable, target:: prma(:)
	integer, pointer:: prm(:)

!	write(*,'(64(1x,i2))') ni,trn(0:ni) !debug
	if(ni.gt.1) then
	 if(ni.gt.max_dim_mem) then; allocate(prma(1:ni)); prm=>prma; else; prm=>prms; endif
	 n=1
	 do while(n.lt.ni)
	  m=n*2
	  do i=1,ni,m
	   k1=i; k2=i+n
	   if(k2.gt.ni) then
	    k2=ni+1; k3=0; k4=0 !no right block, only left block
	   else
	    k3=i+n; k4=min(ni+1,i+m) !right block present
	   endif
	   kf=min(ni+1,i+m)-i; l=0
	   do while(l.lt.kf)
	    if(k3.ge.k4) then !right block is over
	     prm(i+l:i+kf-1)=trn(k1:k2-1); l=kf
	    elseif(k1.ge.k2) then !left block is over
	     prm(i+l:i+kf-1)=trn(k3:k4-1); l=kf
	    else !both blocks present
	     if(slot_gt(keys(trn(k1)),keys(trn(k3)))) then
	      prm(i+l)=trn(k3); k3=k3+1; trn(0)=(1-2*mod(k2-k1,2))*trn(0)
	     else
	      prm(i+l)=trn(k1); k1=k1+1
	     endif
	     l=l+1
	    endif
	   enddo
	  enddo
	  trn(1:ni)=prm(1:ni)
	  n=m
	 enddo
	 nullify(prm); if(ni.gt.max_dim_mem) deallocate(prma)
	endif
	return

	contains

	logical function slot_gt(sl1,sl2) !.TRUE. if SLOT1>SLOT2, .FALSE. otherwise
	 implicit none
	 type(slot), intent(in):: sl1,sl2
	 if(sl1%vert_attr.gt.sl2%vert_attr) then
	  slot_gt=.true.
	 elseif(sl1%vert_attr.lt.sl2%vert_attr) then
	  slot_gt=.false.
	 else
	  if(dble(sl1%conn_typ).gt.dble(sl2%conn_typ)) then
	   slot_gt=.true.
	  elseif(dble(sl1%conn_typ).lt.dble(sl2%conn_typ)) then
	   slot_gt=.false.
	  else
	   if(dimag(sl1%conn_typ).gt.dimag(sl2%conn_typ)) then
	    slot_gt=.true.
	   elseif(dimag(sl1%conn_typ).lt.dimag(sl2%conn_typ)) then
	    slot_gt=.false.
	   else
	    slot_gt=.false.
	   endif
	  endif
	 endif
	 return
	end function slot_gt

	end subroutine sort_slots
!------------------------------------------
	subroutine nullify_sparse_graph(gr)
!This subroutine nullifies an object of type(graph_sparse), gr.
	implicit none
	integer i,j,k,l,m,n,k0,k1,k2,k3,ierr
	type(graph_sparse), intent(inout):: gr

	if(allocated(gr%vertices)) deallocate(gr%vertices)
	if(allocated(gr%descr)) deallocate(gr%descr)
	gr%gcard=0; gr%graph_typ=0; gr%descr_allocated=.false.
	return
	end subroutine nullify_sparse_graph
!--------------------------------------------------
	subroutine copy_sparse_graph(clean,grs,grd)
!This subroutine copies an object of type(graph_sparse), grs, into another object, grd, of the same type.
!If clean=.true. grd will be cleaned before copying (safe mode, slower).
!If clean=.false. the existing allocated structure of grd will be used (faster, user responsibility to avoid segmentation faults).
	implicit none
	integer i,j,k,l,m,n,k0,k1,k2,k3,ierr
	logical, intent(in):: clean
	type(graph_sparse), intent(in):: grs
	type(graph_sparse), intent(inout):: grd

	ierr=0
	if(clean) call nullify_sparse_graph(grd)
	n=grs%gcard
	grd%gcard=grs%gcard; grd%graph_typ=grs%graph_typ; grd%descr_allocated=grs%descr_allocated
	if(n.gt.0) then
	 if(.not.allocated(grd%vertices)) then
	  allocate(grd%vertices(1:n),STAT=ierr)
	  if(ierr.ne.0) then; write(*,*)'ERROR(combinatoric:copy_sparse_graph): allocation 1 failed!'; stop; endif
	 endif
	 do i=1,n
	  grd%vertices(i)%color=grs%vertices(i)%color
	  m=grs%vertices(i)%conn_num; grd%vertices(i)%conn_num=m
	  if(m.gt.0) then
	   if(.not.allocated(grd%vertices(i)%slots)) then
	    allocate(grd%vertices(i)%slots(1:m),STAT=ierr)
	    if(ierr.ne.0) then; write(*,*)'ERROR(combinatoric:copy_sparse_graph): allocation 2 failed!'; stop; endif
	   endif
	   grd%vertices(i)%slots(1:m)=grs%vertices(i)%slots(1:m)
	  endif
	 enddo
	 if(grd%descr_allocated) then
	  if(.not.allocated(grd%descr)) then
	   allocate(grd%descr(1:n),STAT=ierr)
	   if(ierr.ne.0) then; write(*,*)'ERROR(combinatoric:copy_sparse_graph): allocation 3 failed!'; stop; endif
	  endif
	  do i=1,n
	   grd%descr(i)%vertex_class=grs%descr(i)%vertex_class
	   m=grs%descr(i)%conn_num; grd%descr(i)%conn_num=m
	   if(m.gt.0) then
	    if(.not.allocated(grd%descr(i)%slots)) then
	     allocate(grd%descr(i)%slots(1:m),STAT=ierr)
	     if(ierr.ne.0) then; write(*,*)'ERROR(combinatoric:copy_sparse_graph): allocation 4 failed!'; stop; endif
	    endif
	    grd%descr(i)%slots(1:m)=grs%descr(i)%slots(1:m)
	   endif
	  enddo
	 endif
	endif
	return
	end subroutine copy_sparse_graph
!----------------------------------------------------
	integer function compare_descr_slots(sl1,sl2)
!This function compares two objects of type slot, sl1 and sl2.
!Result: -1:sl1<sl2; 0:sl1=sl2; +1:sl1>sl2.
!NOTE: A rigorous (digit-wise) comparison of real numbers is employed!
	implicit none
	type(slot), intent(in):: sl1,sl2
	integer i,j,k,l,m,k0,k1,k2,k3,ks,kf,ierr

	if(sl1%vert_attr.gt.sl2%vert_attr) then
	 compare_descr_slots=+1
	elseif(sl1%vert_attr.lt.sl2%vert_attr) then
	 compare_descr_slots=-1
	else
	 if(dble(sl1%conn_typ).gt.dble(sl2%conn_typ)) then
	  compare_descr_slots=+1
	 elseif(dble(sl1%conn_typ).lt.dble(sl2%conn_typ)) then
	  compare_descr_slots=-1
	 else
	  if(dimag(sl1%conn_typ).gt.dimag(sl2%conn_typ)) then
	   compare_descr_slots=+1
	  elseif(dimag(sl1%conn_typ).lt.dimag(sl2%conn_typ)) then
	   compare_descr_slots=-1
	  else
	   compare_descr_slots=0
	  endif
	 endif
	endif
	return
	end function compare_descr_slots
!---------------------------------------------
	subroutine sort_vertex_connections(gr)
!This subroutine rearranges vertex slots in an ordered form.
	implicit none
	type(graph_sparse), intent(inout):: gr
	integer, parameter:: max_in_mem=1023
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf,ierr
	type(slot), target:: keys(1:max_in_mem)
	type(slot), allocatable, target:: keya(:)
	type(slot), pointer:: key(:)
	integer, target:: trns(0:max_in_mem)
	integer, allocatable, target:: trna(:)
	integer, pointer:: trn(:)

	n=gr%gcard
	if(n.gt.0) then
	 if(n.gt.max_in_mem) then
	  allocate(keya(1:max_in_mem),trna(0:max_in_mem))
	  key=>keya; trn=>trna
	 else
	  key=>keys; trn=>trns
	 endif
	 do i=1,n
	  m=gr%vertices(i)%conn_num
	  if(m.gt.1.and.m.le.n) then
	   trn(0)=+1; do j=1,m; trn(j)=j; enddo !initial permutation
	   key(1:m)=gr%vertices(i)%slots(1:m)
	   call sort_slots(m,key,trn)
	   gr%vertices(i)%slots(1:m)=key(trn(1:m))
	  elseif(m.lt.0.or.m.gt.n) then
	   write(*,*)'ERROR(sort_vertex_connections): invalid graph: abnormal number of vertex connections: ',m,n !stop
	   stop
	  endif
	 enddo
	endif
	if(associated(key)) nullify(key); if(associated(trn)) nullify(trn)
	if(allocated(keya)) deallocate(keya); if(allocated(trna)) deallocate(trna)
	return
	end subroutine sort_vertex_connections
!--------------------------------------------------
	subroutine graph_standard(gr,trn,auto,ierr)
!This subroutine returns one of the standard permutations of an arbitrary (colored, possibly directed) (multi-)graph
!represented by a hermitian complex(8) adjacency matrix. The set of all standard permutations
!constitutes the automorphism group of the graph (the size of the group is returned in AUTO).
!Class Reduction Algorithm is used [D.I.Lyakh, 2011]. Algorithm's worst case is O(N^5*logN*logM),
!where M is the number of connection types (i.e., the number of distinct edge types, 1 for simple graphs).
!The algorithm is based solely on comparisons (hermiticity of the adjacency matrix must be exact!!!).
!INPUT:
! - gr - graph adjacency matrix given in a sparse (packed) form (hermiticity must be strictly imposed!!!);
!OUTPUT:
! - trn(0:*) - standard permutation [1..*] (N2O), where trn(0) is the sign of the permutation;
! - auto - inversed automorphic factor (dimension of the graph automorpishm group); zero means INTEGER*4 overflow;
! - gr with ordered vertex slots;
! - ierr - error code (0 - success).
!NOTES:
! - The vertex numeration begins with 1.
! - The initial size of the largest class (mcl) is not reduced (and not used) in the algorithm.
	implicit none
	type(graph_sparse), intent(inout):: gr
	integer, intent(out):: trn(0:*),auto,ierr
!--------------------------------------------
	integer, parameter:: max_dim_mem=1023 !maximum matrix dimension that does not need an additional memory allocation
	logical, parameter:: no_clean=.false.
	logical, parameter:: clean=.true.
!--------------------------------------------
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf
	integer, target:: iis(0:max_dim_mem)
	integer, allocatable, target:: iia(:)
	integer, pointer:: ii(:)
	integer, target:: vcs(0:max_dim_mem),save_vcs(0:max_dim_mem),save_trns(0:max_dim_mem)
	integer, target:: trnps(0:max_dim_mem),vcls(0:max_dim_mem)
	integer, allocatable, target:: vca(:),save_vca(:),save_trna(:),trnpa(:),vcla(:)
	integer, pointer:: vc(:),save_vc(:),save_trn(:),trnp(:),vcl(:)
	type(slot), target:: sls(1:max_dim_mem)
	type(slot), allocatable, target:: sla(:)
	type(slot), pointer:: sl(:)
	integer ni,ncl,mcl,nclo
	integer save_ncl,save_mcl,save_nclo,checking_stability,class_checked,nclp
	type(graph_sparse) save_gr,grp
	logical add_alloc,context_restored,perturbed

	ierr=0
!	write(*,'("Entered GRAPH_STANDARD:")') !debug
	add_alloc=.false.
	if(gr%gcard.gt.0) then
	 ni=gr%gcard
	 if(ni.gt.max_dim_mem) then
	  allocate(vca(0:ni),STAT=ierr); if(ierr.ne.0) then; ierr=1; goto 999; endif
	  allocate(sla(1:ni),STAT=ierr); if(ierr.ne.0) then; ierr=2; goto 999; endif
	  allocate(iia(0:ni),STAT=ierr); if(ierr.ne.0) then; ierr=3; goto 999; endif
	  allocate(vcla(0:ni),STAT=ierr); if(ierr.ne.0) then; ierr=9; goto 999; endif
	  ii=>iia; vc=>vca; sl=>sla; vcl=>vcla
	 else
	  ii=>iis; vc=>vcs; sl=>sls; vcl=>vcls
	 endif
!Sort vertex connections:
	 call sort_vertex_connections(gr) !order vertex slots
!Sort vertices (diagonal elements of the adjacency matrix) according to (a) their color and (b) amount of connections:
	 auto=1; do i=1,ni; trn(i)=i; enddo; trn(0)=+1 !initial permutation and automorphic factor
	 do i=1,ni
	  sl(i)%vert_attr=0; sl(i)%conn_typ=dcmplx(gr%vertices(i)%color,dble(gr%vertices(i)%conn_num))
	 enddo
	 call sort_slots(ni,sl,trn) !vertices are sorted with respect to (a) their colors and (b) amount of connections
	 ncl=1; mcl=1; m=1; vc(trn(1))=1
	 do i=2,ni
	  if(slot_gt(sl(trn(i)),sl(trn(i-1)))) then
	   ncl=ncl+1; mcl=max(mcl,m); m=1 !new class
	  else
	   m=m+1
	  endif
	  vc(trn(i))=ncl
	 enddo
	 mcl=max(mcl,m); vc(0)=ncl
!Record the initial (preliminary) vcl (vertex classes from the first stable irreducible classification):
	 vcl(0:ni)=vc(0:ni)
!Initialize vertex descriptors:
	 if(allocated(gr%descr)) deallocate(gr%descr)
	 allocate(gr%descr(1:ni),STAT=ierr); if(ierr.ne.0) then; ierr=4; goto 999; endif
	 do i=1,ni !loop over the graph vertices by their IDs (not current positions)
	  gr%descr(i)%vertex_class=vc(i)
	  m=gr%vertices(i)%conn_num
	  gr%descr(i)%conn_num=m
	  if(m.gt.0) then
	   allocate(gr%descr(i)%slots(1:m),STAT=ierr); if(ierr.ne.0) then; ierr=5; goto 999; endif
	  endif
	 enddo
	 gr%descr_allocated=.true.
!Description of the data:
! ni: number of vertices
! ncl: number of classes (=vc(0))
! mcl: the size of the largest class of the initial classification (not used further)
! vc(vertex_ID): vertex class [1..ncl]: all vertices within a class have the same color and the same number of connections (in particular).
! trn(vertex_POSITION): vertex ID [1..ni]; trn(0): sign of the current permutation of vertices.
!	 write(*,'(" Position: ",i2,2x,32(1x,i2))') (j,j=0,ni) !debug
!	 write(*,'(" Item IDs: ",i2,2x,32(1x,i2))') trn(0:ni) !debug
!	 write(*,'(" Classes : ",i2,2x,32(1x,i2))') vc(0),vc(trn(1:ni)); write(*,*) !debug
!Class Reduction Algorithm (D.I.Lyakh, 2011):
	 if(ncl.lt.ni) then
 !Iterations:
	  checking_stability=0; context_restored=.true.; perturbed=.false.
	  nclo=0
	  do while(ncl.lt.ni)
!	   write(*,'("Current automorpism dim = ",i11)') auto !debug
  !Get the irreducible classification (IC):
	   do while(ncl.gt.nclo) !loop over vertex classes
!	    write(*,'(" Begin reduction cycle:")') !debug
!	    write(*,'("  Item IDs: ",i2,2x,32(1x,i2))') trn(0:ni) !debug
!	    write(*,'("  Classes : ",i2,2x,32(1x,i2))') vc(0),vc(trn(1:ni)) !debug
	    nclo=ncl; k1=1
	    do while(k1.le.ni)
	     k2=k1+1; do while(k2.le.ni); if(vc(trn(k2)).ne.vc(trn(k1))) exit; k2=k2+1; enddo ![k1:k2-1]: positions of the current class
   !Construct vertex descriptors:
!	     write(*,'(" Constructing vertex descriptors for positions: ",i4,1x,i4)') k1,k2-1 !debug
	     do i=k1,k2-1 !loop over positions within a class
	      j=trn(i) !vertex ID
	      m=gr%vertices(j)%conn_num !amount of connections of vertex j (which is currently occupies position i)
	      if(m.gt.0) then
	       sl(1:m)%conn_typ=gr%vertices(j)%slots(1:m)%conn_typ
	       sl(1:m)%vert_attr=vc(gr%vertices(j)%slots(1:m)%vert_attr)
	       ii(0)=+1; do k0=1,m; ii(k0)=k0; enddo
	       call sort_slots(m,sl,ii)
	       gr%descr(j)%slots(1:m)=sl(ii(1:m)); gr%descr(j)%vertex_class=vc(j)
	      endif
	     enddo
   !Order vertices according to their descriptors, possibly splitting the class:
!	     write(*,'(" Ordering vertices for positions: ",i4,1x,i4)') k1,k2-1 !debug
	     m=k2-k1 !length of the current class
	     n=gr%vertices(trn(k1))%conn_num !amount of connections for each vertex in the class
	     rloop: do k=1,n !loop over the rows of descriptors of the vertices (row-wise ordering of descriptors)
    !Get a row:
	      do k0=1,m
	       sl(k0)%vert_attr=gr%descr(trn(k1+k0-1))%slots(k)%vert_attr
	       sl(k0)%conn_typ=gr%descr(trn(k1+k0-1))%slots(k)%conn_typ
	      enddo
    !Order the row:
	      ii(0)=+1; do k0=1,m; ii(k0)=k0; enddo
	      call sort_slots(m,sl,ii); ii(1:m)=ii(1:m)+k1-1
	      trn(0)=trn(0)*ii(0); ii(1:m)=trn(ii(1:m)); trn(k1:k2-1)=ii(1:m)
    !Split the class:
	      ks=0; kf=0 !splitting
	      do k0=k1+1,k2-1
	       if(slot_gt(gr%descr(trn(k0))%slots(k),gr%descr(trn(k0-1))%slots(k))) then
	        ks=ks+1; if(kf.eq.0) kf=k0 !splitting occured; 1st splitting occured on element #kf (to be continued from this position)
	       endif
	       vc(trn(k0))=vc(trn(k0))+ks
	      enddo
	      if(ks.gt.0) then; if(k2.le.ni) vc(trn(k2:ni))=vc(trn(k2:ni))+ks; ncl=ncl+ks; vc(0)=ncl; k2=kf; exit rloop; endif
	     enddo rloop
	     k1=k2
	    enddo !k1: beginning position of a class
	   enddo !ncl>nclo
	   vcl(0:ni)=vc(0:ni)
  !Check stability of the irreducible classification:
	   if(ncl.lt.ni.or.checking_stability.gt.0) then
!	    write(*,'(" Checking stability: ",i4)') checking_stability !debug
   !First entrance:
	    if(checking_stability.eq.0) then
	     checking_stability=ni; class_checked=0; nclp=0
    !Allocate additional structures:
	     if(.not.add_alloc) then
	      add_alloc=.true.
	      if(ni.gt.max_dim_mem) then
	       allocate(save_vca(0:ni),STAT=ierr); if(ierr.ne.0) then; ierr=6; goto 999; endif
	       allocate(save_trna(0:ni),STAT=ierr); if(ierr.ne.0) then; ierr=7; goto 999; endif
	       allocate(trnpa(0:ni),STAT=ierr); if(ierr.ne.0) then; ierr=8; goto 999; endif
	       save_vc=>save_vca; save_trn=>save_trna; trnp=>trnpa
	      else
	       save_vc=>save_vcs; save_trn=>save_trns; trnp=>trnps
	      endif
	     endif
    !Save the current context:
	     call save_context !flag context_restored=.true.
	    endif
   !Cross-validation of degenerate elements in each class (perturb a degenerate element, try class reduction, save the result, restore the context):
	    cloop: do while(checking_stability.gt.0) !loop over all degenerate vertices
	     if(save_vc(save_trn(checking_stability)).eq.abs(class_checked)) then !this (degenerate) class has been already entered
	      if(class_checked.lt.0) then !process the result of mutating a vertex of the current (degenerate) class
	       class_checked=-class_checked !set positive class_checked
	       if(nclp.eq.0) then; i=-1; else; i=compare_classifications(ncl,trn(0:ni),gr,nclp,trnp(0:ni),grp); endif !-1:gr<grp; 0:gr=grp; +1:gr>grp
	       if(i.lt.0) then !new leading subclass: mark this element, unmark all previous elements from the same class
	        call copy_sparse_graph(clean,gr,grp); trnp(0:ni)=trn(0:ni); nclp=ncl !set new grp/trnp/nclp, characterizing the new leading subclass
	        save_vc(save_trn(checking_stability))=-save_vc(save_trn(checking_stability)) !negative save_vc marks elements to be separated from the class
	        do j=checking_stability+1,ni !unmark previous elements
	         k1=save_vc(save_trn(j)); if(abs(k1).ne.abs(class_checked)) exit; if(k1.lt.0) save_vc(save_trn(j))=-k1
	        enddo
	       elseif(i.eq.0) then !the same leading subclass: mark this element
	        save_vc(save_trn(checking_stability))=-save_vc(save_trn(checking_stability)) !negative save_vc marks elements to be separated from the class
	       endif
	       call restore_context !context_restored will be .true.
	       checking_stability=checking_stability-1 !to the previous vertex (whichever class it belongs to)
	      else !new vertex of the current (degenerate) class; context_restored=.true. on entrance
	       class_checked=-class_checked !set negative class_checked
	       call mutate_vertex(checking_stability)
	       context_restored=.false.
	       exit cloop !to the class reduction procedure (initiate by a mutation of not the last vertex of a degenerate class)
	      endif
	     else !entering a new (degenerate) class (starting from the last element of the class); context_restored=.true. on entrance
	      if(.not.context_restored) then !trap
	       write(*,*)'ERROR(combinatoric:graph_standard): context not properly restored (1): ',checking_stability
	       stop
	      endif
	      if(checking_stability.gt.1) then
	       if(save_vc(save_trn(checking_stability)).eq.save_vc(save_trn(checking_stability-1))) then !the last element of a degenerate class
	        class_checked=-save_vc(save_trn(checking_stability)) !set negative class_checked
	        call mutate_vertex(checking_stability)
	        nclp=0; context_restored=.false. !begining a new degenerate class
	        exit cloop !to the class reduction procedure (initiated by a mutation of the last vertex of a degenerate class)
	       else !the last element of the class is the only element of it
	        checking_stability=checking_stability-1 !go to the previous element (and class)
	       endif
	      else !checking_stability=1
	       checking_stability=checking_stability-1 !checking_stability=0 (will exit)
	      endif
	     endif
	    enddo cloop !checking_stability>0
	   endif
  !Split classes for unstable classifications or add perturbation for stable irreducible classifications with degeneracy:
	   if(checking_stability.eq.0.and.ncl.lt.ni) then
	    if(.not.context_restored) then !trap
	     write(*,*)'ERROR(combinatoric:graph_standard): context not properly restored (2): ',checking_stability
	     stop
	    endif
!	    write(*,'(" Stability analysis 1:")') !debug
!	    write(*,'("  Item IDs: ",i2,2x,32(1x,i2))') trn(0:ni) !debug
!	    write(*,'("  Classes : ",i2,2x,32(1x,i2))') save_vc(0),save_vc(trn(1:ni)) !debug
   !Unmark (degenerate) classes with all elements mutated (deactivate trivial mutations):
	    k1=0; k2=0; k3=0
	    do i=1,ni
	     j=save_vc(trn(i))
	     if(abs(j).eq.k1+1) then !next class
	      if(k3.eq.0.and.i.gt.k2+1) save_vc(trn(k2:i-1))=-save_vc(trn(k2:i-1)) !reverese trivial mutations in a degenerate class
	      k1=k1+1; k2=i; k3=0 !k2 - begining of the next class
	     elseif(abs(j).lt.k1.or.abs(j).gt.k1+1) then
	      write(*,*)'ERROR(combinatoric:graph_standard): non-monotonic class ordering detected!' !trap
	      stop
	     endif
	     if(j.gt.0) k3=k3+1
	    enddo
	    if(k3.eq.0.and.ni.gt.k2) save_vc(trn(k2:ni))=-save_vc(trn(k2:ni)) !reverese trivial mutations in a degenerate class
   !Now, negative save_vc(:) points to the elements of degenerate classes which are to be separated into new classes:
!	    write(*,'(" Stability analysis 2:")') !debug
!	    write(*,'("  Item IDs: ",i2,2x,32(1x,i2))') trn(0:ni) !debug
!	    write(*,'("  Classes : ",i2,2x,32(1x,i2))') save_vc(0),save_vc(trn(1:ni)) !debug
   !Split classes (if the splitting occured naturally):
	    k1=0; k2=0; k3=0
	    do i=1,ni
	     j=save_vc(trn(i))
	     if(abs(j).eq.k1+1) then !next class
	      k1=k1+1; k2=k2+1+k3; k3=0 !k1 - old class #; k2 - new class #
	     endif
	     if(j.lt.0) then; vc(trn(i))=k2+1; k3=1; else; vc(trn(i))=k2; endif !shift by one the class of the elements separated
	    enddo
	    vc(0)=k2+k3 !new number of classes; vc(1:ni) - new classes
	    if(.not.perturbed) vcl(0:ni)=vc(0:ni)
!	    write(*,'(" Stability analysis 3:")') !debug
!	    write(*,'("  Item IDs: ",i2,2x,32(1x,i2))') trn(0:ni) !debug
!	    write(*,'("  Classes : ",i2,2x,32(1x,i2))') vc(0),vc(trn(1:ni)) !debug
   !Perturb and split:
	    if(vc(0).eq.ncl) then !the degenerate classification was stable: perturb
!	     write(*,'("DEBUG(graph_standard): Vertex Numbers: ",i2,2x,99(1x,i2))') trn(0:ni) !debug
!	     write(*,'("DEBUG(graph_standard): Stable Classes: ",i2,2x,99(1x,i2))') vc(0),vc(trn(1:ni)) !debug
	     perturbed=.true.
	     do i=1,ni-1
	      if(vc(trn(i)).eq.vc(trn(i+1))) then
	       j=2; do while(i+j.le.ni); if(vc(trn(i+j)).ne.vc(trn(i))) exit; j=j+1; enddo
	       if(auto.gt.0) then; if(huge(1_4)/j.ge.auto) then; auto=auto*j; else; auto=0; endif; endif !modify the automorphic factor (0 - out of int*4 rnange)
	       vc(trn(i+1:ni))=vc(trn(i+1:ni))+1; ncl=ncl+1; vc(0)=ncl
!	       write(*,'(" Perturbed position: ",i12)') i !debug
	       exit
	      endif
	     enddo
	    else !the degenerate classification was unstable: next round of class reduction
	     call merge_sort_key_int(ni,vc(1:ni),trn(0:ni))
	     ncl=vc(0) !new number of classes
	    endif
	   endif
	  enddo !ncl<ni
	 else !ncl=ni at the beginning
	  do j=1,ni
	   m=gr%vertices(j)%conn_num !amount of connections of vertex j (which is currently occupies position i)
	   if(m.gt.0) then
	    sl(1:m)%conn_typ=gr%vertices(j)%slots(1:m)%conn_typ
	    sl(1:m)%vert_attr=vc(gr%vertices(j)%slots(1:m)%vert_attr)
	    ii(0)=+1; do k0=1,m; ii(k0)=k0; enddo
	    call sort_slots(m,sl,ii)
	    gr%descr(j)%slots(1:m)=sl(ii(1:m)); gr%descr(j)%vertex_class=vc(j)
	   endif
	  enddo
	 endif !ncl<ni or not
!	 write(*,'(" Automorphic factor = ",i12)') auto !debug
!Restore vertex classes from the first stable irreducible classification:
	 do i=1,ni; if(vc(trn(i)).ne.i) then; ierr=10; goto 999; endif; enddo
!	 gr%descr(1:ni)%vertex_class=vc(1:ni)
	 gr%descr(1:ni)%vertex_class=vcl(1:ni) !true vertex classes
!	 write(*,'(" STANDARD PERMUTATION : ",i2,2x,1024(1x,i2))') trn(0:ni) !debug
!	 write(*,'(" AUTO & VERTEX CLASSES: ",i10,2x,1024(1x,i2))') auto,vcl(trn(1:ni)) !debug
	else
	 ierr=-1; goto 999 !empty graph
	endif
999	if(associated(ii)) nullify(ii); if(associated(sl)) nullify(sl)
        if(associated(vc)) nullify(vc); if(associated(vcl)) nullify(vcl)
	if(allocated(iia)) deallocate(iia); if(allocated(sla)) deallocate(sla)
	if(allocated(vca)) deallocate(vca); if(allocated(vcla)) deallocate(vcla)
	if(add_alloc) then
	 if(associated(save_vc)) nullify(save_vc); if(associated(save_trn)) nullify(save_trn); if(associated(trnp)) nullify(trnp)
	 if(allocated(save_vca)) deallocate(save_vca); if(allocated(save_trna)) deallocate(save_trna)
	 if(allocated(trnpa)) deallocate(trnpa)
	endif
!	write(*,'("Exited GRAPH_STANDARD.")') !debug
	return

	contains

	logical function slot_gt(sl1,sl2) !.TRUE. if SLOT1>SLOT2, .FALSE. otherwise
	 implicit none
	 type(slot), intent(in):: sl1,sl2
	 if(sl1%vert_attr.gt.sl2%vert_attr) then
	  slot_gt=.true.
	 elseif(sl1%vert_attr.lt.sl2%vert_attr) then
	  slot_gt=.false.
	 else
	  if(dble(sl1%conn_typ).gt.dble(sl2%conn_typ)) then
	   slot_gt=.true.
	  elseif(dble(sl1%conn_typ).lt.dble(sl2%conn_typ)) then
	   slot_gt=.false.
	  else
	   if(dimag(sl1%conn_typ).gt.dimag(sl2%conn_typ)) then
	    slot_gt=.true.
	   elseif(dimag(sl1%conn_typ).lt.dimag(sl2%conn_typ)) then
	    slot_gt=.false.
	   else
	    slot_gt=.false.
	   endif
	  endif
	 endif
	 return
	end function slot_gt

	subroutine save_context
	 implicit none
	 save_ncl=ncl; save_mcl=mcl; save_nclo=nclo
	 save_trn(0:ni)=trn(0:ni)
	 save_vc(0:ni)=vc(0:ni)
	 call copy_sparse_graph(clean,gr,save_gr)
	 context_restored=.true.
!	 write(*,'(" Context saved.")') !debug
	 return
	end subroutine save_context

	subroutine restore_context
	 implicit none
	 integer j1
	 ncl=save_ncl; mcl=save_mcl; nclo=save_nclo
	 trn(0:ni)=save_trn(0:ni)
	 do j1=0,ni; vc(j1)=abs(save_vc(j1)); enddo !save_vc is used for class splitting (minus sign marks the elements to be separated)
	 call copy_sparse_graph(no_clean,save_gr,gr)
	 context_restored=.true.
!	 write(*,'(" Context restored.")') !debug
	 return
	end subroutine restore_context

	subroutine mutate_vertex(vert_pos) !assigns a new class (+1) to the vertex at position <vert_pos> and shifts it to the end of the class
	 implicit none
	 integer, intent(in):: vert_pos
	 integer j1,j2
	 j1=vert_pos; j2=trn(vert_pos)
	 do while(j1.lt.ni)
	  if(vc(trn(j1+1)).eq.vc(j2)) then
	   trn(j1)=trn(j1+1); trn(0)=-trn(0); j1=j1+1
	  else
	   exit
	  endif
	 enddo
	 trn(j1)=j2; vc(j2)=vc(j2)+1; vc(trn(j1+1:ni))=vc(trn(j1+1:ni))+1; ncl=ncl+1; vc(0)=ncl
!	 write(*,'(" Mutated vertex: ",i4,1x,i4)') vert_pos,j2 !debug
	 return
	end subroutine mutate_vertex

	integer function compare_classifications(ncl1,trn1,gr1,ncl2,trn2,gr2) !-1:gr1<gr2; 0:gr1=gr2; +1:gr1>gr2
	 implicit none
	 integer, intent(in):: ncl1,ncl2,trn1(0:ni),trn2(0:ni)
	 type(graph_sparse), intent(in):: gr1,gr2
	 integer j1,j2,j3
	 compare_classifications=0
!	 return !debug (not intermediate class reduction)
	 if(ncl1.gt.ncl2) then
	  compare_classifications=-1
	 elseif(ncl1.lt.ncl2) then
	  compare_classifications=+1
	 else
	  do j1=1,ni
	   j2=gr1%descr(trn1(j1))%vertex_class-gr2%descr(trn2(j1))%vertex_class
	   j3=gr1%descr(trn1(j1))%conn_num-gr2%descr(trn2(j1))%conn_num
	   if(j2.ne.0.or.j3.ne.0) exit
	  enddo
	  if(j2.gt.0) then
	   compare_classifications=-1
	  elseif(j2.lt.0) then
	   compare_classifications=+1
	  else
	   if(j3.gt.0) then
	    compare_classifications=-1
	   elseif(j3.lt.0) then
	    compare_classifications=+1
	   else
	    loop1: do j1=1,ni
	     do j2=1,gr1%descr(trn1(j1))%conn_num
	      compare_classifications=compare_descr_slots(gr1%descr(trn1(j1))%slots(j2),gr2%descr(trn2(j1))%slots(j2))
	      if(compare_classifications.ne.0) exit loop1
	     enddo
	    enddo loop1
	   endif
	  endif
	 endif
	 return
	end function compare_classifications

	end subroutine graph_standard
!----------------------------------------------------------------------------------------------------
	subroutine print_graph_status(gr) !prints a detailed allocation status for type(graph_sparse), together with the number of connections for each vertex
	implicit none
	type(graph_sparse), intent(in):: gr
	integer i,j,k,l,m,n,k0,k1,k2,k3,ierr
	n=gr%gcard
	write(*,'(" Printing the status of a graph of cardinality ",i4,1x,l1,l1":")') n,allocated(gr%vertices),allocated(gr%descr)
	if(n.gt.0) then
	 write(*,'(1x,256(1x,i4,1x,l1,1x,l1,1x))') &
	  (gr%vertices(j)%conn_num,allocated(gr%vertices(j)%slots),allocated(gr%descr(j)%slots),j=1,n)
	endif
	return
	end subroutine print_graph_status
!----------------------------------------------------
	logical function graph_isomorphic(g1,g2,ierr)
!This function determines whether two given complex hermitian adjacency matrices are isomorphic or not.
!If they are, the underlying graphs are isomorphic. Hermiticity of the adjacency matrices must be STRICTLY imposed.
!Class Reduction Algorithm is used [D.I.Lyakh, 2011]. Algorithm's worst case is O(N^5*logN*logM), 
!where M is the number of connection types (i.e., the number of distinct edge types, 1 for simple graphs).
!INPUT:
! - g1,g2 - hermitian complex adjacency matrices given in a sparse (packed) form;
!OUTPUT:
! - .true./.false;
! - g1,g2 with ordered vertex slots.
!NOTES:
! - GRAPH_STANDARD rearranges vertex slots in an ordered form (according to COMPARE_DESCR_SLOTS).
	implicit none
	type(graph_sparse), intent(inout):: g1,g2
	integer, parameter:: max_in_mem=1023
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf,ierr
	integer, target:: trns1(0:max_in_mem),trns2(0:max_in_mem)
	integer, allocatable, target:: trna1(:),trna2(:)
	integer, pointer:: trn1(:),trn2(:)
	integer auto1,auto2

	ierr=0; graph_isomorphic=.true.
	n=g1%gcard
	if(n.ge.0) then
	 if(g2%gcard.eq.n) then
	  if(n.gt.0) then
	   if(n.gt.max_in_mem) then
	    allocate(trna1(0:n),trna2(0:n),STAT=ierr); if(ierr.ne.0) then; ierr=-7; goto 999; endif
	    trn1=>trna1; trn2=>trna2
	   else
	    trn1=>trns1; trn2=>trns2
	   endif
 !Canonicalize the graphs (obtain their full invariants):
	   call graph_standard(g1,trn1,auto1,ierr); if(ierr.ne.0) then; graph_isomorphic=.false.; goto 999; endif
	   call graph_standard(g2,trn2,auto2,ierr); if(ierr.ne.0) then; graph_isomorphic=.false.; goto 999; endif
 !Check that the canonical classification is non-degenerate:
	   if(allocated(g1%descr).and.allocated(g2%descr)) then
!	    do i=1,n
!	     if(g1%descr(trn1(i))%vertex_class.ne.i) then; ierr=-3; graph_isomorphic=.false.; goto 999; endif
!	     if(g2%descr(trn2(i))%vertex_class.ne.i) then; ierr=-4; graph_isomorphic=.false.; goto 999; endif
!	    enddo
	   else
	    ierr=-5; graph_isomorphic=.false.
	    goto 999
	   endif
 !Compare vertex descriptors:
	   vloop: do i=1,n
	    m=g1%vertices(trn1(i))%conn_num
	    if(g2%vertices(trn2(i))%conn_num.eq.m.and.g1%vertices(trn1(i))%color.eq.g2%vertices(trn2(i))%color) then
	     do j=1,m
	      if(compare_descr_slots(g1%descr(trn1(i))%slots(j),g2%descr(trn2(i))%slots(j)).ne.0) then
	       graph_isomorphic=.false.; exit vloop
	      endif
	     enddo
	    else
	     graph_isomorphic=.false.
	     exit vloop
	    endif
	   enddo vloop
	  else !n=0 (empty graphs)
	   graph_isomorphic=.true.
	  endif
	 else
	  graph_isomorphic=.false.
	 endif
	else
	 graph_isomorphic=.false.; ierr=-6 !negative cardinality
	endif
999	if(associated(trn1)) nullify(trn1); if(associated(trn2)) nullify(trn2)
	if(allocated(trna1)) deallocate(trna1); if(allocated(trna2)) deallocate(trna2)
	return
	end function graph_isomorphic
!------------------------------------------
	real(8) function dint2frac(dpi,prec)
!This function converts a double-precision integer into a double-precision fractional number less than 1d0
!by reflecting the sequence of decimal digits against the decimal point. For example, 6529450. --> .0549256
!INPUT:
! - dpi - double-precision integer (possible fractional part will be ignored);
! - prec - max number of decimal digits to convert (do not go beyond 12).
!OUTPUT:
! - dint2frac - converted value (pure fractional).
	implicit none
	real(8), intent(in):: dpi
	integer, intent(in):: prec
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf,ierr
	integer(8) i0,i1
	real(8) df,dcheck
	dint2frac=0d0 !clear all digits
	if(prec.gt.0.and.prec.le.12) then
	 i0=sign(int(dpi,8),+1_8); i1=0_8; n=0
	 do while(i0.ne.0_8)
	  i1=i1*10_8+mod(i0,10_8); i0=i0/10_8
	  n=n+1; if(n.eq.prec) exit
	 enddo
	 df=sign(dble(i1)*1d1**(-prec),dpi)
	 dcheck=frac2dint(df,prec)
	 if(dcheck.eq.dint(dpi)) then
	  dint2frac=df
	 else
	  write(*,*)'ERROR(combinatoric:dint2frac): numerical instability detected: ',dpi,df,dcheck
	  stop
	 endif
	elseif(prec.lt.0.or.prec.gt.12) then
	 write(*,*)'ERROR(combinatoric:dint2frac): invalid precision requested: ',prec
	 stop
	endif
	return
	end function dint2frac
!------------------------------------------
	real(8) function frac2dint(dpf,prec)
!This is the inverse function with respect to DINT2FRAC.
!INPUT:
! - dpf - double-precision fractional number with zero integer part (otherwise the integer part is ignored);
! - prec - number of fractional digits to convert (do not use more than 12 for double precision);
!OUTPUT:
! - frac2dint - reflected double-precision integer.
	implicit none
	real(8), intent(in):: dpf
	integer, intent(in):: prec
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf,ierr
	integer(8) i0,i1
	real(8) df,dcheck
	frac2dint=0d0
	if(prec.gt.0.and.prec.le.12) then
	 i0=int(abs(dpf)*1d1**prec,8); i1=0_8
	 do n=1,prec
	  i1=i1*10_8+mod(i0,10_8); i0=i0/10_8
	 enddo
	 df=sign(dble(i1),dpf)
	 dcheck=dint2frac(df,prec)
	 if(dcheck.eq.dpf-dint(dpf)) then
	  frac2dint=df
	 else
	  write(*,*)'ERROR(combinatoric:frac2dint): numerical instability detected: ',dpf,df,dcheck
	  stop
	 endif
	elseif(prec.lt.0.or.prec.gt.12) then
	 write(*,*)'ERROR(combinatoric:frac2dint): invalid precision requested: ',prec
	 stop
	endif
	return
	end function frac2dint
!-----------------------------------------
	subroutine graph2matrix(gr,nv,adj)
!This subroutine converts a graph represented as type(graph_sparse) to a complex hermitian adjacency matrix.
!Vertex numeration starts from 1.
!INPUT:
! - gr - graph type(graph_sparse);
!OUTPUT:
! - nv - number of vertices;
! - adj(:,:) - complex(8) hermitian adjacency matrix (allocated here!).
!NOTES:
! - The adjacency matrix is traversed in a column-wise order (connections of vertex I constitute the Ith column of the matrix).
	implicit none
	type(graph_sparse), intent(in):: gr
	integer, intent(out):: nv
	complex(8), allocatable, intent(inout):: adj(:,:)
	real(8), parameter:: herm_thresh=dp_zero_thresh !hermiticity threshold (real*8 comparison zero threshold)
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf,ierr
	real(8) diff,cl,cf,ct

	if(allocated(adj)) deallocate(adj)
	if(gr%gcard.gt.0) then
	 nv=gr%gcard
	 allocate(adj(1:nv,1:nv),STAT=ierr)
	 if(ierr.ne.0) then; write(*,*)'ERROR(combinatoric:graph2matrix): allocation failed: ',nv; stop; endif
	 adj(:,:)=(0d0,0d0)
	 do i=1,nv !vertex numbers [1..nv]
	  cl=abs(gr%vertices(i)%color)
	  if(cl-dint(cl).le.dp_zero_thresh) then
	   cl=dint(cl); cf=0d0
	   do j=1,gr%vertices(i)%conn_num
	    k=gr%vertices(i)%slots(j)%vert_attr
	    if(k.gt.0.and.k.le.nv) then
	     if(k.ne.i) then !off-diagonal element (a complex*8 connection attribute)
	      adj(k,i)=gr%vertices(i)%slots(j)%conn_typ
	     else !diagonal element, <color.number_of_loops>
	      ct=dble(gr%vertices(i)%slots(j)%conn_typ)
	      if(dabs(dimag(gr%vertices(i)%slots(j)%conn_typ)).le.dp_zero_thresh.and.dabs(ct-dint(ct)).le.dp_zero_thresh) then !no imaginary part, integer real part
	       ct=dint(ct) !real integer
	       if(cf.eq.0d0) then
	        if(ct.ge.0d0) then
	         if(cl.ne.0d0.and.ct.ne.0d0) then !check double-precision overflow when coding a diagonal element
	          if(dlog10(cl)+dlog10(ct).gt.12d0) then !no more than 12 decimal digits total for <vertex_color.number_of_loops>
	           write(*,*)'ERROR(combinatoric:graph2matrix): diagonal overflow: ',cl,ct,i
	           stop
	          endif
	         endif
	         cf=dint2frac(ct,dp_codon_len) !set vertex loops: convert a pure integer into a pure fractional
	        else
	         write(*,*)'ERROR(combinatoric:graph2matrix): negative loop attribute detected: ', &
	          gr%vertices(i)%slots(j)%conn_typ,ct,i,j
	         stop
	        endif
	       else
	        write(*,*)'ERROR(combinatoric:graph2matrix): repeated self-connection (loop) detected: ', &
	         cf,i,j,gr%vertices(i)%slots(j)%conn_typ
	        stop
	       endif
	      else
	       write(*,*)'ERROR(combinatoric:graph2matrix): unable to process non-integer self-connection attribute (loop): ', &
	        gr%vertices(i)%slots(j)%conn_typ,i,j
	       stop
	      endif
	     endif
	    else
	     write(*,*)'ERROR(combinatoric:graph2matrix): invalid vertex number specified in a connection: ',i,j,k,nv
	     stop
	    endif
	   enddo !j
	   adj(i,i)=dcmplx(sign(cl+cf,gr%vertices(i)%color),0d0) !the sign of the diagonal element is determined by the vertex color
	  else
	   write(*,*)'ERROR(combinatoric:graph2matrix): unable to process non-integer vertex colors: ',gr%vertices(i)%color
	   stop
	  endif
	 enddo !i
	 diff=check_hermiticity(nv,adj,.true.)
	 if(diff.gt.herm_thresh) then !trap
	  write(*,*)'ERROR(combinatoric:graph2matrix): hermiticity check failed: ',diff
	  stop
	 endif
	else
	 nv=0 !no allocation is done if there were no vertices
	endif
	return
	end subroutine graph2matrix
!-----------------------------------------
	subroutine matrix2graph(nv,adj,gr)
!This subroutine converts a complex hermitian adjacency matrix into the type(graph_sparse) representation.
!INPUT:
! - nv - number of vertices;
! - adj(1:nv,1:nv) - complex(8) hermitian adjacency matrix.
!OUTPUT:
! - gr - graph of type(graph_sparse).
!NOTES:
! - Empty elements of the adjacency matrix must be exactly zero (0d0,0d0)!
	implicit none
	integer, intent(in):: nv
	complex(8), intent(inout):: adj(1:,1:)
	type(graph_sparse), intent(inout):: gr
	real(8), parameter:: herm_thresh=dp_zero_thresh !hermiticity threshold
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf,ierr
	real(8) diff,de,dl

	if(nv.ge.0) then
	 call nullify_sparse_graph(gr)
	 if(nv.gt.0) then
	  diff=check_hermiticity(nv,adj(1:nv,1:nv),.false.)
	  if(diff.le.herm_thresh) then
	   gr%gcard=nv; allocate(gr%vertices(1:nv))
	   do j=1,nv !columns
	    if(dabs(dimag(adj(j,j))).le.dp_zero_thresh) then !diagonal elements must be real
	     de=dble(adj(j,j))
	     gr%vertices(j)%color=dint(de) !vertex color := integer part of the real part of the diagonal element
	     if(dabs(de-dint(de)).gt.dp_zero_thresh) then !vertex self-connection (loop)
	      dl=frac2dint(de-dint(de),dp_codon_len); m=1
	     else !no loops
	      dl=0d0; m=0
	     endif
	     do i=1,j-1; if(adj(i,j).ne.(0d0,0d0)) m=m+1; enddo; do i=j+1,nv; if(adj(i,j).ne.(0d0,0d0)) m=m+1; enddo !count connections to other vertices
	     allocate(gr%vertices(j)%slots(1:m)); gr%vertices(j)%conn_num=m
	     if(dl.ne.0d0) then; gr%vertices(j)%slots(1)=slot(j,dcmplx(dl,0d0)); m=1; else; m=0; endif !possible self-connection
	     do i=1,j-1
	      if(adj(i,j).ne.(0d0,0d0)) then
	       m=m+1; gr%vertices(j)%slots(m)=slot(i,adj(i,j))
	      endif
	     enddo
	     do i=j+1,nv
	      if(adj(i,j).ne.(0d0,0d0)) then
	       m=m+1; gr%vertices(j)%slots(m)=slot(i,adj(i,j))
	      endif
	     enddo
	    else
	     write(*,*)'ERROR(combinatoric:matrix2graph): imaginary diagonal element detected: ',adj(j,j),j
	     stop
	    endif
	   enddo !j
	  else
	   write(*,*)'ERROR(combinatoric:matrix2graph): hermiticity check failed: ',diff
	   stop
	  endif
	 endif
	else
	 write(*,*)'ERROR(combinatoric:matrix2graph): negative number of vertices: ',nv
	 stop
	endif
	return
	end subroutine matrix2graph
!--------------------------------------------------------
	real(8) function check_hermiticity(nv,adj,correct)
!This function checks hermiticity of a complex(8) matrix, returning the maximal deviation.
!If correct=.true., the hermiticity will be imposed by the upper triangle.
	implicit none
	integer, intent(in):: nv
	complex(8), intent(inout):: adj(1:,1:)
	logical, intent(in):: correct
	integer j0,j1
	check_hermiticity=0d0
!check the diagonal part:
	do j1=1,nv
	 check_hermiticity=max(check_hermiticity,dabs(dimag(adj(j1,j1))))
	enddo
!check the off-diagonal part:
	do j1=2,nv
	 do j0=1,j1-1
	  check_hermiticity=max(check_hermiticity,cdabs(dconjg(adj(j0,j1))-adj(j1,j0)))
	 enddo
	enddo
!correct if requested:
	if(correct) then
	 do j1=1,nv
	  adj(j1,j1)=dcmplx(dble(adj(j1,j1)),0d0)
	  do j0=1,j1-1
	   adj(j1,j0)=dconjg(adj(j0,j1)) !lower trianle := upper triangle
	  enddo
	 enddo
	endif
	return
	end function check_hermiticity
!----------------------------------------------
	subroutine star_graph(nconn,gdim,graph)
!This subroutine generates the adjacency matrix corresponding to a "star" graph of cardinality gdim.
!"Star" graph is a graph in which each vertex is connected exactly with the other NCONN*2 vertices (band adjacency matrix).
!The "star" graph generated here is an undirected mono-color graph without loops.
!INPUT:
! - nconn - half-width of the band;
! - gdim - dimension of the adjacency matrix;
!OUTPUT:
! - graph - band adjacency matrix (complex*8 hermitian).
!NOTES:
! - nconn=0 corresponds to a diagonal matrix;
! - nconn=1 corresponds to a symmetrical polyhedron graph (D_gdim_h symmetry point group);
! - nconn*2=gdim-1 corresponds to a fully-automorphic graph (PG=GP for all P).
! - If the requested graph cardinality is zero, this subroutine does nothing.
	implicit none
	integer, intent(in):: nconn,gdim
	complex(8), intent(out):: graph(1:,1:)
	complex(8),parameter:: diag=(1d0,0d0),offdiag=(1d0,0d0)
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf,ierr

	if(gdim.ge.0.and.nconn.ge.0.and.nconn*2.lt.gdim) then
	 if(gdim.gt.0) then
	  graph(1:gdim,1:gdim)=(0d0,0d0)
!Diagonal:
	  do i=1,gdim; graph(i,i)=diag; enddo
	  if(nconn.gt.0) then
!Off-diagonal:
	   do i=1,gdim
	    do j=i+1,min(i+nconn,gdim); graph(j,i)=offdiag; enddo
	    do j=1,nconn-(gdim-i); graph(j,i)=offdiag; enddo
	    do j=i-1,max(i-nconn,1),-1; graph(j,i)=offdiag; enddo
	    do j=gdim,gdim-(nconn-i),-1; graph(j,i)=offdiag; enddo
	   enddo
	  endif
	 endif
	else
	 write(*,*)'ERROR(combinatoric:star_graph): invalid arguments: ',gdim,nconn
	 stop
	endif
	return
	end subroutine star_graph
!-------------------------------------------------------------------------------
	subroutine random_graph(color_range,valence_range,bond_range,gdim,graph)
!This subroutine generates a random hermitian adjacency matrix of dimension gdim,
!representing an undirected colored multigraph without loops.
!INPUT:
! - color_range - limit for the vertex color: [0..color_range];
! - valence_range - limit for the vertex valence: [0..valence_range];
! - bond_range - limit for the bond multiplicitiy [0..bond_range];
! - gdim - dimension of the adjacency matrix (number of vertices, graph cardinality);
!OUTPUT:
! - graph - complex(8) hermitian adjacency matrix.
!NOTES:
! - If the requested graph dimension is zero, this subroutine does nothing.
	implicit none
	integer, intent(in):: color_range,valence_range,bond_range,gdim
	complex(8), intent(out):: graph(1:,1:)
	integer, parameter:: rnd_chunk=2**12
	integer, parameter:: rnd_chunk2=rnd_chunk/2
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf,ierr
	real(8) rnd_buf(1:rnd_chunk),val

	if(gdim.gt.0.and.color_range.ge.0.and.valence_range.gt.0.and.bond_range.gt.0) then
!clear:
	 graph(1:gdim,1:gdim)=(0d0,0d0)
!set diagonal:
	 do i=1,gdim,rnd_chunk
	  l=min(rnd_chunk,gdim-i+1)
	  call random_number(rnd_buf(1:l))
	  do j=i,i+l-1
	   m=int(rnd_buf(j-i+1)*dble(color_range+1)); if(m.gt.color_range) m=color_range
	   graph(j,j)=dcmplx(dble(m),0d0)
	  enddo
	 enddo
!set off-diagonal:
 !loop over columns:
	 do j=1,gdim-1
 !import and count existing connections in the column j (upper triangle):
	  m=0
	  do i=1,j-1
	   if(graph(j,i).ne.(0d0,0d0)) then; graph(i,j)=dconjg(graph(j,i)); m=m+1; endif
	  enddo
 !set new connections:
	  n=gdim-j !number of free entrees in column j (lower triangle)
	  k0=min(n,valence_range-m) !max number of new connections
	  if(k0.gt.0) then
	   call random_number(val); m=int(val*dble(k0+1)); if(m.gt.k0) m=k0 !actual number of new connections to set
	   do i=1,m,rnd_chunk2
	    l=min(rnd_chunk2,m-i+1)
	    call random_number(rnd_buf(1:l)); rnd_buf(1:l)=rnd_buf(1:l)*dble(n+1) !free entries to occupy
	    call random_number(rnd_buf(l+1:l+l)); rnd_buf(l+1:l+l)=rnd_buf(l+1:l+l)*dble(bond_range)+1d0 !connection multiplicities
	    do k=1,l
	     k1=int(rnd_buf(k)); if(k1.gt.n) k1=n
	     k2=int(rnd_buf(l+k)); if(k2.gt.bond_range) k2=bond_range
	     graph(j+k1,j)=dcmplx(dble(k2),0d0) !set an undirected edge of multiplicity k2: 0<k2<=bond_range
	    enddo !k
	   enddo !i
	  endif
	 enddo !j
	elseif(gdim.lt.0) then
	 write(*,*)'ERROR(combinatoric:random_graph): negative graph cardinality requested: ',gdim
	 stop
	endif
	return
	end subroutine random_graph
!--------------------------------------------------------
	subroutine remove_random_edges(edge_num,gd,graph)
!This subroutine randomly removes EDGE_NUM edges from a graph (as far as possible).
!INPUT:
! - edge_num - number of edges to remove (0 - a random number will be generated);
! - gd - dimension of the adjacency matrix of the graph;
! - graph(1:gd,1:gd) - complex hermitian adjacency matrix of the graph;
!OUTPUT:
! - edge_num - number of edges removed;
! - graph(1:gd,1:gd) - modified adjacency matrix of the graph (without EDGE_NUM edges, if was possible).
!NOTES:
! - No hermiticity check is done here.
! - If the number of edges in the graph is less than the requested amount to delete, all edges will be deleted.
	implicit none
	integer, intent(in):: gd
	integer, intent(inout):: edge_num
	complex(8), intent(inout):: graph(1:,1:)
	integer, parameter:: rnd_chunk=2**12
	integer i,j,k,l,m,n,k0,k1,k2,k3,ks,kf,ierr
	real(8) rn,remove_thresh,rnd_buf(1:rnd_chunk)

	remove_thresh(k,l)=dble(k-l)/dble(k) !k>=l, k!=0: k - amount of objects left; l - amount of objects to still be removed.

	if(gd.ge.0.and.edge_num.ge.0) then
	 m=gd*(gd-1)/2 !number of off-diagonal elements in the upper (or lower) triangle
	 if(m.gt.0) then
!count edges in the upper triangle:
	  n=0 !will be the number of edges in the graph
	  do j=1,gd
	   do i=1,j-1
	    if(graph(i,j).ne.(0d0,0d0)) n=n+1
	   enddo
	  enddo
	  if(edge_num.eq.0) then !get the number of edges to remove (if it was not specified originally)
	   call random_number(rn); edge_num=int(rn*dble(n+1))
	  endif
	  if(edge_num.gt.n) edge_num=n
!remove random edges:
	  if(edge_num.gt.0) then
	   k=n; l=edge_num; k0=0 !initial amounts
	   cloop: do j=1,gd
	    do i=1,j-1
	     if(graph(i,j).ne.(0d0,0d0)) then
	      if(k0.eq.0) then; k0=min(rnd_chunk,k); call random_number(rnd_buf(1:k0)); endif
	      if(rnd_buf(k0).ge.remove_thresh(k,l)) then !remove edge
	       graph(i,j)=(0d0,0d0); graph(j,i)=(0d0,0d0)
	       l=l-1; if(l.eq.0) exit cloop
	      endif
	      k=k-1; k0=k0-1
	     endif
	    enddo
	   enddo cloop
	  endif
	 endif
	else
	 write(*,*)'ERROR(combinatoric:remove_random_edges): invalid arguments passed: ',gd,edge_num
	 stop
	endif
	return
	end subroutine remove_random_edges

	end module combinatoric
