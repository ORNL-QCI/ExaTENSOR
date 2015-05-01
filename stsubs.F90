!Standard subroutines/functions often used by me.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISON: 2015/05/01
	MODULE STSUBS

!Parameters:
        logical, private:: verbose=.false.
	real(8), parameter:: pi=3.14159265358979d0
	real(8), parameter:: BOHR=0.529177249d0    !Bohrs to Angstroms conversion factor

	PUBLIC:: ARRAY2STRING!converts a character*1 array into a string
	PUBLIC:: CAP_ASCII   !makes all English letters capital
	PUBLIC:: CHARNUM     !converts a number given as a string to real*8 and integer numbers.
	PUBLIC:: CREATE_LINE !creates a table line in .txt format with ; separator.
	PUBLIC:: ICHARNUM    !converts a number given as a string to the integer number.
	PUBLIC:: IFCL        !calculates factorial.
	PUBLIC:: IS_IT_NUMBER!checks if the character is an ASCII number
	PUBLIC:: IS_IT_LETTER!checks if the character is an ASCII letter
	PUBLIC:: ITRSIGN     !determines a sign of a given transposition.
	PUBLIC:: LONGNUMCHAR !converts a long integer number to the character representation.
	PUBLIC:: MARKCHF     !counts how many non-blank fields a string contains.
	PUBLIC:: MATMAT      !multiplies matrix on matrix.
	PUBLIC:: MATMATT     !multiplies matrix on transposed matrix.
	PUBLIC:: MATTMAT     !multiplies transposed matrix on matrix.
	PUBLIC:: MATTMATT    !multiplies transposed matrix on transposed matrix.
	PUBLIC:: NAME_HASH   !returns a hash-mask for a given string
	PUBLIC:: NOCOMMENT   !removes comments from a line (!,#)
	PUBLIC:: NORMV       !normalizes a vector.
	PUBLIC:: NOSPACE     !removes blanks from a string.
	PUBLIC:: NOT_A_NUMBER!checks whether a given string solely contains a decimal number
	PUBLIC:: NUMCHAR     !converts an integer number to character representation.
	PUBLIC:: PRINTL      !prints a line of characters
	PUBLIC:: RAND_STR    !returns a random string
	PUBLIC:: ROTS        !rotates an array of points in a 3d-space.
	PUBLIC:: SMALL_ASCII !makes all English letters small
	PUBLIC:: STRING2ARRAY!converts a string into a character*1 array
	PUBLIC:: STRSEARCH   !searches a given fragment in a string.
	PUBLIC:: SYMBOL_TYPE !function which tells you whether the given symbol is a space/tab (0), number (1), letter (2), or other (-1).
	PUBLIC:: TPAUSE      !pause for a given number of seconds.
	PUBLIC:: VALCHAR     !converts a real number to the string of characters.
	PUBLIC:: WAIT_DELAY  !pause for a given time.
	PUBLIC:: WAIT_PRESS  !subroutine makes pause until user presses a key
	PUBLIC:: WR_MAT_IN   !writes a matrix of integers to the screen
	PUBLIC:: WR_MAT_IN8  !writes a matrix of integer8's to the screen
	PUBLIC:: WR_MAT_SP   !writes a matrix of single precision elements to the screen.
	PUBLIC:: WR_MAT_DP   !writes a matrix of double precision elements to the screen.
	PUBLIC:: WR_MAT_DC   !writes a matrix of double complex elements to the screen.
	PUBLIC:: WR_VEC_SP   !writes a vector of single precision elements to the screen.
	PUBLIC:: WR_VEC_DP   !writes a vector of double precision elements to the screen.

	CONTAINS
!------------------------------------------------
	subroutine array2string(str,ar1,arl,ierr)
	implicit none
	character(*), intent(out):: str
	integer, intent(in):: arl
	character(1), intent(in):: ar1(1:arl)
	integer i,j,k,l,m,n,ierr
	ierr=0
	if(arl.le.len(str)) then
	 do i=1,arl; str(i:i)=ar1(i); enddo
	else
	 ierr=1
	endif
	return
	end subroutine array2string
!-------------------------------------------------------------------------------------------
	SUBROUTINE CAP_ASCII(STR) !makes all small Enlgish letters in the string STR capital
	implicit none
	character(*) STR
	integer i,j,k,l,ia1,ia2,ia3

	ia1=iachar('a')
	ia2=iachar('z')
	ia3=iachar('A')
	l=len_trim(STR)
	do k=1,l
	 i=iachar(STR(k:k))
	 if(i.ge.ia1.and.i.le.ia2) STR(k:k)=achar(ia3+i-ia1)
	enddo
	return
	END SUBROUTINE CAP_ASCII
!--------------------------------------------------------------------------------------------------------
	SUBROUTINE CHARNUM(A,DN,N) !converts a real value given as the string A to real and integer forms
	implicit none              !if A is not a number then: int(DN)=-2, N=0 (not a number tag)
	character(*) A             !line that contains a symbolic number
	integer N                  !N - integer output number
	real(8) DN,DSN,DL,DP,DPS   !DN - real output number
	integer I,J,K,K1,K2,L,M,KF,ISN,LA,IA0,IA9

	LA=len_trim(A)       !the length of the treated string without ending blanks
	IA0=iachar('0')
	IA9=iachar('9')
	DSN=+1d0             !default sign is positive
	ISN=+1
	DP=0d0               !power (for number D and E formats)
	DPS=+1d0
	DN=0d0               !output real number
	N=0                  !output integer number
	M=1                  !searcing for integer part of the number first
	DL=0d0               !special
	KF=0                 !flag
!check a sign of the number:
	K1=1
	do while(K1.lt.LA.and.A(K1:K1).eq.' ')  !skip begining blanks
	 K1=K1+1
	enddo
	if(A(K1:K1).eq.'+') then               !positive number
	 if(K1.eq.LA) goto 999                 !not a number case (just (+) appeared)
	 K1=K1+1
	elseif(A(K1:K1).eq.'-') then           !negative number
	 if(K1.eq.LA) goto 999                 !not a number case (just (-) appeared)
	 ISN=-1
	 DSN=-1d0
	 K1=K1+1
	endif
!process the number itself:
	do K=K1,LA
	 I=iachar(A(K:K))
	 if(A(K:K).EQ.'.'.or.A(K:K).EQ.',') then           !decimal dot/comma detected
	  if(M.eq.1) then
	   M=2
	  else
	   goto 999                                        !not a number case (two decimal dots detected)
	  endif
	 elseif(IA0.le.I.and.I.le.IA9) then                !a decimal digit detected
	  KF=1
	  if(M.eq.1) then
	   DN=DN*10d0+DBLE(I-IA0)                          !integer part
	  else
	   DL=DL-1d0
	   DN=DN+DBLE(I-IA0)*10d0**DL                      !only for the real numbers
	  endif
	 elseif(A(K:K).EQ.'D'.or.A(K:K).EQ.'d'.or.A(K:K).EQ.'E'.or.A(K:K).EQ.'e') then  !power symbol
	  L=K+1
	  if(L.gt.LA) goto 500                             !D means D0 by default (the same applies for d,E and e)
	  J=iachar(A(L:L))
	  if(A(L:L).eq.'+') then
	   L=L+1
	  elseif(A(L:L).eq.'-') then
	   DPS=-1d0
	   L=L+1
	  elseif(J.lt.IA0.or.J.gt.IA9) then
	   goto 999                                        !not a number case
	  endif
	  do K2=L,LA
	   J=iachar(A(K2:K2))
	   if(J.lt.IA0.or.J.gt.IA9) goto 999               !not a number case
	   DP=DP*10d0+DBLE(J-IA0)
	  enddo
	  DP=DPS*DP
	  goto 500
	 elseif(A(K:K).EQ.' ') then                        !blank read
	  if(KF.eq.1) goto 500                             !blanks are ignored along with the characters after them
	 else
	  goto 999                                         !not a number case
	 endif
	enddo   !next character A(K:K)
!number succesfully translated:
500	DN=DSN*DN*10d0**DP
	N=int(DN)                                          !must round by neglecting
	return
!not a number case:
999	DN=-2d0                                            !{DN=-2d0 & N=0}: not a number case
	N=0
	return
	END SUBROUTINE CHARNUM
!------------------------------------------------------------------
	SUBROUTINE CREATE_LINE(talen,tabln,mantl,ctrls,numar,symar)
!The subroutine creates a table line in .txt format with some separator.
!INPUT:
! 1) mantl - number of mantissa digits;
! 2) ctrls - control sequence of commands;
! 3) symar - symbolic array of data;
! 4) numar - real*8 array of data;
!OUTPT:
! 1) tabln - table line (string of characters);
! 2) talen - length of the tabln.
	implicit none
!-------------------------------------------------------------
	character(1), parameter:: sep=';'    !default separator
!-------------------------------------------------------------
	integer i,j,k,l,m,n,k1,k2,k3,k4,kf
	integer, intent(in):: mantl
	integer, intent(out):: talen
	character(*), intent(out):: tabln
	integer,intent(in):: ctrls(*)
	real(8), intent(in):: numar(*)
	character(*), intent(in):: symar(*)
	character(256) os

	talen=0
	l=1
	do
	 select case(ctrls(l))
	 case(0)
	  exit
	 case(1)
	  do k=1,ctrls(l+1)
	   talen=talen+1
	   tabln(talen:talen)=sep
	  enddo
	  l=l+2
	 case(2)
	  do k=ctrls(l+1),ctrls(l+2)
	   m=len_trim(symar(k))
	   tabln(talen+1:talen+m+1)=symar(k)(1:m)//sep
	   talen=talen+m+1
	  enddo
	  l=l+3
	 case(3)
	  do k=ctrls(l+1),ctrls(l+2)
	   call valchar(numar(k),mantl,m,os)
	   tabln(talen+1:talen+m+1)=os(1:m)//sep
	   talen=talen+m+1
	  enddo
	  l=l+3
	 case default
	  write(*,*)'ERROR(STSUBS::CREATE_LINE): invalid control command in ctrls:',ctrls(l),l
	  stop
	 end select
	enddo
	return
	END SUBROUTINE CREATE_LINE
!------------------------------------------------------------------------------------------
	INTEGER FUNCTION ICHARNUM(L,OS)  !OS(1:L) - positive or negative number as a string
	 implicit none                   !L=0 --> error
	 integer L,K,M,N,IA0             !returns ICHARNUM as an integer
	 character(*) OS
	 IA0=iachar('0')
	 if(L.gt.0) then
	  if(OS(1:1).eq.'-') then      !negative
	   M=-2
	  elseif(OS(1:1).eq.'+') then  !explicitly positive
	   M=2
	  else                         !implicitly positive
	   M=1
	  endif
	  ICHARNUM=0
	  do K=abs(M),L
	   N=iachar(OS(K:K))-IA0
	   if(N.ge.0.and.N.le.9) then
	    ICHARNUM=ICHARNUM*10 + N
	   else
	    if(verbose) write(*,*)'ERROR(STSUBS::ICHARNUM): invalid string-number: '//OS(1:L)
	    L=0 !error
	    return
	   endif
	  enddo
	  if(M.eq.-2) ICHARNUM=-ICHARNUM
	 else
	  if(verbose) write(*,*)'ERROR(STSUBS::ICHARNUM): string of non-positive length: ',L
	  L=0 !error
	  return
	 endif
	 return
	END FUNCTION ICHARNUM
!----------------------------------------------------------------------------
	INTEGER FUNCTION IFCL(N)   !returns the factorial of a non-negative N
	 implicit none
	 integer N,K
	 IFCL=1
	  do K=2,N
	   IFCL=IFCL*K
	  enddo
	 return
	END FUNCTION IFCL
!----------------------------------------
	INTEGER FUNCTION IS_IT_NUMBER(ch)
	 implicit none
	 character(1), intent(in):: ch
	 if(iachar(ch).ge.iachar('0').and.iachar(ch).le.iachar('9')) then
	  IS_IT_NUMBER=iachar(ch)-iachar('0')
	 else
	  IS_IT_NUMBER=-1
	 endif
	 return
	END FUNCTION IS_IT_NUMBER
!----------------------------------------
	INTEGER FUNCTION IS_IT_LETTER(ch)
	 implicit none
	 character(1), intent(in):: ch
	 if(iachar(ch).ge.iachar('a').and.iachar(ch).le.iachar('z')) then
	  IS_IT_LETTER=1
	 elseif(iachar(ch).ge.iachar('A').and.iachar(ch).le.iachar('Z')) then
	  IS_IT_LETTER=2
	 else
	  IS_IT_LETTER=0
	 endif
	 return
	END FUNCTION IS_IT_LETTER
!------------------------------------------------------------------------------
	SUBROUTINE ITRSIGN(N,ITR) !ITR(0) - the initial sign of the permutation
	 implicit none            !BUBBLE
	 integer ITR(0:*)
	 integer N,K,L
	 K=1
	 do while(K.lt.N)
	  if(ITR(K).gt.ITR(K+1)) then
	   L=ITR(K)
	   ITR(K)=ITR(K+1)
	   ITR(K+1)=L
	   ITR(0)=-ITR(0)
	   if(K.gt.1) then
	    K=K-1
	   else
	    K=2
	   endif
	  else
	   K=K+1
	  endif
	 enddo
	 return !RETURNS an ascending sequence of indices with a sign in ITR(0)
	END SUBROUTINE ITRSIGN
!-----------------------------------------------------------------
	SUBROUTINE LONGNUMCHAR(I,IOSL,OS) !I - long integer number
	implicit none
	integer(8), intent(in):: I
	integer IOSL,L,K1,K2
	integer(8) K
	character(*) OS                   !OS(1:IOSL) - appropriate string
	character(1) A(0:9),CH
	data A/'0','1','2','3','4','5','6','7','8','9'/

	if(I.lt.0) then
	 OS='-'
	 K=-I
	 IOSL=1
	elseif(I.eq.0) then
	 OS='0'
	 IOSL=1
	 return
	else
	 K=I
	 OS=' '
	 IOSL=0
	endif
	L=IOSL
	do while(K.NE.0)
	 K1=mod(K,10)
	 K=K/10
	 L=L+1
	 OS(L:L)=A(K1)
	enddo
	K1=L-IOSL
	if(mod(K1,2).eq.1) K1=K1-1
	K1=K1/2
	do K2=1,K1
	 CH=OS(IOSL+K2:IOSL+K2)
	 OS(IOSL+K2:IOSL+K2)=OS(L+1-K2:L+1-K2)
	 OS(L+1-K2:L+1-K2)=CH
	enddo
	IOSL=L
	return !OS(1:IOSL) - the number as a string (the sign included if the number is negative)
	END SUBROUTINE LONGNUMCHAR
!-----------------------------------------------------------------------------------------
	SUBROUTINE MARKCHF(STR,N,MF)        !It counts number of fields in the string STR.
	 implicit none                      !The field is a sequence of characters without blanks.
	 character(*) STR                   !Blanks split the fields.
	 integer, intent(out):: N           !N - number of fields.
	 integer, intent(out):: MF(2,1:*)   !MF(1:2,K) - start and finish positions of the field#K.
	 integer K,LSTR

	 LSTR=len_trim(STR)
	 N=0
	 K=1
	 do while(K.LE.LSTR)
	  if(STR(K:K).ne.' '.and.iachar(STR(K:K)).ne.9) then !skip TABS also (ASCII code 9)
	   N=N+1
	   MF(1:2,N)=(/K,0/)
	   do while(K.le.LSTR)
	    if(STR(K:K).eq.' '.or.iachar(STR(K:K)).eq.9) exit
	    K=K+1
	   enddo
	   MF(2,N)=K-1
	  else
	   K=K+1
	  endif
	 enddo
	 return
	END SUBROUTINE MARKCHF
!---------------------------------------------------------------------
	SUBROUTINE MATMAT(sa,sb,m,a,b,c) !Matrix-Matrix Multiplication
!C(sa,sb)=A(sa,m)*B(m,sb)+C(sa,sb)
	 implicit none
	 integer, intent(in):: sa,sb,m
	 integer j,k
	 real(8) a(1:sa,1:m),b(1:m,1:sb),c(1:sa,1:sb)

	 do j=1,sb
	  do k=1,m
	   c(1:sa,j)=c(1:sa,j)+a(1:sa,k)*b(k,j)
	  enddo
	 enddo
	 return
	END SUBROUTINE MATMAT
!----------------------------------------------------------------------
	SUBROUTINE MATMATT(sa,sb,m,a,b,c) !Matrix-Matrix multiplication (transposed 2nd argument)
!C(sa,sb)=A(sa,m)*B(sb,m)+C(sa,sb)
	 implicit none
	 integer, intent(in):: sa,sb,m
	 integer j,k
	 real(8) a(1:sa,1:m),b(1:sb,1:m),c(1:sa,1:sb)

	 do j=1,sb
	  do k=1,m
	   c(1:sa,j)=c(1:sa,j)+a(1:sa,k)*b(j,k)
	  enddo
	 enddo
	 return
	END SUBROUTINE MATMATT
!----------------------------------------------------------------------
	SUBROUTINE MATTMAT(sa,sb,m,a,b,c) !Matrix-Matrix multiplication (transposed 1st argument)
!C(sa,sb)=A(m,sa)*B(m,sb)+C(sa,sb)
	 implicit none
	 integer, intent(in):: sa,sb,m
	 integer j,k
	 real(8) a(1:m,1:sa),b(1:m,1:sb),c(1:sa,1:sb)

	 do j=1,sb
	  do k=1,m
	   c(1:sa,j)=c(1:sa,j)+a(k,1:sa)*b(k,j)
	  enddo
	 enddo
	 return
	END SUBROUTINE MATTMAT
!-----------------------------------------------------------------------
	SUBROUTINE MATTMATT(sa,sb,m,a,b,c) !Matrix-Matrix multiplication (both arguments are transposed)
!C(sa,sb)=A(m,sa)*B(sb,m)+C(sa,sb)
	 implicit none
	 integer, intent(in):: sa,sb,m
	 integer j,k
	 real(8) a(1:m,1:sa),b(1:sb,1:m),c(1:sa,1:sb)

	 do j=1,sb
	  do k=1,m
	   c(1:sa,j)=c(1:sa,j)+a(k,1:sa)*b(j,k)
	  enddo
	 enddo
	 return
	END SUBROUTINE MATTMATT
!----------------------------------------------------------------------
        integer function NAME_HASH(max_hash,sl,str)
!This function returns a hash mask (0..max_hash) for the string passed.
!INPUT:
! - max_hash - limit of the hash range (parameter);
! - str(1:sl) - string;
!OUTPUT:
! - name_hash - hash mask [0..max_hash].
        implicit none
        integer, intent(in):: max_hash,sl
        character(*), intent(in):: str
        integer i

        NAME_HASH=0
        do i=1,sl
         NAME_HASH=NAME_HASH+iachar(str(i:i))
        enddo
        NAME_HASH=mod(NAME_HASH,max_hash+1)
        return
        end function NAME_HASH
!-----------------------------------------------------------------------
	subroutine NOCOMMENT(l,str)
!This subroutine removes comments (!,#) from the line given in str(1:l).
!After that l is a modified length.
	implicit none
	integer, intent(inout):: l
	character(*), intent(in):: str
	integer i

	do i=1,l
	 if(str(i:i).eq.'!'.or.str(i:i).eq.'#') then
	  l=i-1
	  exit
	 endif
	enddo
	return
	end subroutine NOCOMMENT
!------------------------------------------------------------------
	SUBROUTINE NORMV(n,vn,v)
!Substitutes the original vector v with its normalized counterpart.
!vn - is an original Euclide length of the vector.
!n - dimension of the vector.
	 implicit none
	 integer, intent(in):: n
	 integer i
	 real(8) vn,v(1:n)

	 vn=0d0
	 do i=1,n
	  vn=vn+v(i)**2
	 enddo
	 if(vn.ne.0d0) then
	  vn=dsqrt(vn)
	  v(1:n)=v(1:n)/vn
	 endif
	 return
	END SUBROUTINE NORMV
!--------------------------------------------------------------
	SUBROUTINE NOSPACE(A,L) !A - a string with a length = L
	 implicit none          !return A without SPACES and TABS. L changes to a reduced length
	 character(*) A
	 integer L,J,M

	 M=len_trim(A)
	 L=0
	 do J=1,M
	  if(A(J:J).NE.' '.and.iachar(A(J:J)).NE.9) then
	   L=L+1
	   A(L:L)=A(J:J)
	  endif
	 enddo
	 return
	END SUBROUTINE NOSPACE
!-----------------------------------------
	LOGICAL FUNCTION NOT_A_NUMBER(str)
	 implicit none
	 character(*), intent(in):: str
	 integer dot_present,exp_present,i,l
	 l=len_trim(str)
	 if(l.gt.0) then
	  NOT_A_NUMBER=.false.; dot_present=0; exp_present=0
	  if(str(1:1).eq.'-'.or.str(1:1).eq.'+') then; i=2; else; i=1; endif
	  do while(i.le.l)
	   if(iachar(str(i:i)).lt.iachar('0').or.iachar(str(i:i)).gt.iachar('9')) then
	    if(str(i:i).eq.'.') then
	     if(dot_present.eq.0) then; dot_present=1; else; NOT_A_NUMBER=.true.; return; endif
	    elseif(str(i:i).eq.'D'.or.str(i:i).eq.'d'.or.str(i:i).eq.'E'.or.str(i:i).eq.'e') then
	     if(exp_present.eq.0) then
	      exp_present=1
	      if(i.lt.l) then; if(str(i+1:i+1).eq.'-'.or.str(i+1:i+1).eq.'+') i=i+1; else; NOT_A_NUMBER=.true.; return; endif
	     else
	      NOT_A_NUMBER=.true.; return
	     endif
	    else
	     NOT_A_NUMBER=.true.; return
	    endif
	    i=i+1
	   else
	    i=i+1
	   endif
	  enddo
	 else !empty string is not a number
	  NOT_A_NUMBER=.true.
	 endif
	 return
	END FUNCTION NOT_A_NUMBER
!--------------------------------------------------------
	SUBROUTINE NUMCHAR(I,IOSL,OS) !I - integer number
	implicit none
	integer, intent(in):: I
	integer IOSL,K,K1,L
	character(*) OS               !OS(1:IOSL) - appropriate string
	character(1) A(0:9),CH
	data A/'0','1','2','3','4','5','6','7','8','9'/

	if(I.lt.0) then
	 OS='-'
	 K=-I
	 IOSL=1
	elseif(I.eq.0) then
	 OS='0'
	 IOSL=1
	 return
	else
	 K=I
	 OS=' '
	 IOSL=0
	endif
	L=IOSL
	do while(K.NE.0)
	 K1=mod(K,10)
	 K=K/10
	 L=L+1
	 OS(L:L)=A(K1)
	enddo
	K1=L-IOSL
	if(mod(K1,2).eq.1) K1=K1-1
	K1=K1/2
	do K=1,K1
	 CH=OS(IOSL+K:IOSL+K)
	 OS(IOSL+K:IOSL+K)=OS(L+1-K:L+1-K)
	 OS(L+1-K:L+1-K)=CH
	enddo
	IOSL=L
	return !OS(1:IOSL) - the number as a string (the sign included if the number is negative)
	END SUBROUTINE NUMCHAR
!---------------------------------------------------------------------------
	SUBROUTINE PRINTL(FH,STR,ADV) !prints a line str(1:*) into a file#FH
	implicit none
	character(32) FRM
	character(*) STR
	integer, intent(in):: FH
	logical, intent(in), optional:: ADV !carriage return control (.false. - no carriage return)
	integer l,k

	l=len_trim(STR)
	if(l.gt.0) then
	 FRM(1:2)='(A'
	 call NUMCHAR(l,k,FRM(3:32))
	 k=3+k; FRM(k:k)=')'
	 if(present(ADV)) then
	  if(ADV) then
	   write(FH,FRM(1:k)) STR(1:l)
	  else
	   write(FH,FRM(1:k),advance='no') STR(1:l)
	  endif
	 else
	  write(FH,FRM(1:k)) STR(1:l)
	 endif
	else
	 if(present(ADV)) then
	  if(ADV) write(FH,*)
	 endif
	endif
	return
	END SUBROUTINE PRINTL
!----------------------------------
	SUBROUTINE RAND_STR(sl,str)
	implicit none
	integer, intent(out):: sl
	character(*), intent(out):: str
	integer, parameter:: max_str_len=8192
	integer, parameter:: buf_size=512
	integer, parameter:: ascii_b=32,ascii_e=126,ascii_range=ascii_e-ascii_b+1
	real(8) rnds(1:buf_size),rn
	integer i,j,k,m,n

	if(len(str).gt.0) then
	 k=min(len(str),max_str_len)
	 call random_number(rn); sl=int(rn*dble(k))+1; if(sl.gt.k) sl=k
	 do i=0,sl-1,buf_size
	  n=min(sl-i,buf_size)
	  call random_number(rnds(1:n))
	  do j=1,n
	   m=int(rnds(j)*dble(ascii_range))+ascii_b; if(m.gt.ascii_e) m=ascii_e
	   str(i+j:i+j)=achar(m)
	  enddo
	 enddo
	else
	 sl=0
	endif
	return
	END SUBROUTINE RAND_STR
!---------------------------------------------------------
	SUBROUTINE ROTS(axis,alpha,n,pnts)
!This subroutine rotates n points placed in the array pnts
!around the axis #axis in a 3d-space.
	implicit none
	integer i,j,k,l,m,n,k1,k2,k3,k4,ks,kf
	integer axis
	real(8) alpha,pnts(3,*),q(3),cosa,sina

	sina=sin(alpha)
	cosa=cos(alpha)
	select case(axis)
	case(1)
!X axis:
	 do k=1,n
	  q(1:3)=pnts(1:3,k)
	  pnts(1,k)=q(1)
	  pnts(2,k)=q(2)*cosa-q(3)*sina
	  pnts(3,k)=q(2)*sina+q(3)*cosa
	 enddo
	case(2)
!Y axis:
	 do k=1,n
	  q(1:3)=pnts(1:3,k)
	  pnts(1,k)=q(1)*cosa-q(3)*sina
	  pnts(2,k)=q(2)
	  pnts(3,k)=q(1)*sina+q(3)*cosa
	 enddo
	case(3)
!Z axis:
	 do k=1,n
	  q(1:3)=pnts(1:3,k)
	  pnts(1,k)=q(1)*cosa-q(2)*sina
	  pnts(2,k)=q(1)*sina+q(2)*cosa
	  pnts(3,k)=q(3)
	 enddo
	case default
	 write(*,*)'ERROR(rots): invalid case encountered: ',axis
	 stop
	end select
	return
	END SUBROUTINE ROTS
!---------------------------------------------------------------------------------------------
	SUBROUTINE SMALL_ASCII(STR) !makes all capital Enlgish letters in the string STR small
	implicit none
	character(*) STR
	integer i,j,k,l,m,n,ia1,ia2,ia3

	ia1=iachar('A')
	ia2=iachar('Z')
	ia3=iachar('a')
	l=len_trim(STR)
	do k=1,l
	 i=iachar(STR(k:k))
	 if(i.ge.ia1.and.i.le.ia2) STR(k:k)=achar(ia3+i-ia1)
	enddo
	return
	END SUBROUTINE SMALL_ASCII
!------------------------------------------------
	subroutine string2array(str,ar1,arl,ierr)
	implicit none
	character(*), intent(in):: str
	character(1), intent(out):: ar1(1:)
	integer, intent(out):: arl
	integer i,j,k,l,m,n,ierr
	ierr=0; arl=len(str)
	if(size(ar1).ge.arl) then
	 do i=1,arl; ar1(i)=str(i:i); enddo
	else
	 ierr=1
	endif
	return
	end subroutine string2array
!---------------------------------------------------------------------------------------------------
	SUBROUTINE STRSEARCH(L,STR,FRG)   !searches exactly the first FRG fragment in the STR string
	 implicit none                    !L - initially the length of STR.
	 integer L                        !if found L = position of FRG in STR
	 character(*) STR,FRG             !if not   L = -1
	 character(1) A
	 integer LFRG,LSTR,K,M

	 LFRG=len(FRG) !spaces are taken into account!
	 if(LFRG.le.0.or.L.le.0) then
	  L=-1
	  return
	 endif
	 LSTR=L
	 M=1
	 A=FRG(1:1)
	 do K=1,LSTR
	  if(STR(K:K).EQ.A) then
	   if(M.eq.LFRG) then
	    L=K+1-LFRG            !the 1st position of the 1st FRG fragment in STR
	    return
	   else
	    M=M+1
	    A=FRG(M:M)
	   endif
	  else
	   M=1
	   A=FRG(1:1)
	  endif
	 enddo
	 L=-1                     !FRG has not been found in STR
	 return
	END SUBROUTINE STRSEARCH
!---------------------------------------
	INTEGER FUNCTION SYMBOL_TYPE(ch)
	 implicit none
	 character(1), intent(in):: ch
	 integer i
	 i=iachar(ch)
	 if(i.eq.9.or.i.eq.32) then !tab or space
	  SYMBOL_TYPE=0
	 elseif(i.ge.iachar('0').and.i.le.iachar('9')) then !number
	  SYMBOL_TYPE=1
	 elseif((i.ge.iachar('a').and.i.le.iachar('z')).or.(i.ge.iachar('A').and.i.le.iachar('Z'))) then
	  SYMBOL_TYPE=2
	 else !other
	  SYMBOL_TYPE=-1
	 endif
	 return
	END FUNCTION SYMBOL_TYPE
!-------------------------------
	SUBROUTINE TPAUSE(ISEC)
	 implicit none
	 integer ISEC
	 real(4) TIM1,TIM2,TIMP
	 TIMP=real(ISEC)
	 call cpu_time(TIM1)
	 call cpu_time(TIM2)
	 do while(TIM2-TIM1.lt.TIMP)
	  call cpu_time(TIM2)
	 enddo
	 return
	END SUBROUTINE TPAUSE
!-----------------------------------------
	SUBROUTINE VALCHAR(val,mantl,l,os)
	 implicit none
	 character(*) os
	 character(1) A(0:9),ch
	 real(8), intent(in):: val
	 integer, intent(in):: mantl
	 integer l,k,m,ls
	 real(8) d
	 data A/'0','1','2','3','4','5','6','7','8','9'/

	 l=0
	 if(val.ge.0d0) then
	  d=val
	 else
	  d=-val
	  os(1:1)='-'
	  l=1
	 endif
	 ls=l+1
	 k=int(d)
	 if(k.gt.0) then
	  do while(k.gt.0)
	   m=mod(k,10)
	   k=k/10
	   l=l+1
	   os(l:l)=A(m)
	  enddo
	  k=0
	  do while(ls+k.lt.l-k)
	   ch=os(ls+k:ls+k)
	   os(ls+k:ls+k)=os(l-k:l-k)
	   os(l-k:l-k)=ch
	   k=k+1
	  enddo
	  l=l+1
	  os(l:l)='.'
	 else
	  os(l+1:l+2)='0.'
	  l=l+2
	 endif
	 d=d-dint(d)
	 k=mantl                       !max length of the mantissa
	 do while(d.gt.0d0.and.k.gt.0)
	  d=d*10d0
	  m=mod(int(d),10)
	  l=l+1
	  os(l:l)=A(m)
	  d=d-dint(d)
	  k=k-1
	 enddo
	 return
	END SUBROUTINE VALCHAR
!---------------------------------
	subroutine WAIT_DELAY(sec)
	implicit none
	real(4), intent(in):: sec
	real(4) sec0,sec1
	call cpu_time(sec0)
	call cpu_time(sec1)
	do while(sec1.lt.sec0+sec)
	 call cpu_time(sec1)
	enddo
	return
	end subroutine WAIT_DELAY
!---------------------------------
	subroutine WAIT_PRESS(msg)
	implicit none
	character(*), intent(in), optional:: msg
	character(1) ch
	if(present(msg)) then
	 call printl(6,msg(1:len_trim(msg)),.false.)
	else
	 write(6,'("Press ENTER to continue ...")',advance='no')
	endif
	read(*,'(A1)') ch
	return
	end subroutine WAIT_PRESS
!----------------------------------
	subroutine WR_MAT_IN(m,n,a)
!This subroutine writes the matrix a(m,n) to the screen.
!The elements of matrix are integers.
	 implicit none
	 integer m,n,i,j,a(1:m,1:n)
	 do i=1,m
	  do j=1,n
	   write(*,'((I10,1x))',advance='no') a(i,j)
	  enddo
	  write(*,*)""
	 enddo
	 return
	end subroutine WR_MAT_IN
!-----------------------------------
	subroutine WR_MAT_IN8(m,n,a)
!This subroutine writes the matrix a(m,n) to the screen.
!The elements of matrix are integer8.
	 implicit none
	 integer m,n,i,j
	 integer(8) a(1:m,1:n)
	 do i=1,m
	  do j=1,n
	   write(*,'((I20,1x))',advance='no') a(i,j)
	  enddo
	  write(*,*)""
	 enddo
	 return
	end subroutine WR_MAT_IN8
!----------------------------------
	subroutine WR_MAT_SP(m,n,a)
!This subroutine writes the matrix a(m,n) to screen.
!The elements of matrix are of single precision real type.
	 implicit none
	 integer m,n,i,j
	 real(4) a(1:m,1:n)
	 do i=1,m
	  do j=1,n
	   write(*,'((F15.7,1x))',advance='no') a(i,j)
	  enddo
	  write(*,*)""
	 enddo
	 return
	end subroutine WR_MAT_SP
!----------------------------------
	subroutine WR_MAT_DP(m,n,a)
!This subroutine writes the matrix a(m,n) to screen.
!The elements of matrix are of double precision real type.
	 implicit none
	 integer m,n,i,j
	 real(8) a(1:m,1:n)
	 do i=1,m
	  do j=1,n
	   write(*,'(D22.14,1x)',advance='no') a(i,j)
	  enddo
	  write(*,*)""
	 enddo
	 return
	end subroutine WR_MAT_DP
!----------------------------------
	subroutine WR_MAT_DC(m,n,a)
!This subroutine writes the matrix a(m,n) to the screen.
!The elements of matrix are of double complex type.
	 implicit none
	 integer m,n,i,j
	 complex(8) a(1:m,1:n)
	 do i=1,m
	  do j=1,n
	   write(*,'(("(",D22.14,",",D22.14,")"))',advance='no') a(i,j)
	  enddo
	  write(*,*)""
	 enddo
	 return
	end subroutine WR_MAT_DC
!--------------------------------
	subroutine WR_VEC_SP(m,a)
	 implicit none
	 integer m,i
	 real(4) a(1:m)
	 do i=1,m
	  write(*,"(F15.7)") a(i)
	 enddo
	 return
	end subroutine WR_VEC_SP
!--------------------------------
	subroutine WR_VEC_DP(m,a)
	 implicit none
	 integer m,i
	 real(8) a(1:m)
	 do i=1,m
	  write(*,"(D22.14)") a(i)
	 enddo
	 return
	end subroutine WR_VEC_DP

	END MODULE STSUBS
