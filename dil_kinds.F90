       module dil_kinds
        use, intrinsic:: ISO_C_BINDING
        implicit none
        public
!BASIC TYPE KINDS:
        integer(C_INT), parameter, public:: INTD=4   !default integer size
        integer(C_INT), parameter, public:: INTL=8   !long integer size
        integer(C_INT), parameter, public:: REALD=8  !default real size
        integer(C_INT), parameter, public:: REALL=16 !long real size
#ifndef NO_PHI
!DIR$ ATTRIBUTES OFFLOAD:mic:: INTD,INTL,REALD,REALL
!DIR$ ATTRIBUTES ALIGN:128:: INTD,INTL,REALD,REALL
#endif

!TENSOR DATA KINDS (keep consistent with tensor_algebra.h):
        integer(C_INT), parameter, public:: NO_TYPE=0 !no type/kind
        integer(C_INT), parameter, public:: R4=4      !float data kind
        integer(C_INT), parameter, public:: R8=8      !double data kind
        integer(C_INT), parameter, public:: C8=16     !double complex data kind
        real(4), parameter, public:: R4_=0.0
        real(8), parameter, public:: R8_=0d0
        complex(8), parameter, public:: C8_=(0d0,0d0)
#ifndef NO_PHI
!DIR$ ATTRIBUTES OFFLOAD:mic:: NO_TYPE,R4,R8,C8,R4_,R8_,C8_
!DIR$ ATTRIBUTES ALIGN:128:: NO_TYPE,R4,R8,C8,R4_,R8_,C8_
#endif

!ALIASES (keep consistent with tensor_algebra.h):
        integer(C_INT), parameter, public:: NOPE=0                      !"NO" answer
        integer(C_INT), parameter, public:: YEP=1                       !"YES" answer
        integer(C_INT), parameter, public:: TRY_LATER=-918273645        !try the action later (resources are currently busy): KEEP THIS UNIQUE!
        integer(C_INT), parameter, public:: NOT_CLEAN=-192837465        !something, like resource release, did not go right, but you can continue: KEEP THIS UNIQUE!
#ifndef NO_PHI
!DIR$ ATTRIBUTES OFFLOAD:mic:: NOPE,YEP,TRY_LATER,NOT_CLEAN
!DIR$ ATTRIBUTES ALIGN:128:: NOPE,YEP,TRY_LATER,NOT_CLEAN
#endif
       end module dil_kinds
