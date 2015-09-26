       module dil_kinds
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
       end module dil_kinds
