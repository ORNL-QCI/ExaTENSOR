!Generic Fortran Containers (GFC): Hash map (unordered map)
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com, liakhdi@ornl.gov
!REVISION: 2016/11/03

!Copyright (C) 2014-2022 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2022 Oak Ridge National Laboratory (UT-Battelle)

!LICENSE: BSD 3-Clause

       module gfc_hash_map
        use gfc_base
        use timers
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !output device
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
        integer(INTD), private:: DEBUG=0    !debugging level (0:none)
!TYPES:
 !

       end module gfc_hash_map
