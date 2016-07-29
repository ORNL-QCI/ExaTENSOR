       program test_components
        use pack_prim_test
        use dil_basic, only: INTD,INTL
        integer:: ierr
        integer(INTD):: errc

!Packing primitives:
        write(*,'("Testing packing primitives ... ")',ADVANCE='NO')
        ierr=test_pack_prim(errc)
        if(ierr.eq.0) then; write(*,'("Passed")'); else; write(*,'("Failed: Error codes ",i11," -> ",i11)') ierr,errc; endif

        stop
       end program test_components
