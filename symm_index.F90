!This module provides infrastructure for symmetric multi-indices.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2014/07/29
!FUNCTIONALITY:
! # i:get_address_table(ndim,lb,ub,ord,handle,iba): obtain the address table;
! # i:multiindex(ndim,mlndx,iba): get a linear offset for a given multiindex;
       module symm_index
!PARAMETERS:
        integer, private:: cons_out=6
        logical, private:: verbose=.true.
        integer, parameter, private:: max_mlndx_length=256
        integer, parameter, private:: max_banks=16
        integer, parameter, private:: tables_per_bank=16384
        integer, parameter, private:: max_addr_tables=max_banks*tables_per_bank
!ALIASES:
        integer, parameter, public:: SYMM_INDEX_EMPTY=-1
        integer, parameter, public:: SYMM_INDEX_NO_ORDER=0
        integer, parameter, public:: SYMM_INDEX_LE_ORDER=1
        integer, parameter, public:: SYMM_INDEX_GE_ORDER=2
        integer, parameter, public:: SYMM_INDEX_LT_ORDER=3
        integer, parameter, public:: SYMM_INDEX_GT_ORDER=4
!TYPES:
 !Address table:
        type, private:: address_table_t
         integer, allocatable, private:: lbounds(:)   !lower index bounds
         integer, allocatable, private:: ubounds(:)   !upper index bounds
         integer, private:: ordering=SYMM_INDEX_EMPTY !type of ordering
         integer, allocatable, private:: incr(:,:)    !table of addressing increaments
        end type address_table_t
 !Bank of address tables:
        type, private:: tables_bank_t
         integer, private:: num_tables=0                           !number of live tables in the bank
         integer, allocatable, private:: free_tab(:)               !free table handles
         type(address_table_t), allocatable, private:: addr_tab(:) !table bank
        end type tables_bank_t
 !Address tables container:
        type, private:: address_tables_t
         integer, private:: banks_in_use=0                    !number of table banks in use
         integer, private:: tables_in_use=0                   !total number of tables in use
         integer, private:: total_table_size=0                !total size of all allocated tables (in integers)
         type(tables_bank_t), private:: tab_bank(1:max_banks) !table banks
        end type address_tables_t
!DATA:
        type(address_tables_t), private:: address_tables
!FUNCTIONS:
        public address_tables_clean
        public delete_address_table
        public get_address_table
       contains
!----------------------------------------
        subroutine address_tables_clean()
        implicit none
        integer j
        do j=1,max_banks
         if(allocated(address_tables%tab_bank(j)%free_tab) deallocate(address_tables%tab_bank(j)%free_tab)
         if(allocated(address_tables%tab_bank(j)%addr_tab) deallocate(address_tables%tab_bank(j)%addr_tab)
         address_tables%tab_bank(j)%num_tables=0
        enddo
        address_tables%banks_in_use=0; address_tables%tables_in_use=0; address_tables%total_table_size=0
        return
        end subroutine address_tables_clean
!------------------------------------------------------------------------
        integer function delete_address_table(handle)
        implicit none
        integer, intent(inout):: handle
        integer i,n
        delete_address_table=0
        if(handle.gt.0.and.handle.le.max_addr_tables) then
         n=((handle-1)/tables_per_bank)+1; i=handle-(n-1)*tables_per_bank
         if(allocated(address_tables%tab_bank(n)%addr_tab.and. &
            allocated(address_tables%tab_bank(n)%free_tab) then
          if(allocated(address_tables%tab_bank(n)%addr_tab(i)%lbounds) &
           deallocate(address_tables%tab_bank(n)%addr_tab(i)%lbounds)
          if(allocated(address_tables%tab_bank(n)%addr_tab(i)%ubounds) &
           deallocate(address_tables%tab_bank(n)%addr_tab(i)%ubounds)
          if(allocated(address_tables%tab_bank(n)%addr_tab(i)%incr) then
           address_tables%total_table_size=address_tables%total_table_size-size(address_tables%tab_bank(n)%addr_tab(i)%incr)
           deallocate(address_tables%tab_bank(n)%addr_tab(i)%incr)
          endif
          address_tables%tab_bank(n)%addr_tab(i)%ordering=SYMM_INDEX_EMPTY
          address_tables%tables_in_use=address_tables%tables_in_use-1
          address_tables%tab_bank(n)%free_tab(address_tables%tab_bank(n)%num_tables)=i
          address_tables%tab_bank(n)%num_tables=address_tables%tab_bank(n)%num_tables-1
          if(address_tables%tab_bank(n)%num_tables.eq.0) then
           deallocate(address_tables%tab_bank(n)%addr_tab,address_tables%tab_bank(n)%free_tab)
           address_tables%banks_in_use=address_tables%banks_in_use-1
          endif
          handle=0
         else
          delete_address_table=-1
         endif
        else
         delete_address_table=-2
        endif
        return
        end function delete_address_table
!------------------------------------------------------------------------
        integer function get_address_table(handle,iba,ndim,lbnd,ubnd,ord)
        implicit none
        integer, intent(inout):: handle            !address table handle (0: empty: to be returned)
        integer, pointer, intent(out):: iba(:,:)   !address table
        integer, intent(in), optional:: ndim       !length of the multi-index
        integer, intent(in), optional:: lbnd(1:*)  !lower bounds of index ranges
        integer, intent(in), optional:: ubnd(1:*)  !upper bounds of index ranges
        integer, intent(in), optional:: ord        !index ordering pattern
        integer i,j,k,l,m,n,ierr
        integer lb(1:max_mlndx_length),ub(1:max_mlndx_length),bank,tab

        get_address_table=0
        if(handle.gt.0.and.handle.le.max_addr_tables) then !handle --> existing table
         n=((handle-1)/tables_per_bank)+1 !bank number
         if(allocated(address_tables%tab_bank(n)%addr_tab) then
          i=handle-(n-1)*tables_per_bank !entry number
          if(allocated(address_tables%tab_bank(n)%addr_tab(i)%incr) then
           iba=>address_tables%tab_bank(n)%addr_tab(i)%incr
          else
           get_address_table=-1
          endif
         else
          get_address_table=-2
         endif
        elseif(handle.le.0) then !empty handle: create a new table
         if(present(ndim).and.present(lb).and.present(ub).and.present(ord)) then
          if(ndim.gt.0.and.ndim.le.max_mlndx_length) then
           do i=1,ndim; if(lbnd(i).gt.ubnd(i)) then; get_address_table=-3; return; endif; enddo !check
           max_range=0; do i=1,ndim; max_range=max(max_range,ubnd(i)-lbnd(i)+1); enddo !get the max index range
           select case(ord)
           case(SYMM_INDEX_NO_ORDER)
            
           case(SYMM_INDEX_LE_ORDER)            
            do i=1,ndim; lb(i)=lbnd(i)-lbnd(1); enddo; do i=1,ndim; ub(i)=ubnd(i)-lbnd(1); enddo !normalize
            do i=1,ndim-1; if(lb(i).gt.lb(i+1)) then; get_address_table=-1; return; endif; enddo !check
            do i=1,ndim-1; if(ub(i).gt.ub(i+1)) then; get_address_table=-1; return; endif; enddo !check
            ierr=get_free_table(bank,tab,handle); if(ierr.ne.0) then; get_address_table=-1; return; endif
            iba=>address_tables%tab_bank(bank)%addr_tab(tab)%incr
            do i=0,ub(1)-lb(1); iba(i,1)=i; enddo
            do m=2,ndim
             l=0
             do i=lb(m),ub(m)
              j=i-lb(m); iba(j,m)=l
              do k=1,m-1; l=l+iba(min(ub(k),i)-lb(k),k); enddo; l=l+1
             enddo
            enddo
           case(SYMM_INDEX_GE_ORDER)
            
           case(SYMM_INDEX_LT_ORDER)
            
           case(SYMM_INDEX_GT_ORDER)
            
           case default
            get_address_table=-4
           end select
          else
           get_address_table=-5
          endif
         else
          get_address_table=-6
         endif
        else
         get_address_table=-7
        endif
        return

        contains

         integer function get_free_table(jb,jt,hndl)
         implicit none
         integer, intent(out):: jb,jt,hndl
         integer:: j0,je
         get_free_table=0
         if(address_tables%tables_in_use.lt.max_addr_tables) then
          jb=0; jt=0; hndl=0
          do j0=1,max_banks
           if(allocated(address_tables%tab_bank(j0)%addr_tab) then
            if(address_tables%tab_bank(j0)%num_tables.lt.tables_per_bank) then
             address_tables%tab_bank(j0)%num_tables=address_tables%tab_bank(j0)%num_tables+1
             jt=address_tables%tab_bank(j0)%free_tab(address_tables%tab_bank(j0)%num_tables); jb=j0
             hndl=jt+(jb-1)*tables_per_bank
             exit
            endif
           else
            allocate(address_tables%tab_bank(j0)%addr_tab(1:tables_per_bank), &
                     address_tables%tab_bank(j0)%free_tab(1:tables_per_bank),STAT=je)
            if(je.ne.0) then; get_free_table=-1; return; endif
            do je=1,tables_per_bank; address_tables%tab_bank(j0)%free_tab(je)=je; enddo
            address_tables%banks_in_use=address_tables%banks_in_use+1
            address_tables%tab_bank(j0)%num_tables=1; jb=j0; jt=1; hndl=jt+(jb-1)*tables_per_bank
            exit
           endif
          enddo
          if(jb.gt.0.and.gt.0) then
           allocate(address_tables%tab_bank(jb)%addr_tab(jt)%lbounds(1:ndim), &
                    address_tables%tab_bank(jb)%addr_tab(jt)%ubounds(1:ndim), &
                    address_tables%tab_bank(jb)%addr_tab(jt)%incr(0:max_range-1,1:ndim),STAT=je)
           if(je.eq.0) then
            address_tables%tab_bank(jb)%addr_tab(jt)%ordering=ord
            address_tables%tab_bank(jb)%addr_tab(jt)%lbounds(1:ndim)=lb(1:ndim)
            address_tables%tab_bank(jb)%addr_tab(jt)%ubounds(1:ndim)=ub(1:ndim)
            address_tables%tables_in_use=address_tables%tables_in_use+1
            address_tables%total_table_size=address_tables%total_table_size+ndim*max_range
           else
            if(allocated(address_tables%tab_bank(jb)%addr_tab(jt)%lbounds)) &
             deallocate(address_tables%tab_bank(jb)%addr_tab(jt)%lbounds)
            if(allocated(address_tables%tab_bank(jb)%addr_tab(jt)%ubounds)) &
             deallocate(address_tables%tab_bank(jb)%addr_tab(jt)%ubounds)
            if(allocated(address_tables%tab_bank(jb)%addr_tab(jt)%incr)) &
             deallocate(address_tables%tab_bank(jb)%addr_tab(jt)%incr)
            get_free_table=-2
           endif
          else
           get_free_table=1
          endif
         else !no free tables left
          get_free_table=2
         endif
         return
         end function get_free_table

        end function get_address_table

       end module symm_index
