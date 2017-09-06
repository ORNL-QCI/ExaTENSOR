!ExaTENSOR: TAVP-Manager (TAVP-MNG) implementation
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/09/06

!Copyright (C) 2014-2017 Dmitry I. Lyakh (Liakh)
!Copyright (C) 2014-2017 Oak Ridge National Laboratory (UT-Battelle)

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

       module tavp_manager
        use virta
        use gfc_base
        use gfc_list
        use gfc_dictionary
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !default output device
        integer(INTD), private:: DEBUG=0    !debugging mode
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
!TYPES:
 !Tensor argument cache entry (TAVP-specific):
        type, extends(tens_cache_entry_t), private:: tens_entry_mng_t
         integer(INTD), private:: owner_id=-1                      !tensor owner id (non-negative)
         contains
          procedure, private:: TensEntryMngCtor                    !ctor
          generic, public:: tens_entry_mng_ctor=>TensEntryMngCtor
          procedure, public:: get_owner_id=>TensEntryMngGetOwnerId !returns the owner id
          final:: tens_entry_mng_dtor
        end type tens_entry_mng_t
 !Tensor operand (encapsulated tensor data processible by a specific TAVP):
        type, extends(ds_oprnd_t), private:: tens_oprnd_t
         class(tens_rcrsv_t), pointer, private:: tensor=>NULL() !non-owning pointer to a persistent recursive tensor
         integer(INTD), private:: owner_id=-1                   !tensor owner id (non-negative)
         contains
          procedure, private:: TensOprndCtor                    !ctor
          generic, public:: tens_oprnd_ctor=>TensOprndCtor
          procedure, public:: get_tensor=>TensOprndGetTensor    !returns a pointer to the tensor
          procedure, public:: get_owner_id=>TensOprndGetOwnerId !returns the tensor owner id
          procedure, public:: set_owner_id=>TensOprndSetOwnerId !sets the tensor owner id
          procedure, public:: is_remote=>TensOprndIsRemote      !returns TRUE if the tensor operand is remote
          procedure, public:: acquire_rsc=>TensOprndAcquireRsc  !explicitly acquires local resources for the tensor operand
          procedure, public:: prefetch=>TensOprndPrefetch       !starts prefetching the remote tensor operand (acquires local resources!)
          procedure, public:: upload=>TensOprndUpload           !starts uploading the tensor operand to its remote location
          procedure, public:: sync=>TensOprndSync               !synchronizes the currently pending communication on the tensor operand
          procedure, public:: release=>TensOprndRelease         !destroys the present local copy of the tensor operand (releases local resources!), but the operand stays defined
          procedure, public:: destruct=>TensOprndDestruct       !performs complete destruction back to an empty state
          final:: tens_oprnd_dtor
        end type tens_oprnd_t
!MODULE DATA:
 !TAVP-MNG microcode (static) table, set by dsvp.configure():
        type(ds_microcode_t), target, private:: microcode(0:TAVP_ISA_SIZE-1)
!VISIBILITY:
 !tens_entry_mng_t:
        private TensEntryMngCtor
        private TensEntryMngGetOwnerId
        public tens_entry_mng_dtor
        private tens_entry_mng_alloc
 !tens_oprnd_t:
        private TensOprndCtor
        private TensOprndGetTensor
        private TensOprndGetOwnerId
        private TensOprndSetOwnerId
        private TensOprndIsRemote
        private TensOprndAcquireRsc
        private TensOprndPrefetch
        private TensOprndUpload
        private TensOprndSync
        private TensOprndRelease
        private TensOprndDestruct
        public tens_oprnd_dtor
!IMPLEMENTATION:
       contains
![tens_entry_mng_t]========================================
        subroutine TensEntryMngCtor(this,tensor,owner,ierr)
!Constructs a <tens_entry_mng_t>. Note move semantics for <tensor>!
         implicit none
         class(tens_entry_mng_t), intent(out):: this          !out: specialized tensor cache entry
         class(tens_rcrsv_t), pointer, intent(inout):: tensor !inout: pointer to an allocated tensor (ownership transfer will occur here!)
         integer(INTD), intent(in):: owner                    !in: tensor owner id
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTD):: errc

         errc=0
         if(associated(tensor)) then
          if(owner.ge.0) then
           call this%set_tensor(tensor,errc)
           if(errc.eq.0) then
            tensor=>NULL() !transfer the ownership
            this%owner_id=owner
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensEntryMngCtor
!------------------------------------------------------------
        function TensEntryMngGetOwnerId(this,ierr) result(id)
         implicit none
         integer(INTD):: id                          !out: tensor owner id
         class(tens_entry_mng_t), intent(in):: this  !in: specialized tensor cache entry
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%is_set(errc)) then
          id=this%owner_id
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensEntryMngGetOwnerId
!-------------------------------------------
        subroutine tens_entry_mng_dtor(this)
         implicit none
         type(tens_entry_mng_t):: this !inout: specialized tensor cache entry

         call this%nullify_tensor(.TRUE.) !deallocate the tensor component (if it is set)
         this%owner_id=-1
         return
        end subroutine tens_entry_mng_dtor
!-------------------------------------------------------------
        function tens_entry_mng_alloc(tens_entry) result(ierr)
!Non-member allocator for tens_entry_mng_t.
         implicit none
         integer(INTD):: ierr
         class(tens_cache_entry_t), allocatable, intent(out):: tens_entry

         allocate(tens_entry_mng_t::tens_entry,STAT=ierr)
         return
        end function tens_entry_mng_alloc
![tens_oprnd_t]=========================================
        subroutine TensOprndCtor(this,tensor,ierr,owner)
!Constructs a tensor operand. The <tensor> must be set.
!The owner id is optional.
         implicit none
         class(tens_oprnd_t), intent(inout):: this        !inout: undefined tensor operand (on entrance)
         class(tens_rcrsv_t), target, intent(in):: tensor !in: defined tensor
         integer(INTD), intent(out), optional:: ierr      !out: error code
         integer(INTD), intent(in), optional:: owner      !in: tensor owner id (no restrictions)
         integer(INTD):: errc

         if(.not.this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(tensor%is_set(errc)) then
            if(errc.eq.TEREC_SUCCESS) then
             this%tensor=>tensor
             if(present(owner)) this%owner_id=owner
             call this%mark_active(errc)
             if(errc.ne.DSVP_SUCCESS) errc=-5
            else
             errc=-4
            endif
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndCtor
!------------------------------------------------------------
        function TensOprndGetTensor(this,ierr) result(tens_p)
!Returns a pointer to the tensor.
         implicit none
         class(tens_rcrsv_t), pointer:: tens_p       !out: pointer to the tensor
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; tens_p=>this%tensor
         if(.not.associated(tens_p)) errc=-1
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetTensor
!---------------------------------------------------------
        function TensOprndGetOwnerId(this,ierr) result(id)
!Returns the tensor owner id.
         implicit none
         integer(INTD):: id                          !out: tensor owner id
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(associated(this%tensor)) then
          id=this%owner_id
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetOwnerId
!------------------------------------------------------
        subroutine TensOprndSetOwnerId(this,owner,ierr)
!Sets the tensor owner id (must be non-negative).
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: active tensor operand
         integer(INTD), intent(in):: owner           !in: tensor owner id (>=0)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(owner.ge.0) then
            this%owner_id=owner
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndSetOwnerId
!--------------------------------------------------------
        function TensOprndIsRemote(this,ierr) result(res)
!Returns TRUE if the tensor operand is remote, FALSE otherwise.
         implicit none
         logical:: res                               !out: result
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,id

         res=.TRUE. !assume remote by default
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           id=this%get_owner_id(errc)
           if(errc.eq.0) then
            res=id.ne.role_rank !tensor is owned by a different TAVP-MNG processor
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensOprndIsRemote
!------------------------------------------------
        subroutine TensOprndAcquireRsc(this,ierr)
!Acquires local resources for the remote tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         !`No local resources are currently needed
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndAcquireRsc
!----------------------------------------------
        subroutine TensOprndPrefetch(this,ierr)
!Starts prefetching the (remote) tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(.not.this%is_present(errc)) then
            if(errc.eq.DSVP_SUCCESS) then
             if(this%get_comm_stat().eq.DS_OPRND_NO_COMM) then
              call this%set_comm_stat(DS_OPRND_FETCHING,errc); if(errc.ne.DSVP_SUCCESS) errc=-6
             else
              errc=-5
             endif
            else
             errc=-4
            endif
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndPrefetch
!----------------------------------------------
        subroutine TensOprndUpload(this,ierr)
!Starts uploading the (remote) tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(this%is_present(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(this%get_comm_stat().eq.DS_OPRND_NO_COMM) then
            call this%set_comm_stat(DS_OPRND_UPLOADING,errc); if(errc.ne.DSVP_SUCCESS) errc=-4
           else
            errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndUpload
!---------------------------------------------------------
        function TensOprndSync(this,ierr,wait) result(res)
!Synchronizes a pending prefetch/upload.
         implicit none
         logical:: res                               !out: TRUE on communication completion, FALSE otherwise
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: TRUE activates WAIT instead of TEST synchronization (default)
         integer(INTD):: errc

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(this%get_comm_stat().ne.DS_OPRND_NO_COMM) then
            call this%mark_delivered(errc); if(errc.ne.DSVP_SUCCESS) errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensOprndSync
!---------------------------------------------
        subroutine TensOprndRelease(this,ierr)
!Releases local resources acquired for the remote tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         logical:: delivered

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(this%get_comm_stat().ne.DS_OPRND_NO_COMM) then
            delivered=this%sync(errc,wait=.TRUE.)
            if((.not.delivered).or.(errc.ne.0)) errc=-3
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         !`No local resources are currently needed
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndRelease
!----------------------------------------------
        subroutine TensOprndDestruct(this,ierr)
!Destructs the tensor operand.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,ier
         logical:: delivered

         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(this%get_comm_stat().ne.DS_OPRND_NO_COMM) then
            delivered=this%sync(errc,wait=.TRUE.)
            if((.not.delivered).or.(errc.ne.0)) errc=-4
           endif
           call this%mark_empty(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-3
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndDestruct
!---------------------------------------
        subroutine tens_oprnd_dtor(this)
         implicit none
         type(tens_oprnd_t):: this
         integer(INTD):: errc

         call this%destruct(errc)
         if(errc.ne.0) call quit(errc,'#FATAL(TAVP-MNG:tens_oprnd_dtor): Destructor failed!')
        end subroutine tens_oprnd_dtor

       end module tavp_manager
