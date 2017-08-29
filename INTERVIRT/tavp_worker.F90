!ExaTENSOR: TAVP-Worker implementation
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/08/29

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

       module tavp_worker
        use virta
        use gfc_base
        use gfc_list
        use gfc_vector
        use gfc_dictionary
        implicit none
        private
!PARAMETERS:
 !Basic:
        integer(INTD), private:: CONS_OUT=6 !default console output
        integer(INTD), private:: DEBUG=0    !debugging mode
        logical, private:: VERBOSE=.TRUE.   !verbosity for errors
 !Distributed memory space:
        integer(INTD), parameter, private:: TAVP_WORKER_NUM_WINS=1 !number of MPI windows in the DDSS distributed space
 !On-node pinned Host memory buffer:
        integer(INTL), protected:: tavp_worker_host_buf_size=1_INTL*(1024_INTL*1024_INTL*1024_INTL) !Host buffer size in bytes
        integer(INTD), protected:: tavp_worker_host_arg_max=0 !is set later: max number of tensors in the pinned Host buffer (initialized by TAL-SH)
 !Elementary tensor instruction granularity classification:
        real(8), protected:: TAVP_WORKER_FLOPS_HEAVY=1d12  !minimal number of Flops to consider the operation as heavy-cost
        real(8), protected:: TAVP_WORKER_FLOPS_MEDIUM=1d10 !minimal number of Flops to consider the operation as medium-cost
        real(8), protected:: TAVP_WORKER_COST_TO_SIZE=1d2  !minimal cost (Flops) to size (Words) ratio to consider the operation arithmetically intensive
!TYPES:
 !Tensor resource (local resource):
        type, extends(ds_resrc_t), private:: tens_resrc_t
         type(C_PTR), private:: base_addr=C_NULL_PTR   !C pointer to a local buffer for tensor body storage
         integer(C_SIZE_T), private:: bytes=0_C_SIZE_T !size of the tensor body storage buffer in bytes
         logical, private:: pinned=.FALSE.             !whether or not the buffer is pinned
         integer(C_INT), private:: dev_id=DEV_NULL     !flat device id where the buffer resides
         integer(C_INT), private:: ref_count=0         !reference count (how many tensor operands are associated with this resource)
         contains
          procedure, public:: is_empty=>TensResrcIsEmpty               !returns TRUE of the tensor resource is empty (unallocated)
          procedure, public:: allocate_buffer=>TensResrcAllocateBuffer !allocates a local buffer for tensor body storage
          procedure, public:: free_buffer=>TensResrcFreeBuffer         !frees the local buffer
          procedure, public:: get_mem_ptr=>TensResrcGetMemPtr          !returns a C pointer to the local memory buffer
          procedure, public:: get_mem_size=>TensResrcGetMemSize        !returns the size of the memory buffer in bytes
          procedure, private:: incr_ref_count=>TensResrcIncrRefCount   !increments the reference count
          procedure, private:: decr_ref_count=>TensResrcDecrRefCount   !decrements the reference count
          final:: tens_resrc_dtor                                      !dtor
        end type tens_resrc_t
#if 0
 !Tensor operand (encapsulated tensor data processible by a specific TAVP):
        type, extends(ds_oprnd_t), private:: tens_oprnd_t
         class(tens_rcrsv_t), pointer, private:: tensor=>NULL()   !non-owning pointer to a persistent recursive tensor
         class(tens_resrc_t), pointer, private:: resource=>NULL() !non-owning pointer to a persistent local tensor resource
         type(talsh_tens_t), private:: talsh_tens                 !TAL-SH tensor object (for performing actual computations)
         contains
          procedure, private:: TensOprndCtor                    !ctor
          generic, public:: tens_oprnd_ctor=>TensOprndCtor
          procedure, public:: set_resource=>TensOprndSetResource!sets the resource component if it has not been set via constructor
          procedure, public:: get_resource=>TensOprndGetResource!returns a pointer to the tensor resource
          procedure, public:: get_tensor=>TensOprndGetTensor    !returns a pointer to the tensor
          procedure, public:: set_talsh_tens=>TensOprndSetTalshTens !sets up the TAL-SH tensor object for further processing with TAL-SH
          procedure, public:: is_remote=>TensOprndIsRemote      !returns TRUE if the tensor operand is remote
          procedure, public:: acquire_rsc=>TensOprndAcquireRsc  !explicitly acquires local resources for the tensor operand
          procedure, public:: prefetch=>TensOprndPrefetch       !starts prefetching the remote tensor operand (acquires local resources!)
          procedure, public:: upload=>TensOprndUpload           !starts uploading the tensor operand to its remote location
          procedure, public:: sync=>TensOprndSync               !synchronizes the currently pending communication on the tensor operand
          procedure, public:: release=>TensOprndRelease         !destroys the present local copy of the tensor operand (releases local resources!), but the operand stays defined
          procedure, public:: destruct=>TensOprndDestruct       !performs complete destruction back to an empty state
          final:: tens_oprnd_dtor                               !dtor
        end type tens_oprnd_t
 !Tensor argument cache entry (TAVP-specific):
        type, extends(tens_cache_entry_t), private:: tens_entry_wrk_t
         type(tens_resrc_t), private:: resource                       !tensor resource
         contains
          procedure, private:: TensEntryWrkCtor
          procedure, public:: tens_entry_wrk_ctor=>TensEntryWrkCtor   !ctor
          procedure, public:: get_resource=>TensEntryWrkGetResource   !returns a pointer to the resource
          final:: tens_entry_wrk_dtor                                 !dtor
        end type tens_entry_wrk_t
 !Tensor instruction (realization of a tensor operation for a specific TAVP):
        type, extends(ds_instr_t), private:: tens_instr_t
         type(talsh_task_t), private:: talsh_task                   !TAL-SH task
         contains
          procedure, private:: TensInstrCtor                        !ctor: constructs a tensor instruction from the specification of a tensor operation
          generic, public:: tens_instr_ctor=>TensInstrCtor
          procedure, public:: encode=>TensInstrEncode               !encoding procedure: Packs the TAVP instruction into a raw byte packet (bytecode)
          final:: tens_instr_dtor                                   !dtor
        end type tens_instr_t
 !TAVP-WRK decoder:
        type, extends(ds_decoder_t), private:: tavp_wrk_decoder_t
         type(tens_cache_t), private:: tens_cache                   !tensor argument cache
         contains
          procedure, public:: start=>TAVPWRKDecoderStart            !starts TAVP-WRK decoder
          procedure, public:: shutdown=>TAVPWRKDecoderShutdown      !shuts down TAVP-WRK decoder
          procedure, public:: decode=>TAVPWRKDecoderDecode          !decodes the DS bytecode into a DS instruction
        end type tavp_wrk_decoder_t
 !TAVP-WRK encoder:
        type, extends(ds_encoder_t), private:: tavp_wrk_encoder_t
         contains
          procedure, public:: start=>TAVPWRKEncoderStart            !starts TAVP-WRK encoder
          procedure, public:: shutdown=>TAVPWRKEncoderShutdown      !shuts down TAVP-WRK encoder
          procedure, public:: encode=>TAVPWRKEncoderEncode          !encodes a DS instruction into the DS bytecode
        end type tavp_wrk_encoder_t
 !TAVP-WRK:
        type, extends(dsvp_t), public:: tavp_wrk_t
         type(tavp_wrk_decoder_t), private:: decoder      !DSVU: decodes incoming DS instructions from the manager
         type(tavp_wrk_encoder_t), private:: retirer      !DSVU: retires processed DS instructions and sends them back to the manager
         contains
          procedure, public:: configure=>TAVPWRKConfigure !configures the TAVP-WRK DSVP
        end type tavp_wrk_t
#endif
!VISIBILITY:
 !non-member test/debug:
        private test_carma
 !non-member control:
        private tavp_worker_set_host_buf_size
 !non-member TAVP microcode implementation:
        private acquire_resource_dummy
        private acquire_resource_basic
        private prefetch_input_dummy
        private prefetch_input_basic
        private sync_prefetch_dummy
        private sync_prefetch_basic
        private upload_output_dummy
        private upload_output_basic
        private sync_upload_dummy
        private sync_upload_basic
        private release_resource_dummy
        private release_resource_basic
        private sync_execution_dummy
        private sync_execution_basic
        private execute_dummy
        private execute_tensor_create
        private execute_tensor_destroy
        private execute_tensor_contract
#if 0
        private init_microcode
#endif
 !tens_resrc_t:
        private TensResrcIsEmpty
        private TensResrcAllocateBuffer
        private TensResrcFreeBuffer
        private TensResrcGetMemPtr
        private TensResrcGetMemSize
        private TensResrcIncrRefCount
        private TensResrcDecrRefCount
        public tens_resrc_dtor
#if 0
 !tens_oprnd_t:
        private TensOprndCtor
        private TensOprndSetResource
        private TensOprndGetResource
        private TensOprndGetTensor
        private TensOprndSetTalshTens
        private TensOprndIsRemote
        private TensOprndAcquireRsc
        private TensOprndPrefetch
        private TensOprndUpload
        private TensOprndSync
        private TensOprndRelease
        private TensOprndDestruct
        public tens_oprnd_dtor
 !tens_entry_wrk_t:
        private TensEntryWrkCtor
        private TensEntryWrkGetResource
        public tens_entry_wrk_dtor
 !tens_instr_t:
        private TensInstrCtor
        private TensInstrEncode
        public tens_instr_dtor
 !tavp_wrk_decoder_t:
        private TAVPWRKDecoderStart
        private TAVPWRKDecoderShutdown
        private TAVPWRKDecoderDecode
 !tavp_wrk_encoder_t:
        private TAVPWRKEncoderStart
        private TAVPWRKEncoderShutdown
        private TAVPWRKEncoderEncode
 !tavp_wrk_t:
        private TAVPWRKConfigure
#endif
!IMPLEMENTATION:
       contains
![non-member:Test/Debug]===========
        subroutine test_carma(ierr)
!DEBUG: Brute-force implementation of a single distributed tensor contraction.
         implicit none
         integer(INTD), intent(out):: ierr
         integer(INTD), parameter:: TENS_RANK=4
         integer(INTD), parameter:: DIM_SEG_SIZE=8
         integer(INTD), parameter:: DIM_NUM_LEVELS=4
         integer(INTD), parameter:: NUM_DIM_SEGS=2**(DIM_NUM_LEVELS-1)
         integer(INTD), parameter:: SPLIT_BASE=2**TENS_RANK
         integer(INTD), parameter:: TOTAL_BLOCKS=SPLIT_BASE**(DIM_NUM_LEVELS-1)
         integer(INTD), parameter:: BLOCK_VOL=DIM_SEG_SIZE**TENS_RANK
         integer(INTD):: num_procs,num_blocks,i,j,k,l,n,lid,rid
         integer(INTD):: dsg(1:TENS_RANK),lsg(1:TENS_RANK),rsg(1:TENS_RANK)
         integer(INTL):: tg
         type(DataDescr_t):: dd,ld,rd
         type(DataDescr_t), allocatable:: ddes(:),ldes(:),rdes(:)
         type(DataDescr_t), allocatable:: ddesa(:),ldesa(:),rdesa(:) !root only
         real(8), allocatable, target:: dtens(:),ltens(:),rtens(:)
         type(pack_env_t):: packenv
         type(obj_pack_t):: packet
         type(comm_handle_t):: chl
         real(8), pointer, contiguous:: block_p(:)
         type(C_PTR):: mem_p,dmem_p,lmem_p,rmem_p
         type(talsh_tens_t):: dtns,ltns,rtns
         type(talsh_task_t):: tsk0
         real(8):: tms,tm,tcs,tc,tts,tt,tas,ta

         ierr=0
         num_procs=role_size
         if(ierr.ne.0) call quit(-1,'Bad CARMA!')
         if(mod(TOTAL_BLOCKS,num_procs).ne.0) call quit(-2,'Bad CARMA!')
         num_blocks=TOTAL_BLOCKS/num_procs !number of tensor blocks per process
!Create <num_blocks> tensor blocks on each process:
 !Destination tensor:
         allocate(dtens(0:BLOCK_VOL*num_blocks-1)); allocate(ddes(0:num_blocks-1))
         dtens(:)=0d0
         do i=0,num_blocks-1
          block_p(0:)=>dtens(BLOCK_VOL*i:BLOCK_VOL*(i+1)-1); mem_p=c_loc(block_p)
          call tavp_addr_space%attach(mem_p,R8,int(BLOCK_VOL,8),ddes(i),ierr)
          if(ierr.ne.0) call quit(-3,'Bad CARMA!')
         enddo
 !Left tensor:
         allocate(ltens(0:BLOCK_VOL*num_blocks-1)); allocate(ldes(0:num_blocks-1))
         ltens(:)=1d-3
         do i=0,num_blocks-1
          block_p(0:)=>ltens(BLOCK_VOL*i:BLOCK_VOL*(i+1)-1); mem_p=c_loc(block_p)
          call tavp_addr_space%attach(mem_p,R8,int(BLOCK_VOL,8),ldes(i),ierr)
          if(ierr.ne.0) call quit(-4,'Bad CARMA!')
         enddo
 !Right tensor:
         allocate(rtens(0:BLOCK_VOL*num_blocks-1)); allocate(rdes(0:num_blocks-1))
         rtens(:)=1d-4
         do i=0,num_blocks-1
          block_p(0:)=>rtens(BLOCK_VOL*i:BLOCK_VOL*(i+1)-1); mem_p=c_loc(block_p)
          call tavp_addr_space%attach(mem_p,R8,int(BLOCK_VOL,8),rdes(i),ierr)
          if(ierr.ne.0) call quit(-5,'Bad CARMA!')
         enddo
!Root collects data descriptors from all other processes:
 !Destination tensor:
         if(role_rank.eq.0) then
          allocate(ddesa(0:TOTAL_BLOCKS-1))
          do i=1,num_procs-1 !remote descriptors
           call packenv%reserve_mem(ierr); if(ierr.ne.0) call quit(-6,'Bad CARMA!')
           do while(.not.packenv%receive(chl,ierr,comm=role_comm))
            if(ierr.ne.0) call quit(-7,'Bad CARMA!')
           enddo
           if(ierr.ne.0) call quit(-8,'Bad CARMA!')
           call chl%wait(ierr); if(ierr.ne.0) call quit(-9,'Bad CARMA!')
           n=packenv%get_num_packets(ierr); if(ierr.ne.0) call quit(-10,'Bad CARMA!')
           do j=1,n
            call packenv%extract_packet(j,packet,ierr,tg,preclean=.TRUE.); if(ierr.ne.0) call quit(-11,'Bad CARMA!')
            call ddesa(tg)%unpack(packet,ierr); if(ierr.ne.0) call quit(-12,'Bad CARMA!')
            call packet%clean(ierr); if(ierr.ne.0) call quit(-13,'Bad CARMA!')
           enddo
           call chl%clean(ierr); if(ierr.ne.0) call quit(-14,'Bad CARMA!')
           call packenv%destroy(ierr); if(ierr.ne.0) call quit(-15,'Bad CARMA!')
          enddo
          do i=0,num_blocks-1 !local descriptors
           ddesa(i)=ddes(i)
          enddo
         else
          call packenv%reserve_mem(ierr); if(ierr.ne.0) call quit(-16,'Bad CARMA!')
          do i=0,num_blocks-1
           call packenv%acquire_packet(packet,ierr); if(ierr.ne.0) call quit(-17,'Bad CARMA!')
           call ddes(i)%pack(packet,ierr); if(ierr.ne.0) call quit(-18,'Bad CARMA!')
           call packenv%seal_packet(ierr,tag=int(role_rank*num_blocks+i,INTL)); if(ierr.ne.0) call quit(-19,'Bad CARMA!')
          enddo
          call packenv%send(0,chl,ierr,comm=role_comm); if(ierr.ne.0) call quit(-20,'Bad CARMA!')
          call chl%wait(ierr); if(ierr.ne.0) call quit(-21,'Bad CARMA!')
          call chl%clean(ierr); if(ierr.ne.0) call quit(-22,'Bad CARMA!')
          call packenv%destroy(ierr); if(ierr.ne.0) call quit(-23,'Bad CARMA!')
         endif
         call role_barrier()
 !Left tensor:
         if(role_rank.eq.0) then
          allocate(ldesa(0:TOTAL_BLOCKS-1))
          do i=1,num_procs-1 !remote descriptors
           call packenv%reserve_mem(ierr); if(ierr.ne.0) call quit(-24,'Bad CARMA!')
           do while(.not.packenv%receive(chl,ierr,comm=role_comm))
            if(ierr.ne.0) call quit(-25,'Bad CARMA!')
           enddo
           if(ierr.ne.0) call quit(-26,'Bad CARMA!')
           call chl%wait(ierr); if(ierr.ne.0) call quit(-27,'Bad CARMA!')
           n=packenv%get_num_packets(ierr); if(ierr.ne.0) call quit(-28,'Bad CARMA!')
           do j=1,n
            call packenv%extract_packet(j,packet,ierr,tg,preclean=.TRUE.); if(ierr.ne.0) call quit(-29,'Bad CARMA!')
            call ldesa(tg)%unpack(packet,ierr); if(ierr.ne.0) call quit(-30,'Bad CARMA!')
            call packet%clean(ierr); if(ierr.ne.0) call quit(-31,'Bad CARMA!')
           enddo
           call chl%clean(ierr); if(ierr.ne.0) call quit(-32,'Bad CARMA!')
           call packenv%destroy(ierr); if(ierr.ne.0) call quit(-33,'Bad CARMA!')
          enddo
          do i=0,num_blocks-1 !local descriptors
           ldesa(i)=ldes(i)
          enddo
         else
          call packenv%reserve_mem(ierr); if(ierr.ne.0) call quit(-34,'Bad CARMA!')
          do i=0,num_blocks-1
           call packenv%acquire_packet(packet,ierr); if(ierr.ne.0) call quit(-35,'Bad CARMA!')
           call ldes(i)%pack(packet,ierr); if(ierr.ne.0) call quit(-36,'Bad CARMA!')
           call packenv%seal_packet(ierr,tag=int(role_rank*num_blocks+i,INTL)); if(ierr.ne.0) call quit(-37,'Bad CARMA!')
          enddo
          call packenv%send(0,chl,ierr,comm=role_comm); if(ierr.ne.0) call quit(-38,'Bad CARMA!')
          call chl%wait(ierr); if(ierr.ne.0) call quit(-39,'Bad CARMA!')
          call chl%clean(ierr); if(ierr.ne.0) call quit(-40,'Bad CARMA!')
          call packenv%destroy(ierr); if(ierr.ne.0) call quit(-41,'Bad CARMA!')
         endif
         call role_barrier()
 !Right tensor:
         if(role_rank.eq.0) then
          allocate(rdesa(0:TOTAL_BLOCKS-1))
          do i=1,num_procs-1 !remote descriptors
           call packenv%reserve_mem(ierr); if(ierr.ne.0) call quit(-42,'Bad CARMA!')
           do while(.not.packenv%receive(chl,ierr,comm=role_comm))
            if(ierr.ne.0) call quit(-43,'Bad CARMA!')
           enddo
           if(ierr.ne.0) call quit(-44,'Bad CARMA!')
           call chl%wait(ierr); if(ierr.ne.0) call quit(-45,'Bad CARMA!')
           n=packenv%get_num_packets(ierr); if(ierr.ne.0) call quit(-46,'Bad CARMA!')
           do j=1,n
            call packenv%extract_packet(j,packet,ierr,tg,preclean=.TRUE.); if(ierr.ne.0) call quit(-47,'Bad CARMA!')
            call rdesa(tg)%unpack(packet,ierr); if(ierr.ne.0) call quit(-48,'Bad CARMA!')
            call packet%clean(ierr); if(ierr.ne.0) call quit(-49,'Bad CARMA!')
           enddo
           call chl%clean(ierr); if(ierr.ne.0) call quit(-50,'Bad CARMA!')
           call packenv%destroy(ierr); if(ierr.ne.0) call quit(-51,'Bad CARMA!')
          enddo
          do i=0,num_blocks-1 !local descriptors
           rdesa(i)=rdes(i)
          enddo
         else
          call packenv%reserve_mem(ierr); if(ierr.ne.0) call quit(-52,'Bad CARMA!')
          do i=0,num_blocks-1
           call packenv%acquire_packet(packet,ierr); if(ierr.ne.0) call quit(-53,'Bad CARMA!')
           call rdes(i)%pack(packet,ierr); if(ierr.ne.0) call quit(-54,'Bad CARMA!')
           call packenv%seal_packet(ierr,tag=int(role_rank*num_blocks+i,INTL)); if(ierr.ne.0) call quit(-55,'Bad CARMA!')
          enddo
          call packenv%send(0,chl,ierr,comm=role_comm); if(ierr.ne.0) call quit(-56,'Bad CARMA!')
          call chl%wait(ierr); if(ierr.ne.0) call quit(-57,'Bad CARMA!')
          call chl%clean(ierr); if(ierr.ne.0) call quit(-58,'Bad CARMA!')
          call packenv%destroy(ierr); if(ierr.ne.0) call quit(-59,'Bad CARMA!')
         endif
         call role_barrier()
!Root creates and sends tasks to all processes:
         if(role_rank.eq.0) then
          do i=1,num_procs-1 !process rank
           call packenv%reserve_mem(ierr,mem_size=16*1048576_INTL,max_packets=65536); if(ierr.ne.0) call quit(-60,'Bad CARMA!')
           do j=i*num_blocks,(i+1)*num_blocks-1 !destination block flat id
            do k=0,NUM_DIM_SEGS-1
             do l=0,NUM_DIM_SEGS-1
              call flat2signa(j,dsg)
              lsg(1:TENS_RANK)=(/dsg(4),k,dsg(2),l/); lid=signa2flat(lsg)
              rsg(1:TENS_RANK)=(/dsg(3),l,dsg(1),k/); rid=signa2flat(rsg)
              call packenv%acquire_packet(packet,ierr); if(ierr.ne.0) call quit(-61,'Bad CARMA!')
              call ddesa(j)%pack(packet,ierr); if(ierr.ne.0) call quit(-62,'Bad CARMA!')
              call ldesa(lid)%pack(packet,ierr); if(ierr.ne.0) call quit(-63,'Bad CARMA!')
              call rdesa(rid)%pack(packet,ierr); if(ierr.ne.0) call quit(-64,'Bad CARMA!')
              call packenv%seal_packet(ierr); if(ierr.ne.0) then; print *,ierr; call quit(-65,'Bad CARMA!'); endif
             enddo
            enddo
           enddo
           call packenv%send(i,chl,ierr,comm=role_comm); if(ierr.ne.0) call quit(-66,'Bad CARMA!')
           call chl%wait(ierr); if(ierr.ne.0) call quit(-67,'Bad CARMA!')
           call chl%clean(ierr); if(ierr.ne.0) call quit(-68,'Bad CARMA!')
           call packenv%destroy(ierr); if(ierr.ne.0) call quit(-69,'Bad CARMA!')
          enddo
          call packenv%reserve_mem(ierr,mem_size=16*1048576_INTL,max_packets=65536); if(ierr.ne.0) call quit(-70,'Bad CARMA!')
          do j=0,num_blocks-1
           do k=0,NUM_DIM_SEGS-1
            do l=0,NUM_DIM_SEGS-1
             call flat2signa(j,dsg)
             lsg(1:TENS_RANK)=(/dsg(4),k,dsg(2),l/); lid=signa2flat(lsg)
             rsg(1:TENS_RANK)=(/dsg(3),l,dsg(1),k/); rid=signa2flat(rsg)
             call packenv%acquire_packet(packet,ierr); if(ierr.ne.0) call quit(-71,'Bad CARMA!')
             call ddesa(j)%pack(packet,ierr); if(ierr.ne.0) call quit(-72,'Bad CARMA!')
             call ldesa(lid)%pack(packet,ierr); if(ierr.ne.0) call quit(-73,'Bad CARMA!')
             call rdesa(rid)%pack(packet,ierr); if(ierr.ne.0) call quit(-74,'Bad CARMA!')
             call packenv%seal_packet(ierr); if(ierr.ne.0) then; print *,ierr; call quit(-75,'Bad CARMA!'); endif
            enddo
           enddo
          enddo
         else
          call packenv%reserve_mem(ierr,mem_size=16*1048576_INTL,max_packets=65536); if(ierr.ne.0) call quit(-76,'Bad CARMA!')
          do while(.not.packenv%receive(chl,ierr,comm=role_comm))
           if(ierr.ne.0) call quit(-77,'Bad CARMA!')
          enddo
          if(ierr.ne.0) call quit(-78,'Bad CARMA!')
          call chl%wait(ierr); if(ierr.ne.0) call quit(-79,'Bad CARMA!')
          call chl%clean(ierr); if(ierr.ne.0) call quit(-80,'Bad CARMA!')
         endif
         call role_barrier()
!All processes synchronize and execute their tasks:
         n=packenv%get_num_packets(ierr); if(ierr.ne.0) call quit(-81,'Bad CARMA!'); print *,role_rank,n
         ierr=mem_allocate(talsh_flat_dev_id(DEV_HOST,0),int(BLOCK_VOL*8,C_SIZE_T),YEP,lmem_p)
         if(ierr.ne.0) call quit(-82,'Bad CARMA!')
         ierr=mem_allocate(talsh_flat_dev_id(DEV_HOST,0),int(BLOCK_VOL*8,C_SIZE_T),YEP,rmem_p)
         if(ierr.ne.0) call quit(-83,'Bad CARMA!')
         call role_barrier()
         tc=0d0; tt=0d0; ta=0d0; tms=thread_wtime()
         do i=1,n
          j=mod((i-1)+role_rank,n)+1
          call packenv%extract_packet(j,packet,ierr,preclean=.TRUE.); if(ierr.ne.0) call quit(-84,'Bad CARMA!')
          call dd%unpack(packet,ierr); if(ierr.ne.0) call quit(-85,'Bad CARMA!')
          call ld%unpack(packet,ierr); if(ierr.ne.0) call quit(-86,'Bad CARMA!')
          call rd%unpack(packet,ierr); if(ierr.ne.0) call quit(-87,'Bad CARMA!')
          tts=thread_wtime()
          call ld%get_data(lmem_p,ierr); if(ierr.ne.0) call quit(-88,'Bad CARMA!')
          call rd%get_data(rmem_p,ierr); if(ierr.ne.0) call quit(-89,'Bad CARMA!')
          tt=tt+thread_wtime(tts)
          tas=thread_wtime()
          dmem_p=dd%get_data_ptr(ierr); if(ierr.ne.0) call quit(-90,'Bad CARMA!')
          ierr=talsh_tensor_construct(dtns,R8,(/DIM_SEG_SIZE,DIM_SEG_SIZE,DIM_SEG_SIZE,DIM_SEG_SIZE/),&
                                     &talsh_flat_dev_id(DEV_HOST,0),dmem_p)
          ierr=talsh_tensor_construct(ltns,R8,(/DIM_SEG_SIZE,DIM_SEG_SIZE,DIM_SEG_SIZE,DIM_SEG_SIZE/),&
                                     &talsh_flat_dev_id(DEV_HOST,0),lmem_p)
          if(ierr.ne.0) call quit(-91,'Bad CARMA!')
          ierr=talsh_tensor_construct(rtns,R8,(/DIM_SEG_SIZE,DIM_SEG_SIZE,DIM_SEG_SIZE,DIM_SEG_SIZE/),&
                                     &talsh_flat_dev_id(DEV_HOST,0),rmem_p)
          if(ierr.ne.0) call quit(-92,'Bad CARMA!')
          ta=ta+thread_wtime(tas)
!          tts=thread_wtime()
!          call ld%flush_data(ierr,local=.FALSE.); if(ierr.ne.0) call quit(-93,'Bad CARMA!')
!          call rd%flush_data(ierr,local=.FALSE.); if(ierr.ne.0) call quit(-94,'Bad CARMA!')
!          tt=tt+thread_wtime(tts)
          tcs=thread_wtime()
          ierr=talsh_tensor_contract('D(a,b,c,d)+=L(d,i,b,j)*R(c,j,a,i)',dtns,ltns,rtns,dev_id=0,dev_kind=DEV_HOST,&
                                    &copy_ctrl=COPY_TTT)
          if(ierr.ne.0) call quit(-95,'Bad CARMA!')
          tc=tc+thread_wtime(tcs)
          tas=thread_wtime()
          ierr=talsh_tensor_destruct(rtns); if(ierr.ne.0) call quit(-96,'Bad CARMA!')
          ierr=talsh_tensor_destruct(ltns); if(ierr.ne.0) call quit(-97,'Bad CARMA!')
          ierr=talsh_tensor_destruct(dtns); if(ierr.ne.0) call quit(-98,'Bad CARMA!')
          ta=ta+thread_wtime(tas)
         enddo
         tm=thread_wtime(tms); print *,'Rank ',role_rank,': Total ',tm,' s: Contract ',tc,' s: Comm ',tt,' s: Alloc ',ta,' s'
         call role_barrier()
         ierr=mem_free(talsh_flat_dev_id(DEV_HOST,0),rmem_p); if(ierr.ne.0) call quit(-99,'Bad CARMA!')
         ierr=mem_free(talsh_flat_dev_id(DEV_HOST,0),lmem_p); if(ierr.ne.0) call quit(-100,'Bad CARMA!')
         call packenv%destroy(ierr); if(ierr.ne.0) call quit(-101,'Bad CARMA!')
!Detach tensor blocks:
 !Right tensor:
         do i=0,num_blocks-1
          call tavp_addr_space%detach(rdes(i),ierr)
          if(ierr.ne.0) call quit(-102,'Bad CARMA!')
         enddo
         if(allocated(rdesa)) deallocate(rdesa)
         deallocate(rdes); deallocate(rtens)
 !Left tensor:
         do i=0,num_blocks-1
          call tavp_addr_space%detach(ldes(i),ierr)
          if(ierr.ne.0) call quit(-103,'Bad CARMA!')
         enddo
         if(allocated(ldesa)) deallocate(ldesa)
         deallocate(ldes); deallocate(ltens)
 !Destination tensor:
         do i=0,num_blocks-1
          call tavp_addr_space%detach(ddes(i),ierr)
          if(ierr.ne.0) call quit(-104,'Bad CARMA!')
         enddo
         if(allocated(ddesa)) deallocate(ddesa)
         deallocate(ddes); deallocate(dtens)
         return

        contains

         subroutine flat2tuple(flat,tuple)
          implicit none
          integer(INTD), intent(in):: flat
          integer(INTD), intent(inout):: tuple(1:*)
          integer(INTD):: i,n

          n=flat
          do i=DIM_NUM_LEVELS-1,1,-1
           tuple(i)=mod(n,SPLIT_BASE)
           n=n/SPLIT_BASE
          enddo
          return
         end subroutine flat2tuple

         function tuple2flat(tuple) result(flat)
          implicit none
          integer(INTD):: flat
          integer(INTD), intent(in):: tuple(1:*)
          integer(INTD):: i

          flat=0
          do i=1,DIM_NUM_LEVELS-1
           flat=flat*SPLIT_BASE+tuple(i)
          enddo
          return
         end function tuple2flat

         function tuple2segment(tuple,dimsn) result(segment)
          implicit none
          integer(INTD):: segment
          integer(INTD), intent(in):: tuple(1:*)
          integer(INTD), intent(in):: dimsn
          integer(INTD):: i,j,m

          segment=0
          do i=1,DIM_NUM_LEVELS-1
           m=tuple(i)
           do j=1,dimsn-1; m=m/2; enddo
           segment=segment*2+mod(m,2)
          enddo
          return
         end function tuple2segment

         function flat2segment(flat,dimsn) result(segment)
          implicit none
          integer(INTD):: segment
          integer(INTD), intent(in):: flat
          integer(INTD), intent(in):: dimsn
          integer(INTD):: tuple(1:DIM_NUM_LEVELS-1)

          call flat2tuple(flat,tuple)
          segment=tuple2segment(tuple,dimsn)
          return
         end function flat2segment

         subroutine tuple2signa(tuple,signa)
          implicit none
          integer(INTD), intent(in):: tuple(1:*)
          integer(INTD), intent(inout):: signa(1:*)
          integer(INTD):: i

          do i=1,TENS_RANK
           signa(i)=tuple2segment(tuple,i)
          enddo
          return
         end subroutine tuple2signa

         subroutine signa2tuple(signa,tuple)
          implicit none
          integer(INTD), intent(in):: signa(1:*)
          integer(INTD), intent(inout):: tuple(1:*)
          integer(INTD):: i,j,m

          tuple(1:DIM_NUM_LEVELS-1)=0
          do i=TENS_RANK,1,-1
           m=signa(i)
           do j=DIM_NUM_LEVELS-1,1,-1
            tuple(j)=tuple(j)*2+mod(m,2)
            m=m/2
           enddo
          enddo
          return
         end subroutine signa2tuple

         subroutine flat2signa(flat,signa)
          implicit none
          integer(INTD), intent(in):: flat
          integer(INTD), intent(inout):: signa(1:*)
          integer(INTD):: tuple(1:DIM_NUM_LEVELS-1)

          call flat2tuple(flat,tuple)
          call tuple2signa(tuple,signa)
          return
         end subroutine flat2signa

         function signa2flat(signa) result(flat)
          implicit none
          integer(INTD):: flat
          integer(INTD), intent(in):: signa(1:*)
          integer(INTD):: tuple(1:DIM_NUM_LEVELS-1)

          call signa2tuple(signa,tuple)
          flat=tuple2flat(tuple)
          return
         end function signa2flat

        end subroutine test_carma
![non-member:Control]=======================================
        subroutine tavp_worker_set_host_buf_size(bytes,ierr)
!Changes the default size of the pinned Host memory buffer.
!This procedure should be called before TAVP is initialized.
         implicit none
         integer(INTL), intent(in):: bytes           !in: size in bytes
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(bytes.gt.0_INTL) then
          tavp_worker_host_buf_size=bytes
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine tavp_worker_set_host_buf_size
![non-member:Microcode Implementation]==============
        subroutine acquire_resource_dummy(this,ierr)
!Dummy procedure for acquiring no resource.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine acquire_resource_dummy
!---------------------------------------------------
        subroutine acquire_resource_basic(this,ierr)
!Acquires resource for each tensor instruction operand.
!If some resources cannot be acquired now, returns TRY_LATER.
!In that case, the successfully acquired resources will be kept,
!unless an error other than TRY_LATER occurred.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code, possibly TRY_LATER
         integer(INTD):: errc,ier,n
         class(ds_oprnd_t), pointer:: oprnd

         n=this%get_num_operands(errc)
         if(errc.eq.DSVP_SUCCESS) then
          do while(n.gt.0)
           oprnd=>this%get_operand(n-1,ier)
           if(ier.eq.DSVP_SUCCESS) then
            call oprnd%acquire_rsc(ier)
            if(ier.ne.0) then
             if(ier.eq.TRY_LATER) then
              errc=ier
             else
              errc=-3; exit
             endif
            endif
           else
            errc=-2; exit
           endif
           n=n-1
          enddo
          if(errc.ne.0.and.errc.ne.TRY_LATER) call this%release_resource(n) !n is error code here
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine acquire_resource_basic
!-------------------------------------------------
        subroutine prefetch_input_dummy(this,ierr)
!Dummy procedure for prefetching no input.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine prefetch_input_dummy
!-------------------------------------------------
        subroutine prefetch_input_basic(this,ierr)
!Starts prefetching input tensor operands.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,ier,n
         class(ds_oprnd_t), pointer:: oprnd

         n=this%get_num_operands(errc)
         if(errc.eq.DSVP_SUCCESS) then
          do while(n.gt.1) !operand 0 is output
           n=n-1
           oprnd=>this%get_operand(n,ier)
           if(ier.eq.DSVP_SUCCESS) then
            call oprnd%prefetch(ier)
            if(ier.ne.0.and.errc.eq.0) errc=-3
           else
            if(errc.eq.0) errc=-2
           endif
          enddo
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine prefetch_input_basic
!---------------------------------------------------------------
        function sync_prefetch_dummy(this,ierr,wait) result(res)
!Dummy procedure for syncing no input prefetch.
         implicit none
         logical:: res                               !out: TRUE if synchronized
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: WAIT or TEST (defaults to WAIT)
         integer(INTD):: errc

         errc=0; res=.TRUE.
         if(present(ierr)) ierr=errc
         return
        end function sync_prefetch_dummy
!---------------------------------------------------------------
        function sync_prefetch_basic(this,ierr,wait) result(res)
!Synchronization on the input prefetch, either WAIT or TEST.
         implicit none
         logical:: res                               !out: TRUE if synchronized (all input operands)
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: WAIT or TEST (defaults to WAIT)
         integer(INTD):: errc,ier,n
         class(ds_oprnd_t), pointer:: oprnd
         logical:: wt

         errc=0; res=.FALSE.
         wt=.TRUE.; if(present(wait)) wt=wait
         n=this%get_num_operands(errc)
         if(errc.eq.DSVP_SUCCESS) then
          res=.TRUE.
          do while(n.gt.1)
           n=n-1
           oprnd=>this%get_operand(n,ier)
           if(ier.eq.DSVP_SUCCESS) then
            res=res.and.oprnd%sync(ier,wt)
            if(ier.ne.0.and.errc.eq.0) then; res=.FALSE.; errc=-3; endif
           else
            res=.FALSE.; if(errc.eq.0) errc=-2
           endif
          enddo
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function sync_prefetch_basic
!------------------------------------------------
        subroutine upload_output_dummy(this,ierr)
!Dummy procedure for uploading no output.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine upload_output_dummy
!------------------------------------------------
        subroutine upload_output_basic(this,ierr)
!Starts uploading the output tensor operand 0.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         class(ds_oprnd_t), pointer:: oprnd

         oprnd=>this%get_operand(0,errc) !output operand
         if(errc.eq.DSVP_SUCCESS) then
          call oprnd%upload(errc); if(errc.ne.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine upload_output_basic
!-------------------------------------------------------------
        function sync_upload_dummy(this,ierr,wait) result(res)
!Dummy procedure for syncing no output upload.
         implicit none
         logical:: res                               !out: TRUE if synchronized
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: WAIT or TEST (defaults to WAIT)
         integer(INTD):: errc

         errc=0; res=.TRUE.
         if(present(ierr)) ierr=errc
         return
        end function sync_upload_dummy
!-------------------------------------------------------------
        function sync_upload_basic(this,ierr,wait) result(res)
!Synchronization on the output upload.
         implicit none
         logical:: res                               !out: TRUE if synchronized
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: WAIT or TEST (defaults to WAIT)
         integer(INTD):: errc
         class(ds_oprnd_t), pointer:: oprnd
         logical:: wt

         errc=0; res=.FALSE.
         wt=.TRUE.; if(present(wait)) wt=wait
         oprnd=>this%get_operand(0,errc)
         if(errc.eq.DSVP_SUCCESS) then
          res=oprnd%sync(errc,wt); if(errc.ne.0) then; res=.FALSE.; errc=-2; endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function sync_upload_basic
!---------------------------------------------------
        subroutine release_resource_dummy(this,ierr)
!Dummy procedure for releasing no resource.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine release_resource_dummy
!---------------------------------------------------
        subroutine release_resource_basic(this,ierr)
!Releases resources occupied by tensor instruction operands,
!but the operands stay defined.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,ier,n
         class(ds_oprnd_t), pointer:: oprnd

         n=this%get_num_operands(errc)
         if(errc.eq.DSVP_SUCCESS) then
          do while(n.gt.0)
           n=n-1
           oprnd=>this%get_operand(n,ier)
           if(ier.eq.DSVP_SUCCESS) then
            call oprnd%release(ier)
            if(ier.ne.0.and.errc.eq.0) errc=-4
           else
            if(errc.eq.0) errc=-3
           endif
          enddo
          select type(this)
          class is(tens_instr_t)
           ier=talsh_task_destruct(this%talsh_task)
           if(ier.ne.TALSH_SUCCESS.and.errc.eq.0) errc=-2
          end select
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine release_resource_basic
!----------------------------------------------------------------
        function sync_execution_dummy(this,ierr,wait) result(res)
!Dummy procedure for syncing no execution.
         implicit none
         logical:: res                               !out: TRUE if synchronized
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: WAIT or TEST (defaults to WAIT)
         integer(INTD):: errc

         errc=0; res=.TRUE.
         if(present(ierr)) ierr=errc
         return
        end function sync_execution_dummy
!----------------------------------------------------------------
        function sync_execution_basic(this,ierr,wait) result(res)
!Synchronization on the tensor instruction execution.
         implicit none
         logical:: res                               !out: TRUE if synchronized
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: WAIT or TEST (defaults to WAIT)
         integer(INTD):: errc,sts,ans
         logical:: wt

         errc=0; res=.FALSE.
         wt=.TRUE.; if(present(wait)) wt=wait
         select type(this)
         class is(tens_instr_t)
          if(wt) then !WAIT
           errc=talsh_task_wait(this%talsh_task,sts)
           res=(errc.eq.TALSH_SUCCESS.and.(sts.eq.TALSH_TASK_COMPLETED.or.sts.eq.TALSH_TASK_ERROR))
          else !TEST
           ans=talsh_task_complete(this%talsh_task,sts,errc)
           res=(ans.eq.YEP.and.errc.eq.0)
          endif
          if(sts.eq.TALSH_TASK_ERROR.and.errc.eq.0) errc=-1
         end select
         if(present(ierr)) ierr=errc
         return
        end function sync_execution_basic
!------------------------------------------
        subroutine execute_dummy(this,ierr)
!Dummy procedure for executing nothing.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine execute_dummy
!--------------------------------------------------
        subroutine execute_tensor_create(this,ierr)
!Executes tensor creation. The tensor layout is assumed already defined.
!The tensor body location is set here via the newly created DDSS descriptor.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,dtk
         integer(INTL):: bytes,vol
         class(ds_oprnd_t), pointer:: oprnd
         class(tens_rcrsv_t), pointer:: tensor
         class(tens_resrc_t), pointer:: resource
         class(tens_body_t), pointer:: tens_body
         class(tens_layout_t), pointer:: tens_layout
         type(DataDescr_t):: descr
         type(C_PTR):: mem_p

         oprnd=>this%get_operand(0,errc)
         if(errc.eq.DSVP_SUCCESS) then
          select type(oprnd)
          class is(tens_oprnd_t)
           resource=>oprnd%get_resource(errc)
           if(errc.eq.0) then
            if(.not.resource%is_empty()) then
             tensor=>oprnd%get_tensor(errc)
             if(errc.eq.0) then
              tens_body=>tensor%get_body(errc)
              if(errc.eq.TEREC_SUCCESS) then
               tens_layout=>tens_body%get_layout(errc)
               if(errc.eq.TEREC_SUCCESS) then
                dtk=tens_layout%get_data_type(errc)
                if(errc.eq.TEREC_SUCCESS) then
                 mem_p=resource%get_mem_ptr(errc,bytes)
                 if(errc.eq.0) then
                  vol=tens_layout%get_volume() !bytes = vol * sizeof(data_kind)
                  call tavp_addr_space%attach(mem_p,dtk,vol,descr,errc)
                  if(errc.eq.0) then
                   call tensor%set_location(descr,errc) !tensor has been located
                   if(errc.ne.TEREC_SUCCESS) errc=-11
                  else
                   errc=-10
                  endif
                 else
                  errc=-9
                 endif
                else
                 errc=-8
                endif
               else
                errc=-7
               endif
              else
               errc=-6
              endif
             else
              errc=-5
             endif
            else
             errc=-4
            endif
           else
            errc=-3
           endif
          class default
           errc=-2
          end select
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine execute_tensor_create
!---------------------------------------------------
        subroutine execute_tensor_destroy(this,ierr)
!Executes tensor destruction.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine execute_tensor_destroy
!----------------------------------------------------
        subroutine execute_tensor_contract(this,ierr)
!Executes a tensor contraction operation.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine execute_tensor_contract
#if 0
!--------------------------------------
        subroutine init_microcode(ierr)
!Initializes the global table of TAVP microcode bindings when configuring TAVP.
!There are three kinds of microcode:
! 1. DUMMY: Does nothing.
! 2. BASIC: Default microcode.
! 3. SPECIAL: Specialized microcode.
         implicit none
         integer(INTD), intent(out), optional:: ierr
         integer(INTD):: errc

         errc=0
 !TENSOR CREATE:
         microcode(TAVP_INSTR_CREATE)%acquire_resource=>acquire_resource_basic
         microcode(TAVP_INSTR_CREATE)%prefetch_input=>prefetch_input_dummy
         microcode(TAVP_INSTR_CREATE)%sync_prefetch=>sync_prefetch_dummy
         microcode(TAVP_INSTR_CREATE)%execute=>execute_tensor_create
         microcode(TAVP_INSTR_CREATE)%sync_execution=>sync_execution_dummy
         microcode(TAVP_INSTR_CREATE)%upload_output=>upload_output_dummy
         microcode(TAVP_INSTR_CREATE)%sync_upload=>sync_upload_dummy
         microcode(TAVP_INSTR_CREATE)%release_resource=>release_resource_dummy
 !TENSOR DESTROY:
         microcode(TAVP_INSTR_DESTROY)%acquire_resource=>acquire_resource_dummy
         microcode(TAVP_INSTR_DESTROY)%prefetch_input=>prefetch_input_dummy
         microcode(TAVP_INSTR_DESTROY)%sync_prefetch=>sync_prefetch_dummy
         microcode(TAVP_INSTR_DESTROY)%execute=>execute_tensor_destroy
         microcode(TAVP_INSTR_DESTROY)%sync_execution=>sync_execution_dummy
         microcode(TAVP_INSTR_DESTROY)%upload_output=>upload_output_dummy
         microcode(TAVP_INSTR_DESTROY)%sync_upload=>sync_upload_dummy
         microcode(TAVP_INSTR_DESTROY)%release_resource=>release_resource_basic
 !TENSOR CONTRACT:
         microcode(TAVP_INSTR_CONTRACT)%acquire_resource=>acquire_resource_basic
         microcode(TAVP_INSTR_CONTRACT)%prefetch_input=>prefetch_input_basic
         microcode(TAVP_INSTR_CONTRACT)%sync_prefetch=>sync_prefetch_basic
         microcode(TAVP_INSTR_CONTRACT)%execute=>execute_tensor_contract
         microcode(TAVP_INSTR_CONTRACT)%sync_execution=>sync_execution_basic
         microcode(TAVP_INSTR_CONTRACT)%upload_output=>upload_output_basic
         microcode(TAVP_INSTR_CONTRACT)%sync_upload=>sync_upload_basic
         microcode(TAVP_INSTR_CONTRACT)%release_resource=>release_resource_basic
 !DONE
         if(present(ierr)) ierr=errc
         return
        end subroutine init_microcode
#endif
!tens_resrc_t]=====================================
        function TensResrcIsEmpty(this) result(ans)
!Returns TRUE if the tensor resource is empty.
         implicit none
         logical:: ans                          !out: answer
         class(tens_resrc_t), intent(in):: this !in: tensor resource

         ans=(this%bytes.le.0_C_SIZE_T)
         return
        end function TensResrcIsEmpty
!---------------------------------------------------------------------------
        subroutine TensResrcAllocateBuffer(this,bytes,ierr,in_buffer,dev_id)
!Allocates local memory either from a system or from a custom buffer.
         implicit none
         class(tens_resrc_t), intent(inout):: this    !inout: tensor resource
         integer(INTL), intent(in):: bytes            !in: size in bytes
         integer(INTD), intent(out), optional:: ierr  !out: error code
         logical, intent(in), optional:: in_buffer    !in: if TRUE the memory will be allocated from a custom buffer, FALSE from the system
         integer(INTD), intent(in), optional:: dev_id !in: flat device id (defaults to Host)
         integer(INTD):: errc
         integer(C_INT):: in_buf,dev
         type(C_PTR):: addr

         errc=0
         if(this%is_empty()) then
          if(bytes.gt.0_INTL) then
           in_buf=NOPE; if(present(in_buffer)) then; if(in_buffer) in_buf=YEP; endif
           dev=talsh_flat_dev_id(DEV_HOST,0); if(present(dev_id)) dev=dev_id
           errc=mem_allocate(dev,int(bytes,C_SIZE_T),in_buf,addr)
           if(errc.eq.0) then
            this%base_addr=addr
            this%bytes=bytes
            this%pinned=(in_buf.ne.NOPE)
            this%dev_id=dev
           else
            errc=-1
           endif
          else
           errc=-2
          endif
         else
          errc=-3
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensResrcAllocateBuffer
!------------------------------------------------
        subroutine TensResrcFreeBuffer(this,ierr)
!Frees the tensor resource buffer.
         implicit none
         class(tens_resrc_t), intent(inout):: this    !inout: tensor resource
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTD):: errc

         errc=0
         if(.not.this%is_empty()) then !free only allocated resources
          if(this%ref_count.eq.0) then
           errc=mem_free(this%dev_id,this%base_addr)
           if(errc.eq.0) then
            this%base_addr=C_NULL_PTR
            this%bytes=0_C_SIZE_T
            this%pinned=.FALSE.
            this%dev_id=DEV_NULL
           else
            errc=-2
           endif
          else
           errc=-1
          endif
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensResrcFreeBuffer
!-----------------------------------------------------------------
        function TensResrcGetMemPtr(this,ierr,bytes) result(mem_p)
!Returns a C pointer to the local memory buffer used by the resource.
         implicit none
         type(C_PTR):: mem_p                          !out: C pointer
         class(tens_resrc_t), intent(in):: this       !in: tensor resource
         integer(INTD), intent(out), optional:: ierr  !out: error code
         integer(INTL), intent(out), optional:: bytes !out: number of bytes
         integer(INTD):: errc

         errc=0; mem_p=C_NULL_PTR
         if(.not.this%is_empty()) then
          mem_p=this%base_addr
          if(present(bytes)) bytes=this%bytes
         else
          if(present(bytes)) bytes=0_INTL
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensResrcGetMemPtr
!------------------------------------------------------------
        function TensResrcGetMemSize(this,ierr) result(bytes)
!Returns the size of the memory buffer in bytes.
         implicit none
         integer(INTL):: bytes                       !out: size in bytes
         class(tens_resrc_t), intent(in):: this      !in: tensor resource
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         if(.not.this%is_empty()) then
          bytes=this%bytes; errc=0
         else
          bytes=0_INTL; errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensResrcGetMemSize
!---------------------------------------------
        subroutine TensResrcIncrRefCount(this)
!Increments the reference count.
         implicit none
         class(tens_resrc_t), intent(inout):: this !inout: tensor resource

         this%ref_count=this%ref_count+1
         return
        end subroutine TensResrcIncrRefCount
!---------------------------------------------
        subroutine TensResrcDecrRefCount(this)
!Decrements the reference count.
         implicit none
         class(tens_resrc_t), intent(inout):: this !inout: tensor resource

         this%ref_count=this%ref_count-1
         return
        end subroutine TensResrcDecrRefCount
!---------------------------------------
        subroutine tens_resrc_dtor(this)
         implicit none
         type(tens_resrc_t):: this
         integer(INTD):: errc

         call this%free_buffer(errc)
         if(errc.ne.0) call quit(errc,'#FATAL(TAVP-WRK:tens_resrc_dtor): Attempt to free a memory buffer which is still in use!')
         return
        end subroutine tens_resrc_dtor
#if 0
![tens_oprnd_t]=================================================
        subroutine TensOprndCtor(this,tensor,ierr,tens_resource)
!Constructs a tensor operand. The <tensor> must be set.
!The associated tensor resource is optional and may still be empty.
         implicit none
         class(tens_oprnd_t), intent(inout):: this                            !inout: tensor operand
         class(tens_rcrsv_t), target, intent(in):: tensor                     !in: tensor
         integer(INTD), intent(out), optional:: ierr                          !out: error code
         class(tens_resrc_t), target, intent(inout), optional:: tens_resource !in: local tensor resource, may still be empty
         integer(INTD):: errc

         errc=0
         if(.not.this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(tensor%is_set(errc)) then
            if(errc.eq.TEREC_SUCCESS) then
             this%tensor=>tensor
             if(present(tens_resource)) then
              call tens_resource%incr_ref_count()
              this%resource=>tens_resource
             endif
             call this%mark_active(errc)
             if(errc.ne.DSVP_SUCCESS) errc=-1
            else
             errc=-2
            endif
           else
            errc=-3
           endif
          else
           errc=-4
          endif
         else
          errc=-5
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndCtor
!---------------------------------------------------------------
        subroutine TensOprndSetResource(this,tens_resource,ierr)
!Sets the resource component if it has not been set via constructor.
         implicit none
         class(tens_oprnd_t), intent(inout):: this                  !inout: active tensor operand
         class(tens_resrc_t), target, intent(inout):: tens_resource !inout: tensor resource (may be empty)
         integer(INTD), intent(out), optional:: ierr                !out: error code
         integer(INTD):: errc

         if(this%is_active(errc)) then
          if(errc.eq.0) then
           if(.not.associated(this%resource)) then
            call tens_resource%incr_ref_count()
            this%resource=>tens_resource
           else
            errc=-1
           endif
          else
           errc=-2
          endif
         else
          errc=-3
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndSetResource
!------------------------------------------------------------------
        function TensOprndGetResource(this,ierr) result(resource_p)
!Returns a pointer to the tensor resource.
         implicit none
         type(tens_resrc_t), pointer:: resource_p    !out: pointer to the tensor resource
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; resource_p=>this%resource
         if(.not.associated(resource_p)) errc=-1
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetResource
!------------------------------------------------------------
        function TensOprndGetTensor(this,ierr) result(tens_p)
!Returns a pointer to the tensor.
         implicit none
         type(tens_rcrsv_t), pointer:: tens_p        !out: pointer to the tensor
         class(tens_oprnd_t), intent(in):: this      !in: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0; tens_p=>this%tensor
         if(.not.associated(tens_p)) errc=-1
         if(present(ierr)) ierr=errc
         return
        end function TensOprndGetTensor
!--------------------------------------------------
        subroutine TensOprndSetTalshTens(this,ierr)
!Sets up the TAL-SH tensor object for further processing with TAL-SH.
!The local tensor body location is imported from the tensor resource.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: active tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,nd,unres,lay,data_kind,dims(1:MAX_TENSOR_RANK)
         integer(INTL):: dim_ext(1:MAX_TENSOR_RANK)
         logical:: shpd,layd,locd
         type(C_PTR):: mem_p
         class(tens_rcrsv_t), pointer:: tens_p
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p
         class(tens_header_t), pointer:: header_p

         if(this%is_active(errc)) then
          if(errc.eq.0) then
           if(talsh_tensor_is_empty(this%talsh_tens)) then
            tens_p=>this%tensor
            if(tens_p%is_set(errc,num_dims=nd,shaped=shpd,unresolved=unres,layed=layd,located=locd)) then
             if((errc.eq.0).and.(unres.eq.0).and.shpd.and.layd) then
              body_p=>tens_p%get_body(errc)
              if(errc.eq.TEREC_SUCCESS) then
               layout_p=>body_p%get_layout(errc)
               if(errc.eq.TEREC_SUCCESS) then
                lay=layout_p%get_layout_kind(errc)
                if(errc.eq.TEREC_SUCCESS.and.lay.eq.TEREC_LAY_FDIMS) then
                 data_kind=layout_p%get_data_type(errc)
                 if(errc.eq.TEREC_SUCCESS) then
                  if(.not.this%resource%is_empty()) then
                   mem_p=this%resource%get_mem_ptr(errc)
                   if(errc.eq.0) then
                    header_p=>tens_p%get_header(errc)
                    if(errc.eq.TEREC_SUCCESS) then
                     call header_p%get_dims(dim_ext,nd,errc)
                     if(errc.eq.TEREC_SUCCESS) then
                      dims(1:nd)=dim_ext(1:nd)
                      errc=talsh_tensor_construct(this%talsh_tens,data_kind,dims(1:nd),ext_mem=mem_p)
                      if(errc.ne.TALSH_SUCCESS) errc=-14
                     else
                      errc=-13
                     endif
                    else
                     errc=-12
                    endif
                    header_p=>NULL()
                   else
                    errc=-11
                   endif
                  else
                   errc=-10
                  endif
                 else
                  errc=-9
                 endif
                else
                 errc=-8
                endif
               else
                errc=-7
               endif
               layout_p=>NULL()
              else
               errc=-6
              endif
              body_p=>NULL()
             else
              errc=-5
             endif
            else
             errc=-4
            endif
            tens_p=>NULL()
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
        end subroutine TensOprndSetTalshTens
!--------------------------------------------------------
        function TensOprndIsRemote(this,ierr) result(res)
!Returns TRUE if the tensor operand is remote, FALSE otherwise.
         implicit none
         logical:: res                               !out: result
         class(tens_oprnd_t), intent(in):: this      !in: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer(INT_MPI):: host_proc_rank,mpi_comm,my_rank
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p
         class(DataDescr_t), pointer:: descr_p

         errc=0
         if(this%is_active()) then
          body_p=>this%tensor%get_body(errc)
          if(errc.eq.TEREC_SUCCESS.and.associated(body_p)) then
           layout_p=>body_p%get_layout(errc)
           if(errc.eq.TEREC_SUCCESS.and.associated(layout_p)) then
            descr_p=>layout_p%get_data_descr(errc)
            if(errc.eq.TEREC_SUCCESS.and.associated(descr_p)) then
             if(descr_p%is_set(errc,host_proc_rank,mpi_comm)) then
              if(errc.eq.0) then
               call MPI_Comm_Rank(mpi_comm,my_rank,errc)
               if(errc.eq.0) then
                res=.not.(host_proc_rank.eq.my_rank)
               else
                errc=-1
               endif
              else
               errc=-2
              endif
             else
              errc=-3
             endif
            else
             errc=-4
            endif
           else
            errc=-5
           endif
          else
           errc=-6
          endif
         else
          errc=-7
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensOprndIsRemote
!------------------------------------------------
        subroutine TensOprndAcquireRsc(this,ierr)
!Acquires local resources for the remote tensor operand.
!If the resources have already been allocated, does nothing.
!If the resource component is not set, an error will be returned.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand (with an associated resource component)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         integer(INTL):: buf_size
         integer(INT_MPI):: host_proc_rank
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p

         errc=0
         if(this%is_active()) then
          if(associated(this%resource)) then
           body_p=>this%tensor%get_body(errc)
           if(errc.eq.TEREC_SUCCESS.and.associated(body_p)) then
            layout_p=>body_p%get_layout(errc)
            if(errc.eq.TEREC_SUCCESS.and.associated(layout_p)) then
             if(layout_p%is_set(errc)) then
              if(errc.eq.TEREC_SUCCESS) then
               buf_size=layout_p%get_body_size(errc)
               if(errc.eq.TEREC_SUCCESS.and.buf_size.gt.0_INTL) then
                if(this%resource%is_empty()) then
                 call this%resource%allocate_buffer(buf_size,errc); if(errc.ne.0) errc=-1
                endif
               else
                errc=-2
               endif
              else
               errc=-3
              endif
             else
              errc=-4
             endif
            else
             errc=-5
            endif
           else
            errc=-6
           endif
          else
           errc=-7
          endif
         else
          errc=-8
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndAcquireRsc
!----------------------------------------------
        subroutine TensOprndPrefetch(this,ierr)
!Starts prefetching the (remote) tensor operand using the local tensor resource.
!If the resource component has not been set, an error will be returned.
!If the local resource has not been allocated, it will be allocated here.
!If the tensor operand has been delivered before, does nothing.
!If there is a pending communication on the tensor operand, returns an error.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand (with an associated resource component)
         integer(INTD), intent(out), optional:: ierr !out: error code, may return TRY_LATER
         integer(INTD):: errc
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p
         class(DataDescr_t), pointer:: descr_p
         type(C_PTR):: cptr

         if(.not.this%is_present(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(associated(this%resource)) then
            body_p=>this%tensor%get_body(errc)
            if(errc.eq.TEREC_SUCCESS.and.associated(body_p)) then
             layout_p=>body_p%get_layout(errc)
             if(errc.eq.TEREC_SUCCESS.and.associated(layout_p)) then
              descr_p=>layout_p%get_data_descr(errc)
              if(errc.eq.TEREC_SUCCESS.and.associated(descr_p)) then
               if(descr_p%is_set(errc)) then
                if(errc.eq.0) then
                 if(this%resource%is_empty()) call this%acquire_rsc(errc)
                 if(errc.eq.0) then
                  if(this%get_comm_stat().eq.DS_OPRND_NO_COMM) then
                   cptr=this%resource%base_addr
                   call descr_p%get_data(cptr,errc,MPI_ASYNC_REQ)
                   if(errc.ne.0.and.errc.ne.TRY_LATER) errc=-1
                   if(errc.eq.0) then
                    call this%set_comm_stat(DS_OPRND_FETCHING,errc); if(errc.ne.DSVP_SUCCESS) errc=-2
                   endif
                  else
                   errc=-3
                  endif
                 else
                  errc=-4
                 endif
                else
                 errc=-5
                endif
               else
                errc=-6
               endif
              else
               errc=-7
              endif
             else
              errc=-8
             endif
            else
             errc=-9
            endif
           else
            errc=-10
           endif
          else
           errc=-11
          endif
         else
          errc=-12
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndPrefetch
!--------------------------------------------
        subroutine TensOprndUpload(this,ierr)
!Starts uploading the (remote) tensor operand from the local tensor resource.
!The tensor operand must be marked as delivered (present), even if it is local.
!If there is a pending communication on the tensor operand, returns an error.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p
         class(DataDescr_t), pointer:: descr_p
         type(C_PTR):: cptr

         if(this%is_present(errc)) then !assumes that the local resource is allocated
          if(errc.eq.DSVP_SUCCESS) then
           body_p=>this%tensor%get_body(errc)
           if(errc.eq.TEREC_SUCCESS.and.associated(body_p)) then
            layout_p=>body_p%get_layout(errc)
            if(errc.eq.TEREC_SUCCESS.and.associated(layout_p)) then
             descr_p=>layout_p%get_data_descr(errc)
             if(errc.eq.TEREC_SUCCESS.and.associated(descr_p)) then
              if(descr_p%is_set(errc)) then
               if(errc.eq.0) then
                if(.not.this%resource%is_empty()) then !trap
                 if(this%get_comm_stat().eq.DS_OPRND_NO_COMM) then
                  cptr=this%resource%base_addr
                  call descr_p%acc_data(cptr,errc,MPI_ASYNC_REQ)
                  if(errc.ne.0.and.errc.ne.TRY_LATER) errc=-1
                  if(errc.eq.0) then
                   call this%set_comm_stat(DS_OPRND_UPLOADING,errc); if(errc.ne.DSVP_SUCCESS) errc=-2
                  endif
                 else
                  errc=-3
                 endif
                else
                 errc=-4
                endif
               else
                errc=-5
               endif
              else
               errc=-6
              endif
             else
              errc=-7
             endif
            else
             errc=-8
            endif
           else
            errc=-9
           endif
          else
           errc=-10
          endif
         else
          errc=-11
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensOprndUpload
!---------------------------------------------------------
        function TensOprndSync(this,ierr,wait) result(res)
!Synchronizes the pending prefetch/upload, either TEST or WAIT.
!A successful synchronization on prefetch will mark the tensor operand
!as delivered (present). A successful synchronization on upload will
!not change the status of the tensor operand (which is present).
         implicit none
         logical:: res                               !out: TRUE on communication completion, FALSE otherwise
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         logical, intent(in), optional:: wait        !in: TRUE activates WAIT instead of TEST synchronization (default)
         integer(INTD):: errc,sts
         class(tens_body_t), pointer:: body_p
         class(tens_layout_t), pointer:: layout_p
         class(DataDescr_t), pointer:: descr_p
         logical:: tw

         errc=0; res=.FALSE.
         tw=.FALSE.; if(present(wait)) tw=wait
         if(this%is_active()) then
          body_p=>this%tensor%get_body(errc)
          if(errc.eq.TEREC_SUCCESS.and.associated(body_p)) then
           layout_p=>body_p%get_layout(errc)
           if(errc.eq.TEREC_SUCCESS.and.associated(layout_p)) then
            descr_p=>layout_p%get_data_descr(errc)
            if(errc.eq.TEREC_SUCCESS.and.associated(descr_p)) then
             if(descr_p%is_set(errc)) then
              if(errc.eq.0) then
               sts=this%get_comm_stat()
               if(sts.ne.DS_OPRND_NO_COMM) then
                if(tw) then
                 call descr_p%wait_data(errc); if(errc.eq.0) then; res=.TRUE.; else; errc=-1; endif
                else
                 res=descr_p%test_data(errc); if(errc.ne.0) errc=-2
                endif
                if(sts.eq.DS_OPRND_FETCHING.and.res) then
                 call this%mark_delivered(errc); if(errc.ne.0) errc=-3
                endif
               else
                errc=-4
               endif
              else
               errc=-5
              endif
             else
              errc=-6
             endif
            else
             errc=-7
            endif
           else
            errc=-8
           endif
          else
           errc=-9
          endif
         else
          errc=-10
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensOprndSync
!---------------------------------------------
        subroutine TensOprndRelease(this,ierr)
!Releases local tensor resources occupied by the tensor operand,
!unless there are other active tensor operands sharing the same resource.
!In the latter case, nothing will be done and no error raised.
!Also, if the resource component is not set, nothing will be done either.
         implicit none
         class(tens_oprnd_t), intent(inout):: this   !inout: tensor operand
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc
         logical:: delivered

         errc=0
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(associated(this%resource)) then
            if(this%get_comm_stat().ne.DS_OPRND_NO_COMM) then
             delivered=this%sync(errc,wait=.TRUE.)
             if(.not.delivered.or.errc.ne.DSVP_SUCCESS) errc=-1
            endif
            if(this%resource%ref_count.eq.1) then !only one (last) tensor operand is associated with this resource
             call this%resource%free_buffer(errc); if(errc.ne.0) errc=-2
            endif
           endif
          else
           errc=-3
          endif
         else
          if(errc.ne.DSVP_SUCCESS) errc=-4
         endif
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

         errc=0
         if(this%is_active(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
           if(this%get_comm_stat().ne.DS_OPRND_NO_COMM) then
            delivered=this%sync(errc,wait=.TRUE.)
            if(.not.delivered.or.errc.ne.DSVP_SUCCESS) errc=-1
           endif
           call this%mark_empty(ier); if(ier.ne.DSVP_SUCCESS.and.errc.eq.0) errc=-2
           if(.not.talsh_tensor_is_empty(this%talsh_tens)) then
            ier=talsh_tensor_destruct(this%talsh_tens)
            if(ier.ne.TALSH_SUCCESS.and.errc.eq.0) errc=-3
           endif
           if(associated(this%resource)) then
            call this%resource%decr_ref_count()
            this%resource=>NULL()
           endif
           this%tensor=>NULL()
          else
           errc=-4
          endif
         else
          if(errc.ne.DSVP_SUCCESS) errc=-5
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
         if(errc.ne.0) call quit(errc,'#FATAL(TAVP-WRK:tens_oprnd_dtor): Destructor failed!')
         return
        end subroutine tens_oprnd_dtor
![tens_entry_wrk_t]==================================
        subroutine TensEntryWrkCtor(this,tensor,ierr)
!Constructs a <tens_entry_wrk_t>.
         implicit none
         class(tens_entry_wrk_t), intent(out):: this          !out: specialized tensor cache entry
         class(tens_rcrsv_t), pointer, intent(inout):: tensor !inout: pointer to an allocated tensor (ownership transfer may occur)
         integer(INTD), intent(out), optional:: ierr          !out: error code
         integer(INTD):: errc

         errc=0
         if(associated(tensor)) then
          this%tensor=>tensor; tensor=>NULL() !transfer the ownership
          allocate(this%resource,STAT=errc); if(errc.ne.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensEntryWrkCtor
!---------------------------------------------------------------------
        function TensEntryWrkGetResource(this,ierr) result(resource_p)
!Returns a pointer to the tensor cache entry resource.
         implicit none
         class(tens_resrc_t), pointer:: resource_p   !out: pointer to the previously allocated resource
         class(tens_entry_wrk_t), intent(in):: this  !in: tensor cache entry
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         resource_p=>NULL()
         if(this%is_set(errc)) then
          if(errc.eq.0) then
           resource_p=>this%resource
           if(.not.associated(resource_p)) errc=-3
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensEntryWrkGetResource
!-------------------------------------------
        subroutine tens_entry_wrk_dtor(this)
         implicit none
         type(tens_entry_wrk_t):: this

         if(associated(this%resource)) deallocate(this%resource)
         if(associated(this%tensor)) deallocate(this%tensor)
         return
        end subroutine tens_entry_wrk_dtor
!-------------------------------------------------------------
        function tens_entry_wrk_alloc(tens_entry) result(ierr)
!Non-member allocator for tens_entry_wrk_t.
         implicit none
         integer(INTD):: ierr
         class(*), allocatable, intent(out):: tens_entry

         allocate(tens_entry_wrk_t::tens_entry,STAT=ierr)
         return
        end function tens_entry_wrk_alloc
![tens_instr_t]============================================
        subroutine TensInstrCtor(this,op_code,ierr,op_spec)
!Constructs a tensor instruction from a given tensor operation.
!The tensor instruction is a realization of a given tensor operation
!for a specific TAVP kind.
         implicit none
         class(tens_instr_t), intent(inout):: this        !out: tensor instruction (must be empty on entrance)
         integer(INTD), intent(in):: op_code              !in: instruction code (see top)
         integer(INTD), intent(out), optional:: ierr      !out: error code
         class(*), intent(in), target, optional:: op_spec !in: formal operation specification
         integer(INTD):: errc

         if(this%is_empty(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
!Construct the instruction:
           select case(op_code)
           case(TAVP_INSTR_NOOP)
           case(TAVP_INSTR_STOP)
           case(TAVP_INSTR_CREATE,TAVP_INSTR_DESTROY)
            call construct_instr_create_destroy(errc); if(errc.ne.0) errc=-6
           case(TAVP_INSTR_CONTRACT)
            call construct_instr_contract(errc); if(errc.ne.0) errc=-5
           case default
            errc=-4 !invalid operation (or not implemented)
           end select
!Activate the instruction:
           if(errc.eq.0) then
            call this%activate(op_code,microcode(op_code),errc); if(errc.ne.0) errc=-3
           else
            call this%set_status(DS_INSTR_RETIRED,errc,-1)
            call tens_instr_dtor(this)
           endif
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return

        contains

         subroutine construct_instr_create_destroy(jerr)
          !CREATE/DESTROY a tensor:
          !op_spec={tens_rcrsv_t}
          integer(INTD), intent(out):: jerr
          class(ds_oprnd_t), pointer:: oprnd
          class(tens_oprnd_t), pointer:: tens_oprnd
          class(tens_rcrsv_t), pointer:: tensor

          jerr=0
          tensor=>NULL()
          select type(op_spec); class is(tens_rcrsv_t); tensor=>op_spec; end select
          if(associated(tensor)) then
           if(tensor%is_set()) then
            call this%alloc_operands(1,jerr)
            if(jerr.eq.DSVP_SUCCESS) then
             allocate(tens_oprnd,STAT=jerr)
             if(jerr.eq.0) then
              call tens_oprnd%tens_oprnd_ctor(tensor,jerr)
              if(jerr.eq.0) then
               oprnd=>tens_oprnd
               call this%set_operand(0,oprnd,jerr)
               if(jerr.ne.DSVP_SUCCESS) jerr=-6
              else
               jerr=-5
              endif
              oprnd=>NULL() !<oprnd> pointer was saved in the tensor instruction and will later be deallocated
              tens_oprnd=>NULL()
             else
              jerr=-4
             endif
            else
             jerr=-3
            endif
           else
            jerr=-2
           endif
           tensor=>NULL() !<tensor> pointed to an external object
          else
           jerr=-1
          endif
          return
         end subroutine construct_instr_create_destroy

         subroutine construct_instr_contract(jerr)
          !CONTRACT two tensors into another tensor:
          !op_spec={tens_contraction_t}
          integer(INTD), intent(out):: jerr
          integer(INTD):: jj
          class(ds_oprnd_t), pointer:: oprnd
          class(ds_instr_ctrl_t), pointer:: instr_ctrl
          class(contr_ptrn_ext_t), pointer:: contr_ptrn
          class(tens_contraction_t), pointer:: tens_contr
          class(ctrl_tens_contr_t), pointer:: tens_contr_ctrl
          class(tens_rcrsv_t), pointer:: tensor
          class(tens_oprnd_t), pointer:: tens_oprnd

          jerr=0
          tens_contr=>NULL()
          select type(op_spec); class is(tens_contraction_t); tens_contr=>op_spec; end select
          if(associated(tens_contr)) then
           if(tens_contr%is_set()) then
            contr_ptrn=>tens_contr%get_ext_contr_ptrn(jerr)
            if(jerr.eq.TEREC_SUCCESS) then
             allocate(tens_contr_ctrl,STAT=jerr)
             if(jerr.eq.0) then
              call tens_contr_ctrl%ctrl_tens_contr_ctor(contr_ptrn,jerr,tens_contr%get_prefactor()) !contraction pattern is cloned by value
              if(jerr.eq.0) then
               instr_ctrl=>tens_contr_ctrl
               call this%set_control(instr_ctrl,jerr)
               if(jerr.eq.DSVP_SUCCESS) then
                call this%alloc_operands(3,jerr)
                if(jerr.eq.DSVP_SUCCESS) then
                 do jj=0,2
                  tensor=>tens_contr%get_argument(jj,jerr); if(jerr.ne.TEREC_SUCCESS) exit
                  allocate(tens_oprnd,STAT=jerr); if(jerr.ne.0) exit
                  call tens_oprnd%tens_oprnd_ctor(tensor,jerr); if(jerr.ne.0) exit
                  oprnd=>tens_oprnd
                  call this%set_operand(jj,oprnd,jerr); if(jerr.ne.DSVP_SUCCESS) exit
                  tensor=>NULL(); tens_oprnd=>NULL(); oprnd=>NULL() !<oprnd> pointer was saved in the tensor instruction and will later be deallocated
                 enddo
                else
                 jerr=-7
                endif
               else
                jerr=-6
               endif
              else
               jerr=-5
              endif
              instr_ctrl=>NULL(); tens_contr_ctrl=>NULL() !<tens_contr_ctrl> pointer was saved in the tensor instruction and will later be deallocated
             else
              jerr=-4
             endif
             contr_ptrn=>NULL()
            else
             jerr=-3
            endif
           else
            jerr=-2
           endif
           tens_contr=>NULL()
          else
           jerr=-1
          endif
          return
         end subroutine construct_instr_contract

        end subroutine TensInstrCtor
!---------------------------------------------------------
        subroutine TensInstrEncode(this,instr_packet,ierr)
!Encodes a tensor instruction into the bytecode packet.
         implicit none
         class(tens_instr_t), intent(in):: this          !in: defined tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet !out: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,op_code

!Pack the instruction code (op_code):
         if(.not.this%is_empty(errc)) then
          op_code=this%get_code(errc)
          if(errc.eq.DSVP_SUCCESS) then
           call pack_builtin(instr_packet,op_code,errc)
           if(errc.eq.0) then
!Pack the instruction body:
            select case(op_code)
            case(TAVP_INSTR_NOOP)
            case(TAVP_INSTR_STOP)
            case(TAVP_INSTR_CREATE,TAVP_INSTR_DESTROY)
             call encode_instr_create_destroy(errc); if(errc.ne.0) errc=-6
            case(TAVP_INSTR_CONTRACT)
             call encode_instr_contract(errc); if(errc.ne.0) errc=-5
            case default
             errc=-4
            end select
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

        contains

         subroutine encode_instr_create_destroy(jerr)
          !CREATE/DESTROY a tensor:
          !Packet format: {op_code|tensor}
          integer(INTD), intent(out):: jerr
          class(ds_oprnd_t), pointer:: oprnd
          class(tens_rcrsv_t), pointer:: tensor

          jerr=0
          oprnd=>this%get_operand(0,jerr)
          if(jerr.eq.DSVP_SUCCESS) then
           select type(oprnd)
           class is(tens_oprnd_t)
            tensor=>oprnd%get_tensor(jerr)
            if(jerr.eq.0) then
             if(tensor%is_set()) then
              call tensor%pack(instr_packet,jerr)
              if(jerr.ne.0) jerr=-5
             else
              jerr=-4
             endif
             tensor=>NULL()
            else
             jerr=-3
            endif
           class default
            jerr=-2
           end select
           oprnd=>NULL()
          else
           jerr=-1
          endif
          return
         end subroutine encode_instr_create_destroy

         subroutine encode_instr_contract(jerr)
          !CONTRACT two tensors: tensor0+=tensor1*tensor2*scalar:
          !Packed format: {op_code|ctrl_tens_contr_t|tensor0,tensor1,tensor2}
          integer(INTD), intent(out):: jerr
          integer(INTD):: jj
          class(ds_oprnd_t), pointer:: oprnd
          class(tens_rcrsv_t), pointer:: tensor
          class(ds_instr_ctrl_t), pointer:: tens_contr_ctrl

          jerr=0
          tens_contr_ctrl=>this%get_control(jerr)
          if(jerr.eq.DSVP_SUCCESS) then
           select type(tens_contr_ctrl)
           class is(ctrl_tens_contr_t)
            call tens_contr_ctrl%pack(instr_packet,jerr)
            if(jerr.ne.0) jerr=-8
           class default
            jerr=-7
           end select
           if(jerr.eq.0) then
            do jj=0,2
             oprnd=>this%get_operand(jj,jerr); if(jerr.ne.DSVP_SUCCESS) then; jerr=-6; exit; endif
             select type(oprnd)
             class is(tens_oprnd_t)
              tensor=>oprnd%get_tensor(jerr); if(jerr.ne.0) then; jerr=-5; exit; endif
              if(.not.tensor%is_set()) then; jerr=-4; exit; endif !trap
              call tensor%pack(instr_packet,jerr); if(jerr.ne.0) then; jerr=-3; exit; endif
             class default
              jerr=-2; exit
             end select
            enddo
            tensor=>NULL(); oprnd=>NULL()
           endif
           tens_contr_ctrl=>NULL()
          else
           jerr=-1
          endif
          return
         end subroutine encode_instr_contract

        end subroutine TensInstrEncode
!---------------------------------------
        subroutine tens_instr_dtor(this)
         implicit none
         type(tens_instr_t):: this !inout: empty or retired tensor instruction
         integer(INTD):: sts,errc

         sts=this%get_status(errc)
         if((sts.eq.DS_INSTR_EMPTY.or.sts.eq.DS_INSTR_RETIRED).and.errc.eq.DSVP_SUCCESS) then
          call this%clean(errc)
          if(errc.ne.0) call quit(errc,'#FATAL(tavp_worker:tens_instr_dtor): Tensor instruction destruction failed!')
         else
          call quit(-1,'#FATAL(tavp_worker:tens_instr_dtor): Attempt to destroy an active TAVP instruction!')
         endif
         return
        end subroutine tens_instr_dtor
![tavp_worker_t]=============================
        subroutine TAVPWorkerStart(this,ierr)
!Initializes TAVP "Worker".
         implicit none
         class(tavp_worker_t), intent(inout):: this  !inout: TAVP "Worker"
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         call init_talsh(errc)
         if(errc.eq.0) then
          call init_distributed_space(errc)
          if(errc.eq.0) then
           call init_microcode(errc)
           if(errc.eq.0) then
            !call test_carma(errc) !debug: initial brute-force distributed benchmark
            call this%shutdown(errc)
           else
            errc=-3
           endif
          else
           errc=-2
           call quit(errc,'FATAL(TAVP-Worker:start): Failed to initialize the distributed address space service!')
          endif
         else
          errc=-1
          call quit(errc,'FATAL(TAVP-Worker:start): Failed to initialize the numerical tensor algebra service!')
         endif
         if(present(ierr)) ierr=errc
         return

        contains

         subroutine init_talsh(jerr)
          implicit none
          integer(INTD), intent(out):: jerr
          integer(INTD):: jj

          jerr=talsh_init(tavp_worker_host_buf_size,tavp_worker_host_arg_max,(/(jj,jj=gpu_start,gpu_start+gpu_count-1)/))
          return
         end subroutine init_talsh

         subroutine init_distributed_space(jerr)
          implicit none
          integer(INTD), intent(out):: jerr

          call tavp_addr_space%create(role_comm,TAVP_WORKER_NUM_WINS,'WorkAddressSpace',jerr)
          return
         end subroutine init_distributed_space

        end subroutine TAVPWorkerStart
!-----------------------------------------------
        subroutine TAVPWorkerShutdown(this,ierr)
!Shuts down TAVP "Worker".
         implicit none
         class(tavp_worker_t), intent(inout):: this  !inout: TAVP "Worker"
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         call null_microcode(errc)
         if(errc.eq.0) then
          call stop_distributed_space(errc)
          if(errc.eq.0) then
           call stop_talsh(errc)
           if(errc.eq.0) then
           else
            errc=-1
           endif
          else
           errc=-2
          endif
         else
          errc=-3
         endif
         if(present(ierr)) ierr=errc
         return

        contains

         subroutine stop_talsh(jerr)
          implicit none
          integer(INTD), intent(out):: jerr

          jerr=talsh_shutdown()
          return
         end subroutine stop_talsh

         subroutine stop_distributed_space(jerr)
          implicit none
          integer(INTD), intent(out):: jerr

          call tavp_addr_space%destroy(jerr)
          return
         end subroutine stop_distributed_space

         subroutine null_microcode(jerr)
          implicit none
          integer(INTD), intent(out):: jerr
          integer(INTD):: jj

          jerr=0
          do jj=0,TAVP_ISA_SIZE-1
           call microcode(jj)%reset()
          enddo
          return
         end subroutine null_microcode

        end subroutine TAVPWorkerShutdown
!------------------------------------------------------------------------------
        subroutine TAVPWorkerDecodeInstruction(this,instr_packet,ds_instr,ierr)
!Decodes a tensor instruction from the bytecode packet.
         implicit none
         class(tavp_worker_t), intent(inout):: this          !inout: TAVP "Worker"
         class(obj_pack_t), intent(inout):: instr_packet     !in: instruction bytecode packet
         class(ds_instr_t), intent(inout), target:: ds_instr !out: tensor instruction (must be empty on entrance)
         integer(INTD), intent(out), optional:: ierr         !out: error code
         integer(INTD):: errc,ier,op_code

         if(ds_instr%is_empty(errc)) then
          if(errc.eq.DSVP_SUCCESS) then
!Extract the instruction op_code:
           call unpack_builtin(instr_packet,op_code,errc)
!Extract the instruction body:
           if(errc.eq.0) then
            select case(op_code)
            case(TAVP_INSTR_NOOP)
            case(TAVP_INSTR_STOP)
            case(TAVP_INSTR_CREATE,TAVP_INSTR_DESTROY)
             call decode_instr_create_destroy(errc); if(errc.ne.0) errc=-7
            case(TAVP_INSTR_CONTRACT)
             call decode_instr_contract(errc); if(errc.ne.0) errc=-6
            case default
             errc=-5 !unknown instruction (or not implemented)
            end select
!Activate the instruction:
            call ds_instr%activate(op_code,microcode(op_code),ier); if(ier.ne.0.and.errc.eq.0) errc=-4
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

         contains

          subroutine decode_instr_create_destroy(jerr)
           !CREATE/DESTROY a tensor
           integer(INTD), intent(out):: jerr
           class(tens_rcrsv_t), pointer:: tensor
           class(tens_entry_t), pointer:: tens_entry
           class(tens_resrc_t), pointer:: tens_resource
           class(tens_oprnd_t), pointer:: tens_oprnd
           class(ds_oprnd_t), pointer:: oprnd
           logical:: res

           jerr=0
           tensor=>NULL(); allocate(tensor,STAT=jerr)
           if(jerr.eq.0) then
            call tensor%tens_rcrsv_ctor(instr_packet,jerr)
            if(jerr.eq.TEREC_SUCCESS) then
             tens_entry=>NULL()
             tens_entry=>this%tens_cache%lookup(tensor,jerr)
             if(jerr.eq.0) then
              select case(op_code)
              case(TAVP_INSTR_CREATE) !CREATE a tensor
               if(.not.associated(tens_entry)) then
                res=this%tens_cache%store(tensor,jerr,tens_entry_p=tens_entry)
                if(res.and.(jerr.eq.0).and.associated(tens_entry)) then
                 tens_resource=>NULL()
                 tens_resource=>tens_entry%get_resource()
                 call ds_instr%alloc_operands(1,jerr)
                 if(jerr.eq.DSVP_SUCCESS) then
                  tens_oprnd=>NULL(); allocate(tens_oprnd,STAT=jerr)
                  if(jerr.eq.0) then
                   call tens_oprnd%tens_oprnd_ctor(tensor,jerr,tens_resource) !tensor and tens_resource are stored in tensor cache
                   if(jerr.eq.0) then
                    oprnd=>tens_oprnd
                    call ds_instr%set_operand(0,oprnd,jerr)
                    if(jerr.ne.DSVP_SUCCESS) then
                     call ds_instr%terminate(TAVP_ERR_GEN_FAILURE,jerr)
                     jerr=-14
                    endif
                   else
                    call ds_instr%terminate(TAVP_ERR_GEN_FAILURE,jerr)
                    jerr=-13
                   endif
                  else
                   call ds_instr%terminate(TAVP_ERR_RSC_UNAVAILABLE,jerr)
                   jerr=-12
                  endif
                 else
                  call ds_instr%terminate(TAVP_ERR_RSC_UNAVAILABLE,jerr)
                  jerr=-11
                 endif
                else
                 call ds_instr%terminate(TAVP_ERR_CHE_FAILURE,jerr)
                 jerr=-10
                endif
               else
                call ds_instr%terminate(TAVP_ERR_ARG_DEFINED,jerr)
                jerr=-9
               endif
              case(TAVP_INSTR_DESTROY) !DESTROY a tensor
               if(associated(tens_entry)) then
                deallocate(tensor); tensor=>NULL() !deallocate the temporary tensor and use
                tensor=>tens_entry%get_tensor()    !the same tensor from the tensor cache
                tens_resource=>NULL()              !as well as its associated resource
                tens_resource=>tens_entry%get_resource()
                call ds_instr%alloc_operands(1,jerr)
                if(jerr.eq.DSVP_SUCCESS) then
                 tens_oprnd=>NULL(); allocate(tens_oprnd,STAT=jerr)
                 if(jerr.eq.0) then
                  call tens_oprnd%tens_oprnd_ctor(tensor,jerr,tens_resource) !tensor and tens_resource are from tensor cache
                  if(jerr.eq.0) then
                   oprnd=>tens_oprnd
                   call ds_instr%set_operand(0,oprnd,jerr)
                   if(jerr.ne.DSVP_SUCCESS) then
                    call ds_instr%terminate(TAVP_ERR_GEN_FAILURE,jerr)
                    jerr=-8
                   endif
                  else
                   call ds_instr%terminate(TAVP_ERR_GEN_FAILURE,jerr)
                   jerr=-7
                  endif
                 else
                  call ds_instr%terminate(TAVP_ERR_RSC_UNAVAILABLE,jerr)
                  jerr=-6
                 endif
                else
                 call ds_instr%terminate(TAVP_ERR_RSC_UNAVAILABLE,jerr)
                 jerr=-5
                endif
               else
                call ds_instr%terminate(TAVP_ERR_ARG_UNDEFINED,jerr)
                jerr=-4
               endif
              end select
             else
              call ds_instr%terminate(TAVP_ERR_CHE_FAILURE,jerr)
              jerr=-3
             endif
            else
             call ds_instr%terminate(TAVP_ERR_BTC_BAD,jerr)
             jerr=-2
            endif
           else
            call ds_instr%terminate(TAVP_ERR_RSC_UNAVAILABLE,jerr)
            jerr=-1
           endif
           return
          end subroutine decode_instr_create_destroy

          subroutine decode_instr_contract(jerr)
           !CONTRACT two tensors
           integer(INTD), intent(out):: jerr
           class(ds_instr_ctrl_t), pointer:: instr_ctrl
           class(ctrl_tens_contr_t), pointer:: tens_contr_ctrl
           class(tens_rcrsv_t), pointer:: tensor
           class(tens_resrc_t), pointer:: tens_resource
           class(tens_entry_t), pointer:: tens_entry
           class(tens_oprnd_t), pointer:: tens_oprnd
           class(ds_oprnd_t), pointer:: oprnd
           integer(INTD):: jj

           jerr=0
           tens_contr_ctrl=>NULL(); allocate(tens_contr_ctrl,STAT=jerr)
           if(jerr.eq.0) then
            call tens_contr_ctrl%unpack(instr_packet,jerr)
            if(jerr.eq.0) then
             instr_ctrl=>tens_contr_ctrl
             call ds_instr%set_control(instr_ctrl,jerr)
             if(jerr.eq.DSVP_SUCCESS) then
              call ds_instr%alloc_operands(3,jerr)
              if(jerr.eq.DSVP_SUCCESS) then
               do jj=0,2
                tensor=>NULL(); allocate(tensor,STAT=jerr)
                if(jerr.ne.0) then; call ds_instr%terminate(TAVP_ERR_RSC_UNAVAILABLE,jerr); jerr=-11; exit; endif
                call tensor%tens_rcrsv_ctor(instr_packet,jerr)
                if(jerr.ne.TEREC_SUCCESS) then; call ds_instr%terminate(TAVP_ERR_BTC_BAD,jerr); jerr=-10; exit; endif
                tens_entry=>this%tens_cache%lookup(tensor,jerr)
                if(jerr.ne.0) then; call ds_instr%terminate(TAVP_ERR_CHE_FAILURE,jerr); jerr=-9; exit; endif
                if(.not.associated(tens_entry)) then; call ds_instr%terminate(TAVP_ERR_ARG_UNDEFINED,jerr); jerr=-8; exit; endif
                deallocate(tensor); tensor=>NULL() !deallocate the temporary tensor and use
                tensor=>tens_entry%get_tensor()    !the same tensor from the tensor cache
                tens_resource=>NULL()              !as well as its associated resource
                tens_resource=>tens_entry%get_resource()
                tens_oprnd=>NULL(); allocate(tens_oprnd,STAT=jerr)
                if(jerr.ne.0) then; call ds_instr%terminate(TAVP_ERR_RSC_UNAVAILABLE,jerr); jerr=-7; exit; endif
                call tens_oprnd%tens_oprnd_ctor(tensor,jerr,tens_resource) !tensor and tens_resource are stored in tensor cache
                if(jerr.ne.0) then; call ds_instr%terminate(TAVP_ERR_GEN_FAILURE,jerr); jerr=-6; exit; endif
                oprnd=>tens_oprnd
                call ds_instr%set_operand(jj,oprnd,jerr)
                if(jerr.ne.DSVP_SUCCESS) then; call ds_instr%terminate(TAVP_ERR_GEN_FAILURE,jerr); jerr=-5; exit; endif
               enddo
              else
               call ds_instr%terminate(TAVP_ERR_RSC_UNAVAILABLE,jerr)
               jerr=-4
              endif
             else
              call ds_instr%terminate(TAVP_ERR_GEN_FAILURE,jerr)
              jerr=-3
             endif
            else
             call ds_instr%terminate(TAVP_ERR_BTC_BAD,jerr)
             jerr=-2
            endif
           else
            call ds_instr%terminate(TAVP_ERR_RSC_UNAVAILABLE,jerr)
            jerr=-1
           endif
           return
          end subroutine decode_instr_contract

        end subroutine TAVPWorkerDecodeInstruction
#endif
       end module tavp_worker
