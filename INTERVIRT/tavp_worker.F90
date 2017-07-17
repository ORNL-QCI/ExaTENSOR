!ExaTENSOR: TAVP "Worker" implementation
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017/07/16

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
        integer(INTD), private:: CONS_OUT=6
        integer(INTD), private:: DEBUG=0
        logical, private:: VERBOSE=.TRUE.
 !Distributed memory space (for workers):
        integer(INTD), parameter, private:: TAVP_WORKER_NUM_WINS=1 !number of MPI windows in the distributed space
 !On-node pinned Host memory buffer:
        integer(INTL), protected:: tavp_worker_host_buf_size=1_INTL*(1024_INTL*1024_INTL*1024_INTL) !size in bytes
        integer(INTD), private:: tavp_worker_host_arg_max=0 !set later: max number of tensors in the pinned Host buffer
 !Elementary tensor instruction granularity classification:
        real(8), public:: EXA_FLOPS_MEDIUM=1d10 !minimal number of Flops to consider the operation as medium-cost
        real(8), public:: EXA_FLOPS_HEAVY=1d12  !minimal number of Flops to consider the operation as heavy-cost
        real(8), public:: EXA_COST_TO_SIZE=1d2  !minimal cost (Flops) to size (Words) ratio to consider the operation arithmetically intensive
!TYPES:
 !Tensor cache entry:
        type, private:: tens_entry_t
         class(tens_rcrsv_t), pointer, private:: tensor=>NULL()   !tensor (either allocated or associated)
         class(tens_resrc_t), pointer, private:: resource=>NULL() !resource (either allocated or associated)
         logical, private:: tens_alloc=.FALSE.                    !TRUE if the tensor pointer is allocated
         logical, private:: res_alloc=.FALSE.                     !TRUE if the resource pointer is allocated
         contains
          procedure, private:: TensEntryCtor                      !ctor
          generic, public:: tens_entry_ctor=>TensEntryCtor
          procedure, public:: is_set=>TensEntryIsSet              !returns TRUE if the tensor cache entry is set
          procedure, public:: get_tensor=>TensEntryGetTensor      !returns a pointer to the tensor
          procedure, public:: get_resource=>TensEntryGetResource  !returns a pointer to the resource
          final:: tens_entry_dtor                                 !dtor
        end type tens_entry_t
 !Tensor cache:
        type, private:: tens_cache_t
         type(dictionary_t), private:: map                        !cache dictionary
         contains
          procedure, public:: lookup=>TensCacheLookup             !looks up a given tensor in the cache
          procedure, public:: store=>TensCacheStore               !stores a given tensor in the cache
          procedure, public:: evict=>TensCacheEvict               !evicts a given tensor from the cache
          procedure, public:: erase=>TensCacheErase               !erases everything from the cache (regardless of MPI communications)
          final:: tens_cache_dtor                                 !dtor
        end type tens_cache_t
 !Tensor instruction (realization of a tensor operation for a specific TAVP):
        type, extends(ds_instr_t), private:: tens_instr_t
         type(talsh_task_t), private:: talsh_task                   !TAL-SH task
         contains
          procedure, private:: TensInstrCtor                        !ctor: constructs a tensor instruction from the specification of a tensor operation
          generic, public:: tens_instr_ctor=>TensInstrCtor
          procedure, public:: decode=>TensInstrDecode               !decoding procedure: Unpacks the raw byte packet (bytecode) and constructs a TAVP instruction
          procedure, public:: encode=>TensInstrEncode               !encoding procedure: Packs the TAVP instruction into a raw byte packet (bytecode)
          procedure, private:: set_microcode=>TensInstrSetMicrocode !sets up instruction dynamic bindings to the corresponding microcode
          procedure, private:: activate=>TensInstrActivate          !activates the instruction
          final:: tens_instr_dtor                                   !dtor
        end type tens_instr_t
 !TAVP specialization "Worker":
        type, extends(dsvp_t), public:: tavp_worker_t
         type(tens_cache_t), private:: tens_cache       !tensor cache (both persistent and temporary tensors)
         type(list_bi_t), private:: instr_queue         !instruction queue
         type(list_iter_t), private:: instr_it          !instruction queue iterator
         type(vector_t), private:: events               !events for dependency enforcement
         type(vector_iter_t), private:: event_it        !events iterator
         contains
          procedure, public:: start=>TAVPWorkerStart                                             !initializes TAVP to an active state
          procedure, public:: shutdown=>TAVPWorkerShutdown                                       !shuts down TAVP
          procedure, public:: fetch_instructions=>TAVPWorkerFetchInstructions                    !fetches a block of tensor instructions from TAVP "Manager"
          procedure, public:: return_retired_instructions=>TAVPWorkerReturnRetiredInstructions   !returns back a block of retired tensor instructions to TAVP "Manager"
          procedure, public:: send_instructions=>TAVPWorkerSendInstructions                      !N/A (dummy)
          procedure, public:: receive_retired_instructions=>TAVPWorkerReceiveRetiredInstructions !N/A (dummy)
        end type tavp_worker_t
!INTERFACES:
!DATA:
 !TAVP instruction microcode bindings (set by the TAVP initialization):
        type(microcode_bind_t), private:: microcode(0:TAVP_ISA_SIZE-1)
 !TAVP "Worker":
        type(tavp_worker_t), protected:: tavpWorker
!VISIBILITY:
 !non-member:
        private test_carma
        private tavp_worker_set_host_buf_size
        private acquire_resource_empty
        private acquire_resource_basic
        private prefetch_input_empty
        private prefetch_input_basic
        private sync_prefetch_empty
        private sync_prefetch_basic
        private upload_output_empty
        private upload_output_basic
        private sync_upload_empty
        private sync_upload_basic
        private release_resource_empty
        private release_resource_basic
        private sync_execution_empty
        private sync_execution_basic
        private execute_empty
        private execute_tensor_create
        private execute_tensor_destroy
        private execute_tensor_contract
 !tens_entry_t:
        private TensEntryCtor
        private TensEntryIsSet
        private TensEntryGetTensor
        private TensEntryGetResource
        private tens_entry_dtor
 !tens_cache_t:
        private TensCacheLookup
        private TensCacheStore
        private TensCacheEvict
        private TensCacheErase
        public tens_cache_dtor
 !tens_instr_t:
        private TensInstrCtor
        private TensInstrDecode
        private TensInstrEncode
        private TensInstrSetMicrocode
        private TensInstrActivate
        private tens_instr_dtor
 !tavp_worker_t:
        private TAVPWorkerStart
        private TAVPWorkerShutdown
        private TAVPWorkerFetchInstructions
        private TAVPWorkerReturnRetiredInstructions
        private TAVPWorkerSendInstructions
        private TAVPWorkerReceiveRetiredInstructions

!IMPLEMENTATION:
       contains
![non-member]======================
        subroutine test_carma(ierr)
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
          j=mod((i-1)+role_rank**2,n)+1
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
!-----------------------------------------------------------
        subroutine tavp_worker_set_host_buf_size(bytes,ierr)
!Changes the default size of the pinned Host buffer.
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
!---------------------------------------------------
        subroutine acquire_resource_empty(this,ierr)
!Dummy procedure for acquiring no resource.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine acquire_resource_empty
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
        subroutine prefetch_input_empty(this,ierr)
!Dummy procedure for prefetching no input.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine prefetch_input_empty
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
        function sync_prefetch_empty(this,ierr,wait) result(res)
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
        end function sync_prefetch_empty
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
        subroutine upload_output_empty(this,ierr)
!Dummy procedure for uploading no output.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine upload_output_empty
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
        function sync_upload_empty(this,ierr,wait) result(res)
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
        end function sync_upload_empty
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
        subroutine release_resource_empty(this,ierr)
!Dummy procedure for releasing no resource.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine release_resource_empty
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
        function sync_execution_empty(this,ierr,wait) result(res)
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
        end function sync_execution_empty
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
        subroutine execute_empty(this,ierr)
!Dummy procedure for executing nothing.
         implicit none
         class(ds_instr_t), intent(inout):: this     !inout: tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine execute_empty
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
![tens_entry_t]========================================================================
        subroutine TensEntryCtor(this,ierr,tensor,resource,tensor_alloc,resource_alloc)
!CTOR: If either <tensor> or <resource> are not present,
!the corresponding components of <this> will be allocated,
!otherwise pointer associated. However, <tensor_alloc> and
!<resource_alloc> flags can be used to mark the corresponding
!associated objects as allocated in case they will need to be
!deallocated later from here.
         implicit none
         class(tens_entry_t), intent(out):: this                       !out: tensor cache entry
         integer(INTD), intent(out), optional:: ierr                   !out: error code
         class(tens_rcrsv_t), intent(in), pointer, optional:: tensor   !in: pointer to the tensor
         class(tens_resrc_t), intent(in), pointer, optional:: resource !in: pointer to the resource
         logical, intent(in), optional:: tensor_alloc                  !in: if TRUE, the associated tensor will be assumed ALLOCATED
         logical, intent(in), optional:: resource_alloc                !in: if TRUE, the associated resource will be assumed ALLOCATED
         integer(INTD):: errc

         errc=0
         if(present(tensor)) then
          if(associated(tensor)) then
           this%tensor=>tensor; this%tens_alloc=.FALSE.
           if(present(tensor_alloc)) this%tens_alloc=tensor_alloc
          else
           errc=-1
          endif
         else
          allocate(this%tensor,STAT=errc)
          if(errc.eq.0) then; this%tens_alloc=.TRUE.; else; errc=-2; endif
         endif
         if(errc.eq.0) then
          if(present(resource)) then
           if(associated(resource)) then
            this%resource=>resource; this%res_alloc=.FALSE.
            if(present(resource_alloc)) this%res_alloc=resource_alloc
           else
            errc=-3
           endif
          else
           allocate(this%resource,STAT=errc)
           if(errc.eq.0) then; this%res_alloc=.TRUE.; else; errc=-4; endif
          endif
         endif
         if(errc.ne.0) call tens_entry_dtor(this)
         if(present(ierr)) ierr=errc
         return
        end subroutine TensEntryCtor
!-----------------------------------------------------
        function TensEntryIsSet(this,ierr) result(ans)
         implicit none
         logical:: ans                               !out: answer
         class(tens_entry_t), intent(in):: this      !in: tensor cache entry
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         ans=associated(this%tensor)
         if(present(ierr)) ierr=errc
         return
        end function TensEntryIsSet
!--------------------------------------------------------------
        function TensEntryGetTensor(this,ierr) result(tensor_p)
         implicit none
         class(tens_rcrsv_t), pointer:: tensor_p     !out: pointer to the tensor
         class(tens_entry_t), intent(in):: this      !in: tensor cache entry
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         tensor_p=>this%tensor
         if(.not.associated(tensor_p)) errc=-1
         if(present(ierr)) ierr=errc
         return
        end function TensEntryGetTensor
!------------------------------------------------------------------
        function TensEntryGetResource(this,ierr) result(resource_p)
         implicit none
         class(tens_resrc_t), pointer:: resource_p   !out: pointer to the resource
         class(tens_entry_t), intent(in):: this      !in: tensor cache entry
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         resource_p=>this%resource
         if(.not.associated(resource_p)) errc=-1
         if(present(ierr)) ierr=errc
         return
        end function TensEntryGetResource
!---------------------------------------
        subroutine tens_entry_dtor(this)
         implicit none
         type(tens_entry_t):: this

         if(associated(this%tensor).and.this%tens_alloc) deallocate(this%tensor)
         if(associated(this%resource).and.this%res_alloc) deallocate(this%resource)
         this%tensor=>NULL(); this%resource=>NULL()
         this%tens_alloc=.FALSE.; this%res_alloc=.FALSE.
         return
        end subroutine tens_entry_dtor
![tens_cache_t]========================================================
        function TensCacheLookup(this,tensor,ierr) result(tens_entry_p)
!Looks up a given tensor in the tensor cache. If found, returns a pointer
!to the corresponding tensor cache entry. If not found, returns NULL.
         implicit none
         class(tens_entry_t), pointer:: tens_entry_p !out: pointer to the tensor cache entry or NULL
         class(tens_cache_t), intent(in):: this      !in: tensor cache
         class(tens_rcrsv_t), intent(in):: tensor    !in: tensor to look up (via its descriptor as the key)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,res
         class(*), pointer:: uptr
         type(dictionary_iter_t):: dit
         type(tens_descr_t):: tens_descr

         tens_entry_p=>NULL()
         errc=dit%init(this%map)
         if(errc.eq.GFC_SUCCESS) then
          uptr=>NULL()
          tens_descr=tensor%get_descriptor(errc)
          if(errc.eq.TEREC_SUCCESS) then
           res=dit%search(GFC_DICT_FETCH_IF_FOUND,cmp_tens_descriptors,tens_descr,value_out=uptr)
           if(res.eq.GFC_FOUND) then
            select type(uptr); class is(tens_entry_t); tens_entry_p=>uptr; end select
            if(.not.associated(tens_entry_p)) errc=-5 !trap
           else
            if(res.ne.GFC_NOT_FOUND) errc=-4
           endif
          else
           errc=-3
          endif
          res=dit%release(); if(res.ne.GFC_SUCCESS.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensCacheLookup
!-------------------------------------------------------------------------------------
        function TensCacheStore(this,tensor,ierr,resource,tens_entry_p) result(stored)
!Stores a tensor (and optionally its associated resource) in the tensor cache,
!unless the tensor is already present in the tensor cache. In any case, an optional
!<tens_entry_p> will point to either existing or newly created tensor cache entry.
!If the resource is absent, and empty resource will be allocated for the newly created
!tensor cache entry. <stored> is TRUE on successful new storage, FALSE otherwise.
!In any case, the <tensor> and <resource> components in the tensor cache entry
!will be marked as ALLOCATED, that is, both should be ALLOCATED pointers!
         implicit none
         logical:: stored                                                   !out: TRUE on successful new store, FALSE otherwise
         class(tens_cache_t), intent(inout):: this                          !inout: tensor cache
         class(tens_rcrsv_t), pointer, intent(in):: tensor                  !in: tensor (must be allocated)
         integer(INTD), intent(out), optional:: ierr                        !out: error code
         class(tens_resrc_t), pointer, intent(in), optional:: resource      !in: resource associated with the tensor (must be allocated)
         class(tens_entry_t), pointer, intent(out), optional:: tens_entry_p !out: tensor cache entry (new or existing)
         integer(INTD):: errc,res
         class(*), pointer:: uptr
         type(dictionary_iter_t):: dit
         type(tens_descr_t):: tens_descr
         type(tens_entry_t):: tens_entry

         stored=.FALSE.
         errc=dit%init(this%map)
         if(errc.eq.GFC_SUCCESS) then
          if(associated(tensor)) then
           tens_descr=tensor%get_descriptor(errc)
           if(errc.eq.TEREC_SUCCESS) then
            if(present(resource)) then
             call tens_entry%tens_entry_ctor(errc,tensor,resource,tensor_alloc=.TRUE.,resource_alloc=.TRUE.)
            else
             call tens_entry%tens_entry_ctor(errc,tensor,tensor_alloc=.TRUE.)
            endif
            if(errc.eq.0) then
             uptr=>NULL()
             res=dit%search(GFC_DICT_ADD_IF_NOT_FOUND,cmp_tens_descriptors,tens_descr,tens_entry,GFC_BY_VAL,value_out=uptr)
             if(res.eq.GFC_NOT_FOUND) then; stored=.TRUE.; else; if(res.ne.GFC_FOUND) errc=-7; endif
             if(present(tens_entry_p).and.errc.eq.0) then
              tens_entry_p=>NULL()
              select type(uptr); class is(tens_entry_t); tens_entry_p=>uptr; end select
              if(.not.associated(tens_entry_p)) errc=-6 !trap
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
          res=dit%release(); if(res.ne.0.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensCacheStore
!----------------------------------------------------------------
        function TensCacheEvict(this,tensor,ierr) result(evicted)
!Evicts a specific tensor cache entry from the tensor cache.
!If the corresponding tensor cache entry is not found,
!no error is risen, but <evicted>=FALSE.
         implicit none
         logical:: evicted                           !TRUE if the tensor cache entry was found and evicted, FALSE otherwise
         class(tens_cache_t), intent(inout):: this   !inout: tensor cache
         class(tens_rcrsv_t), intent(in):: tensor    !in: tensor to find and evict from the cache
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,res
         type(dictionary_iter_t):: dit
         type(tens_descr_t):: tens_descr

         evicted=.FALSE.
         errc=dit%init(this%map)
         if(errc.eq.GFC_SUCCESS) then
          tens_descr=tensor%get_descriptor(errc)
          if(errc.eq.TEREC_SUCCESS) then
           res=dit%search(GFC_DICT_DELETE_IF_FOUND,cmp_tens_descriptors,tens_descr)
           if(res.eq.GFC_FOUND) then; evicted=.TRUE.; else; if(res.ne.GFC_NOT_FOUND) errc=-4; endif
          else
           errc=-3
          endif
          res=dit%release(); if(res.ne.0.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end function TensCacheEvict
!-------------------------------------------
        subroutine TensCacheErase(this,ierr)
!Erases the tensor cache completely and unconditionally.
!In particular, if some resource objects are still in use,
!their deallocation will probably cause a crash.
         implicit none
         class(tens_cache_t), intent(inout):: this   !inout: tensor cache
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,res
         type(dictionary_iter_t):: dit

         errc=dit%init(this%map)
         if(errc.eq.GFC_SUCCESS) then
          call dit%delete_all(errc); if(errc.ne.GFC_SUCCESS) errc=-3
          res=dit%release(); if(res.ne.0.and.errc.eq.0) errc=-2
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensCacheErase
!---------------------------------------
        subroutine tens_cache_dtor(this)
         implicit none
         type(tens_cache_t):: this

         call this%erase()
         return
        end subroutine tens_cache_dtor
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
            call this%activate(op_code,errc); if(errc.ne.0) errc=-3
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
        subroutine TensInstrDecode(this,instr_packet,ierr)
!Decodes a tensor instruction from the bytecode packet.
         implicit none
         class(tens_instr_t), intent(inout):: this       !out: empty tensor instruction
         class(obj_pack_t), intent(inout):: instr_packet !in: instruction bytecode packet
         integer(INTD), intent(out), optional:: ierr     !out: error code
         integer(INTD):: errc,ier,op_code

         if(this%is_empty(errc)) then
!Extract the instruction op_code:
          call unpack_builtin(instr_packet,op_code,errc)
!Extract the instruction body:
          if(errc.eq.0) then
           select case(op_code)
           case(TAVP_INSTR_NOOP)
           case(TAVP_INSTR_STOP)
           case(TAVP_INSTR_CREATE,TAVP_INSTR_DESTROY)
            call decode_instr_create_destroy(errc); if(errc.ne.0) errc=-6
           case(TAVP_INSTR_CONTRACT)
            call decode_instr_contract(errc); if(errc.ne.0) errc=-5
           case default
            errc=-4 !unknown instruction (or not implemented)
           end select
!Activate the instruction:
           call this%activate(op_code,ier); if(ier.ne.0.and.errc.eq.0) errc=-3
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
             tens_entry=>tavpWorker%tens_cache%lookup(tensor,jerr)
             if(jerr.eq.0) then
              select case(op_code)
              case(TAVP_INSTR_CREATE) !CREATE a tensor
               if(.not.associated(tens_entry)) then
                res=tavpWorker%tens_cache%store(tensor,jerr,tens_entry_p=tens_entry)
                if(res.and.(jerr.eq.0).and.associated(tens_entry)) then
                 tens_resource=>NULL()
                 tens_resource=>tens_entry%get_resource()
                 call this%alloc_operands(1,jerr)
                 if(jerr.eq.DSVP_SUCCESS) then
                  tens_oprnd=>NULL(); allocate(tens_oprnd,STAT=jerr)
                  if(jerr.eq.0) then
                   call tens_oprnd%tens_oprnd_ctor(tensor,jerr,tens_resource) !tensor and tens_resource are stored in tensor cache
                   if(jerr.eq.0) then
                    oprnd=>tens_oprnd
                    call this%set_operand(0,oprnd,jerr)
                    if(jerr.ne.DSVP_SUCCESS) then
                     call this%terminate(TAVP_ERR_GEN_FAILURE,jerr)
                     jerr=-14
                    endif
                   else
                    call this%terminate(TAVP_ERR_GEN_FAILURE,jerr)
                    jerr=-13
                   endif
                  else
                   call this%terminate(TAVP_ERR_RSC_UNAVAILABLE,jerr)
                   jerr=-12
                  endif
                 else
                  call this%terminate(TAVP_ERR_RSC_UNAVAILABLE,jerr)
                  jerr=-11
                 endif
                else
                 call this%terminate(TAVP_ERR_CHE_FAILURE,jerr)
                 jerr=-10
                endif
               else
                call this%terminate(TAVP_ERR_ARG_DEFINED,jerr)
                jerr=-9
               endif
              case(TAVP_INSTR_DESTROY) !DESTROY a tensor
               if(associated(tens_entry)) then
                deallocate(tensor); tensor=>NULL() !deallocate the temporary tensor and use
                tensor=>tens_entry%get_tensor()    !the same tensor from the tensor cache
                tens_resource=>NULL()              !as well as its associated resource
                tens_resource=>tens_entry%get_resource()
                call this%alloc_operands(1,jerr)
                if(jerr.eq.DSVP_SUCCESS) then
                 tens_oprnd=>NULL(); allocate(tens_oprnd,STAT=jerr)
                 if(jerr.eq.0) then
                  call tens_oprnd%tens_oprnd_ctor(tensor,jerr,tens_resource) !tensor and tens_resource are from tensor cache
                  if(jerr.eq.0) then
                   oprnd=>tens_oprnd
                   call this%set_operand(0,oprnd,jerr)
                   if(jerr.ne.DSVP_SUCCESS) then
                    call this%terminate(TAVP_ERR_GEN_FAILURE,jerr)
                    jerr=-8
                   endif
                  else
                   call this%terminate(TAVP_ERR_GEN_FAILURE,jerr)
                   jerr=-7
                  endif
                 else
                  call this%terminate(TAVP_ERR_RSC_UNAVAILABLE,jerr)
                  jerr=-6
                 endif
                else
                 call this%terminate(TAVP_ERR_RSC_UNAVAILABLE,jerr)
                 jerr=-5
                endif
               else
                call this%terminate(TAVP_ERR_ARG_UNDEFINED,jerr)
                jerr=-4
               endif
              end select
             else
              call this%terminate(TAVP_ERR_CHE_FAILURE,jerr)
              jerr=-3
             endif
            else
             call this%terminate(TAVP_ERR_BTC_BAD,jerr)
             jerr=-2
            endif
           else
            call this%terminate(TAVP_ERR_RSC_UNAVAILABLE,jerr)
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
             call this%set_control(instr_ctrl,jerr)
             if(jerr.eq.DSVP_SUCCESS) then
              call this%alloc_operands(3,jerr)
              if(jerr.eq.DSVP_SUCCESS) then
               do jj=0,2
                tensor=>NULL(); allocate(tensor,STAT=jerr)
                if(jerr.ne.0) then; call this%terminate(TAVP_ERR_RSC_UNAVAILABLE,jerr); jerr=-11; exit; endif
                call tensor%tens_rcrsv_ctor(instr_packet,jerr)
                if(jerr.ne.TEREC_SUCCESS) then; call this%terminate(TAVP_ERR_BTC_BAD,jerr); jerr=-10; exit; endif
                tens_entry=>tavpWorker%tens_cache%lookup(tensor,jerr)
                if(jerr.ne.0) then; call this%terminate(TAVP_ERR_CHE_FAILURE,jerr); jerr=-9; exit; endif
                if(.not.associated(tens_entry)) then; call this%terminate(TAVP_ERR_ARG_UNDEFINED,jerr); jerr=-8; exit; endif
                deallocate(tensor); tensor=>NULL() !deallocate the temporary tensor and use
                tensor=>tens_entry%get_tensor()    !the same tensor from the tensor cache
                tens_resource=>NULL()              !as well as its associated resource
                tens_resource=>tens_entry%get_resource()
                tens_oprnd=>NULL(); allocate(tens_oprnd,STAT=jerr)
                if(jerr.ne.0) then; call this%terminate(TAVP_ERR_RSC_UNAVAILABLE,jerr); jerr=-7; exit; endif
                call tens_oprnd%tens_oprnd_ctor(tensor,jerr,tens_resource) !tensor and tens_resource are stored in tensor cache
                if(jerr.ne.0) then; call this%terminate(TAVP_ERR_GEN_FAILURE,jerr); jerr=-6; exit; endif
                oprnd=>tens_oprnd
                call this%set_operand(jj,oprnd,jerr)
                if(jerr.ne.DSVP_SUCCESS) then; call this%terminate(TAVP_ERR_GEN_FAILURE,jerr); jerr=-5; exit; endif
               enddo
              else
               call this%terminate(TAVP_ERR_RSC_UNAVAILABLE,jerr)
               jerr=-4
              endif
             else
              call this%terminate(TAVP_ERR_GEN_FAILURE,jerr)
              jerr=-3
             endif
            else
             call this%terminate(TAVP_ERR_BTC_BAD,jerr)
             jerr=-2
            endif
           else
            call this%terminate(TAVP_ERR_RSC_UNAVAILABLE,jerr)
            jerr=-1
           endif
           return
          end subroutine decode_instr_contract

        end subroutine TensInstrDecode
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
!--------------------------------------------------
        subroutine TensInstrSetMicrocode(this,ierr)
!Sets up instruction dynamic bindings to the corresponding microcode.
         implicit none
         class(tens_instr_t), intent(inout):: this   !inout: defined tensor instruction
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc,op_code

         op_code=this%get_code(errc)
         if(errc.eq.DSVP_SUCCESS) then
          if(op_code.ge.0.and.op_code.lt.TAVP_ISA_SIZE) then
           this%acquire_resource=>microcode(op_code)%acquire_resource
           this%prefetch_input=>microcode(op_code)%prefetch_input
           this%sync_prefetch=>microcode(op_code)%sync_prefetch
           this%execute=>microcode(op_code)%execute
           this%sync_execution=>microcode(op_code)%sync_execution
           this%upload_output=>microcode(op_code)%upload_output
           this%sync_upload=>microcode(op_code)%sync_upload
           this%release_resource=>microcode(op_code)%release_resource
          else
           errc=-2
          endif
         else
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensInstrSetMicrocode
!------------------------------------------------------
        subroutine TensInstrActivate(this,op_code,ierr)
!Activates the instruction after it has been constructed or decoded.
         implicit none
         class(tens_instr_t), intent(inout):: this   !inout: defined tensor instruction
         integer(INTD), intent(in):: op_code         !in: instruction code
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         call this%set_code(op_code,errc)
         if(errc.eq.DSVP_SUCCESS) then
          if(this%get_status().eq.DS_INSTR_EMPTY) then
           call this%set_microcode(errc)
           if(errc.eq.0) then
            call this%set_status(DS_INSTR_NEW,errc,DSVP_SUCCESS)
            if(errc.ne.DSVP_SUCCESS) then
             call this%set_status(DS_INSTR_RETIRED,errc,TAVP_ERR_GEN_FAILURE)
             errc=-3
            endif
           else
            call this%set_status(DS_INSTR_RETIRED,errc,TAVP_ERR_BTC_BAD)
            errc=-2
           endif
          endif
         else
          call this%set_status(DS_INSTR_RETIRED,errc,TAVP_ERR_GEN_FAILURE)
          errc=-1
         endif
         if(present(ierr)) ierr=errc
         return
        end subroutine TensInstrActivate
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
            call test_carma(errc) !`remove
            call this%shutdown(errc) !`debug: move
           else
            errc=-3
           endif
          else
           errc=-2
           call quit(errc,'FATAL(TAVP-Worker:start): Failed to initialize distributed address space!')
          endif
         else
          errc=-1
          call quit(errc,'FATAL(TAVP-Worker:start): Failed to initialize TAL-SH!')
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

          jerr=0
          call tavp_addr_space%create(role_comm,TAVP_WORKER_NUM_WINS,'WorkAddressSpace',jerr)
          return
         end subroutine init_distributed_space

         subroutine init_microcode(jerr)
          implicit none
          integer(INTD), intent(out):: jerr

          jerr=0
 !TENSOR CREATE:
          microcode(TAVP_INSTR_CREATE)%acquire_resource=>acquire_resource_basic
          microcode(TAVP_INSTR_CREATE)%prefetch_input=>prefetch_input_empty
          microcode(TAVP_INSTR_CREATE)%sync_prefetch=>sync_prefetch_empty
          microcode(TAVP_INSTR_CREATE)%execute=>execute_tensor_create
          microcode(TAVP_INSTR_CREATE)%sync_execution=>sync_execution_empty
          microcode(TAVP_INSTR_CREATE)%upload_output=>upload_output_empty
          microcode(TAVP_INSTR_CREATE)%sync_upload=>sync_upload_empty
          microcode(TAVP_INSTR_CREATE)%release_resource=>release_resource_empty
 !TENSOR DESTROY:
          microcode(TAVP_INSTR_DESTROY)%acquire_resource=>acquire_resource_empty
          microcode(TAVP_INSTR_DESTROY)%prefetch_input=>prefetch_input_empty
          microcode(TAVP_INSTR_DESTROY)%sync_prefetch=>sync_prefetch_empty
          microcode(TAVP_INSTR_DESTROY)%execute=>execute_tensor_destroy
          microcode(TAVP_INSTR_DESTROY)%sync_execution=>sync_execution_empty
          microcode(TAVP_INSTR_DESTROY)%upload_output=>upload_output_empty
          microcode(TAVP_INSTR_DESTROY)%sync_upload=>sync_upload_empty
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
          return
         end subroutine init_microcode

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
!----------------------------------------------------------------
        subroutine TAVPWorkerFetchInstructions(this,dsvp_id,ierr)
!Fetches a block of tensor instructions from TAVP "Manager".
         implicit none
         class(tavp_worker_t), intent(inout):: this  !inout: TAVP "Worker"
         integer(INTD), intent(in):: dsvp_id         !in: ID of the sender TAVP (or TAVP_ANY_ID wildcard)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWorkerFetchInstructions
!------------------------------------------------------------------------
        subroutine TAVPWorkerReturnRetiredInstructions(this,dsvp_id,ierr)
!Returns a block of retired tensor instructions back to TAVP "Manager".
         implicit none
         class(tavp_worker_t), intent(inout):: this  !inout: TAVP "Worker"
         integer(INTD), intent(in):: dsvp_id         !in: ID of the receiver TAVP
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWorkerReturnRetiredInstructions
!---------------------------------------------------------------
        subroutine TAVPWorkerSendInstructions(this,dsvp_id,ierr) !dummy
!Sends a block of tensor instructions to another TAVP.
         implicit none
         class(tavp_worker_t), intent(inout):: this  !inout: TAVP "Worker"
         integer(INTD), intent(in):: dsvp_id         !in: ID of the receiver TAVP
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWorkerSendInstructions
!-------------------------------------------------------------------------
        subroutine TAVPWorkerReceiveRetiredInstructions(this,dsvp_id,ierr) !dummy
!Receives a block of retired tensor instructions from another TAVP.
         implicit none
         class(tavp_worker_t), intent(inout):: this  !inout: TAVP "Worker"
         integer(INTD), intent(in):: dsvp_id         !in: ID of the sender TAVP (or TAVP_ANY_ID wildcard)
         integer(INTD), intent(out), optional:: ierr !out: error code
         integer(INTD):: errc

         errc=0
         if(present(ierr)) ierr=errc
         return
        end subroutine TAVPWorkerReceiveRetiredInstructions

       end module tavp_worker
