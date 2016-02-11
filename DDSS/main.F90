        module aux
        contains

         function array_norm(arr,arr_vol) result(norm1)
         use service_mpi, only: INT_COUNT
         implicit none
         real(8):: norm1
         real(8), intent(in):: arr(1:*)
         integer(INT_COUNT), intent(in):: arr_vol
         integer(INT_COUNT):: i
         norm1=0d0
!$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i) REDUCTION(+:norm1) SCHEDULE(GUIDED)
         do i=1,arr_vol
          norm1=norm1+abs(arr(i))
         enddo
!$OMP END PARALLEL DO
         return
         end function array_norm

        end module aux
!-------------------------------------------------------------
        program main
        use service_mpi
        use distributed
        use dil_basic
        use timers
        use extern_names, only: c_ptr_value
        use aux
        implicit none

        integer(INT_COUNT), parameter:: MAX_BUF_VOL=40000000   !max volume of the send buffer (for testing)
        integer(INT_MPI), parameter:: NUM_PASSES=8             !number of passes with a different pause length
        real(8), parameter:: PAUS_INC=0.05d0                   !pause length increment in seconds
        integer(INT_MPI), parameter:: NUM_WINS_PER_SPACE=1     !number of MPI windows per distributed space
        integer(INT_MPI), parameter:: MAX_PACK_LEN=1024        !max packet length (internal use)

        real(8), allocatable, target:: send_buf(:),recv_buf(:)
        integer(INT_COUNT):: buf_vol0,buf_vol1,pack_len0,pack_len1
        type(C_PTR):: cptr
        integer(INT_MPI):: i,n,ierr
        type(DistrSpace_t):: dspace0
        type(DataDescr_t):: descr0,descr1
        type(PackCont_t):: dpack0,dpack1
        type(CommHandle_t):: ch0,ch1
        real(8):: paus,rnd,tms,tm,snorm1,snorm2
        integer(ELEM_PACK_SIZE):: packet0(MAX_PACK_LEN),packet1(MAX_PACK_LEN)
        integer:: errc,tmr

!Initialize the MPI infrastructure:
        call dil_process_start(ierr)
        write(jo,*) 'MPI process started(rank,comm,err): ',impir,GLOBAL_MPI_COMM,ierr
        if(ierr.ne.0) call quit(ierr,'ERROR: Failed to start a process!')
!Create a distributed memory space over the Global MPI communicator:
        call dspace0%create(GLOBAL_MPI_COMM,NUM_WINS_PER_SPACE,'My Space',ierr)
        write(jo,*) 'My space created(rank,err): ',impir,ierr
        if(ierr.ne.0) call quit(ierr,'ERROR: Failed to create my space!')
        flush(jo) !jo: output device; impir: current MPI rank.
!Allocate a buffer on each MPI process:
        do i=0,impir; call random_number(rnd); enddo !to make rnd different on different MPI processes
!       buf_vol0=int(dble(MAX_BUF_VOL)*rnd,INT_COUNT)+4 !random send buffer volume
        buf_vol0=MAX_BUF_VOL !let's have them all of the sam size for now
        allocate(send_buf(1:buf_vol0)); cptr=c_loc(send_buf)
!$OMP WORKSHARE
        send_buf(1:buf_vol0)=0d0
!$OMP END WORKSHARE
        write(jo,*) 'Allocated a send buffer(rank,vol,addr): ',impir,buf_vol0,c_ptr_value(cptr)
!Attach the buffer to the distributed memory space as double precision (R8):
        call dspace0%attach(c_loc(send_buf),R8,buf_vol0,descr0,ierr) !send buffer is associated with the data descriptor <descr0>
        write(jo,*) 'Attached a buffer(rank,vol,err): ',impir,buf_vol0,ierr
        if(ierr.ne.0) call quit(ierr,'ERROR: Failed to attach the send buffer!')
!DEBUG: Print the original data descriptor:
        call descr0%print_it(dev_out=jo); flush(jo)
!Pack the data descriptor:
        call descr0%pack(packet0,ierr,pack_len0)
        write(jo,*) 'Packed the original descriptor(rank,len,err): ',impir,pack_len0,ierr
        if(ierr.ne.0) call quit(ierr,'ERROR: Failed to pack a data descriptor!')
!DEBUG: Unpack the packet into another descriptor:
!        pack_len0=0; call descr1%unpack(packet0,ierr,pack_len0)
!        write(jo,*) 'UnPacked into another descriptor(rank,len,err): ',impir,pack_len0,ierr
!        if(ierr.ne.0) call quit(ierr,'ERROR: Failed to unpack a data descriptor!')
!DEBUG: Print the resulting data descriptor:
!        call descr1%print_it(dev_out=jo); flush(jo)
!Exchange data descriptors between processes:
 !Create the data packet container:
        call dpack0%reserve_mem(MAX_PACK_LEN,ierr)
        write(jo,*) 'Created a data packet container 0 (rank,err): ',impir,ierr
        if(ierr.ne.0) call quit(ierr,'ERROR: Failed to create a data packet container!')
        call dpack1%reserve_mem(MAX_PACK_LEN,ierr)
        write(jo,*) 'Created a data packet container 1 (rank,err): ',impir,ierr
        if(ierr.ne.0) call quit(ierr,'ERROR: Failed to create a data packet container!')
 !Append packet 0 to the data packet container:
        call dpack0%append(packet0,ierr)
        write(jo,*) 'Appended a data packet to the data packet container 0 (rank,err): ',impir,ierr
        if(ierr.ne.0) call quit(ierr,'ERROR: Failed to append a data packet to a data packet container!')
 !Send the data packet container to the next MPI process:
        i=mod(impir+1,impis) !next MPI process
        call dpack0%send(ch0,ierr,i) !ch0: communication request handle
        write(jo,*) 'Send initiated for the data packet container (rank_from,rank_to,err): ',impir,i,ierr
        if(ierr.ne.0) call quit(ierr,'ERROR: Failed to initiate a send operation!')
 !Receive a data packet container from the previous MPI process:
        i=mod(impis+impir-1,impis) !previous MPI process
        call dpack1%receive(ch1,ierr,i) !ch1: communication request handle
        write(jo,*) 'Receive initiated for the data packet container (rank_from,rank_to,err): ',i,impir,ierr
        if(ierr.ne.0) call quit(ierr,'ERROR: Failed to initiate a receive operation!')
 !Synchronize the receive:
        call ch1%wait(ierr,ignore_old=.true.)
        write(jo,*) 'Synchronized the receive of the data packet container (rank,err): ',impir,ierr
        if(ierr.ne.0) call quit(ierr,'ERROR: Failed to synchronize the receive operation!')
        write(jo,*) 'Total number of packets in the data packet container = ',dpack1%num_packets()
 !Extract a packet from the data packet container:
        call dpack1%remove(ierr,packet=packet1)
        write(jo,*) 'Extracted a packet from the received data packet container (rank,err): ',impir,ierr
        if(ierr.ne.0) call quit(ierr,'ERROR: Failed to extract a packet from the data packet container!')
 !Get the data descriptor from the packet:
        call descr1%unpack(packet1,ierr,pack_len1)
        buf_vol1=descr1%data_volume()
        write(jo,*) 'Unpacked a data descriptor from the extracted data packet (rank,len,data_vol,err): ',&
        &impir,pack_len1,buf_vol1,ierr
        if(ierr.ne.0.or.buf_vol1.le.0) call quit(ierr,'ERROR: Failed to unpack a data descriptor!')
        call descr1%print_it(dev_out=jo); flush(jo) !print it: DEBUG
!Allocate a receive buffer of appropriate volume:
        allocate(recv_buf(1:buf_vol1))
        write(jo,*) 'Allocated a receive buffer(rank,vol): ',impir,buf_vol1

        paus=0d0
        do n=1,NUM_PASSES

!Fill the send buffer with random numbers:
        call random_number(send_buf); snorm1=array_norm(send_buf,buf_vol0)
        write(jo,*) 'Norm1 of the send buffer = ',snorm1,'; Volume = ',buf_vol0
!Clear receive buffer:
!$OMP WORKSHARE
         recv_buf(1:buf_vol1)=0d0
!$OMP END WORKSHARE
!Sync for clear timing:
         call dil_global_comm_barrier()

!Fetch the remote data into the receive buffer:
 !Initiate fetching:
         tms=thread_wtime()
         call descr1%get_data(c_loc(recv_buf),ierr,MPI_ASYNC_NRM)
         tm=thread_wtime(tms)
         write(jo,*) 'Initiated fetching remote data (rank,time,err): ',impir,tm,ierr
         call ddss_print_stat() !DEBUG
 !Pause (like we are doing some computations now and the MPI message is progressing on the background):
         tms=thread_wtime()
         errc=timer_start(paus,tmr); do while(.not.time_is_off(tmr,errc)); enddo !DEBUG: Pause before FLUSH
         tm=thread_wtime(tms); write(jo,*) 'Pause before FLUSH (sec) = ',tm !DEBUG
 !Now complete fetching (hopefully the MPI message is already here):
         tms=thread_wtime()
         call descr1%flush_data(ierr)
         tm=thread_wtime(tms)
         write(jo,*) 'Completed the fetch after a pause: (rank,time,err): ',impir,tm,ierr
         write(jo,*) 'Fetch bandwidth (GB/s) = ',dble(buf_vol1*8)/(tm*1024d0*1024d0*1024d0)
         call ddss_print_stat() !DEBUG
         write(jo,*) 'Norm1 of the receive buffer = ',array_norm(recv_buf,buf_vol1),'; Volume = ',buf_vol1
!Invert the sign of the receive buffer:
!$OMP WORKSHARE
         recv_buf(1:buf_vol1)=-recv_buf(1:buf_vol1)
!$OMP END WORKSHARE
!Sync for clear timing:
         call dil_global_comm_barrier()

!Accumulate the fetched data back to the target process:
 !Initiate the accumulate:
         tms=thread_wtime()
         call descr1%acc_data(c_loc(recv_buf),ierr,MPI_ASYNC_NRM)
         tm=thread_wtime(tms)
         write(jo,*) 'Initiated a remote accumulate (rank,time,err): ',impir,tm,ierr
         call ddss_print_stat() !DEBUG
 !Pause (like we are doing some computations now and the MPI message is progressing on the background):
         tms=thread_wtime()
         errc=timer_start(paus,tmr); do while(.not.time_is_off(tmr,errc)); enddo !DEBUG: Pause before FLUSH
         tm=thread_wtime(tms); write(jo,*) 'Pause before FLUSH (sec) = ',tm !DEBUG
 !Now complete the accumulate (hopefully it is already completed):
         tms=thread_wtime()
         call descr1%flush_data(ierr)
         tm=thread_wtime(tms)
         write(jo,*) 'Completed the accumulate after a pause: (rank,time,err): ',impir,tm,ierr
         write(jo,*) 'Accumulate bandwidth (GB/s) = ',dble(buf_vol1*8)/(tm*1024d0*1024d0*1024d0)
         call ddss_print_stat() !DEBUG
         call dil_global_comm_barrier()
!Check the resulting array norm:
         snorm2=array_norm(send_buf,buf_vol0)
         write(jo,*) 'Norm1 of the accumulated send buffer = ',snorm2,'; Volume = ',buf_vol0
         flush(jo)

         paus=paus+PAUS_INC
        enddo

        write(jo,*) 'Number of passes done = ',NUM_PASSES
        flush(jo)

!Destroy the data packet container:
        call dpack1%clean(ierr)
        write(jo,*) 'Destroyed the data packet container (rank,ierr): ',impir,ierr
        if(ierr.ne.0) call quit(ierr,'ERROR: Failed to destroy a data packet container!')
        call dpack0%clean(ierr)
        write(jo,*) 'Destroyed the data packet container (rank,ierr): ',impir,ierr
        if(ierr.ne.0) call quit(ierr,'ERROR: Failed to destroy a data packet container!')
!Detach buffers from the distributed memory space:
        call dspace0%detach(descr0,ierr)
        write(jo,*) 'Detached the send buffer(rank,err): ',impir,ierr
        if(ierr.ne.0) call quit(ierr,'ERROR: Failed to detach the send buffer!')
!Deallocate memory buffers:
        deallocate(recv_buf)
        deallocate(send_buf)
!Destroy distributed memory spaces:
        call dspace0%destroy(ierr); write(jo,*) 'My space destroyed(rank,err): ',impir,ierr
        if(ierr.ne.0) call quit(ierr,'ERROR: Failed to destroy my space!')
!Finalize MPI:
        call dil_process_finish(ierr); write(*,*) 'MPI process finished(err): ',ierr
        if(ierr.ne.0) call quit(ierr,'ERROR: Failed to finalize the process!')
        stop
        end program main
