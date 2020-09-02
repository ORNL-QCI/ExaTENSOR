       program main
       use combinatoric
       use stsubs
       implicit none
       integer, parameter:: MAX_OPS=1000000
       real(8), parameter:: FMA_FAC=8d0

       type perf_entry_t
        integer:: id
        real(8):: time
        real(8):: flops_total
        real(8):: flops_process
        real(8):: proc_sec_est
       end type perf_entry_t

       integer:: i,j,k,n,l0,l1,l2,ip,np,trn(0:MAX_OPS)
       integer:: ao_id,ao_dim,ao_est,oc_id,oc_dim,oc_est,vi_id,vi_dim,vi_est
       character(1024):: str0,str1,str2
       real(8):: tms,tmf,flops,flops_est,subspace_dim,subspace_est,fl,prochour
       type(perf_entry_t):: table(MAX_OPS)

       open(10,file='perf_analysis.txt',form='FORMATTED',status='UNKNOWN')
       open(11,file='qforce.log',form='FORMATTED',status='OLD')
 !Read dimensions of spaces:
       open(12,file='space_dims.txt',form='FORMATTED',status='OLD')
       read(12,*) ao_id,ao_dim,ao_est !AO subspace id => AO subspace dimension
       read(12,*) oc_id,oc_dim,oc_est !OCC subspace id => OCC subspace dimension
       read(12,*) vi_id,vi_dim,vi_est !VIRT subspace id => VIRT subspace dimension
       close(12)
 !Get the total number of TAVP-WRK instances:
       read(11,*) np
       write(10,'("#TAVP-WRK processors = ",i7)') np
       n=0
       do
        str0=' '; read(11,'(A1024)',end=100) str0; l0=len_trim(str0)
        if(l0.gt.0) then
         if(index(str0(1:l0),'#MSG(exatensor): New Instruction: CONTRACT TENSORS').gt.0) then
          n=n+1
          str1=' '; read(11,'(A1024)') str1; l1=len_trim(str1)
          str2=' '; read(11,'(A1024)') str2; l2=len_trim(str2)
 !Get instruction id:
          i=index(str0(1:l0),'IP =')+4
          do while(str0(i:i).eq.' '); i=i+1; enddo
          call charnum(str0(i:l0),flops,ip)
 !Get beginning time stamp:
          i=index(str0(1:l0),'[')+1; j=index(str0(1:l0),']')-1
          do while(str0(i:i).eq.' '); i=i+1; enddo
          do while(str0(j:j).eq.' '); j=j-1; enddo
          call charnum(str0(i:j),tms,k)
 !Get ending time stamp:
          i=index(str2(1:l2),'[')+1; j=index(str2(1:l2),']')-1
          do while(str2(i:i).eq.' '); i=i+1; enddo
          do while(str2(j:j).eq.' '); j=j-1; enddo
          call charnum(str2(i:j),tmf,k)
 !Get flop count:
          flops=1d0; flops_est=1d0
          i=index(str1(1:l1),':')+1
          do
           j=index(str1(i:l1),':')+i-1; if(j.lt.i) exit
           i=j+1; do while(is_it_number(str1(i:i)).ge.0); i=i+1; enddo
           call charnum(str1(j+1:i-1),subspace_dim,k)
           if(k.eq.ao_id) then
            subspace_dim=dble(ao_dim); subspace_est=dble(ao_est)
           elseif(k.eq.oc_id) then
            subspace_dim=dble(oc_dim); subspace_est=dble(oc_est)
           elseif(k.eq.vi_id) then
            subspace_dim=dble(vi_dim); subspace_est=dble(vi_est)
           else
            write(*,'("#ERROR: Unknown subspace id: ",i13)') k
            stop
           endif
           flops=flops*subspace_dim !flops performed
           flops_est=flops_est*subspace_est !flops to be performed
           i=j+1
          enddo
          flops=dsqrt(flops)*FMA_FAC
          flops_est=dsqrt(flops_est)*FMA_FAC
 !Save the result:
          if(tmf-tms.gt.0d0) then
           fl=flops/(tmf-tms) !measured flop/s
           table(n)=perf_entry_t(ip,tmf-tms,fl/1d12,fl/1d9/dble(np),flops_est/(fl/dble(np)))
          else
           table(n)=perf_entry_t(ip,tmf-tms,0d0,0d0,0d0)
          endif
         endif
        endif
       enddo
100    close(11)
 !Sort the table and print the results:
       prochour=0d0
       trn(0:n)=(/+1,(j,j=1,n)/)
       call merge_sort_key(n,table(1:n)%time,trn)
       do i=1,n
        j=trn(i)
        write(10,'(i7,": ",F13.4," sec: ",F10.3," Tflop/s: ",F10.3," Gflop/s: Process-Hour est = ",F10.1)')&
        &table(j)%id,table(j)%time,table(j)%flops_total,table(j)%flops_process,table(j)%proc_sec_est/3600d0
        prochour=prochour+table(j)%proc_sec_est/3600d0
       enddo
       write(10,'("Total process-hour estimate = ",F10.1)') prochour
       close(10)
       end program main
