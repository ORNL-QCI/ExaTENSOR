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
       end type perf_entry_t

       integer:: i,j,k,n,l0,l1,l2,ip,np,trn(0:MAX_OPS)
       character(1024):: str0,str1,str2
       real(8):: tms,tmf,flops,subspace
       type(perf_entry_t):: table(MAX_OPS)

       open(10,file='perf_analysis.txt',form='FORMATTED',status='UNKNOWN')
       open(11,file='qforce.log',form='FORMATTED',status='OLD')
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
          flops=1d0
          i=index(str1(1:l1),':')+1
          do
           j=index(str1(i:l1),':')+i-1; if(j.lt.i) exit
           i=j+1; do while(is_it_number(str1(i:i)).ge.0); i=i+1; enddo
           call charnum(str1(j+1:i-1),subspace,k)
           flops=flops*subspace
           i=j+1
          enddo
          flops=dsqrt(flops)*FMA_FAC
 !Save the result:
          if(tmf-tms.gt.0d0) then
           table(n)=perf_entry_t(ip,tmf-tms,flops/1d12/(tmf-tms),flops/1d9/(tmf-tms)/dble(np))
          else
           table(n)=perf_entry_t(ip,tmf-tms,0d0,0d0)
          endif
         endif
        endif
       enddo
100    close(11)
 !Sort the table and print the results:
       trn(0:n)=(/+1,(j,j=1,n)/)
       call merge_sort_key(n,table(1:n)%time,trn)
       do i=1,n
        j=trn(i)
        write(10,'(i7,": ",F13.4," sec: ",F10.3," Tflop/s: ",F10.3," Gflop/s")') table(j)%id,table(j)%time,&
        &table(j)%flops_total,table(j)%flops_process
       enddo
       close(10)
       end program main
