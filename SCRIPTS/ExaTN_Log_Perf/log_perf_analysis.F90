       program main
       use parse_prim
       use combinatoric
       use stsubs
       implicit none
       real(8):: time_start=-1d0  !tracing visualization start time in sec
       real(8):: time_finish=-1d0 !tracing visualization finish time in sec
       !-------------------------
       character(1024):: str
       integer(8):: cnt,perc
       integer:: pred_offset(64),pred_length(64),num_pred,ts,id,opcode,l,n,ierr
       logical:: matched
       real(8):: val,mem_used,mem_free,flop,last_flop,last_time

       open(10,file='exatn_process.json',form='FORMATTED',status='UNKNOWN')
       write(10,'("{ ""traceEvents"": [")')
       open(11,file='exatn_process.txt',form='FORMATTED',status='OLD') !exatn_exec_thread.X.log
       id=-1; opcode=-1; last_flop=-1d0; last_time=-1d0
       do
        str=' '; read(11,'(A1024)',end=100) str; l=len_trim(str)
        if(l.gt.0) then
         matched=match_symb_pattern(str(1:l),"[`](`)[`]: Submitting tensor operation `: Opcode = `: Details:",&
                                   &num_pred,pred_offset,pred_length,ierr)
         if(ierr.ne.0) then
          write(*,'("ERROR: Unable to parse line: ")',ADVANCE='NO'); write(*,*) str(1:l)
          stop
         endif
         if(matched) then !try matching tensor operation submission
          if(num_pred.ne.5) then
           write(*,'("ERROR: Invalid parsing of line: ")',ADVANCE='NO'); write(*,*) str(1:l)
           stop
          endif
          call charnum(str(pred_offset(4):pred_offset(4)+pred_length(4)-1),val,id)
          call charnum(str(pred_offset(5):pred_offset(5)+pred_length(5)-1),val,opcode)
          call charnum(str(pred_offset(1):pred_offset(1)+pred_length(1)-1),val,n)
          if((time_start.lt.0d0.or.(val.ge.time_start)).and.(time_finish.lt.0d0.or.(val.le.time_finish))) then
           ts=anint(val*1d6) !integer time stamp in microseconds
           write(10,'("{ ""cat"": ""exatn"", ""name"": """,i10,""", ""ph"": ""B"", ""pid"": 0, ""tid"": 0, ""ts"": ",i18," },")')&
           &id,ts
          endif
         else !try matching tensor operation synchronization
          matched=match_symb_pattern(str(1:l),": Status = `: Syncing ... Success [`]",&
                                    &num_pred,pred_offset,pred_length,ierr)
          if(ierr.ne.0) then
           write(*,'("ERROR: Unable to parse line: ")',ADVANCE='NO'); write(*,*) str(1:l)
           stop
          endif
          if(matched) then !immediate synchronization
           if(num_pred.ne.2) then
            write(*,'("ERROR: Invalid parsing of line: ")',ADVANCE='NO'); write(*,*) str(1:l)
            stop
           endif
           if(id.ge.0) then
            call charnum(str(pred_offset(2):pred_offset(2)+pred_length(2)-1),val,n)
            if((time_start.lt.0d0.or.(val.ge.time_start)).and.(time_finish.lt.0d0.or.(val.le.time_finish))) then
             ts=anint(val*1d6) !integer time stamp in microseconds
             write(10,'("{ ""cat"": ""exatn"", ""name"": """,i10,""", ""ph"": ""E"", ""pid"": 0, ""tid"": 0, ""ts"": ",i18," },")')&
             &id,ts
            endif
            id=-1; opcode=-1
           else
            write(*,'("ERROR: Invalid input log file!")')
            stop
           endif
          else !deferred synchronization
           matched=match_symb_pattern(str(1:l),"[`](`)[`]: Synced tensor operation `: Opcode = `",&
                                     &num_pred,pred_offset,pred_length,ierr)
           if(ierr.ne.0) then
            write(*,'("ERROR: Unable to parse line: ")',ADVANCE='NO'); write(*,*) str(1:l)
            stop
           endif
           if(matched) then
            if(num_pred.ne.5) then
             write(*,'("ERROR: Invalid parsing of line: ")',ADVANCE='NO'); write(*,*) str(1:l)
             stop
            endif
            call charnum(str(pred_offset(4):pred_offset(4)+pred_length(4)-1),val,id)
            call charnum(str(pred_offset(5):pred_offset(5)+pred_length(5)-1),val,opcode)
            call charnum(str(pred_offset(1):pred_offset(1)+pred_length(1)-1),val,n)
            if((time_start.lt.0d0.or.(val.ge.time_start)).and.(time_finish.lt.0d0.or.(val.le.time_finish))) then
             ts=anint(val*1d6) !integer time stamp in microseconds
             write(10,'("{ ""cat"": ""exatn"", ""name"": """,i10,""", ""ph"": ""E"", ""pid"": 0, ""tid"": 0, ""ts"": ",i18," },")')&
             &id,ts
            endif
            id=-1; opcode=-1
           else !try matching total flop counter
            matched=match_symb_pattern(str(1:l),"[`] Total Flop count = `; Memory usage = `, Free = `",&
                                      &num_pred,pred_offset,pred_length,ierr)
            if(ierr.ne.0) then
             write(*,'("ERROR: Unable to parse line: ")',ADVANCE='NO'); write(*,*) str(1:l)
             stop
            endif
            if(matched) then
             if(num_pred.ne.4) then
              write(*,'("ERROR: Invalid parsing of line: ")',ADVANCE='NO'); write(*,*) str(1:l)
              stop
             endif
             call charnum(str(pred_offset(4):pred_offset(4)+pred_length(4)-1),mem_free,n)
             call charnum(str(pred_offset(3):pred_offset(3)+pred_length(3)-1),mem_used,n)
             call charnum(str(pred_offset(2):pred_offset(2)+pred_length(2)-1),flop,n)
             call charnum(str(pred_offset(1):pred_offset(1)+pred_length(1)-1),val,n)
             if((time_start.lt.0d0.or.(val.ge.time_start)).and.(time_finish.lt.0d0.or.(val.le.time_finish))) then
              if(last_time.ge.0d0) then
               if(val-last_time.ge.1d0) then
                cnt=anint((flop-last_flop)/(val-last_time)/1d9) !GFlop/s
                perc=anint(mem_used*1d2/(mem_used+mem_free)) !percent of memory in use
                ts=anint(val*1d6) !integer time stamp in microseconds
                write(10,'("{ ""cat"": ""exatn"", ""name"": ""performance"", ""ph"": ""C"", ""pid"": 0, ""tid"": 0, ""ts"": ",i9,'&
                &//'", ""args"": {""gflop/s"": ",i18,"} },")') ts,cnt
                write(10,'("{ ""cat"": ""exatn"", ""name"": ""memory"", ""ph"": ""C"", ""pid"": 0, ""tid"": 0, ""ts"": ",i9,'&
                &//'", ""args"": {""percent"": ",i18,"} },")') ts,perc
                last_time=val
                last_flop=flop
               endif
              else
               last_time=val
               last_flop=flop
              endif
             endif
            endif
           endif
          endif
         endif
        endif
       enddo
100    close(11)
       write(10,'("{ ""name"": ""end"", ""ph"": ""I"", ""pid"": 0, ""tid"": 0, ""s"": ""p"", ""ts"": ",i18," }")') (ts+1)
       write(10,'("] }")')
       close(10)
       write(*,'("JSON output generated successfully in exatn_process.json")')
       end program main
