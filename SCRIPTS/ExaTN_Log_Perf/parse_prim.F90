!Useful parsing primitives.
!AUTHOR: Dmitry I. Lyakh (Liakh): quant4me@gmail.com
!REVISION: 2017-02-22
!LICENSE: BSD 3-Clause
       module parse_prim
        use stsubs
        implicit none
        private
!PARAMETERS:
 !Limits:
        integer, parameter, public:: MAX_PARSE_RECUR_LEN=8192
        integer, parameter, public:: MAX_TENSOR_RANK=32
        integer, parameter, public:: MAX_TENSOR_OPERANDS=4
        integer, parameter, public:: MAX_TENSOR_NAME_LEN=512
 !Arguments:
        integer, parameter, public:: ANY_CASE=0     !any case wildcard
        integer, parameter, public:: SMALL_ONLY=1   !small letter filter
        integer, parameter, public:: CAPITAL_ONLY=2 !capital letter filter
        integer, parameter, public:: INDEX_SO=3     !spin-orbital index label filter
        integer, parameter, public:: INDEX_MO=4     !molecular-orbital index label filter
 !Tensor operations:
        integer, parameter, public:: COMMENT_LINE=-1
        integer, parameter, public:: NO_OPERATION=0
        integer, parameter, public:: TENSOR_ASSIGNMENT=1
        integer, parameter, public:: TENSOR_SCALING=2
        integer, parameter, public:: TENSOR_ADDITION=3
        integer, parameter, public:: TENSOR_CONTRACTION=4
 !Auxiliary operations:
        integer, parameter, public:: TENSOR_CREATE=5
        integer, parameter, public:: TENSOR_DELETE=6
        integer, parameter, public:: BARRIER_SYNC=7
 !Internal:
        integer, parameter, private:: ASCII_TAB=9   !ASCII code for TAB
!TYPES:
        type, public:: tensor_info_t
         integer:: tens_rank=-1
         character(MAX_TENSOR_RANK):: tens_shape
         integer:: tens_name_len=0
         character(MAX_TENSOR_NAME_LEN):: tens_name
         real(8):: prefactor
        end type tensor_info_t
!VISIBILITY:
        public is_this_blank
        public skip_blanks
        public remove_blanks
        public is_this_one_of
        public is_this_integer
        public is_this_real_number
        public are_these_letters
        public are_these_alphanumeric
        public begins_with
        public match_symb_pattern
        public is_this_index_label
        public match_index_label
        public is_this_tensor
        public is_this_tensor_operation
        public is_this_auxiliary_operation
        public remove_comments
        public parse_tadl_statement
        public print_tensor_info
        private translate_index_labels

       contains
!-----------------------------------------
        logical function is_this_blank(ch)
!Recognizes a blank (either a space or a TAB).
         implicit none
         character(1), intent(in):: ch !in: character

         if(ch.eq.' '.or.iachar(ch).eq.ASCII_TAB) then !{space,tab}
          is_this_blank=.TRUE.
         else
          is_this_blank=.FALSE.
         endif
         return
        end function is_this_blank
!--------------------------------------------------
        subroutine skip_blanks(str,beg_pos,end_pos)
!Returns the start and the end of the non-blank part of the line.
         implicit none
         character(*), intent(in):: str !in: string
         integer, intent(out):: beg_pos !out: first non-blank position
         integer, intent(out):: end_pos !out: last non-blank position
         integer:: i

         end_pos=len_trim(str); beg_pos=end_pos+1
         if(end_pos.gt.0) then
          do i=1,end_pos
           if(.not.is_this_blank(str(i:i))) then
            beg_pos=i; exit
           endif
          enddo
         endif
         return
        end subroutine skip_blanks
!---------------------------------------------
        subroutine remove_blanks(stri,stro,lo)
!Returns the string after removing all blanks.
         implicit none
         character(*), intent(in):: stri  !in: input string
         character(*), intent(out):: stro !out: output string
         integer, intent(out):: lo        !out: length of the output string
         integer:: i,j,li

         li=len_trim(stri); lo=0
         if(li.gt.0) then
          j=li+1
          do i=1,li
           if(.not.is_this_blank(stri(i:i))) then; j=i; exit; endif
          enddo
          if(j.le.li) then; lo=li-j+1; stro(1:lo)=stri(j:li); endif
         endif
         return
        end subroutine remove_blanks
!----------------------------------------------------
        logical function is_this_one_of(ch,char_list)
!Checks whether a character is one from the given list.
         implicit none
         character(1), intent(in):: ch           !in: input character
         character(1), intent(in):: char_list(:) !in: array of characters to test against
         integer:: i,l

         is_this_one_of=.FALSE.; l=size(char_list)
         do i=1,l
          if(ch.eq.char_list(i)) then; is_this_one_of=.TRUE.; exit; endif
         enddo
         return
        end function is_this_one_of
!----------------------------------------------------
        logical function is_this_integer(str,no_sign)
!Recognizes an integer (signed or unsigned).
         implicit none
         character(*), intent(in):: str          !in: string
         logical, intent(in), optional:: no_sign !in: TRUE does not allow for a sign (defaults to FALSE)
         integer:: i,j,beg_pos,end_pos
         logical:: nsign

         is_this_integer=.FALSE.; nsign=.FALSE.
         if(present(no_sign)) nsign=no_sign
         call skip_blanks(str,beg_pos,end_pos)
         if(beg_pos.le.end_pos) then
          j=beg_pos; if((.not.nsign).and.(str(j:j).eq.'-'.or.str(j:j).eq.'+')) j=j+1 !account for the sign
          if(j.le.end_pos) then
           do i=j,end_pos
            if(iachar(str(i:i)).lt.iachar('0').or.iachar(str(i:i)).gt.iachar('9')) return
           enddo
           is_this_integer=.TRUE.
          endif
         endif
         return
        end function is_this_integer
!----------------------------------------------------
        logical function is_this_real_number(str,val)
!Recognizes a real number in Fortran format.
         implicit none
         character(*), intent(in):: str !in: string
         real(8), intent(out):: val     !out: real value
         integer:: beg_pos,end_pos,n

         is_this_real_number=.FALSE.
         call skip_blanks(str,beg_pos,end_pos)
         if(beg_pos.le.end_pos) then
          call charnum(str(beg_pos:end_pos),val,n)
          if(dabs(val-dble(n)).gt.1.1d0) return
          is_this_real_number=.TRUE.
         endif
         return
        end function is_this_real_number
!----------------------------------------------------------
        logical function are_these_letters(str,letter_case)
!Recognizes whether the string consists of letters only (and blanks).
         implicit none
         character(*), intent(in):: str              !in: string
         integer, intent(in), optional:: letter_case !in: {ANY_CASE,SMALL_ONLY,CAPITAL_ONLY}
         integer:: i,lcase,beg_pos,end_pos

         are_these_letters=.FALSE.; lcase=ANY_CASE
         if(present(letter_case)) lcase=letter_case
         call skip_blanks(str,beg_pos,end_pos)
         if(beg_pos.le.end_pos) then
          select case(lcase)
          case(SMALL_ONLY)
           do i=beg_pos,end_pos
            if(iachar(str(i:i)).lt.iachar('a').or.iachar(str(i:i)).gt.iachar('z')) return
           enddo
           are_these_letters=.TRUE.
          case(CAPITAL_ONLY)
           do i=beg_pos,end_pos
            if(iachar(str(i:i)).lt.iachar('A').or.iachar(str(i:i)).gt.iachar('Z')) return
           enddo
           are_these_letters=.TRUE.
          case default !any case
           do i=beg_pos,end_pos
            if(.not.((iachar(str(i:i)).ge.iachar('a').and.iachar(str(i:i)).le.iachar('z')).or.&
                    &(iachar(str(i:i)).ge.iachar('A').and.iachar(str(i:i)).le.iachar('Z')))) return
           enddo
           are_these_letters=.TRUE.
          end select
         endif
         return
        end function are_these_letters
!---------------------------------------------------------------
        logical function are_these_alphanumeric(str,letter_case)
!Checks whether the string is alphanumeric: Letters, decimal digits, underscores.
         implicit none
         character(*), intent(in):: str              !in: string
         integer, intent(in), optional:: letter_case !in: letter case: {ANY_CASE,SMALL_ONLY,CAPITAL_ONLY}
         integer:: i,beg_pos,end_pos,lcase

         are_these_alphanumeric=.FALSE.; lcase=ANY_CASE
         if(present(letter_case)) lcase=letter_case
         call skip_blanks(str,beg_pos,end_pos)
         if(beg_pos.le.end_pos) then
          do i=beg_pos,end_pos
           if(.not.(is_this_integer(str(i:i),no_sign=.TRUE.).or.&
                   &are_these_letters(str(i:i),letter_case=lcase).or.&
                   &str(i:i).eq.'_')) return
          enddo
          are_these_alphanumeric=.TRUE.
         endif
         return
        end function are_these_alphanumeric
!---------------------------------------------------------------
        logical function begins_with(str,comp_stmt,end_pos,ierr)
!Matches a specific statement consisting of multiple words delimited by spaces,
!starting from the first non-blank position in <str>. That is, it is trying
!to see if <str> looks like this:
! [blanks] <comp_stmt> [other words] [blanks]
!where <comp_stmt> looks like this:
! [blanks] word [word] [word] ... [blanks]
         implicit none
         character(*), intent(in):: str        !in: string
         character(*), intent(in):: comp_stmt  !in: tested composite statement (multiple words)
         integer, intent(out):: end_pos        !out: in case of a match, the last position in <str> (the end of composite statement)
         integer, intent(out), optional:: ierr !out: error code
         integer:: b0,e0,b1,e1,i,j,errc

         begins_with=.FALSE.; end_pos=-1; errc=0
         call skip_blanks(str,b0,e0)
         call skip_blanks(comp_stmt,b1,e1)
         if(b0.le.e0.and.b1.le.e1) then
          i=b0; j=b1
          do while (i.le.e0)
           if(is_this_blank(str(i:i))) then
            if(is_this_blank(comp_stmt(j:j))) then
             i=i+1; do while(is_this_blank(str(i:i))); i=i+1; enddo
             j=j+1; do while(is_this_blank(comp_stmt(j:j))); j=j+1; enddo
            else
             exit
            endif
           else
            if(str(i:i).eq.comp_stmt(j:j)) then
             i=i+1; j=j+1
             if(j.gt.e1) then; begins_with=.TRUE.; exit; endif
            else
             exit
            endif
           endif
          enddo
         else
          errc=1
         endif
         if(present(ierr)) ierr=errc
         return
        end function begins_with
!---------------------------------------------------------------------------------------------------
        function match_symb_pattern(str,pattern,num_pred,pred_offset,pred_length,ierr) result(match)
!Matches a symbolic pattern. ` is a special placeholder symbol for predicates (can be empty).
!The matching is ambiguous in general (the same string may be matched in multiple ways with the same pattern).
!Always the first matching alternative will be returned. The proper matching may then be achieved by more
!elaborated patterns.
!Example:
! Pattern: multiply `(`)`
! Match 1: multiply T12(a1,d1)*S23(a2)*2.0
!          `1 = T12
!          `2 = a1,d1
!          `3 = *S23(a2)*2.0
! Match 2: multiply T12(a1,d1)*S23(a2)*2.0
!          `1 = T12
!          `2 = a1,d1)*S23(a2
!          `3 = *2.0
! Better pattern: multiply `(`)*`(`)`
! Unambiguous:
!          `1 = T12
!          `2 = a1,d1
!          `3 = S23
!          `4 = a2
!          `5 = *2.0
! If we remove "*2.0", that is, "multiply T12(a1,d1)*S23(a2)",
! then `5 = '' (empty).
         implicit none
         logical:: match                          !out: TRUE if matched, FALSE otherwise
         character(*), intent(in):: str           !in: input string
         character(*), intent(in):: pattern       !in: symbolic pattern being matched
         integer, intent(out):: num_pred          !out: number of predicates (each substring of <str> corresponding to each ` in the <pattern>)
         integer, intent(inout):: pred_offset(1:) !out: predicate offsets in str
         integer, intent(inout):: pred_length(1:) !out: predicate lengths in str
         integer, intent(out), optional:: ierr    !out: error code
         integer:: errc,bs,es,bp,ep,sp,pp
         integer:: pred,stp,pos_stk(2,1:MAX_PARSE_RECUR_LEN)
         logical:: over

         errc=0; match=.FALSE.; num_pred=0
         call skip_blanks(str,bs,es)
         call skip_blanks(pattern,bp,ep)
         if(ep.ge.bp.and.es.ge.bs) then
          pred=0; stp=0; sp=bs; pp=bp
          do while(pp.le.ep)
 !Are we at the new predicate boundary?:
           if(pattern(pp:pp).eq.'`') then
            if(pp.lt.ep) then !pattern has more symbols to be matched
             call open_predicate(); pp=pp+1
             if(pattern(pp:pp).eq.'`') then; errc=2; exit; endif !mulitple placeholders ` are not allowed in the pattern
            else !pattern is over
             num_pred=num_pred+1
             if(sp.le.es) then
              do while(is_this_blank(str(sp:sp))); sp=sp+1; enddo
              pred_offset(num_pred)=sp; pred_length(num_pred)=es-sp+1
             else
              pred_offset(num_pred)=es+1; pred_length(num_pred)=0
             endif
             match=.TRUE.; exit
            endif
           endif
           if(is_this_blank(pattern(pp:pp))) then
            if(is_this_blank(str(sp:sp))) then !match: both blanks
             if(pred.gt.0) call close_predicate()
             do while(is_this_blank(str(sp:sp))); sp=sp+1; enddo
             do while(is_this_blank(pattern(pp:pp))); pp=pp+1; enddo
            else
             if(pred.gt.0) then
              sp=sp+1
              if(sp.gt.es) then; call restore_previous(over); if(over) exit; endif
             else
              call restore_previous(over); if(over) exit !no match along this path
             endif
            endif
           else
            if(str(sp:sp).eq.pattern(pp:pp)) then !match: symbols matched
             if(pred.gt.0) call close_predicate()
             pp=pp+1; sp=sp+1
             if(pp.gt.ep) then
              if(sp.gt.es) match=.TRUE.
             else
              if(sp.gt.es.and.(pp.lt.ep.or.pattern(pp:pp).ne.'`')) then
               call restore_previous(over); if(over) exit
              endif
             endif
            else
             if(pred.gt.0) then
              sp=sp+1
              if(sp.gt.es) then; call restore_previous(over); if(over) exit; endif
             else
              call restore_previous(over); if(over) exit !no match along this path
             endif
            endif
           endif
          enddo
         else
          errc=1
         endif
         if(.not.match) num_pred=0
         if(present(ierr)) ierr=errc
         return

         contains

          subroutine open_predicate()
           pred=sp; stp=stp+1; pos_stk(:,stp)=(/pp,sp/)
           return
          end subroutine open_predicate

          subroutine close_predicate()
           num_pred=num_pred+1
           pred_offset(num_pred)=pred
           pred_length(num_pred)=sp-pred
           pred=0
           return
          end subroutine close_predicate

          recursive subroutine restore_previous(done)
           logical, intent(out):: done
           done=.FALSE.
           if(pred.gt.0) then
            pred=0; stp=stp-1; call restore_previous(done)
           else
            if(stp.gt.0) then
             pp=pos_stk(1,stp)+1; pred=pos_stk(2,stp)
             sp=pred_offset(num_pred)+pred_length(num_pred)+1; num_pred=num_pred-1
             if(sp.gt.es) call restore_previous(done)
            else
             done=.TRUE.
            endif
           endif
           return
          end subroutine restore_previous

        end function match_symb_pattern
!---------------------------------------------------------
        logical function is_this_index_label(str,so_or_mo)
!Recognizes an index label.
         implicit none
         character(*), intent(in):: str           !in: string
         integer, intent(in), optional:: so_or_mo !in: {ANY_CASE,INDEX_SO,INDEX_MO}
         integer:: i,l,beg_pos,end_pos,som

         is_this_index_label=.FALSE.; som=ANY_CASE
         if(present(so_or_mo)) som=so_or_mo
         call skip_blanks(str,beg_pos,end_pos)
         if(beg_pos.le.end_pos) then
          i=beg_pos; l=end_pos
          if(are_these_letters(str(i:i))) then !single letter: index kind
           i=i+1
           if(i.le.l) then
            if(som.ne.INDEX_SO) then
             if(str(l:l).eq.'a'.or.str(l:l).eq.'b') l=l-1 !spin-label: {a,b}
            endif
            if(i.le.l) then
             if(is_this_integer(str(i:l),no_sign=.TRUE.)) is_this_index_label=.TRUE. !number
            endif
           endif
          endif
         endif
         return
        end function is_this_index_label
!----------------------------------------------------------------------------
        logical function match_index_label(str,so_or_mo,end_pos,beg_pos,spin)
!Matches an index label starting from the given position.
         implicit none
         character(*), intent(in):: str             !in: string
         integer, intent(in):: so_or_mo             !in: {INDEX_SO,INDEX_MO}
         integer, intent(out):: end_pos             !out: end position of the index label (if any)
         integer, intent(in), optional:: beg_pos    !in: beginning position to start matching from (defaults to 1)
         character(1), intent(out), optional:: spin !out: spin {a,b}
         integer:: i,j,bp,l

         match_index_label=.FALSE.; end_pos=0; bp=1
         if(present(beg_pos)) bp=beg_pos
         l=len_trim(str)
         if(bp.gt.0.and.bp.le.l) then
          if(are_these_letters(str(bp:bp))) then !index kind letter
           bp=bp+1
           if(bp.le.l) then
            j=bp
            do while(is_this_integer(str(bp:bp),no_sign=.TRUE.)) !no sign
             bp=bp+1; if(bp.gt.l) exit
            enddo
            if(bp.gt.j) then !at least one number
             if(so_or_mo.ne.INDEX_SO) then
              if(bp.le.l) then
               if(str(bp:bp).eq.'a'.or.str(bp:bp).eq.'b') then !spin label
                if(present(spin)) spin=str(bp:bp)
                end_pos=bp; match_index_label=.TRUE. !MO index label matched
               endif
              endif
             else
              end_pos=bp-1; match_index_label=.TRUE. !SO index label matched
             endif
            endif
           endif
          endif
         endif
         return
        end function match_index_label
!------------------------------------------------------------------------------------------
        logical function is_this_tensor(str,so_or_mo,tens_name,tnl,tens_shape,tsl,inv_pref)
!Recognizes a tensor.
         implicit none
         character(*), intent(in):: str               !in: string
         integer, intent(in):: so_or_mo               !in: {INDEX_SO,INDEX_MO}
         character(*), intent(out):: tens_name        !out: tensor name: tens_name(1:tnl)
         integer, intent(out):: tnl                   !out: length of the tensor name
         character(*), intent(out):: tens_shape       !out: tensor shape: tens_shape(1:tsl)
         integer, intent(out):: tsl                   !out: length of the tensor shape
         integer(8), intent(out), optional:: inv_pref !out: inversed prefactor for unrestricted summations
         integer:: i,j,na,nb,beg_pos,end_pos,lp,rp,sep
         integer(8):: mult
         character(1):: spin,last_kind

         is_this_tensor=.FALSE.; tnl=0; tsl=0; inv_pref=1
         call skip_blanks(str,beg_pos,end_pos)
         if(beg_pos+2.le.end_pos) then !at least one letter for name and an empty '()'
          lp=index(str(beg_pos:end_pos),'(')+beg_pos-1
          if(lp.gt.beg_pos.and.lp.lt.end_pos) then
           rp=index(str(lp+1:end_pos),')')+lp
           if(rp.eq.end_pos) then
            if(are_these_alphanumeric(str(beg_pos:lp-1))) then
             tnl=lp-beg_pos; tens_name(1:tnl)=str(beg_pos:lp-1) !tensor name
             last_kind=' '; sep=0; na=0; nb=0; mult=1; i=lp+1
             do while(i.le.rp)
              if(match_index_label(str,so_or_mo,j,i,spin)) then
               sep=0; tsl=tsl+1; tens_shape(tsl:tsl)=str(i:i)
               last_kind=str(i:i)
               if(so_or_mo.eq.INDEX_MO) then !alpha indices are capitalized, beta are small
                if(spin.eq.'a') then
                 call cap_ascii(tens_shape(tsl:tsl))
                 na=na+1
                else
                 nb=nb+1
                endif
               else
                na=na+1
               endif
               i=j+1
              else
               if(is_this_one_of(str(i:i),(/',','|',')'/))) then !separator
                if(sep.ne.0) return
                if(last_kind.eq.'l'.or.last_kind.eq.'d') then !contracted general indices
                 if(na.gt.0) mult=mult*ifcl(na)
                 if(nb.gt.0) mult=mult*ifcl(nb)
                endif
                sep=1; na=0; nb=0; i=i+1
               else
                return
               endif
              endif
             enddo
             is_this_tensor=.TRUE.
             if(present(inv_pref)) inv_pref=mult
            endif
           endif
          endif
         endif
         return
        end function is_this_tensor
!---------------------------------------------------------------------------------------------
        logical function is_this_tensor_operation(str,op_kind,num_args,tens_args,unrestricted)
!Recognizes basic tensor operations: assignment, scaling, addition, contraction.
         implicit none
         character(*), intent(in):: str                      !in: string
         integer, intent(out):: op_kind                      !out: operation kind
         integer, intent(out):: num_args                     !out: number of arguments
         type(tensor_info_t), intent(inout):: tens_args(0:*) !out: information on each argument
         logical, intent(in), optional:: unrestricted        !in: if TRUE, an unrestricted prefactor will be supplied with tensors
         integer:: i,j,l,beg_pos,end_pos,eq_pos,lhs_end,rhs_beg,mul_pos,num_pref
         integer(8):: inv_pref
         logical:: compos,unres

         is_this_tensor_operation=.FALSE.; op_kind=NO_OPERATION; num_args=0
         if(present(unrestricted)) then; unres=unrestricted; else; unres=.FALSE.; endif
         call skip_blanks(str,beg_pos,end_pos)
         if(beg_pos.le.end_pos) then
          eq_pos=index(str(beg_pos:end_pos),'=')+beg_pos-1
          if(eq_pos.gt.beg_pos) then
           if(str(eq_pos-1:eq_pos-1).eq.'+'.or.str(eq_pos-1:eq_pos-1).eq.'-'.or.&
             &str(eq_pos-1:eq_pos-1).eq.'*'.or.str(eq_pos-1:eq_pos-1).eq.'/') then
            lhs_end=eq_pos-2; rhs_beg=eq_pos+1; compos=.TRUE.
           else
            lhs_end=eq_pos-1; rhs_beg=eq_pos+1; compos=.FALSE.
           endif
           if(lhs_end.ge.beg_pos.and.rhs_beg.le.end_pos) then
            num_pref=0
            if(is_this_tensor(str(beg_pos:lhs_end),INDEX_MO,tens_args(0)%tens_name,tens_args(0)%tens_name_len,&
                             &tens_args(0)%tens_shape,tens_args(0)%tens_rank,inv_pref)) then !destination tensor
             if(unres) then
              tens_args(num_args)%prefactor=1d0/dble(inv_pref)
             else
              tens_args(num_args)%prefactor=1d0
             endif
             num_args=1; i=rhs_beg
             do while(i.le.end_pos)
              mul_pos=index(str(i:end_pos),'*')+i-1
              if(mul_pos.ge.i) then
               if(mul_pos.eq.end_pos) return !'*' cannot be the last token
               l=mul_pos-1; j=mul_pos+1
              else
               l=end_pos; j=end_pos+1
              endif
              if(is_this_tensor(str(i:l),INDEX_MO,tens_args(num_args)%tens_name,tens_args(num_args)%tens_name_len,&
                               &tens_args(num_args)%tens_shape,tens_args(num_args)%tens_rank,inv_pref)) then !rhs tensor argument
               if(unres) then
                tens_args(num_args)%prefactor=1d0/dble(inv_pref)
               else
                tens_args(num_args)%prefactor=1d0
               endif
               num_args=num_args+1
              elseif(is_this_real_number(str(i:l),tens_args(num_args)%prefactor)) then !prefactor
               tens_args(num_args)%tens_name_len=0; tens_args(num_args)%tens_rank=0
               num_args=num_args+1; num_pref=num_pref+1
              else
               return
              endif
              i=j
             enddo
             if(compos) then
              if(num_args-num_pref.eq.2) then
               op_kind=TENSOR_ADDITION
              elseif(num_args-num_pref.eq.3) then
               op_kind=TENSOR_CONTRACTION
              else
               return
              endif
             else
              if(num_args-num_pref.eq.1) then
               op_kind=TENSOR_ASSIGNMENT
              elseif(num_args-num_pref.eq.2) then
               op_kind=TENSOR_SCALING
              else
               return
              endif
             endif
             do i=0,num_args-1
              call translate_index_labels(tens_args(i)%tens_shape(1:tens_args(i)%tens_rank))
             enddo
             is_this_tensor_operation=.TRUE.
            endif
           endif
          endif
         endif
         return
        end function is_this_tensor_operation
!-----------------------------------------------------------------------------------
        logical function is_this_auxiliary_operation(str,op_kind,num_args,tens_args)
!Recognizes custom auxiliary operations.
         implicit none
         character(*), intent(in):: str                      !in: string
         integer, intent(out):: op_kind                      !out: operation kind
         integer, intent(out):: num_args                     !out: number of arguments
         type(tensor_info_t), intent(inout):: tens_args(0:*) !out: information on each argument
         integer:: beg_pos,end_pos,i

         is_this_auxiliary_operation=.FALSE.; op_kind=NO_OPERATION; num_args=0
         call skip_blanks(str,beg_pos,end_pos)
         if(beg_pos.le.end_pos) then
 !TRY "CREATE_ARRAY":
          i=index(str(beg_pos:end_pos),'CREATE_ARRAY')+beg_pos-1
          if(i.ge.beg_pos) then
           i=i+len('CREATE_ARRAY')
           do while(i.le.end_pos)
            if(.not.is_this_blank(str(i:i))) exit
            i=i+1
           enddo
           if(i.le.end_pos) then
            is_this_auxiliary_operation=is_this_tensor(str(i:end_pos),INDEX_MO,&
                                       &tens_args(num_args)%tens_name,tens_args(num_args)%tens_name_len,&
                                       &tens_args(num_args)%tens_shape,tens_args(num_args)%tens_rank)
            if(is_this_auxiliary_operation) then
             num_args=num_args+1; op_kind=TENSOR_CREATE
            endif
           endif
           return
          endif
 !TRY "DELETE_ARRAY":
          i=index(str(beg_pos:end_pos),'DELETE_ARRAY')+beg_pos-1
          if(i.ge.beg_pos) then
           i=i+len('DELETE_ARRAY')
           do while(i.le.end_pos)
            if(.not.is_this_blank(str(i:i))) exit
            i=i+1
           enddo
           if(i.le.end_pos) then
            is_this_auxiliary_operation=is_this_tensor(str(i:end_pos),INDEX_MO,&
                                       &tens_args(num_args)%tens_name,tens_args(num_args)%tens_name_len,&
                                       &tens_args(num_args)%tens_shape,tens_args(num_args)%tens_rank)
            if(is_this_auxiliary_operation) then
             num_args=num_args+1; op_kind=TENSOR_DELETE
            endif
           endif
           return
          endif
 !TRY "BARRIER"
          i=index(str(beg_pos:end_pos),'BARRIER')+beg_pos-1
          if(i.ge.beg_pos) then
           op_kind=BARRIER_SYNC; is_this_auxiliary_operation=.TRUE.
           return
          endif
         endif
         return
        end function is_this_auxiliary_operation
!---------------------------------------------------------------
        subroutine remove_comments(stri,comm_chars,str_len,stro)
!Removes comments from a string.
         implicit none
         character(*), intent(in):: stri            !in: input string
         character(1), intent(in):: comm_chars(:)   !in: comment activating characters, e.g. '#', '!', etc.: Cannot be blank
         integer, intent(out):: str_len             !out: length of the string without comments
         character(*), intent(out), optional:: stro !out: output string (without comments)
         integer:: beg_pos,end_pos,i,n

         str_len=0
         call skip_blanks(stri,beg_pos,end_pos)
         if(beg_pos.le.end_pos) then
          str_len=end_pos; n=size(comm_chars)
          if(n.gt.0) then
           do i=beg_pos,end_pos
            if(is_this_one_of(stri(i:i),comm_chars)) then
             str_len=i-1
             exit
            endif
           enddo
           do while(str_len.ge.beg_pos)
            if(.not.is_this_blank(stri(str_len:str_len))) exit
            str_len=str_len-1
           enddo
           if(str_len.lt.beg_pos) str_len=0
          endif
          if(str_len.gt.0) then
           if(present(stro)) stro(1:str_len)=stri(1:str_len)
          endif
         endif
         return
        end subroutine remove_comments
!-----------------------------------------------------------------------------------------
        logical function parse_tadl_statement(str,op_kind,num_args,tens_args,unrestricted)
!Parses Tensor Algebra Descriptive Language statements (for DiaGen).
         implicit none
         character(*), intent(in):: str                      !in: string
         integer, intent(out):: op_kind                      !out: operation kind
         integer, intent(out):: num_args                     !out: number of arguments
         type(tensor_info_t), intent(inout):: tens_args(0:*) !out: information on each argument
         logical, intent(in), optional:: unrestricted        !in: if TRUE, an unrestricted summation prefactors will be computed for each tensor
         integer:: beg_pos,end_pos,l
         logical:: unres

         parse_tadl_statement=.FALSE.; op_kind=NO_OPERATION; num_args=0
         if(present(unrestricted)) then; unres=unrestricted; else; unres=.FALSE.; endif
         call skip_blanks(str,beg_pos,end_pos)
         if(beg_pos.le.end_pos) then
          call remove_comments(str,(/'#'/),l)
          if(l.gt.0) then
           parse_tadl_statement=is_this_tensor_operation(str,op_kind,num_args,tens_args,unres)
           if(parse_tadl_statement) return
           op_kind=NO_OPERATION; num_args=0
           parse_tadl_statement=is_this_auxiliary_operation(str,op_kind,num_args,tens_args)
           if(parse_tadl_statement) return
           op_kind=NO_OPERATION; num_args=0
          elseif(l.eq.0) then
           op_kind=COMMENT_LINE; parse_tadl_statement=.TRUE.
          endif
         else
          op_kind=COMMENT_LINE; parse_tadl_statement=.TRUE.
         endif
         return
        end function parse_tadl_statement
!----------------------------------------------
        subroutine print_tensor_info(tens_info)
         implicit none
         type(tensor_info_t), intent(in):: tens_info !in: tensor info

         call printl(6,tens_info%tens_name(1:tens_info%tens_name_len)//': ',ADV=.FALSE.)
         call printl(6,tens_info%tens_shape(1:tens_info%tens_rank)//' * ',ADV=.FALSE.)
         write(*,*) tens_info%prefactor
         return
        end subroutine print_tensor_info
!----------------------------------------------------
        subroutine translate_index_labels(tens_shape)
!Replaces index letters with "o","O","v","V" explanations (for DiaGen).
         implicit none
         character(*), intent(inout):: tens_shape !in: string of index kind labels (single letter per tensor dimension)
         integer:: i

         do i=1,len_trim(tens_shape)
          select case(tens_shape(i:i))
          case('m','n','l') !full occupied range (beta)
           tens_shape(i:i)='o'
          case('M','N','L') !full occupied range (alpha)
           tens_shape(i:i)='O'
          case('e','f','d') !full virtual range (beta)
           tens_shape(i:i)='v'
          case('E','F','D') !full virtual range (alpha)
           tens_shape(i:i)='V'
          end select
         enddo
         return
        end subroutine translate_index_labels

       end module parse_prim
![TESTING]=========================================
       module parse_prim_test
        use stsubs
        use parse_prim
        implicit none

        public test_parse_prim

       contains
!----------------------------------------------
        function test_parse_prim() result(ierr)
         implicit none
         integer:: i,ierr,npred,offs(1:128),lens(1:128)
         character(128):: str,ind1,ind2,ind3
         logical:: match

         ierr=0
         write(*,'("Testing symbolic parsing ... ")',ADVANCE='NO')
         !write(*,'()') !debug
         str='D12(a1,b1)+=L451(k1,b1,l1)*R6(a1,l1,k1)*0.34'
         match=match_symb_pattern(str,'`(`)+=`(`)*`(`)*`',npred,offs,lens,ierr)
         if(match.and.ierr.eq.0) then
          !do i=1,npred; call printl(6,str(offs(i):offs(i)+lens(i)-1)); enddo !debug
          ind1=str(offs(2):offs(2)+lens(2)-1)
          ind2=str(offs(4):offs(4)+lens(4)-1)
          ind3=str(offs(6):offs(6)+lens(6)-1)
          match=match_symb_pattern(ind1,'`,`',npred,offs,lens,ierr)
          if(match.and.ierr.eq.0) then
           !do i=1,npred; call printl(6,ind1(offs(i):offs(i)+lens(i)-1)); enddo !debug
           match=match_symb_pattern(ind2,'`,`',npred,offs,lens,ierr)
           if(match.and.ierr.eq.0) then
            !do i=1,npred; call printl(6,ind2(offs(i):offs(i)+lens(i)-1)); enddo !debug
            match=match_symb_pattern(ind3,'`,`',npred,offs,lens,ierr)
            if(match.and.ierr.eq.0) then
             !do i=1,npred; call printl(6,ind3(offs(i):offs(i)+lens(i)-1)); enddo !debug
             write(*,'("PASSED")')
            else
             write(*,'("FAILED")'); ierr=4
            endif
           else
            write(*,'("FAILED")'); ierr=3
           endif
          else
           write(*,'("FAILED")'); ierr=2
          endif
         else
          write(*,'("FAILED")'); ierr=1
         endif
         return
        end function test_parse_prim

       end module parse_prim_test
