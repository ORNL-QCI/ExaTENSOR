!TEMPLATE:
!This function associates a pointer "tptr" of class T
!with a target pointed to by an unlimited polymorphic
!pointer "uptr" only if the latter points to a target
!of class T, otherwise NULL is returned.
template <typename T>
function T_cast(uptr) result(tptr)
 implicit none
 class(T), pointer:: tptr              !out: typed pointer
 class(*), pointer, intent(in):: uptr  !in: unlimited polymorphic pointer

 tptr=>NULL()
 if(associated(uptr)) then
  select type(uptr)
  class is(T)
   tptr=>uptr
  end select
 endif
 return
end function T_cast
