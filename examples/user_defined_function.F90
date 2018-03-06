!This example illustrates how to write a user-defined (stateless) function
!callable from ExaTENSOR to initialize or transform a tensor block.

function compute_2e_integrals(tens_body,data_kind,tens_shape,tens_base) result(ierr)
 use, intrinsic:: ISO_C_BINDING
 use exatensor
 integer(C_INT):: ierr                             !out: error code (0:success)
 type(C_PTR), value:: tens_body                    !in: C pointer to the tensor body (tensor elements)
 integer(C_INT), value:: data_kind                 !in: data kind: {EXA_DATA_KIND_R4,EXA_DATA_KIND_R8,EXA_DATA_KIND_C4,EXA_DATA_KIND_C8}
 integer(C_INT), intent(in):: tens_shape(1:*)      !in: tensor shape (extent of each tensor dimension)
 integer(C_LONG_LONG), intent(in):: tens_base(1:*) !in: base offsets for each tensor dimension)

 complex(8), pointer:: body(:,:,:,:)

 ierr=0 !set error code to success
 if(data_kind == EXA_DATA_KIND_C8) then
  call c_f_pointer(tens_body,body,tens_shape(1:4)) !map C pointer to 4-dimensional array

  !Here we compute the tensor body:
  !body(1:tens_shape(1),1:tens_shape(2),1:tens_shape(3),1:tens_shape(4))
  !Note that local numeration starts from 1, but globally this is a slice
  !of a larger tensor at position tens_signature(1:4).

 else
  ierr=-1 !invalid data kind requested
 endif
 return
end function compute_2e_integrals
