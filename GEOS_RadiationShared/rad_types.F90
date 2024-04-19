module rad_types

  ! Handy wrappers to pointer arrays.
  ! Used to provide ragged arrays, etc.
  type :: rptr1d_wrap
    real, pointer, dimension(:) :: p
  end type rptr1d_wrap

end module rad_types
