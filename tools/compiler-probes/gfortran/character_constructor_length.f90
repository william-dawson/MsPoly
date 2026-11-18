! Fortran 2008 4.8: the type-spec fixes each element's character length.
program character_constructor
  implicit none
  integer :: i
  character(64) :: header
  do i = 1, 2
     header = "%%MatrixMarket matrix coordinate real general"
     if (i == 2) header = &
          & "%%MatrixMarket matrix coordinate complex general"
     call verify([character(64) :: header, "2 2 1"])
     call verify([character(64) :: &
          & "%%MatrixMarket matrix coordinate " // &
          & trim(merge("real   ", "complex", i == 1)) // " general", &
          & "2 2 1"])
  end do
  print *, "RESULT: all passed"
contains
  subroutine verify(lines)
    character(*), intent(in) :: lines(:)
    if (size(lines) /= 2 .or. len(lines) /= 64) error stop 1
    if (trim(lines(2)) /= "2 2 1") error stop 2
  end subroutine verify
end program character_constructor
