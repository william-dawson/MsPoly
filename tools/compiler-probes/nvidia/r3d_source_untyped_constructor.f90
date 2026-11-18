! Control for r3b: the same allocation without a type-spec parses and runs.
program r3d
  implicit none
  real(8), allocatable :: a(:)
  allocate(a, source=[1d0, 2d0])
  if (size(a) == 2 .and. a(2) == 2d0) then
     write (*, '(A)') "RESULT: all passed"
  else
     write (*, '(A)') "RESULT: FAIL"
     error stop 1
  end if
end program r3d
