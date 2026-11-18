! Control for r3a: same SOURCE= expression, non-polymorphic allocate-object.
program r3b
  implicit none
  real(8), allocatable :: a(:)
  allocate(a, source=[real(8) :: 1, 2])
  if (size(a) == 2 .and. a(2) == 2d0) then
     write (*, '(A)') "RESULT: all passed"
  else
     write (*, '(A)') "RESULT: FAIL"
     error stop 1
  end if
end program r3b
