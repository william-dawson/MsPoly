! As r3b with MOLD= instead of SOURCE=.
program r3e
  implicit none
  real(8), allocatable :: a(:)
  allocate(a, mold=[real(8) :: 1, 2])
  if (size(a) == 2) then
     write (*, '(A)') "RESULT: all passed"
  else
     write (*, '(A)') "RESULT: FAIL"
     error stop 1
  end if
end program r3e
