! Deferred-length character assignment whose right-hand side contains the
! variable itself and a call to ACHAR.
program r6
  implicit none
  character(:), allocatable :: path
  integer :: k
  path = "l1"
  do k = 2, 4
     path = path // ".l" // achar(iachar('0') + k)
  end do
  write (*, '(3A,I0)') "path=[", path, "] len=", len(path)
  if (path == "l1.l2.l3.l4" .and. len(path) == 11) then
     write (*, '(A)') "RESULT: all passed"
  else
     write (*, '(A)') "RESULT: FAIL (expected [l1.l2.l3.l4] len=11)"
     error stop 1
  end if
end program r6
