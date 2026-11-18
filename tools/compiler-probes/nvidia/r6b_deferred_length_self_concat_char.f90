! Control for r6: same loop without ACHAR.
program r6b
  implicit none
  character(:), allocatable :: path
  character(len=1), parameter :: digits(3) = ["2", "3", "4"]
  integer :: k
  path = "l1"
  do k = 1, 3
     path = path // ".l" // digits(k)
  end do
  write (*, '(3A,I0)') "path=[", path, "] len=", len(path)
  if (path == "l1.l2.l3.l4" .and. len(path) == 11) then
     write (*, '(A)') "RESULT: all passed"
  else
     write (*, '(A)') "RESULT: FAIL"
     error stop 1
  end if
end program r6b
