! LEN of ACHAR with a constant and with a run-time argument.
program r6c
  implicit none
  integer :: i
  character(:), allocatable :: s
  i = 50
  s = achar(i)
  write (*, '(A,I0,A,I0,A,I0)') "len(achar(50))=", len(achar(50)), &
       "  len(achar(i)) direct=", len(achar(i)), "  len(s=achar(i))=", len(s)
  if (len(s) == 1 .and. len(achar(i)) == 1) then
     write (*, '(A)') "RESULT: all passed"
  else
     write (*, '(A)') "RESULT: FAIL (ACHAR result must have length 1)"
     error stop 1
  end if
end program r6c
