! TRANSFER with an allocated deferred-length character array as MOLD.
program r9
  use, intrinsic :: iso_fortran_env, only: int8
  implicit none
  integer(int8) :: bytes(6) = int([97, 98, 99, 100, 101, 102], int8)
  character(:), allocatable :: s(:)
  allocate(character(len=3) :: s(2))
  s = transfer(bytes, s, 2)
  if (s(1) == "abc" .and. s(2) == "def") then
     write (*, '(A)') "RESULT: all passed"
  else
     write (*, '(A)') "RESULT: FAIL"
     error stop 1
  end if
end program r9
