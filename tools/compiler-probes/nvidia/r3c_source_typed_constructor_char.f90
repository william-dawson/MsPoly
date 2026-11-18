! As r3a with a character type-spec.
program r3c
  implicit none
  class(*), allocatable :: a(:)
  allocate(a, source=[character(len=5) :: "alpha", "beta"])
  select type (x => a)
  type is (character(*))
     if (size(x) == 2 .and. len(x) == 5 .and. x(2) == "beta") then
        write (*, '(A)') "RESULT: all passed"
        stop
     end if
  end select
  write (*, '(A)') "RESULT: FAIL"
  error stop 1
end program r3c
