! MOVE_ALLOC between two local unlimited polymorphic allocatable scalars.
program r2b
  implicit none
  class(*), allocatable :: a, b
  allocate(a, source=1.5d0)
  call move_alloc(a, b)
  select type (x => b)
  type is (real(8))
     if (x == 1.5d0 .and. .not. allocated(a)) then
        write (*, '(A)') "RESULT: all passed"
        stop
     end if
  end select
  write (*, '(A)') "RESULT: FAIL"
  error stop 1
end program r2b
