! MOVE_ALLOC between two local unlimited polymorphic allocatable arrays.
! Valid Fortran 2008 (13.7.118): FROM and TO may be unlimited polymorphic.
program r2a
  implicit none
  class(*), allocatable :: a(:), b(:)
  real(8) :: z(3) = [1d0, 2d0, 3d0]
  allocate(a, source=z)
  call move_alloc(a, b)
  select type (x => b)
  type is (real(8))
     if (all(x == z) .and. .not. allocated(a)) then
        write (*, '(A)') "RESULT: all passed"
        stop
     end if
  end select
  write (*, '(A)') "RESULT: FAIL"
  error stop 1
end program r2a
