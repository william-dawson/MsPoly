! ALLOCATE with SOURCE= an array constructor with a type-spec, into an
! unlimited polymorphic array. Valid Fortran 2008 (6.7.1, 4.8).
program r3a
  implicit none
  class(*), allocatable :: a(:)
  allocate(a, source=[real(8) :: 1, 2])
  select type (x => a)
  type is (real(8))
     if (size(x) == 2 .and. x(2) == 2d0) then
        write (*, '(A)') "RESULT: all passed"
        stop
     end if
  end select
  write (*, '(A)') "RESULT: FAIL"
  error stop 1
end program r3a
