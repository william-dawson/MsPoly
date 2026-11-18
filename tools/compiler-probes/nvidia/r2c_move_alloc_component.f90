! Same as r2a, but FROM and TO are components of a derived type.
program r2c
  implicit none
  type :: holder_t
     class(*), allocatable :: a(:)
  end type holder_t
  type(holder_t) :: p, q
  real(8) :: z(3) = [1d0, 2d0, 3d0]
  allocate(p%a, source=z)
  call move_alloc(p%a, q%a)
  select type (x => q%a)
  type is (real(8))
     if (all(x == z) .and. .not. allocated(p%a)) then
        write (*, '(A)') "RESULT: all passed"
        stop
     end if
  end select
  write (*, '(A)') "RESULT: FAIL"
  error stop 1
end program r2c
