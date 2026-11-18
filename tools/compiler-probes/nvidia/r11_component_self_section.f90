! Intrinsic assignment of a section of an allocatable array to the array
! itself must reallocate it to the section's shape (Fortran 2008 7.2.1.3).
module r11_mod
  implicit none
  type :: holder_t
     integer, allocatable :: v(:)
  end type holder_t
contains
  subroutine shrink_dummy(v)
    integer, allocatable, intent(inout) :: v(:)
    v = v(1:2)
  end subroutine shrink_dummy
end module r11_mod
program r11
  use r11_mod
  implicit none
  integer, allocatable :: a(:), b(:)
  type(holder_t) :: h
  integer :: failures
  failures = 0
  a = [1, 2, 3, 4]
  a = a(1:2)
  if (size(a) /= 2) failures = failures + 1
  write (*, '(A,I0)') "variable:  size ", size(a)
  h%v = [1, 2, 3, 4]
  h%v = h%v(1:2)
  if (size(h%v) /= 2) failures = failures + 1
  write (*, '(A,I0)') "component: size ", size(h%v)
  b = [1, 2, 3, 4]
  call shrink_dummy(b)
  if (size(b) /= 2) failures = failures + 1
  write (*, '(A,I0)') "dummy:     size ", size(b)
  if (failures == 0) then
     write (*, '(A)') "RESULT: all passed"
  else
     write (*, '(A)') "RESULT: FAIL"
     error stop 1
  end if
end program r11
