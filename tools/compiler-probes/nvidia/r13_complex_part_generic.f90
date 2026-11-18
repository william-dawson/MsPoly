! Generic resolution with a complex part designator as the actual argument.
! z%re for a complex array z is a real array (F2008 6.4.4: the real part
! of a complex designator has the kind and type real), so a generic with a
! real and a complex specific must resolve to the real one.
module pick_mod
  implicit none
  private
  public :: pick

  interface pick
     module procedure pick_r, pick_c
  end interface pick
contains
  integer function pick_r(v)
    real(8), dimension(:), intent(in) :: v
    pick_r = 1
    if (size(v) < 0) pick_r = 0
  end function pick_r

  integer function pick_c(v)
    complex(8), dimension(:), intent(in) :: v
    pick_c = 2
    if (size(v) < 0) pick_c = 0
  end function pick_c
end module pick_mod

program r13_complex_part_generic
  use pick_mod, only: pick
  implicit none
  complex(8), dimension(3) :: z

  z = (1.0d0, 2.0d0)
  if (pick(z%re) == 1 .and. pick(z(1:2)%re) == 1) then
     print '(a)', "RESULT: all passed"
  else
     print '(a)', "RESULT: complex part resolved to the complex specific"
     error stop 1
  end if
end program r13_complex_part_generic
