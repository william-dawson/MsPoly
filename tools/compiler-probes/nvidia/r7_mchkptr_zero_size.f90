! -Mchkptr: a zero-size array constructor passed as an actual argument,
! after another procedure in the program unit passed a non-empty one.
module r7_mod
  implicit none
contains
  subroutine take(x)
    real(8), intent(in) :: x(:)
    write (*, '(A,I0)') "size=", size(x)
  end subroutine take
end module r7_mod

program r7
  use r7_mod
  implicit none
  call first()
  call second()
  write (*, '(A)') "RESULT: all passed"
contains
  subroutine first()
    call take([1d0, 2d0])
  end subroutine first
  subroutine second()
    call take([real(8) ::])
  end subroutine second
end program r7
