! ALLOCATE a class(*) component with SOURCE= the typed associate name of
! SELECT TYPE on another class(*) component. Valid Fortran 2008 (6.7.1.1).
module r10_mod
  implicit none
  type :: entry_t
     class(*), allocatable :: v
  end type entry_t
contains
  subroutine copy_typed(dst, src)
    type(entry_t), intent(inout) :: dst
    type(entry_t), intent(in) :: src
    select type (x => src%v)
    type is (real(8))
       allocate(dst%v, source=x)
    end select
  end subroutine copy_typed
  subroutine copy_default(dst, src)
    type(entry_t), intent(inout) :: dst
    type(entry_t), intent(in) :: src
    select type (x => src%v)
    class default
       allocate(dst%v, source=x)
    end select
  end subroutine copy_default
  logical function holds(e, value)
    type(entry_t), intent(in) :: e
    real(8), intent(in) :: value
    holds = .false.
    if (.not. allocated(e%v)) return
    select type (x => e%v)
    type is (real(8))
       holds = x == value
    end select
  end function holds
end module r10_mod
program r10
  use r10_mod
  implicit none
  type(entry_t) :: a, b, c
  allocate(a%v, source=1.5d0)
  call copy_typed(b, a)
  allocate(c%v, source=1.5d0)
  write (*, '(A,L2,A,L2)') "typed branch:", holds(b, 1.5d0), "  class default:", holds(c, 1.5d0)
  if (holds(b, 1.5d0) .and. holds(c, 1.5d0)) then
     write (*, '(A)') "RESULT: all passed"
  else
     write (*, '(A)') "RESULT: FAIL"
     error stop 1
  end if
end program r10
