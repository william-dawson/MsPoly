! Defined assignment a = a for a type with an allocatable component, where
! the assignment builds the new contents first and then replaces the old
! (Fortran 2008 10.2.1.4: the right-hand side is passed as (a), a value).
module m
  implicit none
  type :: box_t
     integer, allocatable :: v(:)
   contains
     generic :: assignment(=) => assign
     procedure :: assign
  end type
contains
  subroutine assign(lhs, rhs)
    class(box_t), intent(inout) :: lhs
    class(box_t), intent(in) :: rhs
    integer, allocatable :: copied(:)
    copied = rhs%v
    if (allocated(lhs%v)) deallocate(lhs%v)
    call move_alloc(copied, lhs%v)
  end subroutine
end module
program self_assignment
  use m
  implicit none
  type(box_t) :: a
  a%v = [1, 2, 3]
  a = a
  if (allocated(a%v)) then
     if (size(a%v) == 3) then
        if (all(a%v == [1, 2, 3])) then
           write (*, '(A)') "RESULT: all passed"
           stop
        end if
     end if
  end if
  write (*, '(A)') "RESULT: FAIL"
  error stop 1
end program
