! As self_assignment.f90, but the type holds an allocatable array of
! entries, each with a class(*) allocatable component, and the defined
! assignment copies the entries explicitly.
module m
  implicit none
  type :: entry_t
     class(*), allocatable :: v
  end type
  type :: box_t
     type(entry_t), allocatable :: e(:)
   contains
     generic :: assignment(=) => assign
     procedure :: assign
  end type
contains
  subroutine assign(lhs, rhs)
    class(box_t), intent(inout) :: lhs
    class(box_t), intent(in) :: rhs
    type(entry_t), allocatable :: copied(:)
    integer :: k
    allocate(copied(size(rhs%e)))
    do k = 1, size(rhs%e)
       select type (x => rhs%e(k)%v)
       type is (integer)
          allocate(copied(k)%v, source=x)
       end select
    end do
    if (allocated(lhs%e)) deallocate(lhs%e)
    call move_alloc(copied, lhs%e)
  end subroutine
end module
program self_assignment_class_star
  use m
  implicit none
  type(box_t) :: a
  integer :: k
  allocate(a%e(3))
  do k = 1, 3
     allocate(a%e(k)%v, source=k)
  end do
  a = a
  k = 0
  if (allocated(a%e)) then
     select type (x => a%e(3)%v)
     type is (integer)
        k = x
     end select
  end if
  if (k == 3) then
     write (*, '(A)') "RESULT: all passed"
  else
     write (*, '(A)') "RESULT: FAIL"
     error stop 1
  end if
end program
