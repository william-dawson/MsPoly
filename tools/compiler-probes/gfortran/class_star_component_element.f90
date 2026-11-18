! Passing an element of a class(*) character array COMPONENT, reached
! through a derived-type dummy, to a class(*) dummy.
module m
  implicit none
  type :: entry_t
     class(*), allocatable :: a(:)
  end type
contains
  subroutine show(x)
    class(*), intent(in) :: x
    select type (x)
    type is (character(*)); write (*, '(3A,I0)') "  [", x, "] len=", len(x)
    end select
  end subroutine
  subroutine show_all(en)
    type(entry_t), intent(in) :: en
    integer :: k
    do k = 1, size(en%a)
       call show(en%a(k))
    end do
  end subroutine
end module
program p
  use m
  implicit none
  type(entry_t) :: e
  character(len=2) :: names(3) = ["ab", "cd", "ef"]
  integer :: k
  allocate(character(len=2) :: e%a(3))
  select type (x => e%a); type is (character(*)); x = names; end select
  write (*, '(A)') "direct (expect [ab] [cd] [ef]):"
  do k = 1, 3
     call show(e%a(k))
  end do
  write (*, '(A)') "through a type(entry_t) dummy:"
  call show_all(e)
end program
