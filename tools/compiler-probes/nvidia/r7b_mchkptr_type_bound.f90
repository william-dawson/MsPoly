module r7b_mod
  implicit none
  type :: holder_t
     class(*), allocatable :: a(:)
   contains
     generic :: set => set_r1, set_i1
     procedure :: set_r1, set_i1
  end type holder_t
contains
  subroutine set_r1(this, value)
    class(holder_t), intent(inout) :: this
    real(8), dimension(:), intent(in) :: value
    if (allocated(this%a)) deallocate(this%a)
    allocate(real(8) :: this%a(size(value)))
    select type (x => this%a)
    type is (real(8))
       x = value
    end select
  end subroutine set_r1
  subroutine set_i1(this, value)
    class(holder_t), intent(inout) :: this
    integer, dimension(:), intent(in) :: value
    if (allocated(this%a)) deallocate(this%a)
    allocate(integer :: this%a(size(value)))
    select type (x => this%a)
    type is (integer)
       x = value
    end select
  end subroutine set_i1
end module r7b_mod
program r7b
  use r7b_mod
  implicit none
  call first()
  call second()
  write (*, '(A)') "RESULT: all passed"
contains
  subroutine first()
    type(holder_t) :: h
    call h%set([1d0, 2d0])
    call h%set([1])
  end subroutine first
  subroutine second()
    type(holder_t) :: h
    call h%set([real(8) ::])
  end subroutine second
end program r7b
