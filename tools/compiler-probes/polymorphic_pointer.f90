! A class(base) pointer associated with either a local class(base)
! allocatable or a polymorphic allocatable component of a target dummy,
! then used for dispatch.
module ptr_mod
  implicit none
  private
  public :: base_t, real_t, holder_t, total_of
  type, abstract :: base_t
   contains
     procedure(total_i), deferred :: total
  end type base_t
  abstract interface
     real(8) function total_i(this)
       import :: base_t
       class(base_t), intent(in) :: this
     end function total_i
  end interface
  type, extends(base_t) :: real_t
     real(8), allocatable :: v(:)
   contains
     procedure :: total => real_total
  end type real_t
  type :: holder_t
     class(base_t), allocatable :: local
  end type holder_t
contains
  real(8) function real_total(this)
    class(real_t), intent(in) :: this
    real_total = sum(this%v)
  end function real_total
  real(8) function total_of(h, use_copy)
    type(holder_t), intent(in), target :: h
    logical, intent(in) :: use_copy
    class(base_t), allocatable, target :: copy
    class(base_t), pointer :: p
    p => h%local
    if (use_copy) then
       allocate(real_t :: copy)
       select type (copy)
       type is (real_t)
          copy%v = [10d0, 20d0]
       end select
       p => copy
    end if
    total_of = p%total()
  end function total_of
end module ptr_mod
program polymorphic_pointer
  use ptr_mod
  implicit none
  type(holder_t) :: h
  allocate(real_t :: h%local)
  select type (s => h%local)
  type is (real_t)
     s%v = [1d0, 2d0]
  end select
  if (total_of(h, .false.) == 3d0 .and. total_of(h, .true.) == 30d0) then
     write (*, '(A)') "RESULT: all passed"
  else
     write (*, '(A)') "RESULT: FAIL"
     error stop 1
  end if
end program polymorphic_pointer
