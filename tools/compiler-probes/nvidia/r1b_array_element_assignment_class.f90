! As r1, but the component is CLASS(base_t) (not unlimited polymorphic).
module r1b_mod
  implicit none
  type :: base_t
     real(8) :: x = 0
  end type base_t
  type, extends(base_t) :: child_t
     integer :: tag = 0
  end type child_t
  type :: box_t
     class(base_t), allocatable :: v
  end type box_t
end module r1b_mod

program r1b
  use r1b_mod
  implicit none
  type(box_t) :: a
  type(box_t), allocatable :: alloc(:)
  logical :: ok

  allocate(a%v, source=child_t(x=1.5d0, tag=7))
  allocate(alloc(2))
  alloc(1) = a
  ok = .false.
  if (allocated(alloc(1)%v)) then
     select type (c => alloc(1)%v)
     type is (child_t)
        ok = c%x == 1.5d0 .and. c%tag == 7
     end select
  end if
  if (ok) then
     write (*, '(A)') "RESULT: all passed"
  else
     write (*, '(A)') "RESULT: FAIL (allocatable array element = scalar)"
     error stop 1
  end if
end program r1b
