! Intrinsic assignment of a derived type that has an unlimited polymorphic
! allocatable component, into an element of an array, loses the component's
! value (Fortran 2008 10.2.1.3: allocatable components are deep copied).
module r1_mod
  implicit none
  type :: box_t
     class(*), allocatable :: v
  end type box_t
  type :: real_box_t                  ! control: non-polymorphic component
     real(8), allocatable :: v
  end type real_box_t
  integer :: failures = 0
contains
  subroutine check(label, b)
    character(*), intent(in) :: label
    type(box_t), intent(in) :: b
    character(len=60) :: what
    if (.not. allocated(b%v)) then
       what = "FAIL: component not allocated"
    else
       select type (x => b%v)
       type is (real(8))
          if (x == 1.5d0) then
             what = "PASS"
          else
             write (what, '(A,ES10.3)') "FAIL: value is ", x
          end if
       class default
          what = "FAIL: dynamic type lost"
       end select
    end if
    if (what(1:4) /= "PASS") failures = failures + 1
    write (*, '(A,T45,A)') label, trim(what)
  end subroutine check
  subroutine assign_into_dummy(arr, b)
    type(box_t), intent(inout) :: arr(:)
    type(box_t), intent(in) :: b
    arr(1) = b
  end subroutine assign_into_dummy
end module r1_mod

program r1
  use r1_mod
  implicit none
  type(box_t) :: a, scalar, fixed(2)
  type(box_t), allocatable :: alloc(:)
  type(real_box_t) :: ra
  type(real_box_t), allocatable :: ralloc(:)

  allocate(a%v, source=1.5d0)

  scalar = a
  call check("scalar = scalar", scalar)
  fixed(1) = a
  call check("fixed-size array element = scalar", fixed(1))
  allocate(alloc(2))
  alloc(1) = a
  call check("allocatable array element = scalar", alloc(1))
  call assign_into_dummy(alloc, a)
  call check("dummy array element = scalar", alloc(1))

  ra%v = 1.5d0                          ! control without class(*)
  allocate(ralloc(2))
  ralloc(1) = ra
  if (ralloc(1)%v /= 1.5d0) failures = failures + 1
  write (*, '(A,T45,A)') "control: real component, alloc array element", &
       merge("PASS", "FAIL", ralloc(1)%v == 1.5d0)

  if (failures > 0) then
     write (*, '(A,I0,A)') "RESULT: ", failures, " failure(s)"
     error stop 1
  end if
  write (*, '(A)') "RESULT: all passed"
end program r1
