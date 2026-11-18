! Intrinsic assignment from a structure constructor that omits components.
! Omitted components with default initialization must take that value
! (Fortran 2008 4.5.10), so used must become .false. below.
program r4
  implicit none
  type :: entry_t
     character(:), allocatable :: key
     class(*), allocatable :: v
     logical :: used = .false.
  end type entry_t
  type :: plain_t                        ! control: no class(*) component
     character(:), allocatable :: key
     real(8), allocatable :: v
     logical :: used = .false.
  end type plain_t
  type(entry_t) :: e
  type(plain_t) :: p
  integer :: failures = 0

  allocate(e%v, source=1.5d0)
  e%used = .true.
  e = entry_t()
  write (*, '(A,L2,A,L2)') "entry_t(): v allocated =", allocated(e%v), &
       "  used =", e%used
  if (allocated(e%v) .or. e%used) failures = failures + 1

  allocate(p%v, source=1.5d0)
  p%used = .true.
  p = plain_t()
  write (*, '(A,L2,A,L2)') "control plain_t(): v allocated =", &
       allocated(p%v), "  used =", p%used
  if (allocated(p%v) .or. p%used) failures = failures + 1

  if (failures > 0) then
     write (*, '(A,I0,A)') "RESULT: ", failures, " failure(s)"
     error stop 1
  end if
  write (*, '(A)') "RESULT: all passed"
end program r4
