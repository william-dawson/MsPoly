!> Minimal test helpers: count failed checks and fail the program at the end.
module testing_mod
  implicit none
  private

  public :: check, finish

  integer :: checks = 0
  integer :: failures = 0

contains

  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(*), intent(in) :: label

    checks = checks + 1
    if (.not. condition) then
       failures = failures + 1
       write (*, '(A)') "FAIL: " // label
    end if
  end subroutine check

  subroutine finish(name)
    character(*), intent(in) :: name

    write (*, '(A, I0, A, I0, A)') name // ": ", checks - failures, "/", &
         & checks, " checks passed"
    if (failures > 0) error stop 1
  end subroutine finish
end module testing_mod
