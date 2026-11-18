!> Error reporting convention.
!>
!> Procedures that can fail take optional `stat` and `errmsg` arguments, like
!> the intrinsic ALLOCATE. If `stat` is present the error is returned through
!> it; otherwise execution stops with the message.
module msp_error_mod
  use, intrinsic :: iso_fortran_env, only: error_unit
  implicit none
  private

  public :: msp_raise, msp_clear

contains

  !> Report an error through stat/errmsg, or stop if stat is absent.
  subroutine msp_raise(msg, stat, errmsg)
    !> Description of the error.
    character(*), intent(in) :: msg
    !> Set to 1 if present.
    integer, intent(out), optional :: stat
    !> Set to msg if present.
    character(:), allocatable, intent(out), optional :: errmsg

    if (present(stat)) then
       stat = 1
       if (present(errmsg)) errmsg = msg
    else
       write (error_unit, '(A)') msg
       error stop 1
    end if
  end subroutine msp_raise

  !> Mark success on entry to a procedure that can fail.
  subroutine msp_clear(stat)
    !> Set to 0 if present.
    integer, intent(out), optional :: stat

    if (present(stat)) stat = 0
  end subroutine msp_clear
end module msp_error_mod
