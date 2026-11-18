!> Worked example of the calling convention.
!>
!> mock_density (public) -> density_impl -> bounds_impl -> power_impl.
!> The user sets "tolerance" at the top; neither density nor bounds knows
!> about it, but power_impl reads it with scoped lookup. Reports from every
!> level come back as one tree.
module mock_solvers_mod
  use msp_kinds_mod, only: dp
  use msp_dict_mod
  implicit none
  private

  public :: mock_density

contains

  !> Public routine: copy the arguments, open views, run, report unused keys.
  subroutine mock_density(args, info, context)
    type(msp_dict), intent(in) :: args
    type(msp_dict), intent(inout), target :: info
    type(msp_dict), intent(in), target, optional :: context
    type(msp_dict), target :: tracked
    type(msp_args) :: a
    type(msp_info) :: r

    call copy(tracked, args)
    call a%open(tracked, context)
    call r%open(info)
    call density_impl(a, r)
    call a%finish(r)
  end subroutine mock_density

  subroutine density_impl(args, info)
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    type(msp_args) :: inner_args
    type(msp_info) :: inner_info
    character(:), allocatable :: method
    real(dp) :: threshold
    integer :: k, max_iterations

    call args%get("method", method, default="trs2")
    call args%get("threshold", threshold, default=1.0e-8_dp, scoped=.true.)
    call args%get("max_iterations", max_iterations, default=100)

    call args%enter("bounds", inner_args)
    call info%enter("bounds", inner_info)
    call bounds_impl(inner_args, inner_info)

    do k = 1, min(6, max_iterations)
       call info%append("energy_history", -1.0_dp - 0.5_dp**k)
    end do
    call info%set("iterations", min(6, max_iterations))
    call info%set("method", method)
    call info%set("threshold_used", threshold)
    call info%set("converged", .true.)
  end subroutine density_impl

  subroutine bounds_impl(args, info)
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    type(msp_args) :: inner_args
    type(msp_info) :: inner_info
    real(dp) :: threshold
    integer :: max_iterations

    call args%get("threshold", threshold, default=1.0e-8_dp, scoped=.true.)
    call args%get("max_iterations", max_iterations, default=50)

    call args%enter("power", inner_args)
    call info%enter("power", inner_info)
    call power_impl(inner_args, inner_info, max_iterations)

    call info%set("threshold_used", threshold)
    call info%set("max_iterations_used", max_iterations)
  end subroutine bounds_impl

  subroutine power_impl(args, info, max_iterations)
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    integer, intent(in) :: max_iterations
    real(dp) :: tolerance, residual
    integer :: k

    call args%get("tolerance", tolerance, default=1.0e-6_dp, scoped=.true.)
    residual = 1.0_dp
    do k = 1, max_iterations
       residual = residual / 10
       call info%append("residual_history", residual)
       if (residual < tolerance) exit
    end do
    call info%set("tolerance_used", tolerance)
    call info%set("converged", residual < tolerance)
  end subroutine power_impl
end module mock_solvers_mod

program test_nested
  use msp_kinds_mod, only: dp
  use msp_dict_mod
  use mock_solvers_mod, only: mock_density
  use testing_mod, only: check, finish
  implicit none

  type(msp_dict) :: args
  type(msp_dict), target :: info, ctx
  real(dp), allocatable :: hist(:)
  character(:), allocatable :: unused(:)

  call ctx%set("threshold", 1.0e-5_dp)

  call args%set("tolerance", 2.0e-3_dp)
  call args%set("bounds.max_iterations", 20)
  call args%set("bounds.threshold", 1.0e-2_dp)
  call args%set("max_iteratoins", 1000)
  call mock_density(args, info, ctx)

  call check(info%get_real("bounds.power.tolerance_used") == 2.0e-3_dp, &
       & "top-level key reaches a routine two levels down")
  call check(info%get_real("threshold_used") == 1.0e-5_dp, &
       & "context default used by outer routine")
  call check(info%get_real("bounds.threshold_used") == 1.0e-2_dp, &
       & "section overrides for its branch only")
  call check(info%get_int("bounds.max_iterations_used") == 20, &
       & "local key read from its own section")
  call check(info%get_int("iterations") == 6, &
       & "outer max_iterations not set (typo), default used")
  call info%get("bounds.power.residual_history", hist)
  call check(size(hist) == 3, "history from the deepest routine")
  call info%get("unused_arguments", unused)
  call check(size(unused) == 1 .and. unused(1) == "max_iteratoins", &
       & "typo reported as unused")

  write (*, '(A)') "--- info ---"
  call info%write_yaml(6)
  call finish("test_nested")
end program test_nested
