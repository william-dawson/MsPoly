!> Density matrix guesses for a new geometry, from the density matrix of
!> the previous one.
module msp_extrapolate_mod
  use msp_citations_mod, only: msp_cite
  use msp_kinds_mod, only: dp
  use msp_dict_mod, only: msp_dict, msp_args, msp_info, copy
  use msp_matrix_mod
  use msp_sqrt_mod, only: sqrt_impl
  implicit none
  private

  public :: msp_purification_extrapolate, msp_lowdin_extrapolate

contains

  !> D becomes a density matrix for the overlap S, made from the previous
  !> density matrix P by purification in the metric S: D <- D S D when
  !> tr(D S) exceeds args "trace" and D <- 2 D - D S D otherwise, which
  !> drives D S to a projector of that trace. Collective.
  !>
  !> args: trace (required), threshold (scoped, default 0), tolerance
  !> (scoped, default 1e-10: stop when a step changes D by at most this),
  !> max_iterations (default 100).
  !> info: iterations, converged, residual_history.
  subroutine msp_purification_extrapolate(P, S, D, args, info, context)
    type(msp_matrix), intent(in) :: P, S
    type(msp_matrix), intent(inout) :: D
    type(msp_dict), intent(in) :: args
    type(msp_dict), intent(inout), target :: info
    type(msp_dict), intent(in), target, optional :: context
    type(msp_dict), target :: tracked
    type(msp_args) :: a
    type(msp_info) :: r
    type(msp_matrix) :: ds, dsd
    real(dp) :: target, threshold, tolerance, trace, change
    integer :: max_iterations, step, iterations
    logical :: converged

    call copy(tracked, args)
    call a%open(tracked, context)
    call r%open(info)
    call a%get("trace", target)
    call a%get("threshold", threshold, default=0.0_dp, scoped=.true.)
    call a%get("tolerance", tolerance, default=1.0e-10_dp, scoped=.true.)
    call a%get("max_iterations", max_iterations, default=100)
    call msp_cite("niklasson2010trace", r)
    call copy(D, P)
    converged = .false.
    iterations = 0
    do step = 1, max_iterations
       call msp_multiply(D, S, ds, threshold=threshold)
       call msp_multiply(ds, D, dsd, threshold=threshold)
       call msp_dot(D, S, trace)
       if (trace < target) then
          call msp_scale(dsd, -1.0_dp)
          call msp_add(dsd, D, alpha=2.0_dp, threshold=threshold)
       end if
       call msp_add(D, dsd, alpha=-1.0_dp)
       change = msp_norm(D)
       call msp_swap(D, dsd)
       call r%append("residual_history", change)
       iterations = step
       converged = change <= tolerance
       if (converged) exit
    end do
    call r%set("iterations", iterations)
    call r%set("converged", converged)
    call ds%free()
    call dsd%free()
    call a%finish(r)
  end subroutine msp_purification_extrapolate

  !> D becomes the Lowdin extrapolation of the previous density matrix P
  !> from the overlap S_old to S_new:
  !> D = S_new^(-1/2) S_old^(1/2) P S_old^(1/2) S_new^(-1/2), which keeps the
  !> density matrix in the Lowdin orthogonal basis. Collective.
  !>
  !> args: threshold (scoped, default 0), sqrt (args of both square roots;
  !> see sqrt_impl).
  !> info: sqrt (info of the inverse square root of S_new).
  subroutine msp_lowdin_extrapolate(P, S_old, S_new, D, args, info, context)
    type(msp_matrix), intent(in) :: P, S_old, S_new
    type(msp_matrix), intent(inout) :: D
    type(msp_dict), intent(in) :: args
    type(msp_dict), intent(inout), target :: info
    type(msp_dict), intent(in), target, optional :: context
    type(msp_dict), target :: tracked
    type(msp_args) :: a, sqrt_args
    type(msp_info) :: r, sqrt_info
    type(msp_matrix) :: root, inverse_root, work
    real(dp) :: threshold

    call copy(tracked, args)
    call a%open(tracked, context)
    call r%open(info)
    call a%get("threshold", threshold, default=0.0_dp, scoped=.true.)
    call msp_cite("exner2002comparison", r)
    call a%enter("sqrt", sqrt_args)
    call r%enter("sqrt", sqrt_info)
    call sqrt_impl(S_old, sqrt_args, sqrt_info, root=root)
    call sqrt_impl(S_new, sqrt_args, sqrt_info, inverse_root=inverse_root)
    call msp_multiply(inverse_root, root, work, threshold=threshold)
    call msp_multiply(work, P, D, threshold=threshold)
    call msp_multiply(D, root, work, threshold=threshold)
    call msp_multiply(work, inverse_root, D, threshold=threshold)
    call root%free()
    call inverse_root%free()
    call work%free()
    call a%finish(r)
  end subroutine msp_lowdin_extrapolate
end module msp_extrapolate_mod
