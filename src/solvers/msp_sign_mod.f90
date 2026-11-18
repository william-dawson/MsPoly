!> The matrix sign function and the polar decomposition.
module msp_sign_mod
  use msp_citations_mod, only: msp_cite
  use msp_kinds_mod, only: dp
  use msp_error_mod, only: msp_raise
  use msp_dict_mod, only: msp_dict, msp_args, msp_info, copy
  use msp_matrix_mod
  use msp_bounds_mod, only: bounds_impl, operator_norms
  use msp_polynomial_mod, only: identity_like
  use msp_convergence_mod, only: msp_monitor
  use msp_eigen_mod, only: dense_function_impl
  implicit none
  private

  public :: msp_sign, msp_polar, sign_impl

contains

  !> X becomes sign(A) for a nonsingular Hermitian A: the eigenvalues of A
  !> replaced by their signs. Collective. args and info as for sign_impl.
  subroutine msp_sign(A, X, args, info, context)
    type(msp_matrix), intent(in) :: A
    type(msp_matrix), intent(inout) :: X
    type(msp_dict), intent(in) :: args
    type(msp_dict), intent(inout), target :: info
    type(msp_dict), intent(in), target, optional :: context
    type(msp_dict), target :: tracked
    type(msp_args) :: a_args
    type(msp_info) :: r

    call copy(tracked, args)
    call a_args%open(tracked, context)
    call r%open(info)
    call sign_impl(A, a_args, r, X, .false.)
    call a_args%finish(r)
  end subroutine msp_sign

  !> The polar decomposition A = U H of a nonsingular A: U unitary and H
  !> Hermitian positive definite (returned when present). Collective. args
  !> and info as for sign_impl.
  subroutine msp_polar(A, U, args, info, context, H)
    type(msp_matrix), intent(in) :: A
    type(msp_matrix), intent(inout) :: U
    type(msp_dict), intent(in) :: args
    type(msp_dict), intent(inout), target :: info
    type(msp_dict), intent(in), target, optional :: context
    type(msp_matrix), intent(inout), optional :: H
    type(msp_dict), target :: tracked
    type(msp_args) :: a_args
    type(msp_info) :: r
    type(msp_matrix) :: uh
    real(dp) :: threshold

    call copy(tracked, args)
    call a_args%open(tracked, context)
    call r%open(info)
    call sign_impl(A, a_args, r, U, .true.)
    if (present(H)) then
       call a_args%get("threshold", threshold, default=0.0_dp, scoped=.true.)
       call msp_transpose(U, uh, conjugate=.true.)
       call msp_multiply(uh, A, H, threshold=threshold)
       call uh%free()
    end if
    call a_args%finish(r)
  end subroutine msp_polar

  !> Newton-Schulz iteration with Higham's dynamic scaling (Functions of
  !> Matrices, SIAM (2008), chapters 5 and 8), X <- X (3 I - a_k^2 X X) a_k
  !> / 2, with X^H X in place of X X for the polar factor. X starts as A
  !> divided by a bound on its spectral radius (sign: Gershgorin; polar:
  !> sqrt(|A|_1 |A|_inf), at least the largest singular value), so every
  !> eigen- or singular value x lies in (0, 1] in magnitude. Each step
  !> scales the iteration with a_k, the largest factor that still maps the
  !> spectral interval [x_k, 1] into itself, tracked by the scalar
  !> recursion x_{k+1} = a_k x_k (3 - a_k^2 x_k^2) / 2 from the ratio of
  !> the bounds, and capped at 1.69770248526, so badly spread spectra
  !> converge in a few steps where the unscaled iteration would crawl.
  !> The residual |I - X X|_F (or |I - X^H X|_F) need not fall
  !> monotonically under the scaling.
  !>
  !> It stops when the residual is at most the tolerance or when progress
  !> stalls: the residual over the last stagnation_window steps has fallen
  !> by less than stagnation_ratio. The window absorbs the scaling's
  !> non-monotonicity, which a strict non-improvement test would read as
  !> stagnation.
  !>
  !> With method "dense", the sign function of a Hermitian A comes from its
  !> eigendecomposition instead (args and info as for dense_function_impl).
  !> The polar factor then comes from the eigendecomposition of A^H A, U =
  !> A (A^H A)^(-1/2), so its accuracy is that of the squared condition
  !> number of A.
  !>
  !> args: method ("newton_schulz" (default) or "dense"), threshold (scoped,
  !> default 0), tolerance (scoped, default 1e-10), max_iterations (default
  !> 100), stagnation_ratio (default 0.6), stagnation_window (default 4).
  !> info: iterations, converged, stop, residual_history.
  subroutine sign_impl(A, args, info, X, polar)
    type(msp_matrix), intent(in) :: A
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    type(msp_matrix), intent(inout) :: X
    logical, intent(in) :: polar
    type(msp_info) :: bounds_info
    type(msp_matrix) :: xh, square, next, id
    type(msp_monitor) :: mon
    real(dp), dimension(:), allocatable :: residual
    real(dp) :: threshold, tolerance, emin, emax, scale
    real(dp) :: ratio, xk, alpha_k
    integer :: max_iterations, step, window
    character(:), allocatable :: reason, method

    call args%get("method", method, default="newton_schulz")
    if (method == "dense") then
       if (polar) then
          call dense_polar(A, X, args, info)
       else
          call dense_function_impl(A, X, sign_of, args, info)
       end if
       return
    else if (method /= "newton_schulz") then
       call msp_raise(merge("msp_polar", "msp_sign ", polar) // &
            & ": unknown method '" // method // "'")
    end if
    call msp_cite("higham2008functions", info)
    call args%get("threshold", threshold, default=0.0_dp, scoped=.true.)
    call args%get("tolerance", tolerance, default=1.0e-10_dp, scoped=.true.)
    call args%get("max_iterations", max_iterations, default=100)
    call args%get("stagnation_ratio", ratio, default=0.6_dp)
    call args%get("stagnation_window", window, default=4)
    if (polar) then
       call operator_norms(A, emin, emax)
       scale = sqrt(emin * emax)
    else
       call info%enter("bounds", bounds_info)
       call bounds_impl(A, bounds_info, emin, emax)
       scale = max(abs(emin), abs(emax))
    end if
    call identity_like(A, id)
    call mon%init(tolerance, ratio, window)
    call copy(X, A)
    call msp_scale(X, 1 / scale)
    !! The scalar the recursion tracks: the smallest magnitude in the
    !! scaled spectrum, from the same bounds as the scaling.
    xk = min(abs(emin), abs(emax)) / scale
    allocate(residual(0))
    reason = ""
    do step = 1, max_iterations
       if (polar) then
          call msp_transpose(X, xh, conjugate=.true.)
          call msp_multiply(xh, X, square, threshold=threshold)
       else
          call msp_multiply(X, X, square, threshold=threshold)
       end if
       call msp_scale(square, -1.0_dp)
       call msp_add(square, id)
       residual = [residual, msp_norm(square)]
       call info%append("residual_history", residual(step))
       reason = mon%check(residual(step))
       if (reason /= "") exit
       alpha_k = min(sqrt(3.0_dp / (1.0_dp + xk + xk * xk)), &
            & 1.69770248526_dp)
       call msp_scale(square, alpha_k**2)
       call msp_add(square, id, alpha=3.0_dp - alpha_k**2)
       call msp_multiply(X, square, next, alpha=0.5_dp * alpha_k, &
            & threshold=threshold)
       call msp_swap(X, next)
       xk = 0.5_dp * alpha_k * xk * (3.0_dp - alpha_k**2 * xk * xk)
    end do
    if (reason == "") reason = "max_iterations"
    call info%set("iterations", size(residual))
    call info%set("converged", reason /= "max_iterations")
    call info%set("stop", reason)
    call xh%free()
    call square%free()
    call next%free()
    call id%free()
  end subroutine sign_impl

  !> The unitary polar factor of A, U = A (A^H A)^(-1/2). Collective.
  subroutine dense_polar(A, U, args, info)
    type(msp_matrix), intent(in) :: A
    type(msp_matrix), intent(inout) :: U
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    type(msp_matrix) :: ah, ata, isq
    real(dp) :: threshold

    call args%get("threshold", threshold, default=0.0_dp, scoped=.true.)
    call msp_transpose(A, ah, conjugate=.true.)
    call msp_multiply(ah, A, ata, threshold=threshold)
    call dense_function_impl(ata, isq, inverse_root_of, args, info, "positive")
    call msp_multiply(A, isq, U, threshold=threshold)
    call ah%free()
    call ata%free()
    call isq%free()
  end subroutine dense_polar

  real(dp) function inverse_root_of(x)
    real(dp), intent(in) :: x

    inverse_root_of = 1 / sqrt(x)
  end function inverse_root_of

  !> -1, 0 or 1 by the sign of x.
  real(dp) function sign_of(x)
    real(dp), intent(in) :: x

    sign_of = merge(1.0_dp, 0.0_dp, x > 0) - merge(1.0_dp, 0.0_dp, x < 0)
  end function sign_of
end module msp_sign_mod
