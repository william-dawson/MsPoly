!> Square roots and inverse square roots of Hermitian positive definite
!> matrices.
module msp_sqrt_mod
  use msp_citations_mod, only: msp_cite
  use msp_kinds_mod, only: dp
  use msp_error_mod, only: msp_raise
  use msp_dict_mod, only: msp_dict, msp_args, msp_info, copy
  use msp_matrix_mod
  use msp_bounds_mod, only: bounds_impl
  use msp_convergence_mod, only: msp_monitor
  use msp_eigen_mod, only: dense_function_impl
  implicit none
  private

  public :: msp_sqrt, msp_inverse_sqrt, sqrt_impl

contains

  !> X becomes the square root of the Hermitian positive definite matrix S.
  !> Collective. args and info as for sqrt_impl.
  subroutine msp_sqrt(S, X, args, info, context)
    type(msp_matrix), intent(in) :: S
    type(msp_matrix), intent(inout) :: X
    type(msp_dict), intent(in) :: args
    type(msp_dict), intent(inout), target :: info
    type(msp_dict), intent(in), target, optional :: context
    type(msp_dict), target :: tracked
    type(msp_args) :: a
    type(msp_info) :: r

    call copy(tracked, args)
    call a%open(tracked, context)
    call r%open(info)
    call sqrt_impl(S, a, r, root=X)
    call a%finish(r)
  end subroutine msp_sqrt

  !> X becomes the inverse square root of the Hermitian positive definite
  !> matrix S. Collective. args and info as for sqrt_impl.
  subroutine msp_inverse_sqrt(S, X, args, info, context)
    type(msp_matrix), intent(in) :: S
    type(msp_matrix), intent(inout) :: X
    type(msp_dict), intent(in) :: args
    type(msp_dict), intent(inout), target :: info
    type(msp_dict), intent(in), target, optional :: context
    type(msp_dict), target :: tracked
    type(msp_args) :: a
    type(msp_info) :: r

    call copy(tracked, args)
    call a%open(tracked, context)
    call r%open(info)
    call sqrt_impl(S, a, r, inverse_root=X)
    call a%finish(r)
  end subroutine msp_inverse_sqrt

  !> Newton-Schulz iteration. With S scaled by its Gershgorin bound so its
  !> spectrum lies in (0, 1], Y = S and Z = I, and only multiplies are
  !> used, Y drives to the square root and Z to the inverse square root.
  !> The default order 5 follows Jansik et al., J. Chem. Phys. 126, 124104
  !> (2007): from X = ZY - I, Z = p(X) Z and Y = Y p(X), with p the Taylor
  !> polynomial 1 - x/2 + 3x^2/8 - 5x^3/16 + 35x^4/128 of (1 + x)^(-1/2),
  !> so the error contracts with the fifth power at every step. Order 3
  !> keeps the terms through x^2. Order 2 is the coupled iteration
  !> (Higham, Numer. Algorithms 15, 227 (1997)): T = (3I - YZ) / 2,
  !> Y = TY, Z = ZT, converging quadratically.
  !>
  !> It stops when the residual |I - YZ| (Frobenius) is at most the
  !> tolerance (a negative tolerance turns this test off), or when
  !> progress toward the truncation floor stalls: the residual over the
  !> last stagnation_window steps has fallen by less than
  !> stagnation_ratio. The filtered iteration reaches its floor's
  !> neighbourhood in a few steps and then crawls, so the wide-window
  !> test ends it long before a strict non-improvement test would.
  !>
  !> With method "dense", the roots come from the eigendecomposition of S
  !> instead (args and info as for dense_function_impl).
  !>
  !> args: method ("newton_schulz" (default) or "dense"), order (2, 3 or
  !> 5; default 5), threshold (scoped, default 0), tolerance (scoped,
  !> default 1e-10), max_iterations (default 100), stagnation_ratio
  !> (default 0.6), stagnation_window (default 4), bounds (args of the
  !> bounds solver: method "gershgorin" (default) or "power"; power's
  !> estimate is safe here because the final rescaling corrects the
  !> scaling). Work matrices follow the precision of S.
  !> info: iterations, converged, stop ("tolerance", "stagnation" or
  !> "max_iterations"), residual_history, bounds.
  subroutine sqrt_impl(S, args, info, root, inverse_root)
    type(msp_matrix), intent(in) :: S
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    type(msp_matrix), intent(inout), optional :: root, inverse_root
    type(msp_info) :: bounds_info
    type(msp_args) :: bounds_args
    type(msp_monitor) :: mon
    type(msp_matrix) :: y, z, t, xk, id, work
    real(dp) :: residual
    real(dp) :: threshold, tolerance, emin, emax
    real(dp) :: aa, bb, cc, dd, pa, pb, pc, pd
    real(dp) :: ratio
    integer :: max_iterations, order, step, iterations, window
    character(:), allocatable :: reason, method

    call args%get("method", method, default="newton_schulz")
    select case (method)
    case ("dense")
       if (present(root) .and. present(inverse_root)) then
          call dense_function_impl(S, root, square_root, args, info, &
               & "positive", inverse_root, inverse_square_root)
       else if (present(root)) then
          call dense_function_impl(S, root, square_root, args, info, &
               & "nonnegative")
       else if (present(inverse_root)) then
          call dense_function_impl(S, inverse_root, inverse_square_root, &
               & args, info, "positive")
       end if
       return
    case ("newton_schulz")
       call msp_cite("jansik2007linear", info)
    case default
       call msp_raise("msp_sqrt: unknown method '" // method // "'")
    end select
    call args%get("threshold", threshold, default=0.0_dp, scoped=.true.)
    call args%get("tolerance", tolerance, default=1.0e-10_dp, scoped=.true.)
    call args%get("max_iterations", max_iterations, default=100)
    call args%get("order", order, default=5)
    if (order /= 2 .and. order /= 3 .and. order /= 5) &
         & call msp_raise("msp_sqrt: order must be 2, 3 or 5")
    call args%get("stagnation_ratio", ratio, default=0.6_dp)
    call args%get("stagnation_window", window, default=4)
    call mon%init(tolerance, ratio, window)
    call args%enter("bounds", bounds_args)
    call info%enter("bounds", bounds_info)
    call bounds_impl(S, bounds_info, emin, emax, bounds_args)

    call id%init(S%grid, S%n)
    call id%identity()
    if (S%is_single()) call msp_to_single(id)
    if (S%on_device()) call msp_to_device(id)
    call copy(y, S)
    call msp_scale(y, 1 / emax)
    call copy(z, id)
    reason = ""
    iterations = 0
    select case (order)
    case (2)
       do step = 1, max_iterations
          call msp_multiply(y, z, xk, threshold=threshold)
          call copy(t, id)
          call msp_scale(t, 3.0_dp)
          call msp_add(t, xk, alpha=-1.0_dp)
          call msp_scale(t, 0.5_dp)
          call msp_add(xk, id, alpha=-1.0_dp)
          residual = msp_norm(xk)
          reason = mon%check(residual)
          call info%append("residual_history", residual)
          if (reason /= "") exit
          call msp_multiply(t, y, work, threshold=threshold)
          call msp_swap(y, work)
          call msp_multiply(z, t, work, threshold=threshold)
          call msp_swap(z, work)
          iterations = step
       end do
    case default
       !! Taylor polynomials of (1 + x)^(-1/2): from X = Z Y - I, apply
       !! Z = p(X) Z and Y = Y p(X), so the error contracts with the
       !! polynomial's order at every step.
       do step = 1, max_iterations
          call msp_multiply(z, y, xk, threshold=threshold)
          call msp_add(xk, id, alpha=-1.0_dp)
          residual = msp_norm(xk)
          reason = mon%check(residual)
          call info%append("residual_history", residual)
          if (reason /= "") exit
          select case (order)
          case (3)
             !! p(x) = 1 - x/2 + 3x^2/8.
             call msp_multiply(xk, xk, t, threshold=threshold)
             call msp_scale(xk, -0.5_dp)
             call msp_add(xk, id)
             call msp_add(xk, t, alpha=0.375_dp)
          case (5)
             !! p(x) = 1 - x/2 + 3x^2/8 - 5x^3/16 + 35x^4/128, scaled to
             !! x^4 + aa x^3 + bb x^2 + cc x + dd and evaluated by Knuth's
             !! two-product form p = (z + x + pb)(z + pc) + pd with
             !! z = x(x + pa).
             aa = -40.0_dp / 35.0_dp
             bb = 48.0_dp / 35.0_dp
             cc = -64.0_dp / 35.0_dp
             dd = 128.0_dp / 35.0_dp
             pa = (aa - 1.0_dp) / 2.0_dp
             pb = bb * (pa + 1.0_dp) - cc - pa * (pa + 1.0_dp)**2
             pc = bb - pb - pa * (pa + 1.0_dp)
             pd = dd - pb * pc
             call msp_multiply(xk, xk, t, threshold=threshold)
             call msp_add(t, xk, alpha=pa)
             call copy(work, id)
             call msp_scale(work, pb)
             call msp_add(work, xk)
             call msp_add(work, t)
             call msp_add(t, id, alpha=pc)
             call msp_multiply(work, t, xk, threshold=threshold)
             call msp_add(xk, id, alpha=pd)
             call msp_scale(xk, 35.0_dp / 128.0_dp)
          end select
          call msp_multiply(xk, z, work, threshold=threshold)
          call msp_swap(z, work)
          call msp_multiply(y, xk, work, threshold=threshold)
          call msp_swap(y, work)
          iterations = step
       end do
    end select
    if (reason == "") reason = "max_iterations"
    if (present(root)) then
       call msp_swap(root, y)
       call msp_scale(root, sqrt(emax))
    end if
    if (present(inverse_root)) then
       call msp_swap(inverse_root, z)
       call msp_scale(inverse_root, 1 / sqrt(emax))
    end if
    call info%set("iterations", iterations)
    call info%set("converged", reason /= "max_iterations")
    call info%set("stop", reason)
    call y%free()
    call z%free()
    call t%free()
    call xk%free()
    call id%free()
    call work%free()
  end subroutine sqrt_impl

  real(dp) function square_root(x)
    real(dp), intent(in) :: x

    square_root = sqrt(max(x, 0.0_dp))
  end function square_root

  real(dp) function inverse_square_root(x)
    real(dp), intent(in) :: x

    inverse_square_root = 1 / sqrt(x)
  end function inverse_square_root
end module msp_sqrt_mod
