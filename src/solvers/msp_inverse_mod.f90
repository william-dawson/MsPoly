!> Inverses, pseudo-inverses and p-th roots of matrices.
module msp_inverse_mod
  use msp_citations_mod, only: msp_cite
  use msp_kinds_mod, only: dp
  use msp_error_mod, only: msp_raise
  use msp_dict_mod, only: msp_dict, msp_args, msp_info, copy
  use msp_matrix_mod
  use msp_bounds_mod, only: bounds_impl, operator_norms
  use msp_sqrt_mod, only: sqrt_impl
  use msp_polynomial_mod, only: identity_like
  use msp_convergence_mod, only: msp_monitor
  use msp_eigen_mod, only: dense_function_impl, eigen_impl, eigen_compose
  implicit none
  private

  public :: msp_inverse, msp_pseudo_inverse, msp_root, msp_inverse_root, &
       & inverse_impl, root_impl

contains

  !> X becomes the inverse of the nonsingular matrix A. Collective. args and
  !> info as for inverse_impl.
  subroutine msp_inverse(A, X, args, info, context)
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
    call inverse_impl(A, a_args, r, X)
    call a_args%finish(r)
  end subroutine msp_inverse

  !> X becomes the Moore-Penrose pseudo-inverse of A, by the same iteration
  !> as msp_inverse. Collective. args and info as for inverse_impl.
  subroutine msp_pseudo_inverse(A, X, args, info, context)
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
    call inverse_impl(A, a_args, r, X, pseudo=.true.)
    call a_args%finish(r)
  end subroutine msp_pseudo_inverse

  !> X becomes A^(1/p), for a Hermitian positive definite A and the integer
  !> p >= 1 given as args "root". Collective. Other args and info as for
  !> root_impl.
  subroutine msp_root(A, X, args, info, context)
    type(msp_matrix), intent(in) :: A
    type(msp_matrix), intent(inout) :: X
    type(msp_dict), intent(in) :: args
    type(msp_dict), intent(inout), target :: info
    type(msp_dict), intent(in), target, optional :: context

    call root_entry(A, X, args, info, .false., context)
  end subroutine msp_root

  !> X becomes A^(-1/p), for a Hermitian positive definite A and the
  !> integer p >= 1 given as args "root". Collective. Other args and info as
  !> for root_impl.
  subroutine msp_inverse_root(A, X, args, info, context)
    type(msp_matrix), intent(in) :: A
    type(msp_matrix), intent(inout) :: X
    type(msp_dict), intent(in) :: args
    type(msp_dict), intent(inout), target :: info
    type(msp_dict), intent(in), target, optional :: context

    call root_entry(A, X, args, info, .true., context)
  end subroutine msp_inverse_root

  subroutine root_entry(A, X, args, info, inverse, context)
    type(msp_matrix), intent(in) :: A
    type(msp_matrix), intent(inout) :: X
    type(msp_dict), intent(in) :: args
    type(msp_dict), intent(inout), target :: info
    logical, intent(in) :: inverse
    type(msp_dict), intent(in), target, optional :: context
    type(msp_dict), target :: tracked
    type(msp_args) :: a_args
    type(msp_info) :: r
    integer :: p

    call copy(tracked, args)
    call a_args%open(tracked, context)
    call r%open(info)
    call a_args%get("root", p)
    if (p < 1) call msp_raise("msp_root: root must be at least 1")
    call root_impl(A, p, inverse, a_args, r, X)
    call a_args%finish(r)
  end subroutine root_entry

  !> Newton-Schulz (Hotelling) iteration X <- 2X - X A X from
  !> X = A^H / (|A|_1 |A|_inf). Then X A = A^H A / (|A|_1 |A|_inf) is
  !> Hermitian with its spectrum in [0, 1], so the iteration converges
  !> quadratically to the inverse, or, for a singular A, to the
  !> pseudo-inverse, since every iterate keeps the range of A^H (Ben-Israel
  !> and Cohen, SIAM J. Numer. Anal. 3, 410 (1966)).
  !>
  !> It stops when |I - X A|_F is at most the tolerance or no longer falls.
  !> The eigenvalues mu of X A stay in [0, 1] and map to 2 mu - mu^2, so the
  !> residual falls monotonically: to 0 for an inverse, and for a
  !> pseudo-inverse to the square root of the nullity, where it stalls.
  !>
  !> With method "dense", the inverse of a Hermitian A comes from its
  !> eigendecomposition instead (args and info as for dense_function_impl).
  !> For the pseudo-inverse, eigenvalues of magnitude at most cutoff (default
  !> n epsilon) times the largest are taken to be zero, and info reports the
  !> rank.
  !>
  !> args: method ("newton_schulz" (default) or "dense"), threshold (scoped,
  !> default 0), tolerance (scoped, default 1e-10), max_iterations (default
  !> 100), cutoff (dense pseudo-inverse).
  !> info: iterations, converged, stop, residual_history; rank.
  subroutine inverse_impl(A, args, info, X, pseudo)
    type(msp_matrix), intent(in) :: A
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    type(msp_matrix), intent(inout) :: X
    logical, intent(in), optional :: pseudo
    type(msp_matrix) :: xa, next, id
    type(msp_monitor) :: mon
    real(dp), dimension(:), allocatable :: residual
    real(dp) :: threshold, tolerance, one, infinity
    integer :: max_iterations, step
    character(:), allocatable :: reason, method

    call args%get("method", method, default="newton_schulz")
    select case (method)
    case ("dense")
       if (present(pseudo)) then
          if (pseudo) then
             call dense_pseudo_inverse(A, X, args, info)
             return
          end if
       end if
       call dense_function_impl(A, X, reciprocal, args, info, "nonzero")
       return
    case ("newton_schulz")
       call msp_cite("palser1998canonical", info)
    case default
       call msp_raise("msp_inverse: unknown method '" // method // "'")
    end select
    call args%get("threshold", threshold, default=0.0_dp, scoped=.true.)
    call args%get("tolerance", tolerance, default=1.0e-10_dp, scoped=.true.)
    call args%get("max_iterations", max_iterations, default=100)
    call operator_norms(A, one, infinity)
    call msp_transpose(A, X, conjugate=.true.)
    if (one > 0) call msp_scale(X, 1 / (one * infinity))
    call identity_like(A, id)
    call mon%init(tolerance)
    allocate(residual(0))
    reason = ""
    do step = 1, max_iterations
       call msp_multiply(X, A, xa, threshold=threshold)
       call copy(next, id)
       call msp_add(next, xa, alpha=-1.0_dp)
       residual = [residual, msp_norm(next)]
       call info%append("residual_history", residual(step))
       reason = mon%check(residual(step))
       if (reason /= "") exit
       call msp_multiply(xa, X, next, alpha=-1.0_dp, threshold=threshold)
       call msp_add(next, X, alpha=2.0_dp, threshold=threshold)
       call msp_swap(X, next)
    end do
    if (reason == "") reason = "max_iterations"
    call info%set("iterations", size(residual))
    call info%set("converged", reason /= "max_iterations")
    call info%set("stop", reason)
    call xa%free()
    call next%free()
    call id%free()
  end subroutine inverse_impl

  !> X = A^(1/p), or A^(-1/p) with inverse. Square roots come from the
  !> coupled Newton-Schulz iteration (sqrt_impl), the first power from
  !> inverse_impl; other roots from the coupled Newton iteration (Guo and
  !> Higham, SIAM J. Matrix Anal. Appl. 28, 788 (2006)): with c the
  !> Gershgorin upper bound, M = A / c and X = c^(-1/p) I, each step
  !> T = ((p + 1) I - M) / p, X <- X T, M <- T^p M keeps M = X^p A, and as
  !> the spectrum of M lies in (0, 1], M converges to I and X to A^(-1/p).
  !> A^(1/p) is then A X^(p - 1). Even orders reduce through square roots
  !> first: A^(+-1/p) for even p is the half-order root of A^(1/2), and a
  !> square-root iteration costs fewer multiplies than the coupled step's
  !> chain of T powers.
  !>
  !> It stops when |M - I|_F is at most the tolerance or no longer falls.
  !>
  !> With method "dense", every order comes from the eigendecomposition of
  !> the Hermitian A instead (args and info as for dense_function_impl).
  !>
  !> args: method ("newton_schulz" (default) or "dense"), threshold (scoped,
  !> default 0), tolerance (scoped, default 1e-10), max_iterations (default
  !> 100), bounds (args of the bounds solver:
  !> method "gershgorin" (default) or "power"; power's estimate is safe
  !> here because the final rescaling corrects the scaling). The reduction
  !> square roots take their args from the same dictionary; their reports
  !> land under square_root.
  !> info: iterations, converged, stop, residual_history, square_roots (the
  !> number of reduction square roots, when any ran), square_root (args and
  !> info of the last one).
  subroutine root_impl(A, p, inverse, args, info, X)
    type(msp_matrix), intent(in) :: A
    integer, intent(in) :: p
    logical, intent(in) :: inverse
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    type(msp_matrix), intent(inout) :: X
    type(msp_matrix) :: m, t, tp, id, work, base
    type(msp_monitor) :: mon
    type(msp_info) :: bounds_info, inner_info
    type(msp_args) :: bounds_args
    real(dp), dimension(:), allocatable :: residual
    real(dp) :: threshold, tolerance, emin, emax
    integer :: max_iterations, step, k, power, reductions
    character(:), allocatable :: reason, method

    call args%get("method", method, default="newton_schulz")
    if (method == "dense") then
       call dense_function_impl(A, X, root_of, args, info, &
            & trim(merge("positive   ", "nonnegative", inverse)))
       return
    else if (method /= "newton_schulz") then
       call msp_raise("msp_root: unknown method '" // method // "'")
    end if
    if (p == 1 .and. .not. inverse) then
       call copy(X, A)
       return
    else if (p == 1) then
       call inverse_impl(A, args, info, X)
       return
    else if (p == 2 .and. inverse) then
       call sqrt_impl(A, args, info, inverse_root=X)
       return
    else if (p == 2) then
       call sqrt_impl(A, args, info, root=X)
       return
    end if
    call args%get("threshold", threshold, default=0.0_dp, scoped=.true.)
    call args%get("tolerance", tolerance, default=1.0e-10_dp, scoped=.true.)
    call args%get("max_iterations", max_iterations, default=100)
    call msp_cite("higham2008functions", info)
    power = p
    reductions = 0
    call copy(m, A)
    do while (power > 2 .and. mod(power, 2) == 0)
       call info%enter("square_root", inner_info)
       call sqrt_impl(m, args, inner_info, root=work)
       call msp_swap(m, work)
       call work%free()
       power = power / 2
       reductions = reductions + 1
    end do
    if (reductions > 0) call info%set("square_roots", reductions)
    call args%enter("bounds", bounds_args)
    call info%enter("bounds", bounds_info)
    call bounds_impl(m, bounds_info, emin, emax, bounds_args)
    call identity_like(A, id)
    !! The coupled iteration scales m in place; the unscaled reduced root
    !! builds the positive-order result at the end.
    call copy(base, m)
    call msp_scale(m, 1 / emax)
    call copy(X, id)
    call msp_scale(X, emax**(-1.0_dp / power))
    call mon%init(tolerance)
    allocate(residual(0))
    reason = ""
    do step = 1, max_iterations
       call copy(t, id)
       call msp_scale(t, real(power + 1, dp))
       call msp_add(t, m, alpha=-1.0_dp)
       call msp_scale(t, 1.0_dp / power)
       call msp_multiply(X, t, work, threshold=threshold)
       call msp_swap(X, work)
       call copy(tp, t)
       do k = 2, power
          call msp_multiply(tp, t, work, threshold=threshold)
          call msp_swap(tp, work)
       end do
       call msp_multiply(tp, m, work, threshold=threshold)
       call msp_swap(m, work)
       call copy(work, m)
       call msp_add(work, id, alpha=-1.0_dp)
       residual = [residual, msp_norm(work)]
       call info%append("residual_history", residual(step))
       reason = mon%check(residual(step))
       if (reason /= "") exit
    end do
    if (reason == "") reason = "max_iterations"
    if (.not. inverse) then
       call copy(t, base)
       do k = 1, power - 1
          call msp_multiply(t, X, work, threshold=threshold)
          call msp_swap(t, work)
       end do
       call msp_swap(X, t)
    end if
    call info%set("iterations", size(residual))
    call info%set("converged", reason /= "max_iterations")
    call info%set("stop", reason)
    call m%free()
    call t%free()
    call tp%free()
    call id%free()
    call work%free()
    call base%free()

  contains

    !> x^(1/p), or x^(-1/p) with inverse; eigenvalues that rounding left
    !> just below zero count as zero.
    real(dp) function root_of(x)
      real(dp), intent(in) :: x

      root_of = max(x, 0.0_dp)**(merge(-1.0_dp, 1.0_dp, inverse) / p)
    end function root_of
  end subroutine root_impl

  !> X = A^+ for a Hermitian A, from its eigendecomposition. Collective.
  subroutine dense_pseudo_inverse(A, X, args, info)
    type(msp_matrix), intent(in) :: A
    type(msp_matrix), intent(inout) :: X
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    type(msp_matrix) :: v
    real(dp), dimension(:), allocatable :: values, weights
    real(dp) :: cutoff, threshold, limit
    integer :: q

    call msp_cite("higham2008functions", info)
    call args%get("cutoff", cutoff, default=A%n * epsilon(1.0_dp))
    call args%get("threshold", threshold, default=0.0_dp, scoped=.true.)
    call eigen_impl(A, values, args, info, vectors=v)
    limit = cutoff * maxval(abs(values))
    allocate(weights(size(values)))
    do q = 1, size(values)
       weights(q) = 0
       if (abs(values(q)) > limit) weights(q) = 1 / values(q)
    end do
    call eigen_compose(v, weights, X, threshold)
    call info%set("rank", count(weights /= 0))
    call v%free()
  end subroutine dense_pseudo_inverse

  real(dp) function reciprocal(x)
    real(dp), intent(in) :: x

    reciprocal = 1 / x
  end function reciprocal
end module msp_inverse_mod
