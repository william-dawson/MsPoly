!> Exponentials and logarithms of Hermitian matrices.
module msp_exponential_mod
  use msp_citations_mod, only: msp_cite
  use msp_kinds_mod, only: dp
  use msp_error_mod, only: msp_raise
  use msp_dict_mod, only: msp_dict, msp_args, msp_info, copy
  use msp_matrix_mod
  use msp_bounds_mod, only: radius_bound
  use msp_polynomial_mod, only: monomial, chebyshev, identity_like
  use msp_sqrt_mod, only: sqrt_impl
  use msp_linear_mod, only: cg_impl
  use msp_eigen_mod, only: dense_function_impl
  implicit none
  private

  public :: msp_exp, msp_log, exp_impl, log_impl, squarings

contains

  !> X becomes exp(A) for a Hermitian A. Collective. args and info as for
  !> exp_impl.
  subroutine msp_exp(A, X, args, info, context)
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
    call exp_impl(A, a_args, r, X)
    call a_args%finish(r)
  end subroutine msp_exp

  !> X becomes log(A) for a Hermitian positive definite A. Collective. args
  !> and info as for log_impl.
  subroutine msp_log(A, X, args, info, context)
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
    call log_impl(A, a_args, r, X)
    call a_args%finish(r)
  end subroutine msp_log

  !> Scaling and squaring: with A / 2^s of spectral radius at most 1,
  !> exp(A) = exp(A / 2^s)^(2^s). exp on [-1, 1] is evaluated by its Taylor
  !> series to degree 18 (Paterson-Stockmeyer), its Chebyshev series to
  !> degree 17, with coefficients 2 I_k(1) from the modified Bessel
  !> functions, both truncating below 1e-16, or its [6/6] Pade approximant
  !> q(-A)^(-1) q(A), error about 1e-15 on [-1, 1], with the denominator,
  !> positive definite there, solved by conjugate gradients. With method
  !> "dense", exp of a Hermitian A comes from its eigendecomposition
  !> instead (args and info as for dense_function_impl).
  !>
  !> args: method ("taylor" (default), "chebyshev", "pade" or "dense"),
  !> threshold (scoped, default 0), cg (args of the conjugate gradients; see
  !> cg_impl).
  !> info: squarings, cg.
  subroutine exp_impl(A, args, info, X)
    type(msp_matrix), intent(in) :: A
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    type(msp_matrix), intent(inout) :: X
    type(msp_matrix) :: scaled, work
    type(msp_args) :: cg_args
    type(msp_info) :: cg_info
    character(:), allocatable :: method
    real(dp), dimension(0:18) :: c
    real(dp) :: threshold
    integer :: s, k, j

    call args%get("method", method, default="taylor")
    if (method == "dense") then
       call dense_function_impl(A, X, exponential, args, info)
       return
    end if
    call args%get("threshold", threshold, default=0.0_dp, scoped=.true.)
    s = squarings(radius_bound(A))
    call copy(scaled, A)
    call msp_scale(scaled, 0.5_dp**s)
    select case (method)
    case ("taylor")
       call msp_cite("paterson1973number", info)
       c(0) = 1
       do k = 1, 18
          c(k) = c(k - 1) / k
       end do
       call monomial(scaled, c, X, threshold, .true.)
    case ("chebyshev")
       c = 0
       do k = 0, 17
          do j = 0, 20
             c(k) = c(k) + 0.5_dp**(2 * j + k) / (gamma(j + 1.0_dp) * &
                  & gamma(j + k + 1.0_dp))
          end do
       end do
       c(1:) = 2 * c(1:)
       call chebyshev(scaled, c(0:17), X, threshold, .false.)
    case ("pade")
       call args%enter("cg", cg_args)
       call info%enter("cg", cg_info)
       call pade(scaled, threshold, cg_args, cg_info, X)
    case default
       call msp_raise("msp_exp: unknown method '" // method // "'")
    end select
    call msp_cite("higham2008functions", info)
    do k = 1, s
       call msp_multiply(X, X, work, threshold=threshold)
       call msp_swap(X, work)
    end do
    call info%set("squarings", s)
    call scaled%free()
    call work%free()
  end subroutine exp_impl

  !> X = q(-A)^(-1) q(A) for the [6/6] Pade approximant of exp, with
  !> q(x) = sum of c_k x^k, c_k = (12 - k)! 6! / (12! k! (6 - k)!), split
  !> into its even part E(A^2) and odd part A O(A^2): q(+-A) = E +- A O.
  subroutine pade(A, threshold, args, info, X)
    type(msp_matrix), intent(in) :: A
    real(dp), intent(in) :: threshold
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    type(msp_matrix), intent(inout) :: X
    type(msp_matrix) :: a2, even, odd, numerator, denominator
    real(dp), dimension(0:6) :: q
    integer :: k

    do k = 0, 6
       q(k) = gamma(13.0_dp - k) * gamma(7.0_dp) / (gamma(13.0_dp) * &
            & gamma(k + 1.0_dp) * gamma(7.0_dp - k))
    end do
    call msp_multiply(A, A, a2, threshold=threshold)
    call monomial(a2, q(0::2), even, threshold, .false.)
    call monomial(a2, q(1::2), odd, threshold, .false.)
    call msp_multiply(A, odd, numerator, threshold=threshold)
    call copy(denominator, even)
    call msp_add(denominator, numerator, alpha=-1.0_dp, threshold=threshold)
    call msp_add(numerator, even, threshold=threshold)
    call cg_impl(denominator, numerator, args, info, X)
    call a2%free()
    call even%free()
    call odd%free()
    call numerator%free()
    call denominator%free()
  end subroutine pade

  !> Inverse scaling and squaring: square roots bring A to X = A^(1/2^k)
  !> with every eigenvalue within 1/4 of 1, then
  !> log(A) = 2^k log(I + Y), Y = X - I, by the series of log(1 + y) to
  !> degree 26 (Paterson-Stockmeyer), which truncates below 1e-17 there.
  !> Each square root is accurate to its tolerance, which 2^k amplifies.
  !> With method "dense", log of a Hermitian positive definite A comes from
  !> its eigendecomposition instead (args and info as for
  !> dense_function_impl).
  !>
  !> args: method ("inverse_scaling" (default) or "dense"), threshold
  !> (scoped, default 0), max_square_roots (default 64;
  !> a positive definite matrix in double precision never needs as many),
  !> sqrt (args of the square roots; see sqrt_impl).
  !> info: square_roots, sqrt (info of the last square root).
  subroutine log_impl(A, args, info, X)
    type(msp_matrix), intent(in) :: A
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    type(msp_matrix), intent(inout) :: X
    type(msp_args) :: sqrt_args
    type(msp_info) :: sqrt_info
    type(msp_matrix) :: root, y, id
    real(dp), dimension(0:26) :: c
    real(dp) :: threshold
    integer :: k, roots, max_roots
    character(:), allocatable :: method

    call args%get("method", method, default="inverse_scaling")
    select case (method)
    case ("dense")
       call dense_function_impl(A, X, logarithm, args, info, "positive")
       return
    case ("inverse_scaling")
       call msp_cite("higham2008functions", info)
       call msp_cite("paterson1973number", info)
    case default
       call msp_raise("msp_log: unknown method '" // method // "'")
    end select
    call args%get("threshold", threshold, default=0.0_dp, scoped=.true.)
    call args%get("max_square_roots", max_roots, default=64)
    call args%enter("sqrt", sqrt_args)
    call info%enter("sqrt", sqrt_info)
    call identity_like(A, id)
    call copy(root, A)
    roots = 0
    do
       call copy(y, root)
       call msp_add(y, id, alpha=-1.0_dp)
       if (radius_bound(y) <= 0.25_dp) exit
       if (roots == max_roots) call msp_raise("msp_log: more than " // &
            & "max_square_roots square roots needed; is the matrix " // &
            & "positive definite?")
       call sqrt_impl(root, sqrt_args, sqrt_info, root=X)
       call msp_swap(root, X)
       roots = roots + 1
    end do
    c(0) = 0
    do k = 1, 26
       c(k) = (-1)**(k + 1) / real(k, dp)
    end do
    call monomial(y, c, X, threshold, .true.)
    call msp_scale(X, 2.0_dp**roots)
    call info%set("square_roots", roots)
    call root%free()
    call y%free()
    call id%free()
  end subroutine log_impl

  !> The smallest s >= 0 with r / 2^s <= 1.
  integer function squarings(r) result(s)
    real(dp), intent(in) :: r

    s = 0
    do while (r / 2.0_dp**s > 1)
       s = s + 1
    end do
  end function squarings

  real(dp) function exponential(x)
    real(dp), intent(in) :: x

    exponential = exp(x)
  end function exponential

  real(dp) function logarithm(x)
    real(dp), intent(in) :: x

    logarithm = log(x)
  end function logarithm
end module msp_exponential_mod
