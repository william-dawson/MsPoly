!> Sines and cosines of Hermitian matrices.
module msp_trigonometry_mod
  use msp_citations_mod, only: msp_cite
  use msp_kinds_mod, only: dp
  use msp_error_mod, only: msp_raise
  use msp_dict_mod, only: msp_dict, msp_args, msp_info, copy
  use msp_matrix_mod
  use msp_bounds_mod, only: radius_bound
  use msp_polynomial_mod, only: monomial, identity_like
  use msp_exponential_mod, only: squarings
  use msp_eigen_mod, only: dense_function_impl
  implicit none
  private

  public :: msp_sin, msp_cos, trigonometry_impl

contains

  !> X becomes sin(A) for a Hermitian A. Collective. args and info as for
  !> trigonometry_impl.
  subroutine msp_sin(A, X, args, info, context)
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
    call trigonometry_impl(A, a_args, r, sine=X)
    call a_args%finish(r)
  end subroutine msp_sin

  !> X becomes cos(A) for a Hermitian A. Collective. args and info as for
  !> trigonometry_impl.
  subroutine msp_cos(A, X, args, info, context)
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
    call trigonometry_impl(A, a_args, r, cosine=X)
    call a_args%finish(r)
  end subroutine msp_cos

  !> Scaling and doubling: with Y = A / 2^s of spectral radius at most 1,
  !> C = cos(Y) and S = sin(Y) by their Taylor series in Y^2 to degree 20
  !> and 21 (Paterson-Stockmeyer), truncating below 1e-19; then s times
  !> S <- 2 S C and C <- 2 C^2 - I. The cosine alone needs no S. With
  !> method "dense", the functions of a Hermitian A come from its
  !> eigendecomposition instead (args and info as for dense_function_impl).
  !>
  !> args: method ("taylor" (default) or "dense"), threshold (scoped,
  !> default 0).
  !> info: doublings.
  subroutine trigonometry_impl(A, args, info, sine, cosine)
    type(msp_matrix), intent(in) :: A
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    type(msp_matrix), intent(inout), optional :: sine, cosine
    type(msp_matrix) :: y, y2, c, s, id, work
    real(dp), dimension(0:10) :: even, odd
    real(dp) :: threshold
    integer :: doublings, k
    character(:), allocatable :: method

    call args%get("method", method, default="taylor")
    select case (method)
    case ("dense")
       call dense_trigonometry(A, args, info, sine, cosine)
       return
    case ("taylor")
       call msp_cite("higham2003computing", info)
       call msp_cite("paterson1973number", info)
    case default
       call msp_raise("msp_trigonometry: unknown method '" // method // "'")
    end select
    call args%get("threshold", threshold, default=0.0_dp, scoped=.true.)
    doublings = squarings(radius_bound(A))
    call copy(y, A)
    call msp_scale(y, 0.5_dp**doublings)
    call msp_multiply(y, y, y2, threshold=threshold)
    even(0) = 1
    odd(0) = 1
    do k = 1, 10
       even(k) = -even(k - 1) / ((2 * k - 1) * (2 * k))
       odd(k) = -odd(k - 1) / ((2 * k) * (2 * k + 1))
    end do
    call monomial(y2, even, c, threshold, .true.)
    if (present(sine)) then
       call monomial(y2, odd, work, threshold, .true.)
       call msp_multiply(y, work, s, threshold=threshold)
    end if
    call identity_like(A, id)
    do k = 1, doublings
       if (present(sine)) then
          call msp_multiply(s, c, work, alpha=2.0_dp, threshold=threshold)
          call msp_swap(s, work)
       end if
       call msp_multiply(c, c, work, alpha=2.0_dp, threshold=threshold)
       call msp_add(work, id, alpha=-1.0_dp, threshold=threshold)
       call msp_swap(c, work)
    end do
    if (present(sine)) call msp_swap(sine, s)
    if (present(cosine)) call msp_swap(cosine, c)
    call info%set("doublings", doublings)
    call y%free()
    call y2%free()
    call c%free()
    call s%free()
    call id%free()
    call work%free()
  end subroutine trigonometry_impl

  !> sine and cosine of A, whichever are present, from one
  !> eigendecomposition.
  subroutine dense_trigonometry(A, args, info, sine, cosine)
    type(msp_matrix), intent(in) :: A
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    type(msp_matrix), intent(inout), optional :: sine, cosine

    if (present(sine) .and. present(cosine)) then
       call dense_function_impl(A, sine, sine_of, args, info, B2=cosine, &
            & f2=cosine_of)
    else if (present(sine)) then
       call dense_function_impl(A, sine, sine_of, args, info)
    else if (present(cosine)) then
       call dense_function_impl(A, cosine, cosine_of, args, info)
    end if
  end subroutine dense_trigonometry

  real(dp) function sine_of(x)
    real(dp), intent(in) :: x

    sine_of = sin(x)
  end function sine_of

  real(dp) function cosine_of(x)
    real(dp), intent(in) :: x

    cosine_of = cos(x)
  end function cosine_of
end module msp_trigonometry_mod
