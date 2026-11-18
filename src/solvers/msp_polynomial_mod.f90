!> Polynomials of matrices in the monomial, Chebyshev and Hermite bases.
module msp_polynomial_mod
  use msp_citations_mod, only: msp_cite
  use msp_kinds_mod, only: dp
  use msp_error_mod, only: msp_raise
  use msp_dict_mod, only: msp_dict, msp_args, msp_info, copy
  use msp_matrix_mod
  implicit none
  private

  public :: msp_polynomial, polynomial_impl, monomial, chebyshev, identity_like

contains

  !> B = p(A) = sum of c(k + 1) P_k(A) over k = 0, ..., degree, with P_k the
  !> basis polynomials. Collective.
  !>
  !> args: coefficients (required, real array c), basis ("monomial"
  !> (default) x^k, "chebyshev" T_k or "hermite" H_k, the physicists'
  !> Hermite polynomials), method ("horner" (default) or
  !> "paterson_stockmeyer" for monomials; "recurrence" (default) or
  !> "factorized" for Chebyshev), threshold (scoped, default 0).
  !> info: degree.
  subroutine msp_polynomial(A, B, args, info, context)
    type(msp_matrix), intent(in) :: A
    type(msp_matrix), intent(inout) :: B
    type(msp_dict), intent(in) :: args
    type(msp_dict), intent(inout), target :: info
    type(msp_dict), intent(in), target, optional :: context
    type(msp_dict), target :: tracked
    type(msp_args) :: a_args
    type(msp_info) :: r

    call copy(tracked, args)
    call a_args%open(tracked, context)
    call r%open(info)
    call polynomial_impl(A, B, a_args, r)
    call a_args%finish(r)
  end subroutine msp_polynomial

  subroutine polynomial_impl(A, B, args, info)
    type(msp_matrix), intent(in) :: A
    type(msp_matrix), intent(inout) :: B
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    real(dp), dimension(:), allocatable :: c
    character(:), allocatable :: basis, method
    real(dp) :: threshold

    call args%get("coefficients", c)
    call args%get("basis", basis, default="monomial")
    call args%get("threshold", threshold, default=0.0_dp, scoped=.true.)
    select case (basis)
    case ("monomial")
       call args%get("method", method, default="horner")
       call check_method(["horner             ", "paterson_stockmeyer"])
       if (method == "paterson_stockmeyer") &
            & call msp_cite("paterson1973number", info)
       call monomial(A, c, B, threshold, method == "paterson_stockmeyer")
    case ("chebyshev")
       call args%get("method", method, default="recurrence")
       call check_method(["recurrence", "factorized"])
       call chebyshev(A, c, B, threshold, method == "factorized")
    case ("hermite")
       call hermite(A, c, B, threshold)
    case default
       call msp_raise("msp_polynomial: unknown basis '" // basis // "'")
    end select
    call info%set("degree", size(c) - 1)

  contains

    subroutine check_method(known)
      character(*), dimension(:), intent(in) :: known

      if (all(method /= known)) call msp_raise("msp_polynomial: unknown " // &
           & "method '" // method // "' for basis '" // basis // "'")
    end subroutine check_method
  end subroutine polynomial_impl

  !> The identity with the grid, precision and location of A.
  subroutine identity_like(A, id)
    type(msp_matrix), intent(in) :: A
    type(msp_matrix), intent(inout) :: id

    call id%init(A%grid, A%n)
    call id%identity()
    if (A%is_single()) call msp_to_single(id)
    if (A%on_device()) call msp_to_device(id)
  end subroutine identity_like

  !> B = sum of c(k + 1) A^k. Horner's rule takes one product per degree.
  !> Paterson and Stockmeyer (SIAM J. Comput. 2, 60 (1973)) take about
  !> 2 sqrt(degree): with s about sqrt(degree), p is a polynomial in A^s
  !> whose coefficients are polynomials of degree below s in A, built from
  !> the powers I, A, ..., A^s and combined by Horner's rule in A^s.
  subroutine monomial(A, c, B, threshold, stockmeyer)
    type(msp_matrix), intent(in) :: A
    real(dp), dimension(:), intent(in) :: c
    type(msp_matrix), intent(inout) :: B
    real(dp), intent(in) :: threshold
    logical, intent(in) :: stockmeyer
    type(msp_matrix), dimension(:), allocatable :: powers
    type(msp_matrix) :: work, part
    integer :: degree, s, k

    degree = size(c) - 1
    s = 1
    if (stockmeyer) s = max(1, int(sqrt(real(degree, dp))))
    allocate(powers(0:s))
    call identity_like(A, powers(0))
    call copy(powers(1), A)
    do k = 2, s
       call msp_multiply(A, powers(k - 1), powers(k), threshold=threshold)
    end do
    call block(degree / s, B)
    do k = degree / s - 1, 0, -1
       call msp_multiply(powers(s), B, work, threshold=threshold)
       call block(k, part)
       call msp_add(work, part, threshold=threshold)
       call msp_swap(B, work)
    end do
    do k = 0, s
       call powers(k)%free()
    end do
    call work%free()
    call part%free()

  contains

    !> q = sum of c(b s + j + 1) A^j for j = 0, ..., s - 1 within the degree.
    subroutine block(b, q)
      integer, intent(in) :: b
      type(msp_matrix), intent(inout) :: q
      integer :: j

      call copy(q, powers(0))
      call msp_scale(q, c(b * s + 1))
      do j = 1, min(s - 1, degree - b * s)
         call msp_add(q, powers(j), alpha=c(b * s + j + 1), &
              & threshold=threshold)
      end do
    end subroutine block
  end subroutine monomial

  !> B = sum of c(k + 1) T_k(A), with T_k the Chebyshev polynomials. The
  !> recurrence T_(k+1) = 2 A T_k - T_(k-1) takes one product per degree.
  !> The factorized form takes about twice as many products as the degree's
  !> logarithm plus the number of splits: with m the largest power of two
  !> up to the degree d, T_(m+i) = 2 T_m T_i - T_(m-i) for i >= 1 splits p
  !> into L + 2 T_m R, L of degree below m and R of degree d - m (with
  !> c_m T_m in R as c_m / 2 at degree 0), each evaluated the same way from
  !> T_1, T_2, T_4, ... (T_2m = 2 T_m^2 - I).
  subroutine chebyshev(A, c, B, threshold, factorized)
    type(msp_matrix), intent(in) :: A
    real(dp), dimension(:), intent(in) :: c
    type(msp_matrix), intent(inout) :: B
    real(dp), intent(in) :: threshold
    logical, intent(in) :: factorized
    type(msp_matrix), dimension(:), allocatable :: t
    type(msp_matrix) :: id, previous, current, next
    integer :: degree, k, levels

    degree = size(c) - 1
    call identity_like(A, id)
    if (factorized) then
       levels = 0
       do while (2**(levels + 1) <= degree)
          levels = levels + 1
       end do
       allocate(t(0:levels))
       call copy(t(0), A)
       do k = 1, levels
          call msp_multiply(t(k - 1), t(k - 1), t(k), alpha=2.0_dp, &
               & threshold=threshold)
          call msp_add(t(k), id, alpha=-1.0_dp, threshold=threshold)
       end do
       call split(c, B)
       do k = 0, levels
          call t(k)%free()
       end do
       call id%free()
       return
    end if
    call copy(B, id)
    call msp_scale(B, c(1))
    if (degree >= 1) call msp_add(B, A, alpha=c(2), threshold=threshold)
    call copy(previous, id)
    call copy(current, A)
    do k = 2, degree
       call msp_multiply(A, current, next, alpha=2.0_dp, threshold=threshold)
       call msp_add(next, previous, alpha=-1.0_dp, threshold=threshold)
       call msp_add(B, next, alpha=c(k + 1), threshold=threshold)
       call msp_swap(previous, current)
       call msp_swap(current, next)
    end do
    call id%free()
    call previous%free()
    call current%free()
    call next%free()

  contains

    !> out = sum of d(k + 1) T_k(A).
    recursive subroutine split(d, out)
      real(dp), dimension(:), intent(in) :: d
      type(msp_matrix), intent(inout) :: out
      real(dp), dimension(:), allocatable :: left, right
      type(msp_matrix) :: left_part, right_part
      integer :: n, m, j, i

      n = size(d) - 1
      if (n <= 1) then
         call copy(out, id)
         call msp_scale(out, d(1))
         if (n == 1) call msp_add(out, A, alpha=d(2), threshold=threshold)
         return
      end if
      j = 0
      do while (2**(j + 1) <= n)
         j = j + 1
      end do
      m = 2**j
      left = d(1:m)
      do i = 1, n - m
         left(m - i + 1) = left(m - i + 1) - d(m + i + 1)
      end do
      right = d(m + 1:n + 1)
      right(1) = right(1) / 2
      call split(right, right_part)
      call msp_multiply(t(j), right_part, out, alpha=2.0_dp, &
           & threshold=threshold)
      call split(left, left_part)
      call msp_add(out, left_part, threshold=threshold)
      call left_part%free()
      call right_part%free()
    end subroutine split
  end subroutine chebyshev

  !> B = sum of c(k + 1) H_k(A), with H_0 = I, H_1 = 2A and
  !> H_(k+1) = 2 A H_k - 2 k H_(k-1).
  subroutine hermite(A, c, B, threshold)
    type(msp_matrix), intent(in) :: A
    real(dp), dimension(:), intent(in) :: c
    type(msp_matrix), intent(inout) :: B
    real(dp), intent(in) :: threshold
    type(msp_matrix) :: previous, current, next
    integer :: k

    call identity_like(A, previous)
    call copy(B, previous)
    call msp_scale(B, c(1))
    call copy(current, A)
    call msp_scale(current, 2.0_dp)
    if (size(c) > 1) call msp_add(B, current, alpha=c(2), threshold=threshold)
    do k = 2, size(c) - 1
       call msp_multiply(A, current, next, alpha=2.0_dp, threshold=threshold)
       call msp_add(next, previous, alpha=-2.0_dp * (k - 1), &
            & threshold=threshold)
       call msp_add(B, next, alpha=c(k + 1), threshold=threshold)
       call msp_swap(previous, current)
       call msp_swap(current, next)
    end do
    call previous%free()
    call current%free()
    call next%free()
  end subroutine hermite
end module msp_polynomial_mod
