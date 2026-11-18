!> Bounds and estimates of the spectrum of a Hermitian matrix, and operator
!> norms.
module msp_bounds_mod
  use mpi_f08, only: MPI_Allreduce, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, &
       & MPI_SUM, MPI_MIN, MPI_MAX, MPI_DOUBLE_COMPLEX
  use msp_kinds_mod, only: dp
  use msp_error_mod, only: msp_raise
  use msp_dict_mod, only: msp_dict, msp_args, msp_info, copy
  use msp_matrix_mod
  implicit none
  private

  public :: msp_bounds, bounds_impl, gershgorin, radius_bound, power_impl, &
       & operator_norms

contains

  !> Bounds on the spectrum of the Hermitian matrix H. Collective.
  !>
  !> args: method ("gershgorin" (default): every eigenvalue lies in
  !> [min, max]; "power": the eigenvalue of largest magnitude, estimated;
  !> see power_impl).
  !> info: min and max, or largest and iterations.
  subroutine msp_bounds(H, args, info, context)
    type(msp_matrix), intent(in) :: H
    type(msp_dict), intent(in) :: args
    type(msp_dict), intent(inout), target :: info
    type(msp_dict), intent(in), target, optional :: context
    type(msp_dict), target :: tracked
    type(msp_args) :: a
    type(msp_info) :: r
    character(:), allocatable :: method
    real(dp) :: emin, emax

    call copy(tracked, args)
    call a%open(tracked, context)
    call r%open(info)
    call a%get("method", method, default="gershgorin")
    select case (method)
    case ("gershgorin")
       call bounds_impl(H, r, emin, emax)
    case ("power")
       call power_impl(H, a, r, emax)
    case default
       call msp_raise("msp_bounds: unknown method '" // method // "'")
    end select
    call a%finish(r)
  end subroutine msp_bounds

  !> The bounds a solver scales with, from optional args "method" (local):
  !> "gershgorin" (default) bounds the spectrum on both sides; "power"
  !> runs the power method and fills only a magnitude, max as the estimate
  !> of the largest-magnitude eigenvalue and min its negation. The power
  !> estimate is not a certified bound: it suits scalings whose final
  !> rescaling corrects a moderate error (square roots), not mappings that
  !> consume the bounds as exact (initial density guesses). Collective.
  !> info: min and max, or largest, iterations, converged and stop.
  subroutine bounds_impl(H, info, emin, emax, args)
    type(msp_matrix), intent(in) :: H
    type(msp_info), intent(in) :: info
    real(dp), intent(out) :: emin, emax
    type(msp_args), intent(in), optional :: args
    character(:), allocatable :: method

    method = "gershgorin"
    if (present(args)) call args%get("method", method, default="gershgorin")
    select case (method)
    case ("gershgorin")
       call gershgorin(H, emin, emax)
    case ("power")
       call power_impl(H, args, info, emax)
       emin = -emax
       return
    case default
       call msp_raise("msp_bounds: unknown method '" // method // "'")
    end select
    call info%set("min", emin)
    call info%set("max", emax)
  end subroutine bounds_impl

  !> Gershgorin bounds: min over rows of h_ii - r_i and max of h_ii + r_i,
  !> with r_i the sum of |h_ij| over j /= i. They need one pass over the
  !> entries, against many multiplies for a power or Lanczos estimate.
  !> Collective.
  subroutine gershgorin(H, emin, emax)
    type(msp_matrix), intent(in) :: H
    real(dp), intent(out) :: emin, emax
    complex(dp), dimension(:), allocatable :: diag
    real(dp), dimension(:), allocatable :: radius, column

    call sums(H, diag, radius, column)
    emin = minval(real(diag, dp) - radius)
    emax = maxval(real(diag, dp) + radius)
    call MPI_Allreduce(MPI_IN_PLACE, emin, 1, MPI_DOUBLE_PRECISION, MPI_MIN, &
         & H%grid%comm)
    call MPI_Allreduce(MPI_IN_PLACE, emax, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
         & H%grid%comm)
  end subroutine gershgorin

  !> An upper bound on the spectral radius of the Hermitian matrix H: the
  !> smaller of the Gershgorin bound, tight for sparse matrices, and the
  !> Frobenius norm, tight for small dense ones. Collective.
  real(dp) function radius_bound(H) result(r)
    type(msp_matrix), intent(in) :: H
    real(dp) :: emin, emax

    call gershgorin(H, emin, emax)
    r = min(max(abs(emin), abs(emax)), msp_norm(H))
  end function radius_bound

  !> The largest absolute column sum (one) and row sum (infinity) of A: the
  !> operator 1- and infinity-norms. Collective.
  subroutine operator_norms(A, one, infinity)
    type(msp_matrix), intent(in) :: A
    real(dp), intent(out) :: one, infinity
    complex(dp), dimension(:), allocatable :: diag
    real(dp), dimension(:), allocatable :: radius, column

    call sums(A, diag, radius, column)
    infinity = maxval(abs(diag) + radius)
    one = maxval(column)
    call MPI_Allreduce(MPI_IN_PLACE, infinity, 1, MPI_DOUBLE_PRECISION, &
         & MPI_MAX, A%grid%comm)
    call MPI_Allreduce(MPI_IN_PLACE, one, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
         & A%grid%comm)
  end subroutine operator_norms

  !> For the local rows of A: the diagonal, the sum of the
  !> other absolute values; for the local columns: the absolute sum. Summed
  !> over the grid row or column, so every rank holds its rows and columns.
  subroutine sums(A, diag, radius, column)
    type(msp_matrix), intent(in) :: A
    complex(dp), dimension(:), allocatable, intent(out) :: diag
    real(dp), dimension(:), allocatable, intent(out) :: radius, column
    integer, dimension(:), allocatable :: i, j
    complex(dp), dimension(:), allocatable :: v
    integer :: q

    call A%entries(i, j, v)
    allocate(diag(A%rows), radius(A%rows), column(A%cols))
    diag = 0
    radius = 0
    column = 0
    do q = 1, size(i)
       if (i(q) == j(q)) then
          diag(i(q) - A%row0) = diag(i(q) - A%row0) + v(q)
       else
          radius(i(q) - A%row0) = radius(i(q) - A%row0) + abs(v(q))
       end if
       column(j(q) - A%col0) = column(j(q) - A%col0) + abs(v(q))
    end do
    call MPI_Allreduce(MPI_IN_PLACE, diag, A%rows, MPI_DOUBLE_COMPLEX, &
         & MPI_SUM, A%grid%row_comm)
    call MPI_Allreduce(MPI_IN_PLACE, radius, A%rows, MPI_DOUBLE_PRECISION, &
         & MPI_SUM, A%grid%row_comm)
    call MPI_Allreduce(MPI_IN_PLACE, column, A%cols, MPI_DOUBLE_PRECISION, &
         & MPI_SUM, A%grid%col_comm)
  end subroutine sums

  !> Power method for the eigenvalue of largest magnitude of the Hermitian
  !> matrix H. The vector is held as the first column of a matrix, so a
  !> product with H advances it. Rayleigh quotients converge linearly;
  !> Aitken's extrapolation of the last three accelerates them.
  !>
  !> args: max_iterations (default 30), tolerance (scoped, default 1e-6:
  !> stop when the extrapolated estimate changes by at most this).
  !> info: largest, iterations, converged, stop (tolerance, max_iterations,
  !> zero_operator, or empty_start when H annihilates the starting vector).
  subroutine power_impl(H, args, info, largest)
    type(msp_matrix), intent(in) :: H
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    real(dp), intent(out) :: largest
    type(msp_matrix) :: x, y
    real(dp), dimension(3) :: ritz
    real(dp) :: tolerance, quotient, previous, denominator, norm, extrapolated
    integer, dimension(:), allocatable :: rows
    real(dp), dimension(:), allocatable :: values
    integer :: max_iterations, step, k
    logical :: converged
    character(:), allocatable :: reason

    call args%get("max_iterations", max_iterations, default=30)
    call args%get("tolerance", tolerance, default=1.0e-6_dp, scoped=.true.)
    allocate(rows(0), values(0))
    if (H%grid%rank == 0) then
       rows = [(k, k = 1, H%n)]
       values = [(1 / sqrt(real(H%n, dp)), k = 1, H%n)]
    end if
    call x%init(H%grid, H%n)
    call x%fill(rows, [(1, k = 1, size(rows))], values)
    if (H%is_single()) call msp_to_single(x)
    if (H%on_device()) call msp_to_device(x)
    ritz = 0
    largest = 0
    converged = .false.
    reason = "max_iterations"
    do step = 1, max_iterations
       call msp_multiply(H, x, y)
       call msp_dot(x, y, quotient)
       ritz = [ritz(2:3), quotient]
       previous = largest
       largest = quotient
       if (step >= 3) then
          denominator = ritz(3) - 2 * ritz(2) + ritz(1)
          if (abs(denominator) > epsilon(1.0_dp) * abs(ritz(3))) then
             extrapolated = (ritz(3) * ritz(1) - ritz(2)**2) / denominator
             !! Once the quotients have converged the second difference is
             !! rounding noise and the extrapolation can land anywhere;
             !! keep it only when it agrees with the direct estimate.
             if (abs(extrapolated - quotient) <= 1.0e-2_dp * &
                  & abs(quotient)) largest = extrapolated
          end if
       end if
       norm = msp_norm(y)
       if (norm == 0) then
          converged = msp_norm(H) == 0
          if (converged) then
             reason = "zero_operator"
          else
             reason = "empty_start"
          end if
          exit
       end if
       call msp_scale(y, 1 / norm)
       call msp_swap(x, y)
       converged = step > 1 .and. abs(largest - previous) <= tolerance
       if (converged) then
          reason = "tolerance"
          exit
       end if
    end do
    call info%set("largest", largest)
    call info%set("iterations", min(step, max_iterations))
    call info%set("converged", converged)
    call info%set("stop", reason)
    call x%free()
    call y%free()
  end subroutine power_impl
end module msp_bounds_mod
