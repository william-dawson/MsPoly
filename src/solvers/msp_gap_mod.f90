!> Gap edge eigenpairs from particle and hole density-matrix filters.
module msp_gap_mod
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, &
       & ieee_quiet_nan
  use msp_citations_mod, only: msp_cite
  use msp_kinds_mod, only: dp
  use msp_error_mod, only: msp_raise
  use msp_dict_mod, only: msp_dict, msp_args, msp_info, copy
  use msp_matrix_mod
  use msp_density_mod, only: density_impl
  implicit none
  private

  public :: msp_estimate_gap

contains

  !> Estimate the HOMO and LUMO of Hermitian H in an orthogonal basis,
  !> using the filters D^2(I-D), D(I-D)^2 of Truflandier,
  !> arXiv:2603.29849v3, Eq. 20, and fully reorthogonalized Lanczos.
  !> Optional density is a quasi-purified D commuting with H, with
  !> occupied eigenvalues near one and unoccupied eigenvalues near zero.
  !> Otherwise the density section supplies purification arguments,
  !> including trace (method defaults to HPCP, as in the paper); the last
  !> iterate above sampling_tolerance is used. Strongly asymmetric
  !> occupation errors can make both filters select the same edge.
  !> A fully idempotent D has no filter signal and cannot give the edges.
  !>
  !> args: tolerance (scoped, default 1e-8, absolute H residual),
  !> max_iterations (default 100), threshold (scoped, default 0, for
  !> constructing filters), seed (default 1, deterministic start vector),
  !> density (args of density_impl). Lanczos arithmetic uses FP64.
  !> info: homo, lumo, gap, converged; density (when computed), occupied
  !> and unoccupied (iterations, residual, filter_residual, occupation,
  !> converged, stop and histories). Failed edges have converged=false.
  !> Optional vectors have the normalized eigenvector in their first
  !> column; the other columns are zero. Collective.
  subroutine msp_estimate_gap(H, args, info, context, density, &
       & homo_vector, lumo_vector)
    type(msp_matrix), intent(in) :: H
    type(msp_dict), intent(in) :: args
    type(msp_dict), intent(inout), target :: info
    type(msp_dict), intent(in), target, optional :: context
    type(msp_matrix), intent(in), optional :: density
    type(msp_matrix), intent(inout), optional :: homo_vector, lumo_vector
    type(msp_dict), target :: tracked
    type(msp_args) :: a, inner_args
    type(msp_info) :: r, inner_info
    type(msp_matrix) :: d, k, d2, variance, particle, hole, adjoint, hv, lv
    real(dp) :: tolerance, threshold, homo, lumo, target
    integer :: limit, seed
    logical :: hc, lc

    call copy(tracked, args)
    call a%open(tracked, context)
    call r%open(info)
    call a%get("tolerance", tolerance, default=1.0e-8_dp, scoped=.true.)
    call a%get("threshold", threshold, default=0.0_dp, scoped=.true.)
    call a%get("max_iterations", limit, default=100)
    call a%get("seed", seed, default=1)
    if (.not. ieee_is_finite(tolerance) .or. tolerance <= 0) &
         & call msp_raise("msp_estimate_gap: tolerance must be positive")
    if (limit < 1) call msp_raise( &
         & "msp_estimate_gap: max_iterations must be positive")
    call msp_cite("truflandier2026gap", r)
    call msp_cite("lanczos1950iteration", r)
    if (present(density)) then
       if (density%n /= H%n .or. &
            & .not. associated(density%grid, H%grid)) call msp_raise( &
            & "msp_estimate_gap: density differs in grid or dimension")
       call copy(d, density)
    else
       call a%enter("density", inner_args)
       call r%enter("density", inner_info)
       call inner_args%get("trace", target)
       if (target <= 0 .or. target >= H%n .or. &
            & target /= real(nint(target), dp)) call msp_raise( &
            & "msp_estimate_gap: density.trace must be an interior integer")
       call density_impl(H, k, inner_args, inner_info, sample=d, &
            & method_default="hpcp")
       call k%free()
    end if
    if (d%n == 0) then
       call r%set("converged", .false.)
       call r%set("stop", "no_sample")
       call r%set("homo", ieee_value(0.0_dp, ieee_quiet_nan))
       call r%set("lumo", ieee_value(0.0_dp, ieee_quiet_nan))
       call r%set("gap", ieee_value(0.0_dp, ieee_quiet_nan))
       if (present(homo_vector)) call homo_vector%free()
       if (present(lumo_vector)) call lumo_vector%free()
       call a%finish(r)
       return
    end if
    call msp_to_double(d)
    if (H%on_device()) call msp_to_device(d)
    call msp_multiply(d, d, d2, threshold=threshold)
    call copy(variance, d)
    call msp_add(variance, d2, alpha=-1.0_dp)
    call msp_multiply(d, variance, particle, threshold=threshold)
    call copy(hole, variance)
    call msp_add(hole, particle, alpha=-1.0_dp)
    !! Truncation can break commutation; keep each filter Hermitian.
    call msp_transpose(particle, adjoint, conjugate=.true.)
    call msp_add(particle, adjoint)
    call msp_scale(particle, 0.5_dp)
    call msp_transpose(hole, adjoint, conjugate=.true.)
    call msp_add(hole, adjoint)
    call msp_scale(hole, 0.5_dp)
    call r%enter("occupied", inner_info)
    call edge(H, d, particle, limit, tolerance, seed, .true., &
         & inner_info, hv, homo, hc)
    call r%enter("unoccupied", inner_info)
    call edge(H, d, hole, limit, tolerance, seed, .false., &
         & inner_info, lv, lumo, lc)
    call r%set("homo", homo)
    call r%set("lumo", lumo)
    call r%set("gap", lumo - homo)
    if (hc .and. lc) then
       call r%set("converged", lumo >= homo)
    else
       call r%set("converged", .false.)
    end if
    if (present(homo_vector)) call msp_swap(homo_vector, hv)
    if (present(lumo_vector)) call msp_swap(lumo_vector, lv)
    call d%free()
    call d2%free()
    call variance%free()
    call particle%free()
    call hole%free()
    call adjoint%free()
    call hv%free()
    call lv%free()
    call a%finish(r)
  end subroutine msp_estimate_gap

  !> Largest algebraic filter eigenpair; certify it against H and D.
  subroutine edge(H, D, filter, limit, tolerance, seed, occupied, info, &
       & vector, energy, converged)
    type(msp_matrix), intent(in) :: H, D
    type(msp_matrix), intent(inout) :: filter, vector
    integer, intent(in) :: limit, seed
    real(dp), intent(in) :: tolerance
    logical, intent(in) :: occupied
    type(msp_info), intent(in) :: info
    real(dp), intent(out) :: energy
    logical, intent(out) :: converged
    type(msp_matrix), dimension(:), allocatable :: basis
    type(msp_matrix) :: work, tmp, scaled
    real(dp), dimension(:), allocatable :: diagonal, off, weights, start
    integer, dimension(:), allocatable :: rows, cols
    real(dp) :: norm, beta, value, residual, f_residual, occupation
    complex(dp) :: coefficient
    integer :: i, j, pass, count, m
    character(:), allocatable :: reason

    energy = ieee_value(0.0_dp, ieee_quiet_nan)
    converged = .false.
    call vector%free()
    norm = msp_norm(filter)
    call info%set("filter_norm", norm)
    if (norm <= 32 * epsilon(1.0_dp)) then
       call info%set("iterations", 0)
       call info%set("converged", .false.)
       call info%set("stop", "vanished_filter")
       return
    end if
    call msp_scale(filter, 1 / norm)
    m = min(limit, H%n)
    allocate(basis(m), diagonal(m), off(m), weights(m))
    diagonal = 0
    off = 0
    rows = [(i, i = H%grid%rank + 1, H%n, H%grid%size)]
    count = size(rows)
    allocate(cols(count), start(count))
    cols = 1
    start = sin(real(rows, dp) * sqrt(2.0_dp) + real(seed, dp)) + &
         & cos(real(rows, dp) * sqrt(3.0_dp) - real(seed, dp))
    call work%init(H%grid, H%n)
    call work%fill(rows, cols, start)
    if (filter%on_device()) call msp_to_device(work)
    norm = msp_norm(work)
    call msp_multiply(filter, work, basis(1))
    beta = msp_norm(basis(1))
    !! A nonzero filter can annihilate the chosen starting vector.
    if (beta <= 64 * epsilon(1.0_dp) * norm) then
       call info%set("iterations", 0)
       call info%set("converged", .false.)
       call info%set("stop", "empty_start")
       call work%free()
       call basis(1)%free()
       return
    end if
    call msp_scale(basis(1), 1 / beta)
    reason = "max_iterations"
    do i = 1, m
       call msp_multiply(filter, basis(i), work)
       call msp_dot(basis(i), work, diagonal(i))
       call msp_add(work, basis(i), alpha=-diagonal(i))
       if (i > 1) call msp_add(work, basis(i - 1), alpha=-off(i - 1))
       !! Two passes keep the basis orthogonal near an invariant subspace.
       do pass = 1, 2
          do j = 1, i
             call msp_dot(basis(j), work, coefficient)
             call copy(scaled, basis(j))
             call msp_scale(scaled, -coefficient)
             call msp_add(work, scaled)
          end do
       end do
       beta = msp_norm(work)
       off(i) = beta
       call ritz_pair(diagonal(1:i), off(1:i), value, weights(1:i))
       call copy(vector, basis(1))
       call msp_scale(vector, weights(1))
       do j = 2, i
          call msp_add(vector, basis(j), alpha=weights(j))
       end do
       call msp_scale(vector, 1 / msp_norm(vector))
       call msp_multiply(H, vector, tmp)
       call msp_dot(vector, tmp, energy)
       call msp_add(tmp, vector, alpha=-energy)
       residual = msp_norm(tmp)
       call msp_multiply(filter, vector, tmp)
       call msp_add(tmp, vector, alpha=-value)
       f_residual = msp_norm(tmp)
       call info%append("energy_history", energy)
       call info%append("residual_history", residual)
       call info%append("filter_residual_history", f_residual)
       converged = residual <= tolerance .and. &
            & f_residual <= min(tolerance, 1.0e-8_dp)
       if (converged) then
          reason = "tolerance"
          exit
       end if
       if (beta < 64 * epsilon(1.0_dp)) then
          reason = "breakdown"
          exit
       end if
       if (i < m) then
          call copy(basis(i + 1), work)
          call msp_scale(basis(i + 1), 1 / beta)
       end if
    end do
    call msp_multiply(D, vector, tmp)
    call msp_dot(vector, tmp, occupation)
    if ((occupied .and. occupation <= 0.5_dp) .or. &
         & (.not. occupied .and. occupation >= 0.5_dp) .or. value <= 0) then
       converged = .false.
       reason = "wrong_occupation"
    end if
    call info%set("iterations", min(i, m))
    call info%set("converged", converged)
    call info%set("stop", reason)
    call info%set("residual", residual)
    call info%set("filter_residual", f_residual)
    call info%set("occupation", occupation)
    call info%set("ritz", value)
    do j = 1, m
       call basis(j)%free()
    end do
    call work%free()
    call tmp%free()
    call scaled%free()
  end subroutine edge

  !> Largest Ritz pair of the real symmetric Lanczos tridiagonal.
  !> Jacobi rotations solve only this small projected problem.
  subroutine ritz_pair(diagonal, off, value, weights)
    real(dp), dimension(:), intent(in) :: diagonal, off
    real(dp), intent(out) :: value
    real(dp), dimension(:), intent(out) :: weights
    real(dp), dimension(size(diagonal), size(diagonal)) :: a, z
    real(dp) :: tau, t, c, s, ap, aq, maximum
    integer :: n, p, q, j, sweep, largest

    n = size(diagonal)
    a = 0
    z = 0
    do p = 1, n
       a(p, p) = diagonal(p)
       z(p, p) = 1
       if (p == n) cycle
       a(p, p + 1) = off(p)
       a(p + 1, p) = off(p)
    end do
    do sweep = 1, 100
       maximum = 0
       do q = 2, n
          do p = 1, q - 1
             maximum = max(maximum, abs(a(p, q)))
             if (abs(a(p, q)) < epsilon(1.0_dp)) cycle
             tau = (a(q, q) - a(p, p)) / (2 * a(p, q))
             t = sign(1.0_dp, tau) / (abs(tau) + sqrt(1 + tau**2))
             c = 1 / sqrt(1 + t * t)
             s = t * c
             a(p, p) = a(p, p) - t * a(p, q)
             a(q, q) = a(q, q) + t * a(p, q)
             a(p, q) = 0
             a(q, p) = 0
             do j = 1, n
                if (j /= p .and. j /= q) then
                   ap = a(j, p)
                   aq = a(j, q)
                   a(j, p) = c * ap - s * aq
                   a(p, j) = a(j, p)
                   a(j, q) = s * ap + c * aq
                   a(q, j) = a(j, q)
                end if
                ap = z(j, p)
                aq = z(j, q)
                z(j, p) = c * ap - s * aq
                z(j, q) = s * ap + c * aq
             end do
          end do
       end do
       if (maximum < 8 * epsilon(1.0_dp)) exit
    end do
    largest = 1
    do p = 2, n
       if (a(p, p) > a(largest, largest)) largest = p
    end do
    value = a(largest, largest)
    weights = z(:, largest)
  end subroutine ritz_pair
end module msp_gap_mod
