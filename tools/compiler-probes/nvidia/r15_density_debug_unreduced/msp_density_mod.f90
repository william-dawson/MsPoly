!> Density matrices of Hermitian Hamiltonians by purification.
module msp_density_mod
  use msp_citations_mod, only: msp_cite
  use msp_kinds_mod, only: dp
  use msp_error_mod, only: msp_raise
  use msp_dict_mod, only: msp_dict, msp_args, msp_info, copy
  use msp_matrix_mod
  use msp_bounds_mod, only: bounds_impl
  use msp_sqrt_mod, only: sqrt_impl
  implicit none
  private

  public :: msp_density, density_impl, msp_energy_density, msp_mcweeny

  !! The scalar polynomial each step applies to an eigenvalue p of X: p^2,
  !! 2p - p^2, and the steps of TRS4, PM (two forms), HPCP and scale and
  !! fold (two forms), with the step's parameter.
  integer, parameter :: square = 1, fold = 2, trs4_step = 3, pm_high = 4, &
       & pm_low = 5, hpcp_step = 6, shift_square = 7, scaled_fold = 8

contains

  !> K becomes the density matrix of the lowest eigenvectors of H, as many
  !> as args "trace": the projector onto them. Collective. With an overlap
  !> matrix S, these are the generalized eigenvectors of H C = S C E, and
  !> K S is the projector.
  !>
  !> args: trace (required), method ("trs2" (default), "trs4", "pm",
  !> "hpcp" or "scale_and_fold", which also needs homo and lumo, the
  !> highest occupied and lowest unoccupied eigenvalues), threshold (scoped,
  !> default 0), tolerance (scoped, default 1e-6), max_iterations (default
  !> 1000). K has the precision of H: with H in single precision, every
  !> iterate is stored in single precision.
  !> Optional sample keeps the last iterate with |X-X^2|_F above
  !> sampling_tolerance (default 1e-4), in the orthogonal basis when an
  !> overlap is given. Its iteration and error are reported under sample.
  !> info: energy = tr(K H), chemical_potential, iterations, converged,
  !> stop ("tolerance", "stagnation" or "max_iterations"), energy_history,
  !> idempotency_history, bounds, and with an overlap inverse_sqrt (args and
  !> info of the inverse square root of S; see sqrt_impl).
  subroutine msp_density(H, K, args, info, context, overlap, sample)
    type(msp_matrix), intent(in) :: H
    type(msp_matrix), intent(inout) :: K
    type(msp_dict), intent(in) :: args
    type(msp_dict), intent(inout), target :: info
    type(msp_dict), intent(in), target, optional :: context
    type(msp_matrix), intent(in), optional :: overlap
    type(msp_matrix), intent(inout), optional :: sample
    type(msp_dict), target :: tracked
    type(msp_args) :: a
    type(msp_info) :: r

    call copy(tracked, args)
    call a%open(tracked, context)
    call r%open(info)
    call density_impl(H, K, a, r, overlap, sample)
    call a%finish(r)
  end subroutine msp_density

  !> Purification: X starts as a linear function of H, c0 I + c1 H, with the
  !> spectrum of H mapped into [0, 1], and each step applies a polynomial
  !> that pushes every eigenvalue to 0 or 1 while steering the trace to the
  !> target. All methods square X first; the new X is a combination of I, X,
  !> X^2 and, for TRS4, PM and HPCP, one more product.
  !>
  !> It stops when the energy changes by at most the tolerance (a negative
  !> tolerance turns this test off), or when numerical errors dominate,
  !> detected from the order of convergence (Kruchinina, Rudberg and
  !> Rubensson, J. Chem. Theory Comput. 12, 5788 (2016)): with the error
  !> e = |X - X^2|_F and P the polynomial of the last two steps, exact
  !> arithmetic gives e_k <= C e_(k-2)^2, where C is the maximum over
  !> [0, 1] of (P(x) - P(x)^2) / (x - x^2)^2. The expansion stops when
  !> e_k > C e_(k-2)^1.8 (observed order below 1.8) and keeps X as it was
  !> before that step. For two different TRS2 steps C = 4.409, the paper's
  !> constant; here C is computed for whatever pair was applied. Two equal
  !> TRS2 steps converge to first order only at one end, so their C is
  !> unbounded and the test cannot stop the trace-correcting phase early.
  !>
  !> The chemical potential is the eigenvalue that the composed step
  !> polynomials map to 1/2, found by bisection and mapped back.
  !>
  !> With an overlap S, purification runs on S^(-1/2) H S^(-1/2) in the
  !> orthogonal basis, and K = S^(-1/2) X S^(-1/2).
  subroutine density_impl(H, K, args, info, overlap, sample, method_default)
    type(msp_matrix), intent(in) :: H
    type(msp_matrix), intent(inout) :: K
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    type(msp_matrix), intent(in), optional :: overlap
    type(msp_matrix), intent(inout), optional :: sample
    character(*), intent(in), optional :: method_default
    type(msp_args) :: inner_args
    type(msp_info) :: bounds_info, inner_info
    type(msp_matrix) :: x2, x3, f, id, last, isq, work, wh
    real(dp) :: sampling_tolerance
    real(dp) :: target, threshold, tolerance, emin, emax, trace, trace2, &
         & trace3, trace_f, trace_g, energy, previous, homo, lumo
    !! Per step: the idempotency error, the polynomial applied and its
    !! parameter.
    real(dp), dimension(:), allocatable :: e, theta
    integer, dimension(:), allocatable :: steps
    integer :: max_iterations, step, iterations
    character(:), allocatable :: method, reason

    call args%get("trace", target)
    if (present(method_default)) then
       call args%get("method", method, default=method_default)
    else
       call args%get("method", method, default="trs2")
    end if
    call args%get("threshold", threshold, default=0.0_dp, scoped=.true.)
    call args%get("tolerance", tolerance, default=1.0e-6_dp, scoped=.true.)
    call args%get("max_iterations", max_iterations, default=1000)
    if (present(sample)) then
       call sample%free()
       call args%get("sampling_tolerance", sampling_tolerance, &
            & default=1.0e-4_dp)
       if (sampling_tolerance <= 0) call msp_raise( &
            & "msp_density: sampling_tolerance must be positive")
    end if
    if (method == "scale_and_fold") then
       call args%get("homo", homo)
       call args%get("lumo", lumo)
    else if (all(method /= [character(5) :: "trs2", "trs4", "pm", "hpcp"])) &
         & then
       call msp_raise("msp_density: unknown method '" // method // "'")
    end if
    select case (method)
    case ("trs2", "trs4")
       call msp_cite("niklasson2002", info)
    case ("pm")
       call msp_cite("palser1998canonical", info)
    case ("hpcp")
       call msp_cite("truflandier2016communication", info)
    case ("scale_and_fold")
       call msp_cite("rubensson2011nonmonotonic", info)
    end select
    if (.not. present(overlap)) then
       call purify(H)
       return
    end if
    call args%enter("inverse_sqrt", inner_args)
    call info%enter("inverse_sqrt", inner_info)
    call sqrt_impl(overlap, inner_args, inner_info, inverse_root=isq)
    call msp_multiply(isq, H, work, threshold=threshold)
    call msp_multiply(work, isq, wh, threshold=threshold)
    call purify(wh)
    call msp_multiply(isq, K, work, threshold=threshold)
    call msp_multiply(work, isq, K, threshold=threshold)
    call isq%free()
    call work%free()
    call wh%free()

  contains

    !> K becomes the density matrix of hw, which has an orthogonal basis.
    subroutine purify(hw)
      type(msp_matrix), intent(in) :: hw
      real(dp) :: mu_offset, mu_scale, sigma, mean, scale, beta, beta_bar, &
           & alpha
      integer :: n

      call info%enter("bounds", bounds_info)
      call bounds_impl(hw, bounds_info, emin, emax)
      n = hw%n
      call id%init(hw%grid, n)
      call id%identity()
      if (hw%is_single()) call msp_to_single(id)
      if (hw%on_device()) call msp_to_device(id)
      call msp_trace(hw, mean)
      mean = mean / n
      !! X = c0 I + c1 H; the chemical potential is mu_offset +
      !! mu_scale * p for the eigenvalue p of X that ends at 1/2.
      select case (method)
      case ("pm")
         scale = min(target / (emax - mean), (n - target) / (mean - emin))
         call start(hw, (scale * mean + target) / n, -scale / n)
         mu_offset = mean + target / scale
         mu_scale = -n / scale
      case ("hpcp")
         sigma = target / n
         scale = min(sigma / (emax - mean), (1 - sigma) / (mean - emin))
         call start(hw, sigma + scale * mean, -scale)
         mu_offset = mean + sigma / scale
         mu_scale = -1 / scale
      case default
         call start(hw, emax / (emax - emin), -1 / (emax - emin))
         mu_offset = emax
         mu_scale = emin - emax
      end select
      if (method == "scale_and_fold") then
         !! Where the homo and lumo lie in the spectrum of X.
         beta = (emax - lumo) / (emax - emin)
         beta_bar = (emax - homo) / (emax - emin)
      end if

      allocate(steps(max_iterations), theta(max_iterations), &
           & e(max_iterations))
      energy = 0
      reason = "max_iterations"
      iterations = 0
      do step = 1, max_iterations
         call msp_trace(K, trace)
         call msp_multiply(K, K, x2, threshold=threshold)
         call msp_trace(x2, trace2)
         call copy(f, K)
         call msp_add(f, x2, alpha=-1.0_dp)
         e(step) = msp_norm(f)
         call info%append("idempotency_history", e(step))
         if (present(sample)) then
            if (e(step) > sampling_tolerance) then
               call copy(sample, K)
               call info%set("sample.iteration", step - 1)
               call info%set("sample.idempotency", e(step))
            end if
         end if
         if (step > 2) then
            if (step == 3) call msp_cite("kruchinina2016parameterless", info)
            if (e(step - 2) < 1 .and. e(step) > order_constant(steps(step &
                 & - 2:step - 1), theta(step - 2:step - 1)) * &
                 & e(step - 2)**1.8_dp) then
               reason = "stagnation"
               !! Return to X before the step that broke the order.
               call msp_swap(K, last)
               iterations = step - 2
               call msp_dot(K, hw, energy)
               exit
            end if
         end if
         theta(step) = 0
         select case (method)
         case ("trs2")
            steps(step) = merge(square, fold, trace > target)
            if (steps(step) == square) then
               call advance(x2)
            else
               call combine([0.0_dp, 2.0_dp, -1.0_dp, 0.0_dp])
            end if
         case ("trs4")
            !! sigma makes tr(X^2 (4X - 3X^2) + sigma X^2 (I - X)^2) the
            !! target; outside [0, 6] a TRS2 step is safer. Once
            !! tr(X^2 (I - X)^2) is at the rounding level of a trace, X is
            !! idempotent to rounding and sigma would be noise: any sigma in
            !! [0, 6] then converges to second order at 0 and 1, so take 3.
            call copy(f, K)
            call msp_scale(f, 4.0_dp)
            call msp_add(f, x2, alpha=-3.0_dp)
            call msp_dot(x2, f, trace_f)
            call copy(f, id)
            call msp_add(f, K, alpha=-2.0_dp)
            call msp_add(f, x2)
            call msp_dot(x2, f, trace_g)
            sigma = 3
            if (trace_g > n * epsilon(1.0_dp)) sigma = (target - trace_f) / &
                 & trace_g
            if (sigma > 6) then
               steps(step) = fold
               call combine([0.0_dp, 2.0_dp, -1.0_dp, 0.0_dp])
            else if (sigma < 0) then
               steps(step) = square
               call advance(x2)
            else
               steps(step) = trs4_step
               theta(step) = sigma
               call copy(f, K)
               call msp_scale(f, 4 - 2 * sigma)
               call msp_add(f, x2, alpha=sigma - 3)
               call msp_add(f, id, alpha=sigma)
               call msp_multiply(x2, f, x3, threshold=threshold)
               call advance(x3)
            end if
         case ("pm", "hpcp")
            call msp_multiply(K, x2, x3, threshold=threshold)
            call msp_trace(x3, trace3)
            !! A mean of the eigenvalues of X weighted by p (1 - p), so in
            !! [0, 1]; once X is idempotent to rounding it is noise.
            sigma = min(1.0_dp, max(0.0_dp, (trace2 - trace3) / &
                 & max(trace - trace2, tiny(1.0_dp))))
            theta(step) = sigma
            if (method == "hpcp") then
               steps(step) = hpcp_step
               call combine([0.0_dp, 1 - 2 * sigma, 2 + 2 * sigma, -2.0_dp])
            else if (sigma > 0.5_dp) then
               steps(step) = pm_high
               call combine([0.0_dp, 0.0_dp, 1 + 1 / sigma, -1 / sigma])
            else
               steps(step) = pm_low
               call combine([0.0_dp, 1 - 2 * sigma, 1 + sigma, -1.0_dp] / &
                    & (1 - sigma))
            end if
         case default
            if (trace > target) then
               steps(step) = shift_square
               alpha = 2 / (2 - beta)
               call combine([(1 - alpha)**2, 2 * alpha * (1 - alpha), &
                    & alpha**2, 0.0_dp])
            else
               steps(step) = scaled_fold
               alpha = 2 / (1 + beta_bar)
               call combine([0.0_dp, 2 * alpha, -alpha**2, 0.0_dp])
            end if
            theta(step) = alpha
            beta = polynomial(steps(step:step), theta(step:step), beta)
            beta_bar = polynomial(steps(step:step), theta(step:step), &
                 & beta_bar)
         end select
         previous = energy
         call msp_dot(K, hw, energy)
         call info%append("energy_history", energy)
         iterations = step
         if (abs(energy - previous) <= tolerance) then
            reason = "tolerance"
            exit
         end if
      end do
      call info%set("energy", energy)
      call info%set("chemical_potential", mu_offset + mu_scale * &
           & bisect(steps(1:iterations), theta(1:iterations), tolerance))
      call info%set("iterations", iterations)
      call info%set("converged", reason /= "max_iterations")
      call info%set("stop", reason)
      call x2%free()
      call x3%free()
      call f%free()
      call id%free()
      call last%free()
    end subroutine purify

    !> K = c0 I + c1 hw.
    subroutine start(hw, c0, c1)
      type(msp_matrix), intent(in) :: hw
      real(dp), intent(in) :: c0, c1

      call copy(K, id)
      call msp_scale(K, c0)
      call msp_add(K, hw, alpha=c1)
    end subroutine start

    !> K = c(1) I + c(2) K + c(3) X^2 + c(4) X^3, dropping entries below
    !> the threshold.
    subroutine combine(c)
      real(dp), dimension(4), intent(in) :: c

      call msp_scale(x2, c(3))
      call msp_add(x2, K, alpha=c(2), threshold=threshold)
      if (c(4) /= 0) call msp_add(x2, x3, alpha=c(4), threshold=threshold)
      if (c(1) /= 0) call msp_add(x2, id, alpha=c(1), threshold=threshold)
      call advance(x2)
    end subroutine combine

    !> The next X is new; the current one is kept as the last.
    subroutine advance(new)
      type(msp_matrix), intent(inout) :: new

      call msp_swap(last, K)
      call msp_swap(K, new)
    end subroutine advance
  end subroutine density_impl

  !> The maximum over (0, 1) of (P(x) - P(x)^2) / (x - x^2)^2 for the
  !> composition P of the steps, sampled on a fine grid. Steps that are first
  !> order at 0 or 1 give a large value, growing with the grid.
  pure real(dp) function order_constant(steps, theta) result(c)
    integer, dimension(:), intent(in) :: steps
    real(dp), dimension(:), intent(in) :: theta
    integer, parameter :: points = 4000
    real(dp) :: x, p
    integer :: j

    c = 0
    do j = 1, points
       x = real(j, dp) / (points + 1)
       p = polynomial(steps, theta, x)
       c = max(c, (p - p * p) / (x - x * x)**2)
    end do
  end function order_constant

  !> The steps applied in turn to p.
  pure real(dp) function polynomial(steps, theta, p0) result(p)
    integer, dimension(:), intent(in) :: steps
    real(dp), dimension(:), intent(in) :: theta
    real(dp), intent(in) :: p0
    real(dp) :: t
    integer :: k

    p = p0
    do k = 1, size(steps)
       t = theta(k)
       select case (steps(k))
       case (square)
          p = p * p
       case (fold)
          p = 2 * p - p * p
       case (trs4_step)
          p = p * p * (4 * p - 3 * p * p) + t * p * p * (1 - p)**2
       case (pm_high)
          p = ((1 + t) * p * p - p**3) / t
       case (pm_low)
          p = ((1 - 2 * t) * p + (1 + t) * p * p - p**3) / (1 - t)
       case (hpcp_step)
          p = p + 2 * (p * p * (1 - p) - t * p * (1 - p))
       case (shift_square)
          p = (t * p + 1 - t)**2
       case default
          p = 2 * t * p - t * t * p * p
       end select
    end do
  end function polynomial

  !> The point of [0, 1] that the steps map to 1/2, by bisection: the steps
  !> map [0, 1] onto itself, increasing.
  real(dp) function bisect(steps, theta, tolerance) result(middle)
    integer, dimension(:), intent(in) :: steps
    real(dp), dimension(:), intent(in) :: theta
    real(dp), intent(in) :: tolerance
    real(dp) :: low, high, p

    low = 0
    high = 1
    do
       middle = (low + high) / 2
       p = polynomial(steps, theta, middle)
       if (abs(p - 0.5_dp) < tolerance .or. high - low < epsilon(1.0_dp)) exit
       if (p < 0.5_dp) then
          low = middle
       else
          high = middle
       end if
    end do
  end function bisect

  !> ED = D H D, the energy-weighted density matrix. Collective.
  subroutine msp_energy_density(H, D, ED, threshold)
    type(msp_matrix), intent(in) :: H, D
    type(msp_matrix), intent(inout) :: ED
    real(dp), intent(in), optional :: threshold
    type(msp_matrix) :: dh

    call msp_multiply(D, H, dh, threshold=threshold)
    call msp_multiply(dh, D, ED, threshold=threshold)
    call dh%free()
  end subroutine msp_energy_density

  !> One McWeeny purification step, D' = 3 D S D - 2 D S D S D, with S the
  !> identity when no overlap is given. Collective.
  subroutine msp_mcweeny(D, out, overlap, threshold)
    type(msp_matrix), intent(in) :: D
    type(msp_matrix), intent(inout) :: out
    type(msp_matrix), intent(in), optional :: overlap
    real(dp), intent(in), optional :: threshold
    type(msp_matrix) :: ds, dsd

    if (present(overlap)) then
       call msp_multiply(D, overlap, ds, threshold=threshold)
    else
       call copy(ds, D)
    end if
    call msp_multiply(ds, D, dsd, threshold=threshold)
    call msp_multiply(ds, dsd, out, threshold=threshold)
    call msp_scale(out, -2.0_dp)
    call msp_add(out, dsd, alpha=3.0_dp, threshold=threshold)
    call ds%free()
    call dsd%free()
  end subroutine msp_mcweeny
end module msp_density_mod
