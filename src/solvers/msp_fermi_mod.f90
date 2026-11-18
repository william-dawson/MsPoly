!> Density matrices at finite temperature: the Fermi-Dirac function of a
!> Hamiltonian.
module msp_fermi_mod
  use msp_citations_mod, only: msp_cite
  use msp_kinds_mod, only: dp
  use msp_error_mod, only: msp_raise
  use msp_dict_mod, only: msp_dict, msp_args, msp_info, copy
  use msp_matrix_mod
  use msp_sqrt_mod, only: sqrt_impl
  use msp_polynomial_mod, only: identity_like
  use msp_eigen_mod, only: eigen_impl, eigen_compose
  implicit none
  private

  public :: msp_fermi, fermi_impl

contains

  !> K becomes the Fermi-Dirac density matrix 1 / (1 + exp(beta (H - mu)))
  !> of the Hermitian H at inverse temperature beta, with either the
  !> chemical potential mu given (grand canonical) or the trace of K given
  !> (canonical, mu then implied). Collective. With an overlap S, the
  !> function is of S^(-1/2) H S^(-1/2) and K = S^(-1/2) F S^(-1/2).
  !> args and info as for fermi_impl.
  subroutine msp_fermi(H, K, args, info, context, overlap)
    type(msp_matrix), intent(in) :: H
    type(msp_matrix), intent(inout) :: K
    type(msp_dict), intent(in) :: args
    type(msp_dict), intent(inout), target :: info
    type(msp_dict), intent(in), target, optional :: context
    type(msp_matrix), intent(in), optional :: overlap
    type(msp_dict), target :: tracked
    type(msp_args) :: a
    type(msp_info) :: r

    call copy(tracked, args)
    call a%open(tracked, context)
    call r%open(info)
    call fermi_impl(H, K, a, r, overlap)
    call a%finish(r)
  end subroutine msp_fermi

  !> Wave operator minimization (Weber, Head-Gordon et al., J. Chem. Phys.
  !> 156, 224106 (2022)): with K = W^2, W evolves in beta from the infinite
  !> temperature W = I / sqrt(2) (or sqrt(trace / n) I) by
  !> dW/dbeta = -(W - W^3) A / 2, with A = H - mu I; in the canonical form
  !> the component of the right side that would change tr(W^2) is removed,
  !> which keeps the trace. Each step is Heun's method (explicit trapezoid);
  !> the difference from the Euler step estimates the local error, which
  !> sets the step size: steps whose estimate exceeds the step tolerance are
  !> retried smaller, and the next step grows at most twofold, so a nearly
  !> stationary W cannot propose a step that overflows. The equation is
  !> stiff: Heun's method is stable only for steps below about
  !> 4 / max |H - mu|, so a large inverse temperature takes many steps, and
  !> near the zero-temperature limit W hovers within the step tolerance's
  !> reach of stationary; the early stop needs a tolerance above that.
  !>
  !> With method "dense", the occupations come from the eigendecomposition
  !> of H (as NTPoly's ComputeDenseFOE): f(e) = 1 / (1 + exp(beta (e - mu))),
  !> with mu by bisection for a given trace; without an inverse temperature,
  !> the lowest floor(trace) states are filled, the next one fractionally
  !> (or the states below mu, those at mu half), and mu for a trace lies
  !> between the highest filled and lowest empty state.
  !>
  !> args: method ("wom" (default) or "dense"), inverse_temperature
  !> (required, except for method "dense"), chemical_potential (grand
  !> canonical) or trace (canonical), step_tolerance (default 1e-2; the
  !> local error is of order step^2, and the global error of the same order
  !> as this), tolerance (scoped, default 0: stop early once |dW/dbeta|_F is
  !> at most this, W then being stationary), max_steps (default 100000),
  !> threshold (scoped, default 0), inverse_sqrt (args of the inverse square
  !> root of the overlap); for method "dense", threshold (scoped, default
  !> 0) and eigen (args of the eigendecomposition; see eigen_impl).
  !> info: energy = tr(K H), inverse_temperature_reached, steps,
  !> evaluations, converged (the target reached, or W stationary),
  !> inverse_sqrt; for method "dense", energy, chemical_potential, eigen.
  subroutine fermi_impl(H, K, args, info, overlap)
    type(msp_matrix), intent(in) :: H
    type(msp_matrix), intent(inout) :: K
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    type(msp_matrix), intent(in), optional :: overlap
    type(msp_args) :: inner_args
    type(msp_info) :: inner_info
    type(msp_matrix) :: isq, work, wh
    character(:), allocatable :: method

    call args%get("method", method, default="wom")
    if (method /= "wom" .and. method /= "dense") call msp_raise( &
         & "msp_fermi: unknown method '" // method // "'")
    if (.not. present(overlap)) then
       call occupy(H, K, method, args, info)
       return
    end if
    call args%enter("inverse_sqrt", inner_args)
    call info%enter("inverse_sqrt", inner_info)
    call sqrt_impl(overlap, inner_args, inner_info, inverse_root=isq)
    call msp_multiply(isq, H, work)
    call msp_multiply(work, isq, wh)
    call occupy(wh, K, method, args, info)
    call msp_multiply(isq, K, work)
    call msp_multiply(work, isq, K)
    call isq%free()
    call work%free()
    call wh%free()
  end subroutine fermi_impl

  !> K for hw in an orthogonal basis, by the method.
  subroutine occupy(hw, K, method, args, info)
    type(msp_matrix), intent(in) :: hw
    type(msp_matrix), intent(inout) :: K
    character(*), intent(in) :: method
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info

    if (method == "dense") then
       call dense(hw, K, args, info)
    else
       call msp_cite("leamer2024positivity", info)
       call evolve(hw, K, args, info)
    end if
  end subroutine occupy

  !> K = V f(D) V^H from the eigendecomposition hw = V D V^H.
  subroutine dense(hw, K, args, info)
    type(msp_matrix), intent(in) :: hw
    type(msp_matrix), intent(inout) :: K
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    type(msp_args) :: eigen_args
    type(msp_info) :: eigen_info
    type(msp_matrix) :: v
    real(dp), dimension(:), allocatable :: e, f
    real(dp) :: beta, mu, n_occupied, threshold, low, high
    integer :: filled, step
    logical :: finite

    finite = args%has("inverse_temperature")
    if (finite) then
       call args%get("inverse_temperature", beta)
       if (.not. beta > 0) call msp_raise("msp_fermi: " // &
            & "inverse_temperature must be positive")
    end if
    call args%get("threshold", threshold, default=0.0_dp, scoped=.true.)
    call args%enter("eigen", eigen_args)
    call info%enter("eigen", eigen_info)
    call eigen_impl(hw, e, eigen_args, eigen_info, vectors=v)
    allocate(f(size(e)))
    if (args%has("chemical_potential")) then
       call args%get("chemical_potential", mu)
       f = occupations(mu)
    else
       call args%get("trace", n_occupied)
       if (.not. (n_occupied >= 0 .and. n_occupied <= hw%n)) &
            & call msp_raise("msp_fermi: trace must lie between 0 and " // &
            & "the matrix size")
       if (finite) then
          !! The occupation falls below exp(-40) within 40 / beta of mu.
          low = e(1) - 40 / beta
          high = e(size(e)) + 40 / beta
          do step = 1, 200
             mu = (low + high) / 2
             if (sum(occupations(mu)) > n_occupied) then
                high = mu
             else
                low = mu
             end if
          end do
          f = occupations(mu)
       else
          filled = min(floor(n_occupied), size(e))
          f = 0
          f(1:filled) = 1
          if (filled < size(e)) f(filled + 1) = n_occupied - filled
          mu = e(max(filled, 1))
          if (filled < size(e)) mu = mu + (filled + 1 - n_occupied) * &
               & (e(filled + 1) - mu) / 2
       end if
    end if
    call eigen_compose(v, f, K, threshold)
    call info%set("energy", sum(f * e))
    call info%set("chemical_potential", mu)
    call v%free()

  contains

    !> The occupations of the eigenvalues e for the chemical potential.
    function occupations(potential) result(occ)
      real(dp), intent(in) :: potential
      real(dp), dimension(size(e)) :: occ
      real(dp) :: t
      integer :: q

      do q = 1, size(e)
         if (.not. finite) then
            occ(q) = merge(1.0_dp, 0.0_dp, e(q) < potential) + &
                 & merge(0.5_dp, 0.0_dp, e(q) == potential)
         else
            t = beta * (e(q) - potential)
            if (t > 0) then
               occ(q) = exp(-t) / (1 + exp(-t))
            else
               occ(q) = 1 / (1 + exp(t))
            end if
         end if
      end do
    end function occupations
  end subroutine dense

  !> K = W^2 at the inverse temperature, for hw in an orthogonal basis.
  subroutine evolve(hw, K, args, info)
    type(msp_matrix), intent(in) :: hw
    type(msp_matrix), intent(inout) :: K
    type(msp_args), intent(in) :: args
    type(msp_info), intent(in) :: info
    type(msp_matrix) :: a, w, k0, k1, euler, heun, id
    real(dp) :: beta, target, mu, n_occupied, step_tolerance, tolerance, &
         & threshold, step, error, energy
    integer :: steps, evaluations, max_steps
    logical :: canonical, converged

    call args%get("inverse_temperature", target)
    canonical = .not. args%has("chemical_potential")
    call args%get("step_tolerance", step_tolerance, default=1.0e-2_dp)
    call args%get("tolerance", tolerance, default=0.0_dp, scoped=.true.)
    call args%get("threshold", threshold, default=0.0_dp, scoped=.true.)
    call args%get("max_steps", max_steps, default=100000)
    call identity_like(hw, id)
    call copy(a, hw)
    call copy(w, id)
    if (canonical) then
       call args%get("trace", n_occupied)
       call msp_scale(w, sqrt(n_occupied / hw%n))
    else
       call args%get("chemical_potential", mu)
       call msp_add(a, id, alpha=-mu)
       call msp_scale(w, 1 / sqrt(2.0_dp))
    end if
    beta = 0
    step = 1
    steps = 0
    evaluations = 0
    converged = .true.
    do while (beta < target)
       if (steps == max_steps) then
          converged = .false.
          exit
       end if
       step = min(step, target - beta)
       call slope(w, a, id, canonical, threshold, k0)
       evaluations = evaluations + 1
       if (msp_norm(k0) <= tolerance) exit
       do
          call copy(euler, w)
          call msp_add(euler, k0, alpha=step, threshold=threshold)
          call slope(euler, a, id, canonical, threshold, k1)
          evaluations = evaluations + 1
          call copy(heun, w)
          call msp_add(heun, k0, alpha=step / 2, threshold=threshold)
          call msp_add(heun, k1, alpha=step / 2, threshold=threshold)
          call msp_add(euler, heun, alpha=-1.0_dp)
          error = msp_norm(euler)
          if (error <= 1.1_dp * step_tolerance) exit
          step = step * sqrt(step_tolerance / error)
       end do
       call msp_swap(w, heun)
       beta = beta + step
       steps = steps + 1
       step = step * min(2.0_dp, sqrt(step_tolerance / max(error, &
            & tiny(1.0_dp))))
    end do
    call msp_multiply(w, w, K, threshold=threshold)
    call msp_dot(K, hw, energy)
    call info%set("energy", energy)
    call info%set("inverse_temperature_reached", beta)
    call info%set("steps", steps)
    call info%set("evaluations", evaluations)
    call info%set("converged", converged)
    call a%free()
    call w%free()
    call k0%free()
    call k1%free()
    call euler%free()
    call heun%free()
    call id%free()
  end subroutine evolve

  !> s = -(X A - p X) / 2 with X = V - V^3, p = 0 (grand canonical) or
  !> <V, X A> / <X, V>, which makes <V, s> = 0 and keeps tr(V^2).
  subroutine slope(v, a, id, canonical, threshold, s)
    type(msp_matrix), intent(in) :: v, a, id
    logical, intent(in) :: canonical
    real(dp), intent(in) :: threshold
    type(msp_matrix), intent(inout) :: s
    type(msp_matrix) :: x
    real(dp) :: numerator, denominator

    call msp_multiply(v, v, s, threshold=threshold)
    call msp_scale(s, -1.0_dp)
    call msp_add(s, id)
    call msp_multiply(v, s, x, threshold=threshold)
    call msp_multiply(x, a, s, alpha=-0.5_dp, threshold=threshold)
    if (canonical) then
       call msp_dot(v, s, numerator)
       call msp_dot(x, v, denominator)
       call msp_add(s, x, alpha=-numerator / denominator, &
            & threshold=threshold)
    end if
    call x%free()
  end subroutine slope
end module msp_fermi_mod
