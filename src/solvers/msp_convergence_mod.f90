!> Convergence monitoring for the fixed-point iterations (square roots,
!> inverses, roots, sign, purification).
module msp_convergence_mod
  use msp_kinds_mod, only: dp
  use msp_error_mod, only: msp_raise
  implicit none
  private

  public :: msp_monitor

  !> The convergence monitor of one iterative solve. The solver sets the
  !> policy at init and reports every residual to check, which returns
  !> the stop reason: "tolerance" when the residual is at most tolerance
  !> (a negative tolerance turns this test off), or "stagnation" when the
  !> residual has not fallen by the ratio over the window.
  !>
  !> Two policies cover the solvers. The strict policy (ratio 1, window
  !> 1) stops on the first step that does not improve the residual at
  !> all: it protects iterations that their filtering damages
  !> irreparably, where continuing past a non-improving step only
  !> amplifies truncation error. A ratio below one over a wider window
  !> detects the truncation floor of self-correcting iterations against
  !> a fixed target: there the Frobenius residual aggregates the noise
  !> of millions of filtered entries into a smooth crawl instead of the
  !> fluctuations a per-column maximum would show, so a strict
  !> non-improvement test never fires and the iteration grinds on the
  !> floor. The window matters: single-step ratios tight enough to catch
  !> the floor also fire during the slow phases of genuine convergence.
  !>
  !> The monitor is a policy, not a log: it keeps only the residuals of
  !> one window, while the solver reports the full residual history
  !> through its info dictionary.
  type :: msp_monitor
     private
     real(dp), allocatable :: r(:)
     integer :: filled = 0
     real(dp) :: tolerance = 0.0_dp
     real(dp) :: ratio = 1.0_dp
     integer :: window = 1
   contains
     procedure :: init
     procedure :: check
  end type msp_monitor

contains

  !> Set the policy: the residual target, and the progress requirement
  !> that ends the iteration when the target is out of reach. Reusing a
  !> monitor clears its window.
  subroutine init(this, tolerance, ratio, window)
    class(msp_monitor), intent(inout) :: this
    !> Residual target of the iteration; a negative tolerance disables
    !> the test.
    real(dp), intent(in), optional :: tolerance
    !> Required improvement over the window: the last residual must be
    !> below ratio times the residual window steps back. 1 (the default)
    !> ends the iteration on any non-improving step.
    real(dp), intent(in), optional :: ratio
    !> Steps over which the improvement is measured (default 1).
    integer, intent(in), optional :: window

    if (present(tolerance)) this%tolerance = tolerance
    if (present(ratio)) this%ratio = ratio
    if (present(window)) this%window = window
    if (this%ratio <= 0.0_dp) &
         & call msp_raise("msp_monitor: the ratio must be positive")
    if (this%window < 1) &
         & call msp_raise("msp_monitor: the window must be at least 1")
    if (allocated(this%r)) deallocate(this%r)
    allocate(this%r(this%window))
    this%filled = 0
  end subroutine init

  !> Report one residual and evaluate the policy: "" to continue,
  !> "tolerance" or "stagnation" to stop.
  function check(this, value) result(reason)
    class(msp_monitor), intent(inout) :: this
    real(dp), intent(in) :: value
    character(:), allocatable :: reason
    integer :: k

    if (value <= this%tolerance) then
       reason = "tolerance"
       return
    end if
    if (this%filled < this%window) then
       this%filled = this%filled + 1
       this%r(this%filled) = value
    else if (value >= this%ratio * this%r(1)) then
       reason = "stagnation"
       return
    else
       do k = 1, this%window - 1
          this%r(k) = this%r(k + 1)
       end do
       this%r(this%window) = value
    end if
    reason = ""
  end function check
end module msp_convergence_mod
