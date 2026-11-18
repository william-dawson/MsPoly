!> Singular value decompositions through the polar decomposition.
module msp_svd_mod
  use msp_citations_mod, only: msp_cite
  use msp_kinds_mod, only: dp
  use msp_dict_mod, only: msp_dict, msp_args, msp_info, copy
  use msp_matrix_mod
  use msp_sign_mod, only: sign_impl
  use msp_eigen_mod, only: eigen_impl
  implicit none
  private

  public :: msp_svd

contains

  !> The singular values (ascending) of the nonsingular A = L S R^H, and
  !> optionally the left and right singular vectors as the columns of L and
  !> R, as NTPoly's SingularValueDecomposition: with the polar
  !> decomposition A = U H and the eigendecomposition H = R S R^H,
  !> L = U R. Collective.
  !>
  !> args: threshold (scoped, default 0), polar (args of the polar
  !> decomposition; see sign_impl), eigen (args of the eigendecomposition;
  !> see eigen_impl).
  !> info: polar, eigen.
  subroutine msp_svd(A, values, args, info, context, left, right)
    type(msp_matrix), intent(in) :: A
    real(dp), dimension(:), allocatable, intent(out) :: values
    type(msp_dict), intent(in) :: args
    type(msp_dict), intent(inout), target :: info
    type(msp_dict), intent(in), target, optional :: context
    type(msp_matrix), intent(inout), optional :: left, right
    type(msp_dict), target :: tracked
    type(msp_args) :: a_args, inner_args
    type(msp_info) :: r, inner_info
    type(msp_matrix) :: u, uh, h, v
    real(dp) :: threshold

    call copy(tracked, args)
    call a_args%open(tracked, context)
    call r%open(info)
    call msp_cite("higham2008functions", r)
    call a_args%get("threshold", threshold, default=0.0_dp, scoped=.true.)
    call a_args%enter("polar", inner_args)
    call r%enter("polar", inner_info)
    call sign_impl(A, inner_args, inner_info, u, .true.)
    call msp_transpose(u, uh, conjugate=.true.)
    call msp_multiply(uh, A, h, threshold=threshold)
    call a_args%enter("eigen", inner_args)
    call r%enter("eigen", inner_info)
    if (present(left) .or. present(right)) then
       call eigen_impl(h, values, inner_args, inner_info, vectors=v)
       if (present(left)) call msp_multiply(u, v, left, threshold=threshold)
       if (present(right)) call msp_swap(right, v)
    else
       call eigen_impl(h, values, inner_args, inner_info)
    end if
    call u%free()
    call uh%free()
    call h%free()
    call v%free()
    call a_args%finish(r)
  end subroutine msp_svd
end module msp_svd_mod
