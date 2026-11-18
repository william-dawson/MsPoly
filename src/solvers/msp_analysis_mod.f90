!> Reduced descriptions of Hermitian matrices.
module msp_analysis_mod
  use msp_citations_mod, only: msp_cite
  use msp_kinds_mod, only: dp
  use msp_dict_mod, only: msp_dict, msp_args, msp_info, copy
  use msp_matrix_mod
  use msp_density_mod, only: msp_density
  use msp_linear_mod, only: msp_pivoted_cholesky
  implicit none
  private

  public :: msp_reduce_dimension

contains

  !> R becomes the dim x dim matrix L^H H L, where L L^H is the projector
  !> onto the dim lowest eigenvectors of the Hermitian H and L comes from
  !> its pivoted Cholesky factorization (dim columns): H restricted to that
  !> subspace, in the basis of the columns of L. Collective.
  !>
  !> args: dim (required), method (purification method, default "trs4"),
  !> threshold (scoped, default 0), tolerance (scoped, of the
  !> purification).
  !> info: density (info of the purification), cholesky.
  subroutine msp_reduce_dimension(H, R, args, info, context)
    type(msp_matrix), intent(in) :: H
    type(msp_matrix), intent(inout) :: R
    type(msp_dict), intent(in) :: args
    type(msp_dict), intent(inout), target :: info
    type(msp_dict), intent(in), target, optional :: context
    type(msp_dict), target :: tracked
    type(msp_dict) :: inner, inner_info
    type(msp_args) :: a
    type(msp_info) :: r_info
    type(msp_matrix) :: p, l, lh, work
    character(:), allocatable :: method
    real(dp) :: threshold, tolerance
    integer :: dim

    call copy(tracked, args)
    call a%open(tracked, context)
    call r_info%open(info)
    call msp_cite("aquilante2006fast", r_info)
    call a%get("dim", dim)
    call a%get("method", method, default="trs4")
    call a%get("threshold", threshold, default=0.0_dp, scoped=.true.)
    call a%get("tolerance", tolerance, default=1.0e-6_dp, scoped=.true.)
    call inner%set("trace", real(dim, dp))
    call inner%set("method", method)
    call inner%set("threshold", threshold)
    call inner%set("tolerance", tolerance)
    call msp_density(H, p, inner, inner_info)
    call info%set("density", inner_info)
    call inner%clear()
    call inner_info%clear()
    call inner%set("rank", dim)
    call msp_pivoted_cholesky(p, l, inner, inner_info)
    call info%set("cholesky", inner_info)
    call msp_transpose(l, lh, conjugate=.true.)
    call msp_multiply(lh, H, work, threshold=threshold)
    call msp_multiply(work, l, p, threshold=threshold)
    call leading(p, dim, R)
    call p%free()
    call l%free()
    call lh%free()
    call work%free()
    call a%finish(r_info)
  end subroutine msp_reduce_dimension

  !> B becomes the leading m x m block of A, on the same grid, with the
  !> kind, precision and location of A.
  subroutine leading(A, m, B)
    type(msp_matrix), intent(in) :: A
    integer, intent(in) :: m
    type(msp_matrix), intent(inout) :: B
    integer, dimension(:), allocatable :: i, j
    complex(dp), dimension(:), allocatable :: v
    integer :: q, kept

    call A%entries(i, j, v)
    kept = 0
    do q = 1, size(i)
       if (i(q) > m .or. j(q) > m) cycle
       kept = kept + 1
       i(kept) = i(q)
       j(kept) = j(q)
       v(kept) = v(q)
    end do
    call B%init(A%grid, m)
    if (A%is_complex()) then
       call B%fill(i(1:kept), j(1:kept), v(1:kept))
    else
       call B%fill(i(1:kept), j(1:kept), real(v(1:kept), dp))
    end if
    if (A%is_single()) call msp_to_single(B)
    if (A%on_device()) call msp_to_device(B)
  end subroutine leading
end module msp_analysis_mod
