!> Dense references for solver tests: matrices with a known spectrum and
!> eigenvectors, distributed copies of dense matrices, and the largest
!> difference between a distributed and a dense matrix.
module spectral_mod
  use mpi_f08, only: MPI_Allreduce, MPI_IN_PLACE, MPI_DOUBLE_PRECISION, &
       & MPI_MAX
  use msp_kinds_mod, only: dp
  use msp_grid_mod, only: msp_grid
  use msp_matrix_mod, only: msp_matrix, msp_to_device
  implicit none
  private

  public :: spectrum, distribute, deviation

contains

  !> U diag(d) U^H, with U the Householder reflector of the vector with
  !> entries cos(a i) (+ i sin(b i) when complex): a Hermitian matrix with
  !> eigenvalues d and known eigenvectors, so f(U D U^H) = U f(D) U^H.
  function spectrum(d, cplx, a, b) result(m)
    real(dp), dimension(:), intent(in) :: d
    logical, intent(in) :: cplx
    real(dp), intent(in) :: a, b
    complex(dp), dimension(size(d), size(d)) :: m
    complex(dp), dimension(size(d), size(d)) :: u
    complex(dp), dimension(size(d)) :: v
    integer :: i, n

    n = size(d)
    do i = 1, n
       v(i) = cmplx(cos(a * i), 0.0_dp, dp)
       if (cplx) v(i)%im = sin(b * i)
    end do
    u = 0
    do i = 1, n
       u(i, i) = 1
    end do
    u = u - 2 * spread(v, 2, n) * spread(conjg(v), 1, n) / sum(abs(v)**2)
    m = 0
    do i = 1, n
       m(i, i) = d(i)
    end do
    m = matmul(u, matmul(m, conjg(transpose(u))))
  end function spectrum

  !> m becomes the dense d on grid, with complex values when cplx (the real
  !> part otherwise), on the GPU when device. Each rank contributes part of
  !> the entries.
  subroutine distribute(m, d, grid, cplx, device)
    type(msp_matrix), intent(inout) :: m
    complex(dp), dimension(:, :), intent(in) :: d
    type(msp_grid), intent(in), target :: grid
    logical, intent(in) :: cplx, device
    integer, dimension(size(d)) :: i, j
    complex(dp), dimension(size(d)) :: v
    integer :: p, q, count

    count = 0
    do q = 1, size(d, 2)
       do p = 1, size(d, 1)
          if (mod(p + q, grid%size) /= grid%rank) cycle
          count = count + 1
          i(count) = p
          j(count) = q
          v(count) = d(p, q)
       end do
    end do
    call m%init(grid, size(d, 1))
    if (cplx) then
       call m%fill(i(1:count), j(1:count), v(1:count))
    else
       call m%fill(i(1:count), j(1:count), real(v(1:count), dp))
    end if
    if (device) call msp_to_device(m)
  end subroutine distribute

  !> Largest difference over all ranks between m and d.
  real(dp) function deviation(m, d) result(e)
    type(msp_matrix), intent(in) :: m
    complex(dp), dimension(:, :), intent(in) :: d
    complex(dp), dimension(:, :), allocatable :: block
    integer, dimension(:), allocatable :: i, j
    complex(dp), dimension(:), allocatable :: v
    integer :: k

    allocate(block(m%rows, m%cols))
    block = d(m%row0 + 1:m%row0 + m%rows, m%col0 + 1:m%col0 + m%cols)
    call m%entries(i, j, v)
    do k = 1, size(i)
       block(i(k) - m%row0, j(k) - m%col0) = &
            & block(i(k) - m%row0, j(k) - m%col0) - v(k)
    end do
    e = 0
    if (size(block) > 0) e = maxval(abs(block))
    call MPI_Allreduce(MPI_IN_PLACE, e, 1, MPI_DOUBLE_PRECISION, MPI_MAX, &
         & m%grid%comm)
  end function deviation
end module spectral_mod
