!> Matrix Market and binary files: real and complex matrices written and
!> read back, symmetric and Hermitian files expanded, and unreadable files
!> reported on every rank; permutations and entry maps. Run on 1, 2 and 4
!> ranks.
program test_io
  use mpi_f08, only: MPI_Init, MPI_Finalize, MPI_COMM_WORLD, MPI_Barrier
  use msp_kinds_mod, only: dp
  use msp_grid_mod, only: msp_grid
  use msp_matrix_mod
  use msp_io_mod
  use msp_maps_mod
  use spectral_mod, only: spectrum, distribute, deviation
  use testing_mod, only: check, finish
  implicit none

  type(msp_grid), target :: grid
  type(msp_matrix) :: a, b
  complex(dp), dimension(7, 7) :: d
  character(len=32) :: prefix
  character(len=64) :: header
  character(len=64), allocatable :: empty(:)
  character(:), allocatable :: msg
  integer :: stat, k, i
  integer, dimension(7) :: perm
  logical :: cplx, opened

  call MPI_Init()
  call grid%init(MPI_COMM_WORLD)
  write (prefix, '(a, i0, a)') "test_io_", grid%size, "_"

  do k = 1, 2
     cplx = k == 2
     d = spectrum([(real(i, dp), i = -3, 3)], cplx, 0.8_dp, 1.9_dp)
     call distribute(a, d, grid, cplx, .false.)
     call msp_write_matrix_market(a, trim(prefix) // "a.mtx")
     call msp_read_matrix_market(b, grid, trim(prefix) // "a.mtx")
     call check(deviation(b, d) < 1.0e-14_dp .and. (b%is_complex() .eqv. &
          & cplx), merge("complex", "real   ", cplx) // " Matrix Market")
     call msp_write_binary(a, trim(prefix) // "a.bin")
     call msp_read_binary(b, grid, trim(prefix) // "a.bin")
     call check(deviation(b, d) == 0 .and. (b%is_complex() .eqv. cplx), &
          & merge("complex", "real   ", cplx) // " binary")

     perm = msp_random_permutation(7, 424242)
     call msp_permute(a, perm, b)
     call check(deviation(b, permuted(d, perm)) == 0 .and. &
          & all([(count(perm == i) == 1, i = 1, 7)]) .and. &
          & any(perm /= [(i, i = 1, 7)]), &
          & merge("complex", "real   ", cplx) // " permutation")
     call msp_permute(b, perm, a, inverse=.true.)
     call check(deviation(a, d) == 0, &
          & merge("complex", "real   ", cplx) // " inverse permutation")
     call msp_map(a, b, upper)
     call check(deviation(b, merge(2 * d, (0.0_dp, 0.0_dp), &
          & upper_mask(7))) == 0 .and. (b%is_complex() .eqv. cplx), &
          & merge("complex", "real   ", cplx) // " entry map")
  end do

  call text("s.mtx", [character(64) :: &
       & "%%MatrixMarket matrix coordinate real symmetric", &
       & "% a comment                                     ", &
       & "3 3 4                                           ", &
       & "1 1 2.0                                         ", &
       & "2 1 -1.0                                        ", &
       & "3 2 0.5                                         ", &
       & "3 3 1.0                                         "])
  call msp_read_matrix_market(b, grid, trim(prefix) // "s.mtx")
  call check(deviation(b, reshape(cmplx([2.0_dp, -1.0_dp, 0.0_dp, -1.0_dp, &
       & 0.0_dp, 0.5_dp, 0.0_dp, 0.5_dp, 1.0_dp], 0.0_dp, dp), [3, 3])) == 0, &
       & "symmetric file")
  call text("h.mtx", [character(64) :: &
       & "%%MATRIXMARKET MATRIX COORDINATE COMPLEX HERMITIAN", &
       & "2 2 2                                             ", &
       & "1 1 1.0 0.0                                       ", &
       & "2 1 0.0 1.0                                       "])
  call msp_read_matrix_market(b, grid, trim(prefix) // "h.mtx")
  call check(deviation(b, reshape([(1.0_dp, 0.0_dp), (0.0_dp, 1.0_dp), &
       & (0.0_dp, -1.0_dp), (0.0_dp, 0.0_dp)], [2, 2])) == 0, &
       & "Hermitian file")

  call msp_read_matrix_market(b, grid, trim(prefix) // "missing.mtx", stat, &
       & msg)
  call check(stat /= 0 .and. index(msg, "cannot open") > 0, &
       & "missing Matrix Market file")
  call text("x.mtx", [character(64) :: &
       & "%%MatrixMarket matrix array real general", &
       & "2 2                                     "])
  call msp_read_matrix_market(b, grid, trim(prefix) // "x.mtx", stat, msg)
  call check(stat /= 0 .and. index(msg, "not a real or complex") > 0, &
       & "unsupported Matrix Market format")
  call text("r.mtx", [character(64) :: &
       & "%%MatrixMarket matrix coordinate real general", &
       & "2 3 1                                        ", &
       & "1 1 1.0                                      "])
  call msp_read_matrix_market(b, grid, trim(prefix) // "r.mtx", stat, msg)
  call check(stat /= 0 .and. index(msg, "not square") > 0, &
       & "rectangular Matrix Market file")
  call msp_read_binary(b, grid, trim(prefix) // "missing.bin", stat, msg)
  call check(stat /= 0 .and. index(msg, "cannot read") > 0, &
       & "missing binary file")
  if (grid%rank == 0) then
     open (newunit=k, file=trim(prefix) // "r.bin", access='stream', &
          & form='unformatted', status='replace')
     write (k) 2, 3, 0, 0_8
     close (k)
  end if
  call MPI_Barrier(MPI_COMM_WORLD)
  call msp_read_binary(b, grid, trim(prefix) // "r.bin", stat, msg)
  call check(stat /= 0 .and. index(msg, "not square") > 0, &
       & "rectangular binary file")

  inquire(file=trim(prefix) // "r.bin", opened=opened)
  call check(.not. opened, "binary error closes the input file")
  allocate(empty(0))
  call text("empty.mtx", empty)
  call msp_read_matrix_market(b, grid, trim(prefix) // "empty.mtx", stat, msg)
  call check(stat /= 0, "empty Matrix Market file")
  call text("shape.mtx", [character(64) :: &
       & "%%MatrixMarket matrix coordinate real general"])
  call msp_read_matrix_market(b, grid, trim(prefix) // "shape.mtx", stat, msg)
  call check(stat /= 0 .and. index(msg, "shape") > 0, &
       & "missing Matrix Market shape")
  call text("negative.mtx", [character(64) :: &
       & "%%MatrixMarket matrix coordinate real general", "2 2 -1"])
  call msp_read_matrix_market(b, grid, trim(prefix) // "negative.mtx", &
       & stat, msg)
  call check(stat /= 0 .and. index(msg, "invalid shape") > 0, &
       & "negative Matrix Market entry count")
  do i = 1, 2
     header = "%%MatrixMarket matrix coordinate real general"
     if (i == 2) header = &
          & "%%MatrixMarket matrix coordinate complex general"
     call text("short.mtx", [character(64) :: header, "2 2 1"])
     call msp_read_matrix_market(b, grid, trim(prefix) // "short.mtx", &
          & stat, msg)
     call check(stat /= 0 .and. index(msg, "entries") > 0, &
          & "truncated Matrix Market entries")
  end do
  call text("invalid.mtx", [character(64) :: &
       & "%%MatrixMarket matrix coordinate real general", &
       & "2 2 1", "3 1 1.0"])
  call msp_read_matrix_market(b, grid, trim(prefix) // "invalid.mtx", &
       & stat, msg)
  call check(stat /= 0, "Matrix Market fill errors return status")
  do i = 1, 5
     if (grid%rank == 0) then
        open (newunit=k, file=trim(prefix) // "bad.bin", access='stream', &
             & form='unformatted', status='replace')
        select case (i)
        case (1)
           write (k) 2
        case (2)
           write (k) 2, 2, 0, -1_8
        case (3:4)
           write (k) 2, 2, i - 3, 1_8
        case (5)
           write (k) 2, 2, 0, 1_8, 3, 1, 1.0_dp
        end select
        close (k)
     end if
     call MPI_Barrier(MPI_COMM_WORLD)
     call msp_read_binary(b, grid, trim(prefix) // "bad.bin", stat, msg)
     inquire(file=trim(prefix) // "bad.bin", opened=opened)
     call check(stat /= 0 .and. .not. opened, &
          & "malformed binary file returns status and closes the file")
  end do

  call a%free()
  call b%free()
  call grid%free()
  call finish("test_io")
  call MPI_Finalize()

contains

  !> d(perm, perm) placed so that entry (i, j) of d is at (perm(i), perm(j)).
  function permuted(d, perm) result(p)
    complex(dp), dimension(:, :), intent(in) :: d
    integer, dimension(:), intent(in) :: perm
    complex(dp), dimension(size(d, 1), size(d, 2)) :: p

    p(perm, perm) = d
  end function permuted

  !> Keep the upper triangle, doubled.
  logical function upper(i, j, v)
    integer, intent(in) :: i, j
    complex(dp), intent(inout) :: v

    v = 2 * v
    upper = i <= j
  end function upper

  function upper_mask(n) result(mask)
    integer, intent(in) :: n
    logical, dimension(n, n) :: mask
    integer :: p, q

    mask = reshape([((p <= q, p = 1, n), q = 1, n)], [n, n])
  end function upper_mask

  !> The first rank writes the lines to the file prefix // name.
  subroutine text(name, lines)
    character(*), intent(in) :: name
    character(*), dimension(:), intent(in) :: lines
    integer :: unit, k

    if (grid%rank == 0) then
       open (newunit=unit, file=trim(prefix) // name, status='replace')
       do k = 1, size(lines)
          write (unit, '(a)') trim(lines(k))
       end do
       close (unit)
    end if
    call MPI_Barrier(MPI_COMM_WORLD)
  end subroutine text
end program test_io
