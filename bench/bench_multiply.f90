!> Time C = A * A with a threshold: the whole multiply (median and mean of
!> the timed calls, after one warm-up call), then its stages (row band
!> gather, column band gather, local product) separately. A is either a
!> banded n x n matrix with entries exp(-|i - j| / 32) for |i - j| <= width,
!> or read from an NTPoly binary file.
!> Usage: bench_multiply n width threshold repeats [device] [single]
!>        bench_multiply file.bin threshold repeats [device] [single]
!> With single, A is stored in single precision.
program bench_multiply
  use mpi_f08, only: MPI_Init, MPI_Finalize, MPI_COMM_WORLD, MPI_Barrier, &
       & MPI_Wtime
  use msp_kinds_mod, only: dp
  use msp_dict_mod, only: msp_dict
  use msp_grid_mod, only: msp_grid
  use msp_storage_mod, only: msp_storage, msp_free
  use msp_matrix_mod
  use msp_io_mod, only: msp_read_binary
  implicit none

  type(msp_grid), target :: grid
  type(msp_matrix) :: a, c
  type(msp_dict) :: report
  character(len=256) :: arg, file
  integer :: n, width, repeats, k, first
  real(dp) :: threshold, start, trace
  real(dp), dimension(:), allocatable :: t
  real(dp), dimension(3) :: stage
  class(msp_storage), allocatable :: a_band, b_band
  logical :: device, single

  call msp_select_gpu()
  call MPI_Init()
  call grid%init(MPI_COMM_WORLD)
  call get_command_argument(1, arg)
  k = 2
  file = ""
  if (verify(trim(arg), "0123456789") /= 0) then
     file = arg
  else
     read (arg, *) n
     call get_command_argument(2, arg)
     read (arg, *) width
     k = 3
  end if
  call get_command_argument(k, arg)
  read (arg, *) threshold
  call get_command_argument(k + 1, arg)
  read (arg, *) repeats
  device = .false.
  single = .false.
  first = k + 2
  do k = first, command_argument_count()
     call get_command_argument(k, arg)
     device = device .or. arg == "device"
     single = single .or. arg == "single"
  end do

  if (len_trim(file) > 0) then
     call msp_read_binary(a, grid, trim(file))
     n = a%n
  else
     call build(a)
  end if
  if (single) call msp_to_single(a)
  if (device) call msp_to_device(a)
  call msp_multiply(a, a, c, threshold=threshold)
  allocate(t(repeats))
  do k = 1, repeats
     call MPI_Barrier(MPI_COMM_WORLD)
     start = MPI_Wtime()
     call msp_multiply(a, a, c, threshold=threshold)
     call MPI_Barrier(MPI_COMM_WORLD)
     t(k) = MPI_Wtime() - start
  end do
  call msp_trace(c, trace)

  stage = 0
  do k = 1, repeats
     call timed(1)
     call a%local%gather_rows(grid%row_comm, a_band)
     call timed(1)
     call timed(2)
     call a%local%gather_cols(grid%col_comm, b_band)
     call timed(2)
     call timed(3)
     call c%local%multiply(a_band, b_band, (1.0_dp, 0.0_dp), (0.0_dp, 0.0_dp), &
          & threshold)
     call timed(3)
  end do

  call report%set("n", n)
  call report%set("ranks", grid%size)
  call report%set("grid", [grid%rows, grid%cols])
  call report%set("device", device)
  call report%set("single", single)
  call report%set("nnz_a", a%nnz())
  call report%set("nnz_c", c%nnz())
  call report%set("trace_c", trace)
  call report%set("multiply_seconds", sum(t) / repeats)
  call report%set("multiply_median_seconds", median(t))
  if (len_trim(file) > 0) call report%set("input", trim(file))
  call report%set("stage_seconds.gather_rows", stage(1) / repeats)
  call report%set("stage_seconds.gather_cols", stage(2) / repeats)
  call report%set("stage_seconds.local_product", stage(3) / repeats)
  call grid%log(report)
  call msp_free(a_band)
  call msp_free(b_band)
  call a%free()
  call c%free()
  call grid%free()
  call MPI_Finalize()

contains

  real(dp) function median(x)
    real(dp), dimension(:), intent(in) :: x
    real(dp), dimension(size(x)) :: y
    integer :: i, j
    real(dp) :: v

    y = x
    do i = 2, size(y)
       v = y(i)
       j = i - 1
       do while (j >= 1)
          if (y(j) <= v) exit
          y(j + 1) = y(j)
          j = j - 1
       end do
       y(j + 1) = v
    end do
    median = y((size(y) + 1) / 2)
  end function median

  !> First call starts stage s, second call adds its elapsed time.
  subroutine timed(s)
    integer, intent(in) :: s
    real(dp), dimension(3), save :: started = -1

    call MPI_Barrier(MPI_COMM_WORLD)
    if (started(s) < 0) then
       started(s) = MPI_Wtime()
    else
       stage(s) = stage(s) + MPI_Wtime() - started(s)
       started(s) = -1
    end if
  end subroutine timed

  !> Each rank generates the rows of its row band, columns of its column band.
  subroutine build(m)
    type(msp_matrix), intent(inout) :: m
    integer, dimension(:), allocatable :: i, j
    real(dp), dimension(:), allocatable :: v
    integer :: r, col, count

    call m%init(grid, n)
    allocate(i(m%rows * (2 * width + 1)), j(m%rows * (2 * width + 1)), &
         & v(m%rows * (2 * width + 1)))
    count = 0
    do r = m%row0 + 1, m%row0 + m%rows
       do col = max(r - width, m%col0 + 1), min(r + width, m%col0 + m%cols)
          count = count + 1
          i(count) = r
          j(count) = col
          v(count) = exp(-abs(r - col) / 32.0_dp)
       end do
    end do
    call m%fill(i(1:count), j(1:count), v(1:count))
  end subroutine build
end program bench_multiply
