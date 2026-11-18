!> Time dense eigendecompositions of a matrix read from an NTPoly binary
!> file: the first call, then the median of the following calls.
!> Usage: bench_eigen file.bin repeats [device] [backend=name] [values]
!> With values, only the eigenvalues are kept (the eigenvectors are still
!> computed, but not moved back into a sparse matrix).
program bench_eigen
  use mpi_f08, only: MPI_Init, MPI_Finalize, MPI_COMM_WORLD, MPI_Barrier, &
       & MPI_Wtime
  use msp_kinds_mod, only: dp
  use msp_dict_mod, only: msp_dict
  use msp_error_mod, only: msp_raise
  use msp_grid_mod, only: msp_grid
  use msp_matrix_mod
  use msp_eigen_mod, only: msp_eigen, msp_dense_release
  use msp_io_mod, only: msp_read_binary
  implicit none

  type(msp_grid), target :: grid
  type(msp_matrix) :: h, vectors
  type(msp_dict) :: args, info, report
  character(len=256) :: file, arg, backend
  real(dp), dimension(:), allocatable :: values, t
  real(dp) :: start
  integer :: repeats, k, stat
  logical :: device, only_values

  call msp_select_gpu()
  call MPI_Init()
  call grid%init(MPI_COMM_WORLD)
  if (command_argument_count() < 2) call msp_raise( &
       & "usage: bench_eigen file.bin repeats [device] [backend=name] " // &
       & "[values]")
  call get_command_argument(1, file)
  call get_command_argument(2, arg)
  read (arg, *, iostat=stat) repeats
  if (stat /= 0) call msp_raise("bench_eigen: invalid repeats")
  if (repeats < 1) call msp_raise("bench_eigen: repeats must be positive")
  device = .false.
  only_values = .false.
  backend = ""
  do k = 3, command_argument_count()
     call get_command_argument(k, arg)
     device = device .or. arg == "device"
     only_values = only_values .or. arg == "values"
     if (arg(1:8) == "backend=" .and. len_trim(arg) > 8) then
        backend = arg(9:)
     else if (arg /= "device" .and. arg /= "values") then
        call msp_raise("bench_eigen: unknown option '" // trim(arg) // "'")
     end if
  end do
  call msp_read_binary(h, grid, trim(file))
  if (device) call msp_to_device(h)
  if (len_trim(backend) > 0) call args%set("backend", trim(backend))
  allocate(t(0:repeats))
  do k = 0, repeats
     call info%clear()
     call MPI_Barrier(MPI_COMM_WORLD)
     start = MPI_Wtime()
     if (only_values) then
        call msp_eigen(h, values, args, info)
     else
        call msp_eigen(h, values, args, info, vectors=vectors)
     end if
     call MPI_Barrier(MPI_COMM_WORLD)
     t(k) = MPI_Wtime() - start
  end do

  call report%set("n", h%n)
  call report%set("input", trim(file))
  call report%set("repeats", repeats)
  call report%set("only_values", only_values)
  call report%set("grid", [grid%rows, grid%cols])
  call report%set("device", device)
  call report%set("nnz_h", h%nnz())
  if (.not. only_values) call report%set("nnz_vectors", vectors%nnz())
  call report%set("lowest", values(1))
  call report%set("highest", values(size(values)))
  call report%set("sum", sum(values))
  call report%set("sum_squares", sum(values**2))
  call report%set("first_seconds", t(0))
  call report%set("median_seconds", median(t(1:)))
  call report%set("warm_seconds", t(1:))
  call report%set("solver", info)
  call grid%log(report)
  call h%free()
  call vectors%free()
  call msp_dense_release()
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
    if (mod(size(y), 2) == 0) &
         & median = (y(size(y) / 2) + y(size(y) / 2 + 1)) / 2
  end function median
end program bench_eigen
