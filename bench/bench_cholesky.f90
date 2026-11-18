!> Usage: bench_cholesky file.bin repeats [rank] [device].
!> A positive rank selects pivoted Cholesky. Verification is not timed.
program bench_cholesky
  use mpi_f08, only: MPI_Init, MPI_Finalize, MPI_COMM_WORLD, MPI_Barrier, &
       & MPI_Wtime
  use msp_kinds_mod, only: dp
  use msp_dict_mod, only: msp_dict
  use msp_error_mod, only: msp_raise
  use msp_grid_mod, only: msp_grid
  use msp_matrix_mod
  use msp_linear_mod, only: msp_cholesky, msp_pivoted_cholesky
  use msp_io_mod, only: msp_read_binary
  implicit none
  type(msp_grid), target :: grid
  type(msp_matrix) :: a, l, lh, product
  type(msp_dict) :: args, info, report
  character(len=256) :: file, arg
  real(dp), allocatable :: seconds(:)
  real(dp) :: start, residual, norm_a
  integer :: repeats, rank, k, stat
  logical :: device

  call msp_select_gpu()
  call MPI_Init()
  call grid%init(MPI_COMM_WORLD)
  if (command_argument_count() < 2) call msp_raise( &
       & "usage: bench_cholesky file.bin repeats [rank] [device]")
  call get_command_argument(1, file)
  call get_command_argument(2, arg)
  read (arg, *, iostat=stat) repeats
  if (stat /= 0) call msp_raise("bench_cholesky: invalid repeats")
  if (repeats < 1) call msp_raise("bench_cholesky: repeats must be positive")
  rank = 0
  device = .false.
  do k = 3, command_argument_count()
     call get_command_argument(k, arg)
     if (arg == "device") then
        device = .true.
     else
        read (arg, *, iostat=stat) rank
        if (stat /= 0) call msp_raise("bench_cholesky: invalid rank")
        if (rank < 0) call msp_raise("bench_cholesky: negative rank")
     end if
  end do
  call msp_read_binary(a, grid, trim(file))
  if (device) call msp_to_device(a)
  if (rank > 0) call args%set("rank", rank)
  allocate(seconds(0:repeats))
  do k = 0, repeats
     call info%clear()
     call MPI_Barrier(grid%comm)
     start = MPI_Wtime()
     if (rank > 0) then
        call msp_pivoted_cholesky(a, l, args, info)
     else
        call msp_cholesky(a, l, args, info)
     end if
     call MPI_Barrier(grid%comm)
     seconds(k) = MPI_Wtime() - start
  end do
  call msp_transpose(l, lh, conjugate=.true.)
  call msp_multiply(l, lh, product)
  call msp_add(product, a, alpha=-1.0_dp)
  residual = msp_norm(product)
  norm_a = msp_norm(a)
  call report%set("input", trim(file))
  call report%set("n", a%n)
  call report%set("grid", [grid%rows, grid%cols])
  call report%set("device", device)
  call report%set("requested_rank", rank)
  call report%set("nnz_a", a%nnz())
  call report%set("nnz_l", l%nnz())
  call report%set("first_seconds", seconds(0))
  call report%set("warm_seconds", seconds(1:))
  call report%set("relative_residual", residual / norm_a)
  call report%set("solver", info)
  call grid%log(report)
  call a%free()
  call l%free()
  call lh%free()
  call product%free()
  call grid%free()
  call MPI_Finalize()
end program bench_cholesky
