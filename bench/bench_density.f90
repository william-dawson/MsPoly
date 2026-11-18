!> Time TRS2 density matrices of a real Hamiltonian read from an NTPoly
!> binary file: the first call, then the median of the following calls.
!> Usage: bench_density file.bin trace threshold tolerance repeats [device]
!>        [single] [method=name] [out=density.bin]
!> With single, H is stored in single precision (and so the iterates); with
!> out, the last density matrix is written in NTPoly's binary format.
program bench_density
  use mpi_f08, only: MPI_Init, MPI_Finalize, MPI_COMM_WORLD, MPI_Barrier, &
       & MPI_Wtime
  use msp_kinds_mod, only: dp
  use msp_dict_mod, only: msp_dict
  use msp_grid_mod, only: msp_grid
  use msp_matrix_mod
  use msp_density_mod, only: msp_density
  use msp_io_mod, only: msp_read_binary, msp_write_binary
  implicit none

  type(msp_grid), target :: grid
  type(msp_matrix) :: h, density
  type(msp_dict) :: args, info, report
  character(len=256) :: file, arg, out, method
  integer :: repeats, k
  real(dp) :: trace, threshold, tolerance, start
  real(dp), dimension(:), allocatable :: t
  logical :: device, single

  call msp_select_gpu()
  call MPI_Init()
  call grid%init(MPI_COMM_WORLD)
  call get_command_argument(1, file)
  call get_command_argument(2, arg)
  read (arg, *) trace
  call get_command_argument(3, arg)
  read (arg, *) threshold
  call get_command_argument(4, arg)
  read (arg, *) tolerance
  call get_command_argument(5, arg)
  read (arg, *) repeats

  device = .false.
  single = .false.
  out = ""
  method = "trs2"
  do k = 6, command_argument_count()
     call get_command_argument(k, arg)
     device = device .or. arg == "device"
     single = single .or. arg == "single"
     if (arg(1:4) == "out=") out = arg(5:)
     if (arg(1:7) == "method=") method = arg(8:)
  end do
  call msp_read_binary(h, grid, trim(file))
  if (single) call msp_to_single(h)
  if (device) call msp_to_device(h)
  call args%set("trace", trace)
  call args%set("threshold", threshold)
  call args%set("tolerance", tolerance)
  call args%set("method", trim(method))
  allocate(t(0:repeats))
  do k = 0, repeats
     call info%clear()
     call MPI_Barrier(MPI_COMM_WORLD)
     start = MPI_Wtime()
     call msp_density(h, density, args, info)
     call MPI_Barrier(MPI_COMM_WORLD)
     t(k) = MPI_Wtime() - start
  end do

  call report%set("n", h%n)
  call report%set("grid", [grid%rows, grid%cols])
  call report%set("device", device)
  call report%set("single", single)
  call report%set("nnz_h", h%nnz())
  call report%set("nnz_density", density%nnz())
  call report%set("first_seconds", t(0))
  call report%set("median_seconds", median(t(1:)))
  call report%set("solver", info)
  call grid%log(report)
  if (len_trim(out) > 0) call msp_write_binary(density, trim(out))
  call h%free()
  call density%free()
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
end program bench_density
