!> Build the matrices a TRS2 purification multiplies, for multiply
!> benchmarks. A real Hamiltonian H is read from an NTPoly binary file. With
!> Gershgorin bounds [emin, emax], X = (emax I - H) / (emax - emin), then
!> each iteration replaces X by X^2 or 2X - X^2, whichever trace is closer
!> to the target, dropping entries below the threshold. Every iterate is
!> written in the same format as prefix_<k>.bin (k = 0 is the starting X).
!> Usage: make_filled file.bin threshold iterations target prefix [device]
program make_filled
  use mpi_f08, only: MPI_Init, MPI_Finalize, MPI_COMM_WORLD
  use msp_kinds_mod, only: dp
  use msp_dict_mod, only: msp_dict
  use msp_grid_mod, only: msp_grid
  use msp_matrix_mod
  use msp_bounds_mod, only: msp_bounds
  use msp_io_mod, only: msp_read_binary, msp_write_binary
  implicit none

  type(msp_grid), target :: grid
  type(msp_matrix) :: h, x, x2
  type(msp_dict) :: report, args, bounds
  character(len=256) :: file, prefix, arg, name
  integer :: iterations, k
  real(dp) :: threshold, target, emin, emax, tx, tx2

  call msp_select_gpu()
  call MPI_Init()
  call grid%init(MPI_COMM_WORLD)
  call get_command_argument(1, file)
  call get_command_argument(2, arg)
  read (arg, *) threshold
  call get_command_argument(3, arg)
  read (arg, *) iterations
  call get_command_argument(4, arg)
  read (arg, *) target
  call get_command_argument(5, prefix)
  call get_command_argument(6, arg)

  call msp_read_binary(h, grid, trim(file))
  call msp_bounds(h, args, bounds)
  emin = bounds%get_real("min")
  emax = bounds%get_real("max")
  call x%init(grid, h%n)
  call x%identity()
  call msp_scale(x, emax)
  call msp_add(x, h, alpha=-1.0_dp)
  call msp_scale(x, 1 / (emax - emin))
  call h%free()
  if (arg == "device") call msp_to_device(x)
  do k = 0, iterations
     if (k > 0) then
        call msp_multiply(x, x, x2, threshold=threshold)
        call msp_trace(x, tx)
        call msp_trace(x2, tx2)
        if (abs(tx2 - target) >= abs(2 * tx - tx2 - target)) then
           call msp_scale(x2, -1.0_dp)
           call msp_add(x2, x, alpha=2.0_dp, threshold=threshold)
        end if
        call copy(x, x2)
     end if
     write (name, '(a,"_",i0,".bin")') trim(prefix), k
     call msp_write_binary(x, trim(name))
     call msp_trace(x, tx)
     call report%append("nnz", x%nnz())
     call report%append("trace", tx)
  end do
  call report%set("gershgorin", [emin, emax])
  call grid%log(report)
  call x%free()
  call x2%free()
  call grid%free()
  call MPI_Finalize()

end program make_filled
