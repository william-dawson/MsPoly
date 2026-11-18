!> Paths that must stop the program with a message. Each run takes one
!> mode: dict | grid (assignment), raise (error without stat), agree
!> (collective error without stat).
program test_stops
  use mpi_f08, only: MPI_Init, MPI_COMM_WORLD
  use msp_dict_mod, only: msp_dict
  use msp_grid_mod, only: msp_grid
  implicit none

  character(len=8) :: mode
  type(msp_dict) :: a, b
  type(msp_grid) :: g, h
  integer :: i

  call get_command_argument(1, mode)
  select case (mode)
  case ("dict")
     call a%set("x", 1)
     b = a
  case ("raise")
     call a%get("missing", i)
  case ("grid")
     call MPI_Init()
     call g%init(MPI_COMM_WORLD)
     h = g
  case ("agree")
     call MPI_Init()
     call g%init(MPI_COMM_WORLD)
     if (g%rank == g%size - 1) call a%set("message", "boom")
     call g%agree(a)
  end select
  write (*, '(A)') "did not stop"
end program test_stops
