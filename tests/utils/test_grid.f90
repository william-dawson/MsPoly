!> Process grid layout, collective error agreement, dictionary broadcast
!> and logging. Run on 1, 2 and 4 ranks.
program test_grid
  use mpi_f08, only: MPI_Init, MPI_Finalize, MPI_COMM_WORLD, MPI_COMM_NULL, &
       & MPI_Comm_size, MPI_Comm_rank, MPI_Bcast, MPI_INTEGER8
  use msp_kinds_mod, only: dp, i64
  use msp_dict_mod
  use msp_grid_mod
  use testing_mod, only: check, finish
  implicit none

  integer :: nranks

  call MPI_Init()
  call MPI_Comm_size(MPI_COMM_WORLD, nranks)
  call test_layout()
  call test_rows()
  call test_agree()
  call test_bcast()
  call test_log()
  call test_free()
  call finish("test_grid")
  call MPI_Finalize()

contains

  subroutine test_layout()
    type(msp_grid) :: g
    integer :: row_size, col_size, row_rank, col_rank

    call g%init(MPI_COMM_WORLD)
    call check(g%rows * g%cols == nranks .and. g%rows == nranks .and. &
         & g%cols == 1, "default grid is 1D: one row band per rank")
    call MPI_Comm_size(g%row_comm, row_size)
    call MPI_Comm_size(g%col_comm, col_size)
    call MPI_Comm_rank(g%row_comm, row_rank)
    call MPI_Comm_rank(g%col_comm, col_rank)
    call check(row_size == g%cols .and. col_size == g%rows .and. &
         & row_rank == g%col .and. col_rank == g%row, &
         & "row and column communicators")
    call check(g%rank == g%row * g%cols + g%col, "row-major rank layout")
  end subroutine test_layout

  subroutine test_rows()
    type(msp_grid) :: g
    integer :: stat
    character(:), allocatable :: msg

    call g%init(MPI_COMM_WORLD, rows=nranks)
    call check(g%rows == nranks .and. g%cols == 1, "explicit row count")
    call g%init(MPI_COMM_WORLD, rows=nranks + 1, stat=stat, errmsg=msg)
    call check(stat /= 0 .and. index(msg, "do not divide") > 0, &
         & "a row count that does not divide the ranks is an error")
  end subroutine test_rows

  subroutine test_agree()
    type(msp_grid) :: g
    type(msp_dict) :: err
    integer :: stat
    character(:), allocatable :: msg
    character(len=8) :: own

    call g%init(MPI_COMM_WORLD)
    call g%agree(err, stat)
    call check(stat == 0 .and. err%num_keys() == 0, "no failure, no error")

    if (g%rank == g%size - 1) then
       call err%set("message", "bad value")
       call err%set("key", "bounds.tol")
    end if
    call g%agree(err, stat, msg)
    call check(stat /= 0 .and. index(msg, "bad value") > 0 .and. &
         & err%get_int("rank") == g%size - 1 .and. &
         & err%get_string("key") == "bounds.tol" .and. &
         & err%get_int("failed_ranks") == 1, &
         & "every rank receives one rank's error record")

    call err%clear()
    write (own, '(A,I0)') "rank", g%rank
    call err%set("message", trim(own))
    call g%agree(err, stat, msg)
    call check(err%get_string("message") == "rank0" .and. &
         & err%get_int("failed_ranks") == g%size, &
         & "the lowest failing rank's record wins")
  end subroutine test_agree

  subroutine test_bcast()
    type(msp_grid) :: g
    type(msp_dict) :: d
    integer(i64) :: expected
    integer :: root

    call g%init(MPI_COMM_WORLD)
    root = g%size - 1
    if (g%rank == root) then
       call d%set("x", 0.5_dp)
       call d%set("s.list", ["ab", "cd"])
       call d%set("s.flags", [.true., .false.])
       expected = d%hash()
    else
       call d%set("junk", 1)
    end if
    call MPI_Bcast(expected, 1, MPI_INTEGER8, root, MPI_COMM_WORLD)
    call g%bcast(d, root)
    call check(d%hash() == expected .and. .not. d%has("junk"), &
         & "broadcast replaces the dictionary with root's")
  end subroutine test_bcast

  subroutine test_log()
    type(msp_grid) :: g
    type(msp_dict) :: d
    character(len=40) :: line
    integer :: unit, lines, ios

    call g%init(MPI_COMM_WORLD)
    call d%set("iterations", 3)
    call d%set("bounds.converged", .true.)
    open (newunit=unit, status='scratch', action='readwrite')
    call g%log(d, unit)
    rewind (unit)
    lines = 0
    do
       read (unit, '(A)', iostat=ios) line
       if (ios /= 0) exit
       lines = lines + 1
    end do
    close (unit)
    call check(lines == merge(4, 0, g%rank == 0), &
         & "root writes one YAML document; other ranks write nothing")
  end subroutine test_log

  subroutine test_free()
    type(msp_grid) :: g

    call g%init(MPI_COMM_WORLD)
    call g%free()
    call g%free()
    call check(g%comm%MPI_VAL == MPI_COMM_NULL%MPI_VAL .and. g%size == 1, &
         & "free releases the grid and can be repeated")
    call g%init(MPI_COMM_WORLD)
    call check(g%size == nranks, "a freed grid can be built again")
  end subroutine test_free
end program test_grid
