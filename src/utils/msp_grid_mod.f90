!> Two-dimensional process grid, collective error agreement and logging.
module msp_grid_mod
  use, intrinsic :: iso_fortran_env, only: int8, output_unit, error_unit
  use mpi_f08, only: MPI_Comm, MPI_COMM_NULL, MPI_INTEGER, MPI_BYTE, &
       & MPI_SUM, MPI_MIN, MPI_Comm_dup, MPI_Comm_split, MPI_Comm_free, &
       & MPI_Comm_rank, MPI_Comm_size, MPI_Allreduce, MPI_Bcast, &
       & MPI_Finalized, MPI_IN_PLACE, operator(/=)
  use msp_error_mod, only: msp_raise, msp_clear
  use msp_dict_mod, only: msp_dict
  implicit none
  private

  public :: msp_grid

  !> Ranks of a communicator arranged as rows x cols, row-major: rank =
  !> row * cols + col. The grid owns duplicates of its communicators;
  !> release them with free. A grid cannot be copied.
  type :: msp_grid
     !> All ranks of the grid.
     type(MPI_Comm) :: comm = MPI_COMM_NULL
     !> Ranks in the same grid row, ordered by column.
     type(MPI_Comm) :: row_comm = MPI_COMM_NULL
     !> Ranks in the same grid column, ordered by row.
     type(MPI_Comm) :: col_comm = MPI_COMM_NULL
     integer :: rank = 0
     integer :: size = 1
     integer :: rows = 1
     integer :: cols = 1
     integer :: row = 0
     integer :: col = 0
   contains
     procedure :: init_comm
     procedure :: init_int
     generic :: init => init_comm, init_int
     procedure :: free
     procedure :: agree
     procedure :: bcast
     procedure :: log
     generic :: assignment(=) => no_copy
     procedure, private :: no_copy
  end type msp_grid

contains

  !> Build the grid over the ranks of comm (collective). Without rows, the
  !> grid is 1D: one row band per rank. The gather-based multiply balances
  !> band-diagonal matrices only on this decomposition, and the total
  !> multiply work is shape-invariant, so 1D is the right default; pass
  !> rows explicitly (or escalate inside callers, like NTPoly's process
  !> slices) when measurements show gather costs overtaking compute.
  subroutine init_comm(this, comm, rows, stat, errmsg)
    class(msp_grid), intent(inout) :: this
    type(MPI_Comm), intent(in) :: comm
    integer, intent(in), optional :: rows
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg
    integer :: n, r

    call msp_clear(stat)
    call this%free()
    call MPI_Comm_size(comm, n)
    if (present(rows)) then
       r = rows
       if (r < 1 .or. mod(n, max(r, 1)) /= 0) then
          call msp_raise("msp_grid: " // str(r) // " rows do not divide " &
               & // str(n) // " ranks", stat, errmsg)
          return
       end if
    else
       !! One row band per rank.
       r = n
    end if
    call MPI_Comm_dup(comm, this%comm)
    call MPI_Comm_rank(this%comm, this%rank)
    this%size = n
    this%rows = r
    this%cols = n / r
    this%row = this%rank / this%cols
    this%col = mod(this%rank, this%cols)
    call MPI_Comm_split(this%comm, this%row, this%col, this%row_comm)
    call MPI_Comm_split(this%comm, this%col, this%row, this%col_comm)
  end subroutine init_comm

  !> Same as init, for a communicator given as an integer handle in the
  !> numbering of the mpi module, for embedding in codes and middlewares
  !> that do not use mpi_f08 (ELSI, for example). The integer handle and
  !> the mpi_f08 communicator designate the same communicator in the
  !> implementations this supports.
  subroutine init_int(this, comm, rows, stat, errmsg)
    class(msp_grid), intent(inout) :: this
    integer, intent(in) :: comm
    integer, intent(in), optional :: rows
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg
    type(MPI_Comm) :: comm_f08

    comm_f08%MPI_VAL = int(comm, kind=kind(comm_f08%MPI_VAL))
    call init_comm(this, comm_f08, rows, stat, errmsg)
  end subroutine init_int

  !> Release the communicators (collective). Safe to call more than once,
  !> and after MPI_Finalize.
  subroutine free(this)
    class(msp_grid), intent(inout) :: this
    logical :: finalized

    call MPI_Finalized(finalized)
    if (.not. finalized) then
       if (this%row_comm /= MPI_COMM_NULL) call MPI_Comm_free(this%row_comm)
       if (this%col_comm /= MPI_COMM_NULL) call MPI_Comm_free(this%col_comm)
       if (this%comm /= MPI_COMM_NULL) call MPI_Comm_free(this%comm)
    end if
    this%row_comm = MPI_COMM_NULL
    this%col_comm = MPI_COMM_NULL
    this%comm = MPI_COMM_NULL
    this%rank = 0
    this%size = 1
    this%rows = 1
    this%cols = 1
    this%row = 0
    this%col = 0
  end subroutine free

  subroutine no_copy(lhs, rhs)
    class(msp_grid), intent(inout) :: lhs
    class(msp_grid), intent(in) :: rhs

    if (lhs%size < 0 .or. rhs%size < 0) continue
    error stop "msp_grid: a grid cannot be copied"
  end subroutine no_copy

  !> Agree on the outcome of a step (collective). Each rank passes its own
  !> error record in err: empty if the step succeeded, otherwise at least a
  !> "message" and any details. If any rank failed, every rank returns the
  !> record of the lowest failing rank, with "rank" and "failed_ranks"
  !> added. Then stat is set, or, without stat, root writes the record to
  !> standard error and every rank stops.
  subroutine agree(this, err, stat, errmsg)
    class(msp_grid), intent(in) :: this
    type(msp_dict), intent(inout) :: err
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg
    integer :: failed, first

    call msp_clear(stat)
    failed = merge(1, 0, err%num_keys() > 0)
    call MPI_Allreduce(MPI_IN_PLACE, failed, 1, MPI_INTEGER, MPI_SUM, &
         & this%comm)
    if (failed == 0) return
    first = merge(this%rank, this%size, err%num_keys() > 0)
    call MPI_Allreduce(MPI_IN_PLACE, first, 1, MPI_INTEGER, MPI_MIN, &
         & this%comm)
    call this%bcast(err, first)
    call err%set("rank", first)
    call err%set("failed_ranks", failed)
    if (present(stat)) then
       call msp_raise("rank " // str(first) // ": " // &
            & err%get_string("message", "error"), stat, errmsg)
       return
    end if
    if (this%rank == 0) then
       write (error_unit, '(A)') "msp_grid: error"
       call err%write_yaml(error_unit, indent=2)
       flush (error_unit)
    end if
    error stop 1
  end subroutine agree

  !> Replace d on every rank with root's d (collective).
  subroutine bcast(this, d, root)
    class(msp_grid), intent(in) :: this
    type(msp_dict), intent(inout) :: d
    integer, intent(in) :: root
    integer(int8), dimension(:), allocatable :: bytes
    integer :: n

    if (this%rank == root) then
       call d%pack(bytes)
       n = size(bytes)
    end if
    call MPI_Bcast(n, 1, MPI_INTEGER, root, this%comm)
    if (this%rank /= root) allocate(bytes(n))
    call MPI_Bcast(bytes, n, MPI_BYTE, root, this%comm)
    if (this%rank /= root) call d%unpack(bytes)
  end subroutine bcast

  !> Write d on root as one YAML document, then flush, so a log file is
  !> complete up to the last record even if the run is killed.
  subroutine log(this, d, unit)
    class(msp_grid), intent(in) :: this
    type(msp_dict), intent(in) :: d
    !> Default: standard output.
    integer, intent(in), optional :: unit
    integer :: u

    if (this%rank /= 0) return
    u = output_unit
    if (present(unit)) u = unit
    write (u, '(A)') "---"
    call d%write_yaml(u)
    flush (u)
  end subroutine log

  pure function str(i) result(s)
    integer, intent(in) :: i
    character(:), allocatable :: s
    character(len=12) :: buffer

    write (buffer, '(I0)') i
    s = trim(buffer)
  end function str
end module msp_grid_mod
