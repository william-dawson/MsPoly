!> Reading and writing distributed matrices: Matrix Market coordinate files
!> and NTPoly's binary format.
module msp_io_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64
  use mpi_f08, only: MPI_Barrier
  use msp_kinds_mod, only: dp
  use msp_error_mod, only: msp_clear
  use msp_dict_mod, only: msp_dict
  use msp_grid_mod, only: msp_grid
  use msp_matrix_mod, only: msp_matrix
  implicit none
  private

  public :: msp_read_matrix_market, msp_write_matrix_market, &
       & msp_read_binary, msp_write_binary

contains

  !> A becomes the square matrix of the Matrix Market coordinate file name
  !> (real or complex; general, symmetric or Hermitian, the other triangle
  !> filled in). The first rank reads the file. Collective; errors are
  !> agreed on all ranks.
  subroutine msp_read_matrix_market(A, grid, name, stat, errmsg)
    type(msp_matrix), intent(inout) :: A
    type(msp_grid), intent(in), target :: grid
    character(*), intent(in) :: name
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg
    type(msp_dict) :: err
    integer, dimension(:), allocatable :: i, j
    complex(dp), dimension(:), allocatable :: v
    character(len=1024) :: line
    character(len=32), dimension(5) :: word
    integer, dimension(3) :: shape
    integer :: unit, io, k, n, e
    real(dp) :: re, im
    logical :: cplx, mirror, hermitian

    call msp_clear(stat)
    allocate(i(0), j(0), v(0))
    n = 0
    cplx = .false.
    if (grid%rank == 0) then
       open (newunit=unit, file=name, action='read', status='old', &
            & iostat=io)
       if (io /= 0) then
          call err%set("message", "msp_io: cannot open " // name)
       else
          word = ''
          read (unit, '(a)', iostat=io) line
          if (io == 0) read (line, *, iostat=io) word
          word = lowercase(word)
          if (io /= 0 .or. word(1) /= "%%matrixmarket" .or. &
               & word(2) /= "matrix" .or. word(3) /= "coordinate" .or. &
               & all(word(4) /= ["real   ", "complex"]) .or. &
               & all(word(5) /= ["general  ", "symmetric", "hermitian"])) &
               & then
             call err%set("message", "msp_io: " // name // " is not a " // &
                  & "real or complex Matrix Market coordinate matrix")
          else
             cplx = word(4) == "complex"
             mirror = word(5) /= "general"
             hermitian = word(5) == "hermitian"
             do
                read (unit, '(a)', iostat=io) line
                if (io /= 0) exit
                if (len_trim(line) == 0) cycle
                if (line(1:1) /= "%") exit
             end do
             shape = 0
             if (io == 0) read (line, *, iostat=io) shape
             if (io /= 0) then
                call err%set("message", "msp_io: cannot read shape of " // &
                     & name)
             else if (shape(1) /= shape(2)) then
                call err%set("message", "msp_io: " // name // &
                     & " is not square")
             else if (shape(1) < 1 .or. shape(3) < 0) then
                call err%set("message", "msp_io: invalid shape of " // name)
             else
                n = shape(1)
                deallocate(i, j, v)
                allocate(i(2 * shape(3)), j(2 * shape(3)), v(2 * shape(3)))
                k = 0
                do e = 1, shape(3)
                   im = 0
                   k = k + 1
                   if (cplx) then
                      read (unit, *, iostat=io) i(k), j(k), re, im
                   else
                      read (unit, *, iostat=io) i(k), j(k), re
                   end if
                   if (io /= 0) exit
                   v(k) = cmplx(re, im, dp)
                   if (mirror .and. i(k) /= j(k)) then
                      k = k + 1
                      i(k) = j(k - 1)
                      j(k) = i(k - 1)
                      v(k) = v(k - 1)
                      if (hermitian) v(k) = conjg(v(k))
                   end if
                end do
                if (io /= 0) then
                   call err%set("message", "msp_io: cannot read entries of " &
                        & // name)
                else
                   i = i(1:k)
                   j = j(1:k)
                   v = v(1:k)
                end if
             end if
          end if
          close (unit)
       end if
    end if
    call grid%agree(err, stat, errmsg)
    if (err%num_keys() > 0) return
    call broadcast(grid, n, cplx)
    call A%init(grid, n)
    if (cplx) then
       call A%fill(i, j, v, stat=stat, errmsg=errmsg)
    else
       call A%fill(i, j, real(v, dp), stat=stat, errmsg=errmsg)
    end if
  end subroutine msp_read_matrix_market

  !> Write A as a Matrix Market coordinate file (real or complex, general).
  !> Ranks append their entries in turn. Collective.
  subroutine msp_write_matrix_market(A, name)
    type(msp_matrix), intent(in) :: A
    character(*), intent(in) :: name
    integer, dimension(:), allocatable :: i, j
    complex(dp), dimension(:), allocatable :: v
    integer(int64) :: total
    integer :: unit, r, k

    call A%entries(i, j, v)
    total = A%nnz()
    if (A%grid%rank == 0) then
       open (newunit=unit, file=name, action='write', status='replace')
       write (unit, '(a)') "%%MatrixMarket matrix coordinate " // &
            & trim(merge("complex", "real   ", A%is_complex())) // " general"
       write (unit, '(i0, 1x, i0, 1x, i0)') A%n, A%n, total
       close (unit)
    end if
    do r = 0, A%grid%size - 1
       call MPI_Barrier(A%grid%comm)
       if (r /= A%grid%rank) cycle
       open (newunit=unit, file=name, action='write', status='old', &
            & position='append')
       do k = 1, size(i)
          if (A%is_complex()) then
             write (unit, '(i0, 1x, i0, 2(1x, es24.16e3))') i(k), j(k), v(k)
          else
             write (unit, '(i0, 1x, i0, 1x, es24.16e3)') i(k), j(k), v(k)%re
          end if
       end do
       close (unit)
    end do
    call MPI_Barrier(A%grid%comm)
  end subroutine msp_write_matrix_market

  !> A becomes the matrix of the NTPoly binary file name: a header of rows,
  !> columns and a complex flag (4-byte integers) and the entry count
  !> (8 bytes), then triplets of row and column (4 bytes) and value (8 bytes
  !> real or 16 bytes complex), 1-based. Each rank reads an equal share of
  !> the triplets. Collective; errors are agreed on all ranks.
  subroutine msp_read_binary(A, grid, name, stat, errmsg)
    type(msp_matrix), intent(inout) :: A
    type(msp_grid), intent(in), target :: grid
    character(*), intent(in) :: name
    integer, intent(out), optional :: stat
    character(:), allocatable, intent(out), optional :: errmsg
    type(msp_dict) :: err
    integer(int32) :: rows, cols, flag
    integer(int64) :: total, first, count, q, width
    integer, dimension(:), allocatable :: i, j
    real(dp), dimension(:), allocatable :: vr
    complex(dp), dimension(:), allocatable :: vc
    integer :: unit, io
    logical :: opened

    call msp_clear(stat)
    open (newunit=unit, file=name, access='stream', form='unformatted', &
         & action='read', status='old', iostat=io)
    opened = io == 0
    if (opened) read (unit, pos=1, iostat=io) rows, cols, flag, total
    if (io /= 0) then
       call err%set("message", "msp_io: cannot read " // name)
    else if (rows /= cols) then
       call err%set("message", "msp_io: " // name // " is not square")
    else if (rows < 1 .or. total < 0) then
       call err%set("message", "msp_io: invalid shape of " // name)
    end if
    call grid%agree(err, stat, errmsg)
    if (err%num_keys() > 0) then
       if (opened) close (unit)
       return
    end if
    width = merge(24, 16, flag /= 0)
    first = grid%rank * (total / grid%size)
    count = total / grid%size
    if (grid%rank == grid%size - 1) count = total - first
    allocate(i(count), j(count))
    if (flag /= 0) then
       allocate(vc(count))
       if (count > 0) read (unit, pos=21 + width * first, iostat=io) &
            & (i(q), j(q), vc(q), q = 1, count)
    else
       allocate(vr(count))
       if (count > 0) read (unit, pos=21 + width * first, iostat=io) &
            & (i(q), j(q), vr(q), q = 1, count)
    end if
    close (unit)
    if (io /= 0) call err%set("message", &
         & "msp_io: cannot read entries of " // name)
    call grid%agree(err, stat, errmsg)
    if (err%num_keys() > 0) return
    call A%init(grid, int(rows))
    if (flag /= 0) then
       call A%fill(i, j, vc, stat=stat, errmsg=errmsg)
    else
       call A%fill(i, j, vr, stat=stat, errmsg=errmsg)
    end if
  end subroutine msp_read_binary

  !> Write A in NTPoly's binary format (see msp_read_binary). Ranks append
  !> their entries in turn. Collective.
  subroutine msp_write_binary(A, name)
    type(msp_matrix), intent(in) :: A
    character(*), intent(in) :: name
    integer, dimension(:), allocatable :: i, j
    complex(dp), dimension(:), allocatable :: v
    integer(int64) :: q, total
    integer :: unit, r

    call A%entries(i, j, v)
    total = A%nnz()
    if (A%grid%rank == 0) then
       open (newunit=unit, file=name, access='stream', form='unformatted', &
            & action='write', status='replace')
       write (unit) int(A%n, int32), int(A%n, int32), &
            & merge(1_int32, 0_int32, A%is_complex()), total
       close (unit)
    end if
    do r = 0, A%grid%size - 1
       call MPI_Barrier(A%grid%comm)
       if (r /= A%grid%rank) cycle
       open (newunit=unit, file=name, access='stream', form='unformatted', &
            & action='write', status='old', position='append')
       if (A%is_complex()) then
          write (unit) (int(i(q), int32), int(j(q), int32), v(q), &
               & q = 1, size(i, kind=int64))
       else
          write (unit) (int(i(q), int32), int(j(q), int32), real(v(q), dp), &
               & q = 1, size(i, kind=int64))
       end if
       close (unit)
    end do
    call MPI_Barrier(A%grid%comm)
  end subroutine msp_write_binary

  !> n and cplx of the first rank on every rank.
  subroutine broadcast(grid, n, cplx)
    type(msp_grid), intent(in) :: grid
    integer, intent(inout) :: n
    logical, intent(inout) :: cplx
    type(msp_dict) :: d

    call d%set("n", n)
    call d%set("complex", cplx)
    call grid%bcast(d, 0)
    n = d%get_int("n")
    cplx = d%get_logical("complex")
  end subroutine broadcast

  !> The words in lower case.
  elemental function lowercase(word) result(lower)
    character(*), intent(in) :: word
    character(len=len(word)) :: lower
    integer :: k, c

    lower = word
    do k = 1, len(word)
       c = iachar(word(k:k))
       if (c >= iachar("A") .and. c <= iachar("Z")) &
            & lower(k:k) = achar(c + 32)
    end do
  end function lowercase
end module msp_io_mod
