! A logical mask allocated on assignment from an elemental expression, then
! used in pack. F2008 7.2.1.3: an unallocated allocatable variable on the
! left of an intrinsic assignment is allocated with the shape of the
! expression, so keep has the size of i and pack(i, keep) the count of
! true elements. Both an empty and a filled piece are checked.
module entries_mod
  implicit none
  private
  public :: entries, take

contains

  !> Triplets (k, k, k) for k = 1, ..., n.
  subroutine entries(n, i, j, v)
    integer, intent(in) :: n
    integer, dimension(:), allocatable, intent(out) :: i, j
    complex(8), dimension(:), allocatable, intent(out) :: v
    integer :: k

    allocate(i(n), j(n), v(n))
    i = [(k, k = 1, n)]
    j = i
    v = cmplx(i, 0, 8)
  end subroutine entries

  !> The number of triplets handed over.
  integer function take(i, j, v)
    integer, dimension(:), intent(in) :: i, j
    real(8), dimension(:), intent(in) :: v

    take = size(i)
    if (size(j) /= take .or. size(v) /= take) take = -1
  end function take
end module entries_mod

program r14_mask_reallocation_pack
  use entries_mod, only: entries, take
  implicit none
  integer :: failures

  failures = 0
  call leading(0, 6, 0)
  call leading(36, 6, 6)
  if (failures == 0) then
     print '(a)', "RESULT: all passed"
  else
     print '(a)', "RESULT: pack with a reallocated mask gave the wrong count"
     error stop 1
  end if

contains

  !> The entries of an n-entry piece in the leading m x m block.
  subroutine leading(n, m, expected)
    integer, intent(in) :: n, m, expected
    integer, dimension(:), allocatable :: i, j
    complex(8), dimension(:), allocatable :: v
    logical, dimension(:), allocatable :: keep

    call entries(n, i, j, v)
    keep = i <= m .and. j <= m
    if (take(pack(i, keep), pack(j, keep), real(pack(v, keep), 8)) &
         & /= expected) failures = failures + 1
  end subroutine leading
end program r14_mask_reallocation_pack
