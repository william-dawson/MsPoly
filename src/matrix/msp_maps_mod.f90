!> Entry-wise transformations of matrices: symmetric permutations, which
!> spread the entries of a matrix evenly over the ranks, and maps of a
!> function over the entries.
module msp_maps_mod
  use, intrinsic :: iso_fortran_env, only: int64
  use msp_kinds_mod, only: dp
  use msp_matrix_mod, only: msp_matrix, msp_to_single, msp_to_device
  implicit none
  private

  public :: msp_random_permutation, msp_permute, msp_map, msp_entry_map

  abstract interface
     !> Called for each entry (i, j, v) of a matrix, with global indices: may
     !> change v, and returns whether to keep the entry.
     logical function msp_entry_map(i, j, v)
       import :: dp
       integer, intent(in) :: i, j
       complex(dp), intent(inout) :: v
     end function msp_entry_map
  end interface

contains

  !> A random permutation of 1, ..., n from the seed, the same on every rank
  !> (Fisher-Yates shuffle driven by a 64-bit linear congruential
  !> generator).
  function msp_random_permutation(n, seed) result(perm)
    integer, intent(in) :: n, seed
    integer, dimension(n) :: perm
    integer(int64) :: state
    integer :: k, pick, swap

    perm = [(k, k = 1, n)]
    state = seed
    do k = n, 2, -1
       state = state * 6364136223846793005_int64 + 1442695040888963407_int64
       pick = 1 + int(modulo(ishft(state, -33), int(k, int64)))
       swap = perm(k)
       perm(k) = perm(pick)
       perm(pick) = swap
    end do
  end function msp_random_permutation

  !> B becomes P A P^T: entry (i, j) of A moves to (perm(i), perm(j)), or,
  !> with inverse, from (perm(i), perm(j)) back to (i, j). Collective.
  subroutine msp_permute(A, perm, B, inverse)
    type(msp_matrix), intent(in) :: A
    integer, dimension(:), intent(in) :: perm
    type(msp_matrix), intent(inout) :: B
    logical, intent(in), optional :: inverse
    integer, dimension(size(perm)) :: to
    integer :: k

    to = perm
    if (present(inverse)) then
       if (inverse) to(perm) = [(k, k = 1, size(perm))]
    end if
    call transform(A, B, to=to)
  end subroutine msp_permute

  !> B becomes A with f applied to every entry; entries for which f returns
  !> false are dropped. B is real when A is (the real part of each new
  !> value is kept). Collective.
  subroutine msp_map(A, B, f)
    type(msp_matrix), intent(in) :: A
    type(msp_matrix), intent(inout) :: B
    procedure(msp_entry_map) :: f

    call transform(A, B, f=f)
  end subroutine msp_map

  !> B from the entries of A, each passed through f and moved to
  !> (to(i), to(j)) when given, with the kind, precision and location of
  !> A.
  subroutine transform(A, B, f, to)
    type(msp_matrix), intent(in) :: A
    type(msp_matrix), intent(inout) :: B
    procedure(msp_entry_map), optional :: f
    integer, dimension(:), intent(in), optional :: to
    integer, dimension(:), allocatable :: i, j
    complex(dp), dimension(:), allocatable :: v
    logical, dimension(:), allocatable :: keep
    integer :: k

    call A%entries(i, j, v)
    allocate(keep(size(i)))
    keep = .true.
    if (present(f)) then
       do k = 1, size(i)
          keep(k) = f(i(k), j(k), v(k))
       end do
    end if
    if (present(to)) then
       i = to(i)
       j = to(j)
    end if
    call B%init(A%grid, A%n)
    if (A%is_complex()) then
       call B%fill(pack(i, keep), pack(j, keep), pack(v, keep))
    else
       call B%fill(pack(i, keep), pack(j, keep), real(pack(v, keep), dp))
    end if
    if (A%is_single()) call msp_to_single(B)
    if (A%on_device()) call msp_to_device(B)
  end subroutine transform
end module msp_maps_mod
