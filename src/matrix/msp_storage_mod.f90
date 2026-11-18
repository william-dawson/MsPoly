!> Interface of the local piece of a distributed matrix.
!>
!> Each rank owns one rows x cols piece with local indices. The distributed
!> algorithms are written against this interface only; implementations hold
!> real or complex values in double or single precision, on the host or on a
!> device. Arithmetic is done in double precision whatever the storage:
!> single precision values are widened, combined, and rounded once when
!> stored. Binary operations require both operands to have the same
!> implementation: the matrix layer converts operands first. A piece may own
!> resources that Fortran does not manage (device memory): drop pieces with
!> msp_free.
module msp_storage_mod
  use mpi_f08, only: MPI_Comm
  use msp_kinds_mod, only: dp, i64
  implicit none
  private

  public :: msp_storage, msp_free

  type, abstract :: msp_storage
     integer :: rows = 0
     integer :: cols = 0
   contains
     procedure(release_i), deferred :: release
     procedure(is_complex_i), deferred :: is_complex
     procedure(is_complex_i), deferred :: is_single
     procedure(nnz_i), deferred :: nnz
     procedure(copy_i), deferred :: copy_to
     procedure(convert_i), deferred :: convert_to
     procedure(scale_i), deferred :: scale
     procedure(add_i), deferred :: add
     procedure(dot_i), deferred :: dot
     procedure(trace_i), deferred :: trace
     procedure(norm2_i), deferred :: norm2
     procedure(gather_i), deferred :: gather_rows
     procedure(gather_i), deferred :: gather_cols
     procedure(multiply_i), deferred :: multiply
  end type msp_storage

  abstract interface
     !> Free resources that deallocation does not.
     subroutine release_i(this)
       import :: msp_storage
       class(msp_storage), intent(inout) :: this
     end subroutine release_i

     !> Whether values are complex (is_complex) or single precision
     !> (is_single).
     logical function is_complex_i(this)
       import :: msp_storage
       class(msp_storage), intent(in) :: this
     end function is_complex_i

     !> Number of stored entries.
     integer(i64) function nnz_i(this)
       import :: msp_storage, i64
       class(msp_storage), intent(in) :: this
     end function nnz_i

     !> dst becomes a deep copy; a previous dst is freed.
     subroutine copy_i(this, dst)
       import :: msp_storage
       class(msp_storage), intent(in) :: this
       class(msp_storage), allocatable, intent(inout) :: dst
     end subroutine copy_i

     !> dst becomes a copy in the same location with complex (cplx) or real,
     !> single or double precision values; a previous dst is freed. Complex
     !> values cannot become real.
     subroutine convert_i(this, dst, cplx, single)
       import :: msp_storage
       class(msp_storage), intent(in) :: this
       class(msp_storage), allocatable, intent(inout) :: dst
       logical, intent(in) :: cplx, single
     end subroutine convert_i

     !> this = alpha * this. Real storage uses the real part of alpha.
     subroutine scale_i(this, alpha)
       import :: msp_storage, dp
       class(msp_storage), intent(inout) :: this
       complex(dp), intent(in) :: alpha
     end subroutine scale_i

     !> this = alpha * x + this, dropping exact zeros and entries with
     !> |value| < threshold.
     subroutine add_i(this, x, alpha, threshold)
       import :: msp_storage, dp
       class(msp_storage), intent(inout) :: this
       class(msp_storage), intent(in) :: x
       complex(dp), intent(in) :: alpha
       real(dp), intent(in) :: threshold
     end subroutine add_i

     !> Sum of conjg(this) * x over all positions.
     complex(dp) function dot_i(this, x)
       import :: msp_storage, dp
       class(msp_storage), intent(in) :: this, x
     end function dot_i

     !> Sum of the entries on the global diagonal, for a piece whose first
     !> row and column have global offsets row0 and col0.
     complex(dp) function trace_i(this, row0, col0)
       import :: msp_storage, dp
       class(msp_storage), intent(in) :: this
       integer, intent(in) :: row0, col0
     end function trace_i

     !> Sum of |value|**2.
     real(dp) function norm2_i(this)
       import :: msp_storage, dp
       class(msp_storage), intent(in) :: this
     end function norm2_i

     !> Concatenate the pieces of every rank of comm, in rank order:
     !> side by side (gather_rows, a row band) or stacked (gather_cols, a
     !> column band); a previous band is freed. Collective over comm.
     subroutine gather_i(this, comm, band)
       import :: msp_storage, MPI_Comm
       class(msp_storage), intent(in) :: this
       type(MPI_Comm), intent(in) :: comm
       class(msp_storage), allocatable, intent(inout) :: band
     end subroutine gather_i

     !> this = alpha * a * b + beta * this, dropping exact zeros and entries
     !> with |value| < threshold. a is a row band, b a column band.
     subroutine multiply_i(this, a, b, alpha, beta, threshold)
       import :: msp_storage, dp
       class(msp_storage), intent(inout) :: this
       class(msp_storage), intent(in) :: a, b
       complex(dp), intent(in) :: alpha, beta
       real(dp), intent(in) :: threshold
     end subroutine multiply_i
  end interface

contains

  !> Release and deallocate x, if allocated.
  subroutine msp_free(x)
    class(msp_storage), allocatable, intent(inout) :: x

    if (.not. allocated(x)) return
    call x%release()
    deallocate(x)
  end subroutine msp_free
end module msp_storage_mod
