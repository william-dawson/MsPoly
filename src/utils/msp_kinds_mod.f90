!> Kind parameters used throughout MsPoly.
module msp_kinds_mod
  use, intrinsic :: iso_fortran_env, only: int32, int64, real32, real64
  implicit none
  private

  !> Floating point kind for storage and arithmetic.
  integer, parameter, public :: dp = real64
  !> Floating point kind of values sent in single precision.
  integer, parameter, public :: sp = real32
  !> Integer kind for local indices.
  integer, parameter, public :: i32 = int32
  !> Integer kind for sizes, counts and offsets.
  integer, parameter, public :: i64 = int64
end module msp_kinds_mod
