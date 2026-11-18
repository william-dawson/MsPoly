! The polymorphic patterns of the MsPoly matrix layer: an abstract storage
! type with deferred bindings, extensions with allocatable arrays, a
! class(base), allocatable component, dispatch on the other operand,
! explicit copy into an intent(out) polymorphic dummy, move_alloc between
! polymorphic allocatables, and ALLOCATE(MOLD=) from a component.
module storage_mod
  implicit none
  private
  public :: storage_t, real_t, complex_t, holder_t, copy_holder

  type, abstract :: storage_t
     integer :: n = 0
   contains
     procedure(copy_i), deferred :: copy_to
     procedure(dot_i), deferred :: dot
     procedure(scale_i), deferred :: scale
  end type storage_t

  abstract interface
     subroutine copy_i(this, dst)
       import :: storage_t
       class(storage_t), intent(in) :: this
       class(storage_t), allocatable, intent(out) :: dst
     end subroutine copy_i
     complex(8) function dot_i(this, other)
       import :: storage_t
       class(storage_t), intent(in) :: this, other
     end function dot_i
     subroutine scale_i(this, alpha)
       import :: storage_t
       class(storage_t), intent(inout) :: this
       complex(8), intent(in) :: alpha
     end subroutine scale_i
  end interface

  type, extends(storage_t) :: real_t
     real(8), allocatable :: v(:)
   contains
     procedure :: copy_to => copy_real
     procedure :: dot => dot_real
     procedure :: scale => scale_real
  end type real_t

  type, extends(storage_t) :: complex_t
     complex(8), allocatable :: v(:)
   contains
     procedure :: copy_to => copy_complex
     procedure :: dot => dot_complex
     procedure :: scale => scale_complex
  end type complex_t

  type :: holder_t
     integer :: n = 0
     class(storage_t), allocatable :: local
  end type holder_t

contains

  subroutine copy_real(this, dst)
    class(real_t), intent(in) :: this
    class(storage_t), allocatable, intent(out) :: dst
    allocate(real_t :: dst)
    select type (dst)
    type is (real_t)
       dst%n = this%n
       dst%v = this%v
    end select
  end subroutine copy_real

  subroutine copy_complex(this, dst)
    class(complex_t), intent(in) :: this
    class(storage_t), allocatable, intent(out) :: dst
    allocate(complex_t :: dst)
    select type (dst)
    type is (complex_t)
       dst%n = this%n
       dst%v = this%v
    end select
  end subroutine copy_complex

  complex(8) function dot_real(this, other) result(d)
    class(real_t), intent(in) :: this
    class(storage_t), intent(in) :: other
    d = 0
    select type (other)
    type is (real_t)
       d = sum(this%v * other%v)
    end select
  end function dot_real

  complex(8) function dot_complex(this, other) result(d)
    class(complex_t), intent(in) :: this
    class(storage_t), intent(in) :: other
    d = 0
    select type (other)
    type is (complex_t)
       d = sum(conjg(this%v) * other%v)
    end select
  end function dot_complex

  subroutine scale_real(this, alpha)
    class(real_t), intent(inout) :: this
    complex(8), intent(in) :: alpha
    this%v = real(alpha, 8) * this%v
  end subroutine scale_real

  subroutine scale_complex(this, alpha)
    class(complex_t), intent(inout) :: this
    complex(8), intent(in) :: alpha
    this%v = alpha * this%v
  end subroutine scale_complex

  !> Explicit holder copy (no intrinsic assignment of holder_t).
  subroutine copy_holder(dst, src)
    type(holder_t), intent(inout) :: dst
    type(holder_t), intent(in) :: src
    dst%n = src%n
    if (allocated(dst%local)) deallocate(dst%local)
    if (allocated(src%local)) call src%local%copy_to(dst%local)
  end subroutine copy_holder
end module storage_mod

program storage_patterns
  use storage_mod
  implicit none
  type(holder_t) :: a, b
  class(storage_t), allocatable :: promoted, shaped
  integer :: failures = 0
  complex(8) :: d

  ! class(base), allocatable component, typed allocation
  a%n = 3
  allocate(real_t :: a%local)
  select type (s => a%local)
  type is (real_t)
     s%n = 3
     s%v = [1d0, 2d0, 3d0]
  end select

  ! explicit copy into intent(out) polymorphic dummy; independence
  call copy_holder(b, a)
  call b%local%scale((2d0, 0d0))
  d = a%local%dot(a%local)
  call expect(abs(d - 14d0) < 1d-12, "copy is independent of the original")
  d = a%local%dot(b%local)
  call expect(abs(d - 28d0) < 1d-12, "dispatch on the other operand")

  ! promotion real -> complex, swapped in with move_alloc
  select type (s => a%local)
  type is (real_t)
     allocate(complex_t :: promoted)
     select type (p => promoted)
     type is (complex_t)
        p%n = s%n
        p%v = cmplx(s%v, 0d0, 8)
     end select
  end select
  call move_alloc(promoted, a%local)
  call a%local%scale((0d0, 1d0))
  d = a%local%dot(a%local)
  call expect(abs(d - 14d0) < 1d-12 .and. .not. allocated(promoted), &
       "move_alloc of a class(base) allocatable")
  select type (s => a%local)
  type is (complex_t)
     call expect(abs(s%v(2) - (0d0, 2d0)) < 1d-12, "complex scale after swap")
  class default
     call expect(.false., "dynamic type after move_alloc")
  end select

  ! ALLOCATE(MOLD=) from a polymorphic component
  allocate(shaped, mold=b%local)
  select type (shaped)
  type is (real_t)
     call expect(.not. allocated(shaped%v), "mold gives the dynamic type")
  class default
     call expect(.false., "mold dynamic type")
  end select

  ! intent(out) deallocation: copying over an allocated destination
  call copy_holder(b, a)
  select type (s => b%local)
  type is (complex_t)
     call expect(size(s%v) == 3, "copy over an existing destination")
  class default
     call expect(.false., "copy over destination type")
  end select

  if (failures == 0) then
     write (*, '(A)') "RESULT: all passed"
  else
     write (*, '(A,I0,A)') "RESULT: ", failures, " failure(s)"
     error stop 1
  end if
contains
  subroutine expect(ok, label)
    logical, intent(in) :: ok
    character(*), intent(in) :: label
    if (.not. ok) then
       failures = failures + 1
       write (*, '(A)') "FAIL: " // label
    end if
  end subroutine expect
end program storage_patterns
