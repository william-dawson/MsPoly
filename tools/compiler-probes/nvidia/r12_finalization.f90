! Finalization of an extension of an abstract type (whose components have
! default initialization) through class(base) allocatables: deallocate,
! intent(out) dummy, end of scope, move_alloc onto an allocated target, a
! polymorphic component of a holder going out of scope, and reallocation.
! Each finalization releases one "resource"; every case must release
! exactly once.
module final_mod
  implicit none
  private
  public :: base_t, owner_t, holder_t, live, make, take

  integer :: live = 0

  type, abstract :: base_t
     integer :: rows = 0
  end type base_t

  type, extends(base_t) :: owner_t
     logical :: owns = .false.
   contains
     final :: release
  end type owner_t

  type :: holder_t
     class(base_t), allocatable :: local
  end type holder_t
contains
  subroutine release(this)
    type(owner_t), intent(inout) :: this
    if (this%owns) live = live - 1
    this%owns = .false.
  end subroutine release

  subroutine make(x)
    class(base_t), allocatable, intent(out) :: x
    allocate(owner_t :: x)
    select type (x)
    type is (owner_t)
       x%owns = .true.
    end select
    live = live + 1
  end subroutine make

  subroutine take(x)
    class(base_t), allocatable, intent(out) :: x
    if (allocated(x)) live = -100
  end subroutine take
end module final_mod

program finalization_patterns
  use final_mod
  implicit none
  class(base_t), allocatable :: a, b
  type(holder_t), allocatable :: h
  integer :: failures = 0

  call make(a)
  deallocate(a)
  call expect(live == 0, "deallocate finalizes")

  call make(a)
  call take(a)
  call expect(live == 0, "intent(out) finalizes")

  call scoped()
  call expect(live == 0, "a local allocatable is finalized at end of scope")

  call make(a)
  call make(b)
  call move_alloc(a, b)
  call expect(live == 1 .and. .not. allocated(a), "move_alloc finalizes the old target only")
  deallocate(b)
  call expect(live == 0, "moved object finalized once")

  call make(a)
  call make(a)
  call expect(live == 1, "intent(out) in make finalizes the previous object")
  deallocate(a)

  allocate(h)
  call make(h%local)
  deallocate(h)
  call expect(live == 0, "polymorphic component finalized with its holder")

  call holder_scope()
  call expect(live == 0, "component of a local holder finalized at end of scope")

  if (failures == 0) then
     write (*, '(A)') "RESULT: all passed"
  else
     write (*, '(A,I0,A)') "RESULT: ", failures, " failure(s)"
     error stop 1
  end if
contains
  subroutine scoped()
    class(base_t), allocatable :: x
    call make(x)
  end subroutine scoped
  subroutine holder_scope()
    type(holder_t) :: g
    call make(g%local)
  end subroutine holder_scope
  subroutine expect(ok, label)
    logical, intent(in) :: ok
    character(*), intent(in) :: label
    if (.not. ok) then
       failures = failures + 1
       write (*, '(A,I0,A)') "FAIL (live=", live, "): " // label
    end if
  end subroutine expect
end program finalization_patterns
