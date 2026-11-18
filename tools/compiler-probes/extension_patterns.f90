! A concrete extension of an abstract type, extended again with some
! bindings overridden; SELECT TYPE with CLASS IS on the middle type; and an
! overridden copy into an intent(out) polymorphic dummy.
module ext_mod
  implicit none
  private
  public :: base_t, host_t, device_t

  type, abstract :: base_t
   contains
     procedure(name_i), deferred :: name
     procedure(copy_i), deferred :: copy_to
  end type base_t
  abstract interface
     function name_i(this) result(s)
       import :: base_t
       class(base_t), intent(in) :: this
       character(:), allocatable :: s
     end function name_i
     subroutine copy_i(this, dst)
       import :: base_t
       class(base_t), intent(in) :: this
       class(base_t), allocatable, intent(out) :: dst
     end subroutine copy_i
  end interface

  type, extends(base_t) :: host_t
     real(8), allocatable :: v(:)
   contains
     procedure :: name => host_name
     procedure :: copy_to => host_copy
     procedure :: total
  end type host_t

  type, extends(host_t) :: device_t
     integer :: device = -1
   contains
     procedure :: name => device_name
     procedure :: copy_to => device_copy
  end type device_t
contains
  function host_name(this) result(s)
    class(host_t), intent(in) :: this
    character(:), allocatable :: s
    s = "host"
    if (.not. allocated(this%v)) s = "host?"
  end function host_name
  function device_name(this) result(s)
    class(device_t), intent(in) :: this
    character(:), allocatable :: s
    s = "device"
    if (this%device < -1) s = "device?"
  end function device_name
  real(8) function total(this)
    class(host_t), intent(in) :: this
    total = sum(this%v)
  end function total
  subroutine host_copy(this, dst)
    class(host_t), intent(in) :: this
    class(base_t), allocatable, intent(out) :: dst
    allocate(host_t :: dst)
    select type (dst)
    type is (host_t)
       dst%v = this%v
    end select
  end subroutine host_copy
  subroutine device_copy(this, dst)
    class(device_t), intent(in) :: this
    class(base_t), allocatable, intent(out) :: dst
    allocate(device_t :: dst)
    select type (dst)
    type is (device_t)
       dst%v = this%v
       dst%device = this%device
    end select
  end subroutine device_copy
end module ext_mod

program extension_patterns
  use ext_mod
  implicit none
  class(base_t), allocatable :: a, b
  integer :: failures = 0

  allocate(device_t :: a)
  select type (a)
  type is (device_t)
     a%v = [1d0, 2d0]
     a%device = 0
  end select
  call expect(a%name() == "device", "override dispatches")
  select type (a)
  class is (host_t)
     call expect(a%total() == 3d0, "class is (middle type) matches the extension")
  class default
     call expect(.false., "class is (middle type)")
  end select
  call a%copy_to(b)
  call expect(b%name() == "device", "overridden copy keeps the dynamic type")
  select type (b)
  type is (device_t)
     call expect(b%device == 0 .and. size(b%v) == 2, "copied components")
  class default
     call expect(.false., "copy type")
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
end program extension_patterns
