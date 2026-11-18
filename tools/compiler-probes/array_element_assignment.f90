module bis2_mod
  implicit none
  type :: entry_t
     character(:), allocatable :: key
     class(*), allocatable :: v
  end type
contains
  subroutine show(label, x)
    character(*), intent(in) :: label
    class(*), intent(in) :: x
    select type (x)
    type is (real(8)); print '(A,F6.2)', label, x
    class default; print '(A)', label // " wrong type"
    end select
  end subroutine
  subroutine copy_in(dst, src)          ! assignment to dummy array element
    type(entry_t), intent(inout) :: dst(:)
    type(entry_t), intent(in) :: src
    dst(1) = src
  end subroutine
  subroutine copy_out(src, dst)         ! intent(out) dummy
    type(entry_t), intent(in) :: src
    type(entry_t), intent(out) :: dst
    dst = src
  end subroutine
end module
program bis2
  use bis2_mod
  implicit none
  type(entry_t) :: a, b, c
  type(entry_t), allocatable :: arr(:)
  allocate(a%v, source=1.5d0)
  b = a;                         call show("1 scalar = scalar:         ", b%v)
  allocate(arr(2)); arr(1) = a;  call show("2 array elem = scalar:     ", arr(1)%v)
  call copy_in(arr, a);          call show("3 dummy array elem = dummy:", arr(1)%v)
  call copy_out(a, c);           call show("4 intent(out) = dummy:     ", c%v)
  deallocate(arr); allocate(arr(2))
  allocate(arr(1)%v, source=a%v); call show("5 allocate source=:        ", arr(1)%v)
end program
