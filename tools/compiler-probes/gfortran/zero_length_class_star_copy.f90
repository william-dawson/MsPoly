! Copy a class(*) component holding a zero-length character value with
! ALLOCATE(SOURCE=) from the associate name of SELECT TYPE ... CLASS DEFAULT.
module m
  implicit none
  type :: entry_t
     class(*), allocatable :: v
  end type
contains
  subroutine copy_entry(dst, src)
    type(entry_t), intent(inout) :: dst
    type(entry_t), intent(in) :: src
    if (allocated(dst%v)) deallocate(dst%v)
    select type (x => src%v)
    class default
       allocate(dst%v, source=x)
    end select
  end subroutine
end module
program es3
  use m
  implicit none
  type(entry_t) :: a, b
  allocate(a%v, source="ab")
  call copy_entry(b, a)
  select type (y => b%v); type is (character(*)); print *, "len 2 copy: [", y, "]"; end select
  deallocate(a%v)
  allocate(character(len=0) :: a%v)
  call copy_entry(b, a)
  select type (y => b%v); type is (character(*)); print *, "len 0 copy: len", len(y); end select
  print '(A)', "RESULT: all passed"
end program
