! Inside SELECT TYPE on a class(*) array component, pass an element of the
! (now typed) associate name to a class(*) dummy.
module r5_mod
  implicit none
  type :: entry_t
     class(*), allocatable :: a(:)
  end type entry_t
contains
  function describe(x) result(str)
    class(*), intent(in) :: x
    character(:), allocatable :: str
    character(len=30) :: buffer
    select type (x)
    type is (real(8))
       write (buffer, '(F6.2)') x
       str = trim(adjustl(buffer))
    class default
       str = "null"
    end select
  end function describe
  function list(en) result(str)
    type(entry_t), intent(in) :: en
    character(:), allocatable :: str
    integer :: k
    str = ""
    select type (x => en%a)
    type is (real(8))
       do k = 1, size(x)
          str = str // " " // describe(x(k))
       end do
    end select
  end function list
end module r5_mod

program r5
  use r5_mod
  implicit none
  type(entry_t) :: e
  character(:), allocatable :: got
  allocate(real(8) :: e%a(2))
  select type (x => e%a)
  type is (real(8))
     x = [1.5d0, 2.5d0]
  end select
  got = list(e)
  write (*, '(A)') "got:" // got
  if (got == " 1.50 2.50") then
     write (*, '(A)') "RESULT: all passed"
  else
     write (*, '(A)') "RESULT: FAIL (expected ' 1.50 2.50')"
     error stop 1
  end if
end program r5
