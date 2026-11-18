program probe3
  implicit none
  class(*), allocatable :: a(:), b(:), c(:)
  character(len=5) :: names(2) = ["alpha", "beta "]
  allocate(a, source=names)
  call move_alloc(a, b)
  select type (s => b)
  type is (character(*))
     print *, "move_alloc: len=", len(s), " s(1)=[", s(1), "]"
  end select
  allocate(a, source=names)
  c = a     ! intrinsic assignment instead
  select type (s => c)
  type is (character(*))
     print *, "assign:     len=", len(s), " s(1)=[", s(1), "]"
  end select
end program
