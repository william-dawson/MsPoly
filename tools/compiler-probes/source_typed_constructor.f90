program p5
  implicit none
  class(*), allocatable :: a(:)
  allocate(a, source=[character(len=5) :: "alpha", "beta"])
  print *, "ok"
end program
