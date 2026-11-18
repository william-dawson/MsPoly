program p6
  implicit none
  class(*), allocatable :: a(:), b(:)
  real(8), allocatable :: z(:)
  allocate(z(0))
  allocate(a, source=z)
  select type (x => a); type is (real(8)); print *, "source= keeps type, size", size(x)
  class default; print *, "source= LOSES type"; end select
  allocate(real(8) :: b(size(z)))
  select type (x => b); type is (real(8)); x = z; print *, "typed allocate keeps type, size", size(x)
  class default; print *, "typed allocate LOSES type"; end select
end program
