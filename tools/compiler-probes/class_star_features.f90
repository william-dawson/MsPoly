! Features the class(*) dictionary relies on.
module probe_mod
  implicit none
  type :: entry_t
     character(:), allocatable :: key
     class(*), allocatable :: v        ! scalar base type or nested dict
     class(*), allocatable :: a(:)     ! rank-1 array of a base type
  end type
  type :: dict_t
     type(entry_t), allocatable :: e(:)
     integer :: n = 0
  end type
end module
program probe2
  use probe_mod
  implicit none
  type(dict_t) :: d, inner, copy
  type(entry_t), allocatable :: grown(:)
  integer :: fails = 0
  character(len=5) :: names(2) = ["alpha", "beta "]
  real(8) :: hist(3) = [1d0, 2d0, 3d0]
  allocate(d%e(4))
  d%e(1)%key = "threshold"; allocate(d%e(1)%v, source=1.0d-6)
  d%e(2)%key = "method";    allocate(d%e(2)%v, source="trs2")
  d%e(3)%key = "names";     allocate(d%e(3)%a, source=names)
  allocate(inner%e(2)); inner%n = 1
  inner%e(1)%key = "hist"; allocate(inner%e(1)%a, source=hist)
  d%e(4)%key = "bounds";    allocate(d%e(4)%v, source=inner)   ! nested dict in class(*)
  d%n = 4
  ! 1. deep copy by intrinsic assignment, including the nested dict
  copy = d
  select type (x => d%e(4)%v)
  type is (dict_t)
     select type (h => x%e(1)%a); type is (real(8)); h(1) = 99; end select
  end select
  select type (x => copy%e(4)%v)
  type is (dict_t)
     select type (h => x%e(1)%a); type is (real(8)); call chk(h(1) == 1d0, "nested deep copy"); end select
  class default; call chk(.false., "nested type lost")
  end select
  ! 2. strings: scalar and array
  select type (s => copy%e(2)%v); type is (character(*)); call chk(s == "trs2" .and. len(s) == 4, &
       "string scalar")
  class default; call chk(.false., "string scalar type")
  end select
  select type (s => copy%e(3)%a); type is (character(*)); call chk(size(s) == 2 .and. s(2) == "beta " &
       .and. len(s) == 5, "string array")
  class default; call chk(.false., "string array type")
  end select
  ! 3. move_alloc keeps payload when growing the entry array
  allocate(grown(8))
  call move_alloc(d%e(1)%v, grown(1)%v)
  call move_alloc(d%e(3)%a, grown(3)%a)
  select type (r => grown(1)%v); type is (real(8)); call chk(r == 1d-6, "move_alloc scalar"); end select
  select type (s => grown(3)%a); type is (character(*)); call chk(s(1) == "alpha", "move_alloc string array"); end select
  ! 4. wrong-type access is detectable
  select type (r => copy%e(1)%v)
  type is (integer)
     call chk(.false., "wrong type matched")
  class default
     call chk(.true., "x")
  end select
  ! 5. same_type_as / extends_type_of on class(*)
  call chk(same_type_as(copy%e(4)%v, inner), "same_type_as dict")
  ! 6. modify nested in place through select type
  select type (x => copy%e(4)%v)
  type is (dict_t); x%n = 2; x%e(2)%key = "converged"; allocate(x%e(2)%v, source=.true.)
  end select
  select type (x => copy%e(4)%v); type is (dict_t); call chk(x%n == 2, "nested in-place update"); end select
  if (fails == 0) then
     print *, "PROBE OK"
  else
     print *, "PROBE FAILURES:", fails
  end if
contains
  subroutine chk(c, label)
    logical, intent(in) :: c
    character(*), intent(in) :: label
    if (.not. c) then
       fails = fails + 1; print *, "FAIL ", label
    end if
  end subroutine
end program
