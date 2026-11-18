!> Edge cases of msp_dict: aliasing, scale, deep nesting, strict typing,
!> output of every value type.
program test_dict_edges
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_positive_inf, &
       & ieee_negative_inf
  use msp_kinds_mod, only: dp, i64
  use msp_dict_mod
  use testing_mod, only: check, finish
  implicit none

  call test_self_assignment()
  call test_nested_copies_are_independent()
  call test_many_keys()
  call test_deep_nesting()
  call test_key_rules()
  call test_overwrite_rules()
  call test_strict_types()
  call test_hash_types()
  call test_yaml_all_types()
  call finish("test_dict_edges")

contains

  subroutine test_self_assignment()
    type(msp_dict) :: a

    call a%set("x", 1)
    call a%set("s.y", 2.0_dp)
    a = a
    call check(a%num_keys() == 2 .and. a%get_int("x") == 1 .and. &
         & a%get_real("s.y") == 2.0_dp, "a = a keeps contents")
  end subroutine test_self_assignment

  subroutine test_nested_copies_are_independent()
    type(msp_dict) :: a, b, sub

    call a%set("s.t.x", 1)
    b = a
    call b%set("s.t.x", 2)
    call b%set("s.t.y", 3)
    call check(a%get_int("s.t.x") == 1 .and. .not. a%has("s.t.y"), &
         & "assignment copies nested dictionaries")
    call a%get("s", sub)
    call sub%set("t.x", 5)
    call check(a%get_int("s.t.x") == 1, "get of a section returns a copy")
    call a%set("copy", a)
    call check(a%get_int("copy.s.t.x") == 1 .and. a%num_keys() == 2, &
         & "a dictionary can store a copy of itself")
    call check(a%hash() /= b%hash(), "different contents, different hash")
    b = a
    call check(a%hash() == b%hash(), "copy hashes equal")
  end subroutine test_nested_copies_are_independent

  subroutine test_many_keys()
    type(msp_dict) :: d
    integer, parameter :: n = 500
    character(len=16) :: key
    integer :: k, j, sorted_ok, found_ok

    !! Insert in a scrambled order (k*263 mod n visits every value once).
    do k = 0, n - 1
       j = mod(k * 263, n)
       write (key, '(A,I0)') "key_", j
       call d%set(trim(key), j)
    end do
    call check(d%num_keys() == n, "500 keys stored")
    sorted_ok = 0
    do k = 1, n - 1
       if (llt(d%key_at(k), d%key_at(k + 1))) sorted_ok = sorted_ok + 1
    end do
    call check(sorted_ok == n - 1, "keys kept in sorted order")
    found_ok = 0
    do j = 0, n - 1
       write (key, '(A,I0)') "key_", j
       if (d%get_int(trim(key)) == j) found_ok = found_ok + 1
    end do
    call check(found_ok == n, "every key retrievable")
    do j = 0, n - 1, 2
       write (key, '(A,I0)') "key_", j
       call d%remove(trim(key))
    end do
    found_ok = 0
    do j = 1, n - 1, 2
       write (key, '(A,I0)') "key_", j
       if (d%get_int(trim(key)) == j) found_ok = found_ok + 1
    end do
    call check(d%num_keys() == n / 2 .and. found_ok == n / 2, &
         & "removing half keeps the rest")
  end subroutine test_many_keys

  subroutine test_deep_nesting()
    type(msp_dict) :: d
    character(:), allocatable :: path
    integer :: k, unit, ios, lines
    character(len=80) :: line
    character(len=1), parameter :: digits(0:9) = &
         & ["0", "1", "2", "3", "4", "5", "6", "7", "8", "9"]

    !! Not achar(): nvfortran 26.x gives achar(i) length i when assigned
    !! to a deferred-length variable.
    path = "l1"
    do k = 2, 10
       path = path // ".l" // digits(mod(k, 10))
    end do
    call d%set(path, 42)
    call check(d%get_int(path) == 42, "ten levels deep")
    call check(d%type_of("l1.l2.l3") == MSP_DICTIONARY, &
         & "intermediate levels are dictionaries")
    call d%remove(path)
    call check(.not. d%has(path) .and. d%has("l1.l2"), &
         & "removing a deep leaf keeps its sections")

    open (newunit=unit, status='scratch', action='readwrite')
    call d%write_yaml(unit)
    rewind (unit)
    lines = 0
    do
       read (unit, '(A)', iostat=ios) line
       if (ios /= 0) exit
       lines = lines + 1
    end do
    close (unit)
    call check(lines == 9 .and. index(line, "{}") > 0, &
         & "empty innermost section written as {}")
  end subroutine test_deep_nesting

  subroutine test_key_rules()
    type(msp_dict) :: d
    integer :: stat

    call d%set("Mixed_Case_9", 1)
    call d%set("mixed_case_9", 2)
    call check(d%num_keys() == 2 .and. d%get_int("Mixed_Case_9") == 1, &
         & "keys are case-sensitive")
    call d%set("9lives", 3, stat=stat)
    call check(stat == 0, "keys may start with a digit")
    call d%set("a-b", 1, stat=stat)
    call check(stat /= 0, "hyphen rejected")
    call d%set("a.", 1, stat=stat)
    call check(stat /= 0, "trailing dot rejected")
    call d%set(".a", 1, stat=stat)
    call check(stat /= 0, "leading dot rejected")
    call check(.not. d%has("a"), "failed sets leave no value")
  end subroutine test_key_rules

  subroutine test_overwrite_rules()
    type(msp_dict) :: d

    call d%set("bounds.tol", 1.0_dp)
    call d%set("bounds", 5)
    call check(d%type_of("bounds") == MSP_INTEGER .and. &
         & .not. d%has("bounds.tol"), "a value may replace a section")
    call d%set("mu", [1.0_dp, 2.0_dp])
    call d%set("mu", 3.0_dp)
    call check(d%length_of("mu") == 1, "a scalar may replace an array")
  end subroutine test_overwrite_rules

  subroutine test_strict_types()
    type(msp_dict) :: d
    real(dp) :: r
    complex(dp) :: c
    integer :: i, stat
    integer, allocatable :: i1(:)
    logical :: l

    call d%set("z", (1.0_dp, 1.0_dp))
    call d%set("x", 2.0_dp)
    call d%set("flag", .true.)
    call d%set("wide", [1_i64, 9000000000_i64])
    call d%get("z", r, stat=stat)
    call check(stat /= 0, "complex is not read as real")
    call d%get("x", c, stat=stat)
    call check(stat /= 0, "real is not read as complex")
    call d%get("flag", i, stat=stat)
    call check(stat /= 0, "logical is not read as integer")
    call d%get("x", l, default=.false., stat=stat)
    call check(stat /= 0, "a present key of the wrong type is an error " // &
         & "even with a default")
    call d%get("wide", i1, stat=stat)
    call check(stat /= 0, "int64 array overflow on default integer read")
  end subroutine test_strict_types

  subroutine test_hash_types()
    type(msp_dict) :: a, b

    call a%set("v", 1)
    call b%set("v", 1.0_dp)
    call check(a%hash() /= b%hash(), "integer 1 and real 1 hash differently")
    call b%set("v", [1])
    call check(a%hash() /= b%hash(), "scalar and array hash differently")
    call a%clear()
    call b%clear()
    call check(a%hash() == b%hash(), "empty dictionaries hash equally")
  end subroutine test_hash_types

  subroutine test_yaml_all_types()
    type(msp_dict) :: d
    character(len=200), dimension(20) :: lines
    character(len=200) :: expected(8)
    character(len=2) :: names(2) = ["ab", "c "]
    integer :: unit, n, ios

    call d%set("b", [.true., .false.])
    call d%set("i", [1, -2])
    call d%set("e", [real(dp) ::])
    call d%set("c", [(1.0_dp, -1.0_dp)])
    call d%set("s", names)
    call d%set("pinf", ieee_value(1.0_dp, ieee_positive_inf))
    call d%set("ninf", ieee_value(1.0_dp, ieee_negative_inf))
    call d%set("t", "x ")

    open (newunit=unit, status='scratch', action='readwrite')
    call d%write_yaml(unit)
    rewind (unit)
    n = 0
    do
       read (unit, '(A)', iostat=ios) lines(n + 1)
       if (ios /= 0) exit
       n = n + 1
    end do
    close (unit)
    expected(1) = "b: [true, false]"
    expected(2) = "c: [{re: 1.00000000000000000E+000, im: -1.00000000000000000E+000}]"
    expected(3) = "e: []"
    expected(4) = "i: [1, -2]"
    expected(5) = "ninf: -.inf"
    expected(6) = "pinf: .inf"
    expected(7) = 's: ["ab", "c"]'
    expected(8) = 't: "x "'
    call check(n == 8, "yaml line count")
    if (n == 8) then
       call check(all(lines(1:8) == expected), &
            & "yaml for every type (array strings trimmed, scalars exact)")
       if (any(lines(1:8) /= expected)) write (*, '(A)') lines(1:n)
    end if
  end subroutine test_yaml_all_types
end program test_dict_edges
