!> Keys, nesting, ordering, copies, hashing and YAML output of msp_dict.
program test_structure
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
       & ieee_positive_inf, ieee_negative_inf
  use, intrinsic :: iso_fortran_env, only: int8
  use msp_kinds_mod, only: dp
  use msp_dict_mod
  use testing_mod, only: check, finish
  implicit none

  call test_keys()
  call test_nesting()
  call test_order_and_removal()
  call test_copies()
  call test_hash()
  call test_pack()
  call test_yaml()
  call finish("test_structure")

contains

  subroutine test_keys()
    type(msp_dict) :: d
    character(len=8), parameter :: bad(6) = &
         & [character(len=8) :: "", "a..b", ".a", "a.", "a-b", "a b"]
    integer :: k, stat, rejected

    rejected = 0
    do k = 1, size(bad)
       call d%set(trim(bad(k)), 1, stat=stat)
       if (stat /= 0) rejected = rejected + 1
    end do
    call check(rejected == size(bad) .and. d%num_keys() == 0, &
         & "invalid keys are rejected and create nothing")
    call d%set("Key_9", 1)
    call d%set("key_9", 2)
    call d%set("9lives", 3)
    call check(d%num_keys() == 3 .and. d%get_int("Key_9") == 1, &
         & "keys are case-sensitive and may start with a digit")
  end subroutine test_keys

  subroutine test_nesting()
    type(msp_dict) :: d, sub, got
    integer :: stat
    character(:), allocatable :: msg

    call d%set("a.b.c", 1)
    call check(d%num_keys() == 1 .and. d%type_of("a.b") == MSP_DICTIONARY &
         & .and. d%get_int("a.b.c") == 1, "dotted paths create sections")
    call sub%set("x", 2)
    call d%set("a.s", sub)
    call d%get("a", got)
    call check(got%get_int("s.x") == 2 .and. got%get_int("b.c") == 1, &
         & "dictionaries are stored and read back as values")
    call d%set("a.b.c.d", 1, stat=stat, errmsg=msg)
    call check(stat /= 0 .and. index(msg, "'c' holds integer") > 0, &
         & "no key below a value")
    call d%get("a.b.c", got, stat=stat)
    call check(stat /= 0, "a value is not read as a dictionary")
    call d%set("a.b", 5)
    call check(d%type_of("a.b") == MSP_INTEGER .and. d%length_of("a.b") == 1 &
         & .and. .not. d%has("a.b.c"), "a value replaces a section")
    call check(d%type_of("none") == MSP_NONE .and. d%length_of("none") == 0 &
         & .and. msp_type_name(MSP_NONE, .false.) == "nothing", &
         & "absent keys have no type and no length")
  end subroutine test_nesting

  subroutine test_order_and_removal()
    type(msp_dict) :: d
    integer, parameter :: n = 64
    character(len=8) :: key
    integer :: k, ok

    !! k*37 mod 64 visits 0..63 once, in scrambled order.
    do k = 0, n - 1
       write (key, '(A,I0)') "k", mod(k * 37, n)
       call d%set(trim(key), mod(k * 37, n))
    end do
    ok = 0
    do k = 1, n - 1
       if (llt(d%key_at(k), d%key_at(k + 1))) ok = ok + 1
    end do
    call check(d%num_keys() == n .and. ok == n - 1, "keys are kept sorted")
    do k = 0, n - 1, 2
       write (key, '(A,I0)') "k", k
       call d%remove(trim(key))
    end do
    ok = 0
    do k = 1, n - 1, 2
       write (key, '(A,I0)') "k", k
       if (d%get_int(trim(key)) == k) ok = ok + 1
    end do
    call check(d%num_keys() == n / 2 .and. ok == n / 2, &
         & "removal keeps the other keys")
    call d%set("s.p", 1)
    call d%set("s.q.r", 1)
    call d%remove("s.p")
    call d%remove("s.q")
    call d%remove("absent")
    call check(d%has("s") .and. .not. d%has("s.q.r") .and. &
         & d%num_keys() == n / 2 + 1, "remove leaves and whole sections")
  end subroutine test_order_and_removal

  subroutine test_copies()
    type(msp_dict) :: a, b

    call a%set("s.t.x", 1)
    call copy(b, a)
    call b%set("s.t.x", 2)
    call check(a%get_int("s.t.x") == 1 .and. b%get_int("s.t.x") == 2, &
         & "copy is deep")
  end subroutine test_copies

  subroutine test_hash()
    type(msp_dict) :: a, b, empty

    call fill(a)
    call b%set("z.list", ["p", "q"])
    call b%set("z.flags", [.true.])
    call fill(b)
    call check(a%hash() == b%hash(), "hash does not depend on insertion order")
    call copy(b, a)
    call b%set("z.list", ["p", "r"])
    call check(a%hash() /= b%hash(), "hash sees a nested string change")
    call copy(b, a)
    call b%set("n", 1.0_dp)
    call check(a%hash() /= b%hash(), "hash sees a type change")
    call copy(b, a)
    call b%set("n", [1])
    call check(a%hash() /= b%hash(), "hash sees a shape change")
    call check(empty%hash() /= a%hash(), "hash sees contents")
  end subroutine test_hash

  subroutine test_pack()
    type(msp_dict) :: a, b, empty
    integer(int8), allocatable :: bytes(:)
    real(dp), allocatable :: none(:)

    allocate(none(0))
    call fill(a)
    call a%set("edge.empty", empty)
    call a%set("edge.none", none)
    call a%set("edge.blank", "")
    call a%pack(bytes)
    call b%set("junk", 1)
    call b%unpack(bytes)
    call check(a%hash() == b%hash() .and. b%length_of("edge.none") == 0 .and. &
         & b%get_string("name") == "trs2", &
         & "pack and unpack round-trip every type, nesting and empty values")
  end subroutine test_pack

  subroutine fill(d)
    type(msp_dict), intent(inout) :: d

    call d%set("n", 1)
    call d%set("x", 0.5_dp)
    call d%set("c", (1.0_dp, 2.0_dp))
    call d%set("flag", .false.)
    call d%set("name", "trs2")
    call d%set("z.list", ["p", "q"])
    call d%set("z.flags", [.true.])
    call d%set("z.mu", [0.1_dp])
    call d%set("z.c", [(0.0_dp, 1.0_dp)])
    call d%set("z.i", [1, 2])
  end subroutine fill

  subroutine test_yaml()
    type(msp_dict) :: d, empty
    character(len=2) :: names(2) = ["ab", "c "]
    real(dp), allocatable :: none(:)
    character(len=100), dimension(20) :: lines
    character(len=100), dimension(15) :: expected
    integer :: unit, n, ios

    allocate(none(0))
    call d%set("b", [.true., .false.])
    call d%set("c", [(1.0_dp, -1.0_dp)])
    call d%set("e", none)
    call d%set("empty", empty)
    call d%set("i", [1, -2])
    call d%set("nan", ieee_value(1.0_dp, ieee_quiet_nan))
    call d%set("ninf", ieee_value(1.0_dp, ieee_negative_inf))
    call d%set("pinf", ieee_value(1.0_dp, ieee_positive_inf))
    call d%set("q", 'a"b' // achar(92))
    call d%set("s", names)
    call d%set("sec.flag", .true.)
    call d%set("sec.inner.n", 3)
    call d%set("sec.z", (2.0_dp, 0.0_dp))
    call d%set("x", 0.5_dp)

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
    expected(2) = "c: [{re: 1.00000000000000000E+000, " // &
         & "im: -1.00000000000000000E+000}]"
    expected(3) = "e: []"
    expected(4) = "empty: {}"
    expected(5) = "i: [1, -2]"
    expected(6) = "nan: .nan"
    expected(7) = "ninf: -.inf"
    expected(8) = "pinf: .inf"
    expected(9) = 'q: "a' // achar(92) // '"b' // achar(92) // achar(92) // '"'
    expected(10) = 's: ["ab", "c"]'
    expected(11) = "sec:"
    expected(12) = "  flag: true"
    expected(13) = "  inner:"
    expected(14) = "    n: 3"
    expected(15) = "  z: {re: 2.00000000000000000E+000, " // &
         & "im: 0.00000000000000000E+000}"
    call check(n == 16, "yaml line count")
    if (n == 16) then
       call check(all(lines(1:15) == expected) .and. &
            & lines(16) == "x: 5.00000000000000000E-001", &
            & "yaml output of every type, nesting and special values")
       if (any(lines(1:15) /= expected)) write (*, '(A)') lines(1:n)
    end if
  end subroutine test_yaml
end program test_structure
