!> msp_args lookup rules and usage tracking, and msp_info reports.
program test_views
  use msp_kinds_mod, only: dp
  use msp_dict_mod
  use testing_mod, only: check, finish
  implicit none

  call test_lookup()
  call test_context()
  call test_usage()
  call test_strict()
  call test_errors()
  call test_reports()
  call finish("test_views")

contains

  subroutine test_lookup()
    type(msp_dict), target :: root
    type(msp_args) :: top, mid, deep
    real(dp) :: r
    integer :: i

    call root%set("tol", 1.0e-3_dp)
    call root%set("n", 100)
    call root%set("bounds.tol", 1.0e-4_dp)
    call top%open(root)
    call top%enter("bounds", mid)
    call mid%enter("power", deep)
    call deep%get("tol", r, scoped=.true.)
    call check(r == 1.0e-4_dp, "scoped: nearest enclosing scope wins")
    call deep%get("n", i, default=5, scoped=.true.)
    call check(i == 100, "scoped: found at the top")
    call deep%get("n", i, default=5)
    call check(i == 5, "local: enclosing scopes are not searched")
    call check(deep%has("n", scoped=.true.) .and. .not. deep%has("n") .and. &
         & .not. deep%has("missing", scoped=.true.), &
         & "has: looks up keys as get does")
  end subroutine test_lookup

  subroutine test_context()
    type(msp_dict), target :: root, ctx
    type(msp_args) :: top, inner
    real(dp) :: r
    integer :: i

    call ctx%set("tol", 1.0e-6_dp)
    call ctx%set("bounds.n", 7)
    call ctx%set("bounds.tol", 1.0e-7_dp)
    call top%open(root, ctx)
    call top%enter("bounds", inner)
    call inner%get("n", i, default=5)
    call check(i == 7, "local read reaches the context section")
    call root%set("bounds.n", 9)
    call inner%get("n", i, default=5)
    call check(i == 9, "arguments beat the context")
    call root%set("tol", 1.0e-3_dp)
    call inner%get("tol", r, scoped=.true.)
    call check(r == 1.0e-3_dp, &
         & "all argument scopes are searched before the context")
  end subroutine test_context

  subroutine test_usage()
    type(msp_dict) :: caller
    type(msp_dict), target :: root, ctx, report
    type(msp_args) :: top, mid
    type(msp_info) :: info
    character(:), allocatable :: unused(:)
    character(len=8) :: key
    real(dp) :: r
    integer :: k

    call ctx%set("ctx_only", 1)
    call caller%set("tol", 1.0e-3_dp)
    call caller%set("bounds.tol", 1.0e-4_dp)
    do k = 1, 9
       write (key, '(A,I0)') "typo", k
       call caller%set("bounds.sub." // trim(key), k)
    end do
    call copy(root, caller)
    call top%open(root, ctx)
    call info%open(report)
    call top%enter("bounds", mid)
    call mid%get("tol", r, scoped=.true.)
    call mid%get("tol", r, scoped=.true.)
    call top%finish(info)
    call report%get("unused_arguments", unused)
    call check(size(unused) == 10 .and. any(unused == "tol") .and. &
         & any(unused == "bounds.sub.typo9"), &
         & "unread keys reported by path; shadowed and nested keys included")
    call check(.not. any(unused == "ctx_only"), "context keys never reported")
    call caller%unused_keys(unused)
    call check(size(unused) == 11, "the caller's dictionary is not marked")
  end subroutine test_usage

  subroutine test_strict()
    type(msp_dict), target :: root, ctx, report
    type(msp_args) :: top
    type(msp_info) :: info
    integer :: i, stat
    character(:), allocatable :: msg

    call ctx%set("strict_arguments", .true.)
    call root%set("used", 1)
    call top%open(root, ctx)
    call info%open(report)
    call top%get("used", i)
    call top%finish(info, stat)
    call check(stat == 0 .and. .not. report%has("unused_arguments"), &
         & "strict: all arguments used is not an error")
    call root%set("typo", 1)
    call top%finish(info, stat, msg)
    call check(stat /= 0 .and. index(msg, "typo") > 0, &
         & "strict: unused arguments are an error")
  end subroutine test_strict

  subroutine test_errors()
    type(msp_dict), target :: root
    type(msp_args) :: a, b, c
    integer :: i, stat
    character(:), allocatable :: msg

    call root%set("n", 2.5_dp)
    call root%set("x.y.section.inner", 1)
    call a%open(root)
    call a%enter("x", b)
    call b%enter("y", c)
    call c%get("missing", i, stat=stat, errmsg=msg)
    call check(stat /= 0 .and. index(msg, "'x.y.missing'") > 0, &
         & "a missing key is named by its full path")
    call c%get("n", i, scoped=.true., stat=stat, errmsg=msg)
    call check(stat /= 0 .and. index(msg, "'n' holds real") > 0, &
         & "a type error names the key that was found")
    call c%get("section", i, stat=stat, errmsg=msg)
    call check(stat /= 0 .and. index(msg, "holds dictionary") > 0, &
         & "a section is not read as a value")
  end subroutine test_errors

  subroutine test_reports()
    type(msp_dict), target :: report
    type(msp_info) :: info, first, again, deeper

    call info%open(report)
    call info%enter("bounds", first)
    call first%set("iterations", 3)
    call info%enter("bounds", again)
    call again%enter("power", deeper)
    call deeper%set("converged", .true.)
    call info%set("converged", .false.)
    call check(report%get_int("bounds.iterations") == 3 .and. &
         & report%get_logical("bounds.power.converged") .and. &
         & .not. report%get_logical("converged"), &
         & "repeated and nested report views merge into one tree")
  end subroutine test_reports
end program test_views
