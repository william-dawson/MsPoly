!> Logical-set insertion, nested report citations and process history.
program test_citations
  use msp_dict_mod, only: msp_dict, msp_info
  use msp_citations_mod, only: msp_cite, msp_get_citations, &
       & msp_clear_citations
  use testing_mod, only: check, finish
  implicit none

  type(msp_dict), target :: report
  type(msp_dict) :: snapshot, again, section, values
  type(msp_info) :: root, inner
  character(:), allocatable :: message
  integer :: stat, k
  character(len=16) :: mode

  call get_command_argument(1, mode)
  if (mode == "conflict") then
     call report%set("citations", 1)
     call root%open(report)
     call root%cite("niklasson2002")
  end if
  call msp_clear_citations()
  call msp_get_citations(snapshot)
  call check(snapshot%num_keys() == 0, "empty initial history")
  call root%open(report)
  call root%enter("density", inner)
  do k = 1, 1000
     call msp_cite("niklasson2002", inner)
  end do
  call report%get("citations", section)
  call check(section%num_keys() == 1 .and. &
       & section%get_logical("niklasson2002") .and. &
       & .not. report%has("density.citations"), &
       & "nested calls deduplicate at the report root")
  call msp_get_citations(snapshot)
  call report%clear()
  call msp_cite("truflandier2026gap", root)
  call msp_cite("lanczos1950iteration")
  call msp_get_citations(again)
  call check(again%num_keys() == 3 .and. snapshot%num_keys() == 1, &
       & "query owns a snapshot; clearing reports retains history")
  call snapshot%include("custom")
  call msp_get_citations(again)
  call check(.not. again%has("custom"), "snapshot changes are independent")
  call msp_clear_citations()
  call msp_get_citations(again)
  call check(again%num_keys() == 0 .and. &
       & report%get_logical("citations.truflandier2026gap"), &
       & "reset retains existing reports")

  call values%include("nested.a")
  call values%include("nested.b")
  call values%set("nested.a", .false.)
  call values%include("nested.a")
  call check(values%get_logical("nested.a") .and. &
       & values%get_logical("nested.b"), "include sets false to true")
  call values%include("bad..key", stat, message)
  call check(stat /= 0 .and. index(message, "invalid key") > 0 .and. &
       & .not. values%has("bad"), "invalid path does not mutate the set")
  call values%set("scalar", 1)
  call values%include("scalar", stat, message)
  call check(stat /= 0 .and. index(message, "scalar") > 0, &
       & "existing nonlogical scalar is rejected")
  call values%include("scalar.child", stat, message)
  call check(stat /= 0, "cannot insert through a scalar")
  call values%set("array", [.true., .false.])
  call values%include("array", stat, message)
  call check(stat /= 0, "logical arrays are not set members")
  call values%include("array.child", stat, message)
  call check(stat /= 0, "cannot insert through an array")
  call msp_clear_citations()
  call finish("test_citations")
end program test_citations
