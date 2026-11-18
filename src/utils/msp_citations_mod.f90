!> Deduplicated citation keys for calculations performed by this process.
module msp_citations_mod
  use msp_dict_mod, only: msp_dict, msp_info, copy
  implicit none
  private

  public :: msp_cite, msp_get_citations, msp_clear_citations

  type(msp_dict), save :: used

contains

  !> Record a bibliography key once for this process and, when present,
  !> for the calculation's report. Call at method entry, outside loops.
  !> Keys refer to Citations.bib. No MPI collective is performed.
  subroutine msp_cite(key, info)
    character(*), intent(in) :: key
    type(msp_info), intent(in), optional :: info

    call used%include(key)
    if (present(info)) call info%cite(key)
  end subroutine msp_cite

  !> Snapshot of citations recorded so far, as key=true. The caller owns
  !> the result; changing it cannot change the accumulated history.
  !> Query any time, independently on each rank, without a collective.
  subroutine msp_get_citations(citations)
    type(msp_dict), intent(inout) :: citations

    call copy(citations, used)
  end subroutine msp_get_citations

  !> Start a new citation history on this process. Existing snapshots
  !> and per-calculation reports are retained.
  subroutine msp_clear_citations()
    call used%clear()
  end subroutine msp_clear_citations
end module msp_citations_mod
