---
name: calling-convention
description: How public MsPoly routines take arguments and return reports with msp_dict, msp_args and msp_info. Use when adding or changing a public routine, a solver, or any routine that reads parameters or writes results.
---

# Calling convention

## Signatures

- Matrices and other heavy objects are typed arguments. They are never
  dictionary values.
- Everything else comes in through one `args` dictionary and goes out
  through one `info` dictionary. Dictionary values are base types
  (logical, integer, real, complex, string), rank-1 arrays of those, or
  nested dictionaries.
- Low-level operations called in user loops (multiply, axpy, scale) also
  keep a few typed optional scalars (`alpha`, `beta`, `threshold`).

## Public routine pattern

```fortran
subroutine msp_density(H, K, args, info, context)
  type(msp_matrix), intent(in) :: H
  type(msp_matrix), intent(inout) :: K
  type(msp_dict), intent(in) :: args
  type(msp_dict), intent(inout), target :: info
  type(msp_dict), intent(in), target, optional :: context
  type(msp_dict), target :: tracked
  type(msp_args) :: a
  type(msp_info) :: r

  call copy(tracked, args)        ! usage marks go on the copy
  call a%open(tracked, context)
  call r%open(info)
  call density_impl(H, K, a, r)
  call a%finish(r)                ! reports unused_arguments
end subroutine msp_density
```

All the logic lives in `*_impl(..., args, info)`, which takes `msp_args`
and `msp_info` views. Library code calls the `_impl` routines directly.

## Reading parameters

- Read every parameter once, on entry, into typed local variables. Never
  read the dictionary inside an iteration.
- Put the default at the read site:
  `call args%get("max_iterations", n, default=100)`.
- Choose the lookup per key; it is the reader's decision:
  - **local** (the default) reads only this routine's scope. Use it for
    `max_iterations`, `method`, and anything specific to this routine.
  - **scoped** (`scoped=.true.`) also searches each enclosing scope and
    then the context. Use it for thresholds and tolerances, which a user
    should be able to set once at the top.
- Never reject unknown keys. Typos surface as `unused_arguments`.

## Nested calls

```fortran
call args%enter("bounds", inner_args)
call info%enter("bounds", inner_info)
call bounds_impl(H, inner_args, inner_info)
```

Name the section after the nested routine's role, not its implementation.
The nested routine's results appear under that section of the report.

## Errors

A failure that only some ranks can see is agreed on before going on:

```fortran
type(msp_dict) :: err
if (local_problem) then
   call err%set("message", "matrix is not square")
   call err%set("rows", n)
end if
call grid%agree(err, stat, errmsg)   ! same outcome and record on all ranks
```

## Reports

- Write results with `info%set` and per-iteration histories with
  `info%append` (for example `energy_history`).
- Every iterative routine reports `converged`, `iterations` and the
  quantity it converged on. Non-convergence is reported, never silent.
