---
name: testing
description: How MsPoly tests are organised, what a good test looks like, how to measure coverage, and what CI checks. Use when adding or changing code, adding tests, or checking that a change is ready.
---

# Testing

## Shape of the suite

- One directory per library layer: `tests/<layer>/`. Register each test in
  `tests/CMakeLists.txt` with `msp_add_test(<layer> <name> <libraries>)`.
  A `.fypp` test is preprocessed automatically.
- Keep the suite minimal with full coverage:
  - Each check tests one distinct behaviour.
  - Do not repeat a behaviour that another check already covers.
  - Prefer one check over several, as long as its label still says exactly
    what failed.
- Generate tests over value types or kinds with fypp, as
  `tests/utils/test_values.fypp` does, instead of writing one copy per type.
  Every type instance of a generated routine must run.
- Test error paths through `stat`/`errmsg`, and assert on the message
  content that matters (the key or path named).
- Name tests after what they cover (`test_values`, `test_structure`,
  `test_views`). `test_nested` is the worked example of the calling
  convention; keep it readable as documentation.
- Tests use `check(condition, label)` and `finish(name)` from
  `tests/support/testing_mod.f90`.
- MPI tests are registered with `RANKS 1 2 4` and run through `mpiexec`.
  Pass `-DMSPOLY_MPIEXEC_FLAGS=--oversubscribe` where cores are scarce.
- Paths that must stop the program (errors without `stat`, forbidden
  assignments, collective errors) are modes of `tests/utils/test_stops.f90`,
  checked through the message they print.

## Before a change is done

1. **Local:** build and run with gfortran under AddressSanitizer as well
   (`-DCMAKE_Fortran_FLAGS="-fsanitize=address"`); it found bugs plain runs
   missed:
   ```sh
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug && cmake --build build
   ctest --test-dir build --output-on-failure
   ```
2. **Coverage:** run `tools/coverage.sh`.
   - Every reachable line of the changed library sources must be covered.
     The only accepted gaps are gcov artifacts: a line reported unexecuted
     although the next statement in the same block ran, or the end of a
     procedure that always stops.
   - Code that only a GPU build can reach sits inside
     `#:if defined('MSPOLY_GPU')`, so the CPU build has no unreachable
     lines. GPU code is covered by the GPU test runs (label `gpu`) on a GPU
     machine.
   - Code no test can reach is dead: remove it.
3. **All compilers:** CI (`.github/workflows/ci.yml`) runs on every push:
   - gfortran 13 and 14 with Open MPI, each in Debug and Release;
   - Fortran 2008 conformance with `gfortran -std=f2008 -pedantic -Werror`;
   - nvfortran 26.3, 26.5 and 26.9 with HPC-X, on x86-64 and arm64, each in
     Debug and Release. Any nvfortran warning fails the build;
   - coverage, reported in the job summary;
   - the compiler reproducers, for information only.

   Every job must pass. A failure on one compiler only goes through the
   `compiler-bugs` skill. Never weaken a test to make a compiler pass.
