---
name: compiler-bugs
description: Workflow for a suspected compiler defect, from first symptom to a recorded reproducer and coding rule. Use when a test passes on one compiler or optimization level and fails on another, when a compiler rejects code that looks valid, or when results look corrupted.
---

# Handling a suspected compiler bug

1. **Confirm it is not our code.**
   - Build with gfortran `-O0 -fcheck=all` and under AddressSanitizer/UBSan.
   - Build with nvfortran `-g -Mbounds -Mchkptr`.
   - If any of these flags our code, fix our code.
2. **Reduce.** Copy the failing pattern into a standalone program of a few
   dozen lines, then remove pieces until the failure disappears.
   - Add controls: the same program with one feature changed that passes.
     They show the boundary of the defect.
   - Time-box this. If it will not reduce, keep the full failing file.
3. **Make it self-checking.** The program prints `RESULT: all passed`, or
   `RESULT: ...` followed by `error stop 1`. Comments state the expected
   behaviour and the Fortran 2008 clause that requires it.
4. **Run the matrix.** Put nvfortran reproducers in
   `tools/compiler-probes/nvidia/rN_<what>.f90` and gfortran ones in
   `tools/compiler-probes/gfortran/`.
   - `nvidia/run.sh [compiler...]` runs all `r*.f90` at `-O0`, `-O2` and
     the compiler's conformance mode.
   - Run it with every nvfortran and gfortran version available. CI runs it
     with each nvfortran version and reports the result.
   - Save the combined output as `nvidia/results.txt`.
5. **Record it** in `tools/compiler-probes/BUGS.md`: ID (`NVn`/`GFn`), date,
   affected versions, kind (wrong code / rejects valid / crash / false
   positive), reproducer, one-line summary. Unreduced bugs are recorded too,
   marked **unreduced**, with exact steps to trigger them.
6. **Report.** For nvfortran, add a section to `nvidia/REPORT.md`: the
   minimal source, the command, the verbatim output, and the standard
   reference.
7. **Derive the rule.** Add or extend a rule in the `polymorphic-fortran`
   skill so the pattern is avoided everywhere, then fix every occurrence in
   `src/`.

Never describe the bug or the workaround in code comments. The code says
what it does; the history lives in `BUGS.md`.
