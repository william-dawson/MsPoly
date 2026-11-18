# Compiler bug log

A running record of compiler defects found while building MsPoly, with
reproducers and the coding rule each one imposes. Add new entries as they
are found. Keep code comments free of this history.

Compilers in use:
- nvfortran: NVHPC 26.3, 26.5, 26.9 on NVIDIA Grace (aarch64), Ubuntu 24.04
- gfortran: 13.3 (aarch64, Ubuntu 24.04), 16.2 (macOS arm64)

`nvidia/run.sh` runs every `nvidia/r*.f90` with every installed nvfortran
and gfortran 13. `nvidia/REPORT.md` is the draft report to NVIDIA.

## Patterns confirmed to work

`storage_patterns.f90` covers the polymorphic storage design of the matrix
layer. It passes on every compiler listed above, at -O0 and -O2, and under
AddressSanitizer:
- an abstract type with deferred bindings;
- a `class(t), allocatable` component;
- dispatch on a second `class(t)` operand;
- explicit copy into a `class(t), allocatable, intent(out)` dummy;
- `move_alloc` between `class(t)` allocatables;
- `allocate(mold=)` from a component.

Run it again before changing those patterns.

## Hardware measurements

`placement/` measures, on a GB200 node with 4 ranks, the cost of creating
fresh system-memory arrays for the GPU (first touch, prefetch, copy) and
their read bandwidth. Result (2026-11-13): 5 to 10 GB/s to create, and 0.5
to 1.1 TB/s to read, against about 3.4 TB/s for pooled device memory. GPU
data that is created often therefore lives in device memory.

## Coding rules

The rules these bugs impose are in `.claude/skills/polymorphic-fortran/SKILL.md`.
The workflow for a new bug is in `.claude/skills/compiler-bugs/SKILL.md`.

## nvfortran (all of 26.3, 26.5, 26.9 unless noted)

| ID | Found | Kind | Reproducer | Summary |
|---|---|---|---|---|
| NV1 | 2026-11-12 | wrong code | `nvidia/r1_array_element_assignment.f90`, `r1b_...class.f90` | `alloc(i) = x`: polymorphic allocatable component of `x` copied with the right type but value 0. Scalars, fixed-size arrays and dummy arrays are fine. A larger program crashed in jemalloc. |
| NV2 | 2026-11-12 | wrong code | `nvidia/r4_structure_constructor.f90` | `e = entry_t()` ignores default initialization of other components when the type has a `class(*)` component. W-0155 under `-Mstandard`. |
| NV3 | 2026-11-12 | rejects valid | `nvidia/r3a`–`r3c`, `r3e` (control `r3d`) | `allocate(a, source=[real(8) :: 1, 2])` or `mold=` with a type-spec array constructor: S-0034 syntax error, for any allocate-object. |
| NV4 | 2026-11-12 | rejects valid | `nvidia/r2a_move_alloc_local.f90` (controls `r2b`, `r2c`) | `move_alloc` between local `class(*)` arrays: S-0074. Scalars and components work. |
| NV5 | 2026-11-12 | wrong code, -O2 only | `nvidia/r5_typed_element_to_class_star.f90` | Inside `select type (x => en%a)`, passing `x(k)` to a `class(*)` dummy loses the type at `-O2` (`class default` taken). Correct at `-O0`. |
| NV6 | 2026-11-12 | wrong code | `nvidia/r6_deferred_length_self_concat.f90`, `r6c_achar_length.f90` (control `r6b`) | `s = achar(i)` with a run-time `i` and deferred-length `s` gives `len(s) == i`. `len(achar(i))` itself is 1. |
| NV7 | 2026-11-12 | false positive, `-Mchkptr`, **unreduced** | `nvidia/r7_mchkptr_unreduced/test_dict_edges.f90` | "Null pointer for zs$array" when a zero-size array constructor is passed to a type-bound generic, but only after earlier procedures passed non-empty constructors. `r7_mchkptr_zero_size.f90` and `r7b_mchkptr_type_bound.f90` are reduction attempts that do not trigger it. To reproduce: build MsPoly at that commit as Debug with nvfortran, then run that file linked to `libmspoly_core` and `libmsp_testing`. |
| NV8 | 2026-11-12 | documented behaviour | — | `'\'` in a string literal is an escape unless `-Mbackslash`/`-Mstandard`. Not reported. |
| NV9 | 2026-11-13 | wrong code | `nvidia/r8_self_assignment_class_star.f90` (control `r8b`) | Defined assignment `a = a` for a type holding an allocatable array of entries with `class(*)` components, implemented as copy-then-replace: wrong result at -O2 (26.3, 26.5, 26.9), no result at -O0 (26.3, 26.5). Same assignment with a plain allocatable component is correct. |
| NV10 | 2026-11-13 | internal compiler error | `nvidia/r9_transfer_deferred_length_mold.f90` | `s = transfer(bytes, s, n)` with an allocated deferred-length character array `s` as mold: F-0000 "Deferred-length character symbol must have descriptor" (26.3, 26.5, 26.9). |
| NV11 | 2026-11-13 | crash | `nvidia/r10_source_typed_associate.f90` | `allocate(dst%v, source=x)` with `x` the typed associate name of `select type (x => src%v)`, both `class(*)` components: segfault (26.3, 26.9, -O0 and -O2). The `class default` associate name works, the opposite of GF6. |
| NV12 | 2026-11-13 | wrong code | `nvidia/r11_component_self_section.f90` | `h%v = h%v(1:2)` on an allocatable array component does not reallocate: the size stays the same and stale trailing elements remain (26.3, 26.5, 26.9, -O0 and -O2). A plain variable and an allocatable dummy are correct. |
| NV13 | 2026-11-13 | wrong code | `nvidia/r12_finalization.f90` | Final procedures of an extension of an abstract type are not called when `move_alloc` replaces an allocated `class(t)` target, when an allocated `class(t)` actual is passed to an `intent(out)` dummy, or for `class(t)` components of a holder that is deallocated or goes out of scope (26.3, 26.5, 26.9). gfortran 13/16 are correct but warn under `-std=f2008` (f08/0011). |
| NV14 | 2026-11-14 | rejects valid | `nvidia/r13_complex_part_generic.f90` | A complex part designator `z%re` (or `z(1:n)%re`) of a complex array passed to a generic with real and complex specifics: S-0188 "type mismatch" against the complex specific, so the real part is treated as complex (26.3, 26.5, 26.9, x86-64 and arm64, seen in CI). `real(z, 8)` resolves correctly. |
| NV15 | 2026-11-14 | wrong code | `nvidia/r14_mask_reallocation_pack.f90` | A logical mask allocated on assignment (`keep = i <= m .and. j <= m`, `i` and `j` from an `intent(out)` allocatable dummy) and then used in `pack`: `ALLOCATE: 5.6e14 bytes requested` at -O2 in `msp_reduce_dimension` on 2 and 4 ranks (26.3, 26.5, 26.9, x86-64 and arm64, seen in CI); -O0 with `-Mbounds -Mchkptr` is correct. The probe gives a wrong count or exit 127 in some of the CI configurations. Allocating `keep(size(i))` first does not avoid it; selecting the entries with a compacting loop does. |

### NV16: invalid debug metadata in density purification (unreduced)

Confirmed with nvfortran 26.9 on Grace/aarch64, GPU Debug build, at
`-g -O0 -Mbounds -Mstandard -Mchkptr -mp`. Compiling the density module
at 88452d4 emits `!dbg attachment of global variable must be a
DIGlobalVariableExpression` and `warning: ignoring invalid debug info`.
It returns success and all 97 CPU/GPU runtime tests pass. The defect
concerns emitted debug information; no incorrect numerical result was
observed. Release and the CPU CI configurations compile without it.

The retained source and self-checking driver are in
`nvidia/r15_density_debug_unreduced/`. The driver takes a configured
nvfortran GPU Debug build containing the module dependencies. Removing
the call to `msp_cite` from the internal `purify` procedure removes the
diagnostic. Changing the one-line IF into a block IF or naming the string
constant does not. Moving that call into a module-level helper does.
A small standalone reduction using host association, optional derived
arguments and polymorphic saved state does not trigger it on 26.3/26.5/
26.9 at -O0/-O2, or gfortran 16.2 under ASan. Scope beyond the confirmed
full-module configuration remains unknown.

Rule: put this stopping-criterion citation call in a module-level helper;
the internal purification routine calls that helper. Keep the original
case and controls in the outer citation-validation record.

## gfortran

| ID | Found | Versions | Kind | Reproducer | Summary |
|---|---|---|---|---|---|
| GF1 | 2026-11-12 | 13.3, 16.2 | wrong code | `move_alloc_string_array.f90` | `move_alloc` between local `class(*)` arrays holding `character(5)` leaves length 0. |
| GF2 | 2026-11-12 | 13.3, 16.2 | crash | `zero_size_source.f90` | `allocate(a, source=z)` with `class(*)` `a` and zero-size `real(8)` `z` segfaults. |
| GF3 | 2026-11-12 | 13.3, 16.2 at -O2/-O3 | wrong code | `gfortran/r1b_array_element_assignment_class.f90` | `alloc(i) = x`: the `class(base_t)` component's dynamic type degrades from `child_t` to `base_t`. Correct at `-O0`. |
| GF4 | 2026-11-12 | 13.3, 16.2 | wrong code | `gfortran/class_star_component_element.f90` | Passing an element of a `class(*)` character array component to a `class(*)` dummy uses a 1-byte stride: `["ab","cd","ef"]` gives `ab, bc, cd`. |
| GF5 | 2026-11-12 | 13.3 at -O2 | spurious warnings | — | `-Wuninitialized` on deferred-length allocatable locals, and `-Wfunction-elimination` on impure getters. gfortran 16.2 does not warn. |
| GF6 | 2026-11-13 | 13.3, 16.2 | crash | `gfortran/zero_length_class_star_copy.f90` | `allocate(dst%v, source=x)` where `x` is the `class default` associate name of a `class(*)` component holding a zero-length character value: SIGBUS/SIGSEGV. nvfortran 26.9 is correct. |
| GF7 | 2026-11-13 | 16.2 | use-after-free, **unreduced** | `gfortran/gf7_dict_self_assignment_unreduced/test_structure.f90` | `a = a` for `msp_dict`, with a copy-then-replace defined assignment, is a heap use-after-free under AddressSanitizer (freed by the assignment, then touched by the caller). `nvidia/r8_*` and `r8b_*` pass on gfortran. To reproduce: build at that commit with `-fsanitize=address`, then run that file linked to `libmspoly_utils`. |

| GF8 | 2026-11-14 | 16.2 (arm64) | heap-buffer-overflow | `gfortran/character_constructor_length.f90` | A concatenation with run-time `trim(merge(...))` inside a `character(64)` array constructor allocates using the expression length (45) but writes 64 bytes per element at -O0 under ASan. Assigning the header to a fixed-length variable first passes. |
