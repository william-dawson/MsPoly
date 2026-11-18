# nvfortran bug reports (drafted 2026-11-12)

Nine defects in nvfortran, reproduced with NVHPC **26.3, 26.5 and 26.9**
(26.9 is the newest installed), at
`-O0`, `-O2` and `-O0 -Mstandard`. All the test
programs are valid Fortran 2008. gfortran 13.3 and 16.2 compile and run
them correctly at `-O0`; one gfortran `-O2` issue is noted at the end.

- Platform: NVIDIA Grace (aarch64, neoverse-v2), Ubuntu 24.04.5 LTS,
  kernel 7.0.0-1025-nvidia-64k
- Compiler: `nvfortran 26.9-0 linuxarm64 target on aarch64 Linux -tp
  neoverse-v2` (module `nvhpc/26.9`)
- Each reproducer is self-checking: it prints `RESULT: all passed` (exit
  0) or `RESULT: ...failure(s)` (`ERROR STOP 1`).
- `run.sh` runs every reproducer with every installed version, and
  `results.txt` is its output.

| # | Kind | Files | 26.3 | 26.5 | 26.9 |
|---|---|---|---|---|---|
| 1 | Wrong code: polymorphic component lost when assigning into an allocatable array element | `r1_array_element_assignment.f90`, `r1b_array_element_assignment_class.f90` | wrong | wrong | wrong |
| 2 | Wrong code: default initialization ignored in a structure constructor when the type has a `class(*)` component | `r4_structure_constructor.f90` | wrong | wrong | wrong |
| 3 | Rejects valid code: `SOURCE=`/`MOLD=` with a type-spec array constructor | `r3a`, `r3b`, `r3c`, `r3e` (control: `r3d` passes) | error | error | error |
| 4 | Rejects valid code: `MOVE_ALLOC` with local unlimited polymorphic arrays | `r2a_move_alloc_local.f90` (controls: `r2b`, `r2c` pass) | error | error | error |
| 5 | Wrong code at -O2: typed element passed to a `class(*)` dummy inside SELECT TYPE loses its type | `r5_typed_element_to_class_star.f90` | wrong | wrong | wrong |
| 6 | Wrong code: `s = achar(i)` with run-time `i` gives `len(s) == i` | `r6_deferred_length_self_concat.f90`, `r6c_achar_length.f90` (control: `r6b`) | wrong | wrong | wrong |
| 7 | Wrong code: defined assignment `a = a` for a type with `class(*)` components | `r8_self_assignment_class_star.f90` (control: `r8b`) | wrong | wrong | wrong |
| 8 | Internal compiler error: `transfer` with a deferred-length character array mold | `r9_transfer_deferred_length_mold.f90` | ICE | ICE | ICE |
| 9 | Wrong code: `h%v = h%v(1:2)` does not reallocate an allocatable component | `r11_component_self_section.f90` | wrong | wrong | wrong |

---

## 1. Intrinsic assignment into an allocatable array element drops a polymorphic allocatable component

**Severity:** silent wrong results.

For a derived type with a `class(*), allocatable` (or `class(base_t),
allocatable`) component, `alloc(1) = a` leaves the component of `alloc(1)`
allocated with the right dynamic type but the wrong value (`0.0` instead of
`1.5`). All of these work:
- assigning to a scalar;
- assigning to an element of a fixed-size array;
- assigning to an element of an assumed-shape dummy;
- the same type with a non-polymorphic `real(8), allocatable` component.

Fortran 2008 10.2.1.3 requires a deep copy of allocatable components,
including polymorphic ones, so `alloc(1)%v` must have the dynamic type and
value of `a%v`.

```fortran
type :: box_t
   class(*), allocatable :: v
end type box_t
type(box_t) :: a
type(box_t), allocatable :: alloc(:)
allocate(a%v, source=1.5d0)
allocate(alloc(2))
alloc(1) = a          ! alloc(1)%v is real(8) but its value is 0.0
```

```text
$ nvfortran -O0 r1_array_element_assignment.f90 -o r1 && ./r1
scalar = scalar                             PASS
fixed-size array element = scalar           PASS
allocatable array element = scalar          FAIL: value is  0.000E+00
dummy array element = scalar                PASS
control: real component, alloc array elementPASS
RESULT: 1 failure(s)
ERROR STOP 1
```

`r1b` is the same case with a `class(base_t), allocatable` component
holding an extended type `child_t`. nvfortran fails it too.

In a larger program that assigns such derived types into elements of an
allocatable array, the same pattern crashed inside jemalloc (`SIGSEGV` in
`rtree_metadata_read`) instead of only producing wrong values.

Workaround: copy components explicitly (`allocate(dst%v, source=src%v)`,
or a typed `allocate` plus `select type`) instead of intrinsic assignment.

## 2. Structure constructor ignores default initialization when the type has a `class(*)` component

**Severity:** silent wrong results. With `-Mstandard` there is also a
misleading warning: `NVFORTRAN-W-0155-Mismatched data type for member used`.

Components omitted from a structure constructor take their default
initialization (Fortran 2008 4.5.10), so after `e = entry_t()` the logical
`used` must be `.false.`. nvfortran leaves it `.true.`. The same type
without the `class(*)` component behaves correctly.

```fortran
type :: entry_t
   character(:), allocatable :: key
   class(*), allocatable :: v
   logical :: used = .false.
end type entry_t
type(entry_t) :: e
allocate(e%v, source=1.5d0)
e%used = .true.
e = entry_t()         ! e%used stays .true.
```

```text
$ nvfortran -O0 r4_structure_constructor.f90 -o r4 && ./r4
entry_t(): v allocated = F  used = T
control plain_t(): v allocated = F  used = F
RESULT: 1 failure(s)
ERROR STOP 1
```

## 3. `SOURCE=` / `MOLD=` with a type-spec array constructor: syntax error

**Severity:** rejects valid code. This affects any allocate-object, not
only polymorphic ones.

```fortran
real(8), allocatable :: a(:)
allocate(a, source=[real(8) :: 1, 2])   ! NVFORTRAN-S-0034 Syntax error at or near =
allocate(a, mold=[real(8) :: 1, 2])     ! same
allocate(a, source=[1d0, 2d0])          ! accepted (no type-spec)
```

```text
$ nvfortran -O0 r3b_source_typed_constructor_plain.f90
NVFORTRAN-S-0034-Syntax error at or near = (r3b_source_typed_constructor_plain.f90: 5)
  0 inform,   0 warnings,   1 severes, 0 fatal for r3b
```

Fortran 2008 R626/R631 allows any expression as *source-expr*, and R469
allows a type-spec in an array constructor. The polymorphic (`r3a`) and
character (`r3c`) variants fail the same way.

Workaround: assign the constructor to a variable first and pass the
variable to `SOURCE=`.

## 4. `MOVE_ALLOC` rejects local unlimited polymorphic arrays

**Severity:** rejects valid code.

```fortran
class(*), allocatable :: a(:), b(:)
real(8) :: z(3) = [1d0, 2d0, 3d0]
allocate(a, source=z)
call move_alloc(a, b)    ! NVFORTRAN-S-0074
```

```text
$ nvfortran -O0 r2a_move_alloc_local.f90
NVFORTRAN-S-0074-Illegal number or type of arguments to move_alloc - keyword argument to (r2a_move_alloc_local.f90: 8)
  0 inform,   0 warnings,   1 severes, 0 fatal for r2a
```

Fortran 2008 13.7.118 allows FROM and TO to be unlimited polymorphic
(TO must be type compatible with FROM). nvfortran accepts these controls:
- the same call on `class(*), allocatable` scalars (`r2b`);
- the same call on rank-1 `class(*)` arrays that are derived-type
  components (`r2c`).

The diagnostic text ("keyword argument to") also looks truncated.

## 5. Typed array element passed to a `class(*)` dummy inside SELECT TYPE (-O2)

**Severity:** silent wrong results at `-O2`. `-O0` is correct.

```fortran
type :: entry_t
   class(*), allocatable :: a(:)
end type entry_t
...
select type (x => en%a)          ! en is a type(entry_t) dummy
type is (real(8))
   do k = 1, size(x)
      str = str // " " // describe(x(k))   ! describe(x) has class(*), intent(in) :: x
   end do
end select
```

Inside `describe`, `select type` must match `type is (real(8))`. At `-O2`
it takes `class default`:

```text
$ nvfortran -O2 r5_typed_element_to_class_star.f90 -o r5 && ./r5
got: null null
RESULT: FAIL (expected ' 1.50 2.50')
ERROR STOP 1
$ nvfortran -O0 r5_typed_element_to_class_star.f90 -o r5 && ./r5
got: 1.50 2.50
RESULT: all passed
```

## 6. `ACHAR` with a run-time argument assigned to a deferred-length variable

**Severity:** silent wrong results.

```fortran
integer :: i
character(:), allocatable :: s
i = 50
s = achar(i)      ! len(s) is 50; must be 1 (Fortran 2008 13.7.3)
```

```text
$ nvfortran -O0 r6c_achar_length.f90 -o r6c && ./r6c
len(achar(50))=1  len(achar(i)) direct=1  len(s=achar(i))=50
RESULT: FAIL (ACHAR result must have length 1)
ERROR STOP 1
```

In a loop, `path = path // ".l" // achar(iachar('0') + k)` builds
`l1.l2<48 blanks>.l3<49 blanks>...` (`r6`). The same loop using a
`character(len=1)` table instead of `achar` is correct (`r6b`).

## 7. Defined assignment `a = a` with `class(*)` components

**Severity:** silent wrong results at `-O2`. 26.3 and 26.5 give no result at
`-O0`.

The type holds an allocatable array of entries, each with a
`class(*), allocatable` component. Its defined assignment builds the copy,
then deallocates the old array and moves the copy in. After `a = a` the
values are wrong. The same pattern with a plain `integer, allocatable :: v(:)`
component is correct (`r8b`). Fortran 2008 10.2.1.4 passes the right-hand
side as `(a)`, a value, so the result must equal the original.

```text
$ nvfortran -O0 r8_self_assignment_class_star.f90 -o r8 && ./r8
RESULT: FAIL
ERROR STOP 1
```

## 8. `TRANSFER` with a deferred-length character array mold: ICE

**Severity:** rejects valid code with an internal compiler error.

```fortran
integer(int8) :: bytes(6) = int([97, 98, 99, 100, 101, 102], int8)
character(:), allocatable :: s(:)
allocate(character(len=3) :: s(2))
s = transfer(bytes, s, 2)
```

```text
$ nvfortran r9_transfer_deferred_length_mold.f90
NVFORTRAN-F-0000-Internal compiler error. Deferred-length character symbol must have descriptor ...
```

gfortran 13.3 and 16.2 compile it and print `RESULT: all passed`.

## 9. Self-section assignment to an allocatable component

**Severity:** silent wrong results.

```fortran
type :: holder_t
   integer, allocatable :: v(:)
end type holder_t
type(holder_t) :: h
h%v = [1, 2, 3, 4]
h%v = h%v(1:2)          ! size(h%v) must become 2 (Fortran 2008 7.2.1.3)
```

```text
$ nvfortran r11_component_self_section.f90 -o r11 && ./r11
variable:  size 2
component: size 4
dummy:     size 2
RESULT: FAIL
```

The same assignment to a plain allocatable variable and to an allocatable
dummy reallocates correctly.

Also observed, not yet reduced: `-Mchkptr` reports "Null pointer for
zs$array" for a valid zero-size array constructor actual argument in a
larger program. See `../BUGS.md` (NV7).

---

## Complex part designator resolves to the complex specific of a generic

`r13_complex_part_generic.f90` (reduced):

```fortran
interface pick
   module procedure pick_r, pick_c   ! real(8) and complex(8) array dummies
end interface pick
complex(8), dimension(3) :: z
print *, pick(z%re), pick(z(1:2)%re)   ! must call pick_r twice
```

```
$ nvfortran r13_complex_part_generic.f90
NVFORTRAN-S-0188-Argument number 1 to pick_c: type mismatch
```

F2008 6.4.4 (complex parts): a complex part designator `%re` has type real
and the kind of the complex parent, so generic resolution (12.5.5.2) must
select `pick_r`. Observed with 26.3, 26.5 and 26.9 on x86-64 and arm64.
gfortran 13 and 16 resolve it to `pick_r`. Workaround: `real(z, 8)`.

---

## Logical mask allocated on assignment, then used in `pack`

`r14_mask_reallocation_pack.f90` (from MsPoly's `msp_reduce_dimension`):

```fortran
call entries(n, i, j, v)              ! intent(out) allocatable dummies
keep = i <= m .and. j <= m            ! keep unallocated before this
n_kept = take(pack(i, keep), pack(j, keep), real(pack(v, keep), 8))
```

```
$ nvfortran -O2 ...   (MsPoly test_solvers, 2 MPI ranks)
0: ALLOCATE: 562945456924416 bytes requested
```

F2008 7.2.1.3: an unallocated allocatable variable on the left of an
intrinsic assignment is allocated with the shape of the expression, so
`keep` has the size of `i`. Observed in the full program with 26.3, 26.5
and 26.9 on x86-64 and arm64 at -O2, where -O0 is correct. The standalone
probe gives a wrong count (or exits with 127) in some of those
configurations. Allocating `keep` before the assignment does not help;
selecting the entries with a loop does.

---

## NV16: invalid debug metadata in an internal purification procedure

Unreduced source: `r15_density_debug_unreduced/msp_density_mod.f90`,
retained from MsPoly 88452d4. It needs module files from an nvfortran
GPU Debug build with ScaLAPACK and cuSOLVERMp. Reproduce with:

```
python3 r15_density_debug_unreduced/reproduce.py /path/to/build-Debug
```

The driver reuses the build's include paths and flags, compiles into an
isolated temporary directory and checks the compiler diagnostic. The
original command's flags include:

```
nvfortran -g -O0 -Mbounds -Mstandard -Mbounds -Mchkptr -g -mp \
  -I/path/to/build-Debug/include [MPI/CUDA include paths] \
  -c msp_density_mod.f90
```

On Grace/aarch64 with NVHPC 26.9:

```
!dbg attachment of global variable must be a DIGlobalVariableExpression
warning: ignoring invalid debug info in /tmp/nvfortranJWOzPTG0vapZ.ll
```

The compiler returns success, and all 97 Debug CPU/GPU runtime tests pass
with bounds and pointer checking. The warning causes our warning-free
build gate to reject the build. No source error or numerical corruption
was detected by gfortran bounds/ASan checks or these runtime tests.

The trigger is the call to use-associated `msp_cite` from internal
`purify`, passing the parent procedure's `info` dummy by host association
(Fortran 2008, 16.5.1.4). Removing this call eliminates the warning.
Using a block IF or a named string constant does not. A module-level
helper making the same call eliminates it without changing the arguments
or when the citation is registered. A small standalone reduction does
not reproduce it on 26.3/26.5/26.9 at -O0/-O2, so this is retained as an
unreduced full-module case; those versions are not claimed affected based
on that reduction alone.

---

## Appendix: gfortran issues found at the same time (for GCC Bugzilla)

These are not nvfortran issues. They are recorded here because they came up
in the same work. The full list is in `../BUGS.md`.

- `r1b_array_element_assignment_class.f90` fails with gfortran 13.3 and
  16.2 at `-O2`/`-O3`, but passes at `-O0`. The dynamic type of the copied
  `class(base_t)` component degrades from `child_t` to `base_t`.
- `../move_alloc_string_array.f90`: `move_alloc` between local
  `class(*)` arrays holding `character(5)` values leaves the result with
  length 0 (gfortran 13.3, 16.2).
- `../zero_size_source.f90`: `allocate(a, source=z)` with a zero-size
  `real(8)` array `z` and `class(*)` array `a` segfaults (gfortran 13.3,
  16.2).
- `../gfortran/class_star_component_element.f90`: an element of a `class(*)`
  character array component passed to a `class(*)` dummy uses a 1-byte
  stride (gfortran 13.3, 16.2).
