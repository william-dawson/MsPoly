---
name: polymorphic-fortran
description: Mandatory rules for class(*) values, polymorphic components, allocation, copies, strings and allocatable components that nvfortran 26.x and gfortran 13/16 miscompile. Use before writing or reviewing any code that uses class(*), class(t), SELECT TYPE, MOVE_ALLOC, ALLOCATE with SOURCE=/MOLD=, structure constructors, or deferred-length strings.
---

# Polymorphic Fortran rules

Each rule prevents a confirmed compiler defect that silently gives wrong
results or rejects valid code. The defect IDs refer to
`tools/compiler-probes/BUGS.md`. Code that breaks a rule may pass on one
compiler and fail on the other, or fail only at `-O2`.

## 1. No intrinsic assignment or structure constructors for types with polymorphic components

Covers any derived type with a `class(*)` or `class(t)` component,
including types that contain such a type (NV1, NV2, GF3).

```fortran
! no:  this%e(k) = en        ! array element: payload lost (nvfortran)
! no:  en = entry_t()        ! default initialization ignored (nvfortran)
call copy_entry(this%e(k), en)   ! explicit, recursive deep copy
call reset_entry(en)             ! explicit reset
```

Such types are not assignable: give them a type-bound defined assignment
that stops with a message, and a public explicit copy (`copy(dst, src)` for
`msp_dict`). A defined assignment that copies is unsafe for `a = a` on both
compilers (NV9, GF7). Moving whole arrays of such types with `move_alloc` is
fine.

## 2. Never pass a value to a `class(*)` dummy argument

Cast once with `select type` at the point of use. Below that point, write
typed code (NV5, GF4).

```fortran
select type (x => en%a)
type is (real(dp))
   do k = 1, size(x)
      str = str // format_real(x(k))     ! typed helper, not class(*)
   end do
end select
```

## 3. Allocate `class(*)` arrays by type, then assign

```fortran
allocate(real(dp) :: en%a(size(value)))       ! not allocate(en%a, source=value)
select type (x => en%a)
type is (real(dp))
   x = value
end select
allocate(character(len=len(value)) :: en%a(size(value)))   ! strings
```

A scalar `allocate(en%v, source=value)` is fine when `value` is an ordinary
(non-associate) variable of a non-character type. Never use `source=` with
a `select type` associate name: nvfortran crashes on the typed name (NV11),
gfortran on the `class default` name holding an empty string (GF6). To copy
a `class(*)` value, branch on every type and use a typed `allocate` plus
assignment, as for arrays. Strings always use a typed `allocate` with an
explicit length.

## 4. `SOURCE=` and `MOLD=` take a variable

Never pass an array constructor with a type-spec (`[real(dp) :: ...]`).
nvfortran rejects it for every allocate-object (NV3).

## 5. No `move_alloc` on `class(*)` values

Copy them instead (GF1, NV4). Moving arrays of derived types that contain
`class(*)` components is fine.

## 6. Strings

- Never assign `achar(i)` with a run-time `i` to a deferred-length
  variable. Use a `character(len=1)` lookup table or an internal write
  (NV6).
- Write a backslash as a named constant, `achar(92)`. Never use `'\'` in a
  literal (NV8).
- Never use a deferred-length character variable as the `mold` of
  `transfer`. Declare a fixed-length local inside a `block` instead (NV10).

## 7. Shrinking allocatable components

Never assign a section of an allocatable component to itself
(`x%a = x%a(1:n)`). nvfortran keeps the old size (NV12). Copy the section to
a temporary and `move_alloc` it back.

## 8. No final procedures

nvfortran skips finalization in most polymorphic cases (NV13). Resources
other than memory Fortran manages, such as device memory, are released
explicitly through a `release` binding, and a leak counter in the tests
checks that every acquisition is released.

## 9. Zero-size arrays

Pass zero-size arrays as allocated variables. Avoid zero-size constructors
and zero-size explicit-shape declarations: `-Mchkptr` false positive NV7,
warning W-0435.

## Before relying on a new polymorphic pattern

For density purification's stopping-criterion citation, call a module
helper from the internal procedure; keep the call to `msp_cite` in that
helper. The direct call emits invalid debug metadata with nvfortran 26.9
in the GPU Debug configuration (NV16, unreduced). Runtime checks pass;
changing IF syntax or naming the citation string does not remove it.

Write a small self-checking probe (see the `compiler-bugs` skill). Run it
with nvfortran 26.9 at `-O0` and `-O2`, and with gfortran at `-O0` and
`-O2`, before using the pattern in library code.

Character array constructors with a fixed type-spec must use fixed-length
variables for run-time concatenations. Build the concatenation in a
`character(n)` variable first (GF8).
