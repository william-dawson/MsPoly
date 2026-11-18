# MsPoly

Pure Fortran 2008 library for distributed sparse matrix functions (real and
complex, CPU and GPU).

Skills in `.claude/skills/` hold the development rules:

- `fortran-style`: language level, naming, comments, errors, layout
- `polymorphic-fortran`: mandatory rules for `class(*)` and polymorphic
  code
- `calling-convention`: `args`/`info` dictionaries for public routines
- `testing`: what a test suite looks like, coverage, and what CI checks
- `compiler-bugs`: from a suspected compiler bug to a recorded rule

Compiler defects and their reproducers are logged in
`tools/compiler-probes/BUGS.md`.
