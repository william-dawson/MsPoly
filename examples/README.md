tblite_density
================================================================================
Builds the GFN2-xTB core Hamiltonian H and overlap matrix S of caffeine with
[tblite](https://github.com/tblite/tblite) — each rank computes only the columns
its piece of the MsPoly decomposition stores — and purifies the projector onto
the occupied orbitals of H C = S C E with `msp_density`.

> cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DMSPOLY_BUILD_EXAMPLES=ON
> cmake --build build
> mpiexec -n 4 build/examples/tblite_density

Needs tblite 0.8 or newer with mctc-lib, found by pkg-config; the example reads
caffeine.xyz from the working directory. MsPoly and the example are MIT; tblite
(LGPL-3.0-or-later) is linked, never vendored; mctc-lib is Apache-2.0.
