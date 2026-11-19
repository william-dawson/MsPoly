# The example needs tblite (github.com/tblite/tblite) for the Hamiltonian.
# Its Fortran modules and its dependency mctc-lib come with a tblite
# install, which this finds through pkg-config. A static tblite pulls in
# BLAS and LAPACK at link time.
find_package(PkgConfig REQUIRED)
pkg_check_modules(TBLITE REQUIRED IMPORTED_TARGET tblite)
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)

add_executable(tblite_density tblite_density.f90)
target_link_libraries(tblite_density PRIVATE
  MsPoly::solvers MsPoly::matrix MsPoly::utils
  PkgConfig::TBLITE BLAS::BLAS LAPACK::LAPACK)
msp_target_flags(tblite_density)
