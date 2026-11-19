!> The density matrix of a real molecule. tblite builds the GFN2-xTB core
!> Hamiltonian H and the overlap matrix S of a structure, and MsPoly
!> purifies the projector onto the occupied space of the generalized
!> eigenproblem H C = S C E. The distribution of the MsPoly matrices
!> decides which elements each rank computes, so no element is evaluated
!> twice and fill moves no data.
program tblite_density
  use mpi_f08, only: MPI_Init, MPI_Finalize, MPI_COMM_WORLD
  use mctc_env, only: error_type
  use mctc_io, only: structure_type, read_structure
  use tblite_adjlist, only: adjacency_list, new_adjacency_list
  use tblite_basis_type, only: get_cutoff
  use tblite_cutoff, only: get_lattice_points
  use tblite_integral_type, only: integral_type, new_integral
  use tblite_xtb_calculator, only: xtb_calculator
  use tblite_xtb_gfn2, only: new_gfn2_calculator
  use tblite_xtb_h0, only: get_selfenergy, get_hamiltonian
  use msp_kinds_mod, only: dp
  use msp_dict_mod, only: msp_dict
  use msp_grid_mod, only: msp_grid
  use msp_matrix_mod
  use msp_density_mod, only: msp_density
  implicit none

  !> Entries below this size are not stored. Larger structures have more
  !> of them; for caffeine almost nothing is far enough away to drop.
  real(dp), parameter :: drop = 1.0e-10_dp
  !> Accuracy of the integral cutoffs.
  real(dp), parameter :: accuracy = 0.1_dp

  type(msp_grid), target :: grid
  type(msp_matrix) :: h, s, density
  type(msp_dict) :: args, info
  type(structure_type) :: mol
  type(error_type), allocatable :: error
  real(dp) :: energy, mu

  call MPI_Init()
  call grid%init(MPI_COMM_WORLD)

  call read_structure(mol, "caffeine.xyz", error)
  if (allocated(error)) then
     print '(a)', "mctc_io: " // error%message
     error stop 1
  end if

  !> Fetch the Hamiltonian and the overlap from tblite, then solve for the
  !! density matrix.
  call build_hamiltonian(mol, grid, h)
  call build_overlap(mol, grid, s)

  !> The projector onto the occupied orbitals. Caffeine has 74 valence
  !! electrons; the closed-shell projector fills half of them, per orbital.
  call args%set("trace", 74.0_dp / 2)
  call args%set("method", "trs4")
  call args%set("threshold", drop)
  call args%set("tolerance", 1.0e-10_dp)
  call msp_density(h, density, args, info, overlap=s)
  call grid%log(info)

  energy = info%get_real("energy")
  mu = info%get_real("chemical_potential")
  if (grid%rank == 0) then
     print '(a, i0, a, i0, a)', "GFN2-xTB core Hamiltonian: ", mol%nat, &
          & " atoms, ", h%n, " orbitals"
     print '(a, es21.14)', "electronic energy tr(K H) / Eh ", energy
     print '(a, es21.14)', "chemical potential          / Eh ", mu
  end if

  call h%free()
  call s%free()
  call density%free()
  call grid%free()
  call MPI_Finalize()

contains

  !> h becomes the GFN2-xTB core Hamiltonian of mol, distributed over grid:
  !! each rank computes only the columns its piece of the decomposition
  !! stores.
  subroutine build_hamiltonian(mol, grid, h)
    implicit none
    type(structure_type), intent(in) :: mol
    type(msp_grid), intent(in), target :: grid
    type(msp_matrix), intent(inout) :: h

    type(xtb_calculator) :: calc
    type(integral_type) :: ints
    type(adjacency_list) :: list
    type(error_type), allocatable :: error
    integer, allocatable :: columns(:)
    real(dp), allocatable :: cn(:), selfenergy(:), lattr(:, :)
    real(dp) :: cutoff

    call new_gfn2_calculator(calc, mol, error)
    if (allocated(error)) then
       print '(a)', "tblite: " // error%message
       error stop 1
    end if

    ! The coordination numbers shift the atomic levels.
    allocate(cn(mol%nat), selfenergy(calc%bas%nsh))
    call calc%ncoord%get_cn(mol, cn)
    call get_selfenergy(calc%h0, mol%id, calc%bas%ish_at, calc%bas%nsh_id, &
         & cn=cn, selfenergy=selfenergy)

    ! The neighbour list tblite evaluates the integral blocks over.
    cutoff = get_cutoff(calc%bas, accuracy)
    call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
    call new_adjacency_list(list, mol, lattr, cutoff)

    ! The decomposition gives this rank a band of rows and a band of
    ! columns; the Hamiltonian is symmetric, so its row band is also a band
    ! of columns, and that is the range tblite computes for this rank.
    call h%init(grid, calc%bas%nao)
    allocate(columns(2))
    columns = [h%row0 + 1, h%row0 + h%rows]
    call new_integral(ints, calc%bas%nao, columns)
    call get_hamiltonian(mol, lattr, list, calc%bas, calc%h0, selfenergy, &
         & ints%overlap, ints%dipole, ints%quadrupole, ints%hamiltonian, &
         & columns=columns)
    call fill_matrix(h, ints%hamiltonian)
  end subroutine build_hamiltonian

  !> s becomes the overlap matrix of mol, distributed over grid. The
  !! overlap comes out of the same integral evaluation as the Hamiltonian,
  !! whose atomic levels it does not depend on.
  subroutine build_overlap(mol, grid, s)
    implicit none
    type(structure_type), intent(in) :: mol
    type(msp_grid), intent(in), target :: grid
    type(msp_matrix), intent(inout) :: s

    type(xtb_calculator) :: calc
    type(integral_type) :: ints
    type(adjacency_list) :: list
    type(error_type), allocatable :: error
    integer, allocatable :: columns(:)
    real(dp), allocatable :: selfenergy(:), lattr(:, :)
    real(dp) :: cutoff

    call new_gfn2_calculator(calc, mol, error)
    if (allocated(error)) then
       print '(a)', "tblite: " // error%message
       error stop 1
    end if

    allocate(selfenergy(calc%bas%nsh))
    selfenergy = 0

    cutoff = get_cutoff(calc%bas, accuracy)
    call get_lattice_points(mol%periodic, mol%lattice, cutoff, lattr)
    call new_adjacency_list(list, mol, lattr, cutoff)

    call s%init(grid, calc%bas%nao)
    allocate(columns(2))
    columns = [s%row0 + 1, s%row0 + s%rows]
    call new_integral(ints, calc%bas%nao, columns)
    call get_hamiltonian(mol, lattr, list, calc%bas, calc%h0, selfenergy, &
         & ints%overlap, ints%dipole, ints%quadrupole, ints%hamiltonian, &
         & columns=columns)
    call fill_matrix(s, ints%overlap)
  end subroutine build_overlap

  !> Fill matrix from this rank's columns of values: values(:, c) is the
  !! full column matrix%row0 + c of the matrix, and its transpose is what
  !! the local piece stores, so the entries are read transposed. Fill is
  !! collective, but each rank passes exactly the entries its piece stores,
  !! so nothing moves and no entry is counted twice.
  subroutine fill_matrix(matrix, values)
    implicit none
    type(msp_matrix), intent(inout) :: matrix
    real(dp), intent(in) :: values(:, :)

    integer, allocatable :: rows(:), cols(:)
    real(dp), allocatable :: vals(:)
    integer :: i, j, k

    allocate(rows(matrix%rows * matrix%cols), &
         & cols(matrix%rows * matrix%cols), vals(matrix%rows * matrix%cols))
    k = 0
    do i = matrix%row0 + 1, matrix%row0 + matrix%rows
       do j = matrix%col0 + 1, matrix%col0 + matrix%cols
          if (abs(values(j, i - matrix%row0)) >= drop) then
             k = k + 1
             rows(k) = i
             cols(k) = j
             vals(k) = values(j, i - matrix%row0)
          end if
       end do
    end do
    call matrix%fill(rows(:k), cols(:k), vals(:k))
  end subroutine fill_matrix
end program tblite_density
