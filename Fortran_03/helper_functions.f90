!======================================================================
! Helper API: coordinate transforms + evaluate_LRF (external procs)
!   - No module, no saved state
!   - Safe for OpenMP/MPI callers (stateless)
!======================================================================

subroutine coordinate_transformation(new_coordinates, coordinates, coord_format)
  use iso_fortran_env, only : real64
  implicit none
  real(real64), intent(in)  :: coordinates(6)
  character(*),  intent(in) :: coord_format
  real(real64), intent(out) :: new_coordinates(6)

  new_coordinates = coordinates

  if (coord_format == "Euler_ZYZ") then
    new_coordinates(5) = coordinates(5) - 90.0_real64
    new_coordinates(6) = coordinates(6) - 90.0_real64
  else if (coord_format == "Spherical") then
    new_coordinates(5) = 90.0_real64 - coordinates(5)
    new_coordinates(6) = 90.0_real64 - coordinates(6)
  end if
end subroutine coordinate_transformation


subroutine general_coordinates_format(general_coordinates, dim, old_coordinates)
  use iso_fortran_env, only : real64, int32
  implicit none
  integer(int32), intent(in) :: dim
  real(real64),   intent(in) :: old_coordinates(dim)
  real(real64),   intent(out):: general_coordinates(6)

  general_coordinates = 0.0_real64

  select case (dim)
  case (2)
    general_coordinates(1) = old_coordinates(1)  ! R
    general_coordinates(2) = old_coordinates(2)  ! b1
  case (3)
    general_coordinates(1) = old_coordinates(1)  ! R
    general_coordinates(2) = old_coordinates(2)  ! b1
    general_coordinates(5) = old_coordinates(3)  ! c1
  case (4)
    general_coordinates(1:4) = old_coordinates(1:4)  ! R, b1, b2, phi
  case (5)
    general_coordinates(1:5) = old_coordinates(1:5)  ! + c1
  case (6)
    general_coordinates(1:6) = old_coordinates(1:6)  ! + c2
  case default
    ! leave as zeros
  end select
end subroutine general_coordinates_format


subroutine user_coordinates_to_general_coordinates(general_coordinates_ZXZ, xdim, coord_format, user_coordinates)
  use iso_fortran_env, only : real64, int32
  implicit none
  integer(int32), intent(in) :: xdim
  real(real64),   intent(in) :: user_coordinates(xdim)   ! explicit-shape: no explicit interface needed
  character(*),   intent(in) :: coord_format
  real(real64),   intent(out):: general_coordinates_ZXZ(6)

  real(real64) :: general_coordinates(6)

  call general_coordinates_format(general_coordinates, xdim, user_coordinates)
  call coordinate_transformation(general_coordinates_ZXZ, general_coordinates, coord_format)
end subroutine user_coordinates_to_general_coordinates


subroutine evaluate_LRF(total_energy, xdim, coordinates, coord_format, filename)
  use iso_fortran_env,      only : real64, int32
  use Fitting_Constant_v2,  only : get_coeff_zero, get_coeff_index
  use Geometry_Constant_v2, only : get_total_interaction_energy
  implicit none
  real(real64),   intent(out) :: total_energy
  integer(int32), intent(in)  :: xdim
  real(real64),   intent(in)  :: coordinates(xdim)   ! explicit-shape
  character(*),   intent(in)  :: coord_format
  character(*),   intent(in)  :: filename

  real(real64)   :: general_coordinates_ZXZ(6)
  integer(int32) :: i, coeff_index
  real(real64)   :: x1

  coeff_index = get_coeff_index(filename)

  x1 = 0.0_real64
  do i = 1, xdim
    x1 = x1 + abs(coordinates(i))
  end do

  if (x1 <= 1.0e-10_real64) then
    total_energy = get_coeff_zero(coeff_index)
    return
  end if

  call user_coordinates_to_general_coordinates(general_coordinates_ZXZ, xdim, coord_format, coordinates)
  total_energy = get_total_interaction_energy(coeff_index, general_coordinates_ZXZ)
end subroutine evaluate_LRF


!===========================================================
! Batch evaluator:
!  - energies(j) = E( coords(:,j) ) for j=1..n
!  - MPI_flag, OMP_flag are 0/1 selectors
!    (0,0)=serial, (1,0)=MPI, (0,1)=OpenMP, (1,1)=hybrid
!  - MPI combine: ALLGATHERV so EVERY rank ends up with full energies(1:n)
!  - No MPI calls at all unless MPI_flag==1 (safe for OMP-only runs)
!===========================================================
subroutine evaluate_LRF_batch(energies, xdim, coords, n, coord_format, filename, MPI_flag, OMP_flag, ierr)
  use, intrinsic :: iso_fortran_env, only : int32, real64
  use Fitting_Constant_v2,  only : get_coeff_index, get_coeff_zero, ensure_t_index_map_ready, ensure_fact_nn_ready
  use Geometry_Constant_v2, only : get_total_interaction_energy, LMAX
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  integer(int32),               intent(in)  :: xdim, n
  real(real64),                 intent(in)  :: coords(xdim, n)
  character(*),                 intent(in)  :: coord_format, filename
  integer,                      intent(in)  :: MPI_flag, OMP_flag
  real(real64),                 intent(out) :: energies(n)
  integer,            optional, intent(out) :: ierr

  integer :: myerr
  integer(int32) :: j, j1, j2, base, extra
  integer(int32) :: coeff_index
  real(real64)   :: e, x1,zero_e
  real(real64)   :: gen_zxz(6)
  logical :: want_mpi, want_omp
#ifdef USE_MPI
  logical :: mpi_active
  integer :: rank, nproc, mpierr, r
  integer, allocatable :: counts(:), displs(:)
#else
  integer, parameter :: rank  = 0
  integer, parameter :: nproc = 1
#endif

  myerr     = 0
  energies  = 0.0_real64
  want_mpi  = (MPI_flag /= 0)
  want_omp  = (OMP_flag /= 0)

#ifdef USE_MPI
  mpi_active = .false.
  if (want_mpi) then
    call MPI_Initialized(mpi_active, mpierr)
    if (mpi_active) then
      call MPI_Comm_rank(MPI_COMM_WORLD, rank,  mpierr)
      call MPI_Comm_size(MPI_COMM_WORLD, nproc, mpierr)
    else
      ! requested MPI but MPI is not initialized by the caller:
      ! fall back to non-MPI path (no error; just do serial/OMP)
      want_mpi = .false.
    end if
  end if
#endif

  ! Partition
  if (want_mpi) then
    base  = n / nproc
    extra = mod(n, nproc)
    if (rank < extra) then
      j1 = rank*(base+1) + 1
      j2 = j1 + base
    else
      j1 = rank*base + extra + 1
      j2 = j1 + base - 1
    end if
  else
    j1 = 1_int32
    j2 = n
  end if
  if (j2 < j1) then
    j1 = 1_int32
    j2 = 0_int32
  end if

  ! Load coeffs (once per process)
  coeff_index = get_coeff_index(filename)
  if (coeff_index < 1) then
    myerr = -20
    if (present(ierr)) ierr = myerr
    return
  end if

  call ensure_t_index_map_ready(LMAX)
  call ensure_fact_nn_ready(LMAX)

  zero_e = get_coeff_zero(coeff_index)

!$omp parallel do if (want_omp) default(none) &
!$omp& shared(coords, energies, xdim, j1, j2, coeff_index, coord_format) &
!$omp& private(j, x1, gen_zxz, e)
  do j = j1, j2
    x1 = 0.0_real64
    if (xdim >= 1) x1 = sum(abs(coords(1:xdim, j)))
    if (x1 <= 1.0e-10_real64) then
      e = get_coeff_zero(coeff_index)
    else
      call user_coordinates_to_general_coordinates(gen_zxz, xdim, coord_format, coords(:, j))
      e = get_total_interaction_energy(coeff_index, gen_zxz)
    end if
    energies(j) = e
  end do
!$omp end parallel do

#ifdef USE_MPI
  if (want_mpi) then
    base  = n / nproc
    extra = mod(n, nproc)
    allocate(counts(nproc), displs(nproc))
    do r = 0, nproc-1
      if (r < extra) then
        counts(r+1) = base + 1
        displs(r+1) = r * (base + 1)
      else
        counts(r+1) = base
        displs(r+1) = extra * (base + 1) + (r - extra) * base
      end if
    end do
    ! In-place ALLGATHERV: each rank has its local block already in-place
    call MPI_Allgatherv(MPI_IN_PLACE, 0, MPI_DATATYPE_NULL, &
                        energies, counts, displs, MPI_DOUBLE_PRECISION, &
                        MPI_COMM_WORLD, mpierr)
    deallocate(counts, displs)
  end if
#endif

  if (present(ierr)) ierr = myerr
end subroutine evaluate_LRF_batch






