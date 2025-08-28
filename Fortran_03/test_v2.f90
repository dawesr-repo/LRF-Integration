program test_v2
  use, intrinsic :: iso_fortran_env, only : int32, real64
  use testing_v2_mod, only : evaluate_LRF_batch, running_time_performance
#ifdef USE_MPI
  use mpi
#endif
#ifdef _OPENMP
  use omp_lib
#endif
  implicit none

  integer(int32) :: THREADS, NPROC, MPI_flag, OMP_flag, NTEST

  ! ---------------- MPI state (distinct names; avoid clash with NPROC) -----------
  integer :: rank_mpi, nproc_mpi
  logical :: we_started_mpi
#ifdef USE_MPI
  logical :: mpi_active
  integer :: mpierr, provided
#endif

  ! ---------------- Problem constants ----------------
  integer(int32),       parameter :: XDIM = 4_int32
  integer(int32),       parameter :: N    = 5_int32
  character(len=*),     parameter :: COORDINATE_FORMAT   = "Euler_ZYZ"
  character(len=*),     parameter :: PATH_TO_COEFFICIENTS = &
       "../testing_datafiles/coefficients/D_inf_h(1)_Spherical(1)_Coeff.txt"

  ! ---------------- Data ----------------
  real(real64) :: energies(N)
  real(real64) :: coords(XDIM, N)

  ! ====================== MAIN ======================
  call parse_args(THREADS, NPROC, MPI_flag, OMP_flag, NTEST)

  ! Build a tiny batch (vary alpha by +30° each)
  call make_demo_batch(coords)

  we_started_mpi = .false.
  rank_mpi  = 0
  nproc_mpi = 1

#ifdef _OPENMP
  if (THREADS > 0) call omp_set_num_threads(THREADS)
#endif

#ifdef USE_MPI
  if (MPI_flag /= 0) then
    call MPI_Initialized(mpi_active, mpierr)
    if (.not. mpi_active) then
      call MPI_Init_thread(MPI_THREAD_FUNNELED, provided, mpierr)
      we_started_mpi = .true.
    end if
    call MPI_Comm_rank(MPI_COMM_WORLD, rank_mpi,  mpierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc_mpi, mpierr)
  end if
#endif

  ! ---- Batch evaluate (library is MPI-agnostic; only uses MPI if already inited) ----
  call evaluate_LRF_batch(energies, XDIM, coords, N, COORDINATE_FORMAT, PATH_TO_COEFFICIENTS, MPI_flag, OMP_flag)

  if (rank_mpi == 0) then
    write(*,*) "Interaction Energy (N=", N, "): ", energies, " (cm^-1)"
  end if

  ! ---- Performance sweeps (always print once from root) ----
  if (rank_mpi == 0) then
    call running_time_performance(PATH_TO_COEFFICIENTS, NTEST, run_mpi=(MPI_flag /= 0), run_omp=(OMP_flag /= 0))
  end if

#ifdef USE_MPI
  if (MPI_flag /= 0) then
    call MPI_Barrier(MPI_COMM_WORLD, mpierr)
    if (we_started_mpi) call MPI_Finalize(mpierr)
  end if
#endif

contains

  subroutine parse_args(threads, nproc_cli, mpi_f, omp_f, ntest)
    use, intrinsic :: iso_fortran_env, only : int32
    implicit none
    integer(int32), intent(out) :: threads, nproc_cli, mpi_f, omp_f, ntest
    integer :: argc
    character(len=64) :: s

    threads   = 0_int32
    nproc_cli = 1_int32
    mpi_f     = 0_int32
    omp_f     = 0_int32
    ntest     = 1000_int32    ! default if not provided

    argc = command_argument_count()
    if (argc >= 1) then
      call get_command_argument(1, s); read(s, *, err=10) threads
    end if
    if (argc >= 2) then
      call get_command_argument(2, s); read(s, *, err=10) nproc_cli
    end if
    if (argc >= 3) then
      call get_command_argument(3, s); read(s, *, err=10) mpi_f
    end if
    if (argc >= 4) then
      call get_command_argument(4, s); read(s, *, err=10) omp_f
    end if
    if (argc >= 5) then
      call get_command_argument(5, s); read(s, *, err=10) ntest
    end if
10  continue
  end subroutine parse_args

  subroutine make_demo_batch(coords)
    use, intrinsic :: iso_fortran_env, only : real64, int32
    implicit none
    real(real64), intent(out) :: coords(XDIM, N)
    real(real64), parameter   :: PI = dacos(-1.0_real64)
    real(real64) :: R, beta1, beta2, alpha0
    integer(int32) :: k

    R      = 9.224922190454659_real64
    beta1  = dacos(-0.516833742198944_real64) * 180.0_real64 / PI
    beta2  = dacos( 0.761164535894394_real64) * 180.0_real64 / PI
    alpha0 = 0.081548803182827_real64        * 180.0_real64 / PI

    do k = 1, N
      coords(:, k) = [ R, beta1, beta2, alpha0 + 30.0_real64*real(k-1, real64) ]
    end do
  end subroutine make_demo_batch

end program test_v2
