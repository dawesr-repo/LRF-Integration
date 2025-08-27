!===============================================================
! test_v2.f90  — main test driver that your Makefile runs
! - Positional args: THREADS NPROC MPI_flag OMP_flag
! - Safe MPI init/finalize only when MPI_flag==1
! - Prints once (root) because helpers gate output
!===============================================================
program test_v2
  use iso_fortran_env, only : int32, real64
#ifdef _OPENMP
  use omp_lib
#endif
#ifdef USE_MPI
  use mpi
#endif
  implicit none
  interface
    subroutine running_time_performance(coeff_file_name, fileoutput_number)
      use iso_fortran_env, only : int32, real64
      character(len=*), intent(in) :: coeff_file_name
      integer(int32),   optional   :: fileoutput_number
    end subroutine running_time_performance
  end interface

  integer(int32), parameter :: XDIM = 4_int32
  character(len=*), parameter :: COORDINATE_FORMAT = "Euler_ZYZ"
  character(len=*), parameter :: PATH_TO_COEFFICIENTS = &
       "../testing_datafiles/coefficients/D_inf_h(1)_Spherical(1)_Coeff.txt"

  integer(int32) :: nthreads = 0_int32
  integer(int32) :: nproc_in = 1_int32
  integer        :: MPI_flag = 0
  integer        :: OMP_flag = 0
  character(len=64) :: s
  integer :: narg, ios

  ! MPI control
  logical :: want_mpi, we_inited_mpi
#ifdef USE_MPI
  logical :: have_mpi
  integer :: myrank=0, nproc=1, mpierr, provided
#else
  integer, parameter :: myrank = 0, nproc = 1
#endif

  ! ---- parse positional CLI args: THREADS NPROC MPI_flag OMP_flag ----
  narg = command_argument_count()
  if (narg >= 1) then
    call get_command_argument(1, s); read(s,*,iostat=ios) nthreads; if (ios /= 0) nthreads = 0
  end if
  if (narg >= 2) then
    call get_command_argument(2, s); read(s,*,iostat=ios) nproc_in; if (ios /= 0) nproc_in = 1
  end if
  if (narg >= 3) then
    call get_command_argument(3, s); read(s,*,iostat=ios) MPI_flag; if (ios /= 0) MPI_flag = 0
  end if
  if (narg >= 4) then
    call get_command_argument(4, s); read(s,*,iostat=ios) OMP_flag; if (ios /= 0) OMP_flag = 0
  end if

  want_mpi      = (MPI_flag /= 0)
  we_inited_mpi = .false.

#ifdef _OPENMP
  if (OMP_flag /= 0 .and. nthreads > 0) call omp_set_num_threads(nthreads)
#endif

#ifdef USE_MPI
  if (want_mpi) then
    call MPI_Initialized(have_mpi, mpierr)
    if (.not. have_mpi) then
      call MPI_Init_thread(MPI_THREAD_FUNNELED, provided, mpierr)
      we_inited_mpi = .true.
    end if
    call MPI_Comm_rank(MPI_COMM_WORLD, myrank, mpierr)
    call MPI_Comm_size(MPI_COMM_WORLD, nproc, mpierr)
  end if
#endif

  ! ---- Run a small demo that prints a 1xN energy vector (root only) ----
  call run_batch_demo(PATH_TO_COEFFICIENTS, COORDINATE_FORMAT, MPI_flag, OMP_flag, nthreads)

  ! ---- Run the timing suite; prints a summary once (root only) ----
  call running_time_performance(PATH_TO_COEFFICIENTS)

#ifdef USE_MPI
  if (want_mpi .and. we_inited_mpi) call MPI_Finalize(mpierr)
#endif
end program test_v2
