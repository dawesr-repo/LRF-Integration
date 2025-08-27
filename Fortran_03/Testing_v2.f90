!===============================================================
! Testing_v2.f90  — helper routines for the test driver
! - Prints only once (root rank) even under mpirun
! - Never calls MPI unless MPI_flag==1 and code was built with USE_MPI
!===============================================================
subroutine run_batch_demo(coeff_file_name, coord_format, MPI_flag, OMP_flag, nthreads)
  use iso_fortran_env, only : int32, real64
#ifdef _OPENMP
  use omp_lib
#endif
  implicit none
  character(len=*), intent(in) :: coeff_file_name, coord_format
  integer,          intent(in) :: MPI_flag, OMP_flag, nthreads

  integer(int32), parameter :: XDIM = 4_int32
  integer(int32), parameter :: N    = 5_int32
  real(real64),   parameter :: PI   = dacos(-1.0_real64)

  real(real64) :: energies(N)
  real(real64) :: coords(XDIM, N)
  integer      :: j

  logical, external :: am_root_env

  ! ---- explicit interface for batch evaluator (callee has optional ierr) ----
  interface
    subroutine evaluate_LRF_batch(energies, xdim, coords, n, coord_format, filename, MPI_flag, OMP_flag, ierr)
      use iso_fortran_env, only : int32, real64
      integer(int32),               intent(in)  :: xdim, n
      real(real64),                 intent(out) :: energies(n)
      real(real64),                 intent(in)  :: coords(xdim, n)
      character(*),                 intent(in)  :: coord_format, filename
      integer,                      intent(in)  :: MPI_flag, OMP_flag
      integer,            optional, intent(out) :: ierr
    end subroutine
  end interface

  ! ---- small batch by varying alpha (deg) ----
  block
    real(real64) :: R, beta1, beta2, alpha0
    R      = 9.224922190454659_real64
    beta1  = dacos(-0.516833742198944_real64) * 180.0_real64 / PI
    beta2  = dacos( 0.761164535894394_real64) * 180.0_real64 / PI
    alpha0 = 0.081548803182827_real64 * 180.0_real64 / PI
    do j = 1, N
      coords(:, j) = [ R, beta1, beta2, alpha0 + 30.0_real64*real(j-1, real64) ]
    end do
  end block

#ifdef _OPENMP
  if (OMP_flag /= 0 .and. nthreads > 0) call omp_set_num_threads(nthreads)
#endif

  call evaluate_LRF_batch( energies, XDIM, coords, N, coord_format, coeff_file_name, &
                           MPI_flag, OMP_flag )

  if (am_root_env()) then
    write(*,*) 'Interaction Energy (N=', size(energies), '): ', energies, ' (cm^-1)'
  end if
end subroutine run_batch_demo


subroutine running_time_performance(coeff_file_name, fileoutput_number)
  use iso_fortran_env, only : int32, real64
#ifdef _OPENMP
  use omp_lib
#endif
  implicit none
  character(len=*), intent(in) :: coeff_file_name
  integer(int32),   optional   :: fileoutput_number

  character(len=*), parameter :: COORD_FORMAT = "Euler_ZYZ"
  integer(int32),   parameter :: xdim  = 6_int32
  integer(int32),   parameter :: ntest = 1000_int32
  real(real64),     parameter :: coord6(6) = [ 10.27_real64, 30.0_real64, 20.0_real64, &
                                               120.0_real64, 40.0_real64, 50.0_real64 ]

  real(real64) :: energy, t0, t1, t_serial, t_omp, t_mpi, t_hyb
  real(real64), allocatable :: coords(:, :), energies(:)
  integer(int32) :: i
  logical :: have_omp, have_mpi, iam_root
  logical, external :: am_root_env

  ! explicit interfaces for evaluators
  interface
    subroutine evaluate_LRF(energy, xdim, coords, coord_format, filename)
      use iso_fortran_env, only : int32, real64
      real(real64),                 intent(out) :: energy
      integer(int32),               intent(in)  :: xdim
      real(real64),                 intent(in)  :: coords(xdim)
      character(*),                 intent(in)  :: coord_format, filename
    end subroutine
    subroutine evaluate_LRF_batch(energies, xdim, coords, n, coord_format, filename, MPI_flag, OMP_flag, ierr)
      use iso_fortran_env, only : int32, real64
      integer(int32),               intent(in)  :: xdim, n
      real(real64),                 intent(out) :: energies(n)
      real(real64),                 intent(in)  :: coords(xdim, n)
      character(*),                 intent(in)  :: coord_format, filename
      integer,                      intent(in)  :: MPI_flag, OMP_flag
      integer,            optional, intent(out) :: ierr
    end subroutine
  end interface

  call wall_time(t0)
  do i = 1, ntest
    call evaluate_LRF( energy, xdim, coord6, COORD_FORMAT, coeff_file_name )
  end do
  call wall_time(t1)
  t_serial = t1 - t0

  have_omp = .false.
#ifdef _OPENMP
  have_omp = .true.
#endif

  have_mpi = .false.
#ifdef USE_MPI
  have_mpi = .true.
#endif

  t_omp = -1.0_real64
  if (have_omp) then
    allocate(coords(xdim, ntest), energies(ntest))
    coords = spread(coord6, dim=2, ncopies=ntest)
    call wall_time(t0)
    call evaluate_LRF_batch(energies, xdim, coords, ntest, COORD_FORMAT, coeff_file_name, 0, 1)
    call wall_time(t1)
    t_omp = t1 - t0
    deallocate(coords, energies)
  end if

  t_mpi = -1.0_real64
  if (have_mpi) then
    allocate(coords(xdim, ntest), energies(ntest))
    coords = spread(coord6, dim=2, ncopies=ntest)
    call wall_time(t0)
    call evaluate_LRF_batch(energies, xdim, coords, ntest, COORD_FORMAT, coeff_file_name, 1, 0)
    call wall_time(t1)
    t_mpi = t1 - t0
    deallocate(coords, energies)
  end if

  t_hyb = -1.0_real64
  if (have_mpi .and. have_omp) then
    allocate(coords(xdim, ntest), energies(ntest))
    coords = spread(coord6, dim=2, ncopies=ntest)
    call wall_time(t0)
    call evaluate_LRF_batch(energies, xdim, coords, ntest, COORD_FORMAT, coeff_file_name, 1, 1)
    call wall_time(t1)
    t_hyb = t1 - t0
    deallocate(coords, energies)
  end if

  iam_root = am_root_env()
  if (iam_root) then
    write(*,*) "*********************************************************************"
    write(*,'(A,I0,A,F10.4,A)') "* SERIAL PERFORMANCE for 6D (N=", ntest, "): time = ", t_serial, " s *"
    write(*,*) "*********************************************************************"
    if (have_omp) then
      write(*,'(A,I0,A,F10.4,A)') "* OpenMP PERFORMANCE (N=", ntest, "): time = ", t_omp, " s *"
    else
      write(*,*) "* OpenMP PERFORMANCE -- OpenMP not set *"
    end if
    write(*,*) "*********************************************************************"
    if (have_mpi) then
      write(*,'(A,I0,A,F10.4,A)') "* MPI PERFORMANCE     (N=", ntest, "): time = ", t_mpi, " s *"
    else
      write(*,*) "* MPI PERFORMANCE -- MPI not set *"
    end if
    write(*,*) "*********************************************************************"
    if (have_mpi .and. have_omp) then
      write(*,'(A,I0,A,F10.4,A)') "* HYBRID MPI+OpenMP   (N=", ntest, "): time = ", t_hyb, " s *"
    else
      write(*,*) "* HYBRID PERFORMANCE -- MPI and/or OpenMP not set *"
    end if
    write(*,*) "*********************************************************************"
  end if
contains
  subroutine wall_time(t)
    use iso_fortran_env, only : real64
    real(real64), intent(out) :: t
    integer :: c, r, m
    call system_clock(count=c, count_rate=r, count_max=m)
    if (r > 0) then
      t = real(c, real64) / real(r, real64)
    else
      t = 0.0_real64
    end if
  end subroutine wall_time
end subroutine running_time_performance


logical function am_root_env()
  implicit none
  character(len=32) :: s
  integer :: ln, ios, rr
  am_root_env = .true.
  call get_environment_variable("OMPI_COMM_WORLD_RANK", s, length=ln, status=ios)
  if (ios == 0) then
    read(s(1:ln),*,iostat=ios) rr
    if (ios == 0) am_root_env = (rr == 0)
    return
  end if
  call get_environment_variable("PMI_RANK", s, length=ln, status=ios)
  if (ios == 0) then
    read(s(1:ln),*,iostat=ios) rr
    if (ios == 0) am_root_env = (rr == 0)
  end if
end function am_root_env
