module testing_v2_mod
  use, intrinsic :: iso_fortran_env, only : int32, real64
  implicit none
  private

  ! Public API
  public :: am_root_env
  public :: run_batch_demo
  public :: running_time_performance
  public :: serial_6d_timing
  public :: omp_6d_timing
  public :: mpi_6d_timing
  public :: hybrid_6d_timing
  public :: evaluate_LRF, evaluate_LRF_batch

  ! Explicit interfaces for external procedures (implemented elsewhere)
  interface
    subroutine evaluate_LRF(total_energy, xdim, coordinates, coord_format, filename)
      import :: real64, int32
      real(real64),   intent(out) :: total_energy
      integer(int32), intent(in)  :: xdim
      real(real64),   intent(in)  :: coordinates(xdim)
      character(*),   intent(in)  :: coord_format
      character(*),   intent(in)  :: filename
    end subroutine evaluate_LRF

    subroutine evaluate_LRF_batch(energies, xdim, coords, n, coord_format, filename, MPI_flag, OMP_flag, ierr)
      import :: real64, int32
      integer(int32),               intent(in)  :: xdim, n
      real(real64),                 intent(out) :: energies(n)
      real(real64),                 intent(in)  :: coords(xdim, n)
      character(*),                 intent(in)  :: coord_format, filename
      integer,                      intent(in)  :: MPI_flag, OMP_flag
      integer,            optional, intent(out) :: ierr
    end subroutine evaluate_LRF_batch
  end interface

contains

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

  subroutine run_batch_demo(coeff_file_name, coord_format, MPI_flag, OMP_flag, nthreads)
    use, intrinsic :: iso_fortran_env, only : int32, real64
#ifdef _OPENMP
    use omp_lib, only : omp_set_num_threads
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

    call evaluate_LRF_batch(energies, XDIM, coords, N, coord_format, coeff_file_name, MPI_flag, OMP_flag)

    if (am_root_env()) then
      write(*,*) 'Interaction Energy (N=', size(energies), '): ', energies, ' (cm^-1)'
    end if
  end subroutine run_batch_demo

  subroutine running_time_performance(coeff_file_name, ntest, fileoutput_number, run_mpi, run_omp)
    use, intrinsic :: iso_fortran_env, only: real64, int32
    implicit none
    character(len=*), intent(in)           :: coeff_file_name
    integer(int32),   intent(in), optional :: ntest
    integer(int32),   intent(in), optional :: fileoutput_number
    logical,          intent(in), optional :: run_mpi, run_omp

    integer(int32) :: nrep
    logical        :: do_mpi, do_omp

    nrep   = merge(ntest, 1000_int32, present(ntest))
    do_mpi = merge(run_mpi, .false.,    present(run_mpi))
    do_omp = merge(run_omp, .false.,    present(run_omp))

    call serial_6d_timing(coeff_file_name, nrep)

    if (do_omp .and. .not. do_mpi) then
       call omp_6d_timing(coeff_file_name, nrep)
    else
       if (am_root_env()) then
         write(*,*)"*********************************************************************"
         write(*,*)" * OpenMP PERFORMANCE -- OpenMP not set *"
         write(*,*)"*********************************************************************"
       end if
    end if

    if (do_mpi .and. .not. do_omp) then
       call mpi_6d_timing(coeff_file_name, nrep)
    else
       if (am_root_env()) then
         write(*,*)"*********************************************************************"
         write(*,*)" * MPI PERFORMANCE -- MPI not set *"
         write(*,*)"*********************************************************************"
       end if
    end if

    if (do_mpi .and. do_omp) then
       call hybrid_6d_timing(coeff_file_name, nrep)
    else
       if (am_root_env()) then
         write(*,*)"*********************************************************************"
         write(*,*)" * HYBRID PERFORMANCE -- MPI and/or OpenMP not set *"
         write(*,*)"*********************************************************************"
       end if
    end if
  end subroutine running_time_performance

  subroutine serial_6d_timing(coeff_file_name, nrep)
    use, intrinsic :: iso_fortran_env, only : real64, int32
    implicit none
    character(len=*), intent(in) :: coeff_file_name
    integer(int32),   intent(in) :: nrep

    character(len=*), parameter :: COORD_FORMAT = "Euler_ZYZ"
    integer(int32),   parameter :: XDIM = 6_int32

    real(real64) :: coordinates(XDIM)
    real(real64) :: energy, t0, t1
    integer(int32) :: i

    coordinates = [ 10.27_real64, 30.0_real64, 20.0_real64, 120.0_real64, 40.0_real64, 50.0_real64 ]

    call cpu_time(t0)
    do i = 1, nrep
      call evaluate_LRF(energy, XDIM, coordinates, COORD_FORMAT, coeff_file_name)
    end do
    call cpu_time(t1)

    if (am_root_env()) then
      write(*,*)"*********************************************************************"
      write(*,*) "* SERIAL PERFORMANCE for 6D (N=", nrep, "): time = ", t1-t0, " s *"
      write(*,*)"*********************************************************************"
    end if
  end subroutine serial_6d_timing

  subroutine omp_6d_timing(coeff_file_name, nrep)
    use, intrinsic :: iso_fortran_env, only : real64, int32
#ifdef _OPENMP
    use omp_lib
#endif
    implicit none
    character(len=*), intent(in) :: coeff_file_name
    integer(int32),   intent(in) :: nrep

    character(len=*), parameter :: COORD_FORMAT = "Euler_ZYZ"
    integer(int32),   parameter :: XDIM = 6_int32
    integer(int32) :: N
    real(real64), allocatable :: coords(:,:), energies(:)
    real(real64) :: t0, t1
    integer(int32) :: k

    N = nrep
    allocate(coords(XDIM,N), energies(N))

    do k = 1, N
      coords(:,k) = [ 10.27_real64, 30.0_real64, 20.0_real64, &
                      120.0_real64 + 0.05_real64*real(k-1,real64), 40.0_real64, 50.0_real64 ]
    end do

    call cpu_time(t0)
    call evaluate_LRF_batch(energies, XDIM, coords, N, COORD_FORMAT, coeff_file_name, &
                            MPI_flag=0, OMP_flag=1)
    call cpu_time(t1)

    if (am_root_env()) then
      write(*,*)"*********************************************************************"
      write(*,*) "* OpenMP PERFORMANCE  (N=", N, "): time = ", t1-t0, " s *"
      write(*,*)"*********************************************************************"
    end if

    deallocate(coords, energies)
  end subroutine omp_6d_timing

  subroutine mpi_6d_timing(coeff_file_name, nrep)
    use, intrinsic :: iso_fortran_env, only : real64, int32
    implicit none
    character(len=*), intent(in) :: coeff_file_name
    integer(int32),   intent(in) :: nrep

    character(len=*), parameter :: COORD_FORMAT = "Euler_ZYZ"
    integer(int32),   parameter :: XDIM = 6_int32
    integer(int32) :: N
    real(real64), allocatable :: coords(:,:), energies(:)
    real(real64) :: t0, t1
    integer(int32) :: k

    N = nrep
    allocate(coords(XDIM,N), energies(N))

    do k = 1, N
      coords(:,k) = [ 10.27_real64, 30.0_real64, 20.0_real64, &
                      120.0_real64 + 0.05_real64*real(k-1,real64), 40.0_real64, 50.0_real64 ]
    end do

    call cpu_time(t0)
    call evaluate_LRF_batch(energies, XDIM, coords, N, COORD_FORMAT, coeff_file_name, &
                            MPI_flag=1, OMP_flag=0)
    call cpu_time(t1)

    if (am_root_env()) then
      write(*,*)"*********************************************************************"
      write(*,*) "* MPI PERFORMANCE      (N=", N, "): time = ", t1-t0, " s *"
      write(*,*)"*********************************************************************"
    end if

    deallocate(coords, energies)
  end subroutine mpi_6d_timing

  subroutine hybrid_6d_timing(coeff_file_name, nrep)
    use, intrinsic :: iso_fortran_env, only : real64, int32
    implicit none
    character(len=*), intent(in) :: coeff_file_name
    integer(int32),   intent(in) :: nrep

    character(len=*), parameter :: COORD_FORMAT = "Euler_ZYZ"
    integer(int32),   parameter :: XDIM = 6_int32
    integer(int32) :: N
    real(real64), allocatable :: coords(:,:), energies(:)
    real(real64) :: t0, t1
    integer(int32) :: k

    N = nrep
    allocate(coords(XDIM,N), energies(N))

    do k = 1, N
      coords(:,k) = [ 10.27_real64, 30.0_real64, 20.0_real64, &
                      120.0_real64 + 0.05_real64*real(k-1,real64), 40.0_real64, 50.0_real64 ]
    end do

    call cpu_time(t0)
    call evaluate_LRF_batch(energies, XDIM, coords, N, COORD_FORMAT, coeff_file_name, &
                            MPI_flag=1, OMP_flag=1)
    call cpu_time(t1)

    if (am_root_env()) then
      write(*,*)"*********************************************************************"
      write(*,*) "* HYBRID PERFORMANCE   (N=", N, "): time = ", t1-t0, " s *"
      write(*,*)"*********************************************************************"
    end if

    deallocate(coords, energies)
  end subroutine hybrid_6d_timing

end module testing_v2_mod
