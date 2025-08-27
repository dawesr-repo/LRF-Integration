! ======================================================================================
! Minimal usage example for evaluate_LRF_batch
! (no module; evaluate_LRF_batch is an external subroutine)
! ======================================================================================
program min_example_batch
  use iso_fortran_env, only : real64, int32
  implicit none

  integer(int32),       parameter :: XDIM = 4_int32
  integer(int32),       parameter :: N    = 5_int32
  character(len=*),     parameter :: COORDINATE_FORMAT = "Euler_ZYZ"
  character(len=*),     parameter :: PATH_TO_COEFFICIENTS = &
       "../testing_datafiles/coefficients/D_inf_h(1)_Spherical(1)_Coeff.txt"

  real(real64),                  parameter :: PI = dacos(-1.0_real64)
  real(real64)                             :: energies(N)
  real(real64)                             :: coords(XDIM, N)
  integer                                   :: j

  ! Parallel setup
  integer(int32) :: nthreads = 8_int32
  integer(int32) :: nproc_in = 2_int32
  integer        :: MPI_flag = 1 ! ---- defaults (serial) ----
  integer        :: OMP_flag = 1 ! ---- defaults (serial) ----

  interface
    subroutine evaluate_LRF_batch(energies, xdim, coords, n, coord_format, filename, MPI_flag, OMP_flag, ierr)
        use, intrinsic :: iso_fortran_env, only : int32, real64
        implicit none
        integer(int32),               intent(in)  :: xdim, n
        real(real64),                 intent(out) :: energies(n)
        real(real64),                 intent(in)  :: coords(xdim, n)
        character(*),                 intent(in)  :: coord_format, filename
        integer,                      intent(in)  :: MPI_flag, OMP_flag
        integer,            optional, intent(out) :: ierr
    end subroutine evaluate_LRF_batch
   end interface  

  ! build a small batch by varying alpha (degrees)
  block
    real(real64) :: R, beta1, beta2, alpha0
    integer      :: k
    R      = 9.224922190454659_real64
    beta1  = dacos(-0.516833742198944_real64) * 180.0_real64 / PI
    beta2  = dacos( 0.761164535894394_real64) * 180.0_real64 / PI
    alpha0 = 0.081548803182827_real64 * 180.0_real64 / PI
    do k = 1, N
      coords(:, k) = [ R, beta1, beta2, alpha0 + 30.0_real64*real(k-1, real64) ]
    end do
  end block



call evaluate_LRF_batch(    energies,              &
                            XDIM,                  &
                            coords, N,             &
                            COORDINATE_FORMAT,     &
                            PATH_TO_COEFFICIENTS,  &
                            MPI_flag, OMP_flag )



write(*,*) "Interaction Energy : ", energies, " (cm^-1)"


end program min_example_batch
