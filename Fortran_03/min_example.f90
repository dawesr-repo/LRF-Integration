! ======================================================================================
! Minimal usage example for evaluate_LRF (parallelism is handled inside Geometry module)
! ======================================================================================
program min_example
  use iso_fortran_env, only : real64, int32
  use LRF_API,         only : evaluate_LRF
  implicit none

  integer(int32),       parameter :: XDIM = 4_int32
  character(len=*),     parameter :: COORDINATE_FORMAT = "Euler_ZYZ"
  character(len=*),     parameter :: PATH_TO_COEFFICIENTS = "../testing_datafiles/coefficients/D_inf_h(1)_Spherical(1)_Coeff.txt"

  real(real64)                      :: energy
  real(real64), dimension(XDIM)     :: coordinates
  real(real64), parameter           :: PI = acos(-1.0_real64)

  coordinates = [ 9.224922190454659_real64,                                 & ! R
                  acos(-0.516833742198944_real64) * 180.0_real64 / PI,      & ! beta1
                  acos( 0.761164535894394_real64) * 180.0_real64 / PI,      & ! beta2
                  0.081548803182827_real64 * 180.0_real64 / PI ]              ! alpha

  call evaluate_LRF( energy,                 &
                     XDIM,                   &
                     coordinates,            &
                     COORDINATE_FORMAT,      &
                     PATH_TO_COEFFICIENTS )

  write(*,*) "Interaction Energy : ", energy, " (cm^-1)"
end program min_example
