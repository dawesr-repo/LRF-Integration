! ======================================================================================
! Minimal usage example for evaluate_LRF (no module; evaluate_LRF is an external subroutine)
! ======================================================================================
program min_example
  use iso_fortran_env, only : real64, int32
  implicit none

  integer(int32),       parameter :: XDIM = 4_int32
  character(len=*),     parameter :: COORDINATE_FORMAT = "Euler_ZYZ"
  character(len=*),     parameter :: PATH_TO_COEFFICIENTS = "../testing_datafiles/coefficients/D_inf_h(1)_Spherical(1)_Coeff.txt"

  real(real64)                  :: energy
  real(real64)                  :: coordinates(XDIM)
  real(real64), parameter       :: PI = dacos(-1.0_real64)

  ! declare the external subroutine (implemented in your helper_functions.f90)
  external :: evaluate_LRF

  coordinates = [ 9.224922190454659_real64,                                 & ! R
                  dacos(-0.516833742198944_real64) * 180.0_real64 / PI,      & ! beta1
                  dacos( 0.761164535894394_real64) * 180.0_real64 / PI,      & ! beta2
                  0.081548803182827_real64 * 180.0_real64 / PI ]              ! alpha

  call evaluate_LRF( energy,                 &
                     XDIM,                   &
                     coordinates,            &
                     COORDINATE_FORMAT,      &
                     PATH_TO_COEFFICIENTS )

  write(*,*) "Interaction Energy : ", energy, " (cm^-1)"
end program min_example
