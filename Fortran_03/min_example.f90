! version: LRF_Fortran v4
! *********************************************************************************************************************
! Minimal example of how to call evaluate_LR for a given set of coefficients exported by: LRF MATLAB v4.x
! evaluate_LR needs as function parameters:
!  Energy (Output, Real): Interaction energy between the monomers (in cm^-1).
!  XDIM (Input, Integer): Number of degrees of freedom in the system.
!  COORDINATE_FORMAT (Input, Char): Euler convention used to describe the monomers’ orientation.
!  coordinates(XDIM) (Input, Real): Intermolecular distance (Å) followed by the angles (in degrees) describing the
!                                   orientation.
!  PATH_TO_COEFFICIENTS (Input, Char): Path to the coefficients file containing the long-range coefficient expansion.
! *********************************************************************************************************************

program min_example
  use iso_fortran_env, only : real64, int32
  implicit none

  integer(int32),       parameter :: XDIM = 4_int32       ! Coordinates Dimensions
  character(len=*),     parameter :: COORDINATE_FORMAT = "Euler_ZYZ"   ! Coordinate Format
  character(len=*),     parameter :: PATH_TO_COEFFICIENTS = "../testing_datafiles/coefficients/D_inf_h(1)_Spherical(1)_Coeff.txt"

  real(real64)                      :: energy    ! Interaction Energy
  real(real64), dimension(XDIM)     :: coordinates
  real(real64), parameter           :: PI = acos(-1.0_real64)

  coordinates = [ 9.224922190454659_real64,                                 & ! R
                  acos(-0.516833742198944_real64) * 180.0_real64 / PI,      & ! beta1
                  acos( 0.761164535894394_real64) * 180.0_real64 / PI,      & ! beta2
                  0.081548803182827_real64 * 180.0_real64 / PI ]              ! alpha
  ! For XDIM=6 you could extend with gamma1, gamma2 terms.

  ! Evaluate the Potential Energy Surface in the Long-Range region
  call evaluate_LRF( energy,                 &
                     XDIM,                   &
                     coordinates,            &
                     COORDINATE_FORMAT,      &
                     PATH_TO_COEFFICIENTS )

  ! printing the Interaction Energy in the console
  write(*,*) "Interaction Energy : ", energy, " (cm^-1)"
end program min_example