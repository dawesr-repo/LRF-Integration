! version: LRF_Fortran v4
! *********************************************************************************************************************
! Minimal example of how to call evaluate_LR for a given set of coefficients exported by: LRF MATLAB  v4.x
! evaluate_LR needs as function parameters:
!  Energy (Output, Real): Interaction energy between the monomers (in cm^-1).
!  XDIM (Input, Integer): Number of degrees of freedom in the system.
!  COORDINATE_FORMAT (Input, Char): Euler convention used to describe the monomers’ orientation.
!  coordinates(XDIM) (Input, Real): Intermolecular distance (Å) followed by the angles (in degrees) describing the
!                                   orientation.
!  PATH_TO_COEFFICIENTS (Input, Char): Path to the coefficients file containing the long-range coefficient expansion.
! *********************************************************************************************************************

program omp_hello
  use omp_lib
  implicit none
!$omp parallel
  print *, "Hello from thread"
!$omp end parallel
end program










