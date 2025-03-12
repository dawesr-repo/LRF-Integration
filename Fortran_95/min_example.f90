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

program min_example

 implicit none

 integer (kind=4), parameter:: XDIM=6  ! Coordinates Dimensions
 character (len=*,kind=1), parameter:: COORDINATE_FORMAT = "Euler_ZYZ"   ! Coordinate Format
 character (len=*,kind=1), parameter:: PATH_TO_COEFFICIENTS = "../testing_datafiles/coefficients/C1(1)_C1(1)_Coeff.txt"

 real (kind=8):: energy    ! Interaction Energy
 real (kind=8), dimension(XDIM):: coordinates = [10.27d0,& ! R
                                                  30d0,&    ! beta1
                                                  20d0,&    ! beta2
                                                  120d0,&   ! alpha
                                                  0d0,&     ! gamma1
                                                  0d0]     ! gamma2


 ! Evaluate the Potential Energy Surface in the Long-Range region
 call evaluate_LRF( energy,&
                    XDIM,&
                    coordinates,&
                    COORDINATE_FORMAT,&
                    PATH_TO_COEFFICIENTS&
                  )

 ! printing the Interaction Energy in the console
 write(*,*) "Interaction Energy : ", Energy, " (cm^-1)"

end program min_example










