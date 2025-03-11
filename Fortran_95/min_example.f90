! version: LRF_Fortran v4
! *********************************************************************************************************************
! Minimal example of how to call evaluate_LR for a given set of coefficients exported by: LRF MATLAB  v4.x
! evaluate_LR needs as function parameters:
!  Energy (Output, Real): Interaction energy between the monomers.
!  xdim (Input, Integer): Number of degrees of freedom in the system.
!  coordinate_format (Input, Char): Euler convention used to describe the monomers’ orientation.
!  coordinates(xdim) (Input, Real): Intermolecular distance (Å) followed by the angles (in degrees) describing the
!                                   orientation.
!  path_to_coefficients (Input, Char): Path to the coefficients file containing the long-range coefficient expansion.
! *********************************************************************************************************************

program min_example

 implicit none

 integer (kind=4), parameter:: xdim=4  ! Coordinates Dimensions
 character (len=*,kind=1), parameter:: coordinate_format = "Euler_ZYZ"   ! Coordinate Format
 character (len=*,kind=1), parameter:: path_to_coefficients = "../testing_datafiles/coefficients/C1(1)_C1(1)_Coeff.txt"

 real (kind=8) :: Energy    ! Interaction Energy
 real (kind=8), dimension(xdim) :: coordinates = [10.27d0,& ! R
                                                  30d0,&                 ! beta1
                                                  20d0,&                 ! beta2
                                                  120d0]                 ! alpha

 ! Evaluate the Potential Energy Surface in the Long-Range region

 call evaluate_LRF( Energy,&        ! Interaction Energy
                    xdim,&          ! Coordinates Dimensions
                    coordinates,&   ! System of Coordinate
                    coordinate_format,&   ! Coordinate Format
                    path_to_coefficients& ! path to coefficients file
                  )

 !printing Output Energy in the console
 write(*,*) "Interaction Energy : ", Energy, " (cm^-1)"

end program min_example










