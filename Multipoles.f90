
!  Fortran subroutine for 3.1.1 matlab version Optimization in running time
! last test: 26/12/2022


PROGRAM main_subroutine
    IMPLICIT NONE

    real*8 :: E2,E1
    real*8 ,dimension(6):: GeneralCoordenates,coordinates,zero_coordinates
    real*8 :: FourD_coord(4),TwoD_coord(2)
    Character(len = 20) :: coord_format = "Euler_ZYZ"
    INTEGER :: XDim=4,i
    real*8 :: start, finish,pii
    pii=DAcos(-1d0)
    

    if (XDim==3)Then
        coord_format = "Spherical"  
    end if

    FourD_coord(1) = 10.271963668339600d0 !R
    FourD_coord(2) = dacos(0.745608537837982d0)  !b1
    FourD_coord(3) = Dacos(-0.550799772686005d0)  !b2
    FourD_coord(4) = 0.846716722982929*180d0/pii  !Phi

    ! ThreeD_coord(1) = 10d0  !R
    ! ThreeD_coord(2) = 50d0  !b1
    ! ThreeD_coord(3) = 70d0  !c1

    ! TwoD_coord(1) = 10d0  !R
    ! TwoD_coord(2) = 50d0  !b1

    Call General_Coordinates_Format(XDim, FourD_coord, GeneralCoordenates)
    Call Coordinate_Transformation(GeneralCoordenates,coord_format,coordinates)
    CALL Long_Range_Potential(coordinates,E1)

    !call cpu_time(start)

    ! do i=1,100000

    !     Call General_Coordinates_Format(XDim, TwoD_coord, GeneralCoordenates)
    !     Call Coordinate_Transformation(GeneralCoordenates,coord_format,coordinates)
     
    !     !write(*,*)coordinates
    !     CALL Long_Range_Potential(coordinates,E1)
    !     !write(*,*)E1

    ! enddo    
  
    !call cpu_time(finish)
    !print '("Time = ",f6.3," seconds.")',finish-start
    write(*,*)E1

    ! DON't use Coordinate_Transformation function to obtain the zero potential, 
    ! instead pass an array with all of them zero
    zero_coordinates(1) = 0d0  !R
    zero_coordinates(2) = 0d0  !b1
    zero_coordinates(3) = 0d0  !b2
    zero_coordinates(4) = 0d0  !phi
    zero_coordinates(5) = 0d0  !c1
    zero_coordinates(6) = 0d0  !c2


    CALL Long_Range_Potential(zero_coordinates,E2)
  
    write(*,*)E2

    write(*,*)-186812.30755691 - E2,E1,E1-(-186812.30755691 - E2)


    
END PROGRAM main_subroutine








