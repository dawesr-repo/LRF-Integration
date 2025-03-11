SUBROUTINE Coordinate_Transformation(coordinates,coord_format,new_coordinates)
    IMPLICIT NONE
    real*8 ,dimension(6), INTENT(IN):: coordinates

    Character(*), INTENT(IN) :: coord_format
    real*8 , dimension(6), INTENT(INOUT) :: new_coordinates

    new_coordinates=coordinates;


    if (coord_format == "Euler_ZYZ") then
        new_coordinates(5) = coordinates(5) - 90
        new_coordinates(6) = coordinates(6) - 90
    end if
    if (coord_format == "Spherical") then
        new_coordinates(5) = 90 - coordinates(5)
        new_coordinates(6) = 90 -coordinates(6)
    end if

End   SUBROUTINE Coordinate_Transformation

!********************************************************
SUBROUTINE General_Coordinates_Format(Dim, old_coordinates, general_coodinates)
    IMPLICIT NONE

    INTEGER,INTENT(IN) ::  Dim
    real*8 , dimension(6), INTENT(INOUT) ::general_coodinates
    real*8 , dimension (Dim),INTENT(IN) :: old_coordinates


    if (Dim==2) then
        general_coodinates(1) = old_coordinates(1)  !R
        general_coodinates(2) = old_coordinates(2)  !b1
        general_coodinates(3) = 0d0  !b2
        general_coodinates(4) = 0d0  !phi
        general_coodinates(5) = 0d0  !c1
        general_coodinates(6) = 0d0  !c2

    elseif (Dim==3) then
        general_coodinates(1) = old_coordinates(1)  !R
        general_coodinates(2) = old_coordinates(2)  !b1
        general_coodinates(3) = 0d0  !b2
        general_coodinates(4) = 0d0  !phi
        general_coodinates(5) = old_coordinates(3)  !c1
        general_coodinates(6) = 0d0  !c2

    elseif (Dim==4) then
        general_coodinates(1) = old_coordinates(1)  !R
        general_coodinates(2) = old_coordinates(2)  !b1
        general_coodinates(3) = old_coordinates(3)  !b2
        general_coodinates(4) = old_coordinates(4)  !phi
        general_coodinates(5) = 0d0  !c1
        general_coodinates(6) = 0d0  !c2

    elseif (Dim==5) then
        general_coodinates(1) = old_coordinates(1)  !R
        general_coodinates(2) = old_coordinates(2)  !b1
        general_coodinates(3) = old_coordinates(3)  !b2
        general_coodinates(4) = old_coordinates(4)  !phi
        general_coodinates(5) = old_coordinates(5)  !c1
        general_coodinates(6) = 0d0  !c2
    elseif (Dim==6) then
        general_coodinates(1) = old_coordinates(1)  !R
        general_coodinates(2) = old_coordinates(2)  !b1
        general_coodinates(3) = old_coordinates(3)  !b2
        general_coodinates(4) = old_coordinates(4)  !phi
        general_coodinates(5) = old_coordinates(5)  !c1
        general_coodinates(6) = old_coordinates(6)  !c2

    end if

    RETURN
END SUBROUTINE General_Coordinates_Format





SUBROUTINE TotalEnergy_Calc (ind,TotalEnergy)

    implicit none


    Integer, INTENT(IN):: ind ! index of the coefficents
    real*8  , INTENT(INOut) ::TotalEnergy

    real*8   ::EM,ED,EI,term

    call Multipole_Sph3(ind, EM)
    call Induction_Sph3(ind, EI)
    call Dispersion_Sph3(ind, ED)

    !print*, "Electrostatic Energy: ", EM, " Induction Energy: ", EI, " Dispersion Energy: ", ED


    TotalEnergy = EM + ED + EI


end SUBROUTINE TotalEnergy_Calc



! Arg 1 [coordinates] : a coordinate vector [ R , b1, b2, phi] *the angles should be in degrees
! Arg 2 [Coeff_Address] address of the file which contains the longe range expansion coefficients
! Arg 3 [TotalEnergy]   Total Energy calculated


! version 4.0


SUBROUTINE Evaluate_LRF(TotalEnergy,coordinates,coord_format,XDIM,filename)

    use FitConstants, only: Coeff, max_T, Get_Coeff_Index
    use Geometry_Constant_v2, only: init_Tensors_v2

    IMPLICIT NONE

    real*8, INTENT(OUT) :: TotalEnergy
    INTEGER, INTENT(IN) :: XDIM
    real*8 ,dimension(:), INTENT(IN):: coordinates(XDIM)
    Character(*), INTENT(IN) :: coord_format
    Character(*), INTENT(IN) ::  filename
    real*8 ,dimension(6):: GeneralCoordenates,GeneralCoordenates1
    INTEGER :: i,CoeffIndex
    real*8 :: x1

    call Get_Coeff_Index(filename,CoeffIndex)

    x1 = 0d0
    do i=1,XDIM
        x1=x1+dabs(coordinates(i))
    enddo

    if (x1 <= 1d-10) then
        TotalEnergy = Coeff(CoeffIndex)%Zero
        return
    endif

    Call General_Coordinates_Format(XDIM, coordinates, GeneralCoordenates)
    Call Coordinate_Transformation(GeneralCoordenates,coord_format,GeneralCoordenates1)

    Call init_Tensors_v2(max_T,GeneralCoordenates1) ! Initializing in zero the new vectors v2
    Call TotalEnergy_Calc (CoeffIndex,TotalEnergy)

    return

END SUBROUTINE Evaluate_LRF
