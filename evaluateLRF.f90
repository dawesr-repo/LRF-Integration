SUBROUTINE evaluateLR(coordinates,XDIM,E1,filename)
IMPLICIT NONE

real*8, INTENT(OUT) :: E1
INTEGER, INTENT(IN) :: XDIM
real*8 ,dimension(:), INTENT(IN):: coordinates(XDIM)
Character(len = 20) :: coord_format = "Euler_ZYZ" !for Xdim =3, use coord_format ="Spherical" for Autosurf input 
real*8 ,dimension(6):: GeneralCoordenates,GeneralCoordenates1
INTEGER :: i
real*8 :: x1
Character(*), INTENT(IN) ::  filename
    
x1=0d0
do i=1,XDIM
  x1=x1+dabs(coordinates(i))
enddo

if (x1 <= 1d-10) then
  GeneralCoordenates=0d0
  CALL Long_Range_Potential(GeneralCoordenates,E1,filename)
  return
endif

Call General_Coordinates_Format(XDIM, coordinates, GeneralCoordenates)
Call Coordinate_Transformation(GeneralCoordenates,coord_format,GeneralCoordenates1)
CALL Long_Range_Potential(GeneralCoordenates1,E1,filename)
return

END SUBROUTINE evaluateLR
