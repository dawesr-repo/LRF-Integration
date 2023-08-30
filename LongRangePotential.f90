! Arg 1 [coordenates] : a coordenate vector [ R , b1, b2, phi] *the angles should be in degrees  
! Arg 2 [Coeff_Address] address of the file which contains the longe range expansion coefficients   
! Arg 3 [TotalEnergy]   Total Energy calculated 


! version 3.1.1



SUBROUTINE Long_Range_Potential(coordenates,TotalEnergy,filename)

 use Tensors_constant
 use FitConstants

 IMPLICIT NONE

 real*8, INTENT(INOUT)  ::  TotalEnergy
 real*8 ,dimension(6), INTENT(IN)  :: coordenates ! the angles are in degree
 character(len = *), INTENT(IN) :: filename
 real*8::testErr(52)

 real*8 , dimension(3):: Ar 
 real*8 , dimension(3):: Br
 real*8 , dimension(9):: C

 real*8 , dimension(11):: cal_coord
 integer::CoeffIndex
    


    
 call init_Tensors() ! Initializing in zero the new vectors

 call Get_Coeff_Index(filename,CoeffIndex) ! Initializing coefficients Fit for the file named as "filename"


 if (coordenates(1)==0d0 .and. coordenates(2)==0d0 .and. coordenates(3)==0d0 .and. coordenates(4)==0d0 &
     .and. coordenates(5)==0d0 .and. coordenates(6)==0d0) THEN
    TotalEnergy = Coeff(CoeffIndex)%Zero
 else
    Call Generate_Coordenates(coordenates,cal_coord,Ar,Br,C)
    Call TotalEnergy_Calc (cal_coord,Ar,Br,C,CoeffIndex,TotalEnergy,0,testErr)
  end if     

END SUBROUTINE Long_Range_Potential

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
