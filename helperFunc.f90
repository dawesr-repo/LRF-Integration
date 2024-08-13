SUBROUTINE Coordinate_Transformation(coordenates,coord_format,new_coordinates)
    IMPLICIT NONE
    real*8 ,dimension(6), INTENT(IN):: coordenates
    
    Character(len = 20), INTENT(IN) :: coord_format
    real*8 , dimension(6), INTENT(INOUT) :: new_coordinates
    
    new_coordinates=coordenates;

    
    if (coord_format == "Euler_ZYZ") then
        new_coordinates(5) = coordenates(5) - 90
        new_coordinates(6) = coordenates(6) - 90
    end if 
    if (coord_format == "Spherical") then
        new_coordinates(5) = 90 - coordenates(5) 
        new_coordinates(6) = 90 -coordenates(6) 
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
 use FitConstants, only:Coeff,C1,C2,C3
 implicit none

 
    Integer, INTENT(IN):: ind ! index of the coefficents 
    real*8  , INTENT(INOut) ::TotalEnergy

    real*8   ::Ene,EM,ED,EI,term
    Integer :: n ;

    real*8 :: Elect_energy(15),Dispe_energy(10),Induc_energy(12)

     TotalEnergy = 0.d0
     EM  = 0.d0
     ED  = 0.d0
     EI  = 0.d0
     

    call Multipole_Sph3(ind, EM)
    call Induction_Sph3(ind, EI)

    print*, "Electrostatic Energy: ", EM, " Induction Energy: ", EI


    TotalEnergy = EM + ED + EI  


end SUBROUTINE TotalEnergy_Calc



! Arg 1 [coordenates] : a coordenate vector [ R , b1, b2, phi] *the angles should be in degrees  
! Arg 2 [Coeff_Address] address of the file which contains the longe range expansion coefficients   
! Arg 3 [TotalEnergy]   Total Energy calculated 


! version 4.0



SUBROUTINE Long_Range_Potential(coordenates,TotalEnergy,filename)

 use FitConstants
 use Geometry_Constant_v2
 IMPLICIT NONE

 real*8, INTENT(INOUT)  ::  TotalEnergy
 real*8 ,dimension(6), INTENT(IN)  :: coordenates ! the angles are in degree
 character(len = *), INTENT(IN) :: filename
 real*8::testErr(52),T1,T2,r1

 integer::CoeffIndex
 Character(len = 1), dimension(3) :: cpn = ["0","c","s"]


 if (coordenates(1)==0d0 .and. coordenates(2)==0d0 .and. coordenates(3)==0d0 .and. coordenates(4)==0d0 &
     .and. coordenates(5)==0d0 .and. coordenates(6)==0d0) THEN
    TotalEnergy = Coeff(CoeffIndex)%Zero
 else
    
    call Get_Coeff_Index(filename,CoeffIndex) ! Initializing coefficients Fit for the file named as "filename"
    
    call init_Tensors_v2(15,coordenates) ! Initializing in zero the new vectors v2

    Call TotalEnergy_Calc (CoeffIndex,TotalEnergy)
    print*, "Total Energy: ", TotalEnergy
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
