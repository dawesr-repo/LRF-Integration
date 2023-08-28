! Arg 1 [coordenates] : a coordenate vector [ R , b1, b2, phi] *the angles should be in degrees  
! Arg 2 [Coeff_Address] address of the file which contains the longe range expansion coefficients   
! Arg 3 [TotalEnergy]   Total Energy calculated 


! version 3.1.1


MODULE constants
 implicit none
 save
 public
  Integer, dimension (8) :: M_Fit     
  Integer, dimension (3):: D_Fit     
  Integer, dimension (5):: I_Fit     
  Integer, dimension (2):: H_Fit     
  real*8 , dimension(1195):: coeff_arr
  Real*8  ::   Zero  
END MODULE constants


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
    TotalEnergy = Zero
 else
    Call Generate_Coordenates(coordenates,cal_coord,Ar,Br,C)
    call TotalEnergy_Calc(cal_coord,Ar,Br,C,CoeffIndex,TotalEnergy,0,testErr)
  end if     

END SUBROUTINE Long_Range_Potential
