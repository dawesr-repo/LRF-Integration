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

MODULE constants_1
 implicit none
 save
 public
  Integer, dimension (8) :: M_Fit     
  Integer, dimension (3):: D_Fit     
  Integer, dimension (5):: I_Fit     
  Integer, dimension (2):: H_Fit     
  real*8 , dimension(1195):: coeff_arr
  Real*8  ::   Zero  
END MODULE constants_1

MODULE constants_2
 implicit none
 save
 public
  Integer, dimension (8) :: M_Fit     
  Integer, dimension (3):: D_Fit     
  Integer, dimension (5):: I_Fit     
  Integer, dimension (2):: H_Fit     
  real*8 , dimension(1195):: coeff_arr
  Real*8  ::   Zero  
END MODULE constants_2



SUBROUTINE Long_Range_Potential(coordenates,TotalEnergy,filename)
 use Tensors_constant
 use constants
 IMPLICIT NONE

 real*8, INTENT(INOUT)  ::  TotalEnergy
 real*8 ,dimension(6), INTENT(IN)  :: coordenates ! the angles are in degree
 character(len = *), INTENT(IN) :: filename
 real*8::testErr(52)

 real*8 , dimension(3):: Ar 
 real*8 , dimension(3):: Br
 real*8 , dimension(9):: C

 real*8 , dimension(11):: cal_coord
 character(len = 200):: fName
    

 integer :: initflag
 save initflag
 data initflag /1/

 call init_Tensors() ! Initializing in zero the new vectors

 IF(initflag==1)THEN! initialize 
   CALL Prep_Param(filename,coeff_arr,M_Fit ,D_Fit,I_Fit,H_Fit,Zero)
   initflag=2  
 ENDIF

 if (coordenates(1)==0d0 .and. coordenates(2)==0d0 .and. coordenates(3)==0d0 .and. coordenates(4)==0d0 &
     .and. coordenates(5)==0d0 .and. coordenates(6)==0d0) THEN
    TotalEnergy = Zero
 else
    Call Generate_Coordenates(coordenates,cal_coord,Ar,Br,C)
    call TotalEnergy_Calc(cal_coord,Ar,Br,C,coeff_arr, M_Fit ,D_Fit,I_Fit,H_Fit,TotalEnergy,0,testErr)
  end if     

END SUBROUTINE Long_Range_Potential











SUBROUTINE Long_Range_Potential_1(coordenates,TotalEnergy,filename)
 use Tensors_constant
 use constants_1
 IMPLICIT NONE

 real*8, INTENT(INOUT)  ::  TotalEnergy
 real*8 ,dimension(6), INTENT(IN)  :: coordenates ! the angles are in degree
 character(len = *), INTENT(IN) :: filename
 real*8::testErr(52)

 real*8 , dimension(3):: Ar 
 real*8 , dimension(3):: Br
 real*8 , dimension(9):: C

 real*8 , dimension(11):: cal_coord
 character(len = 200):: fName
    
 integer :: initflag1
 save initflag1
 data initflag1 /1/
 
 call init_Tensors() ! Initializing in zero the new vectors

 IF(initflag1==1)THEN! initialize 
   CALL Prep_Param(filename,coeff_arr,M_Fit ,D_Fit,I_Fit,H_Fit,Zero)
   initflag1=2  
 ENDIF

 if (coordenates(1)==0d0 .and. coordenates(2)==0d0 .and. coordenates(3)==0d0 .and. coordenates(4)==0d0 &
     .and. coordenates(5)==0d0 .and. coordenates(6)==0d0) THEN
    TotalEnergy = Zero
 else
    Call Generate_Coordenates(coordenates,cal_coord,Ar,Br,C)
    call TotalEnergy_Calc(cal_coord,Ar,Br,C,coeff_arr, M_Fit ,D_Fit,I_Fit,H_Fit,TotalEnergy,0,testErr)
  end if     

END SUBROUTINE Long_Range_Potential_1


SUBROUTINE Long_Range_Potential_2(coordenates,TotalEnergy,filename)
 use Tensors_constant
 use constants_2
 IMPLICIT NONE

 real*8, INTENT(INOUT)  ::  TotalEnergy
 real*8 ,dimension(6), INTENT(IN)  :: coordenates ! the angles are in degree
 character(len = *), INTENT(IN) :: filename
 real*8::testErr(52)

 real*8 , dimension(3):: Ar 
 real*8 , dimension(3):: Br
 real*8 , dimension(9):: C

 real*8 , dimension(11):: cal_coord
 character(len = 200):: fName
    

 integer :: initflag2
 save initflag2
 data initflag2 /1/

 call init_Tensors() ! Initializing in zero the new vectors

 IF(initflag2==1)THEN! initialize 
   CALL Prep_Param(filename,coeff_arr,M_Fit ,D_Fit,I_Fit,H_Fit,Zero)
   initflag2=2  
 ENDIF

 if (coordenates(1)==0d0 .and. coordenates(2)==0d0 .and. coordenates(3)==0d0 .and. coordenates(4)==0d0 &
     .and. coordenates(5)==0d0 .and. coordenates(6)==0d0) THEN
    TotalEnergy = Zero
 else
    Call Generate_Coordenates(coordenates,cal_coord,Ar,Br,C)
    call TotalEnergy_Calc(cal_coord,Ar,Br,C,coeff_arr, M_Fit ,D_Fit,I_Fit,H_Fit,TotalEnergy,0,testErr)
  end if     

END SUBROUTINE Long_Range_Potential_2

