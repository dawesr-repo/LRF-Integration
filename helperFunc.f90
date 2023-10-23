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





SUBROUTINE TotalEnergy_Calc (ind,TotalEnergy,doTesting,testErr)
 use FitConstants, only:Coeff,C1,C2,C3
 use Geometry_Constant
 use Search, only: Initialize_Search
 implicit none

 
    Integer, INTENT(IN):: ind ! index of the coefficents 
    real*8  , INTENT(INOut) ::TotalEnergy,testErr(52)
    integer, INTENT(IN)::doTesting

    real*8   ::Ene,EM,ED,EH,EI,term,EI_Test
    Integer :: n ;

    real*8 :: Elect_energy(9),Dispe_energy(4),Induc_energy(6),Hyper_energy(3)

    

     call Initialize_Search()
   
  
     Ene = 0.d0
     EM  = 0.d0
     ED  = 0.d0
     EI  = 0.d0
     EH  = 0.d0

     Elect_energy  = 0.d0
     Dispe_energy  = 0.d0
     Induc_energy  = 0.d0
     Hyper_energy  = 0.d0

   
    do n = 1, 8
        EM  = 0.d0
        IF ( Coeff(ind)%M_Fit(n) > 0) THEN

            call Multipole_Sph2(ind,n, EM)
           
            testErr(5 + n) = EM
            Elect_energy(1+n) = EM
            Elect_energy(1) = Elect_energy(1)+EM 
            Ene = Ene+EM
            
         END IF 
         
     end do
     


     do n = 6,8
        IF (Coeff(ind)%D_Fit(n-5) > 0) THEN

            Call Dispersion_Sph2( ind,n, ED)

          
            testErr(20 + n - 5) = ED
            Dispe_energy(n-4) = ED
            Dispe_energy(1) = Dispe_energy(1)+ ED
            Ene = Ene+ED
          
         END IF 
        
     end do


     do n = 4, 8
        ! Call Induction_Sph2( ind,n, EI_Test)  
        ! write(*,*) "EI_Test: ",EI_Test

        IF (Coeff(ind)%I_Fit(n-3) > 0) THEN
            if(n==4)Then
                Call Induction_4_Sph2( Coeff(ind)%A_Mult,Coeff(ind)%B_Mult ,&
                Coeff(ind)%A_Pol,Coeff(ind)%B_Pol ,EI) 
               
            elseif (n==5)Then
                Call Induction_5_Sph2( Coeff(ind)%A_Mult,Coeff(ind)%B_Mult ,&
                Coeff(ind)%A_Pol,Coeff(ind)%B_Pol  ,EI)
            elseif(n==6)Then
                Call Induction_6_Sph2( Coeff(ind)%A_Mult,Coeff(ind)%B_Mult ,&
                Coeff(ind)%A_Pol,Coeff(ind)%B_Pol  ,EI)
                Call Induction_Sph2( ind,n, EI_Test) 
                write(*,*) "EI_6: ",C3*(C1*C2**n)*EI,EI_Test,C3*(C1*C2**n)*EI-EI_Test
            elseif(n==7)Then
                Call Induction_7_Sph2( Coeff(ind)%A_Mult,Coeff(ind)%B_Mult ,&
                Coeff(ind)%A_Pol,Coeff(ind)%B_Pol  ,EI)
                  Call Induction_Sph2( ind,n, EI_Test) 
                  write(*,*) "EI_7: ",C3*(C1*C2**n)*EI,EI_Test,C3*(C1*C2**n)*EI-EI_Test
            elseif(n==8)Then
                Call Induction_8_Sph2( Coeff(ind)%A_Mult,Coeff(ind)%B_Mult ,&
                Coeff(ind)%A_Pol,Coeff(ind)%B_Pol  ,EI) 
                Call Induction_Sph2( ind,n, EI_Test)      
                write(*,*) "EI_8: ",C3*(C1*C2**n)*EI,EI_Test,C3*(C1*C2**n)*EI-EI_Test   
            endif

            term = (C1*C2**n)*EI
            testErr(30 + n-3) = C3*term
            Induc_energy(n-2) = term
            Induc_energy(1) = Induc_energy(1)+term
            Ene = Ene+term
            
         END IF 
        
     end do



    !  do n = 6, 7
    !     IF (Coeff(ind)%H_Fit(n-5) > 0) THEN
    !         if(n==6)Then
    !             Call HyperPolarizability_6_Sph2( Coeff(ind)%A_Mult,&
    !                                 Coeff(ind)%B_Mult ,Coeff(ind)%A_HPol,Coeff(ind)%B_HPol,EH)
    !         elseif(n==7)Then
    !             Call HyperPolarizability_7_Sph2( Coeff(ind)%A_Mult,&
    !                                 Coeff(ind)%B_Mult ,Coeff(ind)%A_HPol,Coeff(ind)%B_HPol,EH)   
    !         endif
    !         term = (C1*C2**n)*EH
    !         testErr(42 + n-5) = C3*term
    !         Hyper_energy(n-4) = term
    !         Hyper_energy(1) = Hyper_energy(1)+term
    !         Ene = Ene+term
        
    !      END IF 
        
    !  end do



    TotalEnergy  = Ene
    if (doTesting>0)Then
        testErr(1) = TotalEnergy
        testErr(2) = Elect_energy(1)
        testErr(3) = Dispe_energy(1)
        testErr(4) = Induc_energy(1)
        testErr(5) = Hyper_energy(1)
    end if 


end SUBROUTINE TotalEnergy_Calc



! Arg 1 [coordenates] : a coordenate vector [ R , b1, b2, phi] *the angles should be in degrees  
! Arg 2 [Coeff_Address] address of the file which contains the longe range expansion coefficients   
! Arg 3 [TotalEnergy]   Total Energy calculated 


! version 3.1.1



SUBROUTINE Long_Range_Potential(coordenates,TotalEnergy,filename)

 use FitConstants
 use Geometry_Constant

 IMPLICIT NONE

 real*8, INTENT(INOUT)  ::  TotalEnergy
 real*8 ,dimension(6), INTENT(IN)  :: coordenates ! the angles are in degree
 character(len = *), INTENT(IN) :: filename
 real*8::testErr(52)

 integer::CoeffIndex
    



 if (coordenates(1)==0d0 .and. coordenates(2)==0d0 .and. coordenates(3)==0d0 .and. coordenates(4)==0d0 &
     .and. coordenates(5)==0d0 .and. coordenates(6)==0d0) THEN
    TotalEnergy = Coeff(CoeffIndex)%Zero
 else
    call init_Tensors() ! Initializing in zero the new vectors

    call Get_Coeff_Index(filename,CoeffIndex) ! Initializing coefficients Fit for the file named as "filename"

    Call Generate_Coordenates(coordenates)
    Call TotalEnergy_Calc (CoeffIndex,TotalEnergy,0,testErr)
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
