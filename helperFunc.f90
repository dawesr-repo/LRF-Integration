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




SUBROUTINE Generate_Coordenates(coordenates,cal_coord,Ar,Br,C)

    real*8 , dimension(11), INTENT(INOUT):: cal_coord
    real*8 ,dimension(6), INTENT(IN)  :: coordenates ! the angles are in degree
    real*8, parameter :: pii = DACOS(-1.d0)  
    real*8 , dimension(3), INTENT(INOUT):: Ar 
    real*8 , dimension(3), INTENT(INOUT):: Br
    real*8 , dimension(9), INTENT(INOUT):: C

    real*8 :: cos_b1,cos_b2,cos_c1,cos_c2,sin_b1,sin_b2,sin_c1,sin_c2,cos_phi,sin_phi

        cal_coord(1) = coordenates(1)
           

        cal_coord(2)  = DCOS(coordenates(2)*pii/180d0)
        cal_coord(3) = DSIN(coordenates(2)*pii/180d0)

        cal_coord(4) = DCOS(coordenates(3)*pii/180d0)
        cal_coord(5) = DSIN(coordenates(3)*pii/180d0)
    
        cal_coord(6) = DCOS(coordenates(4)*pii/180d0)
        cal_coord(7) = DSIN(coordenates(4)*pii/180d0)

        cal_coord(8) = DCOS(coordenates(5)*pii/180d0)
        cal_coord(9) = DSIN(coordenates(5)*pii/180d0)

        cal_coord(10) = DCOS(coordenates(6)*pii/180d0)
        cal_coord(11) = DSIN(coordenates(6)*pii/180d0)

       ! write(*,*)"cal_coord",cal_coord(2),cal_coord(4),coordenates(4)*pii/180d0

        cos_b1 =    cal_coord(2)
        sin_b1 =    cal_coord(3)
        cos_b2 =    cal_coord(4)
        sin_b2 =    cal_coord(5)
        cos_phi =   cal_coord(6)
        sin_phi =   cal_coord(7)
        cos_c1 =    cal_coord(8)
        sin_c1 =    cal_coord(9)
        cos_c2 =    cal_coord(10)
        sin_c2 =    cal_coord(11)

        

        Ar(1) = cos_b1           !Az
        Ar(2) = cos_c1*sin_b1    !Ay
        Ar(3) = sin_b1*sin_c1    !Ax
        
        
        Br(1) = -cos_b2           !Bz
        Br(2) = -cos_c2*sin_b2    !By
        Br(3) = -sin_b2*sin_c2    !Bx
        
        C(1)= cos_b1*cos_b2 + cos_phi*sin_b1*sin_b2 !Czz
        C(2)= cos_b2*sin_b1*sin_c1 - sin_b2 *(cos_c1*sin_phi + cos_phi*cos_b1*sin_c1)  !Cxz
        C(3)= cos_b2*cos_c1*sin_b1 + sin_b2 *(-cos_phi*cos_b1*cos_c1 + sin_phi*sin_c1) !Cyz

        C(4)= cos_c2*sin_phi*sin_b1 + (-cos_phi*cos_b2*sin_b1 + cos_b1*sin_b2)*sin_c2; !Czx
        C(5)= -cos_b1*cos_c2*sin_phi*sin_c1 + (cos_b2*cos_c1*sin_phi + sin_b1*sin_b2*sin_c1)*sin_c2 &
        + cos_phi *(cos_c1*cos_c2 + cos_b1*cos_b2*sin_c1*sin_c2) !Cxx
        C(6)= cos_c1*sin_b1*sin_b2*sin_c2 + cos_b1*cos_c1 *(-cos_c2*sin_phi + cos_phi*cos_b2*sin_c2) - &
        sin_c1 *(cos_phi*cos_c2 + cos_b2*sin_phi*sin_c2)   !Cyx

        C(7)= -cos_phi*cos_b2*cos_c2*sin_b1 + cos_b1*cos_c2*sin_b2 - sin_phi*sin_b1*sin_c2 !Czy
        C(8)= cos_c2*sin_b1*sin_b2*sin_c1 + cos_b2*cos_c2 *(cos_c1*sin_phi + cos_phi*cos_b1*sin_c1) &
        + (-cos_phi*cos_c1 + cos_b1*sin_phi*sin_c1)*sin_c2 !Cxy
        C(9)= -cos_b2*cos_c2*sin_phi*sin_c1 + cos_c1 *(cos_c2*sin_b1*sin_b2 + cos_b1*sin_phi*sin_c2) &
        + cos_phi*(cos_b1*cos_b2*cos_c1*cos_c2 + sin_c1*sin_c2) !Cyy
       

End   SUBROUTINE Generate_Coordenates


SUBROUTINE TotalEnergy_Calc (cal_coord,Ar,Br,C,ind,TotalEnergy,doTesting,testErr)
 use FitConstants
 implicit none

 
    Integer, INTENT(IN):: ind ! index of the coefficents 
    real*8 , dimension(11) , INTENT(IN):: cal_coord
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    real*8  , INTENT(INOut) ::TotalEnergy,testErr(52)
    integer, INTENT(IN)::doTesting

    real*8   ::Ene,EM,ED,EH,EI,cal_coord_temp(11)
    Integer :: n ;

    real*8 :: Elect_energy(9),Dispe_energy(4),Induc_energy(6),Hyper_energy(3),term

    

    cal_coord_temp = cal_coord
   
  
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
        IF ( Coeff(ind)%M_Fit(n) > 0) THEN
            if (n==1)Then
                Call Approx_1_Sph2(cal_coord_temp,Ar,Br,C , Coeff(ind)%A_Mult,Coeff(ind)%B_Mult ,EM)
            elseif(n==2)Then
                Call Approx_2_Sph2(cal_coord_temp,Ar,Br,C , Coeff(ind)%A_Mult,Coeff(ind)%B_Mult ,EM)
            elseif(n==3)Then
                Call Approx_3_Sph2(cal_coord_temp,Ar,Br,C , Coeff(ind)%A_Mult,Coeff(ind)%B_Mult ,EM)
            elseif(n==4)Then
                Call Approx_4_Sph2(cal_coord_temp,Ar,Br,C , Coeff(ind)%A_Mult,Coeff(ind)%B_Mult ,EM)    
            elseif (n==5)Then
                Call Approx_5_Sph2(cal_coord_temp,Ar,Br,C , Coeff(ind)%A_Mult,Coeff(ind)%B_Mult ,EM)
            elseif(n==6)Then
                Call Approx_6_Sph2(cal_coord_temp,Ar,Br,C , Coeff(ind)%A_Mult,Coeff(ind)%B_Mult ,EM)
            elseif(n==7)Then
                Call Approx_7_Sph2(cal_coord_temp,Ar,Br,C , Coeff(ind)%A_Mult,Coeff(ind)%B_Mult ,EM)
            elseif(n==8)Then
                Call Approx_8_Sph2(cal_coord_temp,Ar,Br,C , Coeff(ind)%A_Mult,Coeff(ind)%B_Mult ,EM)         
            endif

            term = (C1*C2**n)*EM
            testErr(5 + n) = C3*term
            Elect_energy(1+n) = term
            Elect_energy(1) = Elect_energy(1)+term 
            Ene = Ene+term
            
         END IF 
         
     end do
     


     do n = 6,8
        IF (Coeff(ind)%D_Fit(n-5) > 0) THEN
            if(n==6)Then
                Call Dispersion_6_Sph2(cal_coord_temp,Ar,Br,C, Coeff(ind)%Disp  ,ED)
            elseif(n==7)Then
                Call Dispersion_7_Sph2(cal_coord_temp,Ar,Br,C, Coeff(ind)%Disp  ,ED)
            elseif(n==8)Then
                Call Dispersion_8_Sph2(cal_coord_temp,Ar,Br,C, Coeff(ind)%Disp  ,ED)         
            endif

            term = (C1*C2**n)*ED
            testErr(20 + n - 5) = C3*term
            Dispe_energy(n-4) = term
            Dispe_energy(1) = Dispe_energy(1)+term
            Ene = Ene+term
          
         END IF 
        
     end do


     do n = 4, 8
        IF (Coeff(ind)%I_Fit(n-3) > 0) THEN
            if(n==4)Then
                Call Induction_4_Sph2(cal_coord_temp,Ar,Br,C , Coeff(ind)%A_Mult,Coeff(ind)%B_Mult ,&
                Coeff(ind)%A_Pol,Coeff(ind)%B_Pol ,EI)    
            elseif (n==5)Then
                Call Induction_5_Sph2(cal_coord_temp,Ar,Br,C , Coeff(ind)%A_Mult,Coeff(ind)%B_Mult ,&
                Coeff(ind)%A_Pol,Coeff(ind)%B_Pol  ,EI)
            elseif(n==6)Then
                Call Induction_6_Sph2(cal_coord_temp,Ar,Br,C , Coeff(ind)%A_Mult,Coeff(ind)%B_Mult ,&
                Coeff(ind)%A_Pol,Coeff(ind)%B_Pol  ,EI)
            elseif(n==7)Then
                Call Induction_7_Sph2(cal_coord_temp,Ar,Br,C , Coeff(ind)%A_Mult,Coeff(ind)%B_Mult ,&
                Coeff(ind)%A_Pol,Coeff(ind)%B_Pol  ,EI)
            elseif(n==8)Then
                Call Induction_8_Sph2(cal_coord_temp,Ar,Br,C , Coeff(ind)%A_Mult,Coeff(ind)%B_Mult ,&
                Coeff(ind)%A_Pol,Coeff(ind)%B_Pol  ,EI)         
            endif

            term = (C1*C2**n)*EI
            testErr(30 + n-3) = C3*term
            Induc_energy(n-2) = term
            Induc_energy(1) = Induc_energy(1)+term
            Ene = Ene+term
            
         END IF 
        
     end do



     do n = 6, 7
        IF (Coeff(ind)%H_Fit(n-5) > 0) THEN
            if(n==6)Then
                Call HyperPolarizability_6_Sph2(cal_coord_temp,Ar,Br,C , Coeff(ind)%A_Mult,&
                                    Coeff(ind)%B_Mult ,Coeff(ind)%A_HPol,Coeff(ind)%B_HPol,EH)
            elseif(n==7)Then
                Call HyperPolarizability_7_Sph2(cal_coord_temp,Ar,Br,C , Coeff(ind)%A_Mult,&
                                    Coeff(ind)%B_Mult ,Coeff(ind)%A_HPol,Coeff(ind)%B_HPol,EH)   
            endif
            term = (C1*C2**n)*EH
            testErr(42 + n-5) = C3*term
            Hyper_energy(n-4) = term
            Hyper_energy(1) = Hyper_energy(1)+term
            Ene = Ene+term
        
         END IF 
        
     end do



    TotalEnergy  = C3*Ene
    if (doTesting>0)Then
        testErr(1) = TotalEnergy
        testErr(2) = C3*Elect_energy(1)
        testErr(3) = C3*Dispe_energy(1)
        testErr(4) = C3*Induc_energy(1)
        testErr(5) = C3*Hyper_energy(1)
    end if 


end SUBROUTINE TotalEnergy_Calc
