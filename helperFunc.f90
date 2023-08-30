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



SUBROUTINE TotalEnergy_Calc (cal_coord,Ar,Br,C,CoeffIndex,TotalEnergy,doTesting,testErr)

    use FitConstants
    implicit none
    real*8, parameter ::  C1=627.5095d0,C2=0.529177249d0,Const=349.757d0


    real*8 , dimension(11) , INTENT(IN):: cal_coord
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C

    real*8, INTENT(INOUT) ::TotalEnergy
    integer,intent(in)::CoeffIndex,doTesting
    real*8,intent(INOUT)::testErr(52)


    real*8   ::Ene,EM,ED,EH,EI,T10,T20,T30,T40,cal_coord_temp(11)
    Integer :: n ;
    real*8::Multipole_Energies(8),Ind_Energ(5),Hyp_Energ(2),Dispersion_Energies(3)

    real*8 , dimension(64) :: A_Mult,B_Mult
    real*8 , dimension(57) :: A_Pol,B_Pol
    real*8 , dimension(40) :: A_HPol,B_HPol
    real*8 , dimension(873):: Disp_AB

    Integer, dimension (8):: M_Fit      
    Integer, dimension (3):: D_Fit     
    Integer, dimension (5):: I_Fit     
    Integer, dimension (2):: H_Fit  

    real*8 :: Elect_energy(9),Dispe_energy(4),Induc_energy(6),Hyper_energy(3),term

    

    cal_coord_temp = cal_coord

    M_Fit = Coeff(CoeffIndex)%M_Fit
    I_Fit = Coeff(CoeffIndex)%I_Fit
    D_Fit = Coeff(CoeffIndex)%D_Fit
    H_Fit = Coeff(CoeffIndex)%H_Fit

    A_Mult = Coeff(CoeffIndex)%A_Mult
    B_Mult = Coeff(CoeffIndex)%A_Mult
    A_Pol  = Coeff(CoeffIndex)%A_Pol
    B_Pol  = Coeff(CoeffIndex)%B_Pol
    A_HPol = Coeff(CoeffIndex)%A_HPol
    B_HPol = Coeff(CoeffIndex)%B_HPol
    Disp_AB= Coeff(CoeffIndex)%Disp


   
  
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
        IF (M_Fit(n) > 0) THEN
            if (n==1)Then
                Call Approx_1_Sph2(cal_coord_temp,Ar,Br,C , A_Mult,B_Mult ,EM)
            elseif(n==2)Then
                Call Approx_2_Sph2(cal_coord_temp,Ar,Br,C , A_Mult,B_Mult ,EM)
            elseif(n==3)Then
                Call Approx_3_Sph2(cal_coord_temp,Ar,Br,C , A_Mult,B_Mult ,EM)
            elseif(n==4)Then
                Call Approx_4_Sph2(cal_coord_temp,Ar,Br,C , A_Mult,B_Mult ,EM)    
            elseif (n==5)Then
                Call Approx_5_Sph2(cal_coord_temp,Ar,Br,C , A_Mult,B_Mult ,EM)
            elseif(n==6)Then
                Call Approx_6_Sph2(cal_coord_temp,Ar,Br,C , A_Mult,B_Mult ,EM)
            elseif(n==7)Then
                Call Approx_7_Sph2(cal_coord_temp,Ar,Br,C , A_Mult,B_Mult ,EM)
            elseif(n==8)Then
                Call Approx_8_Sph2(cal_coord_temp,Ar,Br,C , A_Mult,B_Mult ,EM)         
            endif

            term = (C1*C2**n)*EM
            testErr(5 + n) = Const*term
            Elect_energy(1+n) = term
            Elect_energy(1) = Elect_energy(1)+term 
            Ene = Ene+term
            
         END IF 
         
     end do
     


     do n = 6,8
        IF (D_Fit(n-5) > 0) THEN
            if(n==6)Then
                Call Dispersion_6_Sph2(cal_coord_temp,Ar,Br,C, Disp_AB  ,ED)
            elseif(n==7)Then
                Call Dispersion_7_Sph2(cal_coord_temp,Ar,Br,C, Disp_AB  ,ED)
            elseif(n==8)Then
                Call Dispersion_8_Sph2(cal_coord_temp,Ar,Br,C, Disp_AB  ,ED)         
            endif

            term = (C1*C2**n)*ED
            testErr(20 + n - 5) = Const*term
            Dispe_energy(n-4) = term
            Dispe_energy(1) = Dispe_energy(1)+term
            Ene = Ene+term
         END IF 
        
     end do


     do n = 4, 8
        IF (I_Fit(n-3) > 0) THEN
            if(n==4)Then
                Call Induction_4_Sph2(cal_coord_temp,Ar,Br,C , A_Mult,B_Mult ,A_Pol,B_Pol ,EI)    
            elseif (n==5)Then
                Call Induction_5_Sph2(cal_coord_temp,Ar,Br,C , A_Mult,B_Mult ,A_Pol,B_Pol  ,EI)
            elseif(n==6)Then
                Call Induction_6_Sph2(cal_coord_temp,Ar,Br,C , A_Mult,B_Mult ,A_Pol,B_Pol  ,EI)
            elseif(n==7)Then
                Call Induction_7_Sph2(cal_coord_temp,Ar,Br,C , A_Mult,B_Mult ,A_Pol,B_Pol  ,EI)
            elseif(n==8)Then
                Call Induction_8_Sph2(cal_coord_temp,Ar,Br,C , A_Mult,B_Mult ,A_Pol,B_Pol  ,EI)         
            endif

            term = (C1*C2**n)*EI
            testErr(30 + n-3) = Const*term
            Induc_energy(n-2) = term
            Induc_energy(1) = Induc_energy(1)+term
            Ene = Ene+term
         END IF 
        
     end do



     do n = 6, 7
        IF (H_Fit(n-5) > 0) THEN
            if(n==6)Then
                Call HyperPolarizability_6_Sph2(cal_coord_temp,Ar,Br,C , A_Mult,B_Mult ,A_HPol,B_HPol,EH)
            elseif(n==7)Then
                Call HyperPolarizability_7_Sph2(cal_coord_temp,Ar,Br,C , A_Mult,B_Mult ,A_HPol,B_HPol,EH)   
            endif
            term = (C1*C2**n)*EH
            testErr(42 + n-5) = Const*term
            Hyper_energy(n-4) = term
            Hyper_energy(1) = Hyper_energy(1)+term
            Ene = Ene+term
         END IF 
        
     end do



    TotalEnergy  = Const*Ene
    if (doTesting>0)Then
        testErr(1) = TotalEnergy
        testErr(2) = Const*Elect_energy(1)
        testErr(3) = Const*Dispe_energy(1)
        testErr(4) = Const*Induc_energy(1)
        testErr(5) = Const*Hyper_energy(1)
    end if 


end SUBROUTINE TotalEnergy_Calc



! SUBROUTINE Prep_Param(Coeff_Address, coeff_arr,M_Fit ,D_Fit,I_Fit,H_Fit,Zero)
!     IMPLICIT NONE

!     !Multipoles !
    
!     real*8 , dimension(1)   :: qA,qB !q
!     real*8 , dimension(3)   :: mA,mB !m
!     real*8 , dimension(5)   :: QdA,QdB !Qd 
!     real*8 , dimension(7)   :: OA,OB !O 
!     real*8 , dimension(9)   :: PhiA,PhiB !Phi 
!     real*8 , dimension(11)   :: M5A,M5B !M5
!     real*8 , dimension(13)   :: M6A,M6B !M6 
!     real*8 , dimension(15)   :: M7A,M7B !M7

!     real*8 , dimension(64)   :: A_Mult,B_Mult 

!     !Polarizability!

!     real*8 , dimension(6)   :: mmA,mmB 
!     real*8 , dimension(15)  :: mQdA,mQdB
!     real*8 , dimension(15)  :: QdQdA,QdQdB
!     real*8 , dimension(21)  :: mOA,mOB

!     real*8 , dimension(57)   :: A_Pol,B_Pol 

!     !Hyperpolarizability!

!     real*8 , dimension(10)   :: mmmA,mmmB 
!     real*8 , dimension(30)  :: mmQdA,mmQdB

!     real*8 , dimension(40)   :: A_HPol,B_HPol 

!     !Dispersion!

!     real*8 , dimension(36)   :: mm_mm 
!     real*8 , dimension(90)  ::  mm_mQd,mQd_mm,mm_QdQd,QdQd_mm
!     real*8 , dimension(126)   :: mO_mm,mm_mO 
!     real*8 , dimension(225)   :: mQd_mQd 

!     real*8 , dimension(873)   :: Disp



!     real*8 , dimension(1195), INTENT(INOUT) :: coeff_arr

!     Character(len = *), INTENT(IN)   ::  Coeff_Address

!     Integer, dimension (8),  INTENT(INOUT) :: M_Fit      
!     Integer, dimension (3),  INTENT(INOUT) :: D_Fit     
!     Integer, dimension (5),  INTENT(INOUT) :: I_Fit     
!     Integer, dimension (2),  INTENT(INOUT) :: H_Fit     

!     Real*8 , INTENT(out) :: Zero                         


!     Character(len = 20) :: row
!     ! Integer , dimension(7) :: DataColumn ! R , Cos_b1 , Cos_b2 , alpha ,  Cos_c1 , Cos_c2 , Energy    

!     !write(*,*) 'Reading file'

!     Open( 10, file = Coeff_Address )

!     Read( 10, *) row
!     Read( 10, *) row
!     Read( 10, *) row
!     Read( 10, *) row
!     Read( 10, *) row
!     Read( 10, *) row
!     Read( 10, *) row

!     read(10, *)  M_Fit
!     read(10, *)  D_Fit
!     read(10, *)  I_Fit
!     read(10, *)  H_Fit

!     Read( 10, *) row
!     read(10, *) Zero

!     Read( 10, *) row
!     read(10, *) qA
!     read(10, *) mA
!     read(10, *) QdA
!     read(10, *) OA
!     read(10, *) PhiA
!     read(10, *) M5A
!     read(10, *) M6A
!     read(10, *) M7A

!     Read( 10, *) row
!     read(10, *) qB
!     read(10, *) mB
!     read(10, *) QdB
!     read(10, *) OB
!     read(10, *) PhiB
!     read(10, *) M5B
!     read(10, *) M6B
!     read(10, *) M7B

!     Read( 10, *) row
!     read(10, *) mmA
!     read(10, *) mQdA
!     read(10, *) QdQdA
!     read(10, *) mOA

!     Read( 10, *) row
!     read(10, *) mmB
!     read(10, *) mQdB
!     read(10, *) QdQdB
!     read(10, *) mOB

!     Read( 10, *) row
!     read(10, *) mmmA
!     read(10, *) mmQdA

    
!     Read( 10, *) row
!     read(10, *) mmmB
!     read(10, *) mmQdB

!     Read( 10, *) row
!     read(10, *) mm_mm
!     read(10, *) mQd_mm
!     read(10, *) mm_mQd
!     read(10, *) mO_mm
!     read(10, *) mm_mO
!     read(10, *) QdQd_mm
!     read(10, *) mm_QdQd
!     read(10, *) mQd_mQd


!     close(10) 



!     A_Mult(1:1) = qA
!     A_Mult(2:4) = mA
!     A_Mult(5:9) = QdA
!     A_Mult(10:16) = OA
!     A_Mult(17:25) = PhiA
!     A_Mult(26:36) = M5A
!     A_Mult(37:49) = M6A
!     A_Mult(50:64) = M7A

!     B_Mult(1:1) = qB
!     B_Mult(2:4) = mB
!     B_Mult(5:9) = QdB
!     B_Mult(10:16) = OB
!     B_Mult(17:25) = PhiB
!     B_Mult(26:36) = M5B
!     B_Mult(37:49) = M6B
!     B_Mult(50:64) = M7B


!     A_Pol(1:6) = mmA
!     A_Pol(7:21) = mQdA
!     A_Pol(22:36) = QdQdA
!     A_Pol(37:57) = mOA

    
!     B_Pol(1:6) = mmB
!     B_Pol(7:21) = mQdB
!     B_Pol(22:36) = QdQdB
!     B_Pol(37:57) = mOB

!     A_HPol(1:10) = mmmA
!     A_HPol(11:40) = mmQdA

!     B_HPol(1:10) = mmmB
!     B_HPol(11:40) = mmQdB

!     Disp(1:36) = mm_mm
!     Disp(37:126) = mQd_mm
!     Disp(127:216) = mm_mQd
!     Disp(217:342) = mO_mm
!     Disp(343:468) = mm_mO
!     Disp(469:558) = QdQd_mm
!     Disp(559:648) = mm_QdQd
!     Disp(649:873) = mQd_mQd


!     coeff_arr(1:64)   = A_Mult
!     coeff_arr(65:128) = B_Mult

!     coeff_arr(129:185)   = A_Pol
!     coeff_arr(186:242)   = B_Pol

!     coeff_arr(243:282)   = A_HPol
!     coeff_arr(283:322)   = B_HPol

!     coeff_arr(323:1195)   = Disp


!     RETURN
! END SUBROUTINE Prep_Param



