
!********************************************************
SUBROUTINE Approx_3_Sph2(cal_coord,Ar,Br,C ,A_Multipoles,B_Multipoles , Approx_3_Energy)
    IMPLICIT NONE
    
    !   NEED TO DECLARE ALL THE SUBROUTINE ARGUMENTS and
    !   ANY OTHER VARIABLES LOCAL TO THE SUBROUTINE

    real*8, INTENT(INOUT)  :: Approx_3_Energy
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 :: R 
    real*8 , dimension(11), INTENT(IN):: cal_coord
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    real*8 :: Q20_0,Q20_1,Q11_0

    R =cal_coord(1);
    
    Call Q20(Ar,Br,C ,A_Multipoles,B_Multipoles,0,Q20_0) 
    Call Q20(Ar,Br,C ,A_Multipoles,B_Multipoles,1,Q20_1)

    Call Q11(Ar,Br,C ,A_Multipoles,B_Multipoles,0,Q11_0)



    
    Approx_3_Energy =   (Q20_0 + Q20_1 + Q11_0)/ R**3

    RETURN
END SUBROUTINE Approx_3_Sph2
!********************************************************



SUBROUTINE Q20(Ar,Br,C ,A_Multipoles,B_Multipoles, ind , result)
        IMPLICIT NONE

   
        real*8, INTENT(INOUT) ::result 
        real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
        real*8 , dimension(3), INTENT(IN):: Ar 
        real*8 , dimension(3), INTENT(IN):: Br
        real*8 , dimension(9), INTENT(IN):: C
        Integer, INTENT(IN):: ind
        real*8 :: q
        real*8 , dimension(5):: Qd
        
   
        
        if (ind==0) then
        q =  B_Multipoles(1)
        Qd  = A_Multipoles(5:9) 
        

        else
            q =  A_Multipoles(1)
            Qd  = B_Multipoles(5:9) 
        end if

     

        Call T_l0(Ar,Br,C,q,Qd,ind,2, result)

       
    
        RETURN
END SUBROUTINE Q20

    
SUBROUTINE Q11(Ar,Br,C ,A_Multipoles,B_Multipoles, ind , result)
            IMPLICIT NONE
    
            
            real*8, INTENT(INOUT) ::result 
            real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
            real*8 , dimension(3), INTENT(IN):: Ar 
            real*8 , dimension(3), INTENT(IN):: Br
            real*8 , dimension(9), INTENT(IN):: C
            Integer, INTENT(IN):: ind
            real*8 , dimension(3):: mA,mB
            
           
            
            
            mA =  A_Multipoles(2:4)
            mB  = B_Multipoles(2:4) 
            
    
    
          
    
            Call T_ll(Ar,Br,C,mA,mB,1,1 , result)
        
            RETURN
END SUBROUTINE Q11