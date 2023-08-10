
    !********************************************************
    SUBROUTINE Approx_6_Sph2(cal_coord,Ar,Br,C ,A_Multipoles,B_Multipoles , Approx_6_Energy)
    IMPLICIT NONE
    
    !   NEED TO DECLARE ALL THE SUBROUTINE ARGUMENTS and
    !   ANY OTHER VARIABLES LOCAL TO THE SUBROUTINE

    real*8, INTENT(INOUT)  :: Approx_6_Energy
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 :: R 
    real*8 , dimension(11), INTENT(IN):: cal_coord
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    real*8 :: Q41_0,Q41_1,Q50_0,Q50_1,Q32_0,Q32_1

    R =cal_coord(1);
    
    Call Q41(Ar,Br,C ,A_Multipoles,B_Multipoles,0,Q41_0) 
    Call Q41(Ar,Br,C ,A_Multipoles,B_Multipoles,1,Q41_1)

    Call Q32(Ar,Br,C ,A_Multipoles,B_Multipoles,0,Q32_0)
    Call Q32(Ar,Br,C ,A_Multipoles,B_Multipoles,1,Q32_1)

    Call Q50(Ar,Br,C ,A_Multipoles,B_Multipoles,0,Q50_0) 
    Call Q50(Ar,Br,C ,A_Multipoles,B_Multipoles,1,Q50_1)




    
    Approx_6_Energy =   (Q41_0 + Q41_1 + Q50_0 + Q50_1 + Q32_0 + Q32_1)/ R**6

   

    RETURN
    END SUBROUTINE Approx_6_Sph2
    !********************************************************



    SUBROUTINE Q41(Ar,Br,C ,A_Multipoles,B_Multipoles, ind , result)
        IMPLICIT NONE

        real*8, INTENT(INOUT) ::result 
        real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
        real*8 , dimension(3), INTENT(IN):: Ar 
        real*8 , dimension(3), INTENT(IN):: Br
        real*8 , dimension(9), INTENT(IN):: C
        Integer, INTENT(IN):: ind
 
        real*8 , dimension(3):: m
        real*8 , dimension(9):: Phi
        
   
        
        if (ind==0) then
            m =  B_Multipoles(2:4)
            Phi  = A_Multipoles(17:25) 
            
            call  T_ll(Ar,Br,C,Phi,m,4,1 , result)

        else
            m =  A_Multipoles(2:4)
            Phi  = B_Multipoles(17:25)
            
            call  T_ll(Ar,Br,C,m,Phi,1,4 , result)
        end if


    
        RETURN
        END SUBROUTINE Q41


        SUBROUTINE Q32(Ar,Br,C ,A_Multipoles,B_Multipoles, ind , result)
            IMPLICIT NONE
    
     
            real*8, INTENT(INOUT) ::result 
            real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
            real*8 , dimension(3), INTENT(IN):: Ar 
            real*8 , dimension(3), INTENT(IN):: Br
            real*8 , dimension(9), INTENT(IN):: C
            Integer, INTENT(IN):: ind
      
            real*8 , dimension(5):: Qd
            real*8 , dimension(7):: O
            
       
            
           
            if (ind==0) then
                Qd =  B_Multipoles(5:9)
                O  = A_Multipoles(10:16) 
                
                call  T_ll(Ar,Br,C,O,Qd,3,2 , result)
        
                else
                    Qd =  A_Multipoles(5:9)
                    O  = B_Multipoles(10:16) 

                    call  T_ll(Ar,Br,C,Qd,O,2,3 , result)
                end if
            
    
           
    
                
        
            RETURN
            END SUBROUTINE Q32

        SUBROUTINE Q50(Ar,Br,C ,A_Multipoles,B_Multipoles, ind , result)
            IMPLICIT NONE
    
 
            real*8, INTENT(INOUT) ::result 
            real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
            real*8 , dimension(3), INTENT(IN):: Ar 
            real*8 , dimension(3), INTENT(IN):: Br
            real*8 , dimension(9), INTENT(IN):: C
            Integer, INTENT(IN):: ind
            real*8 :: q
            real*8 , dimension(11):: M5
            
       
            
            if (ind==0) then
            q =  B_Multipoles(1)
            M5  = A_Multipoles(26:36) 
            
    
            else
                q =  A_Multipoles(1)
                M5  = B_Multipoles(26:36)
            end if
    

            Call T_l0(Ar,Br,C,q,M5,ind,5 , result) 

            RETURN
            END SUBROUTINE Q50

    
