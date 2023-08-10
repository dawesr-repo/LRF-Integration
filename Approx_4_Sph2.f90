
    !********************************************************
    SUBROUTINE Approx_4_Sph2(cal_coord,Ar,Br,C ,A_Multipoles,B_Multipoles , Approx_4_Energy)
    IMPLICIT NONE
    
    !   NEED TO DECLARE ALL THE SUBROUTINE ARGUMENTS and
    !   ANY OTHER VARIABLES LOCAL TO THE SUBROUTINE

    real*8, INTENT(INOUT)  :: Approx_4_Energy
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 :: R 
    real*8 , dimension(11), INTENT(IN):: cal_coord
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    real*8 :: Q21_0,Q21_1,Q30_0,Q30_1

    R =cal_coord(1);
    
    Call Q21(Ar,Br,C ,A_Multipoles,B_Multipoles,0,Q21_0) 
    Call Q21(Ar,Br,C ,A_Multipoles,B_Multipoles,1,Q21_1)

    Call Q30(Ar,Br,C ,A_Multipoles,B_Multipoles,0,Q30_0) 
    Call Q30(Ar,Br,C ,A_Multipoles,B_Multipoles,1,Q30_1)


    
    Approx_4_Energy =   (Q21_0 + Q21_1 +Q30_0+Q30_1)/ R**4

    RETURN
    END SUBROUTINE Approx_4_Sph2
    !********************************************************



    SUBROUTINE Q21(Ar,Br,C ,A_Multipoles,B_Multipoles, ind , result)
        IMPLICIT NONE

   
        real*8, INTENT(INOUT) ::result 
        real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
        real*8 , dimension(3), INTENT(IN):: Ar 
        real*8 , dimension(3), INTENT(IN):: Br
        real*8 , dimension(9), INTENT(IN):: C
        Integer, INTENT(IN):: ind
       
        real*8 , dimension(3):: m
        real*8 , dimension(5):: Qd
        
   
        
        if (ind==0) then
            m =  B_Multipoles(2:4)
            Qd  = A_Multipoles(5:9) 
            
            
            Call T_ll(Ar,Br,C,Qd,m,2,1 , result)

        else
            m =  A_Multipoles(2:4)
            Qd  = B_Multipoles(5:9) 

            Call T_ll(Ar,Br,C,m,Qd,1,2 , result)
        end if

        
    
        RETURN
        END SUBROUTINE Q21


        SUBROUTINE Q30(Ar,Br,C ,A_Multipoles,B_Multipoles, ind , result)
            IMPLICIT NONE
    
           
            real*8, INTENT(INOUT) ::result 
            real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
            real*8 , dimension(3), INTENT(IN):: Ar 
            real*8 , dimension(3), INTENT(IN):: Br
            real*8 , dimension(9), INTENT(IN):: C
            Integer, INTENT(IN):: ind
            real*8 :: q
            real*8 , dimension(7):: O
            
       
            
            if (ind==0) then
            q =  B_Multipoles(1)
            O  = A_Multipoles(10:16) 
            
    
            else
                q =  A_Multipoles(1)
                O  = B_Multipoles(10:16) 
            end if
    
            Call T_l0(Ar,Br,C,q,O,ind,3 , result)
        
            RETURN
            END SUBROUTINE Q30

    
