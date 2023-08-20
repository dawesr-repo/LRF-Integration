
    !********************************************************
SUBROUTINE Approx_5_Sph2(cal_coord,Ar,Br,C ,A_Multipoles,B_Multipoles , Approx_5_Energy)
    IMPLICIT NONE
    
    !   NEED TO DECLARE ALL THE SUBROUTINE ARGUMENTS and
    !   ANY OTHER VARIABLES LOCAL TO THE SUBROUTINE

    real*8, INTENT(INOUT)  :: Approx_5_Energy
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 :: R 
    real*8 , dimension(11), INTENT(IN):: cal_coord
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    real*8 :: Q31_0,Q31_1,Q40_0,Q40_1,Q22_0

    R =cal_coord(1);
    
    Call Q31(Ar,Br,C ,A_Multipoles,B_Multipoles,0,Q31_0) 
    Call Q31(Ar,Br,C ,A_Multipoles,B_Multipoles,1,Q31_1)

    Call Q22(Ar,Br,C ,A_Multipoles,B_Multipoles,0,Q22_0)

    Call Q40(Ar,Br,C ,A_Multipoles,B_Multipoles,0,Q40_0) 
    Call Q40(Ar,Br,C ,A_Multipoles,B_Multipoles,1,Q40_1)

    
    Approx_5_Energy =   (Q31_0 + Q31_1+ Q22_0 +Q40_0+Q40_1)/ R**5

    RETURN
END SUBROUTINE Approx_5_Sph2
    !********************************************************



SUBROUTINE Q31(Ar,Br,C ,A_Multipoles,B_Multipoles, ind , result)
    IMPLICIT NONE


    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    Integer, INTENT(IN):: ind
    real*8 , dimension(3):: m
    real*8 , dimension(7):: O
    

    
    if (ind==0) then

        m =  B_Multipoles(2:4)
        O  = A_Multipoles(10:16) 
        call  T_ll(Ar,Br,C,O,m,3,1 , result)

    else
        m =  A_Multipoles(2:4)
        O  = B_Multipoles(10:16) 

        call  T_ll(Ar,Br,C,m,O,1,3 , result)
    end if



        RETURN
END SUBROUTINE Q31


SUBROUTINE Q22(Ar,Br,C ,A_Multipoles,B_Multipoles, ind , result)
        IMPLICIT NONE


        real*8, INTENT(INOUT) ::result 
        real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
        real*8 , dimension(3), INTENT(IN):: Ar 
        real*8 , dimension(3), INTENT(IN):: Br
        real*8 , dimension(9), INTENT(IN):: C
        Integer, INTENT(IN):: ind

        real*8 , dimension(5):: QdA,QdB

        
    
        
        
        QdB =  B_Multipoles(5:9)
        QdA  = A_Multipoles(5:9) 
        

        call  T_ll(Ar,Br,C,QdA,QdB,2,2 , result)
        
        RETURN
END SUBROUTINE Q22

SUBROUTINE Q40(Ar,Br,C ,A_Multipoles,B_Multipoles, ind , result)
        IMPLICIT NONE

        
        real*8, INTENT(INOUT) ::result 
        real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
        real*8 , dimension(3), INTENT(IN):: Ar 
        real*8 , dimension(3), INTENT(IN):: Br
        real*8 , dimension(9), INTENT(IN):: C
        Integer, INTENT(IN):: ind
        real*8 :: q
        real*8 , dimension(9):: Phi
    
        
    
        
        if (ind==0) then
        q =  B_Multipoles(1)
        Phi  = A_Multipoles(17:25) 
        

        else
            q =  A_Multipoles(1)
            Phi  = B_Multipoles(17:25) 
        end if

        Call T_l0(Ar,Br,C,q,Phi,ind,4 , result)    
    
        RETURN
END SUBROUTINE Q40

    
