
!********************************************************
SUBROUTINE Approx_2_Sph2(cal_coord,Ar,Br,C ,A_Multipoles,B_Multipoles , Approx_2_Energy)
    IMPLICIT NONE

    !   NEED TO DECLARE ALL THE SUBROUTINE ARGUMENTS and
    !   ANY OTHER VARIABLES LOCAL TO THE SUBROUTINE

    real*8, INTENT(INOUT)  :: Approx_2_Energy
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 :: R 
    real*8 , dimension(11), INTENT(IN):: cal_coord
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    real*8 :: Q10_0,Q10_1

    R =cal_coord(1);

    Call Q10(Ar,Br,C ,A_Multipoles,B_Multipoles,0,Q10_0) 
    Call Q10(Ar,Br,C ,A_Multipoles,B_Multipoles,1,Q10_1)



    Approx_2_Energy =   (Q10_0 + Q10_1)/ R**2


    RETURN
END SUBROUTINE Approx_2_Sph2
!********************************************************



SUBROUTINE Q10(Ar,Br,C ,A_Multipoles,B_Multipoles, ind , result)
    IMPLICIT NONE

    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    Integer, INTENT(IN):: ind
    real*8 :: q
    real*8 , dimension(3):: m
    
    
    
    if (ind==0) then
    q =  B_Multipoles(1)
    m  = A_Multipoles(2:4) 
    

    else
        q =  A_Multipoles(1)
        m  = B_Multipoles(2:4) 
    end if

    Call T_l0(Ar,Br,C,q,m,ind,1 , result)

        
    RETURN
END SUBROUTINE Q10


        