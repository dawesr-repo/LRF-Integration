
!********************************************************
SUBROUTINE Approx_1_Sph2(cal_coord,Ar,Br,C ,A_Multipoles,B_Multipoles , Approx_1_Energy)
    IMPLICIT NONE

    !   NEED TO DECLARE ALL THE SUBROUTINE ARGUMENTS and
    !   ANY OTHER VARIABLES LOCAL TO THE SUBROUTINE

    real*8, INTENT(INOUT)  :: Approx_1_Energy
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 :: qA,qB
    real*8 :: R 
    real*8 , dimension(11), INTENT(IN):: cal_coord
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C

    R =cal_coord(1);


    qA=A_Multipoles(1)
    qB=B_Multipoles(1)


    Approx_1_Energy =   qA*qB / R

    RETURN
END SUBROUTINE Approx_1_Sph2
!********************************************************