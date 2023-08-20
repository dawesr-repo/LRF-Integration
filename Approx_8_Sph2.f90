
    !********************************************************
SUBROUTINE Approx_8_Sph2(cal_coord,Ar,Br,C ,A_Multipoles,B_Multipoles , Approx_8_Energy)
    IMPLICIT NONE
    
    !   NEED TO DECLARE ALL THE SUBROUTINE ARGUMENTS and
    !   ANY OTHER VARIABLES LOCAL TO THE SUBROUTINE

    real*8, INTENT(INOUT)  :: Approx_8_Energy
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 :: R 
    real*8 , dimension(11), INTENT(IN):: cal_coord
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    real*8 :: Q61_0,Q61_1,Q70_0,Q70_1,Q52_0,Q52_1,Q43_0,Q43_1

    R =cal_coord(1);

    Call Q43(Ar,Br,C ,A_Multipoles,B_Multipoles,0,Q43_0)
    Call Q43(Ar,Br,C ,A_Multipoles,B_Multipoles,1,Q43_1)

    Call Q52(Ar,Br,C ,A_Multipoles,B_Multipoles,0,Q61_0) 
    Call Q52(Ar,Br,C ,A_Multipoles,B_Multipoles,1,Q61_1)

    Call Q70(Ar,Br,C ,A_Multipoles,B_Multipoles,0,Q52_0)
    Call Q70(Ar,Br,C ,A_Multipoles,B_Multipoles,1,Q52_1)

    Call Q61(Ar,Br,C ,A_Multipoles,B_Multipoles,0,Q70_0) 
    Call Q61(Ar,Br,C ,A_Multipoles,B_Multipoles,1,Q70_1)


    
    Approx_8_Energy =   (Q70_0+Q70_1+Q61_0+Q61_1+Q52_0+Q52_1+ Q43_0+Q43_1)/ R**8

   

    RETURN
END SUBROUTINE Approx_8_Sph2
    !********************************************************

SUBROUTINE Q43(Ar,Br,C ,A_Multipoles,B_Multipoles, ind , result)
        IMPLICIT NONE


        real*8, INTENT(INOUT) ::result 
        real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
        real*8 , dimension(3), INTENT(IN):: Ar 
        real*8 , dimension(3), INTENT(IN):: Br
        real*8 , dimension(9), INTENT(IN):: C
        Integer, INTENT(IN):: ind

        real*8 , dimension(7):: O
        real*8 , dimension(9):: Phi

        
   
        
   
        if (ind==0) then
                O    =  B_Multipoles(10:16)
                Phi  =  A_Multipoles(17:25) 
                
                call  T_ll(Ar,Br,C,Phi,O,4,3 , result)
            
            else
                O    =  A_Multipoles(10:16)
                Phi  =  B_Multipoles(17:25) 

                call  T_ll(Ar,Br,C,O,Phi,3,4 , result)
   
            end if
    
        RETURN
END SUBROUTINE Q43


SUBROUTINE Q61(Ar,Br,C ,A_Multipoles,B_Multipoles, ind , result)
        IMPLICIT NONE


        real*8, INTENT(INOUT) ::result 
        real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
        real*8 , dimension(3), INTENT(IN):: Ar 
        real*8 , dimension(3), INTENT(IN):: Br
        real*8 , dimension(9), INTENT(IN):: C
        Integer, INTENT(IN):: ind
      
        real*8 , dimension(3):: m
        real*8 , dimension(13):: M6

        
   
        
        if (ind==0) then
            m =  B_Multipoles(2:4)
            M6  = A_Multipoles(37:49) 
        
            call  T_ll(Ar,Br,C,M6,m,6,1 , result)
        else
            m =  A_Multipoles(2:4)
            M6  = B_Multipoles(37:49)

            call  T_ll(Ar,Br,C,m,M6,1,6 , result)
        end if

       
        RETURN
END SUBROUTINE Q61


SUBROUTINE Q52(Ar,Br,C ,A_Multipoles,B_Multipoles, ind , result)
    IMPLICIT NONE

    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    Integer, INTENT(IN):: ind

    real*8 , dimension(5):: Qd
    real*8 , dimension(11):: M5

    

    
    
    if (ind==0) then
        Qd =  B_Multipoles(5:9)
        M5  = A_Multipoles(26:36) 
    
        call  T_ll(Ar,Br,C,M5,Qd,5,2 , result)
    else
        Qd =  A_Multipoles(5:9)
        M5  = B_Multipoles(26:36)

        call  T_ll(Ar,Br,C,Qd,M5,2,5 , result)
    end if

    RETURN
END SUBROUTINE Q52

SUBROUTINE Q70(Ar,Br,C ,A_Multipoles,B_Multipoles, ind , result)
    IMPLICIT NONE


    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    Integer, INTENT(IN):: ind
    real*8 :: q
    real*8 , dimension(15):: M7
    

    
    if (ind==0) then
        q =  B_Multipoles(1)
        M7  = A_Multipoles(50:64) 

    else
        q =  A_Multipoles(1)
        M7  = B_Multipoles(50:64)
    end if

    Call T_l0(Ar,Br,C,q,M7,ind,7, result) 

    RETURN
END SUBROUTINE Q70

    
