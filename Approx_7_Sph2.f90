
    !********************************************************
    SUBROUTINE Approx_7_Sph2(cal_coord,Ar,Br,C ,A_Multipoles,B_Multipoles , Approx_7_Energy)
    IMPLICIT NONE
    
    !   NEED TO DECLARE ALL THE SUBROUTINE ARGUMENTS and
    !   ANY OTHER VARIABLES LOCAL TO THE SUBROUTINE

    real*8, INTENT(INOUT)  :: Approx_7_Energy
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 :: R 
    real*8 , dimension(11), INTENT(IN):: cal_coord
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    real*8 :: Q51_0,Q51_1,Q60_0,Q60_1,Q42_0,Q42_1,Q33_0

    R =cal_coord(1);

    Call Q33(Ar,Br,C ,A_Multipoles,B_Multipoles,0,Q33_0)

    Call Q51(Ar,Br,C ,A_Multipoles,B_Multipoles,0,Q51_0) 
    Call Q51(Ar,Br,C ,A_Multipoles,B_Multipoles,1,Q51_1)

    Call Q42(Ar,Br,C ,A_Multipoles,B_Multipoles,0,Q42_0)
    Call Q42(Ar,Br,C ,A_Multipoles,B_Multipoles,1,Q42_1)

    Call Q60(Ar,Br,C ,A_Multipoles,B_Multipoles,0,Q60_0) 
    Call Q60(Ar,Br,C ,A_Multipoles,B_Multipoles,1,Q60_1)


    
    Approx_7_Energy =   (Q60_0+Q60_1+Q51_0+Q51_1+Q42_0+Q42_1+ Q33_0)/ R**7

   

    RETURN
    END SUBROUTINE Approx_7_Sph2
    !********************************************************

    SUBROUTINE Q33(Ar,Br,C ,A_Multipoles,B_Multipoles, ind , result)
        IMPLICIT NONE

        real*8, INTENT(INOUT) ::result 
        real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
        real*8 , dimension(3), INTENT(IN):: Ar 
        real*8 , dimension(3), INTENT(IN):: Br
        real*8 , dimension(9), INTENT(IN):: C
        Integer, INTENT(IN):: ind
  
        real*8 , dimension(7):: OA
        real*8 , dimension(7):: OB
        
   
        
   
        OB =  B_Multipoles(10:16)
        OA  = A_Multipoles(10:16)
        

     
        call  T_ll(Ar,Br,C,OA,OB,3,3 , result)
          
    
        RETURN
        END SUBROUTINE Q33


    SUBROUTINE Q51(Ar,Br,C ,A_Multipoles,B_Multipoles, ind , result)
        IMPLICIT NONE

  
        real*8, INTENT(INOUT) ::result 
        real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
        real*8 , dimension(3), INTENT(IN):: Ar 
        real*8 , dimension(3), INTENT(IN):: Br
        real*8 , dimension(9), INTENT(IN):: C
        Integer, INTENT(IN):: ind
      
        real*8 , dimension(3):: m
        real*8 , dimension(11):: M5

        
   
        
        if (ind==0) then
        m =  B_Multipoles(2:4)
        M5  = A_Multipoles(26:36) 
        
        call  T_ll(Ar,Br,C,M5,m,5,1 , result)

        else
            m =  A_Multipoles(2:4)
            M5  = B_Multipoles(26:36)

            call  T_ll(Ar,Br,C,m,M5,1,5 , result)
        end if

          
    
        RETURN
END SUBROUTINE Q51


SUBROUTINE Q42(Ar,Br,C ,A_Multipoles,B_Multipoles, ind , result)
            IMPLICIT NONE
    
            
            real*8, INTENT(INOUT) ::result 
            real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
            real*8 , dimension(3), INTENT(IN):: Ar 
            real*8 , dimension(3), INTENT(IN):: Br
            real*8 , dimension(9), INTENT(IN):: C
            Integer, INTENT(IN):: ind
            
            real*8 , dimension(5):: Qd
            real*8 , dimension(9):: Phi
            
       
            
           
            if (ind==0) then
                    Qd =  B_Multipoles(5:9)
                    Phi  = A_Multipoles(17:25) 
                
                    call  T_ll(Ar,Br,C,Phi,Qd,4,2 , result)
        
                else
                    Qd =  A_Multipoles(5:9)
                    Phi  = B_Multipoles(17:25) 

                    call  T_ll(Ar,Br,C,Qd,Phi,2,4 , result)
                end if
            
    
           
        
            RETURN
END SUBROUTINE Q42

SUBROUTINE Q60(Ar,Br,C ,A_Multipoles,B_Multipoles, ind , result)
            IMPLICIT NONE
    

            real*8, INTENT(INOUT) ::result 
            real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
            real*8 , dimension(3), INTENT(IN):: Ar 
            real*8 , dimension(3), INTENT(IN):: Br
            real*8 , dimension(9), INTENT(IN):: C
            Integer, INTENT(IN):: ind
            real*8 :: q
            real*8 , dimension(13):: M6
            
       
            
            if (ind==0) then
            q =  B_Multipoles(1)
            M6  = A_Multipoles(37:49) 
            
    
            else
                q =  A_Multipoles(1)
                M6  = B_Multipoles(37:49)
            end if
    
            Call T_l0(Ar,Br,C,q,M6,ind,6 , result) 
        
            RETURN
END SUBROUTINE Q60

    
