
    !********************************************************
SUBROUTINE T_ll(Ar,Br,C,M1,M2,lk1,lk2 , result)
    IMPLICIT NONE

    INTEGER ::  i,j
    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    Integer, INTENT(IN):: lk1,lk2
    real*8 :: t_k1
    real*8 , dimension(2*lk1+1):: M1
    real*8 , dimension(2*lk2+1):: M2
    real*8:: eps=EPSILON(result)
    Character :: cp1,cp2


        result = 0d0 

        do i = 1,2*lk1+1
                do j = 1,2*lk2+1

                    if (DABS(M1(i))>eps .and. DABS(M2(j))>eps) Then

                        call Get_Comp(i,cp1)
                        call Get_Comp(j,cp2)

                        Call T_lk(Ar,Br,C,lk1,Floor((i*1d0)/2d0),cp1,lk2,Floor((j*1d0)/2d0),cp2,t_k1)
                        result = result + M1(i)*M2(j)*t_k1 
                    end if 
                
                end do
        end do

    RETURN
END SUBROUTINE T_ll


SUBROUTINE Get_Comp(i,cp)
            !! DECLARACTIONS
            IMPLICIT NONE

            INTEGER, INTENT(IN):: i
            Character, INTENT(INOUT):: cp
            !! EXECUTABLES
            if (i==1) then 
                cp = "0"
            elseif(mod(i,2)==0)then
                cp = "c"
            else
                cp = "s"
            end if
        
            
        RETURN	!!ONLY NEEDED IF WE PLAN TO REACH END FUNCTION ALL OF THE TIME
END SUBROUTINE Get_Comp
    