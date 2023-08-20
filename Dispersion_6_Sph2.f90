
!********************************************************
SUBROUTINE Dispersion_6_Sph2(cal_coord,Ar,Br,C ,Disp_AB, Disp_6_Energy)
    IMPLICIT NONE

    !   NEED TO DECLARE ALL THE SUBROUTINE ARGUMENTS and
    !   ANY OTHER VARIABLES LOCAL TO THE SUBROUTINE

    real*8, INTENT(INOUT)  :: Disp_6_Energy
    real*8 , dimension(216) , INTENT(IN) :: Disp_AB

    real*8 :: R 
    real*8 , dimension(11), INTENT(IN):: cal_coord
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    real*8 :: D_mm_mm

    R =cal_coord(1)

    Call mm_mm_Dispersion(Ar,Br,C ,Disp_AB,D_mm_mm) 




    Disp_6_Energy =   (-1.0d0)*(D_mm_mm)/ R**6

    RETURN
END SUBROUTINE Dispersion_6_Sph2
!********************************************************



SUBROUTINE mm_mm_Dispersion(Ar,Br,C ,Disp_AB, result)
    IMPLICIT NONE

    INTEGER ::  i,j,k,h,g
    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(216) , INTENT(IN) :: Disp_AB
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    real*8 :: r1,r2
    real*8 , dimension(36)::  mm_mm_AB
    real*8:: eps=EPSILON(result)
    Character :: cp_h,cp_i,cp_j,cp_g
        
    result = 0d0  
    
    mm_mm_AB = Disp_AB(1:36)
    

    
        result = 0d0  

        do g = 1,3
            do h = 1,3       
                do i = 1,3
                    do j = 1,3

                    

                        Call GetIndex_mm_mm(g,h,i,j,k)

                        

                        if (DABS(mm_mm_AB(k))>eps) Then
                        
                                ! Call Tmm(Ar,Br,C ,h,j,0,t_jh)
                                
                                ! Call Tmm(Ar,Br,C ,g,i,0,t_ig)


                                ! result = result + mm_mm_AB(k)*t_ig*t_jh

                            call Get_Comp(g,cp_g)
                            call Get_Comp(h,cp_h)
                            call Get_Comp(i,cp_i)
                            call Get_Comp(j,cp_j)

                            Call T_lk(Ar,Br,C,1,Floor((h*1d0)/2d0),cp_h,1,Floor((j*1d0)/2d0),cp_j,r1)

                            Call T_lk(Ar,Br,C,1,Floor((g*1d0)/2d0),cp_g,1,Floor((i*1d0)/2d0),cp_i,r2)
                            
                            
                            result = result + mm_mm_AB(k)*r1*r2

                        end if 
                    end do
                end do
            end do
        end do

        
    RETURN
END SUBROUTINE mm_mm_Dispersion


  