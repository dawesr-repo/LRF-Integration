
    !********************************************************
    SUBROUTINE Dispersion_7_Sph2(cal_coord,Ar,Br,C ,Disp_AB, Disp_7_Energy)
    IMPLICIT NONE
    
    !   NEED TO DECLARE ALL THE SUBROUTINE ARGUMENTS and
    !   ANY OTHER VARIABLES LOCAL TO THE SUBROUTINE

    real*8, INTENT(INOUT)  :: Disp_7_Energy
    real*8 , dimension(216) , INTENT(IN) :: Disp_AB

    real*8 :: R 
    real*8 , dimension(11), INTENT(IN):: cal_coord
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    real*8 :: D_mQd_mm,D_mm_mQd

    R =cal_coord(1)

    Call mQd_mm_Dispersion(Ar,Br,C ,Disp_AB,D_mQd_mm) 
    Call mm_mQd_Dispersion(Ar,Br,C ,Disp_AB,D_mm_mQd) 


    
    Disp_7_Energy =   (-1.0d0)*(D_mQd_mm+D_mm_mQd)/ R**7


    RETURN
    END SUBROUTINE Dispersion_7_Sph2
    !********************************************************



    SUBROUTINE mQd_mm_Dispersion(Ar,Br,C ,Disp_AB, result)
        IMPLICIT NONE

        INTEGER ::  i,j,k,h,g
        real*8, INTENT(INOUT) ::result 
        real*8 , dimension(216) , INTENT(IN) :: Disp_AB
        real*8 , dimension(3), INTENT(IN):: Ar 
        real*8 , dimension(3), INTENT(IN):: Br
        real*8 , dimension(9), INTENT(IN):: C
        real*8 :: r1,r2
        real*8 , dimension(90)::  mQd_mm_AB
        real*8:: eps=EPSILON(result)
        Character :: cp_h,cp_i,cp_j,cp_g
            
        result = 0d0  
        
        mQd_mm_AB = Disp_AB(37:126)
      

        
       
       

            do g = 1,3
                do h = 1,5       
                    do i = 1,3
                        do j = 1,3

                       

                            Call GetIndex_mQd_mm(g,h,i,j,k)

                          
                            

                            if (DABS(mQd_mm_AB(k))>eps) Then
                            
                                    ! Call TQdm(Ar,Br,C ,h,j,0,t_jh)
                                    
                                    ! Call Tmm(Ar,Br,C ,g,i,0,t_ig)


                                    ! result = result + mQd_mm_AB(k)*t_ig*t_jh

                                call Get_Comp(g,cp_g)
                                call Get_Comp(h,cp_h)
                                call Get_Comp(i,cp_i)
                                call Get_Comp(j,cp_j)

                                Call T_lk(Ar,Br,C,1,Floor((g*1d0)/2d0),cp_g,1,Floor((i*1d0)/2d0),cp_i,r1)

                                Call T_lk(Ar,Br,C,2,Floor((h*1d0)/2d0),cp_h,1,Floor((j*1d0)/2d0),cp_j,r2)
                                
                                
                                result = result + mQd_mm_AB(k)*r1*r2

                            end if 
                        end do
                    end do
                end do
            end do
        RETURN
    END SUBROUTINE mQd_mm_Dispersion

    SUBROUTINE mm_mQd_Dispersion(Ar,Br,C ,Disp_AB, result)
        IMPLICIT NONE

        INTEGER ::  i,j,k,h,g
        real*8, INTENT(INOUT) ::result 
        real*8 , dimension(216) , INTENT(IN) :: Disp_AB
        real*8 , dimension(3), INTENT(IN):: Ar 
        real*8 , dimension(3), INTENT(IN):: Br
        real*8 , dimension(9), INTENT(IN):: C
        real*8 :: r1,r2
        real*8 , dimension(90)::  mm_mQd_AB
        real*8:: eps=EPSILON(result)
        Character :: cp_h,cp_i,cp_j,cp_g
            
        result = 0d0  
        
        mm_mQd_AB = Disp_AB(127:216)
      

        

   

            do g = 1,3
                do h = 1,3       
                    do i = 1,3
                        do j = 1,5

                       

                            Call GetIndex_mm_mQd(g,h,i,j,k)

                          

                            if (DABS(mm_mQd_AB(k))>eps) Then
                            
                                    ! Call TQdm(Ar,Br,C ,j,h,1,t_jh)
                                    
                                    ! Call Tmm(Ar,Br,C ,g,i,0,t_ig)


                                    ! result = result + mm_mQd_AB(k)*t_ig*t_jh

                                call Get_Comp(g,cp_g)
                                call Get_Comp(h,cp_h)
                                call Get_Comp(i,cp_i)
                                call Get_Comp(j,cp_j)

                                Call T_lk(Ar,Br,C,1,Floor((g*1d0)/2d0),cp_g,1,Floor((i*1d0)/2d0),cp_i,r1)

                                Call T_lk(Ar,Br,C,1,Floor((h*1d0)/2d0),cp_h,2,Floor((j*1d0)/2d0),cp_j,r2)
                                
             
                               
                                result = result + mm_mQd_AB(k)*r1*r2

                            end if 
                        end do
                    end do
                end do
            end do
        RETURN
    END SUBROUTINE mm_mQd_Dispersion
  