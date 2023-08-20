
    !********************************************************
SUBROUTINE HyperPolarizability_6_Sph2(cal_coord,Ar,Br,C ,A_Multipoles,B_Multipoles ,A_HPol,B_HPol, H_6_Energy)
    IMPLICIT NONE
    
    !   NEED TO DECLARE ALL THE SUBROUTINE ARGUMENTS and
    !   ANY OTHER VARIABLES LOCAL TO THE SUBROUTINE

    real*8, INTENT(INOUT)  :: H_6_Energy
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(40) , INTENT(IN) :: A_HPol,B_HPol
    real*8 :: R 
    real*8 , dimension(11), INTENT(IN):: cal_coord
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    real*8 :: qqq_mmm_0 , qqq_mmm_1

    R =cal_coord(1)
    
    Call qqq_mmm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_HPol,B_HPol,0,qqq_mmm_0) 
    Call qqq_mmm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_HPol,B_HPol,1,qqq_mmm_1)


    
    H_6_Energy =   (1.0d0/6.0d0)*(qqq_mmm_0 + qqq_mmm_1)/ R**6

    RETURN
END SUBROUTINE HyperPolarizability_6_Sph2
    !********************************************************



SUBROUTINE qqq_mmm(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_HPol,B_HPol, ind , result)
    IMPLICIT NONE

    INTEGER ::  i,j,h,k
    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(40) , INTENT(IN) :: A_HPol,B_HPol
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    Integer, INTENT(IN):: ind
    real*8 :: q,r1,r2,r3
    real*8 , dimension(10):: beta
    real*8:: eps=EPSILON(result)
    Character :: cp_h,cp_i,cp_j,cp_g
        
    result = 0d0 
    
    
    
    
    if (ind==0) then
        q =  B_Multipoles(1)
        beta  = A_HPol(1:10) 

        do h = 1,3
            do i = 1,3
                do j = 1,3

                    Call GetIndex_mmm(h,i,j,k)

                    if (DABS(q)>eps .and.  DABS(beta(k))>eps) Then
                            ! Call Tmq(Ar,Br,C ,i,ind,t_i)
                            ! Call Tmq(Ar,Br,C ,j,ind,t_j)
                            ! Call Tmq(Ar,Br,C ,h,ind,t_h)
                            ! result = result + (q**3)*beta(k)*t_i*t_j*t_h

                        
                        call Get_Comp(h,cp_h)
                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp_i,0,0,"0",r1)

                        Call T_lk(Ar,Br,C,1,Floor((j*1d0)/2d0),cp_j,0,0,"0",r2)
                        
                        Call T_lk(Ar,Br,C,1,Floor((h*1d0)/2d0),cp_h,0,0,"0",r3)
                        
                        result = result + (q**3)*beta(k)*r1*r2*r3



                    end if 
                end do
            end do
        end do
    

    else
        q      =  A_Multipoles(1)
        beta  =  B_HPol(1:10) 

        do h = 1,3
            do i = 1,3
                do j = 1,3

                    Call GetIndex_mmm(h,i,j,k)

                    if (DABS(q)>eps .and.  DABS(beta(k))>eps) Then
                            ! Call Tmq(Ar,Br,C ,i,ind,t_i)
                            ! Call Tmq(Ar,Br,C ,j,ind,t_j)
                            ! Call Tmq(Ar,Br,C ,h,ind,t_h)
                            ! result = result + (q**3)*beta(k)*t_i*t_j*t_h

                        
                        call Get_Comp(h,cp_h)
                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(Ar,Br,C,0,0,"0",1,Floor((i*1d0)/2d0),cp_i,r1)

                        Call T_lk(Ar,Br,C,0,0,"0",1,Floor((j*1d0)/2d0),cp_j,r2)
                        
                        Call T_lk(Ar,Br,C,0,0,"0",1,Floor((h*1d0)/2d0),cp_h,r3)
                        
                        result = result + (q**3)*beta(k)*r1*r2*r3



                    end if 
                end do
            end do
        end do
    end if

    

        
    
    RETURN
END SUBROUTINE qqq_mmm


  