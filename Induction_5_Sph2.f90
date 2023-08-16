
    !********************************************************
    SUBROUTINE Induction_5_Sph2(cal_coord,Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, Ind_5_Energy)
    IMPLICIT NONE
    
    !   NEED TO DECLARE ALL THE SUBROUTINE ARGUMENTS and
    !   ANY OTHER VARIABLES LOCAL TO THE SUBROUTINE

    real*8, INTENT(INOUT)  :: Ind_5_Energy
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol
    real*8 :: R 
    real*8 , dimension(11), INTENT(IN):: cal_coord
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    real*8 :: qm_mm_0,qm_mm_1,  qq_mQd_0,qq_mQd_1

    R =cal_coord(1)
    
    Call qm_mm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qm_mm_0) 
    Call qm_mm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qm_mm_1)


    Call qq_mQd(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qq_mQd_0) 
    Call qq_mQd(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qq_mQd_1)

!write(*,*)qm_mm_0,qm_mm_1,qq_mQd_0,qq_mQd_1

    
    Ind_5_Energy =   (-0.5d0)*(2.0d0*qm_mm_0 + 2.0d0*qm_mm_1 +2.0d0*qq_mQd_0+ 2.0d0*qq_mQd_1)/ R**5

    RETURN
    END SUBROUTINE Induction_5_Sph2
    !********************************************************



    SUBROUTINE qm_mm(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
        IMPLICIT NONE

        INTEGER ::  i,j,k,h
        real*8, INTENT(INOUT) ::result 
        real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
        real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol
        real*8 , dimension(3), INTENT(IN):: Ar 
        real*8 , dimension(3), INTENT(IN):: Br
        real*8 , dimension(9), INTENT(IN):: C
        Integer, INTENT(IN):: ind
        real*8 :: q,r1,r2
        real*8 , dimension(3):: m
        real*8 , dimension(6):: alpha
        real*8:: eps=EPSILON(result)

        Character :: cp_h,cp_i,cp_j
        
        
        
        if (ind==0) then
            q =  B_Multipoles(1)
            m =  B_Multipoles(2:4)
            alpha  = A_Pol(1:6) 
        
            result = 0d0  

          
            do h=1,3
                do i = 1,3
                    do j = 1,3

                        

                        if (DABS(q)>eps .and.    DABS(m(h))>eps) Then
                            Call GetIndex_mm(i,j,k)
                            if(DABS(alpha(k))>eps)then

                                call Get_Comp(h,cp_h)
                                call Get_Comp(i,cp_i)
                                call Get_Comp(j,cp_j)

                                Call T_lk(Ar,Br,C,1,Floor((j*1d0)/2d0),cp_j,1,Floor((h*1d0)/2d0),cp_h,r1)

                                Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp_i,0,0,"0",r2)
                                
                               
                                result = result + (q*m(h))*alpha(k)*r1*r2

                            end if 
                        end if 
                    end do
                end do
            end do

        else
            q       =  A_Multipoles(1)
            m       =  A_Multipoles(2:4)
            alpha   =  B_Pol(1:6) 

            result = 0d0  

          
            do h=1,3
                do i = 1,3
                    do j = 1,3

                        

                        if (DABS(q)>eps .and. DABS(m(h))>eps) Then
                            Call GetIndex_mm(i,j,k)
                            if(DABS(alpha(k))>eps)then
                                !Call Tmm(Ar,Br,C ,j,h,ind,t_jh)
                                !Call Tmq(Ar,Br,C ,i,ind,t_i)
                                !result = result + (q*m(h))*alpha(k)*t_i*t_jh

                                call Get_Comp(h,cp_h)
                                call Get_Comp(i,cp_i)
                                call Get_Comp(j,cp_j)

                                Call T_lk(Ar,Br,C,1,Floor((h*1d0)/2d0),cp_h,1,Floor((j*1d0)/2d0),cp_j,r1)

                                Call T_lk(Ar,Br,C,0,0,"0",1,Floor((i*1d0)/2d0),cp_i,r2)
                                
                               
                                result = result + (q*m(h))*alpha(k)*r1*r2

                            end if 
                        end if 
                    end do
                end do
            end do
        end if


    
        RETURN
        END SUBROUTINE qm_mm





        SUBROUTINE qq_mQd(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
            IMPLICIT NONE
    
            INTEGER ::  i,j,k
            real*8, INTENT(INOUT) ::result 
            real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
            real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol
            real*8 , dimension(3), INTENT(IN):: Ar 
            real*8 , dimension(3), INTENT(IN):: Br
            real*8 , dimension(9), INTENT(IN):: C
            Integer, INTENT(IN):: ind
            real*8 :: q,r1,r2
            real*8 , dimension(15):: alpha_mQ
            real*8:: eps=EPSILON(result)
            Character :: cp1,cp2
            
            
            
            
            if (ind==0) then
                q =  B_Multipoles(1)
                alpha_mQ  = A_Pol(7:21) 
            
                result = 0d0  
    
    
                    do i = 1,3
                        do j = 1,5
    
                            if (DABS(q)>eps    ) Then
                                Call GetIndex_mQd(i,j,k)
                                    !Call TQdq(Ar,Br,C ,j,ind,t_j)
                                    !Call Tmq (Ar,Br,C ,i,ind,t_i)
                                if(DABS(alpha_mQ(k))>eps)Then
                                    call Get_Comp(i,cp1)
                                    call Get_Comp(j,cp2)

                                    Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp1,0,0,"0",r1)
                                    Call T_lk(Ar,Br,C,2,Floor((j*1d0)/2d0),cp2,0,0,"0",r2)
                    
                                    result = result + (q**2)*alpha_mQ(k)*r1*r2
                                    ! result = result + (q**2)*alpha_mQ(k)*t_i*t_j
                                end if 
                            end if 
                        end do
                    end do
    
            else
                q       =  A_Multipoles(1)
                alpha_mQ   =  B_Pol(7:21) 

                result = 0d0  
    
    
                    do i = 1,3
                        do j = 1,5
    
                            
    
                            if (DABS(q)>eps  ) Then
                                Call GetIndex_mQd(i,j,k)
                                if(DABS(alpha_mQ(k))>eps)then

                                    call Get_Comp(i,cp1)
                                    call Get_Comp(j,cp2)

                                    Call T_lk(Ar,Br,C,0,0,"0",1,Floor((i*1d0)/2d0),cp1,r1)
                                    Call T_lk(Ar,Br,C,0,0,"0",2,Floor((j*1d0)/2d0),cp2,r2)
                            
                                    result = result + (q**2)*alpha_mQ(k)*r1*r2
                                    ! result = result + (q**2)*alpha_mQ(k)*t_i*t_j
                                end if 
                            end if
                        end do
                    end do
                    
            end if
    
                    

        
            RETURN
            END SUBROUTINE qq_mQd