
    !********************************************************
SUBROUTINE Induction_6_Sph2(A_Multipoles,B_Multipoles ,A_Pol,B_Pol, Ind_6_Energy)
    IMPLICIT NONE
    
    !   NEED TO DECLARE ALL THE SUBROUTINE ARGUMENTS and
    !   ANY OTHER VARIABLES LOCAL TO THE SUBROUTINE

    real*8, INTENT(INOUT)  :: Ind_6_Energy
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol
    real*8 :: R 


    real*8 :: qQd_mm_0,qQd_mm_1,mm_mm_0,mm_mm_1,qm_mQd_0,qm_mQd_1,qm_Qdm_0,qm_Qdm_1
    real*8 :: qq_mO_0,qq_mO_1,qq_QdQd_0,qq_QdQd_1

    R =cal_coord(1)
    
    Call qQd_mm(A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qQd_mm_0) 
    Call qQd_mm(A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qQd_mm_1)

    Call mm_mm(A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,mm_mm_0) 
    Call mm_mm(A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,mm_mm_1)

    Call qm_mQd(A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qm_mQd_0) 
    Call qm_mQd(A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qm_mQd_1)

    Call qm_Qdm(A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qm_Qdm_0) 
    Call qm_Qdm(A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qm_Qdm_1)
    
    Call qq_mO(A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qq_mO_0) 
    Call qq_mO(A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qq_mO_1)

    Call qq_QdQd(A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qq_QdQd_0) 
    Call qq_QdQd(A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qq_QdQd_1)


    Ind_6_Energy =   (-0.5d0)*(  2.0d0*qQd_mm_0 + 2.0d0*qQd_mm_1 &
                                 +mm_mm_0 + mm_mm_1 &
                                 +2.0d0*qm_mQd_0 + 2.0d0*qm_mQd_1 +2.0d0*qm_Qdm_0 + 2.0d0*qm_Qdm_1 &
                                +2.0d0*qq_mO_0  + 2.0d0*qq_mO_1 + qq_QdQd_0 + qq_QdQd_1 &
                                )/ R**6



    RETURN
END SUBROUTINE Induction_6_Sph2
    !********************************************************



SUBROUTINE qQd_mm(A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
    IMPLICIT NONE

    INTEGER ::  i,j,k,h
    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol

    Integer, INTENT(IN):: ind
    real*8 :: q,r1,r2
    real*8 , dimension(5):: Qd
    real*8 , dimension(6):: alpha
    real*8:: eps=EPSILON(result)

    Character :: cp_h,cp_i,cp_j
    
    
    
    if (ind==0) then
        q =  B_Multipoles(1)
        Qd =  B_Multipoles(5:9)
        alpha  = A_Pol(1:6) 
    
        result = 0d0  

        
        do h=1,5
            do i = 1,3
                do j = 1,3

                    

                    if (DABS(q)>eps .and.   DABS(Qd(h))>eps) Then
                            !Call TQdm(h,j,mod(ind+1,2),t_jh)
                            !Call Tmq(i,ind,t_i)
                        Call GetIndex_mm(i,j,k)

                        if(DABS(alpha(k))>eps )Then
                            call Get_Comp(h,cp_h)
                            call Get_Comp(i,cp_i)
                            call Get_Comp(j,cp_j)

                            !result = result + (q*Qd(h))*alpha(k)*t_i*t_jh

                            Call T_lk(1,Floor((i*1d0)/2d0),cp_i,0,0,"0",r1)

                            Call T_lk(1,Floor((j*1d0)/2d0),cp_j,2,Floor((h*1d0)/2d0),cp_h,r2)
                            
                            
                            result = result + (q*Qd(h))*alpha(k)*r1*r2
                        end if
                    end if 
                end do
            end do
        end do

    else
        q       =  A_Multipoles(1)
        Qd       =  A_Multipoles(5:9)
        alpha   =  B_Pol(1:6) 

        result = 0d0  

        
        do h=1,5
            do i = 1,3
                do j = 1,3

                    

                    if (DABS(q)>eps   .and. DABS(Qd(h))>eps) Then
                            !Call TQdm(h,j,mod(ind+1,2),t_jh)
                            !Call Tmq(i,ind,t_i)
                        Call GetIndex_mm(i,j,k)
                        if(DABS(alpha(k))>eps)Then
                            call Get_Comp(h,cp_h)
                            call Get_Comp(i,cp_i)
                            call Get_Comp(j,cp_j)

                            Call T_lk(0,0,"0",1,Floor((i*1d0)/2d0),cp_i,r1)

                            Call T_lk(2,Floor((h*1d0)/2d0),cp_h,1,Floor((j*1d0)/2d0),cp_j,r2)
                            
                            
                            result = result + (q*Qd(h))*alpha(k)*r1*r2
                        end if
                    end if  
                end do
            end do
        end do
    end if


    RETURN
END SUBROUTINE qQd_mm




SUBROUTINE mm_mm(A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
    IMPLICIT NONE

    INTEGER ::  i,j,k,h,g
    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br

    real*8 , dimension(3):: m
    real*8 , dimension(6):: alpha
    real*8:: eps=EPSILON(result)

    Character :: cp_h,cp_i,cp_j,cp_g
    
    result = 0d0  
    
    if (ind==0) then
        
        m =  B_Multipoles(2:4)
        alpha  = A_Pol(1:6) 
        
        do g=1,3 
            do h=1,3
                do i = 1,3
                    do j = 1,3

                        

                        if (DABS(m(g))>eps  .and. DABS(m(h))>eps) Then
                            Call GetIndex_mm(i,j,k)
                            if ( DABS(alpha(k))>eps ) Then
                                ! Call Tmm(j,h,ind,t_jh)
                                ! Call Tmm(i,g,ind,t_ig)
                                ! result = result + (m(g)*m(h))*alpha(k)*t_ig*t_jh

                            call Get_Comp(g,cp_g)
                            call Get_Comp(h,cp_h)
                            call Get_Comp(i,cp_i)
                            call Get_Comp(j,cp_j)

                            Call T_lk(1,Floor((j*1d0)/2d0),cp_j,1,Floor((h*1d0)/2d0),cp_h,r1)

                            Call T_lk(1,Floor((i*1d0)/2d0),cp_i,1,Floor((g*1d0)/2d0),cp_g,r2)
                            
                            
                            result = result + (m(g)*m(h))*alpha(k)*r1*r2
                        end if 
                        end if 
                    end do
                end do
            end do
        end do

    else
        m      =  A_Multipoles(2:4)
        alpha   =  B_Pol(1:6) 


        do g=1,3 
            do h=1,3
                do i = 1,3
                    do j = 1,3

                        

                        if (DABS(m(g))>eps .and. DABS(m(h))>eps) Then
                            Call GetIndex_mm(i,j,k)
                            if ( DABS(alpha(k))>eps ) Then
                                ! Call Tmm(j,h,ind,t_jh)
                                ! Call Tmm(i,g,ind,t_ig)
                                ! result = result + (m(g)*m(h))*alpha(k)*t_ig*t_jh

                            call Get_Comp(g,cp_g)
                            call Get_Comp(h,cp_h)
                            call Get_Comp(i,cp_i)
                            call Get_Comp(j,cp_j)

                            Call T_lk(1,Floor((h*1d0)/2d0),cp_h,1,Floor((j*1d0)/2d0),cp_j,r1)

                            Call T_lk(1,Floor((g*1d0)/2d0),cp_g,1,Floor((i*1d0)/2d0),cp_i,r2)
                            
                            
                            result = result + (m(g)*m(h))*alpha(k)*r1*r2
                        end if 
                        end if 
                    end do
                end do
            end do
        end do


    end if

        

    
    RETURN
END SUBROUTINE mm_mm











SUBROUTINE qm_mQd(A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
    IMPLICIT NONE

    INTEGER ::  i,j,k,h
    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol

    Integer, INTENT(IN):: ind

    real*8:: eps=EPSILON(result)

    
    Character :: cp_h,cp_i,cp_j
    
    result = 0d0  

    if (ind==0) then
        q =  B_Multipoles(1)
        m =  B_Multipoles(2:4)
        alpha_mQ  = A_Pol(7:21) 
        
        

        do h = 1,3
            do i = 1,3
                do j = 1,5

                    

                    if (DABS(m(h))>eps .and. DABS(q)>eps ) Then
                        Call GetIndex_mQd(i,j,k)
                        if (  DABS(alpha_mQ(k))>eps ) Then
                            ! Call TQdq(j,ind,t_j)
                            ! Call Tmm (i,h,ind,t_ih)
                            ! result = result + (q*m(h))*alpha_mQ(k)*t_ih*t_j

                    
                        call Get_Comp(h,cp_h)
                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(1,Floor((i*1d0)/2d0),cp_i,1,Floor((h*1d0)/2d0),cp_h,r1)

                        Call T_lk(2,Floor((j*1d0)/2d0),cp_j,0,0,"0",r2)
                        
                        
                        result = result + (q*m(h))*alpha_mQ(k)*r1*r2

                    end if
                    end if 
                end do
            end do
        end do


    else
        q           =  A_Multipoles(1)
        m           =  A_Multipoles(2:4)
        alpha_mQ    =  B_Pol(7:21) 

        

        do h = 1,3
            do i = 1,3
                do j = 1,5


                    if (DABS(m(h))>eps .and. DABS(q)>eps  ) Then
                        Call GetIndex_mQd(i,j,k)
                        if ( DABS(alpha_mQ(k))>eps ) Then
                            ! Call TQdq(j,ind,t_j)
                            ! Call Tmm (i,h,ind,t_ih)
                            ! result = result + (q*m(h))*alpha_mQ(k)*t_ih*t_j

                    
                        call Get_Comp(h,cp_h)
                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(1,Floor((h*1d0)/2d0),cp_h,1,Floor((i*1d0)/2d0),cp_i,r1)

                        Call T_lk(0,0,"0",2,Floor((j*1d0)/2d0),cp_j,r2)
                        
                        
                        result = result + (q*m(h))*alpha_mQ(k)*r1*r2

                    end if 
                    end if 
                end do
            end do
        end do



    end if

        
    RETURN
END SUBROUTINE qm_mQd


SUBROUTINE qm_Qdm(A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
    IMPLICIT NONE

    INTEGER ::  i,j,k,h
    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol

    Integer, INTENT(IN):: ind
    real*8 :: q,r1,r2
    real*8 , dimension(3):: m

    
    Character :: cp_h,cp_i,cp_j
    
    result = 0d0 

    if (ind==0) then
        q =  B_Multipoles(1)
        m =  B_Multipoles(2:4)
        alpha_mQ  = A_Pol(7:21) 

            

        do h = 1,3
            do i = 1,3
                do j = 1,5

                    

                    if (DABS(m(h))>eps .and. DABS(q)>eps ) Then
                        Call GetIndex_mQd(i,j,k)
                        if ( DABS(alpha_mQ(k))>eps ) Then
                            ! Call TQdm(j,h,ind,t_jh)
                            ! Call Tmq (i,ind,t_i)
                            ! result = result + (q*m(h))*alpha_mQ(k)*t_i*t_jh

                        call Get_Comp(h,cp_h)
                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(2,Floor((j*1d0)/2d0),cp_j,1,Floor((h*1d0)/2d0),cp_h,r1)

                        Call T_lk(1,Floor((i*1d0)/2d0),cp_i,0,0,"0",r2)
                        
                        
                        result = result + (q*m(h))*alpha_mQ(k)*r1*r2
                    end if 
                    end if 
                end do
            end do
        end do
    

    else
        q           =  A_Multipoles(1)
        m           =  A_Multipoles(2:4)
        alpha_mQ    =  B_Pol(7:21) 

    

        do h = 1,3
            do i = 1,3
                do j = 1,5

                    

                    if (DABS(m(h))>eps .and. DABS(q)>eps  ) Then
                        Call GetIndex_mQd(i,j,k)
                        if ( DABS(alpha_mQ(k))>eps ) Then
                            ! Call TQdm(j,h,ind,t_jh)
                            ! Call Tmq (i,ind,t_i)
                            ! result = result + (q*m(h))*alpha_mQ(k)*t_i*t_jh

                        call Get_Comp(h,cp_h)
                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(1,Floor((h*1d0)/2d0),cp_h,2,Floor((j*1d0)/2d0),cp_j,r1)

                        Call T_lk(0,0,"0",1,Floor((i*1d0)/2d0),cp_i,r2)
                        
                        
                        result = result + (q*m(h))*alpha_mQ(k)*r1*r2
                    end if 
                    end if 
                end do
            end do
        end do
    end if

        

    RETURN
END SUBROUTINE qm_Qdm


SUBROUTINE qq_mO(A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
    IMPLICIT NONE

    INTEGER ::  i,j,k
    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol

    Integer, INTENT(IN):: ind
    real*8 :: q,r1,r2
    real*8 , dimension(21):: alpha_mO
    real*8:: eps=EPSILON(result)

    Character :: cp_i,cp_j
    
    
    result = 0d0  

    if (ind==0) then
        q =  B_Multipoles(1)
        alpha_mO  = A_Pol(37:57) 
    
        


            do i = 1,3
                do j = 1,7

                    

                    if ( DABS(q)>eps  ) Then
                        Call GetIndex_mO(i,j,k)
                        if (   DABS(alpha_mO(k))>eps ) Then
                            
                            ! Call TOq(j,ind,t_j)
                            ! Call Tmq (i,ind,t_i)
                            ! result = result + (q**2)*alpha_mO(k)*t_i*t_j

                
                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(3,Floor((j*1d0)/2d0),cp_j,0,0,"0",r1)

                        Call T_lk(1,Floor((i*1d0)/2d0),cp_i,0,0,"0",r2)
                        
                        
                        result = result + (q**2)*alpha_mO(k)*r1*r2
                        end if 
                    end if 
                end do
            end do
    

    else
        q           =  A_Multipoles(1)
        alpha_mO    =  B_Pol(37:57) 


        


            do i = 1,3
                do j = 1,7

                    

                    if ( DABS(q)>eps  ) Then
                        Call GetIndex_mO(i,j,k)
                        if (  DABS(alpha_mO(k))>eps ) Then
                            ! Call TOq(j,ind,t_j)
                            ! Call Tmq (i,ind,t_i)
                            ! result = result + (q**2)*alpha_mO(k)*t_i*t_j

                
                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(0,0,"0",3,Floor((j*1d0)/2d0),cp_j,r1)

                        Call T_lk(0,0,"0",1,Floor((i*1d0)/2d0),cp_i,r2)
                        
                        
                        result = result + (q**2)*alpha_mO(k)*r1*r2
                        end if 
                    end if 
                end do
            end do
    

    end if

            

    RETURN
END SUBROUTINE qq_mO



SUBROUTINE qq_QdQd(A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
    IMPLICIT NONE

    INTEGER ::  i,j,k
    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol

    Integer, INTENT(IN):: ind
    real*8 :: q,r1,r2
    real*8 , dimension(15):: alpha_QdQd
    real*8:: eps=EPSILON(result)

    Character :: cp_i,cp_j
    
    result = 0d0  
    
    if (ind==0) then
        q =  B_Multipoles(1)
        alpha_QdQd  = A_Pol(22:36) 
    
        


            do i = 1,5
                do j = 1,5

                    

                    if ( DABS(q)>eps  ) Then
                        Call GetIndex_QdQd(i,j,k)
                        if (  DABS(alpha_QdQd(k))>eps ) Then
                            ! Call TQdq(j,ind,t_j)
                            ! Call TQdq (i,ind,t_i)
                            ! result = result + (q**2)*alpha_QdQd(k)*t_i*t_j

                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(2,Floor((j*1d0)/2d0),cp_j,0,0,"0",r1)

                        Call T_lk(2,Floor((i*1d0)/2d0),cp_i,0,0,"0",r2)
                        
                        
                        result = result + (q**2)*alpha_QdQd(k)*r1*r2
                    end if
                    end if 
                end do
            end do

    else
        q            =  A_Multipoles(1)
        alpha_QdQd   =  B_Pol(22:36)

        


            do i = 1,5
                do j = 1,5

                    

                    if ( DABS(q)>eps  ) Then
                        Call GetIndex_QdQd(i,j,k)
                        if (  DABS(alpha_QdQd(k))>eps ) Then
                            ! Call TQdq(j,ind,t_j)
                            ! Call TQdq (i,ind,t_i)
                            ! result = result + (q**2)*alpha_QdQd(k)*t_i*t_j

                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(0,0,"0",2,Floor((j*1d0)/2d0),cp_j,r1)

                        Call T_lk(0,0,"0",2,Floor((i*1d0)/2d0),cp_i,r2)
                        
                        
                        result = result + (q**2)*alpha_QdQd(k)*r1*r2

                    end if 
                end if
                end do
            end do

    end if

            
    

    RETURN
END SUBROUTINE qq_QdQd