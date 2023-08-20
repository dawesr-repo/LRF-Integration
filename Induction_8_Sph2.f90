
    !********************************************************
SUBROUTINE Induction_8_Sph2(cal_coord,Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, Ind_8_Energy)
    IMPLICIT NONE
    
    !   NEED TO DECLARE ALL THE SUBROUTINE ARGUMENTS and
    !   ANY OTHER VARIABLES LOCAL TO THE SUBROUTINE

    real*8, INTENT(INOUT)  :: Ind_8_Energy
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol
    real*8 :: R 
    real*8 , dimension(11), INTENT(IN):: cal_coord
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C

    real*8 :: qPhi_mm_0,qPhi_mm_1,mO_mm_0,mO_mm_1,QdQd_mm_0,QdQd_mm_1,qQd_QdQd_0,qQd_QdQd_1
    real*8 :: mm_QdQd_0, mm_QdQd_1,mQd_mQd_0,mQd_mQd_1,mQd_Qdm_0,mQd_Qdm_1,qO_Qdm_0,qO_Qdm_1
    real*8 :: qO_mQd_0,qO_mQd_1,qQd_mO_0,qQd_mO_1,qQd_Om_0,qQd_Om_1,mm_Om_0,mm_Om_1


    R =cal_coord(1)
    
    Call qPhi_mm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qPhi_mm_0) 
    Call qPhi_mm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qPhi_mm_1)

    Call mO_mm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,mO_mm_0) 
    Call mO_mm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,mO_mm_1)

    Call QdQd_mm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,QdQd_mm_0) 
    Call QdQd_mm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,QdQd_mm_1)

    Call qQd_QdQd(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qQd_QdQd_0) 
    Call qQd_QdQd(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qQd_QdQd_1) 

    Call mm_QdQd(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,mm_QdQd_0) 
    Call mm_QdQd(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,mm_QdQd_1) 

    Call mQd_mQd(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,mQd_mQd_0) 
    Call mQd_mQd(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,mQd_mQd_1) 

    Call mQd_Qdm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,mQd_Qdm_0) 
    Call mQd_Qdm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,mQd_Qdm_1) 

    Call qO_Qdm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qO_Qdm_0) 
    Call qO_Qdm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qO_Qdm_1) 

    Call qO_mQd(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qO_mQd_0) 
    Call qO_mQd(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qO_mQd_1) 

    Call qQd_Om(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qQd_Om_0) 
    Call qQd_Om(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qQd_Om_1)

    Call qQd_mO(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qQd_mO_0) 
    Call qQd_mO(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qQd_mO_1) 
    
    Call mm_Om(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,mm_Om_0) 
    Call mm_Om(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,mm_Om_1) 



    Ind_8_Energy =   (-0.5d0)*( 2.0d0*qPhi_mm_0 + 2.0d0*qPhi_mm_1 +2.0d0*mO_mm_0 + 2.0d0*mO_mm_1  &
                                +QdQd_mm_0 + QdQd_mm_1 +2.0d0*qQd_QdQd_0 + 2.0d0*qQd_QdQd_1 &
                                +mm_QdQd_0 + mm_QdQd_1 +2.0d0*mQd_mQd_0 +2.0d0*mQd_mQd_1 &
                                +2.0d0*mQd_Qdm_0 +2.0d0*mQd_Qdm_1 +2.0d0*qO_Qdm_0 +2.0d0*qO_Qdm_1 &
                                +2.0d0*qO_mQd_0 +2.0d0*qO_mQd_1 +2.0d0*qQd_mO_0 +2.0d0*qQd_mO_1 &
                                +2.0d0*qQd_Om_0 +2.0d0*qQd_Om_1 +2.0d0*mm_Om_0 +2.0d0*mm_Om_1 &
                                )/ R**8
               
                                

    


    RETURN
END SUBROUTINE Induction_8_Sph2
    !********************************************************

SUBROUTINE qPhi_mm(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
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
    real*8 , dimension(9):: Phi
    real*8 , dimension(6):: alpha
    real*8:: eps=EPSILON(result)
    Character :: cp_h,cp_i,cp_j
        
    result = 0d0 
    
    
    
    
    if (ind==0) then
        q =  B_Multipoles(1)
        Phi =  B_Multipoles(17:25)
        alpha  = A_Pol(1:6) 

        do h=1,9
            do i = 1,3
                do j = 1,3

                    Call GetIndex_mm(i,j,k)

                    if (DABS(q)>eps .and.  DABS(alpha(k))>eps .and. DABS(Phi(h))>eps) Then
                            ! Call TPhim(Ar,Br,C ,h,j,mod(ind+1,2),t_jh)
                            ! Call Tmq(Ar,Br,C ,i,ind,t_i)
                            ! result = result + (q*Phi(h))*alpha(k)*t_i*t_jh


                            
                            call Get_Comp(h,cp_h)
                            call Get_Comp(i,cp_i)
                            call Get_Comp(j,cp_j)

                            Call T_lk(Ar,Br,C,1,Floor((j*1d0)/2d0),cp_j,4,Floor((h*1d0)/2d0),cp_h,r1)

                            Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp_i,0,0,"0",r2)
                            
                        
                            result = result + (q*Phi(h))*alpha(k)*r1*r2



                    end if 
                end do
            end do
        end do
    

    else
        q       =  A_Multipoles(1)
        Phi     =  A_Multipoles(17:25)
        alpha   =  B_Pol(1:6) 


        do h=1,9
            do i = 1,3
                do j = 1,3

                    Call GetIndex_mm(i,j,k)

                    if (DABS(q)>eps .and.  DABS(alpha(k))>eps .and. DABS(Phi(h))>eps) Then
                            ! Call TPhim(Ar,Br,C ,h,j,mod(ind+1,2),t_jh)
                            ! Call Tmq(Ar,Br,C ,i,ind,t_i)
                            ! result = result + (q*Phi(h))*alpha(k)*t_i*t_jh


                            
                            call Get_Comp(h,cp_h)
                            call Get_Comp(i,cp_i)
                            call Get_Comp(j,cp_j)

                            Call T_lk(Ar,Br,C,4,Floor((h*1d0)/2d0),cp_h,1,Floor((j*1d0)/2d0),cp_j,r1)

                            Call T_lk(Ar,Br,C,0,0,"0",1,Floor((i*1d0)/2d0),cp_i,r2)
                            
                        
                            result = result + (q*Phi(h))*alpha(k)*r1*r2



                    end if 
                end do
            end do
        end do



    end if

        

        
        

    RETURN
END SUBROUTINE qPhi_mm


SUBROUTINE mO_mm(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
    IMPLICIT NONE

    INTEGER ::  i,j,k,h,g
    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    Integer, INTENT(IN):: ind
    real*8 :: r1,r2
    real*8 , dimension(3):: m
    real*8 , dimension(7):: O
    real*8 , dimension(6):: alpha
    real*8:: eps=EPSILON(result)

    Character :: cp_h,cp_i,cp_j,cp_g
        
    result = 0d0  
    
    
    
    if (ind==0) then
        m =    B_Multipoles(2:4)
        O =    B_Multipoles(10:16)
        alpha  = A_Pol(1:6) 
    
        do g=1,3  
            do h=1,7
                do i = 1,3
                    do j = 1,3

                        Call GetIndex_mm(i,j,k)

                        if (DABS(m(g))>eps .and.  DABS(alpha(k))>eps .and. DABS(O(h))>eps) Then
                

                                call Get_Comp(g,cp_g)
                                call Get_Comp(h,cp_h)
                                call Get_Comp(i,cp_i)
                                call Get_Comp(j,cp_j)

                                Call T_lk(Ar,Br,C,1,Floor((j*1d0)/2d0),cp_j,3,Floor((h*1d0)/2d0),cp_h,r1)

                                Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp_i,1,Floor((g*1d0)/2d0),cp_g,r2)
                                
                            
                                result = result + (m(g)*O(h))*alpha(k)*r1*r2


                        end if 
                    end do
                end do
            end do
        end do

    else
        m =    A_Multipoles(2:4)
        O =    A_Multipoles(10:16)
        alpha   =  B_Pol(1:6) 

        do g=1,3  
            do h=1,7
                do i = 1,3
                    do j = 1,3

                        Call GetIndex_mm(i,j,k)

                        if (DABS(m(g))>eps .and.  DABS(alpha(k))>eps .and. DABS(O(h))>eps) Then

                                call Get_Comp(g,cp_g)
                                call Get_Comp(h,cp_h)
                                call Get_Comp(i,cp_i)
                                call Get_Comp(j,cp_j)

                                Call T_lk(Ar,Br,C,3,Floor((h*1d0)/2d0),cp_h,1,Floor((j*1d0)/2d0),cp_j,r1)

                                Call T_lk(Ar,Br,C,1,Floor((g*1d0)/2d0),cp_g,1,Floor((i*1d0)/2d0),cp_i,r2)
                                
                            
                                result = result + (m(g)*O(h))*alpha(k)*r1*r2


                        end if 
                    end do
                end do
            end do
        end do

    end if

        

    

    RETURN
END SUBROUTINE mO_mm

SUBROUTINE QdQd_mm(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
    IMPLICIT NONE

    INTEGER ::  i,j,k,h,g
    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    Integer, INTENT(IN):: ind
    real*8 :: r1,r2
    real*8 , dimension(5):: Qd
    real*8 , dimension(6):: alpha
    real*8:: eps=EPSILON(result)

    Character :: cp_h,cp_i,cp_j,cp_g
        
    result = 0d0  
    
    
    
    if (ind==0) then
        Qd =    B_Multipoles(5:9)
        alpha  = A_Pol(1:6) 
    
        do g=1,5 
            do h=1,5
                do i = 1,3
                    do j = 1,3

                        Call GetIndex_mm(i,j,k)

                        if (DABS(Qd(g))>eps .and.  DABS(alpha(k))>eps .and. DABS(Qd(h))>eps) Then
                                ! Call TQdm(Ar,Br,C ,h,j,mod(ind+1,2),t_jh)
                                ! Call TQdm(Ar,Br,C ,g,i,mod(ind+1,2),t_ig)
                                ! result = result + (Qd(g)*Qd(h))*alpha(k)*t_ig*t_jh
                                call Get_Comp(g,cp_g)
                                call Get_Comp(h,cp_h)
                                call Get_Comp(i,cp_i)
                                call Get_Comp(j,cp_j)

                                Call T_lk(Ar,Br,C,1,Floor((j*1d0)/2d0),cp_j,2,Floor((h*1d0)/2d0),cp_h,r1)

                                Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp_i,2,Floor((g*1d0)/2d0),cp_g,r2)
                                
                            
                                result = result + (Qd(g)*Qd(h))*alpha(k)*r1*r2


                        end if 
                    end do
                end do
            end do
        end do

    else
        Qd =    A_Multipoles(5:9)
        alpha   =  B_Pol(1:6) 

        do g=1,5 
            do h=1,5
                do i = 1,3
                    do j = 1,3

                        Call GetIndex_mm(i,j,k)

                        if (DABS(Qd(g))>eps .and.  DABS(alpha(k))>eps .and. DABS(Qd(h))>eps) Then
                                ! Call TQdm(Ar,Br,C ,h,j,mod(ind+1,2),t_jh)
                                ! Call TQdm(Ar,Br,C ,g,i,mod(ind+1,2),t_ig)
                                ! result = result + (Qd(g)*Qd(h))*alpha(k)*t_ig*t_jh
                                call Get_Comp(g,cp_g)
                                call Get_Comp(h,cp_h)
                                call Get_Comp(i,cp_i)
                                call Get_Comp(j,cp_j)

                                Call T_lk(Ar,Br,C,2,Floor((h*1d0)/2d0),cp_h,1,Floor((j*1d0)/2d0),cp_j,r1)

                                Call T_lk(Ar,Br,C,2,Floor((g*1d0)/2d0),cp_g,1,Floor((i*1d0)/2d0),cp_i,r2)
                                
                            
                                result = result + (Qd(g)*Qd(h))*alpha(k)*r1*r2


                        end if 
                    end do
                end do
            end do
        end do

    end if



    RETURN
END SUBROUTINE QdQd_mm


SUBROUTINE qQd_QdQd(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
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
    real*8 , dimension(5):: Qd
    real*8 , dimension(15):: alpha_QdQd
    real*8:: eps=EPSILON(result)

    Character :: cp_h,cp_i,cp_j,cp_g
        
    result = 0d0  
    
    
    
    if (ind==0) then
        q =  B_Multipoles(1)
        Qd =  B_Multipoles(5:9)
        alpha_QdQd  = A_Pol(22:36) 
    
        do h=1,5
            do i = 1,5
                do j = 1,5

                    Call GetIndex_QdQd(i,j,k)

                    if ( DABS(q)>eps .and. DABS(Qd(h))>eps .and.  DABS(alpha_QdQd(k))>eps ) Then
                            ! Call TQdq(Ar,Br,C ,j,ind,t_j)
                            ! Call TQdQd (Ar,Br,C ,i,h,ind,t_ih)
                            ! result = result + (q*Qd(h))*alpha_QdQd(k)*t_ih*t_j

                        call Get_Comp(h,cp_h)
                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(Ar,Br,C,2,Floor((j*1d0)/2d0),cp_j,0,0,"0",r1)

                        Call T_lk(Ar,Br,C,2,Floor((i*1d0)/2d0),cp_i,2,Floor((h*1d0)/2d0),cp_h,r2)
                        
                        
                        result = result + (q*Qd(h))*alpha_QdQd(k)*r1*r2


                    end if 
                end do
            end do
        end do

    else
        q            =  A_Multipoles(1)
        Qd           =  A_Multipoles(5:9)
        alpha_QdQd   =  B_Pol(22:36)

        do h=1,5
            do i = 1,5
                do j = 1,5

                    Call GetIndex_QdQd(i,j,k)

                    if ( DABS(q)>eps .and. DABS(Qd(h))>eps .and.  DABS(alpha_QdQd(k))>eps ) Then
                            ! Call TQdq(Ar,Br,C ,j,ind,t_j)
                            ! Call TQdQd (Ar,Br,C ,i,h,ind,t_ih)
                            ! result = result + (q*Qd(h))*alpha_QdQd(k)*t_ih*t_j

                        call Get_Comp(h,cp_h)
                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(Ar,Br,C,0,0,"0",2,Floor((j*1d0)/2d0),cp_j,r1)

                        Call T_lk(Ar,Br,C,2,Floor((h*1d0)/2d0),cp_h,2,Floor((i*1d0)/2d0),cp_i,r2)
                        
                        
                        result = result + (q*Qd(h))*alpha_QdQd(k)*r1*r2


                    end if 
                end do
            end do
        end do

    end if

            

        

    RETURN
END SUBROUTINE qQd_QdQd


SUBROUTINE mm_QdQd(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
    IMPLICIT NONE

    INTEGER ::  i,j,k,h,g
    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    Integer, INTENT(IN):: ind
    real*8 :: r1,r2
    real*8 , dimension(3):: m
    real*8 , dimension(15):: alpha_QdQd
    real*8:: eps=EPSILON(result)

    Character :: cp_h,cp_i,cp_j,cp_g
        
    result = 0d0  
    
    
    
    if (ind==0) then
        m=  B_Multipoles(2:4)
        alpha_QdQd  = A_Pol(22:36) 

        do g=1,3
            do h=1,3
                do i = 1,5
                    do j = 1,5

                        Call GetIndex_QdQd(i,j,k)

                        if ( DABS(m(g))>eps .and. DABS(m(h))>eps .and.  DABS(alpha_QdQd(k))>eps ) Then
                                ! Call TQdm(Ar,Br,C ,j,g,ind,t_jg)
                                ! Call TQdm (Ar,Br,C ,i,h,ind,t_ih)
                                ! result = result + (m(g)*m(h))*alpha_QdQd(k)*t_ih*t_jg

                            call Get_Comp(g,cp_g)
                            call Get_Comp(h,cp_h)
                            call Get_Comp(i,cp_i)
                            call Get_Comp(j,cp_j)

                            Call T_lk(Ar,Br,C,2,Floor((j*1d0)/2d0),cp_j,1,Floor((h*1d0)/2d0),cp_h,r1)

                            Call T_lk(Ar,Br,C,2,Floor((i*1d0)/2d0),cp_i,1,Floor((g*1d0)/2d0),cp_g,r2)
                            
                            
                            result = result + (m(g)*m(h))*alpha_QdQd(k)*r1*r2




                        end if 
                    end do
                end do
            end do
        end do


    else
        m           =  A_Multipoles(2:4)
        alpha_QdQd   =  B_Pol(22:36)

        do g=1,3
            do h=1,3
                do i = 1,5
                    do j = 1,5

                        Call GetIndex_QdQd(i,j,k)

                        if ( DABS(m(g))>eps .and. DABS(m(h))>eps .and.  DABS(alpha_QdQd(k))>eps ) Then
                                ! Call TQdm(Ar,Br,C ,j,g,ind,t_jg)
                                ! Call TQdm (Ar,Br,C ,i,h,ind,t_ih)
                                ! result = result + (m(g)*m(h))*alpha_QdQd(k)*t_ih*t_jg

                            call Get_Comp(g,cp_g)
                            call Get_Comp(h,cp_h)
                            call Get_Comp(i,cp_i)
                            call Get_Comp(j,cp_j)

                            Call T_lk(Ar,Br,C,1,Floor((h*1d0)/2d0),cp_h,2,Floor((j*1d0)/2d0),cp_j,r1)

                            Call T_lk(Ar,Br,C,1,Floor((g*1d0)/2d0),cp_g,2,Floor((i*1d0)/2d0),cp_i,r2)
                            
                            
                            result = result + (m(g)*m(h))*alpha_QdQd(k)*r1*r2




                        end if 
                    end do
                end do
            end do
        end do
    end if


    

    RETURN
END SUBROUTINE mm_QdQd

SUBROUTINE mQd_mQd(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
    IMPLICIT NONE

    INTEGER ::  i,j,k,h,g
    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    Integer, INTENT(IN):: ind
    real*8 :: t_ih,t_jg
    real*8 , dimension(3):: m
    real*8 , dimension(5):: Qd
    real*8 , dimension(15):: alpha_mQ
    real*8:: eps=EPSILON(result)
    real*8 :: r1,r2

    Character :: cp_h,cp_i,cp_j,cp_g
        
    result = 0d0  
    
    
    
    if (ind==0) then
        m =  B_Multipoles(2:4)
        Qd =  B_Multipoles(5:9)
        alpha_mQ  = A_Pol(7:21) 
    
        do g = 1,3
            do h = 1,5
                do i = 1,3
                    do j = 1,5

                        Call GetIndex_mQd(i,j,k)

                        if (DABS(Qd(h))>eps .and. DABS(m(g))>eps .and.  DABS(alpha_mQ(k))>eps ) Then
                                ! Call TQdm(Ar,Br,C ,j,g,ind,t_jg)
                                ! Call TQdm (Ar,Br,C ,i,h,mod(ind+1,2),t_ih)
                                ! result = result + (m(g)*Qd(h))*alpha_mQ(k)*t_ih*t_jg

                            call Get_Comp(g,cp_g)
                            call Get_Comp(h,cp_h)
                            call Get_Comp(i,cp_i)
                            call Get_Comp(j,cp_j)

                            Call T_lk(Ar,Br,C,2,Floor((j*1d0)/2d0),cp_j,1,Floor((g*1d0)/2d0),cp_g,r1) ! T1

                            Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp_i,2,Floor((h*1d0)/2d0),cp_h,r2) ! T2
                            
                            
                            result = result + (m(g)*Qd(h))*alpha_mQ(k)*r1*r2


                        end if 
                    end do
                end do
            end do
        end do

    else
        m           =  A_Multipoles(2:4)
        Qd          =  A_Multipoles(5:9)
        alpha_mQ    =  B_Pol(7:21) 

        do g = 1,3
            do h = 1,5
                do i = 1,3
                    do j = 1,5

                        Call GetIndex_mQd(i,j,k)

                        if (DABS(Qd(h))>eps .and. DABS(m(g))>eps .and.  DABS(alpha_mQ(k))>eps ) Then
                                ! Call TQdm(Ar,Br,C ,j,g,ind,t_jg)
                                ! Call TQdm (Ar,Br,C ,i,h,mod(ind+1,2),t_ih)
                                ! result = result + (m(g)*Qd(h))*alpha_mQ(k)*t_ih*t_jg

                            call Get_Comp(g,cp_g)
                            call Get_Comp(h,cp_h)
                            call Get_Comp(i,cp_i)
                            call Get_Comp(j,cp_j)

                            Call T_lk(Ar,Br,C,1,Floor((g*1d0)/2d0),cp_g,2,Floor((j*1d0)/2d0),cp_j,r1)

                            Call T_lk(Ar,Br,C,2,Floor((h*1d0)/2d0),cp_h,1,Floor((i*1d0)/2d0),cp_i,r2)
                            
            
                            
                            result = result + (m(g)*Qd(h))*alpha_mQ(k)*r1*r2


                        end if 
                    end do
                end do
            end do
        end do

    end if

            
    
    RETURN
END SUBROUTINE mQd_mQd

SUBROUTINE mQd_Qdm(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
    IMPLICIT NONE

    INTEGER ::  i,j,k,h,g
    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    Integer, INTENT(IN):: ind
    real*8 :: r1,r2
    real*8 , dimension(3):: m
    real*8 , dimension(5):: Qd
    real*8 , dimension(15):: alpha_mQ
    real*8:: eps=EPSILON(result)

    Character :: cp_h,cp_i,cp_j,cp_g
            
    result = 0d0  
    
    
    
    if (ind==0) then
        m =  B_Multipoles(2:4)
        Qd =  B_Multipoles(5:9)
        alpha_mQ  = A_Pol(7:21) 
    
        do g = 1,3
            do h = 1,5
                do i = 1,3
                    do j = 1,5
    
                        Call GetIndex_mQd(i,j,k)
    
                        if (DABS(Qd(h))>eps .and. DABS(m(g))>eps .and.  DABS(alpha_mQ(k))>eps ) Then
                                ! Call TQdQd(Ar,Br,C ,j,h,ind,t_jh)
                                ! Call Tmm (Ar,Br,C ,i,g,ind,t_ig)
                                ! result = result + (m(g)*Qd(h))*alpha_mQ(k)*t_ig*t_jh
    
                                call Get_Comp(g,cp_g)
                                call Get_Comp(h,cp_h)
                                call Get_Comp(i,cp_i)
                                call Get_Comp(j,cp_j)
    
                                Call T_lk(Ar,Br,C,2,Floor((j*1d0)/2d0),cp_j,2,Floor((h*1d0)/2d0),cp_h,r1)
    
                                Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp_i,1,Floor((g*1d0)/2d0),cp_g,r2)
                                
                            
                                result = result + (m(g)*Qd(h))*alpha_mQ(k)*r1*r2
    
    
                        end if 
                    end do
                end do
            end do
        end do

    else
        m           =  A_Multipoles(2:4)
        Qd          =  A_Multipoles(5:9)
        alpha_mQ    =  B_Pol(7:21) 

        do g = 1,3
            do h = 1,5
                do i = 1,3
                    do j = 1,5
    
                        Call GetIndex_mQd(i,j,k)
    
                        if (DABS(Qd(h))>eps .and. DABS(m(g))>eps .and.  DABS(alpha_mQ(k))>eps ) Then
     
                                call Get_Comp(g,cp_g)
                                call Get_Comp(h,cp_h)
                                call Get_Comp(i,cp_i)
                                call Get_Comp(j,cp_j)
    
                                Call T_lk(Ar,Br,C,2,Floor((h*1d0)/2d0),cp_h,2,Floor((j*1d0)/2d0),cp_j,r1)
    
                                Call T_lk(Ar,Br,C,1,Floor((g*1d0)/2d0),cp_g,1,Floor((i*1d0)/2d0),cp_i,r2)
                                
                            
                                result = result + (m(g)*Qd(h))*alpha_mQ(k)*r1*r2
    
    
                        end if 
                    end do
                end do
            end do
        end do

    end if

         

    RETURN
END SUBROUTINE mQd_Qdm




SUBROUTINE qO_Qdm(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
    IMPLICIT NONE

    INTEGER ::  i,j,k,h,g
    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    Integer, INTENT(IN):: ind
    real*8 ::q,r1,r2
    real*8 , dimension(7):: O
    real*8 , dimension(15):: alpha_mQ
    real*8:: eps=EPSILON(result)

    Character :: cp_h,cp_i,cp_j,cp_g
            
    result = 0d0  
    
    
    
    if (ind==0) then
        q =  B_Multipoles(1)
        O =  B_Multipoles(10:16)
        alpha_mQ  = A_Pol(7:21) 
    
        do h = 1,7
            do i = 1,3
                do j = 1,5

                    Call GetIndex_mQd(i,j,k)

                    if (DABS(O(h))>eps .and. DABS(q)>eps .and.  DABS(alpha_mQ(k))>eps ) Then
                            ! Call TQdq(Ar,Br,C ,j,ind,t_j)
                            ! Call TOm (Ar,Br,C ,i,h,mod(ind+1,2),t_ih)
                            ! result = result + (q*O(h))*alpha_mQ(k)*t_ih*t_j

                        call Get_Comp(g,cp_g)
                        call Get_Comp(h,cp_h)
                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(Ar,Br,C,2,Floor((j*1d0)/2d0),cp_j,0,0,"0",r1)

                        Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp_i,3,Floor((h*1d0)/2d0),cp_h,r2)
                        
                       
                        result = result + (q*O(h))*alpha_mQ(k)*r1*r2


                    end if 
                end do
            end do
        end do

    else
        q =  A_Multipoles(1)
        O =  A_Multipoles(10:16)
        alpha_mQ    =  B_Pol(7:21) 

        do h = 1,7
            do i = 1,3
                do j = 1,5

                    Call GetIndex_mQd(i,j,k)

                    if (DABS(O(h))>eps .and. DABS(q)>eps .and.  DABS(alpha_mQ(k))>eps ) Then
                            ! Call TQdq(Ar,Br,C ,j,ind,t_j)
                            ! Call TOm (Ar,Br,C ,i,h,mod(ind+1,2),t_ih)
                            ! result = result + (q*O(h))*alpha_mQ(k)*t_ih*t_j

                        call Get_Comp(g,cp_g)
                        call Get_Comp(h,cp_h)
                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(Ar,Br,C,0,0,"0",2,Floor((j*1d0)/2d0),cp_j,r1)

                        Call T_lk(Ar,Br,C,3,Floor((h*1d0)/2d0),cp_h,1,Floor((i*1d0)/2d0),cp_i,r2)
                        
                       
                        result = result + (q*O(h))*alpha_mQ(k)*r1*r2


                    end if 
                end do
            end do
        end do

    end if
 
 
        

    RETURN
END SUBROUTINE qO_Qdm

SUBROUTINE qO_mQd(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
    IMPLICIT NONE

    INTEGER ::  i,j,k,h,g
    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    Integer, INTENT(IN):: ind
    real*8 ::q,r1,r2
    real*8 , dimension(7):: O
    real*8 , dimension(15):: alpha_mQ
    real*8:: eps=EPSILON(result)

    Character :: cp_h,cp_i,cp_j,cp_g
            
    result = 0d0  
    
    
    
    if (ind==0) then
        q =  B_Multipoles(1)
        O =  B_Multipoles(10:16)
        alpha_mQ  = A_Pol(7:21) 
    

        do h = 1,7
            do i = 1,3
                do j = 1,5

                    Call GetIndex_mQd(i,j,k)

                    if (DABS(O(h))>eps .and. DABS(q)>eps .and.  DABS(alpha_mQ(k))>eps ) Then
                            ! Call Tmq(Ar,Br,C ,i,ind,t_i)
                            ! Call TOQd (Ar,Br,C ,j,h,mod(ind+1,2),t_jh)
                            ! result = result + (q*O(h))*alpha_mQ(k)*t_i*t_jh

                     
                        call Get_Comp(h,cp_h)
                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp_i,0,0,"0",r1)

                        Call T_lk(Ar,Br,C,2,Floor((j*1d0)/2d0),cp_j,3,Floor((h*1d0)/2d0),cp_h,r2)
                        
                       
                        result = result + (q*O(h))*alpha_mQ(k)*r1*r2

                    end if 
                end do
            end do
        end do


    else
        q =  A_Multipoles(1)
        O =  A_Multipoles(10:16)
        alpha_mQ    =  B_Pol(7:21) 

        
        do h = 1,7
            do i = 1,3
                do j = 1,5

                    Call GetIndex_mQd(i,j,k)

                    if (DABS(O(h))>eps .and. DABS(q)>eps .and.  DABS(alpha_mQ(k))>eps ) Then

                        call Get_Comp(h,cp_h)
                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(Ar,Br,C,0,0,"0",1,Floor((i*1d0)/2d0),cp_i,r1)

                        Call T_lk(Ar,Br,C,3,Floor((h*1d0)/2d0),cp_h,2,Floor((j*1d0)/2d0),cp_j,r2)
                        
                       
                        result = result + (q*O(h))*alpha_mQ(k)*r1*r2

                    end if 
                end do
            end do
        end do


    end if


    RETURN
END SUBROUTINE qO_mQd


SUBROUTINE qQd_Om(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
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
    real*8 , dimension(5):: Qd
    real*8 , dimension(21):: alpha_mO
    real*8:: eps=EPSILON(result)
    Character :: cp_h,cp_i,cp_j,cp_g
            
    result = 0d0  
    
    
    
    
    if (ind==0) then
        q  =  B_Multipoles(1)
        Qd =  B_Multipoles(5:9)
        alpha_mO  = A_Pol(37:57) 
    
        do h = 1,5
            do i = 1,3
                do j = 1,7

                    Call GetIndex_mO(i,j,k)

                    if ( DABS(Qd(h))>eps .and.DABS(q)>eps .and.  DABS(alpha_mO(k))>eps ) Then
                            ! Call TOq(Ar,Br,C ,j,ind,t_j)
                            ! Call TQdm (Ar,Br,C ,i,h,mod(ind+1,2),t_ih)
                            ! result = result + (q*Qd(h))*alpha_mO(k)*t_ih*t_j

                       
                        call Get_Comp(h,cp_h)
                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(Ar,Br,C,3,Floor((j*1d0)/2d0),cp_j,0,0,"0",r1)

                        Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp_i,2,Floor((h*1d0)/2d0),cp_h,r2)
                        
                       
                        result = result + (q*Qd(h))*alpha_mO(k)*r1*r2


                    end if 
                end do
            end do
        end do

    else
        q            =  A_Multipoles(1)
        Qd           =  A_Multipoles(5:9)
        alpha_mO    =  B_Pol(37:57) 


        do h = 1,5
            do i = 1,3
                do j = 1,7

                    Call GetIndex_mO(i,j,k)

                    if ( DABS(Qd(h))>eps .and.DABS(q)>eps .and.  DABS(alpha_mO(k))>eps ) Then
                            ! Call TOq(Ar,Br,C ,j,ind,t_j)
                            ! Call TQdm (Ar,Br,C ,i,h,mod(ind+1,2),t_ih)
                            ! result = result + (q*Qd(h))*alpha_mO(k)*t_ih*t_j

                       
                        call Get_Comp(h,cp_h)
                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(Ar,Br,C,0,0,"0",3,Floor((j*1d0)/2d0),cp_j,r1)

                        Call T_lk(Ar,Br,C,2,Floor((h*1d0)/2d0),cp_h,1,Floor((i*1d0)/2d0),cp_i,r2)
                        
                       
                        result = result + (q*Qd(h))*alpha_mO(k)*r1*r2


                    end if 
                end do
            end do
        end do

    end if

      



    RETURN
END SUBROUTINE qQd_Om




SUBROUTINE qQd_mO(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
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
    real*8 , dimension(5):: Qd
    real*8 , dimension(21):: alpha_mO
    real*8:: eps=EPSILON(result)
    Character :: cp_h,cp_i,cp_j
            
    result = 0d0  
    
    
    
    
    if (ind==0) then
        q  =  B_Multipoles(1)
        Qd =  B_Multipoles(5:9)
        alpha_mO  = A_Pol(37:57) 
    
        do h = 1,5
            do i = 1,3
                do j = 1,7

                    Call GetIndex_mO(i,j,k)

                    if ( DABS(Qd(h))>eps .and.DABS(q)>eps .and.  DABS(alpha_mO(k))>eps ) Then
                            ! Call TOQd(Ar,Br,C ,j,h,ind,t_jh)
                            ! Call Tmq (Ar,Br,C ,i,ind,t_i)
                            ! result = result + (q*Qd(h))*alpha_mO(k)*t_i*t_jh

                        
                        call Get_Comp(h,cp_h)
                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(Ar,Br,C,3,Floor((j*1d0)/2d0),cp_j,2,Floor((h*1d0)/2d0),cp_h,r1)

                        Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp_i,0,0,"0",r2)
                        
                       
                        result = result +  (q*Qd(h))*alpha_mO(k)*r1*r2


                    end if 
                end do
            end do
        end do

    else
        q            =  A_Multipoles(1)
        Qd           =  A_Multipoles(5:9)
        alpha_mO    =  B_Pol(37:57) 

        do h = 1,5
            do i = 1,3
                do j = 1,7

                    Call GetIndex_mO(i,j,k)

                    if ( DABS(Qd(h))>eps .and.DABS(q)>eps .and.  DABS(alpha_mO(k))>eps ) Then
                            ! Call TOQd(Ar,Br,C ,j,h,ind,t_jh)
                            ! Call Tmq (Ar,Br,C ,i,ind,t_i)
                            ! result = result + (q*Qd(h))*alpha_mO(k)*t_i*t_jh

                        
                        call Get_Comp(h,cp_h)
                        call Get_Comp(i,cp_i)
                        call Get_Comp(j,cp_j)

                        Call T_lk(Ar,Br,C,2,Floor((h*1d0)/2d0),cp_h,3,Floor((j*1d0)/2d0),cp_j,r1)

                        Call T_lk(Ar,Br,C,0,0,"0",1,Floor((i*1d0)/2d0),cp_i,r2)
                        
                       
                        result = result +  (q*Qd(h))*alpha_mO(k)*r1*r2


                    end if 
                end do
            end do
        end do

    end if




    RETURN
END SUBROUTINE qQd_mO


SUBROUTINE mm_Om(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
    IMPLICIT NONE

    INTEGER ::  i,j,k,h,g
    real*8, INTENT(INOUT) ::result 
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    Integer, INTENT(IN):: ind
    real*8 :: r1,r2
    real*8 , dimension(3):: m
    real*8 , dimension(21):: alpha_mO
    real*8:: eps=EPSILON(result)
    Character :: cp_h,cp_i,cp_j,cp_g
            
    result = 0d0 
    
    
    
    
    if (ind==0) then
        m =  B_Multipoles(2:4)
        alpha_mO  = A_Pol(37:57) 
    
        do g = 1,3
            do h = 1,3
                do i = 1,3
                    do j = 1,7
    
                        Call GetIndex_mO(i,j,k)
    
                        if ( DABS(m(g))>eps .and. DABS(m(h))>eps  .and.  DABS(alpha_mO(k))>eps ) Then
                                ! Call TOm(Ar,Br,C ,j,h,ind,t_jh)
                                ! Call Tmm (Ar,Br,C ,i,g,ind,t_ig)
                                ! result = result + (m(g)*m(h))*alpha_mO(k)*t_ig*t_jh
    
                                
                            call Get_Comp(g,cp_g)
                            call Get_Comp(h,cp_h)
                            call Get_Comp(i,cp_i)
                            call Get_Comp(j,cp_j)
    
                            Call T_lk(Ar,Br,C,3,Floor((j*1d0)/2d0),cp_j,1,Floor((h*1d0)/2d0),cp_h,r1)
    
                            Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp_i,1,Floor((g*1d0)/2d0),cp_g,r2)
                            
                           
                            result = result + (m(g)*m(h))*alpha_mO(k)*r1*r2
    
    
                        end if 
                    end do
                end do
            end do
        end do

    else
        m          =  A_Multipoles(2:4)
        alpha_mO    =  B_Pol(37:57) 

        do g = 1,3
            do h = 1,3
                do i = 1,3
                    do j = 1,7
    
                        Call GetIndex_mO(i,j,k)
    
                        if ( DABS(m(g))>eps .and. DABS(m(h))>eps  .and.  DABS(alpha_mO(k))>eps ) Then
                                ! Call TOm(Ar,Br,C ,j,h,ind,t_jh)
                                ! Call Tmm (Ar,Br,C ,i,g,ind,t_ig)
                                ! result = result + (m(g)*m(h))*alpha_mO(k)*t_ig*t_jh
    
                                
                            call Get_Comp(g,cp_g)
                            call Get_Comp(h,cp_h)
                            call Get_Comp(i,cp_i)
                            call Get_Comp(j,cp_j)
    
                            Call T_lk(Ar,Br,C,1,Floor((h*1d0)/2d0),cp_h,3,Floor((j*1d0)/2d0),cp_j,r1)
    
                            Call T_lk(Ar,Br,C,1,Floor((g*1d0)/2d0),cp_g,1,Floor((i*1d0)/2d0),cp_i,r2)
                            
                           
                            result = result + (m(g)*m(h))*alpha_mO(k)*r1*r2
    
    
                        end if 
                    end do
                end do
            end do
        end do

    end if

  
    

    RETURN
END SUBROUTINE mm_Om