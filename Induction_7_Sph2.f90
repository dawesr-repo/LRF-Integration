
    !********************************************************
SUBROUTINE Induction_7_Sph2(cal_coord,Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, Ind_7_Energy)
    IMPLICIT NONE
    
    !   NEED TO DECLARE ALL THE SUBROUTINE ARGUMENTS and
    !   ANY OTHER VARIABLES LOCAL TO THE SUBROUTINE

    real*8, INTENT(INOUT)  :: Ind_7_Energy
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol
    real*8 :: R 
    real*8 , dimension(11), INTENT(IN):: cal_coord
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C

    real*8 :: qO_mm_0,qO_mm_1,mQd_mm_0,mQd_mm_1,mm_mQd_0,mm_mQd_1,qQd_mQd_0,qQd_mQd_1
    real*8 :: qQd_Qdm_0,qQd_Qdm_1,qm_QdQd_0,qm_QdQd_1,qm_mO_0,qm_mO_1,qm_Om_0,qm_Om_1
    
    R =cal_coord(1)
    
    Call qO_mm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qO_mm_0) 
    Call qO_mm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qO_mm_1)

    Call mQd_mm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,mQd_mm_0) 
    Call mQd_mm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,mQd_mm_1)

    Call mm_mQd(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,mm_mQd_0) 
    Call mm_mQd(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,mm_mQd_1)

    Call qQd_mQd(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qQd_mQd_0) 
    Call qQd_mQd(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qQd_mQd_1)

    Call qQd_Qdm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qQd_Qdm_0) 
    Call qQd_Qdm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qQd_Qdm_1)

    Call qm_QdQd(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qm_QdQd_0) 
    Call qm_QdQd(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qm_QdQd_1)

    Call qm_mO(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qm_mO_0) 
    Call qm_mO(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qm_mO_1)

    Call qm_Om(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qm_Om_0) 
    Call qm_Om(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qm_Om_1)

    Ind_7_Energy =   (-0.5d0)*(  2.0d0*qO_mm_0 + 2.0d0*qO_mm_1  +2.0d0*mQd_mm_0 +2.0d0*mQd_mm_1 &
                                +2.0d0*mm_mQd_0+ 2.0d0*mm_mQd_1 +2.0d0*qQd_mQd_0+2.0d0*qQd_mQd_1 & 
                                +2.0d0*qQd_Qdm_0+2.0d0*qQd_Qdm_1+2.0d0*qm_QdQd_0+2.0d0*qm_QdQd_1 &
                                +2.0d0*qm_mO_0 + 2.0d0*qm_mO_1 +2.0d0*qm_Om_0 + 2.0d0*qm_Om_1&
                                )/ R**7

    

  
 
    
    
    RETURN
END SUBROUTINE Induction_7_Sph2
    !********************************************************



SUBROUTINE qO_mm(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
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
        real*8 , dimension(7):: O
        real*8 , dimension(6):: alpha
        real*8:: eps=EPSILON(result)

        Character :: cp_h,cp_i,cp_j
            
        result = 0d0  
        
        
        if (ind==0) then
            q =  B_Multipoles(1)
            O =  B_Multipoles(10:16)
            alpha  = A_Pol(1:6) 
          
            do h=1,7
                do i = 1,3
                    do j = 1,3

                        Call GetIndex_mm(i,j,k)

                        if (DABS(q)>eps .and.  DABS(alpha(k))>eps .and. DABS(O(h))>eps) Then


                            call Get_Comp(h,cp_h)
                            call Get_Comp(i,cp_i)
                            call Get_Comp(j,cp_j)

                            Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp_i,3,Floor((h*1d0)/2d0),cp_h,r1)
                            Call T_lk(Ar,Br,C,1,Floor((j*1d0)/2d0),cp_j,0,0,"0",r2)
                                
                               
                            result = result + (q*O(h))*alpha(k)*r1*r2


                        end if 
                    end do
                end do
            end do

        else
            q       =  A_Multipoles(1)
            O       =  A_Multipoles(10:16)
            alpha   =  B_Pol(1:6) 

         

          
            do h=1,7
                do i = 1,3
                    do j = 1,3

                        Call GetIndex_mm(i,j,k)

                        if (DABS(q)>eps .and.  DABS(alpha(k))>eps .and. DABS(O(h))>eps) Then
                   

                            call Get_Comp(h,cp_h)
                            call Get_Comp(i,cp_i)
                            call Get_Comp(j,cp_j)

                            Call T_lk(Ar,Br,C,3,Floor((h*1d0)/2d0),cp_h,1,Floor((i*1d0)/2d0),cp_i,r1)
                            Call T_lk(Ar,Br,C,0,0,"0",1,Floor((j*1d0)/2d0),cp_j,r2)
                                
                               
                            result = result + (q*O(h))*alpha(k)*r1*r2

                        end if 
                    end do
                end do
            end do

        end if

            
    
        RETURN
END SUBROUTINE qO_mm


SUBROUTINE mQd_mm(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
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
            real*8 , dimension(6):: alpha
            real*8:: eps=EPSILON(result)
    
            Character :: cp_g,cp_h,cp_i,cp_j
            
            result = 0d0
            
            
            
            if (ind==0) then
                m        =  B_Multipoles(2:4)
                Qd       =  B_Multipoles(5:9)
                alpha  = A_Pol(1:6) 
            
                do g=1,3
                    do h=1,5
                        do i = 1,3
                            do j = 1,3
        
                                Call GetIndex_mm(i,j,k)
        
                                if (DABS(m(g))>eps .and.  DABS(alpha(k))>eps .and. DABS(Qd(h))>eps) Then
                                        ! Call TQdm(Ar,Br,C ,h,j,mod(ind+1,2),t_jh)
                                        ! Call Tmm(Ar,Br,C ,i,g,ind,t_ig)
                                        ! result = result + (m(g)*Qd(h))*alpha(k)*t_ig*t_jh
    
                                        call Get_Comp(g,cp_g)
                                        call Get_Comp(h,cp_h)
                                        call Get_Comp(i,cp_i)
                                        call Get_Comp(j,cp_j)
    
                                        Call T_lk(Ar,Br,C,1,Floor((j*1d0)/2d0),cp_j,2,Floor((h*1d0)/2d0),cp_h,r1)
    
                                        Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp_i,1,Floor((g*1d0)/2d0),cp_g,r2)
                                        
                                       
                                        result = result + (m(g)*Qd(h))*alpha(k)*r1*r2    
    
                                end if 
                            end do
                        end do
                    end do
                end do
    
            else
                m       =  A_Multipoles(2:4)
                Qd       =  A_Multipoles(5:9)
                alpha   =  B_Pol(1:6) 

                do g=1,3
                    do h=1,5
                        do i = 1,3
                            do j = 1,3
        
                                Call GetIndex_mm(i,j,k)
        
                                if (DABS(m(g))>eps .and.  DABS(alpha(k))>eps .and. DABS(Qd(h))>eps) Then
                                        ! Call TQdm(Ar,Br,C ,h,j,mod(ind+1,2),t_jh)
                                        ! Call Tmm(Ar,Br,C ,i,g,ind,t_ig)
                                        ! result = result + (m(g)*Qd(h))*alpha(k)*t_ig*t_jh
    
                                        call Get_Comp(g,cp_g)
                                        call Get_Comp(h,cp_h)
                                        call Get_Comp(i,cp_i)
                                        call Get_Comp(j,cp_j)
    
                                        Call T_lk(Ar,Br,C,2,Floor((h*1d0)/2d0),cp_h,1,Floor((j*1d0)/2d0),cp_j,r1)
    
                                        Call T_lk(Ar,Br,C,1,Floor((g*1d0)/2d0),cp_g,1,Floor((i*1d0)/2d0),cp_i,r2)
                                        
                                       
                                        result = result + (m(g)*Qd(h))*alpha(k)*r1*r2    
    
                                end if 
                            end do
                        end do
                    end do
                end do

            end if
    


            RETURN
END SUBROUTINE mQd_mm


         
    
SUBROUTINE mm_mQd(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
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
                    real*8 , dimension(15):: alpha_mQ
                    real*8:: eps=EPSILON(result)

                    Character :: cp_h,cp_i,cp_j,cp_g
            
                    result = 0d0 
                    
                    
                    
                    
                    if (ind==0) then
                        m =  B_Multipoles(2:4)
                        alpha_mQ  = A_Pol(7:21) 
                    
                        do g = 1,3
                            do h = 1,3
                                do i = 1,3
                                    do j = 1,5
                
                                        Call GetIndex_mQd(i,j,k)
                
                                        if (DABS(m(h))>eps .and. DABS(m(g))>eps .and.  DABS(alpha_mQ(k))>eps ) Then
                                                ! Call TQdm(Ar,Br,C ,j,h,ind,t_jh)
                                                ! Call Tmm (Ar,Br,C ,i,g,ind,t_ig)
                                                ! result = result + (m(g)*m(h))*alpha_mQ(k)*t_ig*t_jh
                                                
                                                call Get_Comp(g,cp_g)
                                                call Get_Comp(h,cp_h)
                                                call Get_Comp(i,cp_i)
                                                call Get_Comp(j,cp_j)
                
                                                Call T_lk(Ar,Br,C,2,Floor((j*1d0)/2d0),cp_j,1,Floor((h*1d0)/2d0),cp_h,r1)
                
                                                Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp_i,1,Floor((g*1d0)/2d0),cp_g,r2)
                                                
                                            
                                                result = result + (m(g)*m(h))*alpha_mQ(k)*r1*r2
            
            
    
                                        end if 
                                    end do
                                end do
                            end do
                        end do
            
                    else
                        m           =  A_Multipoles(2:4)
                        alpha_mQ    =  B_Pol(7:21) 

                        do g = 1,3
                            do h = 1,3
                                do i = 1,3
                                    do j = 1,5
                
                                        Call GetIndex_mQd(i,j,k)
                
                                        if (DABS(m(h))>eps .and. DABS(m(g))>eps .and.  DABS(alpha_mQ(k))>eps ) Then
                                                ! Call TQdm(Ar,Br,C ,j,h,ind,t_jh)
                                                ! Call Tmm (Ar,Br,C ,i,g,ind,t_ig)
                                                ! result = result + (m(g)*m(h))*alpha_mQ(k)*t_ig*t_jh
                                                
                                                call Get_Comp(g,cp_g)
                                                call Get_Comp(h,cp_h)
                                                call Get_Comp(i,cp_i)
                                                call Get_Comp(j,cp_j)
                
                                                Call T_lk(Ar,Br,C,1,Floor((h*1d0)/2d0),cp_h,2,Floor((j*1d0)/2d0),cp_j,r1)
                
                                                Call T_lk(Ar,Br,C,1,Floor((g*1d0)/2d0),cp_g,1,Floor((i*1d0)/2d0),cp_i,r2)
                                                
                                            
                                                result = result + (m(g)*m(h))*alpha_mQ(k)*r1*r2
            
            
    
                                        end if 
                                    end do
                                end do
                            end do
                        end do


                    end if
            
                     
                    
                    RETURN
END SUBROUTINE mm_mQd
    
    
SUBROUTINE qQd_mQd(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
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
        real*8 , dimension(15):: alpha_mQ
        real*8:: eps=EPSILON(result)

        Character :: cp_h,cp_i,cp_j
            
        result = 0d0 
        
        
        
        if (ind==0) then
            q  =  B_Multipoles(1)
            Qd =  B_Multipoles(5:9)
            alpha_mQ  = A_Pol(7:21) 
        
            do h = 1,5
                do i = 1,3
                    do j = 1,5

                        Call GetIndex_mQd(i,j,k)

                        if (DABS(q)>eps .and. DABS(Qd(h))>eps .and.  DABS(alpha_mQ(k))>eps ) Then
                                ! Call TQdQd(Ar,Br,C ,j,h,ind,t_jh)
                                ! Call Tmq (Ar,Br,C ,i,ind,t_i)
                                ! result = result + (q*Qd(h))*alpha_mQ(k)*t_i*t_jh

                           
                            call Get_Comp(h,cp_h)
                            call Get_Comp(i,cp_i)
                            call Get_Comp(j,cp_j)

                            Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp_i,2,Floor((h*1d0)/2d0),cp_h,r1)

                            Call T_lk(Ar,Br,C,2,Floor((j*1d0)/2d0),cp_j,0,0,"0",r2)
                            
                           
                            result = result + (q*Qd(h))*alpha_mQ(k)*r1*r2


                        end if 
                    end do
                end do
            end do
            

        else
            q           =  A_Multipoles(1)
            Qd           =  A_Multipoles(5:9)
            alpha_mQ    =  B_Pol(7:21) 


            do h = 1,5
                do i = 1,3
                    do j = 1,5

                        Call GetIndex_mQd(i,j,k)

                        if (DABS(q)>eps .and. DABS(Qd(h))>eps .and.  DABS(alpha_mQ(k))>eps ) Then
                           
                            call Get_Comp(h,cp_h)
                            call Get_Comp(i,cp_i)
                            call Get_Comp(j,cp_j)

                            Call T_lk(Ar,Br,C,2,Floor((h*1d0)/2d0),cp_h,1,Floor((i*1d0)/2d0),cp_i,r1)

                            Call T_lk(Ar,Br,C,0,0,"0",2,Floor((j*1d0)/2d0),cp_j,r2)
                            
                           
                            result = result + (q*Qd(h))*alpha_mQ(k)*r1*r2


                        end if 
                    end do
                end do
            end do

        end if

                 
 
            
     
        RETURN
END SUBROUTINE qQd_mQd

SUBROUTINE qQd_Qdm(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
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
        real*8 , dimension(15):: alpha_mQ
        real*8:: eps=EPSILON(result)

        Character :: cp_h,cp_i,cp_j
            
        result = 0d0  


        
        
        if (ind==0) then
            q =  B_Multipoles(1)
            Qd =  B_Multipoles(5:9)
            alpha_mQ  = A_Pol(7:21) 
        
            do h = 1,5
                do i = 1,3
                    do j = 1,5

                        Call GetIndex_mQd(i,j,k)

                        if (DABS(q)>eps .and. DABS(Qd(h))>eps .and.  DABS(alpha_mQ(k))>eps ) Then
                                ! Call TQdq(Ar,Br,C ,j,ind,t_jh)
                                ! Call TQdm (Ar,Br,C ,h,i,mod(ind+1,2),t_i)
                                ! result = result + (q*Qd(h))*alpha_mQ(k)*t_i*t_jh

                                call Get_Comp(h,cp_h)
                                call Get_Comp(i,cp_i)
                                call Get_Comp(j,cp_j)

                                Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp_i,0,0,"0",r1)

                                Call T_lk(Ar,Br,C,2,Floor((j*1d0)/2d0),cp_j,2,Floor((h*1d0)/2d0),cp_h,r2)
                            
                                result = result + (q*Qd(h))*alpha_mQ(k)*r1*r2



                        end if 
                    end do
                end do
            end do

        else
            q           =  A_Multipoles(1)
            Qd          =  A_Multipoles(5:9)
            alpha_mQ    =  B_Pol(7:21) 

            do h = 1,5
                do i = 1,3
                    do j = 1,5

                        Call GetIndex_mQd(i,j,k)


                        if (DABS(q)>eps .and. DABS(Qd(h))>eps .and.  DABS(alpha_mQ(k))>eps ) Then
                                ! Call TQdq(Ar,Br,C ,j,ind,t_jh)
                                ! Call TQdm (Ar,Br,C ,h,i,mod(ind+1,2),t_i)
                                ! result = result + (q*Qd(h))*alpha_mQ(k)*t_i*t_jh

                                call Get_Comp(h,cp_h)
                                call Get_Comp(i,cp_i)
                                call Get_Comp(j,cp_j)

                                Call T_lk(Ar,Br,C,0,0,"0",1,Floor((i*1d0)/2d0),cp_i,r1)

                                Call T_lk(Ar,Br,C,2,Floor((h*1d0)/2d0),cp_h,2,Floor((j*1d0)/2d0),cp_j,r2)
                            
                                result = result + (q*Qd(h))*alpha_mQ(k)*r1*r2



                        end if 
                    end do
                end do
            end do


        end if

                 
 
            
     
        RETURN
END SUBROUTINE qQd_Qdm

SUBROUTINE qm_QdQd(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
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
        real*8 , dimension(15):: alpha_QdQd
        real*8:: eps=EPSILON(result)

        Character :: cp_h,cp_i,cp_j
            
        result = 0d0  
        
        
        
        if (ind==0) then
            q =  B_Multipoles(1)
            m =  B_Multipoles(2:4)
            alpha_QdQd  = A_Pol(22:36) 
            
            do h=1,3
                do i = 1,5
                    do j = 1,5

                        Call GetIndex_QdQd(i,j,k)

                        if ( DABS(q)>eps .and. DABS(m(h))>eps .and.  DABS(alpha_QdQd(k))>eps ) Then
                                ! Call TQdq(Ar,Br,C ,j,ind,t_j)
                                ! Call TQdm (Ar,Br,C ,i,h,ind,t_ih)
                                ! result = result + (q*m(h))*alpha_QdQd(k)*t_ih*t_j

                          
                                call Get_Comp(h,cp_h)
                                call Get_Comp(i,cp_i)
                                call Get_Comp(j,cp_j)

                                Call T_lk(Ar,Br,C,2,Floor((j*1d0)/2d0),cp_j,0,0,"0",r1)

                                Call T_lk(Ar,Br,C,2,Floor((i*1d0)/2d0),cp_i,1,Floor((h*1d0)/2d0),cp_h,r2)
                                
                            
                                result = result + (q*m(h))*alpha_QdQd(k)*r1*r2


                        end if 
                    end do
                end do
            end do

        else
            q            =  A_Multipoles(1)
            m            =  A_Multipoles(2:4)
            alpha_QdQd   =  B_Pol(22:36)

            do h=1,3
                do i = 1,5
                    do j = 1,5

                        Call GetIndex_QdQd(i,j,k)

                        if ( DABS(q)>eps .and. DABS(m(h))>eps .and.  DABS(alpha_QdQd(k))>eps ) Then
                                ! Call TQdq(Ar,Br,C ,j,ind,t_j)
                                ! Call TQdm (Ar,Br,C ,i,h,ind,t_ih)
                                ! result = result + (q*m(h))*alpha_QdQd(k)*t_ih*t_j

                          
                                call Get_Comp(h,cp_h)
                                call Get_Comp(i,cp_i)
                                call Get_Comp(j,cp_j)

                                Call T_lk(Ar,Br,C,0,0,"0",2,Floor((j*1d0)/2d0),cp_j,r1)

                                Call T_lk(Ar,Br,C,1,Floor((h*1d0)/2d0),cp_h,2,Floor((i*1d0)/2d0),cp_i,r2)
                                
                            
                                result = result + (q*m(h))*alpha_QdQd(k)*r1*r2


                        end if 
                    end do
                end do
            end do


        end if


         
           
    
        RETURN
END SUBROUTINE qm_QdQd

SUBROUTINE qm_mO(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
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
        real*8 , dimension(21):: alpha_mO
        real*8:: eps=EPSILON(result)

        Character :: cp_h,cp_i,cp_j
            
        result = 0d0 
        
        
        
        if (ind==0) then
            q =  B_Multipoles(1)
            m =  B_Multipoles(2:4)
            alpha_mO  = A_Pol(37:57) 
        
            do h = 1,3
                do i = 1,3
                    do j = 1,7

                        Call GetIndex_mO(i,j,k)

                        if ( DABS(m(h))>eps .and.DABS(q)>eps .and.  DABS(alpha_mO(k))>eps ) Then
                                ! Call TOm(Ar,Br,C ,j,h,ind,t_jh)
                                ! Call Tmq (Ar,Br,C ,i,ind,t_i)
                                ! result = result + (q*m(h))*alpha_mO(k)*t_i*t_jh

                           
                                call Get_Comp(h,cp_h)
                                call Get_Comp(i,cp_i)
                                call Get_Comp(j,cp_j)

                                Call T_lk(Ar,Br,C,3,Floor((j*1d0)/2d0),cp_j,1,Floor((h*1d0)/2d0),cp_h,r1)

                                Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp_i,0,0,"0",r2)
                                
                            
                                result = result + (q*m(h))*alpha_mO(k)*r1*r2


                        end if 
                    end do
                end do
            end do

        else
            q           =  A_Multipoles(1)
            m           =  A_Multipoles(2:4)
            alpha_mO    =  B_Pol(37:57) 


            do h = 1,3
                do i = 1,3
                    do j = 1,7

                        Call GetIndex_mO(i,j,k)

                        if ( DABS(m(h))>eps .and.DABS(q)>eps .and.  DABS(alpha_mO(k))>eps ) Then
                                ! Call TOm(Ar,Br,C ,j,h,ind,t_jh)
                                ! Call Tmq (Ar,Br,C ,i,ind,t_i)
                                ! result = result + (q*m(h))*alpha_mO(k)*t_i*t_jh

                           
                                call Get_Comp(h,cp_h)
                                call Get_Comp(i,cp_i)
                                call Get_Comp(j,cp_j)

                                Call T_lk(Ar,Br,C,1,Floor((h*1d0)/2d0),cp_h,3,Floor((j*1d0)/2d0),cp_j,r1)

                                Call T_lk(Ar,Br,C,0,0,"0",1,Floor((i*1d0)/2d0),cp_i,r2)
                                
                            
                                result = result + (q*m(h))*alpha_mO(k)*r1*r2


                        end if 
                    end do
                end do
            end do

        end if

              

            
    
        RETURN
END SUBROUTINE qm_mO


SUBROUTINE qm_Om(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
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
            real*8 , dimension(21):: alpha_mO
            real*8:: eps=EPSILON(result)
    
            Character :: cp_h,cp_i,cp_j
            
            result = 0d0 
            
            
            
            if (ind==0) then
                q =  B_Multipoles(1)
                m =  B_Multipoles(2:4)
                alpha_mO  = A_Pol(37:57) 
            
                do h = 1,3
                    do i = 1,3
                        do j = 1,7
    
                            Call GetIndex_mO(i,j,k)
    
                            if ( DABS(m(h))>eps .and.DABS(q)>eps .and.  DABS(alpha_mO(k))>eps ) Then
       
                              
                                call Get_Comp(h,cp_h)
                                call Get_Comp(i,cp_i)
                                call Get_Comp(j,cp_j)

                                Call T_lk(Ar,Br,C,3,Floor((j*1d0)/2d0),cp_j,0,0,"0",r1)

                                Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp_i,1,Floor((h*1d0)/2d0),cp_h,r2)
                                
                               
                                result = result + (q*m(h))*alpha_mO(k)*r1*r2


                            end if 
                        end do
                    end do
                end do
    
            else
                q           =  A_Multipoles(1)
                m           =  A_Multipoles(2:4)
                alpha_mO    =  B_Pol(37:57) 

                do h = 1,3
                    do i = 1,3
                        do j = 1,7
    
                            Call GetIndex_mO(i,j,k)
    
                            if ( DABS(m(h))>eps .and.DABS(q)>eps .and.  DABS(alpha_mO(k))>eps ) Then
                                    ! Call TOq(Ar,Br,C ,j,ind,t_j)
                                    ! Call Tmm (Ar,Br,C ,i,h,ind,t_ih)
                                    ! result = result + (q*m(h))*alpha_mO(k)*t_ih*t_j

                              
                                call Get_Comp(h,cp_h)
                                call Get_Comp(i,cp_i)
                                call Get_Comp(j,cp_j)

                                Call T_lk(Ar,Br,C,0,0,"0",3,Floor((j*1d0)/2d0),cp_j,r1)

                                Call T_lk(Ar,Br,C,1,Floor((h*1d0)/2d0),cp_h,1,Floor((i*1d0)/2d0),cp_i,r2)
                                
                               
                                result = result + (q*m(h))*alpha_mO(k)*r1*r2


                            end if 
                        end do
                    end do
                end do

            end if
    
               
    
                
        
            RETURN
END SUBROUTINE qm_Om


