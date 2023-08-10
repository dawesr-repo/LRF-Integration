
    !********************************************************
           SUBROUTINE T_l0(Ar,Br,C,q,M2,IND,lk , result)
            IMPLICIT NONE
    
            INTEGER ::  i
            real*8, INTENT(INOUT) ::result 
            real*8 , dimension(3), INTENT(IN):: Ar 
            real*8 , dimension(3), INTENT(IN):: Br
            real*8 , dimension(9), INTENT(IN):: C
            Integer, INTENT(IN):: ind,lk
            real*8 :: q,t_k1,t_k2
            real*8 , dimension(2*lk+1):: M2
            real*8:: eps=EPSILON(result)
            
       
    
                result = 0d0 
    
                do i = 0,lk
                

                    if (ind==0) then

                        if (i==0) then
                
                            if (DABS(q)>eps .and. DABS(M2(1))>eps) Then
                                Call T_lk(Ar,Br,C,lk,0,"0",0,0,"0",t_k1)
                             
                                result = result + q*M2(1)*t_k1 
                            end if 
                        else

                            if (DABS(q)>eps .and. DABS(M2(2*i))>eps) Then
                                Call T_lk(Ar,Br,C,lk,i,"c",0,0,"0",t_k1)
                                result = result + q*M2(2*i)*t_k1 
                            end if 
                            if ( DABS(q)>eps .and. DABS(M2(2*i+1))>eps) Then
                                Call T_lk(Ar,Br,C,lk,i,"s",0,0,"0",t_k2)
                                result = result + q*M2(2*i+1)*t_k2
                            end if 
                            
                        end if
                ELSE
                    if (i==0) then

                        if (DABS(q)>eps .and. DABS(M2(1))>eps) Then
                            Call T_lk(Ar,Br,C,0,0,"0",lk,0,"0",t_k1)
                            result = result + q*M2(1)*t_k1 
                        end if 
                    else

                            if (DABS(q)>eps .and. DABS(M2(2*i))>eps) Then
                                Call T_lk(Ar,Br,C,0,0,"0",lk,i,"c",t_k1)
                                result = result + q*M2(2*i)*t_k1 
                            end if 
                            if ( DABS(q)>eps .and. DABS(M2(2*i+1))>eps) Then
                                Call T_lk(Ar,Br,C,0,0,"0",lk,i,"s",t_k2)
                                result = result + q*M2(2*i+1)*t_k2
                            end if 
                            
                        end if
                        
                    end if 
                end do
        
            RETURN
            END SUBROUTINE T_l0


