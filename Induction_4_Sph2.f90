
    !********************************************************
    SUBROUTINE Induction_4_Sph2(cal_coord,Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, Ind_4_Energy)
    IMPLICIT NONE
    
    !   NEED TO DECLARE ALL THE SUBROUTINE ARGUMENTS and
    !   ANY OTHER VARIABLES LOCAL TO THE SUBROUTINE

    real*8, INTENT(INOUT)  :: Ind_4_Energy
    real*8 , dimension(64) , INTENT(IN) :: A_Multipoles,B_Multipoles
    real*8 , dimension(57) , INTENT(IN) :: A_Pol,B_Pol
    real*8 :: R 
    real*8 , dimension(11), INTENT(IN):: cal_coord
    real*8 , dimension(3), INTENT(IN):: Ar 
    real*8 , dimension(3), INTENT(IN):: Br
    real*8 , dimension(9), INTENT(IN):: C
    real*8 :: qq_mm_0,qq_mm_1

    R =cal_coord(1)
    
    Call qq_mm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,0,qq_mm_0) 
    Call qq_mm(Ar,Br,C ,A_Multipoles,B_Multipoles,A_Pol,B_Pol,1,qq_mm_1)


    
    Ind_4_Energy =   (-0.5d0)*(qq_mm_0 + qq_mm_1)/ R**4

    RETURN
    END SUBROUTINE Induction_4_Sph2
    !********************************************************



    SUBROUTINE qq_mm(Ar,Br,C ,A_Multipoles,B_Multipoles ,A_Pol,B_Pol, ind , result)
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
        real*8 , dimension(6):: alpha
        real*8:: eps=EPSILON(result)
        Character :: cp1,cp2

        
        result = 0d0  
        
        !Induction_4_Sph2.o Induction_5_Sph2.o Induction_6_Sph2.o Induction_7_Sph2.o Induction_8_Sph2.o Dispersion_6_Sph2.o Dispersion_7_Sph2.o HyperPolarizability_6_Sph2.o HyperPolarizability_7_Sph2.o
        
        if (ind==0) then
            q =  B_Multipoles(1)
            alpha  = A_Pol(1:6) 
        
            do i = 1,3
                do j = 1,3

                    

                    if (DABS(q)>eps)then
                        
                        Call GetIndex_mm(i,j,k)
                        if(DABS(alpha(k))>eps) Then

 
                            call Get_Comp(i,cp1)
                            call Get_Comp(j,cp2)

                            Call T_lk(Ar,Br,C,1,Floor((i*1d0)/2d0),cp1,0,0,"0",r1)
                            Call T_lk(Ar,Br,C,1,Floor((j*1d0)/2d0),cp2,0,0,"0",r2)
                         
                            result = result + (q**2)*alpha(k)*r1*r2

                    end if 
                    end if 
                end do
            end do

        else
            q      =  A_Multipoles(1)
            alpha  =  B_Pol(1:6) 

            do i = 1,3
                do j = 1,3

                    Call GetIndex_mm(i,j,k)

                    if (DABS(q)>eps .and.  DABS(alpha(k))>eps) Then

                            ! Call Tmq(Ar,Br,C ,i,ind,t_i)
                            ! Call Tmq(Ar,Br,C ,j,ind,t_j)
                            call Get_Comp(i,cp1)
                            call Get_Comp(j,cp2)

                            Call T_lk(Ar,Br,C,0,0,"0",1,Floor((i*1d0)/2d0),cp1,r1)
                            Call T_lk(Ar,Br,C,0,0,"0",1,Floor((j*1d0)/2d0),cp2,r2)
                           
                            result = result + (q**2)*alpha(k)*r1*r2

                    end if 
                end do
            end do
        end if

         

          


    
        RETURN
        END SUBROUTINE qq_mm


   
