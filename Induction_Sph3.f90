! ****************************************************

SUBROUTINE Induction_Sph3(ind,IM)
!INDUCTION Summary of this function goes here
!   Detailed explanation goes here
    use FitConstants, only: Coeff
    IMPLICIT NONE
    integer , INTENT(IN) :: ind
    real*8, INTENT(OUT)  :: Im
    integer :: order
    real*8 :: temp1,temp2

    IM  = 0d0
    temp1 = 0d0
    temp2 = 0d0

    !print * , 'induction Pol',ind,Coeff(ind)%B_Pol(1,1,1:9)
    
    do order=1,15 
        IF ( Coeff(ind)%I_Fit(order) > 0) THEN
            call induction_order(order,ind,1,temp1) ! indB
            call induction_order(order,ind,0,temp2) ! indA
            IM = IM + temp1 + temp2
        END IF  
    end do
end SUBROUTINE Induction_Sph3

SUBROUTINE induction_order(order,ind,index,energy)
    use Geometry_Constant_v2, only: cal_coord_v2
    use FitConstants, only: C1,C2,C3
    IMPLICIT NONE
    real*8, INTENT(OUT)  :: energy
    integer , INTENT(IN) :: order,index,ind
    real*8 :: R,temp
    integer :: l1,l2,i,j,lmin,lmax
    real*8 :: res



    R=cal_coord_v2(1)

    res = 0d0
    temp= 0d0

    do l1=1,order-3
        do l2=1,order-3
            if (l1+l2+2 <= order) Then
             
                do i=0,order-2-l1-l2
                    do j=0,order-2-l1-l2
                        if (i+j+l1+l2+2 == order)Then
                            
                            call induction_ij_l1l2(i,j,l1,l2,ind, index,temp)
                            res  =  res + temp
                        end if
                    end do
                end do
            end if
        end do
    end do
  
    energy =  (-0.5d0*(C3*C1*(C2**order))*res)/(R**order)  

end SUBROUTINE induction_order

SUBROUTINE induction_ij_l1l2(i,j,l1,l2,ind,index,energy)
    use Geometry_Constant_v2, only: T_Tensor_v2
    use FitConstants, only: Coeff
    IMPLICIT NONE
    integer, INTENT(IN) :: i,j,l1,l2,index,ind
    real*8, INTENT(OUT)  :: energy
    real*8 :: Qai,Qbj,comp_a_k1_k2,T_l1_i,T_l2_j,res
    integer :: ci,cj,k1,k2,cpn,ni,nj,nl1,nl2,lmin,lmax
    real*8, allocatable :: Qa_cpn(:),Qb_cpn(:),pol_arr(:)
    real*8 :: EPS = EPSILON(energy)  !,Qa_cpn(2*i+1),Qb_cpn(2*j+1)
  
    res = 0d0
    ni = 2*i+1;
    nj = 2*j+1;
    nl1 = 2*l1+1;
    nl2 = 2*l2+1;
    lmin = minval([l1,l2])
    lmax = maxval([l1,l2])

    allocate(Qa_cpn(ni),Qb_cpn(nj),pol_arr(nl1*nl2))

    if (index==1) then
        Qa_cpn = Coeff(ind)%A_Mult(i**2 + 1:(i+1)**2)
        Qb_cpn = Coeff(ind)%A_Mult(j**2 + 1:(j+1)**2)
        pol_arr = Coeff(ind)%B_Pol(lmin,lmax,1:nl1*nl2)
        
    else
        Qa_cpn = Coeff(ind)%B_Mult(i**2 + 1:(i+1)**2)
        Qb_cpn = Coeff(ind)%B_Mult(j**2 + 1:(j+1)**2)
        pol_arr = Coeff(ind)%A_Pol(lmin,lmax,1:nl1*nl2)
    end if




    do ci = 1,ni
        Qai = Qa_cpn(ci);
        if (Dabs(Qai)>EPS) then
            do cj = 1,nj
                Qbj = Qb_cpn(cj);
                if (Dabs(Qbj)>EPS) then

                    do k1 = 1,nl1
                        do k2 = 1,nl2

                            Call get_IND_cpn(l1,l2,k1,k2,cpn) 
                            comp_a_k1_k2 = pol_arr(cpn)

                            if (Dabs(comp_a_k1_k2)>EPS) then
                                ! index indicate if Im calculating pol over A
                                ! or pol over B
                                if (index==0)   then      
                                    
                                   
                                    res = res + Qai*Qbj*comp_a_k1_k2*(T_Tensor_v2(l1+1,k1,i+1,ci)* &
                                                                      T_Tensor_v2(l2+1,k2,j+1,cj));
                                
                             else
                                    

                                    res = res + Qai*Qbj*comp_a_k1_k2*(T_Tensor_v2(i+1,ci,l1+1,k1)* &
                                                                      T_Tensor_v2(j+1,cj,l2+1,k2));
                                
                              end if

                             
                            end if

                        end do
                    end do
                end if
    
            end do 
        end if
    end do 

    energy = res

    deallocate(Qa_cpn,Qb_cpn,pol_arr)
    
end SUBROUTINE induction_ij_l1l2




SUBROUTINE get_IND_cpn(l1,l2,li,lj,cpn) 
    IMPLICIT NONE
    integer, INTENT(IN) :: l1,l2,li,lj
    integer, INTENT(OUT) :: cpn

    if (l1>l2) then
        cpn = (lj-1)*(2*l1+1) + li;
    else
        cpn = (li-1)*(2*l2+1) + lj;
    end if

end SUBROUTINE get_IND_cpn


























