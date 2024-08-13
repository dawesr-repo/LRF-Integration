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
    
    do order=1,15 
        IF ( Coeff(ind)%I_Fit(order) > 0) THEN
        
            call induction_order(order,Coeff(ind)%A_Mult,Coeff(ind)%B_Pol,1,temp1) ! indB
            call induction_order(order,Coeff(ind)%B_Mult,Coeff(ind)%A_Pol,0,temp2) !indA
    
            IM = IM + temp1 + temp2
        END IF  
    end do
end SUBROUTINE Induction_Sph3

SUBROUTINE induction_order(order,mult_coeff,pol_kk_coeff,index,energy)
    use Geometry_Constant_v2, only: cal_coord_v2
    use FitConstants, only: C1,C2,C3
    IMPLICIT NONE
    real*8, INTENT(OUT)  :: energy
    integer , INTENT(IN) :: order,index
    real*8 :: R, res=0d0,temp=0d0
    integer :: l1,l2,i,j,lmin,lmax
    real*8 :: mult_coeff(225),pol_kk_coeff(6,7,195)

    R=cal_coord_v2(1)

    

    do l1=1,order-3
        do l2=1,order-3
            if (l1+l2+2 <= order) Then
             lmin = minval([l1,l2])
             lmax = maxval([l1,l2])
                do i=0,order-2-2
                    do j=0,order-2-2
                        if (i+j+l1+l2+2 == order)Then
                            call induction_ij_l1l2(i,j,l1,l2,mult_coeff, &
                                                    pol_kk_coeff(lmin,lmax,1:(2*l1+1)*(2*l2+1)), &
                                                    index,temp)
                            res  =  res + temp
                        end if
                    end do
                end do
            end if
        end do
    end do
  
    energy =  (-0.5*(C3*C1*(C2**order))*res)/(R**order)  

end SUBROUTINE induction_order

SUBROUTINE induction_ij_l1l2(i,j,l1,l2,mult,pol_arr,index,energy)
    use Geometry_Constant_v2, only: T_Tensor_v2
    IMPLICIT NONE
    integer, INTENT(IN) :: i,j,l1,l2,index
    real*8, INTENT(IN) :: mult(225),pol_arr(195)
    real*8, INTENT(OUT)  :: energy
    real*8 :: Qai,Qbj,comp_a_k1_k2,T_l1_i,T_l2_j,res
    integer :: ci,cj,k1,k2,cpn
    
    real*8 :: EPS = EPSILON(energy)
  
    res = 0d0

    do ci = 0,2*i
        Qai = mult(i**2 +1+ci)
        if (Dabs(Qai)>EPS) then
            do cj = 0,2*j
                Qbj = mult(j**2 +1+cj)
                if (Dabs(Qbj)>EPS) then

                    do k1 = 0,2*l1
                        do k2 = 0,2*l2

                            Call get_IND_cpn(l1,l2,k1,k2,cpn) 
                            comp_a_k1_k2 = pol_arr(cpn)

                            if (Dabs(comp_a_k1_k2)>EPS) then
                                ! index indicate if Im calculating pol over A
                                ! or pol over B
                                if (index==0)   then      
                                    T_l1_i = T_Tensor_v2(l1+1,k1+1,i+1,ci+1)
                                    T_l2_j = T_Tensor_v2(l2+1,k2+1,j+1,cj+1)
                                else
                                    T_l1_i = T_Tensor_v2(i+1,ci+1,l1+1,k1+1)
                                    T_l2_j = T_Tensor_v2(j+1,cj+1,l2+1,k2+1)
                                end if

                                res = res + Qai*Qbj*comp_a_k1_k2*T_l1_i*T_l2_j

                             
                            end if

                        end do
                    end do
                end if
    
            end do 
        end if
    end do 

    energy = res
    
end SUBROUTINE induction_ij_l1l2




SUBROUTINE get_IND_cpn(l1,l2,li,lj,cpn) 
    IMPLICIT NONE
    integer, INTENT(IN) :: l1,l2,li,lj
    integer, INTENT(OUT) :: cpn

    if (l1>l2) then
        cpn = lj*(2*l1+1) + li + 1
    else
        cpn = li*(2*l2+1) + lj + 1
    end if

end SUBROUTINE get_IND_cpn


























