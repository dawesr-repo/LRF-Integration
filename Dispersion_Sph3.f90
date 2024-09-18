!***********************************************************************
SUBROUTINE Dispersion_Sph3(ind, DM)

    !   Dispersion Summary of this function goes here
    !   Detailed explanation goes here
    use FitConstants, only: Coeff
    IMPLICIT NONE
    integer , INTENT(IN) :: ind
    real*8, INTENT(OUT)  :: DM
    integer :: order
    real*8 :: temp
    
    DM  = 0d0


    do order=1,15 
        IF ( Coeff(ind)%D_Fit(order) > 0) THEN
            call dispersion_order(ind,order,temp) 
            DM = DM + temp
        End IF
    end do

end SUBROUTINE Dispersion_Sph3

SUBROUTINE dispersion_order(ind,order,energy)
    use Geometry_Constant_v2, only: cal_coord_v2
    use FitConstants, only: C1,C2,C3     
    IMPLICIT NONE
    integer , INTENT(IN) :: order,ind
    real*8, INTENT(OUT)  :: energy
    integer :: l1,l2,t1,t2
    real*8 :: res,temp,R
    res = 0d0

    R=cal_coord_v2(1)

    do l1=1,order-2
    do l2=1,order-2-l1
        do t1=1,order-2-l1-l2
        do t2=1,order-2-l1-l2-t1
            
            if (l1+l2+t1+t2+2 == order)Then 

                call dispersion_l1l2_t1t2(ind,l1,l2,t1,t2,temp) 
                res  = res + temp

            end if
        end do
        end do
    end do
    end do
  
    energy =  -((C3*C1*(C2**order))*res)/ (R**order)  

end SUBROUTINE dispersion_order

SUBROUTINE dispersion_l1l2_t1t2(ind,l1,l2,t1,t2,energy)
    use Geometry_Constant_v2, only: T_Tensor_v2
    use FitConstants, only: Coeff
    IMPLICIT NONE
    integer , INTENT(IN) :: l1,l2,t1,t2,ind
    real*8, INTENT(OUT)  :: energy
    integer :: li,lj,ti,tj,cpn
    real*8 :: T_li_ti,T_lj_tj,res,disp_coeff
    real*8 :: disp_arr((2*l1+1)*(2*l2+1)*(2*t1+1)*(2*t2+1))
    Integer::lmin,lmax,tmin,tmax   
    Real*8::eps=EPSILON(res)


    res = 0d0
    lmin = minval([l1,l2])
    lmax = MAXVAL([l1,l2])
    tmin = minval([t1,t2])
    tmax = MAXVAL([t1,t2])

    disp_arr = Coeff(ind)%Disp(lmin,lmax,tmin,tmax, 1:(2*l1+1)*(2*l2+1)*(2*t1+1)*(2*t2+1))

  
    do li = 0,2*l1
        do lj = 0,2*l2
            do ti = 0,2*t1
                do tj = 0,2*t2


                    call get_DISP_cpn(l1,l2,t1,t2,li,lj,ti,tj,cpn)  
                    
                    disp_coeff = disp_arr(cpn)


                    if (Dabs(disp_coeff)>EPS) THEN
                        T_li_ti = T_Tensor_v2(l1+1,li+1,t1+1,ti+1)
                        T_lj_tj = T_Tensor_v2(l2+1,lj+1,t2+1,tj+1)
                        
                        res = res + disp_coeff*T_li_ti*T_lj_tj
                        
                    end if 
                    
                end do
            end do
        end do
    end do

   

    energy = res
end SUBROUTINE dispersion_l1l2_t1t2

SUBROUTINE get_DISP_cpn(l1,l2,t1,t2,li,lj,ti,tj,cpn)

    IMPLICIT NONE
    integer, INTENT(IN) :: l1,l2,t1,t2,li,lj,ti,tj
    integer, INTENT(OUT) :: cpn 

    if (l1>l2 .and. t1>t2) then
        cpn = lj*(2*l1+1)*(2*t2+1)*(2*t1+1) + li*(2*t2+1)*(2*t1+1) +  tj*(2*t1+1)+ ti+1    
    elseif (l1>l2 .and. t1<=t2 ) then  
        cpn = lj*(2*l1+1)*(2*t1+1)*(2*t2+1) + li*(2*t1+1)*(2*t2+1) +  ti*(2*t2+1)+ tj+1
    elseif (l1<=l2 .and. t1>t2 ) then
        cpn = li*(2*l2+1)*(2*t2+1)*(2*t1+1) + lj*(2*t2+1)*(2*t1+1) +  tj*(2*t1+1)+ ti+1
    else
        cpn = li*(2*l2+1)*(2*t1+1)*(2*t2+1) + lj*(2*t1+1)*(2*t2+1) +  ti*(2*t2+1)+ tj+1
    end if 

end SUBROUTINE get_DISP_cpn





























