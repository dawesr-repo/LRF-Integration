
!********************************************************
SUBROUTINE Multipole_Sph2(ind,order, Energy)
    use Geometry_Constant, only: cal_coord
    use FitConstants, only: C1,C2,C3
    IMPLICIT NONE

    real*8, INTENT(OUT)  :: Energy     ! returned Energy
    integer , INTENT(IN) :: ind,order  !Coeff index and multipole order

    integer ::i
    real*8 :: R ,result,M_ij

    result =0d0

    R =cal_coord(1);
    
    do i=0,order-1

                    Call Multipole_Energy(ind,i,order-1-i,M_ij) 
                    result  = result +M_ij
   
    end do
    
    Energy =   (C3*C1*C2**order)*result/ R**order

    RETURN
END SUBROUTINE Multipole_Sph2
!********************************************************

SUBROUTINE Multipole_Energy(ind,i,j,result)
    use Geometry_Constant, only: Get_Comp,T_lk
    use FitConstants, only: Coeff
    IMPLICIT NONE
    integer, INTENT(IN) ::i,j,ind
    real*8, INTENT(OUT)  :: result
    real*8 :: MA(2*i+1),MB(2*j+1) 
    integer::ci,cj
    real*8:: eps=EPSILON(result),result_ci_cj
    Character :: cp1,cp2
   
    Call Coeff(ind)%Get_Mult_Comp("A",i,MA)

    Call Coeff(ind)%Get_Mult_Comp("B",j,MB)


        result = 0d0 

        do ci = 1,2*i+1
                do cj = 1,2*j+1

                    if (DABS(MA(ci)*MB(cj))>eps ) Then

                        call Get_Comp(ci,cp1)
                        call Get_Comp(cj,cp2)

                        Call T_lk(i,Floor((ci*1d0)/2d0),cp1,j,Floor((cj*1d0)/2d0),cp2,result_ci_cj)
                        result = result + MA(ci)*MB(cj)*result_ci_cj
                    end if 
                
                end do
        end do

end SUBROUTINE Multipole_Energy

