SUBROUTINE Dispersion_Sph2( ind,order, Energy)

    use Geometry_Constant, only: cal_coord
    use FitConstants, only: C1,C2,C3

    IMPLICIT NONE
    real*8, INTENT(OUT)  :: Energy
    integer , INTENT(IN) :: ind,order

    integer ::i,j,k,l
    real*8:: R,D_ij,result

    R =cal_coord(1);


    result = 0d0

    do i=1,order-2-3
    do j=i,order-2-2-i
        do k=1,order-2-1-i-j
        do l=k,order-2-i-j-k
                if (i+j+k+l+2 == order)then
                   
                    Call Dispersion_Energy(ind,i,j,k,l ,D_ij) 
                    result  = result +D_ij
                end if
        end do    
        end do
    end do    
    end do
    

    Energy =   -1d0*(C3*C1*C2**order)*result/ R**order

    RETURN
END SUBROUTINE Dispersion_Sph2

SUBROUTINE Dispersion_Energy(ind,i,j,k,l,Energy)
    use Geometry_Constant
    use FitConstants
    use Search
    IMPLICIT NONE
    integer, INTENT(IN) ::i,j,k,l,ind
    real*8, INTENT(OUT)  :: Energy
    real*8, allocatable  :: Darr(:) 
    integer::n
    real*8::result
    real*8:: eps=EPSILON(result), r1,r2
    INTEGER::ai,aj,bk,bl, indK 
    Character :: cp_i,cp_j,cp_k,cp_l

    Call Coeff(ind)%Get_Disp_Comp(i,j,k,l,Darr,n)

        result = 0d0  

        do ai = 1,2*i+1
            do aj = 1,2*j+1      
                do bk = 1,2*k+1
                    do bl = 1,2*l+1


                        !Call GetIndex_mm_mm(ai,aj,bk,bl,indK)
                        Call Get_Disp_Index(i,j,k,l,ai,aj,bk,bl,indK)
                        

                        if (DABS(Darr(indK))>eps) Then
                        

                            call Get_Comp(ai,cp_i)
                            call Get_Comp(aj,cp_j)
                            call Get_Comp(bk,cp_k)
                            call Get_Comp(bl,cp_l)

                            Call T_lk(i,Floor((ai*1d0)/2d0),cp_i,k,Floor((bk*1d0)/2d0),cp_k,r1)

                            Call T_lk(j,Floor((aj*1d0)/2d0),cp_j,l,Floor((bl*1d0)/2d0),cp_l,r2)
                            
                            
                            result = result + Darr(indK)*r1*r2

                        end if 
                    end do
                end do
            end do
        end do




    Energy = result

end SUBROUTINE Dispersion_Energy


























