SUBROUTINE Induction_Sph2( ind,order, Energy)

    use Geometry_Constant, only: cal_coord
    use FitConstants, only: C1,C2,C3

    IMPLICIT NONE
    real*8, INTENT(OUT)  :: Energy
    integer , INTENT(IN) :: ind,order

    integer ::i,j,k,l
    real*8:: R,result,term,fact
    real*8:: I_ij_kl_0,I_ij_kl_1,I_ij_lk_0,I_ij_lk_1

    R =cal_coord(1);


    result = 0d0

    do i=0,order-2-2
    do j=i,order-2-1-i
        do k=1,order-2-i-j
        do l=k,order-2-i-j-k
                if (i+j+k+l+2 == order)then
                   
                    write(*,*) 'order : ',order, i,j,k,l
                    Call Induction_Energy(ind,i,j,k,l ,0,I_ij_kl_0) 
                    Call Induction_Energy(ind,i,j,k,l ,1,I_ij_kl_1) 
                    write(*,*) 'I_ij_kl_0 : ',I_ij_kl_0
                    write(*,*) 'I_ij_kl_1 : ',I_ij_kl_1

                    term  = I_ij_kl_0
                    term  = term + I_ij_kl_1
                    if (i.NE.j .and. k.NE.l)then

                      write(*,*) 'order : ',order, j,i,k,l 
                      
                      Call Induction_Energy(ind,j,i,k,l, 0 ,I_ij_lk_0) 
                      Call Induction_Energy(ind,j,i,k,l, 1 ,I_ij_lk_1)
                      
                      write(*,*) 'I_ij_lk_0 : ',I_ij_lk_0
                      write(*,*) 'I_ij_lk_1 : ',I_ij_lk_1
                      
                      term  = term +I_ij_lk_0
                      term  = term +I_ij_lk_1
                    end if

                    if (i==j .and. k==l)then
                        fact = 1.0d0
                    else 
                        fact = 2.0d0
                    end if

                    result  = result + fact*term

                end if
        end do    
        end do
    end do    
    end do
    

    Energy =   -0.5d0*(C3*C1*C2**order)*result/ R**order

    RETURN
END SUBROUTINE Induction_Sph2


SUBROUTINE Induction_Energy(ind,i,j,k,l,dir,Energy)
    use Geometry_Constant
    use FitConstants
    use Search
    IMPLICIT NONE
    integer, INTENT(IN) ::i,j,k,l,ind,dir
    real*8, INTENT(OUT)  :: Energy
    
    real*8::result
    real*8:: eps=EPSILON(result), r1,r2
    INTEGER::ai,aj,bk,bl, indK 
    Character :: cp_i,cp_j,cp_k,cp_l

    real*8 :: Mi(2*i+1),Mj(2*j+1) 
    real*8, allocatable  :: Ikl(:)
    integer::n

    if (dir==0)then
        Call Coeff(ind)%Get_Mult_Comp("A",i,Mi)
        Call Coeff(ind)%Get_Mult_Comp("A",j,Mj)
        Call Coeff(ind)%Get_Pol_Comp("B",k,l,Ikl,n)
    else
        Call Coeff(ind)%Get_Mult_Comp("B",i,Mi)
        Call Coeff(ind)%Get_Mult_Comp("B",j,Mj)
        Call Coeff(ind)%Get_Pol_Comp("A",k,l,Ikl,n)
    end if 


    
    result = 0d0  

    do ai = 1,2*i+1
        do aj = 1,2*j+1      
            do bk = 1,2*k+1
                do bl = 1,2*l+1
                    Call Get_Ind_Index(k,l,bk,bl,indK)            

                    if (DABS(Ikl(indK))>eps .and. Mi(ai)>eps .and. Mj(aj)>eps) Then
                    

                        call Get_Comp(ai,cp_i)
                        call Get_Comp(aj,cp_j)
                        call Get_Comp(bk,cp_k)
                        call Get_Comp(bl,cp_l)

                        Call T_lk(i,Floor((ai*1d0)/2d0),cp_i,k,Floor((bk*1d0)/2d0),cp_k,r1)

                        Call T_lk(j,Floor((aj*1d0)/2d0),cp_j,l,Floor((bl*1d0)/2d0),cp_l,r2)
                        
                        
                        result = result + Mi(ai)*Mj(aj)*Ikl(indK)*r1*r2

                    end if 
                end do
            end do
        end do
    end do




    Energy = result

end SUBROUTINE Induction_Energy


























