SUBROUTINE HyperPolarizability_Sph2( ind,order, Energy)

    use Geometry_Constant, only: cal_coord
    use FitConstants, only: C1,C2,C3

    IMPLICIT NONE
    real*8, INTENT(OUT)  :: Energy
    integer , INTENT(IN) :: ind,order

    integer ::i,j,k,n,m,p
    real*8:: R,result,term,fact
    real*8:: H0,H1

    R =cal_coord(1);


    result = 0d0

    do i=0,order-3-2
    do j=i,order-3-1-i
    do k=j,order-3-i-j
        do m=1,order-3-i-j-k
        do n=m,order-3-i-j-k-m
        do p=n,order-3-i-j-k-m-n
                if (i+j+k+m+n+p+3 == order)then
                   
                    
                    Call HyperPolarizability_Energy(ind,i,j,k,m,n,p,0,H0) 
                    Call HyperPolarizability_Energy(ind,i,j,k,m,n,p,1,H1) 
                  
                    term  = H0+H1
      
                    ! if (i.NE.j .and. k.NE.l)then
                      
                    !   Call Induction_Energy(ind,j,i,k,l, 0 ,I_ij_lk_0) 
                    !   Call Induction_Energy(ind,j,i,k,l, 1 ,I_ij_lk_1)
                    !   term  = term +I_ij_lk_0+I_ij_lk_1
                    ! end if

                    call HFactor(i,j,k,m,n,p,fact)

                    result  = result + fact*term

                end if
        end do    
        end do
        end do 
    end do    
    end do
    end do 
    

    Energy =   (1.0d0/6.0d0)*(C3*C1*C2**order)*result/ R**order

    RETURN
END SUBROUTINE HyperPolarizability_Sph2

SUBROUTINE HFactor(i,j,k,m,n,p,fact)
IMPLICIT NONE
    integer, INTENT(IN) ::i,j,k,m,n,p
    real*8 , INTENT(OUT):: fact
    real*8 :: fact1, fact2
    fact1 =3d0
    fact2 =3d0

    if (i==j .and. j==k)then
        fact1 =1d0
    elseif (i.NE.j .and. i.NE.k .and. j.NE.k)then  
        fact1 =6d0
    end if
    if (m==n .and. n==p)then
        fact2 =1d0
    elseif (m.NE.n .and. m.NE.p .and. n.NE.p)then  
        fact2 =6d0
    end if
    
    fact = fact1*fact2

end SUBROUTINE HFactor

SUBROUTINE HyperPolarizability_Energy(ind,i,j,k,m,n,p,dir,Energy)
    use Geometry_Constant
    use FitConstants
    use Search
    IMPLICIT NONE
    integer, INTENT(IN) ::i,j,k,m,n,p,ind,dir
    real*8, INTENT(OUT)  :: Energy
    
    real*8::result
    real*8:: eps=EPSILON(result), r1,r2,r3
    INTEGER::ai,aj,ak,bn,bm,bp, indK 
    Character :: cp_i,cp_j,cp_k,cp_m,cp_n,cp_p

    real*8 :: Mi(2*i+1),Mj(2*j+1),Mk(2*k+1)  
    real*8, allocatable  :: Hmnp(:)
    integer::nArr,count
    

    if (dir==1)then
        Call Coeff(ind)%Get_Mult_Comp("A",i,Mi)
        Call Coeff(ind)%Get_Mult_Comp("A",j,Mj)
        Call Coeff(ind)%Get_Mult_Comp("A",k,Mk)
        Call Coeff(ind)%Get_HPol_Comp("B",m,n,p,Hmnp,nArr)
    else
        Call Coeff(ind)%Get_Mult_Comp("B",i,Mi)
        Call Coeff(ind)%Get_Mult_Comp("B",j,Mj)
        Call Coeff(ind)%Get_Mult_Comp("B",k,Mk)
        Call Coeff(ind)%Get_HPol_Comp("A",m,n,p,Hmnp,nArr)
    end if 

    
    result = 0d0  
count=0
    do ai = 1,2*i+1
        do aj = 1,2*j+1
            do ak = 1,2*k+1       
                do bm = 1,2*m+1
                do bn = 1,2*n+1
                do bp = 1,2*p+1
                count = count +1
                  
                    if (DABS(Mi(ai))>eps .and. DABS(Mj(aj))>eps .and. DABS(Mk(ak))>eps) Then
                        Call Get_Hyper_Index(m,n,p,bm,bn,bp,indK)
                         
                        if (DABS(Hmnp(indK))>eps)then
                        
                            call Get_Comp(ai,cp_i)
                            call Get_Comp(aj,cp_j)
                            call Get_Comp(ak,cp_k)
                            call Get_Comp(bm,cp_m)
                            call Get_Comp(bn,cp_n)
                            call Get_Comp(bp,cp_p)

                            if (dir==1)then
                                Call T_lk(i,Floor((ai*1d0)/2d0),cp_i,m,Floor((bm*1d0)/2d0),cp_m,r1)
                                Call T_lk(j,Floor((aj*1d0)/2d0),cp_j,n,Floor((bn*1d0)/2d0),cp_n,r2)
                                Call T_lk(k,Floor((ak*1d0)/2d0),cp_k,p,Floor((bp*1d0)/2d0),cp_p,r3)
                            else 
                                Call T_lk(m,Floor((bm*1d0)/2d0),cp_m,i,Floor((ai*1d0)/2d0),cp_i,r1)
                                Call T_lk(n,Floor((bn*1d0)/2d0),cp_n,j,Floor((aj*1d0)/2d0),cp_j,r2)
                                Call T_lk(p,Floor((bp*1d0)/2d0),cp_p,k,Floor((ak*1d0)/2d0),cp_k,r3)
                            end if
                         
                             result = result + Mi(ai)*Mj(aj)*Mk(ak)*Hmnp(indK)*r1*r2*r3
                        end if 
                    end if 
                end do
                end do
                end do
            end do
        end do
    end do

    Energy = result

end SUBROUTINE HyperPolarizability_Energy


























