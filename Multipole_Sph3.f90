
!********************************************************
SUBROUTINE Multipole_Sph3(ind, Energy)
    use Geometry_Constant_v2, only: cal_coord_v2
    use FitConstants, only: C1,C2,C3,Coeff
    IMPLICIT NONE

    real*8, INTENT(OUT)  :: Energy     ! returned Energy
    integer , INTENT(IN) :: ind        !Coeff index

    integer ::order
    real*8 :: R ,EM_ref
    real*8:: T1,T2

   print*,"ind",ind
   
    R =cal_coord_v2(1);
    Energy = 0d0


    do order = 1, 15
       
        IF ( Coeff(ind)%M_Fit(order) > 0) THEN
            call Multipole_Order(ind,order, EM_ref)
            Energy =  Energy + (C3*C1*(C2**order))*EM_ref/ R**order
            
        END IF 
    end do
    

    RETURN
END SUBROUTINE Multipole_Sph3
!********************************************************

SUBROUTINE Multipole_Order(ind,order,E_order)

    use Geometry_Constant_v2, only: T_Tensor_v2
    use FitConstants, only: Coeff
    integer, INTENT(IN) :: order,ind
    real*8, INTENT(OUT)  :: E_order
    Integer::i,j,ci,cj
    real*8:: eps=EPSILON(E_order)
    real*8:: Qai,Qbj

    E_order = 0d0
   

    do i=0,order-1
    
        j = order-1-i;

    
        do ci = 0,2*i
            Qai = Coeff(ind)%A_Mult(i**2+1+ci);
            if (DABS(Qai)>eps) then
                do cj = 0,2*j
                    Qbj = Coeff(ind)%B_Mult(j**2+1+cj);
                
                    if (DABS(Qbj)>eps) then
                        !print*,"Mult",order,i,j,ci,cj,T_Tensor_v2(i+1,ci+1,j+1,cj+1)
                        E_order = E_order + Qai*Qbj*T_Tensor_v2(i+1,ci+1,j+1,cj+1);
                    end if
        
                end do 
            end if
        end do
    end do
    
    
  
    
END SUBROUTINE Multipole_Order



! SUBROUTINE Multipole_Energy(ind,i,j,result)
!     use Geometry_Constant, only: Get_Comp,T_lk
!     use FitConstants, only: Coeff
!     IMPLICIT NONE
!     integer, INTENT(IN) ::i,j,ind
!     real*8, INTENT(OUT)  :: result
!     real*8 :: MA(2*i+1),MB(2*j+1) 
!     integer::ci,cj
!     real*8:: eps=EPSILON(result),result_ci_cj
!     Character :: cp1,cp2
   
!     Call Coeff(ind)%Get_Mult_Comp("A",i,MA)

!     Call Coeff(ind)%Get_Mult_Comp("B",j,MB)


!         result = 0d0 

!         do ci = 1,2*i+1
!                 do cj = 1,2*j+1

!                     if (DABS(MA(ci)*MB(cj))>eps ) Then

!                         call Get_Comp(ci,cp1)
!                         call Get_Comp(cj,cp2)

!                         Call T_lk(i,Floor((ci*1d0)/2d0),cp1,j,Floor((cj*1d0)/2d0),cp2,result_ci_cj)
!                         result = result + MA(ci)*MB(cj)*result_ci_cj
!                     end if 
                
!                 end do
!         end do

! end SUBROUTINE Multipole_Energy

