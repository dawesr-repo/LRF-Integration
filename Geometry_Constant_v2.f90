
    !********************************************************
module Geometry_Constant_v2  
   implicit none
   
    integer,parameter::Lev = 15
    public
    real*8 :: T_Tensor_v2(Lev,2*Lev+1,Lev,2*Lev+1)
 
    real*8 , dimension(3):: Ar_v2 
    real*8 , dimension(3):: Br_v2
    real*8 , dimension(9):: CC_v2
    real*8 , dimension(11):: cal_coord_v2
  contains

  subroutine init_Tensors_v2(maxlevel,coordenates)  
     implicit none 
     Integer, INTENT(IN) :: maxlevel
     real*8 ,dimension(6), INTENT(IN)  :: coordenates ! the angles are in degree 
     T_Tensor_v2  = 0d0


     call Generate_Coordenates_v2(coordenates)
     call calculate_tensor(maxlevel); 


  end subroutine init_Tensors_v2

  subroutine calculate_tensor(maxlevel) 
            ! % """Calculate all components of t_tensor up to maxlevel
            ! % 
            ! % Ar_v2gs
            ! %     maxlevel (int) will calculate all tensors within the constrain
            ! %                     la + lb < maxlevel
            ! % """

            integer ,INTENT(IN)::maxlevel
            integer :: order,la,ka_,lb,kb_
                

            do order=1,maxlevel
                do la=0, order-1
                    lb = order-la-1
                    do ka_ = 0,2*la
                        do kb_ = 0,2*lb
                            Call t_lk_iter(la, ka_, lb, kb_)
                        end do
                    end do
                end do
            end do
  
    
  end subroutine calculate_tensor

  subroutine t_lk_iter(la, ka_, lb, kb_)
            ! % """Based on the t-tensor recursive relationship with bottom-down
            ! %     approach to improve performance
            ! % Ar_v2gs
            ! % 
            ! %     la   (int)       Multipole order of Molecule A
            ! %     ka_  (int)       Component Order of Molecule A  0 < ka_ < 2*la
            ! %     lb   (int)       Multipole order of Molecule B
            ! %     kb_  (int)       Component Order of Molecule B  0 < kb_ < 2*lb
            ! % 
            ! % """

            Integer, INTENT(IN) :: la, ka_, lb, kb_
            real*8:: EPS=EPSILON(T_Tensor_v2(1,1,1,1))

            real*8:: res,comp_lk,comp_T,prod_comp,fact_prod,la_fact,l2_fact,l3_fact,l4_fact
            real*8:: lb_fact,la2_fact,fact_nn_1,fact_nn_2,cij,const,fact_nn,lb2_fact,m1,m2,m
            real*8:: r_comp 

            
            Integer:: ka1,rka1, kb1,rkb1,  rk1, rk_, rk_i, rk_j, i, j,n
            Character(len = 1), dimension(3):: coord = ["z", "x", "y"]! Cartesian Axis Labels
            Character(len = 1):: str_comp,rka2,rkb2,ka2, kb2,rk2
         
            call get_splitting_componet(ka_,ka2)
            call get_splitting_componet(kb_,kb2)
    
            ka1 =  floor((ka_ + 1.0d0)/2.0d0)
            kb1 =  floor((kb_ + 1.0d0)/2.0d0)

           
            if (la == 0 .and. lb == 0) then
                res = 1d0
            else
                ! reculrsive relation for lb = 0
                if (lb == 0)then
                    ! initializating component
                    comp_lk = 0d0
                    la_fact = (2d0*la-1d0)/(1d0*la)

                    ! loop though every coodinate axis
                    do i=1,3
                        ! new multipole components
                        call N_eta(coord(i), ka1, ka2,  rk1, rk2)
                        call Get_tensor_component(la, rk1, rk2,rk_)
                        call coeff_m(coord(i), ka1, ka2,m)
                        ! coeffcient NN of the recurence
                        call Factorial_nn(la-1, rk1, 0, 0,fact_nn)

                        if  (DABS(m)>EPS &
                            .and. la>=1             &
                            .and. fact_nn>EPS       &
                            .and. rk_ <= 2*(la-1) )then

                            comp_T =  T_Tensor_v2(la-1+1, rk_+1, 1, 1)
                            prod_comp =  Ar_v2(i)*comp_T
                            fact_prod = la_fact*m*fact_nn
                            comp_lk = comp_lk + fact_prod*prod_comp
                        end if
                    end do   

                    if (la >= 2 .and. ka_ <= 2*(la-2) .and. ka_ >=0) then
                        call Factorial_nn(la-2, ka1, 0, 0,fact_nn)
                        la2_fact = (la-1d0)/(1d0*la)
                        comp_lk = comp_lk -la2_fact* fact_nn*T_Tensor_v2(la-2+1, ka_+1, 1, 1)
                    end if

                    call Factorial_nn(la, ka1, lb, kb1,fact_nn)
                    
                    res = comp_lk/fact_nn 
                  
                ! reculrsive relation for la = 0
                elseif (la == 0) then
                    
                    ! initializating component
                    comp_lk = 0d0
                    lb_fact = (2d0*lb-1d0)/(1d0*lb)

                    ! loop though every coodinate axis
                    do i=1,3
                        ! new multipole components
                        call N_eta(coord(i), kb1, kb2,rk1, rk2)
                        call get_tensor_component(lb, rk1, rk2, rk_ )
                        call Coeff_m(coord(i), kb1, kb2,m)
                        call Factorial_nn(0, 0, lb-1, rk1,fact_nn)
                        

                        if (DABS(m) > EPS            &
                            .and. lb >= 1             &
                            .and. fact_nn > EPS       &
                            .and. rk_ <= 2*(lb-1) &
                            .and. rk_ >= 0) then

                            comp_lk = comp_lk+ lb_fact*m*fact_nn* Br_v2(i)* &
                                      T_Tensor_v2( 1, 1, lb-1+1 , rk_+1)
                        end if 
                    end do

                    if (lb >= 2 .and. kb_ <= 2*(lb-2) .and. kb_>=0)   then
                        Call Factorial_nn(0, 0, lb-2, kb1,fact_nn)
                        lb2_fact = (lb-1d0)/(1d0*lb)
                        comp_lk = comp_lk - lb2_fact*fact_nn* & 
                                  T_Tensor_v2( 1, 1, lb-2+1, kb_+1)
                    end if

                    Call Factorial_nn(la, ka1, lb, kb1,fact_nn)
                    res = comp_lk/fact_nn
                    
                !reculrsive relation for lb >0 .and. la>0
                else
                    ! initializating component
                    comp_lk = 0d0

                    if (ka_ <= 2*(la-2))Then
                        Call factorial_nn(la-2, ka1, lb, kb1,fact_nn_1)
                        comp_lk = comp_lk + fact_nn_1*T_Tensor_v2(la-2+1, ka_+1, lb+1, kb_+1)
                    end if 

                    if (kb_ <= 2*(lb-2)) then
                    
                        l2_fact = (2d0*la +lb-1d0)/(1d0*lb)
                        Call Factorial_nn(la, ka1, lb-2, kb1,fact_nn_2)
                        comp_lk = comp_lk-(l2_fact*fact_nn_2)*T_Tensor_v2(la+1, ka_+1, lb-2+1, kb_+1)
                    end if

                    do i=1,3
                        Call N_eta(coord(i), kb1, kb2,rk1, rk2)
                        Call Get_tensor_component(lb, rk1, rk2,rk_i)
                        Call Coeff_m(coord(i), kb1, kb2,m)
                        Call Factorial_nn(la, ka1,lb-1, rk1,fact_nn)
                        l3_fact = (2d0*(la + lb) -1d0)/(lb*1d0)
                        const = l3_fact*m*fact_nn

                        if (DABS(const) > EPS .and. rk_i <= 2*(lb-1)) then
                            comp_lk = comp_lk + const*Br_v2(i)*T_Tensor_v2(la+1, ka_+1, lb-1+1, rk_i+1)        
                        end if
                    end do

                    do i=1,3
                        do j=1,3
                            n = 3*(i-1) + j
                            l4_fact = (2d0*la-1d0)/(1d0*lb)

                            Call N_eta(coord(i), ka1, ka2,rka1, rka2)
                            Call N_eta(coord(j), kb1, kb2,rkb1, rkb2)
                            call Get_tensor_component(la, rka1, rka2, rk_i)
                            call Get_tensor_component(lb, rkb1, rkb2, rk_j)
                            call Coeff_m(coord(i), ka1, ka2, m1)
                            call Coeff_m(coord(j), kb1, kb2, m2)

                            Call Factorial_nn(la-1, rka1, lb-1, rkb1, fact_nn)
                            const = l4_fact*m1*m2*fact_nn
                            if  (DABS(const) > EPS &
                                .and. rk_i <= 2*(la-1) &
                                .and. rk_j <= 2*(lb-1)) then

                                comp_lk = comp_lk + const*CC_v2(n)*T_Tensor_v2(la-1+1, rk_i+1, lb-1+1, rk_j+1)
                            end if
                        enddo
                    enddo

                    call Factorial_nn(la, ka1, lb, kb1,fact_nn)
                    res = comp_lk/fact_nn
                end if
            end if 

          
            T_Tensor_v2( la+1, ka_+1, lb+1, kb_+1) = res
  end subroutine t_lk_iter

  SUBROUTINE Factorial(n,fact)
    IMPLICIT NONE

    integer, intent(in) :: n
    real*8, intent(inout) :: fact
    integer :: i


    fact = 1.0d0
    do i = 2, n
      fact = fact * (i*1d0)
    end do
    
    

  End   SUBROUTINE Factorial

  subroutine factorial_nn(la, ka1, lb, kb1,fn)
            ! % """_summary_
            ! % 
            ! % Ar_v2gs
            ! %     la (int) order of the Mutipole of molecule A
            ! %     ka (int) order of the component pf l-th Mutipole of molecule A
            ! %     lb (int) order of the Mutipole of molecule B
            ! %     kb (int) order of the component pf l-th Mutipole of molecule A
            ! % 
            ! % Returns
            ! %     float value of NN coefficient defined in the T-Tensor recursion
            ! % """

            Integer, INTENT(IN) :: la, ka1, lb, kb1
            Real*8, INTENT(OUT) :: fn
            Real*8 :: nf1, nf2, df1, df2

            if (la < 0 .or. lb < 0 .or. ka1 < 0 .or. kb1 < 0 .or. ka1 > la .or. kb1 > lb) then
                fn =  0d0
            else    
                call Factorial(la + ka1,nf1)
                call Factorial(lb + kb1,nf2)
                call Factorial(la-ka1,df1)
                call Factorial(lb-kb1,df2)

                fn = Dsqrt((nf1/df1)*(nf2/df2))
            end if
  end subroutine factorial_nn

  subroutine N_eta(mu, k1, k2,  ka1, ka2) 
            ! % """Auxiliar function to Calculate the recursive equation of T-Tensors
            ! % 
            ! % Ar_v2gs
            ! %     mu (str) Cartesian Axis "x","y" or "z"
            ! %     k1 (int) integer refering to the component order &&
            ! %     k2 (str) splitting "0","c","s" refer to the spherical components
            ! %               of for that given order Example k =[1,"c"]
            ! % Returns
            ! %     int integer refering to the component order &&
            ! %     str splitting "0","c","s" refer to the spherical components of
            ! %               for that given order Example k =[1,"c"]
            ! % 
            ! % """

            Character(len = 1), INTENT(IN) :: mu,k2
            INTEGER, INTENT(IN) :: k1

            Character(len = 1), INTENT(OUT) :: ka2
            INTEGER, INTENT(OUT) :: ka1
    
            if (mu == "x") then
                if (k1 <= 1) then
                    ka1 = 0
                    ka2 = "0"
                else
                    ka1= k1-1
                    ka2 = k2
                end if    
    
            elseif (mu == "y") then  
                if (k1 <= 1) then
                    ka1 = 0
                    ka2 = "0"
                else
                    if (k2 == "c") then
                        ka1 = k1-1
                        ka2 = "s"
                    else 
                        ka1 = k1-1
                        ka2  = "c"
                    endif 
                endif
            else
                ka1 = k1
                ka2  = k2
            endif
  end subroutine N_eta

  subroutine coeff_m(mu, k1, k2,m)
      ! % """Auxiliar function to Calculate the recursive equation of T-Tensors
      ! % 
      ! % Ar_v2gs
      ! %     mu (str) Cartesian Axis "x","y" or "z"
      ! %     k1 (int) integer refering to the component order
      ! %     k2 (str) splitting "0","c","s" refer to the spherical components
      ! %               of for that given order Example k =[1,"c"]
      ! % Return
      ! %     float    value of M coefficient defined in the T-Tensor recursion
      ! % """


      Character(len = 1), INTENT(IN) :: mu,k2
      Integer, INTENT(IN) :: k1
      real*8, Intent(OUT):: m


      m = 0d0


      if (mu == "x") then
          if (k1 == 1) then
              if (k2 == "c") then
                  m = Dsqrt(2.0d0)
              end if
          else
              m = 1d0*k1
          end if

      elseif (mu == "y") then
          if (k1 == 1) then
              if (k2 == "s") then
                  m = Dsqrt(2.0d0)
              end if
          else
              if (k2 == "s") then
                  m = 1d0*k1  
              else 
                  m = -1d0*k1
              end if
          end if
      else
          m = 1d0
      end if
  end subroutine coeff_m

  subroutine get_splitting_componet(i,str_comp)
            ! % """get the splitting of the Spherical Components given the tensor index
            ! % 
            ! % Ar_v2gs
            ! %     i (int) tensor index. i >= 0
            ! % 
            ! % Returns
            ! %     str the splitting of the Spherical Components
            ! % """
            Integer, INTENT(IN) :: i
            Character(len = 1), INTENT(OUT) :: str_comp

            if (i >= 0)then 
              if (i == 0) then
                  str_comp = "0"
              elseif (mod(i,2)== 1) then
                  str_comp = "c"
              else
                  str_comp = "s"
              end if  
            else
                str_comp = "-1"
            end if   
  end subroutine get_splitting_componet
  
  subroutine get_tensor_component(mult_ord, k1, k2,comp)
      ! % """
      ! % 
      ! % Ar_v2gs
      ! %     mult_ord (int) Multipole order
      ! %     k1 (int)       Component Order
      ! %     k2 (str)       Component splitting
      ! % 
      ! % Returns
      ! %     int tensor index starting by zero
      ! % """

      Integer, INTENT(IN) :: mult_ord, k1
      Character(len = 1), INTENT(IN) :: k2
      Integer, INTENT(OUT) :: comp

      if (k1 < 0 .or. mult_ord < 0 .or. k1 > mult_ord) then
          comp = -1
      elseif (k1 == 0) then
          if (k2 == "0") then
              comp = 0  
          else 
              comp = -1
          end if
      else
          if (k2 == "s") then
              comp = 2*k1
          elseif (k2 == "c") then
              comp = 2*k1-1
          else
              comp = 0
          end if
      end if
  end SUBROUTINE get_tensor_component
  
  SUBROUTINE Generate_Coordenates_v2(coordenates)

 
      real*8 ,dimension(6), INTENT(IN)  :: coordenates ! the angles are in degree
      real*8, parameter :: pii = DACOS(-1.d0)  


      real*8 :: cos_b1,cos_b2,cos_c1,cos_c2,sin_b1,sin_b2,sin_c1,sin_c2,cos_phi,sin_phi

          cal_coord_v2(1) = coordenates(1)
            

          cal_coord_v2(2)  = DCOS(coordenates(2)*pii/180d0)
          cal_coord_v2(3) = DSIN(coordenates(2)*pii/180d0)

          cal_coord_v2(4) = DCOS(coordenates(3)*pii/180d0)
          cal_coord_v2(5) = DSIN(coordenates(3)*pii/180d0)
      
          cal_coord_v2(6) = DCOS(coordenates(4)*pii/180d0)
          cal_coord_v2(7) = DSIN(coordenates(4)*pii/180d0)

          cal_coord_v2(8) = DCOS(coordenates(5)*pii/180d0)
          cal_coord_v2(9) = DSIN(coordenates(5)*pii/180d0)

          cal_coord_v2(10) = DCOS(coordenates(6)*pii/180d0)
          cal_coord_v2(11) = DSIN(coordenates(6)*pii/180d0)


          cos_b1 =    cal_coord_v2(2)
          sin_b1 =    cal_coord_v2(3)
          cos_b2 =    cal_coord_v2(4)
          sin_b2 =    cal_coord_v2(5)
          cos_phi =   cal_coord_v2(6)
          sin_phi =   cal_coord_v2(7)
          cos_c1 =    cal_coord_v2(8)
          sin_c1 =    cal_coord_v2(9)
          cos_c2 =    cal_coord_v2(10)
          sin_c2 =    cal_coord_v2(11)



          

          Ar_v2(1) = cos_b1           !Az
          Ar_v2(2) = sin_b1*sin_c1    !Ax
          Ar_v2(3) = cos_c1*sin_b1    !Ay
          
          
          Br_v2(1) = -cos_b2           !Bz
          Br_v2(2) = -sin_b2*sin_c2    !Bx
          Br_v2(3) = -cos_c2*sin_b2    !By
          
          CC_v2(1)= cos_b1*cos_b2 + cos_phi*sin_b1*sin_b2 !Czz
          CC_v2(2)= cos_c2*sin_phi*sin_b1 + (-cos_phi*cos_b2*sin_b1 + cos_b1*sin_b2)*sin_c2 !Czx
          CC_v2(3)= -cos_phi*cos_b2*cos_c2*sin_b1 + cos_b1*cos_c2*sin_b2 - sin_phi*sin_b1*sin_c2 !Czy


          CC_v2(4)= cos_b2*sin_b1*sin_c1 - sin_b2 *(cos_c1*sin_phi + cos_phi*cos_b1*sin_c1)  !Cxz
          CC_v2(5)= -cos_b1*cos_c2*sin_phi*sin_c1 + (cos_b2*cos_c1*sin_phi + sin_b1*sin_b2*sin_c1)*sin_c2 &
          + cos_phi *(cos_c1*cos_c2 + cos_b1*cos_b2*sin_c1*sin_c2) !Cxx
          CC_v2(6)= cos_c2*sin_b1*sin_b2*sin_c1 + cos_b2*cos_c2 *(cos_c1*sin_phi + cos_phi*cos_b1*sin_c1) &
          + (-cos_phi*cos_c1 + cos_b1*sin_phi*sin_c1)*sin_c2 !Cxy


          CC_v2(7)= cos_b2*cos_c1*sin_b1 + sin_b2 *(-cos_phi*cos_b1*cos_c1 + sin_phi*sin_c1) !Cyz
          CC_v2(8)= cos_c1*sin_b1*sin_b2*sin_c2 + cos_b1*cos_c1 *(-cos_c2*sin_phi + cos_phi*cos_b2*sin_c2) - &
          sin_c1 *(cos_phi*cos_c2 + cos_b2*sin_phi*sin_c2)   !Cyx
          CC_v2(9)= -cos_b2*cos_c2*sin_phi*sin_c1 + cos_c1 *(cos_c2*sin_b1*sin_b2 + cos_b1*sin_phi*sin_c2) &
          + cos_phi*(cos_b1*cos_b2*cos_c1*cos_c2 + sin_c1*sin_c2) !Cyy
        

  End   SUBROUTINE Generate_Coordenates_v2

end module Geometry_Constant_v2
  
