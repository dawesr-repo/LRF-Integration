
!********************************************************
module Geometry_Constant_v2
    implicit none
    real (kind=8), parameter ::  C1=627.5095d0
    real (kind=8), parameter ::  C2=0.529177249d0
    real (kind=8), parameter ::  C3=349.757d0
    integer (kind=4),parameter::L = 15
    real (kind=8) , dimension(3):: ar_v2
    real (kind=8) , dimension(3):: br_v2
    real (kind=8) , dimension(9):: cc_v2
    real (kind=8) , dimension(11):: cal_coord_v2
    real (kind=8) :: t_tensor_v2(L,2*L+1,L,2*L+1)

    private :: L, ar_v2, br_v2, cc_v2, cal_coord_v2, t_tensor_v2 ! variables can not be called out of the module
    private ::  tensors_initialization_v2 ,&
                calculate_tensor,&
                t_lk_iter,&
                factorial,&
                factorial_nn,&
                n_eta,&
                coeff_m,&
                get_splitting_componet,&
                get_tensor_component,&
                generate_coordenates_v2,&
                multipole_order,&
                induction_order,&
                induction_ij_l1l2,&
                get_induction_cpn,&
                dispersion_order,&
                dispersion_l1l2_t1t2,&
                get_dispersion_cpn

    public :: get_total_interaction_energy,&
            ! Testing Only * to be used when debugging individual components
              multipole_sph_v3,&
              induction_sph_v3,&
              dispersion_sph_v3


contains

    subroutine tensors_initialization_v2(maxlevel,coordinates)
        implicit none
        integer (kind=4), intent(in) :: maxlevel
        real (kind=8) ,dimension(6), intent(in)  :: coordinates ! the angles are in degree

        t_tensor_v2  = 0d0
        
        call generate_coordenates_v2(coordinates)
        call calculate_tensor(maxlevel);

    end subroutine tensors_initialization_v2

    subroutine calculate_tensor(maxlevel)
        ! % """Calculate all components of t_tensor up to maxlevel
        ! %
        ! % ar_v2gs
        ! %     maxlevel (int) will calculate all tensors within the constrain
        ! %                     la + lb < maxlevel
        ! % """

        implicit none
        integer (kind=4) ,intent(in)::maxlevel
        integer (kind=4) :: order,la,ka_,lb,kb_


        do order=1,maxlevel
            do la=0, order-1
                lb = order-la-1
                do ka_ = 0,2*la
                    do kb_ = 0,2*lb
                        call t_lk_iter(la, ka_, lb, kb_)
                    end do
                end do
            end do
        end do


    end subroutine calculate_tensor

    subroutine t_lk_iter(la, ka_, lb, kb_)
        ! % """Based on the t-tensor recursive relationship with bottom-down
        ! %     approach to improve performance
        ! % ar_v2gs
        ! %
        ! %     la   (int)       Multipole order of Molecule A
        ! %     ka_  (int)       Component Order of Molecule A  0 < ka_ < 2*la
        ! %     lb   (int)       Multipole order of Molecule B
        ! %     kb_  (int)       Component Order of Molecule B  0 < kb_ < 2*lb
        ! %
        ! % """
        implicit none
        integer (kind=4), intent(in) :: la, ka_, lb, kb_
        real (kind=8),parameter:: EPS = epsilon(t_tensor_v2(1,1,1,1))

        real (kind=8):: res,comp_lk,comp_t,prod_comp,fact_prod,la_fact,l2_fact,l3_fact,l4_fact
        real (kind=8):: lb_fact,la2_fact,fact_nn_1,fact_nn_2,cij,const,fact_nn,lb2_fact,m1,m2,m
        real (kind=8):: r_comp


        integer (kind=4):: ka1,rka1, kb1,rkb1,  rk1, rk_, rk_i, rk_j, i, j,n
        character(len = 1), dimension(3):: coord = ["z", "x", "y"]! Cartesian Axis Labels
        character(len = 1):: str_comp,rka2,rkb2,ka2, kb2,rk2

        ka2 = get_splitting_componet(ka_)
        kb2 = get_splitting_componet(kb_)

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
                    call n_eta(coord(i), ka1, ka2,  rk1, rk2)
                    rk_ = get_tensor_component(la, rk1, rk2)
                    m = coeff_m(coord(i), ka1, ka2)
                    ! coeffcient NN of the recurence
                    fact_nn = factorial_nn(la-1, rk1, 0, 0)

                    if  (dabs(m)>EPS &
                            .and. la>=1             &
                            .and. fact_nn>EPS       &
                            .and. rk_ <= 2*(la-1) )then

                        comp_t =  t_tensor_v2(la-1+1, rk_+1, 1, 1)
                        prod_comp =  ar_v2(i)*comp_t
                        fact_prod = la_fact*m*fact_nn
                        comp_lk = comp_lk + fact_prod*prod_comp
                    end if
                end do

                if (la >= 2 .and. ka_ <= 2*(la-2) .and. ka_ >=0) then

                    la2_fact = (la-1d0)/(1d0*la)
                    comp_lk = comp_lk - la2_fact * factorial_nn(la-2, ka1, 0, 0 ) *&
                                        t_tensor_v2(la-2+1, ka_+1, 1, 1)
                end if

                res = comp_lk/factorial_nn(la, ka1, lb, kb1)

                ! reculrsive relation for la = 0
            elseif (la == 0) then

                ! initializating component
                comp_lk = 0d0
                lb_fact = (2d0*lb-1d0)/(1d0*lb)

                ! loop though every coodinate axis
                do i=1,3
                    ! new multipole components
                    call n_eta(coord(i), kb1, kb2,rk1, rk2)
                    rk_ = get_tensor_component(lb, rk1, rk2)
                    m = coeff_m(coord(i), kb1, kb2)

                    fact_nn = factorial_nn(0, 0, lb-1, rk1)


                    if (dabs(m) > EPS            &
                            .and. lb >= 1             &
                            .and. fact_nn > EPS       &
                            .and. rk_ <= 2*(lb-1) &
                            .and. rk_ >= 0) then

                        comp_lk = comp_lk+ lb_fact*m*fact_nn* br_v2(i)* &
                                t_tensor_v2( 1, 1, lb-1+1 , rk_+1)
                    end if
                end do

                if (lb >= 2 .and. kb_ <= 2*(lb-2) .and. kb_>=0)   then

                    lb2_fact = (lb-1d0)/(1d0*lb)
                    comp_lk = comp_lk - lb2_fact * factorial_nn( 0, 0, lb-2, kb1)* &
                            t_tensor_v2( 1, 1, lb-2+1, kb_+1)
                end if

                res = comp_lk/factorial_nn(la, ka1, lb, kb1)

                !reculrsive relation for lb >0 .and. la>0
            else
                ! initializating component
                comp_lk = 0d0

                if (ka_ <= 2*(la-2))then

                    comp_lk = comp_lk + factorial_nn( la-2, ka1, lb, kb1) *&
                                        t_tensor_v2( la-2+1, ka_+1, lb+1, kb_+1 )
                end if

                if (kb_ <= 2*(lb-2)) then

                    l2_fact = (2d0*la +lb-1d0)/(1d0*lb)

                    comp_lk = comp_lk - (l2_fact * factorial_nn(la, ka1, lb-2, kb1)) *&
                                         t_tensor_v2( la+1, ka_+1, lb-2+1, kb_+1 )
                end if

                do i=1,3
                    call n_eta(coord(i), kb1, kb2,rk1, rk2)
                    rk_i = get_tensor_component(lb, rk1, rk2)
                    m = coeff_m(coord(i), kb1, kb2)
                    l3_fact = (2d0*(la + lb) -1d0)/(lb*1d0)
                    const = l3_fact*m*factorial_nn(la, ka1,lb-1, rk1)

                    if (dabs(const) > EPS .and. rk_i <= 2*(lb-1)) then
                        comp_lk = comp_lk + const*br_v2(i)*t_tensor_v2(la+1, ka_+1, lb-1+1, rk_i+1)
                    end if
                end do

                do i=1,3
                    do j=1,3
                        n = 3*(i-1) + j
                        l4_fact = (2d0*la-1d0)/(1d0*lb)

                        call n_eta(coord(i), ka1, ka2,rka1, rka2)
                        call n_eta(coord(j), kb1, kb2,rkb1, rkb2)
                        rk_i = get_tensor_component(la, rka1, rka2)
                        rk_j = get_tensor_component(lb, rkb1, rkb2)
                        m1 = coeff_m(coord(i), ka1, ka2)
                        m2 = coeff_m(coord(j), kb1, kb2)


                        const = l4_fact*m1*m2*factorial_nn(la-1, rka1, lb-1, rkb1)
                        if  (dabs(const) > EPS &
                                .and. rk_i <= 2*(la-1) &
                                .and. rk_j <= 2*(lb-1)) then

                            comp_lk = comp_lk + const*cc_v2(n)*t_tensor_v2(la-1+1, rk_i+1, lb-1+1, rk_j+1)
                        end if
                    enddo
                enddo

                res = comp_lk/factorial_nn(la, ka1, lb, kb1)
            end if
        end if


        t_tensor_v2( la+1, ka_+1, lb+1, kb_+1) = res
    end subroutine t_lk_iter

    function factorial(n)
        implicit none

        integer (kind=4), intent(in) :: n
        real (kind=8):: factorial
        integer (kind=4) :: i

        if (n>=0)then
            factorial = 1.0d0
            do i = 2, n
                factorial = factorial * (i*1d0)
            end do
        else
            factorial = 0d0
        end if
        
    end  function factorial

    function factorial_nn(la, ka1, lb, kb1)
        ! % """_summary_
        ! %
        ! % ar_v2gs
        ! %     la (int) order of the Mutipole of molecule A
        ! %     ka (int) order of the component pf l-th Mutipole of molecule A
        ! %     lb (int) order of the Mutipole of molecule B
        ! %     kb (int) order of the component pf l-th Mutipole of molecule A
        ! %
        ! % returns
        ! %     float value of NN coefficient defined in the T-Tensor recursion
        ! % """
        implicit none
        integer (kind=4), intent(in) :: la, ka1, lb, kb1
        real (kind=8) :: factorial_nn

        if (la < 0 .or. lb < 0 .or. ka1 < 0 .or. kb1 < 0 .or. ka1 > la .or. kb1 > lb) then
            factorial_nn =  0d0
        else
            factorial_nn = dsqrt((factorial(la + ka1)/factorial(la-ka1))*&
                                 (factorial(lb + kb1)/factorial(lb-kb1)))
        end if
    end function factorial_nn

    subroutine n_eta(mu, k1, k2,  ka1, ka2)
        ! % """Auxiliar function to Calculate the recursive equation of T-Tensors
        ! %
        ! % ar_v2gs
        ! %     mu (str) Cartesian Axis "x","y" or "z"
        ! %     k1 (int) integer (kind=4) refering to the component order &&
        ! %     k2 (str) splitting "0","c","s" refer to the spherical components
        ! %               of for that given order Example k =[1,"c"]
        ! % returns
        ! %     int integer (kind=4) refering to the component order &&
        ! %     str splitting "0","c","s" refer to the spherical components of
        ! %               for that given order Example k =[1,"c"]
        ! %
        ! % """

        implicit none
        character(len = 1), intent(in) :: mu,k2
        integer (kind=4), intent(in) :: k1

        character(len = 1), intent (out) :: ka2
        integer (kind=4), intent (out) :: ka1

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
    end subroutine n_eta

    function coeff_m(mu, k1, k2)
        ! % """Auxiliar function to Calculate the recursive equation of T-Tensors
        ! %
        ! % ar_v2gs
        ! %     mu (str) Cartesian Axis "x","y" or "z"
        ! %     k1 (int) integer (kind=4) refering to the component order
        ! %     k2 (str) splitting "0","c","s" refer to the spherical components
        ! %               of for that given order Example k =[1,"c"]
        ! % return
        ! %     float    value of M coefficient defined in the T-Tensor recursion
        ! % """

        implicit none
        character(len = 1), intent(in) :: mu,k2
        integer (kind=4), intent(in) :: k1
        real (kind=8):: coeff_m

        coeff_m = 0d0

        if (mu == "x") then
            if (k1 == 1) then
                if (k2 == "c") then
                    coeff_m = dsqrt(2.0d0)
                end if
            else
                coeff_m = 1d0*k1
            end if

        elseif (mu == "y") then
            if (k1 == 1) then
                if (k2 == "s") then
                    coeff_m = dsqrt(2.0d0)
                end if
            else
                if (k2 == "s") then
                    coeff_m = 1d0*k1
                else
                    coeff_m = -1d0*k1
                end if
            end if
        else
            coeff_m = 1d0
        end if
    end function coeff_m

    function get_splitting_componet(i)
        ! % """get the splitting of the Spherical Components given the tensor index
        ! %
        ! % ar_v2gs
        ! %     i (int) tensor index. i >= 0
        ! %
        ! % returns
        ! %     str the splitting of the Spherical Components
        ! % """
        implicit none
        integer (kind=4), intent(in) :: i
        character(len = 1):: get_splitting_componet

        if (i >= 0)then
            if (i == 0) then
                get_splitting_componet = "0"
            elseif (mod(i,2)== 1) then
                get_splitting_componet = "c"
            else
                get_splitting_componet = "s"
            end if
        else
            get_splitting_componet = "-1"
        end if
    end function get_splitting_componet

    function get_tensor_component(mult_ord, k1, k2)
        ! % """
        ! %
        ! % ar_v2gs
        ! %     mult_ord (int) Multipole order
        ! %     k1 (int)       Component Order
        ! %     k2 (str)       Component splitting
        ! %
        ! % returns
        ! %     int tensor index starting by zero
        ! % """
        implicit none
        integer (kind=4), intent(in) :: mult_ord, k1
        character(len = 1), intent(in) :: k2
        integer (kind=4):: get_tensor_component

        if (k1 < 0 .or. mult_ord < 0 .or. k1 > mult_ord) then
            get_tensor_component = -1
        elseif (k1 == 0) then
            if (k2 == "0") then
                get_tensor_component = 0
            else
                get_tensor_component = -1
            end if
        else
            if (k2 == "s") then
                get_tensor_component = 2*k1
            elseif (k2 == "c") then
                get_tensor_component = 2*k1-1
            else
                get_tensor_component = 0
            end if
        end if
    end function get_tensor_component

    subroutine generate_coordenates_v2(coordinates)

        implicit none
        real (kind=8) ,dimension(6), intent(in)  :: coordinates ! the angles are in degree
        real (kind=8), parameter :: PII = DACOS(-1.d0)
        real (kind=8) :: cos_b1,cos_b2,cos_c1,cos_c2,sin_b1,sin_b2,sin_c1,sin_c2,cos_phi,sin_phi

        cal_coord_v2(1) = coordinates(1)


        cal_coord_v2(2)  = DCOS(coordinates(2)*PII/180d0)
        cal_coord_v2(3) = DSIN(coordinates(2)*PII/180d0)

        cal_coord_v2(4) = DCOS(coordinates(3)*PII/180d0)
        cal_coord_v2(5) = DSIN(coordinates(3)*PII/180d0)

        cal_coord_v2(6) = DCOS(coordinates(4)*PII/180d0)
        cal_coord_v2(7) = DSIN(coordinates(4)*PII/180d0)

        cal_coord_v2(8) = DCOS(coordinates(5)*PII/180d0)
        cal_coord_v2(9) = DSIN(coordinates(5)*PII/180d0)

        cal_coord_v2(10) = DCOS(coordinates(6)*PII/180d0)
        cal_coord_v2(11) = DSIN(coordinates(6)*PII/180d0)

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

        ar_v2(1) = cos_b1           !Az
        ar_v2(2) = sin_b1*sin_c1    !Ax
        ar_v2(3) = cos_c1*sin_b1    !Ay

        br_v2(1) = -cos_b2           !Bz
        br_v2(2) = -sin_b2*sin_c2    !Bx
        br_v2(3) = -cos_c2*sin_b2    !By

        cc_v2(1)= cos_b1*cos_b2 + cos_phi*sin_b1*sin_b2 !Czz
        cc_v2(2)= cos_c2*sin_phi*sin_b1 + (-cos_phi*cos_b2*sin_b1 + cos_b1*sin_b2)*sin_c2 !Czx
        cc_v2(3)= -cos_phi*cos_b2*cos_c2*sin_b1 + cos_b1*cos_c2*sin_b2 - sin_phi*sin_b1*sin_c2 !Czy


        cc_v2(4)= cos_b2*sin_b1*sin_c1 - sin_b2 *(cos_c1*sin_phi + cos_phi*cos_b1*sin_c1)  !Cxz
        cc_v2(5)= -cos_b1*cos_c2*sin_phi*sin_c1 + (cos_b2*cos_c1*sin_phi + sin_b1*sin_b2*sin_c1)*sin_c2 &
                + cos_phi *(cos_c1*cos_c2 + cos_b1*cos_b2*sin_c1*sin_c2) !Cxx
        cc_v2(6)= cos_c2*sin_b1*sin_b2*sin_c1 + cos_b2*cos_c2 *(cos_c1*sin_phi + cos_phi*cos_b1*sin_c1) &
                + (-cos_phi*cos_c1 + cos_b1*sin_phi*sin_c1)*sin_c2 !Cxy


        cc_v2(7)= cos_b2*cos_c1*sin_b1 + sin_b2 *(-cos_phi*cos_b1*cos_c1 + sin_phi*sin_c1) !Cyz
        cc_v2(8)= cos_c1*sin_b1*sin_b2*sin_c2 + cos_b1*cos_c1 *(-cos_c2*sin_phi + cos_phi*cos_b2*sin_c2) - &
                sin_c1 *(cos_phi*cos_c2 + cos_b2*sin_phi*sin_c2)   !Cyx
        cc_v2(9)= -cos_b2*cos_c2*sin_phi*sin_c1 + cos_c1 *(cos_c2*sin_b1*sin_b2 + cos_b1*sin_phi*sin_c2) &
                + cos_phi*(cos_b1*cos_b2*cos_c1*cos_c2 + sin_c1*sin_c2) !Cyy


    end  subroutine generate_coordenates_v2

    ! ********************  MULTIPOLE INTERACTION  ********************************

    function multipole_sph_v3(ind)
        
        use Fitting_Constant_v2, only: get_coeff_Fit
        
        implicit none

        integer (kind=4), intent(in) :: ind         !coeff index
        real (kind=8) :: multipole_sph_v3           ! returned Energy
        integer (kind=4)::order
        real (kind=8) :: R
        
        R = cal_coord_v2(1)
        multipole_sph_v3 = 0d0
        
        do order = 1, 15

            if ( get_coeff_Fit(ind,order,"M") > 0) then
                multipole_sph_v3 =  multipole_sph_v3 + (C3 * C1 * (C2 ** order)) * multipole_order(ind,order) &
                                                        / R ** order

            end if
        end do

    end function multipole_sph_v3

    function multipole_order(ind,order)

        use Fitting_Constant_v2, only: get_coeff_multipole
        implicit none
        integer (kind=4), intent(in) :: order,ind
        real (kind=8) :: multipole_order
        integer (kind=4)::i,j,ci,cj
        real (kind=8),parameter:: EPS = epsilon(multipole_order)
        real (kind=8):: Qai,Qbj
        real (kind=8),dimension(225)::A_Mult,B_Mult
        multipole_order = 0d0

        A_Mult =  get_coeff_multipole(ind,"A")
        B_Mult =  get_coeff_multipole(ind,"B")

        do i=0,order-1

            j = order-1-i;

            do ci = 0,2*i
                Qai = A_Mult(i**2 +1+ci);
                
                if ( dabs(Qai) > eps ) then
                    do cj = 0,2*j
                        Qbj = B_Mult(j**2 +1+cj);

                        if ( dabs(Qbj) > EPS ) then
                            multipole_order = multipole_order + Qai * Qbj * t_tensor_v2(i+1,ci+1,j+1,cj+1);
                        end if

                    end do
                end if
            end do
        end do

    end function multipole_order

    ! ********************  INDUCTION INTERACTION  ********************************

    function induction_sph_v3(ind)
        !INDUCTION Summary of this function goes here
        !   Detailed explanation goes here
        use Fitting_Constant_v2, only: get_coeff_Fit
        implicit none
        integer(kind=4) , intent(in) :: ind
        real (kind=8) :: induction_sph_v3
        integer(kind=4) :: order


        do order=1,15
            if ( get_coeff_Fit(ind,order,"I") > 0) then
                induction_sph_v3 = induction_sph_v3 + induction_order(order,ind,1)& ! induction of B over A
                                                    + induction_order(order,ind,0)  ! induction of A over B
            end if
        end do
    end function induction_sph_v3

    function induction_order(order,ind,index)
        
        implicit none
        integer (kind=4) , intent(in) :: order,index,ind
        real (kind=8)  :: induction_order
        real (kind=8) :: R
        integer (kind=4) :: l1,l2,i,j,lmin,lmax
        real (kind=8) :: res

        R = cal_coord_v2(1)

        res = 0d0

        do l1=1,order-3
            do l2=1,order-3
                if (l1+l2+2 <= order) then

                    do i=0,order-2-l1-l2
                        do j=0,order-2-l1-l2
                            
                            if (i+j+l1+l2+2 == order)then
                                res  =  res + induction_ij_l1l2(i,j,l1,l2,ind, index)
                            end if
                            
                        end do
                    end do
                end if
            end do
        end do

        induction_order =  ( -0.5d0 * (C3 * C1 * ( C2 ** order )) * res)&
                            /( R ** order)

    end function induction_order

    function induction_ij_l1l2(i,j,l1,l2,ind,index)

        use Fitting_Constant_v2, only:  get_coeff_multipole_by_index, &
                                        get_coeff_polarizability_by_index
        implicit none
        integer (kind=4), intent(in) :: i,j,l1,l2,index,ind
        real (kind=8):: induction_ij_l1l2
        real (kind=8) :: Qai,Qbj,comp_a_k1_k2,T_l1_i,T_l2_j,res
        integer (kind=4) :: ci,cj,k1,k2,cpn,ni,nj,nl1,nl2,lmin,lmax
        real (kind=8), allocatable :: Qa_cpn(:),Qb_cpn(:),pol_arr(:)
        real (kind=8),parameter :: EPS = epsilon(induction_ij_l1l2)

        res = 0d0
        ni = 2*i+1;
        nj = 2*j+1;
        nl1 = 2*l1+1;
        nl2 = 2*l2+1;
        lmin = minval([l1,l2])
        lmax = maxval([l1,l2])

        allocate(Qa_cpn(ni),Qb_cpn(nj),pol_arr(nl1*nl2))

        if ( index == 1 ) then


            Qa_cpn = get_coeff_multipole_by_index(ind,"A",i)
            Qb_cpn = get_coeff_multipole_by_index(ind,"A",j)
            pol_arr = get_coeff_polarizability_by_index(ind,"B",lmin,lmax,nl1*nl2)

        else
            Qa_cpn = get_coeff_multipole_by_index(ind,"B",i)
            Qb_cpn = get_coeff_multipole_by_index(ind,"B",j)
            pol_arr = get_coeff_polarizability_by_index(ind,"A",lmin,lmax,nl1*nl2)
        end if

        do ci = 1,ni
            Qai = Qa_cpn(ci);
            if ( dabs(Qai) > EPS ) then
                do cj = 1,nj
                    Qbj = Qb_cpn(cj);
                    if ( dabs(Qbj) > EPS ) then

                        do k1 = 1,nl1
                            do k2 = 1,nl2

                                cpn = get_induction_cpn(l1,l2,k1,k2)
                                comp_a_k1_k2 = pol_arr(cpn)

                                if ( dabs(comp_a_k1_k2) > EPS ) then
                                    ! index indicate if Im calculating pol over A
                                    ! or pol over B
                                    if ( index == 0 )   then
                                        res = res + Qai * Qbj * comp_a_k1_k2 * &
                                                ( t_tensor_v2(l1+1,k1,i+1,ci)* t_tensor_v2(l2+1,k2,j+1,cj));

                                    else


                                        res = res + Qai * Qbj * comp_a_k1_k2 * &
                                                (t_tensor_v2(i+1,ci,l1+1,k1) * t_tensor_v2(j+1,cj,l2+1,k2));

                                    end if


                                end if

                            end do
                        end do
                    end if

                end do
            end if
        end do

        induction_ij_l1l2 = res

        deallocate(Qa_cpn,Qb_cpn,pol_arr)

    end function induction_ij_l1l2

    function get_induction_cpn(l1,l2,li,lj)
        implicit none
        integer(kind=4), intent(in) :: l1,l2,li,lj
        integer(kind=4) :: get_induction_cpn

        if ( l1 > l2 ) then
            get_induction_cpn = (lj-1) * (2*l1+1) + li;
        else
            get_induction_cpn = (li-1) * (2*l2+1) + lj;
        end if

    end function get_induction_cpn

    !*********************  DISPERSION INTERACTION  *************************************************

    function dispersion_sph_v3(ind)

        !   Dispersion Summary of this function goes here
        !   Detailed explanation goes here
        use Fitting_Constant_v2, only: get_coeff_Fit
        implicit none
        integer (kind=4), intent(in) :: ind
        real (kind=8)  :: dispersion_sph_v3
        integer (kind=4) :: order
        real (kind=8) :: temp

        dispersion_sph_v3  = 0d0
        
        do order=1,15
            if ( get_coeff_Fit(ind,order,"D") > 0) then
                dispersion_sph_v3 = dispersion_sph_v3 + dispersion_order(ind,order)
            end if
        end do

    end function dispersion_sph_v3

    function dispersion_order(ind,order)

        implicit none
        integer (kind=4) , intent(in) :: order,ind
        real (kind=8):: dispersion_order
        integer :: l1,l2,t1,t2
        real (kind=8) :: res,R

        res = 0d0
        R = cal_coord_v2(1)

        do l1=1,order-2
            do l2=1,order-2-l1
                do t1=1,order-2-l1-l2
                    do t2=1,order-2-l1-l2-t1

                        if (l1+l2+t1+t2+2 == order)then

                            res  = res + dispersion_l1l2_t1t2(ind,l1,l2,t1,t2)

                        end if
                    end do
                end do
            end do
        end do

        dispersion_order =  -(( C3 * C1 * ( C2 ** order )) * res)/ (R ** order)

    end function dispersion_order

    function dispersion_l1l2_t1t2(ind,l1,l2,t1,t2)

        use Fitting_Constant_v2, only: get_coeff_dispersion_by_index

        implicit none
        integer (kind=4), intent(in) :: l1,l2,t1,t2,ind
        real (kind=8)  :: dispersion_l1l2_t1t2
        integer (kind=4):: li,lj,ti,tj,cpn
        real (kind=8) :: res,disp_coeff
        real (kind=8) :: disp_arr((2*l1+1)*(2*l2+1)*(2*t1+1)*(2*t2+1))

        real (kind=8) :: EPS = epsilon(res)

        res = 0d0


        disp_arr = get_coeff_dispersion_by_index(ind, l1,l2,t1,t2)

        do li = 0,2*l1
            do lj = 0,2*l2
                do ti = 0,2*t1
                    do tj = 0,2*t2
                        
                        cpn = get_dispersion_cpn(l1,l2,t1,t2,li,lj,ti,tj)
                        disp_coeff = disp_arr(cpn)

                        if ( dabs(disp_coeff) > EPS ) then
                            
                            res = res + disp_coeff * t_tensor_v2(l1+1,li+1,t1+1,ti+1) *&
                                                     t_tensor_v2(l2+1,lj+1,t2+1,tj+1)

                        end if

                    end do
                end do
            end do
        end do

        dispersion_l1l2_t1t2 = res
    end function dispersion_l1l2_t1t2

    function get_dispersion_cpn(l1,l2,t1,t2,li,lj,ti,tj)

        implicit none
        integer (kind=4), intent(in) :: l1,l2,t1,t2,li,lj,ti,tj
        integer (kind=4) :: get_dispersion_cpn

        if (l1>l2 .and. t1>t2) then
            get_dispersion_cpn = lj * (2*l1+1) * (2*t2+1) * (2*t1+1) + li * (2*t2+1) * (2*t1+1) +  tj * (2*t1+1) + ti+1
        elseif (l1>l2 .and. t1<=t2 ) then
            get_dispersion_cpn = lj * (2*l1+1) * (2*t1+1) * (2*t2+1) + li * (2*t1+1) * (2*t2+1) +  ti * (2*t2+1) + tj+1
        elseif (l1<=l2 .and. t1>t2 ) then
            get_dispersion_cpn = li * (2*l2+1) * (2*t2+1) * (2*t1+1) + lj * (2*t2+1) * (2*t1+1) +  tj * (2*t1+1) + ti+1
        else
            get_dispersion_cpn = li * (2*l2+1) * (2*t1+1) * (2*t2+1) + lj * (2*t1+1) * (2*t2+1) +  ti * (2*t2+1)+ tj +1
        end if

    end function get_dispersion_cpn

    !*********************  TOTAL INTERACTION  *************************************************

    function get_total_interaction_energy(coeff_index,general_coordinates_ZXZ)

        use Fitting_Constant_v2, only: get_coeff_max_t_tensor_order

        implicit none
        integer(kind=4), intent(in) :: coeff_index
        real (kind=8) ,dimension(6),intent(in):: general_coordinates_ZXZ
        real (kind=8)  :: get_total_interaction_energy
        integer(kind=4) :: max_t_tensor_order

        max_t_tensor_order = get_coeff_max_t_tensor_order(coeff_index)
        call tensors_initialization_v2(max_t_tensor_order,general_coordinates_ZXZ)

        get_total_interaction_energy =  multipole_sph_v3(coeff_index) +&
                                        induction_sph_v3(coeff_index) +&
                                        dispersion_sph_v3(coeff_index)

    end function get_total_interaction_energy
end module Geometry_Constant_v2