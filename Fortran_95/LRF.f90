!******************************************************************************
!      Compilation Day and Time
!      Month / Day / Year:           11 /           6 /        2024
!      Hr    / Min / Sec :           10 :           4 :           3
!      LRF MATLAB  v4.1.1
!      LRF_Fortran v4.1.1
!******************************************************************************


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

    subroutine init_Tensors_v2(maxlevel,coordinates)
        implicit none
        Integer, INTENT(IN) :: maxlevel
        real*8 ,dimension(6), INTENT(IN)  :: coordinates ! the angles are in degree

        T_Tensor_v2  = 0d0


        call Generate_Coordenates_v2(coordinates)
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

    SUBROUTINE Generate_Coordenates_v2(coordinates)


        real*8 ,dimension(6), INTENT(IN)  :: coordinates ! the angles are in degree
        real*8, parameter :: pii = DACOS(-1.d0)


        real*8 :: cos_b1,cos_b2,cos_c1,cos_c2,sin_b1,sin_b2,sin_c1,sin_c2,cos_phi,sin_phi

        cal_coord_v2(1) = coordinates(1)


        cal_coord_v2(2)  = DCOS(coordinates(2)*pii/180d0)
        cal_coord_v2(3) = DSIN(coordinates(2)*pii/180d0)

        cal_coord_v2(4) = DCOS(coordinates(3)*pii/180d0)
        cal_coord_v2(5) = DSIN(coordinates(3)*pii/180d0)

        cal_coord_v2(6) = DCOS(coordinates(4)*pii/180d0)
        cal_coord_v2(7) = DSIN(coordinates(4)*pii/180d0)

        cal_coord_v2(8) = DCOS(coordinates(5)*pii/180d0)
        cal_coord_v2(9) = DSIN(coordinates(5)*pii/180d0)

        cal_coord_v2(10) = DCOS(coordinates(6)*pii/180d0)
        cal_coord_v2(11) = DSIN(coordinates(6)*pii/180d0)


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

!FittingConstant is the module in charge of all constants related with the fit
!It can handle several coefficient files at the same time


MODULE FitConstants
    save
    public
    real*8, parameter ::  C1=627.5095d0
    real*8, parameter ::  C2=0.529177249d0
    real*8, parameter ::  C3=349.757d0

    INTEGER :: max_T = 0


    TYPE Fit_Contant
        character(:), allocatable :: filename
        real*8 :: Zero
        Integer::initflag

        Integer, dimension (15) :: M_Fit
        Integer, dimension (15) :: D_Fit
        Integer, dimension (15) :: I_Fit


        !Multipoles !

        real*8 , dimension(225)   :: A_Mult,B_Mult

        !Polarizability!

        real*8 , dimension(6,12,195)  :: A_Pol,B_Pol


        !Dispersion!

        real*8 , dimension(5,10,5,10,3087)   :: Disp

    CONTAINS
        PROCEDURE, PASS :: Initializer
        PROCEDURE, PASS :: Read_Parameters

    END TYPE

    Integer,parameter :: NArray = 5
    TYPE(Fit_Contant) :: Coeff(NArray)


CONTAINS

    SUBROUTINE Initializer(this,filename)
        IMPLICIT NONE
        CLASS(Fit_Contant), INTENT(OUT) :: this
        Character(len=*), INTENT(IN) ::filename

        this%initflag = 1
        this%filename = filename

    END SUBROUTINE Initializer

    SUBROUTINE Read_Parameters(this)
        IMPLICIT NONE
        CLASS(Fit_Contant), INTENT(InOut) :: this

        Character(len = 200) :: row
        Integer::iord,mord,dord
        Integer::i,j,l1,l2,t1,t2,ln
        real*8 :: polarr(195)

        if (this%initflag==1)Then


            this%initflag = 2

            Open( 10, file = this%filename )

            Read( 10, *) row
            Read( 10, *) row
            Read( 10, *) row
            Read( 10, *) row
            Read( 10, *) row
            Read( 10, *) row


            Read(10, *)  row,this%Zero
            Read( 10, *) row

            read(10, *)  row,this%M_Fit
            read(10, *)  row,this%I_Fit
            read(10, *)  row,this%D_Fit

            iord = MAXVAL(this%I_Fit);
            mord = MAXVAL(this%M_Fit);
            mord = MAXVAL([mord,iord-3]);
            dord = MAXVAL(this%D_Fit);

            max_T = MAXVAL([max_T,iord-2,mord,dord-3])


            read(10, *)row, this%A_Mult(1:mord**2)  !A_Mult
            read(10, *)row, this%B_Mult(1:mord**2)  !B_Mult



            if (iord>=4) then
                do i=1,iord-3
                    do j=i,iord-3
                        if (i+j<=iord-2) then
                            ln = (2*i+1)*(2*j+1);
                            read(10, *)row, this%A_Pol(i,j,1:ln)  !PA_i-j
                            read(10, *)row,  this%B_Pol(i,j,1:ln) !PB_i-j
                        end if
                    end do
                end do
            end if

            if (dord>=6) then
                do l1=1,dord-5
                    do l2=l1,dord-5
                        do t1=1,dord-5
                            do t2=t1,dord-5

                                if (l1+l2+t1+t2<=dord-2) then
                                    ln = (2*l1+1)*(2*l2+1)*(2*t1+1)*(2*t2+1);
                                    read(10, *)row, this%Disp(l1,l2,t1,t2,1:ln)  !Dispersion l1,l2, t1,t2                                                           disp_kk_vv_coeff{l1,l2,t1,t2});
                                endif
                            enddo
                        enddo
                    enddo
                enddo
            endif

            close(10)


        end if

    END SUBROUTINE Read_Parameters

    SUBROUTINE Find_Coeff_Set(filename,ind)
        IMPLICIT NONE
        Character(*), INTENT(IN) :: filename
        INTEGER, INTENT(OUT) :: ind
        integer:: i
        ind = -1

        if (NArray<1)Then
            Return
        else
            do i=1,NArray
                if (Coeff(i)%filename == filename)then
                    ind = i
                    return
                end if
            end do

        end if

    END SUBROUTINE Find_Coeff_Set

    SUBROUTINE Last_Coeff_Set(lastIndex)
        IMPLICIT NONE
        INTEGER, INTENT(OUT) :: lastIndex
        integer:: i
        lastIndex = 0

        if (NArray<1)Then
            Return
        else

            do i=1,NArray
                if (LEN(Coeff(i)%filename)<1)then
                    lastIndex = i-1
                    return
                end if
            end do

            if (lastIndex==NArray)then
                write(*,*)"The maximun number of coefficients sets is :",NArray,&
                        "to change the maximun go to module MODULE Fit and change NARRAY"
                lastIndex=-10
            end if
        end if


    END SUBROUTINE Last_Coeff_Set

    SUBROUTINE Get_Coeff_Index(filename,indx)
        IMPLICIT NONE
        Character(*), INTENT(IN) :: filename
        Integer, INTENT(OUT) :: indx
        integer::ind,lastIndex

        Call Find_Coeff_Set(filename,ind)


        if (ind < 1)then

            Call Last_Coeff_Set(lastIndex)

            indx = lastIndex + 1


            CALL Coeff(indx)%Initializer(filename)
            Call Coeff(indx)%Read_Parameters()

        else
            indx = ind
            return
        end if

    END SUBROUTINE Get_Coeff_Index

END module FitConstants

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
            Qai = Coeff(ind)%A_Mult(i**2 +1+ci);
            if (DABS(Qai)>eps) then
                do cj = 0,2*j
                    Qbj = Coeff(ind)%B_Mult(j**2 +1+cj);

                    if (DABS(Qbj)>eps) then
                        !print*,"Mult",order,i,j,ci,cj,T_Tensor_v2(i+1,ci+1,j+1,cj+1)
                        E_order = E_order + Qai*Qbj*T_Tensor_v2(i+1,ci+1,j+1,cj+1);
                    end if

                end do
            end if
        end do
    end do




END SUBROUTINE Multipole_Order


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





























SUBROUTINE Coordinate_Transformation(coordinates,coord_format,new_coordinates)
    IMPLICIT NONE
    real*8 ,dimension(6), INTENT(IN):: coordinates

    Character(*), INTENT(IN) :: coord_format
    real*8 , dimension(6), INTENT(INOUT) :: new_coordinates

    new_coordinates=coordinates;


    if (coord_format == "Euler_ZYZ") then
        new_coordinates(5) = coordinates(5) - 90
        new_coordinates(6) = coordinates(6) - 90
    end if
    if (coord_format == "Spherical") then
        new_coordinates(5) = 90 - coordinates(5)
        new_coordinates(6) = 90 -coordinates(6)
    end if

End   SUBROUTINE Coordinate_Transformation

!********************************************************
SUBROUTINE General_Coordinates_Format(Dim, old_coordinates, general_coodinates)
    IMPLICIT NONE

    INTEGER,INTENT(IN) ::  Dim
    real*8 , dimension(6), INTENT(INOUT) ::general_coodinates
    real*8 , dimension (Dim),INTENT(IN) :: old_coordinates


    if (Dim==2) then
        general_coodinates(1) = old_coordinates(1)  !R
        general_coodinates(2) = old_coordinates(2)  !b1
        general_coodinates(3) = 0d0  !b2
        general_coodinates(4) = 0d0  !phi
        general_coodinates(5) = 0d0  !c1
        general_coodinates(6) = 0d0  !c2

    elseif (Dim==3) then
        general_coodinates(1) = old_coordinates(1)  !R
        general_coodinates(2) = old_coordinates(2)  !b1
        general_coodinates(3) = 0d0  !b2
        general_coodinates(4) = 0d0  !phi
        general_coodinates(5) = old_coordinates(3)  !c1
        general_coodinates(6) = 0d0  !c2

    elseif (Dim==4) then
        general_coodinates(1) = old_coordinates(1)  !R
        general_coodinates(2) = old_coordinates(2)  !b1
        general_coodinates(3) = old_coordinates(3)  !b2
        general_coodinates(4) = old_coordinates(4)  !phi
        general_coodinates(5) = 0d0  !c1
        general_coodinates(6) = 0d0  !c2

    elseif (Dim==5) then
        general_coodinates(1) = old_coordinates(1)  !R
        general_coodinates(2) = old_coordinates(2)  !b1
        general_coodinates(3) = old_coordinates(3)  !b2
        general_coodinates(4) = old_coordinates(4)  !phi
        general_coodinates(5) = old_coordinates(5)  !c1
        general_coodinates(6) = 0d0  !c2
    elseif (Dim==6) then
        general_coodinates(1) = old_coordinates(1)  !R
        general_coodinates(2) = old_coordinates(2)  !b1
        general_coodinates(3) = old_coordinates(3)  !b2
        general_coodinates(4) = old_coordinates(4)  !phi
        general_coodinates(5) = old_coordinates(5)  !c1
        general_coodinates(6) = old_coordinates(6)  !c2

    end if

    RETURN
END SUBROUTINE General_Coordinates_Format





SUBROUTINE TotalEnergy_Calc (ind,TotalEnergy)

    implicit none


    Integer, INTENT(IN):: ind ! index of the coefficents
    real*8  , INTENT(INOut) ::TotalEnergy

    real*8   ::EM,ED,EI,term

    call Multipole_Sph3(ind, EM)
    call Induction_Sph3(ind, EI)
    call Dispersion_Sph3(ind, ED)

    !print*, "Electrostatic Energy: ", EM, " Induction Energy: ", EI, " Dispersion Energy: ", ED


    TotalEnergy = EM + ED + EI


end SUBROUTINE TotalEnergy_Calc



! Arg 1 [coordinates] : a coordinate vector [ R , b1, b2, phi] *the angles should be in degrees
! Arg 2 [Coeff_Address] address of the file which contains the longe range expansion coefficients
! Arg 3 [TotalEnergy]   Total Energy calculated


! version 4.0


SUBROUTINE Evaluate_LRF(TotalEnergy,coordinates,coord_format,XDIM,filename)

    use FitConstants, only: Coeff, max_T, Get_Coeff_Index
    use Geometry_Constant_v2, only: init_Tensors_v2

    IMPLICIT NONE

    real*8, INTENT(OUT) :: TotalEnergy
    INTEGER, INTENT(IN) :: XDIM
    real*8 ,dimension(:), INTENT(IN):: coordinates(XDIM)
    Character(*), INTENT(IN) :: coord_format
    Character(*), INTENT(IN) ::  filename
    real*8 ,dimension(6):: GeneralCoordenates,GeneralCoordenates1
    INTEGER :: i,CoeffIndex
    real*8 :: x1

    call Get_Coeff_Index(filename,CoeffIndex)

    x1 = 0d0
    do i=1,XDIM
        x1=x1+dabs(coordinates(i))
    enddo

    if (x1 <= 1d-10) then
        TotalEnergy = Coeff(CoeffIndex)%Zero
        return
    endif

    Call General_Coordinates_Format(XDIM, coordinates, GeneralCoordenates)
    Call Coordinate_Transformation(GeneralCoordenates,coord_format,GeneralCoordenates1)

    Call init_Tensors_v2(max_T,GeneralCoordenates1) ! Initializing in zero the new vectors v2
    Call TotalEnergy_Calc (CoeffIndex,TotalEnergy)

    return

END SUBROUTINE Evaluate_LRF
