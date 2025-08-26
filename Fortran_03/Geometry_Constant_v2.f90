!********************************************************
! Geometry_Constant_v2
!   - Builds T-tensor for a single geometry (serial)
!   - Provides order-level kernels (multipole/induction/dispersion)
!     that read the shared T-tensor and are safe to call from
!     OpenMP parallel loops at a coarse grain (e.g., over 'order').
!********************************************************
module Geometry_Constant_v2
  use, intrinsic :: iso_fortran_env, only: int32, real64
  use Fitting_Constant_v2, only : &
       t_index, ensure_t_index_map_ready, tmap_L, tmap_NNZ, &
       get_coeff_Fit, get_coeff_multipole, get_coeff_multipole_by_index, &
       get_coeff_polarizability_by_index, get_coeff_dispersion_by_index, &
       get_coeff_max_t_tensor_order, &
       ensure_fact_nn_ready, fact_nn
  implicit none
  private

  ! ---- kinds & constants
  integer, parameter :: I4 = int32
  integer, parameter :: RK = real64
  real(RK), parameter :: PI = dacos(-1.0_RK)
  real(RK), parameter :: C1 = 627.5095_RK
  real(RK), parameter :: C2 = 0.529177249_RK
  real(RK), parameter :: C3 = 349.755088236337_RK

  integer(I4), parameter :: LMAX = 15

  ! geometry-dependent helpers
  real(RK)             :: ar_v2(3)
  real(RK)             :: br_v2(3)
  real(RK)             :: cc_v2(9)
  real(RK)             :: cal_coord_v2(11)

  ! Flattened T-tensor (1..NNZ) using global map from Fitting_Constant_v2
  real(RK), allocatable :: t_tensor_v3(:)

  ! ==== Exports ====
  public  :: t_tensor_v3
  public  :: tensors_initialization_v2
  public  :: multipole_sph_v3, induction_sph_v3, dispersion_sph_v3
  public  :: get_total_interaction_energy
  ! Expose order kernels so caller can parallelize over 'order'
  public  :: multipole_order, induction_order, dispersion_order

  ! ==== Local only ====
  private :: calculate_tensor, t_lk_iter
  private :: factorial, factorial_nn, n_eta, coeff_m
  private :: get_splitting_componet, get_tensor_component
  private :: generate_coordenates_v2
  private :: induction_ij_l1l2, get_induction_cpn
  private :: dispersion_l1l2_t1t2, get_dispersion_cpn

contains

  !====================  SETUP (SERIAL)  ====================

  subroutine tensors_initialization_v2(maxlevel, coordinates)
    use Fitting_Constant_v2, only: ensure_t_index_map_ready, ensure_fact_nn_ready
    implicit none
    integer(I4), intent(in) :: maxlevel
    real(RK),    intent(in) :: coordinates(6)   ! angles in degrees

    ! Ensure the global index map is ready (lives in Fitting_Constant_v2)
    call ensure_t_index_map_ready(LMAX)
    call ensure_fact_nn_ready(LMAX)

    if (.not. allocated(t_tensor_v3)) then
      allocate(t_tensor_v3(tmap_NNZ))
    else if (size(t_tensor_v3,1) /= tmap_NNZ) then
      deallocate(t_tensor_v3)
      allocate(t_tensor_v3(tmap_NNZ))
    end if
    t_tensor_v3 = 0.0_RK

    call generate_coordenates_v2(coordinates)

    ! STRICTLY SERIAL build of all T components up to maxlevel
    call calculate_tensor(maxlevel)
  end subroutine tensors_initialization_v2


  subroutine calculate_tensor(maxlevel)
    implicit none
    integer(I4), intent(in) :: maxlevel
    integer(I4) :: order, la, lb, ka_, kb_

    do order = 1, maxlevel
      do la = 0, order-1
        lb = order - la - 1
        do ka_ = 0, 2*la
          do kb_ = 0, 2*lb
            call t_lk_iter(la, ka_, lb, kb_)   ! serial
          end do
        end do
      end do
    end do
  end subroutine calculate_tensor


  subroutine t_lk_iter(la, ka_, lb, kb_)
    implicit none
    integer(I4), intent(in) :: la, ka_, lb, kb_
    real(RK), parameter :: EPS = epsilon(1.0_RK)

    real(RK) :: res, comp_lk, comp_t, prod_comp, fact_prod,fact
    real(RK) :: la_fact, l2_fact, l3_fact, l4_fact
    real(RK) :: lb_fact, la2_fact, fnn, lb2_fact, m1, m2, m, const
    integer(I4) :: ka1, rka1, kb1, rkb1, rk1, rk_, rk_i, rk_j, i, j, n
    character(len=1), parameter :: coord(3) = ["z","x","y"]
    character(len=1) :: rka2, rkb2, ka2, kb2, rk2

    ka2 = get_splitting_componet(ka_)
    kb2 = get_splitting_componet(kb_)
    ka1 = floor((ka_ + 1.0_RK)/2.0_RK)
    kb1 = floor((kb_ + 1.0_RK)/2.0_RK)

    if (la == 0 .and. lb == 0) then
      res = 1.0_RK
    else
      if (lb == 0) then
        comp_lk = 0.0_RK
        la_fact = (2.0_RK*la - 1.0_RK)/(1.0_RK*la)

        do i=1,3
          call n_eta(coord(i), ka1, ka2, rk1, rk2)
          rk_ = get_tensor_component(la, rk1, rk2)
          m   = coeff_m(coord(i), ka1, ka2)
          fact = fact_nn(la, rk1+1, 1, 1)

          if ( dabs(m) > EPS .and. la >= 1 .and. fact > EPS .and. rk_ <= 2*(la-1) ) then
            comp_t    = t_tensor_v3(t_index(la,     rk_+1, 1, 1))   ! la-1+1 = la
            prod_comp = ar_v2(i) * comp_t
            fact_prod = la_fact * m * fact
            comp_lk   = comp_lk + fact_prod * prod_comp
          end if
        end do

        if (la >= 2 .and. ka_ <= 2*(la-2) .and. ka_ >= 0) then
          la2_fact = (la-1.0_RK)/(1.0_RK*la)
          comp_lk = comp_lk - la2_fact * fact_nn(la-1, ka1+1, 1, 1) * &
                    t_tensor_v3(t_index(la-1,  ka_+1, 1, 1))           ! la-2+1 = la-1
        end if

        res = comp_lk / fact_nn(la+1, ka1+1, lb+1, kb1+1)

      elseif (la == 0) then
        comp_lk = 0.0_RK
        lb_fact = (2.0_RK*lb - 1.0_RK)/(1.0_RK*lb)

        do i=1,3
          call n_eta(coord(i), kb1, kb2, rk1, rk2)
          rk_ = get_tensor_component(lb, rk1, rk2)
          m   = coeff_m(coord(i), kb1, kb2)
          fact = fact_nn(1, 1, lb, rk1+1 )

          if ( dabs(m) > EPS .and. lb >= 1 .and. fact > EPS .and. rk_ <= 2*(lb-1) .and. rk_ >= 0) then
            comp_lk = comp_lk + lb_fact * m * fact * br_v2(i) * &
                      t_tensor_v3(t_index(1, 1, lb, rk_+1))         ! lb-1+1 = lb
          end if
        end do

        if (lb >= 2 .and. kb_ <= 2*(lb-2) .and. kb_ >= 0) then
          lb2_fact = (lb-1.0_RK)/(1.0_RK*lb)
          comp_lk = comp_lk - lb2_fact * fact_nn(1, 1, lb-1, kb1+1) * &
                    t_tensor_v3(t_index(1, 1, lb-1, kb_+1))           ! lb-2+1 = lb-1
        end if

        res = comp_lk / fact_nn(la+1, ka1+1, lb+1, kb1+1)

      else
        comp_lk = 0.0_RK

        if (ka_ <= 2*(la-2)) then
          comp_lk = comp_lk + fact_nn(la-1, ka1+1, lb+1, kb1+1) * &
                              t_tensor_v3(t_index(la-1, ka_+1, lb+1, kb_+1))
        end if

        if (kb_ <= 2*(lb-2)) then
          l2_fact = (2.0_RK*la + lb - 1.0_RK)/(1.0_RK*lb)
          comp_lk = comp_lk - (l2_fact * fact_nn(la+1, ka1+1, lb-1, kb1+1)) * &
                              t_tensor_v3(t_index(la+1, ka_+1, lb-1, kb_+1))
        end if

        do i=1,3
          call n_eta(coord(i), kb1, kb2, rk1, rk2)
          rk_i = get_tensor_component(lb, rk1, rk2)
          m    = coeff_m(coord(i), kb1, kb2)
          l3_fact = (2.0_RK*(la + lb) - 1.0_RK)/(lb*1.0_RK)
          const = l3_fact * m * fact_nn(la+1, ka1+1, lb, rk1+1)

          if ( dabs(const) > EPS .and. rk_i <= 2*(lb-1) ) then
            comp_lk = comp_lk + const * br_v2(i) * &
                      t_tensor_v3(t_index(la+1, ka_+1, lb, rk_i+1))  ! lb-1+1 = lb
          end if
        end do

        do i=1,3
          do j=1,3
            n = 3*(i-1) + j
            l4_fact = (2.0_RK*la - 1.0_RK)/(1.0_RK*lb)

            call n_eta(coord(i), ka1, ka2, rka1, rka2)
            call n_eta(coord(j), kb1, kb2, rkb1, rkb2)
            rk_i = get_tensor_component(la, rka1, rka2)
            rk_j = get_tensor_component(lb, rkb1, rkb2)
            m1 = coeff_m(coord(i), ka1, ka2)
            m2 = coeff_m(coord(j), kb1, kb2)

            const = l4_fact * m1 * m2 * fact_nn(la, rka1+1, lb, rkb1+1)
            if ( dabs(const) > EPS .and. rk_i <= 2*(la-1) .and. rk_j <= 2*(lb-1) ) then
              comp_lk = comp_lk + const * cc_v2(n) * &
                        t_tensor_v3(t_index(la, rk_i+1, lb, rk_j+1))  ! la-1+1=la; lb-1+1=lb
            end if
          end do
        end do

        res = comp_lk / fact_nn(la+1, ka1+1, lb+1, kb1+1)
      end if
    end if

    t_tensor_v3(t_index(la+1, ka_+1, lb+1, kb_+1)) = res
  end subroutine t_lk_iter

  !====================  SMALL HELPERS  ====================

  function factorial(n) result(f)
    implicit none
    integer(I4), intent(in) :: n
    real(RK) :: f
    integer(I4) :: i

    if (n >= 0) then
      f = 1.0_RK
      do i = 2, n
        f = f * real(i, RK)
      end do
    else
      f = 0.0_RK
    end if
  end function factorial

  function factorial_nn(la, ka1, lb, kb1) result(fn)
    implicit none
    integer(I4), intent(in) :: la, ka1, lb, kb1
    real(RK) :: fn
    if (la < 0 .or. lb < 0 .or. ka1 < 0 .or. kb1 < 0 .or. ka1 > la .or. kb1 > lb) then
      fn = 0.0_RK
    else
      fn = dsqrt( ( factorial(la + ka1)/factorial(la - ka1) ) * &
                  ( factorial(lb + kb1)/factorial(lb - kb1) ) )
    end if
  end function factorial_nn

  subroutine n_eta(mu, k1, k2, ka1, ka2)
    implicit none
    character(len=1), intent(in)  :: mu, k2
    integer(I4),      intent(in)  :: k1
    character(len=1), intent(out) :: ka2
    integer(I4),      intent(out) :: ka1
    if (mu == "x") then
      if (k1 <= 1) then
        ka1 = 0; ka2 = "0"
      else
        ka1 = k1 - 1; ka2 = k2
      end if
    elseif (mu == "y") then
      if (k1 <= 1) then
        ka1 = 0; ka2 = "0"
      else
        if (k2 == "c") then
          ka1 = k1 - 1; ka2 = "s"
        else
          ka1 = k1 - 1; ka2 = "c"
        end if
      end if
    else
      ka1 = k1; ka2 = k2
    end if
  end subroutine n_eta

  function coeff_m(mu, k1, k2) result(cm)
    implicit none
    character(len=1), intent(in) :: mu, k2
    integer(I4),      intent(in) :: k1
    real(RK) :: cm
    cm = 0.0_RK
    if (mu == "x") then
      if (k1 == 1) then
        if (k2 == "c") cm = dsqrt(2.0_RK)
      else
        cm = real(k1, RK)
      end if
    elseif (mu == "y") then
      if (k1 == 1) then
        if (k2 == "s") cm = dsqrt(2.0_RK)
      else
        if (k2 == "s") then
          cm = real(k1, RK)
        else
          cm = -real(k1, RK)
        end if
      end if
    else
      cm = 1.0_RK
    end if
  end function coeff_m

  function get_splitting_componet(i) result(sp)
    implicit none
    integer(I4), intent(in) :: i
    character(len=1) :: sp
    if (i >= 0) then
      if (i == 0) then
        sp = "0"
      elseif (mod(i,2) == 1) then
        sp = "c"
      else
        sp = "s"
      end if
    else
      sp = "-"
    end if
  end function get_splitting_componet

  function get_tensor_component(mult_ord, k1, k2) result(ic)
    implicit none
    integer(I4), intent(in) :: mult_ord, k1
    character(len=1), intent(in) :: k2
    integer(I4) :: ic
    if (k1 < 0 .or. mult_ord < 0 .or. k1 > mult_ord) then
      ic = -1
    elseif (k1 == 0) then
      if (k2 == "0") then
        ic = 0
      else
        ic = -1
      end if
    else
      if (k2 == "s") then
        ic = 2*k1
      elseif (k2 == "c") then
        ic = 2*k1 - 1
      else
        ic = 0
      end if
    end if
  end function get_tensor_component

  subroutine generate_coordenates_v2(coordinates)
    implicit none
    real(RK), intent(in) :: coordinates(6) ! angles in degree
    real(RK), parameter :: PII = dacos(-1.0_RK)
    real(RK) :: cos_b1,cos_b2,cos_c1,cos_c2,sin_b1,sin_b2,sin_c1,sin_c2,cos_phi,sin_phi

    cal_coord_v2(1)  = coordinates(1)

    cal_coord_v2(2)  = dcos(coordinates(2)*PII/180.0_RK)
    cal_coord_v2(3)  = dsin(coordinates(2)*PII/180.0_RK)

    cal_coord_v2(4)  = dcos(coordinates(3)*PII/180.0_RK)
    cal_coord_v2(5)  = dsin(coordinates(3)*PII/180.0_RK)

    cal_coord_v2(6)  = dcos(coordinates(4)*PII/180.0_RK)
    cal_coord_v2(7)  = dsin(coordinates(4)*PII/180.0_RK)

    cal_coord_v2(8)  = dcos(coordinates(5)*PII/180.0_RK)
    cal_coord_v2(9)  = dsin(coordinates(5)*PII/180.0_RK)

    cal_coord_v2(10) = dcos(coordinates(6)*PII/180.0_RK)
    cal_coord_v2(11) = dsin(coordinates(6)*PII/180.0_RK)

    cos_b1 = cal_coord_v2(2);  sin_b1 = cal_coord_v2(3)
    cos_b2 = cal_coord_v2(4);  sin_b2 = cal_coord_v2(5)
    cos_phi= cal_coord_v2(6);  sin_phi= cal_coord_v2(7)
    cos_c1 = cal_coord_v2(8);  sin_c1 = cal_coord_v2(9)
    cos_c2 = cal_coord_v2(10); sin_c2 = cal_coord_v2(11)

    ar_v2(1) = cos_b1
    ar_v2(2) = sin_b1*sin_c1
    ar_v2(3) = cos_c1*sin_b1

    br_v2(1) = -cos_b2
    br_v2(2) = -sin_b2*sin_c2
    br_v2(3) = -cos_c2*sin_b2

    cc_v2(1)=  cos_b1*cos_b2 + cos_phi*sin_b1*sin_b2
    cc_v2(2)=  cos_c2*sin_phi*sin_b1 + (-cos_phi*cos_b2*sin_b1 + cos_b1*sin_b2)*sin_c2
    cc_v2(3)= -cos_phi*cos_b2*cos_c2*sin_b1 + cos_b1*cos_c2*sin_b2 - sin_phi*sin_b1*sin_c2

    cc_v2(4)=  cos_b2*sin_b1*sin_c1 - sin_b2*(cos_c1*sin_phi + cos_phi*cos_b1*sin_c1)
    cc_v2(5)= -cos_b1*cos_c2*sin_phi*sin_c1 + (cos_b2*cos_c1*sin_phi + sin_b1*sin_b2*sin_c1)*sin_c2 &
              + cos_phi*(cos_c1*cos_c2 + cos_b1*cos_b2*sin_c1*sin_c2)
    cc_v2(6)=  cos_c2*sin_b1*sin_b2*sin_c1 + cos_b2*cos_c2*(cos_c1*sin_phi + cos_phi*cos_b1*sin_c1) &
              + (-cos_phi*cos_c1 + cos_b1*sin_phi*sin_c1)*sin_c2

    cc_v2(7)=  cos_b2*cos_c1*sin_b1 + sin_b2*(-cos_phi*cos_b1*cos_c1 + sin_phi*sin_c1)
    cc_v2(8)=  cos_c1*sin_b1*sin_b2*sin_c2 + cos_b1*cos_c1*(-cos_c2*sin_phi + cos_phi*cos_b2*sin_c2) - &
               sin_c1*(cos_phi*cos_c2 + cos_b2*sin_phi*sin_c2)
    cc_v2(9)= -cos_b2*cos_c2*sin_phi*sin_c1 + cos_c1*(cos_c2*sin_b1*sin_b2 + cos_b1*sin_phi*sin_c2) &
              + cos_phi*(cos_b1*cos_b2*cos_c1*cos_c2 + sin_c1*sin_c2)
  end  subroutine generate_coordenates_v2

  !====================  MULTIPOLE  ====================

  
  function multipole_sph_v3(ind) result(val)
    implicit none
    integer(I4), intent(in) :: ind
    real(RK) :: val
    integer(I4) :: order
    real(RK) :: R
    R = cal_coord_v2(1)
    val = 0.0_RK

    
    do order = 1, LMAX
        if ( get_coeff_Fit(ind,order,"M") > 0 ) then
        val = val + (C3*C1*(C2**order))*multipole_order(ind,order) / (R**order)
        end if
    end do
  
  end function multipole_sph_v3

  function multipole_order(ind,order) result(mo)
    implicit none
    integer(I4), intent(in) :: order, ind
    real(RK) :: mo
    integer(I4) :: i, j, ci, cj
    real(RK), parameter :: EPS = epsilon(1.0_RK)
    real(RK) :: Qai, Qbj
    real(RK) :: A_Mult(225), B_Mult(225)

    mo = 0.0_RK
    A_Mult = get_coeff_multipole(ind,"A")
    B_Mult = get_coeff_multipole(ind,"B")

    

    do i = 0, order-1
      j = order-1-i
      do ci = 0, 2*i
        Qai = A_Mult(i**2 + 1 + ci)
        if ( dabs(Qai) > EPS ) then
          do cj = 0, 2*j
            Qbj = B_Mult(j**2 + 1 + cj)
            if ( dabs(Qbj) > EPS ) then
              mo = mo + Qai*Qbj*t_tensor_v3( t_index(i+1,ci+1,j+1,cj+1) )
            end if
          end do
        end if
      end do
    end do
  end function multipole_order

  !====================  INDUCTION  ====================
  function induction_sph_v3(ind) result(val)
    implicit none
    integer(I4), intent(in) :: ind
    real(RK) :: val
    integer(I4) :: order
    real(RK) :: R
    R = cal_coord_v2(1)
    val = 0.0_RK

    !$omp parallel do reduction(+:val) schedule(static)
    do order=1, LMAX
        if ( get_coeff_Fit(ind,order,"I") > 0 ) then
        val = val + (-0.5_RK*(C3*C1*(C2**order)) *  &
                ( induction_order(order,ind,1) + induction_order(order,ind,0) )) / (R**order)
        end if
    end do
    !$omp end parallel do
  end function induction_sph_v3

  function induction_order(order,ind,index) result(io)
    implicit none
    integer(I4), intent(in) :: order, ind, index
    real(RK) :: io
    integer(I4) :: l1, l2, i, j
    real(RK) :: res
    res = 0.0_RK

    do l1=1,order-3
      do l2=1,order-3
        if (l1 + l2 + 2 <= order) then
          do i=0,order-2-l1-l2
            do j=0,order-2-l1-l2
              if (i + j + l1 + l2 + 2 == order) then
                res = res + induction_ij_l1l2(i,j,l1,l2,ind,index)
              end if
            end do
          end do
        end if
      end do
    end do
    io = res
  end function induction_order

  function induction_ij_l1l2(i,j,l1,l2,ind,index) result(val)
    implicit none
    integer(I4), intent(in) :: i, j, l1, l2, index, ind
    real(RK) :: val
    integer(I4) :: ci, cj, k1, k2, cpn, ni, nj, nl1, nl2, lmin, lmax
    real(RK), allocatable :: Qa_cpn(:), Qb_cpn(:), pol_arr(:)
    real(RK) :: tmpA(225), tmpB(225)
    real(RK), parameter :: EPS = epsilon(1.0_RK)
    real(RK) :: Qai, Qbj, comp_a_k1_k2

    val = 0.0_RK
    ni = 2*i + 1
    nj = 2*j + 1
    nl1 = 2*l1 + 1
    nl2 = 2*l2 + 1
    lmin = minval([l1,l2])
    lmax = maxval([l1,l2])

    allocate(Qa_cpn(ni), Qb_cpn(nj), pol_arr(nl1*nl2))

    if ( index == 1 ) then
      tmpA = get_coeff_multipole_by_index(ind,"A",i)
      tmpB = get_coeff_multipole_by_index(ind,"A",j)
      Qa_cpn(1:ni) = tmpA(1:ni)
      Qb_cpn(1:nj) = tmpB(1:nj)
      pol_arr(1:nl1*nl2) = get_coeff_polarizability_by_index(ind,"B",lmin,lmax,nl1*nl2)
    else
      tmpA = get_coeff_multipole_by_index(ind,"B",i)
      tmpB = get_coeff_multipole_by_index(ind,"B",j)
      Qa_cpn(1:ni) = tmpA(1:ni)
      Qb_cpn(1:nj) = tmpB(1:nj)
      pol_arr(1:nl1*nl2) = get_coeff_polarizability_by_index(ind,"A",lmin,lmax,nl1*nl2)
    end if

    do ci = 1, ni
      Qai = Qa_cpn(ci)
      if ( dabs(Qai) > EPS ) then
        do cj = 1, nj
          Qbj = Qb_cpn(cj)
          if ( dabs(Qbj) > EPS ) then
            do k1 = 1, nl1
              do k2 = 1, nl2
                cpn = get_induction_cpn(l1,l2,k1,k2)
                comp_a_k1_k2 = pol_arr(cpn)
                if ( dabs(comp_a_k1_k2) > EPS ) then
                  if ( index == 0 ) then
                    val = val + Qai*Qbj*comp_a_k1_k2 * &
                          ( t_tensor_v3(t_index(l1+1,k1,i+1,ci)) * t_tensor_v3(t_index(l2+1,k2,j+1,cj)) )
                  else
                    val = val + Qai*Qbj*comp_a_k1_k2 * &
                          ( t_tensor_v3(t_index(i+1,ci,l1+1,k1)) * t_tensor_v3(t_index(j+1,cj,l2+1,k2)) )
                  end if
                end if
              end do
            end do
          end if
        end do
      end if
    end do

    deallocate(Qa_cpn, Qb_cpn, pol_arr)
  end function induction_ij_l1l2

  function get_induction_cpn(l1,l2,li,lj) result(ind)
    implicit none
    integer(I4), intent(in) :: l1,l2,li,lj
    integer(I4) :: ind
    if ( l1 > l2 ) then
      ind = (lj-1) * (2*l1+1) + li
    else
      ind = (li-1) * (2*l2+1) + lj
    end if
  end function get_induction_cpn

  !====================  DISPERSION  ====================
  function dispersion_sph_v3(ind) result(val)
    implicit none
    integer(I4), intent(in) :: ind
    real(RK) :: val
    integer(I4) :: order
    real(RK) :: R
    R = cal_coord_v2(1)
    val = 0.0_RK

    !$omp parallel do reduction(+:val) schedule(static)
    do order=1, LMAX
        if ( get_coeff_Fit(ind,order,"D") > 0 ) then
        val = val + ( - (C3*C1*(C2**order)) * dispersion_order(ind,order) ) / (R**order)
        end if
    end do
    !$omp end parallel do
  end function dispersion_sph_v3

  function dispersion_order(ind,order) result(do_val)
    implicit none
    integer(I4), intent(in) :: order, ind
    real(RK) :: do_val
    integer(I4) :: l1, l2, t1, t2
    real(RK) :: res

    res = 0.0_RK
    do l1=1,order-2
      do l2=1,order-2-l1
        do t1=1,order-2-l1-l2
          do t2=1,order-2-l1-l2-t1
            if (l1 + l2 + t1 + t2 + 2 == order) then
              res = res + dispersion_l1l2_t1t2(ind,l1,l2,t1,t2)
            end if
          end do
        end do
      end do
    end do
    do_val = res
  end function dispersion_order

  function dispersion_l1l2_t1t2(ind,l1,l2,t1,t2) result(val)
    implicit none
    integer(I4), intent(in) :: l1,l2,t1,t2,ind
    real(RK) :: val
    integer(I4) :: li, lj, ti, tj, cpn
    real(RK) :: res, disp_coeff
    real(RK) :: disp_arr((2*l1+1)*(2*l2+1)*(2*t1+1)*(2*t2+1))
    real(RK), parameter :: EPS = epsilon(1.0_RK)

    res = 0.0_RK
    disp_arr = get_coeff_dispersion_by_index(ind, l1,l2,t1,t2)

    do li = 0,2*l1
      do lj = 0,2*l2
        do ti = 0,2*t1
          do tj = 0,2*t2
            cpn = get_dispersion_cpn(l1,l2,t1,t2,li,lj,ti,tj)
            disp_coeff = disp_arr(cpn)
            if ( dabs(disp_coeff) > EPS ) then
              res = res + disp_coeff * t_tensor_v3(t_index(l1+1,li+1,t1+1,ti+1)) * &
                                   t_tensor_v3(t_index(l2+1,lj+1,t2+1,tj+1))
            end if
          end do
        end do
      end do
    end do
    val = res
  end function dispersion_l1l2_t1t2

  function get_dispersion_cpn(l1,l2,t1,t2,li,lj,ti,tj) result(ind)
    implicit none
    integer(I4), intent(in) :: l1,l2,t1,t2,li,lj,ti,tj
    integer(I4) :: ind
    if (l1>l2 .and. t1>t2) then
      ind = lj*(2*l1+1)*(2*t2+1)*(2*t1+1) + li*(2*t2+1)*(2*t1+1) + tj*(2*t1+1) + ti + 1
    elseif (l1>l2 .and. t1<=t2 ) then
      ind = lj*(2*l1+1)*(2*t1+1)*(2*t2+1) + li*(2*t1+1)*(2*t2+1) + ti*(2*t2+1) + tj + 1
    elseif (l1<=l2 .and. t1>t2 ) then
      ind = li*(2*l2+1)*(2*t2+1)*(2*t1+1) + lj*(2*t2+1)*(2*t1+1) + tj*(2*t1+1) + ti + 1
    else
      ind = li*(2*l2+1)*(2*t1+1)*(2*t2+1) + lj*(2*t1+1)*(2*t2+1) + ti*(2*t2+1) + tj + 1
    end if
  end function get_dispersion_cpn

  !====================  TOTAL  ====================

  function get_total_interaction_energy(coeff_index,general_coordinates_ZXZ) result(E)
    implicit none
    integer(I4), intent(in) :: coeff_index
    real(RK),    intent(in) :: general_coordinates_ZXZ(6)
    real(RK) :: E
    integer(I4) :: max_t_tensor_order

    max_t_tensor_order = get_coeff_max_t_tensor_order(coeff_index)
    call tensors_initialization_v2(max_t_tensor_order, general_coordinates_ZXZ)

    E = multipole_sph_v3(coeff_index) + induction_sph_v3(coeff_index) + dispersion_sph_v3(coeff_index)
  end function get_total_interaction_energy

end module Geometry_Constant_v2
