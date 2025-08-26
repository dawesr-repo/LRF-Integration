!********************************************************
module Geometry_Constant_v2
  use, intrinsic :: iso_fortran_env, only: int32, real64
  implicit none
  private

  !-------------------- kinds & constants --------------------
  integer, parameter :: I4 = int32
  integer, parameter :: RK = real64
  real(RK), parameter :: PI  = dacos(-1.0d0)
  real(RK), parameter :: DEG = PI/180.0d0
  real(RK), parameter :: C1  = 627.5095_RK
  real(RK), parameter :: C2  = 0.529177249_RK
  real(RK), parameter :: C3  = 349.755088236337_RK
  integer(I4), parameter :: MAX_L = 15_I4
  real(RK), parameter :: EPS = 1.0e-14_RK

  ! axis tags to avoid string traffic: index order matches ar_v2/br_v2 layout
  integer(I4), parameter :: MU_Z=1_I4, MU_X=2_I4, MU_Y=3_I4
  integer(I4), parameter :: K0=0_I4, KC=1_I4, KS=2_I4   ! spherical split tags: 0, cos, sin

  !-------------------- module storage -----------------------
  real(RK) , dimension(3)              :: ar_v2
  real(RK) , dimension(3)              :: br_v2
  real(RK) , dimension(9)              :: cc_v2
  real(RK) , dimension(11)             :: cal_coord_v2
  real(RK) :: t_tensor_v2(MAX_L, 2*MAX_L+1, MAX_L, 2*MAX_L+1)

  !-------------------- visibility ---------------------------
  public :: t_tensor_v2
  public :: tensors_initialization_v2
  public :: multipole_sph_v3, induction_sph_v3, dispersion_sph_v3
  public :: get_total_interaction_energy

  private :: calculate_tensor, t_lk_iter
  private :: nn_coeff, n_eta, coeff_m, get_splitting_component, get_tensor_component
  private :: generate_coordinates_v2
  private :: multipole_order
  private :: induction_order, induction_ij_l1l2, get_induction_cpn
  private :: dispersion_order, dispersion_l1l2_t1t2, get_dispersion_cpn

contains

  !============================================================
  subroutine tensors_initialization_v2(maxlevel, coordinates)
    !! Initialize coordinates and build all T-tensors up to maxlevel
    integer(I4), intent(in) :: maxlevel
    real(RK)   , intent(in) :: coordinates(6)  ! R, beta1, gamma1, phi, gamma2, beta2 (angles in degrees)

    if (maxlevel < 1 .or. maxlevel > MAX_L) stop 'tensors_initialization_v2: invalid maxlevel'

    t_tensor_v2 = 0.0_RK
    call generate_coordinates_v2(coordinates)
    call calculate_tensor(maxlevel)
  end subroutine tensors_initialization_v2

  !============================================================
  subroutine calculate_tensor(maxlevel)
    !! Fill t_tensor_v2(la,ka,lb,kb) for 0<=la+lb<maxlevel (internally stored 1-based)
    integer(I4), intent(in) :: maxlevel
    integer(I4) :: order, la, ka_, lb, kb_

    do order = 1, maxlevel
      do la = 0, order-1
        lb = order - la - 1
        do ka_ = 0, 2*la
          do kb_ = 0, 2*lb
            call t_lk_iter(la, ka_, lb, kb_)
          end do
        end do
      end do
    end do
  end subroutine calculate_tensor

  !============================================================
  subroutine t_lk_iter(la, ka_, lb, kb_)
    !! Bottom-up evaluation of T-tensor via recursion
    integer(I4), intent(in) :: la, ka_, lb, kb_

    real(RK) :: res, comp_lk, comp_t, prod_comp, fact_prod
    real(RK) :: la_fact, l2_fact, l3_fact, l4_fact, fact_nn, const
    real(RK) :: m, m1, m2
    integer(I4) :: ka1, kb1, rk1, rk2, rka1, rka2, rkb1, rkb2, i, j, n
    integer(I4), parameter :: axis_list(3) = [MU_Z, MU_X, MU_Y]
    integer(I4) :: ka2, kb2, rk_, rk_i, rk_j

    ka2 = get_splitting_component(ka_)
    kb2 = get_splitting_component(kb_)
    ka1 = (ka_ + 1) / 2
    kb1 = (kb_ + 1) / 2

    if (la == 0 .and. lb == 0) then
      res = 1.0_RK

    else if (lb == 0) then
      ! ---- recurrence for lb = 0 ----
      comp_lk = 0.0_RK
      la_fact = (2.0_RK*la - 1.0_RK) / real(la, RK)

      do i = 1, 3
        call n_eta(axis_list(i), ka1, ka2, rk1, rk2)
        rk_ = get_tensor_component(la, rk1, rk2)
        m   = coeff_m(axis_list(i), ka1, ka2)
        fact_nn = nn_coeff(la-1, rk1, 0, 0)

        if (dabs(m) > EPS .and. la >= 1 .and. fact_nn > EPS .and. rk_ <= 2*(la-1)) then
          comp_t    = t_tensor_v2( (la-1)+1, rk_+1, 1, 1 )
          prod_comp = ar_v2(i) * comp_t
          fact_prod = la_fact * m * fact_nn
          comp_lk   = comp_lk + fact_prod * prod_comp
        end if
      end do

      if (la >= 2 .and. ka_ <= 2*(la-2) .and. ka_ >= 0) then
        comp_lk = comp_lk - ((la-1.0_RK)/real(la,RK)) * nn_coeff(la-2, ka1, 0, 0) * &
                             t_tensor_v2( (la-2)+1, ka_+1, 1, 1 )
      end if

      res = comp_lk / nn_coeff(la, ka1, lb, kb1)

    else if (la == 0) then
      ! ---- recurrence for la = 0 ----
      comp_lk = 0.0_RK
      ! lb_fact used implicitly inside constants below as in original form

      do i = 1, 3
        call n_eta(axis_list(i), kb1, kb2, rk1, rk2)
        rk_ = get_tensor_component(lb, rk1, rk2)
        m   = coeff_m(axis_list(i), kb1, kb2)
        fact_nn = nn_coeff(0, 0, lb-1, rk1)

        if (dabs(m) > EPS .and. lb >= 1 .and. fact_nn > EPS .and. rk_ <= 2*(lb-1) .and. rk_ >= 0) then
          comp_lk = comp_lk + ((2.0_RK*lb-1.0_RK)/real(lb,RK)) * m * fact_nn * br_v2(i) * &
                              t_tensor_v2( 1, 1, (lb-1)+1, rk_+1 )
        end if
      end do

      if (lb >= 2 .and. kb_ <= 2*(lb-2) .and. kb_ >= 0) then
        comp_lk = comp_lk - ((lb-1.0_RK)/real(lb,RK)) * nn_coeff(0, 0, lb-2, kb1) * &
                             t_tensor_v2( 1, 1, (lb-2)+1, kb_+1 )
      end if

      res = comp_lk / nn_coeff(la, ka1, lb, kb1)

    else
      ! ---- general recurrence la>0 & lb>0 ----
      comp_lk = 0.0_RK

      if (ka_ <= 2*(la-2)) then
        comp_lk = comp_lk + nn_coeff(la-2, ka1, lb,   kb1) * t_tensor_v2( (la-2)+1, ka_+1, lb+1,    kb_+1 )
      end if

      if (kb_ <= 2*(lb-2)) then
        l2_fact = (2.0_RK*la + lb - 1.0_RK)/real(lb, RK)
        comp_lk = comp_lk - l2_fact * nn_coeff(la,   ka1, lb-2, kb1) * t_tensor_v2( la+1,    ka_+1, (lb-2)+1, kb_+1 )
      end if

      l3_fact = (2.0_RK*(la+lb) - 1.0_RK)/real(lb, RK)
      do i = 1, 3
        call n_eta(axis_list(i), kb1, kb2, rk1, rk2)
        rk_i = get_tensor_component(lb, rk1, rk2)
        m    = coeff_m(axis_list(i), kb1, kb2)
        const = l3_fact * m * nn_coeff(la, ka1, lb-1, rk1)
        if (dabs(const) > EPS .and. rk_i <= 2*(lb-1)) then
          comp_lk = comp_lk + const * br_v2(i) * t_tensor_v2( la+1, ka_+1, (lb-1)+1, rk_i+1 )
        end if
      end do

      l4_fact = (2.0_RK*la - 1.0_RK)/real(lb, RK)
      do i = 1, 3
        do j = 1, 3
          n = 3*(i-1) + j
          call n_eta(axis_list(i), ka1, ka2, rka1, rka2)
          call n_eta(axis_list(j), kb1, kb2, rkb1, rkb2)
          rk_i = get_tensor_component(la, rka1, rka2)
          rk_j = get_tensor_component(lb, rkb1, rkb2)
          m1   = coeff_m(axis_list(i), ka1, ka2)
          m2   = coeff_m(axis_list(j), kb1, kb2)
          const = l4_fact * m1 * m2 * nn_coeff(la-1, rka1, lb-1, rkb1)
          if (dabs(const) > EPS .and. rk_i <= 2*(la-1) .and. rk_j <= 2*(lb-1)) then
            comp_lk = comp_lk + const * cc_v2(n) * t_tensor_v2( (la-1)+1, rk_i+1, (lb-1)+1, rk_j+1 )
          end if
        end do
      end do

      res = comp_lk / nn_coeff(la, ka1, lb, kb1)
    end if

    t_tensor_v2( la+1, ka_+1, lb+1, kb_+1 ) = res
  end subroutine t_lk_iter

  !============================================================
  pure function nn_coeff(la, ka1, lb, kb1) result(val)
    !! Stable evaluation of the NN coefficient using log-gamma ratios
    integer(I4), intent(in) :: la, ka1, lb, kb1
    real(RK) :: val

    if (la<0 .or. lb<0 .or. ka1<0 .or. kb1<0 .or. ka1>la .or. kb1>lb) then
      val = 0.0_RK
    else
      val = dexp( 0.5_RK * ( log_gamma(real(la+ka1+1, RK)) - log_gamma(real(la-ka1+1, RK)) &
                          + log_gamma(real(lb+kb1+1, RK)) - log_gamma(real(lb-kb1+1, RK)) ) )
    end if
  end function nn_coeff

  !============================================================
  pure subroutine n_eta(mu, k1, k2, ka1, ka2)
    !! Helper mapping for recursion: (mu,k1,k2) -> (ka1,ka2)
    integer(I4), intent(in)  :: mu, k1, k2   ! mu in {MU_Z,MU_X,MU_Y}, k2 in {K0,KC,KS}
    integer(I4), intent(out) :: ka1, ka2

    select case (mu)
    case (MU_X)
      if (k1 <= 1) then
        ka1 = 0; ka2 = K0
      else
        ka1 = k1 - 1; ka2 = k2
      end if
    case (MU_Y)
      if (k1 <= 1) then
        ka1 = 0; ka2 = K0
      else
        ka1 = k1 - 1
        if (k2 == KC) then
          ka2 = KS
        else
          ka2 = KC
        end if
      end if
    case default  ! MU_Z
      ka1 = k1; ka2 = k2
    end select
  end subroutine n_eta

  !============================================================
  pure function coeff_m(mu, k1, k2) result(val)
    !! M-coefficient in the T-tensor recursion
    integer(I4), intent(in) :: mu, k1, k2
    real(RK) :: val

    val = 0.0_RK
    select case (mu)
    case (MU_X)
      if (k1 == 1) then
        if (k2 == KC) val = dsqrt(2.0_RK)
      else
        val = real(k1, RK)
      end if
    case (MU_Y)
      if (k1 == 1) then
        if (k2 == KS) val = dsqrt(2.0_RK)
      else
        if (k2 == KS) then
          val = real(k1, RK)
        else
          val = -real(k1, RK)
        end if
      end if
    case default  ! MU_Z
      val = 1.0_RK
    end select
  end function coeff_m

  !============================================================
  pure function get_splitting_component(i) result(tag)
    !! Map tensor index i -> spherical split tag (K0/KC/KS)
    integer(I4), intent(in) :: i
    integer(I4) :: tag

    if (i < 0) then
      tag = -1
    else if (i == 0) then
      tag = K0
    else if (mod(i,2) == 1) then
      tag = KC
    else
      tag = KS
    end if
  end function get_splitting_component

  !============================================================
  pure function get_tensor_component(mult_ord, k1, k2) result(idx)
    !! Return packed tensor component index (0-based) or -1 if invalid
    integer(I4), intent(in) :: mult_ord, k1, k2
    integer(I4) :: idx

    if (k1 < 0 .or. mult_ord < 0 .or. k1 > mult_ord) then
      idx = -1
    else if (k1 == 0) then
      if (k2 == K0) then
        idx = 0
      else
        idx = -1
      end if
    else
      select case (k2)
      case (KS)
        idx = 2*k1
      case (KC)
        idx = 2*k1 - 1
      case default
        idx = 0
      end select
    end if
  end function get_tensor_component

  !============================================================
  subroutine generate_coordinates_v2(coordinates)
    !! Precompute direction cosines and coupling terms from ZXZ Euler angles (degrees)
    real(RK), intent(in) :: coordinates(6)
    real(RK) :: cos_b1,cos_b2,cos_c1,cos_c2,sin_b1,sin_b2,sin_c1,sin_c2,cos_phi,sin_phi

    cal_coord_v2(1) = coordinates(1)

    cal_coord_v2(2)  = dcos(coordinates(2)*DEG)
    cal_coord_v2(3)  = dsin(coordinates(2)*DEG)

    cal_coord_v2(4)  = dcos(coordinates(3)*DEG)
    cal_coord_v2(5)  = dsin(coordinates(3)*DEG)

    cal_coord_v2(6)  = dcos(coordinates(4)*DEG)
    cal_coord_v2(7)  = dsin(coordinates(4)*DEG)

    cal_coord_v2(8)  = dcos(coordinates(5)*DEG)
    cal_coord_v2(9)  = dsin(coordinates(5)*DEG)

    cal_coord_v2(10) = dcos(coordinates(6)*DEG)
    cal_coord_v2(11) = dsin(coordinates(6)*DEG)

    cos_b1 = cal_coord_v2(2);  sin_b1 = cal_coord_v2(3)
    cos_b2 = cal_coord_v2(4);  sin_b2 = cal_coord_v2(5)
    cos_phi= cal_coord_v2(6);  sin_phi= cal_coord_v2(7)
    cos_c1 = cal_coord_v2(8);  sin_c1 = cal_coord_v2(9)
    cos_c2 = cal_coord_v2(10); sin_c2 = cal_coord_v2(11)

    ! ar_v2: A orientation (index mapping: 1->z, 2->x, 3->y)
    ar_v2(1) =  cos_b1             ! Az
    ar_v2(2) =  sin_b1*sin_c1      ! Ax
    ar_v2(3) =  cos_c1*sin_b1      ! Ay

    ! br_v2: B orientation (negative as in original)
    br_v2(1) = -cos_b2             ! Bz
    br_v2(2) = -sin_b2*sin_c2      ! Bx
    br_v2(3) = -cos_c2*sin_b2      ! By

    ! cc_v2: rotation-coupling terms
    cc_v2(1) =  cos_b1*cos_b2 + cos_phi*sin_b1*sin_b2                         ! Czz
    cc_v2(2) =  cos_c2*sin_phi*sin_b1 + (-cos_phi*cos_b2*sin_b1 + cos_b1*sin_b2)*sin_c2 ! Czx
    cc_v2(3) = -cos_phi*cos_b2*cos_c2*sin_b1 + cos_b1*cos_c2*sin_b2 - sin_phi*sin_b1*sin_c2 ! Czy

    cc_v2(4) =  cos_b2*sin_b1*sin_c1 - sin_b2*(cos_c1*sin_phi + cos_phi*cos_b1*sin_c1)      ! Cxz
    cc_v2(5) = -cos_b1*cos_c2*sin_phi*sin_c1 + (cos_b2*cos_c1*sin_phi + sin_b1*sin_b2*sin_c1)*sin_c2 &
               + cos_phi*(cos_c1*cos_c2 + cos_b1*cos_b2*sin_c1*sin_c2)                        ! Cxx
    cc_v2(6) =  cos_c2*sin_b1*sin_b2*sin_c1 + cos_b2*cos_c2*(cos_c1*sin_phi + cos_phi*cos_b1*sin_c1) &
               + (-cos_phi*cos_c1 + cos_b1*sin_phi*sin_c1)*sin_c2                             ! Cxy

    cc_v2(7) =  cos_b2*cos_c1*sin_b1 + sin_b2*(-cos_phi*cos_b1*cos_c1 + sin_phi*sin_c1)      ! Cyz
    cc_v2(8) =  cos_c1*sin_b1*sin_b2*sin_c2 + cos_b1*cos_c1*(-cos_c2*sin_phi + cos_phi*cos_b2*sin_c2) &
               - sin_c1*(cos_phi*cos_c2 + cos_b2*sin_phi*sin_c2)                              ! Cyx
    cc_v2(9) = -cos_b2*cos_c2*sin_phi*sin_c1 + cos_c1*(cos_c2*sin_b1*sin_b2 + cos_b1*sin_phi*sin_c2) &
               + cos_phi*(cos_b1*cos_b2*cos_c1*cos_c2 + sin_c1*sin_c2)                        ! Cyy
  end subroutine generate_coordinates_v2

  !===================== MULTIPOLE ============================
  function multipole_sph_v3(ind) result(E)
    use Fitting_Constant_v2, only: get_coeff_Fit
    integer(I4), intent(in) :: ind
    real(RK) :: E, rinv, rinv_pow, scale
    integer(I4) :: order

    E = 0.0_RK
    rinv = 1.0_RK / cal_coord_v2(1)
    rinv_pow = rinv

    do order = 1, MAX_L
      if (get_coeff_Fit(ind, order, 'M') > 0) then
        scale = (C3*C1) * (C2**order)
        E = E + scale * multipole_order(ind, order) * rinv_pow
      end if
      rinv_pow = rinv_pow * rinv
    end do
  end function multipole_sph_v3

  !------------------------------------------------------------
  function multipole_order(ind, order) result(val)
    use Fitting_Constant_v2, only: get_coeff_multipole
    integer(I4), intent(in) :: ind, order
    real(RK) :: val
    integer(I4) :: i, j, ci, cj
    real(RK) :: Qai, Qbj
    real(RK), parameter :: ZERO = 0.0_RK
    real(RK) :: A_Mult(225), B_Mult(225)

    val = 0.0_RK
    A_Mult = get_coeff_multipole(ind, 'A')
    B_Mult = get_coeff_multipole(ind, 'B')

    do i = 0, order-1
      j = order - 1 - i
      do ci = 0, 2*i
        Qai = A_Mult(i*i + 1 + ci)
        if (dabs(Qai) > EPS) then
          do cj = 0, 2*j
            Qbj = B_Mult(j*j + 1 + cj)
            if (dabs(Qbj) > EPS) then
              val = val + Qai * Qbj * t_tensor_v2(i+1, ci+1, j+1, cj+1)
            end if
          end do
        end if
      end do
    end do
  end function multipole_order

  !===================== INDUCTION ============================
  function induction_sph_v3(ind) result(E)
    use Fitting_Constant_v2, only: get_coeff_Fit
    integer(I4), intent(in) :: ind
    real(RK) :: E
    integer(I4) :: order

    E = 0.0_RK
    do order = 1, MAX_L
      if (get_coeff_Fit(ind, order, 'I') > 0) then
        E = E + induction_order(order, ind, 1_I4) + induction_order(order, ind, 0_I4)
      end if
    end do
  end function induction_sph_v3

  !------------------------------------------------------------
  function induction_order(order, ind, index) result(val)
    integer(I4), intent(in) :: order, ind, index
    real(RK) :: val
    integer(I4) :: l1, l2, i, j
    real(RK) :: res, R

    R = cal_coord_v2(1)
    res = 0.0_RK

    do l1 = 1, order-3
      do l2 = 1, order-3
        if (l1 + l2 + 2 <= order) then
          do i = 0, order-2-l1-l2
            do j = 0, order-2-l1-l2
              if (i + j + l1 + l2 + 2 == order) then
                res = res + induction_ij_l1l2(i, j, l1, l2, ind, index)
              end if
            end do
          end do
        end if
      end do
    end do

    val = (-0.5_RK * (C3*C1) * (C2**order) * res) / (R**order)
  end function induction_order

  !------------------------------------------------------------
  function induction_ij_l1l2(i, j, l1, l2, ind, index) result(val)
    use Fitting_Constant_v2, only: get_coeff_multipole_by_index, get_coeff_polarizability_by_index
    integer(I4), intent(in) :: i, j, l1, l2, ind, index
    real(RK) :: val

    real(RK) :: Qai, Qbj, comp_a_k1_k2
    integer(I4) :: ci, cj, k1, k2, cpn
    integer(I4) :: ni, nj, nl1, nl2, lmin, lmax
    real(RK), allocatable :: Qa_cpn(:), Qb_cpn(:), pol_arr(:)

    val = 0.0_RK

    ni = 2*i + 1; nj = 2*j + 1; nl1 = 2*l1 + 1; nl2 = 2*l2 + 1
    lmin = min(l1, l2); lmax = max(l1, l2)

    allocate(Qa_cpn(ni), Qb_cpn(nj), pol_arr(nl1*nl2))

    if (index == 1) then
      Qa_cpn = get_coeff_multipole_by_index(ind, 'A', i)
      Qb_cpn = get_coeff_multipole_by_index(ind, 'A', j)
      pol_arr = get_coeff_polarizability_by_index(ind, 'B', lmin, lmax, nl1*nl2)
    else
      Qa_cpn = get_coeff_multipole_by_index(ind, 'B', i)
      Qb_cpn = get_coeff_multipole_by_index(ind, 'B', j)
      pol_arr = get_coeff_polarizability_by_index(ind, 'A', lmin, lmax, nl1*nl2)
    end if

    do ci = 1, ni
      Qai = Qa_cpn(ci)
      if (dabs(Qai) > EPS) then
        do cj = 1, nj
          Qbj = Qb_cpn(cj)
          if (dabs(Qbj) > EPS) then
            do k1 = 1, nl1
              do k2 = 1, nl2
                cpn = get_induction_cpn(l1, l2, k1, k2)
                comp_a_k1_k2 = pol_arr(cpn)
                if (dabs(comp_a_k1_k2) > EPS) then
                  if (index == 0) then
                    val = val + Qai * Qbj * comp_a_k1_k2 * &
                                ( t_tensor_v2(l1+1, k1, i+1, ci) * t_tensor_v2(l2+1, k2, j+1, cj) )
                  else
                    val = val + Qai * Qbj * comp_a_k1_k2 * &
                                ( t_tensor_v2(i+1, ci, l1+1, k1) * t_tensor_v2(j+1, cj, l2+1, k2) )
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

  !------------------------------------------------------------
  pure function get_induction_cpn(l1, l2, li, lj) result(idx)
    integer(I4), intent(in) :: l1, l2, li, lj
    integer(I4) :: idx

    if (l1 > l2) then
      idx = (lj-1) * (2*l1+1) + li
    else
      idx = (li-1) * (2*l2+1) + lj
    end if
  end function get_induction_cpn

  !===================== DISPERSION ===========================
  function dispersion_sph_v3(ind) result(E)
    use Fitting_Constant_v2, only: get_coeff_Fit
    integer(I4), intent(in) :: ind
    real(RK) :: E
    integer(I4) :: order

    E = 0.0_RK
    do order = 1, MAX_L
      if (get_coeff_Fit(ind, order, 'D') > 0) then
        E = E + dispersion_order(ind, order)
      end if
    end do
  end function dispersion_sph_v3

  !------------------------------------------------------------
  function dispersion_order(ind, order) result(val)
    integer(I4), intent(in) :: ind, order
    real(RK) :: val
    integer(I4) :: l1, l2, t1, t2
    real(RK) :: res, R

    res = 0.0_RK
    R   = cal_coord_v2(1)

    do l1 = 1, order-2
      do l2 = 1, order-2-l1
        do t1 = 1, order-2-l1-l2
          do t2 = 1, order-2-l1-l2-t1
            if (l1 + l2 + t1 + t2 + 2 == order) then
              res = res + dispersion_l1l2_t1t2(ind, l1, l2, t1, t2)
            end if
          end do
        end do
      end do
    end do

    val = - ( (C3*C1) * (C2**order) * res ) / (R**order)
  end function dispersion_order

  !------------------------------------------------------------
  function dispersion_l1l2_t1t2(ind, l1, l2, t1, t2) result(val)
    use Fitting_Constant_v2, only: get_coeff_dispersion_by_index
    integer(I4), intent(in) :: ind, l1, l2, t1, t2
    real(RK) :: val

    integer(I4) :: li, lj, ti, tj, cpn
    real(RK) :: disp_coeff
    real(RK) :: disp_arr( (2*l1+1)*(2*l2+1)*(2*t1+1)*(2*t2+1) )

    val = 0.0_RK
    disp_arr = get_coeff_dispersion_by_index(ind, l1, l2, t1, t2)

    do li = 0, 2*l1
      do lj = 0, 2*l2
        do ti = 0, 2*t1
          do tj = 0, 2*t2
            cpn = get_dispersion_cpn(l1, l2, t1, t2, li, lj, ti, tj)
            disp_coeff = disp_arr(cpn)
            if (dabs(disp_coeff) > EPS) then
              val = val + disp_coeff * t_tensor_v2(l1+1, li+1, t1+1, ti+1) * &
                                  t_tensor_v2(l2+1, lj+1, t2+1, tj+1)
            end if
          end do
        end do
      end do
    end do
  end function dispersion_l1l2_t1t2

  !------------------------------------------------------------
  pure function get_dispersion_cpn(l1,l2,t1,t2,li,lj,ti,tj) result(idx)
    integer(I4), intent(in) :: l1,l2,t1,t2,li,lj,ti,tj
    integer(I4) :: idx

    if (l1>l2 .and. t1>t2) then
      idx =  lj * (2*l1+1) * (2*t2+1) * (2*t1+1) + li * (2*t2+1) * (2*t1+1) + tj * (2*t1+1) + ti + 1
    else if (l1>l2 .and. t1<=t2) then
      idx =  lj * (2*l1+1) * (2*t1+1) * (2*t2+1) + li * (2*t1+1) * (2*t2+1) + ti * (2*t2+1) + tj + 1
    else if (l1<=l2 .and. t1>t2) then
      idx =  li * (2*l2+1) * (2*t2+1) * (2*t1+1) + lj * (2*t2+1) * (2*t1+1) + tj * (2*t1+1) + ti + 1
    else
      idx =  li * (2*l2+1) * (2*t1+1) * (2*t2+1) + lj * (2*t1+1) * (2*t2+1) + ti * (2*t2+1) + tj + 1
    end if
  end function get_dispersion_cpn

  !====================== TOTAL ===============================
  function get_total_interaction_energy(coeff_index, general_coordinates_ZXZ) result(E)
    use Fitting_Constant_v2, only: get_coeff_max_t_tensor_order
    integer(I4), intent(in) :: coeff_index
    real(RK)  , intent(in) :: general_coordinates_ZXZ(6)
    real(RK) :: E
    integer(I4) :: max_t_tensor_order

    max_t_tensor_order = get_coeff_max_t_tensor_order(coeff_index)
    call tensors_initialization_v2(max_t_tensor_order, general_coordinates_ZXZ)

    E = multipole_sph_v3(coeff_index) + induction_sph_v3(coeff_index) + dispersion_sph_v3(coeff_index)
  end function get_total_interaction_energy

end module Geometry_Constant_v2
