module multipole_impl_omp
  use omp_lib
  implicit none
  integer, parameter :: rk = kind(1.0d0)
  integer, parameter :: maxl = 15
  integer, parameter :: off(0:maxl-1) = (/ (i*i + 1, i=0,maxl-1) /)

  ! If t_tensor_v2 is a pure/external function, give it an explicit interface
  ! and allow the compiler to vectorize calls that appear in simd loops:
  interface
     pure real(rk) function t_tensor_v2(ii,ci,jj,cj)
       integer, intent(in) :: ii, ci, jj, cj
     end function t_tensor_v2
  end interface
  !$omp declare simd(t_tensor_v2) notinbranch

contains

  pure function multipole_order_fast(order, A, B) result(val)
    integer,  intent(in) :: order
    real(rk), intent(in) :: A(:), B(:)
    real(rk) :: val
    integer :: i, j, ci, cj, ia, jb
    real(rk) :: Qai, Qbj, tmp

    val = 0.0_rk
    do i = 0, order-1
       j  = order - 1 - i
       ia = off(i)   ! base index for block i
       jb = off(j)   ! base index for block j

       ! Option B: thread the ci-loop if we are NOT already in a parallel region
       !$omp parallel do default(none) if (.not. omp_in_parallel()) &
       !$omp& shared(i,j,ia,jb,A,B) private(ci,cj,Qai,Qbj,tmp) reduction(+:val)
       do ci = 0, 2*i
          Qai = A(ia + ci)
          if (Qai == 0.0_rk) cycle

          tmp = 0.0_rk
          !$omp simd reduction(+:tmp)
          do cj = 0, 2*j
             Qbj = B(jb + cj)
             if (Qbj /= 0.0_rk) then
                tmp = tmp + Qai * Qbj * t_tensor_v2(i+1,ci+1,j+1,cj+1)
             end if
          end do
          val = val + tmp
       end do
       !$omp end parallel do
    end do
  end function multipole_order_fast


  pure function multipole_sph_v3_fast(R, C1, C2, C3, A_mult, B_mult, Mmask) result(E)
    real(rk), intent(in) :: R, C1, C2, C3
    real(rk), intent(in) :: A_mult(:), B_mult(:)     ! length 225 each
    logical, intent(in)  :: Mmask(maxl)              ! active orders
    real(rk) :: E
    integer  :: order
    real(rk) :: base, ratio, powo

    E     = 0.0_rk
    base  = C3*C1
    ratio = C2 / R
    powo  = ratio

    do order = 1, maxl
       if (Mmask(order)) then
          E = E + base * powo * multipole_order_fast(order, A_mult, B_mult)
       end if
       powo = powo * ratio
    end do
  end function multipole_sph_v3_fast

end module multipole_impl_omp