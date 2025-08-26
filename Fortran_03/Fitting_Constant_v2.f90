!********************************************************
! FittingConstant is the module in charge of all constants related with the fit
! It can handle several coefficient files at the same time
!
! Also provides a global index map (t_index_map) and an indexing function (t_index)
! to flatten t-tensor components (L,2L+1,L,2L+1) -> [1..NNZ] in a stable, precomputed way.
!********************************************************
module Fitting_Constant_v2
    use iso_fortran_env, only : real64, int32
    implicit none
    save

    !========================================================
    ! Types and data for fitted constants
    !========================================================
    type fit_contant
        character(:), allocatable :: filename
        real(real64) :: Zero
        integer(int32) :: initflag
        integer(int32) :: max_t_tensor_order

        integer(int32), dimension(15) :: M_Fit
        integer(int32), dimension(15) :: D_Fit
        integer(int32), dimension(15) :: I_Fit

        ! Multipoles
        real(real64), dimension(225) :: A_Mult, B_Mult

        ! Polarizability
        ! A_Pol(lmin,lmax, (2*lmin+1)*(2*lmax+1)) stored in first ln entries
        real(real64), dimension(6,12,195) :: A_Pol, B_Pol

        ! Dispersion
        ! Disp(lmin,lmax,tmin,tmax, (2*l1+1)*(2*l2+1)*(2*t1+1)*(2*t2+1)) in first ln entries
        real(real64), dimension(5,10,5,10,3087) :: Disp

    contains
        procedure, pass :: initializer
        procedure, pass :: read_parameters
    end type

    integer(int32), parameter :: NARRAY = 5
    type(fit_contant) :: coeff(NARRAY)

    !========================================================
    ! Global t-tensor index map (1-based like the original 4D array)
    ! idx = t_index_map(la+1,ka+1,lb+1,kb+1)
    !========================================================
    integer(int32), parameter :: TMAP_L_CONST = 15_int32

    integer(int32) :: tmap_L     = 0         ! maxlevel used to build map
    integer(int32) :: tmap_NNZ   = 0         ! number of valid components
    logical :: tmap_ready = .false.

    integer(int32), allocatable :: t_index_map(:,:,:,:)  ! (L,2L+1,L,2L+1), 1-based

    !=======================  Precomputed factorial_nn  =======================
    ! Precomputed sqrt-factorial combination map:
    !   fact_nn(la,ka1,lb,kb1) = sqrt( ( (la+ka1)!/(la-ka1)! ) * ( (lb+kb1)!/(lb-kb1)! ) )
    ! Indices are 0-based and valid only when ka1<=la and kb1<=lb; otherwise 0.
    !========================================================
    integer(int32), save :: fnn_L = 0
    logical,        save :: fnn_ready = .false.
    real(real64),   allocatable, save :: fac_tbl(:)                 ! 0..2*L
    real(real64),   allocatable, save :: fact_nn_map(:,:,:,:)       ! (0:L,0:L,0:L,0:L)

    private :: coeff
    public :: find_coeff_set, &
              get_coeff_index, &
              get_coeff_zero, &
              get_coeff_max_t_tensor_order, &
              get_coeff_Fit, &
              get_coeff_multipole, &
              get_coeff_multipole_by_index, &
              get_coeff_polarizability_by_index, &
              get_coeff_dispersion_by_index, &
              ! expose the index infra
              t_index, t_index_map, tmap_L, tmap_NNZ, ensure_t_index_map_ready
    public :: fact_nn, ensure_fact_nn_ready
              
contains

    !========================================================
    ! t-index map helpers
    !========================================================

    subroutine ensure_fact_nn_ready(Lin)
        use iso_fortran_env, only: real64, int32
        implicit none
        integer(int32), intent(in), optional :: Lin
        integer(int32) :: L, la, lb, ka1, kb1, m
        real(real64), allocatable :: fact_arr(:)  ! factorials 0..2L

        L = merge(Lin, TMAP_L_CONST, present(Lin))  ! default to your max (e.g. 15)

        if (fnn_ready .and. fnn_L == L) return

        if (allocated(fact_nn_map)) deallocate(fact_nn_map)
        allocate(fact_nn_map(L+1, L+1, L+1, L+1))  ! store at (la+1,ka1+1,lb+1,kb1+1)

        ! precompute factorials 0..2L in double
        allocate(fact_arr(0:2*L))
        fact_arr(0) = 1.0_real64
        do m = 1, 2*L
            fact_arr(m) = fact_arr(m-1) * real(m, real64)
        end do

        do la = 0, L
            do lb = 0, L
            do ka1 = 0, L
                do kb1 = 0, L
                if (ka1 <= la .and. kb1 <= lb) then
                    fact_nn_map(la+1,ka1+1,lb+1,kb1+1) = &
                    dsqrt( (fact_arr(la+ka1)/fact_arr(la-ka1)) * &
                            (fact_arr(lb+kb1)/fact_arr(lb-kb1)) )
                else
                    fact_nn_map(la+1,ka1+1,lb+1,kb1+1) = 0.0_real64
                end if
                end do
            end do
            end do
        end do

        deallocate(fact_arr)
        fnn_L = L
        fnn_ready = .true.
    end subroutine ensure_fact_nn_ready

    ! Elemental accessor so callers can write:  fn = fact_nn(la+1,ka1+1,lb+1,kb1+1)
    pure real(real64) function fact_nn(la1,ka11,lb1,kb11) result(fn)
        use iso_fortran_env, only: real64, int32
        implicit none
        integer(int32), intent(in) :: la1, ka11, lb1, kb11   ! 1-based indices
        ! ASSUMES ensure_fact_nn_ready() has been called already.
        fn = fact_nn_map(la1, ka11, lb1, kb11)
    end function fact_nn


    subroutine ensure_t_index_map_ready(Lin)
        integer(int32), intent(in), optional :: Lin
        integer(int32) :: Luse
        integer(int32) :: order, la, lb, ka_, kb_, count

        Luse = merge(Lin, TMAP_L_CONST, present(Lin))

        if (tmap_ready .and. tmap_L == Luse) return

        if (allocated(t_index_map)) deallocate(t_index_map)
        allocate(t_index_map(Luse, 2*Luse+1, Luse, 2*Luse+1))
        t_index_map = -1

        count = 1
        do order = 1, Luse
            do la = 0, order-1
                lb = order - la - 1
                do ka_ = 0, 2*la
                    do kb_ = 0, 2*lb
                        t_index_map(la+1, ka_+1, lb+1, kb_+1) = count
                        count = count + 1
                    end do
                end do
            end do
        end do

        tmap_L     = Luse
        tmap_NNZ   = count - 1
        tmap_ready = .true.
    end subroutine ensure_t_index_map_ready

    pure integer(int32) function t_index(i1,i2,i3,i4) result(ind)
        integer(int32), intent(in) :: i1,i2,i3,i4
        ind = t_index_map(i1,i2,i3,i4)
    end function t_index

    !========================================================
    ! Coefficients IO
    !========================================================
    subroutine initializer(this, filename)
        class(fit_contant), intent(out) :: this
        character(len=*), intent(in) :: filename

        this%initflag = 1
        this%filename = filename
    end subroutine initializer

    subroutine read_parameters(this)
        class(fit_contant), intent(inout) :: this

        character(len=200) :: row
        integer(int32) :: iord, mord, dord
        integer(int32) :: i, j, l1, l2, t1, t2, ln
        integer(int32) :: u, ios

        if (this%initflag == 1) then
            this%initflag = 2

            open(newunit=u, file=this%filename, status='old', action='read', iostat=ios)
            if (ios /= 0) then
                write(*,*) 'Error opening coefficients file: ', trim(this%filename), ' iostat=', ios
                return
            end if

            ! Skip header lines as in the original reader
            read(u, *) row
            read(u, *) row
            read(u, *) row
            read(u, *) row
            read(u, *) row
            read(u, *) row

            read(u, *) row, this%Zero
            read(u, *) row

            read(u, *) row, this%M_Fit
            read(u, *) row, this%I_Fit
            read(u, *) row, this%D_Fit

            iord = maxval(this%I_Fit)
            mord = maxval(this%M_Fit)
            mord = maxval([mord, iord - 3])
            dord = maxval(this%D_Fit)

            this%max_t_tensor_order = maxval([iord - 2, mord, dord - 3])

            if (mord > 0) then
                read(u, *) row, this%A_Mult(1:mord**2)  ! A_Mult
                read(u, *) row, this%B_Mult(1:mord**2)  ! B_Mult
            end if

            if (iord >= 4) then
                do i = 1, iord - 3
                    do j = i, iord - 3
                        if (i + j <= iord - 2) then
                            ln = (2*i + 1)*(2*j + 1)
                            read(u, *) row, this%A_Pol(i, j, 1:ln)  ! PA_i-j
                            read(u, *) row, this%B_Pol(i, j, 1:ln)  ! PB_i-j
                        end if
                    end do
                end do
            end if

            if (dord >= 6) then
                do l1 = 1, dord - 5
                    do l2 = l1, dord - 5
                        do t1 = 1, dord - 5
                            do t2 = t1, dord - 5
                                if (l1 + l2 + t1 + t2 <= dord - 2) then
                                    ln = (2*l1 + 1)*(2*l2 + 1)*(2*t1 + 1)*(2*t2 + 1)
                                    read(u, *) row, this%Disp(l1, l2, t1, t2, 1:ln)
                                end if
                            end do
                        end do
                    end do
                end do
            end if

            close(u)
        end if
    end subroutine read_parameters

    !========================================================
    ! Public API for coefficients
    !========================================================
    subroutine find_coeff_set(filename, ind)
        character(*), intent(in) :: filename
        integer(int32), intent(out) :: ind
        integer(int32) :: i

        ind = -1
        if (NARRAY < 1) then
            return
        else
            do i = 1, NARRAY
                if (allocated(coeff(i)%filename)) then
                    if (coeff(i)%filename == filename) then
                        ind = i
                        return
                    end if
                end if
            end do
        end if
    end subroutine find_coeff_set

    function last_coeff_set() result(idx)
        integer(int32) :: idx
        integer(int32) :: i

        idx = 0
        if (NARRAY < 1) then
            return
        else
            do i = 1, NARRAY
                if (.not. allocated(coeff(i)%filename)) then
                    idx = i - 1
                    return
                end if
                if (len(coeff(i)%filename) < 1) then
                    idx = i - 1
                    return
                end if
            end do

            if (idx == NARRAY) then
                write(*,*) "The maximun number of coefficients sets is :", NARRAY, &
                           " to change the maximun go to module Fit and change NARRAY"
                idx = -10
            end if
        end if
    end function last_coeff_set

    function get_coeff_index(filename) result(idx)
        character(*), intent(in) :: filename
        integer(int32) :: idx
        integer(int32) :: ind

        call find_coeff_set(filename, ind)

        if (ind < 1) then
            idx = last_coeff_set() + 1
            if (idx < 1 .or. idx > NARRAY) then
                write(*,*) "No available slot in coeff array."
                idx = -1
                return
            end if
            call coeff(idx)%initializer(filename)
            call coeff(idx)%read_parameters()
        else
            idx = ind
            return
        end if
    end function get_coeff_index

    function get_coeff_zero(coeff_index) result(val)
        integer(int32), intent(in) :: coeff_index
        real(real64) :: val
        val = coeff(coeff_index)%Zero
    end function get_coeff_zero

    function get_coeff_max_t_tensor_order(coeff_index) result(val)
        integer(int32), intent(in) :: coeff_index
        integer(int32) :: val
        val = coeff(coeff_index)%max_t_tensor_order
    end function get_coeff_max_t_tensor_order

    function get_coeff_Fit(coeff_index, order, interaction) result(val)
        integer(int32), intent(in) :: coeff_index, order
        character(len=1), intent(in) :: interaction
        integer(int32) :: val

        select case (interaction)
        case ("M")
            val = coeff(coeff_index)%M_Fit(order)
        case ("I")
            val = coeff(coeff_index)%I_Fit(order)
        case ("D")
            val = coeff(coeff_index)%D_Fit(order)
        case default
            val = 0_int32
        end select
    end function get_coeff_Fit

    function get_coeff_multipole(coeff_index, molecule_label) result(arr)
        integer(int32), intent(in) :: coeff_index
        character(len=1), intent(in) :: molecule_label
        real(real64), dimension(225) :: arr

        if (molecule_label == "A") then
            arr = coeff(coeff_index)%A_Mult
        else
            arr = coeff(coeff_index)%B_Mult
        end if
    end function get_coeff_multipole

    function get_coeff_multipole_by_index(coeff_index, molecule_label, i) result(arr)
        integer(int32), intent(in) :: coeff_index, i
        character(len=1), intent(in) :: molecule_label
        real(real64), dimension(225) :: arr

        if (molecule_label == "A") then
            arr = coeff(coeff_index)%A_Mult(i**2 + 1:(i+1)**2)
        else
            arr = coeff(coeff_index)%B_Mult(i**2 + 1:(i+1)**2)
        end if
    end function get_coeff_multipole_by_index

    function get_coeff_polarizability_by_index(coeff_index, molecule_label, lmin, lmax, nl1l2) result(arr)
        integer(int32), intent(in) :: coeff_index, lmin, lmax, nl1l2
        character(len=1), intent(in) :: molecule_label
        real(real64), dimension(nl1l2) :: arr

        if (molecule_label == "A") then
            arr = coeff(coeff_index)%A_Pol(lmin, lmax, 1:nl1l2)
        else
            arr = coeff(coeff_index)%B_Pol(lmin, lmax, 1:nl1l2)
        end if
    end function get_coeff_polarizability_by_index

    function get_coeff_dispersion_by_index(coeff_index, l1, l2, t1, t2) result(arr)
        integer(int32), intent(in) :: coeff_index, l1, l2, t1, t2
        real(real64), dimension((2*l1+1)*(2*l2+1)*(2*t1+1)*(2*t2+1)) :: arr
        integer(int32) :: lmin, lmax, tmin, tmax

        lmin = minval([l1, l2])
        lmax = maxval([l1, l2])
        tmin = minval([t1, t2])
        tmax = maxval([t1, t2])

        arr = coeff(coeff_index)%Disp(lmin, lmax, tmin, tmax, &
              1:(2*l1+1)*(2*l2+1)*(2*t1+1)*(2*t2+1))
    end function get_coeff_dispersion_by_index

end module Fitting_Constant_v2
