!FittingConstant is the module in charge of all constants related with the fit 
!It can handle several coefficient files at the same time

module Fitting_Constant_v2
    save

    type fit_contant
        character(:), allocatable :: filename
        real (kind=8) :: Zero
        integer (kind=4) ::initflag
        integer (kind=4) ::max_t_tensor_order

        integer (kind=4), dimension (15) :: M_Fit
        integer (kind=4), dimension (15) :: D_Fit
        integer (kind=4), dimension (15) :: I_Fit

        !Multipoles!

        real (kind=8) , dimension(225)   :: A_Mult,B_Mult

        !Polarizability!

        real (kind=8) , dimension(6,12,195)  :: A_Pol,B_Pol

        !Dispersion!

        real (kind=8) , dimension(5,10,5,10,3087)   :: Disp

    contains
        procedure, pass :: initializer
        procedure, pass :: read_parameters

    end type

    integer (kind=4),parameter :: NARRAY = 5
    type(fit_contant) :: coeff(NARRAY)

    private :: coeff
    public::find_coeff_set,&
            get_coeff_index,&
            get_coeff_zero,&
            get_coeff_max_t_tensor_order,&
            get_coeff_Fit,&
            get_coeff_multipole,&
            get_coeff_multipole_by_index,&
            get_coeff_polarizability_by_index,&
            get_coeff_dispersion_by_index

    contains

        subroutine initializer(this,filename)
            implicit none
            class(fit_contant), intent(out) :: this
            character(len=*), intent(in) ::filename

            this%initflag = 1
            this%filename = filename

        end subroutine initializer

        subroutine read_parameters(this)
            implicit none
            class(fit_contant), intent(inout) :: this

            character(len = 200) :: row
            integer (kind=4)::iord,mord,dord
            integer (kind=4)::i,j,l1,l2,t1,t2,ln
            real (kind=8) :: polarr(195)

            if (this%initflag==1)then

                this%initflag = 2

                open( 10, file = this%filename )

                read( 10, *) row
                read( 10, *) row
                read( 10, *) row
                read( 10, *) row
                read( 10, *) row
                read( 10, *) row


                read(10, *)  row,this%Zero
                read( 10, *) row

                read(10, *)  row,this%M_Fit
                read(10, *)  row,this%I_Fit
                read(10, *)  row,this%D_Fit

                iord = maxval(this%I_Fit);
                mord = maxval(this%M_Fit);
                mord = maxval([mord,iord-3]);
                dord = maxval(this%D_Fit);

                this%max_t_tensor_order = maxval([iord-2,mord,dord-3])


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

        end subroutine read_parameters

        subroutine find_coeff_set(filename,ind)
            implicit none
            character(*), intent(in) :: filename
            integer (kind=4), intent(out) :: ind
            integer (kind=4):: i

            ind=-1
            if (NARRAY<1)then
                Return
            else
                do i=1,NARRAY
                    if (coeff(i)%filename == filename)then
                        ind = i
                        return
                    end if
                end do

            end if

        end subroutine find_coeff_set

        function last_coeff_set()
            implicit none
            integer (kind=4) :: last_coeff_set
            integer (kind=4):: i

            last_coeff_set = 0
            if (NARRAY<1)then
                Return
            else

                do i=1,NARRAY
                    if (LEN(coeff(i)%filename)<1)then
                        last_coeff_set = i-1
                        return
                    end if
                end do

                if (last_coeff_set == NARRAY )then
                    write(*,*)"The maximun number of coefficients sets is :",NARRAY,&
                            "to change the maximun go to module module Fit and change NARRAY"
                    last_coeff_set=-10
                end if
            end if


        end function last_coeff_set

        function get_coeff_index(filename)
            implicit none
            character(*), intent(in) :: filename
            integer (kind=4) :: get_coeff_index
            integer (kind=4)::ind,last_index

            call find_coeff_set(filename,ind)


            if (ind < 1)then

                get_coeff_index = last_coeff_set() + 1

                call coeff(get_coeff_index)%initializer(filename)
                call coeff(get_coeff_index)%read_parameters()

            else
                get_coeff_index = ind
                return
            end if

        end function get_coeff_index

        function get_coeff_zero(coeff_index)
            implicit none
            integer (kind=4),intent(in) :: coeff_index
            real (kind=8)::get_coeff_zero

            get_coeff_zero = coeff(coeff_index)%Zero

        end function get_coeff_zero

        function get_coeff_max_t_tensor_order(coeff_index)
            implicit none
            integer (kind=4),intent(in) :: coeff_index
            integer (kind=4)::get_coeff_max_t_tensor_order

            get_coeff_max_t_tensor_order = coeff(coeff_index)%max_t_tensor_order

        end function get_coeff_max_t_tensor_order

        function get_coeff_Fit(coeff_index,order,interaction)
            implicit none
            integer (kind=4),intent(in) :: coeff_index,order
            character(len=1),intent(in)::interaction
            integer (kind=4)::get_coeff_Fit
            if (interaction=="M") then
                get_coeff_Fit = coeff(coeff_index)%M_Fit(order)
            elseif(interaction=="I")  then
                get_coeff_Fit = coeff(coeff_index)%I_Fit(order)
            elseif(interaction=="D") then
                get_coeff_Fit = coeff(coeff_index)%D_Fit(order)
            end if

        end function get_coeff_Fit

        function get_coeff_multipole(coeff_index,molecule_label)
            implicit none
            integer (kind=4),intent(in) :: coeff_index
            character(len=1),intent(in)::molecule_label
            real (kind=8),dimension(225)::get_coeff_multipole

            if (molecule_label == "A")then
                get_coeff_multipole = coeff(coeff_index)%A_Mult
                else
                get_coeff_multipole = coeff(coeff_index)%B_Mult
            end if


        end function get_coeff_multipole

        function get_coeff_multipole_by_index(coeff_index,molecule_label,i)
            implicit none
            integer (kind=4),intent(in) :: coeff_index,i
            character(len=1),intent(in)::molecule_label
            real (kind=8),dimension(225)::get_coeff_multipole_by_index

            if (molecule_label == "A")then
                get_coeff_multipole_by_index = coeff(coeff_index)%A_Mult(i**2 + 1:(i+1)**2)
            else
                get_coeff_multipole_by_index = coeff(coeff_index)%B_Mult(i**2 + 1:(i+1)**2)
            end if


        end function get_coeff_multipole_by_index

        function get_coeff_polarizability_by_index(coeff_index,molecule_label,lmin,lmax,nl1l2)
            implicit none
            integer (kind=4),intent(in) :: coeff_index,lmin,lmax,nl1l2
            character(len=1),intent(in) :: molecule_label
            real (kind=8) , dimension(nl1l2) :: get_coeff_polarizability_by_index

            if (molecule_label == "A")then
                get_coeff_polarizability_by_index = coeff(coeff_index)%A_Pol(lmin,lmax,1:nl1l2)
            else
                get_coeff_polarizability_by_index = coeff(coeff_index)%B_Pol(lmin,lmax,1:nl1l2)
            end if


        end function get_coeff_polarizability_by_index

        function get_coeff_dispersion_by_index(coeff_index,l1,l2,t1,t2)
            implicit none
            integer (kind=4),intent(in) :: coeff_index,l1,l2,t1,t2
            real (kind=8) , dimension((2*l1+1)*(2*l2+1)*(2*t1+1)*(2*t2+1)) :: get_coeff_dispersion_by_index
            integer (kind=4) :: lmin,lmax,tmin,tmax
            lmin = minval([l1,l2])
            lmax = maxval([l1,l2])
            tmin = minval([t1,t2])
            tmax = maxval([t1,t2])

            get_coeff_dispersion_by_index = coeff(coeff_index)%Disp(lmin,lmax,tmin,tmax, &
                                                                    1:(2*l1+1)*(2*l2+1)*(2*t1+1)*(2*t2+1))


        end function get_coeff_dispersion_by_index

end module Fitting_Constant_v2