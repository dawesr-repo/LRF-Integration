

program testing_subroutine

        use Testing_v2, only: file_checking,&
                              running_time_performance,&
                              check_energy_MATLAB

        implicit none
        integer (kind=4):: date_time (8),fileout_number=12
        character (len = 10) big_ben (3)
        integer (kind=4):: n_sys = 49
        integer (kind=4):: xdim_arr(50)
        type :: var_str
                character(len=:), allocatable :: label
        end type var_str

        type(var_str) :: sys(50)
        type(var_str) :: symmetries(15)
        integer (kind=4):: i,k,h,count
        integer (kind=4),dimension(7):: xdim = [3,3,3,3,2,2,0]

        symmetries(1)%label  = "C1(1)"
        symmetries(2)%label  = "T(2)"
        symmetries(3)%label  = "Ih(2)"
        symmetries(4)%label  = "Cs(1)"
        symmetries(5)%label  = "D_inf_h(1)"
        symmetries(6)%label  = "C_inf_v(1)"
        symmetries(7)%label  = "Spherical(1)"

        count = 1
        do k=1,7
                do h=1,7
                        if (xdim(k) >= xdim(h)) then
                                sys(count)%label = symmetries(k)%label//"_"//symmetries(h)%label
                                xdim_arr(count) = xdim(k)+ xdim(h)
                                if (xdim(k)==1 .and. xdim(h)==1) then
                                        xdim_arr(count) = 1
                                end if
                                count = count + 1
                        end if
                end do
        end do

        call date_and_time (big_ben (1), big_ben (2), big_ben (3), date_time)

        call file_checking('../testing_datafiles/output.test.txt',fileout_number)
        rewind(fileout_number)

        write(fileout_number,*)"******************************************************************************"
        write(fileout_number,*)  "Test Day and Time Record"
        write(fileout_number,*)  "Month / Day / Year: ",date_time(2),"/",date_time(3),"/",date_time(1)
        write(fileout_number,*)  "Hr    / Min / Sec : ",date_time(5),":",date_time(6),":",date_time(7)

        call running_time_performance('../testing_datafiles/coefficients/C1(1)_C1(1)_Coeff.txt',fileout_number)

        do i=1,5!count-1
                call check_energy_MATLAB( sys(i)%label,&
                                          xdim_arr(i),&
                                          0,&
                                          fileout_number&
                                        )
        end do

        close(fileout_number)
end program testing_subroutine