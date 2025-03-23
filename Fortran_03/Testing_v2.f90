! Created by albpl on 3/11/2025.

module Testing_v2
contains

    subroutine file_checking(file_name,num)

        implicit none
        character(*),intent(in)::file_name
        Integer,intent(in)::num
        logical :: exist

        inquire(file=file_name, exist=exist)
        if (exist) then
            open(num, file=file_name, status="old", position="append", action="write")
        else
            open(num, file=file_name, status="new", action="write")
        end if

    END subroutine file_checking

    ! TEST Functions
    subroutine running_time_performance(coeff_file_name,fileoutput_number)

        implicit none
        character(len=*),intent(in) :: coeff_file_name
        integer (kind=4),optional:: fileoutput_number
        real (kind=8):: energy
        character (len=9),parameter :: COORD_FORMAT = "Euler_ZYZ"
        integer (kind=4), parameter:: xdim=6
        integer (kind=4):: i,ntest=1000
        real (kind=8):: start, finish
        real (kind=8),dimension(6):: coordinates_set_6D = [ 10.27d0,& !R
                                                            30d0,&  !beta1
                                                            20d0,&  !beta2
                                                            120d0,& !alpha
                                                            40d0,&  !gamma1
                                                            50d0]   !gamma2


        call cpu_time(start)
        
        do i=1,ntest
            call evaluate_LRF( energy,&
                              xdim,&
                              coordinates_set_6D,&
                              COORD_FORMAT,&
                              coeff_file_name)
        end do

        call cpu_time(finish)
        
        write(*,*)"*********************************************************************"
        write(*,*)"* PERFORMANCE for 6D: ",ntest," /time: ",finish-start," *"
        write(*,*)"*********************************************************************"

    end  subroutine running_time_performance


    ! TEST Functions
    subroutine check_energy_MATLAB(system_name,xdim,verbose,fileoutput_number)

        implicit none
        character(len=*),intent(in) :: system_name
        integer (kind=4),intent(in):: xdim,verbose
        integer (kind=4),optional:: fileoutput_number
        integer (kind=4):: i,j,ntest=1000
        real (kind=8):: E0,E1,rmse,Emax,E0_maxval,Erel
        real (kind=8):: coord_from_file(xdim+2)
        real (kind=8), allocatable:: coord(:)
        logical (kind=1):: pass
        real (kind=8), parameter :: PII = DACOS(-1.d0)
        character (len=9),parameter :: COORD_FORMAT = "Euler_ZYZ"
        
        rmse=0d0
        Emax = 0d0
        E0_maxval = 0d0
        Erel = 0d0
        pass = .false.

        open( 17, file = '../testing_datafiles/datasets/'//system_name//'.txt' )
        allocate(coord(xdim))

        do i=1,ntest
            read(17,*)coord_from_file
            E0 = coord_from_file(xdim+2);

            coord(1) = coord_from_file(2)!R
            coord(2) = DACOS(coord_from_file(3))*180d0/PII
            if (xdim==3) then
                coord(3) = coord_from_file(4)*180d0/PII + 90d0
            else
                coord(3) = DACOS(coord_from_file(4))*180d0/PII
                coord(4) = coord_from_file(5)*180d0/PII
                if (xdim>=5) then
                    coord(5) = coord_from_file(6)*180d0/PII + 90d0
                endif
                if (xdim==6) then
                    coord(6) = coord_from_file(7)*180d0/PII + 90d0
                endif
            end if

            call evaluate_LRF(E1,&
                    xdim,&
                    coord,&
                    COORD_FORMAT,&
                    '../testing_datafiles/coefficients/'//system_name//'_Coeff.txt'&
                    )

            rmse = rmse+dabs(E0-E1)**2
            Emax = maxval([Emax,dabs(E0-E1)])
            Erel = Erel+dabs(E0-E1)/dabs(E0)
            E0_maxval = maxval([E0_maxval,dabs(E0)])
        end do

        rmse = dsqrt(rmse/ntest)
        Erel = Erel/ntest


        close(17)

        write(*,*)"*********************************************************************"
        if (verbose == 1) then
            write(*,*)"* System: ",system_name,' - ',xdim," *"
            !write(*,*)"* E0_maxval: ",E0_maxval," *"
            if ((rmse+Erel)/2d0 <= 10d0**(-7)) then
                write(*,*)"* Test: ",char(27)//"[32m"//"Passed!"//char(27)//"[0m"
            else
                write(*,*)"* Test: ",char(27)//"[31m"//"Failure"//char(27)//"[0m"
            end if
            write(*,*)"* Emax: ",Emax," *"
            write(*,*)"* E0_maxval: ",E0_maxval," *"
            write(*,*)"* Erel: ",Erel," *"
            write(*,*)"* rmse: ",rmse," *"
        else
            if ((rmse+Erel)/2d0 <= 10d0**(-7)) then
                write(*,*)"* System: ",system_name," - Test: ",char(27)//"[32m"//"Passed!"//char(27)//"[0m"
            else
                write(*,*)"* System: ",system_name," - Test: ",char(27)//"[31m"//"Failure"//char(27)//"[0m"
            end if

        end if
        deallocate(coord)
    end  subroutine check_energy_MATLAB

end module Testing_v2