! Created by albpl on 3/11/2025.
module Testing_v2
  use iso_fortran_env, only : real64, int32
  use LRF_API,         only : evaluate_LRF, user_coordinates_to_general_coordinates
contains

  subroutine file_checking(file_name, num)
    implicit none
    character(*), intent(in) :: file_name
    integer(int32), intent(in) :: num
    logical :: exist

    inquire(file=file_name, exist=exist)
    if (exist) then
      open(num, file=file_name, status="old", position="append", action="write")
    else
      open(num, file=file_name, status="new", action="write")
    end if
  end subroutine file_checking

  subroutine running_time_performance(coeff_file_name, fileoutput_number)
    implicit none
    character(len=*), intent(in) :: coeff_file_name
    integer(int32), optional     :: fileoutput_number
    real(real64) :: energy
    character(len=*), parameter :: COORD_FORMAT = "Euler_ZYZ"
    integer(int32), parameter :: xdim = 6
    integer(int32) :: i, ntest
    real(real64) :: start, finish
    real(real64), dimension(6) :: coordinates_set_6D

    ntest = 1000_int32
    coordinates_set_6D = [ 10.27_real64, 30.0_real64, 20.0_real64, 120.0_real64, 40.0_real64, 50.0_real64 ]

    call cpu_time(start)
    do i = 1, ntest
      call evaluate_LRF( energy,          &
                         xdim,            &
                         coordinates_set_6D, &
                         COORD_FORMAT,    &
                         coeff_file_name )
    end do
    call cpu_time(finish)

    write(*,*)"*********************************************************************"
    write(*,*)"* PERFORMANCE for 6D: ", ntest, " /time: ", finish - start, " *"
    write(*,*)"*********************************************************************"
  end subroutine running_time_performance

  subroutine t_tensor_test()
    use Geometry_Constant_v2, only : tensors_initialization_v2, t_tensor_v3
    use Fitting_Constant_v2,  only : t_index, ensure_t_index_map_ready
    implicit none
    integer(int32) :: i, j, ntest, cpn, order, la, lb, ka, kb
    real(real64) :: general_coordinates_ZXZ(6), r(6), T(9640)
    real(real64), parameter :: PII = acos(-1.0_real64)

    ntest = 1_int32
    cpn   = 1_int32
    call ensure_t_index_map_ready(15_int32)

    do i = 1, ntest
      call random_number(r)
      general_coordinates_ZXZ(1) = 10.0_real64 + r(1)*10.0_real64
      general_coordinates_ZXZ(2) = r(2)*180.0_real64
      general_coordinates_ZXZ(3) = r(3)*180.0_real64
      general_coordinates_ZXZ(4) = r(4)*360.0_real64
      general_coordinates_ZXZ(5) = r(5)*360.0_real64
      general_coordinates_ZXZ(6) = r(6)*360.0_real64

      call tensors_initialization_v2(15_int32, general_coordinates_ZXZ)

      open(unit=10, file="../testing_datafiles/t_tensors/t_tensors_test.txt", action="write")

      do order = 1, 15
        do la = 0, order - 1
          lb = order - la - 1
          do ka = 0, 2*la
            do kb = 0, 2*lb
              T(cpn) = t_tensor_v3(t_index(la+1,ka+1,lb+1,kb+1))
              cpn = cpn + 1
            end do
          end do
        end do
      end do

      write(10,*) general_coordinates_ZXZ, T
      close(10)
    end do
  end subroutine t_tensor_test

  subroutine check_energy_MATLAB(system_name, xdim, verbose, fileoutput_number)
    implicit none
    character(len=*), intent(in) :: system_name
    integer(int32),   intent(in) :: xdim, verbose
    integer(int32),   optional   :: fileoutput_number
    integer(int32) :: i, ntest
    real(real64)  :: E0, E1, rmse, Emax, E0_maxval, Erel
    real(real64)  :: coord_from_file(xdim+2)
    real(real64), allocatable :: coord(:)
    real(real64), parameter :: PII = acos(-1.0_real64)
    character(len=*), parameter :: COORD_FORMAT = "Euler_ZYZ"

    ntest = 1000_int32
    rmse = 0.0_real64; Emax = 0.0_real64; E0_maxval = 0.0_real64; Erel = 0.0_real64

    open(17, file='../testing_datafiles/datasets/'//system_name//'.txt')
    allocate(coord(xdim))

    do i = 1, ntest
      read(17,*) coord_from_file
      E0 = coord_from_file(xdim+2)

      coord(1) = coord_from_file(2)
      coord(2) = acos(coord_from_file(3))*180.0_real64/PII
      if (xdim == 3) then
        coord(3) = coord_from_file(4)*180.0_real64/PII + 90.0_real64
      else
        coord(3) = acos(coord_from_file(4))*180.0_real64/PII
        coord(4) = coord_from_file(5)*180.0_real64/PII
        if (xdim >= 5) coord(5) = coord_from_file(6)*180.0_real64/PII + 90.0_real64
        if (xdim == 6) coord(6) = coord_from_file(7)*180.0_real64/PII + 90.0_real64
      end if

      call evaluate_LRF( E1, xdim, coord, COORD_FORMAT, '../testing_datafiles/coefficients/'//system_name//'_Coeff.txt' )

      rmse = rmse + abs(E0 - E1)**2
      Emax = max(Emax, abs(E0 - E1))
      Erel = Erel + abs(E0 - E1) / abs(E0)
      E0_maxval = max(E0_maxval, abs(E0))
    end do

    rmse = sqrt(rmse/real(ntest, real64))
    Erel = Erel / real(ntest, real64)

    close(17)

    write(*,*)"*********************************************************************"
    if (verbose == 1) then
      write(*,*)"* System: ", system_name, ' - ', xdim, " *"
      if ((rmse + Erel)/2.0_real64 <= 1.0e-7_real64) then
        write(*,*)"* Test: ", char(27)//"[32m"//"Passed!"//char(27)//"[0m"
      else
        write(*,*)"* Test: ", char(27)//"[31m"//"Failure"//char(27)//"[0m"
      end if
      write(*,*)"* Emax: ", Emax, " *"
      write(*,*)"* E0_maxval: ", E0_maxval, " *"
      write(*,*)"* Erel: ", Erel, " *"
      write(*,*)"* rmse: ", rmse, " *"
    else
      if ((rmse + Erel)/2.0_real64 <= 1.0e-7_real64) then
        write(*,*)"* System: ", system_name, " - Test: ", char(27)//"[32m"//"Passed!"//char(27)//"[0m"
      else
        write(*,*)"* System: ", system_name, " - Test: ", char(27)//"[31m"//"Failure"//char(27)//"[0m"
      end if
    end if

    deallocate(coord)
  end subroutine check_energy_MATLAB

end module Testing_v2
