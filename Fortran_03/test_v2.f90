program testing_subroutine
  use iso_fortran_env, only : int32
  use Testing_v2, only : file_checking, &
                         running_time_performance, &
                         check_energy_MATLAB, &
                         t_tensor_test
  implicit none

  integer(int32)            :: date_time(8), fileout_number
  character(len=10)         :: big_ben(3)
  integer(int32)            :: n_sys
  integer(int32)            :: xdim_arr(50)
  type :: var_str
    character(len=:), allocatable :: label
  end type var_str

  type(var_str)             :: sys(50)
  type(var_str)             :: symmetries(5)
  integer(int32)            :: i, k, h, count
  integer(int32), dimension(5) :: xdim

  ! init scalars/arrays with kind-correct literals
  fileout_number = 12_int32
  n_sys          = 49_int32
  xdim           = [ 3_int32, 3_int32, 2_int32, 2_int32, 0_int32 ]

  symmetries(1)%label  = "C1(1)"
  symmetries(2)%label  = "Cs(1)"
  symmetries(3)%label  = "D_inf_h(1)"
  symmetries(4)%label  = "C_inf_v(1)"
  symmetries(5)%label  = "Spherical(1)"

  count = 1_int32
  do k = 1, 5
    do h = 1, 5
      if (xdim(k) >= xdim(h)) then
        sys(count)%label = symmetries(k)%label // "_" // symmetries(h)%label
        xdim_arr(count)  = xdim(k) + xdim(h)
        if (xdim(k) == 1_int32 .and. xdim(h) == 1_int32) then
          xdim_arr(count) = 1_int32
        end if
        count = count + 1_int32
      end if
    end do
  end do

  call date_and_time(date=big_ben(1), time=big_ben(2), zone=big_ben(3), values=date_time)

  call file_checking('../testing_datafiles/output.test.txt', fileout_number)
  rewind(fileout_number)

  write(fileout_number,*)"******************************************************************************"
  write(fileout_number,*)  "Test Day and Time Record"
  write(fileout_number,*)  "Month / Day / Year: ", date_time(2), "/", date_time(3), "/", date_time(1)
  write(fileout_number,*)  "Hr    / Min / Sec : ", date_time(5), ":", date_time(6), ":", date_time(7)

  ! Performance run (adjust ARGS in Make to pass a different coeff file if needed)
  call running_time_performance('../testing_datafiles/coefficients/C1(1)_C1(1)_Coeff.txt', fileout_number)

  ! Optional: generate one t-tensor dump using the new flattened storage
  ! call t_tensor_test()

  ! Compare energies vs MATLAB datasets for the first 12 systems generated above
  do i = 1, 12
    call check_energy_MATLAB( sys(i)%label, &
                              xdim_arr(i),   &
                              0_int32,       &
                              fileout_number )
  end do

  close(fileout_number)
end program testing_subroutine
