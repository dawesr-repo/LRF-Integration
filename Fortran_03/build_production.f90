subroutine copy_paste(copy_filename, read_unit, writing_unit)
  use iso_fortran_env, only : int32, iostat_end
  implicit none
  character(len=*), intent(in) :: copy_filename
  integer(int32),   intent(in) :: read_unit, writing_unit

  character(len=200) :: command(10000)
  integer(int32)     :: n, i, ios

  ! Open source file for reading on the provided unit
  open(unit=read_unit, file=copy_filename, status='old', action='read', iostat=ios)
  if (ios /= 0) then
    write(*,*) "Error opening file: ", trim(copy_filename), " iostat=", ios
    return
  end if

  n = 0
  do
    read(read_unit, '(A)', iostat=ios) command(n+1)
    if (ios == 0) then
      n = n + 1
    else if (ios == iostat_end) then
      exit
    else
      write(*,*) "Read error on: ", trim(copy_filename), " iostat=", ios
      exit
    end if
  end do

  write(*,*) n, " lines ---- ", trim(copy_filename)

  close(read_unit)

  do i = 1, n
    write(writing_unit, '(A)') trim(command(i))
  end do
end subroutine copy_paste


subroutine file_checking(filename, num)
  use iso_fortran_env, only : int32
  implicit none
  character(*),  intent(in) :: filename
  integer(int32), intent(in) :: num
  logical :: exist

  inquire(file=filename, exist=exist)
  if (exist) then
    open(num, file=filename, status="old", position="append", action="write")
  else
    open(num, file=filename, status="new", action="write")
  end if
end subroutine file_checking


program build_production
  use iso_fortran_env, only : int32
  implicit none

  integer(int32), parameter :: WRITING_UNIT = 101_int32

  type :: string
    character(:), allocatable :: s
  end type string

  integer(int32), parameter :: N = 3_int32
  type(string) :: filenames(N)
  integer(int32) :: i
  integer(int32) :: date_time(8)
  character(len=10) :: big_ben(3)

  call date_and_time (big_ben(1), big_ben(2), big_ben(3), date_time)

  filenames(1)%s = "Fitting_Constant_v2.f90"
  filenames(2)%s = "Geometry_Constant_v2.f90"
  filenames(3)%s = "helper_functions.f90"

  call file_checking('LRF.f90', WRITING_UNIT)
  rewind(WRITING_UNIT)

  write(WRITING_UNIT,*) "!******************************************************************************"
  write(WRITING_UNIT,*) "!      Compilation Day and Time"
  write(WRITING_UNIT,*) "!      Month / Day / Year: ", date_time(2), "/", date_time(3), "/", date_time(1)
  write(WRITING_UNIT,*) "!      Hr    / Min / Sec : ", date_time(5), ":", date_time(6), ":", date_time(7)
  write(WRITING_UNIT,*) "!      LRF MATLAB  v4.1.1"
  write(WRITING_UNIT,*) "!      LRF_Fortran v4.1.1"
  write(WRITING_UNIT,*) "!******************************************************************************"
  write(WRITING_UNIT,*)

  do i = 1, N
    call copy_paste(filenames(i)%s, 1000_int32 + i, WRITING_UNIT)
  end do

  close(WRITING_UNIT)
end program build_production
