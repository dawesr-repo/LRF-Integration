subroutine copy_paste(copy_filename, read_unit, writing_unit)

    implicit none
    character(len=*), intent(in) ::copy_filename
    integer(kind = 4), intent(in) :: read_unit, writing_unit
   
    character(len=200) :: command(10000)
    integer(kind = 4):: n, i,ios

    open(unit=read_unit, file=copy_filename, iostat=ios)
    if ( ios /= 0 ) stop "Error opening file data.dat"

    n = 0

    do
        read(read_unit, '(A)', iostat=ios) command(n+1)
        if (ios /= 0) exit
        n = n + 1
    end do

    print*,  n, "lines"," ---- ",copy_filename

    close(read_unit)

    do i = 1, n
        write(writing_unit,*) TRIM(command(i))
    end do

end subroutine copy_paste


subroutine file_checking(filename,num)
    implicit none
    character(*),intent(in)::filename
    integer(kind = 4), intent(in)::num
    logical (kind = 1):: exist

    inquire(file=filename, exist=exist)
    if (exist) then
        open(num, file=filename, status="old", position="append", action="write")
    else
        open(num, file=filename, status="new", action="write")
    end if

end subroutine file_checking



PROGRAM build_production


    implicit none
    integer(kind = 4), parameter :: WRITING_UNIT = 101

    type string
        character(:), allocatable :: s
    end type

    integer(kind = 4), parameter :: N = 3
    type(string) :: filenames(N)
    integer(kind = 4) :: i,if1,if2,if3,if4,if5
    integer (kind = 4):: date_time (8)
    character (len = 10) big_ben (3)

    call date_and_time (big_ben (1), big_ben (2), big_ben (3), date_time)


    filenames(1)%s = "Fitting_Constant_v2.f90"
    filenames(2)%s = "Geometry_Constant_v2.f90"
    filenames(3)%s = "helper_functions.f90"


    call file_checking('LRF.f90',WRITING_UNIT)
    rewind(WRITING_UNIT)

    write(WRITING_UNIT,*)"!******************************************************************************"
    write(WRITING_UNIT,*)  "!      Compilation Day and Time"
    write(WRITING_UNIT,*)  "!      Month / Day / Year: ",date_time(2),"/",date_time(3),"/",date_time(1)
    write(WRITING_UNIT,*)  "!      Hr    / Min / Sec : ",date_time(5),":",date_time(6),":",date_time(7)
    write(WRITING_UNIT,*)  "!      LRF MATLAB  v4.1.1"
    write(WRITING_UNIT,*)  "!      LRF_Fortran v4.1.1"
    write(WRITING_UNIT,*)"!******************************************************************************"
    write(WRITING_UNIT,*)
    do i = 1, N
        call copy_paste(filenames(i)%s,1000+i,WRITING_UNIT)
    end do

    close(WRITING_UNIT)

end PROGRAM build_production