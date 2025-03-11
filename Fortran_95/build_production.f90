Subroutine Copy_Paste(copy_filename,read_unit,writing_unit)

    implicit none
    character(len=*), intent(IN) ::copy_filename
    integer, intent(IN) :: read_unit,writing_unit

    integer :: ios
    character(len=200) :: command(10000)
    integer :: n, i

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

end Subroutine Copy_Paste


SUBROUTINE FileChecking(fileName,num)
    implicit none
    character(*),INTENT(IN)::fileName
    Integer,INTENT(IN)::num
    logical :: exist

    inquire(file=fileName, exist=exist)
    if (exist) then
        open(num, file=fileName, status="old", position="append", action="write")
    else
        open(num, file=fileName, status="new", action="write")
    end if

END SUBROUTINE FileChecking



PROGRAM build_production


    IMPLICIT NONE
    integer, parameter :: writing_unit = 101

    type string
        character(:), allocatable :: s
    end type

    integer, parameter :: N = 6
    type(string) :: fileNames(N)
    integer :: i,if1,if2,if3,if4,if5
    INTEGER  :: DATE_TIME (8)
    CHARACTER (LEN = 10) BIG_BEN (3)

    CALL DATE_AND_TIME (BIG_BEN (1), BIG_BEN (2), &
            BIG_BEN (3), DATE_TIME)



    fileNames(1)%s = "Geometry_Constant_v2.f90"
    fileNames(2)%s = "FittingConstant.f90"
    fileNames(3)%s = "Multipole_Sph3.f90"
    fileNames(4)%s = "Induction_Sph3.f90"
    fileNames(5)%s = "Dispersion_Sph3.f90"
    fileNames(6)%s = "helperFunc.f90"


    call FileChecking('LRF.f90',writing_unit)
    REWIND(writing_unit)

    write(writing_unit,*)"!******************************************************************************"
    write(writing_unit,*)  "!      Compilation Day and Time"
    write(writing_unit,*)  "!      Month / Day / Year: ",DATE_TIME(2),"/",DATE_TIME(3),"/",DATE_TIME(1)
    write(writing_unit,*)  "!      Hr    / Min / Sec : ",DATE_TIME(5),":",DATE_TIME(6),":",DATE_TIME(7)
    write(writing_unit,*)  "!      LRF MATLAB  v4.1.1"
    write(writing_unit,*)  "!      LRF_Fortran v4.1.1"
    write(writing_unit,*)"!******************************************************************************"
    write(writing_unit,*)
    do i = 1, N
        Call Copy_Paste(fileNames(i)%s,1000+i,writing_unit)
    end do

    close(writing_unit)

END PROGRAM build_production