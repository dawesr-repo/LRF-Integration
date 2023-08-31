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
        write(writing_unit,*) command(i)
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

        integer, parameter :: N = 25
        type(string) :: fileNames(N)
        integer :: i,if1,if2,if3,if4,if5
        INTEGER  :: DATE_TIME (8)
        CHARACTER (LEN = 10) BIG_BEN (3)
        
        CALL DATE_AND_TIME (BIG_BEN (1), BIG_BEN (2), &
        BIG_BEN (3), DATE_TIME)
        

        
        fileNames(1)%s = "T_saved.f90"
        fileNames(2)%s = "FittingConstant.f90"
        fileNames(3)%s = "T_lk.f90"
        fileNames(4)%s = "T_l0.f90"
        fileNames(5)%s = "T_ll.f90"
        fileNames(6)%s = "Index_Searcher.f90"


        if1 = 6
        
        fileNames(if1+1)%s = "Approx_1_Sph2.f90"
        fileNames(if1+2)%s = "Approx_2_Sph2.f90"
        fileNames(if1+3)%s = "Approx_3_Sph2.f90"
        fileNames(if1+4)%s = "Approx_4_Sph2.f90"
        fileNames(if1+5)%s = "Approx_5_Sph2.f90"
        fileNames(if1+6)%s = "Approx_6_Sph2.f90"
        fileNames(if1+7)%s = "Approx_7_Sph2.f90"
        fileNames(if1+8)%s = "Approx_8_Sph2.f90"

        if2 = if1+8

        fileNames(if2+1)%s = "Induction_4_Sph2.f90"
        fileNames(if2+2)%s = "Induction_5_Sph2.f90"
        fileNames(if2+3)%s = "Induction_6_Sph2.f90"
        fileNames(if2+4)%s = "Induction_7_Sph2.f90"
        fileNames(if2+5)%s = "Induction_8_Sph2.f90"

        if3 = if2+5

        fileNames(if3+1)%s = "Dispersion_6_Sph2.f90"
        fileNames(if3+2)%s = "Dispersion_7_Sph2.f90"
        fileNames(if3+3)%s = "Dispersion_8_Sph2.f90"

        if4 = if3+3

        fileNames(if4+1)%s = "HyperPolarizability_6_Sph2.f90"
        fileNames(if4+2)%s = "HyperPolarizability_7_Sph2.f90"

        if5 = if4+2

        fileNames(if5+1)%s = "helperFunc.f90"
     

        call FileChecking('LRF.f90',writing_unit)
        REWIND(writing_unit)

        write(writing_unit,*)"!******************************************************************************"
        write(writing_unit,*)  "!      Compilation Day and Time"
        write(writing_unit,*)  "!      Month / Day / Year: ",DATE_TIME(2),"/",DATE_TIME(3),"/",DATE_TIME(1) 
        write(writing_unit,*)  "!      Hr    / Min / Sec : ",DATE_TIME(5),":",DATE_TIME(6),":",DATE_TIME(7)
        write(writing_unit,*)  "!      LRF MATLAB  v0.7.1"
        write(writing_unit,*)  "!      LRF_Fortran v0.3.1"  
        write(writing_unit,*)"!******************************************************************************"
        write(writing_unit,*)
        do i = 1, N
                Call Copy_Paste(fileNames(i)%s,1000+i,writing_unit) 
        end do

        close(writing_unit)

END PROGRAM build_production






