PROGRAM main_subroutine
    IMPLICIT NONE

    character(len=20) :: bufferc,bufferd
    logical :: exist
    character(len = 22), parameter:: dirData = "./files/test/datasets/"
    character(len = 26), parameter:: dirCoeff="./files/test/coefficients/"
    INTEGER::indCoeff,indData
    

    
    indCoeff = 1
    fileloop: do
            write(bufferc,"(A,I3.3,A)") "coefficients_", indCoeff, ".txt"
            inquire(file= dirCoeff//bufferc, exist=exist)

            if (exist) then

                    write(*,*) "File: '", trim(bufferc), "' found."

                    indData = 1
                    filedloop: do
                            write(bufferd,"(A,I3.3,A,I3.3,A)") "datatest_",indCoeff,"_", indData, ".txt"
                            inquire(file= dirData//bufferd, exist=exist)

                            if (exist) then
                                    write(*,*) "     File: '", trim(bufferd), "' found."
                                    indData = indData + 1
                            else 
                                   if (indData==1) then
                                        write(*,*) "    No Dataset Found for coefficients: '", trim(bufferc)
                                   end if
                                    exit
                            end if
                        
                    end do filedloop
                    indCoeff = indCoeff + 1
            else
                    exit
            end if
        
    end do fileloop
    
END PROGRAM main_subroutine


SUBROUTINE test_subroutine(dir,filename)
    IMPLICIT NONE
    Character(len = *),Intent(IN) :: dir,filename
    real*8 :: E2,E1
    real*8 ,dimension(6):: GeneralCoordenates,coordinates,zero_coordinates
    real*8 :: FourD_coord(4),TwoD_coord(2)
    Character(len = 20) :: coord_format = "Euler_ZXZ"
    INTEGER :: XDim=4,i
    real*8 :: start, finish,pii
    real*8 :: ind, R,cos_b1,cos_b2,phi, E0
    
    pii=DAcos(-1d0)

    zero_coordinates(1) = 0d0  !R
    zero_coordinates(2) = 0d0  !b1
    zero_coordinates(3) = 0d0  !b2
    zero_coordinates(4) = 0d0  !phi
    zero_coordinates(5) = 0d0  !c1
    zero_coordinates(6) = 0d0  !c2


    ! CALL Long_Range_Potential(zero_coordinates,E2)
    ! Open( 20, file = './files/data.dat' )
    
    

    ! do i=1,10
    !     read(20,*) ind, R,cos_b1,cos_b2,phi, E0

    !     if (R>10d0)Then
    !         FourD_coord(1) = R !R
    !         FourD_coord(2) = dacos(cos_b1) *180d0/pii !b1
    !         FourD_coord(3) = Dacos(cos_b2)*180d0/pii  !b2
    !         FourD_coord(4) = phi*180d0/pii  !Phi


    ! Call General_Coordinates_Format(XDim, FourD_coord, GeneralCoordenates)
    ! Call Coordinate_Transformation(GeneralCoordenates,coord_format,coordinates)
    ! CALL Long_Range_Potential(coordinates,E1)

    ! write(*,*)(E0-E2)-E1

    ! !call cpu_time(start)

    ! ! do i=1,100000

    ! !     Call General_Coordinates_Format(XDim, TwoD_coord, GeneralCoordenates)
    ! !     Call Coordinate_Transformation(GeneralCoordenates,coord_format,coordinates)
     
    ! !     !write(*,*)coordinates
    ! !     CALL Long_Range_Potential(coordinates,E1)
    ! !     !write(*,*)E1

    ! ! enddo    
  
    ! !call cpu_time(finish)
    ! !print '("Time = ",f6.3," seconds.")',finish-start

    !     end if
    ! end do
    ! close(20) 



    
END SUBROUTINE test_subroutine



