module Testing_v2
        use Geometry_Constant_v2
        use FitConstants

        contains
        ! Helper functions         
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
        
        ! TEST Functions
        Subroutine RunningTime_Performance(coeff_filename,fileOutputNumber)
                IMPLICIT NONE
                Character(len = *),Intent(IN) :: coeff_filename
                integer,optional:: fileOutputNumber
                real*8 :: E2,E1
                real*8 ,dimension(6):: GeneralCoordenates,coordinates,zero_coordinates
                real*8 :: FourD_coord(4),SixD_coord(6)
                Character(len = 20) :: coord_format = "Euler_ZYZ"
                INTEGER :: XDim=6,i,ntest=1000
                real*8 :: start, finish
                real*8::testErr(52)

                ! FourD_coord(1) = 10.271963668339600d0 !R
                ! FourD_coord(2) = 30d0  !b1
                ! FourD_coord(3) = 20d0  !b2
                ! FourD_coord(4) = 120d0  !Phi
                ! FourD_coord(3) = 20d0  !b2
 

                SixD_coord(1) = 10.271963668339600d0 !R
                SixD_coord(2) = 30d0  !b1
                SixD_coord(3) = 20d0  !b2
                SixD_coord(4) = 120d0 !Phi
                SixD_coord(5) = 40d0  !c1
                SixD_coord(6) = 50d0  !c2
                

                call cpu_time(start)

           
                do i=1,ntest
                        call evaluateLR(SixD_coord,XDIM,E1,coeff_filename)
                end do 
  
                
                call cpu_time(finish)
                write(*,*)"*********************************************************************"
                write(*,*)"* PERFORMANCE for: ",ntest," /time: ",finish-start," *"
                write(*,*)"*********************************************************************"

        end  Subroutine RunningTime_Performance

end module Testing_v2



PROGRAM main_subroutine

        use Testing_v2
        IMPLICIT NONE


        INTEGER  :: DATE_TIME (8),fileOutNumber
        CHARACTER (LEN = 10) BIG_BEN (3)
        integer:: level_init,level_final,nMAX

        nMAX=10
        level_init=1
        level_final=8
        fileOutNumber=12

        CALL DATE_AND_TIME (BIG_BEN (1), BIG_BEN (2), &
        BIG_BEN (3), DATE_TIME)

        call FileChecking('./files/test/output.test.txt',fileOutNumber)
        REWIND(fileOutNumber)

        write(fileOutNumber,*)"******************************************************************************"
        write(fileOutNumber,*)  "Test Day and Time Record"
        write(fileOutNumber,*)  "Month / Day / Year: ",DATE_TIME(2),"/",DATE_TIME(3),"/",DATE_TIME(1) 
        write(fileOutNumber,*)  "Hr    / Min / Sec : ",DATE_TIME(5),":",DATE_TIME(6),":",DATE_TIME(7) 


        call RunningTime_Performance('./files/test/coefficients/C1(1)_C1(1)_Coeff.txt',fileOutNumber)

! this is useful to test a particular tensor component (Debugging)



     
        close(fileOutNumber)
END PROGRAM main_subroutine






