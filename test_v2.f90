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


        ! TEST Functions
        Subroutine Check_MATLAB_ENERGY(SystemName,XDIM,verbose,fileOutputNumber)
                IMPLICIT NONE
                Character(len = *),Intent(IN) :: SystemName
                Integer,Intent(IN):: XDIM,verbose
                integer,optional:: fileOutputNumber
                INTEGER :: i,j,ntest=1000
                real*8 :: E0,E1,RMSE,Emax,E0_MAXVAL,Erel
                real*8 :: coord(XDIM),coord_from_file(XDim+2)
                logical :: pass
                real*8, parameter :: pii = DACOS(-1.d0) 

                RMSE=0d0
                Emax = 0d0   
                E0_MAXVAL = 0d0     
                Erel = 0d0 
                pass = .false.
               
                Open( 17, file = './files/test/datasets/'//SystemName//'.txt' )
           
                do i=1,ntest
                        read(17,*)coord_from_file
                        E0 = coord_from_file(XDim+2);

                        coord(1) = coord_from_file(2)!R
                        coord(2) = DACOS(coord_from_file(3))*180d0/pii
                        if (XDIM==3) then
                                coord(3) = coord_from_file(4)*180d0/pii + 90d0
                        else
                                coord(3) = DACOS(coord_from_file(4))*180d0/pii
                                coord(4) = coord_from_file(5)*180d0/pii
                                if (XDIM>=5) then        
                                        coord(5) = coord_from_file(6)*180d0/pii + 90d0
                                endif
                                if (XDIM==6) then        
                                        coord(6) = coord_from_file(7)*180d0/pii + 90d0
                                endif
                        end if
            
                        call evaluateLR(coord,XDIM,E1,'./files/test/coefficients/'//SystemName//'_Coeff.txt')
                        
                        RMSE = RMSE+dabs(E0-E1)**2
                        Emax = MAXVAL([Emax,dabs(E0-E1)])
                        Erel = Erel+dabs(E0-E1)/dabs(E0)
                        E0_MAXVAL = MAXVAL([E0_MAXVAL,dabs(E0)])
                end do 

                RMSE = Dsqrt(RMSE/ntest)
                Erel = Erel/ntest

             

                


                close(17)
           
                write(*,*)"*********************************************************************"
                if (verbose == 1) then
                        write(*,*)"* System: ",SystemName,' - ',xdim," *"
                        !write(*,*)"* E0_MAXVAL: ",E0_MAXVAL," *"
                        if ((RMSE+Erel)/2d0 <= 10d0**(-7)) then
                                write(*,*)"* Test: ",char(27)//"[32m"//"Passed!"//char(27)//"[0m"
                        else
                                write(*,*)"* Test: ",char(27)//"[31m"//"Failure"//char(27)//"[0m"
                        end if
                        write(*,*)"* Emax: ",Emax," *"
                        write(*,*)"* E0_MAXVAL: ",E0_MAXVAL," *"
                        write(*,*)"* Erel: ",Erel," *"
                        write(*,*)"* RMSE: ",RMSE," *"
                else
                        if ((RMSE+Erel)/2d0 <= 10d0**(-7)) then
                                write(*,*)"* System: ",SystemName," - Test: ",char(27)//"[32m"//"Passed!"//char(27)//"[0m"
                        else
                                write(*,*)"* System: ",SystemName," - Test: ",char(27)//"[31m"//"Failure"//char(27)//"[0m"
                        end if
                        
                end if
            
        end  Subroutine Check_MATLAB_ENERGY

end module Testing_v2





PROGRAM main_subroutine

        use Testing_v2
        IMPLICIT NONE


        INTEGER  :: DATE_TIME (8),fileOutNumber
        CHARACTER (LEN = 10) BIG_BEN (3)
        integer:: n_sys = 7**2
        integer :: xdim_arr(50)
        Type :: varStr
                character(len=:), allocatable :: as_str
        end Type varStr
 
        type(varStr), dimension(1) :: Sys(50)
        type(varStr), dimension(1) :: numStr(15)
        integer :: i,k,h,count,xdim(7)


        
        numStr(1)%as_str  = "C1(1)"
        numStr(2)%as_str  = "T(2)"
        numStr(3)%as_str  = "Ih(2)"
        numStr(4)%as_str  = "Cs(1)"

        numStr(5)%as_str  = "D_inf_h(1)"
        numStr(6)%as_str  = "C_inf_v(1)"
        numStr(7)%as_str  = "Spherical(1)"

        

        xdim= (/3,3,3,3,2,2,0/)

        count = 1
        do k=1,7
           do h=1,7
                if (xdim(k) >= xdim(h)) then
                Sys(count)%as_str = numStr(k)%as_str//"_"//numStr(h)%as_str
                xdim_arr(count) = xdim(k)+ xdim(h)
                if (xdim(k)==1 .and. xdim(h)==1) then
                         xdim_arr(count) = 1
                end if       
                count = count + 1        
                end if
           end do
        end do


        fileOutNumber=12


        CALL DATE_AND_TIME (BIG_BEN (1), BIG_BEN (2), &
        BIG_BEN (3), DATE_TIME)

        call FileChecking('./files/test/output.test.txt',fileOutNumber)
        REWIND(fileOutNumber)

        write(fileOutNumber,*)"******************************************************************************"
        write(fileOutNumber,*)  "Test Day and Time Record"
        write(fileOutNumber,*)  "Month / Day / Year: ",DATE_TIME(2),"/",DATE_TIME(3),"/",DATE_TIME(1) 
        write(fileOutNumber,*)  "Hr    / Min / Sec : ",DATE_TIME(5),":",DATE_TIME(6),":",DATE_TIME(7) 


        ! call RunningTime_Performance('./files/test/coefficients/C1(1)_C1(1)_Coeff.txt',fileOutNumber)
       
        do i=1,count-1
             !print*,"SYSNAME",i,Sys(i)%as_str
             call Check_MATLAB_ENERGY(Sys(i)%as_str,xdim_arr(i),0,fileOutNumber)
        end do

        close(fileOutNumber)
END PROGRAM main_subroutine






