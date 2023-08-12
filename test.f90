module Testing

        contains

        Subroutine READ_DATASET(filename,systName,ntest,XDim,DATASET_,coord_format,M_Fit ,D_Fit,I_Fit,H_Fit)
                
                IMPLICIT NONE
                Character(len = *),Intent(IN) :: filename
                Character(len = 20),Intent(out) :: coord_format,systName
                INTEGER ,Intent(out) :: XDim,ntest    
                real*8, allocatable,INTENT(INOUT) ::DATASET_(:,:)
                
                Integer, dimension (8),  INTENT(INOUT) :: M_Fit      
                Integer, dimension (3),  INTENT(INOUT) :: D_Fit     
                Integer, dimension (5),  INTENT(INOUT) :: I_Fit     
                Integer, dimension (2),  INTENT(INOUT) :: H_Fit     

                integer:: i
                Character(len = 40) :: row

                Open( 100, file = filename )
                Read( 100, *) row
                Read( 100, *) row,row,row,ntest
                Read( 100, *) row
                Read( 100, *) row,row,row,systName
                Read( 100, *) row
                Read( 100, *) row,row,row,coord_format
                Read( 100, *) row,row,XDim
                Read( 100, *) row


                Read( 100, *) row

                read(100, *)  M_Fit
                read(100, *)  D_Fit
                read(100, *)  I_Fit
                read(100, *)  H_Fit

                Read( 100, *) row
                Read( 100, *) row
                Read( 100, *) row
                Read( 100, *) row

                Allocate(DATASET_(ntest,29))
          

                do i=1,ntest
                        Read( 100, *)DATASET_(i,:)
                enddo

                close(10) 
        end Subroutine READ_DATASET

        Subroutine Test_Dataset(coeff_filename, data_filename,fileOutputNumber,ntestMAX)
                IMPLICIT NONE
                Character(len = *),Intent(IN) :: coeff_filename, data_filename
                integer,optional:: fileOutputNumber,ntestMAX

                real*8 :: E_fortran,E_fit
                real*8 ,dimension(6):: GeneralCoordenates,coordinates,zero_coordinates
                real*8, allocatable :: XDim_coord(:),dataset_(:,:)
                Character(len = 20) :: coord_format,systName
                Character(len = 2) :: str
                INTEGER :: XDIM,ntest,nMAX
                real*8 :: start, finish,pii
                real*8 :: ind, R,cos_b1,cos_b2,phi, E0
                Integer, dimension (8) :: M_Fit      
                Integer, dimension (3) :: D_Fit     
                Integer, dimension (5) :: I_Fit     
                Integer, dimension (2) :: H_Fit   

                integer:: i,success_test,fOutputNum,j
                real*8::testArr(52),errArr(52),maxErr,err_tol,diff_comp(52)

                if (present(fileOutputNumber))then
                        fOutputNum = fileOutputNumber
                else
                        fOutputNum = 0
                end if
                
                

                maxErr=0d0
                err_tol = 1d-6
                pii=DAcos(-1d0)
                diff_comp =0d0

                Call READ_DATASET(data_filename,systName,ntest,XDIM,dataset_,coord_format,M_Fit ,D_Fit,I_Fit,H_Fit)

                Allocate(XDim_coord(XDIM))

                success_test = 0

                if (present(ntestMAX))then
                        nMAX = ntestMAX
                else
                        nMAX = ntest
                end if

                do i=1,nMAX

                     if (XDIM==3)Then
                        
                     else
                        XDim_coord = dataset_(i,1:XDIM)
                        E_fit = dataset_(i,7)

                        Call General_Coordinates_Format(XDIM, XDim_coord, GeneralCoordenates)
                        Call Coordinate_Transformation(GeneralCoordenates,coord_format,coordinates)
                        CALL Long_Range_Potential(coordinates,E_fortran,coeff_filename,1,testArr)

                        diff_comp(1) = DABS(E_fortran-E_fit)
                        diff_comp(2) = DABS(dataset_(i,8)-testArr(2))
                        diff_comp(3) = DABS(dataset_(i,10)-testArr(3))
                        diff_comp(4) = DABS(dataset_(i,9)-testArr(4))
                        diff_comp(5) = DABS(dataset_(i,11)-testArr(5))

                        diff_comp(6:13) = DABS(dataset_(i,12:19)-testArr(6:13))  
                        diff_comp(21:23) = DABS(dataset_(i,20:22)-testArr(21:23))  
                        diff_comp(31:35) = DABS(dataset_(i,23:27)-testArr(31:35)) 
                        diff_comp(43:44) = DABS(dataset_(i,28:29)-testArr(43:44)) 

                        write(*,*)testArr(21),dataset_(i,20)

                        if(MAXVAL(diff_comp)>maxErr)then
                            maxErr = MAXVAL(diff_comp)
                        end if

                        if(MAXVAL(diff_comp)>err_tol )then
                                
                                if ( fOutputNum < 0 )then
                                        write(*,*) 
                                        write(*,*) "Row : ",i
                                        write(*,*) "Coordinates : ",XDim_coord
                                        write(*,*) "Energy: ",E_fortran,E_fit,diff_comp(1)
                                        write(*,*) "Elect. : ",dataset_(i,8),testArr(2), diff_comp(2)
                                        write(*,*) "Dispe. : ",dataset_(i,10),testArr(3),diff_comp(3)
                                        write(*,*) "Induc. : ",dataset_(i,9),testArr(4), diff_comp(4)
                                        write(*,*) "Hyper. : ",dataset_(i,11),testArr(5),diff_comp(5)
                                        write(*,*) 


                                elseif ( fOutputNum > 0 )then

                                write(fOutputNum,*) 
                                write(fOutputNum,*) "Row : ",i
                                write(fOutputNum,*) "Coordinates : ",XDim_coord
                                write(fOutputNum,*) "Energy: ",E_fortran,E_fit, diff_comp(1)
                                write(fOutputNum,*) "Elect. : ",dataset_(i,8),testArr(2), diff_comp(2)
                                write(fOutputNum,*) "Dispe. : ",dataset_(i,10),testArr(3), diff_comp(3)
                                write(fOutputNum,*) "Induc. : ",dataset_(i,9),testArr(4),  diff_comp(4)
                                write(fOutputNum,*) "Hyper. : ",dataset_(i,11),testArr(5), diff_comp(5)
                                write(fOutputNum,*) 
                                write(fOutputNum,*) "By Components: "

                                do j=6,52
                                        if (diff_comp(j) > err_tol )then

                                                if (j>5 .and. j<21)then
                                                        write(str,'(A,I1)')"M",j-5
                                                elseif(j>20 .and. j<31)then
                                                        write(str,'(A,I1)')"D",j-20+5
                                                elseif(j>30 .and. j<43)then
                                                        write(str,'(A,I1)')"I",j-30+3  
                                                elseif(j>42)then
                                                        write(str,'(A,I1)')"H",j-42+5               
                                                end if
                                                
                                                
                                                write(fOutputNum,'(A,E15.3 )') str,diff_comp(j)
                                        end if 
                                end do 
                                write(fOutputNum,*) "*****************************************"
                                

                                end if 
                        else
                                success_test =success_test+1
                        end if
                   
                        
                     end if 



                enddo

                Call Print_Test_Results(systName,maxErr,nMAX,success_test,fOutputNum)

        end Subroutine Test_Dataset

        SUBROUTINE Print_Test_Results(testName,maxErr,ntest,success_test,fileOutputNumber)

                IMPLICIT NONE
                character(*),INTENT(IN) :: testName
                integer,INTENT(IN)::ntest,success_test
                character (len =25 ) ::f_res1,f_res2
                character (len =7) :: str,str1
                integer::failed,ptg
                real*8,INTENT(IN)::maxErr
                integer,optional:: fileOutputNumber

                f_res1 = achar(27)//'[32m TEST PASSED '//achar(27)//'[0m'
                f_res2 = achar(27)//'[31m TEST FAILED :('//achar(27)//'[0m'

                ptg = FLOOR(success_test*100d0/ntest)
                failed = ntest-success_test
                write (str1, '(I4,A)') ptg,"%"

                write(*,  *)
                if (ntest > success_test)then
        
                write(*, '(A, I5, A, I1, A, I5, A)') f_res2//"("//achar(27)//'[31m'//str1//achar(27)//'[0m)'//" ----  "&
                        //achar(27)//'[35m'//testName//achar(27)//'[0m'
                else
                write(*,  *) f_res1//"("//achar(27)//'[32m 100% '//achar(27)//'[0m)'//" ----  "&
                //achar(27)//'[35m'//testName//achar(27)//'[0m)'                 
                end if
                write(*,  '(A, I8, A, I8, A, I8)') "                                Num. Test: ",ntest,"     ---> success: "&
                                                ,success_test," / failed : ",failed
                write(*,  '(A, F21.15)') "                                Error Max: ",maxErr                                 
                write(*,  *)

                if (present(fileOutputNumber))then
                         

                        ! write(*,  *)
                        ! if (ntest > success_test)then
                
                        ! write(*, '(A, I5, A, I1, A, I5, A)') f_res2//"("//achar(27)//'[31m'//str1//achar(27)//'[0m)'//" ----  "&
                        !         //achar(27)//'[35m'//testName//achar(27)//'[0m'
                        ! else
                        ! write(*,  *) f_res1//"("//achar(27)//'[32m 100% '//achar(27)//'[0m)'//" ----  "&
                        ! //achar(27)//'[35m'//testName//achar(27)//'[0m)'                 
                        ! end if
                        ! write(*,  '(A, I8, A, I8, A, I8)') "                                Num. Test: ",ntest,"     ---> success: "&
                        !                                 ,success_test," / failed : ",failed
                        ! write(*,  '(A, F21.15)') "                                Error Max: ",maxErr                                 
                        ! write(*,  *)
                end if

        end SUBROUTINE Print_Test_Results

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


        SUBROUTINE Test_All(fileOutNumber)
                IMPLICIT NONE
                INTEGER ,Intent(In)::fileOutNumber
                character(len=20) :: bufferc,bufferd
                logical :: exist
                character(len = 22), parameter:: dirData = "./files/test/datasets/"
                character(len = 26), parameter:: dirCoeff="./files/test/coefficients/"
                INTEGER  ::indCoeff,indData

                 
                indCoeff = 1

                fileloop: do
                        
                        write(bufferc,"(A,I3.3,A)") "coefficients_", indCoeff, ".txt"
                        inquire(file= dirCoeff//bufferc, exist=exist)
                        
                       

                        if (exist) then

                                !write(*,*) "File: '", trim(bufferc), "' found."

                                indData = 1
                                filedloop: do
                                write(bufferd,"(A,I3.3,A,I3.3,A)") "datatest_",indCoeff,"_", indData, ".txt"
                                inquire(file= dirData//bufferd, exist=exist)



                                if (exist) then
                                        !write(*,*) "     File: '", trim(bufferd), "' found."
                                        call Test_Dataset(dirCoeff//bufferc, dirData//bufferd,fileOutNumber)
                                        
                                        
                                        indData = indData + 1
                                else 
                                        if (indData==1) then
                                                write(fileOutNumber,*)
                                                write(fileOutNumber,*) "*******************************************************"
                                                write(fileOutNumber,*) "*    No Dataset Found for coefficients: ", trim(bufferc)
                                                write(fileOutNumber,*) "*******************************************************"
                                                write(fileOutNumber,*)
                                        end if
                                        exit
                                end if
                                        
                                end do filedloop
                                indCoeff = indCoeff + 1
                        else
                                exit
                        end if

                       
               
                end do fileloop
        END SUBROUTINE Test_All

        SUBROUTINE Testing_TTensors()
                IMPLICIT NONE


                call init_Tensors() ! Initializing in zero the new vectors

                 
                
        END SUBROUTINE Testing_TTensors 
end module Testing



PROGRAM main_subroutine

        use Testing
        IMPLICIT NONE


        INTEGER  :: DATE_TIME (8),fileOutNumber
        CHARACTER (LEN = 10) BIG_BEN (3)
    
        fileOutNumber=12

        CALL DATE_AND_TIME (BIG_BEN (1), BIG_BEN (2), &
        BIG_BEN (3), DATE_TIME)

        call FileChecking('./files/test/output.test.txt',fileOutNumber)
        REWIND(fileOutNumber)

        write(fileOutNumber,*)"******************************************************************************"
        write(fileOutNumber,*)  "Test Day and Time Record"
        write(fileOutNumber,*)  "Month / Day / Year: ",DATE_TIME(2),"/",DATE_TIME(3),"/",DATE_TIME(1) 
        write(fileOutNumber,*)  "Hr    / Min / Sec : ",DATE_TIME(5),":",DATE_TIME(6),":",DATE_TIME(7) 


        !call Test_All(fileOutNumber)
        call Test_Dataset("./files/test/coefficients/coefficients_002.txt", "./files/test/datasets/datatest_002_001.txt"&
                        ,fileOutNumber,3)

        close(fileOutNumber)
END PROGRAM main_subroutine






