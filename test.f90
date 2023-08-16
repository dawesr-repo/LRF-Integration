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

                close(100) 
        end Subroutine READ_DATASET

        Subroutine READ_GP(filename,ntest,Ar,Br,C)
                
                IMPLICIT NONE
                Character(len = *),Intent(IN) :: filename
                INTEGER ,Intent(IN) :: ntest    

                real*8 , INTENT(INOUT):: Ar(ntest,3)
                real*8 , INTENT(INOUT):: Br(ntest,3)
                real*8 , INTENT(INOUT):: C(ntest,9)
   

                integer:: i
                Character(len = 40) :: row

                Open( 400, file = filename )
                Read( 400, *) row

                do i=1,ntest
                        Read( 400, *)Ar(i,:),Br(i,:),C(i,:)
                enddo

                close(400) 
        end Subroutine READ_GP

        Subroutine Read_Check_val(filename,Check_val,ntest)
                
                IMPLICIT NONE
                Character(len = *),Intent(IN) :: filename
                INTEGER ,Intent(IN) :: ntest    

                real*8 , INTENT(INOUT):: Check_val(ntest)
            
   

                integer:: i
                Character(len = 40) :: row

                Open( 400, file = filename )
                Read( 400, *) row

                do i=1,ntest
                        Read( 400, *)Check_val(i)
                enddo

                close(400) 
        end Subroutine Read_Check_val

         Subroutine READ_DATASET_TTensor(filename,systName,ntest,DATASET_,coord_format)
                
                IMPLICIT NONE
                Character(len = *),Intent(IN) :: filename
                Character(len = 20),Intent(out) :: coord_format,systName
                INTEGER ,Intent(out) :: ntest    
                real*8, allocatable,INTENT(INOUT) ::DATASET_(:,:)
                

                integer:: i
                Character(len = 40) :: row

                Open( 200, file = filename )
                Read( 200, *) row
                Read( 200, *) row,row,row,ntest
                Read( 200, *) row
                Read( 200, *) row,row,row,systName
                Read( 200, *) row
                Read( 200, *) row,row,row,coord_format
                Read( 200, *) row


                Read( 200, *) row
                Read( 200, *) row

                Allocate(DATASET_(ntest,5))
          

                do i=1,ntest
                        Read( 200, *)DATASET_(i,:)
                enddo

                close(200) 
        end Subroutine READ_DATASET_TTensor

         Subroutine READ_TTensor_Values(ntest,la,lb,DATASET_)
                
                IMPLICIT NONE
                INTEGER ,Intent(IN) :: ntest   ,la,lb 
                real*8 ::DATASET_(ntest,(2*7+1)*(2*7+1))
                Character(len = 23),parameter :: dirTTensor = "./files/test/T_Tensors/"
                logical :: exist
                integer:: i
                Character(len = 40) :: row
                type string
                        character(len=:), allocatable :: s
                end type string
                ! create an array of strings where each element is separately allocatable
                type(string) :: mpArr(13)
       
                mpArr(1)%s = 'q'   
                mpArr(2)%s = 'm'
                mpArr(3)%s = 'Qd'
                mpArr(4)%s = 'O'
                mpArr(5)%s = 'Phi'
                mpArr(6)%s = 'M5'
                mpArr(7)%s = 'M6'
                mpArr(8)%s = 'M7'
                mpArr(9)%s = 'M8'
                mpArr(10)%s = 'M9'
                mpArr(11)%s = 'M10'
                mpArr(12)%s = 'M11'
                mpArr(13)%s = 'M12'


                inquire(file= dirTTensor//mpArr(la+1)%s//mpArr(lb+1)%s//'.txt', exist=exist)
                if (exist)then 
                        Open( 200, file = dirTTensor//mpArr(la+1)%s//mpArr(lb+1)%s//'.txt' )
                        Read( 200, *) row
                        Read( 200, *) row
                        Read( 200, *) row
                        Read( 200, *) row
                        Read( 200, *) row
                        Read( 200, *) row


                        Read( 200, *) row
                        Read( 200, *) row

        

                        do i=1,ntest
                                Read( 200, *)DATASET_(i,1:(2*la+1)*(2*lb+1))
                        enddo

                        close(200) 
                else
                        write(*,*) "file locate at "//dirTTensor//mpArr(la+1)%s//mpArr(lb+1)%s//'.txt' //" DOES NOT EXIST"
                end if 


        end Subroutine READ_TTensor_Values

        Subroutine Test_Dataset(coeff_filename, data_filename,fileOutputNumber,ntestMIN,ntestMAX)
                IMPLICIT NONE
                Character(len = *),Intent(IN) :: coeff_filename, data_filename
                integer,optional:: fileOutputNumber,ntestMIN,ntestMAX

                real*8 :: E_fortran,E_fit
                real*8 ,dimension(6):: GeneralCoordenates,coordinates,zero_coordinates
                real*8, allocatable :: XDim_coord(:),dataset_(:,:)
                Character(len = 20) :: coord_format,systName
                Character(len = 2) :: str
                INTEGER :: XDIM,ntest,nMAX,nMIN
                real*8 :: start, finish,pii
                real*8 :: ind, R,cos_b1,cos_b2,phi, E0
                Integer, dimension (8) :: M_Fit      
                Integer, dimension (3) :: D_Fit     
                Integer, dimension (5) :: I_Fit     
                Integer, dimension (2) :: H_Fit   

                integer:: i,success_test,fOutputNum,j,nj,FailedR(5)
                real*8::testArr(52),errArr(52),maxErr,maxErrp,err_tol,diff_comp(52),Rmin,Rmax

                if (present(fileOutputNumber))then
                        fOutputNum = fileOutputNumber
                else
                        fOutputNum = 0
                end if
                
                

                maxErr=0d0
                maxErrp=0d0
                err_tol = 1d-5
                diff_comp =0d0
                FailedR = 0

                Call READ_DATASET(data_filename,systName,ntest,XDIM,dataset_,coord_format,M_Fit ,D_Fit,I_Fit,H_Fit)

                Allocate(XDim_coord(XDIM))

                Rmin = MINVAL(dataset_(:,1))
                Rmax = MAXVAL(dataset_(:,1))
                success_test = 0

                if (present(ntestMAX))then
                        nMAX = ntestMAX
                else
                        nMAX = ntest
                end if

                if (present(ntestMIN))then
                        nMIN = ntestMIN
                else
                        nMIN = 1
                end if

                do i=nMIN,nMAX

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

                        !write(*,*)testArr(21),dataset_(i,20)

                        if(diff_comp(1)>maxErr)then
                            maxErr = diff_comp(1)
                        end if

                        if(100d0*(diff_comp(1)/E_fit)>maxErrp)then
                            maxErrp = 100d0*(diff_comp(1)/E_fit)
                        end if

                        if(diff_comp(1)>err_tol )then
                           
                               

                               if ( XDim_coord(1) < Rmin+(Rmax-Rmin)*0.2 )then
                               FailedR(1)=FailedR(1)+1
                               else if(  XDim_coord(1) >= Rmin+(Rmax-Rmin)*0.2 .and.&
                                         XDim_coord(1) <  Rmin+(Rmax-Rmin)*0.4 )then
                                FailedR(2)=FailedR(2)+1
                                else if(  XDim_coord(1) >= Rmin+(Rmax-Rmin)*0.4 .and.&
                                         XDim_coord(1) <  Rmin+(Rmax-Rmin)*0.6 )then
                                FailedR(3)=FailedR(3)+1
                                else if(  XDim_coord(1) >= Rmin+(Rmax-Rmin)*0.6 .and.&
                                         XDim_coord(1) <  Rmin+(Rmax-Rmin)*0.8 )then
                                FailedR(4)=FailedR(4)+1
                                else if(  XDim_coord(1) >= Rmin+(Rmax-Rmin)*0.8 .and.&
                                         XDim_coord(1) <=  Rmax )then
                                FailedR(5)=FailedR(5)+1
                               end if 
                                
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
                                        if (diff_comp(j) >0d0)then! err_tol )then

                                                if (j>5 .and. j<21)then
                                                        write(str,'(A,I1)')"M",j-5
                                                        nj = j + 6
                                                elseif(j>20 .and. j<31)then
                                                        write(str,'(A,I1)')"D",j-20+5
                                                        nj = j -1
                                                elseif(j>30 .and. j<43)then
                                                        write(str,'(A,I1)')"I",j-30+3  
                                                        nj = j -8
                                                elseif(j>42)then
                                                        write(str,'(A,I1)')"H",j-42+5      
                                                        nj = j -15        
                                                end if
                                                
                                                
                                                write(fOutputNum,'(A,F15.9,F15.9,E15.3 )') str&
                                                        ,dataset_(i,nj),testArr(j),diff_comp(j)
                                        end if 
                                end do 
                                write(fOutputNum,*) "*****************************************"
                                

                                end if 
                        else
                                success_test =success_test+1
                        end if
                   
                        
                     end if 



                enddo

                Call Print_Test_Results(systName,maxErr,nMAX-nMIN+1,success_test,fOutputNum)
                write(*,*) "FailedR by R: ",Rmin,Rmin+(Rmax-Rmin)*0.05,FailedR(1)
                write(*,*) "FailedR by R: ",Rmin+(Rmax-Rmin)*0.05,Rmin+(Rmax-Rmin)*0.1,FailedR(2)
                write(*,*) "FailedR by R: ",Rmin+(Rmax-Rmin)*0.1,Rmin+(Rmax-Rmin)*0.5,FailedR(3)
                write(*,*) "FailedR by R: ",Rmin+(Rmax-Rmin)*0.5,Rmin+(Rmax-Rmin)*0.8,FailedR(4)
                write(*,*) "FailedR by R: ",Rmin+(Rmax-Rmin)*0.8,Rmax,FailedR(5)
 
                write(*,*) "maxErrp : " ,maxErrp
        end Subroutine Test_Dataset


        Subroutine Test_Component(coeff_filename, data_filename,Comp_filename,&
                                        fileOutputNumber,ntestMIN,ntestMAX)

                use Tensors_constant
                IMPLICIT NONE
                Character(len = *),Intent(IN) :: coeff_filename, data_filename,Comp_filename
                integer,optional:: fileOutputNumber,ntestMIN,ntestMAX

                real*8 :: E_fortran,E_fit
                real*8 ,dimension(6):: GeneralCoordenates,coordenates,zero_coordinates
                real*8, allocatable :: XDim_coord(:),dataset_(:,:),Check_val(:)
                Character(len = 20) :: coord_format,systName
                Character(len = 2) :: str
                INTEGER :: XDIM,ntest,nMAX,nMIN
                real*8 :: start, finish,pii
                real*8 :: ind, R,cos_b1,cos_b2,phi, E0
                Integer, dimension (8) :: M_Fit      
                Integer, dimension (3) :: D_Fit     
                Integer, dimension (5) :: I_Fit     
                Integer, dimension (2) :: H_Fit   

                real*8 , dimension(11):: cal_coord
                real*8 , dimension(1195):: coeff_arr
                Real*8  ::   Zero 

                real*8 , dimension(3):: Ar 
                real*8 , dimension(3):: Br
                real*8 , dimension(9):: C

                integer:: i,success_test,fOutputNum,j,nj
                real*8::testArr(52),errArr(52),maxErr,err_tol,diff_comp(52),Cal_Val,Const

                real*8, parameter ::  C1=627.5095d0,C2=0.529177249d0,C3=349.757d0

                real*8 , dimension(64) :: A_Mult,B_Mult !q, mz, Qz, Oz, Phiz, M5z, M6z, M7z 
                real*8 , dimension(57) :: A_Pol,B_Pol
                real*8 , dimension(40) :: A_HPol,B_HPol
                real*8 , dimension(873) :: Disp_AB

                Const = C3*(C1*C2**2)

                write(*,*)Const


                if (present(fileOutputNumber))then
                        fOutputNum = fileOutputNumber
                else
                        fOutputNum = 0
                end if
                
                

                maxErr=0d0
                err_tol = 1d-6
                diff_comp =0d0

                Call READ_DATASET(data_filename,systName,ntest,XDIM,dataset_,coord_format,M_Fit ,D_Fit,I_Fit,H_Fit)

                Allocate(XDim_coord(XDIM),Check_val(ntest))

                Call Read_Check_val(Comp_filename,Check_val,ntest)



                success_test = 0

                if (present(ntestMAX))then
                        nMAX = ntestMAX
                else
                        nMAX = ntest
                end if

                if (present(ntestMIN))then
                        nMIN = ntestMIN
                else
                        nMIN = 1
                end if

                do i=nMIN,nMAX

                     if (XDIM==3)Then
                        
                     else
                        XDim_coord = dataset_(i,1:XDIM)
                       
                        call init_Tensors() ! Initializing in zero the new vectors

                        CALL Prep_Param(coeff_filename,coeff_arr,M_Fit ,D_Fit,I_Fit,H_Fit,Zero)


                        A_Mult=coeff_arr(1:64)
                        B_Mult=coeff_arr(65:128)
                        A_Pol=coeff_arr(129:185)
                        B_Pol=coeff_arr(186:242)
                        A_HPol=coeff_arr(243:282)
                        B_HPol=coeff_arr(283:322)
                        Disp_AB=coeff_arr(323:1195)

                        Call General_Coordinates_Format(XDIM, XDim_coord, GeneralCoordenates)
                        Call Coordinate_Transformation(GeneralCoordenates,coord_format,coordenates)
                        Call Generate_Coordenates(coordenates,cal_coord,Ar,Br,C)

                        Call Approx_2_Sph2(cal_coord,Ar,Br,C , A_Mult,B_Mult ,Cal_Val)

                        diff_comp(1) = DABS(Const*Cal_Val - Const*Check_val(i))

                      

                        if(MAXVAL(diff_comp)>maxErr)then
                            maxErr = MAXVAL(diff_comp)
                        end if

                        if(diff_comp(1)>err_tol )then
                                
                                if ( fOutputNum < 0 )then
                                        write(*,*) 
                                        write(*,*) "Row : ",i
                                        write(*,*) "Coordinates : ",XDim_coord
                                        write(*,*) "Check_val: ",Const*Cal_Val,Const*Check_val(i),diff_comp(1)
                                        write(*,*) 


                                elseif ( fOutputNum > 0 )then

                                write(fOutputNum,*) 
                                write(fOutputNum,*) "Row : ",i
                                write(fOutputNum,*) "Coordinates : ",XDim_coord
                                write(fOutputNum,*) "Check_val: ",Const*Cal_Val,Const*Check_val(i),diff_comp(1)
                               


                                write(fOutputNum,*) "*****************************************"
                                

                                end if 
                        else
                                success_test =success_test+1
                        end if
                   
                        
                     end if 



                enddo

                Call Print_Test_Results(systName,maxErr,nMAX-nMIN+1,success_test,fOutputNum)

        end Subroutine Test_Component

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
                write(*,  '(A, E15.4)') "                                Error Max: ",maxErr                                 
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

        SUBROUTINE Testing_TTensors(filename,ntestMAX,level_init,level_final)
                use Tensors_constant 
                IMPLICIT NONE
                Character(len = *),Intent(IN) :: filename
                integer,optional:: ntestMAX,level_init,level_final

                real*8 :: E_fortran,E_fit
                real*8 ,dimension(6):: GeneralCoordenates,coordinates,zero_coordinates
                real*8, allocatable :: result(:,:,:,:),dataset_(:,:),res(:,:)
                Character(len = 20) :: coord_format,systName
                Character(len = 2) :: str
                INTEGER :: XDIM,ntest,nMAX
                real*8 :: start, finish,pii
                real*8 , dimension(3):: Ar 
                real*8 , dimension(3):: Br
                real*8 , dimension(9):: C
                Integer la_,level,lb_
                real*8 , dimension(11):: cal_coord  
                character(len=20)::title
                integer:: i,success_test,passed,fOutputNum,j,k,l_init,l_final
                integer,parameter::lmax=20
                real*8::testArr(52),errArr(52),maxErr,err_tol,diff_comp(52),err
                real*8::XDim_coord(6),res_partial(lmax,lmax,(2*lmax+1)*(2*lmax+1))

      

                Call READ_DATASET_TTensor(filename,systName,ntest,DATASET_,coord_format)

                if (present(ntestMAX))then
                        nMAX = ntestMAX
                else
                        nMAX = ntest
                end if

                if (present(level_init) .and. present(level_final))then
                        l_init = level_init 
                        l_final= level_final
                else
                        l_init = 1
                        l_final= 8
                end if


                Allocate(result(ntest,lmax,lmax,(2*lmax+1)*(2*lmax+1)),&
                         res(ntest,(2*lmax+1)*(2*lmax+1)))
            

                result = 0d0
                res = 0d0
                maxErr =0d0

                do level=1,level_final
                        do la_=0,level-1
                         lb_ = level-la_-1
                                Call READ_TTensor_Values(ntest,la_,lb_,res )
                                
                                do k=1,nMAX
                                        result(k,la_+1,lb_+1,1:(2*la_+1)*(2*lb_+1)) = res(k,1:(2*la_+1)*(2*lb_+1))
                                end do 
                        end do
                end do

                success_test = 0



                 do i=1,nMAX
                        
                        call init_Tensors() ! Initializing in zero the new vectors

                        XDim_coord(1) = 1
                        XDim_coord(2:6) = dataset_(i,1:5)
                 

                        Call General_Coordinates_Format(6, XDim_coord, GeneralCoordenates)
                        Call Coordinate_Transformation(GeneralCoordenates,coord_format,coordinates)
                        Call Generate_Coordenates(coordinates,cal_coord,Ar,Br,C)

                        res_partial = 0d0
                        res_partial = result(i,:,:,:) 
                       
                        Call Check_Tensors(Ar,Br,C,res_partial,passed,err,l_init,l_final) 

                        if (err > maxErr)then
                                maxErr = err
                        end if 
                        success_test = success_test + passed

                 end do

                write(title,'(A,I1,A,I1)')"T-Tensors Check: ",level_init,"-",level_final

                Call Print_Test_Results(title,maxErr,nMAX,success_test,fOutputNum)

        END SUBROUTINE Testing_TTensors 

        SUBROUTINE Check_Tensors(Ar,Br,C,result,passed,errMax,linit,lfinal)
                use Tensors_constant 
                IMPLICIT NONE
                
                
                real*8 , dimension(3),Intent(IN):: Ar 
                real*8 , dimension(3),Intent(IN):: Br
                real*8 , dimension(9),Intent(IN):: C
                integer,parameter::lmax=20
                integer,Intent(IN)::linit,lfinal
                integer , intent(out)::passed
                real*8  , intent(out)::errMax
                real*8 , Intent(IN):: result(lmax,lmax,(2*lmax+1)*(2*lmax+1))
                real*8 :: res,check(lmax,lmax,(2*lmax+1)*(2*lmax+1)),err_tol,&
                                diff(lmax,lmax,(2*lmax+1)*(2*lmax+1))
                integer::level,la,lb,ka,kb,mla,mlb,mka,mkb,ind
                character(len=1)::cpa,cpb
     
                 err_tol =1d-6
                 errMax = 0d0
                 check = 0d0
                 passed = 1


                diff=0d0



                do level=linit,lfinal
                        do la=0,level-1
                                lb = level-la-1
                                mla = 2*la+1
                                mlb = 2*lb+1

                                do mka=0,mla-1
                                        do mkb=0,mlb-1
                                          
                                         ! call init_Tensors()

                                          call Get_Comp(mka+1,cpa)
                                          call Get_Comp(mkb+1,cpb)
                                          ka =  (mka+1)/2
                                          kb =  (mkb+1)/2

                                          ind = mka*mlb+mkb+1
                                        
                                          Call T_lk(Ar,Br,C,la,ka,cpa,lb,kb,cpb,check(la+1,lb+1,ind))
                                          diff(la+1,lb+1,ind) = DABS(check(la+1,lb+1,ind)- result(la+1,lb+1,ind))

                                          if (diff(la+1,lb+1,ind) > errMax)then
                                                errMax = diff(la+1,lb+1,ind)
                                          end if 
                                          if (diff(la+1,lb+1,ind) > err_tol)then
                                                passed = 0 
                                                write(*,*)"Problems with ",la,ka,cpa,lb,kb,cpb,check(la+1,lb+1,ind)&
                                                , result(la+1,lb+1,ind),diff(la+1,lb+1,ind)
                                             
                                          end if 
                                        end do
                                end do


                        end do
                end do




                 

        END SUBROUTINE Check_Tensors
        
        Subroutine Testing_GP(data_filename,GP_filename,fileOutputNumber,ntestMIN,ntestMAX)
                IMPLICIT NONE
                Character(len = *),Intent(IN) :: data_filename,GP_filename
                integer,optional:: fileOutputNumber,ntestMIN,ntestMAX

                real*8 :: E_fortran,E_fit
                real*8 ,dimension(6):: GeneralCoordenates,coordenates,zero_coordinates
                real*8, allocatable :: XDim_coord(:),dataset_(:,:)
                Character(len = 20) :: coord_format,systName
                Character(len = 2) :: str
                INTEGER :: XDIM,ntest,nMAX,nMIN
                real*8 :: start, finish,pii
                real*8 :: ind, R,cos_b1,cos_b2,phi, E0
                real*8 , dimension(11):: cal_coord  
                real*8 , allocatable:: Ar(:,:),Br(:,:),C(:,:)
                real*8 , allocatable:: Ar_f(:,:),Br_f(:,:),C_f(:,:)
                integer:: i,success_test,fOutputNum,j,nj
                real*8::testArr(52),errArr(52),maxErr,err_tol,diff_comp(52)
                real*8 , dimension(3):: Ar_
                real*8 , dimension(3):: Br_
                real*8 , dimension(9):: C_
                Integer, dimension (8) :: M_Fit      
                Integer, dimension (3) :: D_Fit     
                Integer, dimension (5) :: I_Fit     
                Integer, dimension (2) :: H_Fit   


                if (present(fileOutputNumber))then
                        fOutputNum = fileOutputNumber
                else
                        fOutputNum = 0
                end if
                
                if (present(ntestMAX))then
                        nMAX = ntestMAX
                else
                        nMAX = ntest
                end if

                if (present(ntestMIN))then
                        nMIN = ntestMIN
                else
                        nMIN = 1
                end if

                maxErr=0d0
                err_tol = 1d-6
                diff_comp =0d0

                Call READ_DATASET(data_filename,systName,ntest,XDIM,dataset_,coord_format,M_Fit ,D_Fit,I_Fit,H_Fit)
                XDIM = 6
                Allocate(XDim_coord(XDIM))
                Allocate(Ar(nMAX,3),Br(nMAX,3),C(nMAX,9))
                Allocate(Ar_f(nMAX,3),Br_f(nMAX,3),C_f(nMAX,9))

                Call READ_GP(GP_filename,nMAX,Ar,Br,C)

            

                success_test = 0

                

                do i=nMIN,nMAX

                     if (XDIM==3)Then
                        
                     else
                        XDim_coord = dataset_(i,1:XDIM)
                        E_fit = dataset_(i,7)

                        Call General_Coordinates_Format(XDIM, XDim_coord, GeneralCoordenates)
                        Call Coordinate_Transformation(GeneralCoordenates,coord_format,coordenates)
                        Call Generate_Coordenates(coordenates,cal_coord,Ar_,Br_,C_)

                        diff_comp(1:3) = DABS(Ar(i,:) - Ar_)
                        diff_comp(4:6) = DABS(Br(i,:) - Br_)
                        diff_comp(7) = DABS(C(i,1) - C_(1))
                        diff_comp(8) = DABS(C(i,2) - C_(3))
                        diff_comp(9) = DABS(C(i,3) - C_(2))
                        diff_comp(10) = DABS(C(i,4) - C_(7))
                        diff_comp(11) = DABS(C(i,5) - C_(9))
                        diff_comp(12) = DABS(C(i,6) - C_(8))   
                        diff_comp(13) = DABS(C(i,7) - C_(4)) 
                        diff_comp(14) = DABS(C(i,8) - C_(6)) 
                        diff_comp(15) = DABS(C(i,9) - C_(5)) 


                        if(MAXVAL(diff_comp)>maxErr)then
                            maxErr = MAXVAL(diff_comp)
                        end if

                        if(MAXVAL(diff_comp)>err_tol )then
                                
                                if ( fOutputNum < 0 )then
                                        write(*,*) 
                                        write(*,*) "Row : ",i
                                        write(*,*) "Coordinates : ",XDim_coord
                                        write(*,*) "Arz: ",Ar(i,1),Ar_(1),diff_comp(1)
                                        write(*,*) "Ary: ",Ar(i,2),Ar_(2),diff_comp(2)
                                        write(*,*) "Arx: ",Ar(i,3),Ar_(3),diff_comp(3)
                                        write(*,*) "Brz: ",Br(i,1),Br_(1),diff_comp(4)
                                        write(*,*) "Bry: ",Br(i,2),Br_(2),diff_comp(5)
                                        write(*,*) "Brx: ",Br(i,3),Br_(3),diff_comp(6)
                                        write(*,*) "Czz: ",C(i,1),C_(1),diff_comp(6+1)
                                        write(*,*) "Cyz: ",C(i,2),C_(2),diff_comp(6+2)
                                        write(*,*) "Cxz: ",C(i,3),C_(3),diff_comp(6+3)
                                        write(*,*) "Czy: ",C(i,4),C_(4),diff_comp(6+4)
                                        write(*,*) "Cyy: ",C(i,5),C_(5),diff_comp(6+5)
                                        write(*,*) "Cxy: ",C(i,6),C_(6),diff_comp(6+6)
                                        write(*,*) "Czx: ",C(i,7),C_(7),diff_comp(6+7)
                                        write(*,*) "Cyx: ",C(i,8),C_(8),diff_comp(6+8)
                                        write(*,*) "Cxx: ",C(i,9),C_(9),diff_comp(6+9)
                                        write(*,*) 


                                elseif ( fOutputNum > 0 )then

                                write(fOutputNum,*) 
                                        write(fOutputNum,*) "Row : ",i
                                        write(fOutputNum,*) "Coordinates : ",XDim_coord
                                        write(fOutputNum,*) "Arz: ",Ar(i,1),Ar_(1),diff_comp(1)
                                        write(fOutputNum,*) "Ary: ",Ar(i,2),Ar_(2),diff_comp(2)
                                        write(fOutputNum,*) "Arx: ",Ar(i,3),Ar_(3),diff_comp(3)
                                        write(fOutputNum,*) "Brz: ",Br(i,1),Br_(1),diff_comp(4)
                                        write(fOutputNum,*) "Bry: ",Br(i,2),Br_(2),diff_comp(5)
                                        write(fOutputNum,*) "Brx: ",Br(i,3),Br_(3),diff_comp(6)
                                        write(fOutputNum,*) "Czz: ",C(i,1),C_(1),diff_comp(6+1)
                                        write(fOutputNum,*) "Cyz: ",C(i,2),C_(2),diff_comp(6+2)
                                        write(fOutputNum,*) "Cxz: ",C(i,3),C_(3),diff_comp(6+3)
                                        write(fOutputNum,*) "Czy: ",C(i,4),C_(4),diff_comp(6+4)
                                        write(fOutputNum,*) "Cyy: ",C(i,5),C_(5),diff_comp(6+5)
                                        write(fOutputNum,*) "Cxy: ",C(i,6),C_(6),diff_comp(6+6)
                                        write(fOutputNum,*) "Czx: ",C(i,7),C_(7),diff_comp(6+7)
                                        write(fOutputNum,*) "Cyx: ",C(i,8),C_(8),diff_comp(6+8)
                                        write(fOutputNum,*) "Cxx: ",C(i,9),C_(9),diff_comp(6+9)
                                        write(fOutputNum,*)  

                                write(fOutputNum,*) "*****************************************"
                                

                                end if 
                        else
                                success_test =success_test+1
                        end if
                   
                        
                      end if 



                enddo

                Call Print_Test_Results(systName,maxErr,nMAX-nMIN+1,success_test,fOutputNum)

        end Subroutine Testing_GP
     
end module Testing



PROGRAM main_subroutine

        use Testing
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


        !call Test_All(fileOutNumber)
        call Test_Dataset("./files/test/coefficients/coefficients_003.txt", "./files/test/datasets/datatest_003_001.txt"&
                       ,fileOutNumber,1,10000)
        ! Call Test_Component("./files/test/coefficients/coefficients_003.txt", &
        !                     "./files/test/datasets/datatest_003_001.txt",&
        !                     "./files/test/M2(CF+_H2+).txt",&
        !                    fileOutNumber,1,10000)
        ! call Testing_GP( "./files/test/datasets/datatest_003_001.txt","./files/test/GPTable(CF+_H2+).txt"&
        !                ,fileOutNumber,1,10000)
        !call Testing_TTensors('./files/test/T_Tensors/datatest.txt',10,level_init,level_final)
     
        close(fileOutNumber)
END PROGRAM main_subroutine






