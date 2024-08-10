
!  Fortran subroutine for 3.1.1 matlab version Optimization in running time
! last test: 26/12/2022


! PROGRAM min_example
 
!  IMPLICIT NONE
!  real*8 :: E1
!  INTEGER :: XDIM=4
!  real*8 :: X1(4)
 
!  X1(1) = 10.271963668339600d0 !R
!  X1(2) = 30d0  !b1
!  X1(3) = 20d0  !b2
!  X1(4) = 120d0  !Phi

!  call evaluateLR(X1,XDIM,E1,'./files/coefficients.txt')

    

!  write(*,*)"Energy: ", E1
 
! END PROGRAM min_example

PROGRAM min_example
 use FitConstants
 use Search
 IMPLICIT NONE
 real*8 :: E1
 INTEGER :: XDIM=4
 real*8 :: X1(4)
 real*8::Energy
 real*8, allocatable:: arrp(:),arrh(:),arrd(:)
 character(len=5),allocatable::arr(:)
 integer ::i,j,k,l,index1,index2,counter

 X1(1) = 10.271963668339600d0 !R
 X1(2) = 30d0  !b1
 X1(3) = 20d0  !b2
 X1(4) = 120d0  !Phi

 
 call evaluateLR(X1,XDIM,E1,'./files/test/coefficients/coefficients_002.txt')
 



!  counter = 0

!  do i=1,3
!     do j=1,3
!     do k=1,3
!                 call Get_Hyper_Index(1,1,1,i,j,k,index1)
!                 Call GetIndex_mmm(i,j,k,index2)
!                  write(*,*)i,j,k,index1,index2,index1-index2
!                 if(ABS(index1-index2)>0)then
!                     if (counter==0)then
                     
! write(*,*)"          i           ","j           ","idnew      ","idold    ","diff"   
! write(*,*)"________________________________________________________________________________________"

!                     end if
!                      write(*,*)i,j,k,index1,index2,index1-index2
!                      counter=counter +1
!                 end if
!     end do 
!     end do 
!  end do 

!  write(*,*)"Index checking ", 1,1,1 , "Counting: ", counter

!   counter = 0

!  do i=1,3
!     do j=1,5
!                 Call Get_Ind_Index(1,2,i,j,index1)
!                 Call GetIndex_mQd(i,j,index2)
!                 if(ABS(index1-index2)>0)then
!                     if (counter==0)then
                     
! write(*,*)"          i           ","j           ","idnew      ","idold    ","diff"   
! write(*,*)"________________________________________________________________________________________"

!                     end if
!                      write(*,*)i,j,index1,index2,index1-index2
!                      counter=counter +1
!                 end if

!     end do 
!  end do 

!  write(*,*)"Index checking ", 1,2 , "Counting: ", counter

!    counter = 0

!  do i=1,5
!     do j=1,5
!                 Call Get_Ind_Index(2,2,i,j,index1)
!                 Call GetIndex_QdQd(i,j,index2)
!                 if(ABS(index1-index2)>0)then
!                     if (counter==0)then
                     
! write(*,*)"          i           ","j           ","idnew      ","idold    ","diff"   
! write(*,*)"________________________________________________________________________________________"

!                     end if
!                      write(*,*)i,j,index1,index2,index1-index2
!                      counter=counter +1
!                 end if

!     end do 
!  end do 

!  write(*,*)"Index checking ", 2,2 , "Counting: ", counter

!    counter = 0

!  do i=1,3
!     do j=1,7
!                 Call Get_Ind_Index(1,3,i,j,index1)
!                 Call GetIndex_mO(i,j,index2)
!                 if(ABS(index1-index2)>0)then
!                     if (counter==0)then
                     
! write(*,*)"          i           ","j           ","idnew      ","idold    ","diff"   
! write(*,*)"________________________________________________________________________________________"

!                     end if
!                      write(*,*)i,j,index1,index2,index1-index2
!                      counter=counter +1
!                 end if

!     end do 
!  end do 

!  write(*,*)"Index checking ", 1,3 , "Counting: ", counter

END PROGRAM min_example










