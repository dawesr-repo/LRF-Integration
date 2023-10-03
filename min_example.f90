
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
 integer ::i,j,k,l,index1,index2

 X1(1) = 10.271963668339600d0 !R
 X1(2) = 30d0  !b1
 X1(3) = 20d0  !b2
 X1(4) = 120d0  !Phi


 call evaluateLR(X1,XDIM,E1,'./files/test/coefficients/coefficients_003.txt')
 
 call Initialize_Search()
 do i=1,3
    do j=1,5
        do k=1,3
            do l=1,3
                Call Get_Disp_Index(1,2,1,1,i,j,k,l,index1)
                Call GetIndex_mQd_mm(i,j,k,l,index2)
                write(*,*)i,j,k,l,index1,index2,index1-index2
            end do 
        end do 
    end do 
 end do 
 write(*,*)"Energy: ", E1


 
END PROGRAM min_example






