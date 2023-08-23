
!  Fortran subroutine for 3.1.1 matlab version Optimization in running time
! last test: 26/12/2022


PROGRAM min_example

 IMPLICIT NONE
 real*8 :: E1,E2
 INTEGER :: XDIM=4
 real*8 :: X1(4)
 
 X1(1) = 10.271963668339600d0 !R
 X1(2) = 30d0  !b1
 X1(3) = 20d0  !b2
 X1(4) = 120d0  !Phi

 call evaluateLR_1(X1,XDIM,E1,'./files/coefficientsAp.txt')
 call evaluateLR_2(X1,XDIM,E2,'./files/coefficients_Adp.txt')

 write(*,*)"Energy: ", E1,E2
 
END PROGRAM min_example








