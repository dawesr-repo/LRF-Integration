PROGRAM min_example

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


 call evaluateLR(X1,XDIM,E1,'./files/test/coefficients/C1(1)_C1(1)_Coeff.txt')
  
END PROGRAM min_example










