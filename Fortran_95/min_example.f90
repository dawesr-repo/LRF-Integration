PROGRAM min_example

 IMPLICIT NONE

 INTEGER :: xdim=4
 real*8 :: coordinates(4)
 real*8 :: Energy
 coordinates(1) = 10.271963668339600d0 !R
 coordinates(2) = 30d0  !b1
 coordinates(3) = 20d0  !b2
 coordinates(4) = 120d0  !Phi


 call Evaluate_LRF(   Energy,&
         coordinates,&
         "Euler_ZYZ",&
         xdim,&
         '../testing_datafiles/coefficients/C1(1)_C1(1)_Coeff.txt')


END PROGRAM min_example










