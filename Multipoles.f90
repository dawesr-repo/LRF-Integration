
!  Fortran subroutine for 3.1.1 matlab version Optimization in running time
! last test: 26/12/2022


PROGRAM main_subroutine

        IMPLICIT NONE
        
        real*8 :: E2,E1,E0
        real*8 ,dimension(6):: GeneralCoordenates,coordinates,zero_coordinates
        real*8 :: FourD_coord(4),TwoD_coord(2)
        Character(len = 20) :: coord_format = "Euler_ZYZ" !for Xdim =3, use coord_format ="Spherical"  
        INTEGER :: XDim=4,i
        real*8 :: start, finish
        real*8::testErr(52)
                FourD_coord(1) = 10.271963668339600d0 !R
                FourD_coord(2) = 30d0  !b1
                FourD_coord(3) = 20d0  !b2
                FourD_coord(4) = 120d0  !Phi

                

                call cpu_time(start)

             

                Call General_Coordinates_Format(XDim, FourD_coord, GeneralCoordenates)
                Call Coordinate_Transformation(GeneralCoordenates,coord_format,coordinates)
        
        
                    CALL Long_Range_Potential(coordinates,E1,&
                "./files/test/coefficients/coefficients_003.txt")
    

            write(*,*)E1
              

        ! DON't use Coordinate_Transformation function to obtain the zero potential, 
        ! instead pass an array with all of them zero
        ! zero_coordinates(1) = 0d0  !R
        ! zero_coordinates(2) = 0d0  !b1
        ! zero_coordinates(3) = 0d0  !b2
        ! zero_coordinates(4) = 0d0  !phi
        ! zero_coordinates(5) = 0d0  !c1
        ! zero_coordinates(6) = 0d0  !c2


        ! CALL Long_Range_Potential(zero_coordinates,E0,"./files/test/coefficients/coefficients_003.txt"&
        ! ,0,testErr)
    
        ! write(*,*)E0




    
END PROGRAM main_subroutine








