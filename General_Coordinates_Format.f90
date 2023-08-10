
    !********************************************************
SUBROUTINE General_Coordinates_Format(Dim, old_coordinates, general_coodinates)
            IMPLICIT NONE
    
            INTEGER,INTENT(IN) ::  Dim
            real*8 , dimension(6), INTENT(INOUT) ::general_coodinates
            real*8 , dimension (Dim),INTENT(IN) :: old_coordinates
       

            if (Dim==2) then
                    general_coodinates(1) = old_coordinates(1)  !R
                    general_coodinates(2) = old_coordinates(2)  !b1
                    general_coodinates(3) = 0d0  !b2
                    general_coodinates(4) = 0d0  !phi
                    general_coodinates(5) = 0d0  !c1
                    general_coodinates(6) = 0d0  !c2

                elseif (Dim==3) then
                    general_coodinates(1) = old_coordinates(1)  !R
                    general_coodinates(2) = old_coordinates(2)  !b1
                    general_coodinates(3) = 0d0  !b2
                    general_coodinates(4) = 0d0  !phi
                    general_coodinates(5) = old_coordinates(3)  !c1
                    general_coodinates(6) = 0d0  !c2

                elseif (Dim==4) then   
                    general_coodinates(1) = old_coordinates(1)  !R
                    general_coodinates(2) = old_coordinates(2)  !b1
                    general_coodinates(3) = old_coordinates(3)  !b2
                    general_coodinates(4) = old_coordinates(4)  !phi
                    general_coodinates(5) = 0d0  !c1
                    general_coodinates(6) = 0d0  !c2 

                elseif (Dim==5) then   
                        general_coodinates(1) = old_coordinates(1)  !R
                        general_coodinates(2) = old_coordinates(2)  !b1
                        general_coodinates(3) = old_coordinates(3)  !b2
                        general_coodinates(4) = old_coordinates(4)  !phi
                        general_coodinates(5) = old_coordinates(5)  !c1
                        general_coodinates(6) = 0d0  !c2     
                elseif (Dim==6) then   
                        general_coodinates(1) = old_coordinates(1)  !R
                        general_coodinates(2) = old_coordinates(2)  !b1
                        general_coodinates(3) = old_coordinates(3)  !b2
                        general_coodinates(4) = old_coordinates(4)  !phi
                        general_coodinates(5) = old_coordinates(5)  !c1
                        general_coodinates(6) = old_coordinates(6)  !c2       

            end if
  
            RETURN
END SUBROUTINE General_Coordinates_Format

