
    SUBROUTINE Coordinate_Transformation(coordenates,coord_format,new_coordinates)
        IMPLICIT NONE
        real*8 ,dimension(6), INTENT(IN):: coordenates
        
        Character(len = 20), INTENT(IN) :: coord_format
        real*8 , dimension(6), INTENT(INOUT) :: new_coordinates
       
        new_coordinates=coordenates;
        
        if (coord_format == "Euler_ZYZ") then
            new_coordinates(5) = coordenates(5) - 90
            new_coordinates(6) = coordenates(6) - 90
        end if 
        if (coord_format == "Spherical") then
            new_coordinates(5) = 90 - coordenates(5) 
            new_coordinates(6) = 90 -coordenates(6) 
        end if 
        
    End   SUBROUTINE Coordinate_Transformation