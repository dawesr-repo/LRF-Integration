subroutine coordinate_transformation(new_coordinates,coordinates,coord_format)
    
    implicit none
    real (kind=8) ,dimension(6), intent(in):: coordinates
    character(*), intent(in) :: coord_format
    real (kind=8) , dimension(6), intent(out) :: new_coordinates

    new_coordinates=coordinates;
    
    if (coord_format == "Euler_ZYZ") then
        new_coordinates(5) = coordinates(5) - 90
        new_coordinates(6) = coordinates(6) - 90
    end if
    if (coord_format == "Spherical") then
        new_coordinates(5) = 90 - coordinates(5)
        new_coordinates(6) = 90 -coordinates(6)
    end if

end  subroutine coordinate_transformation

!********************************************************
subroutine general_coordinates_format(general_coodinates,dim, old_coordinates)
    
    implicit none
    integer (kind=4),intent(in) ::  dim
    real (kind=8) , dimension (dim),intent(in) :: old_coordinates
    real (kind=8) , dimension(6), intent(out) ::general_coodinates

    general_coodinates = 0d0
    if (dim==2) then
        general_coodinates(1) = old_coordinates(1)  !R
        general_coodinates(2) = old_coordinates(2)  !b1
        general_coodinates(3) = 0d0  !b2
        general_coodinates(4) = 0d0  !phi
        general_coodinates(5) = 0d0  !c1
        general_coodinates(6) = 0d0  !c2

    elseif (dim==3) then
        general_coodinates(1) = old_coordinates(1)  !R
        general_coodinates(2) = old_coordinates(2)  !b1
        general_coodinates(3) = 0d0  !b2
        general_coodinates(4) = 0d0  !phi
        general_coodinates(5) = old_coordinates(3)  !c1
        general_coodinates(6) = 0d0  !c2

    elseif (dim==4) then
        general_coodinates(1) = old_coordinates(1)  !R
        general_coodinates(2) = old_coordinates(2)  !b1
        general_coodinates(3) = old_coordinates(3)  !b2
        general_coodinates(4) = old_coordinates(4)  !phi
        general_coodinates(5) = 0d0  !c1
        general_coodinates(6) = 0d0  !c2

    elseif (dim==5) then
        general_coodinates(1) = old_coordinates(1)  !R
        general_coodinates(2) = old_coordinates(2)  !b1
        general_coodinates(3) = old_coordinates(3)  !b2
        general_coodinates(4) = old_coordinates(4)  !phi
        general_coodinates(5) = old_coordinates(5)  !c1
        general_coodinates(6) = 0d0  !c2
    elseif (dim==6) then
        general_coodinates(1) = old_coordinates(1)  !R
        general_coodinates(2) = old_coordinates(2)  !b1
        general_coodinates(3) = old_coordinates(3)  !b2
        general_coodinates(4) = old_coordinates(4)  !phi
        general_coodinates(5) = old_coordinates(5)  !c1
        general_coodinates(6) = old_coordinates(6)  !c2

    end if

    RETURN
end subroutine general_coordinates_format

subroutine user_coordinates_to_general_coordinates(general_coordinates_ZXZ,xdim,coord_format,user_coordinates)

    implicit none
    integer (kind=4),intent(in) ::  xdim
    real (kind=8), dimension (xdim),intent(in) :: user_coordinates
    character(*), intent(in) :: coord_format
    real (kind=8),dimension (6),intent(out):: general_coordinates_ZXZ
    real (kind=8) ,dimension(6):: general_coordenates


    call general_coordinates_format(general_coordenates,xdim, user_coordinates)
    call coordinate_transformation(general_coordinates_ZXZ,general_coordenates,coord_format)

end subroutine user_coordinates_to_general_coordinates
! Arg 1 [coordinates] : a coordinate vector [ R , b1, b2, phi] *the angles should be in degrees
! Arg 2 [coeff_Address] address of the file which contains the longe range expansion coefficients
! Arg 3 [total_energy]   Total Energy calculated


! version 4.0


subroutine evaluate_LRF(total_energy,xdim,coordinates,coord_format,filename)

    use Fitting_Constant_v2, only: get_coeff_zero, get_coeff_index
    use Geometry_Constant_v2, only: get_total_interaction_energy

    implicit none
    real (kind=8), intent(out) :: total_energy
    integer (kind=4), intent(in) :: xdim
    real (kind=8) ,dimension(:), intent(in):: coordinates(xdim)
    character(*), intent(in) :: coord_format
    character(*), intent(in) ::  filename
    
    real (kind=8) ,dimension(6):: general_coordinates_ZXZ
    integer (kind=4) :: i,coeff_index
    real (kind=8) :: x1

    coeff_index = get_coeff_index(filename)

    x1 = 0d0
    do i=1,xdim
        x1=x1+dabs(coordinates(i))
    enddo

    ! if the user coordinate array is set to zero, then
    ! the function will return the dissociation energy (asymptotic energy)
    if (x1 <= 1d-10) then
        total_energy = get_coeff_zero(coeff_index)
        return
    endif

    ! Passing to the user coordinates to the 6D coordinates under Euler-ZXZ convension
    call user_coordinates_to_general_coordinates(general_coordinates_ZXZ,xdim,coord_format,coordinates)

    ! Evaluating the expansion in the general coordinates
    ! total energy is the sum of the individual contibutions
    total_energy = get_total_interaction_energy(coeff_index,general_coordinates_ZXZ)
    return

end subroutine evaluate_LRF


