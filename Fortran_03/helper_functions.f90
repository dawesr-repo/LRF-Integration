!======================================================================
! Helper API: coordinate transforms + evaluate_LRF
!   - No OpenMP here (parallelism is inside Geometry_Constant_v2)
!======================================================================
module LRF_API
  use iso_fortran_env, only : real64, int32
  implicit none
  private
  public :: coordinate_transformation
  public :: general_coordinates_format
  public :: user_coordinates_to_general_coordinates
  public :: evaluate_LRF

contains

  subroutine coordinate_transformation(new_coordinates, coordinates, coord_format)
    use iso_fortran_env, only : real64
    implicit none
    real(real64), intent(in)  :: coordinates(6)
    character(*),  intent(in) :: coord_format
    real(real64), intent(out) :: new_coordinates(6)

    new_coordinates = coordinates

    if (coord_format == "Euler_ZYZ") then
      new_coordinates(5) = coordinates(5) - 90.0_real64
      new_coordinates(6) = coordinates(6) - 90.0_real64
    else if (coord_format == "Spherical") then
      new_coordinates(5) = 90.0_real64 - coordinates(5)
      new_coordinates(6) = 90.0_real64 - coordinates(6)
    end if
  end subroutine coordinate_transformation

  subroutine general_coordinates_format(general_coordinates, dim, old_coordinates)
    use iso_fortran_env, only : real64, int32
    implicit none
    integer(int32), intent(in) :: dim
    real(real64),   intent(in) :: old_coordinates(dim)
    real(real64),   intent(out):: general_coordinates(6)

    general_coordinates = 0.0_real64

    select case (dim)
    case (2)
      general_coordinates(1) = old_coordinates(1)  ! R
      general_coordinates(2) = old_coordinates(2)  ! b1
    case (3)
      general_coordinates(1) = old_coordinates(1)  ! R
      general_coordinates(2) = old_coordinates(2)  ! b1
      general_coordinates(5) = old_coordinates(3)  ! c1
    case (4)
      general_coordinates(1:4) = old_coordinates(1:4)  ! R, b1, b2, phi
    case (5)
      general_coordinates(1:5) = old_coordinates(1:5)  ! + c1
    case (6)
      general_coordinates(1:6) = old_coordinates(1:6)  ! + c2
    case default
      ! leave as zeros
    end select
  end subroutine general_coordinates_format

  subroutine user_coordinates_to_general_coordinates(general_coordinates_ZXZ, xdim, coord_format, user_coordinates)
    use iso_fortran_env, only : real64, int32
    implicit none
    integer(int32), intent(in) :: xdim
    real(real64),   intent(in) :: user_coordinates(xdim)
    character(*),   intent(in) :: coord_format
    real(real64),   intent(out):: general_coordinates_ZXZ(6)

    real(real64) :: general_coordinates(6)

    call general_coordinates_format(general_coordinates, xdim, user_coordinates)
    call coordinate_transformation(general_coordinates_ZXZ, general_coordinates, coord_format)
  end subroutine user_coordinates_to_general_coordinates

  !--------------------------------------------------------------------
  ! MAIN public entry: evaluate_LRF
  !   * Interaction sums may use OpenMP inside Geometry_Constant_v2
  !--------------------------------------------------------------------
  subroutine evaluate_LRF(total_energy, xdim, coordinates, coord_format, filename)
    use iso_fortran_env,      only : real64, int32
    use Fitting_Constant_v2,  only : get_coeff_zero, get_coeff_index
    use Geometry_Constant_v2, only : get_total_interaction_energy
    implicit none
    real(real64),   intent(out) :: total_energy
    integer(int32), intent(in)  :: xdim
    real(real64),   intent(in)  :: coordinates(xdim)
    character(*),   intent(in)  :: coord_format
    character(*),   intent(in)  :: filename

    real(real64)   :: general_coordinates_ZXZ(6)
    integer(int32) :: i, coeff_index
    real(real64)   :: x1

    coeff_index = get_coeff_index(filename)

    x1 = 0.0_real64
    do i = 1, xdim
      x1 = x1 + abs(coordinates(i))
    end do

    ! Return asymptote if all user coords are ~0
    if (x1 <= 1.0e-10_real64) then
      total_energy = get_coeff_zero(coeff_index)
      return
    end if

    ! Convert user coordinates to Euler-ZXZ 6D
    call user_coordinates_to_general_coordinates(general_coordinates_ZXZ, xdim, coord_format, coordinates)

    ! Evaluate contributions (multipole + induction + dispersion)
    total_energy = get_total_interaction_energy(coeff_index, general_coordinates_ZXZ)
  end subroutine evaluate_LRF

end module LRF_API
