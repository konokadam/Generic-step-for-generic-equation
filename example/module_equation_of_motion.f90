 !!!Author: Huseyin Emrah Konokman
 
 module aero_params

    use :: precision, only: rp

    implicit none

contains

    function drag_data_point_mass(shape) result(CD)

        implicit none
        real(rp) :: CD
        integer, intent(in) :: shape

        select case (shape)
        case (1)        !.. x shape
            CD = 1.1_rp     !!!Assuming average CD with no change with Mach number
        !case ...
        case default
            write(*,*) "CD data for fiven shape:", shape, " is not given in the code!"
        end select

    end function drag_data_point_mass

end module aero_params

module equation_of_motion_class

    use :: precision, only: rp
    use :: aero_params, only: drag_data_point_mass
    use :: step_class, only: derivative_abstract_type

    implicit none

    private

    !!!Equation of motion abstract type
    type, abstract, extends(derivative_abstract_type) :: derivative_motion_abstract_type

    end type derivative_motion_abstract_type

    !!!Point mass equations of motion type
    type, extends(derivative_motion_abstract_type), public :: derivative_motion_point_mass_dt_type
        !!!Parameters of the derivatives of the equation set
        real(rp) :: area, mass
        real(rp) :: CD
    contains
        procedure :: derivative => derivative_motion_point_mass_dt
    end type derivative_motion_point_mass_dt_type

    public :: allocate_motion_derivative, params_derivative_motion_point_mass

contains

    !!!Procedure for the allocation of the derivative container according to desired equation of motion type
    subroutine allocate_motion_derivative(derivative_c, eq_motion_type)

        implicit none
        class(derivative_abstract_type), allocatable, intent(inout) :: derivative_c     !!!Derivative type container
        character(*), intent(in) :: eq_motion_type

        if(allocated(derivative_c)) deallocate(derivative_c)
        select case(trim(adjustl(eq_motion_type)))
            case("point_mass_dt")
                allocate(derivative_motion_point_mass_dt_type :: derivative_c)
            case default
                write(*,*) 'Derivative for equations of motion of type "'//trim(adjustl(eq_motion_type)) &
                //'" does not exist in the current version of the code! &
                & (subroutine allocate_motion_derivative(derivative_c, eq_motion_type))'
                !allocate(derivative_motion_point_mass_ds_type :: derivative_c)
        end select

    end subroutine


    subroutine params_derivative_motion_point_mass(self, mass, area, CD)

        implicit none
        class(derivative_motion_point_mass_dt_type), intent(inout) :: self
        real(rp), intent(in) :: mass, area
        real(rp), intent(in) :: CD

        self%mass = mass
        self%area = area
        self%CD = CD

    end subroutine params_derivative_motion_point_mass

    subroutine derivative_motion_point_mass_dt(self, h, x, xdot)

        implicit none
        class(derivative_motion_point_mass_dt_type), intent(inout) :: self
        real(rp), intent(in) :: h
        real(rp), intent(in) :: x(:)
        real(rp), intent(out) :: xdot(size(x))

        real(rp) :: drago
        real(rp) :: V
        real(rp) :: g = 9.80665_rp      !!!Assuming no change with altitude
        real(rp) :: rho = 1.225_rp      !!!Assuming sea level conditions with no change with altitude

        !!!x(1) = Vx
        !!!x(2) = Vy
        !!!x(3) = Vz
        !!!x(4) = x
        !!!x(5) = y
        !!!x(6) = z
        !!!x(7) = s

        V = sqrt(x(1)**2 + x(2)**2 + x(3)**2)
        drago = 0.5_rp * rho * self%area/self%mass * self%CD

        xdot(1) = - (drago * V * x(1))       !!!Vxdot
        xdot(2) = - (drago * V * x(2))       !!!Vydot
        xdot(3) = - (drago * V * x(3) + g)   !!!Vzdot
        xdot(4) = x(1)                       !!!xdot
        xdot(5) = x(2)                       !!!ydot
        xdot(6) = x(3)                       !!!zdot
        xdot(7) = 1._rp                      !!!tdot

    end subroutine derivative_motion_point_mass_dt

end module equation_of_motion_class
