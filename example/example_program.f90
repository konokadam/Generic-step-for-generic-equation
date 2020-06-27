!!!Author: H. Emrah Konokman

module equation_of_motion_class

    use :: precision, only: rp
    use :: step_class, only: derivative_abstract_type

    implicit none

    private

    !!!Point mass equations of motion type
    type, extends(derivative_abstract_type), public :: derivative_motion_point_mass_dt_type
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
    select case(trim(eq_motion_type))
        case("point_mass_dt")
            allocate(derivative_motion_point_mass_dt_type :: derivative_c)
        case default
            write(*,*) trim(eq_motion_type)//" not defined!"
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
        real(rp) :: g = 9.80665_rp
        real(rp) :: rho = 1.225_rp

        !!!x(1) = Vx
        !!!x(2) = Vy
        !!!x(3) = Vz
        !!!x(4) = x
        !!!x(5) = y
        !!!x(6) = z
        !!!x(7) = t

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




program test

    use :: precision, only: rp
    use :: step_class, only: derivative_abstract_type, step_abstract_type, step
    use :: equation_of_motion_class

    implicit none
    class(derivative_abstract_type), allocatable :: eq_of_motion_c   !!!Container for equation of motion
    class(step_abstract_type), allocatable :: step_c        !!!Container for step

    real(rp) :: x, y, z, t, Vx, Vy, Vz, area, mass, CD
    character(len=:), allocatable :: step_type, eq_motion_type
    real(rp) :: step_size, limit

    eq_motion_type = "point_mass_dt"
    step_size = 0.1_rp !!![s]
    limit = 12._rp      !!![s]

    step_type = "Euler"
    call params_init(x, y, z, t, Vx, Vy, Vz, area, mass, CD)

    call allocate_motion_derivative(eq_of_motion_c, eq_motion_type)
    select type(derivative_p => eq_of_motion_c)
        type is (derivative_motion_point_mass_dt_type)
            call params_derivative_motion_point_mass(derivative_p, mass, area, CD)
    end select
    call step(step_c, step_type)

    write(*,*) "Solving by "//step_type//" step with step_size =", step_size
    call motion_point_mass(step_c, eq_of_motion_c, step_size, limit, step_type//".dat")

    step_type = "RK4"
    call params_init(x, y, z, t, Vx, Vy, Vz, area, mass, CD)

    call allocate_motion_derivative(eq_of_motion_c, eq_motion_type)
    select type(derivative_p => eq_of_motion_c)
        type is (derivative_motion_point_mass_dt_type)
            call params_derivative_motion_point_mass(derivative_p, mass, area, CD)
    end select
    call step(step_c, step_type)

    write(*,*) "Solving by "//step_type//" step with step_size =", step_size
    call motion_point_mass(step_c, eq_of_motion_c, step_size, limit, step_type//".dat")

contains

    subroutine params_init(x, y, z, t, Vx, Vy, Vz, area, mass, CD)

        implicit none
        real(rp), intent(out) :: x, y, z, t, Vx, Vy, Vz, area, mass, CD

        x = 0._rp
        y = 0._rp
        z = 0._rp
        t = 0._rp
        Vx = 50._rp
        Vy = 0._rp
        Vz = 50._rp
        area = 0.25_rp * 2._rp* (5._rp*10._rp + 5._rp*15._rp + 10._rp*15._rp) *1e-6_rp
        mass = 7850._rp * (5._rp*10._rp*15._rp) * 1e-9_rp
        CD = 0.3_rp

    end subroutine params_init

    !!!Motion procedure
    subroutine motion_point_mass(step_c, eq_of_motion_c, dt, limit, file)

        implicit none
        class(derivative_abstract_type), intent(in) :: eq_of_motion_c   !!!Container for equation of motion
        class(step_abstract_type), intent(in) :: step_c        !!!Container for step
        real(rp), intent(in) :: dt, limit
        character(*), intent(in) :: file
        real(rp) :: s
        real(rp) :: dx(7)
        integer :: u

        open(newunit=u, file = trim(adjustl(file)))
        write(u,'(20a14)') "x", "y", "z" &
            , "t", "s", "Vx", "Vy", "Vz"
        s = 0._rp
        do while (t <= limit .and. z >= -0.1_rp)
            write(u,'(20f14.4)') x, y, z &
                , t, s, Vx, Vy, Vz
            call step_c%step(eq_of_motion_c, &
                            [Vx, Vy, Vz, x, y, z, t], dt, dx)

            Vx = Vx + dx(1)
            Vy = Vy + dx(2)
            Vz = Vz + dx(3)
            x = x + dx(4)
            y = y + dx(5)
            z = z + dx(6)
            t = t + dx(7)
            s = s + sqrt(dx(4)**2 + dx(5)**2 + dx(6)**2)
        end do
        close(u)

    end subroutine motion_point_mass

end program
