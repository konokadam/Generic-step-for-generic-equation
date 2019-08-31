!!!Author: Huseyin Emrah Konokman

module point_mass_class

    use :: precision, only: rp
    use :: step_class, only: step_abstract_type, derivative_abstract_type, step
    use :: equation_of_motion_class
    use :: aero_params

    implicit none
    type :: point_mass_type
        real(rp) :: area, mass
        real(rp) :: x, y, z, t
        real(rp) :: Vx, Vy, Vz
        real(rp) :: CD
        class(derivative_abstract_type), allocatable :: eq_of_motion_c   !!!Container for equation of motion
        class(step_abstract_type), allocatable :: step_c        !!!Container for step
    contains
        procedure :: motion => point_mass_motion     !!!Procedure for the motion of the point mass
    end type point_mass_type

    public :: point_mass

    !!!Point mass constructor
    interface point_mass
        module procedure point_mass_params
    end interface

contains

    subroutine point_mass_params(point_mass, x, y, z, t, Vx, Vy, Vz, area, mass, shape, step_type, eq_motion_type) !result(point_mass)

        implicit none
        type(point_mass_type) :: point_mass
        real(rp), intent(in) :: x, y, z, t, Vx, Vy, Vz, area, mass
        integer, intent(in) :: shape
        character(*), intent(in) :: step_type, eq_motion_type

        point_mass%x = x
        point_mass%y = y
        point_mass%z = z
        point_mass%z = t
        point_mass%Vx = Vx
        point_mass%Vy = Vy
        point_mass%Vz = Vz
        point_mass%area = area
        point_mass%mass = mass
        call allocate_motion_derivative(point_mass%eq_of_motion_c, eq_motion_type)
        select type(derivative_p => point_mass%eq_of_motion_c)
        type is (derivative_motion_point_mass_dt_type)
            point_mass%CD = drag_data_point_mass(shape)
            call params_derivative_motion_point_mass(derivative_p, point_mass%mass, point_mass%area, point_mass%CD)
        end select
        call allocate_step_type(point_mass%step_c, step_type)

    end subroutine point_mass_params

    !!!Motion procedure
    subroutine point_mass_motion(point_mass, dt, limit, file)

        implicit none
        class(point_mass_type), intent(inout) :: point_mass
        real(rp), intent(in) :: dt, limit
        character(*), intent(in) :: file
        real(rp) :: s
        real(rp) :: dx(7)
        integer :: u

        open(newunit=u, file = trim(adjustl(file)))
        write(u,'(20a14)') "pm%x", "pm%y", "pm%z" &
        , "pm%t", "s", "pm%Vx", "pm%Vy", "pm%Vz"
        s = 0._rp
        do while (point_mass%t <= limit .and. point_mass%z >= 0._rp)
            write(u,'(20f14.4)') point_mass%x, point_mass%y, point_mass%z &
            , point_mass%t, s, point_mass%Vx, point_mass%Vy, point_mass%Vz
            call point_mass%step_c%step(point_mass%eq_of_motion_c &
            , [point_mass%Vx, point_mass%Vy, point_mass%Vz, point_mass%x, point_mass%y, point_mass%z, point_mass%t], dt, dx)

            point_mass%Vx = point_mass%Vx + dx(1)
            point_mass%Vy = point_mass%Vy + dx(2)
            point_mass%Vz = point_mass%Vz + dx(3)
            point_mass%x = point_mass%x + dx(4)
            point_mass%y = point_mass%y + dx(5)
            point_mass%z = point_mass%z + dx(6)
            point_mass%t = point_mass%t + dx(7)
            s = s + sqrt(dx(4)**2 + dx(5)**2 + dx(6)**2)
       end do
       close(u)

    end subroutine point_mass_motion

end module point_mass_class


program test

    use :: precision, only: rp
    use :: point_mass_class

    type(point_mass_type) :: pm
    real(rp) :: x, y, z, t, Vx, Vy, Vz, area, mass
    integer :: shape
    character(len=:), allocatable :: step_type, eq_motion_type
    real(rp) :: step_size, limit

    eq_motion_type = "point_mass_dt"
    step_size = 0.01_rp !!![s]
    limit = 12._rp      !!![s]

    step_type = "Euler"
    call params_init(x, y, z, t, Vx, Vy, Vz, area, mass, shape)
    call point_mass(pm, x, y, z, t, Vx, Vy, Vz, area, mass, shape, step_type, eq_motion_type)

    write(*,*) "Solving by "//step_type//" step with step_size =", step_size
    call pm%motion(step_size, limit, step_type//".dat")

    step_type = "RK4"
    call params_init(x, y, z, t, Vx, Vy, Vz, area, mass, shape)
    call point_mass(pm, x, y, z, t, Vx, Vy, Vz, area, mass, shape, step_type, eq_motion_type)

    write(*,*) "Solving by "//step_type//" step with step_size =", step_size
    call pm%motion(step_size, limit, step_type//".dat")

contains

    subroutine params_init(x, y, z, t, Vx, Vy, Vz, area, mass, shape)

        implicit none
        real(rp), intent(out) :: x, y, z, t, Vx, Vy, Vz, area, mass
        integer, intent(out) :: shape

        x = 0._rp
        y = 0._rp
        z = 0._rp
        t = 0._rp
        Vx = 1000._rp
        Vy = 0._rp
        Vz = 1000._rp
        area = 0.25_rp * 2._rp* (5._rp*10._rp + 5._rp*15._rp + 10._rp*15._rp) *1e-6_rp
        mass = 7850._rp * (5._rp*10._rp*15._rp) * 1e-9_rp
        shape = 1

    end subroutine params_init

end program test
