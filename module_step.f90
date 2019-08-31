module step_class

    !!!Author: Huseyin Emrah Konokman

    use :: precision, only: rp
    
    implicit none
    
    private 
    
    !!!The equation (derivatives of the equation set) to be integrated (to be stepped forward)
    !! The type for the equation set to be stepped forward should be extended from this abstract type
    type, abstract, public :: derivative_abstract_type
        !!!Parameters of the (derivatives of) equation set should be given by the extended type of that equation set
    contains
        procedure(derivative_proc), deferred :: derivative    !!!The procedure of the derivative equation to be stepped forward. This procedure takes x input and output xdot
    end type derivative_abstract_type

    abstract interface
        subroutine derivative_proc(self, h, x, xdot)
            import :: derivative_abstract_type, rp
            class(derivative_abstract_type), intent(inout) :: self
            real(rp),intent(in) :: h                    !!!Step size
            real(rp), intent(in) :: x(:)                !!!Vector of the equation set
            real(rp), intent(out) :: xdot(size(x))      !!!Vector of the derivatives of equation set
        end subroutine
    end interface


    
    !!!Abstract step type
    !! The type for a new step method is extended from this abstract type
    type, abstract, public :: step_abstract_type
        real(rp) :: ds          !!!Step size
    contains
        procedure(step_proc), deferred, nopass :: step
    end type step_abstract_type

    abstract interface
        subroutine step_proc(f, x, h, dx)
            use :: precision, only: rp
            import :: derivative_abstract_type
            class(derivative_abstract_type) :: f        !!!The derivative function come through overloaded deferred procedure f%derivative(...)
            real(rp), intent(in) :: x(:)                !!!Vector of the equation set
            real(rp), intent(in) :: h                   !!!Step size
            real(rp), intent(out) :: dx(:)              !!!Vector of one step of the equation set
        end subroutine step_proc
    end interface

    !!!Euler step
    type, extends(step_abstract_type), public :: step_euler_type

    contains
        procedure, nopass :: step => step_euler
    end type step_euler_type

    !!!RK4 step
    type, extends(step_abstract_type), public :: step_rk4_type

    contains
        procedure, nopass :: step => step_rk4
    end type step_rk4_type
    
    !!!Any new step types can be added here>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    
    
    !!!Any new step types can be added here<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
    
    

    public :: step, allocate_step_type
    
    !!!Constructor
    interface step
        module procedure allocate_step_type
    end interface
    
contains

    !!!Procedure for the allocation of the step container according to desired step type
    subroutine allocate_step_type(step_c, step_type, ds)

        implicit none
        class(step_abstract_type), allocatable, intent(inout) :: step_c     !!!Step type container
        character(*), intent(in) :: step_type
        real(rp), intent(in), optional :: ds

        if(allocated(step_c)) deallocate(step_c)
        select case(trim(adjustl(step_type)))
        case("RK4")
            allocate(step_rk4_type :: step_c)
        case("Euler")
            allocate(step_euler_type :: step_c)
        case default
                allocate(step_rk4_type :: step_c)
        end select
        
        if(present(ds)) step_c%ds = ds

    end subroutine

    !!!Procedures for the steps
    
    !!!Procedure for simple Euler step
    subroutine step_euler(f, x, h, dx)

        implicit none
        class(derivative_abstract_type) :: f
        real(rp), intent(in) :: x(:)
        real(rp), intent(in) :: h
        real(rp), intent(out) :: dx(:)

        real(rp), dimension(size(x)) :: xdot

        call f%derivative(h, x, xdot)
        dx = h * xdot

    end subroutine step_euler

    !!!Procedure for 4th order Runge-Kutta step
    subroutine step_rk4(f, x, h, dx)

        implicit none
        class(derivative_abstract_type) :: f
        real(rp), intent(in) :: x(:)
        real(rp), intent(in) :: h
        real(rp), intent(out) :: dx(:)

        real(rp), dimension(size(x)) :: xdot1, xdot2, xdot3, xdot4
        real(rp) :: ho2

        ho2 = 0.5_rp * h

        call f%derivative(0._rp, x, xdot1)           !!!1st RK4 step
        call f%derivative(ho2, x+ho2*xdot1, xdot2)   !!!2nd RK4 step
        call f%derivative(ho2, x+ho2*xdot2, xdot3)   !!!3rd RK4 step
        call f%derivative(h, x+h*xdot3, xdot4)       !!!4th RK4 step

        dx = h*(xdot1 + 2._rp*xdot2 + 2._rp*xdot3 + xdot4)/6._rp

    end subroutine step_rk4
    
    !!!Procedures for any new step types can be added here>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    
    
    !!!Procedures for any new step types can be added here<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    
end module step_class
