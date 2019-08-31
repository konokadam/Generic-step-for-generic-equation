module precision
    
    private
    !integer, parameter :: real4 = selected_real_kind(6,30)         !.. minimum 6  decimal digits or 10^30,    4 bytes
    !integer, parameter :: real8 = selected_real_kind(14,300)        !.. minimum 14 decimal digits or 10^300,   8 bytes
    integer, parameter, p :: real16 = selected_real_kind(30,4000)      !.. minimum 30 decimal digits or 10^4000, 16 bytes

!    integer, parameter, public :: rp = real4                    !.. working precision for reals
!    integer, parameter, public :: rp = real8                    !.. working precision for reals
    integer, parameter, public :: rp = real16                    !.. working precision for reals
    
end module precision
