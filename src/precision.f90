module precision
    
    private
    integer, parameter :: real16 = selected_real_kind(30,4000)
    integer, parameter, public :: rp = real16
    
end module precision
