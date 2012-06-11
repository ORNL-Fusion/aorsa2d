module timer_mod

implicit none

integer, parameter :: t_dbl = selected_real_kind ( p = 14 )

type, public :: timer
    
    real ( kind = t_dbl ) :: time

end type timer

contains
    subroutine start_timer ( start_time )
    
        implicit none
    
        type ( timer ) :: start_time
    
        integer, dimension ( 8 ) :: value
    
        call date_and_time ( values = value )
    
        start_time%time = 86400.d0 * value ( 3 ) + 3600.d0 * value ( 5 ) &
            + 60.d0 * value ( 6 ) + value ( 7 ) + 0.001d0 * value ( 8 )
    
    end subroutine start_timer


    real function end_timer ( start_time )
    
        implicit none
            
        type ( timer ) :: start_time, current_time
    
        integer, dimension ( 8 ) :: value
        !real ( kind = dbl ) :: curent_time
    
        call date_and_time ( values = value )
        current_time%time   = 86400.d0 * value ( 3 ) + 3600.d0 * value ( 5 ) &
            + 60.d0 * value ( 6 ) + value ( 7 ) + 0.001d0 * value ( 8 )
    
        end_timer   = current_time%time - start_time%time
    
    end function end_timer


end module timer_mod
