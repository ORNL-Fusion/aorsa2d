module check_mod

use netcdf

contains

    subroutine check ( status )

            integer, intent(in) :: status
          
            if(status /= nf90_noerr) then 
                print *, trim(nf90_strerror(status))
                stop "Stopped"
            end if
    
    end subroutine check

end module check_mod
