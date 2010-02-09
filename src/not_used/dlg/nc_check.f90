module nc_check

contains

    subroutine check ( status )
            use netcdf
            integer, intent ( in) :: status
          
            if(status /= nf90_noerr) then 
                print *, trim(nf90_strerror(status))
                stop "Stopped"
            end if
    
    end subroutine check

end module nc_check
