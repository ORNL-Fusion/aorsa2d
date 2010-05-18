module chebyshev_mod

use constants

contains
    function chebT (n,x)

        implicit none

        integer, intent(in) :: n
        real, intent(in) :: x
        real :: chebT
 
        if(abs(x)>1) then
            write(*,*) 'ERROR: chebT arg abs(x)>1', x
        endif

        if(x==1.0)then
            chebT = 1
        elseif(x==-1.0)then
            chebT = (-1)**n 
        else
            !chebT = ( ( x - sqrt ( x**2 - 1 ) )**n &
            !    + ( x + sqrt ( x**2 - 1 ) )**n ) / 2d0
            chebT = cos ( n * acos(x) )
        endif


    end function chebT
   

    function chebU (n,x)
 
        implicit none

        integer, intent(in) :: n
        real, intent(in) :: x
        real(kind=dbl) :: chebU
 
        if(x==1.0) then
            chebU   = n + 1
        elseif(x==-1.0) then
            chebU   = (n+1) * (-1)**n
        else
               
            !chebU   = ( ( x + sqrt ( x**2 - 1 ) )**(n+1) &
            !    - ( x - sqrt ( x**2 - 1 ) )**(n+1) ) &
            !    / ( 2d0 * sqrt ( x**2 - 1 ) )

            chebU   = sin ( (n+1) * acos(x) ) / sin ( acos(x) )

        endif

    
    end function chebU

end module chebyshev_mod
