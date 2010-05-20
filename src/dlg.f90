module dlg

contains
    function dlg_pDeriv ( x, y, array, dir ) result ( dArray )
        implicit none
        
        real, dimension (:,:), allocatable :: dArray
        integer :: nX, nY, i, j
        integer, intent(IN) :: dir
        real, intent(IN) :: array (:,:), x(:), y(:)
        real :: dS

        nX  = size ( array, 1 )
        nY  = size ( array, 2 )

        allocate ( dArray ( nX, nY ) )
        dArray  = 0.0
      
        if ( dir == 2 ) then
            
            do i=1,nX
                do j=1,nY

                    if ( j > 1 .AND. j < nY ) then

                        dS = ((y(j+1)-y(j)) + (y(j)-y(j-1))) / 2.0

                        dArray(i,j) = 1.0 / ( 2.0 * dS ) * &
                            ( array(i,j+1) - array(i,j-1) )

                    elseif ( j == 1 ) then

                        dS = ((y(j+1)-y(j)) + (y(j+2)-y(j+1))) / 2.0

                        dArray(i,j) = 1.0 / ( 2.0 * dS ) * &
                            ( &
                            -3.0 * array(i,j) + 4.0 * array(i,j+1) &
                            - array(i,j+2) )

                    elseif ( j == nY ) then

                        dS = ((y(j-1)-y(j-2)) + (y(j)-y(j-1))) / 2.0

                        dArray(i,j) = 1.0 / ( 2.0 * dS ) * &
                            ( &
                            3.0 * array(i,j) - 4.0 * array(i,j-1) &
                            + array(i,j-2) )

                    endif

                enddo
            enddo
       
        elseif ( dir == 1 ) then
            
            do i=1,nX
                do j=1,nY

                    if ( i > 1 .AND. i < nX ) then

                        dS = ((x(i+1)-x(i)) + (x(i)-x(i-1))) / 2.0

                        dArray(i,j) = 1.0 / ( 2.0 * dS ) * &
                            ( array(i+1,j) - array(i-1,j) )

                    elseif ( i == 1 ) then

                        dS = ((x(i+1)-x(i)) + (x(i+2)-x(i+1))) / 2.0

                        dArray(i,j) = 1.0 / ( 2.0 * dS ) * &
                            ( &
                            -3.0 * array(i,j) + 4.0 * array(i+1,j) &
                            - array(i+2,j) )

                    elseif ( i == nX ) then

                        dS = ((x(i-1)-x(i-2)) + (x(i)-x(i-1))) / 2.0

                        dArray(i,j) = 1.0 / ( 2.0 * dS ) * &
                            ( &
                            3.0 * array(i,j) - 4.0 * array(i-1,j) &
                            + array(i-1,j) )

                    endif

                enddo
            enddo
       
        end if 
        
    end function dlg_pDeriv

end module dlg
