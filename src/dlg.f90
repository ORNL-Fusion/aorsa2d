module dlg

contains
    function dlg_pDeriv ( array, dir, dS ) result ( dArray )
        implicit none
        
        real, dimension (:,:), allocatable :: dArray
        integer :: nX, nY, i, j
        integer, intent(IN) :: dir
        real, intent(IN) :: array (:,:), dS

        nX  = size ( array, 1 )
        nY  = size ( array, 2 )

        allocate ( dArray ( nX, nY ) )
        dArray  = 0.0
      
        if ( dir == 2 ) then
            
            do i=1,nX
                do j=1,nY

                    if ( j > 1 .AND. j < nY ) &
                        dArray(i,j) = 1.0 / ( 2.0 * dS ) * &
                            ( array(i,j+1) - array(i,j-1) )

                    if ( j == 1 ) &
                        dArray(i,j) = 1.0 / ( 2.0 * dS ) * &
                            ( &
                            -3.0 * array(i,j) + 4.0 * array(i,j+1) &
                            - array(i,j+2) )

                    if ( j == nY ) &
                        dArray(i,j) = 1.0 / ( 2.0 * dS ) * &
                            ( &
                            3.0 * array(i,j) - 4.0 * array(i,j-1) &
                            + array(i,j-2) )

                end do
            end do
       
        else if ( dir == 1 ) then
            
            do i=1,nX
                do j=1,nY

                    if ( i > 1 .AND. i < nX ) &
                        dArray(i,j) = 1.0 / ( 2.0 * dS ) * &
                            ( array(i+1,j) - array(i-1,j) )

                    if ( i == 1 ) &
                        dArray(i,j) = 1.0 / ( 2.0 * dS ) * &
                            ( &
                            -3.0 * array(i,j) + 4.0 * array(i+1,j) &
                            - array(i+2,j) )

                    if ( i == nX ) &
                        dArray(i,j) = 1.0 / ( 2.0 * dS ) * &
                            ( &
                            3.0 * array(i,j) - 4.0 * array(i-1,j) &
                            + array(i-1,j) )

                end do
            end do
       
        end if 
        
    end function dlg_pDeriv

end module dlg
