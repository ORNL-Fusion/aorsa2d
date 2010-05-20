module derivatives

implicit none

contains


    subroutine deriv_x ( x, f, i, j, dfdx, d2fdx2 )
    
        implicit none
        
        integer, intent(in) :: i, j
        real, intent(in) :: f(:,:), x(:)
        real, intent(out), optional :: dfdx, d2fdx2
        real :: dx
    
        integer :: nx, ny
    
        nx    = size ( f, dim = 1 )
        ny    = size ( f, dim = 2 )

        ! 2D
        if(nx/=1)then

            if(i .ne. 1 .and. i .ne. nx)then

                dx = ( (x(i)-x(i-1)) + (x(i+1)-x(i)) ) / 2
                dfdx = (f(i+1, j) - f(i-1, j)) / (2.0 * dx)
                d2fdx2 = (f(i+1, j) - 2.0 * f(i,j) + f(i-1, j)) / dx**2

            endif
            
            if(i .eq. 1)then

                dx = ( (x(i+2)-x(i+1)) + (x(i+1)-x(i)) ) / 2
                dfdx = (f(i+1, j) - f(i, j)) / dx
                d2fdx2 = (f(i+2, j) - 2.0 * f(i+1,j) + f(i, j)) / dx**2

            endif
            
            if(i .eq. nx)then

                dx = ( (x(i)-x(i-1)) + (x(i-1)-x(i-2)) ) / 2
                dfdx = (f(i, j) - f(i-1, j)) / dx
                d2fdx2 = (f(i, j) - 2.0 * f(i-1,j) + f(i-2, j)) / dx**2

            endif
        !1D
        else

            dfdx = 0 
            d2fdx2 = 0 

        endif
    
    end subroutine deriv_x

!
!***************************************************************************
!
    subroutine deriv_y ( y, f, i, j, dfdy, d2fdy2 )

        implicit none

        integer, intent(in) :: i, j
        real, intent(in) :: f(:,:), y(:)
        real, intent(out), optional :: dfdy, d2fdy2

        real :: dy
        integer :: nx, ny

        nx    = size ( f, dim = 1 )
        ny    = size ( f, dim = 2 )

        !2D
        if(ny/=1)then
            if(j .ne. 1 .and. j .ne. ny)then

                dy = ( (y(j)-y(j-1)) + (y(j+1)-y(j)) ) / 2
                dfdy = (f(i, j+1) - f(i, j-1)) / (2.0 * dy)
                d2fdy2 = (f(i, j+1) - 2.0 * f(i,j) + f(i, j-1)) / dy**2
            endif

            if(j .eq. 1)then
                dy = ( (y(j+2)-y(j+1)) + (y(j+1)-y(j)) ) / 2
                dfdy = (f(i, j+1) - f(i, j)) / dy
                d2fdy2 = (f(i, j+2) - 2.0 * f(i,j+1) + f(i, j)) / dy**2
            endif


            if(j .eq. ny)then
                dy = ( (y(j)-y(j-1)) + (y(j-1)-y(j-2)) ) / 2
                dfdy = (f(i, j) - f(i, j-1)) / dy
                d2fdy2 = (f(i, j) - 2.0 * f(i,j-1) + f(i, j-2)) / dy**2
            endif

        !1D
        else

            dfdy = 0 
            d2fdy2 = 0 

        endif

      end subroutine deriv_y

!
!***************************************************************************
!
      subroutine deriv_xy ( x, y, f, i, j, d2fdxy)

      implicit none

      integer, intent(in) :: i, j
      real, intent(in) :: f(:,:), x(:), y(:)
      real, intent(out), optional :: d2fdxy

      real :: dx, dy
      integer :: nx, ny

      nx    = size ( f, dim = 1 )
      ny    = size ( f, dim = 2 )

      d2fdxy = 0.0

      if(nx/=1 .and. ny/=1)then

        if (i /= 1 .and. i /= nx .and. &
            j /= 1 .and. j /= ny) then

            dx = ( (x(i)-x(i-1)) + (x(i+1)-x(i)) ) / 2
            dy = ( (y(j)-y(j-1)) + (y(j+1)-y(j)) ) / 2

           d2fdxy = (f(i+1,j+1) - f(i-1,j+1) - f(i+1,j-1) + f(i-1,j-1) ) &
              / (4.0 * dx * dy)

        endif

      endif


      end subroutine deriv_xy



end module derivatives
