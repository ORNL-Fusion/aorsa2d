module derivatives

implicit none

contains


    subroutine deriv_x ( f, i, j, dx, dfdx, d2fdx2 )
    
        implicit none
        
        integer, intent(in) :: i, j
        real, intent(in) :: f(:,:)
        real, intent(out), optional :: dfdx, d2fdx2
        real, intent(in) :: dx
    
        integer :: nx, ny
    
        nx    = size ( f, dim = 1 )
        ny    = size ( f, dim = 2 )

        ! 2D
        if(nx/=1)then

            if(i .ne. 1 .and. i .ne. nx)then

                dfdx = (f(i+1, j) - f(i-1, j)) / (2.0 * dx)
                d2fdx2 = (f(i+1, j) - 2.0 * f(i,j) + f(i-1, j)) / dx**2

            endif
            
            if(i .eq. 1)then

                dfdx = (f(i+1, j) - f(i, j)) / dx
                d2fdx2 = (f(i+2, j) - 2.0 * f(i+1,j) + f(i, j)) / dx**2

            endif
            
            if(i .eq. nx)then

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
    subroutine deriv_y ( f, i, j, dy, dfdy, d2fdy2 )

        implicit none

        integer, intent(in) :: i, j
        real, intent(in) :: f(:,:)
        real, intent(in) :: dy
        real, intent(out), optional :: dfdy, d2fdy2

        integer :: nx, ny

        nx    = size ( f, dim = 1 )
        ny    = size ( f, dim = 2 )

        !2D
        if(ny/=1)then
            if(j .ne. 1 .and. j .ne. ny)then
                dfdy = (f(i, j+1) - f(i, j-1)) / (2.0 * dy)
                d2fdy2 = (f(i, j+1) - 2.0 * f(i,j) + f(i, j-1)) / dy**2
            endif

            if(j .eq. 1)then
                dfdy = (f(i, j+1) - f(i, j)) / dy
                d2fdy2 = (f(i, j+2) - 2.0 * f(i,j+1) + f(i, j)) / dy**2
            endif


            if(j .eq. ny)then
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
      subroutine deriv_xy ( f, i, j, dx, dy, d2fdxy)

      implicit none

      integer, intent(in) :: i, j
      real, intent(in) :: f(:,:)
      real, intent(in) :: dx, dy
      real, intent(out), optional :: d2fdxy

      integer :: nx, ny

      nx    = size ( f, dim = 1 )
      ny    = size ( f, dim = 2 )

      d2fdxy = 0.0

      if(nx/=1 .and. ny/=1)then

        if (i /= 1 .and. i /= nx .and. &
            j /= 1 .and. j /= ny) then

           d2fdxy = (f(i+1,j+1) - f(i-1,j+1) - f(i+1,j-1) + f(i-1,j-1) ) &
              / (4.0 * dx * dy)

        endif

      endif


      end subroutine deriv_xy



end module derivatives
