module generic_biLinearInterp

use constants, only: dbl

implicit none

interface biLinearInterp

    module procedure biLinearInterp_r
    module procedure biLinearInterp_i
    module procedure biLinearInterp_l

end interface biLinearInterp

contains

    function biLinearInterp_r ( r, z, g, f2D )

        use grid, only: gridBlock

        implicit none

        real :: r, z
        type(gridBlock), intent(in) :: g
        real, intent(in) :: f2D(:,:)

        real :: biLinearInterp_r

        integer :: i1,i2, j1,j2
        real :: i, j
        real :: f11, f12, f21, f22
        real :: r1, r2, z1, z2
        real :: denom
        real :: y0,y1,x0,x1,x
        real :: threshold

        i = ((r-g%rMin)/g%rRange*(g%nR-1))+1
        j = ((z-g%zMin)/g%zRange*(g%nZ-1))+1

        ! All the if statements etc are catches for the cases
        ! where you request a point that is at integral values
        ! of i and or j.

        ! No interp
        if(floor(i)==ceiling(i).and.floor(j)==ceiling(j)) then
            biLinearInterp_r = f2D(i,j)
            return

        ! Linear in z
        elseif(floor(i) == ceiling(i)) then 
            j1 = floor(j)
            j2 = ceiling(j)
            y0 = f2D(i,j1)
            y1 = f2D(i,j2)
            x0 = g%z(j1)
            x1 = g%z(j2)
            x  = r
            biLinearInterp_r = &
                y0 + (x-x0)*(y1-y0)/(x1-x0)
            return

        ! Linear in r
        elseif(floor(j) == ceiling(j)) then
            i1 = floor(i)
            i2 = ceiling(i)
            y0 = f2D(i1,j)
            y1 = f2D(i2,j)
            x0 = g%r(i1)
            x1 = g%z(i2)
            x  = z
            biLinearInterp_r = &
                y0 + (x-x0)*(y1-y0)/(x1-x0)
            return

        ! General biLinear
        else 
            i1 = floor(i)
            i2 = ceiling(i)

            j1 = floor(j)
            j2 = ceiling(j)

            f11 = f2D(i1,j1)
            f12 = f2D(i1,j2)
            f21 = f2D(i2,j1)
            f22 = f2D(i2,j2)

            r1 = g%r(i1)
            r2 = g%r(i2)
            z1 = g%z(j1)
            z2 = g%z(j2)

            denom = ((r2-r1)*(z2-z1))

            biLinearInterp_r = &
               f11/denom*(r2-r)*(z2-z) &
              +f21/denom*(r-r1)*(z2-z) &
              +f12/denom*(r2-r)*(z-z1) &
              +f22/denom*(r-r1)*(z-z1)
            return

        endif

    end function biLinearInterp_r

    function biLinearInterp_i ( r, z, g, f2D )

        use grid, only: gridBlock
        use constants, only: dbl

        implicit none

        real :: r, z
        type(gridBlock), intent(in) :: g
        integer, intent(in) :: f2D(:,:)

        real :: biLinearInterp_i

        integer :: i1,i2, j1,j2
        real :: i, j
        real :: f11, f12, f21, f22
        real :: r1, r2, z1, z2
        real :: denom
        real :: y0,y1,x0,x1,x
        real :: threshold
        real, allocatable :: f2D_temp(:,:)

        allocate( f2D_temp(size(f2D,1),size(f2D,2)) )
        f2D_temp = f2D

        i = ((r-g%rMin)/g%rRange*(g%nR-1))+1
        j = ((z-g%zMin)/g%zRange*(g%nZ-1))+1

        ! All the if statements etc are catches for the cases
        ! where you request a point that is at integral values
        ! of i and or j.

        ! No interp
        if(floor(i)==ceiling(i).and.floor(j)==ceiling(j)) then
            biLinearInterp_i = f2D_temp(i,j)
            return

        ! Linear in z
        elseif(floor(i) == ceiling(i)) then 
            j1 = floor(j)
            j2 = ceiling(j)
            y0 = f2D_temp(i,j1)
            y1 = f2D_temp(i,j2)
            x0 = g%z(j1)
            x1 = g%z(j2)
            x  = r
            biLinearInterp_i = &
                y0 + (x-x0)*(y1-y0)/(x1-x0)
            return

        ! Linear in r
        elseif(floor(j) == ceiling(j)) then
            i1 = floor(i)
            i2 = ceiling(i)
            y0 = f2D_temp(i1,j)
            y1 = f2D_temp(i2,j)
            x0 = g%r(i1)
            x1 = g%z(i2)
            x  = z
            biLinearInterp_i = &
                y0 + (x-x0)*(y1-y0)/(x1-x0)
            return

        ! General biLinear
        else 
            i1 = floor(i)
            i2 = ceiling(i)

            j1 = floor(j)
            j2 = ceiling(j)

            f11 = f2D_temp(i1,j1)
            f12 = f2D_temp(i1,j2)
            f21 = f2D_temp(i2,j1)
            f22 = f2D_temp(i2,j2)

            r1 = g%r(i1)
            r2 = g%r(i2)
            z1 = g%z(j1)
            z2 = g%z(j2)

            denom = ((r2-r1)*(z2-z1))

            biLinearInterp_i = &
               f11/denom*(r2-r)*(z2-z) &
              +f21/denom*(r-r1)*(z2-z) &
              +f12/denom*(r2-r)*(z-z1) &
              +f22/denom*(r-r1)*(z-z1)
            return

        endif

    end function biLinearInterp_i

    function biLinearInterp_l ( r, z, g, f2D )

        use grid, only: gridBlock
        use constants, only: dbl

        implicit none

        real :: r, z
        type(gridBlock), intent(in) :: g
        logical, intent(in) :: f2D(:,:)

        real :: biLinearInterp_l

        integer :: i1,i2, j1,j2
        real :: i, j
        real :: f11, f12, f21, f22
        real :: r1, r2, z1, z2
        real :: denom
        real :: y0,y1,x0,x1,x
        real :: threshold
        real, allocatable :: f2D_temp(:,:)
        integer :: ii, jj

        allocate( f2D_temp(size(f2D,1),size(f2D,2)) )
        f2D_temp = 0
        where(f2D) f2D_temp = 1

        i = ((r-g%rMin)/g%rRange*(g%nR-1))+1
        j = ((z-g%zMin)/g%zRange*(g%nZ-1))+1

        ! All the if statements etc are catches for the cases
        ! where you request a point that is at integral values
        ! of i and or j.

        ! No interp
        if(floor(i)==ceiling(i).and.floor(j)==ceiling(j)) then
            biLinearInterp_l = f2D_temp(i,j)
            return

        ! Linear in z
        elseif(floor(i) == ceiling(i)) then 
            j1 = floor(j)
            j2 = ceiling(j)
            y0 = f2D_temp(i,j1)
            y1 = f2D_temp(i,j2)
            x0 = g%z(j1)
            x1 = g%z(j2)
            x  = r
            biLinearInterp_l = &
                y0 + (x-x0)*(y1-y0)/(x1-x0)
            return

        ! Linear in r
        elseif(floor(j) == ceiling(j)) then
            i1 = floor(i)
            i2 = ceiling(i)
            y0 = f2D_temp(i1,j)
            y1 = f2D_temp(i2,j)
            x0 = g%r(i1)
            x1 = g%z(i2)
            x  = z
            biLinearInterp_l = &
                y0 + (x-x0)*(y1-y0)/(x1-x0)
            return

        ! General biLinear
        else 
            i1 = floor(i)
            i2 = ceiling(i)

            j1 = floor(j)
            j2 = ceiling(j)

            f11 = f2D_temp(i1,j1)
            f12 = f2D_temp(i1,j2)
            f21 = f2D_temp(i2,j1)
            f22 = f2D_temp(i2,j2)

            r1 = g%r(i1)
            r2 = g%r(i2)
            z1 = g%z(j1)
            z2 = g%z(j2)

            denom = ((r2-r1)*(z2-z1))

            biLinearInterp_l = &
               f11/denom*(r2-r)*(z2-z) &
              +f21/denom*(r-r1)*(z2-z) &
              +f12/denom*(r2-r)*(z-z1) &
              +f22/denom*(r-r1)*(z-z1)
            return

        endif

    end function biLinearInterp_l

end module generic_biLinearInterp
