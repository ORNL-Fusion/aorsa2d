program test_spline

use spline_dlg
use netcdf
use check_mod

    real, allocatable :: x(:), y(:), x2(:), s(:)
    integer :: n, n2, i
    integer :: nc_id, n_id, x_id, y_id, s_id, x2_id, n2_id
    character(len=100) :: fName
    real, allocatable :: splines(:,:)

    n = 10 
    n2 = 9

    allocate (x(n),y(n))
    allocate (x2(n2),s(n2))
    allocate ( splines(n,4) )

    x = (/ (i*2*3.14/(n-1),i=0,n-1) /)
    y = cos(x)
    x2 = x(1:n-1)+(x(2)-x(1))/2 

    splines = spline ( x, y, 0.0, 0.0 )
    
    do i=1,n2
        s(i) = splines(i,1) &
            + splines(i,2)*(x2(i)-x(i)) &
            + splines(i,3)*(x2(i)-x(i))**2 &
            + splines(i,4)*(x2(i)-x(i))**3
    enddo
    
    
    write(*,*) x
    write(*,*) y
    write(*,*) s

    ! write value to file
    ! -------------------

    fName = 'test_spline.nc'
    call check ( nf90_create ( fName, nf90_clobber, nc_id ) )
    call check ( nf90_def_dim ( nc_id, "n", n, n_id ) )
    call check ( nf90_def_dim ( nc_id, "n2", n2, n2_id ) )

    call check ( nf90_def_var ( nc_id, "x", nf90_real, &
        (/n_id/), x_id ) ) 
    call check ( nf90_def_var ( nc_id, "y", nf90_real, &
        (/n_id/), y_id ) ) 
    call check ( nf90_def_var ( nc_id, "s", nf90_real, &
        (/n2_id/), s_id ) ) 
    call check ( nf90_def_var ( nc_id, "x2", nf90_real, &
        (/n2_id/), x2_id ) ) 
 
    call check ( nf90_enddef ( nc_id ) )

    call check ( nf90_put_var ( nc_id, x_id, x ) )
    call check ( nf90_put_var ( nc_id, y_id, y ) )
    call check ( nf90_put_var ( nc_id, s_id, s ) )
    call check ( nf90_put_var ( nc_id, x2_id, x2 ) )

    call check ( nf90_close ( nc_id ) )
 
end program test_spline
