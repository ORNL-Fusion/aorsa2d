program test_expBesselI

    use bessel_mod
    use constants
    use netcdf
    use check_mod

    implicit none

    integer :: nPts, lMax, lMaxTmp, iErr, i, j, lOut
    integer, allocatable :: l(:)
    complex(kind=dbl), allocatable :: z(:,:), besselI(:), &
        expBesselI(:), expBesselIPrime(:), expBesselIOverGam(:)
    real :: dz, zMax
    complex, allocatable :: res(:,:), expRes(:,:)
    character(len=100) :: fName
    complex(kind=dbl) :: zLarge, zExtraLarge, zExtraSmall

    integer :: nc_id, n_id, res_re_id, res_im_id, &
        z_re_id, z_im_id, expres_re_id, expres_im_id

    lMax = 20
    lOut = 20
    nPts = 201 
    allocate ( z(nPts,nPts), res(nPts,nPts), expRes(nPts,nPts) )
    allocate ( l(-lMax:lMax) )

    zMax = 10.0 
    dz  = zMax * 2.0 / (nPts - 1)
    do i=1,nPts
        do j=1,nPts

            z(i,j)   = complex ( (i-1)*dz - zMax, (j-1)*dz - zMax )

        enddo
    enddo

    ! I_l(z) and exp(-z) * I_l(z) for small z values
    ! ----------------------------------------------

    allocate ( besselI(lMax+1), &
                expBesselI(lMax+1), &
                expBesselIPrime(lMax+1), &
                expBesselIOverGam(lMax+1) )

    do i=1,nPts
        do j=1,nPts

            lMaxTmp = lMax 
            call besic ( z(i,j), lMaxTmp, besselI, iErr)

            if(iErr/=0) then 
                    write(*,*) 'ERROR: iErr'
                    stop
            endif

            res(i,j)    = besselI(lOut+1)
            expRes(i,j) = exp ( -z(i,j) ) * res(i,j)

        enddo
    enddo


    ! I_l(z) and exp(-z) * I_l(z) for large z values
    ! ----------------------------------------------

    zLarge  = complex ( 7.0, 7.0 )
    call besIExp ( zLarge, lMax, expBesselI, expBesselIPrime, expBesselIOverGam )

    write(*,*) zLarge, expBesselI(lOut+1), expBesselIPrime(lOut+1)

    
    ! exp(-z) * I_l(z) for Re(z) < 1e-8
    ! ------------------------------

    zExtraSmall  = complex ( 1d-9, 1d-9 )
    lMaxTmp = lMax
    call bes_expand ( zExtraSmall, lMaxTmp, expBesselI, expBesselIPrime, expBesselIOverGam )

    write(*,*) zExtraSmall, expBesselI(lOut+1), expBesselIPrime(lOut+1), expBesselIOverGam(lOut+1)


    ! exp(-z) * I_l(z) for |z| > 700
    ! ------------------------------

    zExtraLarge  = complex ( 900, 900 )
    call besIExp ( zExtraLarge, lMax, expBesselI, expBesselIPrime, expBesselIOverGam )

    write(*,*) zExtraLarge, expBesselI(lOut+1), expBesselIPrime(lOut+1)




    ! write value to file for plotting in IDL and comparing
    ! to mathematica solutions
    ! -----------------------------------------------------

    fName = 'test_expBesselI.nc'
    call check ( nf90_create ( fName, nf90_clobber, nc_id ) )
    call check ( nf90_def_dim ( nc_id, "nPts", nPts, n_id ) )

    call check ( nf90_def_var ( nc_id, "res_re", nf90_real, &
        (/n_id,n_id/), res_re_id ) ) 
    call check ( nf90_def_var ( nc_id, "res_im", nf90_real, &
        (/n_id,n_id/), res_im_id ) ) 
    call check ( nf90_def_var ( nc_id, "z_re", NF90_REAL, &
        (/n_id,n_id/), z_re_id ) ) 
    call check ( nf90_def_var ( nc_id, "z_im", NF90_REAL, &
        (/n_id,n_id/), z_im_id ) ) 

    call check ( nf90_def_var ( nc_id, "expRes_re", nf90_real, &
        (/n_id,n_id/), expRes_re_id ) ) 
    call check ( nf90_def_var ( nc_id, "expRes_im", nf90_real, &
        (/n_id,n_id/), expRes_im_id ) ) 
 
    call check ( nf90_enddef ( nc_id ) )

    call check ( nf90_put_var ( nc_id, res_re_id, real ( res ) ) )
    call check ( nf90_put_var ( nc_id, res_im_id, aimag ( res ) ) )
    call check ( nf90_put_var ( nc_id, z_re_id, real ( z ) ) )
    call check ( nf90_put_var ( nc_id, z_im_id, aimag ( z ) ) )
    call check ( nf90_put_var ( nc_id, expRes_re_id, real ( expRes ) ) )
    call check ( nf90_put_var ( nc_id, expRes_im_id, aimag ( expRes ) ) )
 
    call check ( nf90_close ( nc_id ) )

 
end program test_expBesselI
