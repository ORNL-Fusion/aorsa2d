module write_data

use netcdf

contains

    subroutine write_runData ( fName )
 
        use aorsa2din_mod, &
        only: nPtsX, nPtsY, nModesX, nModesY, nSpec
        use bField
        use antenna, &
        only: xjx, xjy, xjz
        use grid, &
        only: capR, y, xkxsav, xkysav, xx, yy
        use profiles
        use constants

        implicit none

        character(len=*), intent(in) :: fName 

        integer :: nc_id, nX_id, nY_id, nc_stat
        integer :: nModesX_id, nModesY_id, nSpec_id
        integer :: &
            x_id, y_id, &
            bx_id, by_id, bz_id, bmod_id, &
            jy_re_id, jy_im_id, kx_id, ky_id, &
            dens_id, temp_id
        integer :: xx_re_id, yy_re_id, xx_im_id, yy_im_id

        call check ( nf90_create ( fName, nf90_clobber, nc_id ) )
        call check ( nf90_def_dim ( nc_id, "nPtsX", nPtsX, nX_id ) )
        call check ( nf90_def_dim ( nc_id, "nPtsY", nPtsY, nY_id ) )
        call check ( nf90_def_dim ( nc_id, "nModesX", nModesX, nModesX_id ) )
        call check ( nf90_def_dim ( nc_id, "nModesY", nModesY, nModesY_id ) )
        call check ( nf90_def_dim ( nc_id, "nSpec", nSpec, nSpec_id ) )

        call check ( nf90_def_var ( nc_id, "capR", NF90_REAL, &
            (/nX_id/), x_id ) ) 
        call check ( nf90_def_var ( nc_id, "y", NF90_REAL, &
            (/nY_id/), y_id ) ) 
 
        call check ( nf90_def_var ( nc_id, "bxn", NF90_REAL, &
            (/nX_id,nY_id/), bx_id ) ) 
        call check ( nf90_def_var ( nc_id, "byn", NF90_REAL, &
            (/nX_id,nY_id/), by_id ) ) 
        call check ( nf90_def_var ( nc_id, "bzn", NF90_REAL, &
            (/nX_id,nY_id/), bz_id ) ) 
        call check ( nf90_def_var ( nc_id, "bmod", NF90_REAL, &
            (/nX_id,nY_id/), bmod_id ) ) 
 
        call check ( nf90_def_var ( nc_id, "xjy_re", NF90_REAL, &
            (/nX_id,nY_id/), jy_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "xjy_im", NF90_REAL, &
            (/nX_id,nY_id/), jy_im_id ) ) 

        call check ( nf90_def_var ( nc_id, "xkxsav", NF90_REAL, &
            (/nModesX_id/), kx_id ) ) 
        call check ( nf90_def_var ( nc_id, "xkysav", NF90_REAL, &
            (/nModesY_id/), ky_id ) ) 
 
        call check ( nf90_def_var ( nc_id, "xx_re", NF90_REAL, &
            (/nModesX_id,nX_id/), xx_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "yy_re", NF90_REAL, &
            (/nModesY_id,nY_id/), yy_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "xx_im", NF90_REAL, &
            (/nModesX_id,nX_id/), xx_im_id ) ) 
        call check ( nf90_def_var ( nc_id, "yy_im", NF90_REAL, &
            (/nModesY_id,nY_id/), yy_im_id ) ) 
 

        nc_stat = nf90_def_var ( nc_id, "densitySpec", NF90_REAL, (/nX_id,nY_id,nSpec_id/), dens_id ) 
        nc_stat = nf90_def_var ( nc_id, "tempSpec", NF90_REAL, (/nX_id,nY_id,nSpec_id/), temp_id ) 

        call check ( nf90_enddef ( nc_id ) )
        
        call check ( nf90_put_var ( nc_id, x_id, capR ) )
        call check ( nf90_put_var ( nc_id, y_id, y ) )

        call check ( nf90_put_var ( nc_id, bx_id, bxn ) )
        call check ( nf90_put_var ( nc_id, by_id, byn ) )
        call check ( nf90_put_var ( nc_id, bz_id, bzn ) )
        call check ( nf90_put_var ( nc_id, bmod_id, bmod ) )
        call check ( nf90_put_var ( nc_id, jy_re_id, real(xjy) ) )
        call check ( nf90_put_var ( nc_id, jy_im_id, aimag(xjy) ) )
        call check ( nf90_put_var ( nc_id, kx_id, xkxsav ) )
        call check ( nf90_put_var ( nc_id, ky_id, xkysav ) )
        nc_stat = nf90_put_var ( nc_id, dens_id, densitySpec )
        nc_stat = nf90_put_var ( nc_id, temp_id, ktSpec / q )

        call check ( nf90_put_var ( nc_id, xx_re_id, real ( xx ) ) )
        call check ( nf90_put_var ( nc_id, yy_re_id, real ( yy ) ) )
        call check ( nf90_put_var ( nc_id, xx_im_id, aimag ( xx ) ) )
        call check ( nf90_put_var ( nc_id, yy_im_id, aimag ( yy ) ) )


        call check ( nf90_close ( nc_id ) )

    end subroutine write_runData


    subroutine write_solution ( fName )

        use aorsa2din_mod, &
        only: nPtsX, nPtsY, nModesX, nModesY
        use solve, &
        only: ealpha, ebeta, eB, ealphak, ebetak, eBk

        implicit none

        character(len=*), intent(in) :: fName 

        integer :: nc_id, nX_id, nY_id, nModesX_id, nModesY_id
        integer :: &
            e1_re_id, e1_im_id, &
            e2_re_id, e2_im_id, &
            e3_re_id, e3_im_id
        integer :: &
            e1k_re_id, e1k_im_id, &
            e2k_re_id, e2k_im_id, &
            e3k_re_id, e3k_im_id

        call check ( nf90_create ( fName, nf90_clobber, nc_id ) )
        call check ( nf90_def_dim ( nc_id, "nX", nPtsX, nX_id ) )
        call check ( nf90_def_dim ( nc_id, "nY", nPtsY, nY_id ) )
        call check ( nf90_def_dim ( nc_id, "nModesX", nModesX, nModesX_id ) )
        call check ( nf90_def_dim ( nc_id, "nModesY", nModesY, nModesY_id ) )


        call check ( nf90_def_var ( nc_id, "ealpha_re", NF90_REAL, &
            (/nX_id,nY_id/), e1_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "ealpha_im", NF90_REAL, &
            (/nX_id,nY_id/), e1_im_id ) ) 
        call check ( nf90_def_var ( nc_id, "ebeta_re", NF90_REAL, &
            (/nX_id,nY_id/), e2_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "ebeta_im", NF90_REAL, &
            (/nX_id,nY_id/), e2_im_id ) ) 
        call check ( nf90_def_var ( nc_id, "eB_re", NF90_REAL, &
            (/nX_id,nY_id/), e3_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "eB_im", NF90_REAL, &
            (/nX_id,nY_id/), e3_im_id ) ) 
 
        call check ( nf90_def_var ( nc_id, "ealphak_re", NF90_REAL, &
            (/nModesX_id,nModesY_id/), e1k_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "ealphak_im", NF90_REAL, &
            (/nModesX_id,nModesY_id/), e1k_im_id ) ) 
        call check ( nf90_def_var ( nc_id, "ebetak_re", NF90_REAL, &
            (/nModesX_id,nModesY_id/), e2k_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "ebetak_im", NF90_REAL, &
            (/nModesX_id,nModesY_id/), e2k_im_id ) ) 
        call check ( nf90_def_var ( nc_id, "eBk_re", NF90_REAL, &
            (/nModesX_id,nModesY_id/), e3k_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "eBk_im", NF90_REAL, &
            (/nModesX_id,nModesY_id/), e3k_im_id ) ) 
 
        call check ( nf90_enddef ( nc_id ) )

        call check ( nf90_put_var ( nc_id, e1_re_id, real ( ealpha ) ) )
        call check ( nf90_put_var ( nc_id, e1_im_id, aimag ( ealpha ) ) )
        call check ( nf90_put_var ( nc_id, e2_re_id, real ( ebeta ) ) )
        call check ( nf90_put_var ( nc_id, e2_im_id, aimag ( ebeta ) ) )
        call check ( nf90_put_var ( nc_id, e3_re_id, real ( eB ) ) )
        call check ( nf90_put_var ( nc_id, e3_im_id, aimag ( eB ) ) )

        call check ( nf90_put_var ( nc_id, e1k_re_id, real ( ealphak ) ) )
        call check ( nf90_put_var ( nc_id, e1k_im_id, aimag ( ealphak ) ) )
        call check ( nf90_put_var ( nc_id, e2k_re_id, real ( ebetak ) ) )
        call check ( nf90_put_var ( nc_id, e2k_im_id, aimag ( ebetak ) ) )
        call check ( nf90_put_var ( nc_id, e3k_re_id, real ( eBk ) ) )
        call check ( nf90_put_var ( nc_id, e3k_im_id, aimag ( eBk ) ) )

        call check ( nf90_close ( nc_id ) )

    end subroutine write_solution 

    subroutine write_amat ( fName )

        use mat_fill 
        use grid, &
        only: xx, yy
        use aorsa2din_mod, &
        only: nPtsX, nPtsY, nModesX, nModesY
 
        implicit none

        character(len=*), intent(in) :: fName 

        integer :: nc_id, nX_id, nY_id, scalar_id
        integer :: nX, nY
        integer :: amat_re_id, amat_im_id
        integer :: nPtsX_id, nPtsY_id, nModesX_id, nModesY_id
        integer :: nModesX_dim_id, nModesY_dim_id, nPtsX_dim_id, nPtsY_dim_id
        integer :: xx_re_id, yy_re_id, xx_im_id, yy_im_id

        nX  = size ( aMat, 1 )
        nY  = size ( aMat, 2 )

        call check ( nf90_create ( fName, nf90_clobber, nc_id ) )
        call check ( nf90_def_dim ( nc_id, "nX", nX, nX_id ) )
        call check ( nf90_def_dim ( nc_id, "nY", nY, nY_id ) )
        call check ( nf90_def_dim ( nc_id, "scalar", 1, scalar_id ) )
        call check ( nf90_def_dim ( nc_id, "nModesX", nModesX, nModesX_dim_id ) )
        call check ( nf90_def_dim ( nc_id, "nModesY", nModesY, nModesY_dim_id ) )
        call check ( nf90_def_dim ( nc_id, "nPtsX", nPtsX, nPtsX_dim_id ) )
        call check ( nf90_def_dim ( nc_id, "nPtsY", nPtsY, nPtsY_dim_id ) )

        call check ( nf90_def_var ( nc_id, "amat_re", NF90_REAL, &
            (/nX_id,nY_id/), amat_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "amat_im", NF90_REAL, &
            (/nX_id,nY_id/), amat_im_id ) ) 
        call check ( nf90_def_var ( nc_id, "nPtsX", NF90_INT, &
            scalar_id, nPtsX_id ) ) 
        call check ( nf90_def_var ( nc_id, "nPtsY", NF90_INT, &
            scalar_id, nPtsY_id ) ) 
        call check ( nf90_def_var ( nc_id, "nModesX", NF90_INT, &
            scalar_id, nModesX_id ) ) 
        call check ( nf90_def_var ( nc_id, "nModesY", NF90_INT, &
            scalar_id, nModesY_id ) ) 
        call check ( nf90_def_var ( nc_id, "xx_re", NF90_REAL, &
            (/nModesX_dim_id,nPtsX_dim_id/), xx_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "yy_re", NF90_REAL, &
            (/nModesY_dim_id,nPtsY_dim_id/), yy_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "xx_im", NF90_REAL, &
            (/nModesX_dim_id,nPtsX_dim_id/), xx_im_id ) ) 
        call check ( nf90_def_var ( nc_id, "yy_im", NF90_REAL, &
            (/nModesY_dim_id,nPtsY_dim_id/), yy_im_id ) ) 
 
        call check ( nf90_enddef ( nc_id ) )

        call check ( nf90_put_var ( nc_id, amat_re_id, real ( aMat ) ) )
        call check ( nf90_put_var ( nc_id, amat_im_id, aimag ( aMat ) ) )
        call check ( nf90_put_var ( nc_id, nPtsX_id, nPtsX ) )
        call check ( nf90_put_var ( nc_id, nPtsY_id, nPtsY ) )
        call check ( nf90_put_var ( nc_id, nModesX_id, nModesX ) )
        call check ( nf90_put_var ( nc_id, nModesY_id, nModesY ) )
        call check ( nf90_put_var ( nc_id, xx_re_id, real ( xx ) ) )
        call check ( nf90_put_var ( nc_id, yy_re_id, real ( yy ) ) )
        call check ( nf90_put_var ( nc_id, xx_im_id, aimag ( xx ) ) )
        call check ( nf90_put_var ( nc_id, yy_im_id, aimag ( yy ) ) )


        call check ( nf90_close ( nc_id ) )

    end subroutine write_amat


    subroutine check ( status )

            integer, intent(in) :: status
          
            if(status /= nf90_noerr) then 
                print *, trim(nf90_strerror(status))
                stop "Stopped"
            end if
    
    end subroutine check

end module write_data
