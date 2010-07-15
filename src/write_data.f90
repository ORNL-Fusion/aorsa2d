module write_data

use netcdf
use check_mod

contains

    subroutine write_runData ( g )
 
        use aorsa2din_mod, &
        only: nSpec
        use bField
        use grid
        use constants
        use rotation

        implicit none

        type(gridBlock), intent(in) :: g
        character(len=100) :: fName 

        integer :: nc_id, nX_id, nY_id, nc_stat
        integer :: nModesX_id, nModesY_id, nSpec_id
        integer :: &
            x_id, y_id, &
            bx_id, by_id, bz_id, bmod_id, &
            jy_re_id, jy_im_id, kx_id, ky_id, &
            dens_id, temp_id, omgc_id, omgp2_id
        integer :: &
            xx_re_id, yy_re_id, xx_im_id, yy_im_id, &
            drBfn_re_id, dzBfn_re_id, drbFn_im_id, dzbFn_im_id

        integer :: drUrrid, drUrtid, drUrzid
        integer :: drUtrid, drUttid, drUtzid
        integer :: drUzrid, drUztid, drUzzid

        integer :: dzUrrid, dzUrtid, dzUrzid
        integer :: dzUtrid, dzUttid, dzUtzid
        integer :: dzUzrid, dzUztid, dzUzzid

        integer :: Urrid, Urtid, Urzid, &
            Utrid, Uttid, Utzid, &
            Uzrid, Uztid, Uzzid

        fName = 'runData'//g%fNumber//'.nc'

        call check ( nf90_create ( fName, nf90_clobber, nc_id ) )
        call check ( nf90_def_dim ( nc_id, "nPtsX", g%nR, nX_id ) )
        call check ( nf90_def_dim ( nc_id, "nPtsY", g%nZ, nY_id ) )
        call check ( nf90_def_dim ( nc_id, "nModesX", g%nModesR, nModesX_id ) )
        call check ( nf90_def_dim ( nc_id, "nModesY", g%nModesZ, nModesY_id ) )
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

        !call check ( nf90_def_var ( nc_id, "kxsav", NF90_REAL, &
        !    (/nModesX_id/), kx_id ) ) 
        !call check ( nf90_def_var ( nc_id, "kysav", NF90_REAL, &
        !    (/nModesY_id/), ky_id ) ) 
 
        call check ( nf90_def_var ( nc_id, "xx_re", NF90_REAL, &
            (/nModesX_id,nX_id/), xx_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "yy_re", NF90_REAL, &
            (/nModesY_id,nY_id/), yy_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "xx_im", NF90_REAL, &
            (/nModesX_id,nX_id/), xx_im_id ) ) 
        call check ( nf90_def_var ( nc_id, "yy_im", NF90_REAL, &
            (/nModesY_id,nY_id/), yy_im_id ) ) 
 
        call check ( nf90_def_var ( nc_id, "drBfn_re", NF90_REAL, &
            (/nModesX_id,nX_id/), drBfn_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "dzBfn_re", NF90_REAL, &
            (/nModesY_id,nY_id/), dzBfn_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "drBfn_im", NF90_REAL, &
            (/nModesX_id,nX_id/), drBfn_im_id ) ) 
        call check ( nf90_def_var ( nc_id, "dzBfn_im", NF90_REAL, &
            (/nModesY_id,nY_id/), dzBfn_im_id ) ) 

        nc_stat = nf90_def_var ( nc_id, "densitySpec", NF90_REAL, (/nX_id,nY_id,nSpec_id/), dens_id ) 
        nc_stat = nf90_def_var ( nc_id, "tempSpec", NF90_REAL, (/nX_id,nY_id,nSpec_id/), temp_id ) 
        nc_stat = nf90_def_var ( nc_id, "omgc", NF90_REAL, (/nX_id,nY_id,nSpec_id/), omgc_id ) 
        nc_stat = nf90_def_var ( nc_id, "omgp2", NF90_REAL, (/nX_id,nY_id,nSpec_id/), omgp2_id ) 

        nc_stat = nf90_def_var ( nc_id, "drUrr", NF90_REAL, (/nX_id,nY_id/), drUrrid ) 
        nc_stat = nf90_def_var ( nc_id, "drUrt", NF90_REAL, (/nX_id,nY_id/), drUrtid ) 
        nc_stat = nf90_def_var ( nc_id, "drUrz", NF90_REAL, (/nX_id,nY_id/), drUrzid ) 

        nc_stat = nf90_def_var ( nc_id, "drUtr", NF90_REAL, (/nX_id,nY_id/), drUtrid ) 
        nc_stat = nf90_def_var ( nc_id, "drUtt", NF90_REAL, (/nX_id,nY_id/), drUttid ) 
        nc_stat = nf90_def_var ( nc_id, "drUtz", NF90_REAL, (/nX_id,nY_id/), drUtzid ) 

        nc_stat = nf90_def_var ( nc_id, "drUzr", NF90_REAL, (/nX_id,nY_id/), drUzrid ) 
        nc_stat = nf90_def_var ( nc_id, "drUzt", NF90_REAL, (/nX_id,nY_id/), drUztid ) 
        nc_stat = nf90_def_var ( nc_id, "drUzz", NF90_REAL, (/nX_id,nY_id/), drUzzid ) 

        nc_stat = nf90_def_var ( nc_id, "dzUrr", NF90_REAL, (/nX_id,nY_id/), dzUrrid ) 
        nc_stat = nf90_def_var ( nc_id, "dzUrt", NF90_REAL, (/nX_id,nY_id/), dzUrtid ) 
        nc_stat = nf90_def_var ( nc_id, "dzUrz", NF90_REAL, (/nX_id,nY_id/), dzUrzid ) 

        nc_stat = nf90_def_var ( nc_id, "dzUtr", NF90_REAL, (/nX_id,nY_id/), dzUtrid ) 
        nc_stat = nf90_def_var ( nc_id, "dzUtt", NF90_REAL, (/nX_id,nY_id/), dzUttid ) 
        nc_stat = nf90_def_var ( nc_id, "dzUtz", NF90_REAL, (/nX_id,nY_id/), dzUtzid ) 

        nc_stat = nf90_def_var ( nc_id, "dzUzr", NF90_REAL, (/nX_id,nY_id/), dzUzrid ) 
        nc_stat = nf90_def_var ( nc_id, "dzUzt", NF90_REAL, (/nX_id,nY_id/), dzUztid ) 
        nc_stat = nf90_def_var ( nc_id, "dzUzz", NF90_REAL, (/nX_id,nY_id/), dzUzzid ) 

        nc_stat = nf90_def_var ( nc_id, "Urr", NF90_REAL, (/nX_id,nY_id/), Urrid ) 
        nc_stat = nf90_def_var ( nc_id, "Urt", NF90_REAL, (/nX_id,nY_id/), Urtid ) 
        nc_stat = nf90_def_var ( nc_id, "Urz", NF90_REAL, (/nX_id,nY_id/), Urzid ) 
        
        nc_stat = nf90_def_var ( nc_id, "Utr", NF90_REAL, (/nX_id,nY_id/), Utrid ) 
        nc_stat = nf90_def_var ( nc_id, "Utt", NF90_REAL, (/nX_id,nY_id/), Uttid ) 
        nc_stat = nf90_def_var ( nc_id, "Utz", NF90_REAL, (/nX_id,nY_id/), Utzid ) 

        nc_stat = nf90_def_var ( nc_id, "Uzr", NF90_REAL, (/nX_id,nY_id/), Uzrid ) 
        nc_stat = nf90_def_var ( nc_id, "Uzt", NF90_REAL, (/nX_id,nY_id/), Uztid ) 
        nc_stat = nf90_def_var ( nc_id, "Uzz", NF90_REAL, (/nX_id,nY_id/), Uzzid ) 

        call check ( nf90_enddef ( nc_id ) )
        
        call check ( nf90_put_var ( nc_id, x_id, g%R ) )
        call check ( nf90_put_var ( nc_id, y_id, g%Z ) )

        call check ( nf90_put_var ( nc_id, bx_id, g%bR_unit ) )
        call check ( nf90_put_var ( nc_id, by_id, g%bT_unit ) )
        call check ( nf90_put_var ( nc_id, bz_id, g%bZ_unit ) )
        call check ( nf90_put_var ( nc_id, bmod_id, g%bMag ) )
        call check ( nf90_put_var ( nc_id, jy_re_id, real(g%jZ) ) )
        call check ( nf90_put_var ( nc_id, jy_im_id, aimag(g%jZ) ) )
        !call check ( nf90_put_var ( nc_id, kx_id, kxsav ) )
        !call check ( nf90_put_var ( nc_id, ky_id, kysav ) )
        nc_stat = nf90_put_var ( nc_id, dens_id, g%densitySpec )
        nc_stat = nf90_put_var ( nc_id, temp_id, g%ktSpec / q )
        nc_stat = nf90_put_var ( nc_id, omgc_id, g%omgc )
        nc_stat = nf90_put_var ( nc_id, omgp2_id, g%omgp2 )

        call check ( nf90_put_var ( nc_id, xx_re_id, real ( g%xx ) ) )
        call check ( nf90_put_var ( nc_id, yy_re_id, real ( g%yy ) ) )
        call check ( nf90_put_var ( nc_id, xx_im_id, aimag ( g%xx ) ) )
        call check ( nf90_put_var ( nc_id, yy_im_id, aimag ( g%yy ) ) )

        call check ( nf90_put_var ( nc_id, drBfn_re_id, real ( g%dRbFn_bFn ) ) )
        call check ( nf90_put_var ( nc_id, dzBfn_re_id, real ( g%dZbFn_bFn ) ) )
        call check ( nf90_put_var ( nc_id, drBfn_im_id, aimag ( g%dRbFn_bFn ) ) )
        call check ( nf90_put_var ( nc_id, dzBfn_im_id, aimag ( g%dZbFn_bFn ) ) )

        call check ( nf90_put_var ( nc_id, drUrrid, g%drUrr ) )
        call check ( nf90_put_var ( nc_id, drUrtid, g%drUrt ) )
        call check ( nf90_put_var ( nc_id, drUrzid, g%drUrz ) )

        call check ( nf90_put_var ( nc_id, drUtrid, g%drUtr ) )
        call check ( nf90_put_var ( nc_id, drUttid, g%drUtt ) )
        call check ( nf90_put_var ( nc_id, drUtzid, g%drUtz ) )

        call check ( nf90_put_var ( nc_id, drUzrid, g%drUzr ) )
        call check ( nf90_put_var ( nc_id, drUztid, g%drUzt ) )
        call check ( nf90_put_var ( nc_id, drUzzid, g%drUzz ) )

        call check ( nf90_put_var ( nc_id, dzUrrid, g%dzUrr ) )
        call check ( nf90_put_var ( nc_id, dzUrtid, g%dzUrt ) )
        call check ( nf90_put_var ( nc_id, dzUrzid, g%dzUrz ) )

        call check ( nf90_put_var ( nc_id, dzUtrid, g%dzUtr ) )
        call check ( nf90_put_var ( nc_id, dzUttid, g%dzUtt ) )
        call check ( nf90_put_var ( nc_id, dzUtzid, g%dzUtz ) )

        call check ( nf90_put_var ( nc_id, dzUzrid, g%dzUzr ) )
        call check ( nf90_put_var ( nc_id, dzUztid, g%dzUzt ) )
        call check ( nf90_put_var ( nc_id, dzUzzid, g%dzUzz ) )

        call check ( nf90_put_var ( nc_id, Urrid, g%Urr ) )
        call check ( nf90_put_var ( nc_id, Urtid, g%Urt ) )
        call check ( nf90_put_var ( nc_id, Urzid, g%Urz ) )

        call check ( nf90_put_var ( nc_id, Utrid, g%Utr ) )
        call check ( nf90_put_var ( nc_id, Uttid, g%Utt ) )
        call check ( nf90_put_var ( nc_id, Utzid, g%Utz ) )

        call check ( nf90_put_var ( nc_id, Uzrid, g%Uzr ) )
        call check ( nf90_put_var ( nc_id, Uztid, g%Uzt ) )
        call check ( nf90_put_var ( nc_id, Uzzid, g%Uzz ) )

        call check ( nf90_close ( nc_id ) )

    end subroutine write_runData


    subroutine write_solution ( g )

        use grid
        use aorsa2din_mod, &
        only: nSpec

        implicit none

        type(gridBlock), intent(in) :: g
        character(len=100) :: fName 

        integer :: nc_id, nX_id, nY_id, nModesX_id, nModesY_id, nSpec_id
        integer :: &
            e1_re_id, e1_im_id, &
            e2_re_id, e2_im_id, &
            e3_re_id, e3_im_id
        integer :: &
            e1k_re_id, e1k_im_id, &
            e2k_re_id, e2k_im_id, &
            e3k_re_id, e3k_im_id
        integer :: &
            er_re_id, er_im_id, &
            et_re_id, et_im_id, &
            ez_re_id, ez_im_id
        integer :: &
            j1_re_id, j1_im_id, &
            j2_re_id, j2_im_id, &
            j3_re_id, j3_im_id

        integer :: jouleHeating_id

        fName = 'solution'//g%fNumber//'.nc'

        call check ( nf90_create ( fName, nf90_clobber, nc_id ) )
        call check ( nf90_def_dim ( nc_id, "nX", g%nR, nX_id ) )
        call check ( nf90_def_dim ( nc_id, "nY", g%nZ, nY_id ) )
        call check ( nf90_def_dim ( nc_id, "nModesX", g%nModesR, nModesX_id ) )
        call check ( nf90_def_dim ( nc_id, "nModesY", g%nModesZ, nModesY_id ) )
        call check ( nf90_def_dim ( nc_id, "nSpec", nSpec, nSpec_id ) )

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

        call check ( nf90_def_var ( nc_id, "er_re", NF90_REAL, &
            (/nX_id,nY_id/), er_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "er_im", NF90_REAL, &
            (/nX_id,nY_id/), er_im_id ) ) 
        call check ( nf90_def_var ( nc_id, "et_re", NF90_REAL, &
            (/nX_id,nY_id/), et_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "et_im", NF90_REAL, &
            (/nX_id,nY_id/), et_im_id ) ) 
        call check ( nf90_def_var ( nc_id, "ez_re", NF90_REAL, &
            (/nX_id,nY_id/), ez_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "ez_im", NF90_REAL, &
            (/nX_id,nY_id/), ez_im_id ) ) 

        call check ( nf90_def_var ( nc_id, "jalpha_re", NF90_REAL, &
            (/nX_id,nY_id,nSpec_id/), j1_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "jalpha_im", NF90_REAL, &
            (/nX_id,nY_id,nSpec_id/), j1_im_id ) ) 
        call check ( nf90_def_var ( nc_id, "jbeta_re", NF90_REAL, &
            (/nX_id,nY_id,nSpec_id/), j2_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "jbeta_im", NF90_REAL, &
            (/nX_id,nY_id,nSpec_id/), j2_im_id ) ) 
        call check ( nf90_def_var ( nc_id, "jB_re", NF90_REAL, &
            (/nX_id,nY_id,nSpec_id/), j3_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "jB_im", NF90_REAL, &
            (/nX_id,nY_id,nSpec_id/), j3_im_id ) ) 
 
        call check ( nf90_def_var ( nc_id, "jouleHeating", NF90_REAL, &
            (/nX_id,nY_id,nSpec_id/), jouleHeating_id ) ) 
 

        call check ( nf90_enddef ( nc_id ) )

        call check ( nf90_put_var ( nc_id, e1_re_id, realpart ( g%ealpha ) ) )
        call check ( nf90_put_var ( nc_id, e1_im_id, imagpart ( g%ealpha ) ) )
        call check ( nf90_put_var ( nc_id, e2_re_id, real ( g%ebeta ) ) )
        call check ( nf90_put_var ( nc_id, e2_im_id, aimag ( g%ebeta ) ) )
        call check ( nf90_put_var ( nc_id, e3_re_id, real ( g%eB ) ) )
        call check ( nf90_put_var ( nc_id, e3_im_id, aimag ( g%eB ) ) )

        call check ( nf90_put_var ( nc_id, e1k_re_id, real ( g%ealphak ) ) )
        call check ( nf90_put_var ( nc_id, e1k_im_id, aimag ( g%ealphak ) ) )
        call check ( nf90_put_var ( nc_id, e2k_re_id, real ( g%ebetak ) ) )
        call check ( nf90_put_var ( nc_id, e2k_im_id, aimag ( g%ebetak ) ) )
        call check ( nf90_put_var ( nc_id, e3k_re_id, real ( g%eBk ) ) )
        call check ( nf90_put_var ( nc_id, e3k_im_id, aimag ( g%eBk ) ) )

        call check ( nf90_put_var ( nc_id, er_re_id, real ( g%eR ) ) )
        call check ( nf90_put_var ( nc_id, er_im_id, aimag ( g%eR ) ) )
        call check ( nf90_put_var ( nc_id, et_re_id, real ( g%eTh ) ) )
        call check ( nf90_put_var ( nc_id, et_im_id, aimag ( g%eTh ) ) )
        call check ( nf90_put_var ( nc_id, ez_re_id, real ( g%eZ ) ) )
        call check ( nf90_put_var ( nc_id, ez_im_id, aimag ( g%eZ ) ) )

        call check ( nf90_put_var ( nc_id, j1_re_id, realpart ( g%jalpha ) ) )
        call check ( nf90_put_var ( nc_id, j1_im_id, imagpart ( g%jalpha ) ) )
        call check ( nf90_put_var ( nc_id, j2_re_id, realpart ( g%jbeta ) ) )
        call check ( nf90_put_var ( nc_id, j2_im_id, imagpart ( g%jbeta ) ) )
        call check ( nf90_put_var ( nc_id, j3_re_id, realpart ( g%jB ) ) )
        call check ( nf90_put_var ( nc_id, j3_im_id, imagpart ( g%jB ) ) )

        call check ( nf90_put_var ( nc_id, jouleHeating_id, g%jouleHeating ) )

        call check ( nf90_close ( nc_id ) )

    end subroutine write_solution 

    !subroutine write_amat ( g, fName )

    !    !use mat_fill 
    !    use grid
    !    
    !    implicit none

    !    type(gridBlock), intent(in) :: g
    !    character(len=*), intent(in) :: fName 

    !    integer :: nc_id, nX_id, nY_id, scalar_id
    !    integer :: nX, nY
    !    integer :: amat_re_id, amat_im_id
    !    integer :: nPtsX_id, nPtsY_id, nModesX_id, nModesY_id
    !    integer :: nModesX_dim_id, nModesY_dim_id, nPtsX_dim_id, nPtsY_dim_id
    !    integer :: xx_re_id, yy_re_id, xx_im_id, yy_im_id

    !    nX  = g%nModesR * g%nModesZ * 3 !size ( aMat, 1 )
    !    nY  = g%nR * g%nZ * 3 !size ( aMat, 2 )

    !    call check ( nf90_create ( fName, nf90_clobber, nc_id ) )
    !    call check ( nf90_def_dim ( nc_id, "nX", nX, nX_id ) )
    !    call check ( nf90_def_dim ( nc_id, "nY", nY, nY_id ) )
    !    call check ( nf90_def_dim ( nc_id, "scalar", 1, scalar_id ) )
    !    call check ( nf90_def_dim ( nc_id, "nModesX", g%nModesR, nModesX_dim_id ) )
    !    call check ( nf90_def_dim ( nc_id, "nModesY", g%nModesZ, nModesY_dim_id ) )
    !    call check ( nf90_def_dim ( nc_id, "nPtsX", g%nR, nPtsX_dim_id ) )
    !    call check ( nf90_def_dim ( nc_id, "nPtsY", g%nZ, nPtsY_dim_id ) )

    !    call check ( nf90_def_var ( nc_id, "amat_re", NF90_REAL, &
    !        (/nX_id,nY_id/), amat_re_id ) ) 
    !    call check ( nf90_def_var ( nc_id, "amat_im", NF90_REAL, &
    !        (/nX_id,nY_id/), amat_im_id ) ) 
    !    call check ( nf90_def_var ( nc_id, "nPtsX", NF90_INT, &
    !        scalar_id, nPtsX_id ) ) 
    !    call check ( nf90_def_var ( nc_id, "nPtsY", NF90_INT, &
    !        scalar_id, nPtsY_id ) ) 
    !    call check ( nf90_def_var ( nc_id, "nModesX", NF90_INT, &
    !        scalar_id, nModesX_id ) ) 
    !    call check ( nf90_def_var ( nc_id, "nModesY", NF90_INT, &
    !        scalar_id, nModesY_id ) ) 
    !    call check ( nf90_def_var ( nc_id, "xx_re", NF90_REAL, &
    !        (/nModesX_dim_id,nPtsX_dim_id/), xx_re_id ) ) 
    !    call check ( nf90_def_var ( nc_id, "yy_re", NF90_REAL, &
    !        (/nModesY_dim_id,nPtsY_dim_id/), yy_re_id ) ) 
    !    call check ( nf90_def_var ( nc_id, "xx_im", NF90_REAL, &
    !        (/nModesX_dim_id,nPtsX_dim_id/), xx_im_id ) ) 
    !    call check ( nf90_def_var ( nc_id, "yy_im", NF90_REAL, &
    !        (/nModesY_dim_id,nPtsY_dim_id/), yy_im_id ) ) 
 
    !    call check ( nf90_enddef ( nc_id ) )

    !    call check ( nf90_put_var ( nc_id, amat_re_id, real ( aMat ) ) )
    !    call check ( nf90_put_var ( nc_id, amat_im_id, aimag ( aMat ) ) )
    !    call check ( nf90_put_var ( nc_id, nPtsX_id, g%nR ) )
    !    call check ( nf90_put_var ( nc_id, nPtsY_id, g%nZ ) )
    !    call check ( nf90_put_var ( nc_id, nModesX_id, g%nModesR ) )
    !    call check ( nf90_put_var ( nc_id, nModesY_id, g%nModesZ ) )
    !    call check ( nf90_put_var ( nc_id, xx_re_id, real ( g%xx ) ) )
    !    call check ( nf90_put_var ( nc_id, yy_re_id, real ( g%yy ) ) )
    !    call check ( nf90_put_var ( nc_id, xx_im_id, aimag ( g%xx ) ) )
    !    call check ( nf90_put_var ( nc_id, yy_im_id, aimag ( g%yy ) ) )


    !    call check ( nf90_close ( nc_id ) )

    !end subroutine write_amat


    ! Routines to initialize and incrementally write to the
    ! sigma file
    ! -----------------------------------------------------

    subroutine init_sigma_file &
        ( g, fName, nc_id, sigma_re_id, sigma_im_id, nSpec )
    
        use grid
        
        implicit none

        type(gridBlock), intent(in) :: g
        character(len=*), intent(in) :: fName 
        integer, intent(in) :: nSpec
        integer, intent(out) :: nc_id, sigma_re_id, sigma_im_id

        integer :: nX_id, nY_id, nN_id, nM_id, n3_id, nS_id
        integer :: nX, nY, nN, nM, n3

        nX = g%nR 
        nY = g%nZ 
        nN = g%nModesR
        nM = g%nModesZ
        n3 = 3

        call check ( nf90_create ( fName, nf90_clobber, nc_id ) )

        call check ( nf90_def_dim ( nc_id, "nX", nX, nX_id ) )
        call check ( nf90_def_dim ( nc_id, "nY", nY, nY_id ) )
        call check ( nf90_def_dim ( nc_id, "nN", nN, nN_id ) )
        call check ( nf90_def_dim ( nc_id, "nM", nM, nM_id ) )
        call check ( nf90_def_dim ( nc_id, "n3", n3, n3_id ) )
        call check ( nf90_def_dim ( nc_id, "nS", nSpec, nS_id ) )

        call check ( nf90_def_var ( nc_id, "sigma_re", NF90_REAL, &
            (/nX_id,nY_id,nN_id,nM_id,n3_id,n3_id,nS_id/), sigma_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "sigma_im", NF90_REAL, &
            (/nX_id,nY_id,nN_id,nM_id,n3_id,n3_id,nS_id/), sigma_im_id ) ) 

        call check ( nf90_enddef ( nc_id ) )

    end subroutine init_sigma_file


    subroutine write_sigma_pt ( i,j,sigma, nc_id, sigma_re_id, sigma_im_id, nN, nM, nSpec )

        complex, intent(in) :: sigma(:,:,:,:,:)
        integer, intent(in) :: i, j, nc_id, sigma_re_id, sigma_im_id
        integer, intent(in) :: nN, nM

        call check ( nf90_put_var ( nc_id, sigma_re_id, realpart ( sigma ), &
            start = (/ i, j, 1, 1, 1, 1, 1 /), &
            count = (/ 1, 1, nN, nM, 3, 3, nSpec /) ) )
        call check ( nf90_put_var ( nc_id, sigma_im_id, imagpart ( sigma ), &
            start = (/ i, j, 1, 1, 1, 1, 1 /), &
            count = (/ 1, 1, nN, nM, 3, 3, nSpec /) ) )

    end subroutine write_sigma_pt


    subroutine close_sigma_file ( nc_id )

        integer, intent(in) :: nc_id

        call check ( nf90_close ( nc_id ) )

    end subroutine close_sigma_file


end module write_data
