
module write_data

use netcdf
use check_mod

contains

    subroutine write_runData ( g, rid, rhs )
 
        use aorsaNamelist, &
        only: nSpec, nPhi, freqcy
        use bField
        use grid
        use constants
        use rotation
        use parallel

        implicit none

        type(gridBlock), intent(in) :: g
        integer, intent(in) :: rhs

        character(len=100) :: fName 
        character(len=20) :: rid
        character(len=6) :: rhs_string
        character(len=4) :: nPhi_string

        integer :: nc_id, nX_id, nY_id, nc_stat, NRHS_id
        integer :: nModesX_id, nModesY_id, nSpec_id
        integer :: &
            x_id, y_id, &
            brU_id, btU_id, bzU_id, bmod_id, &
            jr_re_id, jr_im_id, jt_re_id, jt_im_id, &
            jz_re_id, jz_im_id, kx_id, ky_id, &
            dens_id, temp_id, omgc_id, omgp2_id, &
            scalar_id, nPhi_id, freq_id, nuOmg_id
        integer :: &
            xx_re_id, yy_re_id, xx_im_id, yy_im_id, &
            drBfn_re_id, dzBfn_re_id, drbFn_im_id, dzbFn_im_id, &
            LimMask_id

        integer :: drUrrid, drUrtid, drUrzid
        integer :: drUtrid, drUttid, drUtzid
        integer :: drUzrid, drUztid, drUzzid

        integer :: dzUrrid, dzUrtid, dzUrzid
        integer :: dzUtrid, dzUttid, dzUtzid
        integer :: dzUzrid, dzUztid, dzUzzid

        integer :: Urrid, Urtid, Urzid, &
            Utrid, Uttid, Utzid, &
            Uzrid, Uztid, Uzzid

        real, allocatable :: RealTmp(:,:),RealTmp3(:,:,:)
        integer, allocatable :: IntTmp(:,:)
        complex, allocatable :: ComplexTmp2(:,:)
        integer :: p,i,j,s
        integer, allocatable :: Cnt(:,:)

        write(nPhi_string,'(sp,i4.3)'), int(nPhi)
        write(rhs_string,'(i6.6)'), rhs

        fName = trim(rid)//'output/runData_'//g%fNumber//'_'//nPhi_string//'_'//rhs_string//'.nc'

        if(iAm==0)then
            call check ( nf90_create ( fName, nf90_clobber, nc_id ) )
            call check ( nf90_def_dim ( nc_id, "nPtsX", g%nR, nX_id ) )
            call check ( nf90_def_dim ( nc_id, "nPtsY", g%nZ, nY_id ) )
            call check ( nf90_def_dim ( nc_id, "nModesX", g%nModesR, nModesX_id ) )
            call check ( nf90_def_dim ( nc_id, "nModesY", g%nModesZ, nModesY_id ) )
            call check ( nf90_def_dim ( nc_id, "nSpec", nSpec, nSpec_id ) )
            call check ( nf90_def_dim ( nc_id, "scalar", 1, scalar_id ) )
            call check ( nf90_def_dim ( nc_id, "NRHS", NRHS, NRHS_id ) )

            call check ( nf90_def_var ( nc_id, "nPhi", NF90_INT, &
                scalar_id, nPhi_id ) )
            call check ( nf90_def_var ( nc_id, "freq", NF90_REAL, &
                scalar_id, freq_id ) )

            call check ( nf90_def_var ( nc_id, "capR", NF90_REAL, &
                (/nX_id/), x_id ) ) 
            call check ( nf90_def_var ( nc_id, "y", NF90_REAL, &
                (/nY_id/), y_id ) ) 

            call check ( nf90_def_var ( nc_id, "brU", NF90_REAL, &
                (/nX_id,nY_id/), brU_id ) ) 
            call check ( nf90_def_var ( nc_id, "btU", NF90_REAL, &
                (/nX_id,nY_id/), btU_id ) ) 
            call check ( nf90_def_var ( nc_id, "bzU", NF90_REAL, &
                (/nX_id,nY_id/), bzU_id ) ) 
            call check ( nf90_def_var ( nc_id, "bmod", NF90_REAL, &
                (/nX_id,nY_id/), bmod_id ) ) 

            call check ( nf90_def_var ( nc_id, "jr_re", NF90_REAL, &
                (/nX_id,nY_id/), jr_re_id ) ) 
            call check ( nf90_def_var ( nc_id, "jr_im", NF90_REAL, &
                (/nX_id,nY_id/), jr_im_id ) ) 
            call check ( nf90_def_var ( nc_id, "jt_re", NF90_REAL, &
                (/nX_id,nY_id/), jt_re_id ) ) 
            call check ( nf90_def_var ( nc_id, "jt_im", NF90_REAL, &
                (/nX_id,nY_id/), jt_im_id ) ) 
            call check ( nf90_def_var ( nc_id, "jz_re", NF90_REAL, &
                (/nX_id,nY_id/), jz_re_id ) ) 
            call check ( nf90_def_var ( nc_id, "jz_im", NF90_REAL, &
                (/nX_id,nY_id/), jz_im_id ) ) 


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
            nc_stat = nf90_def_var ( nc_id, "nuOmg", NF90_REAL, (/nX_id,nY_id,nSpec_id/), nuOmg_id ) 
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

            nc_stat = nf90_def_var ( nc_id, "LimMask", NF90_INT, (/nX_id,nY_id/), LimMask_id ) 

            call check ( nf90_enddef ( nc_id ) )
        endif

        if(iAm==0)then
            nc_stat = nf90_put_var ( nc_id, nPhi_id, nPhi ) 
            nc_stat = nf90_put_var ( nc_id, freq_id, freqcy ) 
            nc_stat = nf90_put_var ( nc_id, x_id, g%R ) 
            nc_stat = nf90_put_var ( nc_id, y_id, g%Z )    
        endif

        allocate(IntTmp(g%nR,g%nZ)) 
        allocate(Cnt(g%nR,g%nZ))
        IntTmp = 0
        Cnt = 0

        do p=1,size(g%pt)
            i = g%pt(p)%i
            j = g%pt(p)%j
            IntTmp(i,j) = g%isMetal(p)
            Cnt(i,j) = Cnt(i,j)+1
        enddo
        
#ifdef par
        call iGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, IntTmp, g%nR, -1, -1 )
        call iGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, Cnt, g%nR, -1, -1 )
#endif
        IntTmp = IntTmp/Cnt
        IntTmp = abs(IntTmp-1)

        if(iAm==0)then

            nc_stat = nf90_put_var ( nc_id, LimMask_id, IntTmp )
            if(nc_stat.ne.0)then
                    write(*,*) 'ERROR: nc_stat: ', nc_stat
                    stop
            endif
        endif   

        deallocate(IntTmp)

        allocate(RealTmp(g%nR,g%nZ))

        !bR_unit 
        RealTmp = 0 
        Cnt = 0 
        do p=1,size(g%pt)
            i = g%pt(p)%i
            j = g%pt(p)%j
            RealTmp(i,j) = g%bR_unit(p)
            Cnt(i,j) = Cnt(i,j)+1
        enddo
#ifdef par
        call sGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, RealTmp, g%nR, -1, -1 )
        call sGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, Cnt, g%nR, -1, -1 )
#endif
        if(iAm==0)nc_stat = nf90_put_var ( nc_id, brU_id, RealTmp/Cnt )
        RealTmp = 0
        Cnt = 0

        !bT_unit 
        RealTmp = 0 
        Cnt = 0 
        do p=1,size(g%pt)
            i = g%pt(p)%i
            j = g%pt(p)%j
            RealTmp(i,j) = g%bT_unit(p)
            Cnt(i,j) = Cnt(i,j)+1
        enddo
#ifdef par
        call sGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, RealTmp, g%nR, -1, -1 )
        call sGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, Cnt, g%nR, -1, -1 )
#endif
        if(iAm==0)nc_stat = nf90_put_var ( nc_id, btU_id, RealTmp/Cnt )
        RealTmp = 0
        Cnt = 0

        !bZ_unit 
        RealTmp = 0 
        Cnt = 0 
        do p=1,size(g%pt)
            i = g%pt(p)%i
            j = g%pt(p)%j
            RealTmp(i,j) = g%bZ_unit(p)
            Cnt(i,j) = Cnt(i,j)+1
        enddo
#ifdef par
        call sGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, RealTmp, g%nR, -1, -1 )
        call sGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, Cnt, g%nR, -1, -1 )
#endif
        if(iAm==0)nc_stat = nf90_put_var ( nc_id, bzU_id, RealTmp/Cnt )
        RealTmp = 0
        Cnt = 0

        !bMag
        RealTmp = 0 
        Cnt = 0 
        do p=1,size(g%pt)
            i = g%pt(p)%i
            j = g%pt(p)%j
            RealTmp(i,j) = g%bMag(p)
            Cnt(i,j) = Cnt(i,j)+1
        enddo
#ifdef par
        call sGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, RealTmp, g%nR, -1, -1 )
        call sGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, Cnt, g%nR, -1, -1 )
#endif
        if(iAm==0)nc_stat = nf90_put_var ( nc_id, bmod_id, RealTmp/Cnt )
        RealTmp = 0
        Cnt = 0

        allocate(RealTmp3(g%nR,g%nZ,nSpec))

        !Density
        RealTmp3 = 0
        do s=1,nSpec
            RealTmp = 0 
            Cnt = 0 
            do p=1,size(g%pt)
                i = g%pt(p)%i
                j = g%pt(p)%j
                RealTmp(i,j) = g%DensitySpec(p,s)
                Cnt(i,j) = Cnt(i,j)+1
            enddo
#ifdef par
            call sGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, RealTmp, g%nR, -1, -1 )
            call sGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, Cnt, g%nR, -1, -1 )
#endif
            RealTmp3(:,:,s) = RealTmp/Cnt
            RealTmp = 0
            Cnt = 0
        enddo
        if(iAm==0)nc_stat = nf90_put_var ( nc_id, dens_id, RealTmp3 )
        RealTmp3 = 0

        !nuOmg
        RealTmp3 = 0
        do s=1,nSpec
            RealTmp = 0 
            Cnt = 0 
            do p=1,size(g%pt)
                i = g%pt(p)%i
                j = g%pt(p)%j
                RealTmp(i,j) = g%nuOmg(p,s)
                Cnt(i,j) = Cnt(i,j)+1
            enddo
#ifdef par
            call sGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, RealTmp, g%nR, -1, -1 )
            call sGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, Cnt, g%nR, -1, -1 )
#endif
            RealTmp3(:,:,s) = RealTmp/Cnt
            RealTmp = 0
            Cnt = 0
        enddo
        if(iAm==0)nc_stat = nf90_put_var ( nc_id, nuOmg_id, RealTmp3 )
        RealTmp3 = 0

        !Temp_eV
        RealTmp3 = 0
        do s=1,nSpec
            RealTmp = 0 
            Cnt = 0 
            do p=1,size(g%pt)
                i = g%pt(p)%i
                j = g%pt(p)%j
                RealTmp(i,j) = g%ktSpec(p,s)/q
                Cnt(i,j) = Cnt(i,j)+1
            enddo
#ifdef par
            call sGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, RealTmp, g%nR, -1, -1 )
            call sGSUM2D ( iContext, 'All', ' ', g%nR, g%nZ, Cnt, g%nR, -1, -1 )
#endif
            RealTmp3(:,:,s) = RealTmp/Cnt
            RealTmp = 0
            Cnt = 0
        enddo
        if(iAm==0)nc_stat = nf90_put_var ( nc_id, temp_id, RealTmp3 )
        RealTmp3 = 0

        !JAnt
        if(iAm==0)nc_stat = nf90_put_var ( nc_id, jr_re_id, real(g%jR) )
        if(iAm==0)nc_stat = nf90_put_var ( nc_id, jr_im_id, aimag(g%jR) )
        if(iAm==0)nc_stat = nf90_put_var ( nc_id, jt_re_id, real(g%jT) )
        if(iAm==0)nc_stat = nf90_put_var ( nc_id, jt_im_id, aimag(g%jT) )
        if(iAm==0)nc_stat = nf90_put_var ( nc_id, jz_re_id, real(g%jZ) )
        if(iAm==0)nc_stat = nf90_put_var ( nc_id, jz_im_id, aimag(g%jZ) )

        if(iAm==0)call check ( nf90_close ( nc_id ) )

        deallocate(RealTmp,RealTmp3)

    end subroutine write_runData


    subroutine write_solution ( g, rid, rhs )

        use grid
        use aorsaNamelist, &
            only: nSpec, nPhi, freqcy
        use parallel

        implicit none

        type(gridBlock), intent(in) :: g
        integer, intent(in) :: rhs

        character(len=100) :: fName 
        character(len=20) :: rid
        character(len=6) :: rhs_string
        character(len=4) :: nPhi_string

        integer :: nc_id, nX_id, nY_id, nModesX_id, nModesY_id, nSpec_id, NRHS_id
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
        integer :: &
            jP_r_re_id, jP_r_im_id, &
            jP_t_re_id, jP_t_im_id, &
            jP_z_re_id, jP_z_im_id
        integer :: r_id, z_id,nPhi_id,freqcy_id,scalar_id
        integer :: stat

        integer :: jouleHeating_id

        write(nPhi_string,'(sp,i4.3)'), int(nPhi)
        write(rhs_string,'(i6.6)'), rhs

        fName = trim(rid)//'output/solution_'//g%fNumber//'_'//nPhi_string//'_'//rhs_string//'.nc'

        call check ( nf90_create ( fName, nf90_clobber, nc_id ) )
        stat=nf90_def_dim(nc_id,"scalar",1,scalar_id)
        call check ( nf90_def_dim ( nc_id, "nX", g%nR, nX_id ) )
        call check ( nf90_def_dim ( nc_id, "nY", g%nZ, nY_id ) )
        call check ( nf90_def_dim ( nc_id, "nModesX", g%nModesR, nModesX_id ) )
        call check ( nf90_def_dim ( nc_id, "nModesY", g%nModesZ, nModesY_id ) )
        call check ( nf90_def_dim ( nc_id, "nSpec", nSpec, nSpec_id ) )
        call check ( nf90_def_dim ( nc_id, "NRHS", NRHS, NRHS_id ) )

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

        call check ( nf90_def_var ( nc_id, "jP_r_re", NF90_REAL, &
            (/nX_id,nY_id,nSpec_id/), jP_r_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "jP_r_im", NF90_REAL, &
            (/nX_id,nY_id,nSpec_id/), jP_r_im_id ) ) 
        call check ( nf90_def_var ( nc_id, "jP_t_re", NF90_REAL, &
            (/nX_id,nY_id,nSpec_id/), jP_t_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "jP_t_im", NF90_REAL, &
            (/nX_id,nY_id,nSpec_id/), jP_t_im_id ) ) 
        call check ( nf90_def_var ( nc_id, "jP_z_re", NF90_REAL, &
            (/nX_id,nY_id,nSpec_id/), jP_z_re_id ) ) 
        call check ( nf90_def_var ( nc_id, "jP_z_im", NF90_REAL, &
            (/nX_id,nY_id,nSpec_id/), jP_z_im_id ) ) 
 
        call check ( nf90_def_var ( nc_id, "jouleHeating", NF90_REAL, &
            (/nX_id,nY_id,nSpec_id/), jouleHeating_id ) ) 
 
        stat = nf90_def_var(nc_id,"r",NF90_REAL,(/nx_id/),r_id)
        stat = nf90_def_var(nc_id,"z",NF90_REAL,(/ny_id/),z_id)
        stat = nf90_def_var(nc_id,"freqcy",NF90_REAL,(/scalar_id/),freqcy_id)
        stat = nf90_def_var(nc_id,"nPhi",NF90_INT,(/scalar_id/),nPhi_id)

        call check ( nf90_enddef ( nc_id ) )

        call check ( nf90_put_var ( nc_id, e1_re_id, real(real(g%ealpha)) ) )
        call check ( nf90_put_var ( nc_id, e1_im_id, real(aimag(g%ealpha)) ) )
        call check ( nf90_put_var ( nc_id, e2_re_id, real(real(g%ebeta)) ) )
        call check ( nf90_put_var ( nc_id, e2_im_id, real(aimag(g%ebeta)) ) )
        call check ( nf90_put_var ( nc_id, e3_re_id, real(real(g%eB)) ) )
        call check ( nf90_put_var ( nc_id, e3_im_id, real(aimag(g%eB)) ) )

        call check ( nf90_put_var ( nc_id, e1k_re_id, real(real(g%ealphak(:,:,rhs))) ) )
        call check ( nf90_put_var ( nc_id, e1k_im_id, real(aimag(g%ealphak(:,:,rhs))) ) )
        call check ( nf90_put_var ( nc_id, e2k_re_id, real(real(g%ebetak(:,:,rhs))) ) )
        call check ( nf90_put_var ( nc_id, e2k_im_id, real(aimag(g%ebetak(:,:,rhs))) ) )
        call check ( nf90_put_var ( nc_id, e3k_re_id, real(real(g%eBk(:,:,rhs))) ) )
        call check ( nf90_put_var ( nc_id, e3k_im_id, real(aimag(g%eBk(:,:,rhs))) ) )

        call check ( nf90_put_var ( nc_id, er_re_id, real(real(g%eR)) ) )
        call check ( nf90_put_var ( nc_id, er_im_id, real(aimag(g%eR)) ) )
        call check ( nf90_put_var ( nc_id, et_re_id, real(real(g%eTh)) ) )
        call check ( nf90_put_var ( nc_id, et_im_id, real(aimag(g%eTh)) ) )
        call check ( nf90_put_var ( nc_id, ez_re_id, real(real(g%eZ)) ) )
        call check ( nf90_put_var ( nc_id, ez_im_id, real(aimag(g%eZ)) ) )

        call check ( nf90_put_var ( nc_id, j1_re_id, real(real(g%jalpha)) ) )
        call check ( nf90_put_var ( nc_id, j1_im_id, real(aimag(g%jalpha)) ) )
        call check ( nf90_put_var ( nc_id, j2_re_id, real(real(g%jbeta)) ) )
        call check ( nf90_put_var ( nc_id, j2_im_id, real(aimag(g%jbeta)) ) )
        call check ( nf90_put_var ( nc_id, j3_re_id, real(real(g%jB)) ) )
        call check ( nf90_put_var ( nc_id, j3_im_id, real(aimag(g%jB)) ) )

        call check ( nf90_put_var ( nc_id, jP_r_re_id, real(real(g%jP_r)) ) )
        call check ( nf90_put_var ( nc_id, jP_r_im_id, real(aimag(g%jP_r)) ) )
        call check ( nf90_put_var ( nc_id, jP_t_re_id, real(real(g%jP_t)) ) )
        call check ( nf90_put_var ( nc_id, jP_t_im_id, real(aimag(g%jP_t)) ) )
        call check ( nf90_put_var ( nc_id, jP_z_re_id, real(real(g%jP_z)) ) )
        call check ( nf90_put_var ( nc_id, jP_z_im_id, real(aimag(g%jP_z)) ) )

        call check ( nf90_put_var ( nc_id, jouleHeating_id, g%jouleHeating ) )

        stat = nf90_put_var(nc_id,r_id,g%r)
        stat = nf90_put_var(nc_id,z_id,g%z)
        stat = nf90_put_var(nc_id,freqcy_id,freqcy)
        stat = nf90_put_var(nc_id,nPhi_id,nPhi)

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


    subroutine write_sigma_pt ( sigma, nc_id, sigma_re_id, sigma_im_id )!, nN, nM, nSpec )

        implicit none

        complex, intent(in) :: sigma(:,:,:,:,:,:,:)
        integer, intent(in) :: nc_id, sigma_re_id, sigma_im_id!,i, j 
        !integer, intent(in) :: nN, nM, nSpec

        !call check ( nf90_put_var ( nc_id, sigma_re_id, realpart ( sigma ), &
        !    start = (/ i, j, 1, 1, 1, 1, 1 /), &
        !    count = (/ 1, 1, nN, nM, 3, 3, nSpec /) ) )
        !call check ( nf90_put_var ( nc_id, sigma_im_id, imagpart ( sigma ), &
        !    start = (/ i, j, 1, 1, 1, 1, 1 /), &
        !    count = (/ 1, 1, nN, nM, 3, 3, nSpec /) ) )

        call check ( nf90_put_var ( nc_id, sigma_re_id, real(real(sigma)) ) )
        call check ( nf90_put_var ( nc_id, sigma_im_id, real(aimag(sigma)) ) )

    end subroutine write_sigma_pt


    subroutine close_sigma_file ( nc_id )

        implicit none

        integer, intent(in) :: nc_id

        call check ( nf90_close ( nc_id ) )

    end subroutine close_sigma_file

    subroutine WritePerformanceData ( P, rid )

        use Performance

        implicit none

        type(RunPerfData), intent(in) :: P

        integer :: stat
        character(len=100) :: fName 
        character(len=20) :: rid

        integer :: nc_id, scalar_id
        integer :: &
            nProcs_id, &
            nSpatialPoints_id, &
            nRowLocal_id, &
            nColLocal_id, &
            nRowGlobal_id, &
            nColGlobal_id, &
            MatSizeLocal_GB_id, &
            MatSizeGlobal_GB_id, &
            TimeWorkList_id, &
            TimeFill_id, &
            TimeSolve_id, &
            TimeCurrent_id, &
            TimeTotal_id, &
            GflopsFillLocal_id, &
            GflopsFillGlobal_id, &
            GflopsSolveLocal_id, &
            GflopsSolveGlobal_id, &
            GflopsCurrentLocal_id, &
            GflopsCurrentGlobal_id

        fName = trim(rid)//'PerfData.nc'

        stat = nf90_create ( fName, nf90_clobber, nc_id ) 
        stat = nf90_def_dim ( nc_id, "scalar", 1, scalar_id ) 

        stat = nf90_def_var (nc_id,"nProcs",NF90_INT,(/scalar_id/),nProcs_id)
        stat = nf90_def_var (nc_id,"nSpatialPoints",NF90_INT,(/scalar_id/),nSpatialPoints_id)
        stat = nf90_def_var (nc_id,"nRowLocal",NF90_INT,(/scalar_id/),nRowLocal_id)
        stat = nf90_def_var (nc_id,"nColLocal",NF90_INT,(/scalar_id/),nColLocal_id)
        stat = nf90_def_var (nc_id,"nRowGlobal",NF90_INT,(/scalar_id/),nRowGlobal_id)
        stat = nf90_def_var (nc_id,"nColGlobal",NF90_INT,(/scalar_id/),nColGlobal_id)
        stat = nf90_def_var (nc_id,"MatSizeLocal_GB",NF90_REAL,(/scalar_id/),MatSizeLocal_GB_id)
        stat = nf90_def_var (nc_id,"MatSizeGlobal_GB",NF90_REAL,(/scalar_id/),MatSizeGlobal_GB_id)
        stat = nf90_def_var (nc_id,"TimeWorkList",NF90_REAL,(/scalar_id/),TimeWorkList_id)
        stat = nf90_def_var (nc_id,"TimeFill",NF90_REAL,(/scalar_id/),TimeFill_id)
        stat = nf90_def_var (nc_id,"TimeSolve",NF90_REAL,(/scalar_id/),TimeSolve_id)
        stat = nf90_def_var (nc_id,"TimeCurrent",NF90_REAL,(/scalar_id/),TimeCurrent_id)
        stat = nf90_def_var (nc_id,"TimeTotal",NF90_REAL,(/scalar_id/),TimeTotal_id)
        stat = nf90_def_var (nc_id,"GflopsFillLocal",NF90_REAL,(/scalar_id/),GflopsFillLocal_id)
        stat = nf90_def_var (nc_id,"GflopsFillGlobal",NF90_REAL,(/scalar_id/),GflopsFillGlobal_id)
        stat = nf90_def_var (nc_id,"GflopsSolveLocal",NF90_REAL,(/scalar_id/),GflopsSolveLocal_id)
        stat = nf90_def_var (nc_id,"GflopsSolveGlobal",NF90_REAL,(/scalar_id/),GflopsSolveGlobal_id)
        stat = nf90_def_var (nc_id,"GflopsCurrentLocal",NF90_REAL,(/scalar_id/),GflopsCurrentLocal_id)
        stat = nf90_def_var (nc_id,"GflopsCurrentGlobal",NF90_REAL,(/scalar_id/),GflopsCurrentGlobal_id)

        stat = nf90_enddef ( nc_id )

        stat = nf90_put_var ( nc_id, nProcs_id,      P%nProcs )
        stat = nf90_put_var ( nc_id, nSpatialPoints_id,      P%nSpatialPoints )
        stat = nf90_put_var ( nc_id, nRowLocal_id,           P%nRowLocal)          
        stat = nf90_put_var ( nc_id, nColLocal_id,           P%nColLocal)          
        stat = nf90_put_var ( nc_id, nRowGlobal_id,          P%nRowGlobal)         
        stat = nf90_put_var ( nc_id, nColGlobal_id,          P%nColGlobal)         
        stat = nf90_put_var ( nc_id, MatSizeLocal_GB_id,     P%MatSizeLocal_GB)    
        stat = nf90_put_var ( nc_id, MatSizeGlobal_GB_id,    P%MatSizeGlobal_GB)   
        stat = nf90_put_var ( nc_id, TimeWorkList_id,        P%TimeWorkList)       
        stat = nf90_put_var ( nc_id, TimeFill_id,            P%TimeFill)           
        stat = nf90_put_var ( nc_id, TimeSolve_id,           P%TimeSolve)          
        stat = nf90_put_var ( nc_id, TimeCurrent_id,         P%TimeCurrent)        
        stat = nf90_put_var ( nc_id, TimeTotal_id,           P%TimeTotal)          
        stat = nf90_put_var ( nc_id, GflopsFillLocal_id,     P%GflopsFillLocal)    
        stat = nf90_put_var ( nc_id, GflopsFillGlobal_id,    P%GflopsFillGlobal)   
        stat = nf90_put_var ( nc_id, GflopsSolveLocal_id,    P%GflopsSolveLocal)   
        stat = nf90_put_var ( nc_id, GflopsSolveGlobal_id,   P%GflopsSolveGlobal) 
        stat = nf90_put_var ( nc_id, GflopsCurrentLocal_id,  P%GflopsCurrentLocal)  
        stat = nf90_put_var ( nc_id, GflopsCurrentGlobal_id, P%GflopsCurrentGlobal)  

        stat = nf90_close ( nc_id )

    end subroutine WritePerformanceData


end module write_data
