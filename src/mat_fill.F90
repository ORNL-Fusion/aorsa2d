module mat_fill

use constants

implicit none

complex :: bFn, dRBFn, dZBFn, d2RbFn, d2ZbFn, &
    bFnHere, d3RbFn, d4RbFn
complex(kind=dbl) :: &
    sigAlpAlp, sigAlpBet, sigAlpPrl, &
    sigBetAlp, sigBetBet, sigBetPrl, &
    sigPrlAlp, sigPrlBet, sigPrlPrl
complex(kind=dbl) :: &
    sigAlpAlpTmp, sigAlpBetTmp, sigAlpPrlTmp, &
    sigBetAlpTmp, sigBetBetTmp, sigBetPrlTmp, &
    sigPrlAlpTmp, sigPrlBetTmp, sigPrlPrlTmp

complex :: &
    kAlpAlp, kAlpBet, kAlpPrl, &
    kBetAlp, kBetBet, kBetPrl, &
    kPrlAlp, kPrlBet, kPrlPrl
real :: rnx, rny, rnPhi 
complex :: &
    fdk, fek, ffk, &
    fgk, fak, fpk, &
    frk, fqk, fsk
complex, allocatable, dimension(:) :: &
    sss, ttt, qqq

real :: GlobalSizeMB,LocalSizeMB,GlobalSizeGB,LocalSizeGB

#ifndef dblprec
    complex, allocatable :: aMat(:,:)
#else
    complex(kind=dbl), allocatable :: aMat(:,:)
#endif


contains

    subroutine alloc_total_aMat ( nPts_tot )

        use parallel
        use scalapack_mod

        implicit none
        integer, intent(in) :: nPts_tot

#ifdef par

        if (iAm == 0) then
            write(*,*) '    nSpatialPts_tot: ', nPts_tot
            write(*,*) '    nRowLocal: ', LM_A
            write(*,*) '    nColLocal: ', LN_A
            write(*,*) '    nRowGlobal: ', DESCA(M_)
            write(*,*) '    nColGlobal: ', DESCA(N_) 

            LocalSizeMB = LM_A*LN_A*2.0*8.0 / 1024.0**2
            GlobalSizeMB = nPts_tot*3.0*nPts_tot*3.0*2.0*8.0 / 1024.0**2.0
            LocalSizeGB = LocalSizeMB/1024.0
            GlobalSizeGB = GlobalSizeMB/1024.0

            if(GlobalSizeMB<=1)then
                write(*,100), GlobalSizeMB, LocalSizeMB 
                100 format (' Filling aMat [global size: ',f8.1,' MB, local size: ',f8.1' MB]')
            else
                write(*,101), GlobalSizeGB, LocalSizeMB 
                101 format (' Filling aMat [global size: ',f8.1,' GB, local size: ',f8.1' MB]')
            endif
        endif

        allocate ( aMat(LM_A,LN_A) )
#else 
        write(*,*) '    nPts_tot: ', nPts_tot
        write(*,100), &
            nPts_tot*3*nPts_tot*3*2*8.0 / 1024.0**2
        100 format (' Filling aMat [global size: ',f8.1,' MB]')

        allocate ( aMat(nPts_tot*3,nPts_tot*3) )
        write(*,*) '    ',nPts_tot*3,nPts_tot*3  

#endif

        aMat = 0

    end subroutine alloc_total_aMat


    subroutine amat_boundaries ( gAll, nPts_tot )

        use grid
        use antenna
        use aorsaNamelist, &
            only: overlap
        use spline_dlg
        use parallel, only: NRHS

        implicit none

        type(gridBlock), intent(in) :: gAll(:)
        integer, intent(in) :: nPts_tot

        type(gridBlock) :: me, nbr
        integer :: i, j, n, m, iRow, iCol, rhs
        integer(kind=long) :: ii 
        complex(kind=dbl) :: aMatBlock(3,3)
        real :: r, kt
        complex :: bFn_iL, bFn_iR, bFn_iRR
        complex :: bFn_i0, bFn_i1, bFn_i2, bFn_i3, bFn_i4, bFn_i5
        complex :: bFn_R, bFn_Z, bFn_0
        real :: rNorm, zNorm


        real :: h1, h2, alpha, coeffL, coeffR, splines(3,4)
        integer :: iiL, iiR, iiRR, iiArr(1)

        do ii=1,nPts_tot

            if(bndryBlockID(3,ii)/=0)then

                write(*,*) 'BndryType: ', bndryType(ii)
                write(*,*) 'BndryBlockID: ', bndryBlockID(:,ii)

                me = gAll(bndryBlockID(3,ii))

                do n=me%nMin,me%nMax
                    do m=me%mMin,me%mMax

                        aMatBlock = 0

                        if(bndryType(ii)==-999) then

                            ! create E=E @ overlap point

                            i = me%nR-overlap
                            j = bndryBlockID(2,ii)

                            bFn = me%xx(n, i) * me%yy(m, j)

                            aMatBlock(1,1) = -bFn
                            aMatBlock(2,2) = -bFn
                            aMatBlock(3,3) = -bFn

                            !aMatBlock = -99

                        elseif(bndryType(ii)<=-1 .and. bndryType(ii)>=-20) then

                            ! Dirichlet 
                            ! ---------

                            nbr = gAll(bndryBlockID(4,ii))
                            i = bndryBlockID(1,ii)
                            j = bndryBlockID(2,ii)

                            ! linear interpolation
                            ! --------------------

                            ! determine left index for linear interpolation
                            iiArr = maxLoc ( me%R-nbr%R(i), (me%R-nbr%R(i))<0)
                            iiL = iiArr(1)
                            iiR = iiL + 1

                            coeffL = 1 - nbr%R(i) / ( -me%R(iiL)+me%R(iiR) ) &
                                + me%R(iiL) / ( -me%R(iiL)+me%R(iiR) )
                            coeffR = nbr%R(i) / ( -me%R(iiL)+me%R(iiR) ) &
                                - me%R(iiL) / ( -me%R(iiL)+me%R(iiR) )

                            bFn_iL =me%xx(n, iiL) * me%yy(m, j)
                            bFn_iR =me%xx(n, iiR) * me%yy(m, j)

                            aMatBlock(1,1) = -(coeffL*bFn_iL + coeffR*bFn_iR)
                            aMatBlock(2,2) = -(coeffL*bFn_iL + coeffR*bFn_iR)
                            aMatBlock(3,3) = -(coeffL*bFn_iL + coeffR*bFn_iR)

                        endif

                        ! but couple with the neighbour block

                        i = bndryBlockID(1,ii)
                        j = bndryBlockID(2,ii)

                        nbr = gAll(bndryBlockID(4,ii))

                        iRow = (i-1) * 3 * nbr%nZ + (j-1) * 3 + 1
                        iRow = iRow + ( nbr%startRow-1 )

                        iCol = (n-me%nMin) * 3 * me%nModesZ + (m-me%mMin) * 3 + 1
                        iCol = iCol + ( me%startCol-1 )

                        aMat(iRow:iRow+2,iCol:iCol+2) = aMatBlock
                        do rhs=1,NRHS
                            brhs(iRow:iRow+2,rhs) = 0
                        enddo

                    enddo
                enddo

            endif

        enddo

    end subroutine amat_boundaries


    subroutine createWorkList ( g )

        use grid
        use parallel, only : iAm, ICTXT, NRHS, LM_A, LN_A, &
                MYROW, MYCOL, DESCA, DESCB
        use aorsaNamelist, &
            only: npRow, npCol
        use scalapack_mod
 
        implicit none

        type(gridBlock), intent(inout) :: g

        integer :: workListPosition, i, j, n, m
        type(workListEntry) :: thisWorkList
        type(workListEntry), allocatable :: workListTooLong(:),WorkListJustRight(:)

        type(SpatialRow), allocatable ::  SpatialRowsTmp(:)
        integer :: Cnt, row, w

#ifdef par
        integer :: iCol,iRow,ii,jj
        !   scalapack indicies
        !   see http://www.netlib.org/scalapack/slug/node76.html
        integer :: l_sp, m_sp, pr_sp, pc_sp, x_sp, y_sp
        integer :: pr_sp_thisPt(3), pc_sp_thisPt(3)
#endif


        ! Create a work list for this processor
        ! -------------------------------------

#ifndef par
        allocate(g%wl(g%nR*g%nZ*g%nModesR*g%nModesZ))
        workListPosition = 0

        i_workList: &
        do i=1,g%nR
            j_workList: &
            do j=1,g%nZ

                n_workList: &
                do n=g%nMin,g%nMax
                    m_workList: &
                    do m=g%mMin,g%mMax

                        workListPosition = workListPosition + 1
                        g%wl(workListPosition) = workListEntry(i,j,m,n,0)

                    enddo m_workList
                enddo n_workList

            enddo j_workList
        enddo i_workList
#else
        allocate(workListTooLong(LM_A*LN_A))

        ! Touch this piece of code and die!
        WorkListPosition = 0
        do ii=1,LM_A,3
            do jj=1,LN_A,3
       
                iRow = IndxL2G ( ii, DESCA(MB_), MyRow, 0, NpRow )
                iCol = IndxL2G ( jj, DESCA(NB_), MyCol, 0, NpCol )

                i = (iRow-1)/(3*g%nZ)+1
                j = (mod(iRow-1,3*g%nZ)+1-1)/3+1
                n = (iCol-1)/(3*g%nModesZ)+1+g%nMin-1
                m = (mod(iCol-1,3*g%nModesZ)+1-1)/3+1+g%mMin-1
           
                workListPosition = workListPosition + 1
                thisWorkList = workListEntry(i,j,m,n,0)
                workListTooLong(workListPosition) = thisWorkList

            enddo
        enddo

        allocate(g%wl(workListPosition))
        g%wl = workListTooLong(1:workListPosition)

        if(iAM==0)write(*,*) '    MyWork:  ', size(g%wl)
        if(iAM==0)write(*,*) '    AllWork: ', g%nR*g%nZ*g%nModesR*g%nModesZ

#endif

        allocate(SpatialRowsTmp(size(g%wl)))
        SpatialRowsTmp = SpatialRow(0,0,0)
        Cnt = 1

        if(iAm==0)write(*,*) '    Calculating wl -> space mapping array ...'
        do w=1,size(g%wl)
            i = g%wl(w)%i
            j = g%wl(w)%j
            if(i<1.or.i>g%nR.or.j<1.or.j>g%nZ)then
                    write(*,*) 'ERROR: ', i, j
                    stop
            endif
            Row = (i-1)*g%nZ + j
            if(any(SpatialRowsTmp(1:Cnt)%row==Row))then
                g%wl(w)%iPt = minloc(abs(SpatialRowsTmp(1:Cnt)%row-Row),1)
            else
                g%wl(w)%iPt = Cnt
                SpatialRowsTmp(Cnt)%row = Row
                SpatialRowsTmp(Cnt)%i = i 
                SpatialRowsTmp(Cnt)%j = j
                !if(iAm==0)write(*,*) '        Cnt:', cnt, SpatialRowsTmp(Cnt)
                Cnt = Cnt + 1
            endif
        enddo

        Cnt = Cnt - 1

        allocate(g%pt(Cnt))
        g%pt = SpatialRowsTmp(1:Cnt) 
        deallocate(SpatialRowsTmp)

        if(iAm==0)write(*,*) '    DONE'

    end subroutine createWorkList


    subroutine amat_fill ( g )

        use aorsaNamelist, only: &
            delta0, nSpec, &
            iSigma, npRow, npCol, &
            nPhi, square, lsWeightFac, &
            useEqdsk, overlap
        use grid
        use sigma
        use rotation
        use profiles
        use bField
        use parallel, only : iAm, ICTXT, NRHS, DESCA, DESCB
        use chebyshev_mod
        use write_data
        use getMatElements
        use fitPack
        use scalapack_mod

        implicit none

        type(gridBlock), intent(in) :: g

        type(dBfnArg) :: d
        integer :: iRow, iCol, i, j, n, m, p, s, ii, jj, iOL, jOL
        integer :: localRow, localCol
        real :: kr, kt, kz, r, z, kVec_stix(3)

        complex :: mat3by3Block(3,3)

        complex :: bFn_iL, bFn_iR, bFn_i0, bFn_i1, bFn_i2, bFn_i3, bFn_i4, bFn_i5 
        real :: h1, h2, alpha

        complex(kind=dbl) :: sigma_tmp(3,3), sigma_tmp_neg(3,3), sigmaHere(3,3)
        complex, allocatable :: sigma_write(:,:,:,:,:,:,:)

        integer :: nr, nz, iStat
        integer :: sigma_nc_id, sigma_re_id, sigma_im_id


#ifdef par
        !   scalapack indicies
        !   see http://www.netlib.org/scalapack/slug/node76.html

        integer :: l_sp, m_sp, pr_sp, pc_sp, x_sp, y_sp
        integer :: pr_sp_thisPt(3), pc_sp_thisPt(3)
        integer :: Dummy
#endif

        integer :: w
        integer :: aa
        real :: aa_rPts(4), aa_zPts(4)


        if(square) lsWeightFac = 1

#if __sigma__ != 2

        ! Initialize grid sigma file
        ! --------------------------

        if(iAm==0) &
        call init_sigma_file ( g, 'sigma'//g%fNumber//'.nc', &
            sigma_nc_id, sigma_re_id, sigma_im_id, nSpec )

        allocate ( sigma_write(g%nR,g%nZ,g%nMin:g%nMax,g%mMin:g%mMax,3,3,nSpec), stat = iStat )
        if(iStat/=0)then
                write(*,*) 'ERROR src/mat_fill.f90 - allocation failed :('
                stop
        endif

        sigma_write = 0

        ! Calculate sigma seperately outside main loop
        ! --------------------------------------------

        do s=1,nSpec

        if(iAm==0) &
        write(*,*) 'Calculating sigma for species ', s, ' of', nSpec

            do m=g%mMin,g%mMax
                do n=g%nMin,g%nMax
!#ifndef par
!                !   progress indicator
!                !   ------------------
!                do p=1,7 
!                    write(*,'(a)',advance='no') char(8)
!                enddo
!                write(*,'(1x,f5.1,a)',advance='no') &
!                    real((m-g%mMin)*g%nR+(n-g%nMin))/(g%nR*g%nZ)*100, '%'
!#endif
                    do j=1,g%nZ
                        do i=1,g%nR

                            if(chebyshevX) then
                                if(n>1) then
                                    kr = n / sqrt ( sin ( pi * (g%rNorm(i)+1)/2  ) ) * g%normFacR 
                                else
                                    kr = n * g%normFacR
                                endif
                            else
                                kr = n * g%normFacR
                            endif

                            if(chebyshevY) then
                                if(m>1) then
                                    kz = m / sqrt ( sin ( pi * (g%zNorm(j)+1)/2 ) ) * g%normFacZ 
                                else
                                    kz = m * g%normFacZ
                                endif
                            else
                                kz = m * g%normFacZ
                            endif

  
                            hotPlasma:& 
                            if (iSigma==1 .and. (.not. g%isMetal(i,j)) ) then        

                                kVec_stix = matMul( g%U_RTZ_to_ABb(i,j,:,:), (/ kr, g%kPhi(i), kz /) ) 

                                sigma_tmp = sigmaHot_maxwellian&
                                    ( mSpec(s), &
                                    g%ktSpec(i,j,s), g%omgc(i,j,s), g%omgp2(i,j,s), &
                                    kVec_stix, g%R(i), &
                                    omgrf, k0, &
                                    g%k_cutoff, s, &
                                    g%sinTh(i,j), g%bPol(i,j), g%bMag(i,j), g%gradPrlB(i,j), &
                                    g%nuOmg(i,j) )

                            endif hotPlasma


                            coldPlasma: &
                            if (iSigma==0 .and. (.not. g%isMetal(i,j)) ) then 

                                sigma_tmp = sigmaCold_stix &
                                    ( g%omgc(i,j,s), g%omgp2(i,j,s), omgrf, &
                                    g%nuOmg(i,j) )

                            endif coldPlasma

                            sigma_write(i,j,n,m,:,:,s) = sigma_tmp

                        enddo
                    enddo

                enddo
            enddo

        enddo
#ifndef par
        write(*,*)
#endif

        ! Write sigma
        ! -----------
        if(iAm==0) &
        write(*,*) 'Writing sigma ...'
        if(iAm==0) &
        call write_sigma_pt ( sigma_write, &
            sigma_nc_id, sigma_re_id, sigma_im_id )
        if(iAm==0) &
        write(*,*) 'DONE'

        ! Close sigma file
        ! ----------------

        if(iAm==0) &
        call close_sigma_file ( sigma_nc_id )

#endif 
! __sigma__ != 2



        ! Begin loop
        ! ----------

        workListLoop: &
        do w=1,size(g%wl)

            iRow = (g%wl(w)%i-1) * 3 * g%nZ + (g%wl(w)%j-1) * 3 + 1
            iRow = iRow + ( g%startRow-1 )

            iCol = (g%wl(w)%n-g%nMin) * 3 * g%nModesZ + (g%wl(w)%m-g%mMin) * 3 + 1
            iCol = iCol + ( g%startCol-1 )

            bFn = g%xx(g%wl(w)%n, g%wl(w)%i) * g%yy(g%wl(w)%m, g%wl(w)%j)

            interior: &
            if(g%label(g%wl(w)%iPt)==0)then
   
                mat3by3Block = get3by3Block ( g, w)

            endif interior

            ii_loop: &
            do ii=0,2
                jj_loop: &
                do jj=0,2
#ifdef par
                    pr_sp = IndxG2P ( iRow, DESCA(MB_), Dummy, 0, NpRow )
                    pc_sp = IndxG2P ( iCol, DESCA(NB_), Dummy, 0, NpCol )
#if __CheckParallelLocation__==1
                    myProc: &
                    if ( myRow==pr_sp .and. myCol==pc_sp ) then
#endif
                        LocalRow = IndxG2L ( iRow+ii, DESCA(MB_), Dummy, Dummy, NpRow )
                        LocalCol = IndxG2L ( iCol+jj, DESCA(NB_), Dummy, Dummy, NpCol )
#else
                        LocalRow    = iRow+ii
                        LocalCol    = iCol+jj
#endif
                            
                        ! Interior points
                        ! ---------------

                        if(g%label(g%wl(w)%iPt)==0) then 

                            if (ii==0 .and. jj==0) aMat(localRow,localCol) = mat3by3Block(1,1) * bFn  
                            if (ii==0 .and. jj==1) aMat(localRow,localCol) = mat3by3Block(1,2) * bFn   
                            if (ii==0 .and. jj==2) aMat(localRow,localCol) = mat3by3Block(1,3) * bFn   
                                                                                  
                            if (ii==1 .and. jj==0) aMat(localRow,localCol) = mat3by3Block(2,1) * bFn   
                            if (ii==1 .and. jj==1) aMat(localRow,localCol) = mat3by3Block(2,2) * bFn   
                            if (ii==1 .and. jj==2) aMat(localRow,localCol) = mat3by3Block(2,3) * bFn   
                                                                             
                            if (ii==2 .and. jj==0) aMat(localRow,localCol) = mat3by3Block(3,1) * bFn   
                            if (ii==2 .and. jj==1) aMat(localRow,localCol) = mat3by3Block(3,2) * bFn   
                            if (ii==2 .and. jj==2) aMat(localRow,localCol) = mat3by3Block(3,3) * bFn   

                        endif


                        ! Outer boundary points
                        ! ---------------------

                        if (g%label(g%wl(w)%iPt)==888) then
                            
                            if (ii==0 .and. jj==0) aMat(localRow,localCol) = bFn * lsWeightFac 
                            if (ii==0 .and. jj==1) aMat(localRow,localCol) = 0
                            if (ii==0 .and. jj==2) aMat(localRow,localCol) = 0
                        
                            if (ii==1 .and. jj==0) aMat(localRow,localCol) = 0
                            if (ii==1 .and. jj==1) aMat(localRow,localCol) = bFn * lsWeightFac  
                            if (ii==1 .and. jj==2) aMat(localRow,localCol) = 0
                        
                            if (ii==2 .and. jj==0) aMat(localRow,localCol) = 0
                            if (ii==2 .and. jj==1) aMat(localRow,localCol) = 0
                            if (ii==2 .and. jj==2) aMat(localRow,localCol) = bFn * lsWeightFac 

                        endif

                        if (g%label(g%wl(w)%iPt)==999) then

                            if (ii==0 .and. jj==0) aMat(localRow,localCol) = bFn
                            if (ii==0 .and. jj==1) aMat(localRow,localCol) = 0  
                            if (ii==0 .and. jj==2) aMat(localRow,localCol) = 0  
                        
                            if (ii==1 .and. jj==0) aMat(localRow,localCol) = 0  
                            if (ii==1 .and. jj==1) aMat(localRow,localCol) = bFn
                            if (ii==1 .and. jj==2) aMat(localRow,localCol) = 0   
                        
                            if (ii==2 .and. jj==0) aMat(localRow,localCol) = 0   
                            if (ii==2 .and. jj==1) aMat(localRow,localCol) = 0   
                            if (ii==2 .and. jj==2) aMat(localRow,localCol) = bFn

                        endif


                        ! Mesh-Mesh boundary for R bFn 
                        ! ----------------------------

                        if (g%label(g%wl(w)%iPt)>=1 .and. g%label(g%wl(w)%iPt)<=20) then

                            if (ii==0 .and. jj==0) aMat(localRow,localCol) = bFn!Here 
                            if (ii==0 .and. jj==1) aMat(localRow,localCol) = 0  
                            if (ii==0 .and. jj==2) aMat(localRow,localCol) = 0  
                        
                            if (ii==1 .and. jj==0) aMat(localRow,localCol) = 0  
                            if (ii==1 .and. jj==1) aMat(localRow,localCol) = bFn!Here  
                            if (ii==1 .and. jj==2) aMat(localRow,localCol) = 0   
                        
                            if (ii==2 .and. jj==0) aMat(localRow,localCol) = 0   
                            if (ii==2 .and. jj==1) aMat(localRow,localCol) = 0   
                            if (ii==2 .and. jj==2) aMat(localRow,localCol) = bFn!Here 

                        endif

#ifdef par
#if __CheckParallelLocation__==1
                   else
                       write(*,*) 'THE CODE SHOULD NEVER GET HERE!'
                   endif myProc
#endif
#endif
               enddo jj_loop
           enddo ii_loop

        enddo workListLoop

#if __sigma__ != 2
        deallocate ( sigma_write )
#endif

    end subroutine amat_fill

end module mat_fill
