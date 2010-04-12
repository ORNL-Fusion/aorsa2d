module solve

implicit none

complex, allocatable, dimension(:,:) :: &
    ealphak, ebetak, eBk
complex, allocatable, dimension(:,:) :: &
   ealpha, ealphax, ealphay, &
   ebeta, ebetax, ebetay, &
   eB, eBx, eBy


contains

    subroutine solve_lsq ()

        use aorsa2din_mod, &
        only: nPtsX, nPtsY, nModesX, nModesY
        use mat_fill, &
        only: aMat
        use antenna, &
        only: brhs, brhs_global

        implicit none

        integer :: info, nRow, nCol
        integer, allocatable, dimension(:) :: ipiv

        !   lapack lsq solve variables
        
        integer :: M_, N_, NRHS, LDA, LDB, MN_, LWORK, RANK
        real :: RCOND
        real, allocatable :: rWork(:)
        complex, allocatable :: work(:)
        integer, allocatable, dimension(:) :: jpvt
          
        nRow    = nPtsX * nPtsY * 3
        nCol    = nModesX * nModesY * 3

        !allocate ( ipiv ( nRow ), jpvt ( nCol ) )

        !call cgesv ( nRow, 1, aMat, nRow, ipiv, brhs, nRow, info )

        M_  = nRow
        N_  = nCol
        NRHS   = 1
        LDA = maxVal ( (/ 1, M_ /) )
        LDB = maxVal ( (/ 1, M_, N_ /) )
        MN_ = minVal ( (/ M_, N_ /) )
        LWORK   = MN_ + maxVal ( (/ 2 * MN_, N_ + 1, MN_ + NRHS /) )
        RCOND   = 1E-12

        allocate ( WORK ( maxVal ( (/ 1, LWORK /) ) ) )
        allocate ( RWORK ( 2 * N_ ) )
        allocate ( JPVT ( N_ ) ) 

        JPVT    = 0

        write(*,*) shape ( amat ), M_*N_, size(amat)

        call cgelsy ( M_, N_, NRHS, aMat, LDA, brhs, LDB, JPVT, RCOND, RANK, &
                            WORK, LWORK, RWORK, info )

        write(*,*) '    LAPACK status: ', info

        brhs_global = brhs

    end subroutine solve_lsq

#ifdef par

    subroutine solve_lsq_parallel ()

        use aorsa2din_mod, &
        only: nPtsX, nPtsY, nModesX, nModesY, &
            npRow, npCol, square
        use mat_fill, &
        only: aMat
        use antenna, &
        only: brhs, brhs_global
        use parallel

        implicit none

        character :: trans
        integer :: m, n, nrhs
        integer :: ia, ja
        complex :: b
        integer :: ib, jb
        complex, allocatable :: work(:)
        integer :: lWork, info
        integer :: ltau, lwf, lws, MpA0, NqA0, NRHSqB0, MpB0
        integer :: IROFFA, IAROW, IACOL, IROFFB, IBROW, IBCOL
        integer :: mb_a, nb_a, csrc_a, csrc_b, icoffa, mb_b
        integer :: icoffb, nb_b, Npb0, rsrc_a, rsrc_b

        integer, external :: ILCM, NUMROC, INDXG2P
        integer, allocatable :: ipiv(:)

        integer :: ii, gg

        trans   = 'N'

        mb_a    = rowBlockSize
        nb_a    = colBlockSize

        mb_b    = rowBlockSize
        nb_b    = 1

        rsrc_a  = rowStartProc
        csrc_a  = colStartProc
       
        rsrc_b  = rowStartProc 
        csrc_b  = colStartProc

        m   = nRow
        n   = nCol
        nrhs    = 1
        ia  = 1!myRow * rowBlockSize + 1
        ja  = 1!myCol * colBlockSize + 1 
        ib  = 1!myRow * rowBlockSize + 1
        jb  = 1

        iroffa  = mod( ia-1, mb_a )
        icoffa  = mod( ja-1, nb_a )
        iarow   = indxg2p( ia, mb_a, myRow, rsrc_a, npRow )
        iacol   = indxg2p( ja, nb_a, myCol, csrc_a, npCol )
        MpA0    = numroc( m+iroffa, mb_a, myrow, iarow, nprow )
        NqA0    = numroc( n+icoffa, nb_a, mycol, iacol, npcol )

        iroffb  = mod( ib-1, mb_b )
        icoffb  = mod( jb-1, nb_b )
        ibrow   = indxg2p( ib, mb_b, myRow, rsrc_b, npRow )
        ibcol   = indxg2p( jb, nb_b, myCol, csrc_b, npCol )
        Mpb0    = numroc( m+iroffb, mb_b, myRow, ibrow, npRow )
        Npb0    = numroc( n+iroffb, mb_b, myRow, ibrow, npRow )
        nrhsqb0 = numroc( nrhs+icoffb, nb_b, myCol, ibcol, npCol )

        ltau    = numroc( ja+min(m,n)-1, nb_a, myCol, csrc_a, npCol )
        lwf     = nb_a * ( Mpa0 + Nqa0 + nb_a )
        lws     = max( (nb_a*(nb_a-1))/2, (nrhsqb0 + Mpb0)*nb_a ) + &
                    nb_a * nb_a

        lWork   = ltau + max ( lwf, lws )


        ! Least squares or LU decomp for square array
        ! -------------------------------------------

        if(square) then 

            if(iAm==0) write(*,*) &
            "   using LU decomposition for square system"

            allocate ( ipiv ( numroc ( m, mb_a, myRow, rsrc_a, npRow ) + mb_a ) )

            call pcgesv ( n, nrhs, aMat, ia, ja, descriptor_aMat, ipiv, &
                    brhs, ib, jb, descriptor_brhs, info ) 

            if (iAm==0) then 
                write(*,*) '    pcgesv status: ', info
            endif

            deallocate ( ipiv )

        else 

            if(iAm==0) write(*,*) &
            "   using least squares for n x m system"

            allocate ( work ( lWork ) )

            call pcgels ( trans, m, n, nrhs, aMat, ia, ja, descriptor_aMat, &
                        brhs, ib, jb, descriptor_brhs, work, lWork, info ) 

            if (iAm==0) then 
                write(*,*) '    pcgels status: ', info
                write(*,*) '    actual/optimal lwork: ', lwork, real(work(1))
            endif

            deallocate ( work )

        endif 
      
        
        call gather_coeffs ()

    end subroutine solve_lsq_parallel


    subroutine gather_coeffs ()

        use aorsa2din_mod, &
        only: nPtsX, nPtsY, nModesX, nModesY, npRow, npCol
        use antenna, &
        only: brhs, brhs_global
        use parallel

        implicit none

        integer :: ii, gg

        !   Gather the solution vector from all processors by
        !   creating a global solution vector of the full, not
        !   local size, and have each proc fill in its piece.
        !   Then do a global sum on all procs such that all procs
        !   will have a complete copy of the solution.
        !   NOT SURE IF THIS WILL WORK FOR A NON 2x2 grid yet.

        if (myCol==0) then 
            do ii=1,nRowLocal

                gg   =  myRow*rowBlockSize+1 &
                        +mod(ii-1,rowBlockSize) &
                        +(ii-1)/rowBlockSize * npRow*rowBlockSize
                brhs_global(gg) = brhs(ii)

            enddo
        endif

        call cgSum2D ( iContext, 'All', ' ', nRow, 1, &
                brhs_global, nRow, -1, -1 )

    end subroutine gather_coeffs

#endif

    subroutine extract_coeffs ()

        use grid, &
        only: kxL, kxR, kyL, kyR
        use aorsa2din_mod, &
        only: nPtsX, nPtsY, nModesX, nModesY
        use antenna, &
        only: brhs_global
 
        implicit none

        integer :: m, n, iRow

        !   Extract the k coefficents from the solution
        !   for each field component
        !   -------------------------------------------

        allocate ( &
            ealphak(kxL:kxR,kyL:kyR), &
            ebetak(kxL:kxR,kyL:kyR), &
            eBk(kxL:kxR,kyL:kyR) )

        do n=kxL,kxR
            do m=kyL,kyR

                iRow = (m-kyL) * 3 + (n-kxL) * nModesY * 3 + 1
        
                ealphak(n,m)    = brhs_global(iRow)
                ebetak(n,m)     = brhs_global(iRow+1)
                eBk(n,m)        = brhs_global(iRow+2)

            enddo
        enddo

    end subroutine extract_coeffs

end module solve
