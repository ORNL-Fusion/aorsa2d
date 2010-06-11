module solve

use constants

implicit none

contains

    subroutine solve_lsq ()

        use aorsa2din_mod, &
        only: square, magma
        use mat_fill, &
        only: aMat!, aMat_
        use antenna, &
        only: brhs, brhs_global

        implicit none

        integer :: info, nRow, nCol
        integer, allocatable, dimension(:) :: ipiv

        !   lapack lsq solve variables
        
        integer :: M_, N_, NRHS, LDA, LDB, MN_, LWORK, RANK
        real :: RCOND
#ifndef dblprec   
        real, allocatable :: rWork(:)
        complex, allocatable :: work(:)
#else
        real(kind=dbl), allocatable :: rWork(:)
        complex(kind=dbl), allocatable :: work(:)
#endif
        integer, allocatable, dimension(:) :: jpvt

        ! Cuda
        integer(kind=long) :: k1, k2, nb
        integer, external :: magma_get_zgetrf_nb
        integer :: magma_solve, magStat
        integer :: dA_dim

        nRow    = size ( aMat, 2 ) !nPtsX * nPtsY * 3
        nCol    = size ( aMat, 1 ) !nModesX * nModesY * 3

        if(square) then

            allocate ( ipiv(nRow) )
            ipiv = 0
#ifndef dblprec
            call cgesv ( nRow, 1, aMat, nRow, ipiv, brhs, nRow, info )
#else            
            if(magma)then

                !! MAGMA Cuda solve
                !! ----------------

                !nb  = magma_get_zgetrf_nb (nRow)
                !k1  = 32 - mod ( maxVal ( (/nRow,nCol/) ), 32 )
                !k2  = 32 - mod ( nCol, 32 ) 
                !if(k1==32)k1=0
                !if(k2==32)k2=0

                !lWork = nRow * nb
                !dA_dim = ( nRow + k1 )**2 &
                !    + (nCol + k2 ) * nb + 2 * nb**2

                !magStat = magma_solve ( dA_dim, lWork, aMat, ipiv, nRow )

                !call zgetrs ( 'N', nRow, 1, aMat, nRow, ipiv, brhs, nRow, info )

            else 

                ! CPU LAPACK
                ! ----------

                call zgesv ( nRow, 1, aMat, nRow, ipiv, brhs, nRow, info )

            endif

#endif

        else

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
            WORK    = 0
            RWORK   = 0

            write(*,*) shape ( amat ), M_*N_, size(amat)
#ifndef dblprec
            call cgelsy ( M_, N_, NRHS, aMat, LDA, brhs, LDB, JPVT, RCOND, RANK, &
                                WORK, LWORK, RWORK, info )
#else
            call zgelsy ( M_, N_, NRHS, aMat, LDA, brhs, LDB, JPVT, RCOND, RANK, &
                                WORK, LWORK, RWORK, info )
#endif

        endif

        if  (info/=0) then

            write(*,*) '    solve.F90: ERROR, solve did not complete successfully'
            write(*,*) '    info: ', info
            stop

        else
            write(*,*) '    LAPACK status: ', info
        endif 
        

        brhs_global = brhs

    end subroutine solve_lsq

#ifdef par

    subroutine solve_lsq_parallel ()

        use aorsa2din_mod, &
        only: npRow, npCol, square
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
#ifndef dblprec
        complex, allocatable :: work(:)
#else
        complex(kind=dbl), allocatable :: work(:)
#endif
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
            ipiv = 0
#ifndef dblprec
            call pcgesv ( n, nrhs, aMat, ia, ja, descriptor_aMat, ipiv, &
                    brhs, ib, jb, descriptor_brhs, info ) 
#else
            call pzgesv ( n, nrhs, aMat, ia, ja, descriptor_aMat, ipiv, &
                    brhs, ib, jb, descriptor_brhs, info ) 
#endif
            if (iAm==0) then 
                write(*,*) '    pcgesv status: ', info
            endif

            deallocate ( ipiv )

        else 

            if(iAm==0) write(*,*) &
            "   using least squares for n x m system"

            allocate ( work ( lWork ) )
            work = 0

#ifndef dblprec
            call pcgels ( trans, m, n, nrhs, aMat, ia, ja, descriptor_aMat, &
                        brhs, ib, jb, descriptor_brhs, work, lWork, info ) 
#else
            call pzgels ( trans, m, n, nrhs, aMat, ia, ja, descriptor_aMat, &
                        brhs, ib, jb, descriptor_brhs, work, lWork, info ) 
#endif
            if (iAm==0) then 
                write(*,*) '    pcgels status: ', info
                write(*,*) '    actual/optimal lwork: ', lwork, real(work(1))
            endif

            deallocate ( work )

        endif 
     
        if  (info/=0) then

            write(*,*) 'solve.F90: ERROR, solve did not complete successfully'
            stop

        endif 
        
        call gather_coeffs ()

    end subroutine solve_lsq_parallel


    subroutine gather_coeffs ()

        use aorsa2din_mod, &
        only: npRow, npCol
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

                !write(*,*) iAm, brhs(ii), brhs_global(gg)

            enddo
        endif

        call cgSum2D ( iContext, 'All', ' ', nRow, 1, &
                brhs_global, nRow, -1, -1 )

    end subroutine gather_coeffs

#endif

    subroutine extract_coeffs ( g )

        use grid
        use antenna, &
        only: brhs_global
 
        implicit none

        type(gridBlock), intent(inout) :: g

        integer :: m, n, iCol

        !   Extract the k coefficents from the solution
        !   for each field component
        !   -------------------------------------------

        allocate ( &
            g%ealphak(g%nMin:g%nMax,g%mMin:g%mMax), &
            g%ebetak(g%nMin:g%nMax,g%mMin:g%mMax), &
            g%eBk(g%nMin:g%nMax,g%mMin:g%mMax) )

        do n=g%nMin,g%nMax
            do m=g%mMin,g%mMax

                iCol = (m-g%mMin) * 3 + (n-g%nMin) * g%nModesZ * 3 + 1
                iCol = iCol + ( g%StartCol-1 )
       
                g%ealphak(n,m)    = brhs_global(iCol+0)
                g%ebetak(n,m)     = brhs_global(iCol+1)
                g%eBk(n,m)        = brhs_global(iCol+2)

            enddo
        enddo

    end subroutine extract_coeffs

end module solve
