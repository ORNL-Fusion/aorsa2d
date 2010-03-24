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
        only: brhs

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

    end subroutine solve_lsq


    subroutine solve_lsq_parallel ()

        use aorsa2din_mod, &
        only: nPtsX, nPtsY, nModesX, nModesY
        use mat_fill, &
        only: aMat
        use antenna, &
        only: brhs
        use parallel

        implicit none

        character(len=1) :: trans
        integer :: m, n, nrhs
        complex :: a
        integer :: ia, ja
        complex :: b
        integer :: ib, jb
        complex, allocatable :: work(:)
        integer :: lWork, info

        trans   = 'N'

        allocate ( work ( lWork ) )

        call pcgels ( trans, m, n, nrhs, a, ia, ja, descriptor_aMat, &
                    b, ib, jb, descriptor_brhs, work, lWork, info ) 

        deallocate ( work )

    end subroutine solve_lsq_parallel


    subroutine extract_coeffs ()

        use grid, &
        only: kxL, kxR, kyL, kyR
        use aorsa2din_mod, &
        only: nPtsX, nPtsY, nModesX, nModesY
        use antenna, &
        only: brhs
 
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
        
                ealphak(n,m)    = brhs(iRow)
                ebetak(n,m)     = brhs(iRow+1)
                eBk(n,m)        = brhs(iRow+2)

            enddo
        enddo

    end subroutine extract_coeffs

end module solve
