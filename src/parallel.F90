module parallel

use constants, only: long

implicit none

integer :: ICTXT, nProcs, iAm, MYROW, MYCOL
integer :: NRHS

! Global 
integer :: M, N, IA, JA, LLD_A, RSRC, CSRC
integer :: MB, NB, NBRHS, IB, JB, LLD_B

! Local
integer :: LM_A, LN_A
integer :: LM_B, LN_B

! Descriptors
integer :: DESCA(9), DESCB(9)

contains

#ifdef par

    subroutine init_procGrid ( nPts_tot )

        use aorsaNamelist, &
            only: NPROW, NPCOL

        implicit none

        integer, intent(in) :: nPts_tot
        integer :: rowBlockSize, colBlockSize, rhsBlockSize
        integer :: rowStartProc, colStartProc
        integer :: nRow, nCol
        integer :: nRHSLocal
        integer, external :: NUMROC
        integer :: INFO

        if(huge(nPts_tot) < 2147483647)then
                write(*,*) 'ERROR: 32 bit machine?'
                stop
        endif

        nRow    = nPts_tot*3
        nCol    = nPts_tot*3

        M = nRow
        N = nCol

        if(NRHS==0)then
            write(*,*) 'ERROR: NRHS = ', NRHS
            stop
        endif

        ! block sizes

        ! Keep this a multiple of 3 such that each 3 rows of a spatial point
        ! is confined to a single processor

        MB = 3*32
        NB = 3*32
        NBRHS = NB 

        ! proc grid starting locations

        RSRC = 0
        CSRC = 0

        IA = 1
        JA = 1

        IB = 1
        JB = 1

        ! Initialize the processor grid

        call SL_INIT( ICTXT, NPROW, NPCOL )
        call blacs_pINFO ( iAM, NPROCS )
        call blacs_gridINFO ( ICTXT, NPROW, NPCOL, MYROW, MYCOL )

        if ( nProcs /= NPROW*NPCOL ) then 

            write(*,*)
            write(*,*) 'ERROR:'
            write(*,*) '------'
            write(*,*) '    nProcs /= NPROW * NPCOL, please correct and re-run'
            write(*,*) '    nProcs : ', nProcs 
            write(*,*) '    NPROW : ', NPROW 
            write(*,*) '    NPCOL : ', NPCOL 
            return  

        endif

        ! local array sizes as determined by scalapack routine NUMROC

        LM_A = NUMROC( M, MB, MYROW, RSRC, NPROW )
        LN_A = NUMROC( N, NB, MYCOL, CSRC, NPCOL )

        LM_B = NUMROC( N, NB, MYROW, RSRC, NPROW )
        LN_B = NUMROC( NRHS, NBRHS, MYCOL, CSRC, NPCOL )

        LLD_A = max( 1, LM_A )
        LLD_B = max( 1, LM_B )

        call descInit ( DESCA, M, N, MB, NB, RSRC, CSRC, ICTXT, LLD_A, INFO )
        if (INFO.ne.0) then 
            write(*,*) 'ERROR: init desc amat status: ', INFO
            stop
        endif

        call descInit ( DESCB, N, NRHS, NB, NBRHS, RSRC, CSRC, ICTXT, LLD_B, INFO )
        if (INFO.ne.0) then 
            write(*,*) 'ERROR: init desc brhs status: ', INFO
            write(*,*) nCol, NRHS, colBlockSize, rowBlockSize, rowStartProc, colStartProc
            write(*,*) DESCB
            stop
        endif

        if ( MYROW == -1 ) then 

            write(*,*) 'init_procGrid: ERROR - some sort of problem?'
            stop

        endif

        if (iAm==0) then 
            write(*,*) 'Parallel INFO:'
            write(*,*) '    nRow, nCol:      ', M, N 
            write(*,*) '    NPROW, NPCOL:    ', NPROW, NPCOL
            write(*,*) '    MYROW, MYCOL:    ', MYROW, MYCOL
            write(*,*) '    block size:      ', MB, NB 
        endif

    end subroutine init_procGrid


    subroutine release_grid ()

        implicit none

        call blacs_gridExit ( ICTXT )
        call blacs_exit ( 0 )

    end subroutine release_grid

#endif

end module parallel
