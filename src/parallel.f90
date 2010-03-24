module parallel

implicit none

integer :: iContext, nProcs, iAm, myRow, myCol
integer :: descriptor_aMat( 50 ), descriptor_brhs( 50 )
integer :: rowStartProc, colStartProc
integer :: rowBlockSize, colBlockSize
integer :: nRow, nCol, nRowLocal, nColLocal

contains

    subroutine init_procGrid ()

        use aorsa2din_mod, &
        only: npRow, npCol, nPtsX, nPtsY, nModesX, nModesY

        implicit none

        nRow    = nPtsX*nPtsY*3
        nCol    = nModesX*nModesY*3

        ! block sizes

        rowBlockSize  = nRow / npRow
        colBlockSize  = nCol / npCol

        ! proc grid starting locations

        rowStartProc    = 0
        colStartProc    = 0

        ! local array sizes

        nRowLocal   = rowBlockSize + mod ( nRow, rowBlockSize )
        nColLocal   = colBlockSize + mod ( nCol, colBlockSize )

        call blacs_pInfo ( iAm, nProcs )

        write(*,*) 'proc id: ', iAm, nProcs

        if ( nProcs <= 1 ) then 

            write(*,*) 'blacs_pInfo backup needed'

            if ( iAm == 0 ) nProcs = npRow * npCol
            call blacs_setup ( iAm, nProcs )

        endif
        write(*,*) npRow, npCol
        call blacs_get ( -1, 0, iContext )
        call blacs_gridInit ( iContext, 'Row-major', npRow, npCol )
        call blacs_gridInfo ( iContext, npRow, npCol, myRow, myCol )

        if ( myRow == -1 ) then 

            write(*,*) 'init_procGrid: ERROR - some sort of problem?'
            stop

        endif

        write(*,*) 'Parallel info:'
        write(*,*) '--------------'
        write(*,*) 'nRow, nCol: ', nRow, nCol
        write(*,*) 'npRow, npCol: ', npRow, npCol
        write(*,*) 'block size: ', rowBlockSize, colBlockSize 
        write(*,*) 'local amat size: ', nRowLocal, nColLocal

    end subroutine init_procGrid


    subroutine init_parallel_aMat ()

        use aorsa2din_mod, &
        only: npRow, npCol, nPtsX, nPtsY, nModesX, nModesY

        implicit none

        integer :: iRowSrc, iColSrc, info, lld

        iRowSrc = 0
        iColSrc = 0

        lld = nRowLocal 

        call descInit ( descriptor_aMat, &
            nRow, nCol, rowBlockSize, colBlockSize, &
            iRowSrc, iColSrc, iContext, lld, info )

        write(*,*) 'init desc amat status: ', info

        call descInit ( descriptor_brhs, &
            nCol, 1, colBlockSize, 1, iRowSrc, iColSrc, &
            iContext, lld, info )

        write(*,*) 'init desc brhs status: ', info

    end subroutine init_parallel_aMat 


end module parallel
