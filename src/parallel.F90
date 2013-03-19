module parallel

use constants, only: long

implicit none

integer :: iContext, nProcs, iAm, myRow, myCol
integer :: descriptor_aMat( 9 ), descriptor_brhs( 9 )
integer :: rowStartProc, colStartProc
integer :: rowBlockSize, colBlockSize
integer(kind=long) :: nRow, nCol, nRowLocal, nColLocal

contains

#ifdef par

    subroutine init_procGrid ( nPts_tot )

        use aorsaNamelist, &
            only: npRow, npCol

        implicit none

        integer(kind=long), intent(in) :: nPts_tot

        integer, external :: numRoC

        nRow    = nPts_tot*3
        nCol    = nPts_tot*3

        ! block sizes

        rowBlockSize  = nRow / npRow
        colBlockSize  = nCol / npCol

        ! Keep this a multiple of 3 such that each 3 rows of a spatial point
        ! is confined to a single processor

        if ( rowBlockSize >= 3*32 ) rowBlockSize = 3*32 
        if ( colBlockSize >= 3*32 ) colBlockSize = 3*32 

        ! proc grid starting locations

        rowStartProc    = 0
        colStartProc    = 0

        call blacs_pInfo ( iAm, nProcs )

        !write(*,*) 'proc id: ', iAm, nProcs

        if ( nProcs <= 1 ) then 

            write(*,*) 'blacs_pInfo backup needed'

            if ( iAm == 0 ) nProcs = npRow * npCol
            call blacs_setup ( iAm, nProcs )

        endif

        if ( nProcs /= npRow*npCol ) then 

            write(*,*)
            write(*,*)
            write(*,*) 'CONFIG ERROR:'
            write(*,*) '-------------'
            write(*,*) '    nProcs /= npRow * npCol'
            write(*,*) '    Please correct and re-run'
            write(*,*) '    Have a nice day :)'
            write(*,*)
            write(*,*)
 
            return  
                    
        endif

        call blacs_get ( -1, 0, iContext )
        call blacs_gridInit ( iContext, 'Row-major', npRow, npCol )
        call blacs_gridInfo ( iContext, npRow, npCol, myRow, myCol )

        ! local array sizes as determined by scalapack routine 
        ! numRoC

        nRowLocal   = numRoC ( nRow, rowBlockSize, myRow, rowStartProc, npRow )
        nColLocal   = numRoC ( nCol, colBlockSize, myCol, colStartProc, npCol )

        if ( myRow == -1 ) then 

            write(*,*) 'init_procGrid: ERROR - some sort of problem?'
            stop

        endif

        if (iAm==0) then 
            write(*,*) 'Parallel info:'
            write(*,*) '    nRow, nCol:      ', nRow, nCol
            write(*,*) '    npRow, npCol:    ', npRow, npCol
            write(*,*) '    myRow, myCol:    ', myRow, myCol
            write(*,*) '    block size:      ', rowBlockSize, colBlockSize 
            write(*,*) '    local amat size: ', nRowLocal, nColLocal
        endif

    end subroutine init_procGrid


    subroutine init_parallel_aMat (NRHS)

        use aorsaNamelist, &
            only: npRow, npCol

        implicit none

        integer, intent(in) :: NRHS
        integer :: info, lld

        lld = nRowLocal 

        call descInit ( descriptor_aMat, &
            nRow, nCol, rowBlockSize, colBlockSize, &
            rowStartProc, colStartProc, iContext, lld, info )

        !write(*,*) 'init desc amat status: ', info

        call descInit ( descriptor_brhs, &
            nRow, NRHS, rowBlockSize, colBlockSize, rowStartProc, colStartProc, &
            iContext, lld, info )

        !write(*,*) 'init desc brhs status: ', info

    end subroutine init_parallel_aMat 


    subroutine release_grid ()

        implicit none

        call blacs_gridExit ( iContext )
        call blacs_exit ( 0 )

    end subroutine release_grid

#endif

end module parallel
