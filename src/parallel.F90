module parallel

implicit none

integer :: iContext, nProcs, iAm, myRow, myCol
integer :: descriptor_aMat( 9 ), descriptor_brhs( 9 )
integer :: rowStartProc, colStartProc
integer :: rowBlockSize, colBlockSize
integer :: nRow, nCol, nRowLocal, nColLocal

contains

#ifdef par

    subroutine init_procGrid ( nPtsR_tot, nPtsZ_tot, nModesR_tot, nModesZ_tot )

        use aorsa2din_mod, &
        only: npRow, npCol

        implicit none

        integer, intent(in) :: nPtsR_tot, nPtsZ_tot, nModesR_tot, nModesZ_tot

        integer, external :: numRoC

        nRow    = nPtsR_tot*nPtsZ_tot*3
        nCol    = nModesR_tot*nModesZ_tot*3

        ! block sizes

        rowBlockSize  = nRow / npRow
        colBlockSize  = nCol / npCol

        if ( rowBlockSize >= 64 ) rowBlockSize = 64 
        if ( colBlockSize >= 64 ) colBlockSize = 64 

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


    subroutine init_parallel_aMat ()

        use aorsa2din_mod, &
        only: npRow, npCol

        implicit none

        integer :: info, lld

        lld = nRowLocal 

        call descInit ( descriptor_aMat, &
            nRow, nCol, rowBlockSize, colBlockSize, &
            rowStartProc, colStartProc, iContext, lld, info )

        !write(*,*) 'init desc amat status: ', info

        call descInit ( descriptor_brhs, &
            nRow, 1, rowBlockSize, colBlockSize, rowStartProc, colStartProc, &
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
