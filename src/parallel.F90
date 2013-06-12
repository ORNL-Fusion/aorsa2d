module parallel

use constants, only: long

implicit none

integer :: iContext, nProcs, iAm, myRow, myCol
integer :: NRHS

! Global 
integer :: M_A, N_A, MB_A, NB_A, IA, JA, LLD_A, RSRC_A, CSRC_A
integer :: M_B, N_B, MB_B, NB_B, IB, JB, LLD_B, RSRC_B, CSRC_B

! Local
integer :: LM_A, LN_A
integer :: LM_B, LN_B

! Descriptors
integer :: desc_A(9), desc_B(9)

contains

#ifdef par

    subroutine init_procGrid ( nPts_tot )

        use aorsaNamelist, &
            only: npRow, npCol

        implicit none

        integer, intent(in) :: nPts_tot
        integer :: rowBlockSize, colBlockSize, rhsBlockSize
        integer :: rowStartProc, colStartProc
        integer :: nRow, nCol, nRowLocal, nColLocal
        integer :: nRHSLocal
        integer, external :: numRoC
        integer :: info

        if(huge(nPts_tot) < 2147483647)then
                write(*,*) 'ERROR: 32 bit machine?'
                stop
        endif

        nRow    = nPts_tot*3
        nCol    = nPts_tot*3

        M_A = nRow
        N_A = nCol

        M_B = nCol
        N_B = NRHS

        ! block sizes

        ! Keep this a multiple of 3 such that each 3 rows of a spatial point
        ! is confined to a single processor

        rowBlockSize = 3*32
        colBlockSize = 3*32
        rhsBlockSize = 32

        MB_A = rowBlockSize
        NB_A = colBlockSize

        MB_B = colBlockSize
        NB_B = rhsBlocksize

        ! proc grid starting locations

        rowStartProc    = 0
        colStartProc    = 0

        RSRC_A = rowStartProc
        CSRC_A = colStartProc

        RSRC_B = 0
        CSRC_B = 0

        IA = 1
        JA = 1

        IB = 1
        JB = 1

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
        nRHSLocal   = numRoC ( NRHS, rhsBlockSize, myCol, colStartProc, npCol )

        LM_A = NumRoc(M_A,MB_A,myRow,RSRC_A,npRow) 
        LN_A = NumRoc(N_A,NB_A,myCol,CSRC_A,npCol) 

        LM_B = NumRoc(M_B,MB_B,myRow,RSRC_B,npRow) 
        LN_B = NumRoc(N_B,NB_B,myCol,CSRC_B,npCol) 

        LLD_A = max(1,numroc(M_A,MB_A,myRow,RSRC_A,npRow))
        LLD_B = max(1,numroc(M_B,MB_B,myRow,RSRC_B,npRow))

        call descInit ( desc_A, M_A, N_A, MB_A, NB_A, RSRC_A, CSRC_A, iContext, LLD_A, info )
        if (info.ne.0) then 
            write(*,*) 'ERROR: init desc amat status: ', info
            stop
        endif

        call descInit ( desc_B, M_B, N_B, MB_B, NB_B, RSRC_B, CSRC_B, iContext, LLD_B, info )
        if (info.ne.0) then 
            write(*,*) 'ERROR: init desc brhs status: ', info
            write(*,*) nCol, NRHS, colBlockSize, rowBlockSize, rowStartProc, colStartProc
            write(*,*) desc_B
            stop
        endif


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


!    subroutine init_parallel_aMat (NRHS)
!
!        use aorsaNamelist, &
!            only: npRow, npCol
!        use scalapack_mod, only: numroc
!
!        implicit none
!
!        integer, intent(in) :: NRHS
!        integer :: info
!
!        call descInit ( descriptor_aMat, &
!            nRow, nCol, rowBlockSize, colBlockSize, &
!            rowStartProc, colStartProc, iContext, lld, info )
!
!        if (info.ne.0) then 
!            write(*,*) 'ERROR: init desc amat status: ', info
!            stop
!        endif
!
!        call descInit ( descriptor_brhs, &
!            nCol, NRHS, colBlockSize, rhsBlockSize, rowStartProc, colStartProc, &
!            iContext, lld, info )
!
!        if (info.ne.0) then 
!            write(*,*) 'ERROR: init desc brhs status: ', info
!            write(*,*) nCol, NRHS, colBlockSize, rowBlockSize, rowStartProc, colStartProc
!            write(*,*) descriptor_brhs
!            stop
!        endif
!
!
!    end subroutine init_parallel_aMat 


    subroutine release_grid ()

        implicit none

        call blacs_gridExit ( iContext )
        call blacs_exit ( 0 )

    end subroutine release_grid

#endif

end module parallel
