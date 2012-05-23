program aorsa2dMain
   
    use constants
    use eqdsk_dlg
    use aorsaNamelist
    use inv_fourier
    use write_data
    use grid
    use bField
    use interp
    use profiles
    use rotation
    use mat_fill
    use antenna
    use solve
    use parallel
    use timer_mod
    use E_to_lab
    use power
    use setMetal
    !use sigmaInputGeneration
    use ar2Input
    use Performance

    implicit none

#ifdef usepapi
    include "f90papi.h"
#endif
    !include "gptl.inc"

    !   Timer variables

    type ( timer ) :: tFill, tSolve, tTotal, tWorkList, &
        tCurrent, tRotate, tJDotE, tWriteSolution, tCreateGrids, &
        tBfield, tSetMetal, tProfiles, tRotMat, tRHS, tStartup, &
        tOffset, tAllocAMat, tLabel, tSumFlops1, tInitDescriptors, &
        tExtractCoeffs, tInvFourier, tSumFlops2, tSumFlops3, &
        tReleaseGrid

#ifdef usepapi

    !   PAPI variables

    integer :: papi_irc
    real :: papi_rtime, papi_ptime, papi_mflops
    real :: papi_mflops_global
    real :: papi_rtime_zero, papi_ptime_zero, &
        papi_rTime_ProgramZero, papi_pTime_ProgramZero
    real :: papi_rtime_fill, papi_rtime_solve, papi_rtime_current, papi_rtime_total
    real :: papi_ptime_fill, papi_ptime_solve, papi_ptime_current, papi_ptime_total
    integer(kind=long) :: papi_flpins


#endif

    integer(kind=long) :: nPts_tot
    integer :: i
    integer(kind=long) :: nModesR_tot, nModesZ_tot

    !! GPTL vars

    !integer :: stat

    !stat = gptlSetUTR ( GPTLpapitime)
    !write(*,*) 'gptl papiTime stat: ', stat
    !stat = gptlInitialize ()
    !stat = gptlStart ( 'total' )

    call start_timer ( tTotal )
    call start_timer ( tStartup )

    Mem%nAllocations = 0
    Mem%MBytes = 0

#ifdef usepapi
    call PAPIF_flops ( papi_rTime, papi_pTime, papi_flpins, papi_mflops, papi_irc )
    papi_rTime_ProgramZero = papi_rTime
    papi_pTime_ProgramZero = papi_pTime
#endif

    
!   read namelist input data
!   ------------------------

    call read_nameList ()


!   initialise the parallel env
!   ---------------------------

    nPts_tot = sum ( nRAll(1:nGrid)*nZAll(1:nGrid) )

#ifdef par
    call init_procGrid ( nPts_tot )
#endif

#if __noU__==1
    if(iAm==0)write(*,*)'Am using noU ==',__noU__
#else
    if(iAm==0)write(*,*)'Not using noU ==',__noU__
#endif

#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif

    if(iAm==0) write(*,*) '    Time for startup: ', end_timer ( tStartup ),  'seconds'


!   initialise the spatial grid
!   ---------------------------

   
    call start_timer ( tCreateGrids )
    if (iAm==0) &
    write(*,*) 'Creating spatial grid'

    allocate ( allGrids ( nGrid ) )

    do i=1,nGrid 

        allGrids(i) = init_GridBlock ( &
            nRAll(i), nZAll(i), &
            rMinAll(i), rMaxAll(i), &
            zMinAll(i), zMaxAll(i) )

        allGrids(i)%gridNumber = i

    enddo

#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif

    if(iAm==0) write(*,*) '    Time to create grids: ', end_timer ( tCreateGrids ),  'seconds'

    if(iAm==0)then 
        write(*,*) '    nAllocations: ', Mem%nAllocations
        write(*,*) '    Mem Allocated [MBytes]: ', Mem%MBytes
    endif


!   Create workLists
!   ----------------

    call start_timer ( tWorkList )

    do i=1,nGrid

        call createWorkList ( allGrids(i) )

    enddo

#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif
    if (iAm==0) then
        Perf%TimeWorkList = end_timer ( tWorkList )
        write(*,*) '    Time to create WorkList: ', Perf%TimeWorkList
    endif

!   setup magnetic field 
!   --------------------

    call start_timer ( tBfield )

    if (iAm==0) &
    write(*,*) 'Reading magnetic field data'

    if(useEqdsk)then
        !call read_geqdsk ( eqdsk, plot = .false. )
        !call init_interp ()
        stop
    endif

    if(useAr2Input)then
        if(iAm==0) write(*,*) '    Reading ar2 input file ...'
        call ReadAr2Input (AR2InputFileName)
        if(iAm==0) write(*,*) '    DONE'
        if(iAm==0) write(*,*) '    Initializing bField interpolations ...'
        call init_interp ()
        if(iAm==0) write(*,*) '    DONE'
    endif

    do i=1,nGrid

        if (useEqdsk) then
            !call bFieldEqdsk ( allGrids(i) )
            stop
        elseif(useAR2Input)then
            if(iAm==0) write(*,*) '    Interpolating bField to grid ...'
            call bFieldAR2 ( allGrids(i) )
            if(iAm==0) write(*,*) '    DONE' 
        elseif (useSoloviev) then
            !call soloviev ( allGrids(i) )
            stop
        elseif (useCircular) then
            !call bFieldCircular ( allGrids(i) )
            stop
        else
            !call bFieldAnalytical ( allGrids(i) )
            stop
        endif

    enddo
#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif

    if(iAm==0) write(*,*) '    Time to generate bField: ', end_timer ( tBfield ),  'seconds'


!   define the metal regions
!   ------------------------

    call start_timer ( tSetMetal )

    do i=1,nGrid
        if(iAm==0) write(*,*) '    Setting metal regions ...'
        call setMetalRegions( allGrids(i) )
        if(iAm==0) write(*,*) '    DONE'
    enddo
#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif

    if(iAm==0) write(*,*) '    Time to set metal regions: ', end_timer ( tSetMetal ),  'seconds'

!   setup profiles
!   --------------

    call start_timer ( tProfiles )

    if (iAm==0) &
    write(*,*) 'Profile setup'
  
    if(useAR2Input)then 

        do i=1,nGrid
            call init_ar2_profiles( allGrids(i) )
        enddo

    else 

        !call init_profiles (nSpec)

        !do i=1,nGrid

        !    if (useFluxProfiles) then
        !        call flux_profiles ( allGrids(i) )
        !    elseif (useCircularProfiles) then
        !        call circular_profiles ( allGrids(i) )
        !    else
        !        call flat_profiles ( allGrids(i), parabolic = parabolic )
        !    endif

        !enddo
        stop
    endif

#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif

    if(iAm==0) write(*,*) '    Time to generate profiles: ', end_timer ( tProfiles ),  'seconds'


!!   Setup the plasma paraemter splines for sigma input
!!   --------------------------------------------------
!
!    do i=1,nGrid
!
!        call setupSigmaParameterSplines ( allGrids(i) )
!
!    enddo


!   calculate rotation matrix U
!   ---------------------------

    call start_timer ( tRotMat )

    if (iAm==0) &
    write(*,*) 'Building rotation matrix U'

    do i=1,nGrid

        if(iAm==0)write(*,*) '    Initializing rotation matrix ...'
        call init_rotation ( allGrids(i) )
        if(iAm==0)write(*,*) '    DONE'
        if(iAm==0)write(*,*) '    Generating derivatives of the rotation matrix ...'
        call deriv_rotation ( allGrids(i) )
        if(iAm==0)write(*,*) '    DONE'

    enddo

#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif

    if(iAm==0) write(*,*) '    Time to build rotation matrix: ', end_timer ( tRotmat ),  'seconds'


!   Calculate aMat index for first pt in each grid block
!   ----------------------------------------------------

    call start_timer ( tOffset )

    do i=1,nGrid

        if(i==1)then
            allGrids(i)%startRow = 1 
            allGrids(i)%startCol = 1
        else
            allGrids(i)%startRow = sum ( allGrids(1:i-1)%nR * allGrids(1:i-1)%nZ * 3 ) + 1
            allGrids(i)%startCol = sum ( allGrids(1:i-1)%nModesR * allGrids(1:i-1)%nModesZ * 3 ) + 1
        endif

    enddo

#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif

    if(iAm==0) write(*,*) '    Time to generate aMat offsets: ', end_timer ( tOffset ),  'seconds'


!   Antenna current
!   ---------------

    call start_timer ( tRHS )

    if (iAm==0) &
    write(*,*) 'Building antenna current (brhs)'

    call alloc_total_brhs ( nPts_tot )
    do i=1,nGrid
        call init_brhs ( allGrids(i) )
    enddo

#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif
    if(iAm==0) write(*,*) '    Time to build RHS: ', end_timer ( tRHS ),  'seconds'


!!   Write the run input data to disk
!!   --------------------------------
!
!    if (iAm==0) &
!    write(*,*) 'Writing run input data to file'
!
!    if (iAm==0) then
!
!        do i=1,nGrid
!            write(allGrids(i)%fNumber,'(i3.3)') i
!            !call write_runData ( allGrids(i) )
!        enddo
!
!    endif


!   Allocate the matrix
!   -------------------

    call start_timer ( tAllocAMat )
    call alloc_total_aMat ( nPts_tot )
#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif
    if(iAm==0) write(*,*) '    Time to allocate aMat: ', end_timer ( tAllocAMat ),  'seconds'



!   Label each grid blocks boundary points
!   --------------------------------------

    call start_timer ( tLabel )
    call labelPts ( allGrids, nPts_tot )
#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif
    if(iAm==0) write(*,*) '    Time to label points: ', end_timer ( tLabel ),  'seconds'


!   Fill matrix 
!   -----------

    call start_timer ( tFill )

#ifdef usepapi
    call PAPIF_flops ( papi_rTime, papi_pTime, papi_flpins, papi_mflops, papi_irc )
    papi_rTime_zero = papi_rTime
    papi_pTime_zero = papi_pTime
#endif

    do i=1,nGrid

        call aMat_fill ( allGrids(i) )

    enddo

    call amat_boundaries ( allGrids, nPts_tot )

    !where(abs(aMat)/=0)
    !        aMat = 11
    !endwhere
    !write(*,*) shape(amat)
    
    !write(*,'(3x,48(1x,i3.2))') (/ (i,i=1,48) /)
    !write(*,*) '   --------------------------------------------------'
    !do i=1,48
    !write(*,'(i2.2,x,48(1x,i3.2),4x,i2.2)') i,int((aMat(i,:))),int((brhs(i)))
    !enddo

    !do i=1,nGrid
    !    write(*,*) allGrids(i)%label
    !enddo

#ifdef usepapi
    call PAPIF_flops ( papi_rTime, papi_pTime, papi_flpins, papi_mflops, papi_irc )
    papi_rTime_fill = papi_rTime - papi_rTime_zero
    papi_pTime_fill = papi_pTime - papi_pTime_zero
#endif

#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif
    if (iAm==0) then 
        Perf%TimeFill = end_timer ( tFill )
        write(*,*) '    Time to fill aMat: ', Perf%TimeFill
    endif


    call start_timer ( tSumFlops1 )

#ifdef usepapi
    !   sum total mflops from all procs

    papi_mflops_global  = papi_mflops

#ifdef par
    call sgSum2D ( iContext, 'All', ' ', 1, 1, &
                papi_mflops_global, 1, -1, -1 )
#endif

    if (iAm==0) then

        write(*,*) '    PAPI data for fill'
        write(*,'(a30,f12.1)') 'real time: ', papi_rTime_fill
        write(*,'(a30,f12.1)') 'proc time: ', papi_pTime_fill
        write(*,'(a30,f12.2)') 'iAm Gflops/s: ', papi_mflops / 1e3
        write(*,'(a30,f12.2)') 'global Gflops/s: ', papi_mflops_global / 1e3
        write(*,'(a30,i1)') 'status: ', papi_irc

        Perf%GflopsFillLocal = papi_mflops / 1e3
        Perf%GflopsFillGlobal = papi_mflops_global / 1e3

    endif

#endif

#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif
    if(iAm==0) write(*,*) '    Time to sum flops 1: ', end_timer ( tSumFlops1 ),  'seconds'


    !if (iAm==0) &
    !call write_amat ( 'amat.nc' )


!   Solve complex, linear system
!   ----------------------------

    if (iAm==0) &
    write(*,*) 'Solving complex linear system'

    call start_timer ( tInitDescriptors )

#ifdef par
    call init_parallel_aMat ()
#endif

#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif
    if(iAm==0) write(*,*) '    Time to init descriptors: ', end_timer ( tInitDescriptors ),  'seconds'

    call start_timer ( tSolve )

#ifdef usepapi
    call PAPIF_flops ( papi_rTime, papi_pTime, papi_flpins, papi_mflops, papi_irc )
    papi_rTime_zero = papi_rTime
    papi_pTime_zero = papi_pTime
#endif

#ifdef par
    call solve_lsq_parallel ()
#else
    call solve_lsq ()
#endif

#ifdef usepapi
    call PAPIF_flops ( papi_rTime, papi_pTime, papi_flpins, papi_mflops, papi_irc )
    papi_rTime_solve = papi_rTime - papi_rTime_zero
    papi_pTime_solve = papi_pTime - papi_pTime_zero
#endif

#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif
    if (iAm==0) then 
        Perf%TimeSolve = end_timer ( tSolve )
        write(*,*) '    Time to solve: ', Perf%TimeSolve
    endif

    call start_timer ( tSumFlops2 )

#ifdef usepapi
    !   sum total mflops from all procs

    papi_mflops_global  = papi_mflops

#ifdef par
    call sgSum2D ( iContext, 'All', ' ', 1, 1, &
                papi_mflops_global, 1, -1, -1 )
#endif

    if (iAm==0) then

        write(*,*) '    PAPI data for solve'
        write(*,'(a30,f12.1)') 'real time: ', papi_rTime_solve
        write(*,'(a30,f12.1)') 'proc time: ', papi_pTime_solve
        write(*,'(a30,f12.2)') 'iAm Gflops/s: ', papi_mflops / 1e3
        write(*,'(a30,f12.2)') 'global Gflops/s: ', papi_mflops_global / 1e3
        write(*,'(a30,i1)') 'status: ', papi_irc

        Perf%GflopsSolveLocal = papi_mflops / 1e3
        Perf%GflopsSolveGlobal = papi_mflops_global / 1e3

    endif

#endif

#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif
    if(iAm==0) write(*,*) '    Time to sum flops 2: ', end_timer ( tSumFlops2 ),  'seconds'

    call start_timer ( tExtractCoeffs )
    do i=1,nGrid
        call extract_coeffs ( allGrids(i) )    
    enddo
#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif
    if(iAm==0) write(*,*) '    Time to extract coeffs: ', end_timer ( tExtractCoeffs ),  'seconds'


!   Inverse Fourier transform the k coefficents
!   to real space
!   -------------------------------------------

    if (iAm==0) &
    write(*,*) 'Constructin E solution from basis set (alpha,beta,b)'
    call start_timer ( tInvFourier )   
    do i=1,nGrid 
        call sftInv2d ( allGrids(i) )
    enddo

#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif
    if(iAm==0) write(*,*) '    Time to inverse Fourier transform the fields: ', end_timer ( tInvFourier ),  'seconds'


!   Calculate the plasma current
!   ----------------------------

#ifdef usepapi
    call PAPIF_flops ( papi_rTime, papi_pTime, papi_flpins, papi_mflops, papi_irc )
    papi_rTime_zero = papi_rTime
    papi_pTime_zero = papi_pTime
#endif

    call start_timer ( tCurrent )

    if (iAm==0) &
    write(*,*) 'Calculating plasma current'

    !if ( iAm == 0 ) then

        do i=1,nGrid
            call current ( allGrids(i) )
        enddo

    !endif
#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif
    if (iAm==0) then
        Perf%TimeCurrent = end_timer ( tCurrent )
        write(*,*) 'Total to calculate plasma current: ', Perf%TimeCurrent
    endif

    call start_timer ( tSumFlops3 )
#ifdef usepapi
    call PAPIF_flops ( papi_rTime, papi_pTime, papi_flpins, papi_mflops, papi_irc )
    papi_rTime_current = papi_rTime - papi_rTime_zero
    papi_pTime_current = papi_pTime - papi_pTime_zero
#endif

#ifdef usepapi
    !   sum total mflops from all procs

    papi_mflops_global  = papi_mflops

#ifdef par
    call sgSum2D ( iContext, 'All', ' ', 1, 1, &
                papi_mflops_global, 1, -1, -1 )
#endif

    if (iAm==0) then

        write(*,*) '    PAPI data for plasma current'
        write(*,'(a30,f12.1)') 'real time: ', papi_rTime_current
        write(*,'(a30,f12.1)') 'proc time: ', papi_pTime_current
        write(*,'(a30,f12.2)') 'iAm Gflops/s: ', papi_mflops / 1e3
        write(*,'(a30,f12.2)') 'global Gflops/s: ', papi_mflops_global / 1e3
        write(*,'(a30,i1)') 'status: ', papi_irc

        Perf%GflopsCurrentLocal = papi_mflops / 1e3
        Perf%GflopsCurrentGlobal = papi_mflops_global / 1e3

    endif
#endif
#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif
    if(iAm==0) write(*,*) '    Time to sum flops 3: ', end_timer ( tSumFlops3 ),  'seconds'


!   Rotation to Lab frame 
!   ---------------------

    call start_timer ( tRotate )
    if (iAm==0) &
    write(*,*) 'Rotating E solution to lab frame (R,Th,z)'

    if ( iAm == 0 ) then

        do i=1,nGrid 
            call rotate_E_to_lab ( allGrids(i) )
        enddo

    endif
#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif
    if(iAm==0) write(*,*) '    Time to rotate solution to lab frame: ', end_timer ( tCurrent ),  'seconds'

!   Calculate the Joule Heating 
!   ---------------------------

    call start_timer ( tJDotE )
    if (iAm==0) &
    write(*,*) 'Calculating Joule heating' 

    if ( iAm == 0 ) then

        do i=1,nGrid
            call jDotE ( allGrids(i) )
        enddo

    endif
#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif
    if(iAm==0) write(*,*) '    Time to calculate Joule heating: ', end_timer ( tJDotE ),  'seconds'


!   Write soln to file
!   ------------------

    call start_timer ( tWriteSolution )
    if (iAm==0) &
    write(*,*) 'Writing solution to file'

    if ( iAm == 0 ) then

        do i=1,nGrid
            call write_solution ( allGrids(i) )
        enddo

    endif
#ifdef par
    call blacs_barrier ( iContext, 'All' ) 
#endif
    if(iAm==0) write(*,*) '    Time to write solution: ', end_timer ( tWriteSolution ),  'seconds'

#ifdef usepapi
    call PAPIF_flops ( papi_rTime, papi_pTime, papi_flpins, papi_mflops, papi_irc )
    papi_rTime_Total = papi_rTime - papi_rTime_ProgramZero
    papi_pTime_Total = papi_pTime - papi_pTime_ProgramZero
#endif
#ifdef usepapi

    if (iAm==0) then

        write(*,*) '    PAPI data for program'
        write(*,'(a30,f12.1)') 'real time: ', papi_rTime_total
        write(*,'(a30,f12.1)') 'proc time: ', papi_pTime_total
        write(*,'(a30,i1)') 'status: ', papi_irc

    endif
#endif

    call start_timer(tReleaseGrid)
#ifdef par
    call release_grid ()
#endif
    if(iAm==0) write(*,*) '    Time to release mpi grid: ', end_timer ( tReleaseGrid ),  'seconds'

    !stat = gptlStop ('total')
    !stat = gptlpr (0)
    !stat = gptlFinalize ()

    if (iAm==0) then

        Perf%nProcs = npRow*npCol
        Perf%nSpatialPoints = nPts_tot
        Perf%nRowLocal = nRowLocal
        Perf%nColLocal = nColLocal
        Perf%nRowGlobal = nRow
        Perf%nColGlobal = nCol
        Perf%MatSizeLocal_GB = LocalSizeGB
        Perf%MatSizeGlobal_GB = GlobalSizeGB
        Perf%TimeTotal = end_timer ( tTotal )

        write(*,*) '    Total time: ', Perf%TimeTotal

        call WritePerformanceData ( Perf )

    endif

end program aorsa2dMain
