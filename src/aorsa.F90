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
    use ar2Input, ar2_rMin=>rMin, ar2_rMax=>rMax, ar2_zMin=>zMin, ar2_zMax=>zMax
    use AR2SourceLocationsInput, NRHS_FromInputFile=>NRHS
    use Performance
    use read_jp_from_file
    use zfunction_mod

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

    integer :: nPts_tot
    integer :: i, rhs

    integer :: length_rid, status_rid
    character(len=20) :: rid
    character(len=80) :: InputFileName

    logical :: exists

    if(huge(nPts_tot) < 2147483647)then
            write(*,*) 'ERROR: 32 bit machine?'
            stop
    endif


    ! Check if the "output" directory exists
    ! Future: implement the "modFileSys" lib.
    inquire(file='output/.',exist=exists)
    if(.not.exists)then
        call system("mkdir output")
        inquire(file='output/.',exist=exists)
        if(.not.exists)then
            write(*,*) 'ERROR: directory "output/" does not exist'
            stop
        endif
    endif

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

#ifdef USE_GPU
    call cublasinit()
    call profinit()
#endif

!   Get the run id from an env variable
!   -----------------------------------

    call get_environment_variable('AORSA_RID',rid,length_rid,status_rid)

    if(status_rid /= 0)then
        !if(iAm==0)write(*,*) 'AORSA_RID not found with status: ', status_rid
        rid = ''
    endif

!   read namelist input data
!   ------------------------

    InputFileName = 'aorsa2d.in'
    call read_nameList (trim(rid)//trim(InputFileName))

!   read ar2Input file
!   ------------------

    if(useAr2Input)then
        if(iAm==0) write(*,*) '    Reading ar2 input file ...'
        call ReadAr2Input (AR2InputFileName)
        if(iAm==0) write(*,*) '    DONE'
        if(iAm==0) write(*,*) '    Initializing bField interpolations ...'
        call init_interp ()
        if(iAm==0) write(*,*) '    DONE'
    endif


    if(useAR2SourceLocationsFile)then
        call ReadAR2SourceLocations(AR2SourceLocationsFileName) 
        NRHS = NRHS_FromInputFile
    else 
        NRHS = 1
    endif

!   initialise the parallel env
!   ---------------------------

    nPts_tot = sum ( nRAll(1:nGrid)*nZAll(1:nGrid) )

#ifdef par
    call init_procGrid ( nPts_tot )
#endif

    if(iAm==0) &
    write(*,*) 'AORSA_RID: ', rid

#if __noU__==1
    if(iAm==0) &
    write(*,*)'Am using noU ==',__noU__
#else
    if(iAm==0) &
    write(*,*)'Not using noU ==',__noU__
#endif

#ifdef par
    call blacs_barrier ( ICTXT, 'All' ) 
#endif

    if(iAm==0) &
    write(*,*) '    Time for startup: ', end_timer ( tStartup ),  'seconds'


!   initialise the spatial grid
!   ---------------------------

   
    call start_timer ( tCreateGrids )
    if (iAm==0) &
    write(*,*) 'Creating spatial grid'

    allocate ( allGrids ( nGrid ) )

    do i=1,nGrid 

        if(useAR2Input)then
    
            rMinAll(i) = ar2_rMin
            rMaxAll(i) = ar2_rMax
            zMinAll(i) = ar2_zMin
            zMaxAll(i) = ar2_zMax

            if(iAm==0)then
                write(*,*) 'rMin: ', rMinAll(i)
                write(*,*) 'rMax: ', rMaxAll(i)
                write(*,*) 'zMin: ', zMinAll(i)
                write(*,*) 'zMax: ', zMaxAll(i)
            endif

        endif

        allGrids(i) = init_GridBlock ( &
            nRAll(i), nZAll(i), &
            rMinAll(i), rMaxAll(i), &
            zMinAll(i), zMaxAll(i) )

        allGrids(i)%gridNumber = i
        write(allGrids(i)%fNumber,'(i3.3)') i
    enddo

#ifdef par
    call blacs_barrier ( ICTXT, 'All' ) 
#endif

    if(iAm==0) &
    write(*,*) '    Time to create grids: ', end_timer ( tCreateGrids ),  'seconds'

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
    call blacs_barrier ( ICTXT, 'All' ) 
#endif

    Perf%TimeWorkList = end_timer ( tWorkList )
    if(iAm==0) &
    write(*,*) '    Time to create WorkList: ', Perf%TimeWorkList

!   setup magnetic field 
!   --------------------

    call start_timer ( tBfield )

    if (iAm==0) &
    write(*,*) 'Reading magnetic field data'

    do i=1,nGrid

        if(useAR2Input)then
            if(iAm==0) write(*,*) '    Interpolating bField to grid ...'
            call bFieldAR2 ( allGrids(i) )
            if(iAm==0) write(*,*) '    DONE' 
        endif

    enddo
#ifdef par
    call blacs_barrier ( ICTXT, 'All' ) 
#endif

    if(iAm==0) &
    write(*,*) '    Time to generate bField: ', end_timer ( tBfield ),  'seconds'


!   define the metal regions
!   ------------------------

    call start_timer ( tSetMetal )

    do i=1,nGrid
        if(iAm==0) write(*,*) '    Setting metal regions ...'
        call setMetalRegions( allGrids(i) )
        if(iAm==0) write(*,*) '    DONE'
    enddo
#ifdef par
    call blacs_barrier ( ICTXT, 'All' ) 
#endif

    if(iAm==0) &
    write(*,*) '    Time to set metal regions: ', end_timer ( tSetMetal ),  'seconds'

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
        stop
    endif

!   load Z function file 
!   --------------------

    if (iAm==0) &
    write(*,*) 'Loading Z function file'
  
    call z_load_table(zFunctionFileName)

    if(iAm==0) write(*,*) '    DONE'

#ifdef par
    call blacs_barrier ( ICTXT, 'All' ) 
#endif

    if(iAm==0) &
    write(*,*) '    Time to generate profiles: ', end_timer ( tProfiles ),  'seconds'


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
    call blacs_barrier ( ICTXT, 'All' ) 
#endif

    if(iAm==0) &
    write(*,*) '    Time to build rotation matrix: ', end_timer ( tRotmat ),  'seconds'


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
    call blacs_barrier ( ICTXT, 'All' ) 
#endif

    if(iAm==0) &
    write(*,*) '    Time to generate aMat offsets: ', end_timer ( tOffset ),  'seconds'


!   Load plasma current Jp from file if working with Kinetic-J
!   ----------------------------------------------------------

    if(useJpFromFile)then
        call ReadJpFromFile(JpInputFileName)    
        do i=1,nGrid
            call init_JpFromFile( allGrids(i) )
        enddo
    endif


!   Antenna current
!   ---------------

    call start_timer ( tRHS )

    if (iAm==0) &
    write(*,*) 'Building antenna current (brhs)'

    if(iAm==0) &
    write(*,*) 'Initialized jA'

    RHS_PreProcessing_loop: &
    do rhs=1,NRHS

#ifdef par
        call blacs_barrier ( ICTXT, 'All' ) 
#endif
    
        do i=1,nGrid
            call fill_jA ( allGrids(i), rhs )
        enddo

#ifdef par
        call blacs_barrier ( ICTXT, 'All' ) 
#endif
   
    enddo RHS_PreProcessing_loop

    call alloc_total_brhs ( nPts_tot )

    if(iAm==0) &
    write(*,*) 'Allocated B'

    do i=1,nGrid
        call init_brhs ( allGrids(i), nPts_tot )
    enddo

    if(iAm==0.and.NRHS<=1) &
    write(*,*) '    Time to build RHS: ', end_timer ( tRHS ),  'seconds'

    do rhs=1,NRHS

#ifdef par
        call blacs_barrier ( ICTXT, 'All' ) 
#endif
    
        if (iAm==0.and.NRHS<=1) &
        write(*,*) 'Writing run input data to file'
    
        do i=1,nGrid
            call write_runData ( allGrids(i), rid, rhs )
        enddo

#ifdef par
        call blacs_barrier ( ICTXT, 'All' ) 
#endif
   
    enddo 


!   Allocate the matrix
!   -------------------

    call start_timer ( tAllocAMat )
    call alloc_total_aMat ( nPts_tot )
#ifdef par
    call blacs_barrier ( ICTXT, 'All' ) 
#endif
    if(iAm==0) &
    write(*,*) '    Time to allocate aMat: ', end_timer ( tAllocAMat ),  'seconds'


!   Label each grid blocks boundary points
!   --------------------------------------

    call start_timer ( tLabel )
    call labelPts ( allGrids, nPts_tot )
#ifdef par
    call blacs_barrier ( ICTXT, 'All' ) 
#endif
    if(iAm==0) &
    write(*,*) '    Time to label points: ', end_timer ( tLabel ),  'seconds'


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

#ifdef usepapi
    call PAPIF_flops ( papi_rTime, papi_pTime, papi_flpins, papi_mflops, papi_irc )
    papi_rTime_fill = papi_rTime - papi_rTime_zero
    papi_pTime_fill = papi_pTime - papi_pTime_zero
#endif

#ifdef par
    call blacs_barrier ( ICTXT, 'All' ) 
#endif

    Perf%TimeFill = end_timer ( tFill )
    if (iAm==0) &
    write(*,*) '    Time to fill aMat: ', Perf%TimeFill

    call start_timer ( tSumFlops1 )

#ifdef usepapi
    !   sum total mflops from all procs

    papi_mflops_global  = papi_mflops

#ifdef par
    call sgSum2D ( ICTXT, 'All', ' ', 1, 1, &
                papi_mflops_global, 1, -1, -1 )
#endif

    Perf%GflopsFillLocal = papi_mflops / 1e3
    Perf%GflopsFillGlobal = papi_mflops_global / 1e3

    if (iAm==0) then

        write(*,*) '    PAPI data for fill'
        write(*,'(a30,f12.1)') 'real time: ', papi_rTime_fill
        write(*,'(a30,f12.1)') 'proc time: ', papi_pTime_fill
        write(*,'(a30,f12.2)') 'iAm Gflops/s: ', papi_mflops / 1e3
        write(*,'(a30,f12.2)') 'global Gflops/s: ', papi_mflops_global / 1e3
        write(*,'(a30,i1)') 'status: ', papi_irc

    endif

#endif

#ifdef par
    call blacs_barrier ( ICTXT, 'All' ) 
#endif
    if(iAm==0) &
    write(*,*) '    Time to sum flops 1: ', end_timer ( tSumFlops1 ),  'seconds'

    !if (iAm==0) &
    !call write_amat ( 'amat.nc' )


!   Solve complex, linear system
!   ----------------------------

    if (iAm==0) &
    write(*,*) 'Solving complex linear system'

    call start_timer ( tInitDescriptors )

#ifdef par
!    call init_parallel_aMat ( NRHS )
#endif

#ifdef par
    call blacs_barrier ( ICTXT, 'All' ) 
#endif
    if(iAm==0) &
    write(*,*) '    Time to init descriptors: ', end_timer ( tInitDescriptors ),  'seconds'

    call start_timer ( tSolve )
    if(iAm==0) &
    write(*,*) '  tSolve (post init): ', tSolve%time

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
    call blacs_barrier ( ICTXT, 'All' ) 
#endif

    Perf%TimeSolve = end_timer ( tSolve )
    if (iAm==0) &
    write(*,*) '    Time to solve: ', Perf%TimeSolve

    call start_timer ( tSumFlops2 )

#ifdef usepapi
    !   sum total mflops from all procs

    papi_mflops_global  = papi_mflops

#ifdef par
    call sgSum2D ( ICTXT, 'All', ' ', 1, 1, &
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
    call blacs_barrier ( ICTXT, 'All' ) 
#endif
    if(iAm==0) &
    write(*,*) '    Time to sum flops 2: ', end_timer ( tSumFlops2 ),  'seconds'

    call start_timer ( tExtractCoeffs )
    do i=1,nGrid
        call extract_coeffs ( allGrids(i) )    
    enddo
#ifdef par
    call blacs_barrier ( ICTXT, 'All' ) 
#endif
    if(iAm==0) &
    write(*,*) '    Time to extract coeffs: ', end_timer ( tExtractCoeffs ),  'seconds'

    RHS_PostProcessing_loop: &
    do rhs=1,NRHS  

    !   Inverse Fourier transform the k coefficents
    !   to real space
    !   -------------------------------------------
    
        if (iAm==0.and.NRHS<=1) &
        write(*,*) 'Constructin E solution from basis set (alpha,beta,b)'
        call start_timer ( tInvFourier )   

        do i=1,nGrid 
            call sftInv2d ( allGrids(i), rhs )
        enddo

#ifdef par
        call blacs_barrier ( ICTXT, 'All' ) 
#endif
        if(iAm==0.and.NRHS<=1) &
        write(*,*) '    Time to inverse Fourier transform the fields: ', end_timer ( tInvFourier ),  'seconds'


    !   Calculate the plasma current
    !   ----------------------------
    
#ifdef usepapi
        call PAPIF_flops ( papi_rTime, papi_pTime, papi_flpins, papi_mflops, papi_irc )
        papi_rTime_zero = papi_rTime
        papi_pTime_zero = papi_pTime
#endif
    
        call start_timer ( tCurrent )
    
        if (iAm==0.and.NRHS<=1) &
        write(*,*) 'Calculating plasma current'
    
        do i=1,nGrid
            call current ( allGrids(i), rhs )
        enddo
    
#ifdef par
        call blacs_barrier ( ICTXT, 'All' ) 
#endif

        Perf%TimeCurrent = end_timer ( tCurrent )
        if (iAm==0.and.NRHS<=1) &
        write(*,*) 'Total to calculate plasma current: ', Perf%TimeCurrent
    
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
        call sgSum2D ( ICTXT, 'All', ' ', 1, 1, &
                    papi_mflops_global, 1, -1, -1 )
#endif

        Perf%GflopsCurrentLocal = papi_mflops / 1e3
        Perf%GflopsCurrentGlobal = papi_mflops_global / 1e3
    
        if (iAm==0.and.NRHS<=1) then
    
            write(*,*) '    PAPI data for plasma current'
            write(*,'(a30,f12.1)') 'real time: ', papi_rTime_current
            write(*,'(a30,f12.1)') 'proc time: ', papi_pTime_current
            write(*,'(a30,f12.2)') 'iAm Gflops/s: ', papi_mflops / 1e3
            write(*,'(a30,f12.2)') 'global Gflops/s: ', papi_mflops_global / 1e3
            write(*,'(a30,i1)') 'status: ', papi_irc
    
        endif
#endif
#ifdef par
        call blacs_barrier ( ICTXT, 'All' ) 
#endif
        if(iAm==0.and.NRHS<=1) &
        write(*,*) '    Time to sum flops 3: ', end_timer ( tSumFlops3 ),  'seconds'
    
    !   Rotation to Lab frame 
    !   ---------------------
    
        call start_timer ( tRotate )
        if (iAm==0.and.NRHS<=1) &
        write(*,*) 'Rotating E solution to lab frame (R,Th,z)'
    
        if ( iAm == 0 ) then
    
            do i=1,nGrid 
                call rotate_E_to_lab ( allGrids(i), rhs )
            enddo
    
        endif
#ifdef par
        call blacs_barrier ( ICTXT, 'All' ) 
#endif
        if(iAm==0.and.NRHS<=1) &
        write(*,*) '    Time to rotate solution to lab frame: ', end_timer ( tCurrent ),  'seconds'
    
    !   Calculate the Joule Heating 
    !   ---------------------------
    
        call start_timer ( tJDotE )
        if (iAm==0.and.NRHS<=1) &
        write(*,*) 'Calculating Joule heating' 
    
        if ( iAm == 0 ) then
    
            do i=1,nGrid
                call jDotE ( allGrids(i), rhs )
            enddo
    
        endif
#ifdef par
        call blacs_barrier ( ICTXT, 'All' ) 
#endif
        if(iAm==0.and.NRHS<=1) &
        write(*,*) '    Time to calculate Joule heating: ', end_timer ( tJDotE ),  'seconds'
    
    
    !   Write soln to file
    !   ------------------
    
        call start_timer ( tWriteSolution )
        if (iAm==0.and.NRHS<=1) &
        write(*,*) 'Writing solution to file'
    
        if ( iAm == 0 ) then
    
            do i=1,nGrid
                call write_solution ( allGrids(i), rid, rhs )
            enddo
    
        endif
#ifdef par
        call blacs_barrier ( ICTXT, 'All' ) 
#endif
        if(iAm==0.and.NRHS<=1) &
        write(*,*) '    Time to write solution: ', end_timer ( tWriteSolution ),  'seconds'
    
#ifdef usepapi
        call PAPIF_flops ( papi_rTime, papi_pTime, papi_flpins, papi_mflops, papi_irc )
        papi_rTime_Total = papi_rTime - papi_rTime_ProgramZero
        papi_pTime_Total = papi_pTime - papi_pTime_ProgramZero
#endif
#ifdef usepapi
    
        if (iAm==0.and.NRHS<=1) then
    
            write(*,*) '    PAPI data for program'
            write(*,'(a30,f12.1)') 'real time: ', papi_rTime_total
            write(*,'(a30,f12.1)') 'proc time: ', papi_pTime_total
            write(*,'(a30,i1)') 'status: ', papi_irc
    
        endif
#endif


    enddo RHS_PostProcessing_loop

    call start_timer(tReleaseGrid)
#ifdef par
    call release_grid ()
#endif

    if(iAm==0) &
    write(*,*) '    Time to release mpi grid: ', end_timer ( tReleaseGrid ),  'seconds'

    !stat = gptlStop ('total')
    !stat = gptlpr (0)
    !stat = gptlFinalize ()

    if (iAm==0) then

        Perf%nProcs = npRow*npCol
        Perf%nSpatialPoints = nPts_tot
        Perf%nRowLocal = LM_A
        Perf%nColLocal = LN_A
        Perf%nRowGlobal = M
        Perf%nColGlobal = N
        Perf%MatSizeLocal_GB = LocalSizeGB
        Perf%MatSizeGlobal_GB = GlobalSizeGB
        Perf%TimeTotal = end_timer ( tTotal )

        write(*,*) '    Total time: ', Perf%TimeTotal

        call WritePerformanceData ( Perf, rid )

    endif

end program aorsa2dMain
