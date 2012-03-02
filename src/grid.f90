module grid

use constants
use chebyshev_mod
use aorsa2din_mod, &
only: chebyshevX, chebyshevY, cosX, cosY

implicit none

! Define the workList type
! ------------------------

type :: workListEntry
        integer :: i
        integer :: j
        integer :: m
        integer :: n
end type workListEntry


! Define the grid objects
! -----------------------

type :: gridBlock

    ! Grid number
    ! -----------

    integer :: gridNumber
    character(len=3) :: fNumber

    ! Matrix offset
    ! -------------

    integer :: startRow, startCol

    ! Grid itself
    ! -----------
    integer :: nR, nZ, nModesR, nModesZ
    real, allocatable, dimension(:) :: rNorm, zNorm, R, z, kPhi
    real :: rMin, rMax, zMin, zMax, rRange, zRange
    real :: rMinIn, rMaxIn, zMinIn, zMaxIn
    real :: normFacR, normFacZ
    integer, allocatable :: label(:,:)
    integer, allocatable :: neighbr_startRow(:,:),neighbr_startCol(:,:)

    ! Basis functions
    ! ---------------

    integer :: nMin, nMax, mMin, mMax
    real :: k_cutOff
    complex, allocatable, dimension(:,:) :: xx, yy
    complex, allocatable, dimension(:,:) :: &
        drBfn_bfn, dzBfn_bfn, &
        d2rBfn_bfn, d2zBfn_bfn, &
        d3rBfn_bfn, d3zBfn_bfn, &
        d4rBfn_bfn, d4zBfn_bfn

    ! B field, temp, density, frequencies
    ! -----------------------------------
    real, allocatable, dimension(:,:) :: bR_unit, bT_unit, bZ_unit, bMag, rho
    logical, allocatable, dimension(:,:) :: mask
    real, allocatable :: nuOmg(:,:)
    real(kind=dbl), allocatable, dimension(:,:,:) :: densitySpec, ktSpec
    real(kind=dbl), allocatable, dimension(:,:,:) :: omgc, omgp2

    ! Toroidal broadening variables
    ! -----------------------------
    real, allocatable, dimension(:,:,:,:) :: U_RTZ_to_ABb
    real, allocatable :: sinTh(:,:)
    real, allocatable, dimension(:,:) :: gradPrlB
    real, allocatable, dimension(:,:) :: bPol

    ! Antenna Currents
    ! ----------------
    complex, allocatable, dimension(:,:) :: &
        jR, jT, jZ

    ! E field solution
    ! ----------------
    complex, allocatable, dimension(:,:) :: &
        ealphak, ebetak, eBk
    complex, allocatable, dimension(:,:) :: &
       ealpha,  ebeta, eB 
    complex, allocatable, dimension(:,:) :: &
       eR,  eTh, eZ 

    ! Mode integration range
    ! ----------------------
    integer :: nS, nF, mS, mF

    ! Plasma currents per species
    ! ---------------------------
    complex, allocatable, dimension(:,:,:) :: &
        jAlpha, jBeta, jB, &
        jP_r, jP_t, jP_z

    ! Power absorption
    ! ----------------
    real, allocatable, dimension(:,:,:) :: &
        jouleHeating

    ! Rotation matrix and derivatives
    ! -------------------------------
    real, allocatable, dimension(:,:) :: &
        Urr, Urt, Urz, Utr, Utt, Utz, Uzr, Uzt, Uzz

    ! dr first derivatives
    real, allocatable, dimension(:,:) :: &
        drUrr, drUrt, drUrz, &
        drUtr, drUtt, drUtz, &
        drUzr, drUzt, drUzz
    
    ! dz first derivatives
    real, allocatable, dimension(:,:) :: &
        dzUrr, dzUrt, dzUrz, &
        dzUtr, dzUtt, dzUtz, &
        dzUzr, dzUzt, dzUzz
    
    ! drr second derivatives
    real, allocatable, dimension(:,:) :: &
        drrUrr, drrUrt, drrUrz, &
        drrUtr, drrUtt, drrUtz, &
        drrUzr, drrUzt, drrUzz
    
    ! dzz second derivatives
    real, allocatable, dimension(:,:) :: &
        dzzUrr, dzzUrt, dzzUrz, &
        dzzUtr, dzzUtt, dzzUtz, &
        dzzUzr, dzzUzt, dzzUzz
    
    ! drz derivatives
    real, allocatable, dimension(:,:) :: &
        drzUrr, drzUrt, drzUrz, &
        drzUtr, drzUtt, drzUtz, &
        drzUzr, drzUzt, drzUzz

    ! workList
    type(workListEntry), allocatable :: wl(:)

end type gridBlock


! Define bFn derivative argument
! ------------------------------

type :: dBfnArg

    integer :: n, m
    real :: xNorm, yNorm
    real :: normFacX, normFacY    

endtype dBfnArg


type(gridBlock), allocatable :: allGrids(:)

integer(kind=long), allocatable :: bndryBlockID(:,:), bndryType(:)


contains

    function init_gridBlock ( nR, nZ, rMin, rMax, zMin, zMax ) result ( grid )

        use aorsa2din_mod, &
        only : nPhi, xkPerp_cutOff, overlap, &
        rMinAll, rMaxAll, zMinAll, zMaxAll, nGrid, &
        nZ_1D
        use parallel

        implicit none

        type(gridBlock) :: grid
        type(dBfnArg) :: d

        integer, intent(in) :: nR, nZ
        real, intent(in) :: rMin, rMax, zMin, zMax

        integer :: i, j, n, m

        allocate ( grid%rNorm(nR), &
                    grid%zNorm(nZ), &
                    grid%R(nR), &
                    grid%z(nZ), &
                    grid%kPhi(nR) )

            ! Populate grid block parameters
            ! ------------------------------

            grid%nR = nR
            grid%nZ = nZ
            grid%nModesR = nR
            grid%nModesZ = nZ

            grid%rMinIn = rMin
            grid%rMaxIn = rMax
            grid%zMinIn = zMin
            grid%zMaxIn = zMax

            
            ! Create the normalised block grids
            ! for the basis functions.
            ! ---------------------------------

            if(nR>1) then 

                do i = 1, nR
     
                    if(chebyshevX)then
                        ! Gauss-Lobatto grid [-1,1] (Chebyshev basis)
                        grid%rNorm(i)  = -cos(pi*(i-1)/(nR-1))
                    elseif(cosX)then
                        ! Uniform grid [0,pi] (Cos basis)
                        grid%rNorm(i)  = (i-1) * pi / ( nR - 1 )
                    else
                        ! Uniform grid [0,2pi] (Exp basis)
                        grid%rNorm(i)  = (i-1) * 2 * pi / ( nR - 1 )
                    endif
       
                enddo
            else

                grid%rNorm(1) = (rMax-rMin)/2+rMin
                
            endif


            ! Create the real space block grids
            ! taking into account the grid 
            ! block overlap.
            ! ---------------------------------

            grid%rRange  = (rMax - rMin) * &
                (grid%rNorm(nR)-grid%rNorm(1)) / (grid%rNorm(nR-overlap)-grid%rNorm(1+overlap))   

            grid%rMin = rMin - abs(grid%rNorm(1+overlap)-grid%rNorm(1)) / (grid%rNorm(nR)-grid%rNorm(1)) * grid%rRange
            grid%rMax = rMax + abs(grid%rNorm(1+overlap)-grid%rNorm(1)) / (grid%rNorm(nR)-grid%rNorm(1)) * grid%rRange
            grid%zMin = zMin
            grid%zMax = zMax

            ! Trap for non-overlapping end-points
            ! However, this will not work for non-square domains
            ! in 2D so will need fixing.
            ! -----------------------------------

            if(rMin==minVal(rMinAll(1:nGrid))) grid%rMin = rMin
            if(rMax==maxVal(rMaxAll(1:nGrid))) grid%rMax = rMax
            grid%rRange  = (grid%rMax - grid%rMin) 


            if(nR>1)then
                grid%R = (grid%rNorm-grid%rNorm(1)) &
                    / (grid%rNorm(nR)-grid%rNorm(1)) * grid%rRange + grid%rMin
            else
                grid%R(1) = grid%rMin
            endif

            if(chebyshevX)then
                grid%normFacR = 2 / grid%rRange
            elseif(cosX)then
                !Cos
                grid%normFacR = pi / grid%rRange
            else
                !Exp
                grid%normFacR = 2 * pi / grid%rRange
            endif

      
            grid%zRange  = zMax - zMin
            if(nZ>1) then 

                do j = 1, nZ

                    if(chebyshevY) then
                        ! Gauss-Lobatto grid [-1,1] (Chebyshev basis)
                        grid%zNorm(j)  = -cos(pi*(j-1)/(nZ-1))
                    elseif(cosY)then
                        ! Uniform grid [0,pi] (Cos basis)
                        grid%zNorm(j)  = (j-1) * pi / ( nZ - 1 )
                    else
                        ! Uniform grid [0,2pi] (Exp basis)
                        grid%zNorm(j)  = (j-1) * 2 * pi / ( nZ - 1 )
                    endif

                enddo
            else
                grid%z(1)    = (zMax-zMin)/2.0+zMin
                grid%zNorm(1) = 0
            endif

            if(nZ>1) &
            grid%z = (grid%zNorm-grid%zNorm(1)) &
                / (grid%zNorm(nZ)-grid%zNorm(1)) * grid%zRange + zMin

            if(chebyshevY)then
                grid%normFacZ = 2 / grid%zRange
            elseif(cosY)then
                !Cos
                grid%normFacZ = pi / grid%zRange
            else 
                !Exp
                grid%normFacZ = 2 * pi / grid%zRange
            endif


            ! Set the grid block kPhi
            ! -----------------------

            grid%kPhi = nPhi / grid%R


            ! Set the grid block n,m ranges
            ! -----------------------------

            if(chebyshevX)then

                grid%nMin    = 0
                grid%nMax    = nR-1

            elseif(cosX)then

                grid%nMin    = 0
                grid%nMax    = nR-1

            else

                grid%nMin = -nR/2
                grid%nMax =  nR/2

                ! Catch for even number of modes
                ! for producing a square matrix
                ! ------------------------------

                if (mod(nR,2)==0) grid%nMin = grid%nMin+1

            endif


            if(chebyshevY) then

                grid%mMin    = 0
                grid%mMax    = nZ-1

            elseif(cosY)then
                
                grid%mMin    = 0
                grid%mMax    = nZ-1

            else

                grid%mMin = -nZ/2
                grid%mMax =  nZ/2
                if (mod(nZ,2)==0) grid%mMin = grid%mMin+1

            endif

            if(nZ==1)then 
                grid%mMin = nZ_1d
                grid%mMax = nZ_1d
            endif

            if (iAm==0) then
                write(*,*) '    n: ', grid%nMin, grid%nMax
                write(*,*) '    m: ', grid%mMin, grid%mMax
            endif


            ! Set the k above which will be damped
            ! ------------------------------------

            grid%k_cutOff = xkPerp_cutOff * sqrt(&
               (grid%nMax * grid%normFacR)**2+(grid%mMax*grid%normFacZ)**2)


            ! Fill the grid block basis functions
            ! -----------------------------------

            allocate ( &
                grid%xx(grid%nMin:grid%nMax,nR), &
                grid%yy(grid%mMin:grid%mMax,nZ), &
                grid%drBfn_bfn(grid%nMin:grid%nMax,nR), &
                grid%dzBfn_bfn(grid%mMin:grid%mMax,nZ) )


            do i = 1, nR
                do n = grid%nMin, grid%nMax 

                    grid%xx(n,i) = xBasis(n,grid%rNorm(i))

                    d%n = n
                    d%m  = 0
                    d%xNorm = grid%rNorm(i)
                    d%yNorm = 0
                    d%normFacX = grid%normFacR
                    d%normFacY = grid%normFacZ

                    grid%dRbfn_bfn(n,i) = drBfn_bfn(d)

                enddo
            enddo

            do j = 1, nZ
                do m = grid%mMin, grid%mMax 

                    grid%yy(m,j) = yBasis(m,grid%zNorm(j))

                    d%n = 0
                    d%m  = m
                    d%xNorm = 0
                    d%yNorm = grid%zNorm(j)
                    d%normFacX = grid%normFacR
                    d%normFacY = grid%normFacZ

                    grid%dZBfn_bfn(m,j) = dzBfn_bfn(d)

                enddo
            enddo


    end function init_gridBlock


    ! Label each grid point within a block accoring to
    ! interior, outer boundary or inner boundary
    ! ------------------------------------------------

    subroutine labelPts ( gAll, nPts_tot )

        use aorsa2din_mod, &
        only: rMinAll, rMaxAll, zMinAll, zMaxAll, nGrid, overlap

        implicit none

        type(gridBlock), intent(inout) :: gAll(:)
        integer, intent(in) :: nPts_tot

        integer :: iMe, jMe, me, nbr, offSet, label


        ! Boundary conditions
        ! -------------------
        ! 0   interior, Maxwell
        ! 888 exterior, E = 0
        ! 999(-999) E - E = 0
        ! 1(-1)   E - E = 0 @ first overlap pt (some interpolation method)
        ! 2(-2)   E - E = 0     "         " 
        ! 3(-3)   E - E = 0     "         "

        allocate ( bndryBlockID(4,nPts_tot), &
                    bndryType(nPts_tot) )

        bndryBlockID = 0
        bndryType = 0
        offSet = 0

        do me=1,nGrid

            allocate ( gAll(me)%label(gAll(me)%nR,gAll(me)%nZ) )

            ! Interior points
            ! ---------------

            gAll(me)%label = 0

            do iMe=1,gAll(me)%nR
                do jMe=1,gAll(me)%nZ


                    ! Left side of block
                    ! ------------------

                    i_am_an_overlapped_pt_left: &
                    if(iMe<=(1+overlap) .and. gAll(me)%nR>1) then

                        at_mesh_mesh_bndry_rhs: &
                        if(count(gAll(me)%rMinIn==rMaxAll(1:nGrid))>0) then

                            label = abs((iMe-1)-overlap)*2-1 ! odd R derivative terms (1,3,5...)

                            ! special case for E=E 
                            ! (only left boundary, right one will be left zero for Maxwell)

                            if(abs((iMe-1)-overlap)*2==0) label = 999 

                            gAll(me)%label(iMe,jMe) = label

                            ! find the neighbour block
                            do nbr=1,nGrid                         

                                if(gAll(me)%rMinIn==rMaxAll(nbr))then                  

                                    bndryBlockID(:,offSet+(iMe-1)*gAll(me)%nZ+jMe) = (/iMe,jMe,nbr,me/)
                                    bndryType(offSet+(iMe-1)*gAll(me)%nZ+jMe) = -label 

                                endif

                            enddo

                        else ! at domain boundary

                            if(iMe==1) gAll(me)%label(iMe,jMe) = 888 ! outer boundary

                        endif at_mesh_mesh_bndry_rhs

                    endif i_am_an_overlapped_pt_left


                    ! Right side of block
                    ! -------------------

                    i_am_an_overlapped_pt_right: &
                    if(iMe>=gAll(me)%nR-overlap .and. gAll(me)%nR>1) then

                        at_mesh_mesh_bndry_lhs: &
                        if(count(gAll(me)%rMaxIn==rMinAll(1:nGrid))>0) then

                            label = abs((gAll(me)%nR-iMe)-overlap)*2 ! even R derivative terms (2,4,6...)

                            ! if label equals 0 then Maxwells will be used (see
                            ! left block for E=E condition)

                            gAll(me)%label(iMe,jMe) = label ! left part of mesh-mesh boundary bFn 

                            ! find the neighbour block
                            do nbr=1,nGrid                         

                                if(gAll(me)%rMaxIn==rMinAll(nbr))then                  

                                    bndryBlockID(:,offSet+(iMe-1)*gAll(me)%nZ+jMe) = (/iMe,jMe,nbr,me/)
                                    bndryType(offSet+(iMe-1)*gAll(me)%nZ+jMe) = -label ! (left)

                                endif

                            enddo

                        else ! at domain boundary

                            if(iMe==gAll(me)%nR) gAll(me)%label(iMe,jMe) = 888 ! outer boundary

                        endif at_mesh_mesh_bndry_lhs

                    endif i_am_an_overlapped_pt_right


                    ! Bottom block boundary
                    ! ---------------------

                    if(jMe==1 .and. gAll(me)%nZ>1) then
                         if(count(gAll(me)%zMinIn==zMaxAll)>0) then
                            gAll(me)%label(iMe,jMe) = 4 ! inner boundary bot
                        else
                            gAll(me)%label(iMe,jMe) = 888 ! outer boundary
                        endif
                    endif


                    ! Top block boundary
                    ! ------------------

                    if(jMe==gAll(me)%nZ .and. gAll(me)%nZ>1) then
                         if(count(gAll(me)%zMaxIn==zMinAll)>0) then
                            gAll(me)%label(iMe,jMe) = 5 ! inner boundary top
                        else
                            gAll(me)%label(iMe,jMe) = 888 ! outer boundary
                        endif
                    endif


                enddo
            enddo

            offSet = offSet + gAll(me)%nR * gAll(me)%nZ

        enddo

    end subroutine labelPts



    ! Basis function and their derivative routines
    ! for Exp, Chebyshev, Cos
    ! --------------------------------------------

    function xBasis (n,xNorm)

        implicit none
        integer, intent(in) :: n
        real, intent(in) :: xNorm
        complex :: xBasis

        if(chebyshevX) then
            xBasis = chebT ( n, xNorm )
            !xBasis = chebU ( n, xNorm )
        elseif(cosX) then
            xBasis = cos ( n * xNorm )
        else
            xBasis = exp ( zi * n * xNorm )
        endif

    end function xBasis


    function yBasis (m,yNorm)

        implicit none
        integer, intent(in) :: m
        real, intent(in) :: yNorm
        complex :: yBasis

        if (chebyshevY) then
            yBasis = chebT ( m, yNorm )
            !yBasis = chebU ( m, yNorm )
        elseif (cosY) then
            yBasis = cos ( m * yNorm )
        else
            yBasis = exp ( zi * m * yNorm )
        endif

    end function yBasis


    function drBfn_bfn( d )

        ! REMEMBER: this gives drBfn/Bfn 

        implicit none

        type(dBfnArg), intent(in) :: d

        complex :: drBfn_bfn

        if(chebyshevX)then
            ! chebT
            drBfn_bfn = d%n * chebU(-1+d%n,d%xNorm) &
                / chebT(d%n,d%xNorm)
            !! chebU
            !drBfn_bfn = ( (-1-d%n) * chebU(-1+d%n,d%xNorm) &
            !    + d%n * d%xNorm * chebU(d%n,d%xNorm) ) &
            !    / ( (-1+d%xNorm**2)*chebU(d%n,d%xNorm) )
        elseif(cosX)then
            drBfn_bfn = -d%n * tan ( d%n * d%xNorm )
        else
            drBfn_bfn = zi * d%n
        endif

        drBfn_bfn = drBfn_bfn * d%normFacX

    end function drBfn_bfn


    function dzBfn_bfn( d )

        ! REMEMBER: this gives dzBfn/Bfn 

        implicit none

        type(dBfnArg), intent(in) :: d

        complex :: dzBfn_bfn
        
        if(chebyshevY)then
            ! chebT
            dzBfn_bfn = d%m * chebU(-1+d%m,d%yNorm) &
                / chebT(d%m,d%yNorm)
            !! chebU
            !dzBfn_bfn = ( (-1-d%m) * chebU(-1+d%m,d%yNorm) &
            !    + d%m * d%yNorm * chebU(d%m,d%yNorm) ) &
            !    / ( (-1+d%yNorm**2)*chebU(d%m,d%yNorm) )
        elseif(cosY)then
            dzBfn_bfn = -d%m * tan ( d%m * d%yNorm )
        else
            dzBfn_bfn = zi * d%m
        endif

        dzBfn_bfn = dzBfn_bfn * d%normFacY

    end function dzBfn_bfn


    function drrBfn_bfn( d )

        ! REMEMBER: this gives drrBfn/Bfn 

        implicit none

        type(dBfnArg), intent(in) :: d

        complex :: drrBfn_bfn

        if(chebyshevX)then
            ! chebT
            drrBfn_bfn = d%n * &
                ( -d%n * chebU(-2+d%n,d%xNorm) &
                    + (-1+d%n)*d%xNorm * chebU(-1+d%n,d%xNorm )) &
                / ((-1+d%xNorm**2)*chebT(d%n,d%xNorm))
            !! chebU
            !drrBfn_bfn = (d%n*(1+d%n)*chebU(-2+d%n,d%xNorm) &
            !    +(3+d%n-2*d%n**2)*d%xNorm*chebU(-1+d%n,d%xNorm) &
            !    +d%n*(-1+(-1+d%n)*d%xNorm**2)*chebU(d%n,d%xNorm)) &
            !    / ( (-1+d%xNorm**2)**2*chebU(d%n,d%xNorm) )
        elseif(cosX)then
            drrBfn_bfn = -d%n**2
        else
            drrBfn_bfn = -d%n**2
        endif

        drrBfn_bfn = drrBfn_bfn * d%normFacX**2

    end function drrBfn_bfn


    function dzzBfn_bfn( d )

        ! REMEMBER: this gives dzzBfn/Bfn 

        implicit none

        type(dBfnArg), intent(in) :: d

        complex :: dzzBfn_bfn

        if(chebyshevY) then 
            ! chebT
            dzzBfn_bfn = d%m * &
                ( -d%m*chebU(-2+d%m,d%yNorm) &
                    + (-1+d%m)*d%yNorm*chebU(-1+d%m,d%yNorm)) &
                / ( (-1+d%yNorm**2)*chebT(d%m,d%yNorm))
            !! chebU
            !dzzBfn_bfn = (d%m*(1+d%m)*chebU(-2+d%m,d%yNorm) &
            !    +(3+d%m-2*d%m**2)*d%yNorm*chebU(-1+d%m,d%yNorm) &
            !    +d%m*(-1+(-1+d%m)*d%yNorm**2)*chebU(d%m,d%yNorm)) &
            !    / ( (-1+d%yNorm**2)**2*chebU(d%m,d%yNorm) )
        elseif(cosY)then
            dzzBfn_bfn = -d%m**2
        else
            dzzBfn_bfn = -d%m**2
        endif

        dzzBfn_bfn = dzzBfn_bfn * d%normFacY**2

    end function dzzBfn_bfn


    function drzBfn_bfn( d )

        ! REMEMBER: this gives drzBfn/Bfn 

        implicit none

        type(dBfnArg), intent(in) :: d
        complex :: drzBfn_bfn, drBfn_bfn, dzBfn_bfn

        if(chebyshevX)then
            ! chebT
            drBfn_bfn = d%n * chebU(-1+d%n,d%xNorm) &
                / chebT(d%n,d%xNorm)
            !! chebU
            !drBfn_bfn = ( (-1-d%m) * chebU(-1+d%m,d%yNorm) &
            !    + d%m * d%yNorm * chebU(d%m,d%yNorm) ) &
            !    / ( (-1+d%yNorm**2)*chebU(d%m,d%yNorm) )
        elseif(cosX)then
            drBfn_bfn = -d%n * tan( d%n*d%xNorm )
        else
            drBfn_bfn = zi * d%n
        endif

        if(chebyshevY)then
            ! chebT
            dzBfn_bfn = d%m * chebU(-1+d%m,d%yNorm) &
                / chebT(d%m,d%yNorm)
            !! chebU
            !drBfn_bfn = ( (-1-d%m) * chebU(-1+d%m,d%yNorm) &
            !    + d%m * d%yNorm * chebU(d%m,d%yNorm) ) &
            !    / ( (-1+d%yNorm**2)*chebU(d%m,d%yNorm) )
        elseif(cosY)then
            dzBfn_bfn = -d%m * tan( d%m*d%yNorm )
        else
            dzBfn_bfn = zi * d%m
        endif

        drzBfn_bfn = drBfn_bfn * dzBfn_bfn * ( d%normFacX * d%normFacY )

    end function drzBfn_bfn


end module grid
