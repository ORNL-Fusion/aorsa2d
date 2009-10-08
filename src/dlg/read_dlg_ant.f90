module dlg_ant

contains

subroutine read_dlg_ant ( capR, y, dlg_jantx, dlg_janty, ncFileNameIn, gridMatch )

    use netcdf
    use dlg
    implicit none
   
    !   capR, y are the R, z values from AORSA, i.e., capR(i) and y(j)

    integer :: nc_id, j_dim_ids(2), &
        R_nBins, z_nBins, &
        R_binCenters_id, z_binCenters_id, &
        R_binEdges_id, z_binEdges_id, &
        jantx_id, janty_id, ncStat
    integer ::  nR, nz, k, l
    character(len=*), intent(in), optional :: ncFileNameIn
    logical, intent(inout), optional :: gridMatch
    character(len=100) ::  ncFileName
    real, intent(IN) :: capR(:), y(:)

    real, allocatable :: jantx(:,:), janty(:,:), &
        R_bincenters(:), z_binCenters(:), &
        R_binEdges(:), z_binEdges(:)
    integer, allocatable :: R_index(:,:), z_index(:,:)
    integer :: insaneCnt, i, j
    real, intent(inout) :: dlg_jantx(:,:), dlg_janty(:,:)
    real :: minR, minz, RStep, zStep, RRange, zRange

    allocate ( R_index ( size ( capR ), size ( y ) ), &
        z_index ( size ( capR ), size ( y ) ) )

    !   Read particle f netCdf file
   
    if ( present ( ncFileNameIn ) ) then 
        ncFileName = ncFileNameIn
    else 
        ncFileName  = 'dlgAnt.nc'
    endif
 
    if ( .not. present ( gridMatch ) ) gridMatch  = .false. 
    
    ncStat = nf90_open ( path = ncFileName, mode = nf90_nowrite, ncid = nc_id )
    ncStat = nf90_inq_varId ( nc_id, 'jantx', jantx_id )
    ncStat = nf90_inq_varId ( nc_id, 'janty', janty_id )

    ncStat = nf90_inq_varId ( nc_id, 'R_binCenters', R_binCenters_id )
    ncStat = nf90_inq_varId ( nc_id, 'z_binCenters', z_binCenters_id )
    
    ncStat = nf90_inquire_variable ( nc_id, jantx_id, dimIds = j_dim_ids )
    
    ncStat = nf90_inquire_dimension ( nc_id, j_dim_ids(1), len = R_nBins ) 
    ncStat = nf90_inquire_dimension ( nc_id, j_dim_ids(2), len = z_nBins ) 

    allocate ( jantx ( R_nBins, z_nBins ), &
        janty ( R_nBins, z_nBins ), &
        R_binCenters ( R_nBins ), &
        z_binCenters ( z_nBins ) )

    jantx = 0.0
    janty = 0.0
    dlg_jantx   = 0.0
    dlg_janty   = 0.0

    ncStat = nf90_get_var ( nc_id, jantx_id, jantx ) 
    ncStat = nf90_get_var ( nc_id, janty_id, janty ) 

    ncStat = nf90_get_var ( nc_id, R_binCenters_id, R_binCenters ) 
    ncStat = nf90_get_var ( nc_id, z_binCenters_id, z_binCenters ) 
    
    ncStat = nf90_close ( nc_id ) 

    !   Find R,z index of nearest p2f grid point

    RStep   = R_binCenters(2) - R_binCenters(1)
    zStep   = z_binCenters(2) - z_binCenters(1)

    RRange  = ( maxVal ( R_binCenters ) - minVal ( R_binCenters ) ) + RStep
    zRange  = ( maxVal ( z_binCenters ) - minVal ( z_binCenters ) ) + zStep

    minR    = minVal ( R_binCenters ) - RStep / 2.0
    minz    = minVal ( z_binCenters ) - zStep / 2.0

    check_grids: &
    if (gridMatch) then
        
        if ( size ( capR ) /= size ( jantx, 1 ) .or. size ( y ) /= size ( jantx, 2 ) ) then

            write(*,*) 'ERROR: sanity failure on reading antenna current from netcdf file, gridMatch = .true. but check fails'
            write(*,*) size ( capR ), size ( jantx, 1 )
            write(*,*) size ( y ), size ( jantx, 2 )

            stop
     
        endif

    endif check_grids

    do i=1,size(capR) 
        do j=1,size(y)

            exact_aorsa_grid: &
            if (gridMatch) then

                dlg_jantx(i,j)  = jantx(i,j)
                dlg_janty(i,j)  = janty(i,j)

            else

                R_index(i,j) = nint ( ( capR(i) - minR ) / RRange * R_nBins+1 )
                z_index(i,j) = nint ( ( y(j) - minz ) / zRange * z_nBins+1 )

                if ( R_index(i,j) > 0 .and. R_index(i,j) < R_nBins+1 &
                    .and. z_index(i,j) > 0 .and. z_index(i,j) < z_nBins+1 ) then

                    dlg_jantx(i,j)  = jantx(R_index(i,j),z_index(i,j))
                    dlg_janty(i,j)  = janty(R_index(i,j),z_index(i,j))

                endif
            endif exact_aorsa_grid


            sanity_check_x: &
            if ( dlg_jantx(i,j) /= dlg_jantx(i,j) &
                .or. dlg_jantx(i,j) * 0 /= 0 ) then

                write(*,*) 'ERROR: [x] sanity failure on reading antenna current from netcdf file, NaN, Inf detected'
                write(*,*) dlg_jantx(i,j)
                stop
            endif sanity_check_x

            sanity_check_y: &
            if ( dlg_janty(i,j) /= dlg_janty(i,j) &
                .or. dlg_janty(i,j) * 0 /= 0 ) then

                write(*,*) 'ERROR: [y] sanity failure on reading antenna current from netcdf file, NaN, Inf detected'
                write(*,*) dlg_janty(i,j)
                stop
            endif sanity_check_y


        enddo
    enddo

    sanity_check_value: &
    if (count(abs(dlg_janty) > 0) == 0 &
        .and.count(abs(dlg_janty) > 0) == 0) then
 
        write(*,*) 'ERROR: sanity failure on reading antenna current from netcdf file, all zeros'
        write(*,*) maxVal(dlg_janty), minval(dlg_janty),maxVal(dlg_jantx), minval(dlg_jantx), gridMatch
        write(*,*) maxVal(janty), minval(janty),maxVal(jantx), minval(jantx)
        write(*,*) shape(jantx), shape(janty), shape(dlg_jantx), shape(dlg_janty)
        write(*,*) size ( capR ), size ( y )
        write(*,*) maxVal(capR), minVal(capR), maxVal(y), minVal(y)
        stop
    endif sanity_check_value


    deallocate ( jantx, janty, &
        R_binCenters, &
        z_binCenters, &
        R_index, &
        z_index )

    return

end subroutine read_dlg_ant 

end module dlg_ant 
