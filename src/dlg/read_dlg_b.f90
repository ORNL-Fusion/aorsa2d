module dlg_bField

contains

subroutine read_dlg_bField ( capR, y, dlg_bR, dlg_bz, dlg_bPhi, &
                ncFileNameIn, b0, rmaxis, zmaxis, negPol, negTor )

    use netcdf
    use dlg
    implicit none
   
    !   capR, y are the R, z values from AORSA, i.e., capR(i) and y(j)

    integer :: nc_id, b_dim_ids(2), &
        R_nBins, z_nBins, &
        R_id, z_id, &
        bR_id, bz_id, ncStat, &
        bPhi_id
    integer ::  nR, nz, k, l
    character(len=*), intent(in), optional :: ncFileNameIn
    character(len=100) ::  ncFileName
    real, intent(IN) :: capR(:), y(:)
    logical, intent(in), optional :: negPol, negTor

    real, intent(in), optional :: rmaxis, zmaxis
    real, allocatable :: bR(:,:), bz(:,:), bPhi(:,:), &
        R(:), z(:)
    integer, allocatable :: R_index(:,:), z_index(:,:)
    integer :: insaneCnt, i, j
    real, intent(inout) :: dlg_bR(:,:), dlg_bz(:,:), dlg_bPhi(:,:)
    real, intent(out), optional :: b0
    real :: minR, minz, RStep, zStep, RRange, zRange

    allocate ( R_index ( size ( capR ), size ( y ) ), &
        z_index ( size ( capR ), size ( y ) ) )

    !   Read particle f netCdf file
   
    if ( present ( ncFileNameIn ) ) then 
        ncFileName = ncFileNameIn
    else 
        ncFileName  = 'dlg_bField.nc'
    endif
 
    ncStat = nf90_open ( path = ncFileName, mode = nf90_nowrite, ncid = nc_id )
    ncStat = nf90_inq_varId ( nc_id, 'bR', bR_id )
    ncStat = nf90_inq_varId ( nc_id, 'bz', bz_id )
    ncStat = nf90_inq_varId ( nc_id, 'bPhi', bPhi_id )

    ncStat = nf90_inq_varId ( nc_id, 'R', R_id )
    ncStat = nf90_inq_varId ( nc_id, 'z', z_id )
    
    ncStat = nf90_inquire_variable ( nc_id, bR_id, dimIds = b_dim_ids )
    
    ncStat = nf90_inquire_dimension ( nc_id, b_dim_ids(1), len = R_nBins ) 
    ncStat = nf90_inquire_dimension ( nc_id, b_dim_ids(2), len = z_nBins ) 

    allocate ( bR ( R_nBins, z_nBins ), &
        bz ( R_nBins, z_nBins ), &
        bPhi ( R_nBins, z_nBins ), &
        R ( R_nBins ), &
        z ( z_nBins ) )

    bR = 0.0
    bz = 0.0
    bPhi = 0.0
    dlg_bR  = 0.0
    dlg_bz   = 0.0
    dlg_bPhi    = 0.0

    ncStat = nf90_get_var ( nc_id, bR_id, bR ) 
    ncStat = nf90_get_var ( nc_id, bz_id, bz ) 
    ncStat = nf90_get_var ( nc_id, bPhi_id, bPhi ) 

    ncStat = nf90_get_var ( nc_id, R_id, R ) 
    ncStat = nf90_get_var ( nc_id, z_id, z ) 
    
    ncStat = nf90_close ( nc_id ) 

    !   Find R,z index of nearest p2f grid point

    RStep   = R(2) - R(1)
    zStep   = z(2) - z(1)

    RRange  = ( maxVal ( R ) - minVal ( R ) ) + RStep
    zRange  = ( maxVal ( z ) - minVal ( z ) ) + zStep

    minR    = minVal ( R ) - RStep / 2.0
    minz    = minVal ( z ) - zStep / 2.0

    do i=1,size(capR) 
        do j=1,size(y)

                R_index(i,j) = nint ( ( capR(i) - minR ) / RRange * R_nBins+1 )
                z_index(i,j) = nint ( ( y(j) - minz ) / zRange * z_nBins+1 )

                if ( R_index(i,j) > 0 .and. R_index(i,j) < R_nBins+1 &
                    .and. z_index(i,j) > 0 .and. z_index(i,j) < z_nBins+1 ) then

                    dlg_bR(i,j)  = bR(R_index(i,j),z_index(i,j))
                    dlg_bz(i,j)  = bz(R_index(i,j),z_index(i,j))
                    dlg_bPhi(i,j)  = bPhi(R_index(i,j),z_index(i,j))

                endif

            sanity_check_x: &
            if ( dlg_bR(i,j) /= dlg_bR(i,j) &
                .or. dlg_bR(i,j) * 0 /= 0 ) then

                write(*,*) 'ERROR: [x] sanity failure on reading antenna current from netcdf file, NaN, Inf detected'
                write(*,*) dlg_bR(i,j)
                stop
            endif sanity_check_x

            sanity_check_y: &
            if ( dlg_bz(i,j) /= dlg_bz(i,j) &
                .or. dlg_bz(i,j) * 0 /= 0 ) then

                write(*,*) 'ERROR: [y] sanity failure on reading antenna current from netcdf file, NaN, Inf detected'
                write(*,*) dlg_bz(i,j)
                stop
            endif sanity_check_y

            sanity_check_z: &
            if ( dlg_bPhi(i,j) /= dlg_bPhi(i,j) &
                .or. dlg_bPhi(i,j) * 0 /= 0 ) then

                write(*,*) 'ERROR: [z] sanity failure on reading antenna current from netcdf file, NaN, Inf detected'
                write(*,*) dlg_bPhi(i,j)
                stop
            endif sanity_check_z


        enddo
    enddo

    axisField: &
    if ( present ( b0 ) .and. present ( rmaxis ) .and. present ( zmaxis ) ) then

        R_index(1,1) = nint ( ( rmaxis - minR ) / RRange * R_nBins+1 )
        z_index(1,1) = nint ( ( zmaxis - minz ) / zRange * z_nBins+1 )

        b0  = bPhi(R_index(1,1),z_index(1,1))

    endif axisField

    flipSignPol: &
    if ( present ( negPol ) ) then

        dlg_bR  = -dlg_bR
        dlg_bz  = -dlg_bz

    endif flipSignPol

    flipSignTor: &
    if ( present ( negTor ) ) then

        b0  = -b0
        dlg_bPhi = -dlg_bPhi

    endif flipSignTor


    sanity_check_value: &
    if (count(abs(dlg_bR) > 0) == 0 &
        .and.count(abs(dlg_bR) > 0) == 0) then
 
        write(*,*) 'ERROR: sanity failure on reading antenna current from netcdf file, all zeros'
        write(*,*) maxVal(dlg_bR), minval(dlg_bR),maxVal(dlg_bz), minval(dlg_bz)
        write(*,*) maxVal(bR), minval(bR),maxVal(bz), minval(bz)
        write(*,*) shape(bR), shape(bz), shape(dlg_bR), shape(dlg_bz)
        write(*,*) size ( capR ), size ( y )
        write(*,*) maxVal(capR), minVal(capR), maxVal(y), minVal(y)
        stop
    endif sanity_check_value


    deallocate ( bR, bPhi, bz, &
        R, &
        z, &
        R_index, &
        z_index )

    return

end subroutine read_dlg_bField

end module dlg_bField
