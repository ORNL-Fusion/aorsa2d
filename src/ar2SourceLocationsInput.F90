module ar2SourceLocationsInput

use constants, only: dbl

implicit none

integer :: nSources, NRHS
real, allocatable :: CurrentSource_r(:), CurrentSource_z(:), &
        CS_RealFac(:), CS_ImagFac(:)
integer, allocatable :: CurrentSource_ComponentID(:)

contains

subroutine ReadAR2SourceLocations (AR2CurrentSourceFileName)

    use netcdf
    use check_mod

    implicit none

    character(len=*), intent(IN) :: AR2CurrentSourceFileName

    integer :: stat,NRHS_check1,NRHS_check2
    integer :: nc_id,cs_r_id,cs_z_id,cid_id,nSources_id,nRHS_id,realFac_id,imagFac_id

    call check( nf90_open(path=AR2CurrentSourceFileName,mode=nf90_nowrite,ncid=nc_id) )

    call check( nf90_inq_varId(nc_id,'cs_r',cs_r_id) )
    call check( nf90_inq_varId(nc_id,'cs_z',cs_z_id) )
    call check( nf90_inq_varId(nc_id,'realFac',realFac_id) )
    call check( nf90_inq_varId(nc_id,'imagFac',imagFac_id) )

    call check( nf90_inq_varId(nc_id,'component_ident',cid_id) )

    call check( nf90_inq_dimId(nc_id,'nSources',nSources_id) )
    call check( nf90_inquire_dimension(nc_id,nSources_id,len=nSources) )

    call check( nf90_inq_dimId(nc_id,'NRHS',NRHS_id) )
    call check( nf90_inquire_dimension(nc_id,NRHS_id,len=NRHS) )
 
    if(NRHS.lt.1)then
            write(*,*) 'ERROR: NRHS < 1'
            stop
    endif

    allocate(CurrentSource_r(NRHS),CurrentSource_z(NRHS),&
        CurrentSource_ComponentID(NRHS),CS_RealFac(NRHS), CS_ImagFac(NRHS))

    call check( nf90_get_var(nc_id,cs_r_id,CurrentSource_r) )
    call check( nf90_get_var(nc_id,cs_z_id,CurrentSource_z) )
    call check( nf90_get_var(nc_id,cid_id,CurrentSource_ComponentID) )
    call check( nf90_get_var(nc_id,realFac_id,CS_RealFac) )
    call check( nf90_get_var(nc_id,imagFac_id,CS_ImagFac) )

    call check( nf90_close(nc_id) )

end subroutine ReadAR2SourceLocations

end module ar2SourceLocationsInput
