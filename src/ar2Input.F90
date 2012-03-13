module ar2Input

implicit none

integer :: nR,nZ
real, allocatable :: rbbbs(:),zbbbs(:),rlim(:),zlim(:)
real, allocatable :: br(:,:),bt(:,:),bz(:,:),bMag(:,:)
real, allocatable :: r(:),z(:)

contains

subroutine ReadAr2Input (AR2FileName)

    use netcdf
    use check_mod

    implicit none

    character(len=*), intent(IN) :: AR2FileName

    integer :: nc_id,br_id,bt_id,bz_id,nR_id,nZ_id, &
        nbbbs_id,nlim_id,nbbbs,nlim,r_id,z_id

    call check( &
        nf90_open(path=AR2FileName,mode=nf90_nowrite,ncid=nc_id) )
    call check( &
        nf90_inq_varId(nc_id,'r',r_id) )
     call check( &
        nf90_inq_varId(nc_id,'z',z_id) )
    call check( &
        nf90_inq_varId(nc_id,'br',br_id) )
    call check( &
        nf90_inq_varId(nc_id,'bt',bt_id) )
    call check( &
        nf90_inq_varId(nc_id,'bz',bz_id) )
    call check( &
        nf90_inq_dimId(nc_id,'nR',nR_id) )
    call check( &
        nf90_inq_dimId(nc_id,'nZ',nZ_id) )
    call check( &
        nf90_inquire_dimension(nc_id,nR_id,len=nR) )
    call check( &
        nf90_inquire_dimension(nc_id,nZ_id,len=nZ) )
    call check( &
        nf90_inq_dimId(nc_id,'nbbbs',nbbbs_id) )
    call check( &
        nf90_inquire_dimension(nc_id,nbbbs_id,len=nbbbs) )
    call check( &
        nf90_inq_dimId(nc_id,'nlim',nlim_id) )
    call check( &
        nf90_inquire_dimension(nc_id,nlim_id,len=nlim) )
 
    allocate(br(nR,nZ),bt(nR,nZ),bz(nR,nZ))
    allocate(rbbbs(nbbbs),zbbbs(nbbbs),rlim(nlim),zlim(nlim))
    allocate(r(nR),z(nZ))

    call check( &
        nf90_get_var(nc_id,r_id,r) )
    call check( &
        nf90_get_var(nc_id,z_id,z) )
 
    call check( &
        nf90_get_var(nc_id,br_id,br) )
    call check( &
        nf90_get_var(nc_id,bt_id,bt) )
    call check( &
        nf90_get_var(nc_id,bz_id,bz) )

    call check( &
        nf90_close(nc_id) )

end subroutine ReadAr2Input

end module ar2Input
