module ar2Input

use constants, only: dbl

implicit none

integer :: nR,nZ,nS
real, allocatable :: rbbbs(:),zbbbs(:),rlim(:),zlim(:)
real, allocatable :: br(:,:),bt(:,:),bz(:,:),bMag(:,:)
real, allocatable :: r(:),z(:)
real :: rMin,rMax,zMin,zMax
integer, allocatable :: AtomicZ(:), amu(:)
real(kind=dbl), allocatable :: Density_m3(:,:,:), Temp_eV(:,:,:)

contains

subroutine ReadAr2Input (AR2FileName)

    use netcdf
    use check_mod

    implicit none

    character(len=*), intent(IN) :: AR2FileName

    integer :: nc_id,br_id,bt_id,bz_id,nR_id,nZ_id, &
        nbbbs_id,nlim_id,nbbbs,nlim,r_id,z_id, &
        rbbbs_id,zbbbs_id,rLim_id,zLim_id,nS_id, &
        rMin_id,rMax_id,zMin_id,zMax_id,AtomicZ_id, &
        amu_id,Density_id,Temp_id

    call check( nf90_open(path=AR2FileName,mode=nf90_nowrite,ncid=nc_id) )

    call check( nf90_inq_varId(nc_id,'rMin',rMin_id) )
    call check( nf90_inq_varId(nc_id,'rMax',rMax_id) )
    call check( nf90_inq_varId(nc_id,'zMin',zMin_id) )
    call check( nf90_inq_varId(nc_id,'zMax',zMax_id) )
    call check( nf90_inq_varId(nc_id,'AtomicZ',AtomicZ_id) )
    call check( nf90_inq_varId(nc_id,'amu',amu_id) )

    call check( nf90_inq_varId(nc_id,'r',r_id) )
    call check( nf90_inq_varId(nc_id,'z',z_id) )
    call check( nf90_inq_varId(nc_id,'br',br_id) )
    call check( nf90_inq_varId(nc_id,'bt',bt_id) )
    call check( nf90_inq_varId(nc_id,'bz',bz_id) )
    call check( nf90_inq_dimId(nc_id,'nSpec',nS_id) )
    call check( nf90_inquire_dimension(nc_id,nS_id,len=nS) )
    call check( nf90_inq_dimId(nc_id,'nR',nR_id) )
    call check( nf90_inq_dimId(nc_id,'nZ',nZ_id) )
    call check( nf90_inquire_dimension(nc_id,nR_id,len=nR) )
    call check( nf90_inquire_dimension(nc_id,nZ_id,len=nZ) )
    call check( nf90_inq_dimId(nc_id,'nbbbs',nbbbs_id) )
    call check( nf90_inquire_dimension(nc_id,nbbbs_id,len=nbbbs) )
    call check( nf90_inq_dimId(nc_id,'nlim',nlim_id) )
    call check( nf90_inquire_dimension(nc_id,nlim_id,len=nlim) )

    call check( nf90_inq_varId(nc_id,'Density_m3',Density_id) )
    call check( nf90_inq_varId(nc_id,'Temp_eV',Temp_id) )

    call check( nf90_inq_varId(nc_id,'rbbbs',rbbbs_id) )
    call check( nf90_inq_varId(nc_id,'zbbbs',zbbbs_id) )
    call check( nf90_inq_varId(nc_id,'rlim',rLim_id) )
    call check( nf90_inq_varId(nc_id,'zlim',zLim_id) )

    allocate(br(nR,nZ),bt(nR,nZ),bz(nR,nZ))
    allocate(rbbbs(nbbbs),zbbbs(nbbbs),rlim(nlim),zlim(nlim))
    allocate(r(nR),z(nZ))
    allocate(AtomicZ(nS),amu(nS))
    allocate(Density_m3(nR,nZ,nS),Temp_eV(nR,nZ,nS))

    call check( nf90_get_var(nc_id,AtomicZ_id,AtomicZ) )
    call check( nf90_get_var(nc_id,amu_id,amu) )

    call check( nf90_get_var(nc_id,rMin_id,rMin) )
    call check( nf90_get_var(nc_id,rMax_id,rMax) )
    call check( nf90_get_var(nc_id,zMin_id,zMin) )
    call check( nf90_get_var(nc_id,zMax_id,zMax) )

    call check( nf90_get_var(nc_id,r_id,r) )
    call check( nf90_get_var(nc_id,z_id,z) )
 
    call check( nf90_get_var(nc_id,br_id,br) )
    call check( nf90_get_var(nc_id,bt_id,bt) )
    call check( nf90_get_var(nc_id,bz_id,bz) )

    call check( nf90_get_var(nc_id,rbbbs_id,rbbbs) )
    call check( nf90_get_var(nc_id,zbbbs_id,zbbbs) )
    call check( nf90_get_var(nc_id,rLim_id,rLim) )
    call check( nf90_get_var(nc_id,zLim_id,zLim) )

    call check( nf90_get_var(nc_id,Density_id,Density_m3) )
    call check( nf90_get_var(nc_id,Temp_id,Temp_eV) )

    call check( nf90_close(nc_id) )

end subroutine ReadAr2Input

end module ar2Input
