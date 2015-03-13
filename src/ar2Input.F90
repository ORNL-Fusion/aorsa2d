module ar2Input

use constants, only: dbl

implicit none

integer :: nR,nZ,nS
real, allocatable :: br(:,:),bt(:,:),bz(:,:),bMag(:,:)
real, allocatable :: jant_r(:,:),jant_t(:,:),jant_z(:,:),jant(:,:)
real, allocatable :: r(:),z(:)
real :: rMin,rMax,zMin,zMax
real, allocatable :: AtomicZ(:), amu(:)
integer, allocatable :: BbbMask(:,:),LimMask(:,:)
real(kind=dbl), allocatable :: Density_m3(:,:,:), Temp_eV(:,:,:), nuOmg(:,:,:)

contains

subroutine ReadAr2Input (AR2FileName)

    use netcdf
    use check_mod
    use aorsaNamelist, only: noPoloidalField

    implicit none

    character(len=*), intent(IN) :: AR2FileName

    integer :: stat
    integer :: nc_id,br_id,bt_id,bz_id,nR_id,nZ_id, &
        nbbbs_id,nlim_id,nbbbs,nlim,r_id,z_id, &
        rbbbs_id,zbbbs_id,rLim_id,zLim_id,nS_id, &
        rMin_id,rMax_id,zMin_id,zMax_id,AtomicZ_id, &
        amu_id,Density_id,Temp_id,BbbMask_id,LimMask_id, &
        nuOmg_id, jant_id, jant_r_id, jant_t_id, jant_z_id

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

    call check( nf90_inq_varId(nc_id,'jAnt',jant_id) )
    call check( nf90_inq_varId(nc_id,'jAnt_r',jant_r_id) )
    call check( nf90_inq_varId(nc_id,'jAnt_t',jant_t_id) )
    call check( nf90_inq_varId(nc_id,'jAnt_z',jant_z_id) )
 
    call check( nf90_inq_dimId(nc_id,'nSpec',nS_id) )
    call check( nf90_inquire_dimension(nc_id,nS_id,len=nS) )
    call check( nf90_inq_dimId(nc_id,'nR',nR_id) )
    call check( nf90_inq_dimId(nc_id,'nZ',nZ_id) )
    call check( nf90_inquire_dimension(nc_id,nR_id,len=nR) )
    call check( nf90_inquire_dimension(nc_id,nZ_id,len=nZ) )
    call check( nf90_inq_varId(nc_id,'Density_m3',Density_id) )
    call check( nf90_inq_varId(nc_id,'Temp_eV',Temp_id) )
    call check( nf90_inq_varId(nc_id,'nuOmg',nuOmg_id) )

    allocate(br(nR,nZ),bt(nR,nZ),bz(nR,nZ))
    allocate(jant(nR,nZ),jant_r(nR,nZ),jant_t(nR,nZ),jant_z(nR,nZ))
    allocate(r(nR),z(nZ))
    allocate(AtomicZ(nS),amu(nS))
    allocate(Density_m3(nR,nZ,nS),Temp_eV(nR,nZ,nS),nuOmg(nR,nZ,nS))

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

    call check( nf90_get_var(nc_id,jant_r_id,jant_r) )
    call check( nf90_get_var(nc_id,jant_t_id,jant_t) )
    call check( nf90_get_var(nc_id,jant_z_id,jant_z) )
    call check( nf90_get_var(nc_id,jant_id,jant) )

    call check( nf90_get_var(nc_id,Density_id,Density_m3) )
    call check( nf90_get_var(nc_id,Temp_id,Temp_eV) )
    call check( nf90_get_var(nc_id,nuOmg_id,nuOmg) )

    stat = nf90_inq_varId(nc_id,'BbbMask',BbbMask_id)
    stat = nf90_inq_varId(nc_id,'LimMask',LimMask_id)
    allocate(BbbMask(nR,nZ),LimMask(nR,nZ))

    stat = nf90_get_var(nc_id,BbbMask_id,BbbMask)
    stat = nf90_get_var(nc_id,LimMask_id,LimMask)

    call check( nf90_close(nc_id) )

    if(noPoloidalField)then
        br = br*0
        bz = bz*0
    endif

end subroutine ReadAr2Input

end module ar2Input
