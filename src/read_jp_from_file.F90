module read_jp_from_file 

use constants, only: dbl

implicit none

integer :: nR,nZ,nS
complex, allocatable :: Jp_r(:,:),Jp_t(:,:),Jp_z(:,:)
real, allocatable :: r(:),z(:)

contains

subroutine ReadJpFromFile (FileName)

    use netcdf
    use check_mod
    use aorsaNamelist, only: useJpFromFile

    implicit none

    real, allocatable :: Jp_r_re(:,:,:),Jp_t_re(:,:,:),Jp_z_re(:,:,:)
    real, allocatable :: Jp_r_im(:,:,:),Jp_t_im(:,:,:),Jp_z_im(:,:,:)

    character(len=*), intent(IN) :: FileName

    integer :: stat
    integer :: nc_id,br_id,bt_id,bz_id,nR_id,nZ_id, &
        r_id,z_id, nS_id, &
        Jp_r_re_id, Jp_t_re_id, Jp_z_re_id, &
        Jp_r_im_id, Jp_t_im_id, Jp_z_im_id

    call check( nf90_open(path=FileName,mode=nf90_nowrite,ncid=nc_id) )

    call check( nf90_inq_dimId(nc_id,'nR',nR_id) )
    call check( nf90_inq_dimId(nc_id,'nY',nZ_id) )
    call check( nf90_inq_dimId(nc_id,'nSpec',nS_id) )
    call check( nf90_inquire_dimension(nc_id,nR_id,len=nR) )
    call check( nf90_inquire_dimension(nc_id,nZ_id,len=nZ) )
    call check( nf90_inquire_dimension(nc_id,nS_id,len=nS) )

    call check( nf90_inq_varId(nc_id,'r',r_id) )
    call check( nf90_inq_varId(nc_id,'z',z_id) )

    call check( nf90_inq_varId(nc_id,'jP_r_re',Jp_r_re_id) )
    call check( nf90_inq_varId(nc_id,'jP_t_re',Jp_t_re_id) )
    call check( nf90_inq_varId(nc_id,'jP_z_re',Jp_z_re_id) )

    call check( nf90_inq_varId(nc_id,'jP_r_im',Jp_r_im_id) )
    call check( nf90_inq_varId(nc_id,'jP_t_im',Jp_t_im_id) )
    call check( nf90_inq_varId(nc_id,'jP_z_im',Jp_z_im_id) )
 
    allocate(Jp_r_re(nR,nZ,nS),Jp_t_re(nR,nZ,nS),Jp_z_re(nR,nZ,nS))
    allocate(Jp_r_im(nR,nZ,nS),Jp_t_im(nR,nZ,nS),Jp_z_im(nR,nZ,nS))
    allocate(Jp_r(nR,nZ),Jp_t(nR,nZ),Jp_z(nR,nZ))
    allocate(r(nR),z(nZ))

    call check( nf90_get_var(nc_id,r_id,r) )
    call check( nf90_get_var(nc_id,z_id,z) )
 
    call check( nf90_get_var(nc_id,Jp_r_re_id,Jp_r_re) )
    call check( nf90_get_var(nc_id,Jp_t_re_id,Jp_t_re) )
    call check( nf90_get_var(nc_id,Jp_z_re_id,Jp_z_re) )

    call check( nf90_get_var(nc_id,Jp_r_im_id,Jp_r_im) )
    call check( nf90_get_var(nc_id,Jp_t_im_id,Jp_t_im) )
    call check( nf90_get_var(nc_id,Jp_z_im_id,Jp_z_im) )

    call check( nf90_close(nc_id) )

    Jp_r = sum(cmplx(Jp_r_re,Jp_r_im),3)
    Jp_t = sum(cmplx(Jp_t_re,Jp_t_im),3)
    Jp_z = sum(cmplx(Jp_z_re,Jp_z_im),3)

    deallocate(Jp_r_re,Jp_t_re,Jp_z_re)
    deallocate(Jp_r_im,Jp_t_im,Jp_z_im)
 
end subroutine ReadJpFromFile

end module read_jp_from_file
