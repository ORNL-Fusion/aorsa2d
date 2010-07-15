module read_data

use netcdf
use check_mod

contains

subroutine read_sigma ( fName, sigma )
  
    implicit none

    complex, optional, intent(out) :: sigma(:,:,:,:,:,:)

    character(len=*), intent(in) :: fName

    integer :: nc_id, sigma_re_id, sigma_im_id, sigma_dim_ids(6)
    integer :: nR, nZ, nN, nM

    real, allocatable, dimension(:,:,:,:,:,:) :: sigma_re, sigma_im

    write(*,*) '    reading sigma file ', fName

    call check ( nf90_open ( path = fName, mode = nf90_nowrite, ncid = nc_id ) )
    call check ( nf90_inq_varId ( nc_id, 'sigma_re', sigma_re_id ) )
    call check ( nf90_inq_varId ( nc_id, 'sigma_im', sigma_im_id ) )
    
    call check ( nf90_inquire_variable ( &
        nc_id, sigma_re_id, dimIds = sigma_dim_ids ) )
    
    call check ( nf90_inquire_dimension ( nc_id, sigma_dim_ids(1), len = nR ) ) 
    call check ( nf90_inquire_dimension ( nc_id, sigma_dim_ids(2), len = nZ ) ) 
    call check ( nf90_inquire_dimension ( nc_id, sigma_dim_ids(3), len = nN ) ) 
    call check ( nf90_inquire_dimension ( nc_id, sigma_dim_ids(4), len = nM ) ) 

    allocate ( &
        sigma_re ( nR, nZ, nN, nM, 3, 3 ), &
        sigma_im ( nR, nZ, nN, nM, 3, 3 ) )

    call check ( nf90_get_var ( nc_id, sigma_re_id, sigma_re ) )
    call check ( nf90_get_var ( nc_id, sigma_im_id, sigma_im ) )
    
    call check ( nf90_close ( nc_id ) )

    if(present(sigma)) &
        sigma = cmplx ( sigma_re, sigma_im )
    
end subroutine read_sigma

end module read_data
