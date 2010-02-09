      module write_pf
      

      contains

      subroutine init_pf_dlg ( iValue, myId,
     . nNodeX, nNodeY, nuPer, nuPar, 
     . vc_mks, capR, zLoc, uPerp, uPara )

          use dlg_p_check
          use dlg
          use aorsa2din_mod, only: nModesX, nModesY
 
          implicit none

          include 'mpif.h'
#         include <pnetcdf.inc>
      
          ! Parameters
          integer, parameter :: DBL = selected_real_kind ( 13, 300 )
          integer, parameter :: FLT = selected_real_kind ( 6, 37 )
     
          ! Inputs
          integer, intent(in) :: iValue, myId, nNodeX, nNodeY, 
     .         nuPer, nuPar
          real(kind=dbl), intent(in) :: capR(:), zLoc(:), 
     .         uPerp(:), uPara(:)
          real(kind=dbl), intent(in) :: vc_mks
      
          ! Internal
          integer(KIND=MPI_OFFSET_KIND) :: pnc_nR, pnc_nz, 
     .         pnc_nuper, pnc_nupar, pnc_scalar
          character(len=100) :: pncFileName
          integer :: pnR_id, pnz_id, pnuper_id, pnupar_id, scalar_id
          integer :: dimIds(4)
          integer :: pf_R_id, pf_z_id, pf_pScale_id
          integer :: pf_vPer_id, pf_vPar_id, pf_vc_mks_id
          integer(KIND=MPI_OFFSET_KIND) :: pnc_start(4), pnc_cnt(4)
          integer :: pf_dfduper_id, pf_dfdupar_id
          integer :: pnc_f_id
     
          if ( myId == 0 ) write (*,*) 'WRITING output/p_f.nc ...'
          pncFileName = 'output/p_f.nc'
      
          pnc_nR  = nnodex
          pnc_nz  = nnodey
          pnc_nuper   = nuper
          pnc_nupar   = nupar
          pnc_scalar  = 1
          
#ifdef DLGHP 
          call p_check ( nfmpi_create ( 
     .     ivalue, pncFileName, NF_CLOBBER, MPI_INFO_NULL, pnc_f_id ) )
#else
          call p_check ( nfmpi_create ( 
     .     ivalue, pncFileName, ior ( NF_CLOBBER, NF_64BIT_OFFSET ),
     .     MPI_INFO_NULL, pnc_f_id ) )
#endif 
           if ( myId == 0 ) write (*,*) 'Defining variables ...'
     
          call p_check ( nfmpi_def_dim ( 
     . pnc_f_id, "nR", pnc_nR, pnR_id ) )
          call p_check ( nfmpi_def_dim ( 
     . pnc_f_id, "nz", pnc_nz, pnz_id ) )
          call p_check ( nfmpi_def_dim ( 
     . pnc_f_id, "nuper", pnc_nuper, pnuper_id ) )
          call p_check ( nfmpi_def_dim ( 
     . pnc_f_id, "nupar", pnc_nupar, pnupar_id ) )
          call p_check ( nfmpi_def_dim (
     . pnc_f_id, "scalar", pnc_scalar, scalar_id ) )
          
          dimIds  = (/ pnR_id, pnz_id, pnuper_id, pnupar_id /)
      
          call p_check ( nfmpi_def_var ( pnc_f_id, "dfduper", NF_REAL,
     . 4, dimIds, pf_dfduper_id ) )
          call p_check ( nfmpi_def_var ( pnc_f_id, "dfdupar", NF_REAL,
     . 4, dimIds, pf_dfdupar_id ) )
 
          call p_check ( nfmpi_def_var ( pnc_f_id, "R", NF_REAL,
     . 1, (/ pnR_id /), pf_R_id ) )
          call p_check ( nfmpi_def_var ( pnc_f_id, "z", NF_REAL,
     . 1, (/ pnz_id /), pf_z_id ) )
          call p_check ( nfmpi_def_var ( pnc_f_id, "vPer", NF_REAL,
     . 1, (/ pnuper_id /), pf_vPer_id ) )
          call p_check ( nfmpi_def_var ( pnc_f_id, "vPar", NF_REAL,
     . 1, (/ pnupar_id /), pf_vPar_id ) )
          call p_check ( nfmpi_def_var ( pnc_f_id, "vc_mks", NF_REAL,
     . 1, (/ scalar_id /), pf_vc_mks_id ) )
 
          call p_check ( nfmpi_def_var ( pnc_f_id, "pscale", NF_DOUBLE,
     . 1, (/ scalar_id /), pf_pscale_id ) )
      
          call p_check ( nfmpi_enddef ( pnc_f_id ) )
      
              !   Switch to independant data mode to
              !   write the R,z,vPer,vPar coords, i.e., 
              !   so they only get written into the file
              !   once.
              
          if ( myId == 0 ) write (*,*) 'Writing grid in indep mode ...'
     
          call p_check ( nfmpi_begin_indep_data ( pnc_f_id ) ) 
      
          if ( myId == 0 ) then 
     
              call p_check ( nfmpi_put_var_real (
     . pnc_f_id, pf_R_id, real(capR) ) )
              call p_check ( nfmpi_put_var_real ( 
     . pnc_f_id, pf_z_id, real(zLoc) ) )
              call p_check ( nfmpi_put_var_real ( pnc_f_id, 
     . pf_vPer_id, real(uperp * vc_mks) ) )
              call p_check ( nfmpi_put_var_real ( pnc_f_id, 
     . pf_vPar_id, real(upara * vc_mks) ) )
              call p_check ( nfmpi_put_var_real ( pnc_f_id, 
     . pf_vc_mks_id, real(vc_mks) ) )
      
          endif
      
          call p_check ( nfmpi_end_indep_data ( pnc_f_id ) ) 
          call p_check ( nfmpi_close ( pnc_f_id ) )
     
      end subroutine init_pf_dlg

      subroutine write_pf_dlg ( inI, inJ, iValue, 
     . nuPer, nuPar, dfduper, dfdupar )
      
          use dlg_p_check
          use dlg
 
          implicit none

          include 'mpif.h'
#         include <pnetcdf.inc>
      
          ! Parameters
          integer, parameter :: DBL = selected_real_kind ( 13, 300 )
          integer, parameter :: FLT = selected_real_kind ( 6, 37 )
     
          ! Inputs
          integer, intent(in) :: nuPer, nuPar
          integer, intent(in) :: inI, inJ, iValue
          real(kind=dbl), intent(in) :: dfduper(:,:), dfdupar(:,:)
      
          ! Internal

          integer(KIND=MPI_OFFSET_KIND) :: pnc_start(4), pnc_cnt(4)
          character(len=100) :: pncFileName
          integer :: pnc_f_id
          integer :: pf_dfduper_id, pf_dfdupar_id
         
          pncFileName = 'output/p_f.nc'
          write(*,*) 'opening' 

          call p_check ( nfmpi_open ( 
     .       ivalue, pncFileName, NF_WRITE, MPI_INFO_NULL, pnc_f_id ) )
          write(*,*) 'indep' 
          call p_check ( nfmpi_begin_indep_data ( pnc_f_id ) ) 
 
          call p_check ( nfmpi_inq_varid ( 
     .          pnc_f_id, 'dfduper', pf_dfduper_id ) )
           call p_check ( nfmpi_inq_varid ( 
     .          pnc_f_id, 'dfdupar', pf_dfdupar_id ) )
 
          pnc_start(1)    = inI 
          pnc_start(2)    = inJ 
          pnc_start(3)    = 1
          pnc_start(4)    = 1 
      
          pnc_cnt(1) = 1  
          pnc_cnt(2) = 1 
          pnc_cnt(3) = nuper
          pnc_cnt(4) = nupar 
        
         write(*,*) 'write 1' 

          call p_check ( nfmpi_put_vara_real ( 
     . pnc_f_id, pf_dfduper_id, pnc_start, pnc_cnt, real(dfduper) ) )
          write(*,*) 'write 2' 

          call p_check ( nfmpi_put_vara_real ( 
     . pnc_f_id, pf_dfdupar_id, pnc_start, pnc_cnt, real(dfdupar) ) )
            write(*,*) 'done' 

          call p_check ( nfmpi_end_indep_data ( pnc_f_id ) ) 
          call p_check ( nfmpi_close ( pnc_f_id ) )
          write(*,*) 'write_pf_dlg:180 ... closing'
      end subroutine write_pf_dlg

      subroutine write_pf_dlg2 ( inI, inJ, iValue, 
     . nuPer, nuPar, dfduper, dfdupar )
      
          use dlg
          use aorsa2din_mod, only: nModesX, nModesY
          use netcdf
          use nc_check

          implicit none

          include 'mpif.h'
      
          ! Parameters
          integer, parameter :: DBL = selected_real_kind ( 13, 300 )
          integer, parameter :: FLT = selected_real_kind ( 6, 37 )
     
          ! Inputs
          integer, intent(in) :: nuPer, nuPar
          integer, intent(in) :: inI, inJ, iValue
          real(kind=dbl), intent(in) :: dfduper(:,:), dfdupar(:,:)
      
          ! Internal

          integer :: pnc_start(4), pnc_cnt(4)
          character(len=100) :: pncFileName
          integer :: pnc_f_id
          integer :: pf_dfduper_id, pf_dfdupar_id
          integer :: nR, nz, nvper, nvpar, dim_ids(4)
         
          pncFileName = 'output/p_f.nc'

          call check ( nf90_open ( 
     .       pncFileName, NF90_WRITE, pnc_f_id ) )
 
          call check ( nf90_inq_varid ( 
     .          pnc_f_id, 'dfduper', pf_dfduper_id ) )
          call check ( nf90_inq_varid ( 
     .          pnc_f_id, 'dfdupar', pf_dfdupar_id ) )

          call check ( nf90_inquire_variable ( 
     .         pnc_f_id, pf_dfduper_id, dimIds = dim_ids ) )
 
          call check ( nf90_inquire_dimension ( pnc_f_id, dim_ids(1), 
     .         len = nR ) ) 
          call check ( nf90_inquire_dimension ( pnc_f_id, dim_ids(2), 
     .         len = nz ) ) 
          call check ( nf90_inquire_dimension ( pnc_f_id, dim_ids(3), 
     .         len = nvper ) ) 
          call check ( nf90_inquire_dimension ( pnc_f_id, dim_ids(4), 
     .         len = nvpar ) ) 


          pnc_start(1)    = inI 
          pnc_start(2)    = inJ 
          pnc_start(3)    = 1
          pnc_start(4)    = 1 
      
          pnc_cnt(1) = 1  
          pnc_cnt(2) = 1 
          pnc_cnt(3) = nuper
          pnc_cnt(4) = nupar 
        
          call check ( nf90_put_var ( 
     . pnc_f_id, pf_dfduper_id, real(dfduper), pnc_start, pnc_cnt))

          call check ( nf90_put_var ( 
     . pnc_f_id, pf_dfdupar_id, real(dfdupar), pnc_start, pnc_cnt))

          call check ( nf90_close ( pnc_f_id ) )

      end subroutine write_pf_dlg2


      end module write_pf
