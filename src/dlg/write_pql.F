!DLG:   I hate hate hate coding in F77! However
!       here I have no choice since pNetCDF does 
!       not have F90 bindings yet. Serious!? And
!       the whole purpose of this seperate subroutine
!       is just so I can compile it without the freaking
!       -r8 type flag that the rest of AORSA uses. That
!       is just very sloppy programming, comeon guys.

      module dlg_p_check
          
      contains
          subroutine p_check ( status )
#             include <pnetcdf.inc>
              integer, intent ( in) :: status
      
              if(status /= nf_noerr) then 
                  print *, trim(nfmpi_strerror(status))
                  stop "p_check: Stopped"
              end if
      
          end subroutine p_check 
      end module dlg_p_check
 
      module write_pql
      
      contains
      subroutine write_pql_dlg ( iValue, myId,
     . nNodeX, nNodeY, nuPer, nuPar, 
     . vc_mks, capR, zLoc, uPerp, uPara,
     . npRowOut, npColout, myRow, myCol, 
     . bql4d_cyl, eNorm )
      
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
          real(kind=dbl), intent(in) :: vc_mks, eNorm
          integer, intent(in) :: npRowOut, npColOut, myRow, myCol
          real(kind=dbl), intent(in) :: bql4D_cyl(:,:,:,:)
      
          ! Internal
          integer(KIND=MPI_OFFSET_KIND) :: pnc_nR, pnc_nz, 
     .         pnc_nuper, pnc_nupar, pnc_scalar
          integer :: pnc_id
          character(len=100) :: pncFileName
          integer :: pnR_id, pnz_id, pnuper_id, pnupar_id, scalar_id
          integer :: dimIds(4)
          integer :: pql_R_id, pql_z_id, pql_pScale_id
          integer :: pql_vPer_id, pql_vPar_id, pql_vc_mks_id,
     .       pql_eNorm_id
          integer(KIND=MPI_OFFSET_KIND) :: pnc_start(4), pnc_cnt(4)
          integer :: pql_b_id, pql_c_id, pql_e_id, pql_f_id
      
          if ( myId == 0 ) write (*,*) 'WRITING output/p_ql.nc ...'
          pncFileName = 'output/p_ql.nc'
      
          pnc_nR  = nnodex
          pnc_nz  = nnodey
          pnc_nuper   = nuper
          pnc_nupar   = nupar
          pnc_scalar  = 1
          
#ifdef DLGHP 
          call p_check ( nfmpi_create ( 
     .     ivalue, pncFileName, NF_CLOBBER, MPI_INFO_NULL, pnc_id ) )
#else
          call p_check ( nfmpi_create ( 
     .     ivalue, pncFileName, ior ( NF_CLOBBER, NF_64BIT_OFFSET ),
     .     MPI_INFO_NULL, pnc_id ) )
#endif 
           if ( myId == 0 ) write (*,*) 'Defining variables ...'
     
          call p_check ( nfmpi_def_dim ( 
     . pnc_id, "nR", pnc_nR, pnR_id ) )
          call p_check ( nfmpi_def_dim ( 
     . pnc_id, "nz", pnc_nz, pnz_id ) )
          call p_check ( nfmpi_def_dim ( 
     . pnc_id, "nuper", pnc_nuper, pnuper_id ) )
          call p_check ( nfmpi_def_dim ( 
     . pnc_id, "nupar", pnc_nupar, pnupar_id ) )
          call p_check ( nfmpi_def_dim (
     . pnc_id, "scalar", pnc_scalar, scalar_id ) )
          
          dimIds  = (/ pnR_id, pnz_id, pnuper_id, pnupar_id /)
      
          call p_check ( nfmpi_def_var ( pnc_id, "ql_b", NF_REAL,
     . 4, dimIds, pql_b_id ) )
          !call p_check ( nfmpi_def_var ( pnc_id, "ql_c", NF_DOUBLE, 4, dimIds, pql_c_id ) )
          !call p_check ( nfmpi_def_var ( pnc_id, "ql_e", NF_DOUBLE, 4, dimIds, pql_e_id ) )
          !call p_check ( nfmpi_def_var ( pnc_id, "ql_f", NF_DOUBLE, 4, dimIds, pql_f_id ) )
          call p_check ( nfmpi_def_var ( pnc_id, "R", NF_REAL,
     . 1, (/ pnR_id /), pql_R_id ) )
          call p_check ( nfmpi_def_var ( pnc_id, "z", NF_REAL,
     . 1, (/ pnz_id /), pql_z_id ) )
          call p_check ( nfmpi_def_var ( pnc_id, "vPer", NF_REAL,
     . 1, (/ pnuper_id /), pql_vPer_id ) )
          call p_check ( nfmpi_def_var ( pnc_id, "vPar", NF_REAL,
     . 1, (/ pnupar_id /), pql_vPar_id ) )
          call p_check ( nfmpi_def_var ( pnc_id, "vc_mks", NF_REAL,
     . 1, (/ scalar_id /), pql_vc_mks_id ) )
          call p_check ( nfmpi_def_var ( pnc_id, "eNorm", NF_REAL,
     . 1, (/ scalar_id /), pql_eNorm_id ) )
 
          call p_check ( nfmpi_def_var ( pnc_id, "pscale", NF_DOUBLE,
     . 1, (/ scalar_id /), pql_pscale_id ) )
      
          call p_check ( nfmpi_enddef ( pnc_id ) )
      
              !   Switch to independant data mode to
              !   write the R,z,vPer,vPar coords, i.e., 
              !   so they only get written into the file
              !   once.
              
          if ( myId == 0 ) write (*,*) 'Writing grid in indep mode ...'
     
          call p_check ( nfmpi_begin_indep_data ( pnc_id ) ) 
      
          if ( myId == 0 ) then 
     
                !write(*,*) 'write_pql: ', vc_mks
                !write(*,'(65(2x,f5.2))') uperp

              call p_check ( nfmpi_put_var_real (
     . pnc_id, pql_R_id, real(capR) ) )
              call p_check ( nfmpi_put_var_real ( 
     . pnc_id, pql_z_id, real(zLoc) ) )
              call p_check ( nfmpi_put_var_real ( pnc_id, 
     . pql_vPer_id, real(uperp * vc_mks) ) )
              call p_check ( nfmpi_put_var_real ( pnc_id, 
     . pql_vPar_id, real(upara * vc_mks) ) )
              call p_check ( nfmpi_put_var_real ( pnc_id, 
     . pql_vc_mks_id, real(vc_mks) ) )
               call p_check ( nfmpi_put_var_real ( pnc_id, 
     . pql_eNorm_id, real(eNorm) ) )
      
          endif
      
          call p_check ( nfmpi_end_indep_data ( pnc_id ) ) 
      
          !   Switched back to collective data mode
      
          pnc_start(1)    = nModesX/npRowOut * myRow + 1
          pnc_start(2)    = nModesY/npColOut * myCol + 1 
          pnc_start(3)    = 1
          pnc_start(4)    = 1 
      
          pnc_cnt(1) = nModesX/npRowOut  
          pnc_cnt(2) = nModesY/npColOut
          pnc_cnt(3) = nuper
          pnc_cnt(4) = nupar 
     
          if ( myId == 0 ) write (*,*) 'Writing QL values ...'

          call p_check ( nfmpi_put_vara_real_all ( 
     . pnc_id, pql_b_id, pnc_start, pnc_cnt, real(bql4d_cyl) ) )
          !call p_check ( nfmpi_put_vara_double_all ( pnc_id, pql_c_id, pnc_start, pnc_cnt, real(cql4D_cyl) ) )
          !call p_check ( nfmpi_put_vara_double_all ( pnc_id, pql_e_id, pnc_start, pnc_cnt, real(eql4D_cyl) ) )
          !call p_check ( nfmpi_put_vara_double_all ( pnc_id, pql_f_id, pnc_start, pnc_cnt, real(fql4D_cyl) ) )
            
          call p_check ( nfmpi_close ( pnc_id ) )
      
      end subroutine write_pql_dlg
      
      end module write_pql
